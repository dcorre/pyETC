#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from scipy.integrate import dblquad
from scipy.interpolate import interp1d
from skimage.draw import circle
from astropy.io import ascii

from . import photometry as phot


def pixelsize(info_dict):
    """
    Computes the pixel size in m/pixel
    Parameters
    -----------
    info_dict: dictionnary

    Returns
    -------
    pixelsize: float
        size of the pixel in meter/pixel
    """
    # Pixel length (m/pix)
    pixsize = (
        info_dict["cameras"][info_dict["channel"]]["Photocell_SizeX"]
        * info_dict["binning_X"]
    )
    info_dict["pixelSize"] = pixsize
    return info_dict


def pixel_scale(info_dict):
    """
    Returns the pixel scale in arcsec/pixel. Depends on the chosen binning
    Parameters
    -----------
    info_dict: dictionnary

    Return
    -------
    pixel_scale: float
        pixel scale in arcsec/pixel
    """

    # Pixel spatial sampling on axis (arcsec/pix)
    # Field of view is defined is arcmin so *60
    # pixscale = info_dict['FoV_1axis'] /
    # (info_dict['cameras'][info_dict['channel']]['Nphotocell_X']
    # / info_dict['binning_X'])*60
    pixscale_X = (
        info_dict["FoV_axis1"]
        / (
            info_dict["cameras"][info_dict["channel"]]["Nphotocell_X"]
            / info_dict["binning_X"]
        )
        * 60
    )
    pixscale_Y = (
        info_dict["FoV_axis2"]
        / (
            info_dict["cameras"][info_dict["channel"]]["Nphotocell_Y"]
            / info_dict["binning_Y"]
        )
        * 60
    )
    info_dict["pixelScale_X"] = pixscale_X
    info_dict["pixelScale_Y"] = pixscale_Y
    return info_dict


def Focal_length_eff(info_dict):
    """
    Computes the effective focal length give the plate scale and
    the pixel size.
    Parameters
    ----------
    info_dict: dictionary

    Returns
    ---------
    F_eff: float
        effective focal length in m
    """
    # F_eff = 1./ (info_dict['pixelScale'] *
    #              np.pi/(3600*180) / info_dict['pixelSize'])

    if isinstance(info_dict["focal_length"], dict):
        info_dict["foc_len"] = info_dict["focal_length"][info_dict["channel"]]
    else:
        info_dict["foc_len"] = info_dict["focal_length"]
    return info_dict


def field_of_view_tel(info_dict, unit="arcmin"):
    """ Returns the field of view of the telescope
    Parameters
    ----------
    info_dict: dictionary

    unit: string
        unit of the field of view. 'deg' or 'arcmin'
    Returns
    ---------
    FoV: float
        field of view of the in deg2
    """

    fov_axis1 = (
        2
        * np.arctan(
            info_dict["cameras"][info_dict["channel"]]["Nphotocell_X"]
            * info_dict["cameras"][info_dict["channel"]]["Photocell_SizeX"]
            / info_dict["foc_len"]
            / 2.0
        )
        * 180.0
        / np.pi
    )
    fov_axis2 = (
        2
        * np.arctan(
            info_dict["cameras"][info_dict["channel"]]["Nphotocell_Y"]
            * info_dict["cameras"][info_dict["channel"]]["Photocell_SizeY"]
            / info_dict["foc_len"]
            / 2.0
        )
        * 180.0
        / np.pi
    )

    if unit == "arcmin":
        fov_axis1 *= 60
        fov_axis2 *= 60
        fov_tel = fov_axis1 * fov_axis2

    info_dict["FoV_tel"] = fov_tel
    info_dict["FoV_axis1"] = fov_axis1
    info_dict["FoV_axis2"] = fov_axis2
    return info_dict


def obstruction(info_dict):
    """ Obstruction due to M2 """
    # Multiply by 1.05 to take the arms into account
    info_dict["obstruction"] = (
        info_dict["M2_factor"] * info_dict["D_M2"]
    ) ** 2.0 / info_dict["D_M1"] ** 2.0

    return info_dict


def A_tel(info_dict):
    """
    Compute the aperture area of the telescope taking into
    account the obscuration.

    Parameter
    ---------
    info_dict: dictionnary

    Returns
    --------
    A_tel: float
        aperture of the telescope taking into account the obscuration (cm2)
    """
    exposure_area = (
        1e4
        * np.pi
        * ((info_dict["D_M1"]) / 2.0) ** 2.0
        * (1 - info_dict["obstruction"])
    )  # cm2
    info_dict["A_tel"] = exposure_area
    return info_dict


def A_pixel(info_dict):
    """
    Compute the projected area of one pixel on the sky

    Parameter
    ---------
    info_dict: dictionnary

    Returns
    --------
    A_pixel: float
        projected area of one pixel on the sky (arcsec2/pix)
    """
    pix_sky_area = info_dict["pixelScale_X"] * info_dict["pixelScale_Y"]
    info_dict["A_pixel"] = pix_sky_area
    return info_dict


def instrument_background(info_dict):
    """
    Compute the instrument background emission in function
    of the night temperature.
    Parameters
    -----------
    info_dict: dictionary

    Returns:
    ---------
    inst_bg: float
        in electrons/sec/pixel
    """

    try:
        data = ascii.read(
            "%s/local_conditions/instrument_background/%s.dat"
            % (info_dict["MainDirectory"], info_dict["telescope"]),
            "r",
        )
        if info_dict["filter_band"][0] in data.colnames:
            Temp_C = data["Temp"]
            inst_noise = data[info_dict["filter_band"][0]]
            f = interp1d(Temp_C, inst_noise, kind="linear")

            # e-/s/photocells
            inst_bg = f(info_dict["Temp_night"])
            # e-/s/pixel
            inst_bg *= info_dict["binning_X"] * info_dict["binning_Y"]
        else:
            inst_bg = 0  # e-/s/pixel
    except:
        # print ('No file found for computing the instrument'
        #        'background of %s. Set to 0' % info_dict['telescope'])
        inst_bg = 0  # e-/s/pixel

    info_dict["Instrument_bg"] = inst_bg
    return info_dict


def system_response(info_dict):
    """
    Computes the spectral system efficiency of the telescope
    Atmosphere + optics + filters + camera

    Parameters
    ----------
    info_dict: dictionary

    wavelength: array
        wavelength in angstrom
    Returns
    -------
    sys_eff: array
        spectral system efficiency of the telescope, range [0,1]
    """
    if info_dict["detailed_trans"] == 1:
        # Transmission of the atm / optics / filter / camera
        # Trans_atms = local.atmospheric_transmission(info_dict)[1]
        filter_trans = phot.set_filter(info_dict)
        info_dict["Trans_filter"] = filter_trans
        sys_eff = (
            info_dict["Trans_optics_total"]
            * filter_trans
            * info_dict["camera_efficiency"]
        )
    elif info_dict["detailed_trans"] == 0:
        sys_eff = phot.set_filter(info_dict) * info_dict["camera_efficiency"]
    info_dict["system_response"] = sys_eff
    return info_dict


def FWHM(info_dict):
    """ Computes the Full Width at Half Maximum of the PSF

    Parameters
    ----------
    info_dict: dictionary

    wavelength: array
        wavelengths in angstrom

    Returns
    -------
    FWHM: float
        full width at half maximum for each wavelentgh in arcsec

    """
    try:
        fwhm_inst = info_dict["Fwhm_psf_opt"][info_dict["filter_band"]]
    except:
        fwhm_inst = info_dict["Fwhm_psf_opt"]
    # print (band_name,fwhm_inst)
    FWHM_type = info_dict["fwhm_type"]
    if FWHM_type == "seeing":
        # seeing
        # seeing_arcsec = seeing(info_dict)
        # Add instrumental psf
        fwhm_tot = np.sqrt(info_dict["seeing_los_arcsec"]**2.0
                           + fwhm_inst**2.0)

    elif FWHM_type == "diffraction":
        # Diameter of the first mirror in (cm)
        D1 = info_dict["D_M1"] * 1e2
        fwhm_diffraction = (1.22 * info_dict["wavelength_ang"]
                            / D1 * 1e-8 * 206265)
        # Add instrumental psf
        fwhm_tot = np.sqrt(fwhm_diffraction**2.0 + fwhm_inst**2.0)

    elif FWHM_type == "seeing_diff":

        # Diameter of the first mirror in (cm)
        D1 = info_dict["D_M1"] * 1e2
        fwhm_diffraction = (1.22 * info_dict["wavelength_ang"]
                            / D1 * 1e-8 * 206265)

        fwhm_tot = np.sqrt(
            info_dict["seeing_los_arcsec"]**2.0 + fwhm_diffraction**2.0
        )

    elif FWHM_type == "prf":
        # intrinsinc instrumental psf fwhm_int given as an input
        # pixel response function estimated as a gaussian with
        # the plate scale as the pixel width
        fwhm_prf = info_dict["pixelScale_X"]
        fwhm_tot = np.sqrt(
            info_dict["seeing_los_arcsec"]**2. + fwhm_inst**2. + fwhm_prf**2.
        )

    info_dict["FWHM_tot"] = fwhm_tot
    return info_dict


def radius_aperture(info_dict):
    """ Compute the radius of the aperture

    Parameters
    ----------
    info_dict: dictionary

    wavelength: float
        wavelength in angstrom

    Returns
    --------
    aperture: float
        radius of the aperture in arcsec
    """
    if info_dict["psf_type"] == "gaussian":
        aperture = info_dict["radius_int_fwhm"] * info_dict["FWHM_tot"]
    elif info_dict["psf_type"] == "moffat":
        aperture = info_dict["radius_int_fwhm"] * info_dict["FWHM_tot"]
    info_dict["radius_aperture_phot"] = aperture
    return info_dict


def PSF_2D_gaussian(sigx, sigy, x, y):
    """ Compute the PSF of the instrument with a 2D gaussian

    Parameters
    -----------
    sigx: float
        standard deviation along the x axis
    sigy: float
        standard deviation along the y axis
    x: float
        position along the x axis
    y: float
        position along the y axis

    wavelength: float
        effective wavelength of the filter in angstrom
    Returns:
    ---------
    psf: array
        psf of the instrument
    """

    psf = (
        1.0
        / (2.0 * np.pi * sigx * sigy)
        * np.exp(-1.0 / 2 * (x ** 2.0 / sigx ** 2.0 + y ** 2.0 / sigy ** 2.0))
    )
    return psf


def PSF_Moffat(alpha, beta, x, y):
    """ Compute the PSF of the instrument with a Moffat function

    Parameters
    -----------
    alpha: float
        radial parameter
    beta: float
        power indice of the function
    x: float
        position along the x axis
    y: float
        position along the y axis

    wavelength: float
        effective wavelength of the filter in angstrom

    Returns:
    ---------
    psf: array
        psf of the instrument
    """

    psf = (
        (beta - 1.0)
        / (np.pi * alpha * alpha)
        * (1.0 + (x * x + y * y) / (alpha * alpha)) ** (-beta)
    )
    return psf


def Normalisation_factor(info_dict, central_pix):
    """
    Computes the normaliation factor.
    It corresponds to the fraction of the flux we want to keep.
    Number in the range[0,1]
    Parameters
    ----------
    info_dict: dictionnary

    central_pixel: boolean
        if True --> normalisation factor computed for the central pixel
        if False --> normalisation factor computed for the total flux

    wavelength: array
                wavelength in angstrom
    Returns
    -------
    norm_factor: float
                 normalisation factor, number in the range [0,1]
    """
    if info_dict["psf_type"] == "gaussian":
        # computation of the gaussian fraction covered by the brightest pixel
        if central_pix:
            # Max flux at the center of the pixel
            pX = info_dict["pixelScale_X"]  # arcsec
            pY = info_dict["pixelScale_Y"]  # arcsec
            a = pX / 2.0
            b = pY / 2.0
            sigma = info_dict["FWHM_tot"] / (2.0 * np.sqrt(2.0 * np.log(2.0)))

            norm_factor = dblquad(
                lambda x, y: PSF_2D_gaussian(sigma, sigma, x, y),
                -b,
                b,
                lambda x: -a * np.sqrt(1 - x ** 2.0 / b ** 2.0),
                lambda x: a * np.sqrt(1 - x ** 2.0 / b ** 2.0),
            )[0]

        else:  # Flux fraction in the area of radius npix

            # semi X axis
            a = info_dict["radius_aperture_phot"]  # arcsec
            siga = info_dict["FWHM_tot"] / (
                2.0 * np.sqrt(2.0 * np.log(2.0))
            )  # (arcsec)
            # semi Y axis
            b = a
            sigb = siga  # arcsec

            norm_factor = dblquad(
                lambda x, y: PSF_2D_gaussian(siga, sigb, x, y),
                -b,
                b,
                lambda x: -a * np.sqrt(1 - x ** 2.0 / b ** 2.0),
                lambda x: a * np.sqrt(1 - x ** 2.0 / b ** 2.0),
            )[0]

    elif info_dict["psf_type"] == "moffat":
        # computation of the gaussian fraction covered by the brightest pixel
        if central_pix:
            # Max flux at the center of the pixel
            pX = info_dict["pixelScale_X"]  # arcsec
            pY = info_dict["pixelScale_Y"]  # arcsec
            a = pX / 2.0
            b = pY / 2.0
            beta = info_dict["moffat_beta"]
            alpha = info_dict["FWHM_tot"] / (2.0 * np.sqrt(2.0 ** (1.0 / beta)
                                                           - 1.0))

            norm_factor = dblquad(
                lambda x, y: PSF_Moffat(alpha, beta, x, y),
                -b,
                b,
                lambda x: -a * np.sqrt(1 - x ** 2.0 / b ** 2.0),
                lambda x: a * np.sqrt(1 - x ** 2.0 / b ** 2.0),
            )[0]

        else:  # Flux fraction in the area of radius npix

            # semi X axis
            a = info_dict["radius_aperture_phot"]  # arcsec
            b = a
            beta = info_dict["moffat_beta"]
            alpha = info_dict["FWHM_tot"] / (2.0 * np.sqrt(2.0 ** (1.0 / beta)
                                                           - 1.0))

            norm_factor = dblquad(
                lambda x, y: PSF_Moffat(alpha, beta, x, y),
                -b,
                b,
                lambda x: -a * np.sqrt(1 - x ** 2.0 / b ** 2.0),
                lambda x: a * np.sqrt(1 - x ** 2.0 / b ** 2.0),
            )[0]

    return norm_factor


def spatial_binning(info_dict):
    """
    Compute the spatial binning representing the number of pixels
    affected by the light in the spatial direction.

    Parameters
    ----------
    info_dict: dictionnary

    wavelength: float
        wavelength in angstrom
    Returns
    ----------
    spabin: float
        spatial binning in arcsec
    """
    # x2 to get the diameter, /pixel_scale to express it in pixels
    spabin = (2.0 * info_dict["radius_aperture_phot"]
              / info_dict["pixelScale_X"])
    info_dict["spatial_binning"] = spabin
    return info_dict


def nb_pixels(info_dict):
    """ Compute the number of pixels affected by the target light
    Parameters
    ----------
    info_dict: dictionnary

    wavelength: float
        wavelength in angstrom

    Returns
    --------
    npix: float
        number of pixels affected by the target light (ceiled)
    """

    # Radius of the source in pixel
    R_source = (
        info_dict["spatial_binning"] / 2.0
    )
    # R_source = np.ceil(np.pi * (phot_dict['radius_mult_fwhm']
    #                             * Fwhm_psf/2.)**2. / (pixsize1*pixsize2))

    # side length of the image is a square
    # side_length = 2.0 * R_source

    # Computed as in the LAM ETC
    # npix = np.ceil(side_length * side_length)

    # Standard computation
    # npix = np.ceil(np.pi * R_source**2.)

    if info_dict["PSF_position_on_ccd"] == "nogrid":
        npix = np.ceil(
            np.pi
            * info_dict["radius_aperture_phot"] ** 2
            / (info_dict["pixelScale_X"] * info_dict["pixelScale_Y"])
        )

    elif info_dict["PSF_position_on_ccd"] == "center":
        img1 = np.zeros((100, 100), dtype=np.uint8)
        rr1, cc1 = circle(50, 50, R_source + 0.5)
        img1[rr1, cc1] = 1
        npix = np.sum(img1)

    elif info_dict["PSF_position_on_ccd"] == "corner":
        img2 = np.zeros((100, 100), dtype=np.uint8)
        rr2, cc2 = circle(50.5, 50.5, R_source + 0.5)
        img2[rr2, cc2] = 1
        npix = np.sum(img2)

    # print (R_source,npix,npix2,npix3)

    # if info_dict['PSF_position_on_ccd'] == 'nogrid':
    #      info_dict['npix']=npix
    # elif info_dict['PSF_position_on_ccd'] == 'center':
    #      info_dict['npix']=npix2
    # elif info_dict['PSF_position_on_ccd'] == 'corner':
    #      info_dict['npix']=npix3
    # elif info_dict['PSF_position_on_ccd'] == 'mean':
    #      info_dict['npix']=(npix2+npix3)/2
    # print (info_dict['PSF_position_on_ccd'], info_dict['npix'])
    info_dict["npix"] = npix
    return info_dict


def factor_images_averaged(info_dict):
    """
    Computes the factor mutliplying the variance noise in order to
    obtain an unbiaised variance in the case of estimating the
    background noise from several images.
    Parameters
    -----------
    info_dict: dictionary

    Returns
    -------
    factor_averaged: float
        mutltiplying factor for getting an unbiaised variance when
        estimating it with several images

    """
    if info_dict["Nexp"] == 1:
        factor_averaged = 1.0
    else:
        # factor_averaged = 1. + 1./(info_dict['Nexp']-1)
        factor_averaged = 1.0

    return factor_averaged


def t_dithering(info_dict):
    """ Compute the time wasted due to dithering

    Parameters
    ----------
    info_dict: dictionary

    Returns
    -------
    t_dithering: float
        time wasted due to dithering in seconds
    """
    t_dithering = 0.0
    info_dict["T_dithering"] = t_dithering
    return info_dict


def dead_time(info_dict):
    """ Compute the dead time due to the readout and the dithering

    Parameters
    ----------
    info_dict: dictionary

    Returns
    -------
    t_dead: float
        dead time in seconds
    """
    # t_dead = max(t_dithering(), info_dict['readouttime'])
    t_dead = info_dict["T_dithering"]
    info_dict["deadtime_tot"] = t_dead
    return info_dict


def set_preliminary_computation(info_dict):
    """ Set the parameters """
    info_dict = pixelsize(info_dict)
    info_dict = Focal_length_eff(info_dict)
    info_dict = field_of_view_tel(info_dict, unit="arcmin")
    info_dict = pixel_scale(info_dict)
    info_dict = obstruction(info_dict)
    info_dict = A_tel(info_dict)
    info_dict = A_pixel(info_dict)
    info_dict = instrument_background(info_dict)
    info_dict = system_response(info_dict)
    info_dict = FWHM(info_dict)
    info_dict = radius_aperture(info_dict)
    info_dict = spatial_binning(info_dict)
    info_dict = nb_pixels(info_dict)
    info_dict = t_dithering(info_dict)
    info_dict = dead_time(info_dict)

    return info_dict

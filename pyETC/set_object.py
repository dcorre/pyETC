#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from scipy.interpolate import interp1d
from astropy.io import ascii
from . import utils
from . import constants as cc
from . import photometry as phot


def set_object(info_dict):
    """Compute the number of electrons coming from the object per second

    Parameters
    ----------
    info_dict: dictionary

    wavelength : array
        wavelengths in angstrom

    Returns
    ---------
    F_e_s : float
        countrate of the object in e-/s

    """
    if info_dict["object_type"] == "magnitude":
        mag = info_dict["object_magnitude"]
        F_e_s = 10 ** (0.4 * (info_dict["zeropoint"] - mag))  # e/s
        fmag = mag * np.ones(len(info_dict["wavelength_ang"]))
        fJy = phot.mag2Jy(info_dict, fmag)  # Jy
        # erg/s/cm2/A
        flam = utils.fJy_to_flambda(info_dict["wavelength_ang"], fJy)
        # ph/s/cm2/A
        fph = utils.flambda_to_fph(info_dict["wavelength_ang"], flam)

    elif info_dict["object_type"] == "spectrum":

        object_path = (
            info_dict["path"] + info_dict["object_folder"] + info_dict["object_file"]
        )

        inFile = open(object_path, "r")
        lines = inFile.readlines()

        x = []  # wavelength
        y = []  # flux in erg/s/cm2/A
        for line in lines:
            if line[0] != "#" and len(line) > 3:
                bits = line.split()
                y.append(float(bits[1]))
                x.append(float(bits[0]))

        x = np.array(x)
        y = np.array(y, dtype=np.float64)

        f = interp1d(x, y, kind="linear")
        if min(info_dict["wavelength_ang"]) < min(x) or max(
            info_dict["wavelength_ang"]
        ) > max(x):
            print(
                "The wavelength coverage must be smaller or equal to the "
                "one of the input spectrum. Please adapt it in the "
                "configuration file.\nSpectrum wavelength coverage: "
                "%.2f-%.2f Angstroms\n" % (min(x), max(x)),
                "Current chosen wavelength coverage: "
                "%.2f-%.2f Angstroms"
                % (min(info_dict["wavelength_ang"]), max(info_dict["wavelength_ang"])),
            )
        # erg/s/cm2/A
        flam = f(info_dict["wavelength_ang"])
        # ph/s/cm2/A
        fph = utils.flambda_to_fph(info_dict["wavelength_ang"], flam)
        fmag = phot.Jy2Mag(
            info_dict, utils.flambda_to_fJy(info_dict["wavelength_ang"], flam)
        )
        # e/s
        F_e_s = (
            np.trapz(
                fph * info_dict["system_response"] * info_dict["Trans_atmosphere"],
                info_dict["wavelength_ang"],
            )
            * info_dict["A_tel"]
        )

    elif info_dict["object_type"] == "grb_sim":
        # -----------------------------
        # Compute the emission spectra
        # -----------------------------

        if info_dict["grb_model"] != "LightCurve":
            if info_dict["grb_model"] == "gs02":
                try:
                    from pyGRBaglow.synchrotron_model import fireball_afterglow as grb
                except ValueError:
                    print("Package pyGRBaglow not found." "Need to be installed")

                td = info_dict["t_sinceBurst"]  # in days
                DIT = info_dict["exptime"]  # in second
                afterglow = grb(
                    n0=info_dict["n0"],
                    eps_b=info_dict["eps_b"],
                    eps_e=info_dict["eps_e"],
                    E_iso=info_dict["E_iso"],
                    eta=info_dict["eta"],
                    p=info_dict["p"],
                    Y=info_dict["Y"],
                    z=info_dict["grb_redshift"],
                    ism_type=info_dict["ism_type"],
                    disp=0,
                )
                time_grb = np.linspace(
                    td, td + DIT / 86400, 5
                )  # divides exposure time in 5
                frequencies = cc.c_light_m_s / (info_dict["wavelength_ang"] * 1e-10)
                # in mJy
                afterglow_lc = afterglow.light_curve(time_grb, frequencies)
                # afterglow_lc2=afterglow.light_curve(td+DIT/2/86400,frequencies)
                # factor to convert in Jy
                factor_Jy = 1e-3
                factor_time = 86400

            elif info_dict["grb_model"] == "SPL":
                try:
                    from pyGRBaglow.template_models import Templates as grb
                except ValueError:
                    print("Package grb_afterglow not found." "Need to be installed")

                td = float(info_dict["t_sinceBurst"]) * 86400  # in second
                DIT = info_dict["exptime"]  # in second

                time_grb = np.linspace(td, td + DIT, 5)
                afterglow = grb(
                    F0=info_dict["F0"], t0=info_dict["t0"], wvl0=info_dict["wvl0"]
                )
                afterglow_lc = afterglow.light_curve(
                    info_dict["wavelength_ang"],
                    time_grb,
                    [info_dict["alpha"], info_dict["beta"]],
                    model="SPL",
                )  # in Jy
                # conversion factor set to 1 (already in Jy)
                factor_Jy = 1
                factor_time = 1

            elif info_dict["grb_model"] == "BPL":
                try:
                    from pyGRBaglow.template_models import Templates as grb
                except ValueError:
                    print("Package grb_afterglow not found." "Need to be installed")

                td = float(info_dict["t_sinceBurst"]) * 86400  # in second
                DIT = info_dict["exptime"]  # in second

                time_grb = np.linspace(td, td + DIT, 5)

                afterglow = grb(
                    F0=info_dict["F0"], t0=info_dict["t0"], wvl0=info_dict["wvl0"]
                )
                afterglow_lc = afterglow.light_curve(
                    info_dict["wavelength_ang"],
                    time_grb,
                    [
                        info_dict["alpha1"],
                        info_dict["alpha2"],
                        info_dict["beta"],
                        info_dict["s"],
                    ],
                    model="BPL",
                )  # in Jy
                # conversion factor set to 1 (already in Jy)
                factor_Jy = 1
                factor_time = 1

            # Sum over the time
            # dtime=np.diff(time_series)
            sed_stacked = np.zeros(len(info_dict["wavelength_ang"]))
            # sed_stacked2=np.zeros(len(info_dict['wavelength_ang']))
            for i in range(len(info_dict["wavelength_ang"])):
                sed_stacked[i] = np.trapz(afterglow_lc[:, i], time_grb)
                # sed_stacked2[i]=afterglow_lc2[:,i]

            # free memory
            afterglow_lc = None

            # Convert to Jy
            grb_fJy = (
                sed_stacked * factor_Jy / (DIT / factor_time)
            )  # GRB SED in Jy in the observer frame
            """
            # /(DIT/86400.) to go from Jy.day --> Jy because of the
            integration in time due to the time dependancy.
            It has to be divided by the exposure time to recover an e-/s
            unit for the SNR formula.
            It is basically the mean value in Jy in the time interval
            corresponding to the exposure time. This assumption is valid for
            rather short exposures and probably not for long exposures,
            it has to be tested.
            """
            # grb_fJy2 = sed_stacked2*1e-3

            # print (grb_fJy,grb_fJy2)

            # Apply Host galaxy and IGM extinction
            try:
                from pyGRBaglow.reddening_cy import Pei92
            except:
                from pyGRBaglow.reddening import Pei92

            if info_dict["IGM_extinction_model"] == "meiksin":
                try:
                    from pyGRBaglow.igm_cy import meiksin
                except:
                    from pyGRBaglow.igm import meiksin

                grb_fJy *= meiksin(
                    info_dict["wavelength_ang"] / 10.0, info_dict["grb_redshift"]
                )
            elif info_dict["IGM_extinction_model"] == "madau":
                from pyGRBaglow.igm import madau
                grb_fJy *= madau(info_dict["wavelength_ang"], info_dict["grb_redshift"])

            if info_dict["host_extinction_law"] in ["mw", "lmc", "smc"]:
                grb_fJy *= Pei92(
                    info_dict["wavelength_ang"],
                    info_dict["Av_Host"],
                    info_dict["grb_redshift"],
                    ext_law=info_dict["host_extinction_law"],
                    Xcut=True
                )[1]

            if info_dict["galactic_extinction_law"].lower() != "none":
                grb_fJy *= Pei92(
                    info_dict["wavelength_ang"],
                    info_dict["Av_galactic"],
                    0.,
                    ext_law="mw",
                    Xcut=True
                )[1]

            """
            Integration over the exposure time for each wavelength because
            the GRB spectrum varies with the time of observation
            int_grb_sed = []
            def integrand(x,f):
                return afterglow.light_curve(x,f)
            # Flux in mJy.day
            for index in range(len(info_dict['wavelength_ang'])):
                int_grb_sed.append(quad(integrand,
                                        td,
                                        td+DIT/86400.,
                                        args=(frequencies[index]))[0])

            # 1e-3 to convert from mJy --> Jy.
            # GRB SED in Jy in the GRB restframe
            grb_fJy = np.array(int_grb_sed) *1e-3 /(DIT/86400.)
            # /(DIT/86400.) to go from Jy.day --> Jy because of the integration
            in time due to the time dependancy.
            it has to be divided by the exposure time to recover an e-/s unit
            for the SNR formula

            grb_fJy *= meiksin(info_dict['wavelength_ang']/10.,
                               info_dict['grb_redshift'])
            grb_fJy *= reddening(info_dict['wavelength_ang'],
                                 info_dict['grb_redshift'],
                                 Av=info_dict['Av_Host']).Li07(3)
            """

        elif info_dict["grb_model"] == "LightCurve":
            # -----------------------------
            #  Loading light curve
            # -----------------------------
            object_path = (
                info_dict["path"]
                + info_dict["object_folder"]
                + info_dict["object_file"]
            )

            #  Fluxes are assumed to be in miliJansky
            LC_data = ascii.read(object_path)

            #  get wavelength
            wvl_list = []
            for dat in LC_data.group_by(["wvl"]).groups.keys:
                wvl_list.append(dat[0])
            sed_stacked = np.zeros(len(wvl_list))

            td = info_dict["t_sinceBurst"]  # in second
            DIT = info_dict["exptime"]  # in second
            # print (td)

            # Sum each wavelength over the time
            time_grb = np.linspace(td, td + DIT, 5)

            for i, wv in enumerate(wvl_list):
                mask_wl = LC_data["wvl"] == wv

                #  interpolation of light curve for each wavelength
                flux_interp = interp1d(
                    LC_data["Time"][mask_wl], LC_data["flux"][mask_wl]
                )
                fluxes = []
                for t in time_grb:
                    fluxes.append(flux_interp(t))
                sed_stacked[i] = np.trapz(fluxes, time_grb)

            # Resample the wavelength
            flux_interp = interp1d(wvl_list, sed_stacked)
            sed_stacked_resampled = flux_interp(info_dict["wavelength_ang"])

            # print (sed_stacked_resampled)
            factor_Jy = 1e-3
            factor_time = 1

            # Convert to Jy
            grb_fJy = (
                sed_stacked_resampled * factor_Jy / (DIT / factor_time)
            )  # GRB SED in Jy in the observer frame
            """
            # /(DIT/86400.) to go from Jy.day --> Jy because of the
            integration in time due to the time dependancy.
            It has to be divided by the exposure time to recover an e-/s
            unit for the SNR formula.
            It is basically the mean value in Jy in the time interval
            corresponding to the exposure time. This assumption is valid for
            rather short exposures and probably not for long exposures,
            it has to be tested.
            """
            # grb_fJy2 = sed_stacked2*1e-3

            # print (grb_fJy,grb_fJy2)
        # erg/s/cm2/A
        flam = utils.fJy_to_flambda(info_dict["wavelength_ang"], grb_fJy)
        # ph/s/cm2/A
        fph = utils.flambda_to_fph(info_dict["wavelength_ang"], flam)
        # e/s
        F_e_s = (
            np.trapz(
                fph * info_dict["system_response"] * info_dict["Trans_atmosphere"],
                info_dict["wavelength_ang"],
            )
            * info_dict["A_tel"]
        )

    info_dict["Object_fph"] = fph
    info_dict["Object_fes"] = F_e_s
    info_dict["Object_mag"] = phot.Jy2Mag(
        info_dict, utils.flambda_to_fJy(info_dict["wavelength_ang"], flam)
    )
    return info_dict

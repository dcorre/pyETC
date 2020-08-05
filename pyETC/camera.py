#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from scipy.interpolate import interp1d
from . import constants as cc
from pyETC.utils import resample


def digitisation_noise(info_dict):
    """ Computes the digitisation noise in e-/px

    Parameters
    ----------
    info_dict: dictionary

    Returns
    --------
    dig_noise: float
               noise due to the digitisation process in e-/px
    """
    dig_noise = (info_dict["cameras"][info_dict["channel"]]["gain"]
                 / np.sqrt(12.0))
    info_dict["dig_noise"] = dig_noise
    return info_dict


def dark_current(info_dict):
    """ Compute the dark current of the camera
    Parameters
    -----------
    info_dict: dictionary

    Returns:
    ---------
    dark: float
          dark current in electrons/sec/pixel
    """
    if info_dict["camera_type"] == "NIR":
        # Thermal noise
        #  Kelvin
        Temp_K = np.array(
            [50.0, 120.0, 130.0, 143.0, 155.0, 170.0, 183.0, 198.0, 235.0]
        )
        #  e-/s/photocells
        dark_current = np.array([2e-3, 3e-3, 1e-2, 1e-1, 1.0,
                                 10.0, 1e2, 1e3, 1e5])

        f = interp1d(Temp_K, dark_current, kind="linear")

        # e-/s/photocells
        Dark = f(info_dict["Temp_cam"] + cc.zero_celsius)
        # e-/s/pixel
        Dark = Dark * info_dict["bin1"] * info_dict["bin2"]

    elif info_dict["camera_type"] == "VIS1" \
            or info_dict["camera_type"] == "VIS2":
        if info_dict["ccd_type"] == "e2v231_84":
            # Celsius
            Temp_C = np.array([-100.0, -75.0, -55.0])
            # e-/s/photocells
            dark_current = np.array([0.0008, 0.006, 0.05])
        elif info_dict["ccd_type"] == "e2v230_84":
            # Celsius
            Temp_C = np.array([-100.0, -75.0, -55.0])
            # e-/s/photocells
            dark_current = np.array([0.00006, 0.0001, 0.001])
        f = interp1d(Temp_C, dark_current, kind="linear")

        # e-/s/photocells
        info_dict["C_th"] = f(info_dict["Temp_cam"])
        # e-/s/pixel
        Dark = info_dict["C_th"] * info_dict["bin1"] * info_dict["bin2"]
    info_dict["DC"] = Dark
    return info_dict


def camera_efficiency(info_dict):
    """ Compute the camera efficiency over the desired wavelength range

    Parameters
    ----------
    info_dict: dictionary

    wavelength: array
                wavelength in angstrom
    Returns
    -------
    eta: array
         efficiency of the camera, [0,1]
    """
    cam_wavelengths = []
    cam_eta = []
    directory = (
        "%s/transmissions/detectors/" % info_dict["path"]
        + info_dict["cameras"][info_dict["channel"]]["camera_type"]
        + "/"
        + info_dict["cameras"][info_dict["channel"]]["sensor"]
        + ".dat"
    )

    with open(directory, "r") as file:
        for line in file:
            if line[0] != "#" and len(line) > 3:
                b, c = line.split()
                cam_wavelengths.append(float(b))
                cam_eta.append(float(c))

    # angstrom
    cam_wavelengths = np.array(cam_wavelengths)
    cam_eta = np.array(cam_eta)

    # Resampling
    eta = resample(cam_wavelengths, cam_eta,
                   info_dict["wavelength_ang"], 0.0, 1.0)

    info_dict["camera_efficiency"] = eta
    return info_dict


def set_camera(info_dict):

    info_dict = digitisation_noise(info_dict)
    # info_dict=dark_current(info_dict)
    # if info_dict['detailed_trans'] == 1:
    #    info_dict=camera_efficiency(info_dict)
    info_dict = camera_efficiency(info_dict)
    return info_dict

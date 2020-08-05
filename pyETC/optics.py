#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from pyETC.utils import resample


def load_optical_element(info_dict, element, element_type,
                         norm=False, norm_val=1.0):
    """Computes the lense transmission with respect to wavelength

   Parameters
   ----------
   info_dict: dictionary


   element: string
        'mirrors' or 'lenses' or 'dichroics' or 'windows'

   element_type: string
        name of the element to be used
        (see '/transmissions/'element'/' for names)

   norm: boolean
        enables to normalise the values to 'norm_val' (default: False)

   norm_val: float
        value used for normalising the data

   Returns
   ---------
   trans :  array
            transmittance of the lense at a given wavelength  (0-1)
   """
    lense_path = "%s/transmissions/%s/%s.txt" % (
        info_dict["path"],
        element,
        element_type,
    )
    File = open(lense_path, "r")
    lines = File.readlines()

    wvl = []
    trans = []

    for line in lines:
        if line[0] != "#" and len(line) > 3:
            bits = line.split()
            trans.append(float(bits[1]))
            wvl.append(float(bits[0]))

    wvl = np.array(wvl)
    trans = np.array(trans, dtype=np.float64)

    # normalie at norm_val
    if norm:
        trans = trans / max(trans) * norm_val
    # Resample the transmission to the
    trans = resample(wvl, trans, info_dict["wavelength_ang"], 0.0, 1.0)

    return trans


def telescope_efficiency(info_dict):
    """ Computes the total transmission of the telescope

   Parameters
   -----------
   info_dict: dictionnary

   Returns
   --------
   Trans_telescope: array
                    transmission of the telescope only

   """
    Trans_tel = 1
    for optical_element in info_dict["optical_design"]["telescope"]:
        for key, value in info_dict["optical_design"]["telescope"][
            optical_element
        ].items():
            Trans_tel *= load_optical_element(info_dict,
                                              optical_element, key)**value
    return Trans_tel


def instrument_channel_efficiency(info_dict):
    """ Computes the total transmission of the instrument chanel

   Parameters
   -----------
   info_dict: dictionnary

   Returns
   --------
   Trans_inst: array
               transmission of the instrument channel without the filter
   """

    Trans_inst = 1
    for optical_element in info_dict["optical_design"][info_dict["channel"]]:
        for key, value in info_dict["optical_design"][info_dict["channel"]][
            optical_element
        ].items():
            Trans_inst *= load_optical_element(info_dict,
                                               optical_element,
                                               key)**value

    return Trans_inst


def set_total_optics_transmission(info_dict):
    """ Computes the total transmission of the optics without filter

   Parameters
   -----------
   info_dict: dictionnary

   Returns
   --------
   Trans_optics: array
                 transmission of the optics without the filter
   """
    Trans_optics = (telescope_efficiency(info_dict) *
                    instrument_channel_efficiency(info_dict))
    return Trans_optics


def set_optics_transmission(info_dict):

    info_dict["Trans_telescope"] = telescope_efficiency(info_dict)
    info_dict["Trans_inst_channel"] = instrument_channel_efficiency(info_dict)
    info_dict["Trans_optics_total"] = set_total_optics_transmission(info_dict)

    return info_dict

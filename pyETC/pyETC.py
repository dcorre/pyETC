# -*- coding: utf-8 -*-

"""Main module."""

import os
import hjson
import numpy as np
from scipy.interpolate import interp1d

from collections import defaultdict, OrderedDict
from .camera import set_camera
from .local_conditions import set_local_conditions
from .preliminary_computations import set_preliminary_computation
from .photometry import set_photometry
from .solver import etc_computation as etc_comp
from .set_object import set_object
from .optics import set_optics_transmission, load_optical_element
from .write_results import write_results
from .store_results import result_plot

__author__ = 'David Corre'
__version__ = 1.0


class etc():
   """ Exposure Time Calculator """

   def __init__(self,path=os.getenv('pyETC_DIR')+'/pyETC',configFile=os.getenv('pyETC_DIR')+'/pyETC/configFiles/example.hjson',config_type='file',name_telescope='default',scale2Airmass=False):
       """ Class constructor """
  
       self.path = path 
       self.configfile = configFile
       self.scale2Airmass=scale2Airmass
 
       self.information=OrderedDict()
       #update settings with defaults
       self.load_telescope_design(name=name_telescope)

       if config_type == 'file':
           # Load the input file in hjson format into a dictionary
           with open(self.configfile,encoding='utf-8') as f:
               self.information.update(hjson.load(f))
       elif config_type == 'data':
           self.information.update(configFile)

       # Make booleans for verbose, make plots and diplay plots
       if self.information['verbose'].lower() == 'true':
           self.information['verbose'] = True
       else:
           self.information['verbose'] = False
       if self.information['plot'].lower() == 'true':
           self.information['plot'] = True
       else:
           self.information['plot'] = False


       # Add information to dictionary 
       self.information['MainDirectory']=path
       self.information['telescope']=name_telescope
       self.information['scale2Airmass']=scale2Airmass

       # Create spectra bin in microns (observer frame)
       lambdas = np.arange(self.information['lambda_start'],self.information['lambda_end']+self.information['lambda_step'],self.information['lambda_step'])
       # In the following all calculations use wavelength in angstroms
       self.information['wavelength_ang'] = lambdas*1e4   # Angstrom

   def load_telescope_design(self, name='default'):
       """ Load telescope params"""

       with open('%s/telescope_database/%s.hjson' % (self.path,name),encoding='utf-8') as f:
           telescope_params=hjson.load(f)
       #Update parameters
       try:
           self.information.update(telescope_params)
       except ValueError:
           raise KeyError("Unknown telescope name: " + name)

   def set(self, **param_dict):
       """Set parameters of the model by name."""
       for key, val in param_dict.items():
           try:
               i = self.information[key]
           except ValueError:
               raise KeyError("Unknown parameter: " + repr(key))
           self.information[key] = val

   def get(self, name):
       """Get parameter of the model by name."""
       try:
           i = self.information[name]
       except:
           raise KeyError("Model has no parameter " + repr(name))
       return self.information[name]

   def show_all_params(self):
       """ Print the value of all parameters """
       for key, val in self.information.items():
           print ("{} : {}".format(key,val))


   def load_info(self):
       
       #Load design info
       if self.information['detailed_trans'] == 1: self.information.update(set_optics_transmission(self.information))

       #Load camera info
       self.information.update(set_camera(self.information))

       #Load Atmosphere transmission, Airmass, seeing, Sky brightness
       self.information.update(set_local_conditions(self.information))

       #Load some infos (transmissions, PSF, pixel_size...)
       self.information.update(set_preliminary_computation(self.information))

       #Load photometric info (zeropoint, effective wavelength...)
       self.information=set_photometry(self.information)

       #Load Object
       if self.information['etc_type'] in ['snr','time']:
           self.information=set_object(self.information)
       return None

   def sim(self):
       self.load_info()

       # Compute the desired quantity (either SNR, total_exposure_time or the magnitude):
       self.information = etc_comp(self.information)

       #------------------
       # Store Main results
       #--------------------
       if self.information['verbose'] == True: write_results(self.information)
       #print (self.information['etc_plot'])
       if self.information['plot'] == True: 
           
           result_plot(self.information,self.information['wavelength_ang'],self.information['fph'], self.information['Object_mag'])

       return None
 
   def load_element(self,element_type,element_name,norm=False,norm_val=1):

       return load_optical_element(self.information,element_type,element_name,norm=norm,norm_val=norm_val)


   def plot_trans(self,y,filename,title='',ylabel='Transmission',ylim=[0,1],wvl_unit='microns',passband_centered=False):
       """ Plot transmisions curves """

       import matplotlib.pyplot as plt
       plt.figure()
       if passband_centered:
           xlim=np.array([self.information['Passband_cuton']*0.8,self.information['Passband_cutoff']*1.2])
       else: 
           xlim=np.array([self.information['wavelength_ang'][0]*0.95,self.information['wavelength_ang'][-1]*1.05])
       
       if wvl_unit == 'microns': 
           x=self.information['wavelength_ang']*1e-4
           xlim*=1e-4
       elif wvl_unit == 'ang': 
           x=self.information['wavelength_ang']

       plt.plot(x,y)
       if wvl_unit == 'microns':  plt.xlabel(r'$\lambda$ ($\mu m$)',size=14)
       elif wvl_unit == 'ang':  plt.xlabel(r'$\lambdai$ (angstroms)',size=14)

       plt.ylabel('%s' % ylabel,size=14)
 
       ylim=np.array(ylim)
       if max(y) > 1:   plt.ylim(ylim*100)
       elif max(y) <=1: plt.ylim(ylim)

       plt.xlim(xlim)
       plt.title('%s' % title,size=14)
       plt.grid()
       plt.savefig('%s.png' % filename)
       plt.show()


   def write_file_trans(self,y,filename,wvl_unit='microns'):
       """ Write the desired element transmission in a file """
       if wvl_unit == 'microns': factor=1e-4
       elif wvl_unit == 'ang': factor=1
       elif wvl_unit == 'nm': factor=1e-1
 
       x=np.around(self.information['wavelength_ang']*factor,2)
       np.savetxt('%s.txt' % filename,np.array([x,y]).T,fmt='%.2f %.4f')




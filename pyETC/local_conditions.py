import numpy as np
from scipy.interpolate import interp1d
from . import constants as cc
from . import utils
#from . import optics as opt
#from . import camera as cam
#from . import photometry as phot
#from . import preliminary_computations as precomp

#######################################


def sky_brightness(info_dict):
   """Compute the sky brightness according the filter and the moon age.

   Parameters
   ----------
   info_dict: dictionary 
    
   scale2Airmass : boolean
                   whether to scale SB to airmass

   Return
   ---------
   sky_b : array 
           Sky brightness at moon age (mag/A/arcsec2)
   """

   # Load the parameters of the dictionary desvribing local conditions
   moon_age = info_dict['moon_age']   # Age of the moon (days)
   elevation = info_dict['elevation'] # Elevation above the horizon (deg)
   skysite = info_dict['sky_site']     # Sky background site symbol
   lambda1 = info_dict['wavelength_ang'][0]            # First wavelength to study (angstrom)
   lambda2 = info_dict['wavelength_ang'][-1]           # Last wavelength to study (angstrom)

   # load Sky brightness
   File=open('%s/local_conditions/sky_brightness/%s.dat' % (info_dict['MainDirectory'],info_dict['sky_site']), "r")
   lines=File.readlines()
   sb=[]
   for line in lines:
       if line[0] != "#" and len(line) > 3:
           bits=line.split()
           sb.append([i for i in bits])
   sb=np.array(sb)
   # ugly hack to convert array of strings to float. Arrays are homogeneous, list are not, but column manipulations are awful with list and not with array... All strings are thus converted to -1 
   sb[0,1]=-1
   for i in range(len(sb)): sb[i,0]=-1
   sb=sb.astype(np.float)

   # Moon age limits
   if (moon_age < 0.):
        moon_age = 0.
   if (moon_age > 14.):
        moon_age = 14
    
   # Interpolation for moon age
   x = sb[0,2:] 
   xxx = sb[1:,1]
   yyy = [] 
   for k in range(len(sb[1:,2:])):
        y = sb[1+k,2:]
        f = interp1d(x,y,kind='linear')
        yyy.append(f(moon_age))
   yyy = np.array(yyy)

   # Interpolation for wavelengths
   l1 = xxx[0]
   l2 = xxx[-1]
   x=[]
   y=[]
   if lambda1 < l1:       # Extrapolation
        x.append(lambda1)
        y.append(yyy[0])

   for k in range(len(xxx)):
        x.append(xxx[k])
        y.append(yyy[k])

   if lambda2 > l2:
        x.append(lambda2)
        y.append(yyy[-1])
   x = np.array(x)
   y = np.array(y)
  
   f = interp1d(x,y,kind='linear')
   sky_b = f(info_dict['wavelength_ang'])      # Sky brightness at moon age (mag/arcsec2) for the given wavelength
   if info_dict['scale2Airmass']: sky_b=SB_airmass(info_dict,sky_b)
   info_dict['sky_brightness']=sky_b

   #Calculate surface brightness at the effective wavelength
   from .photometry import effectiveWavelength
   info_dict=effectiveWavelength(info_dict)
   info_dict['SB_eff']=f(info_dict['effWavelength'])
    
   return info_dict
#----------------------------------------------------------------------------------------------
def airmass(info_dict):
    """ Compute the airmass for the zenith angle of the target

    Parameters
    ----------
    info_dict : dictionary 

    Return
    -------
    airmass: float
             airmass
    """
    # expression valid up to a zenith angle of 60deg
    #X = 1./np.sin(local_cond_dict['elevation'] * np.pi/180.)

    zenith_angle = 90.- info_dict['elevation']    # in degrees
    #sec_z = 1./np.cos(zenith_angle* np.pi/180.)

    # Bemporad's formula
    # correct (at sea level) to within 0.001 for airmass up to about 6.8
    # and at least 1% accuracy for X <10
    #X = sec_z - 0.0018167*(sec_z-1.) -0.002875*(sec_z-1.)**2. -0.000808*(sec_z - 1.)**3.

    # Kasten and Young's formula
    X = 1./ ( np.cos(zenith_angle* np.pi/180.) + 0.50572 * (96.07995 - zenith_angle)**(-1.6364)  )
    info_dict['airmass']=X
    return X,info_dict

#-----------------------------------------------------------------------------
def SB_airmass(info_dict,SB_ref):
   """ scale the sky brightness with the airmass"""
   X=airmass(info_dict)[0]
   # formula found at http://www.cefca.es/jplusetc/etc_help.html
   # Probably only valid for their site   
   SB = SB_ref*(-0.000278719*X*X*X - 0.0653841*X*X + 1.11979*X - 0.0552132)

   return SB
#-----------------------------------------------------------------------------
def seeing(info_dict):
    """ Computes  the seeing for a given airmass and wavelength
    Parameters
    ----------
    info_dict: dictionnary

    wavelength: array
                wavelength in angstroms
    Returns
    -------
    seeing: float
            seeing in arcsec
    """
    # local_cond_dict['seeing_zenith'] is the seeing at the zenith for an airmass of 1 and at wavelength 500nm
    #seeing_arcsec = local_cond_dict['seeing_zenith'] * X**(3./5)
    seeing_arcsec = info_dict['seeing_zenith'] * info_dict['airmass']**(3./5) * (info_dict['effWavelength']/5000.)**(-0.2)
    info_dict['seeing_los_arcsec']=seeing_arcsec
    return info_dict
#---------------------------------------------------------------------------------------------

def atmospheric_transmission(info_dict):
   """ Compute the atmospheric transmission for the site of San Pedro de Martir. Airmass is also computed 
   Parameters
   ----------
   info_dict : dictionary 

   wavelength : array 
                wavelengths in angstrom

   Return
   -------

   extinction_curve: array
                     extinction in mag/airmass for given wavelengths

   Trans_atm: array
              atmospheric transmission efficiency for the given wavelengths
   """

   if info_dict['atmosphere_type'] == 'data': 
       lambda1 = info_dict['wavelength_ang'][0]     # First wavelength to study (angstrom)
       lambda2 = info_dict['wavelength_ang'][-1]    # Last wavelength to study (angstrom)
       # Load extinction 
       File=open('%s/local_conditions/atmosphere/%s.dat' % (info_dict['MainDirectory'],info_dict['ext_file']), "r")
       lines=File.readlines()
       wvl=[]
       extinction=[]
       for line in lines:
           if line[0] != "#" and len(line) > 3:
               bits=line.split()
               wvl.append(float(bits[0]))
               extinction.append(float(bits[1]))

       wvl=np.array(wvl)
       extinction=np.array(extinction)

       # Interpolation for wavelengths
       l1 = wvl[0]
       l2 = wvl[-1]
       x=[]
       y=[]
       if lambda1 < l1:       # Extrapolation
           x.append(lambda1)
           #y.append(extinction[0])
           y.append((extinction[1]-extinction[0])/(wvl[1] - wvl[0])*lambda1+extinction[0]-(extinction[1]-extinction[0])/(wvl[1] - wvl[0])*wvl[0])
       for k in range(len(wvl)):
           x.append(wvl[k])
           y.append(extinction[k])

       if lambda2 > l2:
           x.append(lambda2)
           y.append(extinction[-1])

       x = np.array(x)
       y = np.array(y)

       f = interp1d(x,y,kind='linear')
       extinction_array = f(info_dict['wavelength_ang'])      # Extinction (mag/airmass) for the given wavelengths

       info_dict['Atmosphere_extinction']=extinction_array

       # alternative: Model (schuster / parrao 2001 the atmospheric extinction of san pedro martir)
       #extinction_array = 0.0254
    
       #-----------------------------------
       # Compute the atmospheric absorption
       trans_atm = 10**(-0.4*extinction_array* airmass(info_dict)[0])


   elif info_dict['atmosphere_type'] == 'file':
       atm_path = '%s/local_conditions/atmosphere/%s.txt' % (info_dict['MainDirectory'],info_dict['atm_file'])
       File=open(atm_path, "r")
       lines=File.readlines()
       wvl=[]
       trans=[]

       for line in lines:
           if line[0] != "#" and len(line) > 3:
               bits=line.split()
               trans.append(float(bits[1]))
               wvl.append(float(bits[0]))

       wvl=np.array(wvl)*1e4
       trans=np.array(trans, dtype=np.float64)

       # Resample the transmission to the 
       trans_atm = utils.resample(wvl,trans,info_dict['wavelength_ang'],0.,1.)

   info_dict['Trans_atmosphere']=trans_atm
   return info_dict
#-----------------------------------------------------------------------------------------------

def sky_countrate(info_dict):
    """Compute the sky brightness according the filter and the moon age.

    Parameters
    ----------
    info_dict : dictionary 

    wavelength : array 
                 wavelengths in angstrom

    Returns
    -------
    Sky_countrate: float
                   sky countrate in e-/s/px

    """
    """
    if info_dict['detailed_trans'] == 1:
         
         optics_trans = opt.set_transmission(info_dict,wavelength)
         filter_trans = phot.set_filter(info_dict,wavelength)
         cam_eta = cam.camera_efficiency(info_dict,wavelength)
         #atm_eta = extinction_atm(info_dict,wavelength)[1]
         system_response = optics_trans *filter_trans*cam_eta
         
    else: 
         system_response=phot.set_filter(info_dict,wavelength) 
    """
    info_dict=sky_brightness(info_dict)
    #Sky_fph = phot.zeromag_to_flux(info_dict,unit='ph',phot_sys='AB')*10.**(-0.4 * info_dict['sky_brightness'])   # (ph/s/cm2/A/arsec2)
    # The magnitudes are given in the AB photometric system
    Sky_fph = utils.flambda_to_fph(info_dict['wavelength_ang'],utils.fJy_to_flambda(info_dict['wavelength_ang'],3631))*10.**(-0.4 * info_dict['sky_brightness'])   # (ph/s/cm2/A/arsec2)
    #Sky_fph = utils.flambda_to_fph(info_dict['wavelength_ang'],utils.fJy_to_flambda(info_dict['wavelength_ang'],info_dict['Flux_zero_Jy']))*10.**(-0.4 * info_dict['SB_eff'])   # (ph/s/cm2/A/arsec2)
    
    #Sky_countrate = np.trapz(Sky_fph * info_dict['Trans_atmosphere']*info_dict['system_response'],info_dict['wavelength_ang']) * info_dict['A_tel'] * info_dict['A_pixel'] # e/s/px
    Sky_countrate = np.trapz(Sky_fph *info_dict['system_response'],info_dict['wavelength_ang']) * info_dict['A_tel'] * info_dict['A_pixel'] # e/s/px
    info_dict['Sky_CountRate']=Sky_countrate
    return info_dict

#-----------------------------------------------------------------------------------------------
def set_local_conditions(info_dict):
    
    info_dict=sky_brightness(info_dict)
    info_dict=airmass(info_dict)[1]
    info_dict=seeing(info_dict)
    info_dict=atmospheric_transmission(info_dict)

    return info_dict

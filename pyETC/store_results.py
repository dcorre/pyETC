import matplotlib.pyplot as plt
from . import camera as cam 
from . import local_conditions as local
from . import optics as opt
from . import photometry as phot

def result_plot(info_dict,wavelength,fph,obj_mag):
    """ Return a multiplot summarising the computation.

    Parameters:
    -----------
    info_dict: dictionary

    wavelength: array
                wavelength in angstrom

    Returns
    -------
    Multiplot saved as resut.png
    """

    # print all values stored in info_dict in 'results/info.txt'
    File = open(info_dict['MainDirectory']+'/results/info.txt','w')
    for key, value in info_dict.items():
         print ('%s: %s' % (key, value), file=File)
    File.close()

    # Load transmissions
    #-------------------
    if info_dict['detailed_trans']==1:
         trans_atms = info_dict['Trans_atmosphere']
         optics_trans = opt.telescope_efficiency(info_dict)*opt.instrument_channel_efficiency(info_dict)
         filter_trans = info_dict['Trans_filter']
         cam_eta = info_dict['camera_efficiency']
         sys_eff = trans_atms * optics_trans * filter_trans * cam_eta
    else:
         cam_eta = info_dict['camera_efficiency']
         optics_trans =phot.set_filter(info_dict)
         sys_eff = optics_trans *cam_eta
    plt.figure(figsize=(15,10))
    
    # Object ph/s/cm2/A
    ax=plt.subplot(3,1,1)
    plt.plot(wavelength,fph,'r-',label='Flux outside atmosphere')
    if info_dict['detailed_trans'] == 1: plt.plot(wavelength,fph*trans_atms*optics_trans*filter_trans,'b-',label='Flux reaching camera')
    else: plt.plot(wavelength,fph*optics_trans,'b-',label='Flux reaching camera')
    plt.xlim(wavelength[0],wavelength[-1])
    plt.xlabel('Wavelength (Angstrom)',size=18)
    plt.ylabel('photons/sec/cm2/A',size=18)
    plt.legend(loc="upper right",fontsize=15)
    plt.tick_params(labelsize=16)
    plt.grid()

    #Object magnitude
    ax=plt.subplot(3,1,2)
    plt.plot(wavelength,obj_mag,'r-')
    #plt.ylim(15,30)
    plt.xlim(wavelength[0],wavelength[-1])
    plt.xlabel('Wavelength (Angstrom)',size=18)
    plt.ylabel('%s mag top atmosphere' % info_dict['photometry_system'],size=18)
    plt.gca().invert_yaxis()
    plt.tick_params(labelsize=16)
    plt.grid()
    # Transmissions
    ax=plt.subplot(3,1,3)
    if info_dict['detailed_trans'] == 1:
         plt.plot(wavelength, trans_atms,'b-.',label='T_atm at elevation')
         plt.plot(wavelength, optics_trans,'b:',label='T_opt (mirror + lenses + obs)')
         plt.plot(wavelength, filter_trans,'b-',label='T_filter')
         plt.plot(wavelength, cam_eta, 'm-', label='QE_ccd')
    else:
         plt.plot(wavelength, optics_trans,'b:',label='T_total_opt')
         plt.plot(wavelength, cam_eta, 'm-', label='QE_ccd')
    plt.plot(wavelength, sys_eff, color='black',label='Trans_tot')
    plt.xlim(wavelength[0],wavelength[-1])
    plt.ylim(0.,1.)
    plt.xlabel('Wavelength (angstrom)',size=18)
    plt.ylabel('Transmission',size=18)
    plt.legend(loc='upper right',fontsize=15)
    #plt.legend(loc='lower center ',bbox_to_anchor=(1., -0.15), ncol=2, fancybox=True, shadow=True)
    plt.grid()
    plt.tick_params(labelsize=16)
    plt.tight_layout()
    plt.savefig(info_dict['MainDirectory']+'/results/plot.png') #% (dir_python,plotname))
    

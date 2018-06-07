import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import local_conditions as local
import sys

wvl = np.arange(0.,25000.,5000.)
local_dict = {'lambda1': wvl[0], 'lambda2': wvl[-1],'elevation': 90.}

extinction_array, trans_atm = local.extinction_atm(local_dict,wvl)

mirror_trans_ag = np.array([0., 0.,  0.65, 0.83,  0.925, 0.978, 0.978, 0.975, 0.975, 0.978, 0.980, 0.979, 0.988, 0.988])
mirror_wavelengths_ag = np.array([0., 170., 200.,  300.,   400.,  500.,    600.,   700.,   800.,   900.,  1000.,   1100.,  2100., 3000.])    # (microns)

f_mirror = interp1d(mirror_wavelengths_ag,mirror_trans_ag,kind='linear')

lense_trans = np.array([0.,0.,0.0, 0.0, 0.06, 0.10, 0.40, 0.60, 0.80, 0.90, 0.91, 0.92, 0.93, 0.94, 0.94, 0.90, 0.93, 0.92, 0.91, 0.88, 0.80, 0.80, 0.77, 0.59,   0.])
lense_wavelengths = np.array([0., 170., 200., 270., 290.,  300.,  315.,  325.,  340.,  380.,  400.,  600.,  800.,  1000., 1350., 1400., 1450., 1600., 1800., 2000., 2200., 2350., 2400., 2600., 2850.])  # (microns)

f_lense = interp1d(lense_wavelengths,lense_trans,kind='linear')

directory = '../detectors/VIS/e2v231_84/BEX2.dat'
cam_wavelengths_vis=[]
cam_eta_vis=[]
with open(directory, 'r') as file:
    for line in file:
         if line[0] != "#" and len(line) > 3:
              b, c = line.split()
              cam_wavelengths_vis.append(float(b))
              cam_eta_vis.append(float(c))

cam_wavelengths_vis = np.array(cam_wavelengths_vis)    # angstrom
cam_eta_vis = np.array(cam_eta_vis)
f_ccd_vis = interp1d(cam_wavelengths_vis,eta_ccd_vis,kind='linear')


directory = '../detectors/'+info_dict['camera_type']+'/hgcdte.dat'
cam_wavelengths_ir=[]
cam_eta_ir=[]
with open(directory, 'r') as file:
    for line in file:
         if line[0] != "#" and len(line) > 3:
              b, c = line.split()
              cam_wavelengths_ir.append(float(b))
              cam_eta_ir.append(float(c))

cam_wavelengths_ir = np.array(cam_wavelengths_ir)    # angstrom
cam_eta_ir = np.array(cam_eta_ir)
f_ccd_ir = interp1d(cam_wavelengths_ir,eta_ccd_ir,kind='linear')

Trans_optic = trans_atm * Trans_optics

if folder == 'bessel' or folder == 'sloan' or folder =='stroemgren':
    trans_ccd=f_ccd_vis(wvl)
    plt.xlim(0.,15000.)
elif folder =='wircam':
    trans_ccd = f_ccd_ir(wvl)
    plt.xlim(0.,25000.)


folder = str(sys.argv[1])
#folder = 'bessel'
#folder = 'stroemgren'
#folder = 'sloan'
if folder == 'bessel':
    bands = ['U','B','V','R','I']
elif folder == 'sloan':
    bands = ['u','g','r','i','z']
elif folder == 'stroemgren':
    bands = ['u','v','b','y']
elif folder == 'wircam':
    bands = ['Y','J' 'H']

for band in bands:
    print (band)
    xx = []
    yy = []
    line = 0
    filename=folder+'_'+band+'.txt'
    with open('%s/%s' % (folder,filename), 'r') as file:
         for row in file:
              if line > 1:
                   a, b = row.split()
                   xx.append(float(a))
                   yy.append(float(b))
              line+=1 

    ang=np.array(xx)    # Angstrom
    Trans=np.array(yy)         # 
    local_dict = {'lambda1': ang[0], 'lambda2': ang[-1],'elevation': 90.}

    extinction_array, trans_atm = local.extinction_atm(local_dict,ang)

    if folder == 'bessel' or folder == 'sloan' or folder =='stroemgren':
         trans_ccd=f_ccd_vis(ang)
         plt.xlim(0.,15000.)
         if design == 1:  # direct imaging
              Trans_optics = f_mirror(ang/10.)**3. * f_lense(ang/10.)**2. *0.98**3.
         elif design == 2:   # didier concept
              Trans_optics = f_mirror(ang/10.)**7. * 0.98**3.

    elif folder =='wircam':
         trans_ccd = f_ccd_ir(ang)
         if design == 1:  # direct imaging
             Trans_optics = f_mirror(ang/10.)**7. * f_lense(ang/10.)* 0.98**3.
         elif design == 2:   # didier concept
             Trans_optics = f_mirror(ang/10.)**7. *0.98**3.
    plt.plot(ang,Trans,label=band)
    plt.plot(ang,Trans*Trans_optic*trans_ccd,'--',label='f')
1

plt.plot(wvl,trans_opt, '--',label='optic + ccd')
#plt.plot(wvl,trans_ccd,'--',color='black',label='ccd')
plt.ylim(0.,1.)
plt.xlabel(r'wavelenghts (nm)')
plt.ylabel('Transmission ')
plt.legend(loc='lower right')
plt.title('%s' % folder)
plt.grid(True)
plt.savefig('filter_%s.png' % folder)
plt.show()



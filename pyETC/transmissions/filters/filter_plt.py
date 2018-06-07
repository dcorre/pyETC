import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
#import local_conditions as local
import sys

def plot_colorfilter(band):
    """ Associate a color to a given filter for nice plotting 
    Parameters 
    ----------
    band: string
          filter band ie 'u','g',...

    Returns
    -------
    band_color: string
                color associated with the band filter ie 'u' with blue
    """

    if band == 'u':
         color_band='purple'
    elif band == 'g':
         color_band='blue'
    elif band == 'r':
         color_band = 'green'
    elif band == 'i':
         color_band = 'orange'
    elif band == 'zs':
         color_band = 'salmon'
    elif  band == 'z':
         color_band = 'salmon'
    elif band == 'y':
         color_band = 'chocolate'
    elif band == 'Y':
         color_band = 'red'
    elif band == 'J':
         color_band = 'maroon'
    elif band == 'H':
         color_band = 'black'

    return color_band


nb_filter_sys=1

folder = str(sys.argv[1])
try:
    folder2 = str(sys.argv[2])
    nb_filter_sys+=1
except:
    print ('')
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
    bands = ['Y','J', 'H']#, 'Ks']
elif folder == 'panstarrs':
    bands = ['g','r','i','z']
elif folder == 'des':
    bands = ['u','g','r','i','z','y']


for band in bands:
    print (band)
    xx = []
    yy = []
    line2 = 1
    filename=folder+'_'+band+'.txt'
    with open('%s/%s' % (folder,filename), 'r') as file:
         for line in file:
              #if line[0] != "#" and len(line) > 3:
              if line2 > 2:
                   a, b = line.split()
                   xx.append(float(a))
                   yy.append(float(b))
              line2+=1 

    ang=np.array(xx)    # Angstrom
    Trans=np.array(yy)
    #print (ang*10,Trans/100.)
    #if folder == 'panstarrs':    Trans=Trans*1e2        # 


    plt.plot(ang*10,Trans/100.,label=band,color=plot_colorfilter(band))


if nb_filter_sys==2:
    if folder2 == 'bessel':
         bands2=['U','B','V','R','I']
    elif folder2 == 'sloan':
         bands2=['u','g','r','i','z']
    elif folder2 == 'stroemgren':
         bands2=['u','v','b','y']
    elif folder2 == 'wircam':
         bands2=['Y','J', 'H']#, 'Ks']
    elif folder2 == 'panstarrs':
         bands2=['g','r','i','z']
    elif folder2 == 'des':
         bands2=['u','g','r','i','z','y']

    for band in bands2:
         print (band)
         xx = []
         yy = []
         line2 = 1
         filename=folder2+'_'+band+'.txt'
         with open('%s/%s' % (folder2,filename), 'r') as file:
              for line in file:
                   #if line[0] != "#" and len(line) > 3:
                   if line2 > 2:
                        a, b = line.split()
                        xx.append(float(a))
                        yy.append(float(b))
                   line2+=1 

         ang=np.array(xx)    # Angstrom
         Trans=np.array(yy)
         #print (ang*10,Trans/100.)
         #if folder == 'panstarrs':    Trans=Trans*1e2        # 
    

         plt.plot(ang*10,Trans/100.,label=band,color=plot_colorfilter(band))



plt.xticks(np.arange(4000,20000,2000))
#plt.yticks(np.arange(0,1.1,0.1))
#plt.yticks(np.arange(0.1, 1.1, 0.2),minor =True)
plt.xlim(3000,22000)
plt.ylim(0.,1.)
plt.xlabel(r'$\lambda$ (Angstroms)',fontsize=15)
plt.ylabel('Transmission ',fontsize=15)
plt.legend(loc='lower right')
plt.grid(True,which='major',lw=1)
plt.grid(True,which='minor',lw=0.3)
plt.minorticks_on()
if nb_filter_sys==1:
    plt.title('%s' % folder)
    plt.savefig('filter_%s.png' % folder)
elif nb_filter_sys==2:
    #plt.title('%s + %s' % (folder,folder2))
    plt.savefig('filter_%s_%s.png' % (folder,folder2))

plt.show()



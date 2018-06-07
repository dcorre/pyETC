import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker

#Create the figure
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
#Set the x and y axis ticks frequency
ax.set_yticks(np.arange(0,1.1,0.1))
ax.set_yticks(np.arange(0, 1.1, 0.05),minor =True)
#Set x and y range
ax.axis([0.3,2, 0,1 ])



files_suffix='sys_'
bands=['g','r','i','z','y','J','H']
files_prefix=['_baseline','_goal']
colors=['purple','blue','green','orange','red','maroon','black']

for f_p in files_prefix:
    for i,band in enumerate(bands):
        f=files_suffix+band+f_p+'.txt'
        data = np.genfromtxt(f,names=['wvl','trans'])

        if f_p == '_baseline': 
           ls='-'
           plt.plot(data['wvl'],data['trans'],ls=ls,color=colors[i],label=band)
        elif f_p == '_goal': 
           ls='--'
           plt.plot(data['wvl'],data['trans'],ls=ls,color=colors[i])


plt.xlabel(r'$\lambda$ ($\mu m$)',fontsize=15)
plt.ylabel('Transmission',fontsize=15)
#plt.title('')
plt.grid(True,which='major',lw=1.3)
plt.grid(True,which='minor',lw=0.3)
plt.legend(loc='upper right')
plt.savefig('sys_eff_all.png')



plt.show()

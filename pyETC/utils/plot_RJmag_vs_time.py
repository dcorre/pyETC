import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker



#Create the figure
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
#Set the x and y axis ticks frequency
ax.set_yticks(np.arange(4,28,2))
ax.set_yticks(np.arange(3, 29, 2),minor =True)
#Set x and y range
ax.axis([0.07,3600 , 28, 14])
# Set x axis to log
ax.set_xscale('log')

ax.axvline(300, color='blue', linestyle='--',lw=1.5)
ax.axhline(22, color='green', linestyle='--',lw=1.5)
ax.axhline(20, color='red', linestyle='--',lw=1.5)
plt.text(300,14.3,'300s',rotation=0,color='blue',fontsize=15)
plt.text(0.2,22.15,'R=22.0',rotation=0,color='green',fontsize=15)
plt.text(0.2,20.15,'J=20.0',rotation=0,color='red',fontsize=15)

texp=np.array([0.1,0.5,1,2,5,10,20,30,45,60,90,120,150,180,210,240,270,300,360,420,480,540,600,720,840,960,1080,1300,1500,1800,2100,2500,3000,3600])
R_baseline=np.array([15.88,17.60,18.31,18.99,19.81,20.37,20.86,21.13,21.38,21.56,21.80,21.97,22.10,22.20,22.29,22.36,22.43,22.49,22.45,22.55,22.64,22.71,22.78,22.86,22.96,23.04,23.11,23.22,23.30,23.41,23.49,23.59,23.69,23.79])
R_goal=np.array([16.12,17.82,18.53,19.19,19.99,20.53,21.01,21.27,21.52,21.69,21.93,22.09,22.22,22.32,22.41,22.49,22.55,22.61,22.58,22.68,22.77,22.84,22.90,22.98,23.08,23.16,23.24,23.35,23.43,23.53,23.62,23.71,23.81,23.91])
J_baseline=np.array([14.57,16.22,16.87,17.46,18.15,18.61,19.04,19.27,19.51,19.67,19.90,20.06,20.18,20.28,20.37,20.44,20.50,20.56,20.55,20.65,20.73,20.81,20.87,20.94,21.03,21.12,21.19,21.30,21.38,21.48,21.56,21.66,21.76,21.86])
J_goal=np.array([14.89,16.51,17.14,17.70,18.36,18.80,19.22,19.45,19.68,19.84,20.07,20.23,20.35,20.45,20.54,20.61,20.67,20.73,20.72,20.82,20.91,20.98,21.04,21.11,21.20,21.29,21.36,21.47,21.55,21.65,21.73,21.83,21.93,22.03])

print (len(texp),len(R_baseline))
plt.plot(texp,R_baseline,color='green',label='R baseline')
plt.plot(texp,R_goal,color='green',ls='--',label='R goal')
plt.plot(texp,J_baseline,color='red',label='J baseline')
plt.plot(texp,J_goal,color='red',ls='--',label='J goal')


plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.gca().invert_yaxis()
plt.xlabel('Exposure time (s)',fontsize=15)
plt.ylabel('AB magnitude',fontsize=15)
#plt.title('')
plt.grid(True,which='major',lw=1.3)
plt.grid(True,which='minor',lw=0.3)
plt.legend(loc='upper left')
plt.savefig('RJ_mag_vs_time.png')
plt.show()

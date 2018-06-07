import numpy as np
import sys
import matplotlib.pyplot as plt

filename = sys.argv[1] #alpha_lyr_stis_006
start = float(sys.argv[2])
end = float(sys.argv[3])


xx = []
yy = []
with open(filename+'.txt', 'r') as file:
    for row in file:
         a, b = row.split()
         xx.append(float(a))
         yy.append(float(b))

microns=np.array(xx)*1e-4    # Angstrom to microns
Flux=np.array(yy)         # erg/s/cm2/A

plt.plot(microns,Flux)
plt.xlim(start,end)
plt.xlabel(r'wavelenghts ($\mu m$)')
plt.ylabel('Flux (erg/s/cm2/A) ')
plt.show()

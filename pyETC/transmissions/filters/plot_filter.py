import numpy as np
import matplotlib.pyplot as plt
import sys


filename = sys.argv[1]


xx = []
yy = []
line = 0
with open(filename, 'r') as file:
    for row in file:
         if line > 1:
              a, b = row.split()
              xx.append(float(a))
              yy.append(float(b))
         line+=1 
ang=np.array(xx)    # Angstrom
Trans=np.array(yy)         # %

plt.plot(ang,Trans)
#plt.xlim(start,end)
plt.xlabel(r'wavelenghts (Angstroms)')
plt.ylabel('Transmission ')
plt.savefig('filter_%s.png' % filename)
plt.show()



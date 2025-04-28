"""
Created in July 2021

@author: Jean Kempf, Francois Rincon

The test is inspired from the following paper:
Parrish, Ian J., et al. "The effects of anisotropic viscosity on turbulence and heat transport in the intracluster medium." Monthly Notices of the Royal Astronomical Society 422.1 (2012): 704-718.
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

fid=open('./output/timevol.dat',"r")
# read the first line to get data names
varnames=fid.readline().split()
fid.close()
# load the bulk of the file
data=np.loadtxt('./output/timevol.dat',skiprows=1)
# store this in our data structure
V={}
i=0
for name in varnames:
    V[name]=data[:,i]
    i=i+1

plt.figure()
plt.plot(V['t'], V['kinx'])
plt.plot(V['t'], V['kiny'])
plt.savefig('kin.png')
plt.show()


'''
This Code calculates density for the deprojected data from normalization using the formula:
norm = (10^-14 / 4*pi*D^2) * Integral(n_e*n_H dV)
Where D is the distance to the source (cm) => Angular Diameter Distance.
n_e and n_H are electron and H densities (cm^-3)
It is integarted over each shell.
'''

import numpy as np
from header import *
import math

r = np.loadtxt("radius.dat",unpack=True)
norm = np.loadtxt("norm_DSDEPROJ.dat",unpack=True)

r = r.tolist()
r = r[::-1]
r.append(0)
r = r[::-1]
r = np.array(r)*kpc

# Shell Volume Calculation

Vol = []
for i in range(len(norm)):
	temp = (4./3.)*np.pi*(r[i+1]**3-r[i]**3)
	Vol.append(temp)

Vol = np.array(Vol)

# Norm to Density conversion

nn = np.array((norm*(4*math.pi*(Da*(1+zz))**2)/(Vol*(1e-14)))) # Da is the angular diameter distance
n = np.sqrt(nn*mu_e*mu_i/(mu**2))
n_e = mu*np.array(n)/mu_e

np.savetxt("density_DSDEPROJ.dat",n_e)

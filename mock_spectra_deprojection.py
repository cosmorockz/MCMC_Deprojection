'''
This code runs the spec_deproj script (Deprojection Script) for every energy channel. Then it saves all the deprojected spectra in the folder /deprojected.
'''

import numpy as np
from header import *
from spec_deproj import spectra_deprojection
import math
from astropy.io import fits
import os
import shutil

r = np.loadtxt("radius.dat",unpack =True)
r = r*kpc # The radius was in kpc
number_of_channels = 1980

path = os.getcwd()
newpath = path + str('/deprojected')
if os.path.exists(newpath) == True:
	shutil.rmtree(newpath)
if not os.path.exists(newpath):
	os.makedirs(newpath)

r = r.tolist()
r = r[::-1]
r.append(0)
r = r[::-1]
Vol = []

# Volume of each of the shells

for i in range(len(r)-1):
	temp = (4./3.)*np.pi*(r[i+1]**3-r[i]**3)
	Vol.append(temp)

# Deprojecting the spectra for every energy range

for i in range(number_of_channels):
	C = spectra_deprojection("spectra1/channel"+str(i+1)+".dat",r) # Gives the deprojected count rate
	spec = np.array(C)*np.array(Vol) # Multiplies the count rate with volume to get the total counts
	np.savetxt("deprojected/channel"+str(i+1)+".dat",spec)



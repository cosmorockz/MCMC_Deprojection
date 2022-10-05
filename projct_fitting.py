'''
This code fits the non-deprojected Data and saves the Temperature and Normalization.
'''

import numpy as np
from header import *
from xspec import *
import os
from astropy.io import fits

r = np.loadtxt("radius.dat",unpack=True)
n_annulus = len(r)

f = open("script_pha","w")
f.write("#!/bin/bash\n")
f.write("ls *.pha >File2\n")
f.close()

os.system("chmod +x script_pha")
os.system("./script_pha")

f = open("File2")
ff = f.readlines()
files = [i.strip('\n') for i in ff]
norm,T = [],[]

for i in range(n_annulus):
	Xset.abund = 'aspl' # The model for incorporating metallicity
	s = Spectrum(str(i+1)+".pha") # Loading the spectrum file
	m = Model("mekal") # Setting the model to mekal
	m.setPars(2.0,nH,Z,zz,0,1e-4)
	'''
	par1	 = 	plasma temperature in keV => our initial guess 2.0

	par2	 = 	hydrogen density in cm^-3 => neutral hydrogen density. we freeze it to nH=1e-06


	par3	 = 	Metal abundances (He fixed at cosmic). The elements included are


			C, N, O, Ne, Na, Mg, Al, Si, S, Ar, Ca, Fe, Ni. Abundances are


			set by the abund command. => we fix it to Z=0.2


	par4	 = 	redshift, zz=0.001


	par5	 =      0 -> calculate, 1 -> interpolate => our choice 0


	par6	 =	normalization => we calculate the density from it => initial guess = 1e-04
	'''
	Fit.method = "leven 100 1e-60" # Levenberg Marquadt method with tolerance 1e-60 for fitting
	Fit.perform()
	norm.append(m(6).values[0]) # Parameter 6 for normalization => density
	T.append(m(1).values[0]) # Parameter 1 for temperature
	AllData.clear()


np.savetxt("norm_PROJCT.dat",norm)
np.savetxt("temp_PROJCT.dat",T)
	



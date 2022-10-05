'''
This code creates the .pha files from the .dat files in /deprojected1 folder and then fits them using the MEKAL plasma model in XSpec and saves the temperature and normalization.
'''

import numpy as np
from header import *
from xspec import *
import os
from astropy.io import fits

n_channels = 1980

# Bash Script to copy all the .pha files and their response file to /deprojected1 folder

f = open("script_copy","w")
f.write("#!/bin/bash\n")
f.write("cp *.pha deprojected1\n")
f.write("cp *.rsp deprojected1\n")
f.write("ls *.pha >File2\n")
f.write("\n")
f.close()
os.system("cd spectra")
os.system("chmod +x script_copy")
os.system("./script_copy")

f = open("File2")
ff = f.readlines()
files = [i.strip('\n') for i in ff]
n_annuli = len(files)

# Takes the .pha files for the non-deprojected spectra and replaces the CHANNEL,RATE and STAT_ERR column. Keeps the CHANNEL and STAT_ERR column same as the actual .pha file. Just changes the RATE column with the deprojected RATE

location = "deprojected1"
path = os.getcwd()
for j in range(n_annuli):
	newpath = path+str('/')+str(location)+str('/annulus')+str(j+1)+str('.pha')
	if os.path.exists(newpath):
		os.remove(newpath)
	spec = np.loadtxt(str(location)+"/annulus"+str(j+1)+".dat",unpack=True)
	hdu = fits.open(files[j])
	channel = hdu[1].data.CHANNEL
	rate = spec
	stat_err = hdu[1].data.STAT_ERR
	data = [(channel[i], rate[i], stat_err[i])for i in range(n_channels)]
	data1 = np.array(data, dtype = [('CHANNEL' ,'>i4'),('RATE','>f4'),('STAT_ERR','>i2')])
	hdu[1].data = data1
	hdu.writeto(str(location)+"/annulus"+str(j+1)+".pha")

# Fitting Using XSpec

norm,T = [],[]
for i in range(n_annuli):
	s = Spectrum("deprojected1/"+"annulus"+str(i+1)+".pha") # Loading the spectrum file
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
	Xset.abund = 'aspl' # The model for incorporating metallicity
	#Fit.statMethod = "cstat"
	#Fit.statTest = "pchi"
	Fit.method = "leven 100 1e-60" # Levenberg Marquadt method with tolerance 1e-60 for fitting
	Fit.perform()
	norm.append(m(6).values[0]) # Parameter 6 for normalization => density
	T.append(m(1).values[0]) # Parameter 1 for temperature
	AllData.clear()


np.savetxt("norm_DSDEPROJ.dat",norm)
np.savetxt("temp_DSDEPROJ.dat",T)



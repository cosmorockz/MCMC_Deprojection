'''
Extract the count RATE column from the .pha files, stores them .dat files. Saves the .dat files in a new folder /spectra
'''

import numpy as np
from astropy.io import fits
import os
import shutil

n_channels = 1980 # Number of energy channels

# Bash Script that lists the name of all the files in the folder in a single file.
f = open("script_pha","w")
f.write("#!/bin/bash\n")
f.write("ls *.pha >File2\n")
f.close()
os.system("chmod +x script_pha")
os.system("./script_pha")
f = open("File2")

ff = f.readlines()
files = [i.strip('\n') for i in ff]
n_annulus = len(files)

def dat_generator(n_annuli):
	path = os.getcwd()
	newpath = path +"/spectra" 
	# Creating a new folder to save the .dat files
	if os.path.exists(newpath) == True:
		shutil.rmtree(newpath)
	if not os.path.exists(newpath):
		os.makedirs(newpath)
	for spectra in files:
		newpath1 = spectra
		newpath2 = newpath+'/'+spectra.strip('.pha')+str('.dat')
		if os.path.exists(newpath2):
			os.remove(newpath2)
		hdu = fits.open(newpath1)
		spec = hdu[1].data.RATE # Extracting the count RATE column
		np.savetxt(newpath2,spec)

dat_generator(n_annulus)


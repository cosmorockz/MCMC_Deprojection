'''
Creates .pha files using XSpec from .spec files.
'''

import numpy as np
from astropy.io import fits
import os
from xspec import *

# Lists all the spectral file names in a single file

f = open("script1","w")
f.write("#!/bin/bash\n")
f.write("ls *.spec >File1\n")
f.close()
os.system("chmod +x script1")
os.system("./script1")

n_channels=1980
f = open("File1")
ff = f.readlines()
files = [i.strip('\n') for i in ff]

for spectra in files:
	e,counts = np.loadtxt(spectra,usecols=(0,1),unpack=True)
	lower = e-4e-3
	upper = e+4e-3
	counts = np.array(counts)*1e+3
	error = [0 for i in range(len(counts))]
	filename = spectra.strip('.spec') + str('.dat')
	np.savetxt(filename,np.c_[lower,upper,counts,error])
	pha_name = filename.strip('.dat')+str('.pha') # .pha filename
	resp_name = filename.strip('.dat')+str('.rsp') # response filename
	f = open("script2","w")
	f.write("#!/bin/bash\n")
	f.write('flx2xsp '+ str(filename)+' '+ str(pha_name)+' '+ str(resp_name) +' yunit=ph/cm^2/s/MeV clobber=yes') # Using flx2xsp command to create the data and response file from .spec files
	f.close()
	os.system("chmod +x script2")
	os.system("./script2")




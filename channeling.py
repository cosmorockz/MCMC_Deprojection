'''
Takes the .dat files from /spectra folder. Converts from "specrta (energy range vs counts) of each annuli" to "counts vs annuli" for each energy range. Saves the resulting files /spectra1 folder.
'''

import numpy as np
import math
import os
import shutil
import sys

number_of_channels = 1980

# Bash Script to list all the .dat files in /spectra folder in a single file.

f = open("script_pha","w")
f.write("#!/bin/bash\n")
f.write("cd spectra\n")
f.write("ls *.dat >File2\n")
f.write("cp File2 ..\n")
f.close()
os.system("chmod +x script_pha")
os.system("./script_pha")

f = open("File2")
ff = f.readlines()
files = [i.strip('\n') for i in ff]
n_annulus = len(files)
os.system("cd ..")

def annuli_to_channel(n_annuli,statement):
	path = os.getcwd()
	newpath1 = path +"/"+str(statement)
	newpath2 = newpath1+"1"

	# Creating /spectra1 folder

	if os.path.exists(newpath2) == True:
		shutil.rmtree(newpath2)
	if not os.path.exists(newpath2):
		os.makedirs(newpath2)
	
	# Storing the spectra for all the annuli in a matrix "spec"
	
	spec = []
	for i in range(n_annulus):
		temp = np.loadtxt(newpath1+"/"+str(i+1)+".dat",unpack=True)
		spec.append(temp)

	# Transposing the matrix spec to have "count vs annuli"
	
	for i in range(number_of_channels):
		temp = [spec[j][i] for j in range(n_annulus)]
		np.savetxt(str(newpath2)+"/channel"+str(i+1)+".dat",temp)


annuli_to_channel(n_annulus,'spectra')




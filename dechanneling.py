'''
This code takes the deprojected channeled data (count vs radius for each energy range) as input. It makes count vs energy channel files (for each annuli) and saves them in the folder /deprojected1.
'''
import numpy as np
import math
import os
import shutil
import sys

r = np.loadtxt("radius.dat",unpack=True)
n_annulus = len(r)
number_of_channels = 1980

# Bash Script that lists the name of all the .dat files in the folder in a single file.

f = open("script_pha","w")
f.write("#!/bin/bash\n")
f.write("cd deprojected\n")
f.write("ls *.dat >File2\n")
f.write("cp File2 ..\n")
f.close()
os.system("chmod +x script_pha")
os.system("./script_pha")

f = open("File2")
ff = f.readlines()
files = [i.strip('\n') for i in ff]
n_channels = len(files)
os.system("cd ..")

def channel_to_annuli(n_annuli,statement):
	
	# Craeting a new folder /deprojected1
	
	path = os.getcwd()
	newpath1 = path +"/"+str(statement)
	newpath2 = newpath1+"1"
	if os.path.exists(newpath2) == True:
		shutil.rmtree(newpath2)
	if not os.path.exists(newpath2):
		os.makedirs(newpath2)
	
	# Storing the spectra for all the energy channels in a matrix "spec"
	
	spec = []
	for i in range(len(files)):
		temp = np.loadtxt(newpath1+'/'+"channel"+str(i+1)+".dat",unpack=True)
		spec.append(temp)
	
	# Transposing the matrix spec to have "count vs energy channel"
	
	for i in range(n_annulus):
		temp = [spec[j][i] for j in range(n_channels)]
		np.savetxt(str(newpath2)+"/annulus"+str(i+1)+".dat",temp)


channel_to_annuli(n_channels,'deprojected')




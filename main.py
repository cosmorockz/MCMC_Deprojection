'''
Inputs: spectral data, grid.txt
Outputs: density_PROJCT.dat,density_DSDSEPROJ.dat, temp_PROJCT.dat, temp_DSDEPROJ.dat
'''

'''
This code outputs the temperature and density of the spectral data(Using both PROJCT and DSDEPROJ).
The steps for this process are: 
1. This code takes the spectral data as .dat/.txt files. 
2. Using the spectral data, it creates the .pha (by XSpec fakeit command). 
3. Then it converts the .pha files in .dat files.  
4. Then separates the data into energy channles (Channeling).
5. The separate energy channels are then deprojected.
6. We then combine the energy channels (Dechanneling).
7. We fit the data using XSpec. We obtain the temperature and normalization.
8. We convert the normalization to density.
'''

import os
import numpy as np

# Bash Script to execute rename_file. It renames the spectral files in the ascending order of redius.

f = open("script_deproject_and_fitting","w")
f.write('#!/bin/bash\n')
f.write('cd ..\n')
f.write("data=$(ls *.txt)\n")
f.write('cp $data \spectra \n')
f.write('cd \spectra\n')
f.write('gcc rename_file.c -o rename_file\n') # rename_file.c needs grid.txt
f.write('./rename_file\n')
f.close()
os.system("chmod +x script_deproject_and_fitting")
os.system("./script_deproject_and_fitting")

# Bash Script to extract the radius of the annuli from grid.txt. It creates radius.dat file.

f = open("script_deproject_and_fitting2","w")
f.write('#!/bin/bash\n')
f.write('echo "import numpy as np\n')
f.write("r1,r2=np.loadtxt('grid.txt',usecols=(0,1),unpack=True)\n")
f.write("np.savetxt('radius.dat',np.c_[r2])\n")
f.write('exit()" | python')
f.close()
os.system("chmod +x script_deproject_and_fitting2")
os.system("./script_deproject_and_fitting2")

# Bash Script to create data_original.dat. It stores the volume_integrated density and temperature.

f = open("temp","w")
f.write('#!/bin/bash\n')
f.write('cd ..\n')
f.write('data=$(ls data*.tab)\n')
f.write('cp $data \spectra\n')
f.write('cd \spectra\n')
f.write('mv $data data_original.dat\n')
f.close()
os.system("chmod +x temp")
os.system("./temp")


r = np.loadtxt('radius.dat',unpack=True)
n_annulus = len(r)


import pha_file_maker
import pha_to_dat
import channeling
import mock_spectra_deprojection
import dechanneling
import fit
import norm_to_density
import projct_fitting
import norm_to_density_projct

# Bash Script to combine the density and temperature to a single file.

f = open("script_deproject_and_fitting3","w")
f.write('#!/bin/bash\n')
f.write('gcc combine.c -o combine\n') # combine.c needs grid.txt, density_DSDEPROJ.dat and temp_DSDEPROJ.dat
f.write('./combine\n')
f.close()
os.system("chmod +x script_deproject_and_fitting3")
os.system("./script_deproject_and_fitting3")

import plot



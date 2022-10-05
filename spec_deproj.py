'''
This code deprojects the spectra of a given channel.
'''

import xspec
import math
import numpy as np


def spectra_deprojection(fName,r):
	'''
	fName: Name of the channeled spectra
	r: radius of the corresponding annuli
	This function follows the algorithm DSDEPROJ
	Returns the count rate for each shell
	'''
	S = np.loadtxt(fName, unpack =True)
	r = r[::-1]
	number_of_annuli = len(r) - 1
	Vol = []
	S = S[::-1]
	for i in range(number_of_annuli+1):
		Vol.append([])
		for j in range(number_of_annuli+1):
			if j > i:
				temp =(4/3.)*math.pi*(r[i]**2 - r[j]**2)**(3/2.)
			if j > i:
				Vol[i].append(temp)
			else:
				Vol[i].append(0)

	A = []
	for i in range(number_of_annuli):
		temp = math.pi*(r[i]**2-r[i+1]**2)
		A.append(temp)

	C = []
	for m in range(number_of_annuli):
		temp1 = S[m]
		temp2 = 0
		for i in range(m):
			temp2 += C[i]*((Vol[i][m+1]-Vol[i+1][m+1])-(Vol[i][m]-Vol[i+1][m]))
		temp3 = temp1 - temp2
		temp = temp3/((Vol[m][m+1]-Vol[m+1][m+1])-(Vol[m][m]-Vol[m+1][m]))
		C.append(temp)

	C = C[::-1]
	return C









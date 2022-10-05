'''
Plots the deprojected and non-deprojected density and temperature and compares it with the volume integrated density and temperature
'''
import numpy as np
import matplotlib.pyplot as plt

n = np.loadtxt("density_DSDEPROJ.dat",unpack=True)
T = np.loadtxt("temp_DSDEPROJ.dat",unpack=True)
TT = np.loadtxt("temp_PROJCT.dat",unpack=True)
nn = np.loadtxt("density_PROJCT.dat",unpack=True)

r = np.loadtxt("radius.dat",unpack=True)
rr = r.tolist()
rr = rr[::-1]
rr.append(0)
rr = rr[::-1]
r_data = [(rr[i]+rr[i+1])/2.0 for i in range(len(r))]
r_error = [(rr[i+1]-rr[i])/2.0 for i in range(len(r))]
#r1,rho,e = np.loadtxt("de_vs_re0.dat",usecols=(0,1,2),unpack=True)
r1,rho,e = np.loadtxt("data_original.dat",usecols=(0,2,4),unpack=True)
r1 = np.array(r1)*100.0

plt.figure(1)
plt.xlabel("Radius(kpc)")
plt.ylabel("Density")
plt.xscale('log')
plt.yscale('log')
plt.grid()
plt.scatter(r_data,n,label='Deprojected')
plt.scatter(r_data,nn,label='Projected')
ne = rho/(1.151)
plt.plot(r1,ne,label='Simulated Density Profile',color='green')
plt.legend()
plt.show()

plt.figure(2)
plt.xlabel("Radius(kpc)")
plt.ylabel("Temperature")
plt.xscale('log')
plt.grid()
plt.scatter(r_data,T,label='Deprojected')
plt.scatter(r_data,TT,label='Projected')
T = e*0.5987*1.67e-8/(1.38e-16 * np.array(rho))
T = T/1.16e+7
plt.plot(r1,T,label='Simulated Temperature Profile',color='green')
plt.legend()
plt.show()

''' 
Uncomment the commented part when we have the volume integrated density and temperature.
'''




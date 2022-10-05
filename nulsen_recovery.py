import numpy as np
from scipy.integrate import quad
import corner
import matplotlib.pyplot as plt

G = 6.67e-8
kpc = 3.086e+21
k_B = 1.380658e-16
m_H = 1.6733e-24
rho_crit = 8.638e-30 #'''in gm/cm^3 '''
mu = 0.59876
mue = 1.151
keV = 1.16e7
v_c = 3.5e+7

r11, r22, ne, T = np.loadtxt("data.dat",usecols=(0,1,2,3), unpack=True)
ne = np.array(ne)

r1 = r11
r1 = r1.tolist()
r1.append(r22[len(r22)-1])

#ne = rho
T  = np.array(T)*keV
r1 = np.array(r1)*kpc
r_data = (np.array(r11)+np.array(r22) )/2.0*kpc

c_act = 4.7
m200_act = 0.5e+48
v_c_act = 350e+5
rho1_act = 1.27e-01
rii = 0.2*kpc

def g(r,c,m200,r200):
        temp1 = G*m200/(np.log(1+c) - c/(1+c))
        r   = np.sqrt(r*r+rii*rii)
	x = r/r200
	temp2 = - c/(1.+c*x) + np.log(1+c*x)/x
	temp3 = v_c**2/r
        temp = temp1*temp2*x/r**2 + temp3
	return temp

def rho_calculator(rho1,c,m200):
	r200 = (3.0*m200/(4.*np.pi*200.*rho_crit))**(1./3.)
	ne_midpoint = []
	rho_model = []

	for i in range(len(r1)-1):
		rho_model.append(rho1)

		Dphi	= quad(lambda r: g(r,c,m200,r200), r1[i], r_data[i])
		ne_mid  = rho1*np.exp(- Dphi[0]*mu*m_H/(k_B*T[i]))
		ne_midpoint.append( ne_mid )
		
		if (i == len(r1)-2):
#			print "endpoint reached"
			continue
		Dphi	= quad(lambda r: g(r,c,m200,r200), r1[i], r1[i+1])
		ne_ip1  = rho1*np.exp(- Dphi[0]*mu*m_H/(k_B*T[i]))
		rho1	= T[i]*ne_ip1/T[i+1]	

	return np.array(ne_midpoint)

rho_1,c_1,m200_1 = np.loadtxt("chain.dat",usecols=(0,1,2),unpack=True)
T = np.loadtxt("Temperature.dat",unpack=True)
rho1_mean = np.mean(rho_1[2000:])
c1_mean = np.mean(c_1[2000:])
m200_mean = np.mean(m200_1[2000:])
#vc_mean = np.mean(v_c_1[2000:])
rho1_std = np.std(rho_1[2000:])
c1_std = np.std(c_1[2000:])
m200_std = np.std(m200_1[2000:])
#vc_std = np.std(v_c_1[2000:])

values = [rho1_mean,c1_mean,m200_mean]#,vc_mean]
uncertainties = [rho1_std,c1_std,m200_std]#,vc_std]

np.savetxt("Recov_by_MCMC.params",np.c_[values,uncertainties])
rho_model_1 = rho_calculator(rho1_mean,c1_mean,m200_mean)#,vc_mean)
n = mue/mu*np.array(rho_model_1)
p = n*k_B*np.array(T)
np.savetxt("MC_fit_nulsen.out", np.c_[r_data, rho_model_1, p])  # r, ne, p


#rho_model = rho_calculator(rho1_act, c_act, m200_act)
#plt.figure(1)
#plt.xscale('log')
#plt.yscale('log')
#plt.plot(r_data,(ne),label="Given Data",color='green')
#plt.plot(r_data,(rho_model),label="Simulated Data with actual values",color='red')
#plt.plot(r_data,rho_model_1,label="Recovered Data using MCMC",color='blue')
#plt.legend()
#plt.show()

#plt.figure(2)
#plt.xscale('log')
#plt.yscale('log')
#plt.plot(r_data,e,label="Given Data",color='green')
#plt.plot(r_data,p,label="Recovered Data using MCMC",color='blue')
#plt.legend()
#plt.show()



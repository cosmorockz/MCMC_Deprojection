import numpy as np
from emcee import *
from scipy.integrate import quad
import corner

G = 6.67e-8
kpc = 3.086e+21
rho_crit = 8.638e-30
k_B = 1.380658e-16
m_H = 1.6733e-24
mue = 1.157
mu = 0.59876
v_c = 3.5e+7 #in cm/s

r11,r22,rho,T = np.loadtxt("data.dat",usecols=(0,1,2,3),unpack=True)
r1 = r11
r1 = r1.tolist()
r1.append(r22[len(r22)-1])
rho_err = np.array(rho)*1e-2*np.random.normal(0,1)
rho = rho+rho_err
r_data = ((np.array(r11)+np.array(r22))/2.0)
T = np.array(T)*1.16e+7
np.savetxt("Temperature.dat",np.c_[T])

r1 = np.array(r1) * kpc
r_data = np.array(r_data) * kpc

i = 0
j = 0
rii = 0.2*kpc

def g(r,c,m200,r200):
        temp1 = G*m200/(np.log(1+c) - c/(1+c))
        r   = np.sqrt(r*r+rii*rii)
        x = r/r200
        temp2 = - c/(1.+c*x) + np.log(1+c*x)/x
        temp3 = (v_c**2)/r
        temp = temp1*temp2*x/r**2 + temp3
        return temp

def rho_calculator(rho1,c,m200):
        r200 = (3.0*m200/(4.*np.pi*200.*rho_crit))**(1./3.)
        ne_midpoint = []
        rho_model = []

        for i in range(len(r1)-1):
                rho_model.append(rho1)

                Dphi    = quad(lambda r: g(r,c,m200,r200), r1[i], r_data[i])
                ne_mid  = rho1*np.exp(- Dphi[0]*mu*m_H/(k_B*T[i]))
                ne_midpoint.append( ne_mid )

                if (i == len(r1)-2):
                        continue
                Dphi    = quad(lambda r: g(r,c,m200,r200), r1[i], r1[i+1])
                ne_ip1  = rho1*np.exp(- Dphi[0]*mu*m_H/(k_B*T[i]))
                rho1    = T[i]*ne_ip1/T[i+1]

        return np.array(ne_midpoint)


def lnprior(theta):
	#global i
	#i += 1
	#print i
	rho1,c,m200 = theta
	r200 = (m200/(200.*4.*np.pi*rho_crit/3.0))**(1./3.)
	if rho1 <=1e-3 or c <=2.0 or m200 <=5e+46 :
		return -np.inf
	if rho1 >=1.0 or c >=8.0 or m200 >= 5e+48 :
		return -np.inf
	'''if rho1 <=0 or c <=0 or m200 <=0 or v_c <= 0:
		return -np.inf'''
	return 0.0

def lnprob(theta,rho_err,r_data):
	lp = lnprior(theta)
	if not np.isfinite(lp):
		return -np.inf
	return lp + lnlike(theta,rho_err)

def lnlike(theta,rho_err):
	global j
	j += 1

	rho1,c,m200 = theta
	rho_model = rho_calculator(rho1, c, m200)
	inv_sigma2 = 1./rho_err**2.0

	chisq = -0.5*np.sum((rho-rho_model)**2.0 *inv_sigma2 -np.log(inv_sigma2))

	return chisq

ndim,nwalkers = 3,100
posi = [[max(1.010e-3,9.775e-02 + (4e-3*np.random.randn())),max(2.1,4.7 + 1e-1*np.random.randn()), max(2e+46,5e+47 + (2e+45*np.random.randn()))] for i in range(nwalkers)]

sampler = EnsembleSampler(nwalkers,ndim,lnprob,args=(rho_err,r_data),threads=16)
sampler.run_mcmc(posi,1000)
samples = sampler.chain[:,50:,:].reshape((-1, ndim))
np.savetxt("chain.dat",np.c_[samples])
fig = corner.corner(samples)
fig.savefig('params.png')
pp = samples.T

f = open("nulsen-fit.out","w")
f.write("%3.4e %3.4e " % (np.mean(pp[0]), np.std(pp[0])))
f.write("%3.4e %3.4e " % (np.mean(pp[1]), np.std(pp[1])))
f.write("%3.4e %3.4e " % (np.mean(pp[2]), np.std(pp[2])))
f.write("%3.4e %3.4e " % (v_c, 0.0))
f.close()


print "rho1 :",np.mean(pp[0]),"+-",np.std(pp[0])
print "c :",np.mean(pp[1]),"+-",np.std(pp[1])
print "m200 :",np.mean(pp[2]),"+-",np.std(pp[2])
print "v_c :",v_c,"+-",0.0

''' To run this code keep data.dat in (r1 r2 rho T(in KeV) format) in the same folder.
In the chain.dat the columns are rho1 c m200 v_c take their means. Then generate rho using nulsen's method. '''




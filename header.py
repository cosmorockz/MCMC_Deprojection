'''
This file stores all the parameters for the deprojection code.
'''
from astropy.cosmology import WMAP9 as cosmo

pc = 3.08657e+18 # Parsec to Centimeter
kpc = 3.08657e+21 # Kiloparsec to Centimeter
mpc = 3.08657e+24 # Megaparsec to Centimeter
nH = 1e-06 # Neutral Hydrogen Density
Z = 0.2 # Metallicity
zz = 0.001 # Redshift
Da_obj = cosmo.angular_diameter_distance(zz) # Angular Diameter Distance (as object)
Da_Mpc = Da_obj.value # Angular Diameter Distance in Megaparsec
Da = Da_Mpc * mpc # Angular Diameter Distance in Centimeter

# Relation between electron and ionic densities

mu = 0.59876
mu_e = 1.15068
mu_i = 1./((1./mu)-(1./mu_e))


#include "head.h"

#include<stdio.h>
#include<stdlib.h>

int WithinCavity(double, double, double);

/* *************************************************** */
double AvgCoolingTime(double **x, double **y, int i)
{
 double r, ne, ne_avg, T, M, Surf, dS, theta, dtheta, phi, dphi, x1, y1, z1; 

 dtheta	= 0.01;
 dphi	= 0.01;

 r	= (x[i][0]+x[i][1])/2.0/kpc;
 ne	= y[i][0];
 M	= 0.0;
 Surf	= 0.0;

 theta	= 0.0;
 while (theta <= PI){
     phi	= 0.0;
     while (phi <= 2*PI){
         x1	= r*sin(theta)*cos(phi);
         y1	= r*sin(theta)*sin(phi);
         z1	= r*cos(theta);
         dS	= sin(theta)*dtheta*dphi;
         Surf	+= dS; 

         M	+= dS*ne;
         phi 	+= dphi;
     }
     theta	+= dtheta;
 }
 ne_avg	= M/Surf;
 T	= y[i][2];

 return CoolingTime(ne_avg, T);
}


/* ********************************************** */
int WithinCavity(double x, double y, double z)
{
 double r0 = 40.0;

 if ( (x*x+y*y+(z-r0)*(z-r0) <= r0*r0) || (x*x+y*y+(z+r0)*(z+r0) <= r0*r0) ){
     return YES;
 }else{
     return NO;
 }
}

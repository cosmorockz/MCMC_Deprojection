
#include "head.h"

#include<stdio.h>
#include<stdlib.h>

/* *************************************************** */
double FreeFallTime(double r)
{
/*
 static int first_time = 1;
 static double M_200, c, r_200, H_0, vc;
 double x, g_nfw, g_bcg, g;
 
 if (first_time){
     M_200	= 1.4e48*gm;
     c		= 4.7;
     H_0	= 70.*km/sec/Mpc;
     r_200	= pow(3.*M_200/(4.*PI*Delta*rho_crit), 1./3.);
     vc		= 350.*km/sec;
     first_time	= 0;
 }


 x	= r/r_200;
 g_nfw	= G*M_200/f(c) * x/r/r * (-c/(1.+c*x) + log(1.+c*x)/x);
 g_bcg	= vc*vc/r;
 g	= fabs(g_nfw) + g_bcg + gbh(r, 1e10*Msun);
*/

 double g;      
 g 	= fabs(gnfw(r, 7.e14*Msun, 4.7)) + gbcg(r, 3.5e7) + gbh(r, 1e10*Msun);

 return sqrt(2.*r/g);
}

/* ******************************************** */
double CoolingTime(double ne, double T)
{
 static int first_call = 1;
 static double mu_e, mu_i;
 double n, ni, rho_by_mp;

 if (first_call){
     mu_e	= mue;
     mu_i	= 1.248;
     first_call = 0;
 }
 
 rho_by_mp	= ne*mu_e;
 n		= rho_by_mp/mu;
 ni		= rho_by_mp/mu_i;

 return	3./2.*n*kB*T/(ne*ni*Lambda(T));
}

/* *********************************************** */

double Lambda (double T)
{
 int    klo, khi, kmid;
 double   Tmid, dT, tlo, thi, llo, lhi, pllo, plhi, dz, metallicity = 0.2;
 static int  ntab1, ntab3;
 static double  *L_tab1, *L_tab3;
 static double  *T_tab1, *T_tab3; 

 FILE *fcool1, *fcool3 ;

// if (T<1e4 || T> 3e8) return 0.0;

/* -------------------------------------------
        Read tabulated cooling function
   ------------------------------------------- */
 if (T_tab1 == NULL){
     fcool1 = fopen("cooltables/sd_0_tef.dat","r");
     if (fcool1 == NULL){
         printf ("! cooltables/sd_0_tef.dat does not exists\n");
         exit(1);
     }
     L_tab1 = calloc(100, sizeof(double));
     T_tab1 = calloc(100, sizeof(double));

     ntab1 = 0;
     while (fscanf(fcool1, "%lf %lf\n", T_tab1 + ntab1, 
                                      L_tab1 + ntab1)!=EOF) {
     ntab1++;
     }
 }

 if (T_tab3 == NULL){
     fcool3 = fopen("cooltables/sd_m_1_tef.dat","r");
     if (fcool3 == NULL){
         printf ("! cooltables/sd_m_1_tef.dat does not exists\n");
         exit(1);
     }
     L_tab3 = calloc(100, sizeof(double));
     T_tab3 = calloc(100, sizeof(double));

     ntab3 = 0;
     while (fscanf(fcool3, "%lf %lf\n", T_tab3 + ntab3, 
                                      L_tab3 + ntab3)!=EOF) {
     ntab3++;
     }
 }

/* ----------------------------------------------
        Table lookup by binary search  
   ---------------------------------------------- */

 klo = 0;
 khi = ntab1 - 1;
 while (klo != (khi - 1)){
     kmid = (klo + khi)/2;
     Tmid = pow(10.0, T_tab1[kmid]);
     if (T <= Tmid){
         khi = kmid;
     }else if (T > Tmid){
         klo = kmid;
     }
 }
  thi	= pow(10.0 , T_tab1[khi]);	tlo	= pow( 10.0 , T_tab1[klo] );

/* *****************extrapolating the cooling data ******** */
  dz 	= log10( metallicity ) ;
  plhi	= L_tab1[khi] + dz/1.0*(L_tab1[khi] - L_tab3[khi]) ;	
  pllo	= L_tab1[klo] + dz/1.0*(L_tab1[klo] - L_tab3[klo]) ;
  lhi	= pow(10.0 , plhi);		llo = pow(10.0, pllo) ;

  if (T > thi || T < tlo){
    printf (" ! T out of range   T = %12.6e\n",T);
    exit(1);
  }

  dT      = thi - tlo;
  return llo*( thi - T )/dT + lhi*( T - tlo )/dT;
}

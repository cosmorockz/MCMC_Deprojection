#include<stdio.h>
#include<stdlib.h>

#include "codes/head.h"

double dpdr(double, double *);
double g_mc(double, double *);
void GetNLparams(double [], double []);
void GetMCparams(double [], double []);
void GetneFromMCMC(double *);

int main()
{
 cons = pow(3.0/(4.0*PI*Delta*rho_crit), 1./3.);

 int i, stat;
 double **x, **y, **sig;
 double par_nl[MA], dpar_nl[MA], minch_nl, par_mc[3], dpar_mc[3], minch_mc;
 double r, ne, rho, temp, g_nl, tff_nl, tff_mc, tff_act, tc_act, tc_nl, tc_mc;
 double *ne_mc, g_act, gmc;
 char infile[50] = "data.dat", outfile[50] ="out.out";
 FILE *fout = fopen(outfile, "w");

 x      = Array2D(NPT, 2, sizeof(double));
 y      = Array2D(NPT, 3, sizeof(double));
 sig    = Array2D(NPT, 3, sizeof(double));
 ne_mc	= Array1D(NPT, sizeof(double));

 printf("going to read data.dat \n");
 ReadFile(infile, x, y, sig);
 printf("data.dat read\n");

// if( nl_fit(x, y, sig, par_nl, dpar_nl, &minch_nl) != 0){
//     printf("Non-linear fitting could not be done.\n");
// }else{
//     printf("Non-linear fit performed successfully.\n");
// }

 GetNLparams(par_nl, dpar_nl);
 par_nl[1]	*= kpc;
 dpar_nl[1] 	*= kpc;

// if ( system("python nulsen.py") != 0){
//     printf("MC_fit failed.\n");
// }else{
//     printf("Nulsen's method applied succesfully.\n");
// }

 GetMCparams(par_mc, dpar_mc);
 GetneFromMCMC(ne_mc);

 fprintf(fout, "# [1]r, [2]tc_nl/Myr, [3]tc_mc/Myr, [4]tff_nl/Myr, [5]tff_mc/Myr, [6]g_nl, [7]gmc [8]tff_act, [9]g_act\n");
 for (i=0; i<NPT; i++){
     r		= (x[i][0]+x[i][1])/2.0;
     ne		= y[i][0];
     rho	= ne*mue*mp;
     temp	= y[i][2];

//     tc_act	= AvgCoolingTime(x, y, i);
//     g_act	= gnfw(r, 7e14*Msun, 4.7) + gbcg(r, 3.5e7) + gbh(r, 1e10*Msun);
//     tff_act	= FreeFallTime(r);

     g_nl	= -1./rho*dpdr(r, par_nl);
     tff_nl	= sqrt(2.*r/g_nl);
     tc_nl	= CoolingTime(ne, temp);

     gmc	= g_mc(r, par_mc);
     tff_mc	= sqrt(2.*r/gmc);
     tc_mc	= CoolingTime(ne_mc[i], temp);

     fprintf(fout, "%3.4e %3.4e %3.4e %3.4e %3.4e %3.4e %3.4e \n",
                    r, tc_nl/Myr, tc_mc/Myr, tff_nl/Myr, tff_mc/Myr, g_nl, gmc );

 }

 printf("Deprojected points are done.\n");
 CalculateAvgValues("../../","data.0084.dbl");

 fclose(fout); 
 FreeArray2D(NPT, (void **)x);
 FreeArray2D(NPT, (void **)y);
 FreeArray2D(NPT, (void **)sig);

 return 0;
}


/* ************************************************** */
double g_mc(double r, double *par)
{
 static int This_is_first_time = 1;
 static double Mvir, c, vc, r200;
 double x;

 if (This_is_first_time){
     Mvir = par[1];
     c	  = par[2];
     vc	  = 3.5e7; //par[3];
     r200 = cons*pow(Mvir, 1./3.);
     This_is_first_time = 0;
 }
 
 return gnfw(r, Mvir, c) + gbcg(r, vc);
}
/* ******************************************** */
double dpdr (double r, double *par)
{
 double P0, a, a1, a2, ra, denom;
 P0 = par[0]; a = par[1];
 a1 = par[2]; a2 = par[3];

 ra = r/a;
 denom = pow(ra, a1) + pow(ra, a2);

 return -P0/(denom*denom)*(a1/r*pow(ra, a1) + a2/r*pow(ra, a2));
}

/* ******************************************************* */
void GetMCparams(double par[], double dpar[])
{
 double x11, x12, x21, x22, x31, x32, x41, x42;

 FILE *fin = fopen("nulsen-fit.out","r");
 fscanf(fin, "%lf %lf %lf %lf %lf %lf", &x11, &x12, &x21, &x22, &x31, &x32);
 fclose(fin);

 par[0]	= x11; dpar[0] = x12;
 par[1]	= x31; dpar[1] = x32;
 par[2]	= x21; dpar[2] = x22;
// par[3]	= x41; dpar[3] = x42;
 
}

/* ******************************************************* */
void GetNLparams(double par_nl[], double dpar_nl[])
{
 double x1, x2, x3, x4;
 FILE *fp = fopen("nlfit.param","r");
 fscanf(fp, "%lf %lf %lf %lf", &x1, &x2, &x3, &x4);
 par_nl[0] = x1;
 par_nl[1] = x2;
 par_nl[2] = x3;
 par_nl[3] = x4;
 fscanf(fp, "%lf %lf %lf %lf", &x1, &x2, &x3, &x4);
 dpar_nl[0] = x1;
 dpar_nl[1] = x2;
 dpar_nl[2] = x3;
 dpar_nl[3] = x4;

 fclose(fp);
}


/* **************************************************** */
void GetneFromMCMC(double *ne)
{
 if ( system("python nulsen_recovery.py") != 0) printf("recovery of data from  Nulsen's fit failed.\n");

 int i = 0;
 double x1, x2, x3;
 FILE *fp = fopen ("MC_fit_nulsen.out","r");
 
 while (fscanf(fp, "%lf %lf %lf", &x1, &x2, &x3) != EOF){
     ne[i] = x2;
     i++ ;
 }

 fclose(fp);
}

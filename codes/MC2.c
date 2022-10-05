#include "head.h"
#include<stdio.h>
#include<stdlib.h>

#define NP MA
#define rii	(0.2*kpc)

/* ******************************************************************** */
int MC_fit(double **x, double **y, double **sig, double par[], double dpar[], double *minch)
{

 cons = pow(3.0/(4.0*PI*Delta*rho_crit), 1./3.);

 int i, j, k, l, m, Na[NP], iter, maxiter=0;
 double a[NP], al[NP], ah[NP], da[NP], pda[NP], ****chisq;
 double minchisq, bfa[NP], acc;
 char outfile[50] = "MC_fit.out";
 FILE *fout = fopen(outfile,"w");

 al[0] = 1.e-2;		ah[0] = 1.0;		Na[0] = 40;
 al[1] = 1.e47;		ah[1] = 1.e48;		Na[1] = 40;
 al[2] = 2.0;		ah[2] = 8.0;		Na[2] = 20;
 al[3] = 200.*km;	ah[3] = 500.*km;	Na[3] = 20;
 
 for (m=0; m<2; m++) da[m] = (log10(ah[m])-log10(al[m]))/Na[m];
 for (m=2; m<NP; m++) da[m] = (ah[m]-al[m])/Na[m];
 
 chisq	= Array4D(Na[3], Na[2], Na[1], Na[0], sizeof(dtype)); 
 FindMinimum(x, y, sig, a, al, Na, da, chisq, &minchisq, bfa);
 for (iter=0; iter<maxiter; iter++){
     for (m=2; m<NP; m++){
	 al[m] = bfa[m] - da[m];
	 ah[m] = bfa[m] + da[m];
	 da[m] = (4.*da[m])/(Na[m]*1.0);
     }
     for (m=0; m<2; m++){
         al[m]	= bfa[m]*pow(10., -da[m]);
         ah[m]	= bfa[m]*pow(10., da[m]);
         da[m]	= 4.0*da[m]/(Na[m]*1.0);
     }
     FindMinimum(x, y, sig, a, al, Na, da, chisq, &minchisq, bfa);     
 }

 fprintf(fout, "# Minimum Chisq = %3.8e at\n", minchisq); 
 for (m=0; m<2; m++){
     fprintf(fout, "# \ta[%d]=%3.8e(*or/)%3.8e\t(%3.6e%%)\n ", m, bfa[m], pow(10.,da[m]), (pow(10.,da[m])-1.0)*100.0);
     par[m]	= bfa[m];
     dpar[m]	= pow(10.,da[m]);
 }
 for (m=2; m<NP; m++){
     fprintf(fout, "# \ta[%d]=%3.8e+/-%3.8e\t(%3.6e%%)\n ", m, bfa[m], da[m], da[m]*100.0/bfa[m]);
     par[m]     = bfa[m];
     dpar[m]    = da[m];
 }
 *minch = minchisq;

 double r, Ti, *ne;
 ne	= Array1D(NPT, sizeof(double)); 
 Get_ne(x, y, bfa, ne);
 
 for (i=0; i<NPT; i++){
     Ti         = y[i][2];
     r          = (x[i][0]+x[i][1])/2.0;
     fprintf(fout, "%3.6e\t %3.4e\t %3.4e\t %3.4e\n", r, ne[i], ne[i]*mue/mu*kB*Ti, Ti);
 }
 fclose(fout);

 FreeArray1D(ne);
 FreeArray4D(Na[3], Na[2], Na[1], (void ****)chisq);

 return 0;
}

/* **************************************************** */
void Get_ne(double **x, double **y, double bfa[], double *ne)
{
 int i;
 double r, ne_r, ne1, ne_rip1, Ti, Tip1;
 ne1 = bfa[0];
 for (i=0; i<NPT; i++){
     Ti         = y[i][2];
     r          = (x[i][0]+x[i][1])/2.0;
     ne[i]      = model(r, x[i][0], ne1, Ti, bfa);
     if (i == NPT-1) break;
     ne_rip1    = model(x[i][1], x[i][0], ne1, Ti, bfa);
     Tip1       = y[i+1][2];
     ne1    	= ne_rip1*Ti/Tip1;
 }

}

/* ************************************************************** */
void FindMinimum(double **x, double **y, double **sig, double a[], double al[], 
                 int Na[], double da[], double ****chisq, double *minchisq, double bfa[])
{
 int i, j, k, l;
 double minch = *minchisq, cs;

 for (l=0; l<Na[3]; l++){
     a[3]	= al[3] + l*da[3];
     for (k=0; k<Na[2]; k++){
         a[2]	= al[2] + k*da[2];
 	 for (j=0; j<Na[1]; j++){
             a[1]	= al[1]*pow(10., j*da[1]);
 	     for (i=0; i<Na[0]; i++){
	         a[0]	= al[0]* pow(10., i*da[0]);
    		 chisq[l][k][j][i]	=  GetChisq(x, y, sig, a);
		 FindMinChisq(chisq[l][k][j][i], &minch, a, da, bfa);
//		 printf("vc=%3.4e, c=%3.4e, Mvir=%3.4e, ne0=%3.4e, chisq=%3.4e\n", a[3], a[2], a[1], a[0],chisq[l][k][j][i]);
 	      }
	  }

     }

 printf("> task done : %3.2f\n", l*100.0/(Na[3]*1.0));
 }
 *minchisq = minch;

}

/* *****************************/
void FindMinChisq(double newchisq, double *minchisq, double a[], double da[], double bfa[])
{
 int m;
 static int first_time = 1;
 double ar2, dar2, ar, dar;
 if (first_time){
     *minchisq = newchisq;
     for (m=0; m<NP; m++) bfa[m] = a[m];
     first_time = 0;
     return;
 }
 
 if (newchisq < *minchisq){
     *minchisq = newchisq;
/*     ar2 = 0.0;
     for (m=0; m<NP; m++){
 	 ar2 += (a[m]-bfa[m])*(a[m]-bfa[m]);
	 dar2 += da[m]*da[m];
     }
     ar = sqrt(ar2); dar = sqrt(dar2);
     if (ar > 16.0*dar ){
	 printf("Going to a second far away peak at a distance from the current centre %3.8e\n", ar);
         for (m=0; m<NP; m++) printf("a[%d] = %3.4f, ", m, a[m]);
	 printf("\n");
     }
*/
     for (m=0; m<NP; m++) bfa[m] = a[m];
 }

}

/* ************************************************* */
double GetChisq(double **x, double **y, double **sigma, double *a)
{
 int i;
 double cs, r, ne_i, ne_r, ne_rip1, Ti, Tip1;

 cs = 0.0;
 ne_i = a[0];
 for (i=0; i<NPT; i++){
     Ti		= y[i][2];
     r	 	= (x[i][0]+x[i][1])/2.0;
     ne_r 	= model(r, x[i][0], ne_i, Ti, a);
     cs 	+= pow( (y[i][0]- ne_r)/sigma[i][0], 2.0);
     if (i == NPT-1) break;
     ne_rip1	= model(x[i][1], x[i][0], ne_i, Ti, a);
     Tip1	= y[i+1][2];
     ne_i	= ne_rip1*Ti/Tip1;
 }
 return cs;
}
/* ************************************************* */
double model(double r, double ri, double nei, double Ti, double *a)
{
 double Mvir, c, vc;
 Mvir = a[1]; c = a[2]; vc = a[3];

 return nei*exp(-KELVIN/Ti*DPhi(r, ri, Mvir, c, vc));
}

/* ********************************************** */
double DPhi(double r, double ri, double Mvir, double c, double vc)
{
 double rd = ri, dr = (r-ri)/100.0, ret = 0.0;
 while(rd <= r){
  ret += ( gnfw(rd, Mvir, c) + gbcg(rd, vc) )*dr;
  rd	+= dr;
 }
 return ret;

}

/* *********************************************** */
double gnfw(double r, double Mvir, double c)
{
 double x, r200;

 r200	= cons*pow(Mvir, 1./3.);
 r	= sqrt(r*r+rii*rii);
 x 	= r/r200;

 return G*Mvir/f(c)*x/(r*r)*(log(1.+c*x)/x - c/(1.+c*x));
}
double gbcg(double r, double vc)
{
 r = sqrt(r*r+rii*rii);
 return vc*vc/r;
}
double gbh (double r, double Mbh)
{
 r = sqrt(r*r+0.02*0.02*kpc*kpc);
 return G*Mbh/(r*r);
}

/* ******************************************* */
double f(double c)
{
 return log(1.+c)-c/(1.+c);
}

/* *********************************************** */
void ReadFile(char *file, double **x, double **y, double **sigma)
{
 int i;
 double x1, x2, x3, x4;
 FILE *fin = fopen(file, "r");
 i = 0;
 printf ("reading data.dat \n");
 while (fscanf(fin, "%lf %lf %lf %lf", &x1, &x2, &x3, &x4) != EOF){
     x[i][0]	= x1*kpc;
     x[i][1]	= x2*kpc;
     y[i][0]	= x3;				// ne
     y[i][2]	= x4*keV;				// T in K
     y[i][1]	= y[i][0]*mue/mu*kB*y[i][2];	// Prs = ne*mue/mu*kB*T
     sigma[i][0]= 0.01*y[i][0]*rand()/(RAND_MAX*1.0);
     sigma[i][1]= 0.01*y[i][1]*rand()/(RAND_MAX*1.0);
     sigma[i][2]= 0.01*y[i][2]*rand()/(RAND_MAX*1.0);
//     printf("%3.4e %3.4e %3.4e %3.4e %3.4e\n", x[i][0], x[i][1], y[i][0], y[i][1], y[i][2]);
     i++;
 }
 

 fclose(fin);
}

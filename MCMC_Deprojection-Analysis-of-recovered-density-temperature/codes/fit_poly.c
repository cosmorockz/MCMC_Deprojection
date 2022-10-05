#include <stdio.h>
#include <math.h>
#define NRANSI
#include "head.h"
#include "nr.h"
#include "nrutil.h"


/* *********************************************************************** */
int nl_fit(double **x1, double **x2, double **sigma, double param[], double d_param[], double *minch)
{
 long idum=(-911);
 int i,*ia,iter,itst,j,k,mfit=MA, count;
 float alamda,chisq,ochisq,*x,*y,*sig,**covar,**alpha;
 static float a[MA+1];
 static float gues[MA+1]={0.0, 1.e-9, 20., 0.5, 1.0};
 float *dumm, val;

 ia	= ivector(1,MA);
 x	= vector(1,NPT);
 y	= vector(1,NPT);
 sig	= vector(1,NPT);
 covar	= matrix(1,MA,1,MA);
 alpha	= matrix(1,MA,1,MA);
 dumm	= vector(1, NPT);

 FILE *fout	= fopen("nlfit.out", "w");
 for (i=1; i<= NPT; i++){
    x[i]	= (float )((x1[i-1][0]+x1[i-1][1])/2.0/kpc);		// reading in kpc
    y[i]	= (float )x2[i-1][1];			// taking pressure
    sig[i]	= (float )sigma[i-1][1];
 }

 for (i=1;i<=mfit;i++) ia[i]=1;
 for (i=1;i<=MA;i++) a[i]=gues[i];
 for (iter=1;iter<2;iter++){
     alamda = -1;
     mrqmin(x,y,sig,NPT,a,ia,MA,covar,alpha,&chisq,fbrpl,&alamda);
     k=1;
     itst=0;
     for (;;) {
         fprintf(fout, "# %s %2d %17s %10.4e %10s %9.2e\n","Iteration #",k,
		"chi-squared:",chisq,"alamda:",alamda);
	 k++;
	 ochisq=chisq;
	 mrqmin(x,y,sig,NPT,a,ia,MA,covar,alpha,&chisq,fbrpl,&alamda);
	 if (chisq > ochisq)
	     itst=0;
	 else if (fabs(ochisq-chisq) < 0.1)
	     itst++;
	 if (itst < 10) continue;
	 alamda=0.0;
	 mrqmin(x,y,sig,NPT,a,ia,MA,covar,alpha,&chisq,fbrpl,&alamda);
	 break;
     }
 }

 *minch = (double)chisq;

 fprintf(fout,"# Best fit parameters\n");
 for (i=1;i<=MA;i++){
     fprintf(fout, "# a[%d] = %6.6e +- %6.6e (%3.2f%%)\n", i, a[i], sqrt(covar[i][i]), sqrt(covar[i][i])/a[i]*100.0);
     param[i-1] = (double)(a[i]);
     d_param[i-1] = (double)(sqrt(covar[i][i]));
 }
 fprintf(fout, "\n");

 for (i=1; i<= NPT; i++){
     fbrpl(x[i], a, &val, dumm, MA);
     fprintf(fout, "%6.6f \t %6.6e \t %3.3e\n", x[i], y[i], val);
 }
 fclose(fout);

 free_matrix(alpha,1,MA,1,MA);
 free_matrix(covar,1,MA,1,MA);
 free_vector(sig,1,NPT);
 free_vector(y,1,NPT);
 free_vector(x,1,NPT);
 free_ivector(ia,1,MA);

 return 0;
}

/* ******************************************************** */
void poly(float x, float a[], float *y, float dyda[], int na)
{
 int i;
 *y = a[5]*x*x*x*x + a[4]*x*x*x + a[3]*x*x + a[2]*x + a[1];
 dyda[1] = 1.0;
 dyda[2] = x;
 dyda[3] = x*x;
 dyda[4] = x*x*x;
 dyda[5] = x*x*x*x;
}

/* ******************************************************** */
void fbrpl(float x, float a[], float *y, float dyda[], int na)
{
 float denom, denom2, P0, ab, a1, a2;
 P0 = a[1];
 ab = a[2];
 a1 = a[3];
 a2 = a[4];

 denom = pow(x/ab, a1) + pow(x/ab, a2);
 denom2 = denom*denom;
 *y = (P0/denom);
 
 dyda[1] =  1./denom;
 dyda[2] =  P0/denom2*( a1*pow(x, a1)/pow(ab, a1+1.) + a2*pow(x, a2)/pow(ab, a2+1.) ); 
 dyda[3] = -P0/denom2*( pow(x/ab, a1)*log(x/ab) );
 dyda[4] = -P0/denom2*( pow(x/ab, a2)*log(x/ab) );
// printf("x=%3.4e\tP0=%3.4e\tab=%3.4e\ta1=%3.4e\ta2=%3.4e\n", x, P0, ab, a1, a2);
// printf("y=%3.4e\tdyda[1]=%3.4e\tdyda[2]=%3.4e\tdyda[3]=%3.4e\tdyda[4]=%3.4e\n", *y, dyda[1], dyda[2], dyda[3], dyda[4]);

}

#undef NRANSI

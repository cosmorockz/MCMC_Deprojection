
#include "head.h"

#include<stdio.h>
#include<stdlib.h>

int WithinCavity(double, double, double);

/* *************************************************** */
void CalculateAllValues()
/*!
 * Calculates the tcool and tff at each grid point but saves them only
 * as a function of r. The reason is to calculate all the tcool and tff 
 * at any particular r. This is not neccessary in 1D case of 3D-without 
 * perturbation case. 
 * ************************************************************** */
{
 int nx=256, ny=128, nz=32, nvar=6;
 int i, j, k, cnt;
 double x1, x2, x3, mui;
 double r, ne, n, n_act, T, M, Prs, Surf, dS, theta, dtheta, phi, dphi, n_avg, x, y, z; 
 double *d, *r1, *r2, *th1, *th2, *phi1, *phi2;
 double Tcool, Tff; 
 FILE *fout = fopen("avg_values.txt","w");
 mui = 1./(1./mu - 1./mue);

 d	= calloc (nx*ny*nz*nvar, sizeof(double));
 r1	= calloc (nx, sizeof(double));
 r2	= calloc (nx, sizeof(double));
 th1	= calloc (ny, sizeof(double));
 th2	= calloc (ny, sizeof(double));
 phi1	= calloc (nz, sizeof(double));
 phi2	= calloc (nz, sizeof(double));

 FILE *fdata = fopen("../../data.0000.dbl","rb");
 fseek(fdata, SEEK_SET, 0);
 fread(d, sizeof(double), nx*ny*nz*nvar, fdata);
 fclose(fdata);
 printf("data file read.\n");

 FILE *fgrid = fopen("../../grid1.out","r");
 cnt	= 0;
 while (fscanf(fgrid,"%lf %lf", &x1, &x2) != EOF){
     if (cnt < nx){
         r1[cnt]	= x1*100.0;
         r2[cnt]	= x2*100.0;
     }else if(cnt < nx+ny){
         th1[cnt-nx]	= x1;
         th2[cnt-nx]	= x2;
     }else{
         phi1[cnt-nx-ny]	= x1;
         phi2[cnt-nx-ny]	= x2;
     }
     cnt	+= 1;
 }
 fclose(fgrid);
 printf("grid file read.\n");

 for (i=0; i<nx; i++){
     r		= (r1[i]+r2[i])/2.0;
     for (j=0; j<ny; j++){
	 for (k=0; k<nz; k++){
	     n_act	= d[i+j*nx+k*ny*nx]/mu;
     	     Prs	= d[i+4*nx*ny*nz]*1.673e-8;
             T		= Prs/(n_act*kB);
             ne		= n_act*mu/mue;

             Tcool	= CoolingTime(ne, T);
             Tff	= FreeFallTime(r*kpc);
             fprintf(fout,"%2.4f %3.4e %3.4e \n", r, Tcool/Myr, Tff/Myr);
	 }
     }
//     fprintf(fout,"\n");
 }

 fclose(fout);
 free(d); free(r1); free(r2); free(th1); free(th2); free(phi1); free(phi2);

}


/* *************************************************** */
void CalculateAvgValues(char *datafolder, char *datafilename)
/* ***************************************************
 * Calculates angle averaged quantities for a single r
 * ***************************************************** */
{
 printf("calculating avg.\n");
 int nx=256, ny=128, nz=32, nvar=6;
 int i, j, k, cnt;
 double x1, x2, x3, x, y, z, mui;
 double r, ne, ni, n, T, theta, dtheta, phi, dphi, n_avg, T_avg, dV, vol; 
 double *d, *r1, *r2, *th1, *th2, *phi1, *phi2, dr;
 double Tcool, Tff, emis, prs; 
 char datafile[128], gridfile[128];
 FILE *fout = fopen("avg_values.txt","w");
 mui = 1./(1./mu-1./mue);

 d	= calloc (nx*ny*nz*nvar, sizeof(double));
 r1	= calloc (nx, sizeof(double));
 r2	= calloc (nx, sizeof(double));
 th1	= calloc (ny, sizeof(double));
 th2	= calloc (ny, sizeof(double));
 phi1	= calloc (nz, sizeof(double));
 phi2	= calloc (nz, sizeof(double));

 sprintf (datafile,"%s/%s",datafolder,datafilename);
 FILE *fdata = fopen(datafile,"rb");
 fseek(fdata, SEEK_SET, 0);
 fread(d, sizeof(double), nx*ny*nz*nvar, fdata);
 fclose(fdata);

 sprintf(gridfile,"%s/grid1.out",datafolder);
 FILE *fgrid = fopen(gridfile,"r");
 cnt	= 0;
 while (fscanf(fgrid,"%lf %lf", &x1, &x2) != EOF){
     if (cnt < nx){
         r1[cnt]	= x1;
         r2[cnt]	= x2;
     }else if(cnt < nx+ny){
         th1[cnt-nx]	= x1;
         th2[cnt-nx]	= x2;
     }else{
         phi1[cnt-nx-ny]	= x1;
         phi2[cnt-nx-ny]	= x2;
     }
     cnt	+= 1;
 }
 fclose(fgrid);

 for (i=0; i<nx; i++){
     r		= (r1[i]+r2[i])/2.0;
     dr		= r2[i]-r1[i];

     n_avg	= 0.0;
     T_avg	= 0.0;
     vol	= 0.0;

     for (j=0; j<ny; j++){
	 theta	= (th1[j]+th2[j])/2.0;
	 dtheta	= th2[j]-th1[j];
	 for (k=0; k<nz; k++){
//	     phi	= (phi1[k]+phi2[k])/2.0;
	     dphi	= phi2[k]-phi1[k];
	     dV	= r*r*dr*sin(theta)*dtheta*dphi;

	     n		= d[i+j*nx+k*ny*nx]/(mu*mp);
	     prs	= d[i+j*nx+k*ny*nx + 4*nx*ny*nz];
	     T		= prs/(n*kB);
  	     if (T <= 1.e4 || T > 3.e8) continue;

	     ne		= n*mu/mue;
	     ni		= n*mu/mui;
//	     if (T< 5.8e6) emis = 0.0;
	     if (T< 5.8e6) emis = 0.0;
 	     else        emis = ne*ni*Lambda(T);

	     vol	+= dV*emis; 
	     n_avg	+= n*emis*dV;
	     T_avg 	+= T*emis*dV;
	 }
     }
     n_avg	= n_avg/vol;
     T_avg	= T_avg/vol;
     ne		= n_avg*mu/mue;
     Tcool	= CoolingTime(ne, T_avg);
     Tff	= FreeFallTime(r);
     fprintf(fout,"%2.4f %3.4e %3.4e %3.4e %3.4e\n", r/kpc, Tcool/Myr, Tff/Myr, ne, T_avg);
 }

 fclose(fout);
 free(d); free(r1); free(r2); free(th1); free(th2); free(phi1); free(phi2);
 printf("Avg calculated.\n");
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

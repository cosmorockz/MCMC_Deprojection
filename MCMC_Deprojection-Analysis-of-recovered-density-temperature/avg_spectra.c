#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include <sys/types.h>
#include <sys/stat.h>

#define dH	(4285.*1e3)	// size of the horizon [kpc]; assuming h=0.7 [standard]

#define dE	0.008		// keV				[standard - DO NOT CHANGE]
#define nEbins	1980		// number of bins in the spectra [ standard - DO NOT CHANGE]
#define Ebeg	0.106		// keV				[standard - DO NOT CHANGE]
#define dcl	500.0		// size of the cluster [kpc]
#define min_count1	2e3	// 5e7
#define min_count2	2e3	//5e7
#define Aeff	100.0		// cm^2; effective area of the instrument
#define t_exp	1e5		// sec; exposure time

#define	nx	128
#define ny	128

float z, dA;			// z=redshift, dA=angular diameter distance
char *spec_dir;
	
void SetTotPhotonsToZero(float *);
void ReadEnergy(char [], float *);
void ReadSpectra(float *, char []);
void FindTotalSpectra(float *, float *);
void FindAvgSpectra(float *, int, float *);
void FindTotalCounts(float *, float *, float, float, float*);
void PrintSpectra(float, float, float *, float *);

int main(int argc, char *argv[])
/*!
 * Defines the annuli so that each annuli contains a minimum 'min_count' number of counnts in them
 * ans finds out averaged spectrum for that annulus.
 * Input spectral files are taken from pass_1.3/ (located at ~/plutoworks/Cluster/)
 * ***************************************************************************** */
{
 if (argc < 3){
     printf("> Rule is ...$ ./avg_spectra z spec_dir\n> Exiting.\n");
     exit(1);
 }
 z	= atof(argv[1]);
 spec_dir = argv[2]; 
 
 int i, j, k, s, m, n;
 long int file_count;
 float total_count, min_count, r2, r1, r1p, dm, dn, r_m, r_n, dx_initial;
 float r, dr, x, y;
 float *E, *photons, *tot_photons, *avg_spectra;
 char specfile[100];
 FILE *fgrid = fopen("grid.txt", "w");
 FILE *flog  = fopen("logfile","w");

 E		= calloc(nEbins, sizeof(float));
 photons	= calloc(nEbins, sizeof(float));
 tot_photons	= calloc(nEbins, sizeof(float));
 avg_spectra	= calloc(nEbins, sizeof(float));


 dA	= z*dH;			// angular diameter distance
 dm	= (dcl/nx);		// physical size each grid represents
 dn	= (dcl/ny);
 dx_initial	= 0.5;
 dr	= dx_initial;		// differential annuli
 r_m	= 1.02613831;
 r_n	= r_m;
 r1	= 0.0;			// starting annulus
 r2	= 2.*dr ;


 while(r2 < dcl){
     SetTotPhotonsToZero(tot_photons);
     SetTotPhotonsToZero(avg_spectra);
     file_count	= 0;
     total_count = 0.0;
     r1p	= r1;

     min_count	= (r2 < 20.)? min_count1 : min_count2 ; 
     while (total_count < min_count){
	 for (m=-nx+1; m<nx; m++){	// finds the spectral files and adds counts for the annulus r1p-r2
	 for (n=-ny+1; n<ny; n++){
             x	= (m<0)? -dx_initial*(1.-pow(r_m, abs(m)))/(1.-r_m) : dx_initial*(1.-pow(r_m, m+1))/(1.-r_m);
             y	= (n<0)? -dx_initial*(1.-pow(r_n, abs(n)))/(1.-r_n) : dx_initial*(1.-pow(r_n, n+1))/(1.-r_n);
     	     r	= sqrt(x*x+y*y);
	     if ((r < r1p) || (r >= r2)) continue;
	     if (r2 >= dcl){
//		 fprintf(flog, "!Warning: r2 has exceeded cluster radius.\n exiting now.\n"); 
		 break;
	     }
 	     sprintf(specfile,"%s/x%d_y%d.mdl", spec_dir, m, n);
	     ReadEnergy(spec_dir, E);
 	     ReadSpectra(photons, specfile);
	     file_count	+= 1;
	     FindTotalSpectra(photons, tot_photons);
	 }}

	 FindAvgSpectra(tot_photons, file_count, avg_spectra);
 	 FindTotalCounts(E, avg_spectra, r1, r2, &total_count);
	 fprintf(flog, "file_count=%ld, total_count=%2.4e\n", file_count, total_count);		// for test. goes into logfile
	 if (r2 >= dcl){
//	     fprintf(flog, "!Warning: r2 has exceeded cluster radius.\n exiting now.\n"); 
	     break;
	 }

	 r1p	= r2;
	 r2	+= dr;
     }
     fprintf(flog, "%2.4f %2.4f %2.4e\n", r1, r2, total_count);		// for test purpose. goes into logfile
     fprintf(fgrid, "%2.4f %2.4f %2.4e\n", r1, r2, total_count);
     PrintSpectra(r1, r2, E, avg_spectra);

     r1	= r2;
     r2	= r1+dr;
 }

 fclose(fgrid);
 fclose(flog);
 free(E); free(photons); free(tot_photons); free(avg_spectra);
}

/* =================================================================
 *       S U B R O U T I N E S
 * ================================================================= */

/* ****************************************************** */
void ReadEnergy(char spec_dir[], float *E)
{
 int i;
 char eng_file[128];
 static int first_time = 1;

 if (first_time){
     sprintf(eng_file, "%s/energy.mdl", spec_dir);
     FILE *feng = fopen(eng_file, "rb");
     if (feng == NULL){
	 printf("> %s does not exist! Check spectra directory.\n Exiting now.\n", eng_file);
	 exit(1);
     }     

     fseek(feng, 0, SEEK_SET);
     fread(E, sizeof(float), nEbins, feng);
     for (i=0; i<nEbins; i++) E[i]	= E[i]/(1.0+z);
     fclose(feng);
     first_time = 0;
 }


}

/* **************************************************** */
void ReadSpectra(float *photons, char specfile[])
/*!
 * Reads the energy spectra from a given binary spectra file.
 * E [out]		= Energy in keV
 * photons [out]	= photons/s/cm^2/Sr/keV
 * specfile [In]	= spectral file [binary] to be read.
 * ********************************************************** */
{
 int i;
 FILE *fspectra = fopen(specfile ,"rb");
 if (fspectra == NULL){
     printf("> specfile %s does not exist! Check the file path carefully.\n Exiting now.\n", specfile);
     exit(1);
 }

 fseek(fspectra, 0, SEEK_SET);
 fread(photons, sizeof(float), nEbins, fspectra);
 fclose(fspectra);


}

/* ************************************************************* */
void FindTotalCounts(float *E, float *photons, float r1, float r2, float *total_count)
/*!
 * Finds out total counts for the current annulus r1 - r2.
 * total_count considers the Solid angle (dOmega) of the annulus, the effective area
 * of the instrument and the exposure time for the observation.
 * 
 * E [I]	= energy array [keV]
 * photons[I]	= spectra [photons/s/cm^2/keV/Sr]
 * r1, r2 [I]	= inner and outer radius of the annulus [kpc]
 * total_count[O] = total count
 * *********************************************************************************** */
{
 int i;
 float dOmega = 3.14*(r2*r2-r1*r1)/(dA*dA);

 *total_count = 0.0;
 for (i=0; i<nEbins; i++){
     *total_count	+= photons[i]*dE*dOmega*Aeff*t_exp;
     if (E[i] > 10.0) break;
 }

}

/* ********************************************************* */
void FindTotalSpectra(float *photons, float *tot_photons)
{
 int i;
 for (i=0; i<nEbins; i++)
     tot_photons[i]	+= photons[i];
}

/* **************************************************** */
void SetTotPhotonsToZero(float *A)
/*!
 * Initialises the arrays to zero.
 * ******************************** */
{
 int i;
 for (i=0; i<nEbins; i++) A[i] = 0.0;
}

/* ****************************************************** */
void FindAvgSpectra(float *tot_photons, int file_counts, float *avg_spectra)
{
 int i;
 for (i=0; i<nEbins; i++) avg_spectra[i] = tot_photons[i]/file_counts;
}

/* ************************************************************************* */
void PrintSpectra(float r1, float r2, float *E, float *avg_spectra)
/*!
 * Prints the averaged spectra [photons/s/cm^2/keV] of a particular annulus to ascii files
 * 
 * r1, r2 [I]	= inner and outer radius of the annulus.
 * E [I]	= Energy
 * avg_spectra [I] = averaged spectra [photons/s/cm^2/keV/Sr]
 * 
 * ********************************************************************** */
{
 int i; 
 char outspecfile[100];
 float dOmega = 3.14*(r2*r2-r1*r1)/(dA*dA);
 int errdir = mkdir("spectra", 0770);				// creating directory if it does not exist

 sprintf(outspecfile,"spectra/S_%2.1f_%2.1f.spec", r1, r2);
 FILE *fout = fopen(outspecfile, "w");

 for (i=0; i<nEbins; i++) fprintf(fout, "%2.4f  %2.4e\n", E[i], avg_spectra[i]*dOmega);

 fclose(fout);
}

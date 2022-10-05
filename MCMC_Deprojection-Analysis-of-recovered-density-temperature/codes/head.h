#ifndef HEAD_H
#define HEAD_H

//#include<stdio.h>
//#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<time.h>

#define YES	1
#define NO	0

#define dtype   double
#define Delta	200.0
#define G 	6.67259e-8
#define kB      1.380658e-16 //1.38065e-16
#define PI 	3.14159265359
#define pc      3.086e18 //3.0857e18
#define kpc     (1.e3*pc)
#define Mpc	(1.e3*kpc)
#define gm	1.0
#define Msun    1.989e33
#define mp      1.6733e-24
#define km      1.e5
#define sec	1.0
#define yr      3.1536e7
#define Myr     (1.e6*yr)
#define keV	1.16e7

#define rho_crit	8.8638e-30
#define mu		0.59876
#define mue		1.151
#define KELVIN		(mu*mp/kB)

#define NPT 74
#define MA 4

double cons;

dtype     *Array1D (int, size_t);
dtype    **Array2D (int, int, size_t);
dtype   ***Array3D (int, int, int, size_t);
dtype  ****Array4D (int, int, int, int, size_t);
dtype *****Array5D (int, int, int, int, int, size_t);
void FreeArray1D (void *);
void FreeArray2D (int, void **);
void FreeArray3D (int, int, void ***);
void FreeArray4D (int, int, int, void ****);
void FreeArray5D (int, int, int, int, void *****);
dtype **NullArray2D (dtype, dtype, int);

void ReadFile(char *, double **, double **, double **);
double GetChisq(double **, double **, double **, double *);
double model(double, double, double, double, double *);
void FindMinChisq(double, double *, double [], double [], double []);
void FindMinimum(double **, double **, double **, double [], double [],
                 int [], double [], double ****, double *, double []);
int MC_fit(double **, double **, double **, double [], double [], double *);
void Get_ne(double **, double **, double [], double *);

double DPhi(double, double, double, double, double);
double gnfw(double, double, double);
double gbcg(double, double);
double gbh (double, double);
double f(double);

void poly(float, float[], float *, float[], int);
void fbrpl(float, float[], float *, float[], int);
int nl_fit(double **, double **, double **, double *, double *, double *);

double CoolingTime(double, double);
double AvgCoolingTime(double **, double **, int);
void CalculateAvgValues(char *, char*);
void CalculateAllValues();
double FreeFallTime(double);
double Lambda (double);


#endif

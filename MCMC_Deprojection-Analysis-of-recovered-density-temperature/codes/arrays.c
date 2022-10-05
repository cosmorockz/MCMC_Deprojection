#include "head.h"
#include<stdlib.h>
#include<stdio.h>


dtype *Array1D( int NX , size_t dsize )
/* 1D array definition */ 
{
 int i;
 dtype *v;

 v = malloc(NX*dsize);
 
 return v;
}


/* ****************************************************** */
dtype **Array2D ( int NY, int NX, size_t dsize )
/* 2D array definition */
{
 int j;
 dtype **v;
 v = malloc(NY*sizeof(dtype *));

 for (j=0; j<NY; j++) 
     v[j] = Array1D(NX, dsize);

 return v;
}

/* ******************************************************** */
dtype ***Array3D (int NZ, int NY, int NX, size_t dsize)
/* 3D array definition */
{
 int k;
 dtype ***v;
 v = malloc(NZ*sizeof(dtype **));

 for (k=0; k<NZ; k++)
      v[k] = Array2D(NY, NX, dsize);

 return v;
}


/* ********************************************************* */
dtype ****Array4D (int PX, int NZ, int NY, int NX, size_t dsize)
/* 4D array definition */
{
 int l;
 dtype ****v;
 v = malloc(PX*sizeof(dtype ***));

 for (l=0; l<PX; l++)
      v[l] = Array3D(NZ, NY, NX, dsize);

 return v;
}

/* ********************************************************* */
dtype *****Array5D (int PY, int PX, int NZ, int NY, int NX, size_t dsize)
/* 5D array definition */
{
 int l;
 dtype *****v;
 v = malloc(PY*sizeof(dtype ****));

 for (l=0; l<PY; l++)
      v[l] = Array4D(PX, NZ, NY, NX, dsize);

 return v;
}

dtype **NullArray2D (dtype NY, dtype NX, int dsize)
{
 int i, j, k;
 dtype **v;
 v	= Array2D(NY, NX, dsize);
 for (j=0; j<NY; j++ )
 for (i=0; i<NX; i++ )
     v[j][i]	= 0.0;

 return v;

}


/* ********************************************************** *
 * Functions to free higher dimensional arrays 
 * ********************************************************** */
void FreeArray1D (void *v)
{
 free(v);
}

/* ********************************************************* */
void FreeArray2D (int NY, void **v)
{
 int i;
 for (i=0; i<NY; i++) free( v[i] );
 free(v);
}

/* ******************************************************* */
void FreeArray3D (int NZ, int NY, void ***v)
{
 int i;
 for(i=0; i<NZ; i++) FreeArray2D (NY, (void **)v[i]);
 free(v);
}


/* ********************************************************* */
void FreeArray4D (int PX, int NZ, int NY, void ****v)
{
 int i;
 for (i=0; i<PX; i++) FreeArray3D (NZ, NY, (void ***)v[i]);
 free(v);
}
/* ********************************************************* */
void FreeArray5D (int PY, int PX, int NZ, int NY, void *****v)
{
 int i;
 for (i=0; i<PY; i++) FreeArray4D (PX, NZ, NY, (void ****)v[i]);
 free(v);
}

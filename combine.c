// Outputs a single file with radius, density and temperture. The output data.dat can be used for the next steps of the analysis

#include<stdio.h>
#include<stdlib.h>

int main()
{
 double x1, x2, x3;
 FILE *fgrid = fopen("grid.txt","r");
 FILE *fden = fopen("density_DSDEPROJ.dat","r");
 FILE *ftemp = fopen("temp_DSDEPROJ.dat","r");
 FILE *fout = fopen("./data.dat","w");
 
 while( fscanf(fgrid, "%lf %lf %lf", &x1, &x2, &x3) != EOF){
     fprintf(fout, "%3.3e\t %3.3e\t ", x1, x2);
     fscanf(fden, "%lf", &x1);
     fscanf(ftemp, "%lf", &x2);
     fprintf(fout, "%3.3e\t %3.3e \n", x1, x2);
 }

 fclose(fgrid); fclose(fden); fclose(ftemp); fclose(fout);
 return 0;
}

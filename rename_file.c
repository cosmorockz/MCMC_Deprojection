// Renames the .spec files to n.dat where n is the number of annuli in ascending order
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include <sys/types.h>
#include <sys/stat.h>

int main()
{
 int count;
 double x1, x2, x3;
 char infile[50], outfile[50];

 FILE *fgrid = fopen("grid.txt", "r");
 FILE *fin;

 count = 1;
 while(fscanf(fgrid, "%lf %lf %lf",&x1, &x2, &x3) != EOF ){
     sprintf(infile, "S_%3.1f_%3.1f.spec", x1, x2);
     //sprintf(infile,"%d.dat", count);
     sprintf(outfile,"%d.spec", count);
     rename(infile,outfile);
//     printf("filename: %s\n", infile);
//     fin = fopen(infile,"r");
//     if (fin == NULL) printf("file %s could not be found.\n", infile);
     count++;
 }


}

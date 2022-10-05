#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main()
{
 int i;
 double x, dx, y, a[5], sigma;
 a[0] = 1.0;
 a[1] = 0.5;
 a[2] = 0.25;
 a[3] = 0.1;
 a[4] = 0.01;
 x = 1.0;
 dx = 0.1;
 srand(5.0);
 FILE *fdata = fopen("datafile.txt","w");
 while (x < 10.0){
     y		= a[4]*x*x*x*x + a[3]*x*x*x + a[2]*x*x + a[1]*x + a[0];
     sigma	= 0.1*rand()/RAND_MAX;
     fprintf(fdata, "%3.3f %3.3lf %3.3lf\n",  x, y, y*sigma);
     x	+= dx;
 }
 
fclose(fdata);

}

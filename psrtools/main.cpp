#include <stdio.h>

#include "polifitgsl.h"

#define NP 11
double x[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
double y[] = {9.02, 1.01, 2.01, 10.02, 3.98, 5.97, 7.96, 7.94, 2.96, 11.96, 12.97, 0.97, 3.94, 6.94, 5.94, 4.92, 12.96, 13.9, 1.85, 5.9};

#define DEGREE 2
double coeff[DEGREE];

int main()
{
  int i;
  float* resid;
  resid = (float *)malloc(NP * sizeof(float));

  resid = polynomialfit(NP, DEGREE, x, y, coeff);
  for(i=0; i < NP; i++) {
	//printf("%lf\n", coeff[i]);
	printf("%lf %lf %f\n",x[i],y[i],resid[i]);
	
  }
  return 0;
}

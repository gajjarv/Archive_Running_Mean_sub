#ifndef _POLIFITGSL_H
#define _POLIFITGSL_H
#include <gsl/gsl_multifit.h>
#include <stdbool.h>
#include <math.h>
float* polynomialfit(int obs, int degree, 
           float *dx, float *dy, double *store); /* n, p */
#endif


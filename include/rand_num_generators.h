#ifndef _RAND_NUM_GENERATORS_H
#define _RAND_NUM_GENERATORS_H

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

#ifndef M_PI
#define M_PI (3.1415926535)
#endif

#include "ranlxd.h"
#include "common.h"

void       init_rand(int seed);
void       free_rand();

double     rand_double(double a, double b);                                                 //Random real number in the interval [a; b]
int        rand_int(int imax);                                                              //Random integer number in 0 .. imax - 1
double     rand_expbcos(double beta, double chi, int nbins, int *ntrials);                  //Random number distributed with the weight~ exp(B*cos(phi)), generated using the biased Metropolis algorithm with nbins bins
int        rand_poisson(double lambda);                                                     //Generate random integer according to Poisson distribution
double     rand_lorentz(double x0, double gamma);                                           //Generate random double according to Lorentz distribution
int        rand_choice(double *probs, int n);                                               //Chooses among n choices with probability probs[i]
double     rand_gauss_double(double offset, double sigma);                                  //Gaussian number with mean <offset> and dispersion sigma^2
t_complex  rand_gauss_complex(t_complex offset, double sigma);                              //Gaussian complex number with mean <offset> and dispersion sigma^2 (gaussian profile in the complex plane with dispersion sigma and offset <offset>)
void       rand_gauss_double_array(    double *A, int n, double offset, double sigma);      //Same as rand_gauss_double, but fills array of n elements with random numbers
void       rand_gauss_complex_array(t_complex *A, int n, t_complex offset, double sigma);   //Same as rand_gauss_complex, but fills array of n elements

#endif

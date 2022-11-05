#include "rand_num_generators.h"

void    init_rand(int seed)
{
 rlxd_init(2, seed);
}

void    free_rand()
{

}

double rand_double(double a, double b) /* Random real number in the interval [a; b] */
{
 double q;
 ranlxd(&q, 1);      
 return a + (b - a)*q;   
}

int rand_int(int imax) /* Random integer number in 0 .. imax - 1 */
{
 double q;
 ranlxd(&q, 1);
 return (int)(floor((double)imax*q));
}

double rand_expbcos(double beta, double chi, int nbins, int *ntrials)
{
 int i, itrials;
 double sprobs, s, alpha, phi_next;
 double q[3];    //Storage for random numbers
 //First generating the angle according to the "rough" distribution, approximating the real distribution by a histogram with nbins bins
 double *probs = (double *)malloc(nbins*sizeof(double));
 sprobs = 0.0;
 for(i=0; i<nbins; i++)
 {
  probs[i] = exp( beta*cos(M_PI*(double)i/(double)nbins) );
  sprobs += probs[i];
 };
 itrials = 0;
 do
 {
  ranlxd(q, 3);
  s = probs[0];
  i = 0;
  while((q[0]*sprobs>s)&&(i<nbins-1))
  {
   i++;
   s += probs[i];
  };
  phi_next  =  M_PI*((double)i + q[1])/(double)nbins;
  alpha     =  exp(beta*cos(phi_next))/probs[i];
  itrials ++;
  //Metropolis accept/reject step
 }while(q[2]>=alpha);
 free(probs);
 if(ntrials)
  (*ntrials) += itrials;
 //Now transforming phi_next so that it is a symmetric distribution on a unit circle centered around chi  
 ranlxd(q, 1);
 if(q[0] < 0.5)
  phi_next *= -1.0;
 phi_next += chi;
 if(phi_next>M_PI) 
  phi_next -= 2.0*M_PI;
 if(phi_next<-M_PI) 
  phi_next += 2.0*M_PI; 
 return phi_next;
}

int rand_choice(double *probs, int n)
{
 double q;
 ranlxd(&q, 1);
 double s  = probs[0];
 int i = 0;
 while((q>s)&&(i<n-1))
 {
  i++;         
  s += probs[i];          
 };                
 return i;
}

int rand_poisson(double lambda)
{
 double l = exp(-lambda);
 int    k = 0;
 double p = 1.0;
 double q;
 do
 {
  ranlxd(&q, 1);
  k++;
  p *= q;
 }while(p > l);
 return (k-1);   
}

double rand_lorentz(double x0, double gamma)
{
  double q = rand_double(0, 1);
  double x = x0 + gamma * tan(M_PI * (q - 0.5));
  return x;
}

double rand_gauss_double(double offset, double sigma)
{
 double r[2];
 ranlxd(r, 2);
 return (offset + sigma*sqrt(-2.0*log(r[0]))*cos(2.0*M_PI*r[1]));
}

t_complex rand_gauss_complex(t_complex offset, double sigma)
{
 double r[2];
 ranlxd(r, 2);
 return (offset + sigma*sqrt(-log(r[0]))*exp(2.0*M_PI*t_complex(0, 1)*r[1]));
}

void       rand_gauss_double_array(    double *A, int n, double offset, double sigma)   //Same as rand_gauss_double, but fills array of n elements with random numbers
{
 int i;
 for(i=0; i<n; i++)
  A[i] = rand_gauss_double(offset, sigma);
}

void       rand_gauss_complex_array(t_complex *A, int n, t_complex offset, double sigma)   //Same as rand_gauss_complex, but fills array of n elements
{
 int i;
 for(i=0; i<n; i++)
  A[i] = rand_gauss_complex(offset, sigma);
}

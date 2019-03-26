#ifndef CBLAS_TEMPLATES_H_
#define CBLAS_TEMPLATES_H_

#include "common.h"
#include <lapacke.h>
#include <cblas.h>

// TODO: Implement for types other than real and complex numbers
template<typename T> struct CblasHelper
{
  static void cblas_copy(const int N, const void *X, const int incX, void *Y, const int incY);
  static void cblas_swap(const int N, void *X, const int incX, void *Y, const int incY);
  static void cblas_scal(const int N, const void *alpha, void *X, const int incX);
  static void cblas_axpy(const int N, const void *alpha, const void *X, const int incX, void *Y, const int incY);
  static double cblas_nrm2(const int N, const void *X, const int incX);
  static void cblas_dot_sub(const int N, const void *X, const int incX, const void *Y, const int incY, void *dotc);
};

template<>
struct CblasHelper<float>
{
  static void cblas_copy(const int N, const void *X, const int incX, void *Y, const int incY)
  {
    cblas_scopy(N, (const float *)X, incX, (float *)Y, incY);
  }
  static void cblas_swap(const int N, void *X, const int incX, void *Y, const int incY)
  {
    cblas_sswap(N, (float *)X, incX, (float *)Y, incY);
  }
  static void cblas_scal(const int N, const float *alpha, void *X, const int incX)
  {
    cblas_sscal(N, *alpha, (float *)X, incX);
  }
  static void cblas_axpy(const int N, const float *alpha, const void *X, const int incX, void *Y, const int incY)
  {
    cblas_saxpy(N, *alpha, (const float *)X, incX, (float *)Y, incY);
  }
  static float cblas_nrm2(const int N, const void *X, const int incX)
  {
    return cblas_snrm2(N, (const float *)X, incX);
  }
  static void cblas_dot_sub(const int N, const void *X, const int incX, const void *Y, const int incY, void *dot)
  {
    float *val = (float *) dot;
    (*val) = cblas_sdot(N, (const float *)X, incX, (const float *)Y, incY);
  }
};

template<>
struct CblasHelper<double>
{
  static void cblas_copy(const int N, const void *X, const int incX, void *Y, const int incY)
  {
    cblas_dcopy(N, (const double *)X, incX, (double *)Y, incY);
  }
  static void cblas_swap(const int N, void *X, const int incX, void *Y, const int incY)
  {
    cblas_dswap(N, (double *)X, incX, (double *)Y, incY);
  }
  static void cblas_scal(const int N, const double *alpha, void *X, const int incX)
  {
    cblas_dscal(N, *alpha, (double *)X, incX);
  }
  static void cblas_axpy(const int N, const double *alpha, const void *X, const int incX, void *Y, const int incY)
  {
    cblas_daxpy(N, *alpha, (const double *)X, incX, (double *)Y, incY);
  }
  static double cblas_nrm2(const int N, const void *X, const int incX)
  {
    return cblas_dnrm2(N, (const double *)X, incX);
  }
  static void cblas_dot_sub(const int N, const void *X, const int incX, const void *Y, const int incY, void *dot)
  {
    double *val = (double *) dot;
    (*val) = cblas_ddot(N, (const double *)X, incX, (const double *)Y, incY);
  }
};

template<>
struct CblasHelper<t_complex_float>
{
  static void cblas_copy(const int N, const void *X, const int incX, void *Y, const int incY)
  {
    cblas_ccopy(N, (const float *)X, incX, (float *)Y, incY);
  }
  static void cblas_swap(const int N, void *X, const int incX, void *Y, const int incY)
  {
    cblas_cswap(N, (float *)X, incX, (float *)Y, incY);
  }
  static void cblas_scal(const int N, const void *alpha, void *X, const int incX)
  {
    cblas_cscal(N, (const float *)alpha, (float *)X, incX);
  }
  static void cblas_axpy(const int N, const void *alpha, const void *X, const int incX, void *Y, const int incY)
  {
    cblas_caxpy(N, (const float *)alpha, (const float *)X, incX, (float *)Y, incY);
  }
  static float cblas_nrm2(const int N, const void *X, const int incX)
  {
    return cblas_scnrm2(N, (const float *)X, incX);
  }
  static void cblas_dot_sub(const int N, const void *X, const int incX, const void *Y, const int incY, void *dotc)
  {
    cblas_cdotc_sub(N, (const float *)X, incX, (const float *)Y, incY, (openblas_complex_float *)dotc);
  }
};

template<>
struct CblasHelper<t_complex>
{
  static void cblas_copy(const int N, const void *X, const int incX, void *Y, const int incY)
  {
    cblas_zcopy(N, (const double *)X, incX, (double *)Y, incY);
  }
  static void cblas_swap(const int N, void *X, const int incX, void *Y, const int incY)
  {
    cblas_zswap(N, (double *)X, incX, (double *)Y, incY);
  }
  static void cblas_scal(const int N, const void *alpha, void *X, const int incX)
  {
    cblas_zscal(N, (const double *)alpha, (double *)X, incX);
  }
  static void cblas_axpy(const int N, const void *alpha, const void *X, const int incX, void *Y, const int incY)
  {
    cblas_zaxpy(N, (const double *)alpha, (const double *)X, incX, (double *)Y, incY);
  }
  static double cblas_nrm2(const int N, const void *X, const int incX)
  {
    return cblas_dznrm2(N, (const double *)X, incX);
  }
  static void cblas_dot_sub(const int N, const void *X, const int incX, const void *Y, const int incY, void *dotc)
  {
    cblas_zdotc_sub(N, (const double *)X, incX, (const double *)Y, incY, (openblas_complex_double *)dotc);
  }
};

#endif

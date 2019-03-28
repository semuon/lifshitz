#include "linalg.h"
#include <lapacke.h>
#include "profiler.h"

#define CHECK_LAPACK_INFO(info)  \
  if ((info) != 0)         \
  {    \
    THROW_EXCEPTION(LapackFailure, "LAPACK routine returned with a error code " + TO_STRING(info)); \
  }

void Linalg::LapackInvert(t_complex *A, int n)
{
  int info;
  VECTOR<int> piv(n);

  info = LAPACKE_zgetrf(LAPACK_ROW_MAJOR, n, n, A, n, piv.data());
  CHECK_LAPACK_INFO(info);

  info = LAPACKE_zgetri(LAPACK_ROW_MAJOR, n, A, n, piv.data());
  CHECK_LAPACK_INFO(info);
}

void Linalg::LapackLinearSolve(t_complex *x, t_complex *A, int n)
{
  int info;
  VECTOR<int> piv(n);

  info = LAPACKE_zgetrf(LAPACK_ROW_MAJOR, n, n, A, n, piv.data());
  CHECK_LAPACK_INFO(info);

  info = LAPACKE_zgetrs(LAPACK_ROW_MAJOR, 'N', n, 1, A, n, piv.data(), x, 1);
  CHECK_LAPACK_INFO(info);
}

void Linalg::LapackLinearSolve(double *x, double *A, int n)
{
  int info;
  VECTOR<int> piv(n);

  info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, A, n, piv.data());
  CHECK_LAPACK_INFO(info);

  info = LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N', n, 1, A, n, piv.data(), x, 1);
  CHECK_LAPACK_INFO(info);
}

void Linalg::LapackLinearSolveSymmetric(t_complex *x, t_complex *A, int n)
{
  int info;
  int lwork;
  t_complex work_query;
  VECTOR<int> piv(n);
  VECTOR<t_complex> work;

  info = LAPACKE_zsytrf_work(LAPACK_ROW_MAJOR, 'U', n, NULL, n, piv.data(), &work_query, -1);
  CHECK_LAPACK_INFO(info);

  lwork = LAPACK_Z2INT(work_query);
  work.resize(lwork);

  info = LAPACKE_zsytrf_work(LAPACK_ROW_MAJOR, 'U', n, A, n, piv.data(), work.data(), lwork);
  CHECK_LAPACK_INFO(info);

  info = LAPACKE_zsytrs(LAPACK_ROW_MAJOR, 'U', n, 1, A, n, piv.data(), x, 1);
  CHECK_LAPACK_INFO(info);
}

void Linalg::LapackEigensystem(t_complex *A, t_complex *evals, t_complex *evecs, int n)
{
  int info;

  info = LAPACKE_zgeev(LAPACK_ROW_MAJOR, 'N', 'V', n, A, n, evals, NULL, n, evecs, n);
  CHECK_LAPACK_INFO(info);
}

void Linalg::LapackHermitianEigensystem(t_complex *A, double *evals, int n)
{
  int info;
  int lwork;
  t_complex work_query;
  VECTOR<double> rwork;
  VECTOR<t_complex> work;

  rwork.resize(MAX(1 , 3 * n - 2));
  info = LAPACKE_zheev_work( LAPACK_ROW_MAJOR, 'V', 'U', n, NULL, n, NULL, &work_query, -1, rwork.data());
  CHECK_LAPACK_INFO(info);

  lwork = LAPACK_Z2INT(work_query);
  work.resize(lwork);

  info = LAPACKE_zheev_work(LAPACK_ROW_MAJOR, 'V', 'U', n, A, n, evals, work.data(), lwork, rwork.data());
  CHECK_LAPACK_INFO(info);
}

void Linalg::LapackCholeskyDecomposition(t_complex *A, int n)
{
  int info;

  info = LAPACKE_zpotrf(LAPACK_ROW_MAJOR, 'L', n, A, n);
  CHECK_LAPACK_INFO(info);
}

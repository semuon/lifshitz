#include "linalg.h"
#include <lapacke.h>
#include "profiler.h"

#define CHECK_ARPACK_INFO(info)  \
  if ((info) != 0)         \
  {    \
    THROW_EXCEPTION(ArnoldiFailure, "Arnoldi returned with a error code " + TO_STRING(info)); \
  }

//Interface to Fortran code of ARPACK
extern "C"
{
  extern void znaupd_(int*, char*, int*, char*, int*, double*, t_complex*, 
                      int*, t_complex*, int*, int*, int*, t_complex*, t_complex*,
                      int*, double* ,int*, int, int);
  extern void zneupd_(int*, char*, int*, t_complex*, t_complex*, int*, t_complex*,
                      t_complex*, 
                      char*, int*, char*, int*, double*, t_complex*, 
                      int*, t_complex*, int*, int*, int*, t_complex*, t_complex*,
                      int*, double* ,int*, int, int, int);
}

void Linalg::Arnoldi(MatrixByVector F, void *args, int n, VECTOR<t_complex> &evals,
                     int nev, double tol, int max_iter, ArnoldiMode mode)
{
  VECTOR<BaseLinearVector> evecs;
  bool compute_evecs = false;

  return Arnoldi(F, args, n, evals, compute_evecs, evecs, nev, tol, max_iter, mode,
                 c_default_arnoldi_ncv_factor_nom, c_default_arnoldi_ncv_factor_denom);
}

void Linalg::Arnoldi(MatrixByVector F, void *args, int n, VECTOR<t_complex> &evals,
                     bool compute_evecs, VECTOR<BaseLinearVector> &evecs, int nev, double tol, int max_iter,
                     ArnoldiMode mode, int arnoldi_ncv_factor_nom, int arnoldi_ncv_factor_denom)
{
  ASSERT(evals.size() == (uint)nev);

  if (compute_evecs)
    ASSERT(evecs.size() == (uint)nev);

  pGlobalProfiler.StartTimer("ARNOLDI");

  int i, ido, ncv;                  
  char bmat = 'I';            // standart eigenvalue problem 
  char which[3];
  VECTOR<t_complex> resid;
  VECTOR<t_complex> v;
  VECTOR<t_complex> workd;
  VECTOR<t_complex> workl;
  VECTOR<double>    rwork;
  VECTOR<BaseLinearVector> sf_workd;
  int          ldv   = n;
  int          ldz   = n;
  int   iparam[11];
  int    ipntr[14];
  int       lworkl;
  int         info   = 0;
  int  count;
  bool cont;

  // neupd variables 
  int         rvec; // 0 - eigenvalues, 1 - eigenvectors + eigenvalues
  char howmny = 'A';
  t_complex sigma;
  VECTOR<t_complex> d;
  VECTOR<t_complex> workev;
  VECTOR<int> select;

  switch(mode)
  {
    case SmallestAbs:
      sprintf(which, "SM");
      break;
    case LargestAbs:
      sprintf(which, "LM");
      break;
    case SmallestReal:
      sprintf(which, "SR");
      break;
    case LargestReal:
      sprintf(which, "LR");
      break;
    case SmallestImag:
      sprintf(which, "SI");
      break;
    case LargestImag:
      sprintf(which, "LI");
      break;     
    default:
      throw AssertFailure("Uknown Arnoldi mode: " + TO_STRING(mode));
    break; 
  }; 

  ncv = (arnoldi_ncv_factor_nom * nev) / arnoldi_ncv_factor_denom;
  if(ncv - nev < 2)
    ncv = 2 + nev;
  if(ncv > n)
    ncv = n;
  lworkl = 3 * ncv * ncv + 5 * ncv; 

  pGlobalProfiler.StartTimer("ARNOLDI: memory allocations");

  // naupd arguments
  resid.resize(n);
  v.resize(n * ncv);
  workd.resize(3 * n);
  workl.resize(lworkl);
  rwork.resize(ncv);
  sf_workd.reserve(3);
  for(i = 0; i < 3; i++)
  {
    // Construct 3 vectors upon workd memory
    // We shall pass this to F
    sf_workd.emplace_back(n, &workd[i * n]);
  }

  pGlobalProfiler.StopTimer("ARNOLDI: memory allocations");

  for(i = 0; i < 11; i++)
  {
    iparam[i] = 0;
  }
  iparam[0] = 1;
  iparam[2] = max_iter;
  iparam[3] = 1;
  iparam[6] = 1;

  // ARPACK arnoldi routine
  cont  = true;
  count = 0;
  ido   = 0; // first call to reverse communication interface

  while(cont)
  {
    znaupd_(&ido, &bmat, &n, which, &nev, &tol, resid.data(), &ncv, v.data(),
           &ldv, iparam, ipntr, workd.data(), workl.data(), &lworkl, rwork.data(), &info, 1, 2);

    if(abs(ido) == 1)
    {
      int idx_in = (ipntr[0] - 1) / n;
      int idx_out = (ipntr[1] - 1) / n;

      F(sf_workd[idx_out], sf_workd[idx_in], args);
    }
    else
    {
      cont = false;
    }

    count ++;
  };

  CHECK_ARPACK_INFO(info);

  // get eigenvalues
  rvec = (compute_evecs ? 1 : 0);
  howmny = 'A';

  select.resize(ncv);
  d.resize(nev + 1);
  workev.resize(2 * ncv);

  // we can use already allocated array v to store eigenvectors if
  // we are not interested in Schur basis. Read neupd.f for details 
  sigma = 0.0 + I * 0.0;

  // Calling neupd to finally calculate the eigenspectrum
  zneupd_(&rvec, &howmny, select.data(), d.data(), v.data(), &ldz, &sigma, workev.data(), &bmat, &n, which,
          &nev, &tol, resid.data(), &ncv, v.data(), &ldv, iparam, ipntr, workd.data(), workl.data(),
          &lworkl, rwork.data(), &info, 1, 2, 1);

  CHECK_ARPACK_INFO(info);

  // Save eigenvalues
  for(i = 0; i < nev; i++)
  {
    evals[i] = d[i];
  }

  // Save eigenvectors if requested
  if(compute_evecs)
  {
    for(i = 0; i < nev; i++)
    {
      t_complex *evec_ptr = &v[i * n];
      memcpy(evecs[i].DataPtr(), evec_ptr, sizeof(t_complex) * n);
    }
  }

  pGlobalProfiler.StopTimer("ARNOLDI");
}
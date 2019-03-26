#include "linalg.h"
#include <lapacke.h>
#include "profiler.h"

// ---------------------------------------------------------------- //
// ---------------------------------------------------------------- //
// ---------------------------------------------------------------- //
//                             LAPACK                               //
// ---------------------------------------------------------------- //
// ---------------------------------------------------------------- //
// ---------------------------------------------------------------- //

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

#undef CHECK_LAPACK_INFO

// ---------------------------------------------------------------- //
// ---------------------------------------------------------------- //
// ---------------------------------------------------------------- //
//                               CG                                 //
// ---------------------------------------------------------------- //
// ---------------------------------------------------------------- //
// ---------------------------------------------------------------- //

// This is algorithm from Jegerlehner, B. , http://arxiv.org/abs/hep-lat/9612014, indexes in last line changed i->i+1
// F|solution > = |source > , |solution > = F^-1 |source>
void Linalg::MultishiftCG(Linalg::MatrixByVectorFunc F, void *args, const Lattice &lat, const SpinorField &source,
                          int nsys, VECTOR<SpinorField> &solution, VECTOR<double> sigma,
                          double tol, int imax, double &max_err, int &num_iter)
{

// This macro prints out all variables to stdout
#define PRINT_CG_STATE \
      {\
        std::cout << "==================================================" << std::endl; \
        std::cout << "=================== MCG STATE ====================" << std::endl; \
        std::cout << "==================================================" << std::endl; \
        std::cout << "    requested tolerance: " << tol << std::endl; \
        std::cout << "    step: " << counter << std::endl; \
        std::cout << "    c: " << c << std::endl; \
        std::cout << "    d: " << d << std::endl << std::endl; \
        std::cout << "    a: " << a << std::endl; \
        std::cout << "    a_next: " << a_next << std::endl << std::endl; \
        std::cout << "    b_prev: " << b_prev << std::endl; \
        std::cout << "    b: " << b << std::endl; \
        std::cout << "    b_next: " << b_next << std::endl << std::endl; \
        for (int j = 0; j < nsys; j++) \
          std::cout << "    a_s_prev[" << j << "]: " << a_s_prev[j] << std::endl; \
        std::cout << std::endl; \
        for (int j = 0; j < nsys; j++) \
          std::cout << "    a_s[" << j << "]: " << a_s[j] << std::endl; \
        std::cout << std::endl; \
        for (int j = 0; j < nsys; j++) \
          std::cout << "    a_s_next[" << j << "]: " << a_s_next[j] << std::endl; \
        std::cout << std::endl; \
        for (int j = 0; j < nsys; j++) \
          std::cout << "    b_s_prev[" << j << "]: " << b_s_prev[j] << std::endl; \
        std::cout << std::endl; \
        for (int j = 0; j < nsys; j++) \
          std::cout << "    b_s[" << j << "]: " << b_s[j] << std::endl; \
        std::cout << std::endl; \
        for (int j = 0; j < nsys; j++) \
          std::cout << "    b_s_next[" << j << "]: " << b_s_next[j] << std::endl; \
        std::cout << std::endl; \
        for (int j = 0; j < nsys; j++) \
          std::cout << "    z_s_prev[" << j << "]: " << z_s_prev[j] << std::endl; \
        std::cout << std::endl; \
        for (int j = 0; j < nsys; j++) \
          std::cout << "    z_s[" << j << "]: " << z_s[j] << std::endl; \
        std::cout << std::endl; \
        for (int j = 0; j < nsys; j++) \
          std::cout << "    z_s_next[" << j << "]: " << z_s_next[j] << std::endl; \
        std::cout << std::endl; \
      }

// This macro tests condition x and breaks algorithm when test fails
#define CG_ASSERT(x) \
{ \
  if (!(x)) \
  { \
    WRITE_ERROR("MCG ASSERT FAILED"); \
    PRINT_CG_STATE; \
    THROW_EXCEPTION(CGFailure, "MCG assertion failed"); \
  } \
}

  ASSERT(nsys > 0);
  ASSERT(solution.size() == (uint)nsys);
  ASSERT(sigma.size() == (uint)nsys);

  pGlobalProfiler.StartTimer("MCG");

  int counter = 0, i, minnumber;
  double err = 0.0, min;
  double b, b_prev, b_next = 0.0, c, d;
  double err_sys;
  double a, a_next;

  VECTOR<bool> sys_conv(nsys, false);

  VECTOR<double> z_s_prev(nsys);
  VECTOR<double> z_s     (nsys);
  VECTOR<double> z_s_next(nsys);

  VECTOR<double> a_s_prev(nsys);
  VECTOR<double> a_s     (nsys);
  VECTOR<double> a_s_next(nsys);

  VECTOR<double> b_s_prev(nsys);
  VECTOR<double> b_s     (nsys);
  VECTOR<double> b_s_next(nsys);
  VECTOR<double> err_s(nsys);

  pGlobalProfiler.StartTimer("MCG: memory allocations");

  SpinorField chck(lat);
  SpinorField r   (lat);
  SpinorField Ar  (lat);
  VECTOR<SpinorField> x_s;
  VECTOR<SpinorField> p_s;
  VECTOR<SpinorField> r_s;
  VECTOR<SpinorField> q_s;

  x_s.reserve(nsys);
  p_s.reserve(nsys);
  r_s.reserve(nsys);
  q_s.reserve(nsys);
  for(i = 0; i < nsys; i++)
  {
    x_s.emplace_back(lat);
    p_s.emplace_back(lat);
    r_s.emplace_back(lat);
    q_s.emplace_back(lat);
  }

  pGlobalProfiler.StopTimer("MCG: memory allocations");

  max_err = 0;

  min = fabs(sigma[0]);
  minnumber = 0;
  for (i = 0; i < nsys; i++)
  {
    if (min > fabs(sigma[i]))
    {
      min = fabs(sigma[i]);
      minnumber = i;
    }
  }

  // Initial step (I)
  for (i = 0; i < nsys; i++)
    x_s[i].Clear();

  // r_0 = source
  r = source;
  // p_s = source , p = source
  // q_s = A * source
  F(Ar, source, args);
  Ar.Add_bB(sigma[minnumber], source);

  for (i = 0; i < nsys; i++)
  {
    p_s[i] = source;
    r_s[i] = source;
    q_s[i] = Ar;
  }

  // b_prev = z_s_prev = z_s = 1
  // a_s = 0
  b_prev = 1.0;
  a = 0.0;
  a_s.assign(nsys, 0);
  z_s_prev.assign(nsys, 1.0);
  z_s.assign(nsys, 1.0);

  // Iteration step (II)
  for(counter = 0; counter < imax; counter++)
  {
    // b = -(r  , r ) / (p_s  , F(p_s) )
    c = r.NormSqr();
    d = real(p_s[minnumber] * q_s[minnumber]);
    b = -c / d;
    CG_ASSERT(std::isfinite(b));

    for (i = 0; i < nsys; i++)
    {
      if(sys_conv[i])
        continue;

      // z_s_next = z_s * z_s_prev * b_prev/ (b*a (z_s_prev - z_s) + (z_s_prev*b_prev)*(1-sigma*b))
      z_s_next[i] = z_s[i] * z_s_prev[i] * b_prev / (b * a * (z_s_prev[i] - z_s[i]) + z_s_prev[i] * b_prev * (1.0 - (sigma[i] - sigma[minnumber]) * b));
      CG_ASSERT(std::isfinite(z_s_next[i]));
    }

    for (i = 0; i < nsys; i++)
    {
      if(sys_conv[i])
        continue;

      // b_s = b * z_s_next / z_s
      b_s[i] = b * z_s_next[i] / z_s[i];
      CG_ASSERT(std::isfinite(b_s[i]));

      // x_s = x_s - b_s * p_s
      x_s[i].Add_bB(-1.0 * b_s[i], p_s[i]);
    }

    // r =  r + b * F(p_s)
    r.Add_bB(b, q_s[minnumber]);
    err = r.NormSqr();

    for (i = 0; i < nsys; i++)
    {
      if(sys_conv[i])
        continue;

      r_s[i].Add_bB(b_s[i], q_s[i]);
      r_s[i].Add_bB(b_s[i] * (sigma[i] - sigma[minnumber]), p_s[i]);
      err_s[i] = r_s[i].Norm();

      // Test system convergence
      if (err_s[i] < tol)
      {
        F(chck, x_s[i], args);
        chck.Add_bB(sigma[i], x_s[i]);
        chck.Add_bB(-1.0, source);
        err_sys = chck.Norm();

        if(err_sys < tol)
        {
          sys_conv[i] = true;

          if (max_err < err_sys)
            max_err = err_sys;
        }
      }
    }

    // a_next = ( r_next , r_next ) /( r_old, r_old )
    a_next = err / c;
    CG_ASSERT(std::isfinite(a_next));

    F(Ar, r, args);
    Ar.Add_bB(sigma[minnumber], r);

    for (i = 0; i < nsys; i++)
    {
      if(sys_conv[i])
        continue;

      // a_s_next = a_next * z_s_next * b_s / (z_s * b)
      a_s_next[i] = a_next * z_s_next[i] * b_s[i] / (z_s[i] * b);
      CG_ASSERT(std::isfinite(a_s_next[i]));

      // p_s = a_s * p_s + z_s * r
      p_s[i].Assign_aA_plus_bB(a_s_next[i], z_s_next[i], r);
      // q_s = a_s * q_s + z_s * A * r
      q_s[i].Assign_aA_plus_bB(a_s_next[i], z_s_next[i], Ar);
    }

    a = a_next;
    b_prev = b;
    b = b_next;

    for (i = 0; i < nsys; i++)
    {
      if(sys_conv[i])
        continue;

      a_s_prev[i] = a_s[i];
      a_s[i] = a_s_next[i];

      b_s_prev[i] = b_s[i];
      b_s[i] = b_s_next[i];

      z_s_prev[i] = z_s[i];
      z_s[i] = z_s_next[i];
    }

    // Check that all systems converged
    bool all_converged = true;
    for(i = 0; i < nsys; i++)
    {
      all_converged = all_converged && sys_conv[i];
    }

    if(all_converged)
      break;
  }

  num_iter = counter;

  bool success = true;
  for (i = 0; i < nsys; i++)
  {
    solution[i] = x_s[i];

    if(!sys_conv[i])
      success = false;
  }

  pGlobalProfiler.StopTimer("MCG");

  if (!success)
  {
    THROW_EXCEPTION(CGFailure, "Multishift Conjugate Gradient was not able to find solution within given number of iterations");
  }

#undef CG_ASSERT
#undef PRINT_CG_STATE
}

// ---------------------------------------------------------------- //
// ---------------------------------------------------------------- //
// ---------------------------------------------------------------- //
//                             ARNOLDI                              //
// ---------------------------------------------------------------- //
// ---------------------------------------------------------------- //
// ---------------------------------------------------------------- //

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

void Linalg::Arnoldi(MatrixByVectorFunc F, void *args, const Lattice &lat, VECTOR<t_complex> &evals,
                     int nev, double tol, int max_iter, ArnoldiMode mode)
{
  VECTOR<SpinorField> evecs;
  bool compute_evecs = false;

  return Arnoldi(F, args, lat, evals, compute_evecs, evecs, nev, tol, max_iter, mode,
                 c_default_arnoldi_ncv_factor_nom, c_default_arnoldi_ncv_factor_denom);
}

void Linalg::Arnoldi(MatrixByVectorFunc F, void *args, const Lattice &lat, VECTOR<t_complex> &evals,
                     bool compute_evecs, VECTOR<SpinorField> &evecs, int nev, double tol, int max_iter,
                     ArnoldiMode mode, int arnoldi_ncv_factor_nom, int arnoldi_ncv_factor_denom)
{

#define CHECK_ARPACK_INFO(info)  \
  if ((info) != 0)         \
  {    \
    THROW_EXCEPTION(ArnoldiFailure, "Arnoldi returned with a error code " + TO_STRING(info)); \
  }

  ASSERT(evals.size() == (uint)nev);

  if (compute_evecs)
    ASSERT(evecs.size() == (uint)nev);

  pGlobalProfiler.StartTimer("ARNOLDI");

  int n = lat.Volume() * lat.Ndirac() * lat.SUNrank();

  int i, ido, ncv;                  
  char bmat = 'I';            // standart eigenvalue problem 
  char which[3];
  VECTOR<t_complex> resid;
  VECTOR<t_complex> v;
  VECTOR<t_complex> workd;
  VECTOR<t_complex> workl;
  VECTOR<double>    rwork;
  VECTOR<SpinorField> sf_workd;
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
    // Construct 3 SpinorFields upon workd memory
    // We shall pass this to F
    sf_workd.emplace_back(lat, &workd[i * n]);
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

#undef CHECK_ARPACK_INFO
}

// Auxiliary struct to keep eigenvalues and eigenvectors in pairs
struct eigen { 
  int i;
  t_complex val;
  Linalg::SortingMode mode;

  eigen(int i, t_complex val, Linalg::SortingMode mode) :
    i(i), val(val), mode(mode) {}

  bool operator<(eigen const &other) const
  {
    switch(mode)
    {
      case Linalg::AscendingAbs:
        return abs(val) < abs(other.val);
      case Linalg::DescendingAbs:
        return abs(val) > abs(other.val);
      case Linalg::AscendingReal:
        return real(val) < real(other.val);
      case Linalg::DescendingReal:
        return real(val) > real(other.val);
      case Linalg::AscendingImag:
        return imag(val) < imag(other.val);
      case Linalg::DescendingImag:
        return imag(val) > imag(other.val);
      default:
        THROW_EXCEPTION(AssertFailure, "Unknown sorting mode");
    }
  }
};
 
void Linalg::SortEigensystem(VECTOR<int> &sorted_idx, const VECTOR<t_complex> &evals, const SortingMode mode)
{
  int n = evals.size();
  VECTOR<eigen> pairs;

  pairs.reserve(n);
  for(int i = 0; i < n; i++)
    pairs.emplace_back(i, evals[i], mode);

  std::sort(pairs.begin(), pairs.end());

  sorted_idx.resize(n);
  for(int i = 0; i < n; i++)
    sorted_idx[i] = pairs[i].i;
}

void Linalg::TestOrthogonality(VECTOR<SpinorField> &evecs, const VECTOR<int> idx_to_test, const double tol)
{
 ASSERT(idx_to_test.size() <= evecs.size());
  
 bool should_check = true;
 bool first_check = true; 
 
 while(should_check)
  {
    bool success = true;
  
   for(uint iv = 0; iv < idx_to_test.size(); iv++)
   {
     int i = idx_to_test.at(iv);
      VECTOR<int> vec_for_orth;

      for(uint jv = iv + 1; jv < idx_to_test.size(); jv++)
      {
        int j = idx_to_test.at(jv);

        t_complex sp = evecs.at(i) * evecs.at(j);
        double spnorm = abs(sp);
        if (spnorm > tol)
        {
          success = false;

          vec_for_orth.push_back(j);
        }
      }

      // If there are some non-orthogonal vectors, do Gramm-Schmidt
      if (vec_for_orth.size() > 0)
      {
        vec_for_orth.push_back(i);

        Linalg::GrammSchmidt(vec_for_orth, evecs);
      }
    }

    if (!first_check && !success)
    {
      THROW_EXCEPTION(EigenTestFailure,
            "Eigenvectors are not orthogonal (FAILED TO FIX)");
    }

    should_check = first_check && !success;
    first_check = false;
  };
}

void Linalg::TestEigensystem(MatrixByVectorFunc F, void *args,
                             const VECTOR<t_complex> &evals, VECTOR<SpinorField> &evecs,
                             const VECTOR<int> idx_to_test, const double tol)
{
  ASSERT(evals.size() == evecs.size());
  ASSERT(idx_to_test.size() <= evecs.size());

  if (evecs.size() == 0)
    return;

  SpinorField tmp(evecs.at(0));

  for(uint iv = 0; iv < idx_to_test.size(); iv++)
  {
    int i = idx_to_test.at(iv);

    F(tmp, evecs.at(i), args);
    tmp.Add_bB(-1.0 * evals.at(i), evecs.at(i));

    double resid = tmp.Norm();
    if (resid > tol)
    {
       THROW_EXCEPTION(EigenTestFailure,
         "Eigensystem " + TO_STRING(i) + " is not accurate: |H * v - lambda * v| = " + TO_STRING(resid) + " > tol = " + TO_STRING(tol) + "");
    }

    resid = fabs(evecs.at(i).Norm() - 1);
    if (resid > tol)
    {
      evecs.at(i).Normalize();
      // THROW_EXCEPTION(EigenTestFailure,
      //     "Eigenvector norm is not 1: |v_" + TO_STRING(i) + "| - 1 = " + TO_STRING(resid) + " > tol = " + TO_STRING(tol) + "");
    }
  }

}

void Linalg::GrammSchmidt(const VECTOR<int> &vec_idx, VECTOR<SpinorField> &vecs)
{
  ASSERT(vec_idx.size() <= vecs.size());

  if (vec_idx.size() == 0)
    return;

  const int n = vec_idx.size();

  for(int i = 0; i < n; i++)
  {
    int idx_i = vec_idx.at(i);

    SpinorField &v_i = vecs.at(idx_i);

    for(int j = 0; j < i; j++)
    {
      int idx_j = vec_idx.at(j);

      SpinorField &u_j = vecs.at(idx_j);

      u_j.Add_bB(-1.0 * (v_i * u_j) / v_i.NormSqr(), v_i);
    }
  }

  for(int i = 0; i < n; i++)
  {
    int idx_i = vec_idx.at(i);
    vecs.at(idx_i).Normalize();
  }
}

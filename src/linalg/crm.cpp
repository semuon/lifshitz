#include "linalg.h"
#include <lapacke.h>
#include "profiler.h"

// This macro prints out all variables to stdout
#define PRINT_CR_STATE \
      {\
        std::cout << "==================================================" << std::endl; \
        std::cout << "=================== CRM STATE ====================" << std::endl; \
        std::cout << "==================================================" << std::endl; \
        std::cout << "    requested tolerance: " << tol << std::endl; \
        std::cout << "    step: " << counter << std::endl; \
        std::cout << "    err: " << err << std::endl; \
        std::cout << "    a: " << a << std::endl; \
        std::cout << "    b: " << b << std::endl; \
        std::cout << "    c: " << c << std::endl; \
        std::cout << "    d: " << d << std::endl; \
      }

// This macro tests condition x and breaks algorithm when test fails
#define CR_ASSERT(x) \
{ \
  if (!(x)) \
  { \
    WRITE_ERROR("CRM ASSERT FAILED"); \
    PRINT_CR_STATE; \
    THROW_EXCEPTION(CGFailure, "CRM assertion failed"); \
  } \
}

// Standard Conjugate Residual Method
void Linalg::CRM(Linalg::MatrixByVector F, void *args,
                 const LinearVector &source, LinearVector &solution,
                 double tol, int imax, double &max_err, int &num_iter)
{
  ASSERT(source.Count() == solution.Count());

  uint vol = source.Count();
  const Lattice &lat = solution.GetLattice();

  int counter;
  double a, b, c, d;

  double err;
  bool is_converged = false;

  pGlobalProfiler.StartTimer("CRM");

  pGlobalProfiler.StartTimer("CRM: memory allocations");

  LinearVector x_k(lat, vol);
  LinearVector p_k(lat, vol);
  LinearVector r_k(lat, vol);
  LinearVector Ar_k(lat, vol);
  LinearVector Ap_k(lat, vol);
  LinearVector chck(lat, vol);

  pGlobalProfiler.StopTimer("CRM: memory allocations");

  // Initial step (k = 0)
  x_k.Clear();

  // r_0 = source - A x_0
  r_k = source;
  F(Ar_k, x_k, args);
  r_k.Add_bB(-1.0, Ar_k);
  // p_0 = r_0
  p_k = r_k;

  // Iteration steps (k = 0, 1, ...)
  for(counter = 0; counter < imax; counter++)
  {
    F(Ar_k, r_k, args);
    F(Ap_k, p_k, args);

    // a = (r_k  , F(r_k) ) / (F(p_k)  , F(p_k) )
    c = real(r_k * Ar_k);
    d = Ap_k.NormSqr();
    a = c / d;
    CR_ASSERT(std::isfinite(a));

    x_k.Assign_aA_plus_bB(1.0, a, p_k);
    r_k.Assign_aA_plus_bB(1.0, -a, Ap_k);

    err = r_k.Norm();

    if (err < tol)
    {
      F(chck, x_k, args);
      chck.Add_bB(-1.0, source);
      err = chck.Norm();

      if(err < tol)
      {
        is_converged = true;

        solution = x_k;

        if (max_err < err)
          max_err = err;

        break;
      }
    }

    // b = (r_{k+1} , F(r_{k+1}) ) / (r_k  , F(r_k) )
    F(Ar_k, r_k, args);
    d = real(r_k * Ar_k);
    b = d / c;
    CR_ASSERT(std::isfinite(b));

    p_k.Assign_aA_plus_bB(b, 1.0, r_k);
  }

  num_iter = counter;

  pGlobalProfiler.StopTimer("CRM");

  if (!is_converged)
  {
    PRINT_CR_STATE;
    THROW_EXCEPTION(CGFailure, "Conjugate Residual Method was not able to find solution within given number of iterations");
  }
}
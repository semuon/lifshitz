#include "linalg.h"
#include <lapacke.h>
#include "profiler.h"

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

// This is algorithm from Jegerlehner, B. , http://arxiv.org/abs/hep-lat/9612014, indexes in last line changed i->i+1
// F|solution > = |source > , |solution > = F^-1 |source>
void Linalg::MultishiftCG(Linalg::MatrixByVectorFunc F, void *args, const Lattice &lat, const SpinorField &source,
                          int nsys, VECTOR<SpinorField> &solution, VECTOR<double> sigma,
                          double tol, int imax, double &max_err, int &num_iter)
{
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
}
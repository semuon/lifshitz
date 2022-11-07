#include "common.h"
#include "option_parser.h"
#include "utils.h"
#include "parameters.h"
#include "vector_field.h"
#include "spinor_field.h"
#include "scalar_field.h"
#include "scalar_model.h"
#include "matrix.h"
#include "lattice.h"
#include "formats.h"
#include "linalg.h"
#include <time.h>
#include <mpi_module.h>

#include "finite_diff.h"

// Don't use entire namespace std
using std::cout;
using std::endl;
using std::string;

typedef struct PhysicalParams_struct
{
  double m2;
  double invM2;
  double Z;
  double lambda;

  SHARED_PTR<FiniteDifference<int>> laplace_ptr;
  SHARED_PTR<FiniteDifference<int>> laplace_sqr_ptr;
} tPhysicalParams;

template <typename T> bool main_IsFinite(const T value)
{
  return std::isfinite(value);
}

template <> bool main_IsFinite<t_complex>(const t_complex value)
{
  return std::isfinite(value.real()) && std::isfinite(value.imag());
}

bool main_IsReal(const t_complex value)
{
  const double c_treshold = 1e-11;
  const double im = value.imag();

  return std::isfinite(im) && (std::abs(im) < c_treshold);
}

void main_CreateLatticeOperators(tPhysicalParams &params, const uint ndim, const int n_stencil_points)
{
  const uint op_dim = 1;

  auto fwd1d = FiniteDifference<int>::MakeOneSidedDiff(op_dim, n_stencil_points);
  auto bwd1d = FiniteDifference<int>::MakeOneSidedDiff(op_dim, -n_stencil_points);

  params.laplace_ptr = FiniteDifference<int>::MakeLaplacian(ndim, *fwd1d, *bwd1d);
  params.laplace_sqr_ptr = FiniteDifference<int>::ComposeOperators(*params.laplace_ptr, *params.laplace_ptr);
}

template <typename T> T main_LatticeOp(const PhysicalParams_struct &params, const VECTOR<double> &p, const T epsilon)
{
  THROW_EXCEPTION_VERB(std::invalid_argument, "METHOD IS NOT IMPLEMENTED");

  return 0;
}

template <> t_complex main_LatticeOp(const PhysicalParams_struct &params, const VECTOR<double> &p, const t_complex epsilon)
{
  double m2 = params.m2;
  double invM2 = params.invM2;
  double Z = params.Z;

  auto lap = params.laplace_ptr;
  auto lap2 = params.laplace_sqr_ptr;

  t_complex val = invM2 * lap2->Eval(p) - Z * lap->Eval(p);

  val += m2 + 2.0 * epsilon;

  return val;
}

template <> double main_LatticeOp(const PhysicalParams_struct &params, const VECTOR<double> &p, const double epsilon)
{
  double m2 = params.m2;
  double invM2 = params.invM2;
  double Z = params.Z;

  auto lap = params.laplace_ptr;
  auto lap2 = params.laplace_sqr_ptr;

  t_complex c_val = invM2 * lap2->Eval(p) - Z * lap->Eval(p);

  if (!main_IsReal(c_val))
  {
    pStdLogs.Write("IM: %2.15le\n", c_val.imag());
    THROW_EXCEPTION_VERB(std::invalid_argument, "main_LatticeOp<doule>: stencil operator is not real");
  }

  double val = c_val.real() + m2 + 2.0 * epsilon;

  return val;
}

template <typename T> void main_LatticePropAndDerivatives(VECTOR<T> &res, const uint nderiv, const PhysicalParams_struct &params, const Lattice &lat, const T epsilon)
{
  uint ndim = lat.Dim();
  uint vol = lat.Volume();

  res.clear();
  res.resize(nderiv, 0);

  VECTOR<double> p(ndim);
  VECTOR<int> x(ndim);

  bool stop = false;

  for(uint idx = 0; idx < vol; idx++)
  {
    lat.SiteCoordinates(x, idx);

    for(uint i = 0; i < ndim; i++)
    {
      uint L = lat.LatticeSize(i);
      p[i] = 2.0 * M_PI * (double) x[i] / (double) L;
    }

    T latop = main_LatticeOp(params, p, epsilon);

    int num = 1;
    T denum = 1.0;
    for(uint i = 0; i < nderiv; i++)
    {
      denum *= latop;

      res[i] += (double)num / (denum * (double)vol);

      if ( !main_IsFinite(res[i]) ) stop = true;

      num *= -2 * (i + 1);
    }

    if (stop) break;
  }
}

// This always returns 3 derivatives
template <typename T> void main_ContinuumPropAndDerivatives(VECTOR<T> &res, const PhysicalParams_struct &params, const T epsilon)
{
  res.clear();
  res.resize(3);

  double m2 = params.m2;
  double M = sqrt(1 / params.invM2);
  double Z = params.Z;

  T mstar = sqrt(m2 + 2.0 * epsilon);
  T mmodul = sqrt(M * (2.0 * mstar + Z * M));

  res[0] = M * M / (4.0 * M_PI * mmodul);
  res[1] = -1.0 * M * M * M / (4.0 * M_PI * mmodul * mmodul * mmodul * mstar);
  res[2] = 3.0 * M * M * M * M / (4.0 * M_PI * mmodul * mmodul * mmodul * mmodul * mmodul * mstar * mstar) + M * M * M / (4.0 * M_PI * mmodul * mmodul * mmodul * mstar * mstar * mstar);
}

// mplus, mminus, mmodul, kmodul
template <typename T> void main_ContinuumPoles(VECTOR<T> &res, const PhysicalParams_struct &params, const T epsilon)
{
  res.clear();
  res.resize(4);

  double m2 = params.m2;
  double M = sqrt(1 / params.invM2);
  double Z = params.Z;

  double mstar = sqrt(m2 + 2.0 * epsilon);
  double gamma = sqrt(1.0 - 4.0 * mstar * mstar / (Z * Z * M * M));

  double mp = sqrt(0.5 * Z * M * M * (1.0 + gamma)); // m+
  double mm = sqrt(0.5 * Z * M * M * (1.0 - gamma)); // m-
  double m0 = 0.5 * sqrt(M * (2.0 * mstar + Z * M)); // m0
  double k0 = 0.5 * sqrt(M * (2.0 * mstar - Z * M)); // k0

  bool is_symm_phase = (4.0 * mstar * mstar < Z * Z * M * M);

  res[0] = (is_symm_phase) ? mp : 0;
  res[1] = (is_symm_phase) ? mm : 0;
  res[2] = (is_symm_phase) ? 0 : m0;
  res[3] = (is_symm_phase) ? 0 : k0;
}

bool main_ContinuumBisection(double &epsilon, uint &iters, const PhysicalParams_struct &params, const double tol, const uint maxiter)
{
  double m2 = params.m2;
  double M = sqrt(1 / params.invM2);
  double Z = params.Z;
  double lambda = params.lambda;

  VECTOR<double> derivs(3);

  double fval0;
  double fval1;
  double fval2;

  double eps0;
  double eps1;
  double eps2;

  bool converged = false;
  bool stop = false;

  double eps_min = -m2 / 2.0;
  if (Z < 0)
  {
    double eps_min2 = (Z * Z * M * M - 4.0 * m2) / 8.0;
    eps_min = MAX(eps_min, eps_min2);
  }

  eps0 = eps_min + 10 * tol;
  eps1 = MAX(eps0, 1.0);
  eps2 = MAX(eps0, 1.0);

  // Find brackets
  for(uint iter = 0; iter < maxiter; iter++)
  {
    main_ContinuumPropAndDerivatives(derivs, params, eps1);
    fval1 = lambda * derivs[0] / 2.0;
    main_ContinuumPropAndDerivatives(derivs, params, eps2);
    fval2 = lambda * derivs[0] / 2.0;

    if ( !main_IsFinite(fval1) || !main_IsFinite(fval2) )
    {
      stop = true;
      break;
    }

    if (fval1 > eps1 && fval2 < eps2)
      break;

    if (fval1 < eps1)
      eps1 -= (eps1 - eps_min) / 2.0;

    if (fval2 > eps2)
      eps2 *= 2.0;
  }

  if (!stop)
  {
    for(uint iter = 0; iter < maxiter; iter++)
    {
      iters = iter + 1;

      eps0 = eps1 + (eps2 - eps1) / 2.0;

      main_ContinuumPropAndDerivatives(derivs, params, eps0);
      fval0 = lambda * derivs[0] / 2.0;
      main_ContinuumPropAndDerivatives(derivs, params, eps1);
      fval1 = lambda * derivs[0] / 2.0;
      main_ContinuumPropAndDerivatives(derivs, params, eps2);
      fval2 = lambda * derivs[0] / 2.0;

      if ( !main_IsFinite(fval0) || !main_IsFinite(fval1) || !main_IsFinite(fval2) )
        break;

      if (std::abs(eps2 - eps1) < tol && std::abs(fval1 - eps1) < tol)
      {
        epsilon = eps1;
        converged = true;
        break;
      }

      if (fval0 < eps0)
        eps2 = eps0;
      else
        eps1 = eps0;
    }
  }

  return converged;
}

void main_LatticeDispersionDerivatives(VECTOR<double> &res, const double p, const uint ndim, const PhysicalParams_struct &params)
{
  res.clear();
  res.resize(3);

  double invM2 = params.invM2;
  double Z = params.Z;

  // double K1 = 4.0 * invM2;
  // double K2 = -8.0 * ((double)ndim * invM2 + Z / 3.0);
  // double K3 = Z / 6.0;
  // double K4 = 4.0 * (double)ndim * (double)ndim * invM2 + 5.0 * (double)ndim * Z / 2.0;

  // res[0] = K1 * (cos(p) + 2.0) * (cos(p) + 2.0) + K2 * (cos(p) + 2.0) + K3 * (cos(2.0 * p) + 2.0) + K4;
  // res[1] = -2.0 * K1 * (cos(p) + 2.0) * sin(p) - K2 * sin(p) - 2.0 * K3 * sin(2.0 * p);
  // res[2] = 2.0 * K1 * sin(p) * sin(p) - 2.0 * K1 * (cos(p) + 2.0) * cos(p) - K2 * cos(p) - 4.0 * K3 * cos(2.0 * p);

  VECTOR<uint> order(ndim, 0);
  VECTOR<double> pvec(ndim, 0);

  pvec[0] = p;

  auto lap = params.laplace_ptr;
  auto lap2 = params.laplace_sqr_ptr;

  t_complex val = invM2 * lap2->EvalDerivative(order, pvec) - Z * lap->EvalDerivative(order, pvec);
  if (!main_IsReal(val))
  {
    pStdLogs.Write("IM: %2.15le\n", val.imag());
    THROW_EXCEPTION_VERB(std::invalid_argument, "main_LatticeDispersionDerivatives: stencil operator is not real");
  }
  res[0] = val.real();

  order[0] = 1;
  val = invM2 * lap2->EvalDerivative(order, pvec) - Z * lap->EvalDerivative(order, pvec);
  if (!main_IsReal(val))
  {
    pStdLogs.Write("IM: %2.15le\n", val.imag());
    THROW_EXCEPTION_VERB(std::invalid_argument, "main_LatticeDispersionDerivatives: 1st derivative stencil operator is not real");
  }
  res[1] = val.real();

  order[0] = 2;
  val = invM2 * lap2->EvalDerivative(order, pvec) - Z * lap->EvalDerivative(order, pvec);
  if (!main_IsReal(val))
  {
    pStdLogs.Write("IM: %2.15le\n", val.imag());
    THROW_EXCEPTION_VERB(std::invalid_argument, "main_LatticeDispersionDerivatives: 2nd derivative stencil operator is not real");
  }
  res[2] = val.real();

  //double K1 = 1.0 / (M * M) + Z / 12.0;
  //double K2 = -Z;

  //res[0] = 2.0 * K1 * (cos(2.0 * p) - 4.0 * cos(p) + 3.0) + 2.0 * K2 * (cos(p) - 1.0);
  //res[1] = 2.0 * K1 * (-2.0 * sin(2.0 * p) + 4.0 * sin(p)) - 2.0 * K2 * sin(p);
  //res[2] = 2.0 * K1 * (-4.0 * cos(2.0 * p) + 4.0 * cos(p)) - 2.0 * K2 * cos(p);
}

bool main_LatticeDispersionMinimum(VECTOR<double> &res, const double tol, const uint maxiter, const uint ndim, const PhysicalParams_struct &params)
{
  res.clear();
  res.resize(2);

  double p0;
  double p1 = 0;
  double p2 = M_PI;

  double pmin = 0;
  double fmin = 0;

  double e0;
  double e1;
  double e2;

  double f1;
  double f2;

  bool stop = false;
  bool converged = false;

  VECTOR<double> derivs;

  main_LatticeDispersionDerivatives(derivs, p1, ndim, params);
  f1 = derivs[0];
  e1 = derivs[2];
  main_LatticeDispersionDerivatives(derivs, p2, ndim, params);
  f2 = derivs[0];
  e2 = derivs[2];

  if (e1 >= 0 && e2 < 0)
  {
    pmin = p1;
    fmin = f1;
    converged = true;
  }
  else if (e1 < 0 && e2 >= 0)
  {
    pmin = p2;
    fmin = f2;
    converged = true;
  }
  else
  {
    // Find brackets
    p1 = M_PI / 2.0;
    p2 = M_PI / 2.0;
    stop = true;
    for(uint iter = 0; iter < maxiter; iter++)
    {
      main_LatticeDispersionDerivatives(derivs, p1, ndim, params);
      e1 = derivs[1];
      main_LatticeDispersionDerivatives(derivs, p2, ndim, params);
      e2 = derivs[1];

      if (e1 * e2 < 0)
      {
        stop = false;
        break;
      }
      else
      {
        p1 = p1 / 2.0;
        p2 = M_PI - p1;
      }
    }

    if (!stop)
    {
      for(uint iter = 0; iter < maxiter; iter++)
      {
        p0 = p1 + (p2 - p1) / 2.0;

        main_LatticeDispersionDerivatives(derivs, p0, ndim, params);
        e0 = derivs[1];
        main_LatticeDispersionDerivatives(derivs, p1, ndim, params);
        f1 = derivs[0];
        e1 = derivs[1];
        main_LatticeDispersionDerivatives(derivs, p2, ndim, params);
        e2 = derivs[1];

        if (std::abs(p2 - p1) < tol && std::abs(e1) < tol)
        {
          pmin = p1;
          fmin = f1;
          converged = true;
          break;
        }

        if (e0 * e1 < 0)
          p2 = p0;
        else
          p1 = p0;
      }
    }
  }

  if (converged)
  {
    res[0] = pmin;
    res[1] = fmin;
  }
  else
  {
    res[0] = 0;
    res[0] = 0;
  }

  return converged;
}

bool main_LatticeBisection(double &epsilon, double &pmin, uint &iters, const PhysicalParams_struct &params, const Lattice &lat, const double tol, const uint maxiter)
{
  uint ndim = lat.Dim();

  double m2 = params.m2;
  //double M = sqrt(1 / params.invM2);
  //double Z = params.Z;
  double lambda = params.lambda;

  const uint nderivs = 1;
  VECTOR<double> derivs(nderivs);

  double fval0;
  double fval1;
  double fval2;

  double eps0;
  double eps1;
  double eps2;

  bool converged = false;
  bool stop = false;

  iters = 0;

  VECTOR<double> disper_min(2);
  stop = !main_LatticeDispersionMinimum(disper_min, tol, maxiter, ndim, params);

  if (stop)
    return false;

  double eps_min = -1.0 * (m2 + disper_min[1]) / 2.0;

  pmin = disper_min[0];

  eps0 = eps_min + 10 * tol;
  eps1 = MAX(eps0, 1.0);
  eps2 = MAX(eps0, 1.0);

  // Find brackets
  for(uint iter = 0; iter < maxiter; iter++)
  {
    main_LatticePropAndDerivatives(derivs, nderivs, params, lat, eps1);
    fval1 = lambda * derivs[0] / 2.0;
    main_LatticePropAndDerivatives(derivs, nderivs, params, lat, eps2);
    fval2 = lambda * derivs[0] / 2.0;

    if ( !main_IsFinite(fval1) || !main_IsFinite(fval2) )
    {
      stop = true;
      break;
    }

    if (fval1 > eps1 && fval2 < eps2)
      break;

    if (fval1 < eps1)
      eps1 -= (eps1 - eps_min) / 2.0;

    if (fval2 > eps2)
      eps2 *= 2.0;
  }

  if (!stop)
  {
    for(uint iter = 0; iter < maxiter; iter++)
    {
      iters = iter + 1;

      eps0 = eps1 + (eps2 - eps1) / 2.0;

      main_LatticePropAndDerivatives(derivs, nderivs, params, lat, eps0);
      fval0 = lambda * derivs[0] / 2.0;
      main_LatticePropAndDerivatives(derivs, nderivs, params, lat, eps1);
      fval1 = lambda * derivs[0] / 2.0;
      main_LatticePropAndDerivatives(derivs, nderivs, params, lat, eps2);
      fval2 = lambda * derivs[0] / 2.0;

      if ( !main_IsFinite(fval0) || !main_IsFinite(fval1) || !main_IsFinite(fval2) )
        break;

      if (std::abs(eps2 - eps1) < tol && std::abs(fval1 - eps1) < tol)
      {
        epsilon = eps1;
        converged = true;
        break;
      }

      if (fval0 < eps0)
        eps2 = eps0;
      else
        eps1 = eps0;
    }
  }

  return converged;
}

void main_LatticeCorrelator(VECTOR<t_complex> &corr, const PhysicalParams_struct &params, const Lattice &lat, const t_complex epsilon)
{
  uint ndim = lat.Dim();
  uint vol = lat.Volume();

  uint Lmin = lat.LatticeSize(0);
  for(uint i = 1; i < ndim; i++)
  {
    if (Lmin > lat.LatticeSize(i))
      Lmin = lat.LatticeSize(i);
  }

  corr.clear();
  corr.resize(Lmin, 0);

  VECTOR<double> p(ndim);
  VECTOR<int> x(ndim);

  bool stop = false;

  for(uint idx = 0; idx < vol; idx++)
  {
    lat.SiteCoordinates(x, idx);

    for(uint i = 0; i < ndim; i++)
    {
      uint L = lat.LatticeSize(i);
      p[i] = 2.0 * M_PI * (double) x[i] / (double) L;
    }

    t_complex latop = main_LatticeOp(params, p, epsilon);

    for(uint i = 0; i < Lmin; i++)
    {
      corr[i] += exp(t_complex(0, 1) * (double)i * p[0]) / (latop * (double)vol);

      if ( !main_IsFinite(corr[i]) )
      {
        stop = true;
        break;
      }
    }

    if (stop) break;
  }
}

template <typename T> bool main_Newton(T &epsilon, uint &iters, const PhysicalParams_struct &params, const Lattice &lat, const T epsilon0, const tNewtonMethod method, const double tol, const uint maxiter)
{
  epsilon = 0;
  iters = 0;

  double lambda = params.lambda;

  T epsilon1 = epsilon0;
  T epsilon2 = 0;

  T fval = 0;
  T d1fval = 0;
  T d2fval = 0;

  uint nderivs;
  VECTOR<T> derivs;

  bool converged = false;
  bool stop = false;

  switch(method)
  {
    case ITERATIONS_NEWTON:
      nderivs = 2;
      break;
    case ITERATIONS_HALLEY:
      nderivs = 3;
      break;
    default:
      THROW_EXCEPTION(AssertFailure, "Unknown method");
  }

  for(uint iter = 0; iter < maxiter; iter++)
  {
    iters = iter + 1;

    main_LatticePropAndDerivatives(derivs, nderivs, params, lat, epsilon1);

    for(uint i = 0; i < nderivs; i++)
    {
      if ( !main_IsFinite(derivs[i]) )
      {
        stop = true;
        break;
      }
    }

    if (stop) break;

    switch(method)
    {
      case ITERATIONS_NEWTON:
        fval = (lambda / 2.0) * derivs[0] - epsilon1;
        d1fval = (lambda / 2.0) * derivs[1] - 1.0;
        epsilon2 = epsilon1 - fval / d1fval;
        break;
      case ITERATIONS_HALLEY:
        fval = (lambda / 2.0) * derivs[0] - epsilon1;
        d1fval = (lambda / 2.0) * derivs[1] - 1.0;
        d2fval = (lambda / 2.0) * derivs[2];
        epsilon2 = epsilon1 - 2.0 * fval * d1fval / (2.0 * d1fval * d1fval - fval * d2fval);
        break;
      default:
        THROW_EXCEPTION(AssertFailure, "Unknown method");
    }

    if ( !main_IsFinite(epsilon2) ) break;

    if ( std::abs(fval) < tol && std::abs(epsilon2 - epsilon1) < tol )
    {
      converged = true;
      break;
    }

    epsilon1 = epsilon2;
  }

  epsilon = epsilon1;
  return converged;
}

template <typename T> bool main_NewtonContinuum(T &epsilon, uint &iters, const PhysicalParams_struct &params, const T epsilon0, const tNewtonMethod method, const double tol, const uint maxiter)
{
  const uint nderivs = 3;

  epsilon = 0;
  iters = 0;

  double lambda = params.lambda;

  T epsilon1 = epsilon0;
  T epsilon2 = 0;

  VECTOR<T> derivs(3);

  T fval = 0;
  T d1fval = 0;
  T d2fval = 0;

  bool converged = false;
  bool stop = false;

  for(uint iter = 0; iter < maxiter; iter++)
  {
    iters = iter + 1;

    main_ContinuumPropAndDerivatives(derivs, params, epsilon1);

    for(uint i = 0; i < nderivs; i++)
    {
      if ( !main_IsFinite(derivs[i]) )
      {
        stop = true;
        break;
      }
    }

    if (stop) break;

    switch(method)
    {
      case ITERATIONS_NEWTON:
        fval = (lambda / 2.0) * derivs[0] - epsilon1;
        d1fval = (lambda / 2.0) * derivs[1] - 1.0;
        epsilon2 = epsilon1 - fval / d1fval;
        break;
      case ITERATIONS_HALLEY:
        fval = (lambda / 2.0) * derivs[0] - epsilon1;
        d1fval = (lambda / 2.0) * derivs[1] - 1.0;
        d2fval = (lambda / 2.0) * derivs[2];
        epsilon2 = epsilon1 - 2.0 * fval * d1fval / (2.0 * d1fval * d1fval - fval * d2fval);
        break;
      default:
        THROW_EXCEPTION(AssertFailure, "Unknown method");
    }

    if ( !main_IsFinite(epsilon2) ) break;

    if ( std::abs(fval) < tol && std::abs(epsilon2 - epsilon1) < tol )
    {
      converged = true;
      break;
    }

    epsilon1 = epsilon2;
  }

  epsilon = epsilon1;
  return converged;
}

int main(int argc, char **argv)
{
  const string f_bin_attr = "wb";
  const string f_txt_attr = "w";

  common_AppInit(argc, argv, "Lifshitz regime. Large N homogeneous epsilon");
  int64_t begin = Utils::GetTimeMs64();

  const uint nd = 1;
  const uint nc = 1;

  VECTOR<uint> latdims(pL);
  uint ndim = pDim;

  Lattice lat(latdims, ndim, nc, nd);
  //uint vol = lat.Volume();

  uint Lmin = lat.LatticeSize(0);
  for(uint i = 1; i < ndim; i++)
  {
    if (Lmin > lat.LatticeSize(i))
      Lmin = lat.LatticeSize(i);
  }

  cout << "LMIN: " << Lmin << endl;

  double invM2 = pInvM2;
  double m2 = pm2;
  double Z = pZ;
  double lambda = pLambdaN;
  double n_stencil_pts = pNStencilPts;

  tPhysicalParams params;
  params.lambda = lambda;
  params.invM2 = invM2;
  params.m2 = m2;
  params.Z = Z;
  
  main_CreateLatticeOperators(params, ndim, n_stencil_pts);

  int n_iters = pNiters;
  int n_tries = pNtries;
  double tolerance = pTolerance;
  //double relax_alpha = pRelaxAlpha;
  //double random_range = pRandomRange;
  //tNewtonMethod method = pMethod;

  //t_complex epsilon0;
  //t_complex epsilon1;
  double epsilon1;
  double pmin;
  double epsilon_cont;

  VECTOR<t_complex> corr;
  VECTOR<double> poles;

  const string history_name = "history.txt";
 
  FILE *f_history = NULL;

  if (ModuleMPI::IsMasterNode())
  {
    f_history = pDataDir.OpenFile(history_name, f_txt_attr);
  }

  // // Testing HMC action
  // const uint n = 1;
  // tScalarModelParams hmc_params;

  // const double eps0 = 0;

  // hmc_params.N = n;
  // hmc_params.Z = Z;
  // hmc_params.invM2 = invM2;
  // hmc_params.m2 = m2;
  // hmc_params.kappa = 0;
  // hmc_params.lambdaN = 0;

  // params.lambda = 0;

  // RealScalarFieldN phi_x(lat, n);
  // ScalarFieldN phi_k(lat, n);

  // for(uint x = 0; x < vol; x++)
  // {
  //   const uint i = 0;

  //   phi_x(x, i) = rand_double(-1.0, 1.0);
  // }

  // VECTOR<int> kvec(ndim);
  // VECTOR<int> xvec(ndim);
  // VECTOR<double> d_kvec(ndim);

  // for(uint k = 0; k < vol; k++)
  // {
  //   const uint i = 0;

  //   lat.SiteCoordinates(kvec, k);

  //   phi_k(k, i) = 0;

  //   for(uint x = 0; x < vol; x++)
  //   {
  //     lat.SiteCoordinates(xvec, x);

  //     double phase = 0;
  //     for(uint mu = 0; mu < ndim; mu++)
  //       phase += 2.0 * M_PI * (double)kvec[mu] * (double)xvec[mu] / (double)lat.LatticeSize(mu);

  //     phi_k(k, i) += phi_x(x, i) * exp(I * phase) / sqrt(vol);
  //   }
  // }

  // // Momentum space action
  // double action_k = 0;
  // for(uint k = 0; k < vol; k++)
  // {
  //   const uint i = 0;

  //   double abs_phi = std::abs(phi_k(k, i));

  //   lat.SiteCoordinates(kvec, k);

  //   for(uint mu = 0; mu < ndim; mu++)
  //     d_kvec[mu] = 2.0 * M_PI * (double)kvec[mu] / (double)lat.LatticeSize(mu);

  //   action_k += 0.5 * abs_phi * abs_phi * main_LatticeOp(params, d_kvec, eps0);
  // }

  // cout << "K SPACE ACTION: " << action_k << endl;

  // // Coordinate space action
  // double action_x = ScalarModel::Action(hmc_params, phi_x);

  // cout << "X SPACE ACTION: " << action_x << endl;
  // cout << "DIFF: " << action_x - action_k << endl;

  // exit(0);

  //t_complex eps = 0.133;
  //VECTOR<double> ppp(3);
  //ppp[0] = 1.0;
  //ppp[1] = 0.3;
  //ppp[2] = -0.15;
  //t_complex opop = main_LatticeOp(params, ppp, eps);

  //cout << "OP: " << opop << endl;

  //exit(0);

  // VECTOR<double> derivs;
  // uint npts = 3000;
  // FILE *f_tmp = pDataDir.OpenFile("litter.txt", "w");
  // bool converged = main_LatticeDispersionMinimum(derivs, tolerance, n_iters, ndim, params);
  // pmin = derivs[0];
  // double fmin = derivs[1];
  // double epsmin = -(params.m2 + fmin) / 2.0;
  // //double epsmin = -3.0;
  // for(uint i = 0; i < npts + 1; i++)
  // {
  //   double eps = epsmin+ 2.0 * (double) i / (double)npts;

  //   main_LatticePropAndDerivatives(derivs, 1, params, lat, eps);
  //   SAFE_FPRINTF(f_tmp, "%2.15le\t%2.15le\n", eps, derivs[0]);
  // }
  // fclose(f_tmp);

  // exit(0);

  for(int i_try = 0; i_try < n_tries; i_try++)
  {
    uint iters;
    uint iters_cont;

    //epsilon0 = rand_double(-random_range, random_range) + I * rand_double(-random_range, random_range);

    //bool is_converged = main_Newton(epsilon1, iters, params, lat, epsilon0, method, tolerance, n_iters);
    bool is_converged = main_LatticeBisection(epsilon1, pmin, iters, params, lat, tolerance, n_iters);
    //bool is_converged_cont = main_NewtonContinuum(epsilon_cont, iters_cont, params, epsilon0.real(), method, tolerance, n_iters);
    bool is_converged_cont = main_ContinuumBisection(epsilon_cont, iters_cont, params, tolerance, n_iters);

    //t_complex action = 0;

    if (is_converged)
    {
      main_LatticeCorrelator(corr, params, lat, epsilon1);
    }
    else
    {
      corr.clear();
      corr.resize(Lmin, 0);
    }

    if (is_converged_cont)
    {
      main_ContinuumPoles(poles, params, epsilon_cont);
    }
    else
    {
      poles.clear();
      poles.resize(4, 0);
    }

    SAFE_FPRINTF(f_history, "%s\t%d\t%2.15le\t%2.15le\t%2.15le\t%2.15le", (is_converged) ? "Y" : "N", iters, epsilon1, 0.0, pmin, 0.0);
    SAFE_FPRINTF(f_history, "\t%s\t%d\t%2.15le\t%2.15le", (is_converged_cont) ? "Y" : "N", iters_cont, epsilon_cont, 0.0);
    for(uint i = 0; i < 4; i++)
    {
      SAFE_FPRINTF(f_history, "\t%2.15le", poles[i]);
    }
    for(uint i = 0; i < Lmin; i++)
    {
      SAFE_FPRINTF(f_history, "\t%2.15le\t%2.15le", corr[i].real(), corr[i].imag());
    }
    SAFE_FPRINTF(f_history, "\n");

    if (is_converged)
    {
      pStdLogs.Write("Converged at try %d within %d iterations\n", i_try, iters);

      //for(uint x = 0; x < vol; x++)
      //  action0 -= epsilon0(x) * epsilon0(x) / (double) vol;
      //for(uint n = 0; n < vol; n++)
      //  action1 += lambda * log((t_complex)evals[n]) / 2.0;
      //action = action0 + action1;
    }
    else
    {
      pStdLogs.Write("Failed to converge at try %d\n", i_try);
    }
    
  }

  if (ModuleMPI::IsMasterNode())
  {
    fclose(f_history);
  }

  int64_t end = Utils::GetTimeMs64();

  cout << endl;
  cout << "Program realized in " << (1.0 * end - 1.0 * begin) / 1000. << " s" << endl;

  common_AppFin();

  return 0;
}

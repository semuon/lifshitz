#include "common.h"
#include "option_parser.h"
#include "utils.h"
#include "parameters.h"
#include "vector_field.h"
#include "spinor_field.h"
#include "scalar_field.h"
#include "matrix.h"
#include "lattice.h"
#include "formats.h"
#include "linalg.h"
#include <time.h>
#include <mpi_module.h>

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
} tPhysicalParams;

template <typename T> bool main_IsFinite(const T value)
{
  return std::isfinite(value);
}

template <> bool main_IsFinite<t_complex>(const t_complex value)
{
  return std::isfinite(value.real()) && std::isfinite(value.imag());
}

template <typename T> T main_LatticeOp(const PhysicalParams_struct &params, const VECTOR<double> &p, const T epsilon)
{
  uint ndim = p.size();

  double m2 = params.m2;
  double invM2 = params.invM2;
  double Z = params.Z;

  T val = 0;

  for(uint i = 0; i < ndim; i++)
  {
    double pi = p[i];

    val += 2.0 * (invM2 + Z / 12.0) * ( cos(2.0 * pi) - 4.0 * cos(pi) + 3.0 ) - 
           2.0 * Z * ( cos(pi) - 1.0 );
  }

  val += m2 + 2.0 * epsilon;

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

  double mstar = sqrt(m2 + 2.0 * epsilon.real());
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
      corr[i] += exp(I * (double)i * p[0]) / (latop * (double)vol);

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
  uint vol = lat.Volume();

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

  tPhysicalParams params;
  params.lambda = lambda;
  params.invM2 = invM2;
  params.m2 = m2;
  params.Z = Z;

  int n_iters = pNiters;
  int n_tries = pNtries;
  double tolerance = pTolerance;
  double relax_alpha = pRelaxAlpha;
  double random_range = pRandomRange;
  tNewtonMethod method = pMethod;

  t_complex epsilon0;
  t_complex epsilon1;
  t_complex epsilon_cont;

  VECTOR<t_complex> corr;
  VECTOR<t_complex> poles;

  const string history_name = "history.txt";
 
  FILE *f_history = NULL;

  if (ModuleMPI::IsMasterNode())
  {
    f_history = pDataDir.OpenFile(history_name, f_txt_attr);
  }

  for(int i_try = 0; i_try < n_tries; i_try++)
  {
    uint iters;
    uint iters_cont;

    epsilon0 = rand_double(-random_range, random_range) + I * rand_double(-random_range, random_range);

    double is_converged = main_Newton(epsilon1, iters, params, lat, epsilon0, method, tolerance, n_iters);
    double is_converged_cont = main_NewtonContinuum(epsilon_cont, iters_cont, params, epsilon0, method, tolerance, n_iters);

    t_complex action = 0;

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

    SAFE_FPRINTF(f_history, "%s\t%d\t%2.15le\t%2.15le\t%2.15le\t%2.15le", (is_converged) ? "Y" : "N", iters, epsilon1.real(), epsilon1.imag(), action.real(), action.imag());
    SAFE_FPRINTF(f_history, "\t%s\t%d\t%2.15le\t%2.15le", (is_converged_cont) ? "Y" : "N", iters_cont, epsilon_cont.real(), epsilon_cont.imag());
    for(uint i = 0; i < 4; i++)
    {
      SAFE_FPRINTF(f_history, "\t%2.15le", poles[i].real());
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

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

typedef enum NewtonMethod_enum
{
  ITERATIONS_NEWTON,
  ITERATIONS_HALLEY
} tNewtonMethod;

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

  t_complex solution_action = 0;
  t_complex solution_action0 = 0;
  t_complex solution_action1 = 0;

  uint nderiv = 4;
  VECTOR<double> p(3);
  VECTOR<double> res;
  p[0] = 1.0;
  p[1] = 0.3;
  p[2] = -0.15;
  double eee = 0.133;
  cout << "" << "ACTION: " << main_LatticeOp(params, p, eee) << endl;
  main_LatticePropAndDerivatives(res, nderiv, params, lat, eee);
  for(uint i = 0; i < nderiv; i++)
  {
    cout << "" << "N = " << i << " : " << res[i] << endl;
  }

  t_complex e1 = 0.5 + I * 0.33;
  t_complex e2;
  double uuu;
  uint iters;
  bool converged = main_Newton(e2, iters, params, lat, e1, ITERATIONS_HALLEY, 1e-8, 100);

  cout << "CONVERGED: " << converged << endl;
  cout << "ITERS: " << iters << endl;
  cout << "EPSILON: " << e2 << endl;

  exit(0);

  // ScalarField nonsense(lat);
  // ScalarField nonsense2(lat);
  // ScalarField nonsense3(lat);
  // ScalarField nonsense4(lat);
  // ScalarField nonsense5(lat);
  // nonsense3.Clear();
  // nonsense4.Clear();
  // nonsense5.Clear();
  // for(uint x = 0; x < vol; x++)
  // {
  //   nonsense(x) = rand_double(-10,10);
  //   nonsense2(x) = rand_double(-10,10);
  // }

  // pGlobalProfiler.StartSection("TEST");
  // int crm_iter;
  // double crm_err;
  // params.epsilon = &nonsense;

  // pGlobalProfiler.StartTimer("LAPACK");
  // main_Loperator(Lop, lat, nonsense, params);
  // Linalg::LapackHermitianEigensystem(Lop, evals.data(), vol);
  // //Linalg::LapackInvert(Lop, vol);
  // nonsense3 = Lop * nonsense2;
  // pGlobalProfiler.StopTimer("LAPACK");

  // double minnnnnn = evals[0];
  // for(uint x = 1; x < vol; x++)
  //   if (fabs(evals[x]) < fabs(minnnnnn))
  //     minnnnnn = evals[x];
  // cout << "LAPACK MIN ABS EVAL = " << minnnnnn << endl;

  // try
  // {
  //   pGlobalProfiler.StartTimer("WITH DEFLATION");
  //   main_InvertLoperatorWithDeflation(nonsense4, nonsense2, vol, params, 16, 16, 1000000);
  //   pGlobalProfiler.StopTimer("WITH DEFLATION");
  //   pGlobalProfiler.StartTimer("WITHOUT DEFLATION");
  //   Linalg::CRM(main_LoperatorTimesVec, (void *)&params, vol, nonsense2, nonsense4, 1e-8, 1000000, crm_err, crm_iter);
  //   cout << "CRM ERR = " << crm_err << ", ITER = " << crm_iter << endl;
  //   //main_InvertLoperatorWithDeflation(nonsense4, nonsense2, vol, params, 0, 16, 100000);
  //   pGlobalProfiler.StopTimer("WITH DEFLATION");
  //   //Linalg::CRM(main_LoperatorTimesVec, (void *)&params, vol, nonsense2, nonsense4, 1e-8, 100000, crm_err, crm_iter);
  // }
  // catch(std::exception &exc)
  // {
  //   cout  << "CRM failed to converge: " << exc.what() << endl;
  // }

  // //main_LoperatorTimesVec(nonsense4, nonsense2, (void*)&params);
  // nonsense5 = nonsense4 - nonsense3;
  // cout << endl;
  // cout << "DISCREPANCY = " << nonsense5.Norm() << endl;
  // //cout << "CRM ERR = " << crm_err << ", ITER = " << crm_iter << endl;
  // cout << endl;
  // pGlobalProfiler.EndSection("TEST");

  // pGlobalProfiler.PrintStatistics();

  const string history_name = "history.bin";
 
  FILE *f_history = NULL;

  if (ModuleMPI::IsMasterNode())
  {
    f_history = pDataDir.OpenFile(history_name, f_bin_attr);
  }

  bool is_solution_found = false;

  for(int i_try = 0; i_try < n_tries; i_try++)
  {
    bool is_converged = false;

    for(uint x = 0; x < vol; x++)
    {
      //epsilon0(x) = rand_double(-random_range, random_range);
      //epsilon0(x) = 1.0 * x;
    }

    for(int i_iter = 0; i_iter < n_iters; i_iter++)
    {
      //main_Loperator(Lop, lat, epsilon0, params);
      //FILE *fff = pDataDir.OpenFile("M.txt", "w");
      //Formats::PrintMatrix(fff, Lop);
      //fclose(fff);
      //return 0;

      //Linalg::LapackHermitianEigensystem(Lop, evals.data(), vol);
      //Lop.Transpose();

      for(uint n = 0; n < vol; n++)
      {
        //if (fabs(evals[n]) < c_flt_epsilon)
        {
          //pStdLogs.Write("Zero eigenvalue detected at try #%d, iteration #%d: e[%u] = %2.15le\n", i_try, i_iter, n, evals[n]);
          break;
        }
      }

      for(uint x = 0; x < vol; x++)
      {
        //epsilon1(x) = 0;

        //for(uint n = 0; n < vol; n++)
        //  epsilon1(x) += (double)vol * lambda * norm(Lop(n, x)) / (2.0 * evals[n]);
      }

      //depsilon = epsilon1 - epsilon0;
      //epsilon0 = (1.0 - relax_alpha) * epsilon0 + relax_alpha * epsilon1;

      //double de = depsilon.Norm();
      double de = 0;

      if (isinf(de) || isnan(de))
      {
        pStdLogs.Write("Encoutered numerical nonsense at iteration: %d\n", i_iter);
        break;
      }

      if (de <= tolerance)
      {
        is_converged = true;
        break;
      }
    }

    if (is_converged)
    {
      pStdLogs.Write("Converged at try %d\n", i_try);

      t_complex action = 0;
      t_complex action0 = 0;
      t_complex action1 = 0;

      //for(uint x = 0; x < vol; x++)
      //  action0 -= epsilon0(x) * epsilon0(x) / (double) vol;

      //main_Loperator(Lop, lat, epsilon0, params);
      //Linalg::LapackHermitianEigensystem(Lop, evals.data(), vol);
      //Lop.Transpose();

      //for(uint n = 0; n < vol; n++)
      //  action1 += lambda * log((t_complex)evals[n]) / 2.0;

      //action = action0 + action1;

      //bool is_homogeneous = main_IsHomogeneousVector(epsilon0, c_flt_epsilon);

      //pStdLogs.Write("Found %s solution with the action: (%2.15le, %2.15le)\n", (is_homogeneous) ? "homogeneous": "non-homogeneous", i_try, real(action), imag(action));

      //if (is_homogeneous)
      {
        //t_complex e0 = epsilon0[0];
        //pStdLogs.Write("e0 = (%2.15le, %2.15le), m^2_eff = (%2.15le, %2.15le)\n", real(e0), imag(e0), real(2.0 * e0 + m2), imag(2.0 * e0 + m2));
      }

      SAFE_FWRITE(&action, sizeof(t_complex), 1, f_history);
      SAFE_FWRITE(&action0, sizeof(t_complex), 1, f_history);
      SAFE_FWRITE(&action1, sizeof(t_complex), 1, f_history);
      //Formats::DumpBinary(f_history, epsilon0);

      if (!is_solution_found || (real(solution_action) > real(action)))
      {
        pStdLogs.Write("This solution has a lower action\n");

        //solution = epsilon0;
        //solution_evecs = Lop;
        //solution_evals = evals;

        solution_action = action;
        solution_action0 = action0;
        solution_action1 = action1;

        is_solution_found = true;
      }
    }
    else
    {
      pStdLogs.Write("Failed to converge at try %d\n", i_try);
    }
    
  }

  const string loperator_name = "Loperator.bin";
  const string solution_name = "solution.bin";
 
  FILE *f_loperator = NULL;
  FILE *f_solution = NULL;

  if (ModuleMPI::IsMasterNode())
  {
    fclose(f_history);

    f_loperator = pDataDir.OpenFile(loperator_name, f_bin_attr);
    f_solution = pDataDir.OpenFile(solution_name, f_bin_attr);

    //Formats::DumpBinary(f_loperator, solution_evals);
    //Formats::DumpBinary(f_loperator, solution_evecs);

    SAFE_FWRITE(&solution_action, sizeof(t_complex), 1, f_solution);
    SAFE_FWRITE(&solution_action0, sizeof(t_complex), 1, f_solution);
    SAFE_FWRITE(&solution_action1, sizeof(t_complex), 1, f_solution);
    //Formats::DumpBinary(f_solution, solution);

    fclose(f_loperator);
    fclose(f_solution);
  }

  int64_t end = Utils::GetTimeMs64();

  cout << endl;
  cout << "Program realized in " << (1.0 * end - 1.0 * begin) / 1000. << " s" << endl;

  common_AppFin();

  return 0;
}
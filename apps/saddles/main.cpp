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
  ScalarField *epsilon;
} tPhysicalParams;

void main_Loperator(Matrix &Lop, const Lattice &lat, const ScalarField &epsilon, const tPhysicalParams &params)
{
  int ndim = lat.Dim();
  uint vol = lat.Volume();

  double m2 = params.m2;
  double invM2 = params.invM2;
  double Z = params.Z;

  ASSERT(Lop.Ncols() == vol);
  ASSERT(Lop.IsSquareMatrix());
  ASSERT(epsilon.Count() == vol);

  Lop.Clear();

  for(uint x = 0; x < vol; x++)
  {
    Lop(x, x) += m2 + 2.0 * epsilon(x) - 2.0 * ndim * Z + 4.0 * ndim * ndim * invM2;

    for(int mu = 0; mu < ndim; mu++)
    {
      uint xmu_fwd = lat.SiteIndexForward(x, mu);
      uint xmu_bwd = lat.SiteIndexBackward(x, mu);

      double val = Z - 4.0 * ndim * invM2;
      Lop(xmu_fwd, x) += val;
      Lop(xmu_bwd, x) += val;

      for(int nu = 0; nu < ndim; nu++)
      {
        uint xmu_fwd_nu_fwd = lat.SiteIndexForward(xmu_fwd, nu);
        uint xmu_fwd_nu_bwd = lat.SiteIndexBackward(xmu_fwd, nu);
        uint xmu_bwd_nu_fwd = lat.SiteIndexForward(xmu_bwd, nu);
        uint xmu_bwd_nu_bwd = lat.SiteIndexBackward(xmu_bwd, nu);

        Lop(xmu_fwd_nu_fwd, x) += invM2;
        Lop(xmu_fwd_nu_bwd, x) += invM2;
        Lop(xmu_bwd_nu_fwd, x) += invM2;
        Lop(xmu_bwd_nu_bwd, x) += invM2;
      }
    }
  }
}

void main_LoperatorTimesVec(BaseLinearVector &out, const BaseLinearVector &in, void *args)
{
  ASSERT(args != NULL);

  tPhysicalParams &params = *((tPhysicalParams *)args);

  double m2 = params.m2;
  double invM2 = params.invM2;
  double Z = params.Z;
  ScalarField &epsilon = *params.epsilon;

  const Lattice &lat = epsilon.GetLattice();

  int ndim = lat.Dim();
  uint vol = lat.Volume();

  ASSERT(out.Count() == in.Count());
  ASSERT(epsilon.Count() == in.Count());
  ASSERT(in.Count() == vol);

  out.Clear();

  for(uint x = 0; x < vol; x++)
  {
    out[x] += (m2 + 2.0 * epsilon(x) - 2.0 * ndim * Z + 4.0 * ndim * ndim * invM2) * in[x];

    for(int mu = 0; mu < ndim; mu++)
    {
      uint xmu_fwd = lat.SiteIndexForward(x, mu);
      uint xmu_bwd = lat.SiteIndexBackward(x, mu);

      double val = Z - 4.0 * ndim * invM2;
      out[x] += val * (in[xmu_fwd] + in[xmu_bwd]);

      for(int nu = 0; nu < ndim; nu++)
      {
        uint xmu_fwd_nu_fwd = lat.SiteIndexForward(xmu_fwd, nu);
        uint xmu_fwd_nu_bwd = lat.SiteIndexBackward(xmu_fwd, nu);
        uint xmu_bwd_nu_fwd = lat.SiteIndexForward(xmu_bwd, nu);
        uint xmu_bwd_nu_bwd = lat.SiteIndexBackward(xmu_bwd, nu);

        out[x] += invM2 * (in[xmu_fwd_nu_fwd] + in[xmu_fwd_nu_bwd] + in[xmu_bwd_nu_fwd] + in[xmu_bwd_nu_bwd]);
      }
    }
  }
}

int main(int argc, char **argv)
{
  common_AppInit(argc, argv, "Lifshitz regime");
  int64_t begin = Utils::GetTimeMs64();

  const uint nd = 1;
  const uint nc = 1;

  VECTOR<uint> latdims(pL);
  uint ndim = pDim;

  Lattice lat(latdims, ndim, nc, nd);
  uint vol = lat.Volume();

  Matrix Lop(vol);
  Matrix solution_evecs(vol);

  double invM2 = pInvM2;
  double m2 = pm2;
  double Z = pZ;
  double lambda = pLambda;

  ScalarField epsilon0(lat);
  ScalarField epsilon1(lat);
  ScalarField solution(lat);
  ScalarField depsilon(lat);

  tPhysicalParams params;
  params.lambda = lambda;
  params.invM2 = invM2;
  params.m2 = m2;
  params.Z = Z;
  params.epsilon = &epsilon0;

  VECTOR<double> evals(vol);
  VECTOR<double> solution_evals(vol);

  Lop.Clear();
  solution_evecs.Clear();
  epsilon0.Clear();
  epsilon1.Clear();
  solution.Clear();
  depsilon.Clear();

  int n_iters = pNiters;
  int n_tries = pNtries;
  double tolerance = pTolerance;
  double relax_alpha = pRelaxAlpha;
  double random_range = pRandomRange;

  t_complex solution_action = 0;
  t_complex solution_action0 = 0;
  t_complex solution_action1 = 0;

  ScalarField nonsense(lat);
  ScalarField nonsense2(lat);
  ScalarField nonsense3(lat);
  ScalarField nonsense4(lat);
  ScalarField nonsense5(lat);
  nonsense3.Clear();
  nonsense4.Clear();
  nonsense5.Clear();
  for(uint x = 0; x < vol; x++)
  {
    nonsense(x) = rand_double(-10,10);
    nonsense2(x) = rand_double(-10,10);
  }

  pGlobalProfiler.StartSection("TEST");
  int crm_iter;
  double crm_err;
  params.epsilon = &nonsense;

  pGlobalProfiler.StartTimer("LAPACK");
  main_Loperator(Lop, lat, nonsense, params);
  Linalg::LapackInvert(Lop, vol);
  nonsense3 = Lop * nonsense2;
  pGlobalProfiler.StopTimer("LAPACK");

  try
  {
    Linalg::CRM(main_LoperatorTimesVec, (void *)&params, vol, nonsense2, nonsense4, 1e-8, 100000, crm_err, crm_iter);
  }
  catch(std::exception &exc)
  {
    cout  << "CRM failed to converge: " << exc.what() << endl;
  }

  //main_LoperatorTimesVec(nonsense4, nonsense2, (void*)&params);
  nonsense5 = nonsense4 - nonsense3;
  cout << endl;
  cout << "DISCREPANCY = " << nonsense5.Norm() << endl;
  cout << "CRM ERR = " << crm_err << ", ITER = " << crm_iter << endl;
  cout << endl;
  pGlobalProfiler.EndSection("TEST");

  pGlobalProfiler.PrintStatistics();

  bool is_solution_found = false;

  for(int i_try = 0; i_try < n_tries; i_try++)
  {
    bool is_converged = false;

    for(uint x = 0; x < vol; x++)
    {
      epsilon0(x) = rand_double(-random_range, random_range);
      //epsilon0(x) = 1.0 * x;
    }

    for(int i_iter = 0; i_iter < n_iters; i_iter++)
    {
      main_Loperator(Lop, lat, epsilon0, params);
      Linalg::LapackHermitianEigensystem(Lop, evals.data(), vol);
      Lop.Transpose();

      for(uint n = 0; n < vol; n++)
      {
        if (abs(evals[n]) < __FLT_EPSILON__)
        {
          pStdLogs.Write("Zero eigenvalue detected at try #%d, iteration #%d: e[%u] = %2.15le\n", i_try, i_iter, n, evals[n]);
          break;
        }
      }

      for(uint x = 0; x < vol; x++)
      {
        epsilon1(x) = 0;

        for(uint n = 0; n < vol; n++)
          epsilon1(x) += lambda * norm(Lop(n, x)) / (2.0 * evals[n]);
      }

      depsilon = epsilon1 - epsilon0;
      epsilon0 = (1.0 - relax_alpha) * epsilon0 + relax_alpha * epsilon1;

      double de = depsilon.Norm();

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
      t_complex action = 0;
      t_complex action0 = 0;
      t_complex action1 = 0;

      for(uint x = 0; x < vol; x++)
        action0 -= epsilon0(x) * epsilon0(x) / (double) vol;

      main_Loperator(Lop, lat, epsilon0, params);
      Linalg::LapackHermitianEigensystem(Lop, evals.data(), vol);
      Lop.Transpose();

      for(uint n = 0; n < vol; n++)
        action1 += lambda * log((t_complex)evals[n]) / 2.0;

      action = action0 + action1;

      if (!is_solution_found || (real(solution_action) > real(action)))
      {
        solution = epsilon0;
        solution_evecs = Lop;
        solution_evals = evals;

        solution_action = action;
        solution_action0 = action0;
        solution_action1 = action1;

        pStdLogs.Write("Found solution at try #%d with the action: (%2.15le, %2.15le)\n", i_try, real(action), imag(action));

        is_solution_found = true;
      }
    }
  }

  const string loperator_name = "Loperator.bin";
  const string solution_name = "solution.bin";
 
  FILE *f_loperator = NULL;
  FILE *f_solution = NULL;

  if (ModuleMPI::IsMasterNode())
  {
    const string f_bin_attr = "wb";
    const string f_txt_attr = "w";

    f_loperator = pDataDir.OpenFile(loperator_name, f_bin_attr);
    f_solution = pDataDir.OpenFile(solution_name, f_bin_attr);

    Formats::DumpBinary(f_loperator, solution_evals);
    Formats::DumpBinary(f_loperator, solution_evecs);

    SAFE_FWRITE(&solution_action, sizeof(t_complex), 1, f_solution);
    SAFE_FWRITE(&solution_action0, sizeof(t_complex), 1, f_solution);
    SAFE_FWRITE(&solution_action1, sizeof(t_complex), 1, f_solution);
    Formats::DumpBinary(f_solution, solution);

    fclose(f_loperator);
    fclose(f_solution);
  }

  int64_t end = Utils::GetTimeMs64();

  cout << endl;
  cout << "Program realized in " << (1.0 * end - 1.0 * begin) / 1000. << " s" << endl;

  common_AppFin();

  return 0;
}

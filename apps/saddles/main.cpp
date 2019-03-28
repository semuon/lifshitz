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
  double M2;
  double Z;
  double lambda;
} tPhysicalParams;

void main_Loperator(Matrix &Lop, const Lattice &lat, const ScalarField &epsilon, const tPhysicalParams &params)
{
  int ndim = lat.Dim();
  uint vol = lat.Volume();

  double m2 = params.m2;
  double M2 = params.M2;
  double Z = params.Z;

  ASSERT(Lop.Ncols() == vol);
  ASSERT(Lop.IsSquareMatrix());
  ASSERT(epsilon.Count() == vol);

  Lop.Clear();

  for(uint x = 0; x < vol; x++)
  {
    Lop(x, x) += m2 + 2.0 * epsilon(x) - 2.0 * ndim * Z + 4.0 * ndim * ndim * M2;

    for(int mu = 0; mu < ndim; mu++)
    {
      uint xmu_fwd = lat.SiteIndexForward(x, mu);
      uint xmu_bwd = lat.SiteIndexBackward(x, mu);

      double val = Z - 4.0 * ndim * M2;
      Lop(xmu_fwd, x) += val;
      Lop(xmu_bwd, x) += val;

      for(int nu = 0; nu < ndim; nu++)
      {
        uint xmu_fwd_nu_fwd = lat.SiteIndexForward(xmu_fwd, nu);
        uint xmu_fwd_nu_bwd = lat.SiteIndexBackward(xmu_fwd, nu);
        uint xmu_bwd_nu_fwd = lat.SiteIndexForward(xmu_bwd, nu);
        uint xmu_bwd_nu_bwd = lat.SiteIndexBackward(xmu_bwd, nu);

        Lop(xmu_fwd_nu_fwd, x) += M2;
        Lop(xmu_fwd_nu_bwd, x) += M2;
        Lop(xmu_bwd_nu_fwd, x) += M2;
        Lop(xmu_bwd_nu_bwd, x) += M2;
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

  double M2 = pM2;
  double m2 = pm2;
  double Z = pZ;
  double lambda = pLambda;

  tPhysicalParams params;
  params.lambda = lambda;
  params.M2 = M2;
  params.m2 = m2;
  params.Z = Z;

  ScalarField epsilon0(lat);
  ScalarField epsilon1(lat);
  ScalarField solution(lat);
  ScalarField depsilon(lat);

  VECTOR<double> evals(vol);
  VECTOR<double> solution_evals(vol);

  Lop.Clear();
  solution_evecs.Clear();
  epsilon0.Clear();
  epsilon1.Clear();
  solution.Clear();
  depsilon.Clear();

  int n_iters = 1;
  int n_tries = 1;
  double tolerance = 10^-8;
  double relax_alpha = 0.5;
  double random_range = 10.0;

  t_complex solution_action = 0;
  t_complex solution_action0 = 0;
  t_complex solution_action1 = 0;

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
      epsilon0 = epsilon1;

      double de = depsilon.Norm();

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
        action1 += lambda * log(evals[n]) / 2.0;

      action = action0 + action1;

      if (!is_solution_found || (real(solution_action) > real(action)))
      {
        solution = epsilon0;
        solution_evecs = Lop;
        solution_evals = evals;

        solution_action = action;
        solution_action0 = action0;
        solution_action1 = action1;

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

    Formats::DumpBinary(f_loperator, evals);
    Formats::DumpBinary(f_loperator, Lop);

    SAFE_FWRITE(&solution_action, sizeof(t_complex), 1, f_solution);
    SAFE_FWRITE(&solution_action0, sizeof(t_complex), 1, f_solution);
    SAFE_FWRITE(&solution_action1, sizeof(t_complex), 1, f_solution);
    Formats::DumpBinary(f_solution, solution);

    fclose(f_loperator);
    fclose(f_solution);
  }

  int64_t end = Utils::GetTimeMs64();

  cout << "Program realized in " << (1.0 * end - 1.0 * begin) / 1000. << " s" << endl;

  common_AppFin();

  return 0;
}

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
    Lop(x, x) += m2 - 2.0 * I * epsilon(x) - 2.0 * ndim * Z + 4.0 * ndim * ndim * M2;

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

  tPhysicalParams params;
  params.lambda = pLambda;
  params.M2 = pM2;
  params.m2 = pm2;
  params.Z = pZ;

  ScalarField epsilon(lat);

  Lop.Clear();
  epsilon.Clear();

  main_Loperator(Lop, lat, epsilon, params);

  const string vector_current_name = "Loperator.txt";
 
  FILE *f_loperator = NULL;

  if (ModuleMPI::IsMasterNode())
  {
    const string f_bin_attr = "wb";
    const string f_txt_attr = "w";

    f_loperator = pDataDir.OpenFile(vector_current_name, f_txt_attr);

    Formats::PrintMatrix(f_loperator, Lop);
  }

  if (ModuleMPI::IsMasterNode())
  {
    fclose(f_loperator);
  }

  int64_t end = Utils::GetTimeMs64();

  cout << "Program realized in " << (1.0 * end - 1.0 * begin) / 1000. << " s" << endl;

  common_AppFin();

  return 0;
}

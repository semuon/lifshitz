#include "common.h"
#include "option_parser.h"
#include "utils.h"
#include "parameters.h"
#include "scalar_field_n.h"
#include "lattice.h"
#include "scalar_model.h"
#include "formats.h"
#include "linalg.h"
#include <time.h>
#include <mpi_module.h>

// Don't use entire namespace std
using std::cout;
using std::endl;
using std::string;

double main_VectorMean(const VECTOR<double> &vec)
{
  double mean = 0;
  int vec_size = vec.size();

  for(int i = 0; i < vec_size; i++)
    mean += vec[i] / vec_size;

  return mean;
}

double main_VectorSigma(const VECTOR<double> &vec)
{
  double sigma = 0;
  int vec_size = vec.size();

  if (vec_size > 1)
  {
    double mean = main_VectorMean(vec);

    for(int i = 0; i < vec_size; i++)
      sigma += (vec[i] - mean) * (vec[i] - mean) / (vec_size - 1);
  }

  return sqrt(sigma);
}

void main_LoadConfAt(FILE *f, off_t file_size, uint conf_idx, RealScalarFieldN &field)
{
  off_t conf_size = sizeof(double) * field.Count();

  off_t num_confs = file_size / conf_size;
  off_t remaining = file_size % conf_size;

  if (remaining != 0)
  {
    THROW_EXCEPTION(FileFailure, "Configuration file contains unknown extra data");
  }

  if (num_confs == 0)
  {
    THROW_EXCEPTION(FileFailure, "Configuration file does not contain any data");
  }

  off_t conf_begin = conf_size * conf_idx;

  SAFE_FSEEK(f, conf_begin, SEEK_SET);
  SAFE_FREAD(field.DataPtr(), sizeof(double), field.Count(), f);
}

int main_GetNumConfsFile(off_t file_size, uint filed_size)
{
  off_t conf_size = sizeof(double) * filed_size;

  off_t num_confs = file_size / conf_size;
  off_t remaining = file_size % conf_size;

  if (remaining != 0)
  {
    THROW_EXCEPTION(FileFailure, "Configuration file contains unknown extra data");
  }

  if (num_confs == 0)
  {
    THROW_EXCEPTION(FileFailure, "Configuration file does not contain any data");
  }

  return num_confs;
}

int main(int argc, char **argv)
{
  const string f_bin_write_attr = "wb";
  const string f_bin_append_attr = "ab";
  const string f_bin_read_attr = "rb";
  const string f_txt_write_attr = "w";
  const string f_txt_read_attr = "r";

  common_AppInit(argc, argv, "Lifshitz regime: Calculations");
  int64_t begin = Utils::GetTimeMs64();

  const uint nd = 1;
  const uint nc = 1;

  VECTOR<uint> latdims(pL);
  uint ndim = pDim;

  Lattice lat(latdims, ndim, nc, nd);
  uint vol = lat.Volume();

  uint corr_ls = latdims[0];
  for(uint mu = 0; mu < ndim; mu++)
  {
    if (latdims[mu] < corr_ls) 
    {
      corr_ls = latdims[mu];
    }
  }

  int n = pN;

  tScalarModelParams params;
  params.lambdaN = pLambdaN;
  params.invM2 = pInvM2;
  params.m2 = pm2;
  params.Z = pZ;
  params.N = pN;
  params.kappa = pKappa;

  if (pIsLatticeParamsSet)
  {
    tLatticeScalarModelParams lattice_params;
    lattice_params.k1 = pLatK1;
    lattice_params.k2 = pLatK2;
    lattice_params.lambda = pLatLambda;
    lattice_params.kappa = pLatKappa;

    ScalarModel::ConvertCouplings(lattice_params, ndim, params);

    pStdLogs.Write("\nEffective couplings are:\n");
    pStdLogs.Write("  kappa:                                      % -2.15le\n", params.kappa);
    pStdLogs.Write("  m^2:                                        % -2.15le\n", params.m2);
    pStdLogs.Write("  lambda*N:                                   % -2.15le\n", params.lambdaN);
    pStdLogs.Write("  m^2:                                        % -2.15le\n", params.m2);
    pStdLogs.Write("  1/M^2:                                      % -2.15le\n", params.invM2);
    pStdLogs.Write("  Z:                                          % -2.15le\n\n", params.Z);
  }

  int num_skip_first = pNumSkipFirst;
  int num_skip_last = pNumSkipLast;
  int conf_step = pConfStep;

  bool corr_vol_avg = pIsVolAvgCorr;

  const std::string fname_confs = pFnameConfs;

  const std::string fname_corr = "correlator.bin";

  FILE *f_confs = pDataDir.OpenFile(fname_confs, f_bin_read_attr);

  VECTOR<FILE *> f_corrs(ndim);
  for(uint mu = 0; mu < ndim; mu++)
    f_corrs[mu] = NULL;

  if (pIsComputeCorr)
  {
    for(uint mu = 0; mu < ndim; mu++)
    {
      std::string fname_corr_mu = TO_STRING(mu) + "." + fname_corr;

      f_corrs[mu] = pDataDir.OpenFile(fname_corr_mu, f_bin_write_attr);
    }
  }

  RealScalarFieldN phi_field(lat, n);

  VECTOR<double> corr(corr_ls);

  off_t conf_file_size = Utils::GetFileSize(f_confs);
  int num_file_confs = main_GetNumConfsFile(conf_file_size, phi_field.Count());

  int num_confs = (num_file_confs - num_skip_first - num_skip_last) / conf_step;

  if (num_confs <= 0)
  {
    pStdLogs.Write("\nNOTHING TO DO\n\n");
  }
  else
  {
    pStdLogs.Write("\nCALC BEGIN\n\n");

    pStdLogs.Write("\nNumber of confs in the file: %d\n", num_file_confs);
    pStdLogs.Write("Number of confs to process: %d\n", num_confs);
    pStdLogs.Write("Correlator max. distance: %d\n\n", corr_ls);

    for(int conf_idx = num_skip_first; conf_idx < num_file_confs - num_skip_last; conf_idx += conf_step)
    {
      main_LoadConfAt(f_confs, conf_file_size, conf_idx, phi_field);

      pStdLogs.Write("\nConf. id = %d\n", conf_idx + 1);
      pStdLogs.Write("ACTION = %2.15le\n", ScalarModel::Action(params, phi_field)/(double)vol);

      if (pIsComputeCorr)
      {
        for(uint mu = 0; mu < ndim; mu++)
        {
          ScalarModel::CorrelationFunction(phi_field, corr_vol_avg, mu, corr);

          Formats::DumpBinary(f_corrs[mu], corr);
        }
      }
    }

    pGlobalProfiler.PrintStatistics();

    pStdLogs.Write("\nCALC COMPLETED\n\n");
  }

  if (pIsComputeCorr)
  {
    for(uint mu = 0; mu < ndim; mu++)
      fclose(f_corrs[mu]);
  }

  int64_t end = Utils::GetTimeMs64();

  cout << endl;
  cout << "Program realized in " << (1.0 * end - 1.0 * begin) / 1000. << " s" << endl;

  common_AppFin();

  return 0;
}

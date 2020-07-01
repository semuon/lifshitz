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

void main_LoadLastConf(FILE *f, RealScalarFieldN &field)
{
  off_t file_size = Utils::GetFileSize(f);
  off_t conf_size = sizeof(double) * field.Count();

  off_t num_confs = file_size / conf_size;
  off_t remaining = file_size % conf_size;

  if (remaining != 0)
  {
    THROW_EXCEPTION(FileFailure, "Configuration file contains unknown extra data");
    //TERMINATE("Please check configuration file:\nsize = %lld\nconf. size = %lld\nnum. of confs = %lld\nremaining data = %lld\n",
    //  file_size, conf_size, num_confs, remaining);
  }

  if (num_confs == 0)
  {
    THROW_EXCEPTION(FileFailure, "Configuration file does not contain any data");
    //TERMINATE("Configuration file does not contain any data\n");
  }

  off_t conf_begin = conf_size * (num_confs - 1);

  SAFE_FSEEK(f, conf_begin, SEEK_SET);
  SAFE_FREAD(field.DataPtr(), sizeof(double), field.Count(), f);
}

void main_HMCUpdatePhi(const tScalarModelParams &params, const double dt,
                       RealScalarFieldN &phi, const RealScalarFieldN &pi)
{
  for(uint i = 0; i < phi.Count(); i++)
    phi[i] += pi[i] * dt;
}

void main_HMCUpdatePi(const tScalarModelParams &params, const double dt,
                      const RealScalarFieldN &phi, RealScalarFieldN &pi, RealScalarFieldN &force)
{
  ScalarModel::HMCforce(params, phi, force);

  for(uint i = 0; i < phi.Count(); i++)
    pi[i] -= force[i] * dt;
}

void main_HMCOmelyan(const tScalarModelParams &params, const double dt, const int num_steps,
                     RealScalarFieldN &phi, RealScalarFieldN &pi)
{
  const double xi = 0.1931833;

  pGlobalProfiler.StartTimer("HMC Omelyan");

  double n = params.N;

  ASSERT(n == phi.N());
  ASSERT(n == pi.N());
  ASSERT(phi.GetLattice() == pi.GetLattice());

  const Lattice &lat = phi.GetLattice();

  RealScalarFieldN force(lat, n);

  // First step
  main_HMCUpdatePhi(params, dt * xi, phi, pi);
  main_HMCUpdatePi(params, dt * 0.5, phi, pi, force);
  main_HMCUpdatePhi(params, dt * (1 - 2.0 * xi), phi, pi);
  main_HMCUpdatePi(params, dt * 0.5, phi, pi, force);

  // Intermediate steps
  for(int step_idx = 0; step_idx < num_steps - 1; step_idx++)
  {
    main_HMCUpdatePhi(params, dt * 2.0 * xi, phi, pi);
    main_HMCUpdatePi(params, dt * 0.5, phi, pi, force);
    main_HMCUpdatePhi(params, dt * (1 - 2.0 * xi), phi, pi);
    main_HMCUpdatePi(params, dt * 0.5, phi, pi, force);
  }

  // Last step
  main_HMCUpdatePhi(params, dt * xi, phi, pi);

  pGlobalProfiler.StopTimer("HMC Omelyan");
}

void main_HMCLeapfrog(const tScalarModelParams &params, const double dt, const int num_steps,
                      RealScalarFieldN &phi, RealScalarFieldN &pi)
{
  pGlobalProfiler.StartTimer("HMC Leapfrog");

  double n = params.N;

  ASSERT(n == phi.N());
  ASSERT(n == pi.N());
  ASSERT(phi.GetLattice() == pi.GetLattice());

  const Lattice &lat = phi.GetLattice();

  RealScalarFieldN force(lat, n);

  // First step of leapfrog
  main_HMCUpdatePhi(params, dt / 2.0, phi, pi);
  main_HMCUpdatePi(params, dt, phi, pi, force);

  // Intermediate steps of leapfrog
  for(int step_idx = 0; step_idx < num_steps - 1; step_idx++)
  {
    main_HMCUpdatePhi(params, dt, phi, pi);
    main_HMCUpdatePi(params, dt, phi, pi, force);
  }

  // Last step of leapfrog
  main_HMCUpdatePhi(params, dt / 2.0, phi, pi);

  pGlobalProfiler.StopTimer("HMC Leapfrog");
}

int main(int argc, char **argv)
{
  const string f_bin_write_attr = "wb";
  const string f_bin_append_attr = "ab";
  const string f_bin_read_attr = "rb";
  const string f_txt_write_attr = "w";
  const string f_txt_read_attr = "r";

  common_AppInit(argc, argv, "Lifshitz regime: HMC");
  int64_t begin = Utils::GetTimeMs64();

  const uint nd = 1;
  const uint nc = 1;

  VECTOR<uint> latdims(pL);
  uint ndim = pDim;

  Lattice lat(latdims, ndim, nc, nd);
  uint vol = lat.Volume();

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

  double hmc_dt = pHmcDt;
  int hmc_num_steps = pHmcNumSteps;
  int hmc_num_conf = pHmcNumConf;
  int hmc_num_conf_step = pHmcNumConfStep;

  int hmc_num_accepted = 0;
  int hmc_num_saved = 0;
  double hmc_accept_rate = 0;

  VECTOR<double> action_history;
  VECTOR<double> dh_history;
  VECTOR<double> exp_dh_history;
  VECTOR<double> magnetization_abs;
  VECTOR<double> magnetization_pwr_2;

  VECTOR<double> magnetization_i(n);

  const double auto_tune_dh = 0.1;
  const double auto_tune_K = pAutoTuneK;
  VECTOR<double> dh_all_history;
  VECTOR<double> hmc_dt_history;
  VECTOR<int> acceptance_history;

  action_history.reserve(hmc_num_conf);
  dh_history.reserve(hmc_num_conf);
  exp_dh_history.reserve(hmc_num_conf);
  magnetization_abs.reserve(hmc_num_conf);
  magnetization_pwr_2.reserve(hmc_num_conf);

  const std::string fname_load_conf = pFnameStartConf;
  tStartConfigurationType start_type = pStartType;

  tIntegratorType integrator_type = pIntegratorType;

  const std::string fname_hmc_stat = "hmc_stat.txt";
  const std::string fname_simple_observables = "simple_observables.txt";
  const std::string fname_magnetization = "magnetization.bin";
  const std::string fname_confs = "confs.bin";
  const std::string fname_auto_tune = "auto_tune.txt";

  FILE *f_hmc_stat = pDataDir.OpenFile(fname_hmc_stat, f_txt_write_attr);
  FILE *f_confs = pDataDir.OpenFile(fname_confs, f_bin_write_attr);
  FILE *f_simple_observables = pDataDir.OpenFile(fname_simple_observables, f_txt_write_attr);
  FILE *f_magnetization = pDataDir.OpenFile(fname_magnetization, f_bin_write_attr);

  FILE *f_auto_tune = NULL;
  if (pAutoTune)
  {
    f_auto_tune = pDataDir.OpenFile(fname_auto_tune, f_txt_write_attr);
  }

  RealScalarFieldN pi_field(lat, n);
  RealScalarFieldN phi_field_0(lat, n);
  RealScalarFieldN phi_field_1(lat, n);
  RealScalarFieldN phi_force(lat, n);

  switch(start_type)
  {
    case START_CONFIGURATION_ZERO:
    {
      for(uint i = 0; i < n * vol; i++)
        phi_field_0[i] = 0;
    }
    break;
    case START_CONFIGURATION_LOAD:
    {
      FILE *f_load_conf = pDataDir.OpenFile(fname_load_conf, f_bin_read_attr);
      main_LoadLastConf(f_load_conf, phi_field_0);
      fclose(f_load_conf);
    }
    break;
    {
    case START_CONFIGURATION_RANDOM:
    default:
      for(uint i = 0; i < n * vol; i++)
        phi_field_0[i] = rand_double(-1, 1);
    }
    break;
  }

  pStdLogs.Write("\n\nHMC BEGIN\n\n");

  double start_action = ScalarModel::Action(params, phi_field_0) / (double)vol;
  pStdLogs.Write("Action on start configuration: %2.15le\n", start_action);
  pStdLogs.Write("Norm of start configuration: %2.15le\n\n", phi_field_0.Norm());

  int conf_idx;
  for(conf_idx = 0; conf_idx < hmc_num_conf; conf_idx++)
  {
    double h0 = 0;
    double h1 = 0;
    double conf_action = 0;

    pGlobalProfiler.StartSection("HMC Trajectory");

    // Initialize hmc and compute initial hamiltonian h0
    phi_field_1 = phi_field_0;

    for(uint i = 0; i < n * vol; i++)
    {
      pi_field[i] = rand_gauss_double(0, 1);

      h0 += pi_field[i] * pi_field[i] / 2.0;
    }

    h0 += ScalarModel::Action(params, phi_field_1);

    switch(integrator_type)
    {
      case INTEGRATOR_OMELYAN:
        main_HMCOmelyan(params, hmc_dt, hmc_num_steps, phi_field_1, pi_field);
      break;
      case INTEGRATOR_LEAPFROG:
      default:
        main_HMCLeapfrog(params, hmc_dt, hmc_num_steps, phi_field_1, pi_field);
      break;
    }

    for(uint i = 0; i < n * vol; i++)
      h1 += pi_field[i] * pi_field[i] / 2.0;

    conf_action = ScalarModel::Action(params, phi_field_1);
    h1 += conf_action;

    // Accept/reject step
    double hmc_dh = h1 - h0;
    bool accepted = false;

    // Check if energy was finite
    if (!isfinite(hmc_dh))
    {
      pStdLogs.Write("\n\nValue of dh is not finite.\nExiting...\n\n");
      pStdLogs.WriteError(__FILE__, __LINE__, "Value of dh is not finite");
      break;
    }

    if (hmc_dh <= 0)
    {
      accepted = true;
    }
    else
    {
      double alpha = exp(-hmc_dh);
      double dice = rand_double(0, 1);

      accepted = (dice < alpha);
    }

    dh_all_history.push_back(hmc_dh);
    hmc_dt_history.push_back(hmc_dt);
    acceptance_history.push_back( (accepted) ? 1 : 0 );

    // Auto-tune
    if (pAutoTune)
    {
      SAFE_FPRINTF(f_auto_tune, "%2.15le\t%2.15le\t%d\n", hmc_dh, hmc_dt, (accepted) ? 1 : 0);
      fflush(f_auto_tune);

      double R = auto_tune_dh - hmc_dh;
      double F;

      switch(integrator_type)
      {
        case INTEGRATOR_OMELYAN:
          F = auto_tune_K * R / (2.0 * hmc_dt);
        break;
        case INTEGRATOR_LEAPFROG:
          F = auto_tune_K * R;
        default:
          F = 0;
        break;
      }

      if ((!accepted || fabs(R) <= auto_tune_dh) && isfinite(F))
      {
        double new_dt = hmc_dt + F;

        if (new_dt < 0.5 * hmc_dt)
          hmc_dt = hmc_dt / 2;
        else if (new_dt > 2 * hmc_dt)
          hmc_dt = hmc_dt * 2;
        else
          hmc_dt = new_dt;
      }
    }

    pStdLogs.Write("conf = %d\n", conf_idx + 1);

    if(accepted)
    {
      pStdLogs.Write("status: ACCEPTED\n");

      phi_field_0 = phi_field_1;

      // Measurements
      if ((hmc_num_accepted + 1) % hmc_num_conf_step == 0)
      {
        pStdLogs.Write("MEASUREMENTS\n");

        // Condensate
        double m_pwr_2 = 0;
        double m_abs = 0;

        for(int i = 0; i < n; i++)
        {
          double cond = 0;

          for(uint x = 0; x < vol; x++)
            cond += phi_field_0(x, i);

          magnetization_i[i] = cond;
          m_pwr_2 += cond * cond;
          m_abs += fabs(cond);
        }

        Formats::DumpBinary(f_magnetization, magnetization_i);
        magnetization_abs.push_back(m_abs);
        magnetization_pwr_2.push_back(m_pwr_2);
        action_history.push_back(conf_action / vol);
        dh_history.push_back(hmc_dh);
        exp_dh_history.push_back(exp(-hmc_dh));

        // Save configuration
        Formats::DumpBinary(f_confs, phi_field_0);

        hmc_num_saved++;
      }

      hmc_num_accepted++;

      pStdLogs.Write("action = %2.15le\n", conf_action / vol);
    }
    else
    {
      pStdLogs.Write("status: REJECTED\n");
    }

    pStdLogs.Write("dh = %2.15le\n\n", hmc_dh);

    pGlobalProfiler.EndSection("HMC Trajectory");
  }

  pGlobalProfiler.PrintStatistics();

  pStdLogs.Write("\n\nHMC COMPLETED\n\n");

  hmc_accept_rate = (conf_idx > 0) ? (double)hmc_num_accepted / (double)conf_idx : 0;

  pStdLogs.Write("Num. accepted: %d / %d\n", hmc_num_accepted, conf_idx);
  pStdLogs.Write("Acceptance rate: %2.2lf\n", hmc_accept_rate);
  pStdLogs.Write("<exp(dh)> = %2.4lf +/- %2.4lf\n", main_VectorMean(exp_dh_history), main_VectorSigma(exp_dh_history));
  pStdLogs.Write("<dh> = %2.4lf +/- %2.4lf\n\n", main_VectorMean(dh_history), main_VectorSigma(dh_history));

  for(int i = 0; i < hmc_num_saved; i++)
  {
    SAFE_FPRINTF(f_hmc_stat, "%2.15le\t%2.15le\t%2.15le\n", action_history[i], exp_dh_history[i], dh_history[i]);
    SAFE_FPRINTF(f_simple_observables, "%2.15le\t%2.15le\t%2.15le\n", action_history[i], magnetization_abs[i], magnetization_pwr_2[i]);
  }

  // for(int i = 0; i < hmc_dt_history.size(); i++)
  // {
  //   SAFE_FPRINTF(f_auto_tune, "%2.15le\t%2.15le\t%2.15le\n", dh_all_history[i], hmc_dt_history[i], acceptance_history[i]);
  // }

  fclose(f_hmc_stat);
  fclose(f_confs);
  fclose(f_simple_observables);
  fclose(f_magnetization);

  if (pAutoTune)
  {
    fclose(f_auto_tune);
  }

  int64_t end = Utils::GetTimeMs64();

  cout << endl;
  cout << "Program realized in " << (1.0 * end - 1.0 * begin) / 1000. << " s" << endl;

  common_AppFin();

  return 0;
}

#include "common.h"
#include "option_parser.h"
#include "utils.h"
#include "parameters.h"
#include "scalar_field_n.h"
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
  double lambdaN;
  double N;
  double kappa;
} tPhysicalParams;

double main_Action(const tPhysicalParams &params, const RealScalarFieldN &phi)
{
  pGlobalProfiler.StartTimer("Action");

  double m2 = params.m2;
  double invM2 = params.invM2;
  double Z = params.Z;
  double lambdaN = params.lambdaN;
  double n = params.N;
  double kappa = params.kappa;

  ASSERT(n == phi.N());

  const Lattice &lat = phi.GetLattice();

  int ndim = lat.Dim();
  uint vol = lat.Volume();

  double res = 0;

  // Kinetic terms
  for(uint i = 0; i < n; i++)
  {
    for(uint x = 0; x < vol; x++)
    {
      res += (m2 - 2.0 * ndim * Z + 6.0 * ndim * invM2) * phi(x, i) * phi(x, i) / 2.0;

      for(int mu = 0; mu < ndim; mu++)
      {
        uint xmu_fwd = lat.SiteIndexForward(x, mu);
        uint xmu_bwd = lat.SiteIndexBackward(x, mu);
        uint x2mu_fwd = lat.SiteIndexForward(xmu_fwd, mu);
        uint x2mu_bwd = lat.SiteIndexBackward(xmu_bwd, mu);

        res += phi(x, i) * invM2 * (phi(x2mu_fwd, i) + phi(x2mu_bwd, i)) / 2.0;
        res += phi(x, i) * (Z - 4.0 * invM2) * (phi(xmu_fwd, i) + phi(xmu_bwd, i)) / 2.0;
      }
    }
  }

  // Interaction terms
  for(uint x = 0; x < vol; x++)
  {
    double phi2 = 0;

    for(uint i = 0; i < n; i++)
      phi2 += phi(x, i) * phi(x, i);

    res += lambdaN * phi2 * phi2 / (4.0 * n);
    res += kappa * phi2 * phi2 * phi2 / 6.0;
  }

  pGlobalProfiler.StopTimer("Action");

  return res;
}

void main_HMCforce(const tPhysicalParams &params, const RealScalarFieldN &phi, RealScalarFieldN &force)
{
  pGlobalProfiler.StartTimer("HMC Force");

  double m2 = params.m2;
  double invM2 = params.invM2;
  double Z = params.Z;
  double lambdaN = params.lambdaN;
  double n = params.N;
  double kappa = params.kappa;

  ASSERT(n == phi.N());
  ASSERT(n == force.N());

  const Lattice &lat = phi.GetLattice();

  int ndim = lat.Dim();
  uint vol = lat.Volume();

  // Kinetic terms
  for(uint i = 0; i < n; i++)
  {
    for(uint x = 0; x < vol; x++)
    {
      force(x, i) = (m2 - 2.0 * ndim * Z + 6.0 * ndim * invM2) * phi(x, i);

      for(int mu = 0; mu < ndim; mu++)
      {
        uint xmu_fwd = lat.SiteIndexForward(x, mu);
        uint xmu_bwd = lat.SiteIndexBackward(x, mu);
        uint x2mu_fwd = lat.SiteIndexForward(xmu_fwd, mu);
        uint x2mu_bwd = lat.SiteIndexBackward(xmu_bwd, mu);

        force(x, i) += invM2 * (phi(x2mu_fwd, i) + phi(x2mu_bwd, i));
        force(x, i) += (Z - 4.0 * invM2) * (phi(xmu_fwd, i) + phi(xmu_bwd, i));
      }
    }
  }

  // Interaction terms
  for(uint x = 0; x < vol; x++)
  {
    double phi2 = 0;

    for(uint i = 0; i < n; i++)
      phi2 += phi(x, i) * phi(x, i);

    for(uint i = 0; i < n; i++)
    {
      force(x, i) += lambdaN * phi2 * phi(x, i) / n;
      force(x, i) += kappa * phi2 * phi2 * phi(x, i);
    }
  }

  pGlobalProfiler.StopTimer("HMC Force");
}

int main(int argc, char **argv)
{
  const string f_bin_attr = "wb";
  const string f_txt_attr = "w";

  common_AppInit(argc, argv, "Lifshitz regime: HMC");
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
  double lambdaN = pLambdaN;
  double n = 1;
  double kappa = 0;

  double hmc_dt = 0.05;
  int hmc_num_steps = 20;
  int hmc_num_traj = 20;

  double hmc_accept_rate = 0;
  double hmc_avg_dh = 0;
  double hmc_avg_exp_dh = 0;

  tStartConfigurationType start_type = START_CONFIGURATION_RANDOM;

  string fname_load_conf = "init.conf";

  tPhysicalParams params;
  params.lambdaN = lambdaN;
  params.invM2 = invM2;
  params.m2 = m2;
  params.Z = Z;
  params.N = n;
  params.kappa = kappa;

  RealScalarFieldN phi(lat, n);
  RealScalarFieldN force(lat, n);

  for(uint idx = 0; idx < n * vol; idx++)
    phi[idx] = (double) idx;

  double action = main_Action(params, phi);
  main_HMCforce(params, phi, force);

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
      FILE *f_load_conf = pDataDir.OpenFile(fname_load_conf, f_bin_attr);
      SAFE_FREAD(phi_field_0.DataPtr(), sizeof(double), vol * n, f_load_conf);
      fclose(f_load_conf);
    }
    break;
    {
    case START_CONFIGURATION_RANDOM:
    default:
      for(uint i = 0; i < n * vol; i++)
        phi_field_0[i] = rand_gauss_double(0, 1);
    }
    break;
  }

  pStdLogs.Write("\n\nHMC BEGIN\n\n");

  for(int traj_idx = 0; traj_idx < hmc_num_traj; traj_idx++)
  {
    double h0 = 0;
    double h1 = 0;
    double traj_action = 0;

    pGlobalProfiler.StartSection("HMC Trajectory");

    // Initialize hmc and compute initial hamiltonian h0
    for(uint i = 0; i < n * vol; i++)
    {
      pi_field[i] = rand_gauss_double(0, 1);

      h0 += pi_field[i] * pi_field[i] / 2.0;
    }

    phi_field_1 = phi_field_0;

    h0 += main_Action(params, phi_field_1);

    // First step of leapfrog
    main_HMCforce(params, phi_field_1, phi_force);

    for(uint i = 0; i < n * vol; i++)
    {
      pi_field[i] -= phi_force[i] * hmc_dt / 2.0;
      phi_field_1[i] += pi_field[i] * hmc_dt;
    }

    // Intermediate steps of leapfrog
    for(int step_idx = 0; step_idx < hmc_num_steps - 1; step_idx++)
    {
      main_HMCforce(params, phi_field_1, phi_force);

      for(uint i = 0; i < n * vol; i++)
      {
        pi_field[i] -= phi_force[i] * hmc_dt;
        phi_field_1[i] += pi_field[i] * hmc_dt;
      }
    }

    // Last step of leapfrog and computation of final hamiltonian h1
    main_HMCforce(params, phi_field_1, phi_force);

    for(uint i = 0; i < n * vol; i++)
    {
      pi_field[i] -= phi_force[i] * hmc_dt / 2.0;

      h1 += pi_field[i] * pi_field[i] / 2.0;
    }

    traj_action = main_Action(params, phi_field_1);
    h1 += traj_action;

    // Accept/reject step
    double hmc_dh = h1 - h0;
    bool accepted = false;

    if (hmc_dh <= 0)
    {
      accepted = true;
    }
    else
    {
      double alpha = exp(-hmc_dh);
      double dice = rand_double(0, 1);

      accepted = (alpha > dice);
    }

    pStdLogs.Write("step = %d\n", traj_idx + 1);

    if(accepted)
    {
      phi_field_0 = phi_field_1;
      hmc_accept_rate += 1.0 / hmc_num_traj;

      hmc_avg_dh += hmc_dh / hmc_num_traj;
      hmc_avg_exp_dh += exp(-hmc_dh) / hmc_num_traj;

      pStdLogs.Write("status: ACCEPTED\n");
      pStdLogs.Write("action = %2.15le\n", traj_action);
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
  pStdLogs.Write("<exp(dh)> = %2.15le\n", hmc_avg_exp_dh);
  pStdLogs.Write("<dh> = %2.15le\n\n", hmc_avg_dh);

  pStdLogs.Write("Action: %2.15le", action);
  FILE *f = pDataDir.OpenFile("test.bin", "wb");
  Formats::DumpBinary(f, force);
  fclose(f);

  int64_t end = Utils::GetTimeMs64();

  cout << endl;
  cout << "Program realized in " << (1.0 * end - 1.0 * begin) / 1000. << " s" << endl;

  common_AppFin();

  return 0;
}

#include "scalar_model.h"
#include "profiler.h"

void ScalarModel::ConvertCouplings(const tLatticeScalarModelParams &lattice_params, const int ndim, tScalarModelParams &phys_params)
{
  phys_params.invM2 = 2.0 * lattice_params.k2;
  phys_params.Z = 2.0 * lattice_params.k1 - 8.0 * lattice_params.k2;
  phys_params.m2 = 2.0 + 4.0 * ndim * lattice_params.k2 - 4.0 * ndim * lattice_params.k1 - 4.0 * lattice_params.lambda;
  phys_params.lambdaN = 4.0 * phys_params.N * lattice_params.lambda;
  phys_params.kappa = 6.0 * lattice_params.kappa;
}

void ScalarModel::CorrelationFunction(const RealScalarFieldN &phi, const bool vol_avg, VECTOR<double> &corr)
{
  const Lattice &lat = phi.GetLattice();
  const uint ls = corr.size();

  uint ndim = lat.Dim();
  uint vol = lat.Volume();
  uint n = phi.N();

  for(uint mu = 0; mu < ndim; mu++)
    ASSERT(lat.LatticeSize(mu) >= ls);

  for(uint i = 0; i < ls; i++)
    corr[i] = 0;

  uint xmax = (vol_avg) ? vol : 1;

  for(uint x = 0; x < xmax; x++)
  for(uint mu = 0; mu < ndim; mu++)
  {
    uint x1 = x;

    for(uint s = 0; s < ls; s++)
    {
      for(uint i = 0; i < n; i++)
      {
        corr[s] += phi(x, i) * phi(x1, i) / (n * xmax * ndim);
      }

      x1 = lat.SiteIndexForward(x1, mu);
    }
  }
}

double ScalarModel::Action(const tScalarModelParams &params, const RealScalarFieldN &phi)
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
      res += (m2 + invM2 * 4.0 * (double)(ndim * ndim) + Z * 5.0 * (double)ndim / 2.0) * phi(x, i) * phi(x, i);

      for(int mu = 0; mu < ndim; mu++)
      {
        uint xmu_fwd = lat.SiteIndexForward(x, mu);
        uint xmu_bwd = lat.SiteIndexBackward(x, mu);
        uint x2mu_fwd = lat.SiteIndexForward(xmu_fwd, mu);
        uint x2mu_bwd = lat.SiteIndexBackward(xmu_bwd, mu);

        for(int nu = 0; nu < ndim; nu++)
        {
          uint xmu_fwd_nu_fwd = lat.SiteIndexForward(xmu_fwd, nu);
          uint xmu_fwd_nu_bwd = lat.SiteIndexBackward(xmu_fwd, nu);
          uint xmu_bwd_nu_fwd = lat.SiteIndexForward(xmu_bwd, nu);
          uint xmu_bwd_nu_bwd = lat.SiteIndexBackward(xmu_bwd, nu);

          res += phi(x, i) * invM2 * (phi(xmu_fwd_nu_fwd, i) + phi(xmu_fwd_nu_bwd, i) + phi(xmu_bwd_nu_fwd, i) + phi(xmu_bwd_nu_bwd, i));
        }

        res -= phi(x, i) * 4.0 * (double)ndim * invM2 * (phi(xmu_fwd, i) + phi(xmu_bwd, i));

        res -= phi(x, i) * Z * (phi(xmu_fwd, i) + phi(xmu_bwd, i));
        res += phi(x, i) * (Z / 12.0) * (phi(x2mu_fwd, i) + phi(x2mu_bwd, i) - 4.0 * phi(xmu_fwd, i) - 4.0 * phi(xmu_bwd, i));
      }
    }
  }

  res *= 1.0 /2.0;

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

void ScalarModel::HMCforce(const tScalarModelParams &params, const RealScalarFieldN &phi, RealScalarFieldN &force)
{
  pGlobalProfiler.StartTimer("HMC Force");

  double m2 = params.m2;
  double invM2 = (params.invM2 + params.Z / 12.0);
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
      force(x, i) = (m2 + 2.0 * ndim * Z + 6.0 * ndim * invM2) * phi(x, i);

      for(int mu = 0; mu < ndim; mu++)
      {
        uint xmu_fwd = lat.SiteIndexForward(x, mu);
        uint xmu_bwd = lat.SiteIndexBackward(x, mu);
        uint x2mu_fwd = lat.SiteIndexForward(xmu_fwd, mu);
        uint x2mu_bwd = lat.SiteIndexBackward(xmu_bwd, mu);

        force(x, i) += invM2 * (phi(x2mu_fwd, i) + phi(x2mu_bwd, i));
        force(x, i) += (-Z - 4.0 * invM2) * (phi(xmu_fwd, i) + phi(xmu_bwd, i));
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

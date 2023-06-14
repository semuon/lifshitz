#include "scalar_model.h"
#include "profiler.h"

void ScalarModel::CreateLatticeOperators(tScalarModelParams &params, const uint ndim, const int n_stencil_points)
{
  const uint op_dim = 1;

  auto fwd1d = FiniteDifference<int64_t>::MakeOneSidedDiff(op_dim, n_stencil_points);
  auto bwd1d = FiniteDifference<int64_t>::MakeOneSidedDiff(op_dim, -n_stencil_points);

  params.laplace_ptr = FiniteDifference<int64_t>::MakeLaplacian(ndim, *fwd1d, *bwd1d);
  params.laplace_sqr_ptr = FiniteDifference<int64_t>::ComposeOperators(*params.laplace_ptr, *params.laplace_ptr);
}

void ScalarModel::ConvertCouplings(const tLatticeScalarModelParams &lattice_params, const int ndim, tScalarModelParams &phys_params)
{
  // NOT IMPLEMENTED
  ASSERT(false);

  phys_params.invM2 = 2.0 * lattice_params.k2;
  phys_params.Z = 2.0 * lattice_params.k1 - 8.0 * lattice_params.k2;
  phys_params.m2 = 2.0 + 4.0 * ndim * lattice_params.k2 - 4.0 * ndim * lattice_params.k1 - 4.0 * lattice_params.lambda;
  phys_params.lambdaN = 4.0 * phys_params.N * lattice_params.lambda;
  phys_params.kappa = 6.0 * lattice_params.kappa;
}

void ScalarModel::CorrelationFunction(const RealScalarFieldN &phi, const bool vol_avg, const uint mu, VECTOR<double> &corr)
{
  const Lattice &lat = phi.GetLattice();

  uint vol = lat.Volume();
  uint n = phi.N();
  uint ls = lat.LatticeSize(mu);

  corr.resize(ls);

  for(uint i = 0; i < ls; i++)
    corr[i] = 0;

  uint xmax = (vol_avg) ? vol : 1;

  for(uint x = 0; x < xmax; x++)
  {
    uint x1 = x;

    for(uint s = 0; s < ls; s++)
    {
      for(uint i = 0; i < n; i++)
      {
        corr[s] += phi(x, i) * phi(x1, i) / (n * xmax);
      }

      x1 = lat.SiteIndexForward(x1, mu);
    }
  }
}

uint ScalarModel::SiteIndexByOffset(const Lattice &lat, uint x, const AuxVector<int> &offset)
{
  uint ndim = lat.Dim();
  uint res = x;

  for(uint mu = 0; mu < ndim; mu++)
  {
    int dx = offset.GetComponent(mu);

    if (dx > 0)
    {
      for(int i = 0; i < dx; i++)
        res = lat.SiteIndexForward(res, mu);
    }
    else if (dx < 0)
    {
      for(int i = 0; i < -dx; i++)
        res = lat.SiteIndexBackward(res, mu);
    }
  }

  return res;
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

  // Finite differences stencils
  const VECTOR<StencilPoint<int, int64_t>> &lap = params.laplace_ptr->GetStencil();
  const VECTOR<StencilPoint<int, int64_t>> &lap_sqr = params.laplace_sqr_ptr->GetStencil();

  double res = 0;

  for(uint i = 0; i < n; i++)
  {
    for(uint x = 0; x < vol; x++)
    {
      // Laplacian
      for(uint si = 0; si < lap.size(); si++)
      {
        uint y = SiteIndexByOffset(lat, x, lap[si].offset);
        boost::rational<int64_t> coef = lap[si].coef;

        res += -0.5 * Z * boost::rational_cast<double>(coef) * phi(x, i) * phi(y, i);
      }

      // Laplacian squared
      for(uint si = 0; si < lap_sqr.size(); si++)
      {
        uint y = SiteIndexByOffset(lat, x, lap_sqr[si].offset);
        boost::rational<int64_t> coef = lap_sqr[si].coef;

        res += 0.5 * invM2 * boost::rational_cast<double>(coef) * phi(x, i) * phi(y, i);
      }

      // Mass
      res += 0.5 * m2 * phi(x, i) * phi(x, i);
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

void ScalarModel::HMCforce(const tScalarModelParams &params, const RealScalarFieldN &phi, RealScalarFieldN &force)
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
      force(x, i) = (m2 + invM2 * 4.0 * (double)(ndim * ndim) + Z * 5.0 * (double)ndim / 2.0) * phi(x, i);

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

          force(x, i) += invM2 * (phi(xmu_fwd_nu_fwd, i) + phi(xmu_fwd_nu_bwd, i) + phi(xmu_bwd_nu_fwd, i) + phi(xmu_bwd_nu_bwd, i));
        }

        force(x, i) -= 4.0 * (double)ndim * invM2 * (phi(xmu_fwd, i) + phi(xmu_bwd, i));

        force(x, i) -= Z * (phi(xmu_fwd, i) + phi(xmu_bwd, i));
        force(x, i) += (Z / 12.0) * (phi(x2mu_fwd, i) + phi(x2mu_bwd, i) - 4.0 * phi(xmu_fwd, i) - 4.0 * phi(xmu_bwd, i));
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

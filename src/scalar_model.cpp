#include "scalar_model.h"
#include "profiler.h"

// Gives a linear index dimension which is N*(N+1)/2
uint dim_upper_triangular_index(const uint N)
{
  return N * (N + 1) / 2;
}

// Converts a pair 0 <= i <= j < N into a linear index 0 <= k < N*(N+1)/2
uint pack_upper_triangular_index(const uint N, const uint i, const uint j)
{
  return N * (N - 1) / 2 - (N - i) * (N - i -1) / 2 + j;
}

// Converts a linear index 0 <= k < N*(N+1)/2 into pair 0 <= i <= j < N 
void unpack_upper_triangular_index(const uint N, const uint k, uint &i, uint &j)
{
  i = N - 1 - (int)std::floor(std::sqrt(-8 * k + 4 * N * (N + 1) - 7)/2.0 - 0.5);
  j = k - N * (N - 1) / 2 + (N - i) * ((N - i) - 1) / 2;
}

void ScalarModel::CreateLatticeOperators(tScalarModelParams &params, const uint ndim, const int n_stencil_points)
{
  const uint op_dim = 2;
  auto diff1d = FiniteDifference<int64_t>::MakeSymmetricDiff(op_dim, (uint)n_stencil_points);;

  params.laplace_ptr = FiniteDifference<int64_t>::MakeLaplacian(ndim, *diff1d);
  params.laplace_sqr_ptr = FiniteDifference<int64_t>::ComposeOperators(*params.laplace_ptr, *params.laplace_ptr);
}

void ScalarModel::CreateHoppings(tScalarModelParams &params, const Lattice &lat)
{
  uint vol = lat.Volume();
  uint ndim = lat.Dim();

  // Finite differences stencils
  const VECTOR<StencilPoint<int, int64_t>> &lap = params.laplace_ptr->GetStencil();
  const VECTOR<StencilPoint<int, int64_t>> &lap_sqr = params.laplace_sqr_ptr->GetStencil();
  AuxVector<int> offset_bwd(ndim);

  uint lap_npts = lap.size();
  uint lap_sqr_npts = lap_sqr.size();

  params.lap_hoppings.resize(vol * lap_npts * 2);
  params.lap_sqr_hoppings.resize(vol * lap_sqr_npts * 2);

  for(uint x = 0; x < vol; x++)
  {
    // Laplacian
    for(uint si = 0; si < lap_npts; si++)
    {
      for(uint mu = 0; mu < ndim; mu++) { offset_bwd[mu] = -lap[si].offset.GetComponent(mu); }

      uint yfwd = SiteIndexByOffset(lat, x, lap[si].offset);
      uint ybwd = SiteIndexByOffset(lat, x, offset_bwd);

      params.lap_hoppings[x * lap_npts * 2 + si * 2 + 0] = yfwd;
      params.lap_hoppings[x * lap_npts * 2 + si * 2 + 1] = ybwd;
    }

    // Laplacian squared
    for(uint si = 0; si < lap_sqr_npts; si++)
    {
      for(uint mu = 0; mu < ndim; mu++) { offset_bwd[mu] = -lap_sqr[si].offset.GetComponent(mu); }

      uint yfwd = SiteIndexByOffset(lat, x, lap_sqr[si].offset);
      uint ybwd = SiteIndexByOffset(lat, x, offset_bwd);

      params.lap_sqr_hoppings[x * lap_sqr_npts * 2 + si * 2 + 0] = yfwd;
      params.lap_sqr_hoppings[x * lap_sqr_npts * 2 + si * 2 + 1] = ybwd;
    }
  }
}

void ScalarModel::ConvertCouplings(const tLatticeScalarModelParams &lattice_params, const int ndim, tScalarModelParams &phys_params)
{
  // NOT IMPLEMENTED
  // ASSERT(false);

  phys_params.invM2 = 2.0 * lattice_params.k2;
  phys_params.Z = 2.0 * lattice_params.k1 - 8.0 * lattice_params.k2;
  phys_params.m2 = 2.0 + 4.0 * ndim * lattice_params.k2 - 4.0 * ndim * lattice_params.k1 - 4.0 * lattice_params.lambda;
  phys_params.lambdaN = 4.0 * phys_params.N * lattice_params.lambda;
  phys_params.kappa = 6.0 * lattice_params.kappa;
}

void ScalarModel::ExternalField(RealScalarFieldN &h, const double h0, const double k0, const double sigma0)
{
  const Lattice &lat = h.GetLattice();

  uint vol = lat.Volume();
  uint ndim = lat.Dim();
  uint n = h.N();

  const int mu = 0;
  const uint l = lat.LatticeSize((uint)mu);

  VECTOR<int> xs(ndim);

  h.Clear();

  for(uint x = 0; x < vol; x++)
  {
    lat.SiteCoordinates(xs, x);

    for(uint k = 0; k < l; k++)
    {
      double p = 2 * M_PI * k / (double)l;
      double prefactor = h0 * exp(-(p - k0) * (p - k0) / (2 * sigma0)) / ((double)l * sigma0 * sqrt(2 * M_PI));

      for (uint i = 0; i < n; i++)
      {
        if (i == 0)
        {
          h(x, i) += prefactor * cos(p * (double)xs[mu]);
        }
        else if (i == 1)
        {
          h(x, i) += prefactor * sin(p * (double)xs[mu]);
        }
      }
    }
  }
}

void ScalarModel::CorrelationFunctionOLD(const RealScalarFieldN &phi, const bool vol_avg, const uint mu, VECTOR<double> &corr)
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

void ScalarModel::CorrelationMatrix(const RealScalarFieldN &phi, const bool vol_avg, const uint mu, VECTOR<double> &corr)
{
  const Lattice &lat = phi.GetLattice();

  uint vol = lat.Volume();
  uint n = phi.N();
  uint ls = lat.LatticeSize(mu);

  corr.resize(ls * n * n);

  for(uint i = 0; i < corr.size(); i++)
    corr[i] = 0;

  uint xmax = (vol_avg) ? vol : 1;

  for(uint x = 0; x < xmax; x++)
  {
    uint x1 = x;

    for(uint s = 0; s < ls; s++)
    {
      for(uint i = 0; i < n; i++)
      for(uint j = i; j < n; j++) // I'm lazy
      {
        corr[s * n * n + i * n + j] += phi(x, i) * phi(x1, j) / (xmax);
      }

      x1 = lat.SiteIndexForward(x1, mu);
    }
  }
}

void ScalarModel::FullTwoPointFunction(const RealScalarFieldN &phi, VECTOR<double> &corr)
{
  const Lattice &lat = phi.GetLattice();

  uint vol = lat.Volume();
  uint n = phi.N();

  corr.resize(vol * (vol + 1) * n * (n + 1) / 4);

  for(uint i = 0; i < corr.size(); i++)
    corr[i] = 0;

  uint idx = 0;
  for(uint x = 0; x < vol; x++)
  for(uint y = x; y < vol; y++)
  {
      for(uint i = 0; i < n; i++)
      for(uint j = i; j < n; j++)
      {
        corr[idx] += phi(x, i) * phi(y, j);
        idx++;
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
  const RealScalarFieldN &ext_h = *params.h_ptr;

  uint vol = lat.Volume();

  // Finite differences stencils
  const VECTOR<StencilPoint<int, int64_t>> &lap = params.laplace_ptr->GetStencil();
  const VECTOR<StencilPoint<int, int64_t>> &lap_sqr = params.laplace_sqr_ptr->GetStencil();

  uint lap_npts = lap.size();
  uint lap_sqr_npts = lap_sqr.size();

  double res = 0;

  for(uint i = 0; i < n; i++)
  {
    for(uint x = 0; x < vol; x++)
    {
      // Laplacian
      for(uint si = 0; si < lap_npts; si++)
      {
        uint hop_idx = x * lap_npts * 2 + si * 2 + 0;

        uint y = params.lap_hoppings[hop_idx];
        boost::rational<int64_t> coef = lap[si].coef;

        res += -0.5 * Z * boost::rational_cast<double>(coef) * phi(x, i) * phi(y, i);
      }

      // Laplacian squared
      for(uint si = 0; si < lap_sqr_npts; si++)
      {
        uint hop_idx = x * lap_sqr_npts * 2 + si * 2 + 0;

        uint y = params.lap_sqr_hoppings[hop_idx];
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

  // External field
  for(uint i = 0; i < n; i++)
  {
    for(uint x = 0; x < vol; x++)
    {
      res += -ext_h(x, i) * phi(x, i);
    }
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
  const RealScalarFieldN &ext_h = *params.h_ptr;

  uint vol = lat.Volume();

  // Finite differences stencils
  const VECTOR<StencilPoint<int, int64_t>> &lap = params.laplace_ptr->GetStencil();
  const VECTOR<StencilPoint<int, int64_t>> &lap_sqr = params.laplace_sqr_ptr->GetStencil();

  uint lap_npts = lap.size();
  uint lap_sqr_npts = lap_sqr.size();

  for(uint i = 0; i < n; i++)
  {
    for(uint x = 0; x < vol; x++)
    {
      // Mass
      force(x, i) = m2 * phi(x, i);

      // Laplacian
      for(uint si = 0; si < lap_npts; si++)
      {
        uint hop_idx = x * lap_npts * 2 + si * 2;

        uint yfwd = params.lap_hoppings[hop_idx + 0];
        uint ybwd = params.lap_hoppings[hop_idx + 1];

        boost::rational<int64_t> coef = lap[si].coef;

        force(x, i) += -0.5 * Z * boost::rational_cast<double>(coef) * (phi(yfwd, i) + phi(ybwd, i));
      }

      // Laplacian squared
      for(uint si = 0; si < lap_sqr_npts; si++)
      {
        uint hop_idx = x * lap_sqr_npts * 2 + si * 2;

        uint yfwd = params.lap_sqr_hoppings[hop_idx + 0];
        uint ybwd = params.lap_sqr_hoppings[hop_idx + 1];

        boost::rational<int64_t> coef = lap_sqr[si].coef;

        force(x, i) += 0.5 * invM2 * boost::rational_cast<double>(coef) * (phi(yfwd, i) + phi(ybwd, i));
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

  // External field
  for(uint i = 0; i < n; i++)
  {
    for(uint x = 0; x < vol; x++)
    {
      force(x, i) += -ext_h(x, i);
    }
  }

  pGlobalProfiler.StopTimer("HMC Force");
}

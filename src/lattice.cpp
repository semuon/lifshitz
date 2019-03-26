#include "lattice.h"

// ================================================= //
// =================  Constructors  ================ //
// ================================================= //

Lattice::Lattice(const VECTOR<uint> &L, uint ndim, uint sunrank, uint ndirac) :
  vol(1), sun(sunrank), nd(ndirac), dim(ndim), lat_sizes(dim)
{
  ASSERTVERB(L.size() == ndim, "Lattice constructor: ndim should be equal to dimension of L");

  for(uint i = 0; i < dim; i++)
  {
    lat_sizes[i] = L[i];
    vol *= L[i];
  }
}

Lattice::Lattice(const Lattice &other, uint sunrank, uint ndirac) :
  vol(other.vol), sun(sunrank), nd(ndirac), dim(other.dim), lat_sizes(other.lat_sizes)
{
}

Lattice::Lattice(const Lattice &other) :
  vol(other.vol), sun(other.sun), nd(other.nd), dim(other.dim), lat_sizes(other.lat_sizes)
{
}

Lattice::~Lattice()
{
}

// ================================================= //
// ============    Private members     ============= //
// ================================================= //

// ================================================= //
// ========  Memeber Operator overloading  ========= //
// ================================================= //

bool Lattice::operator==(const Lattice &rhs)
{
  bool res = (dim == rhs.dim) && (sun == rhs.sun) && (nd == rhs.nd);

  if (res == true)
  {
    for(uint i = 0; i < dim; i++)
      res = res && (lat_sizes[i] == rhs.lat_sizes[i]);
  }

  return res;
}

bool Lattice::operator!=(const Lattice &rhs)
{
  return !( (*this) == rhs );
}

// ================================================= //
// ===============  Public members  ================ //
// ================================================= //

uint Lattice::Volume() const
{
  return vol;
}

uint Lattice::Dim() const
{
  return dim;
}

uint Lattice::SUNrank() const
{
  return sun;
}

uint Lattice::Ndirac() const
{
  return nd;
}

void Lattice::LatticeSizes(VECTOR<uint> &L) const
{
  DEBUG_ASSERT(L.size() == dim);

  for(uint i = 0; i < dim; i++)
    L[i] = lat_sizes[i];
}

uint Lattice::LatticeSize(const uint mu) const
{
  DEBUG_ASSERT(mu < dim);

  return lat_sizes[mu];
}

uint Lattice::SiteIndex(const VECTOR<int> &x) const
{
  uint idx = 0;
  uint dimsize = 1;

  for(uint i = 0; i < dim; i++)
  {
    int div = x[i] / (int)lat_sizes[i];
    int rem = x[i] - div * (int)lat_sizes[i];
    uint lat_x = (uint)( (rem < 0) ? (rem + (int)lat_sizes[i]) : rem );

    idx += lat_x * dimsize;
    dimsize *= lat_sizes[i];
  }

  return idx;
}

void Lattice::SiteCoordinates(VECTOR<int> &x, const uint idx) const
{
  uint midx = idx;
  uint dimsize = vol;

  for(uint j = 0; j < dim; j++)
  {
    uint i = dim - 1 - j;
    dimsize /= lat_sizes[i];
    x[i] = midx / dimsize;
    midx = midx % dimsize;
  }
}

uint Lattice::LinkIndex(const uint sidx, const uint mu) const
{
  return (mu + dim * sidx);
}

uint Lattice::LinkIndex(const VECTOR<int> &x, const uint mu) const
{
  uint idx = SiteIndex(x);
  return LinkIndex(idx, mu);
}


void Lattice::LinkCoordinates(uint &sidx, uint &mu, const uint idx) const
{
  sidx = idx / dim;
  mu = idx % dim;
}

void Lattice::LinkCoordinates(VECTOR<int> &x, uint &mu, const uint idx) const
{
  uint site_idx;

  LinkCoordinates(site_idx, mu, idx);
  SiteCoordinates(x, site_idx);
}

uint Lattice::SiteIndexForward(const uint idx, const uint mu) const
{
  DEBUG_ASSERT(mu < dim);

  uint dimsize = 1;
  uint x_mu = 0;
  uint lower_idx = 0;
  uint upper_idx = 0;
  uint res = 0;

  for(uint i = 0; i < mu; i++)
    dimsize *= lat_sizes[i];

  lower_idx = idx % dimsize;
  upper_idx = idx / dimsize;
  x_mu = upper_idx % lat_sizes[mu];

  upper_idx -= x_mu;
  x_mu = (x_mu + 1) % lat_sizes[mu];
  upper_idx += x_mu;

  res = lower_idx + dimsize * upper_idx;

  return res;
}

uint Lattice::SiteIndexBackward(const uint idx, const uint mu) const
{
  DEBUG_ASSERT(mu < dim);

  uint dimsize = 1;
  uint x_mu = 0;
  uint lower_idx = 0;
  uint upper_idx = 0;
  uint res = 0;

  for(uint i = 0; i < mu; i++)
    dimsize *= lat_sizes[i];

  lower_idx = idx % dimsize;
  upper_idx = idx / dimsize;
  x_mu = upper_idx % lat_sizes[mu];

  upper_idx -= x_mu;
  x_mu = (x_mu == 0) ? lat_sizes[mu] - 1 : x_mu - 1;
  upper_idx += x_mu;

  res = lower_idx + dimsize * upper_idx;

  return res;
}

uint Lattice::LinkIndexForward(const uint idx, const uint mu) const
{
  DEBUG_ASSERT(mu < dim);

  uint site_idx;
  uint link_mu;
  uint res;

  LinkCoordinates(site_idx, link_mu, idx);
  site_idx = SiteIndexForward(site_idx, mu);
  res = LinkIndex(site_idx, link_mu);

  return res;
}

uint Lattice::LinkIndexBackward(const uint idx, const uint mu) const
{
  DEBUG_ASSERT(mu < dim);

  uint site_idx;
  uint link_mu;
  uint res;

  LinkCoordinates(site_idx, link_mu, idx);
  site_idx = SiteIndexBackward(site_idx, mu);
  res = LinkIndex(site_idx, link_mu);

  return res;
}

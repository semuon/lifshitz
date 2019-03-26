#ifndef _LATTICE_NEW_H_
#define _LATTICE_NEW_H_

#include "common.h"

class Lattice
{
private:
  uint vol;
  uint sun;
  uint nd;
  uint dim;
  VECTOR<uint> lat_sizes;

public:
  Lattice(const VECTOR<uint> &L, uint ndim, uint sunrank, uint ndirac);
  Lattice(const Lattice &other, uint sunrank, uint ndirac);
  Lattice(const Lattice &other);
  ~Lattice();

  bool operator==(const Lattice &rhs);
  bool operator!=(const Lattice &rhs);

  uint Volume() const;
  uint Dim() const;
  uint SUNrank() const;
  uint Ndirac() const;
  void LatticeSizes(VECTOR<uint> &L) const;
  uint LatticeSize(const uint mu) const;

  uint SiteIndex(const VECTOR<int> &x) const;
  void SiteCoordinates(VECTOR<int> &x, const uint idx) const;

  uint LinkIndex(const uint sidx, const uint mu) const;
  uint LinkIndex(const VECTOR<int> &x, const uint mu) const;
  void LinkCoordinates(uint &sidx, uint &mu, const uint idx) const;
  void LinkCoordinates(VECTOR<int> &x, uint &mu, const uint idx) const;

  uint SiteIndexForward(const uint idx, const uint mu) const;
  uint SiteIndexBackward(const uint idx, const uint mu) const;

  uint LinkIndexForward(const uint idx, const uint mu) const;
  uint LinkIndexBackward(const uint idx, const uint mu) const;
};

#endif
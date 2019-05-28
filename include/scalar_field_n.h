#ifndef SCALAR_FIELD_N_H_
#define SCALAR_FIELD_N_H_

#include <common.h>
#include <lattice.h>
#include <linear_vector.h>

template <class T> class TScalarFieldN : public TLinearVector<T>
{
protected:
  uint vol;
  uint n;

public:
  TScalarFieldN(const Lattice &lat, uint n)
  : TLinearVector<T>(lat, n * lat.Volume())
  , vol(lat.Volume()), n(n)
  { }

  TScalarFieldN(const Lattice &lat, uint n, T* memptr)
  : TLinearVector<T>(lat, n * lat.Volume(), memptr)
  , vol(lat.Volume()), n(n)
  { }

  TScalarFieldN(const TScalarFieldN<T> &other)
  : TLinearVector<T>(other)
  , vol(other.vol), n(other.n)
  { }

  ~TScalarFieldN() { }

  T& operator()(const VECTOR<int> &x, const uint i) const
  {
    return (*this)[ElementIndex(x, i)];
  }

  T& operator()(const uint idx, const uint i) const
  {
    return (*this)[ElementIndex(idx, i)];
  }

  TScalarFieldN<T>& operator=(const TScalarFieldN<T>& rhs)
  {
    TLinearVector<T>::operator=(rhs);
    return *this;
  }

  uint ElementIndex(const uint idx, const uint i) const
  {
    DEBUG_ASSERT(i < n);
    return TScalarFieldN<T>::ElementIndex(idx, i, this->lat);
  }

  uint ElementIndex(const VECTOR<int> &x, const uint i) const
  {
    DEBUG_ASSERT(i < n);
    return TScalarFieldN<T>::ElementIndex(x, i, this->lat);
  }

  static uint ElementIndex(const uint idx, const uint i, const Lattice &lat)
  {
    return (idx + lat.Volume() * i);
  }

  static uint ElementIndex(const VECTOR<int> &x, const uint i, const Lattice &lat)
  {
    uint idx = lat.SiteIndex(x);
    return TScalarFieldN<T>::ElementIndex(idx, i, lat);
  }

  uint Volume() const { return vol; };
  uint N() const { return n; };
};

typedef TScalarFieldN<t_complex> ScalarFieldN;
typedef TScalarFieldN<double> RealScalaFieldN;

#endif /* SCALAR_FIELD_N_H_ */ 

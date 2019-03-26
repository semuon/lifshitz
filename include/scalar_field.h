#ifndef SCALAR_FIELD_H_
#define SCALAR_FIELD_H_

#include <common.h>
#include <lattice.h>
#include <linear_vector.h>

template <class T> class TScalarField : public TLinearVector<T>
{
protected:
  uint vol;

public:
  TScalarField(const Lattice &lat)
  : TLinearVector<T>(lat, lat.Volume())
  , vol(lat.Volume())
  { }

  TScalarField(const Lattice &lat, T* memptr)
  : TLinearVector<T>(lat, lat.Volume(), memptr)
  , vol(lat.Volume())
  { }

  TScalarField(const TScalarField<T> &other)
  : TLinearVector<T>(other)
  , vol(other.vol)
  { }

  ~TScalarField() { }

  T& operator()(const uint idx) const
  {
    DEBUG_ASSERT(idx < vol);
    return (*this)[idx];
  }

  T& operator()(const VECTOR<int> &x) const
  {
    uint idx = TLinearVector<T>::lat.SiteIndex(x);
    return (*this)[idx];
  }

  TScalarField<T>& operator=(const TScalarField<T>& rhs)
  {
    TLinearVector<T>::operator=(rhs);
    return *this;
  }

  uint Volume() const { return vol; };
};

typedef TScalarField<t_complex> ScalarField;

#endif /* SCALAR_FIELD_H_ */ 

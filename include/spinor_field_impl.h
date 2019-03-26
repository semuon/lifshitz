#ifndef SPINOR_FIELD_IMPL_H_
#define SPINOR_FIELD_IMPL_H_

#include "profiler.h"

// ================================================= //
// =================  Constructors  ================ //
// ================================================= //

template <class T> TSpinorField<T>::TSpinorField(const Lattice &lat) 
  : TLinearVector<T>(lat, lat.SUNrank() * lat.Ndirac() * lat.Volume())
  , vol(lat.Volume()), nd(lat.Ndirac()), sun(lat.SUNrank())
{
  pGlobalProfiler.IncrementCounter("TSpinorField: constructor with internal memory");
}

template <class T> TSpinorField<T>::TSpinorField(const Lattice &lat, T* memptr)
  : TLinearVector<T>(lat, lat.SUNrank() * lat.Ndirac() * lat.Volume(), memptr)
  , vol(lat.Volume()), nd(lat.Ndirac()), sun(lat.SUNrank())
{
  pGlobalProfiler.IncrementCounter("TSpinorField: constructor with external memory");
}

template <class T> TSpinorField<T>::TSpinorField(const TSpinorField<T> &other)
  : TLinearVector<T>(other)
  , vol(other.vol), nd(other.nd), sun(other.sun)
{
  pGlobalProfiler.IncrementCounter("TSpinorField: Copy constructor");
}

template <class T> TSpinorField<T>::~TSpinorField()
{
}

// ================================================= //
// ============    Private members     ============= //
// ================================================= //

// ================================================= //
// =========  Member Operator overloading  ========= //
// ================================================= //

template <class T> T& TSpinorField<T>::operator()(const uint idx, const uint d, const uint c) const
{
  uint i = ElementIndex(idx, d, c);

  return TBaseLinearVector<T>::mem[i];
}

template <class T> T& TSpinorField<T>::operator()(const VECTOR<int> &x, const uint d, const uint c) const
{
  uint i = ElementIndex(x, d, c);

  return TBaseLinearVector<T>::mem[i];
}

template <class T> TSpinorField<T>& TSpinorField<T>::operator=(const TSpinorField<T> &rhs)
{
  TLinearVector<T>::operator=(rhs);

  return *this;
}

// ================================================= //
// ======  Non-Memeber Operator overloading  ======= //
// ================================================= //

// ================================================= //
// ===============  Public members  ================ //
// ================================================= //

template <class T> uint TSpinorField<T>::ElementIndex(const uint idx, const uint d, const uint c) const
{
  return (c + sun * d + sun * nd * idx);
}

template <class T> uint TSpinorField<T>::ElementIndex(const VECTOR<int> &x, const uint d, const uint c) const
{
  uint idx = TLinearVector<T>::lat.SiteIndex(x);
  return ElementIndex(idx, d, c);
}

template <class T> uint TSpinorField<T>::ElementIndex(const uint idx, const uint d, const uint c, const Lattice &lat)
{
  return (c + lat.SUNrank() * d + lat.SUNrank() * lat.Ndirac() * idx);
}

template <class T> uint TSpinorField<T>::ElementIndex(const VECTOR<int> &x, const uint d, const uint c, const Lattice &lat)
{
  uint idx = lat.SiteIndex(x);
  return ElementIndex(idx, d, c);
}

// ================================================= //
// ==========  Operator specializations  =========== //
// ================================================= //

#endif /* VECTOR_FIELD_IMPL_H_ */ 

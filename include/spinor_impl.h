#ifndef SPINOR_IMPL_H_
#define SPINOR_IMPL_H_

// ================================================= //
// =================  Constructors  ================ //
// ================================================= //

template <class T> TSpinor<T>::TSpinor(const Lattice &lat) 
  : TLinearVector<T>(lat, lat.SUNrank() * lat.Ndirac())
  , nd(lat.Ndirac()), sun(lat.SUNrank()), color_vecs(NULL)
{
  InitClass();
}

template <class T> TSpinor<T>::TSpinor(const Lattice &lat, T* memptr)
  : TLinearVector<T>(lat, lat.SUNrank() * lat.Ndirac(), memptr)
  , nd(lat.Ndirac()), sun(lat.SUNrank()), color_vecs(NULL)
{
  InitClass();
}

template <class T> TSpinor<T>::TSpinor(const TSpinor<T> &other)
  : TLinearVector<T>(other)
  , nd(other.nd), sun(other.sun), color_vecs(NULL)
{
  InitClass();
}

template <class T> TSpinor<T>::~TSpinor()
{
  FinilizeClass();
}

// ================================================= //
// ============    Private members     ============= //
// ================================================= //

template <class T> void TSpinor<T>::InitClass()
{
  char *color_mem = new char [sizeof(TColorVector<T>) * nd];
  color_vecs = reinterpret_cast<TColorVector<T>*>(color_mem);

  for(uint d = 0; d < nd; d++)
  {
    T* color_ptr = TBaseLinearVector<T>::mem + d * sun;
    new (color_vecs + d) TColorVector<T>(TLinearVector<T>::lat, color_ptr);
  }
}

template <class T> void TSpinor<T>::FinilizeClass()
{
  for(uint d = 0; d < nd; d++)
    color_vecs[d].~TColorVector();

  delete[] reinterpret_cast<char *>(color_vecs);
}

// ================================================= //
// =========  Member Operator overloading  ========= //
// ================================================= //

template <class T> TColorVector<T>& TSpinor<T>::operator()(const uint d) const
{
  DEBUG_ASSERT(d < nd);

  return color_vecs[d];
}

template <class T> T& TSpinor<T>::operator()(const uint d, const uint c) const
{
  DEBUG_ASSERT(d < nd && c < sun);

  uint i = ElementIndex(d, c);

  return TBaseLinearVector<T>::mem[i];
}

template <class T> TSpinor<T>& TSpinor<T>::operator=(const TSpinor<T>& rhs)
{
  DEBUG_ASSERT(sun == rhs.sun && nd == rhs.nd);

  TBaseLinearVector<T>::operator=(rhs);

  return *this;
}

// ================================================= //
// ======  Non-Memeber Operator overloading  ======= //
// ================================================= //

// ================================================= //
// ===============  Public members  ================ //
// ================================================= //

template <class T> uint TSpinor<T>::ElementIndex(const uint d, const uint c) const
{
  return (c + sun * d);
}

template <class T> uint TSpinor<T>::ElementIndex(const uint d, const uint c, const Lattice &lat)
{
  return (c + lat.SUNrank() * d);
}

// ================================================= //
// ==========  Operator specializations  =========== //
// ================================================= //


#endif /* SPINOR_IMPL_H_ */ 

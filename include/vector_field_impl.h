#ifndef VECTOR_FIELD_IMPL_H_
#define VECTOR_FIELD_IMPL_H_

// ================================================= //
// =================  Constructors  ================ //
// ================================================= //

template <class T> TVectorField<T>::TVectorField(const Lattice &lat) 
  : TLinearVector<T>(lat, lat.Dim() * lat.Volume())
  , vol(lat.Volume()), dim(lat.Dim()), space_vecs(NULL)
{
  InitClass();
}

template <class T> TVectorField<T>::TVectorField(const Lattice &lat, T* memptr)
  : TLinearVector<T>(lat, lat.Dim() * lat.Volume(), memptr)
  , vol(lat.Volume()), dim(lat.Dim()), space_vecs(NULL)
{
  InitClass();
}

template <class T> TVectorField<T>::TVectorField(const TVectorField<T> &other)
  : TLinearVector<T>(other)
  , vol(other.vol), dim(other.dim), space_vecs(NULL)
{
  InitClass();
}

template <class T> TVectorField<T>::~TVectorField()
{
  FinilizeClass();
}

// ================================================= //
// ============    Private members     ============= //
// ================================================= //

template <class T> void TVectorField<T>::InitClass()
{
  // Allocate vectors
  char *space_mem = new char [sizeof(TSpaceVector<T>) * vol];
  space_vecs = reinterpret_cast<TSpaceVector<T>*>(space_mem);

  for(uint i = 0; i < vol; i++)
  {
    T* space_ptr = TBaseLinearVector<T>::mem + i * dim;
    new (space_vecs + i) TSpaceVector<T>(TLinearVector<T>::lat, space_ptr);
  }
}

template <class T> void TVectorField<T>::FinilizeClass()
{
  // Delete vectors
  for(uint i = 0; i < vol; i++)
    space_vecs[i].~TSpaceVector();

  delete[] reinterpret_cast<char *>(space_vecs);
}

// ================================================= //
// ========  Memeber Operator overloading  ========= //
// ================================================= //

template <class T> TSpaceVector<T>& TVectorField<T>::operator()(const uint idx) const
{
  return space_vecs[idx];
}

template <class T> TSpaceVector<T>& TVectorField<T>::operator()(const VECTOR<int> &x) const
{
  uint idx = TLinearVector<T>::lat.SiteIndex(x);

  return space_vecs[idx];
}

template <class T> T& TVectorField<T>::operator()(const uint idx, const uint mu) const
{
  uint i = ElementIndex(idx, mu);

  return TBaseLinearVector<T>::mem[i];
}

template <class T> T& TVectorField<T>::operator()(const VECTOR<int> &x, const uint mu) const
{
  uint i = ElementIndex(x, mu);

  return TBaseLinearVector<T>::mem[i];
}

template <class T> TVectorField<T>& TVectorField<T>::operator=(const TVectorField<T> &rhs)
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

template <class T> uint TVectorField<T>::ElementIndex(const uint idx, const uint mu) const
{
  return TLinearVector<T>::lat.LinkIndex(idx, mu);
}

template <class T> uint TVectorField<T>::ElementIndex(const VECTOR<int> &x, const uint mu) const
{
  return TLinearVector<T>::lat.LinkIndex(x, mu);
}

template <class T> uint TVectorField<T>::ElementIndex(const uint idx, const uint mu, const Lattice &lat)
{
  return lat.LinkIndex(idx, mu);
}

template <class T> uint TVectorField<T>::ElementIndex(const VECTOR<int> &x, const uint mu, const Lattice &lat)
{
  return lat.LinkIndex(x, mu);
}

// ================================================= //
// ==========  Operator specializations  =========== //
// ================================================= //

#endif /* VECTOR_FIELD_IMPL_H_ */ 

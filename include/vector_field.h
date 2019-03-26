#ifndef VECTOR_FIELD_H_
#define VECTOR_FIELD_H_

#include <common.h>
#include <lattice.h>
#include <linear_vector.h>

template <class T> class TSpaceVector : public TLinearVector<T>
{
protected:
  uint dim;

public:
  TSpaceVector(const Lattice &lat)
  : TLinearVector<T>(lat, lat.Dim())
  , dim(lat.Dim())
  { }

  TSpaceVector(const Lattice &lat, T* memptr)
  : TLinearVector<T>(lat, lat.Dim(), memptr)
  , dim(lat.Dim())
  { }

  TSpaceVector(const TSpaceVector<T> &other)
  : TLinearVector<T>(other)
  , dim(other.dim)
  { }

  ~TSpaceVector() { }

  T& operator()(const uint mu) const
  {
    DEBUG_ASSERT(mu < dim);
    return (*this)[mu];
  }

  TSpaceVector<T>& operator=(const TSpaceVector<T>& rhs)
  {
    DEBUG_ASSERT(dim = rhs.dim);
    TLinearVector<T>::operator=(rhs);
    return *this;
  }

  uint Dim() const { return dim; }
};

template <class T> class TVectorField : public TLinearVector<T>
{
protected:
  uint vol;
  uint dim;
  TSpaceVector<T> *space_vecs;

  void InitClass();
  void FinilizeClass();

public:
  TVectorField(const Lattice &lat);
  TVectorField(const Lattice &lat, T* memptr);
  TVectorField(const TVectorField<T> &other);
  ~TVectorField();

  TSpaceVector<T>& operator()(const uint idx)                           const;
  TSpaceVector<T>& operator()(const VECTOR<int> &x)                const;
                T& operator()(const uint idx, const uint mu)            const;
                T& operator()(const VECTOR<int> &x, const uint mu) const;

  TVectorField<T>& operator=(const TVectorField<T> &rhs);

         uint ElementIndex(const uint idx, const uint mu)               const;
         uint ElementIndex(const VECTOR<int> &x, const uint mu)    const;

  static uint ElementIndex(const uint idx, const uint mu, const Lattice &lat);
  static uint ElementIndex(const VECTOR<int> &x, const uint mu, const Lattice &lat);

  uint Volume() const { return vol; }
  uint Dim()    const { return dim; }
};

#include <vector_field_impl.h>

typedef TSpaceVector<int>  Coordinates;

typedef TSpaceVector<t_complex> SpaceVector;
typedef TVectorField<t_complex> VectorField;

#endif /* VECTOR_FIELD_H_ */ 

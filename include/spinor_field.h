#ifndef SPINOR_FIELD_H_
#define SPINOR_FIELD_H_

#include <common.h>
#include <lattice.h>
#include <linear_vector.h>

template <class T> class TColorVector : public TLinearVector<T>
{
protected:
  uint sun;

public:
  TColorVector(const Lattice &lat)
  : TLinearVector<T>(lat, lat.SUNrank())
  , sun(lat.SUNrank())
  { }

  TColorVector(const Lattice &lat, T* memptr)
  : TLinearVector<T>(lat, lat.SUNrank(), memptr)
  , sun(lat.SUNrank())
  { }

  TColorVector(const TColorVector<T> &other)
  : TLinearVector<T>(other)
  , sun(other.sun)
  { }

  ~TColorVector() { }

  T& operator()(const uint c) const
  {
    DEBUG_ASSERT(c < sun);
    return (*this)[c];
  }

  TColorVector<T>& operator=(const TColorVector<T>& rhs)
  {
    DEBUG_ASSERT(sun == rhs.sun);
    TBaseLinearVector<T>::operator=(rhs);
    return *this;
  }

  uint SUNrank() const { return sun; }
};

template <class T> class TSpinor : public TLinearVector<T>
{
protected:
  uint nd;
  uint sun;
  TColorVector<T> *color_vecs;

  void InitClass();
  void FinilizeClass();

public:
  TSpinor(const Lattice &lat);
  TSpinor(const Lattice &lat, T* memptr);
  TSpinor(const TSpinor<T> &other);
  ~TSpinor();

  TColorVector<T>& operator()(const uint d) const;
                T& operator()(const uint d, const uint c) const;

  TSpinor<T>& operator=(const TSpinor<T>& rhs);

         uint ElementIndex(const uint d, const uint c) const;
  static uint ElementIndex(const uint d, const uint c, const Lattice &lat);

  uint SUNrank() const { return sun; }
  uint Ndirac()  const { return nd;  }
};

template <class T> class TSpinorField : public TLinearVector<T>
{
protected:
  const uint vol;
  const uint nd;
  const uint sun;

public:
  TSpinorField(const Lattice &lat);
  TSpinorField(const Lattice &lat, T* memptr);
  TSpinorField(const TSpinorField<T> &other);
  ~TSpinorField();

  T& operator()(const uint idx, const uint d, const uint c)            const;
  T& operator()(const VECTOR<int> &x, const uint d, const uint c) const;

  TSpinorField<T>& operator=(const TSpinorField<T> &rhs);

         uint ElementIndex(const uint idx, const uint d, const uint c)               const;
         uint ElementIndex(const VECTOR<int> &x, const uint d, const uint c)    const;

  static uint ElementIndex(const uint idx, const uint d, const uint c, const Lattice &lat);
  static uint ElementIndex(const VECTOR<int> &x, const uint d, const uint c, const Lattice &lat);

  uint SUNrank() const { return sun; }
  uint Ndirac()  const { return nd;  }
  uint Volume()  const { return vol; }
};

#include <spinor_impl.h>
#include <spinor_field_impl.h>

typedef TColorVector<t_complex> ColorVector;
typedef TSpinor<t_complex>      Spinor;
typedef TSpinorField<t_complex> SpinorField;

#endif /* VECTOR_FIELD_H_ */ 

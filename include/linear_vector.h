#ifndef LINEAR_VECTOR_H_
#define LINEAR_VECTOR_H_

#include <common.h>
#include <lattice.h>
#include <cblas_templates.h>

// ====================================== //
// ==== Linear vector template class ==== //
// ====================================== //

// Deduce the return type for Norm functions of TBaseLinearVector
// For real-valued vectors Norm is real
// For complex-valued vectors - also real
template<class T>
struct TBaseLinearVector_NormType{ typedef T type; };
template<>
struct TBaseLinearVector_NormType<t_complex>{ typedef double type; };
template<>
struct TBaseLinearVector_NormType<t_complex_float>{ typedef float type; };

template <typename T> class TBaseLinearVector
{
protected:
  T *mem;
  uint mem_dim;
  bool flag_self_mem;

  void AllocateStorage(const uint n, T *memptr);

public:
  TBaseLinearVector(const uint n);
  TBaseLinearVector(const uint n, T *memptr);
  TBaseLinearVector(const TBaseLinearVector<T> &other);
  ~TBaseLinearVector();

  T& operator[](const uint i) const;
  operator T*() const;

  TBaseLinearVector<T>& operator=(const TBaseLinearVector<T> &rhs);

  TBaseLinearVector<T>& operator+=(const TBaseLinearVector<T>& rhs);
  TBaseLinearVector<T>& operator-=(const TBaseLinearVector<T>& rhs);
  TBaseLinearVector<T>& operator*=(const T& rhs);

  uint Count() const;
  bool IsMemoryExternal() const;

  T* DataPtr() const { return mem; }

  // ====================================== //
  // ====== Linear vector operations ====== //
  // ====================================== //

  // ------------------ //
  // A <= 0             //
  // ------------------ //
  void Clear();

  // ------------------ //
  // A <-> B            //
  // ------------------ //
  void Swap(TBaseLinearVector<T>& B);

  // ------------------ //
  // Returns |A|        //
  // ------------------ //
  typename TBaseLinearVector_NormType<T>::type Norm() const;

  // ------------------ //
  // Returns |A|^2      //
  // ------------------ //
  typename TBaseLinearVector_NormType<T>::type NormSqr() const;

  // ------------------ //
  // Make |A| = 1       //
  // ------------------ //
  void Normalize();

  // ------------------ //
  // Returns \sum v_x   //
  // ------------------ //
  T Total() const;

  // ------------------ //
  // A <= A + B         //
  // ------------------ //
  void Add(const TBaseLinearVector<T>& B);

  // ------------------ //
  // A <= A - B         //
  // ------------------ //
  void Substract(const TBaseLinearVector<T>& B);

  // ------------------ //
  // A <= b * B + c * C //
  // ------------------ //
  void Assign_bB_plus_cC(const T b, const TBaseLinearVector<T> &B, const T c, const TBaseLinearVector<T> &C);

  // ------------------ //
  // A <= a * A + b * B //
  // ------------------ //
  void Assign_aA_plus_bB(const T a, const T b, const TBaseLinearVector<T> &B);

  // ------------------ //
  // A <= b * B         //
  // ------------------ //
  void Assign_bB(const T b, const TBaseLinearVector<T> &B);

  // ------------------ //
  // A += b * B         //
  // ------------------ //
  void Add_bB(const T b, const TBaseLinearVector<T> &B);

  // ------------------ //
  // A *= a             //
  // ------------------ //
  void MultiplyBy(const T a);

  // ------------------ //
  // < A, B >           //
  // ------------------ //
  T ScalarProduct(const TBaseLinearVector<T> &B) const;
};

template <typename T, template <typename Z> class X> X<T> operator+(X<T> lhs, const X<T>& rhs);
template <typename T, template <typename Z> class X> X<T> operator-(X<T> lhs, const X<T>& rhs);
template <typename T, template <typename Z> class X> X<T> operator*(X<T> lhs, const T& rhs);
template <typename T, template <typename Z> class X> X<T> operator*(const T& lhs, X<T> rhs);

template <typename T> T operator*(const TBaseLinearVector<T> &lhs, const TBaseLinearVector<T> &rhs);

template <template <class t_complex> class X> X<t_complex> operator*(X<t_complex> rhs, const double& lhs);
template <template <class t_complex> class X> X<t_complex> operator*(const double& lhs, X<t_complex> rhs);

#include <linear_vector_impl.h>

template <class T> class TLinearVector : public TBaseLinearVector<T>
{
protected:
  Lattice lat;

public:
  TLinearVector(const Lattice &lat, const uint n)
  : TBaseLinearVector<T>(n)
  , lat(lat)
  {}

  TLinearVector(const Lattice &lat, const uint n, T *memptr)
  : TBaseLinearVector<T>(n, memptr)
  , lat(lat)
  {}

  TLinearVector(const TLinearVector<T> &other)
  : TBaseLinearVector<T>(other)
  , lat(other.lat)
  {}

  ~TLinearVector() {}

  TLinearVector<T>& operator=(const TLinearVector<T>& rhs)
  {
    DEBUG_ASSERT(lat == rhs.lat);
    TBaseLinearVector<T>::operator=(rhs);
    return *this;
  }

  Lattice GetLattice() const { return lat; }
};

typedef TBaseLinearVector<t_complex> BaseLinearVector;
typedef TLinearVector<t_complex>     LinearVector;

#endif /* LINEAR_VECTOR_H_ */ 

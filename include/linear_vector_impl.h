#ifndef BASE_LINEAR_VECTOR_IMPL_H_
#define BASE_LINEAR_VECTOR_IMPL_H_


// ================================================= //
// =================  Constructors  ================ //
// ================================================= //

template <typename T> TBaseLinearVector<T>::TBaseLinearVector(const uint n)
  : mem(NULL), mem_dim(0), flag_self_mem(false)
{
  AllocateStorage(n, NULL);
}

template <typename T> TBaseLinearVector<T>::TBaseLinearVector(const uint n, T *memptr)
  : mem(NULL), mem_dim(0), flag_self_mem(false)
{
  AllocateStorage(n, memptr);
}

template <typename T> TBaseLinearVector<T>::TBaseLinearVector(const TBaseLinearVector<T> &other)
  : mem(NULL), mem_dim(0), flag_self_mem(other.flag_self_mem)
{
  AllocateStorage(other.mem_dim, NULL);

  for(uint i = 0; i < mem_dim; i++)
    mem[i] = other.mem[i];
}

template <typename T>TBaseLinearVector<T>::~TBaseLinearVector()
{
  if (flag_self_mem)
    delete[] mem;
}

// ================================================= //
// ============    Private members     ============= //
// ================================================= //

template <typename T> void TBaseLinearVector<T>::AllocateStorage(const uint n, T *memptr)
{
  mem_dim = n;

  if (flag_self_mem && mem != NULL)
    delete[] mem;

  if (memptr == NULL)
  {
    mem = new T[mem_dim];
    flag_self_mem = true;
  }
  else
  {
    mem = memptr;
    flag_self_mem = false;
  }
}

// ================================================= //
// ========  Memeber Operator overloading  ========= //
// ================================================= //

template <typename T> T& TBaseLinearVector<T>::operator[](const uint i) const
{
  DEBUG_ASSERT(i < mem_dim);

  return mem[i];
}

template <typename T> TBaseLinearVector<T>::operator T*() const
{
  return mem;
}

template <typename T> TBaseLinearVector<T>& TBaseLinearVector<T>::operator=(const TBaseLinearVector<T>& rhs)
{
  DEBUG_ASSERT( (this != &rhs) && (rhs.mem_dim == mem_dim) );

  for(uint i = 0; i < mem_dim; i++)
    mem[i] = rhs.mem[i];

  return *this;
}

template <typename T> TBaseLinearVector<T>& TBaseLinearVector<T>::operator+=(const TBaseLinearVector<T> &rhs)
{
  Add(rhs);
  return *this;
}

template <typename T> TBaseLinearVector<T>& TBaseLinearVector<T>::operator-=(const TBaseLinearVector<T> &rhs)
{
  Substract(rhs);
  return *this;
}

template <typename T> TBaseLinearVector<T>& TBaseLinearVector<T>::operator*=(const T &rhs)
{
  MultiplyBy(rhs);
  return *this;
}

// ================================================= //
// ======  Non-Memeber Operator overloading  ======= //
// ================================================= //

template <typename T, template <typename Z> class X> X<T> operator+(X<T> lhs, const X<T>& rhs)
{
  *(static_cast<TBaseLinearVector<T>*>(&lhs)) += static_cast<TBaseLinearVector<T> >(rhs);
  return lhs;
}

template <typename T, template <typename Z> class X> X<T> operator-(X<T> lhs, const X<T>& rhs)
{
  *(static_cast<TBaseLinearVector<T>*>(&lhs)) -= static_cast<TBaseLinearVector<T> >(rhs);
  return lhs;
}

template <typename T, template <typename Z> class X> X<T> operator*(X<T> lhs, const T& rhs)
{
  *(static_cast<TBaseLinearVector<T>*>(&lhs)) *= rhs;
  return lhs;
}

template <typename T, template <typename Z> class X> X<T> operator*(const T& lhs, X<T> rhs)
{
  *(static_cast<TBaseLinearVector<T>*>(&rhs)) *= lhs;
  return rhs;
}

template <typename T> T operator*(const TBaseLinearVector<T> &lhs, const TBaseLinearVector<T> &rhs)
{
  return lhs.ScalarProduct(rhs);
}

template <template <class t_complex> class X> X<t_complex> operator*(X<t_complex> lhs, const double& rhs)
{
  *(static_cast<TBaseLinearVector<t_complex>*>(&lhs)) *= (t_complex)rhs;
  return lhs;
}

template <template <class t_complex> class X> X<t_complex> operator*(const double& lhs, X<t_complex> rhs)
{
  *(static_cast<TBaseLinearVector<t_complex>*>(&rhs)) *= (t_complex)lhs;
  return rhs;
}

// ================================================= //
// ===============  Public members  ================ //
// ================================================= //

template <typename T> void TBaseLinearVector<T>::Clear()
{
  for(uint i = 0; i < mem_dim; i++)
    mem[i] = 0;
}

template <typename T> void TBaseLinearVector<T>::Swap(TBaseLinearVector<T>& B)
{
  CblasHelper<T>::cblas_swap(mem_dim, mem, 1, B.mem, 1);
}

template <typename T> uint TBaseLinearVector<T>::Count() const
{
  return mem_dim;
}

template <typename T> bool TBaseLinearVector<T>::IsMemoryExternal() const
{
  return !flag_self_mem;
}

// Returns norm |A|
template <typename T> typename TBaseLinearVector_NormType<T>::type TBaseLinearVector<T>::Norm() const
{
  return CblasHelper<T>::cblas_nrm2(mem_dim, mem, 1);
}

// Returns square norm |A|^2
template <typename T> typename TBaseLinearVector_NormType<T>::type TBaseLinearVector<T>::NormSqr() const
{
  typename TBaseLinearVector_NormType<T>::type res = Norm();
  return res * res;
}

// Make |A| = 1
template <typename T> void TBaseLinearVector<T>::Normalize()
{
  typename TBaseLinearVector_NormType<T>::type norm = Norm();
  MultiplyBy(1 / norm);
}

// Returns sum of all elements
template <typename T> T TBaseLinearVector<T>::Total() const
{
  T sum = 0;

  for(uint i = 0; i < mem_dim; i++)
    sum += mem[i];

  return sum;
}

// ------------------ //
// A <= A + B         //
// ------------------ //
template <typename T> void TBaseLinearVector<T>::Add(const TBaseLinearVector<T>& B)
{
  DEBUG_ASSERT(mem_dim == B.mem_dim);

  const T b = 1;
  CblasHelper<T>::cblas_axpy(mem_dim, &b, B.mem, 1, mem, 1);
}

// ------------------ //
// A <= A - B         //
// ------------------ //
template <typename T> void TBaseLinearVector<T>::Substract(const TBaseLinearVector<T>& B)
{
  DEBUG_ASSERT(mem_dim == B.mem_dim);

  const T b = -1;
  CblasHelper<T>::cblas_axpy(mem_dim, &b, B.mem, 1, mem, 1);
}

// ------------------ //
// A <= b * B + c * C //
// ------------------ //
template <typename T> void TBaseLinearVector<T>::Assign_bB_plus_cC(const T b, const TBaseLinearVector<T> &B, const T c, const TBaseLinearVector<T> &C)
{
  DEBUG_ASSERT(mem_dim == B.mem_dim && mem_dim == C.mem_dim);

  CblasHelper<T>::cblas_copy(mem_dim, B.mem, 1, mem, 1);     // A = B
  CblasHelper<T>::cblas_scal(mem_dim, &b, mem, 1);           // A = b * B
  CblasHelper<T>::cblas_axpy(mem_dim, &c, C.mem, 1, mem, 1); // A = b * B + c * C
}

// ------------------ //
// A <= a * A + b * B //
// ------------------ //
template <typename T> void TBaseLinearVector<T>::Assign_aA_plus_bB(const T a, const T b, const TBaseLinearVector<T> &B)
{
  DEBUG_ASSERT(mem_dim == B.mem_dim);

  CblasHelper<T>::cblas_scal(mem_dim, &a, mem, 1);           // A = a * A
  CblasHelper<T>::cblas_axpy(mem_dim, &b, B.mem, 1, mem, 1); // A = a * A + b * B
}

// ------------------ //
// A <= b * B         //
// ------------------ //
template <typename T> void TBaseLinearVector<T>::Assign_bB(const T b, const TBaseLinearVector<T> &B)
{
  DEBUG_ASSERT(mem_dim == B.mem_dim);

  CblasHelper<T>::cblas_copy(mem_dim, B.mem, 1, mem, 1);    // A = B
  CblasHelper<T>::cblas_scal(mem_dim, &b, mem, 1);          // A = b * B
}

// ------------------ //
// A += b * B         //
// ------------------ //
template <typename T> void TBaseLinearVector<T>::Add_bB(const T b, const TBaseLinearVector<T> &B)
{
  DEBUG_ASSERT(mem_dim == B.mem_dim);

  CblasHelper<T>::cblas_axpy(mem_dim, &b, B.mem, 1, mem, 1); // A = A + b * B
}

// ------------------ //
// A *= a             //
// ------------------ //
template <typename T> void TBaseLinearVector<T>::MultiplyBy(const T a)
{
  CblasHelper<T>::cblas_scal(mem_dim, &a, mem, 1);
}

// ------------------ //
// < A, B >           //
// ------------------ //
template <typename T> T TBaseLinearVector<T>::ScalarProduct(const TBaseLinearVector<T> &B) const
{
  DEBUG_ASSERT(mem_dim == B.mem_dim);

  T res;
  CblasHelper<T>::cblas_dot_sub(mem_dim, mem, 1, B.mem, 1, &res);
  return res;
}

#endif /* BASE_LINEAR_VECTOR_IMPL_H_ */
#ifndef MATRIX_IMPL_H_
#define MATRIX_IMPL_H_


// ================================================= //
// =================  Constructors  ================ //
// ================================================= //

template <class T> TMatrix<T>::TMatrix(const uint n)
  : TBaseLinearVector<T>(n * n)
  , nr(n), nc(n), is_square(true)
{
}

template <class T> TMatrix<T>::TMatrix(const uint nrows, const uint ncols)
  : TBaseLinearVector<T>(nrows * ncols)
  , nr(nrows), nc(ncols), is_square(nrows == ncols)
{
}

template <class T> TMatrix<T>::TMatrix(const uint n, T *memptr)
  : TBaseLinearVector<T>(n * n, memptr)
  , nr(n), nc(n), is_square(true)
{
}

template <class T> TMatrix<T>::TMatrix(const uint nrows, const uint ncols, T *memptr)
  : TBaseLinearVector<T>(nrows * ncols, memptr)
  , nr(nrows), nc(ncols), is_square(nrows == ncols)
{
}

template <class T> TMatrix<T>::TMatrix(const TMatrix<T> &other)
  : TBaseLinearVector<T>(other)
  , nr(other.nr), nc(other.nc), is_square(other.is_square)
{
}

template <class T> TMatrix<T>::~TMatrix()
{
}

// ================================================= //
// ============   Protected members    ============= //
// ================================================= //

template <class T> void TMatrix<T>::MatrixProduct(T* out_mem, const T* A_mem, const T* B_mem,
                                                  const uint A_nrows, const uint A_ncols, const uint B_ncols)
{
  for(uint a = 0; a < A_nrows; a++)
  for(uint b = 0; b < B_ncols; b++)
  {
    uint i1 = TMatrix<T>::ElementIndex(a, b, A_nrows, B_ncols);

    out_mem[i1] = 0;

    for(uint c = 0; c < A_ncols; c++)
    {
      uint i2 = TMatrix<T>::ElementIndex(a, c, A_nrows, A_ncols);
      uint i3 = TMatrix<T>::ElementIndex(c, b, A_ncols, B_ncols);

      out_mem[i1] += A_mem[i2] * B_mem[i3];
    }
  }
}

// ================================================= //
// ========  Memeber Operator overloading  ========= //
// ================================================= //

template <class T> T& TMatrix<T>::operator()(const uint row, const uint col) const
{
  DEBUG_ASSERT(row < nr && col < nc);

  uint idx = ElementIndex(row, col);

  return TBaseLinearVector<T>::mem[idx];
}

template <class T> TMatrix<T>& TMatrix<T>::operator=(const TMatrix<T> &rhs)
{
  DEBUG_ASSERT(nr == rhs.nr && nc == rhs.nc);

  TBaseLinearVector<T>::operator=(rhs);

  return *this;
}

// ================================================= //
// ======  Non-Memeber Operator overloading  ======= //
// ================================================= //

template <class T> TMatrix<T> operator*(const TMatrix<T>& lhs, const TMatrix<T>& rhs)
{
  DEBUG_ASSERT(lhs.nc == rhs.nr);

  const uint A_nr = lhs.nr;
  const uint A_nc = lhs.nc;
  const uint B_nc = rhs.nc;

  TMatrix<T> res(A_nr, B_nc);

  TMatrix<T>::MatrixProduct(res.DataPtr(), lhs.DataPtr(), rhs.DataPtr(), A_nr, A_nc, B_nc);

  return res;
}

template <class T> TBaseLinearVector<T> operator*(const TMatrix<T>& lhs, const TBaseLinearVector<T>& rhs)
{
  DEBUG_ASSERT(lhs.nc == rhs.Count());

  const uint A_nr = lhs.nr;
  const uint A_nc = lhs.nc;
  const uint B_nc = 1;

  TBaseLinearVector<T> res(A_nr);

  TMatrix<T>::MatrixProduct(res.DataPtr(), lhs.DataPtr(), rhs.DataPtr(), A_nr, A_nc, B_nc);

  return res;
}

template <class T> TLinearVector<T> operator*(const TMatrix<T>& lhs, const TLinearVector<T>& rhs)
{
  DEBUG_ASSERT(lhs.nc == rhs.Count());

  const uint A_nr = lhs.nr;
  const uint A_nc = lhs.nc;
  const uint B_nc = 1;

  TLinearVector<T> res(rhs.GetLattice(), A_nr);

  TMatrix<T>::MatrixProduct(res.DataPtr(), lhs.DataPtr(), rhs.DataPtr(), A_nr, A_nc, B_nc);

  return res;
}

template <class T> TColorVector<T> operator*(const TMatrix<T>& lhs, const TColorVector<T>& rhs)
{
  DEBUG_ASSERT(lhs.IsSquareMatrix() && lhs.nc == rhs.Count());

  const uint A_nr = lhs.nr;
  const uint A_nc = lhs.nc;
  const uint B_nc = 1;

  TColorVector<T> res(rhs.GetLattice());

  TMatrix<T>::MatrixProduct(res.DataPtr(), lhs.DataPtr(), rhs.DataPtr(), A_nr, A_nc, B_nc);

  return res;
}

template <class T> TSpinor<T> operator*(const TMatrix<T>& lhs, const TSpinor<T>& rhs)
{
  DEBUG_ASSERT(lhs.IsSquareMatrix() && lhs.nc == rhs.Count());

  const uint A_nr = lhs.nr;
  const uint A_nc = lhs.nc;
  const uint B_nc = 1;

  TSpinor<T> res(rhs.GetLattice());

  TMatrix<T>::MatrixProduct(res.DataPtr(), lhs.DataPtr(), rhs.DataPtr(), A_nr, A_nc, B_nc);

  return res;
}

template <class T> TSpaceVector<T> operator*(const TMatrix<T>& lhs, const TSpaceVector<T>& rhs)
{
  DEBUG_ASSERT(lhs.IsSquareMatrix() && lhs.nc == rhs.Count());

  const uint A_nr = lhs.nr;
  const uint A_nc = lhs.nc;
  const uint B_nc = 1;

  TSpaceVector<T> res(rhs.GetLattice());

  TMatrix<T>::MatrixProduct(res.DataPtr(), lhs.DataPtr(), rhs.DataPtr(), A_nr, A_nc, B_nc);

  return res;
}

template <class T> TSpinorField<T> operator*(const TMatrix<T>& lhs, const TSpinorField<T>& rhs)
{
  DEBUG_ASSERT(lhs.IsSquareMatrix() && lhs.nc == rhs.Count());

  const uint A_nr = lhs.nr;
  const uint A_nc = lhs.nc;
  const uint B_nc = 1;

  TSpinorField<T> res(rhs.GetLattice());

  TMatrix<T>::MatrixProduct(res.DataPtr(), lhs.DataPtr(), rhs.DataPtr(), A_nr, A_nc, B_nc);

  return res;
}

template <class T> TVectorField<T> operator*(const TMatrix<T>& lhs, const TVectorField<T>& rhs)
{
  DEBUG_ASSERT(lhs.IsSquareMatrix() && lhs.nc == rhs.Count());

  const uint A_nr = lhs.nr;
  const uint A_nc = lhs.nc;
  const uint B_nc = 1;

  TVectorField<T> res(rhs.GetLattice());

  TMatrix<T>::MatrixProduct(res.DataPtr(), lhs.DataPtr(), rhs.DataPtr(), A_nr, A_nc, B_nc);

  return res;
}

template <class T> TScalarField<T> operator*(const TMatrix<T>& lhs, const TScalarField<T>& rhs)
{
  DEBUG_ASSERT(lhs.IsSquareMatrix() && lhs.nc == rhs.Count());

  const uint A_nr = lhs.nr;
  const uint A_nc = lhs.nc;
  const uint B_nc = 1;

  TScalarField<T> res(rhs.GetLattice());

  TMatrix<T>::MatrixProduct(res.DataPtr(), lhs.DataPtr(), rhs.DataPtr(), A_nr, A_nc, B_nc);

  return res;
}

// ================================================= //
// ===============  Public members  ================ //
// ================================================= //

template <class T> uint TMatrix<T>::ElementIndex(const uint row, const uint col) const
{
  return (row * nc + col);
}

template <class T> uint TMatrix<T>::ElementIndex(const uint row, const uint col, const uint nrows, const uint ncols)
{
  return (row * ncols + col);
}

template <class T> uint TMatrix<T>::ElementIndex(const uint row, const uint col, const uint n)
{
  return (row * n + col);
}

template <class T> uint TMatrix<T>::Nrows() const
{
  return nr;
}

template <class T> uint TMatrix<T>::Ncols() const
{
  return nc;
}

template <class T> bool TMatrix<T>::IsSquareMatrix() const
{
  return is_square;
}

template <class T> void TMatrix<T>::Transpose()
{
  DEBUG_ASSERT(is_square);

  for(uint a = 0; a < nr; a++)
  for(uint b = 0; b < a; b++)
  {
    T elem = (*this)(a, b);
    (*this)(a, b) = (*this)(b, a);
    (*this)(b, a) = elem;
  }
}

template <> void TMatrix<t_complex>::ConjugateTranspose();
template <> void TMatrix<double>::ConjugateTranspose();

template <class T> void TMatrix<T>::ConjugateTranspose()
{
  ASSERTVERB(true, "Conjugate transpose is not implemented for the type T");
}

template <class T> TMatrix<T> TMatrix<T>::Transpose(TMatrix<T>& m)
{
  TMatrix<T> res(m);
  res.Transpose();
  return res;
}

template <class T> TMatrix<T> TMatrix<T>::ConjugateTranspose(TMatrix<T>& m)
{
  TMatrix<T> res(m);
  res.ConjugateTranspose();
  return res;
}

template <class T> T TMatrix<T>::Trace() const
{
  DEBUG_ASSERT(is_square);

  T res = 0;

  for(uint i = 0; i < nr; i++)
    res += (*this)(i, i);

  return res;
}

// ================================================= //
// ==========  Operator specializations  =========== //
// ================================================= //

#endif /* MATRIX_IMPL_H_ */
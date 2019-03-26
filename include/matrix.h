#ifndef MATRIX_H_
#define MATRIX_H_

#include <common.h>
#include <linear_vector.h>
#include <spinor_field.h>
#include <vector_field.h>
#include <scalar_field.h>

template <class T> class TMatrix;

template <class T> TMatrix<T>           operator*(const TMatrix<T>& lhs, const TMatrix<T>& rhs          );
template <class T> TBaseLinearVector<T> operator*(const TMatrix<T>& lhs, const TBaseLinearVector<T>& rhs);
template <class T> TLinearVector<T>     operator*(const TMatrix<T>& lhs, const TLinearVector<T>& rhs    );
template <class T> TColorVector<T>      operator*(const TMatrix<T>& lhs, const TColorVector<T>& rhs     );
template <class T> TSpinor<T>           operator*(const TMatrix<T>& lhs, const TSpinor<T>& rhs          );
template <class T> TSpaceVector<T>      operator*(const TMatrix<T>& lhs, const TSpaceVector<T>& rhs     );
template <class T> TSpinorField<T>      operator*(const TMatrix<T>& lhs, const TSpinorField<T>& rhs     );
template <class T> TVectorField<T>      operator*(const TMatrix<T>& lhs, const TVectorField<T>& rhs     );
template <class T> TScalarField<T>      operator*(const TMatrix<T>& lhs, const TScalarField<T>& rhs     );

template <class T> class TMatrix : public TBaseLinearVector<T>
{
protected:
  uint nr;
  uint nc;
  bool is_square;

  static void MatrixProduct(T* out_mem, const T* A_mem, const T* B_mem,
                            const uint A_nrows, const uint A_ncols, const uint B_ncols);

public:
  TMatrix(const uint n);
  TMatrix(const uint nrows, const uint ncols);
  TMatrix(const uint n, T *memptr);
  TMatrix(const uint nrows, const uint ncols, T *memptr);
  TMatrix(const TMatrix<T> &other);
  ~TMatrix();

  T& operator()(const uint row, const uint col) const;

  TMatrix<T>& operator=(const TMatrix<T> &rhs);

  friend TMatrix<T>           operator*<>(const TMatrix<T>& lhs, const TMatrix<T>& rhs          );
  friend TBaseLinearVector<T> operator*<>(const TMatrix<T>& lhs, const TBaseLinearVector<T>& rhs);
  friend TLinearVector<T>     operator*<>(const TMatrix<T>& lhs, const TLinearVector<T>& rhs    );
  friend TColorVector<T>      operator*<>(const TMatrix<T>& lhs, const TColorVector<T>& rhs     );
  friend TSpinor<T>           operator*<>(const TMatrix<T>& lhs, const TSpinor<T>& rhs          );
  friend TSpaceVector<T>      operator*<>(const TMatrix<T>& lhs, const TSpaceVector<T>& rhs     );
  friend TSpinorField<T>      operator*<>(const TMatrix<T>& lhs, const TSpinorField<T>& rhs     );
  friend TVectorField<T>      operator*<>(const TMatrix<T>& lhs, const TVectorField<T>& rhs     );
  friend TScalarField<T>      operator*<>(const TMatrix<T>& lhs, const TScalarField<T>& rhs     );

         uint ElementIndex(const uint row, const uint col) const;
  static uint ElementIndex(const uint row, const uint col, const uint nrows, const uint ncols);
  static uint ElementIndex(const uint row, const uint col, const uint n);

  uint Nrows() const;
  uint Ncols() const;
  bool IsSquareMatrix() const;

         void Transpose();
         void ConjugateTranspose();
  static TMatrix<T> Transpose(TMatrix<T>& m);
  static TMatrix<T> ConjugateTranspose(TMatrix<T>& m);

  T Trace() const;
};

#include <matrix_impl.h>

typedef TMatrix<t_complex> Matrix;

#endif /* MATRIX_H_ */ 

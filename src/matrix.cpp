#include <matrix.h>

// ================================================= //
// ==========  Operator specializations  =========== //
// ================================================= //

template <> void TMatrix<t_complex>::ConjugateTranspose()
{
  DEBUG_ASSERT(is_square);

  for(uint a = 0; a < nr; a++)
  for(uint b = 0; b < a; b++)
  {
    t_complex elem = (*this)(a, b);
    (*this)(a, b) = conj( (*this)(b, a) );
    (*this)(b, a) = conj(elem);
  }

  for(uint a = 0; a < nr; a++)
  {
    (*this)(a, a) = conj( (*this)(a, a) );
  }
}

template <> void TMatrix<double>::ConjugateTranspose()
{
  Transpose();
}

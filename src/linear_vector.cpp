#include <linear_vector.h>
#include <linalg.h>

// =================================================================== //
// ======  Non-Memeber Operator overloading (Specializations)  ======= //
// =================================================================== //

// ================================================= //
// ============  std::complex<double>  ============= //
// ================================================= //

template <> TBaseLinearVector<t_complex>& TBaseLinearVector<t_complex>::operator=(const TBaseLinearVector<t_complex>& rhs)
{
  DEBUG_ASSERT( (this != &rhs) && (rhs.mem_dim == mem_dim) );
  memcpy(mem, rhs.mem, sizeof(t_complex) * mem_dim);
  return *this;
}

template <> void TBaseLinearVector<t_complex>::Clear()
{
  memset(mem, 0, sizeof(t_complex) * mem_dim);
}

// ================================================= //
// =============  std::complex<float>  ============= //
// ================================================= //

template <> TBaseLinearVector<t_complex_float>& TBaseLinearVector<t_complex_float>::operator=(const TBaseLinearVector<t_complex_float>& rhs)
{
  DEBUG_ASSERT( (this != &rhs) && (rhs.mem_dim == mem_dim) );
  memcpy(mem, rhs.mem, sizeof(t_complex_float) * mem_dim);
  return *this;
}

template <> void TBaseLinearVector<t_complex_float>::Clear()
{
  memset(mem, 0, sizeof(t_complex_float) * mem_dim);
}

// ================================================= //
// ============        double          ============= //
// ================================================= //

template <> TBaseLinearVector<double>& TBaseLinearVector<double>::operator=(const TBaseLinearVector<double>& rhs)
{
  DEBUG_ASSERT( (this != &rhs) && (rhs.mem_dim == mem_dim) );
  memcpy(mem, rhs.mem, sizeof(double) * mem_dim);
  return *this;
}

template <> void TBaseLinearVector<double>::Clear()
{
  memset(mem, 0, sizeof(double) * mem_dim);
}

// ================================================= //
// =============        float          ============= //
// ================================================= //

template <> TBaseLinearVector<float>& TBaseLinearVector<float>::operator=(const TBaseLinearVector<float>& rhs)
{
  DEBUG_ASSERT( (this != &rhs) && (rhs.mem_dim == mem_dim) );
  memcpy(mem, rhs.mem, sizeof(float) * mem_dim);
  return *this;
}

template <> void TBaseLinearVector<float>::Clear()
{
  memset(mem, 0, sizeof(float) * mem_dim);
}

#ifndef FINITE_DIFF_H_
#define FINITE_DIFF_H_

#include "common.h"
#include "boost/rational.hpp"

using boost::rational;

template <typename T> class AuxVector
{
private:
  uint ndim;
  VECTOR<T> xs;

public:
  AuxVector(uint dim) : ndim(dim), xs(dim, 0) {}
  AuxVector(const AuxVector &other) : ndim(other.ndim), xs(other.xs) {}

  uint Dim() const  { return ndim; }
  VECTOR<T> Components() const { return xs; }

  T GetComponent(const uint i) const  { return xs[i]; }
  void SetComponent(const uint i, T val)  { xs[i] = val; }

  T& operator[](const uint i) { return xs[i]; }

  bool operator==(const AuxVector<T>& rhs) const
  {
    if (ndim != rhs.ndim) return false;
    bool res = true;
    for(uint i = 0; i < ndim; i++)  res = res && xs[i] == rhs[i];
    return res;
  }

  AuxVector<T>& operator=(const AuxVector<T> &rhs)  { ASSERT(ndim == rhs.ndim); for(uint i = 0; i < ndim; i++) xs[i] = rhs.xs[i]; }
  AuxVector<T>& operator+=(const AuxVector<T>& rhs) { ASSERT(ndim == rhs.ndim); for(uint i = 0; i < ndim; i++) xs[i] += rhs.xs[i]; }
  AuxVector<T>& operator-=(const AuxVector<T>& rhs) { ASSERT(ndim == rhs.ndim); for(uint i = 0; i < ndim; i++) xs[i] -= rhs.xs[i]; }  
  AuxVector<T>& operator*=(const T& rhs)            {                           for(uint i = 0; i < ndim; i++) xs[i] *= rhs; }
};

template <typename xT, typename coefT> class StencilPoint
{
public:
  AuxVector<xT> offset;
  rational<coefT> coef;

  StencilPoint(uint dim) : offset(dim), coef(0, 1) {}
  StencilPoint(uint dim, rational<coefT> c) : offset(dim), coef(c) {}
  StencilPoint(const AuxVector<xT> &off) : offset(off), coef(0, 1) {}
  StencilPoint(const AuxVector<xT> &off, rational<coefT> c) : offset(off), coef(c) {}
  StencilPoint(const StencilPoint<xT, coefT> &other) : offset(other.offset), coef(other.coef) {}
};

template <typename T> class FiniteDifference
{
private:
  uint dim;
  VECTOR<StencilPoint<int, T>> stencil;

  static T ipow(T x, uint p)
  {
    T i = 1;
    for (uint j = 1; j <= p; j++)  i *= x;
    return i;
  }

  static T factorial(T n)
  {
    T i = 1;
    for (T j = 1; j <= n; j++)  i *= j;
    return i;
  }

  static void SolveWithGaussianElimination(uint dim, VECTOR<VECTOR<rational<T>>> &system)
  {
    if (dim == 0) return;

    // Check if system has appropriate dimensions
    ASSERT(dim == system.size())
    for(uint i = 0; i < dim; i++)
    {
      ASSERT(dim + 1 == system[i].size());
    }

    // Forward propagation
    for(uint i = 0; i < dim; i++)
    {
      rational<T> c = system[i][i];
      if (c == 0) THROW_EXCEPTION(std::invalid_argument, "System has no unique solution")

      for(uint j = 0; j < dim; j++)
      {
        if (j == i) continue;

        rational<T> mult = system[j][i] / c;

        for(uint k = i; k < dim + 1; k++)
        {
          system[j][k] -= mult * system[i][k];
        }
      }
    }

    // Rescaling
    for(uint i = 0; i < dim; i++)
    {
      system[i][dim] *= 1 / system[i][i];
      system[i][i] = 1;
    }
  }

public:
  FiniteDifference(uint ndim) : dim(ndim) {}

  uint Dim() const { return dim; }

  VECTOR<StencilPoint<int, T>> GetStencil() const { return stencil; }

  void PrintStencil()
  {
    for(uint i = 0; i < stencil.size(); i++)
    {
      pStdLogs.Write("(");
      for(uint j = 0; j + 1 < dim; j++)
        pStdLogs.Write("% -d, ", stencil[i].offset[j]);
      pStdLogs.Write("% -d):\t", stencil[i].offset[dim - 1]);

      pStdLogs.Write("% -d/%d\n", stencil[i].coef.numerator(), stencil[i].coef.denominator());
    }
  }

  static SHARED_PTR<FiniteDifference<T>> MakeOneSidedDiff(uint order, int n_points)
  {
    ASSERT(order < std::abs(n_points));

    const uint dim = 1;
    SHARED_PTR<FiniteDifference<T>> res_ptr = MAKE_SHARED<FiniteDifference<T>>(dim);

    const uint npts = std::abs(n_points);
    const int direction = (n_points < 0) ? -1 : 1;

    VECTOR<VECTOR<rational<T>>> system(npts);
    for(uint i = 0; i < npts; i++)
    {
      system[i].resize(npts + 1);

      for(uint j = 0; j < npts; j++)
      {
        rational<T> a_ij(ipow((T)direction * j, i), factorial(i));
        system[i][j] = a_ij;
      }

      system[i][npts] = (i == order) ? 1 : 0;
    }

    SolveWithGaussianElimination(npts, system);

    res_ptr->stencil.reserve(npts);
    for(uint i = 0; i < npts; i++)
    {
      AuxVector<int> offset(dim);
      offset[0] = direction * (int)i;

      StencilPoint<int, T> stencil_pt(offset, system[i][npts]);

      res_ptr->stencil.push_back(stencil_pt);
      //pStdLogs.Write("%d : %d\n", system[i][npts].numerator(), system[i][npts].denominator());
    }

    return res_ptr;
  }
};

#endif /* FINITE_DIFF_H_ */

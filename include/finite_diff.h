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
    for(uint i = 0; i < ndim; i++)  res = res && xs[i] == rhs.GetComponent(i);
    return res;
  }

  AuxVector<T>& operator=(const AuxVector<T> &rhs)  { ASSERT(ndim == rhs.ndim); for(uint i = 0; i < ndim; i++) xs[i] = rhs.GetComponent(i); return *this; }
  AuxVector<T>& operator+=(const AuxVector<T>& rhs) { ASSERT(ndim == rhs.ndim); for(uint i = 0; i < ndim; i++) xs[i] += rhs.GetComponent(i); return *this; }
  AuxVector<T>& operator-=(const AuxVector<T>& rhs) { ASSERT(ndim == rhs.ndim); for(uint i = 0; i < ndim; i++) xs[i] -= rhs.GetComponent(i); return *this; }  
  AuxVector<T>& operator*=(const T& rhs)            {                           for(uint i = 0; i < ndim; i++) xs[i] *= rhs; return *this; }
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
  uint NumPoints() const { return stencil.size(); }

  VECTOR<StencilPoint<int, T>> GetStencil() const { return stencil; }

  void PrintStencil()
  {
    for(uint i = 0; i < stencil.size(); i++)
    {
      pStdLogs.Write("(");
      for(uint j = 0; j + 1 < dim; j++)
        pStdLogs.Write("% -d, ", stencil[i].offset[j]);
      pStdLogs.Write("% -d):\t", stencil[i].offset[dim - 1]);

      if (stencil[i].coef.denominator() != 1)
        pStdLogs.Write("% -d/%d\n", stencil[i].coef.numerator(), stencil[i].coef.denominator());
      else
        pStdLogs.Write("% -d\n", stencil[i].coef.numerator());
    }
  }

  static SHARED_PTR<FiniteDifference<T>> ComposeOperators(
    const FiniteDifference<T> &A, const FiniteDifference<T> &B)
  {
    ASSERT(A.dim == B.dim);

    const uint dim = A.dim;
    SHARED_PTR<FiniteDifference<T>> res_ptr = MAKE_SHARED<FiniteDifference<T>>(dim);
    FiniteDifference<T> &res = *res_ptr;

    const uint nptA = A.NumPoints();
    const uint nptB = B.NumPoints();

    AuxVector<int> aux(dim);
    rational<T> c;

    for(uint i = 0 ; i < nptA; i++)
    for(uint j = 0 ; j < nptB; j++)
    {
      const StencilPoint<int, T> &ptA = A.stencil[i];
      const StencilPoint<int, T> &ptB = B.stencil[j];

      aux = ptA.offset;
      aux += ptB.offset;
      c = ptA.coef * ptB.coef;

      StencilPoint<int, T> *_ptr = nullptr;
      for(uint k = 0; k < res.NumPoints(); k++)
      {
        if (res.stencil[k].offset == aux)
        {
          _ptr = &res.stencil[k];
          break;
        }
      }

      if (_ptr != nullptr)
      {
        _ptr->coef += c;
      }
      else
      {
        StencilPoint<int, T> new_pt(aux, c);
        res.stencil.push_back(new_pt);
      }
    }

    return res_ptr;
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
    }

    return res_ptr;
  }

  static SHARED_PTR<FiniteDifference<T>> MakeDiffMarc()
  {
    const uint dim = 1;
    const uint npts = 5;
    SHARED_PTR<FiniteDifference<T>> res_ptr = MAKE_SHARED<FiniteDifference<T>>(dim);

    AuxVector<int> aux(dim);

    res_ptr->stencil.reserve(npts);

    aux[0] = -2;
    res_ptr->stencil.emplace_back(aux, rational<T>(-1,12));
    aux[0] = -1;
    res_ptr->stencil.emplace_back(aux, rational<T>(16,12));
    aux[0] = 0;
    res_ptr->stencil.emplace_back(aux, rational<T>(-30,12));
    aux[0] = 1;
    res_ptr->stencil.emplace_back(aux, rational<T>(16,12));
    aux[0] = 2;
    res_ptr->stencil.emplace_back(aux, rational<T>(-1,12));

    return res_ptr;
  }

  static SHARED_PTR<FiniteDifference<T>> MakeLaplacian(
    uint dim, const FiniteDifference<T> &diff1d)
  {
    ASSERT(dim > 0);
    ASSERT(diff1d.Dim() == 1);

    auto res_ptr = MAKE_SHARED<FiniteDifference<T>>(dim);
    auto &res = *res_ptr;

    auto &stencil1d = diff1d.stencil;
    uint npt1d = stencil1d.size();

    AuxVector<int> aux(dim);

    for(uint i = 0; i < dim * npt1d; i++)
    {
      const uint d = (uint)floor((double)i / npt1d);
      const uint sidx = i % npt1d;

      for(uint j = 0; j < dim; j++)
        aux[j] = (j == d) ? stencil1d[sidx].offset.GetComponent(0) : 0;

      bool exists = false;
      for(uint j = 0; j < res.stencil.size(); j++)
      {
        if (aux == res.stencil[j].offset)
        {
          res.stencil[j].coef += stencil1d[sidx].coef;
          exists = true;
          break;
        }
      }

      if (!exists)
        res.stencil.emplace_back(aux, stencil1d[sidx].coef);
    }

    return res_ptr;
  }

  static t_complex EvalDerivative(const VECTOR<uint> &order, const VECTOR<double> &p, const FiniteDifference<T> &A)
  {
    ASSERT(p.size() == A.Dim());
    ASSERT(order.size() == A.Dim());

    t_complex res = 0;

    uint dim = A.Dim();
    auto &stencil = A.stencil;

    for(uint i = 0; i < stencil.size(); i++)
    {
      auto &offset = stencil[i].offset;
      t_complex coef = boost::rational_cast<double>(stencil[i].coef);

      double phase = 0;
      for(uint j = 0; j < dim; j++)
      {
        phase += p[j] * offset.GetComponent(j);
        
        if (order[j] > 0)
          coef *= pow(t_complex(0, 1) * (double)offset.GetComponent(j), order[j]);
      }

      res += coef * exp(t_complex(0, 1) * phase);
    }

    return res;
  }

  static t_complex Eval(const VECTOR<double> &p, const FiniteDifference<T> &A)
  {
    ASSERT(p.size() == A.Dim());

    t_complex res = 0;

    uint dim = A.Dim();
    auto &stencil = A.stencil;

    for(uint i = 0; i < stencil.size(); i++)
    {
      auto &offset = stencil[i].offset;
      auto coef = stencil[i].coef;

      double phase = 0;
      for(uint j = 0; j < dim; j++)
        phase += p[j] * offset.GetComponent(j);

      res += boost::rational_cast<double>(coef) * exp(t_complex(0, 1) * phase);
    }

    return res;
  }

  t_complex EvalDerivative(const VECTOR<uint> &order, const VECTOR<double> &p)
  {
    return FiniteDifference<T>::EvalDerivative(order, p, *this);
  }

  t_complex Eval(const VECTOR<double> &p) const
  {
    return FiniteDifference<T>::Eval(p, *this);
  }
};

#endif /* FINITE_DIFF_H_ */

#ifndef FINITE_DIFF_H_
#define FINITE_DIFF_H_

#include "common.h"
#include "boost/rational.hpp"

using boost::rational;

typedef VECTOR<int> tOffset;
typedef std::tuple<tOffset, rational<int>> tStencilPt;

class FiniteDifference
{
private:
  uint dim;
  VECTOR<tStencilPt> stencil;

  static uint ipow(uint x, uint p)
  {
    uint i = 1;
    for (uint j = 1; j <= p; j++)  i *= x;
    return i;
  }

  static uint factorial(uint n)
  {
    uint i = 1;
    for (uint j = 1; j <= n; j++)  i *= j;
    return i;
  }

  static void SolveWithGaussianElimination(uint dim, VECTOR<VECTOR<rational<int>>> &system)
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
      rational<int> c = system[i][i];
      if (c == 0) THROW_EXCEPTION(std::invalid_argument, "System has no unique solution")

      for(uint j = 0; j < dim; j++)
      {
        if (j == i) continue;

        rational<int> mult = system[j][i] / c;

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
  FiniteDifference(uint ndim) {}

  uint Dim() { return dim; }

  VECTOR<tStencilPt> GetStencil() { return stencil; }

  static SHARED_PTR<FiniteDifference> MakeOneSidedDiff(uint order, int n_points)
  {
    ASSERT(order < std::abs(n_points));

    const uint dim = 1;
    SHARED_PTR<FiniteDifference> res_ptr = MAKE_SHARED<FiniteDifference>(dim);

    const uint npts = std::abs(n_points);
    const int direction = (n_points < 0) ? -1 : 1;

    VECTOR<VECTOR<rational<int>>> system(npts);
    for(uint i = 0; i < npts; i++)
    {
      system[i].resize(npts + 1);

      for(uint j = 0; j < npts; j++)
      {
        rational<int> a_ij(ipow(direction * j, i), factorial(i));
        system[i][j] = a_ij;
      }

      system[i][npts] = (i == order) ? 1 : 0;
    }

    SolveWithGaussianElimination(npts, system);

    res_ptr->stencil.reserve(npts);
    for(uint i = 0; i < npts; i++)
    {
      tOffset offset = { direction * (int)i };
      tStencilPt stencil_pt = std::make_tuple(offset, system[i][npts]);

      res_ptr->stencil.push_back(stencil_pt);
      //pStdLogs.Write("%d : %d\n", system[i][npts].numerator(), system[i][npts].denominator());
    }

    return res_ptr;
  }
};

#endif /* FINITE_DIFF_H_ */

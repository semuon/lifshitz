#include "linalg.h"
#include <lapacke.h>
#include "profiler.h"

// Auxiliary struct to keep eigenvalues and eigenvectors in pairs
struct eigen { 
  int i;
  t_complex val;
  Linalg::SortingMode mode;

  eigen(int i, t_complex val, Linalg::SortingMode mode) :
    i(i), val(val), mode(mode) {}

  bool operator<(eigen const &other) const
  {
    switch(mode)
    {
      case Linalg::AscendingAbs:
        return abs(val) < abs(other.val);
      case Linalg::DescendingAbs:
        return abs(val) > abs(other.val);
      case Linalg::AscendingReal:
        return real(val) < real(other.val);
      case Linalg::DescendingReal:
        return real(val) > real(other.val);
      case Linalg::AscendingImag:
        return imag(val) < imag(other.val);
      case Linalg::DescendingImag:
        return imag(val) > imag(other.val);
      default:
        THROW_EXCEPTION(AssertFailure, "Unknown sorting mode");
    }
  }
};
 
void Linalg::SortEigensystem(VECTOR<int> &sorted_idx, const VECTOR<t_complex> &evals, const SortingMode mode)
{
  int n = evals.size();
  VECTOR<eigen> pairs;

  pairs.reserve(n);
  for(int i = 0; i < n; i++)
    pairs.emplace_back(i, evals[i], mode);

  std::sort(pairs.begin(), pairs.end());

  sorted_idx.resize(n);
  for(int i = 0; i < n; i++)
    sorted_idx[i] = pairs[i].i;
}

void Linalg::TestOrthogonality(VECTOR<SpinorField> &evecs, const VECTOR<int> idx_to_test, const double tol)
{
 ASSERT(idx_to_test.size() <= evecs.size());
  
 bool should_check = true;
 bool first_check = true; 
 
 while(should_check)
  {
    bool success = true;
  
   for(uint iv = 0; iv < idx_to_test.size(); iv++)
   {
     int i = idx_to_test.at(iv);
      VECTOR<int> vec_for_orth;

      for(uint jv = iv + 1; jv < idx_to_test.size(); jv++)
      {
        int j = idx_to_test.at(jv);

        t_complex sp = evecs.at(i) * evecs.at(j);
        double spnorm = abs(sp);
        if (spnorm > tol)
        {
          success = false;

          vec_for_orth.push_back(j);
        }
      }

      // If there are some non-orthogonal vectors, do Gramm-Schmidt
      if (vec_for_orth.size() > 0)
      {
        vec_for_orth.push_back(i);

        Linalg::GrammSchmidt(vec_for_orth, evecs);
      }
    }

    if (!first_check && !success)
    {
      THROW_EXCEPTION(EigenTestFailure,
            "Eigenvectors are not orthogonal (FAILED TO FIX)");
    }

    should_check = first_check && !success;
    first_check = false;
  };
}

void Linalg::TestEigensystem(MatrixByVectorFunc F, void *args,
                             const VECTOR<t_complex> &evals, VECTOR<SpinorField> &evecs,
                             const VECTOR<int> idx_to_test, const double tol)
{
  ASSERT(evals.size() == evecs.size());
  ASSERT(idx_to_test.size() <= evecs.size());

  if (evecs.size() == 0)
    return;

  SpinorField tmp(evecs.at(0));

  for(uint iv = 0; iv < idx_to_test.size(); iv++)
  {
    int i = idx_to_test.at(iv);

    F(tmp, evecs.at(i), args);
    tmp.Add_bB(-1.0 * evals.at(i), evecs.at(i));

    double resid = tmp.Norm();
    if (resid > tol)
    {
       THROW_EXCEPTION(EigenTestFailure,
         "Eigensystem " + TO_STRING(i) + " is not accurate: |H * v - lambda * v| = " + TO_STRING(resid) + " > tol = " + TO_STRING(tol) + "");
    }

    resid = fabs(evecs.at(i).Norm() - 1);
    if (resid > tol)
    {
      evecs.at(i).Normalize();
      // THROW_EXCEPTION(EigenTestFailure,
      //     "Eigenvector norm is not 1: |v_" + TO_STRING(i) + "| - 1 = " + TO_STRING(resid) + " > tol = " + TO_STRING(tol) + "");
    }
  }

}

void Linalg::GrammSchmidt(const VECTOR<int> &vec_idx, VECTOR<SpinorField> &vecs)
{
  ASSERT(vec_idx.size() <= vecs.size());

  if (vec_idx.size() == 0)
    return;

  const int n = vec_idx.size();

  for(int i = 0; i < n; i++)
  {
    int idx_i = vec_idx.at(i);

    SpinorField &v_i = vecs.at(idx_i);

    for(int j = 0; j < i; j++)
    {
      int idx_j = vec_idx.at(j);

      SpinorField &u_j = vecs.at(idx_j);

      u_j.Add_bB(-1.0 * (v_i * u_j) / v_i.NormSqr(), v_i);
    }
  }

  for(int i = 0; i < n; i++)
  {
    int idx_i = vec_idx.at(i);
    vecs.at(idx_i).Normalize();
  }
}

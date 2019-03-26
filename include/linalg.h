#ifndef _LINALG_H
#define _LINALG_H

#include <common.h>
#include <linear_vector.h>
#include <spinor_field.h>

// Exception of Lapack routines
class LapackFailure : public std::runtime_error
{
public:
  LapackFailure(const std::string &msg) : std::runtime_error(msg) {}
  LapackFailure(const char *msg) : std::runtime_error(msg) {}
};

// Exception of Conjugate Gradient iterative solver
class CGFailure : public std::runtime_error
{
public:
  CGFailure(const std::string &msg) : std::runtime_error(msg) {}
  CGFailure(const char *msg) : std::runtime_error(msg) {}
};

// Exception of Arnoldi
class ArnoldiFailure : public std::runtime_error
{
public:
  ArnoldiFailure(const std::string &msg) : std::runtime_error(msg) {}
  ArnoldiFailure(const char *msg) : std::runtime_error(msg) {}
};

// Exception of eigensystem test
class EigenTestFailure : public std::runtime_error
{
public:
  EigenTestFailure(const std::string &msg) : std::runtime_error(msg) {}
  EigenTestFailure(const char *msg) : std::runtime_error(msg) {}
};

class Linalg
{
private:
  Linalg() {}

public:
  // Callback function for matrix by vector operations
  typedef void (*MatrixByVectorFunc)(SpinorField &out, const SpinorField &in, void *args);

  // Interface to LAPACK functions
  static void LapackInvert(t_complex *A, int n);

  static void LapackLinearSolve(t_complex *x, t_complex *A, int n);
  static void LapackLinearSolve(double *x, double *A, int n);
  static void LapackLinearSolveSymmetric(t_complex *x, t_complex *A, int n);

  static void LapackEigensystem(t_complex *A, t_complex *evals, t_complex *evecs, int n);
  static void LapackHermitianEigensystem(t_complex *A, double *evals, int n);

  static void LapackCholeskyDecomposition(t_complex *A, int n);

  // Iterative algorithms

  // This is algorithm from Jegerlehner, B. , http://arxiv.org/abs/hep-lat/9612014
  // F|solution > = |source >, |solution > = F^-1 |source>
  static void MultishiftCG(MatrixByVectorFunc F, void *args, const Lattice &lat, const SpinorField &source,
                           int nsys, VECTOR<SpinorField> &solution, VECTOR<double> sigma,
                           double tol, int imax, double &max_err, int &num_iter);

/**************************************************************************************/
/***************** Front-end interface to ARPACK **************************************/
/**************************************************************************************/

  typedef enum enumArnoldiMode
  {
    SmallestAbs  = 0,
    LargestAbs   = 1,
    SmallestReal = 2,
    LargestReal  = 3,
    SmallestImag = 4,
    LargestImag  = 5
  } ArnoldiMode;

/*! Some internal parameters which control arnoldi performance and which potentially require tuning. Refer to ARPACK documentation for their meaning */
  static const int c_default_arnoldi_ncv_factor_nom   = 12;
  static const int c_default_arnoldi_ncv_factor_denom = 1;
//!< See #arnoldi_ncv_factor_nom

/*!  \brief Front-end interface to ARPACK
 *
 * \param[in] F          is the function which implements the multiplication of the vector by operator
 * \param[in] n          is the size of the linear space
 * \param[out] **evals   is the pointer to the 1D array (size nev) of t_complex where the eigenvalues 
 *                       will be stored
 * \param[out] **evecs   is the pointer to the 1D array (size nev x n) of t_complex where the eigenvectors
 *                        will be stored
 * \param[in] nev        is the number of eigenstates to be calculated
 * \param[in] tol        is the required precision
 * \param[in] max_iter   is the maximal allowed number of Arnoldi iterations
 * \param[in] mode       defines which eigenvalues will be found:
 */
  static void Arnoldi(MatrixByVectorFunc F, void *args, const Lattice &lat, VECTOR<t_complex> &evals,
                      int nev, double tol, int max_iter, ArnoldiMode mode);

  static void Arnoldi(MatrixByVectorFunc F, void *args, const Lattice &lat, VECTOR<t_complex> &evals,
                      bool compute_evecs, VECTOR<SpinorField> &evecs, int nev, double tol, int max_iter,
                      ArnoldiMode mode, int arnoldi_ncv_factor_nom = c_default_arnoldi_ncv_factor_nom,
                      int arnoldi_ncv_factor_denom = c_default_arnoldi_ncv_factor_denom);

  typedef enum enumSortingMode
  {
    AscendingAbs    = 0,
    DescendingAbs   = 1,
    AscendingReal   = 2,
    DescendingReal  = 3,
    AscendingImag   = 4,
    DescendingImag  = 5
  } SortingMode;

  static void SortEigensystem(VECTOR<int> &sorted_idx, const VECTOR<t_complex> &evals, const SortingMode mode);
  
  static void TestOrthogonality(VECTOR<SpinorField> &evecs, const VECTOR<int> idx_to_test, const double tol);

  static void TestEigensystem(MatrixByVectorFunc F, void *args,
                              const VECTOR<t_complex> &evals, VECTOR<SpinorField> &evecs,
                              const VECTOR<int> idx_to_test, const double tol);

  static void GrammSchmidt(const VECTOR<int> &vec_idx, VECTOR<SpinorField> &vecs);
};

#endif

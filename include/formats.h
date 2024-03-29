#ifndef FORMATS_H_
#define FORMATS_H_

#include "common.h"
#include "vector_field.h"
#include "scalar_field.h"
#include "matrix.h"
#include "lattice.h"

class Formats
{
private:
  Formats() {}

public:
  template <class T> static void DumpBinary(FILE *f, const TBaseLinearVector<T> &field)
  {
    SAFE_FWRITE(field.DataPtr(), sizeof(T), field.Count(), f);
  }

  template <class T> static void DumpBinary(FILE *f, const VECTOR<T> &field)
  {
    SAFE_FWRITE(field.data(), sizeof(T), field.size(), f);
  }

  static void PrintMatrix(FILE *f, const Matrix &m);

  static void PrintVectorField(FILE *f, const double t, const VectorField &vec, const Lattice &lat);
  static void PrintVector(FILE *f, const double t, const TBaseLinearVector<double> &vec);
  static void PrintVector(FILE *f, const double t, const TBaseLinearVector<t_complex> &vec);

  static void PrintScalar(FILE *f, const double t, const double val);
  static void PrintScalarField(FILE *f, const double t, const ScalarField  &val, const Lattice &lat);
  static void PrintScalar(FILE *f, const double t, const t_complex val);

  static void ReadVectorField( const std::string name, const double t, VectorField &vec, const double rescale, const Lattice &lat);

};

#endif

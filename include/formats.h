#ifndef FORMATS_H_
#define FORMATS_H_

#include "common.h"
#include "vector_field.h"
#include "scalar_field.h"
#include "lattice.h"

class Formats
{
private:
  Formats() {}

public:
  template <class T> static void DumpBinary(FILE *f, const double t, const TBaseLinearVector<T> &field)
  {
    SAFE_FWRITE(&t, sizeof(double), 1, f);
    SAFE_FWRITE(field.DataPtr(), sizeof(T), field.Count(), f);
  }

  static void PrintVectorField(FILE *f, const double t, const VectorField &vec, const Lattice &lat);
  static void PrintVector(FILE *f, const double t, const TBaseLinearVector<double> &vec);
  static void PrintVector(FILE *f, const double t, const TBaseLinearVector<t_complex> &vec);

  static void PrintScalar(FILE *f, const double t, const double val);
  static void PrintScalarField(FILE *f, const double t, const ScalarField  &val, const Lattice &lat);
  static void PrintScalar(FILE *f, const double t, const t_complex val);

  static void ReadVectorField( const std::string name, const double t, VectorField &vec, const double rescale, const Lattice &lat);

  
  static void PrintOccupationNumber(FILE *f, const double t, const VECTOR< VECTOR< double > > &occ_num, const Lattice &lat, const Lattice &klat);
  static void PrintTotOccupationNumber(FILE *f, const double t, const  VECTOR< double >  &tot_occ_num,  const Lattice &klat); 

};

#endif

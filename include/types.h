#ifndef TYPES_H_
#define TYPES_H_

#include <complex>

// Complex number of double precision
typedef std::complex<double> t_complex;
typedef std::complex<float> t_complex_float;

#define I t_complex(0.0, 1.0)

#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>

typedef unsigned int uint;
typedef unsigned long ulong;

#endif /* TYPES_H_ */

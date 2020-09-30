#ifndef PHI4_H_
#define PHI4_H_

#include <common.h>
#include "scalar_field_n.h"

typedef struct ScalarModelParams_struct
{
  double m2;
  double invM2;
  double Z;
  double lambdaN;
  double N;
  double kappa;
} tScalarModelParams;

typedef struct LatticeScalarModelParams_struct
{
  double k1;
  double k2;
  double lambda;
  double kappa;
} tLatticeScalarModelParams;

class ScalarModel
{
private:
  ScalarModel() {}

public:
  static void ConvertCouplings(const tLatticeScalarModelParams &lattice_params, const int ndim, tScalarModelParams &phys_params);

  static double Action(const tScalarModelParams &params, const RealScalarFieldN &phi);
  static void HMCforce(const tScalarModelParams &params, const RealScalarFieldN &phi, RealScalarFieldN &force);

  static void CorrelationFunction(const RealScalarFieldN &phi, const bool vol_avg, const uint mu, VECTOR<double> &corr);
};

#endif

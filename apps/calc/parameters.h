#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <common.h>

typedef enum EnumStartConfigurationType
{
  START_CONFIGURATION_ZERO = 0,
  START_CONFIGURATION_RANDOM = 1,
  START_CONFIGURATION_LOAD = 2
}tStartConfigurationType;

typedef enum EnumIntegratorType
{
  INTEGRATOR_LEAPFROG = 0,
  INTEGRATOR_OMELYAN = 1
}tIntegratorType;

extern double pLatK1;
extern double pLatK2;
extern double pLatLambda;
extern double pLatKappa;

extern bool pIsLatticeParamsSet;

extern bool pIsComputeCorr;
extern bool pIsFullCorr;
extern bool pIsVolAvgCorr;

extern std::string pFnameConfs;

extern int pNumSkipFirst;
extern int pNumSkipLast;
extern int pConfStep;

#endif

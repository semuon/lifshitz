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

extern int pN;
extern double pKappa;

extern double pLatK1;
extern double pLatK2;
extern double pLatLambda;
extern double pLatKappa;

extern bool pIsLatticeParamsSet;

extern tStartConfigurationType pStartType;
extern tIntegratorType pIntegratorType;
extern std::string pFnameStartConf;

extern int pHmcNumSteps;
extern double pHmcDt;
extern int pHmcNumConf;
extern int pHmcNumConfStep;

extern bool pAutoTune;
extern double pAutoTuneK;

#endif

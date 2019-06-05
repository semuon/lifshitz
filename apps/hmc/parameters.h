#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <common.h>

typedef enum EnumStartConfigurationType
{
  START_CONFIGURATION_ZERO = 0,
  START_CONFIGURATION_RANDOM = 1,
  START_CONFIGURATION_LOAD = 2
}tStartConfigurationType;

extern int pN;
extern double pKappa;

extern tStartConfigurationType pStartType;
extern std::string pFnameStartConf;

extern int pHmcNumSteps;
extern double pHmcDt;
extern int pHmcNumConf;
extern int pHmcNumConfStep;

#endif

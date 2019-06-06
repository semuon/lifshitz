#include "parameters.h"

int pN = 1;
double pKappa = 0;

double pLatK1 = 1.0;
double pLatK2 = 0.0;
double pLatLambda = 1.0;
double pLatKappa = 0.0;

bool pIsLatticeParamsSet = false;

tIntegratorType pIntegratorType = INTEGRATOR_LEAPFROG;

tStartConfigurationType pStartType = START_CONFIGURATION_RANDOM;
std::string pFnameStartConf = "init.conf";

int pHmcNumSteps = 20;
double pHmcDt = 0.05;
int pHmcNumConf = 200;
int pHmcNumConfStep = 1;

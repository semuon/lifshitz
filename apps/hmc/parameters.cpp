#include "parameters.h"

int pN = 1;
double pKappa = 0;

tStartConfigurationType pStartType = START_CONFIGURATION_RANDOM;
std::string pFnameStartConf = "init.conf";

int pHmcNumSteps = 20;
double pHmcDt = 0.05;
int pHmcNumConf = 200;
int pHmcNumConfStep = 1;

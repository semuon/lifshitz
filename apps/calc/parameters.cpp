#include "parameters.h"

int pN = 1;
double pKappa = 0;

double pLatK1 = 1.0;
double pLatK2 = 0.0;
double pLatLambda = 1.0;
double pLatKappa = 0.0;

bool pIsLatticeParamsSet = false;

bool pIsComputeCorr = false;
bool pIsVolAvgCorr = false;

std::string pFnameConfs = "confs.bin";

int pNumSkipFirst = 0;
int pNumSkipLast = 0;
int pConfStep = 1;

#include "common_parameters.h"

uint pDim = 3;
uint pNumDiagDims = 0;
VECTOR<uint> pL(pDim, 1);

double pLambdaN = 1.0;
double pInvM2 = 0.0;
double pm2 = 0.0;
double pZ = 1.0;

double pExtH0 = 0;
double pExtK0 = 0.2 * M_PI;
double pExtSigma0 = 0.1;

int pNStencilPts = c_Winstel_stencil;

std::string pLogsDirPath;
std::string pDataDirPath;
std::string pPrefix;
bool pOverwriteLogs = true;

bool pSetSeed = false;
int pSeed = 0;

int pNthreads = 1;

DataDir pPrivateMPILogsDir;
DataDir pPrivateMPIDataDir;
DataDir pLogsDir;
DataDir pDataDir;

#ifdef WITH_MPI_
int pMPIRank = 0;
int pMPISize = 0;
std::string pMPIProcessorName;
#endif

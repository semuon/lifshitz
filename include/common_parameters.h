#ifndef COMMON_PARAMETERS_H_
#define COMMON_PARAMETERS_H_

#include "common.h"

const double c_flt_epsilon = __FLT_EPSILON__;

extern uint pDim;
extern uint pNumDiagDims;
extern VECTOR<uint> pL;

extern double pLambda;
extern double pInvM2;
extern double pm2;
extern double pZ;

extern std::string pLogsDirPath;
extern std::string pDataDirPath;
extern std::string pPrefix;
extern bool pOverwriteLogs;

extern bool pSetSeed;
extern int pSeed;

extern int pNthreads;

extern DataDir pPrivateMPILogsDir;
extern DataDir pPrivateMPIDataDir;
extern DataDir pLogsDir;
extern DataDir pDataDir;

#ifdef WITH_MPI_
const int pConstMPIMaster = 0;

extern int pMPIRank;
extern int pMPISize;
extern std::string pMPIProcessorName;
#endif

#endif /* COMMON_PARAMETERS_H_ */

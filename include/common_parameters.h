#ifndef COMMON_PARAMETERS_H_
#define COMMON_PARAMETERS_H_

#include "common.h"

extern uint pDim;
extern uint pNumDiagDims;
extern VECTOR<uint> pL;

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

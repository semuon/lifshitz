#ifndef MPI_MODULE_H_
#define MPI_MODULE_H_

#include <common.h>

class ModuleMPI
{
private:
  ModuleMPI() {}

public:
  static bool IsMasterNode();

  static int MasterRank();
  static int Rank();
  static int Size();

  static void ReduceRealVector(void *vec, void *buf, const uint count);

  static void AssignModes(const uint xmodes, const uint kmodes,
                          uint &knum, VECTOR<uint> &k, VECTOR<VECTOR<uint> >&n);
  static void AssignModes(const int mpi_size, const int mpi_rank, const uint xmodes, const uint kmodes,
                          uint &knum, VECTOR<uint> &k, VECTOR<VECTOR<uint> >&n);
};

#endif

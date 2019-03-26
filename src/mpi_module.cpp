#include <mpi_module.h>

bool ModuleMPI::IsMasterNode()
{
  bool is_master_node = true;

#ifdef WITH_MPI_
  is_master_node = (pMPIRank == pConstMPIMaster);
#endif

  return is_master_node;
}

int ModuleMPI::MasterRank()
{
  int master_rank = 0;

#ifdef WITH_MPI_
  master_rank = pConstMPIMaster;
#endif

  return master_rank;
}

int ModuleMPI::Rank()
{
  int rank = 0;

#ifdef WITH_MPI_
  rank = pMPIRank;
#endif

  return rank;
}

int ModuleMPI::Size()
{
  int size = 1;

#ifdef WITH_MPI_
  size = pMPISize;
#endif

  return size;
}

void ModuleMPI::ReduceRealVector(void *vec, void *buf, const uint count)
{
#ifdef WITH_MPI_
  SAFE_MPI_CALL(MPI_Reduce(vec, buf, count, MPI_DOUBLE, MPI_SUM,
                           pConstMPIMaster, MPI_COMM_WORLD
                          ));

  if (ModuleMPI::IsMasterNode())
    std::memcpy(vec, buf, count * sizeof(double));
#endif
}

void ModuleMPI::AssignModes(const uint xmodes, const uint kmodes,
                            uint &knum, VECTOR<uint> &k, VECTOR<VECTOR<uint> >&n)
{
  int mpi_size = 1;
  int mpi_rank = 0;

#ifdef WITH_MPI_
  mpi_size = pMPISize;
  mpi_rank = pMPIRank;
#endif

  AssignModes(mpi_size, mpi_rank, xmodes, kmodes, knum, k, n);
}

void ModuleMPI::AssignModes(const int mpi_size, const int mpi_rank, const uint xmodes, const uint kmodes,
                            uint &knum, VECTOR<uint> &k, VECTOR<VECTOR<uint> >&n)
{
  n.clear();
  k.clear();

  // Number of fermionic modes \psi_n(k) is always Vol * Nd * Nc
  uint nmodes_tot = xmodes * kmodes;
  uint nmodes_per_node = MAX(nmodes_tot / mpi_size, 1);

  if(nmodes_per_node * mpi_size < nmodes_tot)
    nmodes_per_node++;

  uint nmodes_start = mpi_rank * nmodes_per_node;        // Range of modes specific for a given node
  uint nmodes_end = nmodes_start + nmodes_per_node;

  // Check that we don't run out the range on the last node
  nmodes_end = (nmodes_end > nmodes_tot) ? nmodes_tot : nmodes_end;

  // Check if number of MPI nodes is large than number of modes,
  // Then for some nodes no computations should be performed
  nmodes_start = (nmodes_start < nmodes_tot) ? nmodes_start : nmodes_end;

  uint nmodes = nmodes_end - nmodes_start;

  // Initialize indexes
  VECTOR<uint> all_n(nmodes_tot);
  VECTOR<uint> all_k(nmodes_tot);
  for(uint i = 0; i < kmodes; i++)
  for(uint j = 0; j < xmodes; j++)
  {
    uint idx = xmodes * i + j;
    all_n[idx] = j;
    all_k[idx] = i;
  }

  // Take the range of indexes which corresponds to this node
  VECTOR<uint> this_n(all_n.begin() + nmodes_start, all_n.begin() + nmodes_end);
  VECTOR<uint> this_k(all_k.begin() + nmodes_start, all_k.begin() + nmodes_end);
  VECTOR<uint> unique_k(this_k);

  // Obtain soreted list of unique values of k on this node
  std::sort(unique_k.begin(), unique_k.end());
  unique_k.erase( std::unique(unique_k.begin(), unique_k.end()), unique_k.end() );

  knum = unique_k.size();
  n.resize(knum);

  for(uint i = 0; i < knum; i++)
  {
    uint myk = unique_k[i];
    k.push_back(myk);
    for(uint j = 0; j < nmodes; j++)
    {
      if (this_k[j] == myk)
        n[i].push_back(this_n[j]);
    }
  }
}

#include "common.h"
#include "utils.h"
#include <time.h>
#include "common_option_parser.h"
#include <signal.h>

using std::string;
using std::cout;
using std::endl;

static
void sigterm_handler(int signum)
{
  std::cout << std::endl << std::endl;
  std::cerr << std::endl << std::endl;

  switch (signum)
  {
    case SIGTERM:
      std::cout << "SIGTERM";
      std::cerr << "SIGTERM";
      break;
    case SIGINT:
      std::cout << "SIGINT";
      std::cerr << "SIGINT";
      break;
  }

  std::cout << " was caught, exiting." << std::endl;
  std::cerr << " was caught, exiting." << std::endl;

  common_AppFin();

  exit(0);
}

void
common_AppInit(int argc, char **argv, const std::string &app_name)
{
  // Handle SIGTERM and SIGINT signals
  struct sigaction handler_action;
  memset(&handler_action, 0, sizeof(struct sigaction));
  handler_action.sa_handler = sigterm_handler;
  sigaction(SIGTERM, &handler_action, NULL);
  sigaction(SIGINT,  &handler_action, NULL);

#ifdef WITH_MPI_
  VECTOR<char> processor_name(MPI_MAX_PROCESSOR_NAME);
  int p_name_length;

  SAFE_MPI_CALL(MPI_Init(&argc, &argv));
  SAFE_MPI_CALL(MPI_Comm_size(MPI_COMM_WORLD, &pMPISize));
  SAFE_MPI_CALL(MPI_Comm_rank(MPI_COMM_WORLD, &pMPIRank));
  SAFE_MPI_CALL(MPI_Get_processor_name(processor_name.data(), &p_name_length));

  pMPIProcessorName = std::string(processor_name.data());
#endif

  pSeed = time(NULL);

  common_option_parser_Parse(argc, argv, app_name);

  string pSuffix;
#ifdef WITH_MPI_
  pSuffix = "mpi_" + TO_STRING(pMPIRank);
#else
  pSuffix = "";
#endif

  if (!pLogsDirPath.empty())
  {
    // Logs without suffix
    pLogsDir.SetDirectory(pLogsDirPath);
    pLogsDir.SetPrefix(pPrefix);
    pLogsDir.SetSuffix("");

    // Private mpi logs
    pPrivateMPILogsDir.SetDirectory(pLogsDirPath);
    pPrivateMPILogsDir.SetPrefix(pPrefix);
    pPrivateMPILogsDir.SetSuffix(pSuffix);

    // Redirect stdout and stderr to private mpi logfiles
    string main_log = pPrivateMPILogsDir.GetFullPath("main_log.txt");
    string error_log = pPrivateMPILogsDir.GetFullPath("error_log.txt");
    string logs_attr = (pOverwriteLogs) ? "w" : "a";
    pStdLogs.SetStdOut(main_log, logs_attr);
    pStdLogs.SetStdErr(error_log, logs_attr);
  }

  pPrivateMPIDataDir.SetDirectory(pDataDirPath);
  pPrivateMPIDataDir.SetPrefix(pPrefix);
  pPrivateMPIDataDir.SetSuffix(pSuffix);

  pDataDir.SetDirectory(pDataDirPath);
  pDataDir.SetPrefix(pPrefix);
  pDataDir.SetSuffix("");

  omp_set_dynamic(0);
  int n_threads_max = omp_get_max_threads();
  if ((pNthreads < 1) || (pNthreads > n_threads_max))
    pNthreads = n_threads_max;

  omp_set_num_threads(pNthreads);

  Utils::PrintGreetings();

  cout << "Program description: " << app_name << endl;
  cout << "Running with " << pNthreads << " parallel OpenMP threads (maximum available " << n_threads_max << ")" << endl << endl;

#ifdef WITH_MPI_
  cout << "MPI is enabled" << endl;
  cout << "  Size: " << pMPISize << endl;
  cout << "  Rank: " << pMPIRank << endl;
  cout << "  Processor name: " << pMPIProcessorName << endl;
#else
  cout << "MPI is disabled" << endl;
#endif

  cout << endl;

  common_option_parser_PrintParameters(argc, argv);

  init_rand(pSeed);
}

void
common_AppFin()
{
  Utils::PrintEnd();
  pStdLogs.CloseAll();

#ifdef WITH_MPI_
  MPI_Finalize();
#endif
}

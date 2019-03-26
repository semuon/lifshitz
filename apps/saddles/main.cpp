#include "common.h"
#include "option_parser.h"
#include "utils.h"
#include "parameters.h"
#include "vector_field.h"
#include "spinor_field.h"
#include "scalar_field.h"
#include "lattice.h"
#include "formats.h"
#include "linalg.h"
#include <time.h>
#include <mpi_module.h>

// Don't use entire namespace std
using std::cout;
using std::endl;
using std::string;

int main(int argc, char **argv)
{
  common_AppInit(argc, argv, "Lifshitz regime");
  int64_t begin = Utils::GetTimeMs64();

  const string vector_current_name = "vector_current.txt";
 
  FILE *f_vector_current = NULL;

  if (ModuleMPI::IsMasterNode())
  {
    const string f_attr = "wb";
    const string f_tot_attr = "w";

    f_vector_current = pDataDir.OpenFile(vector_current_name, f_tot_attr);
  }

  if (ModuleMPI::IsMasterNode())
  {
    fclose(f_vector_current);
  }

  int64_t end = Utils::GetTimeMs64();

  cout << "Program realized in " << (1.0 * end - 1.0 * begin) / 1000. << " s" << endl;

  common_AppFin();

  return 0;
}

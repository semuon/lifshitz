#include "option_parser.h"
#include "parameters.h"
#include "utils.h"
#include "common_option_parser.h"

static void option_parser_CheckParameters()
{ 
}

void option_parser_Parse(TCLAP::CmdLine &cmd, int argc, char **argv)
{
  TCLAP::ValueArg<double> argTol("", "tol", "Tolerance.", false, pTolerance, "n", cmd);
  TCLAP::ValueArg<int> argNiters("", "n-iters", "Number of iterations.", false, pNiters, "n", cmd);
  TCLAP::ValueArg<int> argNtries("", "n-tries", "Number of distinct random initial conditions for iterations.", false, pNtries, "n", cmd);
  TCLAP::ValueArg<double> argRandRange("", "rand-range", "Initial values of epsilon field will be in the range [-<n>, +<n>].", false, pRandomRange, "n", cmd);
  TCLAP::ValueArg<double> argRelaxAlpha("", "relax-alpha", "Value of relaxation parameter.", false, pRelaxAlpha, "n", cmd);
  
  cmd.parse(argc, argv);

  pTolerance = argTol.getValue();
  pRandomRange = argRandRange.getValue();
  pRelaxAlpha = argRelaxAlpha.getValue();
  pNiters = argNiters.getValue();
  pNtries = argNtries.getValue();

  option_parser_CheckParameters();
}

void option_parser_PrintParameters()
{
  pStdLogs.Write("\nParameters of saddle point search:\n");

  pStdLogs.Write("  Tolerance:                                   %2.4le\n", pTolerance);
  pStdLogs.Write("  Relaxation parameter:                        %2.4le\n", pRelaxAlpha);
  pStdLogs.Write("  Max. iterations:                             %d\n", pNiters);
  pStdLogs.Write("  Num. of tries:                               %d\n", pNtries);
  pStdLogs.Write("  Rand. range:                                 %2.4le\n", pRandomRange);

  pStdLogs.Write("\n");
}

#include "option_parser.h"
#include "parameters.h"
#include "utils.h"
#include "common_option_parser.h"

static void option_parser_CheckParameters()
{
}

void option_parser_Parse(TCLAP::CmdLine &cmd, int argc, char **argv)
{
  TCLAP::ValueArg<int> argHmcNumSteps("", "hmc-num-steps", "Number of steps in one HMC trajectory.", false, pHmcNumSteps, "n", cmd);
  TCLAP::ValueArg<int> argHmcNumConf("", "hmc-num-conf", "Number of HMC configurations.", false, pHmcNumConf, "n", cmd);
  TCLAP::ValueArg<double> argHmcNumDt("", "hmc-dt", "Time step of HMC integrator.", false, pHmcDt, "n", cmd);

  cmd.parse(argc, argv);

  pHmcNumSteps = argHmcNumSteps.getValue();
  pHmcNumConf = argHmcNumConf.getValue();
  pHmcDt = argHmcNumDt.getValue();

  option_parser_CheckParameters();
}

void option_parser_PrintParameters()
{
  pStdLogs.Write("\nParameters of HMC:\n");

  pStdLogs.Write("  Integrator time step:                        %2.4le\n", pHmcDt);
  pStdLogs.Write("  Number of integration steps:                 %d\n", pHmcNumSteps);
  pStdLogs.Write("  Number of configurations:                    %d\n", pHmcNumConf);

  pStdLogs.Write("\n");
}

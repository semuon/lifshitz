#include "option_parser.h"
#include "parameters.h"
#include "utils.h"
#include "common_option_parser.h"

static void option_parser_CheckParameters()
{ 
}

void option_parser_Parse(TCLAP::CmdLine &cmd, int argc, char **argv)
{
  /*TCLAP::ValueArg<double> argDt("", "dt", "Integrator time step. By default is 0.01.", false, 0.01, "n", cmd);
  TCLAP::ValueArg<double> argStartT("", "start-t", "Initial moment of time t0. By default is 0.", false, 0.0, "n", cmd);
  TCLAP::ValueArg<double> argEndT("", "end-t", "Final moment of time. By default is 1.", false, 1.0, "n", cmd);
  TCLAP::ValueArg<uint>   argNmode("", "n-modes", "Number of initial electromagnetic modes. By default is 0.", false, 0, "n", cmd);
  
  ArgProvidedVisitor argNmodefixedVisitor;
  TCLAP::ValueArg<uint>   argNmodefixed("", "n-modes_fixed", "Number of initial electromagnetic modes in direction 1. By default is random.", false, 0, "n", cmd, &argNmodefixedVisitor);
  
  TCLAP::ValueArg<double> argAmplf("", "modes-amplitude", "Amplitude of initial electromagnetic modes. By default is 0.", false, 0, "n", cmd);
  TCLAP::ValueArg<double> argMua("", "mu5", "Chiral chemical potential. By default is 0.", false, 0, "n", cmd);
  TCLAP::ValueArg<double> argFreq("","freq", "Frequence for Sin Bump ext Field. By default 0.", false, 0, "n", cmd);
  TCLAP::ValueArg<int> argKsin("","ksin", "Mode number for Sin Bump ext Field. By default 0.", false, 0, "n", cmd);
  
  TCLAP::SwitchArg switchLoadBackup("", "load-backup", "Load backup files and start simulation from that point. Complete list of backup files should be specified", cmd, false);
  TCLAP::UnlabeledMultiArg<std::string> multiLoadBackupFiles("backupfiles", "Complete list of backup files if --load-backup is specified", false, "backup files list", cmd);
  
  TCLAP::SwitchArg switchNoBackReaction("", "nobackreaction", "Backreaction or not. By default false", cmd, false);
*/
  cmd.parse(argc, argv);
/*
  pMua = argMua.getValue();
  
  if (argNmodefixedVisitor.IsFlagSet())
  {
    pnmodefixed= argNmodefixed.getValue();
    pSetNmodefixed=1;
  }

*/

  option_parser_CheckParameters();
}

void option_parser_PrintParameters()
{
  pStdLogs.Write("\nParameters of saddle point search:\n");

  //pStdLogs.Write("  dt:                                          %2.4le\n", pDt);

  pStdLogs.Write("\n");
}

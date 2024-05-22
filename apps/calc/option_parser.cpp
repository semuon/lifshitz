#include "option_parser.h"
#include "parameters.h"
#include "utils.h"
#include "common_option_parser.h"

static void option_parser_CheckParameters()
{
}

void option_parser_Parse(TCLAP::CmdLine &cmd, int argc, char **argv)
{
  TCLAP::ValueArg<int> argNumSkipFirst("", "num-skip-first", "Skip first <num-skip-first> configurations.", false, pNumSkipFirst, "n", cmd);
  TCLAP::ValueArg<int> argNumSkipLast("", "num-skip-last", "Skip last <num-skip-last> configurations.", false, pNumSkipLast, "n", cmd);
  TCLAP::ValueArg<int> argConfStep("", "conf-step", "Process every <conf-step> configuration.", false, pConfStep, "n", cmd);

  TCLAP::ValueArg<std::string> argFnameConf("", "confs-file", "File name of start configuration. Overloads --start-type option if provided.", true, pFnameConfs, "file name", cmd);

  ArgProvidedVisitor argLatCouplingsVisitor;
  TCLAP::ValueArg<double> argLatK1("", "lat-k1", "Dimensionless lattice nearest neighbour coupling. This will be converted into standard couplings.", false, pLatK1, "n", cmd, &argLatCouplingsVisitor);
  TCLAP::ValueArg<double> argLatK2("", "lat-k2", "Dimensionless lattice next-to-nearest neighbour coupling. This will be converted into standard couplings.", false, pLatK2, "n", cmd, &argLatCouplingsVisitor);
  TCLAP::ValueArg<double> argLatLambda("", "lat-lambda", "Dimensionless quartic coupling. This will be converted into standard couplings.", false, pLatLambda, "n", cmd, &argLatCouplingsVisitor);
  TCLAP::ValueArg<double> argLatKappa("", "lat-kappa", "Dimensionless sextic coupling. This will be converted into standard couplings.", false, pLatKappa, "n", cmd, &argLatCouplingsVisitor);

  TCLAP::SwitchArg switchIsComputeCorr("", "compute-corr", "Compute correlation function", cmd, false);
  TCLAP::SwitchArg switchIsVolAvgCorr("", "vol-avg-corr", "Average correlation function over volume", cmd, false);
  TCLAP::SwitchArg switchIsFullCorr("", "full-corr", "Compute full two-pont function C(x,y)", cmd, false);

  TCLAP::ValueArg<int> argBlockSize("", "block-size", "Size of block average.", false, pConfStep, "n", cmd);


  cmd.parse(argc, argv);

  pFnameConfs = argFnameConf.getValue();
  
  pNumSkipFirst = argNumSkipFirst.getValue();
  pNumSkipLast = argNumSkipLast.getValue();
  pConfStep = argConfStep.getValue();
  
  pLatK1 = argLatK1.getValue();
  pLatK2 = argLatK2.getValue();
  pLatLambda = argLatLambda.getValue();
  pLatKappa = argLatKappa.getValue();

  pIsLatticeParamsSet = argLatCouplingsVisitor.IsFlagSet();

  pIsComputeCorr = switchIsComputeCorr.getValue();
  pIsFullCorr = switchIsFullCorr.getValue();
  pIsVolAvgCorr = switchIsVolAvgCorr.getValue();

  pBlockSize = argBlockSize.getValue();

  option_parser_CheckParameters();
}

void option_parser_PrintParameters()
{
  if (pIsLatticeParamsSet)
  {
    pStdLogs.Write("\nLattice couplings are set. The above couplings will be discarded.\n");
    pStdLogs.Write("  k1:                                         % -2.4le\n", pLatK1);
    pStdLogs.Write("  k2:                                         % -2.4le\n", pLatK2);
    pStdLogs.Write("  lambda:                                     % -2.4le\n", pLatLambda);
    pStdLogs.Write("  kappa:                                      % -2.4le\n", pLatKappa);
  }

  pStdLogs.Write("\nParameters of calculation:\n");

  pStdLogs.Write("  Configurations file:                         ");
  std::cout << pFnameConfs << std::endl;

  pStdLogs.Write("  Skip first #confs:                          % -d\n", pNumSkipFirst);
  pStdLogs.Write("  Skip last #confs:                           % -d\n", pNumSkipLast); 
  pStdLogs.Write("  Process each #conf:                         % -d\n", pConfStep);

  pStdLogs.Write("\n  Compute correlation function:                %s\n", (pIsComputeCorr) ? "YES" : "NO");
  if (pIsComputeCorr)
  {
    pStdLogs.Write("    Average over volume:                       %s\n", (pIsVolAvgCorr) ? "YES" : "NO");
    pStdLogs.Write("    Full C(x,y):                               %s\n", (pIsFullCorr) ? "YES" : "NO");
    pStdLogs.Write("    Block size:                               % -d\n", pBlockSize);
  }

  pStdLogs.Write("\n");
}

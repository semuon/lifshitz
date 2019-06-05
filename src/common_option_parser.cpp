#include "common_option_parser.h"

static void common_option_parser_ParseLatticeSizes(std::string &str);
static void common_option_parser_CheckParameters();

void common_option_parser_Parse(int argc, char **argv, const std::string &app_name)
{
  try
  {
    TCLAP::CmdLine cmd(app_name.c_str(), ' ', "1.0");

    ArgProvidedVisitor argSeedVisitor;
    TCLAP::ValueArg<uint> argSeed("", "seed", "Seed for random number generator. By default current time is used.", false, 0, "path", cmd, &argSeedVisitor);

    TCLAP::ValueArg<std::string> argPrefix("", "prefix", "Prefix which will appear in all filenames: prefix.filename.", false, "", "prefix", cmd);
    TCLAP::ValueArg<int> argNthreads("", "n-threads", "Number of OpenMP threads. By default is 1. \
                                                       Use argument 0 in order to determine the value automatically (maximal available).", 
                                                      false, 1, "n", cmd);

    TCLAP::SwitchArg switchOverwriteLogs("", "append-logs", "Append log files if they exist", cmd, true);
    TCLAP::ValueArg<std::string> argLogsDir("", "logs-dir", "Directory for log files. By default stdout and stderr are used.", false, "", "path", cmd);
    TCLAP::ValueArg<std::string> argDataDir("", "data-dir", "Directory for data files. By default current directory is used.", false, ".", "path", cmd);
    TCLAP::ValueArg<uint> argDim("", "dim", "Number of spatial dimensions. By default is 3.", false, 3, "dimension", cmd);
    TCLAP::ValueArg<std::string> argL("", "L", "Lattice size", true, "5,5,5", "L1,L2,L3,...", cmd);

    TCLAP::ValueArg<double> argInvM2("", "inv-M2", "Inverse of heavy mass scale M^2", false, 0, "real number", cmd);
    TCLAP::ValueArg<double> argm2("", "m2", "Mass (squared)", false, 1.0, "real number", cmd);
    TCLAP::ValueArg<double> argZ("", "Z", "Spatial derivative factor", false, 1.0, "real number", cmd);
    TCLAP::ValueArg<double> argLambda("", "lambda-N", "t'Hooft coupling lambda * N", false, 1.0, "real number", cmd);

    try
    {
      option_parser_Parse(cmd, argc, argv);
    }
    catch(WrongCmdLineOption &ex)
    {
      std::cerr << "Wrong option value: " << ex.what() << std::endl;
      std::exit(-1);
    }

    pInvM2 = argInvM2.getValue();
    pm2 = argm2.getValue();
    pZ = argZ.getValue();
    pLambdaN = argLambda.getValue();

    pLogsDirPath = argLogsDir.getValue();
    pDataDirPath = argDataDir.getValue();
    pOverwriteLogs = switchOverwriteLogs.getValue();
    std::string strL = argL.getValue();
    pDim = argDim.getValue();
    pL.resize(pDim);
    common_option_parser_ParseLatticeSizes(strL);
    
    pNthreads = argNthreads.getValue();
    pPrefix = argPrefix.getValue();
    
    if (argSeedVisitor.IsFlagSet()) pSeed = argSeed.getValue();

    try
    {
      common_option_parser_CheckParameters();
    }
    catch(WrongCmdLineOption &ex)
    {
      std::cerr << "Wrong option value: " << ex.what() << std::endl;
      std::exit(-1);
    }

  }
  catch (TCLAP::ArgException &e)
  {
    std::cerr << "Option error: " << e.error() << " for " << e.argId() << std::endl;
    std::exit(-1);
  }
}

static void
common_option_parser_ParseLatticeSizes(std::string &str)
{
  pL.clear();

  VECTOR<int> vect;
  std::stringstream ss(str);
  int num;

  while (ss >> num)
  {
    vect.push_back(num);

    if (ss.peek() == ',')
        ss.ignore();
  }

  if (vect.size() != pDim)
    THROW_EXCEPTION(WrongCmdLineOption, "Number of lattice sizes should be equal to number of dimesnions");

  for (uint i = 0; i < vect.size(); i++)
    pL.push_back(vect.at(i));
}

static void
common_option_parser_CheckParameters()
{
  //if (pT < 0)
  //  THROW_EXCEPTION(WrongCmdLineOption, "Temperature cannot be negative");
}

void
common_option_parser_PrintParameters(int argc, char **argv)
{
  pStdLogs.Write("Command line:\n");

  for (int i = 0; i < argc; i++)
    pStdLogs.Write("%s ", argv[i]);

  pStdLogs.Write("\n\n");

  common_option_parser_CheckParameters();

  if (!pPrefix.empty())
  {
    pStdLogs.Write("Prefix:                                        %s\n", pPrefix.c_str());
  }
  if (!pLogsDirPath.empty())
  {
    pStdLogs.Write("Logs output directory:                         %s\n", pLogsDirPath.c_str());
    pStdLogs.Write("  Overwrite logs:                              %s\n", (pOverwriteLogs) ? "YES" : "NO");
  }

  pStdLogs.Write("Data output directory:                         %s\n", pDataDir.GetDirectory().c_str());

  if (pSetSeed)
    pStdLogs.Write("Seed for ranlxd:                              % -d\n", pSeed);
  else
    pStdLogs.Write("Seed for ranlxd (automatic):                  % -d\n", pSeed);

  pStdLogs.Write("\n");
  pStdLogs.Write("Run parameters:\n");

  pStdLogs.Write("  L:                                           ");
  for(uint dir = 0; dir < pDim - 1; dir++)
    pStdLogs.Write("%-3d x   ", pL[dir]);
  pStdLogs.Write("%-3d\n", pL[pDim - 1]);

  pStdLogs.Write("  Dim:                                        % -d\n", pDim);

  pStdLogs.Write("  lambda*N:                                   % -2.4le\n", pLambdaN);
  pStdLogs.Write("  m^2:                                        % -2.4le\n", pm2);
  pStdLogs.Write("  1/M^2:                                      % -2.4le\n", pInvM2);
  pStdLogs.Write("  Z:                                          % -2.4le\n", pZ);
  
  option_parser_PrintParameters();
}

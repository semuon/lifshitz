#include "option_parser.h"
#include "parameters.h"
#include "utils.h"
#include "common_option_parser.h"

// Auxiliary structures for --start-type option
struct StartTypeParser
{
  typedef TCLAP::ValueLike ValueCategory;
  tStartConfigurationType type;

  StartTypeParser() : type(START_CONFIGURATION_RANDOM) {}
  StartTypeParser(const tStartConfigurationType &val) : type(val) {}

  bool operator==(const StartTypeParser &other) const { return (type == other.type); }
};

std::ostream& operator<<(std::ostream& os, const StartTypeParser &val)
{
  switch (val.type)
  {
    case START_CONFIGURATION_ZERO:       os << "zero"; break;
    case START_CONFIGURATION_RANDOM:     os << "random"; break;
    case START_CONFIGURATION_LOAD:       os << "load"; break;
    default:                             os.setstate(std::ios::failbit);
  }

  return os;
}

std::istream& operator>>(std::istream& is, StartTypeParser &val)
{
  std::string arg;
  std::getline(is, arg);

  if      (arg == "zero")                   val.type = START_CONFIGURATION_ZERO;
  else if (arg == "random")                 val.type = START_CONFIGURATION_RANDOM;
  else if (arg == "load")                   val.type = START_CONFIGURATION_LOAD;
  else                                      is.setstate(std::ios::failbit);

  return is;
}

// Auxiliary structures for --integrator option
struct IntegratorTypeParser
{
  typedef TCLAP::ValueLike ValueCategory;
  tIntegratorType type;

  IntegratorTypeParser() : type(INTEGRATOR_LEAPFROG) {}
  IntegratorTypeParser(const tIntegratorType &val) : type(val) {}

  bool operator==(const IntegratorTypeParser &other) const { return (type == other.type); }
};

std::ostream& operator<<(std::ostream& os, const IntegratorTypeParser &val)
{
  switch (val.type)
  {
    case INTEGRATOR_LEAPFROG:       os << "leapfrog"; break;
    case INTEGRATOR_OMELYAN:        os << "omelyan"; break;
    default:                             os.setstate(std::ios::failbit);
  }

  return os;
}

std::istream& operator>>(std::istream& is, IntegratorTypeParser &val)
{
  std::string arg;
  std::getline(is, arg);

  if      (arg == "leapfrog")               val.type = INTEGRATOR_LEAPFROG;
  else if (arg == "omelyan")                val.type = INTEGRATOR_OMELYAN;
  else                                      is.setstate(std::ios::failbit);

  return is;
}

static void option_parser_CheckParameters()
{
}

void option_parser_Parse(TCLAP::CmdLine &cmd, int argc, char **argv)
{
  TCLAP::ValueArg<double> argKappa("", "kappa", "Sextic coupling constant: (kappa * phi^6) / 6", false, pKappa, "n", cmd);
  TCLAP::ValueArg<int> argN("", "N", "Rank of the symmetry group O(N).", false, pN, "n", cmd);

  TCLAP::ValueArg<int> argHmcNumSteps("", "hmc-num-steps", "Number of steps in one HMC trajectory.", false, pHmcNumSteps, "n", cmd);
  TCLAP::ValueArg<int> argHmcNumConf("", "hmc-num-conf", "Number of HMC configurations.", false, pHmcNumConf, "n", cmd);
  TCLAP::ValueArg<int> argHmcNumConfStep("", "hmc-num-conf-step", "Save every <hmc-num-conf-step> configuration.", false, pHmcNumConfStep, "n", cmd);
  TCLAP::ValueArg<double> argHmcNumDt("", "hmc-dt", "Time step of HMC integrator.", false, pHmcDt, "n", cmd);

  TCLAP::SwitchArg switchAutoTune("", "auto-tune", "Enable automatic adjustment of HMC MD time step.", cmd, false);
  TCLAP::ValueArg<double> argAutoTuneK("", "auto-tune-K", "Damping coefficient for automatic adjustment of HMC MD time step.", false, pAutoTuneK, "n", cmd);

  StartTypeParser default_start_type;
  std::vector<StartTypeParser> allowed_start_types;
  allowed_start_types.push_back(START_CONFIGURATION_ZERO);
  allowed_start_types.push_back(START_CONFIGURATION_RANDOM);
  TCLAP::ValuesConstraint<StartTypeParser> allowedVals( allowed_start_types );
  TCLAP::ValueArg<StartTypeParser> argStartType("", "start-type", "Type of start configuration", false, default_start_type, &allowedVals, cmd);

  IntegratorTypeParser default_integrator_type;
  std::vector<IntegratorTypeParser> allowed_integrator_types;
  allowed_integrator_types.push_back(INTEGRATOR_LEAPFROG);
  allowed_integrator_types.push_back(INTEGRATOR_OMELYAN);
  TCLAP::ValuesConstraint<IntegratorTypeParser> allowedIntegratorVals( allowed_integrator_types );
  TCLAP::ValueArg<IntegratorTypeParser> argIntegratorType("", "integrator", "Type of HMC integrator", false, default_integrator_type, &allowedIntegratorVals, cmd);

  ArgProvidedVisitor argFnameConfVisitor;
  TCLAP::ValueArg<std::string> argFnameConf("", "start-conf", "File name of start configuration. Overloads --start-type option if provided.", false, pFnameStartConf, "file name", cmd, &argFnameConfVisitor);

  ArgProvidedVisitor argLatCouplingsVisitor;
  TCLAP::ValueArg<double> argLatK1("", "lat-k1", "Dimensionless lattice nearest neighbour coupling. This will be converted into standard couplings.", false, pLatK1, "n", cmd, &argLatCouplingsVisitor);
  TCLAP::ValueArg<double> argLatK2("", "lat-k2", "Dimensionless lattice next-to-nearest neighbour coupling. This will be converted into standard couplings.", false, pLatK2, "n", cmd, &argLatCouplingsVisitor);
  TCLAP::ValueArg<double> argLatLambda("", "lat-lambda", "Dimensionless quartic coupling. This will be converted into standard couplings.", false, pLatLambda, "n", cmd, &argLatCouplingsVisitor);
  TCLAP::ValueArg<double> argLatKappa("", "lat-kappa", "Dimensionless sextic coupling. This will be converted into standard couplings.", false, pLatKappa, "n", cmd, &argLatCouplingsVisitor);

  cmd.parse(argc, argv);

  pN = argN.getValue();
  pKappa = argKappa.getValue();

  pFnameStartConf = argFnameConf.getValue();
  pStartType = (argFnameConfVisitor.IsFlagSet()) ? START_CONFIGURATION_LOAD : argStartType.getValue().type;

  pHmcNumSteps = argHmcNumSteps.getValue();
  pHmcNumConf = argHmcNumConf.getValue();
  pHmcNumConfStep = argHmcNumConfStep.getValue();
  pHmcDt = argHmcNumDt.getValue();

  pLatK1 = argLatK1.getValue();
  pLatK2 = argLatK2.getValue();
  pLatLambda = argLatLambda.getValue();
  pLatKappa = argLatKappa.getValue();

  pIsLatticeParamsSet = argLatCouplingsVisitor.IsFlagSet();

  pIntegratorType = argIntegratorType.getValue().type;

  pAutoTune = switchAutoTune.getValue();
  pAutoTuneK = argAutoTuneK.getValue();

  option_parser_CheckParameters();
}

void option_parser_PrintParameters()
{
  pStdLogs.Write("  kappa:                                      % -2.4le\n", pKappa);
  pStdLogs.Write("  N:                                          % -d\n", pN);

  if (pIsLatticeParamsSet)
  {
    pStdLogs.Write("\nLattice couplings are set. The above couplings will be discarded.\n");
    pStdLogs.Write("  k1:                                         % -2.4le\n", pLatK1);
    pStdLogs.Write("  k2:                                         % -2.4le\n", pLatK2);
    pStdLogs.Write("  lambda:                                     % -2.4le\n", pLatLambda);
    pStdLogs.Write("  kappa:                                      % -2.4le\n", pLatKappa);
  }

  pStdLogs.Write("\nParameters of HMC:\n");

  pStdLogs.Write("  Start configuration:                         ");
  switch (pStartType)
  {
    case START_CONFIGURATION_ZERO:        pStdLogs.Write("zero\n"); break;
    case START_CONFIGURATION_RANDOM:      pStdLogs.Write("random\n"); break;
    case START_CONFIGURATION_LOAD:        std::cout << pFnameStartConf << std::endl; break;
    default:                              pStdLogs.Write("UNKNOWN\n"); break;;
  }

  pStdLogs.Write("  Integrator:                                  ");
  switch (pIntegratorType)
  {
    case INTEGRATOR_LEAPFROG:             pStdLogs.Write("leapfrog\n"); break;
    case INTEGRATOR_OMELYAN:              pStdLogs.Write("omelyan\n"); break;
    default:                              pStdLogs.Write("UNKNOWN\n"); break;;
  }

  pStdLogs.Write("  Integrator time step:                       % -2.4le\n", pHmcDt);
  pStdLogs.Write("  Number of integration steps:                % -d\n", pHmcNumSteps);
  pStdLogs.Write("  Number of configurations:                   % -d\n", pHmcNumConf);
  pStdLogs.Write("  Save configuration on each step:            % -d\n", pHmcNumConfStep);

  pStdLogs.Write("  Auto-adjustment of time step:                %s\n", (pAutoTune) ? "YES" : "NO");
  if (pAutoTune)
  {
    pStdLogs.Write("    Damping coefficient K:                    % -2.4le\n", pAutoTuneK);
  }

  pStdLogs.Write("\n");
}

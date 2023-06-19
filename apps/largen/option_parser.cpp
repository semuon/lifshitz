#include "option_parser.h"
#include "parameters.h"
#include "utils.h"
#include "common_option_parser.h"

// Auxiliary structures for --method option
struct NewtonMethodParser
{
  typedef TCLAP::ValueLike ValueCategory;
  tNewtonMethod type;

  NewtonMethodParser() : type(ITERATIONS_NEWTON) {}
  NewtonMethodParser(const tNewtonMethod &val) : type(val) {}

  bool operator==(const NewtonMethodParser &other) const { return (type == other.type); }
};

std::ostream& operator<<(std::ostream& os, const NewtonMethodParser &val)
{
  switch (val.type)
  {
    case ITERATIONS_NEWTON:       os << "newton"; break;
    case ITERATIONS_HALLEY:       os << "halley"; break;
    default:                             os.setstate(std::ios::failbit);
  }

  return os;
}

std::istream& operator>>(std::istream& is, NewtonMethodParser &val)
{
  std::string arg;
  std::getline(is, arg);

  if      (arg == "newton")                 val.type = ITERATIONS_NEWTON;
  else if (arg == "halley")                 val.type = ITERATIONS_HALLEY;
  else                                      is.setstate(std::ios::failbit);

  return is;
}

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
  //TCLAP::ValueArg<int> argNStencilPts("", "n-stencil", "Number of points of finite difference stencil, by default is 2", false, pNStencilPts, "n", cmd);

  NewtonMethodParser default_method;
  std::vector<NewtonMethodParser> allowed_methods;
  allowed_methods.push_back(ITERATIONS_NEWTON);
  allowed_methods.push_back(ITERATIONS_HALLEY);
  TCLAP::ValuesConstraint<NewtonMethodParser> allowedVals( allowed_methods );
  TCLAP::ValueArg<NewtonMethodParser> argMethod("", "method", "Iterations method", false, default_method, &allowedVals, cmd);

  cmd.parse(argc, argv);

  pTolerance = argTol.getValue();
  pRandomRange = argRandRange.getValue();
  pRelaxAlpha = argRelaxAlpha.getValue();
  pNiters = argNiters.getValue();
  pNtries = argNtries.getValue();
  pMethod = argMethod.getValue().type;
  //pNStencilPts = argNStencilPts.getValue();

  option_parser_CheckParameters();
}

void option_parser_PrintParameters()
{
  pStdLogs.Write("\nParameters of saddle point search:\n");

  pStdLogs.Write("  Iterations:                                  ");
  switch (pMethod)
  {
    case ITERATIONS_NEWTON:               pStdLogs.Write("Newton\n"); break;
    case ITERATIONS_HALLEY:               pStdLogs.Write("Halley\n"); break;
    default:                              pStdLogs.Write("UNKNOWN\n"); break;;
  }

  pStdLogs.Write("  Tolerance:                                   %2.4le\n", pTolerance);
  pStdLogs.Write("  Relaxation parameter:                        %2.4le\n", pRelaxAlpha);
  pStdLogs.Write("  Max. iterations:                             %d\n", pNiters);
  pStdLogs.Write("  Num. of tries:                               %d\n", pNtries);
  pStdLogs.Write("  Rand. range:                                 %2.4le\n", pRandomRange);
  //pStdLogs.Write("  Number of stencil points:                    %d\n", pNStencilPts);

  pStdLogs.Write("\n");
}

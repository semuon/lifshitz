#ifndef COMMON_OPTION_PARSER_H_
#define COMMON_OPTION_PARSER_H_ 

#include "common.h"
#include "tclap/CmdLine.h"

// Main function to parse command line options
void common_option_parser_Parse(int argc, char **argv, const std::string &app_name);
void common_option_parser_PrintParameters(int argc, char **argv);

// ==========================================================
// Application MUST provide the following functions:
// ==========================================================

// This function should add to cmd all needed additional arguments
// and then call cmd.parse(argc, argv).
// All exceptions will be caugth at calling function, no need for try/catch.
// If value of some parameter is wrong then WrongCmdLineOption should be thrown (see below)
extern void option_parser_Parse(TCLAP::CmdLine &cmd, int argc, char **argv);

// This function should print out all application specific parameters to main log  
extern void option_parser_PrintParameters();

// ==========================================================
// Auxiliary classes
// ==========================================================
// Exception which should be thrown if some parameter is wrong
class WrongCmdLineOption : public std::runtime_error
{
public:
  WrongCmdLineOption(const std::string &msg) : std::runtime_error(msg) {}
  WrongCmdLineOption(const char *msg) : std::runtime_error(msg) {}
};

// Visitor class which helps to detect if the option have been provided
// See http://tclap.sourceforge.net/manual.html#VISITORS for details
class ArgProvidedVisitor : public TCLAP::Visitor
{
private:
  bool flag_set;
public:
  ArgProvidedVisitor() : flag_set(false) {}
  virtual void visit() { flag_set = true; }

  bool IsFlagSet() const { return flag_set; }
};


#endif /* COMMON_OPTION_PARSER_H_ */
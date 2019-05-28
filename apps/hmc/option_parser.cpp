#include "option_parser.h"
#include "parameters.h"
#include "utils.h"
#include "common_option_parser.h"

static void option_parser_CheckParameters()
{ 
}

void option_parser_Parse(TCLAP::CmdLine &cmd, int argc, char **argv)
{
  cmd.parse(argc, argv);

  option_parser_CheckParameters();
}

void option_parser_PrintParameters()
{
}

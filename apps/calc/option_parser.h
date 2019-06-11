#ifndef OPTION_PARSER_H_
#define OPTION_PARSER_H_

#include "common.h"
#include "common_option_parser.h"

// Functions as needed by common_option_parser.h

void option_parser_Parse(TCLAP::CmdLine &cmd, int argc, char **argv); 
void option_parser_PrintParameters();

#endif /* OPTION_PARSER_H_ */

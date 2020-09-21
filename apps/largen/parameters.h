#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <common.h>

typedef enum EnumNewtonMethod
{
  ITERATIONS_NEWTON,
  ITERATIONS_HALLEY
} tNewtonMethod;

extern int pNiters;
extern int pNtries;
extern double pTolerance;
extern double pRelaxAlpha;
extern double pRandomRange;
extern tNewtonMethod pMethod;

#endif

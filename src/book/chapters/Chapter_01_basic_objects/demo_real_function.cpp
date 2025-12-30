#include "MMLBase.h"

#include "mml/base/Function.h"

using namespace MML;

// two simple ways to create RealFunction - function pointer and lambda

Real already_existing_func(Real x) {
  return sin(x) * (1 + x * x / 2);
}

void Real_functions_case_1_usage()
{
	// creating a function object from an already existing (standalone) function
	RealFunction fReal(already_existing_func);

	// we can also use lambda directly
	RealFunction fReal2( [](Real x) { return sin(x) * (1 + x * x / 2); } );
}
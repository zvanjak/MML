#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "core/Function.h"
#include "core/InterpolatedFunction.h"

#include "algorithms/FunctionAnalyzers.h"
#endif

#include "../test_data/real_functions_test_bed.h"


using namespace MML;

void Readme_function_analyzer()
{
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "****                  README - function analyzer                   ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	auto fTan = TestBeds::RealFunctionsTestBed::getTestFunctionReal("Tan");
	RealFunctionAnalyzer anTan(fTan._func, "tan(x)");
	anTan.PrintIntervalAnalysis(-5.0, 5.0, 50, 1e-4);

	auto fExp = TestBeds::RealFunctionsTestBed::getTestFunctionReal("Exp");
	RealFunctionAnalyzer anExp(fExp._func, "exp(x)");
	anExp.PrintIntervalAnalysis(-5.0, 5.0, 50, 1e-4);

	RealFunction stepFunc([](Real x) { 
			if( x < 0) return 0.0;
			else if( x > 0) return 1.0;
			else return 0.5;
		});
	RealFunctionAnalyzer anStep(stepFunc, "step(x)");
	anStep.PrintIntervalAnalysis(-5.0, 5.0, 50, 1e-4);

	RealFunction test1([](Real x) { return 1 / (x - 1); });
	RealFunctionAnalyzer an(test1, "1 / (x - 1)");
	an.PrintIntervalAnalysis(-5.0, 5.0, 50, 1e-4);

/* OUTPUT
f(x) = tan(x) - Function analysis in interval [-5.00000000, 5.00000000] with 50 points:
  Defined    : yes
  Continuous : yes
  Monotonic  : no
  Min        : -34.23253274
  Max        : 34.23253274
f(x) = exp(x) - Function analysis in interval [-5.00000000, 5.00000000] with 50 points:
  Defined    : yes
  Continuous : yes
  Monotonic  : yes
  Min        : 0.00673795
  Max        : 121.51041752
f(x) = step(x) - Function analysis in interval [-5.00000000, 5.00000000] with 50 points:
  Defined    : yes
  Continuous : no  Not continuous at points: 0.00000000
  Monotonic  : no
  Min        : 0.00000000
  Max        : 1.00000000
f(x) = 1 / (x - 1) - Function analysis in interval [-5.00000000, 5.00000000] with 50 points:
  Defined    : no  Not defined at points: 1.00000000
  Continuous : yes
  Monotonic  : no
  Min        : -5.00000000
  Max        : inf
*/
}
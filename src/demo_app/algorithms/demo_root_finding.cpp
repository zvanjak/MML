#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Function.h"

#include "algorithms/RootFinding.h"
#endif

using namespace std;
using namespace MML;

void Demo_Root_finding()
{
	std::cout << endl;
	std::cout << "***********************************************************************" << endl;
	std::cout << "****                       ROOT FINDING                            ****" << endl;
	std::cout << "***********************************************************************" << endl;

	RealFunction f([](Real x) { return x * x - 2; });
	Real a = 0.0;
	Real b = 2.0;
	bool isBracketed = RootFinding::BracketRoot(f, a, b);
}
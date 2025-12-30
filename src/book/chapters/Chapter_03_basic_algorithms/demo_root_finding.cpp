#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "mml/base/Function.h"

#include "mml/algorithms/RootFinding.h"

#include "mml/tools/Visualizer.h"
#endif

using namespace MML;
using namespace MML::RootFinding;


void Demo_Root_finding()
{
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "****                   EXAMPLE 2 - Root finding                      ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	// Finding the root of the function f(x) = x^2 - 2 
	// (with known root sqrt(2) = 1.414213562373095)
	RealFunction f{ [](Real x) { return x * x - 2; } };
	Real a = 0.0;
	Real b = 10.0;
	
	Real root = FindRootBisection(f, a, b, 0.0001);
	
	std::cout << "Root of the function using the bisection method: " << root << std::endl;
	std::cout << "Exact root: " << sqrt(2) << std::endl;
	std::cout << "Error: " << std::abs(root - sqrt(2)) << std::endl;

	// example for function with multiple roots
	Real x1 = -20.0;
	Real x2 = 20.0;
	RealFunction f2{ [](Real x) { return x * x * sin(x) * exp(2-Abs(x/2)); } };

	Vector<Real> root_brack_x1(10), root_brack_x2(10);
	int	numFoundRoots = FindRootBrackets(f2, x1, x2, 100, root_brack_x1, root_brack_x2);

	std::cout << "Number of found roots: " << numFoundRoots << std::endl;

	Vector<Real> roots(numFoundRoots);
	for (int i = 0; i < numFoundRoots; i++)
	{
		roots[i] = FindRootBisection(f2, root_brack_x1[i], root_brack_x2[i], 1e-7);

		std::cout << "Root " << i << " : " << roots[i] << std::endl;
	}

	// visualize the function
  Visualizer::VisualizeRealFunction(f2, "f2", x1, x2, 500, "Demo_root_finding.txt");
}

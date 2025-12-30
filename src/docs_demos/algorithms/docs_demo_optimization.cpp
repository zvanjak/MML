#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Function.h"
#include "algorithms/Optimization.h"
#include "algorithms/Optimization/OptimizationMultidim.h"
#endif

using namespace MML;

///////////////////////////////////////////////////////////////////////////////////////////
///                    BRACKETING A MINIMUM                                             ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_Optimization_Bracketing()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Bracketing a Minimum\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nBefore minimizing, we need to bracket the minimum: find (a,b,c)\n";
	std::cout << "such that f(b) < f(a) and f(b) < f(c), guaranteeing min in [a,c].\n";
	
	// Simple parabola with minimum at x=2
	RealFunctionFromStdFunc func([](Real x) { return (x - 2.0) * (x - 2.0) + 1.0; });
	
	std::cout << "\n--- f(x) = (x-2)² + 1  (minimum at x=2) ---\n";
	
	// Find bracket starting from [0, 1]
	MinimumBracket bracket = Minimization::BracketMinimum(func, 0.0, 1.0);
	
	std::cout << "\nBracket found:\n";
	std::cout << "  ax = " << bracket.ax << ", f(ax) = " << bracket.fa << std::endl;
	std::cout << "  bx = " << bracket.bx << ", f(bx) = " << bracket.fb << std::endl;
	std::cout << "  cx = " << bracket.cx << ", f(cx) = " << bracket.fc << std::endl;
	std::cout << "  Valid: " << (bracket.valid ? "yes" : "no") << std::endl;
	
	std::cout << "\nNote: bx is the lowest point, guaranteeing min in [ax, cx]\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                    GOLDEN SECTION SEARCH                                            ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_Optimization_GoldenSection()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Golden Section Search\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nGolden section: simple, robust, linear convergence rate.\n";
	std::cout << "Uses golden ratio φ = (1+√5)/2 ≈ 1.618 for optimal bracketing.\n";
	
	// Function with minimum at x = 3
	RealFunctionFromStdFunc func([](Real x) { 
		return (x - 3.0) * (x - 3.0) + 2.0; 
	});
	
	std::cout << "\n--- f(x) = (x-3)² + 2  (minimum at x=3, f=2) ---\n";
	
	// Bracket first
	MinimumBracket bracket = Minimization::BracketMinimum(func, 0.0, 1.0);
	
	// Golden section search
	MinimizationResult result = Minimization::GoldenSectionSearch(func, bracket);
	
	std::cout << "\nGolden Section Result:\n";
	std::cout << "  x_min = " << result.xmin << " (exact: 3.0)\n";
	std::cout << "  f_min = " << result.fmin << " (exact: 2.0)\n";
	std::cout << "  Iterations: " << result.iterations << std::endl;
	std::cout << "  Converged: " << (result.converged ? "yes" : "no") << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                    BRENT'S METHOD                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_Optimization_Brent()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Brent's Method (Parabolic Interpolation)\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nBrent's method: combines golden section with parabolic interpolation.\n";
	std::cout << "Superlinear convergence for smooth functions, never worse than golden.\n";
	
	// More complex function
	RealFunctionFromStdFunc func([](Real x) { 
		return std::sin(x) + 0.1 * x * x - 2.0;  // min around x ≈ -0.45
	});
	
	std::cout << "\n--- f(x) = sin(x) + 0.1*x² - 2 ---\n";
	
	// Bracket in [-2, 2]
	MinimumBracket bracket = Minimization::BracketMinimum(func, -2.0, 2.0);
	
	// Brent's method
	MinimizationResult result = Minimization::BrentMinimize(func, bracket);
	
	std::cout << "\nBrent's Method Result:\n";
	std::cout << "  x_min = " << result.xmin << std::endl;
	std::cout << "  f_min = " << result.fmin << std::endl;
	std::cout << "  Iterations: " << result.iterations << std::endl;
	
	// Verify
	std::cout << "\nVerification: f'(x_min) ≈ " 
	          << (func(result.xmin + 1e-6) - func(result.xmin - 1e-6)) / 2e-6 
	          << " (should be ~0)\n";
}

void Docs_Demo_Optimization_Brent_WithDeriv()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Brent's Method with Derivatives\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nWhen derivatives are available, convergence is faster.\n";
	
	// Function and its derivative
	RealFunctionFromStdFunc func([](Real x) { 
		return x * x * x * x - 3.0 * x * x + 2.0;  // x⁴ - 3x² + 2
	});
	RealFunctionFromStdFunc dfunc([](Real x) { 
		return 4.0 * x * x * x - 6.0 * x;  // 4x³ - 6x
	});
	
	std::cout << "\n--- f(x) = x⁴ - 3x² + 2, f'(x) = 4x³ - 6x ---\n";
	std::cout << "Minima at x = ±√(3/2) ≈ ±1.225\n";
	
	// Bracket around positive minimum
	MinimumBracket bracket = Minimization::BracketMinimum(func, 0.5, 2.0);
	
	// Brent with derivatives
	MinimizationResult result = Minimization::BrentMinimizeWithDeriv(func, dfunc, bracket);
	
	std::cout << "\nBrent with Derivatives Result:\n";
	std::cout << "  x_min = " << result.xmin << " (exact: " << std::sqrt(1.5) << ")\n";
	std::cout << "  f_min = " << result.fmin << std::endl;
	std::cout << "  f'(x_min) = " << dfunc(result.xmin) << " (should be ~0)\n";
	std::cout << "  Iterations: " << result.iterations << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                    1D MINIMIZATION WITH BOUNDS                                      ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_Optimization_Bounded()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Bounded Minimization\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nFind minimum within specified bounds [a, b].\n";
	
	// Function with global minimum outside bounds
	RealFunctionFromStdFunc func([](Real x) { 
		return (x - 5.0) * (x - 5.0);  // min at x=5
	});
	
	std::cout << "\n--- f(x) = (x-5)²  (global min at x=5) ---\n";
	std::cout << "Searching in [0, 3] - minimum is outside this range!\n";
	
	// Minimize in [0, 3]
	MinimumBracket bracket;
	bracket.ax = 0.0;
	bracket.cx = 3.0;
	bracket.bx = 1.5;
	bracket.fa = func(bracket.ax);
	bracket.fb = func(bracket.bx);
	bracket.fc = func(bracket.cx);
	bracket.valid = true;
	
	MinimizationResult result = Minimization::GoldenSectionSearch(func, bracket, 1e-8);
	
	std::cout << "\nConstrained Result:\n";
	std::cout << "  x_min = " << result.xmin << " (at right boundary ≈ 3)\n";
	std::cout << "  f_min = " << result.fmin << " (compared to global min = 0 at x=5)\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                    COMPARING METHODS                                                ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_Optimization_Comparison()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Method Comparison\n";
	std::cout << "==========================================================================\n";
	
	// Rosenbrock-like 1D function (steep valley)
	RealFunctionFromStdFunc func([](Real x) { 
		return 100.0 * std::pow(1.0 - x, 2) + std::pow(x - 1.0, 4);
	});
	
	std::cout << "\n--- Rosenbrock-like: 100*(1-x)² + (x-1)⁴ ---\n";
	std::cout << "Minimum at x = 1\n";
	
	MinimumBracket bracket = Minimization::BracketMinimum(func, -5.0, 5.0);
	
	// Golden section
	MinimizationResult golden = Minimization::GoldenSectionSearch(func, bracket);
	
	// Brent
	MinimizationResult brent = Minimization::BrentMinimize(func, bracket);
	
	std::cout << "\nComparison:\n";
	std::cout << "                Golden Section     Brent\n";
	std::cout << "  x_min:        " << golden.xmin << "              " << brent.xmin << std::endl;
	std::cout << "  f_min:        " << golden.fmin << "           " << brent.fmin << std::endl;
	std::cout << "  Iterations:   " << golden.iterations << "                    " << brent.iterations << std::endl;
	
	std::cout << "\nBrent typically needs fewer iterations due to parabolic acceleration.\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                    MULTIDIMENSIONAL OPTIMIZATION (PREVIEW)                          ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_Optimization_Multidim()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Multidimensional Optimization Overview\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nFor N-dimensional minimization, see Optimization/OptimizationMultidim.h:\n";
	std::cout << "  - Coordinate descent (simple, dimension-by-dimension)\n";
	std::cout << "  - Powell's method (conjugate directions)\n";
	std::cout << "  - Nelder-Mead (simplex, derivative-free)\n";
	std::cout << "  - Gradient descent (with line search)\n";
	std::cout << "  - BFGS (quasi-Newton, fast convergence)\n";
	
	std::cout << "\nExample setup:\n";
	std::cout << "  ScalarFunction<2> f([](const VectorN<Real,2>& x) {\n";
	std::cout << "      return (x[0]-1)*(x[0]-1) + (x[1]-2)*(x[1]-2);\n";
	std::cout << "  });\n";
	std::cout << "  VectorN<Real,2> x0 = {0, 0};\n";
	std::cout << "  auto result = MultidimMinimization::NelderMead(f, x0);\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         MAIN DEMO FUNCTION                                          ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_Optimization()
{
	std::cout << "\n##########################################################################\n";
	std::cout << "#                     OPTIMIZATION DEMOS                                 #\n";
	std::cout << "##########################################################################\n";
	
	Docs_Demo_Optimization_Bracketing();
	Docs_Demo_Optimization_GoldenSection();
	Docs_Demo_Optimization_Brent();
	Docs_Demo_Optimization_Brent_WithDeriv();
	Docs_Demo_Optimization_Bounded();
	Docs_Demo_Optimization_Comparison();
	Docs_Demo_Optimization_Multidim();
	
	std::cout << "\n=== All Optimization Demos Complete ===\n";
}

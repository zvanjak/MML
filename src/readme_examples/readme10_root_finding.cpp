///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        readme10_root_finding.cpp                                           ///
///  Description: README example - Root Finding Algorithms                            ///
///               Demonstrates Bisection, Newton, Brent, Ridders methods              ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Function.h"
#include "algorithms/RootFinding.h"
#endif

#include <iostream>
#include <iomanip>

using namespace MML;

void Readme_RootFinding()
{
    std::cout << std::endl;
    std::cout << "=== Root Finding Algorithms ===" << std::endl;

    // Find root of f(x) = x³ - 2x - 5 (has root near x ≈ 2.0945)
    RealFunction f{[](Real x) { return x*x*x - 2*x - 5; }};

    std::cout << "Finding root of f(x) = x³ - 2x - 5 in [2, 3]" << std::endl;
    std::cout << std::setprecision(15);

    // Compare different methods
    Real root_bisect = RootFinding::FindRootBisection(f, 2.0, 3.0, 1e-12);
    Real root_brent  = RootFinding::FindRootBrent(f, 2.0, 3.0, 1e-12);
    Real root_newton = RootFinding::FindRootNewton(f, 2.0, 3.0, 1e-12);
    Real root_ridder = RootFinding::FindRootRidders(f, 2.0, 3.0, 1e-12);

    std::cout << std::endl << "Method comparison:" << std::endl;
    std::cout << "  Bisection: " << root_bisect << std::endl;
    std::cout << "  Brent:     " << root_brent << std::endl;
    std::cout << "  Newton:    " << root_newton << std::endl;
    std::cout << "  Ridders:   " << root_ridder << std::endl;

    // Get detailed convergence info using config
    RootFinding::RootFindingConfig config;
    config.tolerance = 1e-14;
    config.max_iterations = 100;
    
    auto result = RootFinding::FindRootBrent(f, 2.0, 3.0, config);
    std::cout << std::endl << "Brent with detailed config:" << std::endl;
    std::cout << "  Root:       " << result.root << std::endl;
    std::cout << "  Iterations: " << result.iterations_used << std::endl;
    std::cout << "  f(root) =   " << result.function_value << std::endl;
    std::cout << "  Converged:  " << (result.converged ? "yes" : "no") << std::endl;

    // Find multiple roots by bracketing
    RealFunction g{[](Real x) { return sin(x); }};  // Has roots at 0, π, 2π, ...
    Vector<Real> brackets_lo, brackets_hi;
    int numRoots = RootFinding::FindRootBrackets(g, -1.0, 10.0, 100, brackets_lo, brackets_hi);
    
    std::cout << std::endl << "Finding roots of sin(x) in [-1, 10]:" << std::endl;
    std::cout << "  Found " << numRoots << " root bracket(s)" << std::endl;
    
    for (int i = 0; i < numRoots; i++) {
        Real root = RootFinding::FindRootBisection(g, brackets_lo[i], brackets_hi[i], 1e-10);
        std::cout << "    Root " << (i+1) << ": x = " << root << std::endl;
    }

    std::cout << std::endl;
}

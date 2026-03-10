///////////////////////////////////////////////////////////////////////////////////////////
///  File:        docs_demo_richardson.cpp                                            ///
///  Description: Brief demonstration of RichardsonExtrapolation.h                    ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#include "core/RichardsonExtrapolation.h"
#endif

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

using namespace MML;

void Docs_Demo_Richardson()
{
    std::cout << "=== Richardson Extrapolation Demo ===" << std::endl;
    std::cout << std::fixed << std::setprecision(10);
    
    // Example: Approximate derivative of sin(x) at x=1 using decreasing step sizes
    // f'(x) ≈ (f(x+h) - f(x-h)) / (2h)
    // True value: cos(1) ≈ 0.5403023059
    
    Real x = 1.0;
    std::vector<Real> values;
    std::vector<Real> h_squared;
    
    std::cout << "Central difference approximations for d/dx sin(1):" << std::endl;
    Real h = 0.1;
    for (int i = 0; i < 5; i++) {
        Real approx = (std::sin(x + h) - std::sin(x - h)) / (2 * h);
        values.push_back(approx);
        h_squared.push_back(h * h);  // Error is O(h²)
        std::cout << "  h=" << std::setw(6) << h << ": " << approx << std::endl;
        h /= 2;
    }
    
    // Apply Neville extrapolation
    RichardsonResult result = Richardson::NevilleExtrapolate(values, h_squared);
    
    std::cout << "\nNeville extrapolation to h=0:" << std::endl;
    std::cout << "  Extrapolated: " << result.value << std::endl;
    std::cout << "  Exact cos(1): " << std::cos(1.0) << std::endl;
    std::cout << "  Error est:    " << result.error_estimate << std::endl;
    std::cout << "  Actual error: " << std::abs(result.value - std::cos(1.0)) << std::endl;
    
    std::cout << "=== Demo Complete ===" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///  File:        docs_demo_function_helpers.cpp                                      ///
///  Description: Brief demonstration of FunctionHelpers.h - function utilities       ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#include "base/Function.h"
#include "core/FunctionHelpers.h"
#endif

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace MML;

void Docs_Demo_FunctionHelpers()
{
    std::cout << "=== FunctionHelpers Demo ===" << std::endl;
    std::cout << std::fixed << std::setprecision(6);
    
    // Taylor series expansion
    RealFunctionFromStdFunc sinFunc([](Real x) { return std::sin(x); });
    
    Real a = 0.0;  // Expand around x=0
    PolynomRealFunc taylor2 = TaylorSeries2(sinFunc, a);
    PolynomRealFunc taylor3 = TaylorSeries3(sinFunc, a);
    
    std::cout << "Taylor expansion of sin(x) around x=0:" << std::endl;
    
    Real x = 0.5;
    std::cout << "At x = " << x << ":" << std::endl;
    std::cout << "  sin(x) exact     = " << std::sin(x) << std::endl;
    std::cout << "  Taylor (order 2) = " << taylor2(x) << std::endl;
    std::cout << "  Taylor (order 3) = " << taylor3(x) << std::endl;
    
    // Derivative wrapper
    RealFunctionFromStdFunc cosFunc([](Real x) { return std::cos(x); });
    RealFuncDerived6 cosDerivative(cosFunc);  // d/dx cos(x) = -sin(x)
    
    std::cout << "\nDerivative of cos(x) at x=0.5:" << std::endl;
    std::cout << "  Computed: " << cosDerivative(0.5) << std::endl;
    std::cout << "  Exact (-sin): " << -std::sin(0.5) << std::endl;
    
    // Second derivative wrapper
    RealFuncSecondDerived4 cosSec(cosFunc);  // d²/dx² cos(x) = -cos(x)
    std::cout << "\nSecond derivative of cos(x) at x=0.5:" << std::endl;
    std::cout << "  Computed: " << cosSec(0.5) << std::endl;
    std::cout << "  Exact (-cos): " << -std::cos(0.5) << std::endl;
    
    std::cout << "=== Demo Complete ===" << std::endl;
}

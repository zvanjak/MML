///////////////////////////////////////////////////////////////////////////////////////////
///  File:        docs_demo_numeric_validation.cpp                                    ///
///  Description: Brief demonstration of NumericValidation.h - input validation       ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#include "core/NumericValidation.h"
#endif

#include <iostream>
#include <cmath>
#include <limits>

using namespace MML;

void Docs_Demo_NumericValidation()
{
    std::cout << "=== NumericValidation Demo ===" << std::endl;
    
    // Check if values are finite
    Real good = 3.14;
    Real nan = std::numeric_limits<Real>::quiet_NaN();
    Real inf = std::numeric_limits<Real>::infinity();
    
    std::cout << "IsFinite checks:" << std::endl;
    std::cout << "  IsFinite(3.14) = " << (IsFinite(good) ? "true" : "false") << std::endl;
    std::cout << "  IsFinite(NaN)  = " << (IsFinite(nan) ? "true" : "false") << std::endl;
    std::cout << "  IsFinite(Inf)  = " << (IsFinite(inf) ? "true" : "false") << std::endl;
    
    // Validate finite values
    try {
        ValidateFinite(good, "initial guess");
        std::cout << "\nValidateFinite(3.14) passed" << std::endl;
    } catch (const NumericInputError& e) {
        std::cout << "Error: " << e.what() << std::endl;
    }
    
    try {
        ValidateFinite(nan, "tolerance");
    } catch (const NumericInputError& e) {
        std::cout << "ValidateFinite(NaN) threw: " << e.what() << std::endl;
    }
    
    // Validate bounds
    try {
        ValidateBounds(0.0, 1.0, "bisection");
        std::cout << "ValidateBounds(0, 1) passed" << std::endl;
    } catch (const NumericInputError& e) {
        std::cout << "Error: " << e.what() << std::endl;
    }
    
    // Validate tolerance
    try {
        ValidateTolerance(1e-10, "Newton");
        std::cout << "ValidateTolerance(1e-10) passed" << std::endl;
    } catch (const NumericInputError& e) {
        std::cout << "Error: " << e.what() << std::endl;
    }
    
    std::cout << "=== Demo Complete ===" << std::endl;
}

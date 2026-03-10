///////////////////////////////////////////////////////////////////////////////////////////
///  File:        docs_demo_standard_functions.cpp                                    ///
///  Description: Brief demonstration of StandardFunctions.h - special math functions ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#include "base/StandardFunctions.h"
#endif

#include <iostream>
#include <iomanip>

using namespace MML;

void Docs_Demo_StandardFunctions()
{
    std::cout << "=== StandardFunctions Demo ===" << std::endl;
    
    Real x = 0.5;
    
    // Basic trig functions
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Sin(" << x << ") = " << Functions::Sin(x) << std::endl;
    std::cout << "Cos(" << x << ") = " << Functions::Cos(x) << std::endl;
    std::cout << "Tan(" << x << ") = " << Functions::Tan(x) << std::endl;
    
    // Hyperbolic functions
    std::cout << "Sinh(" << x << ") = " << Functions::Sinh(x) << std::endl;
    std::cout << "Cosh(" << x << ") = " << Functions::Cosh(x) << std::endl;
    
    // Exponential and log
    std::cout << "Exp(" << x << ") = " << Functions::Exp(x) << std::endl;
    std::cout << "Log(" << x << ") = " << Functions::Log(x) << std::endl;
    
    // Special functions
    std::cout << "Erf(" << x << ") = " << Functions::Erf(x) << std::endl;
    std::cout << "TGamma(5) = " << Functions::TGamma(5.0) << " (= 4! = 24)" << std::endl;
    
#if MML_HAS_STD_SPECIAL_FUNCTIONS
    std::cout << "Bessel J0(1) = " << Functions::CylBesselJ(0, 1.0) << std::endl;
    std::cout << "Legendre P2(0.5) = " << Functions::Legendre(2, 0.5) << std::endl;
#else
    std::cout << "(Special functions not available on this platform)" << std::endl;
#endif
    
    std::cout << "=== Demo Complete ===" << std::endl;
}

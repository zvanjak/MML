///////////////////////////////////////////////////////////////////////////////////////////
///  File:        docs_demo_singularity.cpp                                           ///
///  Description: Brief demonstration of SingularityHandling.h                        ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#include "core/SingularityHandling.h"
#endif

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace MML;

void Docs_Demo_Singularity()
{
    std::cout << "=== Singularity Handling Demo ===" << std::endl;
    std::cout << std::fixed << std::setprecision(6);
    
    // Singularity policies
    std::cout << "Available singularity policies:" << std::endl;
    std::cout << "  Throw     - " << SingularityPolicyToString(SingularityPolicy::Throw) << std::endl;
    std::cout << "  ReturnNaN - " << SingularityPolicyToString(SingularityPolicy::ReturnNaN) << std::endl;
    std::cout << "  ReturnInf - " << SingularityPolicyToString(SingularityPolicy::ReturnInf) << std::endl;
    std::cout << "  Clamp     - " << SingularityPolicyToString(SingularityPolicy::Clamp) << std::endl;
    std::cout << "  ReturnZero- " << SingularityPolicyToString(SingularityPolicy::ReturnZero) << std::endl;
    
    // Safe division with different policies
    Real numerator = 1.0;
    Real denominator = 0.0;
    
    std::cout << "\nSafeDivide(1.0, 0.0) with different policies:" << std::endl;
    
    Real result = Singularity::SafeDivide(numerator, denominator, SingularityPolicy::ReturnNaN);
    std::cout << "  ReturnNaN: " << result << " (isnan=" << (std::isnan(result) ? "true" : "false") << ")" << std::endl;
    
    result = Singularity::SafeDivide(numerator, denominator, SingularityPolicy::ReturnInf);
    std::cout << "  ReturnInf: " << result << std::endl;
    
    result = Singularity::SafeDivide(numerator, denominator, SingularityPolicy::Clamp);
    std::cout << "  Clamp:     " << result << " (large but finite)" << std::endl;
    
    result = Singularity::SafeDivide(numerator, denominator, SingularityPolicy::ReturnZero);
    std::cout << "  ReturnZero: " << result << std::endl;
    
    // Check for singularities
    std::cout << "\nSingularity detection:" << std::endl;
    Real r_small = 1e-12;
    Real r_normal = 1.0;
    std::cout << "  IsNearZero(1e-12) = " << (Singularity::IsNearZero(r_small) ? "true" : "false") << std::endl;
    std::cout << "  IsNearZero(1.0)   = " << (Singularity::IsNearZero(r_normal) ? "true" : "false") << std::endl;
    
    std::cout << "=== Demo Complete ===" << std::endl;
}

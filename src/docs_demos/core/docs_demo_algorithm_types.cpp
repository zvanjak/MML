///////////////////////////////////////////////////////////////////////////////////////////
///  File:        docs_demo_algorithm_types.cpp                                       ///
///  Description: Brief demonstration of AlgorithmTypes.h - algorithm base types      ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#include "core/AlgorithmTypes.h"
#endif

#include <iostream>

using namespace MML;

void Docs_Demo_AlgorithmTypes()
{
    std::cout << "=== AlgorithmTypes Demo ===" << std::endl;
    
    // Algorithm status codes
    std::cout << "Algorithm status codes:" << std::endl;
    std::cout << "  Success:              " << ToString(AlgorithmStatus::Success) << std::endl;
    std::cout << "  MaxIterationsExceeded:" << ToString(AlgorithmStatus::MaxIterationsExceeded) << std::endl;
    std::cout << "  NumericalInstability: " << ToString(AlgorithmStatus::NumericalInstability) << std::endl;
    std::cout << "  SingularMatrix:       " << ToString(AlgorithmStatus::SingularMatrix) << std::endl;
    
    // Base algorithm configuration (ConfigBase)
    ConfigBase config;
    config.tolerance = 1e-10;
    config.max_iterations = 100;
    config.verbose = false;
    
    std::cout << "\nBase algorithm config:" << std::endl;
    std::cout << "  Tolerance: " << config.tolerance << std::endl;
    std::cout << "  Max iterations: " << config.max_iterations << std::endl;
    
    // Base algorithm result (IterativeResultBase)
    IterativeResultBase result;
    result.status = AlgorithmStatus::Success;
    result.iterations_used = 15;
    result.achieved_tolerance = 1.5e-11;
    result.converged = true;
    
    std::cout << "\nAlgorithm result:" << std::endl;
    std::cout << "  Status: " << ToString(result.status) << std::endl;
    std::cout << "  Iterations: " << result.iterations_used << std::endl;
    std::cout << "  Final error: " << result.achieved_tolerance << std::endl;
    std::cout << "  Converged: " << (result.converged ? "Yes" : "No") << std::endl;
    
    std::cout << "=== Demo Complete ===" << std::endl;
}

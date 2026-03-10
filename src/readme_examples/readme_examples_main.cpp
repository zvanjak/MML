///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        readme_examples_main.cpp                                            ///
///  Description: Main entry point for README examples                                ///
///               Demonstrates code exactly as shown in README.md                     ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Vector/Vector.h"
#include "base/Matrix/Matrix.h"
#endif

using namespace MML;

// Function declarations - Flagship demonstration
void Readme_FundamentalTheorems();   // Gauss, Stokes, Green's theorems

// Function declarations - existing sections
void Readme_QuickStart();
void Readme_VectorsMatrices();
void Readme_LinearSystems();
void Readme_FunctionsInterpolation();
void Readme_NumericalCalculus();
void Readme_FieldOperations();
void Readme_ODESolvers();

// Function declarations - new sections
void Readme_Polynomials();
void Readme_RootFinding();
void Readme_CoordTransforms();
void Readme_ParametricCurves();
void Readme_PathIntegrals();
void Readme_FunctionAnalysis();
void Readme_DynamicalSystems();  // Crown Jewel!

int main()
{
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                        README EXAMPLES                        ****" << std::endl;
    std::cout << "****            (Code matches README.md exactly)                   ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    // FLAGSHIP: Fundamental Theorems of Vector Calculus
    Readme_FundamentalTheorems();

    // Original examples
    Readme_QuickStart();
    Readme_VectorsMatrices();
    Readme_LinearSystems();
    Readme_FunctionsInterpolation();
    Readme_NumericalCalculus();
    Readme_FieldOperations();
    Readme_ODESolvers();

    // New examples
    Readme_Polynomials();
    Readme_RootFinding();
    Readme_CoordTransforms();
    Readme_ParametricCurves();
    Readme_PathIntegrals();
    Readme_FunctionAnalysis();

    // Crown Jewel - Dynamical Systems
    Readme_DynamicalSystems();

    return 0;
}

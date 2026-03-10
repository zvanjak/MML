///////////////////////////////////////////////////////////////////////////////////////////
///  File:        docs_demo_field_analyzers.cpp                                       ///
///  Description: Brief demonstration of FieldAnalyzers.h - field analysis tools      ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#include "base/Function.h"
#include "algorithms/FieldAnalyzers.h"
#endif

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace MML;

void Docs_Demo_FieldAnalyzers()
{
    std::cout << "=== FieldAnalyzers Demo ===" << std::endl;
    std::cout << std::fixed << std::setprecision(4);
    
    // Define a scalar field: f(x,y,z) = x² + y² + z² (distance squared)
    ScalarFunctionFromStdFunc<3> distSqField([](const VectorN<Real,3>& p) {
        return p[0]*p[0] + p[1]*p[1] + p[2]*p[2];
    });
    
    ScalarFieldAnalyzer scalarAnalyzer(distSqField);
    
    Vec3Cart point(1.0, 2.0, 3.0);
    std::cout << "Scalar field f = x² + y² + z² at (1,2,3):" << std::endl;
    
    Vec3Cart grad = scalarAnalyzer.GradientAt(point);
    std::cout << "  Gradient: (" << grad.X() << ", " << grad.Y() << ", " << grad.Z() << ")" << std::endl;
    std::cout << "  Expected: (2, 4, 6)" << std::endl;
    
    Real laplacian = scalarAnalyzer.LaplacianAt(point);
    std::cout << "  Laplacian: " << laplacian << " (expected: 6)" << std::endl;
    
    // Check for critical points (where gradient = 0)
    Vec3Cart origin(0.0, 0.0, 0.0);
    std::cout << "  Critical at origin: " << (scalarAnalyzer.IsCriticalPoint(origin) ? "Yes" : "No") << std::endl;
    
    // Define a vector field: F(x,y,z) = (y, -x, 0) (rotation around z-axis)
    VectorFunctionFromStdFunc<3> rotationField([](const VectorN<Real,3>& p) {
        return VectorN<Real,3>{p[1], -p[0], 0.0};
    });
    
    VectorFieldAnalyzer vectorAnalyzer(rotationField);
    
    std::cout << "\nVector field F = (y, -x, 0) at (1,2,3):" << std::endl;
    Real div = vectorAnalyzer.DivergenceAt(point);
    std::cout << "  Divergence: " << div << " (expected: 0)" << std::endl;
    
    Vec3Cart curl = vectorAnalyzer.CurlAt(point);
    std::cout << "  Curl: (" << curl.X() << ", " << curl.Y() << ", " << curl.Z() << ")" << std::endl;
    std::cout << "  Expected: (0, 0, -2)" << std::endl;
    
    std::cout << "=== Demo Complete ===" << std::endl;
}

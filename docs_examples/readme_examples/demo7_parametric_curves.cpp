#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "core/Function.h"

#include "core/Curves.h"

#include "algorithms/DiffGeometryAlgorithms.h"
#endif

#include "../test_data/parametric_curves_test_bed.h"

using namespace MML;

void Readme_parametric_curves()
{
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                  README - parametric curves                   ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    // creating curve directly with lambda
    ParametricCurve<3>        test_curve1( [](double t) -> VectorN<Real, 3> { return VectorN<Real, 3>{t, t*t, t*t*t}; } );
    
    // using predefined curve
    Curves::HelixCurve        helix(2.0, 2.0);
    
    // using curve from TestData
    const ParametricCurve<3> &test_curve = TestBeds::ParametricCurvesTestBed::getTestCurve(0)._curve;

    double t = 0.5;
    auto tangent   = Vector3Cartesian( DiffGeometry::getTangent(test_curve, t) );
    auto unit_tang = Vector3Cartesian( DiffGeometry::getTangentUnit(test_curve, t) );
    auto normal    = Vector3Cartesian( DiffGeometry::getNormal(test_curve, t) );
    auto unit_norm = Vector3Cartesian( DiffGeometry::getNormalUnit(test_curve, t) );
    auto binormal  = VectorProd(unit_tang, unit_norm);
    
    auto curv_vec   = DiffGeometry::getCurvatureVector(test_curve, t);
    auto curvature  = DiffGeometry::getCurvature(test_curve, t);
    auto curvature3 = DiffGeometry::getCurvature3(test_curve, t);
}
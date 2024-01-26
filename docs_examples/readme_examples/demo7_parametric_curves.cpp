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
    Curves::LemniscateCurve   lemniscate;
    Curves::ToroidalSpiralCurve torus(3, 2.0);
    
    // using curve from TestData
    const ParametricCurve<3> &test_curve = TestBeds::ParametricCurvesTestBed::getTestCurve("Helix")._curve;

    double t = 0.5;
    auto tangent   = DiffGeometry::getTangent(test_curve, t);
    auto unit_tang = DiffGeometry::getTangentUnit(test_curve, t);
    auto normal    = DiffGeometry::getNormal(test_curve, t);
    auto unit_norm = DiffGeometry::getNormalUnit(test_curve, t);
    auto binormal  = VectorProd(Vector3Cartesian(unit_tang), Vector3Cartesian(unit_norm));
    
    auto curv_vec   = DiffGeometry::getCurvatureVector(test_curve, t);
    auto curvature  = DiffGeometry::getCurvature(test_curve, t);
    auto curvature3 = DiffGeometry::getCurvature3(test_curve, t);
    // TODO 0.8 - vizualizirati neku krivulju, i u jednoj točki vizualizirati (World view) 3 vektora tangente, normale i binormale, te vektore zakrivljenosti
}
#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "core/Function.h"
#include "core/Curves.h"

#include "algorithms/ParametricCurveAnalyzer.h"
#endif

#include "../test_data/parametric_curves_test_bed.h"

using namespace MML;

void Readme_parametric_curves()
{
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                  README - parametric curves                   ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    // creating curve directly with lambda
    ParametricCurve<3>        test_curve1( [](Real t) -> VectorN<Real, 3> { return VectorN<Real, 3>{t, t*t, t*t*t}; } );
    
    // using predefined curve
    Curves2D::LemniscateCurve     lemniscate;
    Curves3D::ToroidalSpiralCurve torus(3, 2.0);
    
    // using curve from TestData
    const ParametricCurve<3> &test_curve = TestBeds::ParametricCurvesTestBed::getTestCurve("Helix")._curve;

    // interpolated curve

    // as a result of OE solution

    double t = 0.5;
    auto tangent   = ParametricCurveAnalyzer::getTangent(test_curve, t);
    auto unit_tang = ParametricCurveAnalyzer::getTangentUnit(test_curve, t);
    auto normal    = ParametricCurveAnalyzer::getNormal(test_curve, t);
    auto unit_norm = ParametricCurveAnalyzer::getNormalUnit(test_curve, t);
    auto binormal  = VectorProd(Vector3Cartesian(unit_tang), Vector3Cartesian(unit_norm));
    
    auto curv_vec   = ParametricCurveAnalyzer::getCurvatureVector(test_curve, t);
    auto curvature  = ParametricCurveAnalyzer::getCurvature(test_curve, t);
    auto curvature3 = ParametricCurveAnalyzer::getCurvature3(test_curve, t);
    
    // TODO 1.0 - vizualizirati neku krivulju, i u jednoj točki vizualizirati (World view) 3 vektora tangente, normale i binormale, te vektore zakrivljenosti
}
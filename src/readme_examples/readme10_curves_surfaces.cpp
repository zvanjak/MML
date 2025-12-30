#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Function.h"
#include "core/Curves.h"
#endif

#include "../test_data/parametric_curves_test_bed.h"

using namespace MML;

void Readme_curves_surfaces()
{
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                  README - parametric curves                   ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    // creating curve directly with lambda
    Curves::CurveCartesian3D test_curve1([](Real t) -> VectorN<Real, 3> { return VectorN<Real, 3>{t, t*t, t*t*t}; });
    
    // using predefined curve
    Curves::LemniscateCurve     lemniscate;
    Curves::ToroidalSpiralCurve torus(3, 2.0);
    
    // using curve from TestData
    const Curves::CurveCartesian3D &test_curve = TestBeds::ParametricCurvesTestBed::getTestCurve("Helix")._curve;

    // interpolated curve

    // as a result of OE solution

    double t = 0.5;
    auto tangent   = test_curve1.getTangent(t);
    auto unit_tang = test_curve1.getTangentUnit(t);
    auto normal    = test_curve1.getNormal(t);
    auto unit_norm = test_curve1.getNormalUnit(t);
    auto binormal  = test_curve1.getBinormal(t);
    
    auto curv_vec   = test_curve1.getCurvatureVector(t);
    auto curvature  = test_curve1.getCurvature(t);
    
    // TODO 1.0 - vizualizirati neku krivulju, i u jednoj točki vizualizirati (World view) 3 vektora tangente, normale i binormale, te vektore zakrivljenosti
}
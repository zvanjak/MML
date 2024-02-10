#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/VectorN.h"
#include "core/Function.h"
#include "core/Derivation.h"

#include "core/Curves.h"
#include "base/Geometry3D.h"

#include "algorithms/ParametricCurveAnalyzer.h"
#endif

#include "../test_data/parametric_curves_test_bed.h"
#include "../test_data/parametric_surfaces_test_bed.h"

using namespace MML;

void Demo_curves()
{
    std::cout << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                    DIFFERENTIAL GEOMETRY                      ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    ParametricCurve<3>        test_curve1( [](Real t) -> VectorN<Real, 3> { return VectorN<Real, 3>{t, t*t, t*t*t}; } );
    Curves3D::HelixCurve      helix(2.0, 2.0);
    const ParametricCurve<3> &test_curve = TestBeds::ParametricCurvesTestBed::getTestCurve(0)._curve;

    std::cout << "          Tangent                   Tangent unit                   Normal                  Normal unit                    Binormal                   Curv.vec.                Curv.vec.norm.        Curvature\n";

    std::cout << std::fixed;

    for(double t=0.0; t<2*Constants::PI; t+=0.4)
    {
        auto tangent   = Vector3Cartesian( ParametricCurveAnalyzer::getTangent(test_curve, t) );
        auto unit_tang = Vector3Cartesian( ParametricCurveAnalyzer::getTangentUnit(test_curve, t) );
        auto normal    = Vector3Cartesian( ParametricCurveAnalyzer::getNormal(test_curve, t) );
        auto unit_norm = Vector3Cartesian( ParametricCurveAnalyzer::getNormalUnit(test_curve, t) );
        auto binormal  = VectorProd(unit_tang, unit_norm);
        
        auto curv_vec   = ParametricCurveAnalyzer::getCurvatureVector(test_curve, t);
        auto curvature  = ParametricCurveAnalyzer::getCurvature(test_curve, t);
        auto curvature3 = ParametricCurveAnalyzer::getCurvature3(test_curve, t);

        tangent.Print(std::cout, 7, 3); std::cout << " ";
        unit_tang.Print(std::cout, 7, 3); std::cout << " ";
        normal.Print(std::cout, 7, 3); std::cout << " ";
        unit_norm.Print(std::cout, 7, 3); std::cout << " ";
        binormal.Print(std::cout, 7, 3); std::cout << " ";
        
        curv_vec.Print(std::cout, 7, 3); std::cout << " ";
        auto curv_vec_norm = curv_vec / curv_vec.NormL2();
        curv_vec_norm.Print(std::cout, 7, 3);
        std::cout << "   " << curvature << "   " << std::endl;
    }

    bool b = ParametricCurveAnalyzer::isArcLengthParametrized(test_curve, 0.0, 2*Constants::PI);
    std::cout << "Is arc length parametrized : " << b << std::endl;
}

void Demo_surfaces()
{
    //const ParametricCurve<3> &test_curve = MML::TestData::ParametricCurvesTestBed::getTestCurve(0)._curve;

}

void Demo_Diff_geometry()
{
    Demo_curves();
    Demo_surfaces();
}
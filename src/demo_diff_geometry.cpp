#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include <iostream>
#include <iomanip>
#include <cmath>

#include "core/Constants.h"
#include "core/VectorN.h"

#include "basic_types/Function.h"
#include "basic_types/Curves.h"
#include "basic_types/Geometry3D.h"

#include "algorithms/DiffGeometryAlgorithms.h"
#include "algorithms/Derivation.h"
#endif

#include "../test_data/parametric_curves_test_bed.h"

using namespace MML;

void Demo_Diff_geometry()
{
    std::cout << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                    DIFFERENTIAL GEOMETRY                      ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    ParametricCurve<3>   test_curve1( [](double t) -> VectorN<Real, 3> { return VectorN<Real, 3>{t, t*t, t*t*t}; } );
    
    MML::Curves::HelixCurve     helix(2.0, 2.0);
    //IParametricCurve<3>  &test_curve = helix;
    const ParametricCurve<3>  &test_curve = MML::TestData::ParametricCurvesTestBed::_listCurves[0]._curve;

    std::cout << "          Tangent                   Tangent unit                   Normal                  Normal unit                    Binormal                   Curv.vec.                Curv.vec.norm.        Curvature\n";

    std::cout << std::fixed;

    for(double t=0.0; t<2*Constants::PI; t+=0.4)
    {
        auto tangent   = Vector3Cartesian( DiffGeometry::getTangent(test_curve, t) );
        auto unit_tang = Vector3Cartesian( DiffGeometry::getTangentUnit(test_curve, t) );
        auto normal    = Vector3Cartesian( DiffGeometry::getNormal(test_curve, t) );
        auto unit_norm = Vector3Cartesian( DiffGeometry::getNormalUnit(test_curve, t) );
        auto binormal  = VectorProd(unit_tang, unit_norm);
        
        auto curv_vec   = DiffGeometry::getCurvatureVector(test_curve, t);
        auto curvature  = DiffGeometry::getCurvature(test_curve, t);
        auto curvature3 = DiffGeometry::getCurvature3(test_curve, t);

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

    bool b = DiffGeometry::isArcLengthParametrized(test_curve, 0.0, 2*Constants::PI);
    std::cout << "Is arc length parametrized : " << b << std::endl;
}
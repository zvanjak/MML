#include "../catch/catch.hpp"

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

TEST_CASE("Test_Helix_Curvature", "[simple]") 
{
    const ParametricCurve<3>  &schaum_curve = MML::TestData::ParametricCurvesTestBed::_listCurves[2]._curve;
    const RealFunction  &schaum_curve_curv = MML::TestData::ParametricCurvesTestBed::_listCurves[2]._curvatureFunc;

    double curv = DiffGeometry::getCurvature(schaum_curve, 0.5);
    REQUIRE( schaum_curve_curv(0.5) == Approx(curv) );
}
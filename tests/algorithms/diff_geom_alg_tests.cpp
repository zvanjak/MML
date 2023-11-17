#include "../catch/catch.hpp"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include <iostream>
#include <iomanip>
#include <cmath>

#include "utilities/Constants.h"
#include "core/VectorN.h"

#include "core/Function.h"
#include "basic_types/Curves.h"
#include "basic_types/Geometry3D.h"

#include "algorithms/DiffGeometryAlgorithms.h"
#include "core/Derivation.h"
#endif

#include "../test_data/parametric_curves_test_bed.h"

using namespace MML;

// TODO - BIG!!! - implement these tests
// za krivulje iz test beda u točkama verificirati preciznost izračuna tangente, normale, binormale, vektora krivine, krivine
TEST_CASE("Test_Helix_Curvature", "[simple]") 
{
    const ParametricCurve<3>  &schaum_curve = TestBeds::ParametricCurvesTestBed::getTestCurve(2)._curve;
    const RealFunction  &schaum_curve_curv = TestBeds::ParametricCurvesTestBed::getTestCurve(2)._curvatureFunc;

    double curv = DiffGeometry::getCurvature(schaum_curve, 0.5);
    REQUIRE( schaum_curve_curv(0.5) == Approx(curv) );
}
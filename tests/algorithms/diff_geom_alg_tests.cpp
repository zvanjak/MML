#include "../catch/catch.hpp"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/VectorN.h"
#include "base/Geometry3D.h"

#include "core/Function.h"
#include "core/Curves.h"
#include "core/Derivation.h"

#include "algorithms/ParametricCurveAnalyzer.h"
#endif

#include "../test_data/parametric_curves_test_bed.h"

using namespace MML;

// TODO - HIGH, EMPTY!!! - implement these tests
// za krivulje iz test beda u točkama verificirati preciznost izračuna tangente, normale, binormale, vektora krivine, krivine
// za dvije varijante Helixa, arc-param i 3x, provesti testove i usporediti s točnim teorijskim vrijednostima
TEST_CASE("Test_Helix_getTangent", "[simple]") 
{

}

TEST_CASE("Test_Helix_getTangentUnit", "[simple]") 
{

}
TEST_CASE("Test_Helix_getNormal", "[simple]") 
{

}
TEST_CASE("Test_Helix_getNormalUnit", "[simple]") 
{

}
TEST_CASE("Test_Helix_getPrincipalNormal", "[simple]") 
{

}
TEST_CASE("Test_Helix_getBinormal", "[simple]") 
{

}
TEST_CASE("Test_Helix_getCurvatureVector", "[simple]") 
{

}

TEST_CASE("Test_Helix_Curvature", "[simple]") 
{
    const ParametricCurve<3>  &schaum_curve = TestBeds::ParametricCurvesTestBed::getTestCurve(2)._curve;
    const RealFunction  &schaum_curve_curv = TestBeds::ParametricCurvesTestBed::getTestCurve(2)._curvatureFunc;

    double curv = ParametricCurveAnalyzer::getCurvature(schaum_curve, 0.5);
    REQUIRE( schaum_curve_curv(0.5) == Approx(curv) );
}


TEST_CASE("Test_Helix_getCurvature3", "[simple]") 
{

}
TEST_CASE("Test_Helix_getTorsion3", "[simple]") 
{

}
TEST_CASE("Test_Helix_getOsculationPlane", "[simple]") 
{

}
TEST_CASE("Test_Helix_getNormalPlane", "[simple]") 
{

}
TEST_CASE("Test_Helix_getRectifyingPlane", "[simple]") 
{

}
TEST_CASE("Test_Helix_isArcLengthParametrized", "[simple]") 
{

}

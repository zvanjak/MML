#include "../../test_data/parametric_curves_test_bed.h"
#include "../../mml/core/Derivation/DerivationParametricCurve.h"

#include <string>

#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

using namespace MML;
using namespace MML::Testing;
using namespace MML::Derivation;
using namespace MML::TestBeds;

/********************************************************************************************************************/
/**                         Test curves with known analytical derivatives                                         **/
/********************************************************************************************************************/

// Simple line curve: r(t) = {t, 2t, 3t}
// r'(t) = {1, 2, 3}
// r''(t) = {0, 0, 0}
class TestLineCurve : public IParametricCurve<3>
{
public:
    Real getMinT() const override { return -REAL(10.0); }
    Real getMaxT() const override { return REAL(10.0); }
    
    VectorN<Real, 3> operator()(Real t) const override {
        return VectorN<Real, 3>{ t, 2*t, 3*t };
    }
    
    VectorN<Real, 3> analytical_first(Real t) const {
        return VectorN<Real, 3>{ REAL(1.0), REAL(2.0), REAL(3.0) };
    }
    
    VectorN<Real, 3> analytical_second(Real t) const {
        return VectorN<Real, 3>{ REAL(0.0), REAL(0.0), REAL(0.0) };
    }
    
    VectorN<Real, 3> analytical_third(Real t) const {
        return VectorN<Real, 3>{ REAL(0.0), REAL(0.0), REAL(0.0) };
    }
};

// Parabolic curve: r(t) = {t, t², t³}
// r'(t) = {1, 2t, 3t²}
// r''(t) = {0, 2, 6t}
// r'''(t) = {0, 0, 6}
class TestParabolicCurve : public IParametricCurve<3>
{
public:
    Real getMinT() const override { return -REAL(10.0); }
    Real getMaxT() const override { return REAL(10.0); }
    
    VectorN<Real, 3> operator()(Real t) const override {
        return VectorN<Real, 3>{ t, t*t, t*t*t };
    }
    
    VectorN<Real, 3> analytical_first(Real t) const {
        return VectorN<Real, 3>{ REAL(1.0), 2*t, 3*t*t };
    }
    
    VectorN<Real, 3> analytical_second(Real t) const {
        return VectorN<Real, 3>{ REAL(0.0), REAL(2.0), 6*t };
    }
    
    VectorN<Real, 3> analytical_third(Real t) const {
        return VectorN<Real, 3>{ REAL(0.0), REAL(0.0), REAL(6.0) };
    }
};

// Circle curve: r(t) = {cos(t), sin(t), 0}
// r'(t) = {-sin(t), cos(t), 0}
// r''(t) = {-cos(t), -sin(t), 0}
// r'''(t) = {sin(t), -cos(t), 0}
class TestCircleCurve : public IParametricCurve<3>
{
public:
    Real getMinT() const override { return REAL(0.0); }
    Real getMaxT() const override { return REAL(2.0) * Constants::PI; }
    
    VectorN<Real, 3> operator()(Real t) const override {
        return VectorN<Real, 3>{ cos(t), sin(t), REAL(0.0) };
    }
    
    VectorN<Real, 3> analytical_first(Real t) const {
        return VectorN<Real, 3>{ -sin(t), cos(t), REAL(0.0) };
    }
    
    VectorN<Real, 3> analytical_second(Real t) const {
        return VectorN<Real, 3>{ -cos(t), -sin(t), REAL(0.0) };
    }
    
    VectorN<Real, 3> analytical_third(Real t) const {
        return VectorN<Real, 3>{ sin(t), -cos(t), REAL(0.0) };
    }
};

// Helix curve: r(t) = {cos(t), sin(t), t}
// r'(t) = {-sin(t), cos(t), 1}
// r''(t) = {-cos(t), -sin(t), 0}
// r'''(t) = {sin(t), -cos(t), 0}
class TestHelixCurve : public IParametricCurve<3>
{
public:
    Real getMinT() const override { return REAL(0.0); }
    Real getMaxT() const override { return REAL(4.0) * Constants::PI; }
    
    VectorN<Real, 3> operator()(Real t) const override {
        return VectorN<Real, 3>{ cos(t), sin(t), t };
    }
    
    VectorN<Real, 3> analytical_first(Real t) const {
        return VectorN<Real, 3>{ -sin(t), cos(t), REAL(1.0) };
    }
    
    VectorN<Real, 3> analytical_second(Real t) const {
        return VectorN<Real, 3>{ -cos(t), -sin(t), REAL(0.0) };
    }
    
    VectorN<Real, 3> analytical_third(Real t) const {
        return VectorN<Real, 3>{ sin(t), -cos(t), REAL(0.0) };
    }
};

/********************************************************************************************************************/
/**                                         FIRST DERIVATIVE TESTS                                                **/
/********************************************************************************************************************/

TEST_CASE("Parametric curve - NDer1 (first order forward difference)", "[parametric_curve][derivative][NDer1]")
{
    TestLineCurve line;
    TestParabolicCurve para;
    TestCircleCurve circle;
    TestHelixCurve helix;
    
    Real t = REAL(1.0);
    Real tol = 1e-4;  // NDer1 is less accurate
    
    // Line curve
    VectorN<Real, 3> num_deriv = NDer1(line, t);
    VectorN<Real, 3> analytical = line.analytical_first(t);
    for(int i = 0; i < 3; i++)
        REQUIRE(std::abs(num_deriv[i] - analytical[i]) < tol);
    
    // Parabolic curve
    num_deriv = NDer1(para, t);
    analytical = para.analytical_first(t);
    for(int i = 0; i < 3; i++)
        REQUIRE(std::abs(num_deriv[i] - analytical[i]) < tol);
    
    // Circle
    num_deriv = NDer1(circle, t);
    analytical = circle.analytical_first(t);
    for(int i = 0; i < 3; i++)
        REQUIRE(std::abs(num_deriv[i] - analytical[i]) < tol);
    
    // Helix
    num_deriv = NDer1(helix, t);
    analytical = helix.analytical_first(t);
    for(int i = 0; i < 3; i++)
        REQUIRE(std::abs(num_deriv[i] - analytical[i]) < tol);
}

TEST_CASE("Parametric curve - NDer2 (second order central difference)", "[parametric_curve][derivative][NDer2]")
{
    TestLineCurve line;
    TestParabolicCurve para;
    TestCircleCurve circle;
    TestHelixCurve helix;
    
    Real t = REAL(1.0);
    Real tol = 1e-7;
    
    // Line curve - exact for linear functions
    VectorN<Real, 3> num_deriv = NDer2(line, t);
    VectorN<Real, 3> analytical = line.analytical_first(t);
    for(int i = 0; i < 3; i++)
        REQUIRE(std::abs(num_deriv[i] - analytical[i]) < tol);
    
    // Parabolic curve
    num_deriv = NDer2(para, t);
    analytical = para.analytical_first(t);
    for(int i = 0; i < 3; i++)
        REQUIRE(std::abs(num_deriv[i] - analytical[i]) < tol);
    
    // Circle
    num_deriv = NDer2(circle, t);
    analytical = circle.analytical_first(t);
    for(int i = 0; i < 3; i++)
        REQUIRE(std::abs(num_deriv[i] - analytical[i]) < tol);
    
    // Helix
    num_deriv = NDer2(helix, t);
    analytical = helix.analytical_first(t);
    for(int i = 0; i < 3; i++)
        REQUIRE(std::abs(num_deriv[i] - analytical[i]) < tol);
}

TEST_CASE("Parametric curve - NDer4 (fourth order)", "[parametric_curve][derivative][NDer4]")
{
    TestParabolicCurve para;
    TestCircleCurve circle;
    TestHelixCurve helix;
    
    Real t = REAL(1.0);
    Real tol = 1e-9;
    
    // Parabolic curve
    VectorN<Real, 3> num_deriv = NDer4(para, t);
    VectorN<Real, 3> analytical = para.analytical_first(t);
    for(int i = 0; i < 3; i++)
        REQUIRE(std::abs(num_deriv[i] - analytical[i]) < tol);
    
    // Circle
    num_deriv = NDer4(circle, t);
    analytical = circle.analytical_first(t);
    for(int i = 0; i < 3; i++)
        REQUIRE(std::abs(num_deriv[i] - analytical[i]) < tol);
    
    // Helix
    num_deriv = NDer4(helix, t);
    analytical = helix.analytical_first(t);
    for(int i = 0; i < 3; i++)
        REQUIRE(std::abs(num_deriv[i] - analytical[i]) < tol);
}

TEST_CASE("Parametric curve - NDer6 (sixth order)", "[parametric_curve][derivative][NDer6]")
{
    TestParabolicCurve para;
    TestCircleCurve circle;
    TestHelixCurve helix;
    
    Real t = REAL(1.0);
    Real tol = 1e-10;
    
    // Parabolic curve
    VectorN<Real, 3> num_deriv = NDer6(para, t);
    VectorN<Real, 3> analytical = para.analytical_first(t);
    for(int i = 0; i < 3; i++)
        REQUIRE(std::abs(num_deriv[i] - analytical[i]) < tol);
    
    // Circle
    num_deriv = NDer6(circle, t);
    analytical = circle.analytical_first(t);
    for(int i = 0; i < 3; i++)
        REQUIRE(std::abs(num_deriv[i] - analytical[i]) < tol);
    
    // Helix
    num_deriv = NDer6(helix, t);
    analytical = helix.analytical_first(t);
    for(int i = 0; i < 3; i++)
        REQUIRE(std::abs(num_deriv[i] - analytical[i]) < tol);
}

TEST_CASE("Parametric curve - NDer8 (eighth order)", "[parametric_curve][derivative][NDer8]")
{
    TestParabolicCurve para;
    TestCircleCurve circle;
    TestHelixCurve helix;
    
    Real t = REAL(1.0);
    Real tol = 1e-10;
    
    // Parabolic curve
    VectorN<Real, 3> num_deriv = NDer8(para, t);
    VectorN<Real, 3> analytical = para.analytical_first(t);
    for(int i = 0; i < 3; i++)
        REQUIRE(std::abs(num_deriv[i] - analytical[i]) < tol);
    
    // Circle
    num_deriv = NDer8(circle, t);
    analytical = circle.analytical_first(t);
    for(int i = 0; i < 3; i++)
        REQUIRE(std::abs(num_deriv[i] - analytical[i]) < tol);
    
    // Helix
    num_deriv = NDer8(helix, t);
    analytical = helix.analytical_first(t);
    for(int i = 0; i < 3; i++)
        REQUIRE(std::abs(num_deriv[i] - analytical[i]) < tol);
}

/********************************************************************************************************************/
/**                                         SECOND DERIVATIVE TESTS                                               **/
/********************************************************************************************************************/

TEST_CASE("Parametric curve - NSecDer2 (second order second derivative)", "[parametric_curve][derivative][NSecDer2]")
{
    TestLineCurve line;
    TestParabolicCurve para;
    TestCircleCurve circle;
    TestHelixCurve helix;
    
    Real t = REAL(1.0);
    Real tol_zero = 1e-5;
    Real tol = 1e-5;
    
    // Line curve - second derivative is zero
    VectorN<Real, 3> num_deriv = NSecDer2(line, t);
    VectorN<Real, 3> analytical = line.analytical_second(t);
    for(int i = 0; i < 3; i++)
        REQUIRE(std::abs(num_deriv[i] - analytical[i]) < tol_zero);
    
    // Parabolic curve
    num_deriv = NSecDer2(para, t);
    analytical = para.analytical_second(t);
    for(int i = 0; i < 3; i++)
        REQUIRE(std::abs(num_deriv[i] - analytical[i]) < tol);
    
    // Circle
    num_deriv = NSecDer2(circle, t);
    analytical = circle.analytical_second(t);
    for(int i = 0; i < 3; i++)
        REQUIRE(std::abs(num_deriv[i] - analytical[i]) < tol);
    
    // Helix
    num_deriv = NSecDer2(helix, t);
    analytical = helix.analytical_second(t);
    for(int i = 0; i < 3; i++)
        REQUIRE(std::abs(num_deriv[i] - analytical[i]) < tol);
}

TEST_CASE("Parametric curve - NSecDer4 (fourth order second derivative)", "[parametric_curve][derivative][NSecDer4]")
{
    TestParabolicCurve para;
    TestCircleCurve circle;
    TestHelixCurve helix;
    
    Real t = REAL(1.0);
    Real tol = 1e-7;
    
    // Parabolic curve
    VectorN<Real, 3> num_deriv = NSecDer4(para, t);
    VectorN<Real, 3> analytical = para.analytical_second(t);
    for(int i = 0; i < 3; i++)
        REQUIRE(std::abs(num_deriv[i] - analytical[i]) < tol);
    
    // Circle
    num_deriv = NSecDer4(circle, t);
    analytical = circle.analytical_second(t);
    for(int i = 0; i < 3; i++)
        REQUIRE(std::abs(num_deriv[i] - analytical[i]) < tol);
    
    // Helix
    num_deriv = NSecDer4(helix, t);
    analytical = helix.analytical_second(t);
    for(int i = 0; i < 3; i++)
        REQUIRE(std::abs(num_deriv[i] - analytical[i]) < tol);
}

/********************************************************************************************************************/
/**                                         THIRD DERIVATIVE TESTS                                                **/
/********************************************************************************************************************/

TEST_CASE("Parametric curve - NThirdDer2 (second order third derivative)", "[parametric_curve][derivative][NThirdDer2]")
{
    TestParabolicCurve para;
    TestCircleCurve circle;
    TestHelixCurve helix;
    
    Real t = REAL(1.0);
    Real tol = 1e-4;
    
    // Parabolic curve
    VectorN<Real, 3> num_deriv = NThirdDer2(para, t);
    VectorN<Real, 3> analytical = para.analytical_third(t);
    for(int i = 0; i < 3; i++)
        REQUIRE(std::abs(num_deriv[i] - analytical[i]) < tol);
    
    // Circle
    num_deriv = NThirdDer2(circle, t);
    analytical = circle.analytical_third(t);
    for(int i = 0; i < 3; i++)
        REQUIRE(std::abs(num_deriv[i] - analytical[i]) < tol);
    
    // Helix
    num_deriv = NThirdDer2(helix, t);
    analytical = helix.analytical_third(t);
    for(int i = 0; i < 3; i++)
        REQUIRE(std::abs(num_deriv[i] - analytical[i]) < tol);
}

TEST_CASE("Parametric curve - NThirdDer4 (fourth order third derivative)", "[parametric_curve][derivative][NThirdDer4]")
{
    TestParabolicCurve para;
    TestCircleCurve circle;
    TestHelixCurve helix;
    
    Real t = REAL(1.0);
    Real tol = 1e-5;
    
    // Parabolic curve
    VectorN<Real, 3> num_deriv = NThirdDer4(para, t);
    VectorN<Real, 3> analytical = para.analytical_third(t);
    for(int i = 0; i < 3; i++)
        REQUIRE(std::abs(num_deriv[i] - analytical[i]) < tol);
    
    // Circle
    num_deriv = NThirdDer4(circle, t);
    analytical = circle.analytical_third(t);
    for(int i = 0; i < 3; i++)
        REQUIRE(std::abs(num_deriv[i] - analytical[i]) < tol);
    
    // Helix
    num_deriv = NThirdDer4(helix, t);
    analytical = helix.analytical_third(t);
    for(int i = 0; i < 3; i++)
        REQUIRE(std::abs(num_deriv[i] - analytical[i]) < tol);
}

/********************************************************************************************************************/
/**                                    ACCURACY PROGRESSION TESTS                                                 **/
/********************************************************************************************************************/

TEST_CASE("Parametric curve - First derivative accuracy progression", "[parametric_curve][derivative][accuracy]")
{
    TestHelixCurve helix;
    Real t = Constants::PI / 4;
    VectorN<Real, 3> analytical = helix.analytical_first(t);
    
    // Compute errors for each order
    Real error1 = (NDer1(helix, t) - analytical).NormL2();
    Real error2 = (NDer2(helix, t) - analytical).NormL2();
    Real error4 = (NDer4(helix, t) - analytical).NormL2();
    Real error6 = (NDer6(helix, t) - analytical).NormL2();
    Real error8 = (NDer8(helix, t) - analytical).NormL2();
    
    // Verify progressive improvement (up to roundoff)
    REQUIRE(error1 > error2);
    REQUIRE(error2 > error4);
    REQUIRE(error4 > error6);
    // NDer8 may not beat NDer6 due to roundoff
    REQUIRE(error8 < 1e-9);
}

TEST_CASE("Parametric curve - Second derivative accuracy progression", "[parametric_curve][derivative][accuracy]")
{
    TestHelixCurve helix;
    Real t = Constants::PI / 4;
    VectorN<Real, 3> analytical = helix.analytical_second(t);
    
    // Compute errors for each order - now using direct formulas (NSecDer2 and NSecDer4 only)
    Real error2 = (NSecDer2(helix, t) - analytical).NormL2();
    Real error4 = (NSecDer4(helix, t) - analytical).NormL2();
    
    // Verify NSecDer4 is more accurate than NSecDer2
    REQUIRE(error2 > error4);
    REQUIRE(error4 < 1e-7);
}

TEST_CASE("Parametric curve - Third derivative accuracy progression", "[parametric_curve][derivative][accuracy]")
{
    TestHelixCurve helix;
    Real t = Constants::PI / 4;
    VectorN<Real, 3> analytical = helix.analytical_third(t);
    
    // Compute errors for each order - now using direct formulas (NThirdDer2 and NThirdDer4 only)
    Real error2 = (NThirdDer2(helix, t) - analytical).NormL2();
    Real error4 = (NThirdDer4(helix, t) - analytical).NormL2();
    
    // Verify NThirdDer4 is more accurate than NThirdDer2
    REQUIRE(error2 > error4);
    REQUIRE(error4 < 1e-5);
}

/********************************************************************************************************************/
/**                                    MULTIPLE TEST POINTS                                                       **/
/********************************************************************************************************************/

TEST_CASE("Parametric curve - Multiple test points for robustness", "[parametric_curve][derivative][robustness]")
{
    TestHelixCurve helix;
    Real tol = 1e-9;
    
    std::vector<Real> test_points = {REAL(0.0), Constants::PI/6, Constants::PI/4, Constants::PI/3, Constants::PI/2, Constants::PI};
    
    for(Real t : test_points)
    {
        // First derivative
        VectorN<Real, 3> num_deriv = NDer4(helix, t);
        VectorN<Real, 3> analytical = helix.analytical_first(t);
        for(int i = 0; i < 3; i++)
            REQUIRE(std::abs(num_deriv[i] - analytical[i]) < tol);
        
        // Second derivative
        num_deriv = NSecDer4(helix, t);
        analytical = helix.analytical_second(t);
        for(int i = 0; i < 3; i++)
            REQUIRE(std::abs(num_deriv[i] - analytical[i]) < 1e-7);
    }
}

/********************************************************************************************************************/
/**                                    TEST BED CURVES VERIFICATION                                               **/
/********************************************************************************************************************/

TEST_CASE("Parametric curve - Test bed curves first derivatives", "[parametric_curve][derivative][testbed]")
{
    Real tol = 1e-8;
    
    // Test with Helix from test bed
    auto helix = ParametricCurvesTestBed::getTestCurve("Helix");
    Real t = Constants::PI / 4;
    
    VectorN<Real, 3> num_deriv = NDer4(helix._curve, t);
    VectorN<Real, 3> analytical = helix._curveDerived(t);
    
    for(int i = 0; i < 3; i++)
        REQUIRE(std::abs(num_deriv[i] - analytical[i]) < tol);
    
    // Test with Circle3DXY from test bed
    auto circle = ParametricCurvesTestBed::getTestCurve("Circle3DXY");
    t = Constants::PI / 3;
    
    num_deriv = NDer4(circle._curve, t);
    analytical = circle._curveDerived(t);
    
    for(int i = 0; i < 3; i++)
        REQUIRE(std::abs(num_deriv[i] - analytical[i]) < tol);
}

TEST_CASE("Parametric curve - Test bed curves second derivatives", "[parametric_curve][derivative][testbed]")
{
    Real tol = 1e-6;
    
    // Test with Helix from test bed
    auto helix = ParametricCurvesTestBed::getTestCurve("Helix");
    Real t = Constants::PI / 4;
    
    VectorN<Real, 3> num_deriv = NSecDer4(helix._curve, t);
    VectorN<Real, 3> analytical = helix._curveDerSecond(t);
    
    for(int i = 0; i < 3; i++)
        REQUIRE(std::abs(num_deriv[i] - analytical[i]) < tol);
    
    // Test with TwistedCubic from test bed
    // auto twisted = ParametricCurvesTestBed::getTestCurve("TwistedCubic");
    // t = REAL(0.5);
    
    // num_deriv = NSecDer4(twisted._curve, t);
    // analytical = twisted._curveDerSecond(t);
    
    // for(int i = 0; i < 3; i++)
    //     REQUIRE(std::abs(num_deriv[i] - analytical[i]) < tol);
}

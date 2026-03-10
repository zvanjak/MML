/**
 * @file singularity_handling_tests.cpp
 * @brief Comprehensive tests for SingularityHandling.h
 * 
 * Tests cover:
 * - SingularityPolicy enum values
 * - SingularityPolicyToString conversion
 * - Detection functions (IsNearZero, IsAtSphericalOrigin, IsAtSphericalPole, 
 *   IsAtCylindricalAxis, IsAtSphericalSingularity, DescribeSphericalSingularity)
 * - Safe arithmetic operations (SafeDivide, SafeInverseR, SafeInverseR2,
 *   SafeInverseRSinTheta, SafeInverseR2Sin2Theta, SafeCotThetaOverR2)
 */

#include <catch2/catch_all.hpp>
#include <cmath>
#include <limits>
#include <string>

#include "core/SingularityHandling.h"

using namespace MML;
using namespace MML::Singularity;

namespace MML::Tests::Core::SingularityHandlingTests {

///////////////////////////////////////////////////////////////////////////////
//                        SINGULARITY POLICY ENUM                            //
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("SingularityPolicy - Enum values exist", "[Singularity][Policy]") {
    // Verify all policy values are defined
    REQUIRE(static_cast<int>(SingularityPolicy::Throw) >= 0);
    REQUIRE(static_cast<int>(SingularityPolicy::ReturnNaN) >= 0);
    REQUIRE(static_cast<int>(SingularityPolicy::ReturnInf) >= 0);
    REQUIRE(static_cast<int>(SingularityPolicy::Clamp) >= 0);
    REQUIRE(static_cast<int>(SingularityPolicy::ReturnZero) >= 0);
}

TEST_CASE("SingularityPolicy - Enum values are distinct", "[Singularity][Policy]") {
    // All policies should have different values
    REQUIRE(static_cast<int>(SingularityPolicy::Throw) != static_cast<int>(SingularityPolicy::ReturnNaN));
    REQUIRE(static_cast<int>(SingularityPolicy::Throw) != static_cast<int>(SingularityPolicy::ReturnInf));
    REQUIRE(static_cast<int>(SingularityPolicy::Throw) != static_cast<int>(SingularityPolicy::Clamp));
    REQUIRE(static_cast<int>(SingularityPolicy::Throw) != static_cast<int>(SingularityPolicy::ReturnZero));
    REQUIRE(static_cast<int>(SingularityPolicy::ReturnNaN) != static_cast<int>(SingularityPolicy::ReturnInf));
    REQUIRE(static_cast<int>(SingularityPolicy::ReturnNaN) != static_cast<int>(SingularityPolicy::Clamp));
    REQUIRE(static_cast<int>(SingularityPolicy::ReturnNaN) != static_cast<int>(SingularityPolicy::ReturnZero));
    REQUIRE(static_cast<int>(SingularityPolicy::ReturnInf) != static_cast<int>(SingularityPolicy::Clamp));
    REQUIRE(static_cast<int>(SingularityPolicy::ReturnInf) != static_cast<int>(SingularityPolicy::ReturnZero));
    REQUIRE(static_cast<int>(SingularityPolicy::Clamp) != static_cast<int>(SingularityPolicy::ReturnZero));
}

///////////////////////////////////////////////////////////////////////////////
//                       POLICY TO STRING CONVERSION                         //
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("SingularityPolicyToString - All policies have names", "[Singularity][Policy]") {
    REQUIRE(std::string(SingularityPolicyToString(SingularityPolicy::Throw)) == "Throw");
    REQUIRE(std::string(SingularityPolicyToString(SingularityPolicy::ReturnNaN)) == "ReturnNaN");
    REQUIRE(std::string(SingularityPolicyToString(SingularityPolicy::ReturnInf)) == "ReturnInf");
    REQUIRE(std::string(SingularityPolicyToString(SingularityPolicy::Clamp)) == "Clamp");
    REQUIRE(std::string(SingularityPolicyToString(SingularityPolicy::ReturnZero)) == "ReturnZero");
}

TEST_CASE("SingularityPolicyToString - Unknown policy", "[Singularity][Policy]") {
    // Cast an invalid value - should return "Unknown"
    auto unknown = static_cast<SingularityPolicy>(999);
    REQUIRE(std::string(SingularityPolicyToString(unknown)) == "Unknown");
}

///////////////////////////////////////////////////////////////////////////////
//                           DEFAULT CONSTANTS                               //
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("Singularity - Default constants are reasonable", "[Singularity][Constants]") {
    REQUIRE(DEFAULT_SINGULARITY_TOL > 0.0);
    REQUIRE(DEFAULT_SINGULARITY_TOL < 1e-6);  // Should be very small
    REQUIRE(DEFAULT_SINGULARITY_TOL == Catch::Approx(1e-12));
    REQUIRE(DEFAULT_POLICY == SingularityPolicy::Throw);  // Safe default
}

///////////////////////////////////////////////////////////////////////////////
//                            IS NEAR ZERO                                   //
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("IsNearZero - Zero value", "[Singularity][Detection]") {
    REQUIRE(IsNearZero(0.0) == true);
    REQUIRE(IsNearZero(0.0, 1e-10) == true);
    REQUIRE(IsNearZero(0.0, 1e-20) == true);
}

TEST_CASE("IsNearZero - Very small positive values", "[Singularity][Detection]") {
    REQUIRE(IsNearZero(1e-15) == true);   // Below default tolerance
    REQUIRE(IsNearZero(1e-13) == true);   // Below default tolerance
    REQUIRE(IsNearZero(1e-11) == false);  // Above default tolerance
    REQUIRE(IsNearZero(1e-6) == false);   // Well above tolerance
}

TEST_CASE("IsNearZero - Very small negative values", "[Singularity][Detection]") {
    REQUIRE(IsNearZero(-1e-15) == true);   // Below default tolerance
    REQUIRE(IsNearZero(-1e-13) == true);   // Below default tolerance
    REQUIRE(IsNearZero(-1e-11) == false);  // Above default tolerance
    REQUIRE(IsNearZero(-1e-6) == false);   // Well above tolerance
}

TEST_CASE("IsNearZero - Custom tolerance", "[Singularity][Detection]") {
    double value = 0.001;
    REQUIRE(IsNearZero(value, 1e-12) == false);  // Stricter tolerance
    REQUIRE(IsNearZero(value, 0.01) == true);    // Looser tolerance
    // Note: Uses strict < comparison, so exactly at tolerance returns false
    REQUIRE(IsNearZero(value, 0.001) == false);  // Exactly at tolerance
    REQUIRE(IsNearZero(value, 0.0011) == true);  // Just above tolerance
}

TEST_CASE("IsNearZero - Normal values are not near zero", "[Singularity][Detection]") {
    REQUIRE(IsNearZero(1.0) == false);
    REQUIRE(IsNearZero(-1.0) == false);
    REQUIRE(IsNearZero(0.5) == false);
    REQUIRE(IsNearZero(1000.0) == false);
}

///////////////////////////////////////////////////////////////////////////////
//                        SPHERICAL ORIGIN DETECTION                         //
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("IsAtSphericalOrigin - Zero radius", "[Singularity][Detection]") {
    REQUIRE(IsAtSphericalOrigin(0.0) == true);
    REQUIRE(IsAtSphericalOrigin(0.0, 1e-10) == true);
}

TEST_CASE("IsAtSphericalOrigin - Near-zero radius", "[Singularity][Detection]") {
    REQUIRE(IsAtSphericalOrigin(1e-15) == true);
    REQUIRE(IsAtSphericalOrigin(1e-13) == true);
    REQUIRE(IsAtSphericalOrigin(1e-11) == false);
}

TEST_CASE("IsAtSphericalOrigin - Positive radius", "[Singularity][Detection]") {
    REQUIRE(IsAtSphericalOrigin(1.0) == false);
    REQUIRE(IsAtSphericalOrigin(0.1) == false);
    REQUIRE(IsAtSphericalOrigin(1e-6) == false);
}

TEST_CASE("IsAtSphericalOrigin - Negative radius (unusual but checked)", "[Singularity][Detection]") {
    // Radii should typically be positive
    // Implementation uses r < tol (no abs), so negative values are always < tol
    REQUIRE(IsAtSphericalOrigin(-1e-15) == true);
    REQUIRE(IsAtSphericalOrigin(-1.0) == true);  // Any negative is < tol
}

///////////////////////////////////////////////////////////////////////////////
//                         SPHERICAL POLE DETECTION                          //
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("IsAtSphericalPole - At theta = 0 (north pole)", "[Singularity][Detection]") {
    REQUIRE(IsAtSphericalPole(0.0) == true);
    REQUIRE(IsAtSphericalPole(1e-15) == true);
    REQUIRE(IsAtSphericalPole(-1e-15) == true);  // Near zero from negative side
}

TEST_CASE("IsAtSphericalPole - At theta = pi (south pole)", "[Singularity][Detection]") {
    const double pi = 3.14159265358979323846;
    REQUIRE(IsAtSphericalPole(pi) == true);
    REQUIRE(IsAtSphericalPole(pi - 1e-15) == true);
    REQUIRE(IsAtSphericalPole(pi + 1e-15) == true);  // Slightly past pi
}

TEST_CASE("IsAtSphericalPole - Away from poles", "[Singularity][Detection]") {
    const double pi = 3.14159265358979323846;
    REQUIRE(IsAtSphericalPole(pi / 2) == false);     // Equator
    REQUIRE(IsAtSphericalPole(pi / 4) == false);     // 45 degrees
    REQUIRE(IsAtSphericalPole(3 * pi / 4) == false); // 135 degrees
    REQUIRE(IsAtSphericalPole(0.1) == false);
    REQUIRE(IsAtSphericalPole(pi - 0.1) == false);
}

TEST_CASE("IsAtSphericalPole - Custom tolerance", "[Singularity][Detection]") {
    const double pi = 3.14159265358979323846;
    // With looser tolerance
    REQUIRE(IsAtSphericalPole(0.001, 0.01) == true);
    REQUIRE(IsAtSphericalPole(pi - 0.001, 0.01) == true);
    // With stricter tolerance
    REQUIRE(IsAtSphericalPole(0.001, 1e-6) == false);
}

///////////////////////////////////////////////////////////////////////////////
//                       CYLINDRICAL AXIS DETECTION                          //
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("IsAtCylindricalAxis - Zero radius", "[Singularity][Detection]") {
    REQUIRE(IsAtCylindricalAxis(0.0) == true);
    REQUIRE(IsAtCylindricalAxis(1e-15) == true);
    REQUIRE(IsAtCylindricalAxis(-1e-15) == true);
}

TEST_CASE("IsAtCylindricalAxis - Away from axis", "[Singularity][Detection]") {
    REQUIRE(IsAtCylindricalAxis(1.0) == false);
    REQUIRE(IsAtCylindricalAxis(0.001) == false);
    REQUIRE(IsAtCylindricalAxis(1e-6) == false);
}

TEST_CASE("IsAtCylindricalAxis - Custom tolerance", "[Singularity][Detection]") {
    REQUIRE(IsAtCylindricalAxis(0.0001, 0.001) == true);
    REQUIRE(IsAtCylindricalAxis(0.0001, 1e-6) == false);
}

///////////////////////////////////////////////////////////////////////////////
//                    COMBINED SPHERICAL SINGULARITY CHECK                   //
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("IsAtSphericalSingularity - Origin only", "[Singularity][Detection]") {
    const double pi = 3.14159265358979323846;
    // r=0, theta not at pole
    REQUIRE(IsAtSphericalSingularity(0.0, pi / 2) == true);
    REQUIRE(IsAtSphericalSingularity(1e-15, pi / 4) == true);
}

TEST_CASE("IsAtSphericalSingularity - Pole only", "[Singularity][Detection]") {
    const double pi = 3.14159265358979323846;
    // r nonzero, theta at pole
    REQUIRE(IsAtSphericalSingularity(1.0, 0.0) == true);
    REQUIRE(IsAtSphericalSingularity(1.0, pi) == true);
    REQUIRE(IsAtSphericalSingularity(1.0, 1e-15) == true);
}

TEST_CASE("IsAtSphericalSingularity - Both origin and pole", "[Singularity][Detection]") {
    // r=0 and theta=0
    REQUIRE(IsAtSphericalSingularity(0.0, 0.0) == true);
}

TEST_CASE("IsAtSphericalSingularity - No singularity", "[Singularity][Detection]") {
    const double pi = 3.14159265358979323846;
    REQUIRE(IsAtSphericalSingularity(1.0, pi / 2) == false);
    REQUIRE(IsAtSphericalSingularity(0.5, pi / 4) == false);
    REQUIRE(IsAtSphericalSingularity(10.0, 1.0) == false);
}

///////////////////////////////////////////////////////////////////////////////
//                    DESCRIBE SPHERICAL SINGULARITY                         //
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("DescribeSphericalSingularity - Origin", "[Singularity][Detection]") {
    const double pi = 3.14159265358979323846;
    std::string desc = DescribeSphericalSingularity(0.0, pi / 2);
    REQUIRE(desc == "origin (r=0)");
    
    desc = DescribeSphericalSingularity(1e-15, pi / 4);
    REQUIRE(desc == "origin (r=0)");
}

TEST_CASE("DescribeSphericalSingularity - Pole", "[Singularity][Detection]") {
    const double pi = 3.14159265358979323846;
    std::string desc = DescribeSphericalSingularity(1.0, 0.0);
    REQUIRE(desc == "pole (sin(θ)=0)");
    
    desc = DescribeSphericalSingularity(1.0, pi);
    REQUIRE(desc == "pole (sin(θ)=0)");
}

TEST_CASE("DescribeSphericalSingularity - Origin takes precedence", "[Singularity][Detection]") {
    // When both origin and pole, origin is returned (checked first)
    std::string desc = DescribeSphericalSingularity(0.0, 0.0);
    REQUIRE(desc == "origin (r=0)");
}

TEST_CASE("DescribeSphericalSingularity - No singularity", "[Singularity][Detection]") {
    const double pi = 3.14159265358979323846;
    std::string desc = DescribeSphericalSingularity(1.0, pi / 2);
    REQUIRE(desc == "none");
}

///////////////////////////////////////////////////////////////////////////////
//                            SAFE DIVIDE                                    //
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("SafeDivide - Normal division", "[Singularity][SafeArithmetic]") {
    REQUIRE(SafeDivide(10.0, 2.0) == Catch::Approx(5.0));
    REQUIRE(SafeDivide(1.0, 4.0) == Catch::Approx(0.25));
    REQUIRE(SafeDivide(-6.0, 3.0) == Catch::Approx(-2.0));
    REQUIRE(SafeDivide(0.0, 5.0) == Catch::Approx(0.0));
}

TEST_CASE("SafeDivide - Policy::Throw on singular", "[Singularity][SafeArithmetic]") {
    REQUIRE_THROWS_AS(
        SafeDivide(1.0, 0.0, SingularityPolicy::Throw),
        DomainError
    );
    REQUIRE_THROWS_AS(
        SafeDivide(1.0, 1e-15, SingularityPolicy::Throw),
        DomainError
    );
}

TEST_CASE("SafeDivide - Policy::Throw with context", "[Singularity][SafeArithmetic]") {
    try {
        SafeDivide(1.0, 0.0, SingularityPolicy::Throw, "test context");
        FAIL("Should have thrown");
    } catch (const DomainError& e) {
        std::string msg = e.what();
        REQUIRE(msg.find("test context") != std::string::npos);
    }
}

TEST_CASE("SafeDivide - Policy::ReturnNaN", "[Singularity][SafeArithmetic]") {
    Real result = SafeDivide(1.0, 0.0, SingularityPolicy::ReturnNaN);
    REQUIRE(std::isnan(result));
    
    result = SafeDivide(-5.0, 1e-15, SingularityPolicy::ReturnNaN);
    REQUIRE(std::isnan(result));
}

TEST_CASE("SafeDivide - Policy::ReturnInf positive numerator", "[Singularity][SafeArithmetic]") {
    Real result = SafeDivide(1.0, 0.0, SingularityPolicy::ReturnInf);
    REQUIRE(std::isinf(result));
    REQUIRE(result > 0);  // Positive infinity
}

TEST_CASE("SafeDivide - Policy::ReturnInf negative numerator", "[Singularity][SafeArithmetic]") {
    Real result = SafeDivide(-1.0, 0.0, SingularityPolicy::ReturnInf);
    REQUIRE(std::isinf(result));
    REQUIRE(result < 0);  // Negative infinity
}

TEST_CASE("SafeDivide - Policy::ReturnInf zero numerator", "[Singularity][SafeArithmetic]") {
    Real result = SafeDivide(0.0, 0.0, SingularityPolicy::ReturnInf);
    REQUIRE(std::isinf(result));
    REQUIRE(result > 0);  // 0 >= 0, so positive infinity
}

TEST_CASE("SafeDivide - Policy::Clamp positive denominator", "[Singularity][SafeArithmetic]") {
    // Clamp denominator to tolerance
    Real result = SafeDivide(1.0, 0.0, SingularityPolicy::Clamp, nullptr, 1e-12);
    REQUIRE(std::isfinite(result));
    REQUIRE(result == Catch::Approx(1.0 / 1e-12));
}

TEST_CASE("SafeDivide - Policy::Clamp negative denominator", "[Singularity][SafeArithmetic]") {
    // Clamp to -tolerance when denominator is negative
    Real result = SafeDivide(1.0, -1e-15, SingularityPolicy::Clamp, nullptr, 1e-12);
    REQUIRE(std::isfinite(result));
    REQUIRE(result == Catch::Approx(1.0 / -1e-12));
}

TEST_CASE("SafeDivide - Policy::ReturnZero", "[Singularity][SafeArithmetic]") {
    Real result = SafeDivide(1.0, 0.0, SingularityPolicy::ReturnZero);
    REQUIRE(result == Catch::Approx(0.0));
    
    result = SafeDivide(1000.0, 1e-15, SingularityPolicy::ReturnZero);
    REQUIRE(result == Catch::Approx(0.0));
}

///////////////////////////////////////////////////////////////////////////////
//                           SAFE INVERSE R                                  //
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("SafeInverseR - Normal operation", "[Singularity][SafeArithmetic]") {
    REQUIRE(SafeInverseR(2.0) == Catch::Approx(0.5));
    REQUIRE(SafeInverseR(0.5) == Catch::Approx(2.0));
    REQUIRE(SafeInverseR(10.0) == Catch::Approx(0.1));
}

TEST_CASE("SafeInverseR - At singularity with policies", "[Singularity][SafeArithmetic]") {
    REQUIRE_THROWS_AS(SafeInverseR(0.0, SingularityPolicy::Throw), DomainError);
    REQUIRE(std::isnan(SafeInverseR(0.0, SingularityPolicy::ReturnNaN)));
    REQUIRE(std::isinf(SafeInverseR(0.0, SingularityPolicy::ReturnInf)));
    REQUIRE(SafeInverseR(0.0, SingularityPolicy::ReturnZero) == Catch::Approx(0.0));
}

TEST_CASE("SafeInverseR - Context in error message", "[Singularity][SafeArithmetic]") {
    try {
        SafeInverseR(0.0, SingularityPolicy::Throw, "custom 1/r context");
        FAIL("Should have thrown");
    } catch (const DomainError& e) {
        std::string msg = e.what();
        REQUIRE(msg.find("custom 1/r context") != std::string::npos);
    }
}

///////////////////////////////////////////////////////////////////////////////
//                          SAFE INVERSE R²                                  //
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("SafeInverseR2 - Normal operation", "[Singularity][SafeArithmetic]") {
    REQUIRE(SafeInverseR2(2.0) == Catch::Approx(0.25));
    REQUIRE(SafeInverseR2(0.5) == Catch::Approx(4.0));
    REQUIRE(SafeInverseR2(10.0) == Catch::Approx(0.01));
}

TEST_CASE("SafeInverseR2 - At singularity with policies", "[Singularity][SafeArithmetic]") {
    REQUIRE_THROWS_AS(SafeInverseR2(0.0, SingularityPolicy::Throw), DomainError);
    REQUIRE(std::isnan(SafeInverseR2(0.0, SingularityPolicy::ReturnNaN)));
    REQUIRE(std::isinf(SafeInverseR2(0.0, SingularityPolicy::ReturnInf)));
    REQUIRE(SafeInverseR2(0.0, SingularityPolicy::ReturnZero) == Catch::Approx(0.0));
}

TEST_CASE("SafeInverseR2 - Near-zero r triggers singularity", "[Singularity][SafeArithmetic]") {
    // r² will be even smaller, so near-zero r triggers singularity
    REQUIRE_THROWS_AS(SafeInverseR2(1e-7, SingularityPolicy::Throw), DomainError);
}

///////////////////////////////////////////////////////////////////////////////
//                      SAFE INVERSE R*SIN(THETA)                            //
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("SafeInverseRSinTheta - Normal operation", "[Singularity][SafeArithmetic]") {
    const double pi = 3.14159265358979323846;
    // At equator (theta = pi/2), sin(theta) = 1
    REQUIRE(SafeInverseRSinTheta(2.0, pi / 2) == Catch::Approx(0.5));
    
    // At theta = pi/6, sin = 0.5
    REQUIRE(SafeInverseRSinTheta(1.0, pi / 6) == Catch::Approx(2.0));
}

TEST_CASE("SafeInverseRSinTheta - Singularity at r=0", "[Singularity][SafeArithmetic]") {
    const double pi = 3.14159265358979323846;
    REQUIRE_THROWS_AS(
        SafeInverseRSinTheta(0.0, pi / 2, SingularityPolicy::Throw),
        DomainError
    );
}

TEST_CASE("SafeInverseRSinTheta - Singularity at pole", "[Singularity][SafeArithmetic]") {
    const double pi = 3.14159265358979323846;
    // At pole, sin(theta) = 0
    REQUIRE_THROWS_AS(
        SafeInverseRSinTheta(1.0, 0.0, SingularityPolicy::Throw),
        DomainError
    );
    REQUIRE_THROWS_AS(
        SafeInverseRSinTheta(1.0, pi, SingularityPolicy::Throw),
        DomainError
    );
}

TEST_CASE("SafeInverseRSinTheta - ReturnNaN policy", "[Singularity][SafeArithmetic]") {
    REQUIRE(std::isnan(SafeInverseRSinTheta(0.0, 1.0, SingularityPolicy::ReturnNaN)));
}

///////////////////////////////////////////////////////////////////////////////
//                     SAFE INVERSE R²*SIN²(THETA)                           //
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("SafeInverseR2Sin2Theta - Normal operation", "[Singularity][SafeArithmetic]") {
    const double pi = 3.14159265358979323846;
    // At equator, sin²(pi/2) = 1
    REQUIRE(SafeInverseR2Sin2Theta(2.0, pi / 2) == Catch::Approx(0.25));
    
    // At theta = pi/6, sin² = 0.25, so 1/(r²*0.25) = 4/r²
    REQUIRE(SafeInverseR2Sin2Theta(1.0, pi / 6) == Catch::Approx(4.0));
}

TEST_CASE("SafeInverseR2Sin2Theta - Singularity at origin", "[Singularity][SafeArithmetic]") {
    const double pi = 3.14159265358979323846;
    REQUIRE_THROWS_AS(
        SafeInverseR2Sin2Theta(0.0, pi / 2, SingularityPolicy::Throw),
        DomainError
    );
}

TEST_CASE("SafeInverseR2Sin2Theta - Singularity at pole", "[Singularity][SafeArithmetic]") {
    REQUIRE_THROWS_AS(
        SafeInverseR2Sin2Theta(1.0, 0.0, SingularityPolicy::Throw),
        DomainError
    );
}

TEST_CASE("SafeInverseR2Sin2Theta - ReturnInf policy", "[Singularity][SafeArithmetic]") {
    const double pi = 3.14159265358979323846;
    Real result = SafeInverseR2Sin2Theta(0.0, pi / 2, SingularityPolicy::ReturnInf);
    REQUIRE(std::isinf(result));
}

///////////////////////////////////////////////////////////////////////////////
//                       SAFE COT(THETA) / R²                                //
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("SafeCotThetaOverR2 - Normal operation", "[Singularity][SafeArithmetic]") {
    const double pi = 3.14159265358979323846;
    // At theta = pi/4, cot = 1, so result = 1/r²
    REQUIRE(SafeCotThetaOverR2(2.0, pi / 4) == Catch::Approx(0.25));
    
    // At equator (theta = pi/2), cot = 0, so result ≈ 0 (with floating-point tolerance)
    REQUIRE(SafeCotThetaOverR2(1.0, pi / 2) == Catch::Approx(0.0).margin(1e-14));
}

TEST_CASE("SafeCotThetaOverR2 - Singularity at origin", "[Singularity][SafeArithmetic]") {
    const double pi = 3.14159265358979323846;
    REQUIRE_THROWS_AS(
        SafeCotThetaOverR2(0.0, pi / 4, SingularityPolicy::Throw),
        DomainError
    );
}

TEST_CASE("SafeCotThetaOverR2 - Singularity at pole", "[Singularity][SafeArithmetic]") {
    // At pole, sin(theta) = 0, so we divide by zero
    REQUIRE_THROWS_AS(
        SafeCotThetaOverR2(1.0, 0.0, SingularityPolicy::Throw),
        DomainError
    );
}

TEST_CASE("SafeCotThetaOverR2 - Negative cotangent", "[Singularity][SafeArithmetic]") {
    const double pi = 3.14159265358979323846;
    // For theta > pi/2, cot is negative
    Real result = SafeCotThetaOverR2(1.0, 3 * pi / 4);
    REQUIRE(result < 0);
    REQUIRE(result == Catch::Approx(-1.0));  // cot(3π/4) = -1
}

TEST_CASE("SafeCotThetaOverR2 - ReturnZero policy", "[Singularity][SafeArithmetic]") {
    Real result = SafeCotThetaOverR2(0.0, 1.0, SingularityPolicy::ReturnZero);
    REQUIRE(result == Catch::Approx(0.0));
}

///////////////////////////////////////////////////////////////////////////////
//                        EDGE CASES AND INTEGRATION                         //
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("Singularity - Negative tolerance is handled", "[Singularity][EdgeCases]") {
    // Negative tolerance: abs(1e-15) < -1e-12 is false
    // The implementation does NOT use abs on tolerance
    REQUIRE(IsNearZero(1e-15, -1e-12) == false);  // 1e-15 < -1e-12 is false
}

TEST_CASE("Singularity - Very large values", "[Singularity][EdgeCases]") {
    double large = 1e100;
    REQUIRE(IsNearZero(large) == false);
    REQUIRE(SafeDivide(large, 2.0) == Catch::Approx(large / 2.0));
}

TEST_CASE("Singularity - Denormalized numbers", "[Singularity][EdgeCases]") {
    double denorm = std::numeric_limits<double>::denorm_min();
    REQUIRE(IsNearZero(denorm) == true);  // Denormalized is very close to zero
    REQUIRE(IsNearZero(denorm, 0.0) == false);  // But not exactly zero with tol=0
}

TEST_CASE("Singularity - Special IEEE values", "[Singularity][EdgeCases]") {
    double inf = std::numeric_limits<double>::infinity();
    double nan = std::numeric_limits<double>::quiet_NaN();
    
    // Infinity is not near zero
    REQUIRE(IsNearZero(inf) == false);
    REQUIRE(IsNearZero(-inf) == false);
    
    // NaN comparisons are false
    REQUIRE(IsNearZero(nan) == false);
}

TEST_CASE("Singularity - Tolerance boundary cases", "[Singularity][EdgeCases]") {
    double tol = 1e-6;
    
    // Uses strict < comparison, so exactly at tolerance returns false
    REQUIRE(IsNearZero(tol, tol) == false);       // |tol| < tol is false
    
    // Just above tolerance
    REQUIRE(IsNearZero(tol * 1.01, tol) == false);
    
    // Just below tolerance
    REQUIRE(IsNearZero(tol * 0.99, tol) == true);
}

} // namespace MML::Tests::Core::SingularityHandlingTests

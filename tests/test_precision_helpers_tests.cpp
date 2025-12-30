#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "TestPrecision.h"
#include "TestMatchers.h"

using namespace MML::Testing;
using namespace MML::Testing::Matchers;

TEST_CASE("TestPrecision - Precision Detection", "[test-infrastructure]") {
    TEST_PRECISION_INFO();
    
    SECTION("Detect Real type") {
        // Use runtime detection instead of preprocessor defines
        RealPrecision detected = GetRealPrecision();
        const char* name = GetRealPrecisionName();
        
        if (detected == RealPrecision::Float) {
            REQUIRE(std::string(name) == "float");
            REQUIRE(sizeof(Real) == 4);
        } else if (detected == RealPrecision::LongDouble) {
            REQUIRE(std::string(name) == "long double");
        } else {
            REQUIRE(detected == RealPrecision::Double);
            REQUIRE(std::string(name) == "double");
            REQUIRE(sizeof(Real) == 8);
        }
    }
    
    SECTION("Precision digits") {
        int digits = GetRealDigits();
        REQUIRE(digits > 0);
        
        // Use runtime detection instead of preprocessor defines
        if (GetRealPrecision() == RealPrecision::Float) {
            REQUIRE(digits >= 6);  // float has ~7 digits
        } else if (GetRealPrecision() == RealPrecision::LongDouble) {
            REQUIRE(digits >= 18); // long double has 18-33 digits
        } else {
            REQUIRE(digits >= 15); // double has ~15 digits
        }
    }
}

TEST_CASE("TestPrecision - Tolerance Scaling", "[test-infrastructure]") {
    TEST_PRECISION_INFO();
    
    SECTION("ScaleTolerance adjusts for precision") {
        Real baseline = REAL(1e-9);
        Real scaled = ScaleTolerance(baseline);
        
        // Use runtime detection instead of preprocessor defines
        if (GetRealPrecision() == RealPrecision::Float) {
            // Float should scale UP (looser tolerance)
            REQUIRE(scaled > baseline);
            REQUIRE_THAT(scaled, RealApprox(REAL(1e-7))); // 100x larger
        } else if (GetRealPrecision() == RealPrecision::LongDouble) {
            // Long double should scale DOWN (tighter tolerance)
            REQUIRE(scaled < baseline);
            REQUIRE_THAT(scaled, RealApprox(REAL(1e-11))); // 100x smaller
        } else {
            // Double should stay the same
            REQUIRE(scaled == baseline);
        }
    }
    
    SECTION("Predefined tolerances are accessible") {
        REQUIRE(Tolerance::Loose > REAL(0));
        REQUIRE(Tolerance::Standard > REAL(0));
        REQUIRE(Tolerance::Strict > REAL(0));
        REQUIRE(Tolerance::Loose > Tolerance::Standard);
        REQUIRE(Tolerance::Standard > Tolerance::Strict);
    }
    
    SECTION("Operation-specific tolerances") {
        REQUIRE(Tolerance::Integration::Standard > REAL(0));
        REQUIRE(Tolerance::Differentiation::Standard > REAL(0));
        REQUIRE(Tolerance::LinearAlgebra::Standard > REAL(0));
        REQUIRE(Tolerance::RootFinding::Standard > REAL(0));
        REQUIRE(Tolerance::Statistics::Standard > REAL(0));
    }
}

TEST_CASE("TestPrecision - ApproxEqual Helper", "[test-infrastructure]") {
    TEST_PRECISION_INFO();
    
    SECTION("Exact equality") {
        REQUIRE(ApproxEqual(REAL(1.0), REAL(1.0)));
        REQUIRE(ApproxEqual(REAL(0.0), REAL(0.0)));
        REQUIRE(ApproxEqual(REAL(-5.5), REAL(-5.5)));
    }
    
    SECTION("Close values within tolerance") {
        Real a = REAL(1.0);
        Real b = REAL(1.0) + Tolerance::Standard * REAL(0.5);
        REQUIRE(ApproxEqual(a, b, Tolerance::Standard));
    }
    
    SECTION("Different values outside tolerance") {
        Real a = REAL(1.0);
        Real b = REAL(1.1);
        REQUIRE_FALSE(ApproxEqual(a, b, Tolerance::Strict));
    }
    
    SECTION("ApproxZero") {
        REQUIRE(ApproxZero(REAL(0.0)));
        REQUIRE(ApproxZero(Tolerance::Standard * REAL(0.1)));
        REQUIRE_FALSE(ApproxZero(REAL(1.0)));
    }
}

TEST_CASE("TestPrecision - Machine Constants", "[test-infrastructure]") {
    TEST_PRECISION_INFO();
    
    SECTION("Machine epsilon") {
        Real eps = MachineEpsilon();
        REQUIRE(eps > REAL(0));
        REQUIRE(eps < REAL(1.0));
        
        // 1 + eps should be > 1
        REQUIRE(REAL(1.0) + eps > REAL(1.0));
    }
    
    SECTION("Min and Max values") {
        Real minVal = MinPositive();
        Real maxVal = MaxFinite();
        
        REQUIRE(minVal > REAL(0));
        REQUIRE(maxVal > REAL(1.0));
        REQUIRE(maxVal > minVal);
    }
}

TEST_CASE("TestPrecision - REAL macro", "[test-infrastructure]") {
    TEST_PRECISION_INFO();
    
    SECTION("REAL macro prevents narrowing") {
        // These should compile without warnings regardless of Real type
        Real a = REAL(3.14159265358979323846);
        Real b = REAL(1e-15);
        Real c = REAL(0.0);
        Real d = REAL(1.0);
        
        REQUIRE(a > REAL(3.0));
        REQUIRE(b > REAL(0.0));
        REQUIRE(c == REAL(0.0));
        REQUIRE(d == REAL(1.0));
    }
}

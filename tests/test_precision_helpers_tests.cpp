#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "TestPrecision.h"
#include "TestMatchers.h"

using namespace MML::Testing;
using namespace MML::Testing::Matchers;

namespace MML::Tests::TestInfrastructure::PrecisionHelpersTests {

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
        } else if (GetRealPrecision() == RealPrecision::LongDouble && LongDoubleHasExtraPrecision) {
            REQUIRE(digits >= 18); // long double has 18-33 digits (GCC/Linux 80-bit)
        } else {
            REQUIRE(digits >= 15); // double (or MSVC long double) has ~15 digits
        }
    }
}

TEST_CASE("TestPrecision - Tolerance Scaling", "[test-infrastructure]") {
    TEST_PRECISION_INFO();
    
    SECTION("ScaleTolerance adjusts for precision") {
        Real baseline = static_cast<MML::Real>(1e-9);
        Real scaled = ScaleTolerance(baseline);
        
        // Use runtime detection instead of preprocessor defines
        if (GetRealPrecision() == RealPrecision::Float) {
            // Float should scale UP (looser tolerance)
            REQUIRE(scaled > baseline);
            REQUIRE_THAT(scaled, RealApprox(static_cast<MML::Real>(1e-7))); // 100x larger
        } else if (GetRealPrecision() == RealPrecision::LongDouble && LongDoubleHasExtraPrecision) {
            // Long double should scale DOWN (tighter tolerance) - only when genuinely 80-bit
            REQUIRE(scaled < baseline);
            REQUIRE_THAT(scaled, RealApprox(static_cast<MML::Real>(1e-11))); // 100x smaller
        } else {
            // Double (or MSVC long double with no extra precision) should stay the same
            REQUIRE(scaled == baseline);
        }
    }
    
    SECTION("Predefined tolerances are accessible") {
        REQUIRE(Tolerance::Loose > static_cast<MML::Real>(0));
        REQUIRE(Tolerance::Standard > static_cast<MML::Real>(0));
        REQUIRE(Tolerance::Strict > static_cast<MML::Real>(0));
        REQUIRE(Tolerance::Loose > Tolerance::Standard);
        REQUIRE(Tolerance::Standard > Tolerance::Strict);
    }
    
    SECTION("Operation-specific tolerances") {
        REQUIRE(Tolerance::Integration::Standard > static_cast<MML::Real>(0));
        REQUIRE(Tolerance::Differentiation::Standard > static_cast<MML::Real>(0));
        REQUIRE(Tolerance::LinearAlgebra::Standard > static_cast<MML::Real>(0));
        REQUIRE(Tolerance::RootFinding::Standard > static_cast<MML::Real>(0));
        REQUIRE(Tolerance::Statistics::Standard > static_cast<MML::Real>(0));
    }
}

TEST_CASE("TestPrecision - ApproxEqual Helper", "[test-infrastructure]") {
    TEST_PRECISION_INFO();
    
    SECTION("Exact equality") {
        REQUIRE(ApproxEqual(static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0)));
        REQUIRE(ApproxEqual(static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0)));
        REQUIRE(ApproxEqual(static_cast<MML::Real>(-5.5), static_cast<MML::Real>(-5.5)));
    }
    
    SECTION("Close values within tolerance") {
        Real a = static_cast<MML::Real>(1.0);
        Real b = static_cast<MML::Real>(1.0) + Tolerance::Standard * static_cast<MML::Real>(0.5);
        REQUIRE(ApproxEqual(a, b, Tolerance::Standard));
    }
    
    SECTION("Different values outside tolerance") {
        Real a = static_cast<MML::Real>(1.0);
        Real b = static_cast<MML::Real>(1.1);
        REQUIRE_FALSE(ApproxEqual(a, b, Tolerance::Strict));
    }
    
    SECTION("ApproxZero") {
        REQUIRE(ApproxZero(static_cast<MML::Real>(0.0)));
        REQUIRE(ApproxZero(Tolerance::Standard * static_cast<MML::Real>(0.1)));
        REQUIRE_FALSE(ApproxZero(static_cast<MML::Real>(1.0)));
    }
}

TEST_CASE("TestPrecision - Machine Constants", "[test-infrastructure]") {
    TEST_PRECISION_INFO();
    
    SECTION("Machine epsilon") {
        Real eps = MachineEpsilon();
        REQUIRE(eps > static_cast<MML::Real>(0));
        REQUIRE(eps < static_cast<MML::Real>(1.0));
        
        // 1 + eps should be > 1
        REQUIRE(static_cast<MML::Real>(1.0) + eps > static_cast<MML::Real>(1.0));
    }
    
    SECTION("Min and Max values") {
        Real minVal = MinPositive();
        Real maxVal = MaxFinite();
        
        REQUIRE(minVal > static_cast<MML::Real>(0));
        REQUIRE(maxVal > static_cast<MML::Real>(1.0));
        REQUIRE(maxVal > minVal);
    }
}

TEST_CASE("TestPrecision - REAL macro", "[test-infrastructure]") {
    TEST_PRECISION_INFO();
    
    SECTION("REAL macro prevents narrowing") {
        // These should compile without warnings regardless of Real type
        Real a = static_cast<MML::Real>(3.14159265358979323846);
        Real b = static_cast<MML::Real>(1e-15);
        Real c = static_cast<MML::Real>(0.0);
        Real d = static_cast<MML::Real>(1.0);
        
        REQUIRE(a > static_cast<MML::Real>(3.0));
        REQUIRE(b > static_cast<MML::Real>(0.0));
        REQUIRE(c == static_cast<MML::Real>(0.0));
        REQUIRE(d == static_cast<MML::Real>(1.0));
    }
}

} // namespace MML::Tests::TestInfrastructure::PrecisionHelpersTests

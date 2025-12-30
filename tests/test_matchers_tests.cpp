#include <catch2/catch_test_macros.hpp>
#include "TestMatchers.h"

// No using namespace - avoiding potential macro conflicts

TEST_CASE("TestMatchers - RealEquals", "[test-infrastructure][matchers]") {
    TEST_PRECISION_INFO();
    
    SECTION("Exact equality") {
        REQUIRE_THAT(REAL(1.0), MML::Testing::Matchers::RealEquals(REAL(1.0)));
        REQUIRE_THAT(REAL(0.0), MML::Testing::Matchers::RealEquals(REAL(0.0)));
        REQUIRE_THAT(REAL(-5.5), MML::Testing::Matchers::RealEquals(REAL(-5.5)));
    }
    
    SECTION("Not equal") {
        Real a = REAL(1.0);
        Real b = REAL(1.0001);
        REQUIRE_FALSE(MML::Testing::Matchers::RealEquals(a).match(b));
    }
    
    SECTION("NaN handling") {
        Real nan1 = std::numeric_limits<Real>::quiet_NaN();
        Real nan2 = std::numeric_limits<Real>::quiet_NaN();
        REQUIRE_THAT(nan1, MML::Testing::Matchers::RealEquals(nan2));  // NaN == NaN in this matcher
    }
}

TEST_CASE("TestMatchers - RealApprox", "[test-infrastructure][matchers]") {
    TEST_PRECISION_INFO();
    
    using namespace MML::Testing;
    using namespace MML::Testing::Matchers;
    
    SECTION("Default tolerance") {
        Real a = REAL(1.0);
        Real b = REAL(1.0) + Tolerance::Standard * REAL(0.1);
        REQUIRE_THAT(a, RealApprox(b));
    }
    
    SECTION("Custom epsilon") {
        Real a = REAL(1.0);
        Real b = REAL(1.00001);
        REQUIRE_THAT(a, RealApprox(b).epsilon(REAL(0.0001)));
    }
    
    SECTION("Custom margin") {
        Real a = REAL(1.0);
        Real b = REAL(1.01);
        REQUIRE_THAT(a, RealApprox(b).margin(REAL(0.02)));
    }
    
    SECTION("Values outside tolerance fail") {
        Real a = REAL(1.0);
        Real b = REAL(2.0);
        REQUIRE_FALSE(RealApprox(b).match(a));
    }
}

TEST_CASE("TestMatchers - RealWithinRel", "[test-infrastructure][matchers]") {
    TEST_PRECISION_INFO();
    
    using namespace MML::Testing;
    using namespace MML::Testing::Matchers;
    
    SECTION("Within relative tolerance") {
        Real expected = REAL(100.0);
        Real actual = REAL(100.1);
        
        // 0.1% error = 0.001 relative
        REQUIRE_THAT(actual, RealWithinRel(expected, REAL(0.01)));
    }
    
    SECTION("Tolerance auto-scales with precision") {
        Real expected = REAL(1.0);
        Real actual = REAL(1.0) + ScaleRelTolerance(REAL(1e-10)) * REAL(0.5);
        
        // Should pass with auto-scaled tolerance
        REQUIRE_THAT(actual, RealWithinRel(expected, REAL(1e-10)));
    }
    
    SECTION("No scaling version") {
        Real expected = REAL(1.0);
        Real actual = REAL(1.0) + REAL(1e-10) * REAL(0.5);
        
        REQUIRE_THAT(actual, RealWithinRelNoScale(expected, REAL(1e-9)));
    }
}

TEST_CASE("TestMatchers - RealWithinAbs", "[test-infrastructure][matchers]") {
    TEST_PRECISION_INFO();
    
    using namespace MML::Testing;
    using namespace MML::Testing::Matchers;
    
    SECTION("Within absolute tolerance") {
        Real expected = REAL(1.0);
        Real actual = REAL(1.05);
        
        REQUIRE_THAT(actual, RealWithinAbs(expected, REAL(0.1)));
    }
    
    SECTION("Tolerance auto-scales") {
        Real expected = REAL(1.0);
        Real actual = REAL(1.0) + ScaleTolerance(REAL(1e-10)) * REAL(0.5);
        
        REQUIRE_THAT(actual, RealWithinAbs(expected, REAL(1e-10)));
    }
}

TEST_CASE("TestMatchers - RealIsZero", "[test-infrastructure][matchers]") {
    TEST_PRECISION_INFO();
    
    using namespace MML::Testing;
    using namespace MML::Testing::Matchers;
    
    SECTION("Exact zero") {
        REQUIRE_THAT(REAL(0.0), RealIsZero());
    }
    
    SECTION("Approximately zero") {
        Real tiny = Tolerance::Standard * REAL(0.1);
        REQUIRE_THAT(tiny, RealIsZero());
    }
    
    SECTION("Not zero") {
        REQUIRE_FALSE(RealIsZero().match(REAL(1.0)));
    }
}

TEST_CASE("TestMatchers - RealIsPositive/Negative", "[test-infrastructure][matchers]") {
    TEST_PRECISION_INFO();
    
    using namespace MML::Testing::Matchers;
    
    SECTION("Positive values") {
        REQUIRE_THAT(REAL(1.0), RealIsPositive());
        REQUIRE_THAT(REAL(0.0001), RealIsPositive());
        REQUIRE_FALSE(RealIsPositive().match(REAL(0.0)));
        REQUIRE_FALSE(RealIsPositive().match(REAL(-1.0)));
    }
    
    SECTION("Negative values") {
        REQUIRE_THAT(REAL(-1.0), RealIsNegative());
        REQUIRE_THAT(REAL(-0.0001), RealIsNegative());
        REQUIRE_FALSE(RealIsNegative().match(REAL(0.0)));
        REQUIRE_FALSE(RealIsNegative().match(REAL(1.0)));
    }
}

TEST_CASE("TestMatchers - Real-world usage examples", "[test-infrastructure][matchers]") {
    TEST_PRECISION_INFO();
    
    using namespace MML::Testing;
    using namespace MML::Testing::Matchers;
    
    SECTION("Integration result") {
        // Simulating integral result that might have small errors
        Real integral_result = REAL(2.0) + Tolerance::Integration::Standard * REAL(0.5);
        Real expected = REAL(2.0);
        
        REQUIRE_THAT(integral_result, RealWithinAbs(expected, Tolerance::Integration::Standard));
    }
    
    SECTION("Derivative result") {
        // Derivative with numerical error
        Real deriv_result = REAL(3.14159) + Tolerance::Differentiation::Standard * REAL(0.3);
        Real expected = REAL(3.14159);
        
        REQUIRE_THAT(deriv_result, RealApprox(expected).epsilon(Tolerance::Differentiation::Standard));
    }
    
    SECTION("Root finding result") {
        // Root near zero
        Real root = Tolerance::RootFinding::Standard * REAL(0.1);
        
        REQUIRE_THAT(root, RealIsZero(Tolerance::RootFinding::Standard));
    }
    
    SECTION("Statistical mean") {
        Real mean = REAL(100.0) + Tolerance::Statistics::Standard * REAL(0.8);
        Real expected = REAL(100.0);
        
        REQUIRE_THAT(mean, RealWithinRel(expected, Tolerance::Statistics::Standard));
    }
}

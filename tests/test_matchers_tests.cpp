#include <catch2/catch_test_macros.hpp>
#include "TestMatchers.h"

// No using namespace - avoiding potential macro conflicts

namespace MML::Tests::TestInfrastructure::MatchersTests {

TEST_CASE("TestMatchers - RealEquals", "[test-infrastructure][matchers]") {
    TEST_PRECISION_INFO();
    
    SECTION("Exact equality") {
        REQUIRE_THAT(static_cast<MML::Real>(1.0), MML::Testing::Matchers::RealEquals(static_cast<MML::Real>(1.0)));
        REQUIRE_THAT(static_cast<MML::Real>(0.0), MML::Testing::Matchers::RealEquals(static_cast<MML::Real>(0.0)));
        REQUIRE_THAT(static_cast<MML::Real>(-5.5), MML::Testing::Matchers::RealEquals(static_cast<MML::Real>(-5.5)));
    }
    
    SECTION("Not equal") {
        Real a = static_cast<MML::Real>(1.0);
        Real b = static_cast<MML::Real>(1.0001);
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
        Real a = static_cast<MML::Real>(1.0);
        Real b = static_cast<MML::Real>(1.0) + Tolerance::Standard * static_cast<MML::Real>(0.1);
        REQUIRE_THAT(a, RealApprox(b));
    }
    
    SECTION("Custom epsilon") {
        Real a = static_cast<MML::Real>(1.0);
        Real b = static_cast<MML::Real>(1.00001);
        REQUIRE_THAT(a, RealApprox(b).epsilon(static_cast<MML::Real>(0.0001)));
    }
    
    SECTION("Custom margin") {
        Real a = static_cast<MML::Real>(1.0);
        Real b = static_cast<MML::Real>(1.01);
        REQUIRE_THAT(a, RealApprox(b).margin(static_cast<MML::Real>(0.02)));
    }
    
    SECTION("Values outside tolerance fail") {
        Real a = static_cast<MML::Real>(1.0);
        Real b = static_cast<MML::Real>(2.0);
        REQUIRE_FALSE(RealApprox(b).match(a));
    }
}

TEST_CASE("TestMatchers - RealWithinRel", "[test-infrastructure][matchers]") {
    TEST_PRECISION_INFO();
    
    using namespace MML::Testing;
    using namespace MML::Testing::Matchers;
    
    SECTION("Within relative tolerance") {
        Real expected = static_cast<MML::Real>(100.0);
        Real actual = static_cast<MML::Real>(100.1);
        
        // 0.1% error = 0.001 relative
        REQUIRE_THAT(actual, RealWithinRel(expected, static_cast<MML::Real>(0.01)));
    }
    
    SECTION("Tolerance respects precision") {
        Real expected = static_cast<MML::Real>(1.0);
        Real actual = static_cast<MML::Real>(1.0) + static_cast<MML::Real>(1e-10) * static_cast<MML::Real>(0.5);
        
        // Tolerance is used as-is (no auto-scaling)
        REQUIRE_THAT(actual, RealWithinRel(expected, static_cast<MML::Real>(1e-9)));
    }
    
    SECTION("No scaling version") {
        Real expected = static_cast<MML::Real>(1.0);
        Real actual = static_cast<MML::Real>(1.0) + static_cast<MML::Real>(1e-10) * static_cast<MML::Real>(0.5);
        
        REQUIRE_THAT(actual, RealWithinRelNoScale(expected, static_cast<MML::Real>(1e-9)));
    }
}

TEST_CASE("TestMatchers - RealWithinAbs", "[test-infrastructure][matchers]") {
    TEST_PRECISION_INFO();
    
    using namespace MML::Testing;
    using namespace MML::Testing::Matchers;
    
    SECTION("Within absolute tolerance") {
        Real expected = static_cast<MML::Real>(1.0);
        Real actual = static_cast<MML::Real>(1.05);
        
        REQUIRE_THAT(actual, RealWithinAbs(expected, static_cast<MML::Real>(0.1)));
    }
    
    SECTION("Tolerance respects precision") {
        Real expected = static_cast<MML::Real>(1.0);
        Real actual = static_cast<MML::Real>(1.0) + static_cast<MML::Real>(1e-10) * static_cast<MML::Real>(0.5);
        
        REQUIRE_THAT(actual, RealWithinAbs(expected, static_cast<MML::Real>(1e-9)));
    }
}

TEST_CASE("TestMatchers - RealIsZero", "[test-infrastructure][matchers]") {
    TEST_PRECISION_INFO();
    
    using namespace MML::Testing;
    using namespace MML::Testing::Matchers;
    
    SECTION("Exact zero") {
        REQUIRE_THAT(static_cast<MML::Real>(0.0), RealIsZero());
    }
    
    SECTION("Approximately zero") {
        Real tiny = Tolerance::Standard * static_cast<MML::Real>(0.1);
        REQUIRE_THAT(tiny, RealIsZero());
    }
    
    SECTION("Not zero") {
        REQUIRE_FALSE(RealIsZero().match(static_cast<MML::Real>(1.0)));
    }
}

TEST_CASE("TestMatchers - RealIsPositive/Negative", "[test-infrastructure][matchers]") {
    TEST_PRECISION_INFO();
    
    using namespace MML::Testing::Matchers;
    
    SECTION("Positive values") {
        REQUIRE_THAT(static_cast<MML::Real>(1.0), RealIsPositive());
        REQUIRE_THAT(static_cast<MML::Real>(0.0001), RealIsPositive());
        REQUIRE_FALSE(RealIsPositive().match(static_cast<MML::Real>(0.0)));
        REQUIRE_FALSE(RealIsPositive().match(static_cast<MML::Real>(-1.0)));
    }
    
    SECTION("Negative values") {
        REQUIRE_THAT(static_cast<MML::Real>(-1.0), RealIsNegative());
        REQUIRE_THAT(static_cast<MML::Real>(-0.0001), RealIsNegative());
        REQUIRE_FALSE(RealIsNegative().match(static_cast<MML::Real>(0.0)));
        REQUIRE_FALSE(RealIsNegative().match(static_cast<MML::Real>(1.0)));
    }
}

TEST_CASE("TestMatchers - Real-world usage examples", "[test-infrastructure][matchers]") {
    TEST_PRECISION_INFO();
    
    using namespace MML::Testing;
    using namespace MML::Testing::Matchers;
    
    SECTION("Integration result") {
        // Simulating integral result that might have small errors
        Real integral_result = static_cast<MML::Real>(2.0) + Tolerance::Integration::Standard * static_cast<MML::Real>(0.5);
        Real expected = static_cast<MML::Real>(2.0);
        
        REQUIRE_THAT(integral_result, RealWithinAbs(expected, Tolerance::Integration::Standard));
    }
    
    SECTION("Derivative result") {
        // Derivative with numerical error
        Real deriv_result = static_cast<MML::Real>(3.14159) + Tolerance::Differentiation::Standard * static_cast<MML::Real>(0.3);
        Real expected = static_cast<MML::Real>(3.14159);
        
        REQUIRE_THAT(deriv_result, RealApprox(expected).epsilon(Tolerance::Differentiation::Standard));
    }
    
    SECTION("Root finding result") {
        // Root near zero
        Real root = Tolerance::RootFinding::Standard * static_cast<MML::Real>(0.1);
        
        REQUIRE_THAT(root, RealIsZero(Tolerance::RootFinding::Standard));
    }
    
    SECTION("Statistical mean") {
        Real mean = static_cast<MML::Real>(100.0) + Tolerance::Statistics::Standard * static_cast<MML::Real>(0.8);
        Real expected = static_cast<MML::Real>(100.0);
        
        REQUIRE_THAT(mean, RealWithinAbs(expected, Tolerance::Statistics::Standard));
    }
}

} // namespace MML::Tests::TestInfrastructure::MatchersTests

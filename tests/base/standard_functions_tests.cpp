///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML) Tests                            ///
///                                                                                   ///
///  File:        standard_functions_tests.cpp                                        ///
///  Description: Tests for StandardFunctions.h                                       ///
///               Focus on MML implementations (factorial functions)                  ///
///               Wrapper functions around std:: are not tested (testing stdlib)      ///
///                                                                                   ///
///  Coverage:    Functions::Factorial                                                ///
///               Functions::FactorialInt                                             ///
///               Functions::FactorialStirling                                        ///
///////////////////////////////////////////////////////////////////////////////////////////

#include "../TestPrecision.h"
#include "../../mml/base/StandardFunctions.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <cmath>
#include <limits>

using namespace MML;
using namespace MML::Functions;
using Catch::Approx;

namespace MML::Tests::Base::StandardFunctionsTests {

///////////////////////////////////////////////////////////////////////////////////////////
///                         Factorial Tests                                            ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Functions::Factorial - Known values", "[StandardFunctions][Factorial]")
{
    SECTION("Factorial of 0 is 1")
    {
        REQUIRE(Factorial(0) == Approx(1.0));
    }

    SECTION("Factorial of 1 is 1")
    {
        REQUIRE(Factorial(1) == Approx(1.0));
    }

    SECTION("Factorial of small integers")
    {
        REQUIRE(Factorial(2) == Approx(2.0));
        REQUIRE(Factorial(3) == Approx(6.0));
        REQUIRE(Factorial(4) == Approx(24.0));
        REQUIRE(Factorial(5) == Approx(120.0));
        REQUIRE(Factorial(6) == Approx(720.0));
        REQUIRE(Factorial(7) == Approx(5040.0));
        REQUIRE(Factorial(8) == Approx(40320.0));
        REQUIRE(Factorial(9) == Approx(362880.0));
        REQUIRE(Factorial(10) == Approx(3628800.0));
    }

    SECTION("Factorial of larger integers")
    {
        REQUIRE(Factorial(12) == Approx(479001600.0));
        REQUIRE(Factorial(15) == Approx(1307674368000.0));
        REQUIRE(Factorial(20) == Approx(2432902008176640000.0).epsilon(1e-10));
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         FactorialInt Tests                                         ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Functions::FactorialInt - Known values", "[StandardFunctions][FactorialInt]")
{
    SECTION("FactorialInt of 0 is 1")
    {
        REQUIRE(FactorialInt(0) == 1);
    }

    SECTION("FactorialInt of 1 is 1")
    {
        REQUIRE(FactorialInt(1) == 1);
    }

    SECTION("FactorialInt of small integers - exact values")
    {
        REQUIRE(FactorialInt(2) == 2);
        REQUIRE(FactorialInt(3) == 6);
        REQUIRE(FactorialInt(4) == 24);
        REQUIRE(FactorialInt(5) == 120);
        REQUIRE(FactorialInt(6) == 720);
        REQUIRE(FactorialInt(7) == 5040);
        REQUIRE(FactorialInt(8) == 40320);
        REQUIRE(FactorialInt(9) == 362880);
        REQUIRE(FactorialInt(10) == 3628800);
    }

    SECTION("FactorialInt of larger integers - exact values")
    {
        REQUIRE(FactorialInt(12) == 479001600LL);
        REQUIRE(FactorialInt(15) == 1307674368000LL);
        REQUIRE(FactorialInt(20) == 2432902008176640000LL);
    }

    SECTION("FactorialInt matches Factorial for small values")
    {
        for (int n = 0; n <= 12; ++n)
        {
            REQUIRE(static_cast<Real>(FactorialInt(n)) == Approx(Factorial(n)));
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         FactorialStirling Tests                                    ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Functions::FactorialStirling - Approximation quality", "[StandardFunctions][FactorialStirling]")
{
    SECTION("FactorialStirling of 0 returns 1 (special case)")
    {
        REQUIRE(FactorialStirling(0) == Approx(1.0));
    }

    SECTION("FactorialStirling of 1 returns 1 (special case)")
    {
        REQUIRE(FactorialStirling(1) == Approx(1.0));
    }

    SECTION("FactorialStirling approximates well for n >= 10")
    {
        // Stirling's approximation: n! ≈ sqrt(2πn) * (n/e)^n
        // Error is O(1/n), so accuracy improves for larger n
        
        // For n=10, error should be < 1%
        Real exact_10 = Factorial(10);
        Real stirling_10 = FactorialStirling(10);
        Real rel_error_10 = std::abs(stirling_10 - exact_10) / exact_10;
        REQUIRE(rel_error_10 < 0.01);  // < 1% error

        // For n=20, error should be < 0.5%
        Real exact_20 = Factorial(20);
        Real stirling_20 = FactorialStirling(20);
        Real rel_error_20 = std::abs(stirling_20 - exact_20) / exact_20;
        REQUIRE(rel_error_20 < 0.005);  // < 0.5% error

        // For n=50, error should be < 0.2%
        Real exact_50 = std::tgamma(51.0);  // Use gamma function for large factorials
        Real stirling_50 = FactorialStirling(50);
        Real rel_error_50 = std::abs(stirling_50 - exact_50) / exact_50;
        REQUIRE(rel_error_50 < 0.002);  // < 0.2% error
    }

    SECTION("FactorialStirling is always less than exact factorial")
    {
        // Stirling's approximation is always an underestimate
        for (int n = 2; n <= 20; ++n)
        {
            Real exact = Factorial(n);
            Real stirling = FactorialStirling(n);
            INFO("n = " << n << ", exact = " << exact << ", stirling = " << stirling);
            REQUIRE(stirling < exact);
        }
    }

    SECTION("FactorialStirling throws on negative input")
    {
        REQUIRE_THROWS_AS(FactorialStirling(-1), std::domain_error);
        REQUIRE_THROWS_AS(FactorialStirling(-10), std::domain_error);
    }
}

TEST_CASE("Functions::Factorial - Negative input throws", "[StandardFunctions][Factorial]")
{
    REQUIRE_THROWS_AS(Factorial(-1), std::domain_error);
    REQUIRE_THROWS_AS(Factorial(-5), std::domain_error);
    REQUIRE_THROWS_AS(FactorialInt(-1), std::domain_error);
    REQUIRE_THROWS_AS(FactorialInt(-5), std::domain_error);
    REQUIRE_THROWS_AS(FactorialStirling(-1), std::domain_error);
}

TEST_CASE("Functions::FactorialStirling - Convergence to exact", "[StandardFunctions][FactorialStirling]")
{
    SECTION("Relative error decreases as n increases")
    {
        Real prev_rel_error = 1.0;
        
        for (int n = 5; n <= 30; n += 5)
        {
            Real exact = std::tgamma(n + 1.0);  // n! = Gamma(n+1)
            Real stirling = FactorialStirling(n);
            Real rel_error = std::abs(stirling - exact) / exact;
            
            INFO("n = " << n << ", rel_error = " << rel_error << ", prev = " << prev_rel_error);
            REQUIRE(rel_error < prev_rel_error);  // Error should decrease
            
            prev_rel_error = rel_error;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
///                Reciprocal Trig Singularity Tests                                  ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Functions::Reciprocal trig - singularity guards", "[StandardFunctions][Trig]")
{
    const Real pi = Constants::PI;

    SECTION("Sec throws at cos(x)=0")
    {
        REQUIRE_THROWS_AS(Sec(pi / 2), std::domain_error);
        REQUIRE_THROWS_AS(Sec(-pi / 2), std::domain_error);
        REQUIRE_NOTHROW(Sec(0.0));  // cos(0) = 1
    }

    SECTION("Csc throws at sin(x)=0")
    {
        REQUIRE_THROWS_AS(Csc(0.0), std::domain_error);
        REQUIRE_THROWS_AS(Csc(pi), std::domain_error);
        REQUIRE_NOTHROW(Csc(pi / 2));  // sin(pi/2) = 1
    }

    SECTION("Ctg throws at tan(x)=0")
    {
        REQUIRE_THROWS_AS(Ctg(0.0), std::domain_error);
        REQUIRE_THROWS_AS(Ctg(pi), std::domain_error);
        REQUIRE_NOTHROW(Ctg(pi / 4));  // tan(pi/4) = 1
    }

    SECTION("Csch throws at sinh(x)=0")
    {
        REQUIRE_THROWS_AS(Csch(0.0), std::domain_error);
        REQUIRE_NOTHROW(Csch(1.0));  // sinh(1) != 0
    }

    SECTION("Ctgh throws at tanh(x)=0")
    {
        REQUIRE_THROWS_AS(Ctgh(0.0), std::domain_error);
        REQUIRE_NOTHROW(Ctgh(1.0));  // tanh(1) != 0
    }
}

} // namespace MML::Tests::Base::StandardFunctionsTests

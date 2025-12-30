#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"
#include <cmath>

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "algorithms/ChebyshevApproximation.h"
#include "core/OrthogonalBasis/ChebyshevBasis.h"
#include "base/Function.h"
#endif

using namespace MML;
using namespace MML::Testing;

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace MML::Tests::Base::ChebyshevTests
{
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///                         CHEBYSHEV POLYNOMIAL T_n TESTS                              ///
    ///////////////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("ChebyshevT - Base cases T_0 and T_1", "[chebyshev][polynomial]")
    {
        // T_0(x) = 1 for all x
        REQUIRE_THAT(ChebyshevT(0, REAL(0.0)), WithinAbs(REAL(1.0), REAL(1e-14)));
        REQUIRE_THAT(ChebyshevT(0, REAL(0.5)), WithinAbs(REAL(1.0), REAL(1e-14)));
        REQUIRE_THAT(ChebyshevT(0, REAL(1.0)), WithinAbs(REAL(1.0), REAL(1e-14)));
        REQUIRE_THAT(ChebyshevT(0, -REAL(1.0)), WithinAbs(REAL(1.0), REAL(1e-14)));

        // T_1(x) = x
        REQUIRE_THAT(ChebyshevT(1, REAL(0.0)), WithinAbs(REAL(0.0), REAL(1e-14)));
        REQUIRE_THAT(ChebyshevT(1, REAL(0.5)), WithinAbs(REAL(0.5), REAL(1e-14)));
        REQUIRE_THAT(ChebyshevT(1, REAL(1.0)), WithinAbs(REAL(1.0), REAL(1e-14)));
        REQUIRE_THAT(ChebyshevT(1, -REAL(1.0)), WithinAbs(-REAL(1.0), REAL(1e-14)));
    }

    TEST_CASE("ChebyshevT - Known polynomial values", "[chebyshev][polynomial]")
    {
        // T_2(x) = 2x² - 1
        REQUIRE_THAT(ChebyshevT(2, REAL(0.0)), WithinAbs(-REAL(1.0), REAL(1e-14)));
        REQUIRE_THAT(ChebyshevT(2, REAL(0.5)), WithinAbs(2*REAL(0.25) - 1, REAL(1e-14)));
        REQUIRE_THAT(ChebyshevT(2, REAL(1.0)), WithinAbs(REAL(1.0), REAL(1e-14)));

        // T_3(x) = 4x³ - 3x
        REQUIRE_THAT(ChebyshevT(3, REAL(0.0)), WithinAbs(REAL(0.0), REAL(1e-14)));
        REQUIRE_THAT(ChebyshevT(3, REAL(0.5)), WithinAbs(4*REAL(0.125) - REAL(1.5), REAL(1e-14)));
        REQUIRE_THAT(ChebyshevT(3, REAL(1.0)), WithinAbs(REAL(1.0), REAL(1e-14)));

        // T_4(x) = 8x⁴ - 8x² + 1
        REQUIRE_THAT(ChebyshevT(4, REAL(0.0)), WithinAbs(REAL(1.0), REAL(1e-14)));
        REQUIRE_THAT(ChebyshevT(4, REAL(1.0)), WithinAbs(REAL(1.0), REAL(1e-14)));
        Real x = REAL(0.5);
        Real expected_T4 = 8*x*x*x*x - 8*x*x + 1;
        REQUIRE_THAT(ChebyshevT(4, x), WithinAbs(expected_T4, REAL(1e-14)));
    }

    TEST_CASE("ChebyshevT - Property: T_n(1) = 1", "[chebyshev][polynomial]")
    {
        for (int n = 0; n <= 20; n++)
        {
            REQUIRE_THAT(ChebyshevT(n, REAL(1.0)), WithinAbs(REAL(1.0), REAL(1e-12)));
        }
    }

    TEST_CASE("ChebyshevT - Property: T_n(-1) = (-1)^n", "[chebyshev][polynomial]")
    {
        for (int n = 0; n <= 20; n++)
        {
            Real expected = (n % 2 == 0) ? REAL(1.0) : -REAL(1.0);
            REQUIRE_THAT(ChebyshevT(n, -REAL(1.0)), WithinAbs(expected, REAL(1e-12)));
        }
    }

    TEST_CASE("ChebyshevT - Property: T_n(cos(theta)) = cos(n*theta)", "[chebyshev][polynomial]")
    {
        const Real pi = Constants::PI;
        
        for (int n = 0; n <= 10; n++)
        {
            for (Real theta : {REAL(0.0), pi/6, pi/4, pi/3, pi/2, 2*pi/3, pi})
            {
                Real x = std::cos(theta);
                Real T_n = ChebyshevT(n, x);
                Real expected = std::cos(n * theta);
                REQUIRE_THAT(T_n, WithinAbs(expected, REAL(1e-12)));
            }
        }
    }

    TEST_CASE("ChebyshevT - Property: |T_n(x)| <= 1 for x in [-1,1]", "[chebyshev][polynomial]")
    {
        for (int n = 0; n <= 15; n++)
        {
            for (int i = 0; i <= 100; i++)
            {
                Real x = -REAL(1.0) + REAL(2.0) * i / REAL(100.0);
                Real T_n = ChebyshevT(n, x);
                REQUIRE(std::abs(T_n) <= REAL(1.0) + 1e-12);
            }
        }
    }

    TEST_CASE("ChebyshevT - Invalid input", "[chebyshev][polynomial]")
    {
        REQUIRE_THROWS_AS(ChebyshevT(-1, REAL(0.5)), std::invalid_argument);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    ///                         CHEBYSHEV POLYNOMIAL U_n TESTS                              ///
    ///////////////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("ChebyshevU - Base cases U_0 and U_1", "[chebyshev][polynomial]")
    {
        // U_0(x) = 1
        REQUIRE_THAT(ChebyshevU(0, REAL(0.0)), WithinAbs(REAL(1.0), REAL(1e-14)));
        REQUIRE_THAT(ChebyshevU(0, REAL(0.5)), WithinAbs(REAL(1.0), REAL(1e-14)));
        REQUIRE_THAT(ChebyshevU(0, REAL(1.0)), WithinAbs(REAL(1.0), REAL(1e-14)));

        // U_1(x) = 2x
        REQUIRE_THAT(ChebyshevU(1, REAL(0.0)), WithinAbs(REAL(0.0), REAL(1e-14)));
        REQUIRE_THAT(ChebyshevU(1, REAL(0.5)), WithinAbs(REAL(1.0), REAL(1e-14)));
        REQUIRE_THAT(ChebyshevU(1, REAL(1.0)), WithinAbs(REAL(2.0), REAL(1e-14)));
    }

    TEST_CASE("ChebyshevU - Known polynomial values", "[chebyshev][polynomial]")
    {
        // U_2(x) = 4x² - 1
        REQUIRE_THAT(ChebyshevU(2, REAL(0.0)), WithinAbs(-REAL(1.0), REAL(1e-14)));
        REQUIRE_THAT(ChebyshevU(2, REAL(0.5)), WithinAbs(4*REAL(0.25) - 1, REAL(1e-14)));
        REQUIRE_THAT(ChebyshevU(2, REAL(1.0)), WithinAbs(REAL(3.0), REAL(1e-14)));

        // U_3(x) = 8x³ - 4x
        REQUIRE_THAT(ChebyshevU(3, REAL(0.0)), WithinAbs(REAL(0.0), REAL(1e-14)));
        REQUIRE_THAT(ChebyshevU(3, REAL(1.0)), WithinAbs(REAL(4.0), REAL(1e-14)));
    }

    TEST_CASE("ChebyshevU - Property: U_n(1) = n + 1", "[chebyshev][polynomial]")
    {
        for (int n = 0; n <= 20; n++)
        {
            REQUIRE_THAT(ChebyshevU(n, REAL(1.0)), WithinAbs(static_cast<Real>(n + 1), 1e-10));
        }
    }

    TEST_CASE("ChebyshevU - Property: U_n(-1) = (-1)^n * (n + 1)", "[chebyshev][polynomial]")
    {
        for (int n = 0; n <= 20; n++)
        {
            Real expected = ((n % 2 == 0) ? REAL(1.0) : -REAL(1.0)) * (n + 1);
            REQUIRE_THAT(ChebyshevU(n, -REAL(1.0)), WithinAbs(expected, REAL(1e-10)));
        }
    }

    TEST_CASE("ChebyshevU - Invalid input", "[chebyshev][polynomial]")
    {
        REQUIRE_THROWS_AS(ChebyshevU(-1, REAL(0.5)), std::invalid_argument);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    ///                         CHEBYSHEV ROOTS AND EXTREMA                                 ///
    ///////////////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("ChebyshevRoots - T_n(root) = 0", "[chebyshev][roots]")
    {
        for (int n = 1; n <= 10; n++)
        {
            Vector<Real> roots = ChebyshevRoots(n);
            REQUIRE(roots.size() == n);
            
            for (int k = 0; k < n; k++)
            {
                Real root = roots[k];
                REQUIRE(root > -REAL(1.0));
                REQUIRE(root < REAL(1.0));
                REQUIRE_THAT(ChebyshevT(n, root), WithinAbs(REAL(0.0), REAL(1e-12)));
            }
        }
    }

    TEST_CASE("ChebyshevExtrema - |T_n(extremum)| = 1", "[chebyshev][extrema]")
    {
        for (int n = 1; n <= 10; n++)
        {
            Vector<Real> extrema = ChebyshevExtrema(n);
            REQUIRE(extrema.size() == n + 1);
            
            // Check endpoints
            REQUIRE_THAT(extrema[0], WithinAbs(REAL(1.0), REAL(1e-14)));
            REQUIRE_THAT(extrema[n], WithinAbs(-REAL(1.0), REAL(1e-14)));
            
            for (int k = 0; k <= n; k++)
            {
                Real ext = extrema[k];
                Real T_n = ChebyshevT(n, ext);
                REQUIRE_THAT(std::abs(T_n), WithinAbs(REAL(1.0), REAL(1e-12)));
            }
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    ///                       CHEBYSHEV APPROXIMATION TESTS                                 ///
    ///////////////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("ChebyshevApproximation - Constant function", "[chebyshev][approximation]")
    {
        // Approximate f(x) = 5 on [0, 1]
        auto f = [](Real) { return REAL(5.0); };
        ChebyshevApproximation approx(f, REAL(0.0), REAL(1.0), 10);

        // Only c_0 should be significant
        REQUIRE_THAT(approx.Coefficient(0), WithinAbs(REAL(10.0), REAL(1e-12))); // c_0 = 2*f (due to normalization)
        for (int j = 1; j < 10; j++)
        {
            REQUIRE_THAT(approx.Coefficient(j), WithinAbs(REAL(0.0), REAL(1e-12)));
        }

        // Evaluation should return 5
        REQUIRE_THAT(approx(REAL(0.0)), WithinAbs(REAL(5.0), REAL(1e-12)));
        REQUIRE_THAT(approx(REAL(0.5)), WithinAbs(REAL(5.0), REAL(1e-12)));
        REQUIRE_THAT(approx(REAL(1.0)), WithinAbs(REAL(5.0), REAL(1e-12)));
    }

    TEST_CASE("ChebyshevApproximation - Linear function", "[chebyshev][approximation]")
    {
        // Approximate f(x) = 2x + 3 on [-1, 1]
        auto f = [](Real x) { return REAL(2.0)*x + REAL(3.0); };
        ChebyshevApproximation approx(f, -REAL(1.0), REAL(1.0), 10);

        // On [-1,1], f(x) = 2x + 3 has Chebyshev expansion c_0 = 6, c_1 = 2
        REQUIRE_THAT(approx.Coefficient(0), WithinAbs(REAL(6.0), REAL(1e-12)));
        REQUIRE_THAT(approx.Coefficient(1), WithinAbs(REAL(2.0), REAL(1e-12)));
        for (int j = 2; j < 10; j++)
        {
            REQUIRE_THAT(approx.Coefficient(j), WithinAbs(REAL(0.0), REAL(1e-12)));
        }

        // Test evaluation
        REQUIRE_THAT(approx(-REAL(1.0)), WithinAbs(REAL(1.0), REAL(1e-12)));
        REQUIRE_THAT(approx(REAL(0.0)), WithinAbs(REAL(3.0), REAL(1e-12)));
        REQUIRE_THAT(approx(REAL(1.0)), WithinAbs(REAL(5.0), REAL(1e-12)));
    }

    TEST_CASE("ChebyshevApproximation - sin(x) on [0, pi]", "[chebyshev][approximation]")
    {
        const Real pi = Constants::PI;
        auto f = [](Real x) { return std::sin(x); };
        ChebyshevApproximation approx(f, REAL(0.0), pi, 20);

        // Test accuracy at various points
        REQUIRE_THAT(approx(REAL(0.0)), WithinAbs(REAL(0.0), REAL(1e-10)));
        REQUIRE_THAT(approx(pi/2), WithinAbs(REAL(1.0), REAL(1e-10)));
        REQUIRE_THAT(approx(pi), WithinAbs(REAL(0.0), REAL(1e-10)));

        // Maximum error should be small
        Real maxErr = approx.MaxError(f);
        REQUIRE(maxErr < 1e-10);
    }

    TEST_CASE("ChebyshevApproximation - exp(x) on [-1, 1]", "[chebyshev][approximation]")
    {
        auto f = [](Real x) { return std::exp(x); };
        ChebyshevApproximation approx(f, -REAL(1.0), REAL(1.0), 20);

        // Test accuracy
        REQUIRE_THAT(approx(-REAL(1.0)), WithinAbs(std::exp(-REAL(1.0)), 1e-10));
        REQUIRE_THAT(approx(REAL(0.0)), WithinAbs(REAL(1.0), REAL(1e-10)));
        REQUIRE_THAT(approx(REAL(1.0)), WithinAbs(std::exp(REAL(1.0)), 1e-10));

        // Coefficients should decay rapidly for analytic functions
        Real maxErr = approx.MaxError(f);
        REQUIRE(maxErr < 1e-12);
    }

    TEST_CASE("ChebyshevApproximation - Runge function (Chebyshev advantage)", "[chebyshev][approximation]")
    {
        // Runge function: 1/(1 + 25x²)
        // This function is notorious for Runge's phenomenon with equidistant interpolation,
        // but Chebyshev approximation handles it well with enough terms.
        auto runge = [](Real x) { return REAL(1.0) / (REAL(1.0) + REAL(25.0)*x*x); };
        
        // Use more terms for this challenging function
        ChebyshevApproximation approx(runge, -REAL(1.0), REAL(1.0), 80);

        // Maximum error should be reasonably small
        Real maxErr = approx.MaxError(runge);
        REQUIRE(maxErr < 1e-6);  // Much better than equidistant polynomial interpolation

        // Test at critical points
        REQUIRE_THAT(approx(REAL(0.0)), WithinAbs(REAL(1.0), REAL(1e-6)));
        REQUIRE_THAT(approx(REAL(0.2)), WithinAbs(runge(REAL(0.2)), REAL(1e-6)));
    }

    TEST_CASE("ChebyshevApproximation - Domain mapping", "[chebyshev][approximation]")
    {
        // Test on different domains
        auto f = [](Real x) { return x*x; };  // f(x) = x²

        // On [0, 2]
        ChebyshevApproximation approx1(f, REAL(0.0), REAL(2.0), 10);
        REQUIRE_THAT(approx1(REAL(0.0)), WithinAbs(REAL(0.0), REAL(1e-12)));
        REQUIRE_THAT(approx1(REAL(1.0)), WithinAbs(REAL(1.0), REAL(1e-12)));
        REQUIRE_THAT(approx1(REAL(2.0)), WithinAbs(REAL(4.0), REAL(1e-12)));

        // On [10, 20]
        ChebyshevApproximation approx2(f, REAL(10.0), REAL(20.0), 10);
        REQUIRE_THAT(approx2(REAL(10.0)), WithinAbs(REAL(100.0), REAL(1e-10)));
        REQUIRE_THAT(approx2(REAL(15.0)), WithinAbs(REAL(225.0), REAL(1e-10)));
        REQUIRE_THAT(approx2(REAL(20.0)), WithinAbs(REAL(400.0), REAL(1e-10)));
    }

    TEST_CASE("ChebyshevApproximation - Coefficient convergence", "[chebyshev][approximation]")
    {
        // For analytic functions, coefficients decay exponentially
        auto f = [](Real x) { return std::exp(x); };
        ChebyshevApproximation approx(f, -REAL(1.0), REAL(1.0), 30);

        // Check that later coefficients are much smaller
        Real c0 = std::abs(approx.Coefficient(0));
        Real c10 = std::abs(approx.Coefficient(10));
        Real c20 = std::abs(approx.Coefficient(20));

        REQUIRE(c10 < c0 * 1e-3);
        REQUIRE(c20 < c10 * 1e-3);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    ///                       DERIVATIVE TESTS                                              ///
    ///////////////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("ChebyshevApproximation::Derivative - sin -> cos", "[chebyshev][derivative]")
    {
        const Real pi = Constants::PI;
        auto sin_f = [](Real x) { return std::sin(x); };
        auto cos_f = [](Real x) { return std::cos(x); };

        ChebyshevApproximation sin_approx(sin_f, REAL(0.0), 2*pi, 30);
        ChebyshevApproximation deriv = sin_approx.Derivative();

        // Derivative of sin should approximate cos
        for (Real x = REAL(0.1); x < 2*pi - REAL(0.1); x += REAL(0.2))
        {
            REQUIRE_THAT(deriv(x), WithinAbs(cos_f(x), REAL(1e-8)));
        }
    }

    TEST_CASE("ChebyshevApproximation::Derivative - exp -> exp", "[chebyshev][derivative]")
    {
        auto exp_f = [](Real x) { return std::exp(x); };

        ChebyshevApproximation exp_approx(exp_f, -REAL(1.0), REAL(1.0), 25);
        ChebyshevApproximation deriv = exp_approx.Derivative();

        // Derivative of exp should be exp
        for (Real x = -REAL(0.9); x < REAL(0.9); x += REAL(0.1))
        {
            REQUIRE_THAT(deriv(x), WithinAbs(exp_f(x), REAL(1e-10)));
        }
    }

    TEST_CASE("ChebyshevApproximation::Derivative - polynomial", "[chebyshev][derivative]")
    {
        // f(x) = x³ - 2x + 1, f'(x) = 3x² - 2
        auto f = [](Real x) { return x*x*x - 2*x + 1; };
        auto df = [](Real x) { return 3*x*x - 2; };

        ChebyshevApproximation approx(f, -REAL(1.0), REAL(1.0), 10);
        ChebyshevApproximation deriv = approx.Derivative();

        for (Real x = -REAL(0.9); x < REAL(0.9); x += REAL(0.1))
        {
            REQUIRE_THAT(deriv(x), WithinAbs(df(x), REAL(1e-10)));
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    ///                       INTEGRAL TESTS                                                ///
    ///////////////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("ChebyshevApproximation::Integral - cos -> sin", "[chebyshev][integral]")
    {
        const Real pi = Constants::PI;
        auto cos_f = [](Real x) { return std::cos(x); };

        ChebyshevApproximation cos_approx(cos_f, REAL(0.0), pi, 30);
        ChebyshevApproximation integ = cos_approx.Integral();

        // Integral of cos is sin (plus constant)
        // At x = 0, integral should be 0 (by construction)
        Real C = integ(REAL(0.0));  // Should be close to 0
        
        for (Real x = REAL(0.1); x < pi - REAL(0.1); x += REAL(0.2))
        {
            // sin(x) - sin(0) = sin(x)
            REQUIRE_THAT(integ(x) - C, WithinAbs(std::sin(x), 1e-8));
        }
    }

    TEST_CASE("ChebyshevApproximation::Integral then Derivative = original", "[chebyshev][calculus]")
    {
        auto f = [](Real x) { return std::sin(x) * std::exp(-x/2); };

        ChebyshevApproximation approx(f, REAL(0.0), REAL(5.0), 40);
        ChebyshevApproximation integ = approx.Integral();
        ChebyshevApproximation recovered = integ.Derivative();

        // Derivative of integral should recover original
        for (Real x = REAL(0.1); x < REAL(4.9); x += REAL(0.2))
        {
            REQUIRE_THAT(recovered(x), WithinAbs(f(x), REAL(1e-7)));
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    ///                       POLYNOMIAL CONVERSION TESTS                                   ///
    ///////////////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("ChebyshevApproximation::ToPolynomial - T_n expansions", "[chebyshev][conversion]")
    {
        // T_2(x) = 2x² - 1, coefficients: c_0 = -1, c_1 = 0, c_2 = 2
        // But as Chebyshev approximation of T_2, we have c_0 = 0, c_1 = 0, c_2 = 1 (approximately)
        
        // Let's test converting a polynomial through Chebyshev
        auto f = [](Real x) { return 2*x*x - 1; };  // This is T_2
        ChebyshevApproximation approx(f, -REAL(1.0), REAL(1.0), 5);
        
        Polynom<Real> poly = approx.ToPolynomial();

        // The polynomial should evaluate to the same values as original
        for (Real x = -REAL(1.0); x <= REAL(1.0); x += REAL(0.1))
        {
            REQUIRE_THAT(poly(x), WithinAbs(f(x), REAL(1e-10)));
        }
    }

    TEST_CASE("ChebyshevApproximation::ToPolynomial - accuracy", "[chebyshev][conversion]")
    {
        auto f = [](Real x) { return std::exp(x); };
        ChebyshevApproximation approx(f, -REAL(1.0), REAL(1.0), 15);
        
        Polynom<Real> poly = approx.ToPolynomial();

        // Polynomial should give similar results to Chebyshev evaluation
        for (Real x = -REAL(0.9); x <= REAL(0.9); x += REAL(0.1))
        {
            REQUIRE_THAT(poly(x), WithinAbs(approx(x), REAL(1e-10)));
        }
    }

    TEST_CASE("ChebyshevApproximation::FromPolynomial - round trip", "[chebyshev][conversion]")
    {
        // Create a polynomial
        Polynom<Real> p({REAL(1.0), -REAL(2.0), REAL(3.0), -REAL(1.0)});  // 1 - 2x + 3x² - x³

        // Convert to Chebyshev and back
        ChebyshevApproximation approx = ChebyshevApproximation::FromPolynomial(p, -REAL(1.0), REAL(1.0));
        Polynom<Real> p_back = approx.ToPolynomial();

        // Should be close to original
        for (Real x = -REAL(1.0); x <= REAL(1.0); x += REAL(0.1))
        {
            REQUIRE_THAT(p_back(x), WithinAbs(p(x), REAL(1e-10)));
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    ///                       TRUNCATION TESTS                                              ///
    ///////////////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("ChebyshevApproximation::Truncate", "[chebyshev][truncation]")
    {
        auto f = [](Real x) { return std::exp(x); };
        ChebyshevApproximation approx(f, -REAL(1.0), REAL(1.0), 50);

        // Truncate based on threshold
        int original_terms = approx.NumTerms();
        approx.Truncate(1e-14);
        int truncated_terms = approx.NumTerms();

        // Should have fewer terms
        REQUIRE(truncated_terms <= original_terms);
        REQUIRE(truncated_terms > 0);

        // Should still be accurate
        Real maxErr = approx.MaxError(f);
        REQUIRE(maxErr < 1e-12);
    }

    TEST_CASE("ChebyshevApproximation::SetNumTerms", "[chebyshev][truncation]")
    {
        auto f = [](Real x) { return std::sin(x); };
        ChebyshevApproximation approx(f, REAL(0.0), Constants::PI, 30);

        // Use fewer terms
        approx.SetNumTerms(10);
        REQUIRE(approx.NumTerms() == 10);

        // Should still approximate reasonably
        REQUIRE_THAT(approx(Constants::PI/2), WithinAbs(REAL(1.0), REAL(1e-5)));
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    ///                       CONSTRUCTOR AND VALIDATION TESTS                              ///
    ///////////////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("ChebyshevApproximation - Invalid construction", "[chebyshev][validation]")
    {
        auto f = [](Real x) { return x; };

        // Invalid n
        REQUIRE_THROWS_AS(ChebyshevApproximation(f, REAL(0.0), REAL(1.0), 0), std::invalid_argument);
        REQUIRE_THROWS_AS(ChebyshevApproximation(f, REAL(0.0), REAL(1.0), -5), std::invalid_argument);

        // Invalid domain (a >= b)
        REQUIRE_THROWS_AS(ChebyshevApproximation(f, REAL(1.0), REAL(1.0), 10), std::invalid_argument);
        REQUIRE_THROWS_AS(ChebyshevApproximation(f, REAL(2.0), REAL(1.0), 10), std::invalid_argument);

        // Empty coefficients
        Vector<Real> empty;
        REQUIRE_THROWS_AS(ChebyshevApproximation(empty, REAL(0.0), REAL(1.0)), std::invalid_argument);
    }

    TEST_CASE("ChebyshevApproximation - Out of range evaluation", "[chebyshev][validation]")
    {
        auto f = [](Real x) { return x; };
        ChebyshevApproximation approx(f, REAL(0.0), REAL(1.0), 10);

        // Evaluation outside domain should throw
        REQUIRE_THROWS_AS(approx(-REAL(0.1)), std::domain_error);
        REQUIRE_THROWS_AS(approx(REAL(1.1)), std::domain_error);

        // Edge cases should be fine
        REQUIRE_NOTHROW(approx(REAL(0.0)));
        REQUIRE_NOTHROW(approx(REAL(1.0)));
    }

    TEST_CASE("ChebyshevApproximation - Accessors", "[chebyshev][accessors]")
    {
        auto f = [](Real x) { return x*x; };
        ChebyshevApproximation approx(f, -REAL(2.0), REAL(3.0), 15);

        REQUIRE(approx.Degree() == 14);
        REQUIRE(approx.NumCoefficients() == 15);
        REQUIRE(approx.NumTerms() == 15);
        REQUIRE_THAT(approx.DomainMin(), WithinAbs(-REAL(2.0), REAL(1e-14)));
        REQUIRE_THAT(approx.DomainMax(), WithinAbs(REAL(3.0), REAL(1e-14)));

        // Coefficient access
        REQUIRE_NOTHROW(approx.Coefficient(0));
        REQUIRE_NOTHROW(approx.Coefficient(14));
        REQUIRE_THROWS_AS(approx.Coefficient(-1), std::out_of_range);
        REQUIRE_THROWS_AS(approx.Coefficient(15), std::out_of_range);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    ///                       COPY AND MOVE TESTS                                           ///
    ///////////////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("ChebyshevApproximation - Copy construction", "[chebyshev][copy]")
    {
        auto f = [](Real x) { return std::sin(x); };
        ChebyshevApproximation original(f, REAL(0.0), Constants::PI, 20);

        ChebyshevApproximation copy(original);

        REQUIRE(copy.NumCoefficients() == original.NumCoefficients());
        REQUIRE_THAT(copy(REAL(1.0)), WithinAbs(original(REAL(1.0)), REAL(1e-14)));
    }

    TEST_CASE("ChebyshevApproximation - Move construction", "[chebyshev][move]")
    {
        auto f = [](Real x) { return std::sin(x); };
        ChebyshevApproximation original(f, REAL(0.0), Constants::PI, 20);
        Real val_at_1 = original(REAL(1.0));

        ChebyshevApproximation moved(std::move(original));

        REQUIRE_THAT(moved(REAL(1.0)), WithinAbs(val_at_1, REAL(1e-14)));
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    ///                       CHEBYSHEV BASIS TESTS (First Kind)                            ///
    ///////////////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("ChebyshevBasis - Evaluate matches ChebyshevT", "[chebyshev][basis]")
    {
        ChebyshevBasis basis;
        
        for (int n = 0; n <= 10; n++)
        {
            for (Real x : {-REAL(1.0), -REAL(0.5), REAL(0.0), REAL(0.5), REAL(1.0)})
            {
                REQUIRE_THAT(basis.Evaluate(n, x), WithinAbs(ChebyshevT(n, x), REAL(1e-14)));
            }
        }
    }

    TEST_CASE("ChebyshevBasis - Weight function", "[chebyshev][basis]")
    {
        ChebyshevBasis basis;
        
        // w(x) = 1/sqrt(1 - x^2)
        REQUIRE_THAT(basis.WeightFunction(REAL(0.0)), WithinAbs(REAL(1.0), REAL(1e-14)));
        REQUIRE_THAT(basis.WeightFunction(REAL(0.5)), WithinAbs(REAL(1.0) / std::sqrt(REAL(0.75)), REAL(1e-12)));
        
        // At x = 0.8, w(0.8) = 1/sqrt(1 - 0.64) = 1/sqrt(0.36) = 1/0.6 ≈ 1.667
        REQUIRE_THAT(basis.WeightFunction(REAL(0.8)), WithinAbs(REAL(1.0) / REAL(0.6), REAL(1e-12)));
    }

    TEST_CASE("ChebyshevBasis - Normalization", "[chebyshev][basis]")
    {
        ChebyshevBasis basis;
        const Real pi = Constants::PI;
        
        // ||T_0||^2 = pi
        REQUIRE_THAT(basis.Normalization(0), WithinAbs(pi, REAL(1e-14)));
        
        // ||T_n||^2 = pi/2 for n > 0
        for (int n = 1; n <= 10; n++)
        {
            REQUIRE_THAT(basis.Normalization(n), WithinAbs(pi / REAL(2.0), REAL(1e-14)));
        }
    }

    TEST_CASE("ChebyshevBasis - Domain", "[chebyshev][basis]")
    {
        ChebyshevBasis basis;
        
        REQUIRE_THAT(basis.DomainMin(), WithinAbs(-REAL(1.0), REAL(1e-14)));
        REQUIRE_THAT(basis.DomainMax(), WithinAbs(REAL(1.0), REAL(1e-14)));
    }

    TEST_CASE("ChebyshevBasis - Recurrence coefficients", "[chebyshev][basis]")
    {
        ChebyshevBasis basis;
        Real a, b, c;
        
        // T_{n+1} = 2x*T_n - T_{n-1}, so a=2, b=0, c=-1
        for (int n = 0; n <= 10; n++)
        {
            basis.RecurrenceCoefficients(n, a, b, c);
            REQUIRE_THAT(a, WithinAbs(REAL(2.0), REAL(1e-14)));
            REQUIRE_THAT(b, WithinAbs(REAL(0.0), REAL(1e-14)));
            REQUIRE_THAT(c, WithinAbs(-REAL(1.0), REAL(1e-14)));
        }
    }

    TEST_CASE("ChebyshevBasis - GetExtrema", "[chebyshev][basis]")
    {
        ChebyshevBasis basis;
        
        // For T_n, extrema are at x_k = cos(k*pi/n) for k=0,1,...,n
        // At these points |T_n| = 1
        for (int n = 1; n <= 5; n++)
        {
            auto extrema = basis.GetExtrema(n);
            REQUIRE(extrema.size() == n + 1);
            
            for (int k = 0; k <= n; k++)
            {
                Real expected = std::cos(k * Constants::PI / n);
                REQUIRE_THAT(extrema[k], WithinAbs(expected, REAL(1e-14)));
                
                // |T_n(extremum)| = 1
                REQUIRE_THAT(std::abs(ChebyshevT(n, extrema[k])), WithinAbs(REAL(1.0), REAL(1e-12)));
            }
        }
    }

    TEST_CASE("ChebyshevBasis - GetZeros", "[chebyshev][basis]")
    {
        ChebyshevBasis basis;
        
        // Zeros of T_n: x_k = cos((2k-1)*pi/(2n)) for k=1,...,n
        for (int n = 1; n <= 5; n++)
        {
            auto zeros = basis.GetZeros(n);
            REQUIRE(zeros.size() == n);
            
            for (int k = 0; k < n; k++)
            {
                // Verify T_n(zero) = 0
                REQUIRE_THAT(ChebyshevT(n, zeros[k]), WithinAbs(REAL(0.0), REAL(1e-12)));
            }
        }
    }

    TEST_CASE("ChebyshevBasis - Invalid inputs", "[chebyshev][basis][validation]")
    {
        ChebyshevBasis basis;
        Real a, b, c;
        
        REQUIRE_THROWS_AS(basis.Evaluate(-1, REAL(0.5)), std::invalid_argument);
        REQUIRE_THROWS_AS(basis.Normalization(-1), std::invalid_argument);
        REQUIRE_THROWS_AS(basis.RecurrenceCoefficients(-1, a, b, c), std::invalid_argument);
        REQUIRE_THROWS_AS(basis.GetExtrema(-1), std::invalid_argument);
        REQUIRE_THROWS_AS(basis.GetZeros(0), std::invalid_argument);
        REQUIRE_THROWS_AS(basis.GetZeros(-1), std::invalid_argument);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    ///                       CHEBYSHEV BASIS TESTS (Second Kind)                           ///
    ///////////////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("ChebyshevBasisSecondKind - Evaluate matches ChebyshevU", "[chebyshev][basis]")
    {
        ChebyshevBasisSecondKind basis;
        
        for (int n = 0; n <= 10; n++)
        {
            for (Real x : {-REAL(1.0), -REAL(0.5), REAL(0.0), REAL(0.5), REAL(1.0)})
            {
                REQUIRE_THAT(basis.Evaluate(n, x), WithinAbs(ChebyshevU(n, x), REAL(1e-12)));
            }
        }
    }

    TEST_CASE("ChebyshevBasisSecondKind - Weight function", "[chebyshev][basis]")
    {
        ChebyshevBasisSecondKind basis;
        
        // w(x) = sqrt(1 - x^2)
        REQUIRE_THAT(basis.WeightFunction(REAL(0.0)), WithinAbs(REAL(1.0), REAL(1e-14)));
        REQUIRE_THAT(basis.WeightFunction(REAL(0.5)), WithinAbs(std::sqrt(REAL(0.75)), REAL(1e-12)));
        
        // At boundaries w = 0
        REQUIRE_THAT(basis.WeightFunction(REAL(1.0)), WithinAbs(REAL(0.0), REAL(1e-14)));
        REQUIRE_THAT(basis.WeightFunction(-REAL(1.0)), WithinAbs(REAL(0.0), REAL(1e-14)));
    }

    TEST_CASE("ChebyshevBasisSecondKind - Normalization", "[chebyshev][basis]")
    {
        ChebyshevBasisSecondKind basis;
        const Real pi = Constants::PI;
        
        // ||U_n||^2 = pi/2 for all n >= 0
        for (int n = 0; n <= 10; n++)
        {
            REQUIRE_THAT(basis.Normalization(n), WithinAbs(pi / REAL(2.0), REAL(1e-14)));
        }
    }

    TEST_CASE("ChebyshevBasisSecondKind - Domain", "[chebyshev][basis]")
    {
        ChebyshevBasisSecondKind basis;
        
        REQUIRE_THAT(basis.DomainMin(), WithinAbs(-REAL(1.0), REAL(1e-14)));
        REQUIRE_THAT(basis.DomainMax(), WithinAbs(REAL(1.0), REAL(1e-14)));
    }

    TEST_CASE("ChebyshevBasisSecondKind - Recurrence coefficients", "[chebyshev][basis]")
    {
        ChebyshevBasisSecondKind basis;
        Real a, b, c;
        
        // U_{n+1} = 2x*U_n - U_{n-1}, so a=2, b=0, c=-1 (same as T_n)
        for (int n = 0; n <= 10; n++)
        {
            basis.RecurrenceCoefficients(n, a, b, c);
            REQUIRE_THAT(a, WithinAbs(REAL(2.0), REAL(1e-14)));
            REQUIRE_THAT(b, WithinAbs(REAL(0.0), REAL(1e-14)));
            REQUIRE_THAT(c, WithinAbs(-REAL(1.0), REAL(1e-14)));
        }
    }

    TEST_CASE("ChebyshevBasisSecondKind - GetZeros", "[chebyshev][basis]")
    {
        ChebyshevBasisSecondKind basis;
        
        // Zeros of U_n: x_k = cos(k*pi/(n+1)) for k=1,...,n
        for (int n = 1; n <= 5; n++)
        {
            auto zeros = basis.GetZeros(n);
            REQUIRE(zeros.size() == n);
            
            for (int k = 0; k < n; k++)
            {
                // Verify U_n(zero) = 0
                REQUIRE_THAT(ChebyshevU(n, zeros[k]), WithinAbs(REAL(0.0), REAL(1e-12)));
            }
        }
    }

    TEST_CASE("ChebyshevBasisSecondKind - Invalid inputs", "[chebyshev][basis][validation]")
    {
        ChebyshevBasisSecondKind basis;
        Real a, b, c;
        
        REQUIRE_THROWS_AS(basis.Evaluate(-1, REAL(0.5)), std::invalid_argument);
        REQUIRE_THROWS_AS(basis.Normalization(-1), std::invalid_argument);
        REQUIRE_THROWS_AS(basis.RecurrenceCoefficients(-1, a, b, c), std::invalid_argument);
        REQUIRE_THROWS_AS(basis.GetZeros(0), std::invalid_argument);
        REQUIRE_THROWS_AS(basis.GetZeros(-1), std::invalid_argument);
    }

} // namespace MML::Tests::Base::ChebyshevTests

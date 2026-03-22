///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML) Tests                            ///
///                                                                                   ///
///  File:        function_helpers_tests.cpp                                          ///
///  Description: Tests for FunctionHelpers.h - function composition, transformations,///
///               arithmetic operations, and derivative wrappers                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "core/FunctionHelpers.h"
#endif

using namespace MML;
using namespace MML::Testing;

namespace MML::Tests::Core::FunctionHelpersTests
{
    // Test functions for use in tests
    class LinearFunc : public IRealFunction {
        Real _a, _b;  // f(x) = a*x + b
    public:
        LinearFunc(Real a, Real b) : _a(a), _b(b) {}
        Real operator()(Real x) const override { return _a * x + _b; }
    };

    class QuadraticFunc : public IRealFunction {
        Real _a, _b, _c;  // f(x) = a*x² + b*x + c
    public:
        QuadraticFunc(Real a, Real b, Real c) : _a(a), _b(b), _c(c) {}
        Real operator()(Real x) const override { return _a * x * x + _b * x + _c; }
    };

    class SinFunc : public IRealFunction {
    public:
        Real operator()(Real x) const override { return std::sin(x); }
    };

    class CosFunc : public IRealFunction {
    public:
        Real operator()(Real x) const override { return std::cos(x); }
    };

    class ExpFunc : public IRealFunction {
    public:
        Real operator()(Real x) const override { return std::exp(x); }
    };

    ///////////////////////////////////////////////////////////////////////////////////////////
    ///                        TAYLOR SERIES EXPANSION TESTS                                ///
    ///////////////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("FunctionHelpers::TaylorSeries2 - quadratic function", "[taylor][exact]")
    {
        // For a quadratic function, Taylor series of degree 2 should be exact
        QuadraticFunc f(2.0, -3.0, 1.0);  // f(x) = 2x^2 - 3x + 1
        Real a = 0.0;  // Expand around 0

        auto taylor = TaylorSeries2(f, a);

        // Taylor approximation should be exact for quadratic
        // f(0) = 1, f(1) = 0, f(-1) = 6, f(0.5) = -0.25
        REQUIRE_THAT(taylor(0.0), RealWithinAbs(f(0.0), 1e-6));
        REQUIRE_THAT(taylor(1.0), RealWithinAbs(f(1.0), 1e-6));    // f(1) = 0, use Abs
        REQUIRE_THAT(taylor(-1.0), RealWithinRel(f(-1.0), 1e-6));
        REQUIRE_THAT(taylor(0.5), RealWithinAbs(f(0.5), 1e-6));
    }

    TEST_CASE("FunctionHelpers::TaylorSeries2 - exp function", "[taylor]")
    {
        ExpFunc f;
        Real a = 0.0;  // Expand around 0

        auto taylor = TaylorSeries2(f, a);

        // e^x ≈ 1 + x + x²/2 near x=0
        REQUIRE_THAT(taylor(0.0), RealWithinRel(f(0.0), 1e-8));
        REQUIRE_THAT(taylor(0.1), RealWithinRel(f(0.1), 1e-3));  // Good near expansion point
        REQUIRE_THAT(taylor(0.5), RealWithinRel(f(0.5), 0.05));  // Acceptable further out
    }

    TEST_CASE("FunctionHelpers::TaylorSeries3 - cubic function", "[taylor][exact]")
    {
        // For a cubic function, Taylor series of degree 3 should be exact
        auto cubic = [](Real x) { return x * x * x - 2 * x * x + 3 * x - 1; };
        class CubicFunc : public IRealFunction {
        public:
            Real operator()(Real x) const override { return x * x * x - 2 * x * x + 3 * x - 1; }
        };
        CubicFunc f;
        Real a = 0.0;

        auto taylor = TaylorSeries3(f, a);

        REQUIRE_THAT(taylor(0.0), RealWithinRel(f(0.0), 1e-8));
        REQUIRE_THAT(taylor(1.0), RealWithinRel(f(1.0), 1e-4));
        REQUIRE_THAT(taylor(-0.5), RealWithinRel(f(-0.5), 1e-4));
    }

    TEST_CASE("FunctionHelpers::TaylorSeries3 - sin function", "[taylor]")
    {
        SinFunc f;
        Real a = 0.0;

        auto taylor = TaylorSeries3(f, a);

        // sin(x) ≈ x - x³/6 near x=0
        REQUIRE_THAT(taylor(0.0), RealWithinAbs(0.0, 1e-10));
        REQUIRE_THAT(taylor(0.1), RealWithinRel(f(0.1), 1e-5));
        REQUIRE_THAT(taylor(0.5), RealWithinRel(f(0.5), 1e-2));
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    ///                      FIRST DERIVATIVE WRAPPER TESTS                                 ///
    ///////////////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("FunctionHelpers::RealFuncDerived - linear function", "[derivative]")
    {
        // Derivative of 3x + 2 should be 3
        LinearFunc f(3.0, 2.0);

        RealFuncDerived1 d1(f);
        RealFuncDerived2 d2(f);
        RealFuncDerived4 d4(f);
        RealFuncDerived6 d6(f);
        RealFuncDerived8 d8(f);

        REQUIRE_THAT(d1(0.0), RealWithinRel(3.0, 1e-4));
        REQUIRE_THAT(d2(0.0), RealWithinRel(3.0, 1e-8));
        REQUIRE_THAT(d4(0.0), RealWithinRel(3.0, 1e-10));
        REQUIRE_THAT(d6(0.0), RealWithinRel(3.0, 1e-10));
        REQUIRE_THAT(d8(0.0), RealWithinRel(3.0, 1e-10));
    }

    TEST_CASE("FunctionHelpers::RealFuncDerived - sin function", "[derivative]")
    {
        // Derivative of sin(x) should be cos(x)
        SinFunc f;
        CosFunc expected;

        RealFuncDerived2 d2(f);
        RealFuncDerived4 d4(f);
        RealFuncDerived6 d6(f);

        Real x = 1.0;
        REQUIRE_THAT(d2(x), RealWithinRel(expected(x), 1e-8));
        REQUIRE_THAT(d4(x), RealWithinRel(expected(x), 1e-10));
        REQUIRE_THAT(d6(x), RealWithinRel(expected(x), 1e-10));
    }

    TEST_CASE("FunctionHelpers::RealFuncDerived - custom step size", "[derivative]")
    {
        QuadraticFunc f(1.0, 0.0, 0.0);  // f(x) = x²

        RealFuncDerived2 d2_auto(f);
        RealFuncDerived2 d2_custom(f, 1e-4);

        Real x = 2.0;  // f'(2) = 4
        REQUIRE_THAT(d2_auto(x), RealWithinRel(4.0, 1e-8));
        REQUIRE_THAT(d2_custom(x), RealWithinRel(4.0, 1e-6));
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    ///                     SECOND DERIVATIVE WRAPPER TESTS                                 ///
    ///////////////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("FunctionHelpers::RealFuncSecondDerived - quadratic", "[second-derivative]")
    {
        // f(x) = 3x², f''(x) = 6
        QuadraticFunc f(3.0, 0.0, 0.0);

        RealFuncSecondDerived2 d2(f);
        RealFuncSecondDerived4 d4(f);
        RealFuncSecondDerived6 d6(f);

        REQUIRE_THAT(d2(0.0), RealWithinRel(6.0, 1e-4));
        REQUIRE_THAT(d4(0.0), RealWithinRel(6.0, 1e-6));
        REQUIRE_THAT(d6(0.0), RealWithinRel(6.0, 1e-6));
    }

    TEST_CASE("FunctionHelpers::RealFuncSecondDerived - sin function", "[second-derivative]")
    {
        // f(x) = sin(x), f''(x) = -sin(x)
        SinFunc f;

        RealFuncSecondDerived4 d4(f);

        Real x = 1.0;
        REQUIRE_THAT(d4(x), RealWithinRel(-std::sin(x), 1e-6));
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    ///                      FUNCTION ARITHMETIC TESTS                                      ///
    ///////////////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("FunctionHelpers::RealFuncSum", "[arithmetic]")
    {
        LinearFunc f1(2.0, 1.0);   // 2x + 1
        LinearFunc f2(3.0, -2.0);  // 3x - 2

        RealFuncSum sum(f1, f2);   // 5x - 1

        REQUIRE_THAT(sum(0.0), RealWithinAbs(-1.0, 1e-10));
        REQUIRE_THAT(sum(1.0), RealWithinAbs(4.0, 1e-10));
        REQUIRE_THAT(sum(2.0), RealWithinAbs(9.0, 1e-10));
    }

    TEST_CASE("FunctionHelpers::RealFuncDiff", "[arithmetic]")
    {
        LinearFunc f1(5.0, 3.0);   // 5x + 3
        LinearFunc f2(2.0, 1.0);   // 2x + 1

        RealFuncDiff diff(f1, f2); // 3x + 2

        REQUIRE_THAT(diff(0.0), RealWithinAbs(2.0, 1e-10));
        REQUIRE_THAT(diff(1.0), RealWithinAbs(5.0, 1e-10));
        REQUIRE_THAT(diff(2.0), RealWithinAbs(8.0, 1e-10));
    }

    TEST_CASE("FunctionHelpers::RealFuncProduct", "[arithmetic]")
    {
        LinearFunc f1(2.0, 0.0);   // 2x
        LinearFunc f2(3.0, 0.0);   // 3x

        RealFuncProduct prod(f1, f2); // 6x²

        REQUIRE_THAT(prod(0.0), RealWithinAbs(0.0, 1e-10));
        REQUIRE_THAT(prod(1.0), RealWithinAbs(6.0, 1e-10));
        REQUIRE_THAT(prod(2.0), RealWithinAbs(24.0, 1e-10));
    }

    TEST_CASE("FunctionHelpers::RealFuncQuotient", "[arithmetic]")
    {
        QuadraticFunc f1(1.0, 0.0, 0.0);  // x²
        LinearFunc f2(1.0, 0.0);          // x

        RealFuncQuotient quot(f1, f2);    // x

        REQUIRE_THAT(quot(1.0), RealWithinAbs(1.0, 1e-10));
        REQUIRE_THAT(quot(2.0), RealWithinAbs(2.0, 1e-10));
        REQUIRE_THAT(quot(5.0), RealWithinAbs(5.0, 1e-10));
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    ///                      FUNCTION COMPOSITION TESTS                                     ///
    ///////////////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("FunctionHelpers::RealFuncCompose - linear composition", "[composition]")
    {
        LinearFunc f(2.0, 1.0);   // f(x) = 2x + 1
        LinearFunc g(3.0, 0.0);   // g(x) = 3x

        RealFuncCompose fog(f, g); // f(g(x)) = 2(3x) + 1 = 6x + 1

        REQUIRE_THAT(fog(0.0), RealWithinAbs(1.0, 1e-10));
        REQUIRE_THAT(fog(1.0), RealWithinAbs(7.0, 1e-10));
        REQUIRE_THAT(fog(2.0), RealWithinAbs(13.0, 1e-10));
    }

    TEST_CASE("FunctionHelpers::RealFuncCompose - sin(x^2)", "[composition]")
    {
        SinFunc f;
        QuadraticFunc g(1.0, 0.0, 0.0);  // x^2

        RealFuncCompose fog(f, g);  // sin(x^2)

        REQUIRE_THAT(fog(0.0), RealWithinAbs(0.0, 1e-10));
        REQUIRE_THAT(fog(1.0), RealWithinAbs(std::sin(1.0), 1e-10));
        REQUIRE_THAT(fog(2.0), RealWithinAbs(std::sin(4.0), 1e-10));
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    ///                      SCALAR TRANSFORMATION TESTS                                    ///
    ///////////////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("FunctionHelpers::RealFuncScale", "[transform]")
    {
        LinearFunc f(1.0, 0.0);  // f(x) = x

        RealFuncScale scaled(f, 3.0);  // 3x

        REQUIRE_THAT(scaled(0.0), RealWithinAbs(0.0, 1e-10));
        REQUIRE_THAT(scaled(2.0), RealWithinAbs(6.0, 1e-10));
        REQUIRE_THAT(scaled(-1.0), RealWithinAbs(-3.0, 1e-10));
    }

    TEST_CASE("FunctionHelpers::RealFuncShift", "[transform]")
    {
        LinearFunc f(2.0, 0.0);  // f(x) = 2x

        RealFuncShift shifted(f, 5.0);  // 2x + 5

        REQUIRE_THAT(shifted(0.0), RealWithinAbs(5.0, 1e-10));
        REQUIRE_THAT(shifted(1.0), RealWithinAbs(7.0, 1e-10));
        REQUIRE_THAT(shifted(-2.0), RealWithinAbs(1.0, 1e-10));
    }

    TEST_CASE("FunctionHelpers::RealFuncNegate", "[transform]")
    {
        LinearFunc f(3.0, 2.0);  // f(x) = 3x + 2

        RealFuncNegate neg(f);   // -(3x + 2)

        REQUIRE_THAT(neg(0.0), RealWithinAbs(-2.0, 1e-10));
        REQUIRE_THAT(neg(1.0), RealWithinAbs(-5.0, 1e-10));
        REQUIRE_THAT(neg(-1.0), RealWithinAbs(1.0, 1e-10));
    }

    TEST_CASE("FunctionHelpers::RealFuncAbs", "[transform]")
    {
        LinearFunc f(1.0, -2.0);  // f(x) = x - 2

        RealFuncAbs absf(f);      // |x - 2|

        REQUIRE_THAT(absf(0.0), RealWithinAbs(2.0, 1e-10));
        REQUIRE_THAT(absf(2.0), RealWithinAbs(0.0, 1e-10));
        REQUIRE_THAT(absf(4.0), RealWithinAbs(2.0, 1e-10));
    }

    TEST_CASE("FunctionHelpers::RealFuncPow", "[transform]")
    {
        LinearFunc f(1.0, 0.0);  // f(x) = x

        RealFuncPow pow2(f, 2.0);  // x²
        RealFuncPow pow3(f, 3.0);  // x³

        REQUIRE_THAT(pow2(2.0), RealWithinAbs(4.0, 1e-10));
        REQUIRE_THAT(pow2(3.0), RealWithinAbs(9.0, 1e-10));
        REQUIRE_THAT(pow3(2.0), RealWithinAbs(8.0, 1e-10));
        REQUIRE_THAT(pow3(3.0), RealWithinAbs(27.0, 1e-10));
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    ///                      COMPARISON / DISTANCE TESTS                                    ///
    ///////////////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("FunctionHelpers::RealFuncAbsDiff", "[comparison]")
    {
        LinearFunc f1(2.0, 0.0);  // 2x
        LinearFunc f2(1.0, 0.0);  // x

        RealFuncAbsDiff absdiff(f1, f2);  // |2x - x| = |x|

        REQUIRE_THAT(absdiff(0.0), RealWithinAbs(0.0, 1e-10));
        REQUIRE_THAT(absdiff(3.0), RealWithinAbs(3.0, 1e-10));
        REQUIRE_THAT(absdiff(-3.0), RealWithinAbs(3.0, 1e-10));
    }

    TEST_CASE("FunctionHelpers::RealFuncDiffSqr", "[comparison]")
    {
        LinearFunc f1(3.0, 0.0);  // 3x
        LinearFunc f2(1.0, 0.0);  // x

        RealFuncDiffSqr diffsqr(f1, f2);  // (3x - x)² = 4x²

        REQUIRE_THAT(diffsqr(0.0), RealWithinAbs(0.0, 1e-10));
        REQUIRE_THAT(diffsqr(2.0), RealWithinAbs(16.0, 1e-10));
        REQUIRE_THAT(diffsqr(-2.0), RealWithinAbs(16.0, 1e-10));
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    ///                      COMBINED / COMPLEX TESTS                                       ///
    ///////////////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("FunctionHelpers - nested composition", "[composition][complex]")
    {
        // Compose three functions: h(g(f(x)))
        LinearFunc f(2.0, 0.0);  // f(x) = 2x
        LinearFunc g(1.0, 1.0);  // g(x) = x + 1
        LinearFunc h(3.0, 0.0);  // h(x) = 3x

        RealFuncCompose gof(g, f);   // g(f(x)) = 2x + 1
        RealFuncCompose hogof(h, gof); // h(g(f(x))) = 3(2x + 1) = 6x + 3

        REQUIRE_THAT(hogof(0.0), RealWithinAbs(3.0, 1e-10));
        REQUIRE_THAT(hogof(1.0), RealWithinAbs(9.0, 1e-10));
        REQUIRE_THAT(hogof(2.0), RealWithinAbs(15.0, 1e-10));
    }

    TEST_CASE("FunctionHelpers - derivative of composition", "[derivative][composition]")
    {
        // d/dx[sin(x²)] = 2x*cos(x²)
        SinFunc sinf;
        QuadraticFunc sqr(1.0, 0.0, 0.0);

        RealFuncCompose f(sinf, sqr);  // sin(x²)
        RealFuncDerived4 df(f);

        auto expected = [](Real x) { return 2.0 * x * std::cos(x * x); };

        Real x = 1.0;
        REQUIRE_THAT(df(x), RealWithinRel(expected(x), 1e-6));

        x = 0.5;
        REQUIRE_THAT(df(x), RealWithinRel(expected(x), 1e-6));
    }

    TEST_CASE("FunctionHelpers - arithmetic then derivative", "[arithmetic][derivative]")
    {
        // d/dx[sin(x) + cos(x)] = cos(x) - sin(x)
        SinFunc sinf;
        CosFunc cosf;

        RealFuncSum sum(sinf, cosf);
        RealFuncDerived4 dsum(sum);

        auto expected = [](Real x) { return std::cos(x) - std::sin(x); };

        Real x = 1.0;
        REQUIRE_THAT(dsum(x), RealWithinRel(expected(x), 1e-8));
    }

    TEST_CASE("FunctionHelpers - product rule verification", "[arithmetic][derivative]")
    {
        // d/dx[x * sin(x)] = sin(x) + x*cos(x)
        LinearFunc id(1.0, 0.0);  // x
        SinFunc sinf;

        RealFuncProduct prod(id, sinf);
        RealFuncDerived4 dprod(prod);

        auto expected = [](Real x) { return std::sin(x) + x * std::cos(x); };

        Real x = 1.0;
        REQUIRE_THAT(dprod(x), RealWithinRel(expected(x), 1e-6));
    }

    TEST_CASE("FunctionHelpers - function algebra classes", "[FunctionHelpers]")
    {
        LinearFunc f1(2.0, 0.0);
        LinearFunc f2(1.0, 0.0);

        RealFuncDiff diff(f1, f2);
        RealFuncAbsDiff absdiff(f1, f2);
        RealFuncDiffSqr diffsqr(f1, f2);

        REQUIRE_THAT(diff(2.0), RealWithinAbs(2.0, 1e-10));
        REQUIRE_THAT(absdiff(2.0), RealWithinAbs(2.0, 1e-10));
        REQUIRE_THAT(diffsqr(2.0), RealWithinAbs(4.0, 1e-10));
    }

} // namespace

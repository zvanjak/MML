#include <catch2/catch_all.hpp>
#include "../../TestPrecision.h"
#include "../../TestMatchers.h"

#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/MMLBase.h>
#include <mml/core/Integration/MonteCarloIntegration.h>
#include <mml/base/Function.h>
#endif

#include <stdexcept>
#include <vector>

using namespace MML;
using namespace MML::Testing;

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace MML::Tests::Algorithms::MonteCarloIntegrationTests
{
    ///////////////////////////    TEST FUNCTIONS    ///////////////////////////

    // f(x) = 1 (constant) - integral over [0,1]^N = static_cast<MML::Real>(1.0)
    template<int N>
    class ConstantOne : public IScalarFunction<N>
    {
    public:
        Real operator()(const VectorN<Real, N>& x) const override { return static_cast<MML::Real>(1.0); }
    };

    // f(x) = x[0] (linear) - integral over [0,1]^N = static_cast<MML::Real>(0.5)
    template<int N>
    class LinearX0 : public IScalarFunction<N>
    {
    public:
        Real operator()(const VectorN<Real, N>& x) const override { return x[0]; }
    };

    // f(x) = x^2 (1D) - integral over [0,1] = 1/3
    class XSquared1D : public IScalarFunction<1>
    {
    public:
        Real operator()(const VectorN<Real, 1>& x) const override { return x[0] * x[0]; }
    };

    // f(x,y) = x*y - integral over [0,1]^2 = static_cast<MML::Real>(0.25)
    class ProductXY : public IScalarFunction<2>
    {
    public:
        Real operator()(const VectorN<Real, 2>& x) const override { return x[0] * x[1]; }
    };

    // f(x,y,z) = x*y*z - integral over [0,1]^3 = static_cast<MML::Real>(0.125)
    class ProductXYZ : public IScalarFunction<3>
    {
    public:
        Real operator()(const VectorN<Real, 3>& x) const override { return x[0] * x[1] * x[2]; }
    };

    // f(x) = sin(x) - integral over [0,π] = static_cast<MML::Real>(2.0)
    class SinFunc1D : public IScalarFunction<1>
    {
    public:
        Real operator()(const VectorN<Real, 1>& x) const override { return std::sin(x[0]); }
    };

    // f(x) = exp(-x^2) - integral over [-inf,inf] = sqrt(π)
    // But we integrate over finite domain [a,b]
    class Gaussian1D : public IScalarFunction<1>
    {
    public:
        Real operator()(const VectorN<Real, 1>& x) const override { 
            return std::exp(-x[0] * x[0]); 
        }
    };

    // f(x,y) = exp(-(x^2+y^2)) - 2D Gaussian
    class Gaussian2D : public IScalarFunction<2>
    {
    public:
        Real operator()(const VectorN<Real, 2>& x) const override {
            return std::exp(-(x[0]*x[0] + x[1]*x[1]));
        }
    };

    // Sum of coordinates: f(x) = sum(x_i)
    template<int N>
    class SumOfCoords : public IScalarFunction<N>
    {
    public:
        Real operator()(const VectorN<Real, N>& x) const override {
            Real sum = static_cast<MML::Real>(0.0);
            for (int i = 0; i < N; ++i) sum += x[i];
            return sum;
        }
    };

    ///////////////////////////    BASIC FUNCTIONALITY TESTS    ///////////////////////////

    TEST_CASE("MonteCarlo_1D_Constant", "[montecarlo][integration][1d]")
    {
        ConstantOne<1> func;
        MonteCarloIntegrator<1> mc(42);  // Fixed seed for reproducibility
        
        VectorN<Real, 1> lower{static_cast<MML::Real>(0.0)}, upper{static_cast<MML::Real>(1.0)};
        
        auto result = mc.integrate(func, lower, upper, MonteCarloConfig().samples(50000));
        
        REQUIRE_THAT(result.value, WithinAbs(static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.01)));
        REQUIRE(result.samples_used == 50000);
    }

    TEST_CASE("MonteCarlo_1D_Linear", "[montecarlo][integration][1d]")
    {
        LinearX0<1> func;
        MonteCarloIntegrator<1> mc(42);
        
        VectorN<Real, 1> lower{static_cast<MML::Real>(0.0)}, upper{static_cast<MML::Real>(1.0)};
        
        // Integral of x from 0 to 1 = static_cast<MML::Real>(0.5)
        auto result = mc.integrate(func, lower, upper, MonteCarloConfig().samples(100000));
        
        REQUIRE_THAT(result.value, WithinAbs(static_cast<MML::Real>(0.5), static_cast<MML::Real>(0.01)));
    }

    TEST_CASE("MonteCarlo_1D_XSquared", "[montecarlo][integration][1d]")
    {
        XSquared1D func;
        MonteCarloIntegrator<1> mc(42);
        
        VectorN<Real, 1> lower{static_cast<MML::Real>(0.0)}, upper{static_cast<MML::Real>(1.0)};
        
        // Integral of x^2 from 0 to 1 = 1/3
        auto result = mc.integrate(func, lower, upper, MonteCarloConfig().samples(100000));
        
        REQUIRE_THAT(result.value, WithinAbs(static_cast<MML::Real>(1.0)/static_cast<MML::Real>(3.0), static_cast<MML::Real>(0.01)));
    }

    TEST_CASE("MonteCarlo_1D_Sin_0_to_Pi", "[montecarlo][integration][1d]")
    {
        SinFunc1D func;
        MonteCarloIntegrator<1> mc(42);
        
        VectorN<Real, 1> lower{static_cast<MML::Real>(0.0)}, upper{static_cast<MML::Real>(MML::Constants::PI)};
        
        // Integral of sin(x) from 0 to π = static_cast<MML::Real>(2.0)
        auto result = mc.integrate(func, lower, upper, MonteCarloConfig().samples(100000));
        
        REQUIRE_THAT(result.value, WithinAbs(static_cast<MML::Real>(2.0), static_cast<MML::Real>(0.02)));
    }

    ///////////////////////////    2D INTEGRATION TESTS    ///////////////////////////

    TEST_CASE("MonteCarlo_2D_Constant", "[montecarlo][integration][2d]")
    {
        ConstantOne<2> func;
        MonteCarloIntegrator<2> mc(42);
        
        VectorN<Real, 2> lower{static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0)}, upper{static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0)};
        
        auto result = mc.integrate(func, lower, upper, MonteCarloConfig().samples(50000));
        
        REQUIRE_THAT(result.value, WithinAbs(static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.01)));
    }

    TEST_CASE("MonteCarlo_2D_ProductXY", "[montecarlo][integration][2d]")
    {
        ProductXY func;
        MonteCarloIntegrator<2> mc(42);
        
        VectorN<Real, 2> lower{static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0)}, upper{static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0)};
        
        // Integral of x*y over [0,1]^2 = static_cast<MML::Real>(0.25)
        auto result = mc.integrate(func, lower, upper, MonteCarloConfig().samples(100000));
        
        REQUIRE_THAT(result.value, WithinAbs(static_cast<MML::Real>(0.25), static_cast<MML::Real>(0.01)));
    }

    TEST_CASE("MonteCarlo_2D_RectangularDomain", "[montecarlo][integration][2d]")
    {
        ConstantOne<2> func;
        MonteCarloIntegrator<2> mc(42);
        
        VectorN<Real, 2> lower{static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0)}, upper{static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0)};
        
        // Volume = 2 * 3 = 6
        auto result = mc.integrate(func, lower, upper, MonteCarloConfig().samples(50000));
        
        REQUIRE_THAT(result.value, WithinAbs(static_cast<MML::Real>(6.0), static_cast<MML::Real>(0.05)));
    }

    ///////////////////////////    3D INTEGRATION TESTS    ///////////////////////////

    TEST_CASE("MonteCarlo_3D_Constant_UnitCube", "[montecarlo][integration][3d]")
    {
        ConstantOne<3> func;
        MonteCarloIntegrator<3> mc(42);
        
        VectorN<Real, 3> lower{static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0)}, upper{static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0)};
        
        auto result = mc.integrate(func, lower, upper, MonteCarloConfig().samples(50000));
        
        REQUIRE_THAT(result.value, WithinAbs(static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.01)));
    }

    TEST_CASE("MonteCarlo_3D_ProductXYZ", "[montecarlo][integration][3d]")
    {
        ProductXYZ func;
        MonteCarloIntegrator<3> mc(42);
        
        VectorN<Real, 3> lower{static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0)}, upper{static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0)};
        
        // Integral of x*y*z over [0,1]^3 = static_cast<MML::Real>(0.125)
        auto result = mc.integrate(func, lower, upper, MonteCarloConfig().samples(100000));
        
        REQUIRE_THAT(result.value, WithinAbs(static_cast<MML::Real>(0.125), static_cast<MML::Real>(0.005)));
    }

    ///////////////////////////    HIGHER DIMENSIONS    ///////////////////////////

    TEST_CASE("MonteCarlo_4D_SumOfCoords", "[montecarlo][integration][4d]")
    {
        SumOfCoords<4> func;
        MonteCarloIntegrator<4> mc(42);
        
        VectorN<Real, 4> lower{static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0)}, upper{static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0)};
        
        // Integral of (x1+x2+x3+x4) over [0,1]^4 = 4 * static_cast<MML::Real>(0.5) = static_cast<MML::Real>(2.0)
        auto result = mc.integrate(func, lower, upper, MonteCarloConfig().samples(100000));
        
        REQUIRE_THAT(result.value, WithinAbs(static_cast<MML::Real>(2.0), static_cast<MML::Real>(0.02)));
    }

    TEST_CASE("MonteCarlo_5D_Constant", "[montecarlo][integration][5d]")
    {
        ConstantOne<5> func;
        MonteCarloIntegrator<5> mc(42);
        
        VectorN<Real, 5> lower, upper;
        for (int i = 0; i < 5; ++i) {
            lower[i] = static_cast<MML::Real>(0.0);
            upper[i] = static_cast<MML::Real>(1.0);
        }
        
        auto result = mc.integrate(func, lower, upper, MonteCarloConfig().samples(50000));
        
        REQUIRE_THAT(result.value, WithinAbs(static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.02)));
    }

    ///////////////////////////    PI ESTIMATION (CLASSIC EXAMPLE)    ///////////////////////////

    TEST_CASE("MonteCarlo_EstimatePi", "[montecarlo][integration][pi][classic]")
    {
        // Classic textbook example: estimate π using hit-or-miss
        auto result = EstimatePi(500000, 42);
        
        // π ≈ static_cast<MML::Real>(3.14159)...
        REQUIRE_THAT(result.value, WithinAbs(static_cast<MML::Real>(MML::Constants::PI), static_cast<MML::Real>(0.01)));
        
        // Error estimate should be reasonable
        REQUIRE(result.error_estimate > static_cast<MML::Real>(0.0));
        REQUIRE(result.error_estimate < static_cast<MML::Real>(0.02));
    }

    ///////////////////////////    UNIT BALL VOLUME    ///////////////////////////

    TEST_CASE("MonteCarlo_UnitBallVolume_2D", "[montecarlo][integration][ball]")
    {
        // 2D "ball" = circle, volume = π
        auto result = EstimateUnitBallVolume<2>(500000, 42);
        
        REQUIRE_THAT(result.value, WithinAbs(static_cast<MML::Real>(MML::Constants::PI), static_cast<MML::Real>(0.02)));
    }

    TEST_CASE("MonteCarlo_UnitBallVolume_3D", "[montecarlo][integration][ball]")
    {
        // 3D ball volume = (4/3)π ≈ static_cast<MML::Real>(4.189)
        auto result = EstimateUnitBallVolume<3>(500000, 42);
        
        Real expected = (static_cast<MML::Real>(4.0)/static_cast<MML::Real>(3.0)) * static_cast<MML::Real>(MML::Constants::PI);
        REQUIRE_THAT(result.value, WithinAbs(expected, static_cast<MML::Real>(0.05)));
    }

    TEST_CASE("MonteCarlo_UnitBallVolume_4D", "[montecarlo][integration][ball]")
    {
        // 4D ball volume = (π²/2) ≈ static_cast<MML::Real>(4.935)
        auto result = EstimateUnitBallVolume<4>(500000, 42);
        
        Real expected = (static_cast<MML::Real>(MML::Constants::PI) * static_cast<MML::Real>(MML::Constants::PI)) / static_cast<MML::Real>(2.0);
        REQUIRE_THAT(result.value, WithinAbs(expected, static_cast<MML::Real>(0.08)));
    }

    ///////////////////////////    HIT-OR-MISS INTEGRATOR    ///////////////////////////

    TEST_CASE("HitOrMiss_CircleArea", "[montecarlo][hitorrmiss]")
    {
        HitOrMissIntegrator<2> integrator(42);
        
        VectorN<Real, 2> lower{static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0)}, upper{static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0)};
        
        // Quarter circle in first quadrant: x² + y² ≤ 1
        auto inside_quarter_circle = [](const VectorN<Real, 2>& p) {
            return p[0]*p[0] + p[1]*p[1] <= static_cast<MML::Real>(1.0);
        };
        
        auto result = integrator.estimateVolume(inside_quarter_circle, lower, upper, 500000);
        
        // Area of quarter circle = π/4 ≈ static_cast<MML::Real>(0.785)
        REQUIRE_THAT(result.value, WithinAbs(static_cast<MML::Real>(MML::Constants::PI)/static_cast<MML::Real>(4.0), static_cast<MML::Real>(0.01)));
    }

    TEST_CASE("HitOrMiss_TriangleArea", "[montecarlo][hitorrmiss]")
    {
        HitOrMissIntegrator<2> integrator(42);
        
        VectorN<Real, 2> lower{static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0)}, upper{static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0)};
        
        // Triangle: x + y ≤ 1
        auto inside_triangle = [](const VectorN<Real, 2>& p) {
            return p[0] + p[1] <= static_cast<MML::Real>(1.0);
        };
        
        auto result = integrator.estimateVolume(inside_triangle, lower, upper, 100000);
        
        // Area of triangle = static_cast<MML::Real>(0.5)
        REQUIRE_THAT(result.value, WithinAbs(static_cast<MML::Real>(0.5), static_cast<MML::Real>(0.01)));
    }

    ///////////////////////////    ANTITHETIC VARIATES    ///////////////////////////

    TEST_CASE("MonteCarlo_AntitheticVariates", "[montecarlo][variance_reduction]")
    {
        LinearX0<1> func;
        
        VectorN<Real, 1> lower{static_cast<MML::Real>(0.0)}, upper{static_cast<MML::Real>(1.0)};
        
        // Without antithetic
        MonteCarloIntegrator<1> mc1(42);
        auto result1 = mc1.integrate(func, lower, upper, 
                                     MonteCarloConfig().samples(10000).antithetic(false));
        
        // With antithetic variates (should have lower variance for monotonic functions)
        MonteCarloIntegrator<1> mc2(42);
        auto result2 = mc2.integrate(func, lower, upper, 
                                     MonteCarloConfig().samples(10000).antithetic(true));
        
        // Both should give correct answer
        REQUIRE_THAT(result1.value, WithinAbs(static_cast<MML::Real>(0.5), static_cast<MML::Real>(0.02)));
        REQUIRE_THAT(result2.value, WithinAbs(static_cast<MML::Real>(0.5), static_cast<MML::Real>(0.02)));
        
        // Note: antithetic typically has lower variance for monotonic functions
        // but this is a statistical property, not guaranteed for every run
    }

    TEST_CASE("MonteCarlo_RejectsInvalidConfiguration", "[montecarlo][validation]")
    {
        ConstantOne<1> func;
        MonteCarloIntegrator<1> mc(42);

        VectorN<Real, 1> lower{static_cast<MML::Real>(0.0)}, upper{static_cast<MML::Real>(1.0)};

        REQUIRE_THROWS_AS(mc.integrate(func, lower, upper, MonteCarloConfig().samples(0)), ArgumentError);
        REQUIRE_THROWS_AS(mc.integrate(func, lower, upper, MonteCarloConfig().samples(1).antithetic(true)), ArgumentError);
    }

    TEST_CASE("MonteCarlo_RejectsInvalidBounds", "[montecarlo][validation]")
    {
        ConstantOne<2> func;
        MonteCarloIntegrator<2> mc(42);

        VectorN<Real, 2> lower{static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0)}, upper{static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0)};
        VectorN<Real, 2> reversed{static_cast<MML::Real>(-1.0), static_cast<MML::Real>(1.0)};
        VectorN<Real, 2> nonFinite{static_cast<MML::Real>(1.0), std::numeric_limits<Real>::infinity()};

        REQUIRE_THROWS_AS(mc.integrate(func, upper, reversed, MonteCarloConfig().samples(10)), ArgumentError);
        REQUIRE_THROWS_AS(mc.integrate(func, lower, nonFinite, MonteCarloConfig().samples(10)), ArgumentError);
    }

    ///////////////////////////    STRATIFIED SAMPLING    ///////////////////////////

    TEST_CASE("Stratified_2D_ProductXY", "[montecarlo][stratified]")
    {
        ProductXY func;
        StratifiedMonteCarloIntegrator<2> mc(42);
        
        VectorN<Real, 2> lower{static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0)}, upper{static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0)};
        
        // 10 strata per dimension, 10 samples per stratum = 1000 total
        auto result = mc.integrate(func, lower, upper, 10, 10, 42);
        
        // Integral of x*y over [0,1]^2 = static_cast<MML::Real>(0.25)
        REQUIRE_THAT(result.value, WithinAbs(static_cast<MML::Real>(0.25), static_cast<MML::Real>(0.01)));
    }

    TEST_CASE("Stratified_3D_Constant", "[montecarlo][stratified]")
    {
        ConstantOne<3> func;
        StratifiedMonteCarloIntegrator<3> mc(42);
        
        VectorN<Real, 3> lower{static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0)}, upper{static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0)};
        
        // 5 strata per dimension, 5 samples per stratum = 625 total
        auto result = mc.integrate(func, lower, upper, 5, 5, 42);
        
        REQUIRE_THAT(result.value, WithinAbs(static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.02)));
    }

            TEST_CASE("Stratified_RejectsInvalidConfiguration", "[montecarlo][stratified][validation]")
            {
                ConstantOne<2> func;
                StratifiedMonteCarloIntegrator<2> mc(42);

                VectorN<Real, 2> lower{static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0)}, upper{static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0)};

                REQUIRE_THROWS_AS(mc.integrate(func, lower, upper, 0, 10), ArgumentError);
                REQUIRE_THROWS_AS(mc.integrate(func, lower, upper, 10, 0), ArgumentError);
                REQUIRE_THROWS_AS(mc.integrate(func, upper, lower, 10, 10), ArgumentError);
            }

            TEST_CASE("HitOrMiss_RejectsInvalidConfiguration", "[montecarlo][hitorrmiss][validation]")
            {
                HitOrMissIntegrator<2> integrator(42);
                auto indicator = [](const VectorN<Real, 2>&) { return true; };

                VectorN<Real, 2> lower{static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0)}, upper{static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0)};

                REQUIRE_THROWS_AS(integrator.estimateVolume(indicator, lower, upper, 0), ArgumentError);
                REQUIRE_THROWS_AS(integrator.estimateVolume(indicator, upper, lower, 10), ArgumentError);
            }

    ///////////////////////////    CONVENIENCE FUNCTION TESTS    ///////////////////////////

    TEST_CASE("IntegrateMonteCarlo1D_Convenience", "[montecarlo][1d][convenience]")
    {
        // Test the convenience wrapper for 1D IRealFunction
        RealFunction func([](Real x) { return x * x; });
        
        auto result = IntegrateMonteCarlo1D(func, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0),
                                            MonteCarloConfig().samples(100000).randomSeed(42));
        
        // Integral of x^2 from 0 to 1 = 1/3
        REQUIRE_THAT(result.value, WithinAbs(static_cast<MML::Real>(1.0)/static_cast<MML::Real>(3.0), static_cast<MML::Real>(0.01)));
    }

    ///////////////////////////    ERROR ESTIMATION TESTS    ///////////////////////////

    TEST_CASE("MonteCarlo_ErrorEstimate_Reasonable", "[montecarlo][error]")
    {
        ConstantOne<2> func;
        MonteCarloIntegrator<2> mc(42);
        
        VectorN<Real, 2> lower{static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0)}, upper{static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0)};
        
        auto result = mc.integrate(func, lower, upper, MonteCarloConfig().samples(10000));
        
        // Error should be positive and small for constant function
        REQUIRE(result.error_estimate >= static_cast<MML::Real>(0.0));
        
        // The actual error should be within a few standard errors
        Real actual_error = std::abs(result.value - static_cast<MML::Real>(1.0));
        REQUIRE(actual_error < 3 * result.error_estimate + static_cast<MML::Real>(0.001));  // 3-sigma + small buffer
    }

    TEST_CASE("MonteCarlo_MoreSamples_LowerError", "[montecarlo][convergence]")
    {
        ProductXY func;
        
        VectorN<Real, 2> lower{static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0)}, upper{static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0)};
        
        // Run with different sample sizes
        MonteCarloIntegrator<2> mc1(42);
        auto result1 = mc1.integrate(func, lower, upper, MonteCarloConfig().samples(1000));
        
        MonteCarloIntegrator<2> mc2(42);
        auto result2 = mc2.integrate(func, lower, upper, MonteCarloConfig().samples(100000));
        
        // More samples should give lower error estimate
        REQUIRE(result2.error_estimate < result1.error_estimate);
    }

    ///////////////////////////    GAUSSIAN INTEGRAL TESTS    ///////////////////////////

    TEST_CASE("MonteCarlo_Gaussian_FiniteDomain", "[montecarlo][gaussian]")
    {
        Gaussian1D func;
        MonteCarloIntegrator<1> mc(42);
        
        // Integrate exp(-x²) from -3 to 3
        // This captures most of the Gaussian (erf(3) ≈ static_cast<MML::Real>(0.9999779))
        VectorN<Real, 1> lower{-static_cast<MML::Real>(3.0)}, upper{static_cast<MML::Real>(3.0)};
        
        auto result = mc.integrate(func, lower, upper, MonteCarloConfig().samples(100000));
        
        // Full integral = √π ≈ static_cast<MML::Real>(1.7724538509), but we get slightly less
        Real expected = std::sqrt(static_cast<MML::Real>(MML::Constants::PI)) * std::erf(static_cast<MML::Real>(3.0));
        REQUIRE_THAT(result.value, WithinAbs(expected, static_cast<MML::Real>(0.02)));
    }

    TEST_CASE("MonteCarlo_2DGaussian", "[montecarlo][gaussian][2d]")
    {
        Gaussian2D func;
        MonteCarloIntegrator<2> mc(42);
        
        // Integrate exp(-(x²+y²)) over [-3,3]²
        VectorN<Real, 2> lower{-static_cast<MML::Real>(3.0), -static_cast<MML::Real>(3.0)}, upper{static_cast<MML::Real>(3.0), static_cast<MML::Real>(3.0)};
        
        auto result = mc.integrate(func, lower, upper, MonteCarloConfig().samples(200000));
        
        // Full 2D Gaussian integral = π, truncated integral ≈ π * erf(3)²
        Real expected = static_cast<MML::Real>(MML::Constants::PI) * std::erf(static_cast<MML::Real>(3.0)) * std::erf(static_cast<MML::Real>(3.0));
        REQUIRE_THAT(result.value, WithinAbs(expected, static_cast<MML::Real>(0.05)));
    }

    ///////////////////////////    ALGORITHM STATUS TESTS    ///////////////////////////

    TEST_CASE("MonteCarlo_AlgorithmStatus_Success", "[montecarlo][status]")
    {
        ConstantOne<2> func;
        MonteCarloIntegrator<2> mc(42);
        
        VectorN<Real, 2> lower{static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0)}, upper{static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0)};
        auto result = mc.integrate(func, lower, upper, MonteCarloConfig().samples(1000));
        
        REQUIRE(result.status == AlgorithmStatus::Success);
        REQUIRE(result.IsSuccess());
        REQUIRE(result.algorithm_name == "MonteCarloIntegrator");
        REQUIRE(result.elapsed_time_ms >= static_cast<MML::Real>(0.0));
        REQUIRE(result.function_evaluations == 1000);
    }

    TEST_CASE("Stratified_AlgorithmStatus_Success", "[montecarlo][stratified][status]")
    {
        ConstantOne<2> func;
        StratifiedMonteCarloIntegrator<2> smc(42);
        
        VectorN<Real, 2> lower{static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0)}, upper{static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0)};
        auto result = smc.integrate(func, lower, upper, 10, 10);
        
        REQUIRE(result.status == AlgorithmStatus::Success);
        REQUIRE(result.IsSuccess());
        REQUIRE(result.algorithm_name == "StratifiedMonteCarloIntegrator");
        REQUIRE(result.elapsed_time_ms >= static_cast<MML::Real>(0.0));
        REQUIRE(result.function_evaluations > 0);
    }

    TEST_CASE("HitOrMiss_AlgorithmStatus_Completed", "[montecarlo][hitorrmiss][status]")
    {
        auto indicator = [](const VectorN<Real, 2>& x) -> bool {
            return x[0]*x[0] + x[1]*x[1] <= static_cast<MML::Real>(1.0);
        };
        
        HitOrMissIntegrator<2> hom(42);
        
        VectorN<Real, 2> lower{-static_cast<MML::Real>(1.0), -static_cast<MML::Real>(1.0)}, upper{static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0)};
        auto result = hom.estimateVolume(indicator, lower, upper, 10000);
        
        // Fixed sample count: no tolerance was checked -> Completed, not Success (RESULT_POLICY)
        REQUIRE(result.status == AlgorithmStatus::Completed);
        REQUIRE(result.IsSuccess());
        REQUIRE_FALSE(result.converged);
        REQUIRE(HasEstimate(result.error_estimate));                 // the statistical estimate is real
        REQUIRE(result.algorithm_name == "HitOrMissIntegrator");
        REQUIRE(result.elapsed_time_ms >= static_cast<MML::Real>(0.0));
        REQUIRE(result.function_evaluations == 10000);
    }


    TEST_CASE("MonteCarlo variance survives a large constant offset (Welford)", "[montecarlo][variance][wewi.1.12]")
    {
        // Resolve x in the samples while making E[f^2] - E[f]^2 lose the variance.
        class OffsetLinear : public IScalarFunction<1> {
            Real _offset;
        public:
            explicit OffsetLinear(Real offset) : _offset(offset) {}
            Real operator()(const VectorN<Real, 1>& p) const override { return _offset + p[0]; }
        };
        const Real offset = TOL(1e8, std::ldexp(Real(1), std::numeric_limits<Real>::digits / 2));
        OffsetLinear func(offset);
        REQUIRE(func(VectorN<Real, 1>{0}) != func(VectorN<Real, 1>{1}));
        MonteCarloIntegrator<1> mc(42);
        VectorN<Real, 1> lower{static_cast<MML::Real>(0.0)}, upper{static_cast<MML::Real>(1.0)};

        auto result = mc.integrate(func, lower, upper, MonteCarloConfig().samples(20000));
        INFO("variance " << result.variance << " error " << result.error_estimate);
        REQUIRE_THAT(result.variance, WithinRel(static_cast<MML::Real>(1.0) / static_cast<MML::Real>(12.0), static_cast<MML::Real>(0.1)));
        REQUIRE(result.error_estimate > static_cast<MML::Real>(0.0));
        REQUIRE_THAT(result.value, WithinAbs(offset + Real(0.5), Real(0.05)));

        // Welford accumulator unit check against the closed form for 1..n
        Detail::WelfordAccumulator acc;
        const Real sequenceOffset = TOL(1e9, std::ldexp(Real(1), std::numeric_limits<Real>::digits - 3));
        for (int i = 1; i <= 100; ++i) acc.add(Real(i) + sequenceOffset);
        REQUIRE_THAT(acc.mean, WithinRel(sequenceOffset + Real(50.5), TOL(1e-12, Real(2) * Precision::Eps)));
        REQUIRE_THAT(acc.sampleVariance(), WithinRel(Real(841.6666666666666), TOL(1e-9, Real(4) * Precision::Eps)));
    }

    TEST_CASE("Stratified MC uncertainty matches independent within-stratum statistics", "[montecarlo][stratified][kko3]")
    {
        for (int strata : {1, 10}) {
            const int samples = 1000 / strata;
            std::vector<Real> values;
            ScalarFunctionFromStdFunc<1> linear([&](const VectorN<Real, 1>& x) {
                values.push_back(x[0]);
                return x[0];
            });
            StratifiedMonteCarloIntegrator<1> integrator;
            const auto result = integrator.integrate(linear, VectorN<Real, 1>{0}, VectorN<Real, 1>{1}, strata, samples, 42);
            REQUIRE(values.size() == 1000);
            long double expectedMean = 0;
            long double expectedVariance = 0;
            for (int s = 0; s < strata; ++s) {
                long double sum = 0;
                for (int i = 0; i < samples; ++i)
                    sum += values[s * samples + i];
                const long double mean = sum / samples;
                long double squaredDeviations = 0;
                for (int i = 0; i < samples; ++i) {
                    const long double delta = static_cast<long double>(values[s * samples + i]) - mean;
                    squaredDeviations += delta * delta;
                }
                expectedMean += mean / strata;
                expectedVariance += squaredDeviations / (samples - 1) / samples / strata / strata;
            }
            CHECK_THAT(result.value, WithinRel(static_cast<Real>(expectedMean), Precision::EqualityTolerance));
            CHECK_THAT(result.variance, WithinRel(static_cast<Real>(expectedVariance), Precision::EqualityTolerance));
            CHECK_THAT(result.error_estimate, WithinRel(static_cast<Real>(std::sqrt(expectedVariance)), Precision::EqualityTolerance));
            CHECK(result.error_estimate > 0);
            CHECK(std::abs(result.value - Real(0.5)) <= Real(6) * result.error_estimate);
            CHECK(result.samples_used == values.size());
            CHECK(result.function_evaluations == static_cast<int>(values.size()));
            CHECK(result.error_estimated);
            CHECK_FALSE(result.converged);
            CHECK(result.status == AlgorithmStatus::MaxIterationsExceeded);
            CHECK_FALSE(result.IsUsable());
            CHECK_FALSE(result.error_message.empty());
        }
    }

    TEST_CASE("Stratified MC excludes between-stratum variation from uncertainty", "[montecarlo][stratified][kko3]")
    {
        int calls = 0;
        ScalarFunctionFromStdFunc<1> piecewise([&](const VectorN<Real, 1>& x) {
            ++calls;
            return x[0] < Real(0.5) ? Real(-1) : Real(1);
        });
        const auto result = StratifiedMonteCarloIntegrator<1>().integrate(
            piecewise, VectorN<Real, 1>{0}, VectorN<Real, 1>{1}, 2, 10, 42);
        CHECK(result.value == 0);
        CHECK(result.variance == 0);
        CHECK(result.error_estimate == 0);
        CHECK(result.error_estimated);
        CHECK(result.converged);
        CHECK(result.status == AlgorithmStatus::Success);
        CHECK(result.IsUsable());
        CHECK(result.error_message.empty());
        CHECK(calls == 20);
        CHECK(result.samples_used == 20);
        CHECK(result.function_evaluations == calls);
    }

    TEST_CASE("Stratified MC zero and tiny means cannot bypass uncertainty", "[montecarlo][stratified][kko3]")
    {
        for (Real amplitude : {Real(1), Precision::DivisionSafetyThreshold / Real(2)}) {
            int calls = 0;
            ScalarFunctionFromStdFunc<1> alternating([&](const VectorN<Real, 1>&) {
                return ++calls % 2 ? -amplitude : amplitude;
            });
            const auto result = StratifiedMonteCarloIntegrator<1>().integrate(
                alternating, VectorN<Real, 1>{0}, VectorN<Real, 1>{1}, 1, 2, 42);
            CHECK(result.value == 0);
            CHECK_THAT(result.error_estimate, WithinRel(amplitude, Precision::EqualityTolerance));
            CHECK_FALSE(result.converged);
            CHECK(result.status == AlgorithmStatus::MaxIterationsExceeded);
            CHECK_FALSE(result.error_message.empty());
            CHECK(calls == 2);
            CHECK(result.function_evaluations == calls);
            CHECK(result.samples_used == 2);
        }
    }

    TEST_CASE("Stratified MC single observations retain values without invented uncertainty", "[montecarlo][stratified][kko3]")
    {
        for (int strata : {1, 4}) {
            std::vector<Real> values;
            ScalarFunctionFromStdFunc<1> linear([&](const VectorN<Real, 1>& x) {
                values.push_back(x[0]);
                return x[0];
            });
            const auto result = StratifiedMonteCarloIntegrator<1>().integrate(
                linear, VectorN<Real, 1>{0}, VectorN<Real, 1>{1}, strata, 1, 42);
            Real expected = 0;
            for (Real value : values)
                expected += value / static_cast<Real>(strata);
            CHECK_THAT(result.value, WithinRel(expected, Precision::EqualityTolerance));
            CHECK(result.status == AlgorithmStatus::Completed);
            CHECK_FALSE(result.converged);
            CHECK_FALSE(result.error_estimated);
            CHECK_FALSE(HasEstimate(result.error_estimate));
            CHECK_FALSE(HasEstimate(result.variance));
            CHECK(result.IsUsable());
            CHECK(result.error_message.empty());
            CHECK(values.size() == static_cast<std::size_t>(strata));
            CHECK(result.samples_used == values.size());
            CHECK(result.function_evaluations == strata);
        }
    }

    TEST_CASE("Stratified MC checks grid and diagnostic budgets before callbacks", "[montecarlo][stratified][kko3]")
    {
        int calls = 0;
        ScalarFunctionFromStdFunc<64> guard64([&](const VectorN<Real, 64>&) -> Real {
            ++calls;
            throw std::logic_error("oversized grid reached the callback");
        });
        VectorN<Real, 64> lower64, upper64;
        for (int d = 0; d < 64; ++d) {
            lower64[d] = 0;
            upper64[d] = 1;
        }
        CHECK_THROWS_AS(StratifiedMonteCarloIntegrator<64>().integrate(guard64, lower64, upper64, 2, 2), ResourceLimitError);
        CHECK(calls == 0);

        ScalarFunctionFromStdFunc<2> guard2([&](const VectorN<Real, 2>&) -> Real {
            ++calls;
            throw std::logic_error("oversized work reached the callback");
        });
        CHECK_THROWS_AS(StratifiedMonteCarloIntegrator<2>().integrate(
            guard2, VectorN<Real, 2>{0, 0}, VectorN<Real, 2>{1, 1}, 46341, 1), ResourceLimitError);
        CHECK(calls == 0);
        ScalarFunctionFromStdFunc<1> guard1([&](const VectorN<Real, 1>&) -> Real {
            ++calls;
            throw std::logic_error("counter boundary reached the callback");
        });
        CHECK_THROWS_AS(StratifiedMonteCarloIntegrator<1>().integrate(
            guard1, VectorN<Real, 1>{0}, VectorN<Real, 1>{1}, 2, std::numeric_limits<int>::max()), ResourceLimitError);
        CHECK(calls == 0);
        calls = 0;
        CHECK_THROWS_AS(StratifiedMonteCarloIntegrator<1>().integrate(
            guard1, VectorN<Real, 1>{0}, VectorN<Real, 1>{1}, 1, std::numeric_limits<int>::max()), std::logic_error);
        CHECK(calls == 1);
    }

    TEST_CASE("Stratified MC scales values and standard errors by the domain volume", "[montecarlo][stratified][kko3]")
    {
        std::vector<Real> unitValues, scaledValues;
        auto run = [](Real width, std::vector<Real>& values) {
            ScalarFunctionFromStdFunc<2> linear([&](const VectorN<Real, 2>& x) {
                const Real value = x[0] / width + Real(2) * x[1];
                values.push_back(value);
                return value;
            });
            return StratifiedMonteCarloIntegrator<2>().integrate(
                linear, VectorN<Real, 2>{0, 0}, VectorN<Real, 2>{width, 1}, 2, 8, 123);
        };
        const auto unit = run(1, unitValues);
        const auto scaled = run(4, scaledValues);
        CHECK(unitValues == scaledValues);
        CHECK(scaled.value == Real(4) * unit.value);
        CHECK(scaled.error_estimate == Real(4) * unit.error_estimate);
        CHECK(scaled.variance == unit.variance);
        CHECK(unit.samples_used == 32);
        CHECK(scaled.samples_used == 32);
        CHECK(unit.function_evaluations == 32);
        CHECK(scaled.function_evaluations == 32);
        CHECK(unit.status == AlgorithmStatus::MaxIterationsExceeded);
        CHECK(scaled.status == unit.status);
    }

    TEST_CASE("Stratified MC retains finite large constants and rejects nonfinite samples", "[montecarlo][stratified][kko3]")
    {
        const Real large = std::numeric_limits<Real>::max() * Real(0.75);
        int calls = 0;
        ScalarFunctionFromStdFunc<1> constant([&](const VectorN<Real, 1>&) {
            ++calls;
            return large;
        });
        const auto result = StratifiedMonteCarloIntegrator<1>().integrate(
            constant, VectorN<Real, 1>{0}, VectorN<Real, 1>{Real(1) / large}, 2, 4, 42);
        CHECK_THAT(result.value, WithinAbs(Real(1), Precision::EqualityTolerance));
        CHECK(result.error_estimate == 0);
        CHECK(result.status == AlgorithmStatus::Success);
        CHECK(result.converged);
        CHECK(result.samples_used == 8);
        CHECK(result.function_evaluations == calls);
        for (Real invalid : {std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity()}) {
            calls = 0;
            ScalarFunctionFromStdFunc<1> nonfinite([&](const VectorN<Real, 1>&) {
                ++calls;
                return invalid;
            });
            CHECK_THROWS_AS(StratifiedMonteCarloIntegrator<1>().integrate(
                nonfinite, VectorN<Real, 1>{0}, VectorN<Real, 1>{1}, 2, 4), NumericalMethodError);
            CHECK(calls == 1);
        }
    }

    TEST_CASE("Serial antithetic averaging does not overflow a representable integral", "[montecarlo][antithetic][kko3]")
    {
        const Real large = std::numeric_limits<Real>::max() * Real(0.75);
        int calls = 0;
        ScalarFunctionFromStdFunc<1> constant([&](const VectorN<Real, 1>&) {
            ++calls;
            return large;
        });
        MonteCarloConfig config;
        config.num_samples = 4;
        config.use_antithetic = true;
        config.seed = 42;
        const auto result = MonteCarloIntegrator<1>().integrate(
            constant, VectorN<Real, 1>{0}, VectorN<Real, 1>{Real(1) / large}, config);
        CHECK_THAT(result.value, WithinAbs(Real(1), Precision::EqualityTolerance));
        CHECK(result.error_estimate == 0);
        CHECK(result.variance == 0);
        CHECK(result.converged);
        CHECK(result.status == AlgorithmStatus::Success);
        CHECK(calls == 4);
        CHECK(result.samples_used == 4);
        CHECK(result.function_evaluations == calls);
    }

    TEST_CASE("Serial antithetic reflections remain inside extreme finite bounds", "[montecarlo][antithetic][kko3]")
    {
        for (Real sign : {Real(-1), Real(1)}) {
            const Real first = sign * std::numeric_limits<Real>::max() * Real(0.6);
            const Real second = sign * std::numeric_limits<Real>::max() * Real(0.7);
            const Real lower = std::min(first, second);
            const Real upper = std::max(first, second);
            const Real width = upper - lower;
            std::vector<Real> points;
            ScalarFunctionFromStdFunc<1> constant([&](const VectorN<Real, 1>& x) {
                points.push_back(x[0]);
                CHECK(std::isfinite(x[0]));
                CHECK(x[0] >= lower);
                CHECK(x[0] <= upper);
                return Real(1) / width;
            });
            MonteCarloConfig config;
            config.num_samples = 4;
            config.use_antithetic = true;
            config.seed = 42;
            const auto result = MonteCarloIntegrator<1>().integrate(
                constant, VectorN<Real, 1>{lower}, VectorN<Real, 1>{upper}, config);
            REQUIRE(points.size() == 4);
            for (std::size_t i = 0; i < points.size(); i += 2) {
                const Real reflectedSum = (points[i] - lower) / width + (points[i + 1] - lower) / width;
                CHECK_THAT(reflectedSum, WithinAbs(Real(1), Precision::EqualityTolerance));
            }
            CHECK_THAT(result.value, WithinAbs(Real(1), Precision::EqualityTolerance));
            CHECK(result.status == AlgorithmStatus::Success);
            CHECK(result.function_evaluations == 4);
        }
    }

    TEST_CASE("Stratified MC zero-volume domains perform no sampling", "[montecarlo][stratified][kko3]")
    {
        int calls = 0;
        ScalarFunctionFromStdFunc<2> guard([&](const VectorN<Real, 2>&) -> Real {
            ++calls;
            throw std::logic_error("zero-volume domain reached the callback");
        });
        const auto result = StratifiedMonteCarloIntegrator<2>().integrate(
            guard, VectorN<Real, 2>{0, 0}, VectorN<Real, 2>{0, 1}, 2, 4, 42);
        CHECK(calls == 0);
        CHECK(result.value == 0);
        CHECK(result.error_estimate == 0);
        CHECK(result.error_estimated);
        CHECK_FALSE(result.converged);
        CHECK(result.status == AlgorithmStatus::Completed);
        CHECK(result.samples_used == 0);
        CHECK(result.function_evaluations == 0);
    }

    TEST_CASE("Stratified MC rejects overflowing scaled payloads and uncertainties", "[montecarlo][stratified][kko3]")
    {
        int calls = 0;
        ScalarFunctionFromStdFunc<1> constant([&](const VectorN<Real, 1>&) {
            ++calls;
            return std::numeric_limits<Real>::max() * Real(0.75);
        });
        CHECK_THROWS_AS(StratifiedMonteCarloIntegrator<1>().integrate(
            constant, VectorN<Real, 1>{0}, VectorN<Real, 1>{2}, 1, 2, 42), NumericalMethodError);
        CHECK(calls == 2);
        calls = 0;
        ScalarFunctionFromStdFunc<1> alternating([&](const VectorN<Real, 1>&) {
            const Real amplitude = std::sqrt(std::numeric_limits<Real>::max()) / Real(4);
            return ++calls % 2 ? -amplitude : amplitude;
        });
        CHECK_THROWS_AS(StratifiedMonteCarloIntegrator<1>().integrate(
            alternating, VectorN<Real, 1>{0}, VectorN<Real, 1>{std::numeric_limits<Real>::max() / Real(4)}, 1, 2, 42),
            NumericalMethodError);
        CHECK(calls == 2);
    }

}  // namespace

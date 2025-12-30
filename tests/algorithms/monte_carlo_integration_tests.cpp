#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#include "algorithms/MonteCarloIntegration.h"
#include "base/Function.h"
#endif

using namespace MML;
using namespace MML::Testing;

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace MML::Tests::Algorithms::MonteCarloIntegrationTests
{
    ///////////////////////////    TEST FUNCTIONS    ///////////////////////////

    // f(x) = 1 (constant) - integral over [0,1]^N = REAL(1.0)
    template<int N>
    class ConstantOne : public IScalarFunction<N>
    {
    public:
        Real operator()(const VectorN<Real, N>& x) const override { return REAL(1.0); }
    };

    // f(x) = x[0] (linear) - integral over [0,1]^N = REAL(0.5)
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

    // f(x,y) = x*y - integral over [0,1]^2 = REAL(0.25)
    class ProductXY : public IScalarFunction<2>
    {
    public:
        Real operator()(const VectorN<Real, 2>& x) const override { return x[0] * x[1]; }
    };

    // f(x,y,z) = x*y*z - integral over [0,1]^3 = REAL(0.125)
    class ProductXYZ : public IScalarFunction<3>
    {
    public:
        Real operator()(const VectorN<Real, 3>& x) const override { return x[0] * x[1] * x[2]; }
    };

    // f(x) = sin(x) - integral over [0,π] = REAL(2.0)
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
            Real sum = REAL(0.0);
            for (int i = 0; i < N; ++i) sum += x[i];
            return sum;
        }
    };

    ///////////////////////////    BASIC FUNCTIONALITY TESTS    ///////////////////////////

    TEST_CASE("MonteCarlo_1D_Constant", "[montecarlo][integration][1d]")
    {
        ConstantOne<1> func;
        MonteCarloIntegrator<1> mc(42);  // Fixed seed for reproducibility
        
        VectorN<Real, 1> lower{REAL(0.0)}, upper{REAL(1.0)};
        
        auto result = mc.integrate(func, lower, upper, MonteCarloConfig().samples(50000));
        
        REQUIRE_THAT(result.value, WithinAbs(REAL(1.0), REAL(0.01)));
        REQUIRE(result.samples_used == 50000);
    }

    TEST_CASE("MonteCarlo_1D_Linear", "[montecarlo][integration][1d]")
    {
        LinearX0<1> func;
        MonteCarloIntegrator<1> mc(42);
        
        VectorN<Real, 1> lower{REAL(0.0)}, upper{REAL(1.0)};
        
        // Integral of x from 0 to 1 = REAL(0.5)
        auto result = mc.integrate(func, lower, upper, MonteCarloConfig().samples(100000));
        
        REQUIRE_THAT(result.value, WithinAbs(REAL(0.5), REAL(0.01)));
    }

    TEST_CASE("MonteCarlo_1D_XSquared", "[montecarlo][integration][1d]")
    {
        XSquared1D func;
        MonteCarloIntegrator<1> mc(42);
        
        VectorN<Real, 1> lower{REAL(0.0)}, upper{REAL(1.0)};
        
        // Integral of x^2 from 0 to 1 = 1/3
        auto result = mc.integrate(func, lower, upper, MonteCarloConfig().samples(100000));
        
        REQUIRE_THAT(result.value, WithinAbs(REAL(1.0)/REAL(3.0), REAL(0.01)));
    }

    TEST_CASE("MonteCarlo_1D_Sin_0_to_Pi", "[montecarlo][integration][1d]")
    {
        SinFunc1D func;
        MonteCarloIntegrator<1> mc(42);
        
        VectorN<Real, 1> lower{REAL(0.0)}, upper{REAL(MML::Constants::PI)};
        
        // Integral of sin(x) from 0 to π = REAL(2.0)
        auto result = mc.integrate(func, lower, upper, MonteCarloConfig().samples(100000));
        
        REQUIRE_THAT(result.value, WithinAbs(REAL(2.0), REAL(0.02)));
    }

    ///////////////////////////    2D INTEGRATION TESTS    ///////////////////////////

    TEST_CASE("MonteCarlo_2D_Constant", "[montecarlo][integration][2d]")
    {
        ConstantOne<2> func;
        MonteCarloIntegrator<2> mc(42);
        
        VectorN<Real, 2> lower{REAL(0.0), REAL(0.0)}, upper{REAL(1.0), REAL(1.0)};
        
        auto result = mc.integrate(func, lower, upper, MonteCarloConfig().samples(50000));
        
        REQUIRE_THAT(result.value, WithinAbs(REAL(1.0), REAL(0.01)));
    }

    TEST_CASE("MonteCarlo_2D_ProductXY", "[montecarlo][integration][2d]")
    {
        ProductXY func;
        MonteCarloIntegrator<2> mc(42);
        
        VectorN<Real, 2> lower{REAL(0.0), REAL(0.0)}, upper{REAL(1.0), REAL(1.0)};
        
        // Integral of x*y over [0,1]^2 = REAL(0.25)
        auto result = mc.integrate(func, lower, upper, MonteCarloConfig().samples(100000));
        
        REQUIRE_THAT(result.value, WithinAbs(REAL(0.25), REAL(0.01)));
    }

    TEST_CASE("MonteCarlo_2D_RectangularDomain", "[montecarlo][integration][2d]")
    {
        ConstantOne<2> func;
        MonteCarloIntegrator<2> mc(42);
        
        VectorN<Real, 2> lower{REAL(0.0), REAL(0.0)}, upper{REAL(2.0), REAL(3.0)};
        
        // Volume = 2 * 3 = 6
        auto result = mc.integrate(func, lower, upper, MonteCarloConfig().samples(50000));
        
        REQUIRE_THAT(result.value, WithinAbs(REAL(6.0), REAL(0.05)));
    }

    ///////////////////////////    3D INTEGRATION TESTS    ///////////////////////////

    TEST_CASE("MonteCarlo_3D_Constant_UnitCube", "[montecarlo][integration][3d]")
    {
        ConstantOne<3> func;
        MonteCarloIntegrator<3> mc(42);
        
        VectorN<Real, 3> lower{REAL(0.0), REAL(0.0), REAL(0.0)}, upper{REAL(1.0), REAL(1.0), REAL(1.0)};
        
        auto result = mc.integrate(func, lower, upper, MonteCarloConfig().samples(50000));
        
        REQUIRE_THAT(result.value, WithinAbs(REAL(1.0), REAL(0.01)));
    }

    TEST_CASE("MonteCarlo_3D_ProductXYZ", "[montecarlo][integration][3d]")
    {
        ProductXYZ func;
        MonteCarloIntegrator<3> mc(42);
        
        VectorN<Real, 3> lower{REAL(0.0), REAL(0.0), REAL(0.0)}, upper{REAL(1.0), REAL(1.0), REAL(1.0)};
        
        // Integral of x*y*z over [0,1]^3 = REAL(0.125)
        auto result = mc.integrate(func, lower, upper, MonteCarloConfig().samples(100000));
        
        REQUIRE_THAT(result.value, WithinAbs(REAL(0.125), REAL(0.005)));
    }

    ///////////////////////////    HIGHER DIMENSIONS    ///////////////////////////

    TEST_CASE("MonteCarlo_4D_SumOfCoords", "[montecarlo][integration][4d]")
    {
        SumOfCoords<4> func;
        MonteCarloIntegrator<4> mc(42);
        
        VectorN<Real, 4> lower{REAL(0.0), REAL(0.0), REAL(0.0), REAL(0.0)}, upper{REAL(1.0), REAL(1.0), REAL(1.0), REAL(1.0)};
        
        // Integral of (x1+x2+x3+x4) over [0,1]^4 = 4 * REAL(0.5) = REAL(2.0)
        auto result = mc.integrate(func, lower, upper, MonteCarloConfig().samples(100000));
        
        REQUIRE_THAT(result.value, WithinAbs(REAL(2.0), REAL(0.02)));
    }

    TEST_CASE("MonteCarlo_5D_Constant", "[montecarlo][integration][5d]")
    {
        ConstantOne<5> func;
        MonteCarloIntegrator<5> mc(42);
        
        VectorN<Real, 5> lower, upper;
        for (int i = 0; i < 5; ++i) {
            lower[i] = REAL(0.0);
            upper[i] = REAL(1.0);
        }
        
        auto result = mc.integrate(func, lower, upper, MonteCarloConfig().samples(50000));
        
        REQUIRE_THAT(result.value, WithinAbs(REAL(1.0), REAL(0.02)));
    }

    ///////////////////////////    PI ESTIMATION (CLASSIC EXAMPLE)    ///////////////////////////

    TEST_CASE("MonteCarlo_EstimatePi", "[montecarlo][integration][pi][classic]")
    {
        // Classic textbook example: estimate π using hit-or-miss
        auto result = EstimatePi(500000, 42);
        
        // π ≈ REAL(3.14159)...
        REQUIRE_THAT(result.value, WithinAbs(REAL(MML::Constants::PI), REAL(0.01)));
        
        // Error estimate should be reasonable
        REQUIRE(result.error_estimate > REAL(0.0));
        REQUIRE(result.error_estimate < REAL(0.02));
    }

    ///////////////////////////    UNIT BALL VOLUME    ///////////////////////////

    TEST_CASE("MonteCarlo_UnitBallVolume_2D", "[montecarlo][integration][ball]")
    {
        // 2D "ball" = circle, volume = π
        auto result = EstimateUnitBallVolume<2>(500000, 42);
        
        REQUIRE_THAT(result.value, WithinAbs(REAL(MML::Constants::PI), REAL(0.02)));
    }

    TEST_CASE("MonteCarlo_UnitBallVolume_3D", "[montecarlo][integration][ball]")
    {
        // 3D ball volume = (4/3)π ≈ REAL(4.189)
        auto result = EstimateUnitBallVolume<3>(500000, 42);
        
        Real expected = (REAL(4.0)/REAL(3.0)) * REAL(MML::Constants::PI);
        REQUIRE_THAT(result.value, WithinAbs(expected, REAL(0.05)));
    }

    TEST_CASE("MonteCarlo_UnitBallVolume_4D", "[montecarlo][integration][ball]")
    {
        // 4D ball volume = (π²/2) ≈ REAL(4.935)
        auto result = EstimateUnitBallVolume<4>(500000, 42);
        
        Real expected = (REAL(MML::Constants::PI) * REAL(MML::Constants::PI)) / REAL(2.0);
        REQUIRE_THAT(result.value, WithinAbs(expected, REAL(0.08)));
    }

    ///////////////////////////    HIT-OR-MISS INTEGRATOR    ///////////////////////////

    TEST_CASE("HitOrMiss_CircleArea", "[montecarlo][hitorrmiss]")
    {
        HitOrMissIntegrator<2> integrator(42);
        
        VectorN<Real, 2> lower{REAL(0.0), REAL(0.0)}, upper{REAL(1.0), REAL(1.0)};
        
        // Quarter circle in first quadrant: x² + y² ≤ 1
        auto inside_quarter_circle = [](const VectorN<Real, 2>& p) {
            return p[0]*p[0] + p[1]*p[1] <= REAL(1.0);
        };
        
        auto result = integrator.estimateVolume(inside_quarter_circle, lower, upper, 500000);
        
        // Area of quarter circle = π/4 ≈ REAL(0.785)
        REQUIRE_THAT(result.value, WithinAbs(REAL(MML::Constants::PI)/REAL(4.0), REAL(0.01)));
    }

    TEST_CASE("HitOrMiss_TriangleArea", "[montecarlo][hitorrmiss]")
    {
        HitOrMissIntegrator<2> integrator(42);
        
        VectorN<Real, 2> lower{REAL(0.0), REAL(0.0)}, upper{REAL(1.0), REAL(1.0)};
        
        // Triangle: x + y ≤ 1
        auto inside_triangle = [](const VectorN<Real, 2>& p) {
            return p[0] + p[1] <= REAL(1.0);
        };
        
        auto result = integrator.estimateVolume(inside_triangle, lower, upper, 100000);
        
        // Area of triangle = REAL(0.5)
        REQUIRE_THAT(result.value, WithinAbs(REAL(0.5), REAL(0.01)));
    }

    ///////////////////////////    ANTITHETIC VARIATES    ///////////////////////////

    TEST_CASE("MonteCarlo_AntitheticVariates", "[montecarlo][variance_reduction]")
    {
        LinearX0<1> func;
        
        VectorN<Real, 1> lower{REAL(0.0)}, upper{REAL(1.0)};
        
        // Without antithetic
        MonteCarloIntegrator<1> mc1(42);
        auto result1 = mc1.integrate(func, lower, upper, 
                                     MonteCarloConfig().samples(10000).antithetic(false));
        
        // With antithetic variates (should have lower variance for monotonic functions)
        MonteCarloIntegrator<1> mc2(42);
        auto result2 = mc2.integrate(func, lower, upper, 
                                     MonteCarloConfig().samples(10000).antithetic(true));
        
        // Both should give correct answer
        REQUIRE_THAT(result1.value, WithinAbs(REAL(0.5), REAL(0.02)));
        REQUIRE_THAT(result2.value, WithinAbs(REAL(0.5), REAL(0.02)));
        
        // Note: antithetic typically has lower variance for monotonic functions
        // but this is a statistical property, not guaranteed for every run
    }

    ///////////////////////////    STRATIFIED SAMPLING    ///////////////////////////

    TEST_CASE("Stratified_2D_ProductXY", "[montecarlo][stratified]")
    {
        ProductXY func;
        StratifiedMonteCarloIntegrator<2> mc(42);
        
        VectorN<Real, 2> lower{REAL(0.0), REAL(0.0)}, upper{REAL(1.0), REAL(1.0)};
        
        // 10 strata per dimension, 10 samples per stratum = 1000 total
        auto result = mc.integrate(func, lower, upper, 10, 10, 42);
        
        // Integral of x*y over [0,1]^2 = REAL(0.25)
        REQUIRE_THAT(result.value, WithinAbs(REAL(0.25), REAL(0.01)));
    }

    TEST_CASE("Stratified_3D_Constant", "[montecarlo][stratified]")
    {
        ConstantOne<3> func;
        StratifiedMonteCarloIntegrator<3> mc(42);
        
        VectorN<Real, 3> lower{REAL(0.0), REAL(0.0), REAL(0.0)}, upper{REAL(1.0), REAL(1.0), REAL(1.0)};
        
        // 5 strata per dimension, 5 samples per stratum = 625 total
        auto result = mc.integrate(func, lower, upper, 5, 5, 42);
        
        REQUIRE_THAT(result.value, WithinAbs(REAL(1.0), REAL(0.02)));
    }

    ///////////////////////////    CONVENIENCE FUNCTION TESTS    ///////////////////////////

    TEST_CASE("IntegrateMonteCarlo1D_Convenience", "[montecarlo][1d][convenience]")
    {
        // Test the convenience wrapper for 1D IRealFunction
        RealFunction func([](Real x) { return x * x; });
        
        auto result = IntegrateMonteCarlo1D(func, REAL(0.0), REAL(1.0), 
                                            MonteCarloConfig().samples(100000).randomSeed(42));
        
        // Integral of x^2 from 0 to 1 = 1/3
        REQUIRE_THAT(result.value, WithinAbs(REAL(1.0)/REAL(3.0), REAL(0.01)));
    }

    ///////////////////////////    ERROR ESTIMATION TESTS    ///////////////////////////

    TEST_CASE("MonteCarlo_ErrorEstimate_Reasonable", "[montecarlo][error]")
    {
        ConstantOne<2> func;
        MonteCarloIntegrator<2> mc(42);
        
        VectorN<Real, 2> lower{REAL(0.0), REAL(0.0)}, upper{REAL(1.0), REAL(1.0)};
        
        auto result = mc.integrate(func, lower, upper, MonteCarloConfig().samples(10000));
        
        // Error should be positive and small for constant function
        REQUIRE(result.error_estimate >= REAL(0.0));
        
        // The actual error should be within a few standard errors
        Real actual_error = std::abs(result.value - REAL(1.0));
        REQUIRE(actual_error < 3 * result.error_estimate + REAL(0.001));  // 3-sigma + small buffer
    }

    TEST_CASE("MonteCarlo_MoreSamples_LowerError", "[montecarlo][convergence]")
    {
        ProductXY func;
        
        VectorN<Real, 2> lower{REAL(0.0), REAL(0.0)}, upper{REAL(1.0), REAL(1.0)};
        
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
        // This captures most of the Gaussian (erf(3) ≈ REAL(0.9999779))
        VectorN<Real, 1> lower{-REAL(3.0)}, upper{REAL(3.0)};
        
        auto result = mc.integrate(func, lower, upper, MonteCarloConfig().samples(100000));
        
        // Full integral = √π ≈ REAL(1.7724538509), but we get slightly less
        Real expected = std::sqrt(REAL(MML::Constants::PI)) * std::erf(REAL(3.0));
        REQUIRE_THAT(result.value, WithinAbs(expected, REAL(0.02)));
    }

    TEST_CASE("MonteCarlo_2DGaussian", "[montecarlo][gaussian][2d]")
    {
        Gaussian2D func;
        MonteCarloIntegrator<2> mc(42);
        
        // Integrate exp(-(x²+y²)) over [-3,3]²
        VectorN<Real, 2> lower{-REAL(3.0), -REAL(3.0)}, upper{REAL(3.0), REAL(3.0)};
        
        auto result = mc.integrate(func, lower, upper, MonteCarloConfig().samples(200000));
        
        // Full 2D Gaussian integral = π, truncated integral ≈ π * erf(3)²
        Real expected = REAL(MML::Constants::PI) * std::erf(REAL(3.0)) * std::erf(REAL(3.0));
        REQUIRE_THAT(result.value, WithinAbs(expected, REAL(0.05)));
    }

}  // namespace

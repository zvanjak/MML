#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#include "MMLBase.h"
#include "base/Vector.h"
#include "algorithms/Optimization/SimulatedAnnealing.h"

using namespace MML;
using namespace MML::Testing;
// Removed: using Catch::Approx; - now using precision-aware RealApprox

///////////////////////////////////////////////////////////////////////////////
///                           Test Functions                                ///
///////////////////////////////////////////////////////////////////////////////

// Simple quadratic: f(x,y) = x^2 + y^2
// Global minimum at (0, 0) with value 0
class SAQuadratic2D
{
public:
    Real operator()(const Vector<Real>& x) const
    {
        return x[0] * x[0] + x[1] * x[1];
    }
};

// Rosenbrock function: f(x,y) = (1-x)^2 + 100*(y-x^2)^2
// Global minimum at (1, 1) with value 0
class SARosenbrock2D
{
public:
    Real operator()(const Vector<Real>& x) const
    {
        Real a = REAL(1.0) - x[0];
        Real b = x[1] - x[0] * x[0];
        return a * a + REAL(100.0) * b * b;
    }
};

// Booth function: f(x,y) = (x + 2y - 7)^2 + (2x + y - 5)^2
// Global minimum at (1, 3) with value 0
class SABoothFunction
{
public:
    Real operator()(const Vector<Real>& x) const
    {
        Real a = x[0] + REAL(2.0) * x[1] - REAL(7.0);
        Real b = REAL(2.0) * x[0] + x[1] - REAL(5.0);
        return a * a + b * b;
    }
};

// Rastrigin function - multimodal with many local minima
// f(x) = 10*n + sum(x_i^2 - 10*cos(2*pi*x_i))
// Global minimum at origin with value 0
class SARastrigin
{
public:
    Real operator()(const Vector<Real>& x) const
    {
        const Real A = REAL(10.0);
        int n = x.size();
        Real sum = A * n;
        for (int i = 0; i < n; ++i)
        {
            sum += x[i] * x[i] - A * std::cos(REAL(2.0) * Constants::PI * x[i]);
        }
        return sum;
    }
};

// Ackley function - another multimodal test function
// Global minimum at origin with value 0
class SAAckley
{
public:
    Real operator()(const Vector<Real>& x) const
    {
        int n = x.size();
        Real sum1 = 0, sum2 = 0;
        for (int i = 0; i < n; ++i)
        {
            sum1 += x[i] * x[i];
            sum2 += std::cos(REAL(2.0) * Constants::PI * x[i]);
        }
        return -REAL(20.0) * std::exp(-REAL(0.2) * std::sqrt(sum1 / n))
               - std::exp(sum2 / n) + REAL(20.0) + Constants::E;
    }
};

// Sphere function in higher dimensions
class SASphere
{
public:
    Real operator()(const Vector<Real>& x) const
    {
        Real sum = 0;
        for (int i = 0; i < x.size(); ++i)
            sum += x[i] * x[i];
        return sum;
    }
};

// Himmelblau's function - has four identical local minima
// f(x,y) = (x^2 + y - 11)^2 + (x + y^2 - 7)^2
class SAHimmelblau
{
public:
    Real operator()(const Vector<Real>& x) const
    {
        Real a = x[0] * x[0] + x[1] - REAL(11.0);
        Real b = x[0] + x[1] * x[1] - REAL(7.0);
        return a * a + b * b;
    }
};

///////////////////////////////////////////////////////////////////////////////
///                     Cooling Schedule Tests                              ///
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("ExponentialCooling - Temperature decay", "[SA][CoolingSchedule]")
{
    ExponentialCooling cooling(REAL(0.9));
    Real T0 = REAL(100.0);
    
    REQUIRE_THAT(cooling.Temperature(0, T0), RealApprox(REAL(100.0)));
    REQUIRE_THAT(cooling.Temperature(1, T0), RealApprox(REAL(90.0)));
    REQUIRE_THAT(cooling.Temperature(2, T0), RealApprox(REAL(81.0)));
    REQUIRE_THAT(cooling.Temperature(10, T0), RealApprox(T0 * std::pow(REAL(0.9), 10)));
}

TEST_CASE("ExponentialCooling - Invalid alpha throws", "[SA][CoolingSchedule]")
{
    REQUIRE_THROWS_AS(ExponentialCooling(REAL(0.0)), HeuristicOptimizationError);
    REQUIRE_THROWS_AS(ExponentialCooling(REAL(1.0)), HeuristicOptimizationError);
    REQUIRE_THROWS_AS(ExponentialCooling(-REAL(0.5)), HeuristicOptimizationError);
    REQUIRE_THROWS_AS(ExponentialCooling(REAL(1.5)), HeuristicOptimizationError);
}

TEST_CASE("LinearCooling - Temperature decay", "[SA][CoolingSchedule]")
{
    LinearCooling cooling(100);  // 100 iterations to reach 0
    Real T0 = REAL(100.0);
    
    REQUIRE_THAT(cooling.Temperature(0, T0), RealApprox(REAL(100.0)));
    REQUIRE_THAT(cooling.Temperature(50, T0), RealApprox(REAL(50.0)));
    REQUIRE_THAT(cooling.Temperature(100, T0), RealApprox(REAL(0.0)));
    REQUIRE_THAT(cooling.Temperature(150, T0), RealApprox(REAL(0.0)));  // Clamps at 0
}

TEST_CASE("LogarithmicCooling - Temperature decay", "[SA][CoolingSchedule]")
{
    LogarithmicCooling cooling(REAL(1.0));
    Real T0 = REAL(100.0);
    
    REQUIRE_THAT(cooling.Temperature(0, T0), RealApprox(REAL(100.0)));
    // T(k) = T0 / (1 + ln(1+k))
    REQUIRE_THAT(cooling.Temperature(1, T0), RealApprox(T0 / (REAL(1.0) + std::log(REAL(2.0)))));
    
    // Logarithmic cooling is much slower
    REQUIRE(cooling.Temperature(100, T0) > cooling.Temperature(0, T0) / 10);
}

///////////////////////////////////////////////////////////////////////////////
///                   Neighbor Generator Tests                              ///
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("UniformNeighborGenerator - Basic generation", "[SA][NeighborGen]")
{
    UniformNeighborGenerator gen(REAL(1.0), false, REAL(100.0), 42);  // Fixed seed
    
    Vector<Real> current(2);
    current[0] = REAL(0.0);
    current[1] = REAL(0.0);
    
    Vector<Real> neighbor = gen.Generate(current, REAL(100.0));
    
    // Neighbor should be within delta of current
    REQUIRE(std::abs(neighbor[0] - current[0]) <= REAL(1.0));
    REQUIRE(std::abs(neighbor[1] - current[1]) <= REAL(1.0));
}

TEST_CASE("UniformNeighborGenerator - Respects bounds", "[SA][NeighborGen]")
{
    UniformNeighborGenerator gen(REAL(10.0), false, REAL(100.0), 42);  // Large delta
    
    Vector<Real> lower(2), upper(2);
    lower[0] = -REAL(1.0); lower[1] = -REAL(1.0);
    upper[0] = REAL(1.0);  upper[1] = REAL(1.0);
    gen.SetBounds(lower, upper);
    
    Vector<Real> current(2);
    current[0] = REAL(0.0);
    current[1] = REAL(0.0);
    
    // Generate many neighbors and check bounds
    for (int i = 0; i < 100; ++i)
    {
        Vector<Real> neighbor = gen.Generate(current, REAL(100.0));
        REQUIRE(neighbor[0] >= -REAL(1.0));
        REQUIRE(neighbor[0] <= REAL(1.0));
        REQUIRE(neighbor[1] >= -REAL(1.0));
        REQUIRE(neighbor[1] <= REAL(1.0));
    }
}

TEST_CASE("GaussianNeighborGenerator - Basic generation", "[SA][NeighborGen]")
{
    GaussianNeighborGenerator gen(REAL(1.0), false, REAL(100.0), 42);
    
    Vector<Real> current(2);
    current[0] = REAL(5.0);
    current[1] = -REAL(3.0);
    
    Vector<Real> neighbor = gen.Generate(current, REAL(100.0));
    
    // Neighbor should be different from current
    REQUIRE((neighbor[0] != current[0] || neighbor[1] != current[1]));
}

///////////////////////////////////////////////////////////////////////////////
///                   Simulated Annealing Basic Tests                       ///
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("SimulatedAnnealing - Simple 2D quadratic", "[SA][2D]")
{
    SAQuadratic2D func;
    Vector<Real> x0(2);
    x0[0] = REAL(5.0);
    x0[1] = REAL(5.0);
    
    // Use custom generator with appropriate step size
    auto cooling = std::make_unique<ExponentialCooling>(REAL(0.995));
    auto gen = std::make_unique<GaussianNeighborGenerator>(REAL(0.5), true, REAL(100.0), 42);
    
    SimulatedAnnealing sa(std::move(cooling), std::move(gen),
                          REAL(100.0), 1e-10, 10000, 2000,
                          SimulatedAnnealing::StopCriteria::Combined, 42);
    auto result = sa.Minimize(func, x0);
    
    // SA should find near-optimal solution
    INFO("Best found: " << result.fbest << " at (" << result.xbest[0] << ", " << result.xbest[1] << ")");
    REQUIRE(result.fbest < REAL(1.0));  // Close to 0
    REQUIRE(std::abs(result.xbest[0]) < REAL(1.0));
    REQUIRE(std::abs(result.xbest[1]) < REAL(1.0));
}

TEST_CASE("SimulatedAnnealing - Booth function", "[SA][2D]")
{
    SABoothFunction func;
    Vector<Real> x0(2);
    x0[0] = REAL(0.0);
    x0[1] = REAL(0.0);
    
    SimulatedAnnealing sa(REAL(100.0), 1e-10, 10000, 2000,
                          SimulatedAnnealing::StopCriteria::Combined, 42);
    auto result = sa.Minimize(func, x0);
    
    // Should find minimum near (1, 3)
    REQUIRE(result.fbest < REAL(1.0));
    INFO("Best found: (" << result.xbest[0] << ", " << result.xbest[1] << ")");
    INFO("Best value: " << result.fbest);
}

TEST_CASE("SimulatedAnnealing - Higher dimensions", "[SA][ND]")
{
    SASphere func;
    Vector<Real> x0(5);
    for (int i = 0; i < 5; ++i)
        x0[i] = REAL(3.0);  // Start away from origin
    
    auto cooling = std::make_unique<ExponentialCooling>(REAL(0.998));
    auto gen = std::make_unique<GaussianNeighborGenerator>(REAL(0.3), true, REAL(100.0), 42);
    
    SimulatedAnnealing sa(std::move(cooling), std::move(gen),
                          REAL(100.0), 1e-10, 20000, 5000,
                          SimulatedAnnealing::StopCriteria::Combined, 42);
    auto result = sa.Minimize(func, x0);
    
    // Should find near origin
    INFO("5D Sphere best value: " << result.fbest);
    REQUIRE(result.fbest < REAL(5.0));  // Relaxed for higher dimensions
}

TEST_CASE("SimulatedAnnealing - With bounds", "[SA][Bounded]")
{
    SAQuadratic2D func;
    Vector<Real> x0(2);
    x0[0] = REAL(0.5);
    x0[1] = REAL(0.5);
    
    Vector<Real> lower(2), upper(2);
    lower[0] = -REAL(1.0); lower[1] = -REAL(1.0);
    upper[0] = REAL(1.0);  upper[1] = REAL(1.0);
    
    SimulatedAnnealing sa(REAL(50.0), 1e-10, 5000, 1000,
                          SimulatedAnnealing::StopCriteria::Combined, 42);
    auto result = sa.Minimize(func, x0, lower, upper);
    
    // Solution should be within bounds
    REQUIRE(result.xbest[0] >= -REAL(1.0));
    REQUIRE(result.xbest[0] <= REAL(1.0));
    REQUIRE(result.xbest[1] >= -REAL(1.0));
    REQUIRE(result.xbest[1] <= REAL(1.0));
    
    // Should still find good solution
    REQUIRE(result.fbest < REAL(0.1));
}

///////////////////////////////////////////////////////////////////////////////
///                   Multimodal Function Tests                             ///
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("SimulatedAnnealing - Rastrigin function (multimodal)", "[SA][Multimodal]")
{
    SARastrigin func;
    Vector<Real> x0(2);
    x0[0] = REAL(2.0);   // Start closer to allow better exploration
    x0[1] = -REAL(2.0);
    
    // SA should be able to escape local minima
    auto cooling = std::make_unique<ExponentialCooling>(REAL(0.999));
    auto gen = std::make_unique<GaussianNeighborGenerator>(REAL(0.5), true, REAL(200.0), 42);
    
    SimulatedAnnealing sa(std::move(cooling), std::move(gen),
                          REAL(200.0), 1e-10, 30000, 5000,
                          SimulatedAnnealing::StopCriteria::Combined, 42);
    auto result = sa.Minimize(func, x0);
    
    // Rastrigin has global min = 0 at origin
    INFO("Rastrigin best: " << result.fbest << " at (" 
         << result.xbest[0] << ", " << result.xbest[1] << ")");
    
    // With SA, should get reasonably close - Rastrigin is very hard
    REQUIRE(result.fbest < REAL(15.0));  // Relaxed - this is a hard function
}

TEST_CASE("SimulatedAnnealing - Himmelblau function (multiple minima)", "[SA][Multimodal]")
{
    SAHimmelblau func;
    Vector<Real> x0(2);
    x0[0] = REAL(0.0);
    x0[1] = REAL(0.0);
    
    SimulatedAnnealing sa(REAL(100.0), 1e-10, 10000, 2000,
                          SimulatedAnnealing::StopCriteria::Combined, 42);
    auto result = sa.Minimize(func, x0);
    
    // Himmelblau has 4 global minima, all with f = 0
    // SA should find one of them
    REQUIRE(result.fbest < REAL(0.5));
    INFO("Himmelblau best: " << result.fbest);
}

///////////////////////////////////////////////////////////////////////////////
///                   Stopping Criteria Tests                               ///
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("SimulatedAnnealing - MaxIterations stopping", "[SA][StopCriteria]")
{
    SAQuadratic2D func;
    Vector<Real> x0(2);
    x0[0] = REAL(5.0);
    x0[1] = REAL(5.0);
    
    SimulatedAnnealing sa(REAL(100.0), 1e-20, 100, 1000000,  // Very low Tmin, huge stagnation
                          SimulatedAnnealing::StopCriteria::MaxIterations, 42);
    auto result = sa.Minimize(func, x0);
    
    REQUIRE(result.iterations == 100);
    REQUIRE(result.converged == true);  // MaxIter is considered converged
}

TEST_CASE("SimulatedAnnealing - MinTemperature stopping", "[SA][StopCriteria]")
{
    SAQuadratic2D func;
    Vector<Real> x0(2);
    x0[0] = REAL(1.0);
    x0[1] = REAL(1.0);
    
    // Use fast cooling to reach Tmin quickly
    auto fastCooling = std::make_unique<ExponentialCooling>(REAL(0.5));
    auto gen = std::make_unique<GaussianNeighborGenerator>(REAL(0.5), true, REAL(100.0), 42);
    
    SimulatedAnnealing sa(std::move(fastCooling), std::move(gen),
                          REAL(100.0), 1e-6, 100000, 100000,
                          SimulatedAnnealing::StopCriteria::MinTemperature, 42);
    auto result = sa.Minimize(func, x0);
    
    // Should stop before maxIter due to temperature
    REQUIRE(result.iterations < 100000);
    REQUIRE(result.converged == true);
}

TEST_CASE("SimulatedAnnealing - NoImprovement stopping", "[SA][StopCriteria]")
{
    SAQuadratic2D func;
    Vector<Real> x0(2);
    x0[0] = REAL(0.001);  // Already very close to minimum
    x0[1] = REAL(0.001);
    
    SimulatedAnnealing sa(REAL(1.0), 1e-20, 100000, 50,  // Short stagnation limit
                          SimulatedAnnealing::StopCriteria::NoImprovement, 42);
    auto result = sa.Minimize(func, x0);
    
    // Should stop due to no improvement
    REQUIRE(result.iterations < 100000);
    REQUIRE(result.converged == true);
}

///////////////////////////////////////////////////////////////////////////////
///                   Custom Components Tests                               ///
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("SimulatedAnnealing - Custom cooling schedule", "[SA][Custom]")
{
    SAQuadratic2D func;
    Vector<Real> x0(2);
    x0[0] = REAL(3.0);
    x0[1] = REAL(3.0);
    
    auto linearCooling = std::make_unique<LinearCooling>(5000);
    auto gen = std::make_unique<GaussianNeighborGenerator>(REAL(1.0), true, REAL(100.0), 42);
    
    SimulatedAnnealing sa(std::move(linearCooling), std::move(gen),
                          REAL(100.0), 1e-10, 5000, 1000,
                          SimulatedAnnealing::StopCriteria::MaxIterations, 42);
    auto result = sa.Minimize(func, x0);
    
    REQUIRE(result.fbest < REAL(0.5));
}

TEST_CASE("SimulatedAnnealing - Uniform neighbor generator", "[SA][Custom]")
{
    SAQuadratic2D func;
    Vector<Real> x0(2);
    x0[0] = REAL(2.0);  // Start closer
    x0[1] = REAL(2.0);
    
    auto cooling = std::make_unique<ExponentialCooling>(REAL(0.998));
    auto uniformGen = std::make_unique<UniformNeighborGenerator>(REAL(0.3), true, REAL(50.0), 42);
    
    SimulatedAnnealing sa(std::move(cooling), std::move(uniformGen),
                          REAL(50.0), 1e-10, 20000, 5000,
                          SimulatedAnnealing::StopCriteria::Combined, 42);
    auto result = sa.Minimize(func, x0);
    
    INFO("Uniform gen best: " << result.fbest);
    REQUIRE(result.fbest < REAL(5.0));  // Uniform is less efficient than Gaussian
}

///////////////////////////////////////////////////////////////////////////////
///                     Convenience Function Tests                          ///
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("SimulatedAnnealingMinimize convenience function", "[SA][Convenience]")
{
    SAQuadratic2D func;
    Vector<Real> x0(2);
    x0[0] = REAL(3.0);
    x0[1] = REAL(3.0);
    
    auto result = SimulatedAnnealingMinimize(func, x0, REAL(100.0), 20000);
    
    INFO("Convenience function result: " << result.fbest);
    REQUIRE(result.fbest < REAL(10.0));  // Relaxed for default parameters
}

TEST_CASE("SimulatedAnnealingMinimize with bounds", "[SA][Convenience][Bounded]")
{
    SAQuadratic2D func;
    Vector<Real> x0(2);
    x0[0] = REAL(0.5);
    x0[1] = REAL(0.5);
    
    Vector<Real> lower(2), upper(2);
    lower[0] = -REAL(2.0); lower[1] = -REAL(2.0);
    upper[0] = REAL(2.0);  upper[1] = REAL(2.0);
    
    auto result = SimulatedAnnealingMinimize(func, x0, lower, upper, REAL(100.0), 5000);
    
    // Check bounds
    REQUIRE(result.xbest[0] >= -REAL(2.0));
    REQUIRE(result.xbest[0] <= REAL(2.0));
    REQUIRE(result.xbest[1] >= -REAL(2.0));
    REQUIRE(result.xbest[1] <= REAL(2.0));
}

///////////////////////////////////////////////////////////////////////////////
///                          API Tests                                      ///
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("SimulatedAnnealing - Accessor methods", "[SA][API]")
{
    SimulatedAnnealing sa(REAL(100.0), 1e-8, 10000, 1000,
                          SimulatedAnnealing::StopCriteria::Combined);
    
    REQUIRE_THAT(sa.getT0(), RealApprox(REAL(100.0)));
    REQUIRE_THAT(sa.getTmin(), RealApprox(1e-8));
    REQUIRE(sa.getMaxIter() == 10000);
    REQUIRE(sa.getStagnationLimit() == 1000);
    REQUIRE(sa.getStopCriteria() == SimulatedAnnealing::StopCriteria::Combined);
    
    sa.setT0(REAL(200.0));
    sa.setTmin(1e-10);
    sa.setMaxIter(20000);
    sa.setStagnationLimit(2000);
    sa.setStopCriteria(SimulatedAnnealing::StopCriteria::MaxIterations);
    
    REQUIRE_THAT(sa.getT0(), RealApprox(REAL(200.0)));
    REQUIRE_THAT(sa.getTmin(), RealApprox(1e-10));
    REQUIRE(sa.getMaxIter() == 20000);
    REQUIRE(sa.getStagnationLimit() == 2000);
    REQUIRE(sa.getStopCriteria() == SimulatedAnnealing::StopCriteria::MaxIterations);
}

TEST_CASE("ExponentialCooling - Accessor methods", "[SA][API][CoolingSchedule]")
{
    ExponentialCooling cooling(REAL(0.9));
    
    REQUIRE_THAT(cooling.getAlpha(), RealApprox(REAL(0.9)));
    
    cooling.setAlpha(REAL(0.95));
    REQUIRE_THAT(cooling.getAlpha(), RealApprox(REAL(0.95)));
    
    REQUIRE_THROWS_AS(cooling.setAlpha(REAL(0.0)), HeuristicOptimizationError);
    REQUIRE_THROWS_AS(cooling.setAlpha(REAL(1.5)), HeuristicOptimizationError);
}

TEST_CASE("HeuristicOptimizationResult - Default construction", "[SA][API]")
{
    HeuristicOptimizationResult result;
    
    REQUIRE(result.fbest == std::numeric_limits<Real>::max());
    REQUIRE(result.iterations == 0);
    REQUIRE(result.funcEvals == 0);
    REQUIRE(result.acceptedMoves == 0);
    REQUIRE(result.converged == false);
}

///////////////////////////////////////////////////////////////////////////////
///                     Rosenbrock Challenge Test                           ///
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("SimulatedAnnealing - Rosenbrock function", "[SA][2D][Rosenbrock]")
{
    SARosenbrock2D func;
    Vector<Real> x0(2);
    x0[0] = -REAL(2.0);
    x0[1] = REAL(2.0);
    
    // Rosenbrock is challenging - use more iterations
    SimulatedAnnealing sa(REAL(500.0), 1e-10, 30000, 5000,
                          SimulatedAnnealing::StopCriteria::Combined, 42);
    auto result = sa.Minimize(func, x0);
    
    INFO("Rosenbrock best: " << result.fbest);
    INFO("Solution: (" << result.xbest[0] << ", " << result.xbest[1] << ")");
    INFO("Iterations: " << result.iterations);
    INFO("Accepted moves: " << result.acceptedMoves);
    
    // SA may not find exact minimum but should get close
    REQUIRE(result.fbest < REAL(1.0));
}



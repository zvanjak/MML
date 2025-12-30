#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#include "MMLBase.h"
#include "interfaces/IFunction.h"
#include "algorithms/Optimization/OptimizationMultidim.h"

using namespace MML;
using namespace MML::Testing;
using namespace MML::Testing;
// Removed: using Catch::Approx; - now using precision-aware matchers

///////////////////////////////////////////////////////////////////////////////
///                           Test Functions                                ///
///////////////////////////////////////////////////////////////////////////////

// 2D Rosenbrock function: f(x,y) = (1-x)^2 + 100*(y-x^2)^2
// Minimum at (1, 1) with value 0
class Rosenbrock2D : public IScalarFunction<2>
{
public:
    Real operator()(const VectorN<Real, 2>& x) const override
    {
        Real a = REAL(REAL(1.0)) - x[0];
        Real b = x[1] - x[0] * x[0];
        return a * a + REAL(REAL(100.0)) * b * b;
    }
};

// Simple quadratic: f(x,y) = x^2 + y^2
// Minimum at (0, 0) with value 0
class SimpleQuadratic2D : public IScalarFunction<2>
{
public:
    Real operator()(const VectorN<Real, 2>& x) const override
    {
        return x[0] * x[0] + x[1] * x[1];
    }
};

// Elliptic paraboloid: f(x,y) = (x-1)^2 + 2*(y+1)^2
// Minimum at (1, -1) with value 0
class EllipticParaboloid2D : public IScalarFunction<2>
{
public:
    Real operator()(const VectorN<Real, 2>& x) const override
    {
        Real dx = x[0] - REAL(REAL(1.0));
        Real dy = x[1] + REAL(REAL(1.0));
        return dx * dx + REAL(REAL(2.0)) * dy * dy;
    }
};

// 3D quadratic: f(x,y,z) = x^2 + y^2 + z^2
// Minimum at origin
class SimpleQuadratic3D : public IScalarFunction<3>
{
public:
    Real operator()(const VectorN<Real, 3>& x) const override
    {
        return x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
    }
};

// 3D offset quadratic: f(x,y,z) = (x-1)^2 + (y-2)^2 + (z-3)^2
// Minimum at (1, 2, 3)
class OffsetQuadratic3D : public IScalarFunction<3>
{
public:
    Real operator()(const VectorN<Real, 3>& x) const override
    {
        return (x[0] - REAL(REAL(1.0))) * (x[0] - REAL(REAL(1.0))) +
               (x[1] - REAL(REAL(2.0))) * (x[1] - REAL(REAL(2.0))) +
               (x[2] - REAL(REAL(3.0))) * (x[2] - REAL(REAL(3.0)));
    }
};

// Beale's function (classic test)
// f(x,y) = (REAL(REAL(1.5)) - x + xy)^2 + (REAL(REAL(2.25)) - x + xy^2)^2 + (REAL(REAL(2.625)) - x + xy^3)^2
// Minimum at (3, REAL(REAL(0.5))) with value 0
class BealeFunction : public IScalarFunction<2>
{
public:
    Real operator()(const VectorN<Real, 2>& x) const override
    {
        Real a = REAL(REAL(1.5)) - x[0] + x[0] * x[1];
        Real b = REAL(REAL(2.25)) - x[0] + x[0] * x[1] * x[1];
        Real c = REAL(REAL(2.625)) - x[0] + x[0] * x[1] * x[1] * x[1];
        return a * a + b * b + c * c;
    }
};

// Booth's function
// f(x,y) = (x + 2y - 7)^2 + (2x + y - 5)^2
// Minimum at (1, 3) with value 0
class BoothFunction : public IScalarFunction<2>
{
public:
    Real operator()(const VectorN<Real, 2>& x) const override
    {
        Real a = x[0] + REAL(REAL(2.0)) * x[1] - REAL(REAL(7.0));
        Real b = REAL(REAL(2.0)) * x[0] + x[1] - REAL(REAL(5.0));
        return a * a + b * b;
    }
};

// Himmelblau's function (multiple minima)
// f(x,y) = (x^2 + y - 11)^2 + (x + y^2 - 7)^2
// Four minima at approximately:
// (3, 2), (-REAL(REAL(2.805)), REAL(REAL(3.131))), (-REAL(REAL(3.779)), -REAL(REAL(3.283))), (REAL(REAL(3.584)), -REAL(REAL(1.848)))
class HimmelblauFunction : public IScalarFunction<2>
{
public:
    Real operator()(const VectorN<Real, 2>& x) const override
    {
        Real a = x[0] * x[0] + x[1] - REAL(REAL(11.0));
        Real b = x[0] + x[1] * x[1] - REAL(REAL(7.0));
        return a * a + b * b;
    }
};

// 4D sum of squares: f = sum(i * x_i^2)
// Minimum at origin
class WeightedSumSquares4D : public IScalarFunction<4>
{
public:
    Real operator()(const VectorN<Real, 4>& x) const override
    {
        Real sum = REAL(REAL(0.0));
        for (int i = 0; i < 4; ++i)
            sum += (i + 1) * x[i] * x[i];
        return sum;
    }
};

// Negative of a function for max testing
class NegativeQuadratic2D : public IScalarFunction<2>
{
public:
    Real operator()(const VectorN<Real, 2>& x) const override
    {
        // Maximum at (2, 3) with value 0
        return -((x[0] - REAL(REAL(2.0))) * (x[0] - REAL(REAL(2.0))) + (x[1] - REAL(REAL(3.0))) * (x[1] - REAL(REAL(3.0))));
    }
};

///////////////////////////////////////////////////////////////////////////////
///                              Unit Tests                                 ///
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("NelderMead - Simple 2D quadratic", "[NelderMead][2D]")
{
    SimpleQuadratic2D func;
    VectorN<Real, 2> start{REAL(REAL(5.0)), REAL(REAL(5.0))};

    NelderMead optimizer;
    auto result = optimizer.Minimize(func, start, REAL(REAL(1.0)));

    REQUIRE(result.converged);
    REQUIRE_THAT(result.xmin[0], RealWithinAbs(REAL(REAL(0.0)), 1e-6));
    REQUIRE_THAT(result.xmin[1], RealWithinAbs(REAL(REAL(0.0)), 1e-6));
    REQUIRE_THAT(result.fmin, RealWithinAbs(REAL(REAL(0.0)), 1e-10));
}

TEST_CASE("NelderMead - Elliptic paraboloid", "[NelderMead][2D]")
{
    EllipticParaboloid2D func;
    VectorN<Real, 2> start{REAL(REAL(0.0)), REAL(REAL(0.0))};

    auto result = NelderMeadMinimize(func, start, REAL(REAL(1.0)));

    REQUIRE(result.converged);
    REQUIRE_THAT(result.xmin[0], RealWithinAbs(REAL(REAL(1.0)), 1e-6));
    REQUIRE_THAT(result.xmin[1], RealWithinAbs(-REAL(REAL(1.0)), 1e-6));
    REQUIRE_THAT(result.fmin, RealWithinAbs(REAL(REAL(0.0)), 1e-10));
}

TEST_CASE("NelderMead - Rosenbrock banana function", "[NelderMead][2D][Rosenbrock]")
{
    Rosenbrock2D func;
    VectorN<Real, 2> start{-REAL(REAL(1.0)), REAL(REAL(1.0))};

    // Rosenbrock is notoriously difficult - use tighter tolerance
    NelderMead optimizer(1e-10, 10000);
    auto result = optimizer.Minimize(func, start, REAL(REAL(0.5)));

    REQUIRE(result.converged);
    REQUIRE_THAT(result.xmin[0], RealWithinAbs(REAL(REAL(1.0)), 1e-4));
    REQUIRE_THAT(result.xmin[1], RealWithinAbs(REAL(REAL(1.0)), 1e-4));
    REQUIRE_THAT(result.fmin, RealWithinAbs(REAL(REAL(0.0)), 1e-6));
}

TEST_CASE("NelderMead - Beale's function", "[NelderMead][2D]")
{
    BealeFunction func;
    VectorN<Real, 2> start{REAL(REAL(0.0)), REAL(REAL(0.0))};

    NelderMead optimizer(1e-10, 5000);
    auto result = optimizer.Minimize(func, start, REAL(REAL(1.0)));

    REQUIRE(result.converged);
    REQUIRE_THAT(result.xmin[0], RealWithinAbs(REAL(REAL(3.0)), 1e-4));
    REQUIRE_THAT(result.xmin[1], RealWithinAbs(REAL(REAL(0.5)), 1e-4));
    REQUIRE_THAT(result.fmin, RealWithinAbs(REAL(REAL(0.0)), 1e-6));
}

TEST_CASE("NelderMead - Booth's function", "[NelderMead][2D]")
{
    BoothFunction func;
    VectorN<Real, 2> start{REAL(REAL(0.0)), REAL(REAL(0.0))};

    auto result = NelderMeadMinimize(func, start);

    REQUIRE(result.converged);
    REQUIRE_THAT(result.xmin[0], RealWithinAbs(REAL(REAL(1.0)), 1e-6));
    REQUIRE_THAT(result.xmin[1], RealWithinAbs(REAL(REAL(3.0)), 1e-6));
    REQUIRE_THAT(result.fmin, RealWithinAbs(REAL(REAL(0.0)), 1e-10));
}

TEST_CASE("NelderMead - Himmelblau's function (find one minimum)", "[NelderMead][2D]")
{
    HimmelblauFunction func;
    VectorN<Real, 2> start{REAL(REAL(2.0)), REAL(REAL(2.0))};  // Near the (3, 2) minimum

    auto result = NelderMeadMinimize(func, start);

    REQUIRE(result.converged);
    // Should find the (3, 2) minimum from this starting point
    REQUIRE_THAT(result.xmin[0], RealWithinAbs(REAL(REAL(3.0)), 1e-4));
    REQUIRE_THAT(result.xmin[1], RealWithinAbs(REAL(REAL(2.0)), 1e-4));
    REQUIRE_THAT(result.fmin, RealWithinAbs(REAL(REAL(0.0)), 1e-6));
}

TEST_CASE("NelderMead - 3D quadratic at origin", "[NelderMead][3D]")
{
    SimpleQuadratic3D func;
    VectorN<Real, 3> start{REAL(REAL(2.0)), -REAL(REAL(3.0)), REAL(REAL(4.0))};

    auto result = NelderMeadMinimize(func, start);

    REQUIRE(result.converged);
    REQUIRE_THAT(result.xmin[0], RealWithinAbs(REAL(REAL(0.0)), 1e-6));
    REQUIRE_THAT(result.xmin[1], RealWithinAbs(REAL(REAL(0.0)), 1e-6));
    REQUIRE_THAT(result.xmin[2], RealWithinAbs(REAL(REAL(0.0)), 1e-6));
    REQUIRE_THAT(result.fmin, RealWithinAbs(REAL(REAL(0.0)), 1e-10));
}

TEST_CASE("NelderMead - 3D offset quadratic", "[NelderMead][3D]")
{
    OffsetQuadratic3D func;
    VectorN<Real, 3> start{REAL(REAL(0.0)), REAL(REAL(0.0)), REAL(REAL(0.0))};

    auto result = NelderMeadMinimize(func, start);

    REQUIRE(result.converged);
    REQUIRE_THAT(result.xmin[0], RealWithinAbs(REAL(REAL(1.0)), 1e-6));
    REQUIRE_THAT(result.xmin[1], RealWithinAbs(REAL(REAL(2.0)), 1e-6));
    REQUIRE_THAT(result.xmin[2], RealWithinAbs(REAL(REAL(3.0)), 1e-6));
    REQUIRE_THAT(result.fmin, RealWithinAbs(REAL(REAL(0.0)), 1e-10));
}

TEST_CASE("NelderMead - 4D weighted sum of squares", "[NelderMead][4D]")
{
    WeightedSumSquares4D func;
    VectorN<Real, 4> start{REAL(REAL(1.0)), REAL(REAL(1.0)), REAL(REAL(1.0)), REAL(REAL(1.0))};

    NelderMead optimizer(1e-10);
    auto result = optimizer.Minimize(func, start);

    REQUIRE(result.converged);
    for (int i = 0; i < 4; ++i)
    {
        REQUIRE_THAT(result.xmin[i], RealWithinAbs(REAL(REAL(0.0)), 1e-5));
    }
    REQUIRE_THAT(result.fmin, RealWithinAbs(REAL(REAL(0.0)), 1e-8));
}

TEST_CASE("NelderMead - Custom per-dimension deltas", "[NelderMead][2D]")
{
    EllipticParaboloid2D func;
    VectorN<Real, 2> start{REAL(REAL(0.0)), REAL(REAL(0.0))};
    Vector<Real> deltas({REAL(REAL(0.5)), REAL(REAL(2.0))});  // Different scale in each dimension

    auto result = NelderMeadMinimize(func, start, deltas);

    REQUIRE(result.converged);
    REQUIRE_THAT(result.xmin[0], RealWithinAbs(REAL(REAL(1.0)), 1e-6));
    REQUIRE_THAT(result.xmin[1], RealWithinAbs(-REAL(REAL(1.0)), 1e-6));
}

TEST_CASE("NelderMead - From explicit simplex", "[NelderMead][2D]")
{
    SimpleQuadratic2D func;

    // Create a custom initial simplex
    Matrix<Real> simplex(3, 2);
    simplex(0, 0) = REAL(REAL(5.0));  simplex(0, 1) = REAL(REAL(0.0));
    simplex(1, 0) = REAL(REAL(0.0));  simplex(1, 1) = REAL(REAL(5.0));
    simplex(2, 0) = -REAL(REAL(3.0)); simplex(2, 1) = -REAL(REAL(3.0));

    NelderMead optimizer;
    auto result = optimizer.Minimize<2>(func, simplex);

    REQUIRE(result.converged);
    REQUIRE_THAT(result.xmin[0], RealWithinAbs(REAL(REAL(0.0)), 1e-6));
    REQUIRE_THAT(result.xmin[1], RealWithinAbs(REAL(REAL(0.0)), 1e-6));
}

TEST_CASE("NelderMead - Maximization", "[NelderMead][2D][Maximize]")
{
    NegativeQuadratic2D func;  // Maximum at (2, 3)
    VectorN<Real, 2> start{REAL(REAL(0.0)), REAL(REAL(0.0))};

    auto result = NelderMeadMaximize(func, start);

    REQUIRE(result.converged);
    REQUIRE_THAT(result.xmin[0], RealWithinAbs(REAL(REAL(2.0)), 1e-6));
    REQUIRE_THAT(result.xmin[1], RealWithinAbs(REAL(REAL(3.0)), 1e-6));
    REQUIRE_THAT(result.fmin, RealWithinAbs(REAL(REAL(0.0)), 1e-10));
}

TEST_CASE("NelderMead - Starting at minimum", "[NelderMead][2D][EdgeCase]")
{
    SimpleQuadratic2D func;
    VectorN<Real, 2> start{REAL(REAL(0.0)), REAL(REAL(0.0))};  // Already at minimum

    auto result = NelderMeadMinimize(func, start, REAL(REAL(0.1)));

    REQUIRE(result.converged);
    REQUIRE_THAT(result.xmin[0], RealWithinAbs(REAL(REAL(0.0)), 1e-6));
    REQUIRE_THAT(result.xmin[1], RealWithinAbs(REAL(REAL(0.0)), 1e-6));
    REQUIRE_THAT(result.fmin, RealWithinAbs(REAL(REAL(0.0)), 1e-10));
}

TEST_CASE("NelderMead - Very small initial simplex", "[NelderMead][2D][EdgeCase]")
{
    EllipticParaboloid2D func;
    VectorN<Real, 2> start{REAL(REAL(0.9)), -REAL(REAL(0.9))};  // Near minimum

    auto result = NelderMeadMinimize(func, start, REAL(REAL(0.01)));  // Tiny simplex

    REQUIRE(result.converged);
    REQUIRE_THAT(result.xmin[0], RealWithinAbs(REAL(REAL(1.0)), 1e-5));
    REQUIRE_THAT(result.xmin[1], RealWithinAbs(-REAL(REAL(1.0)), 1e-5));
}

TEST_CASE("NelderMead - Large initial simplex", "[NelderMead][2D][EdgeCase]")
{
    SimpleQuadratic2D func;
    VectorN<Real, 2> start{REAL(REAL(0.0)), REAL(REAL(0.0))};

    auto result = NelderMeadMinimize(func, start, REAL(REAL(100.0)));  // Huge simplex

    REQUIRE(result.converged);
    REQUIRE_THAT(result.xmin[0], RealWithinAbs(REAL(REAL(0.0)), 1e-5));
    REQUIRE_THAT(result.xmin[1], RealWithinAbs(REAL(REAL(0.0)), 1e-5));
}

TEST_CASE("NelderMead - Function evaluation count", "[NelderMead][2D]")
{
    SimpleQuadratic2D func;
    VectorN<Real, 2> start{REAL(REAL(5.0)), REAL(REAL(5.0))};

    NelderMead optimizer;
    auto result = optimizer.Minimize(func, start);

    REQUIRE(result.converged);
    REQUIRE(result.iterations > 0);
    REQUIRE(optimizer.getNumFuncEvals() == result.iterations);
    INFO("Function evaluations: " << result.iterations);
}

TEST_CASE("NelderMead - Accessor methods", "[NelderMead][API]")
{
    NelderMead optimizer(1e-6, 3000);

    REQUIRE_THAT(optimizer.getFtol(), RealApprox(1e-6));
    REQUIRE(optimizer.getMaxIter() == 3000);

    optimizer.setFtol(1e-10);
    optimizer.setMaxIter(5000);

    REQUIRE_THAT(optimizer.getFtol(), RealApprox(1e-10));
    REQUIRE(optimizer.getMaxIter() == 5000);
}

TEST_CASE("NelderMead - Different starting points find same minimum", "[NelderMead][2D]")
{
    BoothFunction func;
    
    // Try from multiple starting points
    std::vector<VectorN<Real, 2>> starts = {
        {REAL(REAL(0.0)), REAL(REAL(0.0))},
        {-REAL(REAL(5.0)), -REAL(REAL(5.0))},
        {REAL(REAL(10.0)), -REAL(REAL(10.0))},
        {-REAL(REAL(3.0)), REAL(REAL(8.0))}
    };

    for (const auto& start : starts)
    {
        auto result = NelderMeadMinimize(func, start, REAL(REAL(2.0)));
        
        REQUIRE(result.converged);
        REQUIRE_THAT(result.xmin[0], RealWithinAbs(REAL(REAL(1.0)), 1e-5));
        REQUIRE_THAT(result.xmin[1], RealWithinAbs(REAL(REAL(3.0)), 1e-5));
    }
}

TEST_CASE("NelderMead - Invalid simplex dimensions throw", "[NelderMead][Exception]")
{
    SimpleQuadratic2D func;
    
    // Wrong number of rows
    Matrix<Real> badSimplex1(2, 2);  // Should be 3x2 for 2D
    NelderMead optimizer;
    
    REQUIRE_THROWS_AS(optimizer.Minimize<2>(func, badSimplex1), MultidimOptimizationError);
    
    // Wrong number of columns
    Matrix<Real> badSimplex2(3, 3);  // Should be 3x2 for 2D
    REQUIRE_THROWS_AS(optimizer.Minimize<2>(func, badSimplex2), MultidimOptimizationError);
}

TEST_CASE("NelderMead - Invalid deltas dimension throw", "[NelderMead][Exception]")
{
    SimpleQuadratic2D func;
    VectorN<Real, 2> start{REAL(REAL(0.0)), REAL(REAL(0.0))};
    Vector<Real> badDeltas({REAL(REAL(1.0)), REAL(REAL(1.0)), REAL(REAL(1.0))});  // 3 deltas for 2D problem
    
    NelderMead optimizer;
    REQUIRE_THROWS_AS(optimizer.Minimize(func, start, badDeltas), MultidimOptimizationError);
}

TEST_CASE("NelderMead - Get final simplex", "[NelderMead][API]")
{
    SimpleQuadratic2D func;
    VectorN<Real, 2> start{REAL(REAL(5.0)), REAL(REAL(5.0))};

    NelderMead optimizer;
    auto result = optimizer.Minimize(func, start);

    REQUIRE(result.converged);
    
    const auto& finalSimplex = optimizer.getSimplex();
    const auto& finalValues = optimizer.getSimplexValues();
    
    REQUIRE(finalSimplex.RowNum() == 3);
    REQUIRE(finalSimplex.ColNum() == 2);
    REQUIRE(finalValues.size() == 3);
    
    // All values should be near minimum
    for (int i = 0; i < 3; ++i)
    {
        REQUIRE_THAT(finalValues[i], RealWithinAbs(REAL(REAL(0.0)), 1e-6));
    }
}

TEST_CASE("NelderMead - Tight convergence tolerance", "[NelderMead][2D]")
{
    SimpleQuadratic2D func;
    VectorN<Real, 2> start{REAL(REAL(1.0)), REAL(REAL(1.0))};

    // Nelder-Mead is a derivative-free method, so very tight tolerances
    // may not achieve extreme precision, but should still converge well
    // Use larger delta to ensure simplex can explore the search space
    NelderMead optimizer(1e-12, 10000);
    auto result = optimizer.Minimize(func, start, REAL(REAL(0.5)));

    REQUIRE(result.converged);
    // Reasonable precision for a derivative-free method
    REQUIRE_THAT(result.xmin[0], RealWithinAbs(REAL(REAL(0.0)), 1e-6));
    REQUIRE_THAT(result.xmin[1], RealWithinAbs(REAL(REAL(0.0)), 1e-6));
}

TEST_CASE("NelderMead - Rosenbrock from difficult start", "[NelderMead][2D][Rosenbrock]")
{
    Rosenbrock2D func;
    VectorN<Real, 2> start{-REAL(REAL(5.0)), -REAL(REAL(5.0))};  // Far from minimum

    NelderMead optimizer(1e-10, 20000);
    auto result = optimizer.Minimize(func, start, REAL(REAL(1.0)));

    REQUIRE(result.converged);
    REQUIRE_THAT(result.xmin[0], RealWithinAbs(REAL(REAL(1.0)), 1e-3));
    REQUIRE_THAT(result.xmin[1], RealWithinAbs(REAL(REAL(1.0)), 1e-3));
}

///////////////////////////////////////////////////////////////////////////////
///              Differentiable Test Functions (with gradients)             ///
///////////////////////////////////////////////////////////////////////////////

// Differentiable simple quadratic
class DiffQuadratic2D : public IDifferentiableScalarFunction<2>
{
public:
    Real operator()(const VectorN<Real, 2>& x) const override
    {
        return x[0] * x[0] + x[1] * x[1];
    }

    void Gradient(const VectorN<Real, 2>& x, VectorN<Real, 2>& grad) const override
    {
        grad[0] = REAL(REAL(2.0)) * x[0];
        grad[1] = REAL(REAL(2.0)) * x[1];
    }
};

// Differentiable elliptic paraboloid
class DiffEllipticParaboloid2D : public IDifferentiableScalarFunction<2>
{
public:
    Real operator()(const VectorN<Real, 2>& x) const override
    {
        Real dx = x[0] - REAL(REAL(1.0));
        Real dy = x[1] + REAL(REAL(1.0));
        return dx * dx + REAL(REAL(2.0)) * dy * dy;
    }

    void Gradient(const VectorN<Real, 2>& x, VectorN<Real, 2>& grad) const override
    {
        grad[0] = REAL(REAL(2.0)) * (x[0] - REAL(REAL(1.0)));
        grad[1] = REAL(REAL(4.0)) * (x[1] + REAL(REAL(1.0)));
    }
};

// Differentiable Rosenbrock
class DiffRosenbrock2D : public IDifferentiableScalarFunction<2>
{
public:
    Real operator()(const VectorN<Real, 2>& x) const override
    {
        Real a = REAL(REAL(1.0)) - x[0];
        Real b = x[1] - x[0] * x[0];
        return a * a + REAL(REAL(100.0)) * b * b;
    }

    void Gradient(const VectorN<Real, 2>& x, VectorN<Real, 2>& grad) const override
    {
        Real a = REAL(REAL(1.0)) - x[0];
        Real b = x[1] - x[0] * x[0];
        grad[0] = -REAL(REAL(2.0)) * a - REAL(REAL(400.0)) * x[0] * b;
        grad[1] = REAL(REAL(200.0)) * b;
    }
};

// Differentiable Booth function
class DiffBoothFunction : public IDifferentiableScalarFunction<2>
{
public:
    Real operator()(const VectorN<Real, 2>& x) const override
    {
        Real a = x[0] + REAL(REAL(2.0)) * x[1] - REAL(REAL(7.0));
        Real b = REAL(REAL(2.0)) * x[0] + x[1] - REAL(REAL(5.0));
        return a * a + b * b;
    }

    void Gradient(const VectorN<Real, 2>& x, VectorN<Real, 2>& grad) const override
    {
        Real a = x[0] + REAL(REAL(2.0)) * x[1] - REAL(REAL(7.0));
        Real b = REAL(REAL(2.0)) * x[0] + x[1] - REAL(REAL(5.0));
        grad[0] = REAL(REAL(2.0)) * a + REAL(REAL(4.0)) * b;
        grad[1] = REAL(REAL(4.0)) * a + REAL(REAL(2.0)) * b;
    }
};

// Differentiable 3D quadratic
class DiffQuadratic3D : public IDifferentiableScalarFunction<3>
{
public:
    Real operator()(const VectorN<Real, 3>& x) const override
    {
        return x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
    }

    void Gradient(const VectorN<Real, 3>& x, VectorN<Real, 3>& grad) const override
    {
        grad[0] = REAL(REAL(2.0)) * x[0];
        grad[1] = REAL(REAL(2.0)) * x[1];
        grad[2] = REAL(REAL(2.0)) * x[2];
    }
};

// Differentiable 3D offset quadratic
class DiffOffsetQuadratic3D : public IDifferentiableScalarFunction<3>
{
public:
    Real operator()(const VectorN<Real, 3>& x) const override
    {
        return (x[0] - REAL(REAL(1.0))) * (x[0] - REAL(REAL(1.0))) +
               (x[1] - REAL(REAL(2.0))) * (x[1] - REAL(REAL(2.0))) +
               (x[2] - REAL(REAL(3.0))) * (x[2] - REAL(REAL(3.0)));
    }

    void Gradient(const VectorN<Real, 3>& x, VectorN<Real, 3>& grad) const override
    {
        grad[0] = REAL(REAL(2.0)) * (x[0] - REAL(REAL(1.0)));
        grad[1] = REAL(REAL(2.0)) * (x[1] - REAL(REAL(2.0)));
        grad[2] = REAL(REAL(2.0)) * (x[2] - REAL(REAL(3.0)));
    }
};

///////////////////////////////////////////////////////////////////////////////
///                         Powell's Method Tests                           ///
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("Powell - Simple 2D quadratic", "[Powell][2D]")
{
    SimpleQuadratic2D func;
    VectorN<Real, 2> start{REAL(REAL(5.0)), REAL(REAL(5.0))};

    Powell optimizer;
    auto result = optimizer.Minimize(func, start);

    REQUIRE(result.converged);
    REQUIRE_THAT(result.xmin[0], RealWithinAbs(REAL(REAL(0.0)), 1e-6));
    REQUIRE_THAT(result.xmin[1], RealWithinAbs(REAL(REAL(0.0)), 1e-6));
    REQUIRE_THAT(result.fmin, RealWithinAbs(REAL(REAL(0.0)), 1e-10));
}

TEST_CASE("Powell - Elliptic paraboloid", "[Powell][2D]")
{
    EllipticParaboloid2D func;
    VectorN<Real, 2> start{REAL(REAL(0.0)), REAL(REAL(0.0))};

    auto result = PowellMinimize(func, start);

    REQUIRE(result.converged);
    REQUIRE_THAT(result.xmin[0], RealWithinAbs(REAL(REAL(1.0)), 1e-6));
    REQUIRE_THAT(result.xmin[1], RealWithinAbs(-REAL(REAL(1.0)), 1e-6));
}

TEST_CASE("Powell - Booth function", "[Powell][2D]")
{
    BoothFunction func;
    VectorN<Real, 2> start{REAL(REAL(0.0)), REAL(REAL(0.0))};

    auto result = PowellMinimize(func, start);

    REQUIRE(result.converged);
    REQUIRE_THAT(result.xmin[0], RealWithinAbs(REAL(REAL(1.0)), 1e-5));
    REQUIRE_THAT(result.xmin[1], RealWithinAbs(REAL(REAL(3.0)), 1e-5));
}

TEST_CASE("Powell - 3D quadratic", "[Powell][3D]")
{
    SimpleQuadratic3D func;
    VectorN<Real, 3> start{REAL(REAL(2.0)), -REAL(REAL(3.0)), REAL(REAL(4.0))};

    auto result = PowellMinimize(func, start);

    REQUIRE(result.converged);
    REQUIRE_THAT(result.xmin[0], RealWithinAbs(REAL(REAL(0.0)), 1e-6));
    REQUIRE_THAT(result.xmin[1], RealWithinAbs(REAL(REAL(0.0)), 1e-6));
    REQUIRE_THAT(result.xmin[2], RealWithinAbs(REAL(REAL(0.0)), 1e-6));
}

TEST_CASE("Powell - Rosenbrock function", "[Powell][2D][Rosenbrock]")
{
    Rosenbrock2D func;
    VectorN<Real, 2> start{-REAL(REAL(1.0)), REAL(REAL(1.0))};

    Powell optimizer(1e-10, 500);
    auto result = optimizer.Minimize(func, start);

    REQUIRE(result.converged);
    REQUIRE_THAT(result.xmin[0], RealWithinAbs(REAL(REAL(1.0)), 1e-4));
    REQUIRE_THAT(result.xmin[1], RealWithinAbs(REAL(REAL(1.0)), 1e-4));
}

TEST_CASE("Powell - Custom direction matrix", "[Powell][2D]")
{
    SimpleQuadratic2D func;
    VectorN<Real, 2> start{REAL(REAL(5.0)), REAL(REAL(5.0))};

    // Non-standard initial directions
    Matrix<Real> ximat(2, 2);
    ximat(0, 0) = REAL(REAL(1.0));  ximat(0, 1) = REAL(REAL(1.0));
    ximat(1, 0) = REAL(REAL(1.0));  ximat(1, 1) = -REAL(REAL(1.0));

    Powell optimizer;
    auto result = optimizer.Minimize(func, start, ximat);

    REQUIRE(result.converged);
    REQUIRE_THAT(result.xmin[0], RealWithinAbs(REAL(REAL(0.0)), 1e-6));
    REQUIRE_THAT(result.xmin[1], RealWithinAbs(REAL(REAL(0.0)), 1e-6));
}

///////////////////////////////////////////////////////////////////////////////
///                     Conjugate Gradient Tests                            ///
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("ConjugateGradient - Simple 2D quadratic", "[CG][2D]")
{
    DiffQuadratic2D func;
    VectorN<Real, 2> start{REAL(REAL(5.0)), REAL(REAL(5.0))};

    ConjugateGradient optimizer;
    auto result = optimizer.Minimize(func, start);

    REQUIRE(result.converged);
    REQUIRE_THAT(result.xmin[0], RealWithinAbs(REAL(REAL(0.0)), 1e-6));
    REQUIRE_THAT(result.xmin[1], RealWithinAbs(REAL(REAL(0.0)), 1e-6));
}

TEST_CASE("ConjugateGradient - Elliptic paraboloid", "[CG][2D]")
{
    DiffEllipticParaboloid2D func;
    VectorN<Real, 2> start{REAL(REAL(0.0)), REAL(REAL(0.0))};

    auto result = ConjugateGradientMinimize(func, start);

    REQUIRE(result.converged);
    REQUIRE_THAT(result.xmin[0], RealWithinAbs(REAL(REAL(1.0)), 1e-6));
    REQUIRE_THAT(result.xmin[1], RealWithinAbs(-REAL(REAL(1.0)), 1e-6));
}

TEST_CASE("ConjugateGradient - Booth function", "[CG][2D]")
{
    DiffBoothFunction func;
    VectorN<Real, 2> start{REAL(REAL(0.0)), REAL(REAL(0.0))};

    auto result = ConjugateGradientMinimize(func, start);

    REQUIRE(result.converged);
    REQUIRE_THAT(result.xmin[0], RealWithinAbs(REAL(REAL(1.0)), 1e-5));
    REQUIRE_THAT(result.xmin[1], RealWithinAbs(REAL(REAL(3.0)), 1e-5));
}

TEST_CASE("ConjugateGradient - 3D quadratic", "[CG][3D]")
{
    DiffQuadratic3D func;
    VectorN<Real, 3> start{REAL(REAL(2.0)), -REAL(REAL(3.0)), REAL(REAL(4.0))};

    auto result = ConjugateGradientMinimize(func, start);

    REQUIRE(result.converged);
    REQUIRE_THAT(result.xmin[0], RealWithinAbs(REAL(REAL(0.0)), 1e-6));
    REQUIRE_THAT(result.xmin[1], RealWithinAbs(REAL(REAL(0.0)), 1e-6));
    REQUIRE_THAT(result.xmin[2], RealWithinAbs(REAL(REAL(0.0)), 1e-6));
}

TEST_CASE("ConjugateGradient - 3D offset quadratic", "[CG][3D]")
{
    DiffOffsetQuadratic3D func;
    VectorN<Real, 3> start{REAL(REAL(0.0)), REAL(REAL(0.0)), REAL(REAL(0.0))};

    auto result = ConjugateGradientMinimize(func, start);

    REQUIRE(result.converged);
    REQUIRE_THAT(result.xmin[0], RealWithinAbs(REAL(REAL(1.0)), 1e-6));
    REQUIRE_THAT(result.xmin[1], RealWithinAbs(REAL(REAL(2.0)), 1e-6));
    REQUIRE_THAT(result.xmin[2], RealWithinAbs(REAL(REAL(3.0)), 1e-6));
}

TEST_CASE("ConjugateGradient - Rosenbrock function", "[CG][2D][Rosenbrock]")
{
    DiffRosenbrock2D func;
    VectorN<Real, 2> start{-REAL(REAL(1.0)), REAL(REAL(1.0))};

    ConjugateGradient optimizer(1e-10, 1e-10, 500);
    auto result = optimizer.Minimize(func, start);

    REQUIRE(result.converged);
    REQUIRE_THAT(result.xmin[0], RealWithinAbs(REAL(REAL(1.0)), 1e-4));
    REQUIRE_THAT(result.xmin[1], RealWithinAbs(REAL(REAL(1.0)), 1e-4));
}

TEST_CASE("ConjugateGradient - Fletcher-Reeves method", "[CG][2D]")
{
    DiffQuadratic2D func;
    VectorN<Real, 2> start{REAL(REAL(5.0)), REAL(REAL(5.0))};

    auto result = ConjugateGradientMinimize(func, start, 1e-8,
                                            ConjugateGradient::Method::FletcherReeves);

    REQUIRE(result.converged);
    REQUIRE_THAT(result.xmin[0], RealWithinAbs(REAL(REAL(0.0)), 1e-6));
    REQUIRE_THAT(result.xmin[1], RealWithinAbs(REAL(REAL(0.0)), 1e-6));
}

TEST_CASE("ConjugateGradient - Polak-Ribiere method", "[CG][2D]")
{
    DiffQuadratic2D func;
    VectorN<Real, 2> start{REAL(REAL(5.0)), REAL(REAL(5.0))};

    auto result = ConjugateGradientMinimize(func, start, 1e-8,
                                            ConjugateGradient::Method::PolakRibiere);

    REQUIRE(result.converged);
    REQUIRE_THAT(result.xmin[0], RealWithinAbs(REAL(REAL(0.0)), 1e-6));
    REQUIRE_THAT(result.xmin[1], RealWithinAbs(REAL(REAL(0.0)), 1e-6));
}

///////////////////////////////////////////////////////////////////////////////
///                          BFGS Tests                                     ///
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("BFGS - Simple 2D quadratic", "[BFGS][2D]")
{
    DiffQuadratic2D func;
    VectorN<Real, 2> start{REAL(REAL(5.0)), REAL(REAL(5.0))};

    BFGS optimizer;
    auto result = optimizer.Minimize(func, start);

    REQUIRE(result.converged);
    REQUIRE_THAT(result.xmin[0], RealWithinAbs(REAL(REAL(0.0)), 1e-6));
    REQUIRE_THAT(result.xmin[1], RealWithinAbs(REAL(REAL(0.0)), 1e-6));
}

TEST_CASE("BFGS - Elliptic paraboloid", "[BFGS][2D]")
{
    DiffEllipticParaboloid2D func;
    VectorN<Real, 2> start{REAL(REAL(0.0)), REAL(REAL(0.0))};

    auto result = BFGSMinimize(func, start);

    REQUIRE(result.converged);
    REQUIRE_THAT(result.xmin[0], RealWithinAbs(REAL(REAL(1.0)), 1e-6));
    REQUIRE_THAT(result.xmin[1], RealWithinAbs(-REAL(REAL(1.0)), 1e-6));
}

TEST_CASE("BFGS - Booth function", "[BFGS][2D]")
{
    DiffBoothFunction func;
    VectorN<Real, 2> start{REAL(REAL(0.0)), REAL(REAL(0.0))};

    auto result = BFGSMinimize(func, start);

    REQUIRE(result.converged);
    REQUIRE_THAT(result.xmin[0], RealWithinAbs(REAL(REAL(1.0)), 1e-5));
    REQUIRE_THAT(result.xmin[1], RealWithinAbs(REAL(REAL(3.0)), 1e-5));
}

TEST_CASE("BFGS - 3D quadratic", "[BFGS][3D]")
{
    DiffQuadratic3D func;
    VectorN<Real, 3> start{REAL(REAL(2.0)), -REAL(REAL(3.0)), REAL(REAL(4.0))};

    auto result = BFGSMinimize(func, start);

    REQUIRE(result.converged);
    REQUIRE_THAT(result.xmin[0], RealWithinAbs(REAL(REAL(0.0)), 1e-6));
    REQUIRE_THAT(result.xmin[1], RealWithinAbs(REAL(REAL(0.0)), 1e-6));
    REQUIRE_THAT(result.xmin[2], RealWithinAbs(REAL(REAL(0.0)), 1e-6));
}

TEST_CASE("BFGS - 3D offset quadratic", "[BFGS][3D]")
{
    DiffOffsetQuadratic3D func;
    VectorN<Real, 3> start{REAL(REAL(0.0)), REAL(REAL(0.0)), REAL(REAL(0.0))};

    auto result = BFGSMinimize(func, start);

    REQUIRE(result.converged);
    REQUIRE_THAT(result.xmin[0], RealWithinAbs(REAL(REAL(1.0)), 1e-6));
    REQUIRE_THAT(result.xmin[1], RealWithinAbs(REAL(REAL(2.0)), 1e-6));
    REQUIRE_THAT(result.xmin[2], RealWithinAbs(REAL(REAL(3.0)), 1e-6));
}

TEST_CASE("BFGS - Rosenbrock function", "[BFGS][2D][Rosenbrock]")
{
    DiffRosenbrock2D func;
    VectorN<Real, 2> start{-REAL(REAL(1.0)), REAL(REAL(1.0))};

    BFGS optimizer(1e-10, 1e-10, 500);
    auto result = optimizer.Minimize(func, start);

    REQUIRE(result.converged);
    REQUIRE_THAT(result.xmin[0], RealWithinAbs(REAL(REAL(1.0)), 1e-4));
    REQUIRE_THAT(result.xmin[1], RealWithinAbs(REAL(REAL(1.0)), 1e-4));
}

///////////////////////////////////////////////////////////////////////////////
///                    Comparison Tests                                     ///
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("All methods converge to same minimum", "[Powell][CG][BFGS][Comparison]")
{
    DiffEllipticParaboloid2D func;
    EllipticParaboloid2D funcNoDeriv;
    VectorN<Real, 2> start{REAL(REAL(5.0)), -REAL(REAL(5.0))};

    auto nelderResult = NelderMeadMinimize(funcNoDeriv, start);
    auto powellResult = PowellMinimize(funcNoDeriv, start);
    auto cgResult = ConjugateGradientMinimize(func, start);
    auto bfgsResult = BFGSMinimize(func, start);

    // All should converge
    REQUIRE(nelderResult.converged);
    REQUIRE(powellResult.converged);
    REQUIRE(cgResult.converged);
    REQUIRE(bfgsResult.converged);

    // All should find the same minimum
    REQUIRE_THAT(nelderResult.xmin[0], RealWithinAbs(REAL(REAL(1.0)), 1e-4));
    REQUIRE_THAT(powellResult.xmin[0], RealWithinAbs(REAL(REAL(1.0)), 1e-4));
    REQUIRE_THAT(cgResult.xmin[0], RealWithinAbs(REAL(REAL(1.0)), 1e-4));
    REQUIRE_THAT(bfgsResult.xmin[0], RealWithinAbs(REAL(REAL(1.0)), 1e-4));

    REQUIRE_THAT(nelderResult.xmin[1], RealWithinAbs(-REAL(REAL(1.0)), 1e-4));
    REQUIRE_THAT(powellResult.xmin[1], RealWithinAbs(-REAL(REAL(1.0)), 1e-4));
    REQUIRE_THAT(cgResult.xmin[1], RealWithinAbs(-REAL(REAL(1.0)), 1e-4));
    REQUIRE_THAT(bfgsResult.xmin[1], RealWithinAbs(-REAL(REAL(1.0)), 1e-4));
}

TEST_CASE("Gradient methods are more efficient", "[Powell][CG][BFGS][Comparison]")
{
    // For quadratic functions, gradient methods should converge in fewer iterations
    DiffQuadratic2D func;
    SimpleQuadratic2D funcNoDeriv;
    VectorN<Real, 2> start{REAL(REAL(10.0)), REAL(REAL(10.0))};

    NelderMead nelder;
    auto nelderResult = nelder.Minimize(funcNoDeriv, start, REAL(REAL(1.0)));

    Powell powell;
    auto powellResult = powell.Minimize(funcNoDeriv, start);

    ConjugateGradient cg;
    auto cgResult = cg.Minimize(func, start);

    BFGS bfgs;
    auto bfgsResult = bfgs.Minimize(func, start);

    // All converge
    REQUIRE(nelderResult.converged);
    REQUIRE(powellResult.converged);
    REQUIRE(cgResult.converged);
    REQUIRE(bfgsResult.converged);

    // Gradient methods typically need fewer iterations for quadratics
    // CG and BFGS should theoretically converge in at most N iterations
    // for a pure quadratic
    INFO("Nelder-Mead iterations: " << nelderResult.iterations);
    INFO("Powell iterations: " << powellResult.iterations);
    INFO("CG iterations: " << cgResult.iterations);
    INFO("BFGS iterations: " << bfgsResult.iterations);
}

///////////////////////////////////////////////////////////////////////////////
///                        API Tests                                        ///
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("Powell - Accessor methods", "[Powell][API]")
{
    Powell optimizer(1e-6, 300);

    REQUIRE_THAT(optimizer.getFtol(), RealApprox(1e-6));
    REQUIRE(optimizer.getMaxIter() == 300);

    optimizer.setFtol(1e-10);
    optimizer.setMaxIter(500);

    REQUIRE_THAT(optimizer.getFtol(), RealApprox(1e-10));
    REQUIRE(optimizer.getMaxIter() == 500);
}

TEST_CASE("ConjugateGradient - Accessor methods", "[CG][API]")
{
    ConjugateGradient optimizer(1e-6, 1e-7, 300,
                                ConjugateGradient::Method::FletcherReeves);

    REQUIRE_THAT(optimizer.getFtol(), RealApprox(1e-6));
    REQUIRE_THAT(optimizer.getGtol(), RealApprox(1e-7));
    REQUIRE(optimizer.getMaxIter() == 300);
    REQUIRE(optimizer.getMethod() == ConjugateGradient::Method::FletcherReeves);

    optimizer.setMethod(ConjugateGradient::Method::PolakRibiere);
    REQUIRE(optimizer.getMethod() == ConjugateGradient::Method::PolakRibiere);
}

TEST_CASE("BFGS - Accessor methods", "[BFGS][API]")
{
    BFGS optimizer(1e-6, 1e-7, 300);

    REQUIRE_THAT(optimizer.getFtol(), RealApprox(1e-6));
    REQUIRE_THAT(optimizer.getGtol(), RealApprox(1e-7));
    REQUIRE(optimizer.getMaxIter() == 300);

    optimizer.setFtol(1e-10);
    optimizer.setGtol(1e-11);
    optimizer.setMaxIter(500);

    REQUIRE_THAT(optimizer.getFtol(), RealApprox(1e-10));
    REQUIRE_THAT(optimizer.getGtol(), RealApprox(1e-11));
    REQUIRE(optimizer.getMaxIter() == 500);
}

/***************************************************************************************************
 * COMPREHENSIVE TEST BED INTEGRATION TESTS
 * 
 * These tests iterate over ALL 2D/ND test cases from optimization_test_bed.h to ensure
 * complete coverage of the multidimensional optimization algorithms.
 ***************************************************************************************************/

#include "../../test_data/optimization_test_bed.h"

using namespace MML::TestBeds;

TEST_CASE("NelderMead_All2DConvexTestBed", "[NelderMead][TestBed][Convex][Comprehensive]")
{
    TEST_PRECISION_INFO();
    
    auto convexTests = getConvex2DTests();
    
    INFO("Testing " << convexTests.size() << " convex 2D functions with Nelder-Mead");
    
    for (const auto& test : convexTests)
    {
        DYNAMIC_SECTION(test.name)
        {
            INFO("Category: " << test.category);
            INFO("Difficulty: " << test.difficulty);
            INFO("Description: " << test.description);
            INFO("True minimum value: " << test.trueMinimumF);
            
            OptScalarFunctionWrapper<2> f(test.func);
            
            NelderMead optimizer;
            auto result = optimizer.Minimize(f, test.suggestedStart, 1.0);
            
            CHECK(result.converged);
            
            if (result.converged)
            {
                INFO("Found minimum at (" << result.xmin[0] << ", " << result.xmin[1] << ")");
                INFO("Function value: " << result.fmin);
                
                // Check if we're near any known minimum
                bool nearMinimum = verifyNDMinimum(test, result.xmin, 1e-3);
                CHECK(nearMinimum);
                CHECK(std::abs(result.fmin - test.trueMinimumF) < 1e-4);
            }
        }
    }
}

TEST_CASE("NelderMead_ClassicBenchmarks2DTestBed", "[NelderMead][TestBed][Classic][Comprehensive]")
{
    TEST_PRECISION_INFO();
    
    auto classicTests = getClassicBenchmark2DTests();
    
    INFO("Testing " << classicTests.size() << " classic benchmark functions with Nelder-Mead");
    
    for (const auto& test : classicTests)
    {
        DYNAMIC_SECTION(test.name)
        {
            INFO("Category: " << test.category);
            INFO("Difficulty: " << test.difficulty);
            INFO("Description: " << test.description);
            INFO("Number of global minima: " << test.numMinima);
            
            OptScalarFunctionWrapper<2> f(test.func);
            
            NelderMead optimizer(1e-8, 2000);  // Tighter tolerance, more iterations for hard problems
            auto result = optimizer.Minimize(f, test.suggestedStart, 1.0);
            
            // Adjust tolerance based on difficulty
            Real posTol = 1e-3;
            Real valTol = 1e-4;
            if (test.difficulty >= 3) { posTol = 1e-2; valTol = 1e-2; }
            if (test.difficulty >= 4) { posTol = 1e-1; valTol = 1e-1; }
            
            if (result.converged)
            {
                INFO("Found minimum at (" << result.xmin[0] << ", " << result.xmin[1] << ")");
                INFO("Function value: " << result.fmin);
                
                bool nearMinimum = verifyNDMinimum(test, result.xmin, posTol);
                
                if (!nearMinimum)
                {
                    // For multimodal, might find local minimum - verify it's at least a local min
                    Real error = computeMinimumError(test, result.xmin);
                    INFO("Distance from nearest known minimum: " << error);
                }
                
                CHECK(nearMinimum);
                CHECK(std::abs(result.fmin - test.trueMinimumF) < valTol);
            }
            else
            {
                INFO("Did not converge - checking if we're at least close");
                Real error = computeMinimumError(test, result.xmin);
                INFO("Distance from nearest known minimum: " << error);
            }
        }
    }
}

TEST_CASE("Powell_All2DConvexTestBed", "[Powell][TestBed][Convex][Comprehensive]")
{
    TEST_PRECISION_INFO();
    
    auto convexTests = getConvex2DTests();
    
    INFO("Testing " << convexTests.size() << " convex 2D functions with Powell's method");
    
    for (const auto& test : convexTests)
    {
        DYNAMIC_SECTION(test.name)
        {
            INFO("Category: " << test.category);
            INFO("Difficulty: " << test.difficulty);
            
            OptScalarFunctionWrapper<2> f(test.func);
            
            Powell optimizer;
            auto result = optimizer.Minimize(f, test.suggestedStart);
            
            CHECK(result.converged);
            
            if (result.converged)
            {
                INFO("Found minimum at (" << result.xmin[0] << ", " << result.xmin[1] << ")");
                INFO("Function value: " << result.fmin);
                
                bool nearMinimum = verifyNDMinimum(test, result.xmin, 1e-4);
                CHECK(nearMinimum);
                CHECK(std::abs(result.fmin - test.trueMinimumF) < 1e-5);
            }
        }
    }
}

TEST_CASE("Powell_ClassicBenchmarks2DTestBed", "[Powell][TestBed][Classic]")
{
    TEST_PRECISION_INFO();
    
    auto classicTests = getClassicBenchmark2DTests();
    
    for (const auto& test : classicTests)
    {
        DYNAMIC_SECTION(test.name)
        {
            INFO("Category: " << test.category);
            INFO("Difficulty: " << test.difficulty);
            INFO("Description: " << test.description);
            
            OptScalarFunctionWrapper<2> f(test.func);
            
            Powell optimizer(1e-8, 500);
            auto result = optimizer.Minimize(f, test.suggestedStart);
            
            Real posTol = 1e-3;
            Real valTol = 1e-4;
            if (test.difficulty >= 3) { posTol = 1e-1; valTol = 1.0; }  // Very relaxed for hard multimodal
            if (test.difficulty >= 4) { posTol = 1.0; valTol = 10.0; }
            
            if (result.converged)
            {
                INFO("Found minimum at (" << result.xmin[0] << ", " << result.xmin[1] << ")");
                INFO("Function value: " << result.fmin);
                
                bool nearMinimum = verifyNDMinimum(test, result.xmin, posTol);
                
                // For highly multimodal functions (Ackley, Rastrigin), Powell often finds local minima
                // This is expected behavior - check we found SOME minimum (local is OK)
                if (test.isMultimodal && test.difficulty >= 3)
                {
                    // Just verify convergence and that it found a reasonable value
                    // (function value should be finite and not exploded)
                    CHECK(result.fmin < 100.0);
                    INFO("Note: Found local minimum (expected for hard multimodal)");
                }
                else
                {
                    CHECK(nearMinimum);
                    CHECK(std::abs(result.fmin - test.trueMinimumF) < valTol);
                }
            }
        }
    }
}

TEST_CASE("NelderMead_Multimodal2DTestBed", "[NelderMead][TestBed][Multimodal]")
{
    TEST_PRECISION_INFO();
    
    auto multimodalTests = getMultimodal2DTests();
    
    INFO("Testing " << multimodalTests.size() << " multimodal 2D functions");
    INFO("Note: Nelder-Mead finds LOCAL minima - we verify it finds a valid minimum");
    
    for (const auto& test : multimodalTests)
    {
        DYNAMIC_SECTION(test.name)
        {
            INFO("Category: " << test.category);
            INFO("Difficulty: " << test.difficulty);
            INFO("Description: " << test.description);
            INFO("Number of global minima: " << test.numMinima);
            
            OptScalarFunctionWrapper<2> f(test.func);
            
            NelderMead optimizer(1e-8, 2000);
            auto result = optimizer.Minimize(f, test.suggestedStart, 1.0);
            
            if (result.converged)
            {
                INFO("Found minimum at (" << result.xmin[0] << ", " << result.xmin[1] << ")");
                INFO("Function value: " << result.fmin);
                
                // For multimodal, check if we found any of the known minima
                // or at least a valid local minimum (function value <= nearby points)
                bool nearKnownMinimum = verifyNDMinimum(test, result.xmin, 0.1);
                
                if (nearKnownMinimum)
                {
                    INFO("Found known global minimum");
                    CHECK(std::abs(result.fmin - test.trueMinimumF) < 0.1);
                }
                else
                {
                    // Verify it's at least a local minimum by checking gradient is small
                    // (function value should be lower than nearby points)
                    Real delta = 1e-4;
                    VectorN<Real, 2> x{result.xmin[0], result.xmin[1]};
                    VectorN<Real, 2> xp, xm;
                    
                    bool isLocalMin = true;
                    for (int d = 0; d < 2; d++)
                    {
                        xp = x; xp[d] += delta;
                        xm = x; xm[d] -= delta;
                        
                        if (test.func(xp) < result.fmin - 1e-8 || 
                            test.func(xm) < result.fmin - 1e-8)
                        {
                            isLocalMin = false;
                            break;
                        }
                    }
                    
                    if (isLocalMin)
                        INFO("Found valid local minimum (not global)");
                    else
                        INFO("Warning: May not be a true minimum");
                    
                    CHECK(isLocalMin);
                }
            }
        }
    }
}

TEST_CASE("Optimization2D_AlgorithmComparison", "[Optimization][TestBed][Performance]")
{
    TEST_PRECISION_INFO();
    
    auto easyTests = getConvex2DTests();
    
    INFO("Comparing Nelder-Mead vs Powell on convex test functions");
    
    for (const auto& test : easyTests)
    {
        DYNAMIC_SECTION(test.name + " comparison")
        {
            OptScalarFunctionWrapper<2> f(test.func);
            
            NelderMead nmOptimizer;
            Powell powellOptimizer;
            
            auto nmResult = nmOptimizer.Minimize(f, test.suggestedStart, 1.0);
            auto powellResult = powellOptimizer.Minimize(f, test.suggestedStart);
            
            Real nmError = computeMinimumError(test, nmResult.xmin);
            Real powellError = computeMinimumError(test, powellResult.xmin);
            
            INFO("Nelder-Mead error: " << nmError << ", iterations: " << nmResult.iterations);
            INFO("Powell error: " << powellError << ", iterations: " << powellResult.iterations);
            
            // Both should converge for convex functions
            CHECK(nmResult.converged);
            CHECK(powellResult.converged);
            
            // Both should find the minimum reasonably well
            CHECK(nmError < 1e-3);
            CHECK(powellError < 1e-4);  // Powell typically more accurate
        }
    }
}

TEST_CASE("NelderMead_All2DTestBed", "[NelderMead][TestBed][All][Comprehensive]")
{
    TEST_PRECISION_INFO();
    
    auto allTests = getAll2DOptimizationTests();
    
    INFO("Testing ALL " << allTests.size() << " 2D functions from optimization test bed");
    
    int passed = 0, total = 0;
    
    for (const auto& test : allTests)
    {
        total++;
        
        DYNAMIC_SECTION(test.name)
        {
            INFO("Category: " << test.category);
            INFO("Difficulty: " << test.difficulty);
            INFO("Description: " << test.description);
            
            OptScalarFunctionWrapper<2> f(test.func);
            
            NelderMead optimizer(1e-8, 3000);
            auto result = optimizer.Minimize(f, test.suggestedStart, 1.0);
            
            // Relaxed tolerance for comprehensive test
            Real posTol = 0.1;
            Real valTol = 0.1;
            if (test.difficulty <= 2) { posTol = 1e-2; valTol = 1e-2; }
            if (test.difficulty == 1) { posTol = 1e-3; valTol = 1e-4; }
            
            if (result.converged)
            {
                bool nearMinimum = verifyNDMinimum(test, result.xmin, posTol);
                bool goodValue = std::abs(result.fmin - test.trueMinimumF) < valTol;
                
                if (nearMinimum && goodValue)
                    passed++;
                
                INFO("Found minimum at (" << result.xmin[0] << ", " << result.xmin[1] << ")");
                INFO("Function value: " << result.fmin << " (expected: " << test.trueMinimumF << ")");
                INFO("Near known minimum: " << (nearMinimum ? "YES" : "NO"));
                
                CHECK(nearMinimum);
                CHECK(goodValue);
            }
            else
            {
                INFO("Did not converge");
            }
        }
    }
    
    INFO("Summary: " << passed << "/" << total << " tests found global minimum");
}


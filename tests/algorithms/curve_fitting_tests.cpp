#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../TestPrecision.h"
#include "../TestMatchers.h"

#include "../../mml/algorithms/CurveFitting.h"

using namespace MML;
using namespace MML::Testing;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

TEST_CASE("LinearLeastSquares - Exact fit through two points", "[CurveFitting][LinearLeastSquares]")
{
    Vector<Real> x_data({0.0, 1.0});
    Vector<Real> y_data({1.0, 3.0});
    
    Real a, b;
    Real residual = LinearLeastSquares(x_data, y_data, a, b);
    
    REQUIRE_THAT(a, WithinAbs(REAL(2.0), REAL(1e-10)));
    REQUIRE_THAT(b, WithinAbs(REAL(1.0), REAL(1e-10)));
    REQUIRE_THAT(residual, WithinAbs(REAL(0.0), REAL(1e-10)));
}

TEST_CASE("LinearLeastSquares - Horizontal line (zero slope)", "[CurveFitting][LinearLeastSquares]")
{
    Vector<Real> x_data({-2.0, -1.0, 0.0, 1.0, 2.0});
    Vector<Real> y_data({5.0, 5.0, 5.0, 5.0, 5.0});
    
    Real a, b;
    Real residual = LinearLeastSquares(x_data, y_data, a, b);
    
    REQUIRE_THAT(a, WithinAbs(REAL(0.0), REAL(1e-10)));
    REQUIRE_THAT(b, WithinAbs(REAL(5.0), REAL(1e-10)));
    REQUIRE_THAT(residual, WithinAbs(REAL(0.0), REAL(1e-10)));
}

TEST_CASE("LinearLeastSquares - Exact linear data y = 2x + 3", "[CurveFitting][LinearLeastSquares]")
{
    Vector<Real> x_data({-2.0, -1.0, 0.0, 1.0, 2.0, 3.0});
    Vector<Real> y_data({-1.0, 1.0, 3.0, 5.0, 7.0, 9.0});
    
    Real a, b;
    Real residual = LinearLeastSquares(x_data, y_data, a, b);
    
    REQUIRE_THAT(a, WithinAbs(REAL(2.0), REAL(1e-10)));
    REQUIRE_THAT(b, WithinAbs(REAL(3.0), REAL(1e-10)));
    REQUIRE_THAT(residual, WithinAbs(REAL(0.0), REAL(1e-10)));
}

TEST_CASE("LinearLeastSquares - Noisy data approximation", "[CurveFitting][LinearLeastSquares]")
{
    // True line: y = 0.5x + 2.0, with small noise
    Vector<Real> x_data({1.0, 2.0, 3.0, 4.0, 5.0});
    Vector<Real> y_data({2.5, 3.1, 3.4, 4.2, 4.5});  // y = 0.5x + 2 + noise
    
    Real a, b;
    Real residual = LinearLeastSquares(x_data, y_data, a, b);
    
    // Should be close to true values but not exact due to noise
    REQUIRE_THAT(a, WithinAbs(REAL(0.5), REAL(0.15)));   // Within 30% of true value
    REQUIRE_THAT(b, WithinAbs(REAL(2.0), REAL(0.3)));
    REQUIRE(residual > 0.0);  // Should have non-zero residual
    REQUIRE(residual < 1.0);  // But should be small
}

TEST_CASE("LinearLeastSquares - Textbook example", "[CurveFitting][LinearLeastSquares]")
{
    // Classic example from numerical analysis textbooks
    // Data points roughly following y = 2x - 1
    Vector<Real> x_data({0.0, 1.0, 2.0, 3.0, 4.0});
    Vector<Real> y_data({-1.1, 0.9, 3.2, 5.1, 6.8});
    
    Real a, b;
    Real residual = LinearLeastSquares(x_data, y_data, a, b);
    
    // Expected: a ≈ 2.0, b ≈ -1.0
    REQUIRE_THAT(a, WithinAbs(REAL(2.0), REAL(0.1)));
    REQUIRE_THAT(b, WithinAbs(REAL(-1.0), REAL(0.2)));
    
    // Verify predictions are reasonable
    Real y_pred_0 = a * 0.0 + b;
    Real y_pred_4 = a * 4.0 + b;
    
    REQUIRE_THAT(y_pred_0, WithinAbs(REAL(-1.0), REAL(0.5)));
    REQUIRE_THAT(y_pred_4, WithinAbs(REAL(7.0), REAL(0.5)));
}

TEST_CASE("LinearLeastSquares - Negative slope", "[CurveFitting][LinearLeastSquares]")
{
    // y = -3x + 10
    Vector<Real> x_data({0.0, 1.0, 2.0, 3.0});
    Vector<Real> y_data({10.0, 7.0, 4.0, 1.0});
    
    Real a, b;
    Real residual = LinearLeastSquares(x_data, y_data, a, b);
    
    REQUIRE_THAT(a, WithinAbs(REAL(-3.0), REAL(1e-10)));
    REQUIRE_THAT(b, WithinAbs(REAL(10.0), REAL(1e-10)));
    REQUIRE_THAT(residual, WithinAbs(REAL(0.0), REAL(1e-10)));
}

TEST_CASE("LinearLeastSquares - Large dataset", "[CurveFitting][LinearLeastSquares]")
{
    // Generate 100 points following y = 1.5x - 2.5
    Vector<Real> x_data(100);
    Vector<Real> y_data(100);
    
    for (int i = 0; i < 100; i++) {
        x_data[i] = -10.0 + i * 0.2;  // x from -10 to 9.8
        y_data[i] = 1.5 * x_data[i] - 2.5;
    }
    
    Real a, b;
    Real residual = LinearLeastSquares(x_data, y_data, a, b);
    
    REQUIRE_THAT(a, WithinAbs(REAL(1.5), REAL(1e-10)));
    REQUIRE_THAT(b, WithinAbs(REAL(-2.5), REAL(1e-10)));
    REQUIRE_THAT(residual, WithinAbs(REAL(0.0), REAL(1e-8)));
}

TEST_CASE("LinearLeastSquares - Single point (degenerate case)", "[CurveFitting][LinearLeastSquares]")
{
    Vector<Real> x_data({3.5});
    Vector<Real> y_data({7.2});
    
    Real a, b;
    Real residual = LinearLeastSquares(x_data, y_data, a, b);
    
    // Should return horizontal line through the point
    REQUIRE_THAT(a, WithinAbs(REAL(0.0), REAL(1e-10)));
    REQUIRE_THAT(b, WithinAbs(REAL(7.2), REAL(1e-10)));
    REQUIRE_THAT(residual, WithinAbs(REAL(0.0), REAL(1e-10)));
}

TEST_CASE("LinearLeastSquares - Error handling: empty data", "[CurveFitting][LinearLeastSquares]")
{
    Vector<Real> x_data;
    Vector<Real> y_data;
    
    Real a, b;
    REQUIRE_THROWS_AS(LinearLeastSquares(x_data, y_data, a, b), std::invalid_argument);
}

TEST_CASE("LinearLeastSquares - Error handling: size mismatch", "[CurveFitting][LinearLeastSquares]")
{
    Vector<Real> x_data({1.0, 2.0, 3.0});
    Vector<Real> y_data({1.0, 2.0});  // Different size!
    
    Real a, b;
    REQUIRE_THROWS_AS(LinearLeastSquares(x_data, y_data, a, b), std::invalid_argument);
}

TEST_CASE("LinearLeastSquares - Error handling: identical x values (singular)", "[CurveFitting][LinearLeastSquares]")
{
    Vector<Real> x_data({5.0, 5.0, 5.0, 5.0});  // All x values identical
    Vector<Real> y_data({1.0, 2.0, 3.0, 4.0});
    
    Real a, b;
    REQUIRE_THROWS_AS(LinearLeastSquares(x_data, y_data, a, b), SingularMatrixError);
}

TEST_CASE("LinearLeastSquares - Vertical points (nearly singular)", "[CurveFitting][LinearLeastSquares]")
{
    Vector<Real> x_data({2.0, 2.0});  // Two identical x values
    Vector<Real> y_data({1.0, 5.0});
    
    Real a, b;
    Real residual = LinearLeastSquares(x_data, y_data, a, b);
    
    // Should return horizontal line at average y
    REQUIRE_THAT(a, WithinAbs(REAL(0.0), REAL(1e-10)));
    REQUIRE_THAT(b, WithinAbs(REAL(3.0), REAL(1e-10)));  // (1 + 5) / 2
    // residual = sqrt((1-3)^2 + (5-3)^2) = sqrt(4 + 4) = sqrt(8) = 2.828...
    REQUIRE_THAT(residual, WithinAbs(std::sqrt(8.0), REAL(1e-10)));
}


// Tests for detailed statistics version
TEST_CASE("LinearLeastSquaresDetailed - Perfect fit R^2 = 1", "[CurveFitting][LinearLeastSquaresDetailed]")
{
    Vector<Real> x_data({1.0, 2.0, 3.0, 4.0, 5.0});
    Vector<Real> y_data({2.0, 4.0, 6.0, 8.0, 10.0});  // y = 2x
    
    auto result = LinearLeastSquaresDetailed(x_data, y_data);
    
    REQUIRE_THAT(result.a, WithinAbs(REAL(2.0), REAL(1e-10)));
    REQUIRE_THAT(result.b, WithinAbs(REAL(0.0), REAL(1e-10)));
    REQUIRE_THAT(result.residual_norm, WithinAbs(REAL(0.0), REAL(1e-10)));
    REQUIRE_THAT(result.r_squared, WithinAbs(REAL(1.0), REAL(1e-10)));
    REQUIRE_THAT(result.mean_squared_error, WithinAbs(REAL(0.0), REAL(1e-10)));
}

TEST_CASE("LinearLeastSquaresDetailed - Noisy data R^2 < 1", "[CurveFitting][LinearLeastSquaresDetailed]")
{
    Vector<Real> x_data({1.0, 2.0, 3.0, 4.0, 5.0});
    Vector<Real> y_data({2.5, 3.8, 6.2, 7.9, 10.1});  // y ≈ 2x with noise
    
    auto result = LinearLeastSquaresDetailed(x_data, y_data);
    
    REQUIRE_THAT(result.a, WithinAbs(REAL(2.0), REAL(0.15)));
    REQUIRE(result.r_squared > 0.95);  // Should be very good fit
    REQUIRE(result.r_squared < 1.0);   // But not perfect
    REQUIRE(result.mean_squared_error > 0.0);
    REQUIRE(result.mean_squared_error < 0.5);
}

TEST_CASE("LinearLeastSquaresDetailed - Constant y values R^2 handling", "[CurveFitting][LinearLeastSquaresDetailed]")
{
    Vector<Real> x_data({1.0, 2.0, 3.0, 4.0});
    Vector<Real> y_data({5.0, 5.0, 5.0, 5.0});  // All y identical
    
    auto result = LinearLeastSquaresDetailed(x_data, y_data);
    
    REQUIRE_THAT(result.a, WithinAbs(REAL(0.0), REAL(1e-10)));
    REQUIRE_THAT(result.b, WithinAbs(REAL(5.0), REAL(1e-10)));
    REQUIRE_THAT(result.r_squared, WithinAbs(REAL(1.0), REAL(1e-10)));  // Perfect fit to horizontal line
}

TEST_CASE("LinearLeastSquaresDetailed - Poor fit example", "[CurveFitting][LinearLeastSquaresDetailed]")
{
    // Data that doesn't follow linear pattern well
    Vector<Real> x_data({1.0, 2.0, 3.0, 4.0, 5.0});
    Vector<Real> y_data({1.0, 4.0, 2.0, 8.0, 3.0});  // Very scattered
    
    auto result = LinearLeastSquaresDetailed(x_data, y_data);
    
    // Should have low R^2
    REQUIRE(result.r_squared < 0.5);  // Poor fit
    REQUIRE(result.mean_squared_error > 1.0);  // High error
}


// =============================================================================
// General Linear Least Squares Tests
// =============================================================================

// Simple helper classes for basis functions
class ConstantFunction : public IRealFunction {
public:
    Real operator()(Real) const override { return REAL(1.0); }
};

class LinearFunction : public IRealFunction {
public:
    Real operator()(Real x) const override { return x; }
};

class QuadraticFunction : public IRealFunction {
public:
    Real operator()(Real x) const override { return x * x; }
};

class CubicFunction : public IRealFunction {
public:
    Real operator()(Real x) const override { return x * x * x; }
};

class SineFunction : public IRealFunction {
public:
    Real operator()(Real x) const override { return std::sin(x); }
};

class CosineFunction : public IRealFunction {
public:
    Real operator()(Real x) const override { return std::cos(x); }
};

class ExponentialFunction : public IRealFunction {
public:
    Real operator()(Real x) const override { return std::exp(x); }
};


TEST_CASE("GeneralLinearLeastSquares - Constant function fit", "[CurveFitting][GeneralLLS]")
{
    Vector<Real> x_data({0.0, 1.0, 2.0, 3.0, 4.0});
    Vector<Real> y_data({5.0, 5.0, 5.0, 5.0, 5.0});  // Constant y = 5
    
    ConstantFunction f0;
    Vector<const IRealFunction*> basis({&f0});
    
    auto result = GeneralLinearLeastSquares(x_data, y_data, basis);
    
    REQUIRE(result.coefficients.size() == 1);
    REQUIRE_THAT(result.coefficients[0], WithinAbs(REAL(5.0), REAL(1e-10)));
    REQUIRE_THAT(result.residual_norm, WithinAbs(REAL(0.0), REAL(1e-10)));
    REQUIRE_THAT(result.r_squared, WithinAbs(REAL(1.0), REAL(1e-10)));
}

TEST_CASE("GeneralLinearLeastSquares - Linear fit with IRealFunction", "[CurveFitting][GeneralLLS]")
{
    // y = 2x + 3
    Vector<Real> x_data({0.0, 1.0, 2.0, 3.0, 4.0});
    Vector<Real> y_data({3.0, 5.0, 7.0, 9.0, 11.0});
    
    ConstantFunction f0;  // 1
    LinearFunction f1;    // x
    Vector<const IRealFunction*> basis({&f0, &f1});
    
    auto result = GeneralLinearLeastSquares(x_data, y_data, basis);
    
    REQUIRE(result.coefficients.size() == 2);
    REQUIRE_THAT(result.coefficients[0], WithinAbs(REAL(3.0), REAL(1e-10)));  // constant term
    REQUIRE_THAT(result.coefficients[1], WithinAbs(REAL(2.0), REAL(1e-10)));  // linear term
    REQUIRE_THAT(result.residual_norm, WithinAbs(REAL(0.0), REAL(1e-10)));
    REQUIRE_THAT(result.r_squared, WithinAbs(REAL(1.0), REAL(1e-10)));
}

TEST_CASE("GeneralLinearLeastSquares - Quadratic fit", "[CurveFitting][GeneralLLS]")
{
    // y = x^2 - 2x + 1 = (x-1)^2
    Vector<Real> x_data({-1.0, 0.0, 1.0, 2.0, 3.0});
    Vector<Real> y_data({4.0, 1.0, 0.0, 1.0, 4.0});  // parabola
    
    ConstantFunction f0;
    LinearFunction f1;
    QuadraticFunction f2;
    Vector<const IRealFunction*> basis({&f0, &f1, &f2});
    
    auto result = GeneralLinearLeastSquares(x_data, y_data, basis);
    
    REQUIRE(result.coefficients.size() == 3);
    REQUIRE_THAT(result.coefficients[0], WithinAbs(REAL(1.0), REAL(1e-9)));   // constant
    REQUIRE_THAT(result.coefficients[1], WithinAbs(REAL(-2.0), REAL(1e-9)));  // linear
    REQUIRE_THAT(result.coefficients[2], WithinAbs(REAL(1.0), REAL(1e-9)));   // quadratic
    REQUIRE_THAT(result.residual_norm, WithinAbs(REAL(0.0), REAL(1e-9)));
}

TEST_CASE("GeneralLinearLeastSquares - Cubic fit with overdetermined system", "[CurveFitting][GeneralLLS]")
{
    // y = x^3 - x
    Vector<Real> x_data({-2.0, -1.0, 0.0, 1.0, 2.0, 3.0});
    Vector<Real> y_data({-6.0, 0.0, 0.0, 0.0, 6.0, 24.0});
    
    ConstantFunction f0;
    LinearFunction f1;
    QuadraticFunction f2;
    CubicFunction f3;
    Vector<const IRealFunction*> basis({&f0, &f1, &f2, &f3});
    
    auto result = GeneralLinearLeastSquares(x_data, y_data, basis);
    
    REQUIRE(result.coefficients.size() == 4);
    REQUIRE_THAT(result.coefficients[0], WithinAbs(REAL(0.0), REAL(1e-9)));   // constant = 0
    REQUIRE_THAT(result.coefficients[1], WithinAbs(REAL(-1.0), REAL(1e-9)));  // x term = -1
    REQUIRE_THAT(result.coefficients[2], WithinAbs(REAL(0.0), REAL(1e-9)));   // x^2 term = 0
    REQUIRE_THAT(result.coefficients[3], WithinAbs(REAL(1.0), REAL(1e-9)));   // x^3 term = 1
    REQUIRE_THAT(result.residual_norm, WithinAbs(REAL(0.0), REAL(1e-9)));
}

TEST_CASE("GeneralLinearLeastSquares - Trigonometric fit", "[CurveFitting][GeneralLLS]")
{
    // y = 2 + 3*sin(x) - cos(x)
    const Real PI = REAL(3.14159265358979323846);
    Vector<Real> x_data(20);
    Vector<Real> y_data(20);
    
    for (int i = 0; i < 20; i++) {
        x_data[i] = -PI + i * (REAL(2) * PI / REAL(19.0));
        y_data[i] = REAL(2.0) + REAL(3.0) * std::sin(x_data[i]) - std::cos(x_data[i]);
    }
    
    ConstantFunction f0;
    SineFunction f1;
    CosineFunction f2;
    Vector<const IRealFunction*> basis({&f0, &f1, &f2});
    
    auto result = GeneralLinearLeastSquares(x_data, y_data, basis);
    
    REQUIRE(result.coefficients.size() == 3);
    REQUIRE_THAT(result.coefficients[0], WithinAbs(REAL(2.0), REAL(1e-9)));   // constant
    REQUIRE_THAT(result.coefficients[1], WithinAbs(REAL(3.0), REAL(1e-9)));   // sin term
    REQUIRE_THAT(result.coefficients[2], WithinAbs(REAL(-1.0), REAL(1e-9)));  // cos term
    REQUIRE_THAT(result.r_squared, WithinAbs(REAL(1.0), REAL(1e-10)));
}

TEST_CASE("GeneralLinearLeastSquares - Lambda basis functions", "[CurveFitting][GeneralLLS]")
{
    // y = 3*exp(x) + 2*x
    Vector<Real> x_data({0.0, 0.5, 1.0, 1.5, 2.0});
    Vector<Real> y_data(5);
    for (int i = 0; i < 5; i++) {
        y_data[i] = 3.0 * std::exp(x_data[i]) + 2.0 * x_data[i];
    }
    
    Vector<std::function<Real(Real)>> basis({
        [](Real x) -> Real { return std::exp(x); },
        [](Real x) -> Real { return x; }
    });
    
    auto result = GeneralLinearLeastSquares(x_data, y_data, basis);
    
    REQUIRE(result.coefficients.size() == 2);
    REQUIRE_THAT(result.coefficients[0], WithinAbs(REAL(3.0), REAL(1e-8)));  // exp term
    REQUIRE_THAT(result.coefficients[1], WithinAbs(REAL(2.0), REAL(1e-8)));  // linear term
    REQUIRE_THAT(result.r_squared, WithinAbs(REAL(1.0), REAL(1e-10)));
}

TEST_CASE("GeneralLinearLeastSquares - Noisy data", "[CurveFitting][GeneralLLS]")
{
    // Noisy quadratic data: y = x^2 + noise
    Vector<Real> x_data({-2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0});
    Vector<Real> y_data({4.1, 1.2, -0.1, 0.9, 4.2, 8.8, 16.1});  // x^2 + small noise
    
    Vector<std::function<Real(Real)>> basis({
        [](double) { return 1.0; },
        [](Real x) -> Real { return x; },
        [](Real x) -> Real { return x * x; }
    });
    
    auto result = GeneralLinearLeastSquares(x_data, y_data, basis);
    
    REQUIRE(result.coefficients.size() == 3);
    REQUIRE_THAT(result.coefficients[0], WithinAbs(REAL(0.0), REAL(0.5)));  // constant ≈ 0
    REQUIRE_THAT(result.coefficients[1], WithinAbs(REAL(0.0), REAL(0.5)));  // linear ≈ 0
    REQUIRE_THAT(result.coefficients[2], WithinAbs(REAL(1.0), REAL(0.2)));  // quadratic ≈ 1
    REQUIRE(result.r_squared > 0.99);  // Good fit despite noise
    REQUIRE(result.residual_norm > 0.0);  // Non-zero due to noise
}

TEST_CASE("GeneralLinearLeastSquares - Statistics validation", "[CurveFitting][GeneralLLS]")
{
    Vector<Real> x_data({0.0, 1.0, 2.0, 3.0, 4.0});
    Vector<Real> y_data({1.0, 2.0, 3.0, 4.0, 5.0});  // y = x + 1
    
    Vector<std::function<Real(Real)>> basis({
        [](double) { return 1.0; },
        [](Real x) -> Real { return x; }
    });
    
    auto result = GeneralLinearLeastSquares(x_data, y_data, basis);
    
    REQUIRE(result.num_data_points == 5);
    REQUIRE(result.num_basis_functions == 2);
    REQUIRE(result.effective_rank == 2);  // Full rank
    REQUIRE(result.condition_number < 100.0);  // Well-conditioned
    REQUIRE(result.adjusted_r_squared <= result.r_squared);  // Adjusted is always <= R²
}

TEST_CASE("GeneralLinearLeastSquares - Evaluate fitted function", "[CurveFitting][GeneralLLS]")
{
    // y = 2 + 3x - x^2
    Vector<Real> x_data({0.0, 1.0, 2.0, 3.0, 4.0});
    Vector<Real> y_data({2.0, 4.0, 4.0, 2.0, -2.0});
    
    ConstantFunction f0;
    LinearFunction f1;
    QuadraticFunction f2;
    Vector<const IRealFunction*> basis({&f0, &f1, &f2});
    
    auto result = GeneralLinearLeastSquares(x_data, y_data, basis);
    
    // Verify coefficients
    REQUIRE_THAT(result.coefficients[0], WithinAbs(REAL(2.0), REAL(1e-9)));
    REQUIRE_THAT(result.coefficients[1], WithinAbs(REAL(3.0), REAL(1e-9)));
    REQUIRE_THAT(result.coefficients[2], WithinAbs(REAL(-1.0), REAL(1e-9)));
    
    // Test evaluate function at data points
    for (int i = 0; i < 5; i++) {
        Real y_pred = result.evaluate(x_data[i], basis);
        REQUIRE_THAT(y_pred, WithinAbs(y_data[i], REAL(1e-9)));
    }
    
    // Test at a new point x = 2.5: y = 2 + 3*2.5 - 2.5^2 = 2 + 7.5 - 6.25 = 3.25
    REQUIRE_THAT(result.evaluate(2.5, basis), WithinAbs(REAL(3.25), REAL(1e-9)));
}


// =============================================================================
// Polynomial Fitting Tests
// =============================================================================

TEST_CASE("PolynomialFit - Degree 0 (constant)", "[CurveFitting][PolynomialFit]")
{
    Vector<Real> x_data({1.0, 2.0, 3.0, 4.0, 5.0});
    Vector<Real> y_data({10.0, 10.0, 10.0, 10.0, 10.0});
    
    auto result = PolynomialFit(x_data, y_data, 0);
    
    REQUIRE(result.coefficients.size() == 1);
    REQUIRE_THAT(result.coefficients[0], WithinAbs(REAL(10.0), REAL(1e-10)));
}

TEST_CASE("PolynomialFit - Degree 1 (linear)", "[CurveFitting][PolynomialFit]")
{
    // y = 2 + 3x
    Vector<Real> x_data({0.0, 1.0, 2.0, 3.0, 4.0});
    Vector<Real> y_data({2.0, 5.0, 8.0, 11.0, 14.0});
    
    auto result = PolynomialFit(x_data, y_data, 1);
    
    REQUIRE(result.coefficients.size() == 2);
    REQUIRE_THAT(result.coefficients[0], WithinAbs(REAL(2.0), REAL(1e-10)));  // c_0
    REQUIRE_THAT(result.coefficients[1], WithinAbs(REAL(3.0), REAL(1e-10)));  // c_1
}

TEST_CASE("PolynomialFit - Degree 2 (quadratic)", "[CurveFitting][PolynomialFit]")
{
    // y = 1 - 2x + 3x^2
    Vector<Real> x_data({-1.0, 0.0, 1.0, 2.0, 3.0});
    Vector<Real> y_data({6.0, 1.0, 2.0, 9.0, 22.0});
    
    auto result = PolynomialFit(x_data, y_data, 2);
    
    REQUIRE(result.coefficients.size() == 3);
    REQUIRE_THAT(result.coefficients[0], WithinAbs(REAL(1.0), REAL(1e-9)));   // constant
    REQUIRE_THAT(result.coefficients[1], WithinAbs(REAL(-2.0), REAL(1e-9)));  // x
    REQUIRE_THAT(result.coefficients[2], WithinAbs(REAL(3.0), REAL(1e-9)));   // x^2
}

TEST_CASE("PolynomialFit - High degree fitting", "[CurveFitting][PolynomialFit]")
{
    // y = x^4 - 2x^3 + x^2 - x + 1
    Vector<Real> x_data({-1.0, 0.0, 0.5, 1.0, 1.5, 2.0});
    Vector<Real> y_data(6);
    
    for (int i = 0; i < 6; i++) {
        Real x = x_data[i];
        y_data[i] = x*x*x*x - 2*x*x*x + x*x - x + 1;
    }
    
    auto result = PolynomialFit(x_data, y_data, 4);
    
    REQUIRE(result.coefficients.size() == 5);
    REQUIRE_THAT(result.coefficients[0], WithinAbs(REAL(1.0), REAL(1e-8)));   // constant
    REQUIRE_THAT(result.coefficients[1], WithinAbs(REAL(-1.0), REAL(1e-8)));  // x
    REQUIRE_THAT(result.coefficients[2], WithinAbs(REAL(1.0), REAL(1e-8)));   // x^2
    REQUIRE_THAT(result.coefficients[3], WithinAbs(REAL(-2.0), REAL(1e-8)));  // x^3
    REQUIRE_THAT(result.coefficients[4], WithinAbs(REAL(1.0), REAL(1e-8)));   // x^4
}

TEST_CASE("PolynomialFit - Interpolation (degree = n-1)", "[CurveFitting][PolynomialFit]")
{
    // 4 points -> degree 3 polynomial (exact interpolation)
    Vector<Real> x_data({0.0, 1.0, 2.0, 3.0});
    Vector<Real> y_data({1.0, 2.0, 9.0, 28.0});  // Arbitrary points
    
    auto result = PolynomialFit(x_data, y_data, 3);
    
    // Should pass through all points exactly
    REQUIRE_THAT(result.residual_norm, WithinAbs(REAL(0.0), REAL(1e-10)));
    REQUIRE_THAT(result.r_squared, WithinAbs(REAL(1.0), REAL(1e-10)));
}


// =============================================================================
// Evaluate Polynomial Tests
// =============================================================================

TEST_CASE("EvaluatePolynomial - Constant polynomial", "[CurveFitting][EvaluatePolynomial]")
{
    Vector<Real> coeffs({5.0});  // p(x) = 5
    
    REQUIRE_THAT(EvaluatePolynomial(REAL(0.0), coeffs), WithinAbs(REAL(5.0), REAL(1e-10)));
    REQUIRE_THAT(EvaluatePolynomial(REAL(10.0), coeffs), WithinAbs(REAL(5.0), REAL(1e-10)));
    REQUIRE_THAT(EvaluatePolynomial(REAL(-100.0), coeffs), WithinAbs(REAL(5.0), REAL(1e-10)));
}

TEST_CASE("EvaluatePolynomial - Linear polynomial", "[CurveFitting][EvaluatePolynomial]")
{
    Vector<Real> coeffs({2.0, 3.0});  // p(x) = 2 + 3x
    
    REQUIRE_THAT(EvaluatePolynomial(REAL(0.0), coeffs), WithinAbs(REAL(2.0), REAL(1e-10)));
    REQUIRE_THAT(EvaluatePolynomial(REAL(1.0), coeffs), WithinAbs(REAL(5.0), REAL(1e-10)));
    REQUIRE_THAT(EvaluatePolynomial(REAL(2.0), coeffs), WithinAbs(REAL(8.0), REAL(1e-10)));
}

TEST_CASE("EvaluatePolynomial - Quadratic polynomial", "[CurveFitting][EvaluatePolynomial]")
{
    Vector<Real> coeffs({1.0, -2.0, 1.0});  // p(x) = 1 - 2x + x^2 = (x-1)^2
    
    REQUIRE_THAT(EvaluatePolynomial(REAL(0.0), coeffs), WithinAbs(REAL(1.0), REAL(1e-10)));
    REQUIRE_THAT(EvaluatePolynomial(REAL(1.0), coeffs), WithinAbs(REAL(0.0), REAL(1e-10)));  // Root
    REQUIRE_THAT(EvaluatePolynomial(REAL(2.0), coeffs), WithinAbs(REAL(1.0), REAL(1e-10)));
    REQUIRE_THAT(EvaluatePolynomial(REAL(3.0), coeffs), WithinAbs(REAL(4.0), REAL(1e-10)));
}

TEST_CASE("EvaluatePolynomial - Higher order (Horner's method stability)", "[CurveFitting][EvaluatePolynomial]")
{
    // p(x) = 1 + 2x + 3x^2 + 4x^3 + 5x^4
    Vector<Real> coeffs({1.0, 2.0, 3.0, 4.0, 5.0});
    
    Real x = 2.0;
    Real expected = 1 + 2*2 + 3*4 + 4*8 + 5*16;  // = 1 + 4 + 12 + 32 + 80 = 129
    
    REQUIRE_THAT(EvaluatePolynomial(x, coeffs), WithinAbs(expected, REAL(1e-10)));
}

TEST_CASE("EvaluatePolynomial - Empty coefficients", "[CurveFitting][EvaluatePolynomial]")
{
    Vector<Real> coeffs;
    
    REQUIRE_THAT(EvaluatePolynomial(REAL(5.0), coeffs), WithinAbs(REAL(0.0), REAL(1e-10)));
}


// =============================================================================
// Error Handling Tests for GeneralLinearLeastSquares
// =============================================================================

TEST_CASE("GeneralLinearLeastSquares - Error: empty data", "[CurveFitting][GeneralLLS][Errors]")
{
    Vector<Real> x_data;
    Vector<Real> y_data;
    
    Vector<std::function<Real(Real)>> basis({
        [](double) { return 1.0; }
    });
    
    REQUIRE_THROWS_AS(GeneralLinearLeastSquares(x_data, y_data, basis), std::invalid_argument);
}

TEST_CASE("GeneralLinearLeastSquares - Error: empty basis functions", "[CurveFitting][GeneralLLS][Errors]")
{
    Vector<Real> x_data({1.0, 2.0, 3.0});
    Vector<Real> y_data({1.0, 2.0, 3.0});
    
    Vector<std::function<Real(Real)>> basis;  // Empty!
    
    REQUIRE_THROWS_AS(GeneralLinearLeastSquares(x_data, y_data, basis), std::invalid_argument);
}

TEST_CASE("GeneralLinearLeastSquares - Error: size mismatch", "[CurveFitting][GeneralLLS][Errors]")
{
    Vector<Real> x_data({1.0, 2.0, 3.0});
    Vector<Real> y_data({1.0, 2.0});  // Different size
    
    Vector<std::function<Real(Real)>> basis({
        [](double) { return 1.0; }
    });
    
    REQUIRE_THROWS_AS(GeneralLinearLeastSquares(x_data, y_data, basis), std::invalid_argument);
}

TEST_CASE("GeneralLinearLeastSquares - Error: more basis functions than data points", "[CurveFitting][GeneralLLS][Errors]")
{
    Vector<Real> x_data({1.0, 2.0});  // Only 2 points
    Vector<Real> y_data({1.0, 4.0});
    
    Vector<std::function<Real(Real)>> basis({
        [](double) { return 1.0; },
        [](Real x) -> Real { return x; },
        [](Real x) -> Real { return x * x; }  // 3 basis functions!
    });
    
    REQUIRE_THROWS_AS(GeneralLinearLeastSquares(x_data, y_data, basis), std::invalid_argument);
}

TEST_CASE("PolynomialFit - Error: negative degree", "[CurveFitting][PolynomialFit][Errors]")
{
    Vector<Real> x_data({1.0, 2.0, 3.0});
    Vector<Real> y_data({1.0, 2.0, 3.0});
    
    REQUIRE_THROWS_AS(PolynomialFit(x_data, y_data, -1), std::invalid_argument);
}

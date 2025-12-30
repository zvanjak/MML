///////////////////////////////////////////////////////////////////////////////////////////
// fourier_series_tests.cpp - Tests for Fourier Series decomposition
///////////////////////////////////////////////////////////////////////////////////////////

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../../mml/algorithms/Fourier/FourierSeries.h"
#include "../../mml/base/Function.h"
#include <cmath>

using namespace MML;

///////////////////////////////////////////////////////////////////////////////////////////
// Test helpers
///////////////////////////////////////////////////////////////////////////////////////////

// Simple periodic functions for testing
static Real SquareWave(Real x, Real L)
{
    // Square wave with period 2L: +1 for x ∈ (-L, 0), -1 for x ∈ (0, L)
    Real x_mod = std::fmod(x + L, 2.0 * L);
    if (x_mod < 0) x_mod += 2.0 * L;
    return (x_mod < L) ? 1.0 : -1.0;
}

static Real SawtoothWave(Real x, Real L)
{
    // Sawtooth wave with period 2L: linear from -L to L
    Real x_mod = std::fmod(x + L, 2.0 * L);
    if (x_mod < 0) x_mod += 2.0 * L;
    return (x_mod / L) - 1.0;
}

static Real TriangleWave(Real x, Real L)
{
    // Triangle wave with period 2L
    Real x_mod = std::fmod(x + L, 2.0 * L);
    if (x_mod < 0) x_mod += 2.0 * L;
    
    if (x_mod < L)
        return -1.0 + 2.0 * (x_mod / L);
    else
        return 3.0 - 2.0 * (x_mod / L);
}

///////////////////////////////////////////////////////////////////////////////////////////
// Basic Construction and Evaluation Tests
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("FourierSeries - Construction from coefficients", "[fourier][series][construction]")
{
    SECTION("Constant function (N=0)")
    {
        Vector<Real> a{2.0};  // a₀ = 2
        Vector<Real> b{};     // No sine terms
        FourierSeries fs(a, b, Constants::PI);
        
        REQUIRE(fs.NumTerms() == 0);
        REQUIRE(std::abs(fs(0.0) - 1.0) < 1e-10);  // f(x) = a₀/2 = 1
        REQUIRE(std::abs(fs(1.0) - 1.0) < 1e-10);
    }

    SECTION("Single cosine term")
    {
        Vector<Real> a{0.0, 1.0};  // a₀=0, a₁=1
        Vector<Real> b{0.0};       // b₁=0
        FourierSeries fs(a, b, Constants::PI);
        
        REQUIRE(fs.NumTerms() == 1);
        REQUIRE(std::abs(fs(0.0) - 1.0) < 1e-10);     // cos(0) = 1
        REQUIRE(std::abs(fs(Constants::PI) - (-1.0)) < 1e-10);  // cos(π) = -1
    }

    SECTION("Single sine term")
    {
        Vector<Real> a{0.0, 0.0};  // a₀=0, a₁=0
        Vector<Real> b{1.0};       // b₁=1
        FourierSeries fs(a, b, Constants::PI);
        
        REQUIRE(fs.NumTerms() == 1);
        REQUIRE(std::abs(fs(0.0)) < 1e-10);                    // sin(0) = 0
        REQUIRE(std::abs(fs(Constants::PI / 2.0) - 1.0) < 1e-10);  // sin(π/2) = 1
    }

    SECTION("Invalid coefficient sizes throw")
    {
        Vector<Real> a{1.0, 2.0};
        Vector<Real> b{1.0, 2.0};  // Wrong size (should be size 1)
        
        REQUIRE_THROWS_AS(FourierSeries(a, b, Constants::PI), std::invalid_argument);
    }
}

TEST_CASE("FourierSeries - Construction from function", "[fourier][series][construction]")
{
    SECTION("Approximate sin(x) on [-π, π]")
    {
        RealFunctionFromStdFunc sin_func([](Real x) { return std::sin(x); });
        FourierSeries fs(sin_func, Constants::PI, 10);
        
        REQUIRE(fs.NumTerms() == 10);
        
        // sin(x) should have only b₁ = 1, all other coefficients ≈ 0
        const auto& a = fs.CosineCoefficients();
        const auto& b = fs.SineCoefficients();
        
        // All cosine coefficients should be near zero
        for (int i = 0; i <= 10; i++)
        {
            REQUIRE(std::abs(a[i]) < 0.01);
        }
        
        // b₁ should be near 1
        REQUIRE(std::abs(b[0] - 1.0) < 0.01);
        
        // Other sine coefficients should be near zero
        for (int i = 1; i < 10; i++)
        {
            REQUIRE(std::abs(b[i]) < 0.01);
        }
    }

    SECTION("Approximate cos(x) on [-π, π]")
    {
        RealFunctionFromStdFunc cos_func([](Real x) { return std::cos(x); });
        FourierSeries fs(cos_func, Constants::PI, 10);
        
        // cos(x) should have only a₁ = 1, all other coefficients ≈ 0
        const auto& a = fs.CosineCoefficients();
        const auto& b = fs.SineCoefficients();
        
        REQUIRE(std::abs(a[0]) < 0.01);  // a₀ ≈ 0
        REQUIRE(std::abs(a[1] - 1.0) < 0.01);  // a₁ ≈ 1
        
        for (int i = 2; i <= 10; i++)
        {
            REQUIRE(std::abs(a[i]) < 0.01);
        }
        
        for (int i = 0; i < 10; i++)
        {
            REQUIRE(std::abs(b[i]) < 0.01);
        }
    }

    SECTION("Approximate x² on [-1, 1]")
    {
        RealFunctionFromStdFunc x_squared([](Real x) { return x * x; });
        FourierSeries fs(x_squared, 1.0, 10);
        
        // Even function → only cosine terms
        const auto& b = fs.SineCoefficients();
        
        for (int i = 0; i < 10; i++)
        {
            REQUIRE(std::abs(b[i]) < 0.01);
        }
        
        // Check approximation quality
        REQUIRE(std::abs(fs(0.0) - 0.0) < 0.05);
        REQUIRE(std::abs(fs(0.5) - 0.25) < 0.05);
        REQUIRE(std::abs(fs(1.0) - 1.0) < 0.1);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// Derivative and Integral Tests
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("FourierSeries - Derivative", "[fourier][series][calculus]")
{
    SECTION("Derivative of sin(x) is cos(x)")
    {
        // sin(x) = 0·a₀ + 0·cos(x) + 1·sin(x)
        Vector<Real> a{0.0, 0.0};  // a₀=0, a₁=0
        Vector<Real> b{1.0};       // b₁=1
        FourierSeries sin_series(a, b, Constants::PI);
        
        FourierSeries deriv = sin_series.Derivative();
        
        // d/dx[sin(x)] = cos(x) should have a₁ = 1
        const auto& a_deriv = deriv.CosineCoefficients();
        const auto& b_deriv = deriv.SineCoefficients();
        
        REQUIRE(std::abs(a_deriv[0]) < 1e-10);  // a₀ = 0
        REQUIRE(std::abs(a_deriv[1] - 1.0) < 1e-10);  // a₁ = 1 (derivative factor included)
        REQUIRE(std::abs(b_deriv[0]) < 1e-10);  // b₁ ≈ 0
        
        // Verify evaluation
        REQUIRE(std::abs(deriv(0.0) - 1.0) < 0.01);  // cos(0) = 1
        REQUIRE(std::abs(deriv(Constants::PI / 2.0)) < 0.01);  // cos(π/2) = 0
    }

    SECTION("Derivative of cos(x) is -sin(x)")
    {
        // cos(x) = 0·a₀ + 1·cos(x) + 0·sin(x)
        Vector<Real> a{0.0, 1.0};
        Vector<Real> b{0.0};
        FourierSeries cos_series(a, b, Constants::PI);
        
        FourierSeries deriv = cos_series.Derivative();
        
        const auto& a_deriv = deriv.CosineCoefficients();
        const auto& b_deriv = deriv.SineCoefficients();
        
        REQUIRE(std::abs(a_deriv[0]) < 1e-10);
        REQUIRE(std::abs(a_deriv[1]) < 1e-10);
        REQUIRE(std::abs(b_deriv[0] - (-1.0)) < 1e-10);  // b₁ = -1
        
        REQUIRE(std::abs(deriv(0.0)) < 0.01);  // -sin(0) = 0
        REQUIRE(std::abs(deriv(Constants::PI / 2.0) - (-1.0)) < 0.01);  // -sin(π/2) = -1
    }

    SECTION("Derivative of constant is zero")
    {
        Vector<Real> a{4.0};  // a₀ = 4
        Vector<Real> b{};
        FourierSeries const_series(a, b, Constants::PI);
        
        FourierSeries deriv = const_series.Derivative();
        
        REQUIRE(std::abs(deriv(0.0)) < 1e-10);
        REQUIRE(std::abs(deriv(1.0)) < 1e-10);
        REQUIRE(std::abs(deriv(2.0)) < 1e-10);
    }
}

TEST_CASE("FourierSeries - Integral", "[fourier][series][calculus]")
{
    SECTION("Integral of cos(x) is sin(x)")
    {
        Vector<Real> a{0.0, 1.0};
        Vector<Real> b{0.0};
        FourierSeries cos_series(a, b, Constants::PI);
        
        FourierSeries integral = cos_series.Integral();
        
        const auto& b_int = integral.SineCoefficients();
        
        // ∫cos(x)dx = sin(x) should have b₁ = 1
        REQUIRE(std::abs(b_int[0] - 1.0) < 1e-10);
    }

    SECTION("Integral of sin(x) is -cos(x)")
    {
        Vector<Real> a{0.0, 0.0};
        Vector<Real> b{1.0};
        FourierSeries sin_series(a, b, Constants::PI);
        
        FourierSeries integral = sin_series.Integral();
        
        const auto& a_int = integral.CosineCoefficients();
        
        // ∫sin(x)dx = -cos(x) should have a₁ = -1
        REQUIRE(std::abs(a_int[1] - (-1.0)) < 1e-10);
    }

    SECTION("Second derivative returns to original (with scaling)")
    {
        Vector<Real> a{0.0, 0.0, 1.0};  // cos(2πx/L)
        Vector<Real> b{0.0, 0.0};
        FourierSeries series(a, b, 1.0);
        
        FourierSeries deriv1 = series.Derivative();
        FourierSeries deriv2 = deriv1.Derivative();
        
        // d²/dx²[cos(2πx)] = -4π²cos(2πx)
        Real factor = -4.0 * Constants::PI * Constants::PI;
        
        // Check that second derivative has correct magnitude
        REQUIRE(std::abs(deriv2(0.0) - factor * 1.0) < 0.1);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// Complex Coefficients Tests
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("FourierSeries - Complex coefficients", "[fourier][series][complex]")
{
    SECTION("Convert real to complex and back")
    {
        Vector<Real> a{1.0, 2.0, 3.0};  // a₀, a₁, a₂
        Vector<Real> b{4.0, 5.0};       // b₁, b₂
        FourierSeries fs(a, b, Constants::PI);
        
        Vector<Complex> c = fs.ComplexCoefficients();
        
        REQUIRE(c.size() == 5);  // c₋₂, c₋₁, c₀, c₁, c₂
        
        // c₀ = a₀/2
        REQUIRE(std::abs(c[2].real() - 0.5) < 1e-10);
        REQUIRE(std::abs(c[2].imag()) < 1e-10);
        
        // Reconstruct
        FourierSeries fs2 = FourierSeries::FromComplexCoefficients(c, Constants::PI);
        
        const auto& a2 = fs2.CosineCoefficients();
        const auto& b2 = fs2.SineCoefficients();
        
        for (int i = 0; i < 3; i++)
        {
            REQUIRE(std::abs(a[i] - a2[i]) < 1e-10);
        }
        
        for (int i = 0; i < 2; i++)
        {
            REQUIRE(std::abs(b[i] - b2[i]) < 1e-10);
        }
    }

    SECTION("Complex coefficients satisfy conjugate symmetry")
    {
        Vector<Real> a{0.0, 1.0, 2.0};
        Vector<Real> b{3.0, 4.0};
        FourierSeries fs(a, b, Constants::PI);
        
        Vector<Complex> c = fs.ComplexCoefficients();
        
        // c₋ₙ = c̄ₙ (complex conjugate)
        REQUIRE(std::abs(c[0] - std::conj(c[4])) < 1e-10);  // c₋₂ = c̄₂
        REQUIRE(std::abs(c[1] - std::conj(c[3])) < 1e-10);  // c₋₁ = c̄₁
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// Energy and Norms (Parseval's Theorem)
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("FourierSeries - Energy and Parseval's theorem", "[fourier][series][energy]")
{
    SECTION("Energy of constant function")
    {
        Vector<Real> a{2.0};  // a₀ = 2 → f(x) = 1
        Vector<Real> b{};
        FourierSeries fs(a, b, 1.0);  // L = 1
        
        // Energy = ∫[-1,1] 1² dx = 2
        Real energy = fs.Energy();
        REQUIRE(std::abs(energy - 2.0) < 0.01);
    }

    SECTION("Energy of sin(πx)")
    {
        Vector<Real> a{0.0, 0.0};
        Vector<Real> b{1.0};  // sin(πx) with L=1
        FourierSeries fs(a, b, 1.0);
        
        // Energy = ∫[-1,1] sin²(πx) dx = 1
        Real energy = fs.Energy();
        REQUIRE(std::abs(energy - 1.0) < 0.01);
    }

    SECTION("L2 norm consistency")
    {
        Vector<Real> a{2.0, 1.0};
        Vector<Real> b{1.0};
        FourierSeries fs(a, b, Constants::PI);
        
        Real norm = fs.L2Norm();
        Real energy = fs.Energy();
        
        // ||f||₂ = √(Energy/2L)
        Real expected_norm = std::sqrt(energy / (2.0 * Constants::PI));
        REQUIRE(std::abs(norm - expected_norm) < 1e-10);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// Series Operations Tests
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("FourierSeries - Arithmetic operations", "[fourier][series][operations]")
{
    SECTION("Addition")
    {
        Vector<Real> a1{1.0, 2.0};
        Vector<Real> b1{3.0};
        FourierSeries fs1(a1, b1, Constants::PI);
        
        Vector<Real> a2{4.0, 5.0};
        Vector<Real> b2{6.0};
        FourierSeries fs2(a2, b2, Constants::PI);
        
        FourierSeries sum = fs1 + fs2;
        
        const auto& a_sum = sum.CosineCoefficients();
        const auto& b_sum = sum.SineCoefficients();
        
        REQUIRE(std::abs(a_sum[0] - 5.0) < 1e-10);  // 1+4
        REQUIRE(std::abs(a_sum[1] - 7.0) < 1e-10);  // 2+5
        REQUIRE(std::abs(b_sum[0] - 9.0) < 1e-10);  // 3+6
        
        // Verify evaluation
        Real x = 0.5;
        REQUIRE(std::abs(sum(x) - (fs1(x) + fs2(x))) < 1e-10);
    }

    SECTION("Subtraction")
    {
        Vector<Real> a1{10.0, 20.0};
        Vector<Real> b1{30.0};
        FourierSeries fs1(a1, b1, Constants::PI);
        
        Vector<Real> a2{4.0, 5.0};
        Vector<Real> b2{6.0};
        FourierSeries fs2(a2, b2, Constants::PI);
        
        FourierSeries diff = fs1 - fs2;
        
        const auto& a_diff = diff.CosineCoefficients();
        const auto& b_diff = diff.SineCoefficients();
        
        REQUIRE(std::abs(a_diff[0] - 6.0) < 1e-10);   // 10-4
        REQUIRE(std::abs(a_diff[1] - 15.0) < 1e-10);  // 20-5
        REQUIRE(std::abs(b_diff[0] - 24.0) < 1e-10);  // 30-6
    }

    SECTION("Scalar multiplication")
    {
        Vector<Real> a{1.0, 2.0, 3.0};
        Vector<Real> b{4.0, 5.0};
        FourierSeries fs(a, b, Constants::PI);
        
        FourierSeries scaled = fs * 3.0;
        
        const auto& a_scaled = scaled.CosineCoefficients();
        const auto& b_scaled = scaled.SineCoefficients();
        
        REQUIRE(std::abs(a_scaled[0] - 3.0) < 1e-10);
        REQUIRE(std::abs(a_scaled[1] - 6.0) < 1e-10);
        REQUIRE(std::abs(a_scaled[2] - 9.0) < 1e-10);
        REQUIRE(std::abs(b_scaled[0] - 12.0) < 1e-10);
        REQUIRE(std::abs(b_scaled[1] - 15.0) < 1e-10);
        
        // Check evaluation
        Real x = 1.0;
        REQUIRE(std::abs(scaled(x) - 3.0 * fs(x)) < 1e-10);
    }

    SECTION("Truncation")
    {
        Vector<Real> a{1.0, 2.0, 3.0, 4.0};  // N=3
        Vector<Real> b{5.0, 6.0, 7.0};
        FourierSeries fs(a, b, Constants::PI);
        
        FourierSeries truncated = fs.Truncate(2);  // Keep only first 2 terms
        
        REQUIRE(truncated.NumTerms() == 2);
        
        const auto& a_trunc = truncated.CosineCoefficients();
        REQUIRE(a_trunc.size() == 3);  // a₀, a₁, a₂
        REQUIRE(std::abs(a_trunc[0] - 1.0) < 1e-10);
        REQUIRE(std::abs(a_trunc[1] - 2.0) < 1e-10);
        REQUIRE(std::abs(a_trunc[2] - 3.0) < 1e-10);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// Periodic Function Approximation Tests
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("FourierSeries - Square wave approximation", "[fourier][series][periodic]")
{
    SECTION("Converges to square wave")
    {
        Real L = Constants::PI;
        RealFunctionFromStdFunc square([L](Real x) { return SquareWave(x, L); });
        
        // Square wave Fourier series: f(x) = (4/π)Σ[sin((2k-1)x)/(2k-1)]
        // Uses only odd harmonics
        FourierSeries fs(square, L, 50);
        
        const auto& a = fs.CosineCoefficients();
        const auto& b = fs.SineCoefficients();
        
        // All cosine coefficients should be near zero (odd function)
        for (int i = 0; i <= 50; i++)
        {
            REQUIRE(std::abs(a[i]) < 0.1);
        }
        
        // Odd sine coefficients should dominate
        REQUIRE(std::abs(b[0]) > 1.0);  // b₁ is large
        
        // Check approximation at center
        REQUIRE(std::abs(fs(-L/2.0) - 1.0) < 0.3);  // Gibbs phenomenon
        REQUIRE(std::abs(fs(L/2.0) - (-1.0)) < 0.3);
    }
}

TEST_CASE("FourierSeries - Sawtooth wave approximation", "[fourier][series][periodic]")
{
    SECTION("Approximate sawtooth")
    {
        Real L = 1.0;
        RealFunctionFromStdFunc sawtooth([L](Real x) { return SawtoothWave(x, L); });
        
        FourierSeries fs(sawtooth, L, 30);
        
        // Sawtooth is odd → only sine terms
        const auto& a = fs.CosineCoefficients();
        
        for (int i = 0; i <= 30; i++)
        {
            REQUIRE(std::abs(a[i]) < 0.1);
        }
        
        // Check some approximation points
        REQUIRE(std::abs(fs(0.0)) < 0.2);
        REQUIRE(std::abs(fs(-0.5) - (-0.5)) < 0.2);
        REQUIRE(std::abs(fs(0.5) - 0.5) < 0.2);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// Convergence Analysis Tests
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("FourierSeries - Convergence analysis", "[fourier][series][convergence]")
{
    SECTION("Smooth function has fast decay")
    {
        // sin(x) is smooth → exponential coefficient decay
        RealFunctionFromStdFunc sin_func([](Real x) { return std::sin(x); });
        FourierSeries fs(sin_func, Constants::PI, 20);
        
        // Check coefficient magnitudes decay
        Real mag1 = fs.CoefficientMagnitude(1);
        Real mag5 = fs.CoefficientMagnitude(5);
        Real mag10 = fs.CoefficientMagnitude(10);
        
        // Should have: mag1 >> mag5 >> mag10 for smooth functions
        REQUIRE(mag1 > 10.0 * mag5);
        REQUIRE(mag5 > 5.0 * mag10);
    }

    SECTION("MaxCoefficientMagnitude")
    {
        Vector<Real> a{0.5, 1.0, 0.3, 2.5};
        Vector<Real> b{0.2, 0.8, 0.1};
        FourierSeries fs(a, b, Constants::PI);
        
        Real max_mag = fs.MaxCoefficientMagnitude();
        
        // Should be approximately |a₃| = 2.5 (largest coefficient)
        REQUIRE(max_mag > 2.4);
        REQUIRE(max_mag < 2.6);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// Edge Cases and Error Handling
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("FourierSeries - Edge cases", "[fourier][series][edge]")
{
    SECTION("Zero series evaluates to zero")
    {
        Vector<Real> a{0.0};
        Vector<Real> b{};
        FourierSeries fs(a, b, Constants::PI);
        
        REQUIRE(std::abs(fs(0.0)) < 1e-10);
        REQUIRE(std::abs(fs(1.0)) < 1e-10);
        REQUIRE(std::abs(fs(Constants::PI)) < 1e-10);
    }

    SECTION("Different periods throw on operations")
    {
        Vector<Real> a1{1.0, 2.0};
        Vector<Real> b1{3.0};
        FourierSeries fs1(a1, b1, 1.0);
        
        Vector<Real> a2{4.0, 5.0};
        Vector<Real> b2{6.0};
        FourierSeries fs2(a2, b2, 2.0);  // Different period!
        
        REQUIRE_THROWS_AS(fs1 + fs2, std::invalid_argument);
        REQUIRE_THROWS_AS(fs1 - fs2, std::invalid_argument);
    }

    SECTION("Negative L throws")
    {
        RealFunctionFromStdFunc func([](Real x) { return x; });
        REQUIRE_THROWS_AS(FourierSeries(func, -1.0, 10), std::invalid_argument);
    }

    SECTION("Negative N throws")
    {
        RealFunctionFromStdFunc func([](Real x) { return x; });
        REQUIRE_THROWS_AS(FourierSeries(func, 1.0, -5), std::invalid_argument);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// FourierBasis Unit Tests
// Tests for IOrthogonalBasis interface and FourierBasis implementation
///////////////////////////////////////////////////////////////////////////////////////////

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../../mml/algorithms/Fourier/FourierBasis.h"
#include "../../mml/interfaces/IFunction.h"
#include <cmath>

using namespace MML;

// Use anonymous namespace to avoid conflicts with other test files
namespace {

///////////////////////////////////////////////////////////////////////////////////////////
// Test Tags:
// [fourier][basis] - All FourierBasis tests
// [orthogonal]     - Orthogonality verification
// [expansion]      - Function expansion tests
// [complex]        - Complex exponential basis tests
///////////////////////////////////////////////////////////////////////////////////////////

// Helper function classes (local to this file)
class SinFunc : public IRealFunction {
    Real _k;  // Frequency multiplier
public:
    explicit SinFunc(Real k = 1.0) : _k(k) {}
    Real operator()(Real x) const override { return std::sin(_k * x); }
};

class CosFunc : public IRealFunction {
    Real _k;
public:
    explicit CosFunc(Real k = 1.0) : _k(k) {}
    Real operator()(Real x) const override { return std::cos(_k * x); }
};

class ConstFunc : public IRealFunction {
    Real _c;
public:
    explicit ConstFunc(Real c = 1.0) : _c(c) {}
    Real operator()(Real x) const override { return _c; }
};

class SquareWaveFunc : public IRealFunction {
public:
    Real operator()(Real x) const override {
        // Square wave: +1 for x > 0, -1 for x < 0
        return (x >= 0) ? 1.0 : -1.0;
    }
};

} // end anonymous namespace

///////////////////////////////////////////////////////////////////////////////////////////
// Basic Construction Tests
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("FourierBasis construction", "[fourier][basis]")
{
    SECTION("Default construction (L = π)")
    {
        FourierBasis basis;
        REQUIRE(std::abs(basis.HalfPeriod() - Constants::PI) < 1e-10);
        REQUIRE(std::abs(basis.Period() - 2.0 * Constants::PI) < 1e-10);
        REQUIRE(std::abs(basis.DomainMin() - (-Constants::PI)) < 1e-10);
        REQUIRE(std::abs(basis.DomainMax() - Constants::PI) < 1e-10);
    }

    SECTION("Custom period")
    {
        FourierBasis basis(2.0);
        REQUIRE(std::abs(basis.HalfPeriod() - 2.0) < 1e-10);
        REQUIRE(std::abs(basis.Period() - 4.0) < 1e-10);
        REQUIRE(std::abs(basis.DomainMin() - (-2.0)) < 1e-10);
        REQUIRE(std::abs(basis.DomainMax() - 2.0) < 1e-10);
    }

    SECTION("Invalid period throws")
    {
        REQUIRE_THROWS_AS(FourierBasis(0.0), std::invalid_argument);
        REQUIRE_THROWS_AS(FourierBasis(-1.0), std::invalid_argument);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// Basis Function Evaluation Tests
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("FourierBasis evaluation", "[fourier][basis]")
{
    FourierBasis basis(Constants::PI);
    
    SECTION("Constant basis function (n=0)")
    {
        // φ₀(x) = 1
        REQUIRE(std::abs(basis.Evaluate(0, 0.0) - 1.0) < 1e-10);
        REQUIRE(std::abs(basis.Evaluate(0, 1.0) - 1.0) < 1e-10);
        REQUIRE(std::abs(basis.Evaluate(0, -1.0) - 1.0) < 1e-10);
    }

    SECTION("First cosine basis (n=1)")
    {
        // φ₁(x) = cos(πx/L) = cos(x) when L = π
        REQUIRE(std::abs(basis.Evaluate(1, 0.0) - 1.0) < 1e-10);               // cos(0)
        REQUIRE(std::abs(basis.Evaluate(1, Constants::PI/2) - 0.0) < 1e-10);   // cos(π/2)
        REQUIRE(std::abs(basis.Evaluate(1, Constants::PI) - (-1.0)) < 1e-10);  // cos(π)
    }

    SECTION("First sine basis (n=2)")
    {
        // φ₂(x) = sin(πx/L) = sin(x) when L = π
        REQUIRE(std::abs(basis.Evaluate(2, 0.0) - 0.0) < 1e-10);               // sin(0)
        REQUIRE(std::abs(basis.Evaluate(2, Constants::PI/2) - 1.0) < 1e-10);   // sin(π/2)
        REQUIRE(std::abs(basis.Evaluate(2, Constants::PI) - 0.0) < 1e-10);     // sin(π)
    }

    SECTION("Higher harmonics")
    {
        // φ₃(x) = cos(2x), φ₄(x) = sin(2x), etc.
        Real x = Constants::PI / 4;
        
        // cos(2x) at x = π/4 → cos(π/2) = 0
        REQUIRE(std::abs(basis.Evaluate(3, x) - 0.0) < 1e-10);
        
        // sin(2x) at x = π/4 → sin(π/2) = 1
        REQUIRE(std::abs(basis.Evaluate(4, x) - 1.0) < 1e-10);
    }

    SECTION("Negative index throws")
    {
        REQUIRE_THROWS_AS(basis.Evaluate(-1, 0.0), std::invalid_argument);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// Index Conversion Tests
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("FourierBasis index conversions", "[fourier][basis]")
{
    SECTION("Harmonic number")
    {
        REQUIRE(FourierBasis::HarmonicNumber(0) == 0);
        REQUIRE(FourierBasis::HarmonicNumber(1) == 1);  // cos(1·πx/L)
        REQUIRE(FourierBasis::HarmonicNumber(2) == 1);  // sin(1·πx/L)
        REQUIRE(FourierBasis::HarmonicNumber(3) == 2);  // cos(2·πx/L)
        REQUIRE(FourierBasis::HarmonicNumber(4) == 2);  // sin(2·πx/L)
        REQUIRE(FourierBasis::HarmonicNumber(5) == 3);
        REQUIRE(FourierBasis::HarmonicNumber(6) == 3);
    }

    SECTION("Cosine/Sine detection")
    {
        REQUIRE(FourierBasis::IsCosine(0) == true);   // constant (technically a₀)
        REQUIRE(FourierBasis::IsCosine(1) == true);   // a₁
        REQUIRE(FourierBasis::IsCosine(2) == false);  // b₁
        REQUIRE(FourierBasis::IsCosine(3) == true);   // a₂
        REQUIRE(FourierBasis::IsCosine(4) == false);  // b₂
        
        REQUIRE(FourierBasis::IsSine(0) == false);
        REQUIRE(FourierBasis::IsSine(1) == false);
        REQUIRE(FourierBasis::IsSine(2) == true);
        REQUIRE(FourierBasis::IsSine(3) == false);
        REQUIRE(FourierBasis::IsSine(4) == true);
    }

    SECTION("Traditional index conversion")
    {
        auto [type0, k0] = FourierBasis::ToTraditionalIndex(0);
        REQUIRE(type0 == 'a');
        REQUIRE(k0 == 0);
        
        auto [type1, k1] = FourierBasis::ToTraditionalIndex(1);
        REQUIRE(type1 == 'a');
        REQUIRE(k1 == 1);
        
        auto [type2, k2] = FourierBasis::ToTraditionalIndex(2);
        REQUIRE(type2 == 'b');
        REQUIRE(k2 == 1);
        
        auto [type3, k3] = FourierBasis::ToTraditionalIndex(3);
        REQUIRE(type3 == 'a');
        REQUIRE(k3 == 2);
    }

    SECTION("Round-trip index conversion")
    {
        for (int n = 0; n < 10; n++)
        {
            auto [type, k] = FourierBasis::ToTraditionalIndex(n);
            int recovered = FourierBasis::FromTraditionalIndex(type, k);
            REQUIRE(recovered == n);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// Normalization Tests
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("FourierBasis normalization", "[fourier][basis]")
{
    FourierBasis basis(Constants::PI);
    
    SECTION("Constant basis normalization")
    {
        // ||φ₀||² = ∫[-π,π] 1² dx = 2π
        REQUIRE(std::abs(basis.Normalization(0) - 2.0 * Constants::PI) < 1e-10);
    }

    SECTION("Trigonometric basis normalization")
    {
        // ||φₙ||² = ∫[-π,π] cos²(kx) dx = π (for n > 0)
        // Same for sin²
        for (int n = 1; n < 10; n++)
        {
            REQUIRE(std::abs(basis.Normalization(n) - Constants::PI) < 1e-10);
        }
    }

    SECTION("Weight function is uniform")
    {
        REQUIRE(std::abs(basis.WeightFunction(0.0) - 1.0) < 1e-10);
        REQUIRE(std::abs(basis.WeightFunction(1.0) - 1.0) < 1e-10);
        REQUIRE(std::abs(basis.WeightFunction(-1.0) - 1.0) < 1e-10);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// Orthogonality Tests
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("FourierBasis orthogonality", "[fourier][basis][orthogonal]")
{
    FourierBasis basis(Constants::PI);
    
    SECTION("Numerical orthogonality verification")
    {
        // Verify ⟨φₘ, φₙ⟩ ≈ 0 for m ≠ n
        // ⟨φₘ, φₙ⟩ = ∫[-π,π] φₘ(x)φₙ(x) dx
        
        auto inner_product = [&basis](int m, int n) {
            class Integrand : public IRealFunction {
                const FourierBasis& _b;
                int _m, _n;
            public:
                Integrand(const FourierBasis& b, int m, int n) : _b(b), _m(m), _n(n) {}
                Real operator()(Real x) const override {
                    return _b.Evaluate(_m, x) * _b.Evaluate(_n, x);
                }
            };
            
            Integrand integrand(basis, m, n);
            return IntegrateTrap(integrand, basis.DomainMin(), basis.DomainMax(), 1e-10).value;
        };
        
        // Test orthogonality for first several basis functions
        for (int m = 0; m < 5; m++)
        {
            for (int n = m + 1; n < 5; n++)
            {
                Real ip = inner_product(m, n);
                REQUIRE(std::abs(ip) < 1e-8);
            }
        }
    }

    SECTION("Normalization verification")
    {
        // Verify ⟨φₙ, φₙ⟩ = ||φₙ||²
        auto norm_sq = [&basis](int n) {
            class Integrand : public IRealFunction {
                const FourierBasis& _b;
                int _n;
            public:
                Integrand(const FourierBasis& b, int n) : _b(b), _n(n) {}
                Real operator()(Real x) const override {
                    Real val = _b.Evaluate(_n, x);
                    return val * val;
                }
            };
            
            Integrand integrand(basis, n);
            return IntegrateTrap(integrand, basis.DomainMin(), basis.DomainMax(), 1e-10).value;
        };
        
        for (int n = 0; n < 5; n++)
        {
            Real computed = norm_sq(n);
            Real expected = basis.Normalization(n);
            REQUIRE(std::abs(computed - expected) < 1e-6);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// Coefficient Computation Tests
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("FourierBasis coefficient computation", "[fourier][basis][expansion]")
{
    FourierBasis basis(Constants::PI);
    
    SECTION("Constant function integration verification")
    {
        // Verify constant function integrates correctly
        ConstFunc f(3.0);
        
        // Verify the function returns 3.0
        REQUIRE(std::abs(f(0.0) - 3.0) < 1e-10);
        REQUIRE(std::abs(f(1.0) - 3.0) < 1e-10);
        REQUIRE(std::abs(f(-1.0) - 3.0) < 1e-10);
        
        // Direct integration of f(x) over [-π, π]
        auto integral_result = IntegrateTrap(f, -Constants::PI, Constants::PI, 1e-10);
        Real expected_integral = 6.0 * Constants::PI;  // ∫3 dx from -π to π = 6π
        REQUIRE(std::abs(integral_result.value - expected_integral) < 1e-6);
    }
    
    SECTION("Constant function expansion")
    {
        ConstFunc f(3.0);  // f(x) = 3
        
        // c₀ = (1/||φ₀||²) ∫ f(x)·φ₀(x) dx = (1/2π) ∫ 3·1 dx = (1/2π)·6π = 3
        Real c0 = basis.ComputeCoefficient(f, 0);
        REQUIRE(std::abs(c0 - 3.0) < 1e-6);
        
        // All other coefficients should be zero
        for (int n = 1; n < 5; n++)
        {
            Real cn = basis.ComputeCoefficient(f, n);
            REQUIRE(std::abs(cn) < 1e-6);
        }
    }

    SECTION("Cosine function expansion")
    {
        CosFunc f(1.0);  // f(x) = cos(x) = cos(πx/π) = φ₁
        
        // c₁ should be 1, all others should be 0
        for (int n = 0; n < 5; n++)
        {
            Real cn = basis.ComputeCoefficient(f, n);
            Real expected = (n == 1) ? 1.0 : 0.0;
            REQUIRE(std::abs(cn - expected) < 1e-6);
        }
    }

    SECTION("Sine function expansion")
    {
        SinFunc f(1.0);  // f(x) = sin(x) = sin(πx/π) = φ₂
        
        // c₂ should be 1, all others should be 0
        for (int n = 0; n < 5; n++)
        {
            Real cn = basis.ComputeCoefficient(f, n);
            Real expected = (n == 2) ? 1.0 : 0.0;
            REQUIRE(std::abs(cn - expected) < 1e-6);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// Function Expansion Tests
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("FourierBasis Expand method", "[fourier][basis][expansion]")
{
    FourierBasis basis(Constants::PI);
    
    SECTION("Expand constant function")
    {
        ConstFunc f(2.0);
        Vector<Real> coeffs = basis.Expand(f, 3);  // 2*3+1 = 7 terms
        
        REQUIRE(coeffs.size() == 7);
        REQUIRE(std::abs(coeffs[0] - 2.0) < 1e-6);  // c₀ = 2
        for (int n = 1; n < 7; n++)
        {
            REQUIRE(std::abs(coeffs[n]) < 1e-6);  // All others zero
        }
    }

    SECTION("Expand pure harmonic")
    {
        SinFunc f(2.0);  // f(x) = sin(2x) = φ₄
        Vector<Real> coeffs = basis.Expand(f, 3);
        
        // Only c₄ should be nonzero
        REQUIRE(std::abs(coeffs[4] - 1.0) < 1e-6);
        for (int n = 0; n < 7; n++)
        {
            if (n != 4)
                REQUIRE(std::abs(coeffs[n]) < 1e-6);
        }
    }

    SECTION("Evaluate expansion")
    {
        CosFunc f(1.0);
        Vector<Real> coeffs = basis.Expand(f, 5);
        
        // Expansion should reproduce original function
        for (Real x = -Constants::PI; x <= Constants::PI; x += 0.1)
        {
            Real original = f(x);
            Real expanded = basis.EvaluateExpansion(coeffs, x);
            REQUIRE(std::abs(original - expanded) < 1e-6);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// FourierSeries Conversion Tests
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("FourierBasis to FourierSeries conversion", "[fourier][basis][expansion]")
{
    FourierBasis basis(Constants::PI);
    
    SECTION("Convert pure cosine")
    {
        CosFunc f(1.0);
        Vector<Real> coeffs = basis.Expand(f, 3);
        FourierSeries fs = basis.ToFourierSeries(coeffs);
        
        // Verify FourierSeries reproduces function
        for (Real x = -Constants::PI; x <= Constants::PI; x += 0.2)
        {
            REQUIRE(std::abs(fs(x) - std::cos(x)) < 1e-6);
        }
    }

    SECTION("Convert mixed function")
    {
        // f(x) = 2 + 3cos(x) + 4sin(2x)
        class MixedFunction : public IRealFunction {
        public:
            Real operator()(Real x) const override {
                return 2.0 + 3.0 * std::cos(x) + 4.0 * std::sin(2.0 * x);
            }
        };
        
        MixedFunction f;
        Vector<Real> coeffs = basis.Expand(f, 5);
        FourierSeries fs = basis.ToFourierSeries(coeffs);
        
        for (Real x = -Constants::PI; x <= Constants::PI; x += 0.2)
        {
            REQUIRE(std::abs(fs(x) - f(x)) < 1e-5);
        }
    }

    SECTION("Round-trip conversion")
    {
        CosFunc f(2.0);  // cos(2x)
        Vector<Real> coeffs1 = basis.Expand(f, 5);
        FourierSeries fs = basis.ToFourierSeries(coeffs1);
        Vector<Real> coeffs2 = basis.FromFourierSeries(fs);
        
        REQUIRE(coeffs1.size() == coeffs2.size());
        for (size_t i = 0; i < coeffs1.size(); i++)
        {
            REQUIRE(std::abs(coeffs1[i] - coeffs2[i]) < 1e-10);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// Complex Fourier Basis Tests
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("ComplexFourierBasis construction", "[fourier][basis][complex]")
{
    SECTION("Default construction")
    {
        ComplexFourierBasis basis;
        REQUIRE(std::abs(basis.HalfPeriod() - Constants::PI) < 1e-10);
    }

    SECTION("Invalid period throws")
    {
        REQUIRE_THROWS_AS(ComplexFourierBasis(0.0), std::invalid_argument);
        REQUIRE_THROWS_AS(ComplexFourierBasis(-1.0), std::invalid_argument);
    }
}

TEST_CASE("ComplexFourierBasis evaluation", "[fourier][basis][complex]")
{
    ComplexFourierBasis basis(Constants::PI);
    
    SECTION("Basis function at zero")
    {
        // φₙ(0) = e^0 = 1 for all n
        for (int n = -5; n <= 5; n++)
        {
            Complex val = basis.Evaluate(n, 0.0);
            REQUIRE(std::abs(val.real() - 1.0) < 1e-10);
            REQUIRE(std::abs(val.imag()) < 1e-10);
        }
    }

    SECTION("Positive and negative frequencies")
    {
        Real x = Constants::PI / 4;
        
        // φ₁(x) = e^(ix) and φ₋₁(x) = e^(-ix) should be conjugates
        Complex pos = basis.Evaluate(1, x);
        Complex neg = basis.Evaluate(-1, x);
        
        REQUIRE(std::abs(pos.real() - neg.real()) < 1e-10);
        REQUIRE(std::abs(pos.imag() + neg.imag()) < 1e-10);  // Conjugate
    }
}

TEST_CASE("ComplexFourierBasis expansion", "[fourier][basis][complex][expansion]")
{
    ComplexFourierBasis basis(Constants::PI);
    
    SECTION("Expand real constant")
    {
        ConstFunc f(5.0);
        Vector<Complex> coeffs = basis.Expand(f, 3);  // n = -3 to 3
        
        REQUIRE(coeffs.size() == 7);
        
        // c₀ should be 5.0, all others ~0
        int center = 3;  // Index of c₀
        REQUIRE(std::abs(coeffs[center].real() - 5.0) < 1e-5);
        REQUIRE(std::abs(coeffs[center].imag()) < 1e-5);
        
        for (size_t i = 0; i < coeffs.size(); i++)
        {
            if (i != center)
            {
                REQUIRE(std::abs(coeffs[i]) < 1e-5);
            }
        }
    }

    SECTION("Expand cosine (real even function)")
    {
        CosFunc f(1.0);  // cos(x) = (e^(ix) + e^(-ix))/2
        Vector<Complex> coeffs = basis.Expand(f, 3);
        
        // For cos(x): c₁ = c₋₁ = 0.5, others = 0
        int center = 3;
        
        // c₁ (index center+1) ≈ 0.5
        REQUIRE(std::abs(coeffs[center + 1].real() - 0.5) < 1e-5);
        REQUIRE(std::abs(coeffs[center + 1].imag()) < 1e-5);
        
        // c₋₁ (index center-1) ≈ 0.5
        REQUIRE(std::abs(coeffs[center - 1].real() - 0.5) < 1e-5);
        REQUIRE(std::abs(coeffs[center - 1].imag()) < 1e-5);
    }

    SECTION("Expand sine (real odd function)")
    {
        SinFunc f(1.0);  // sin(x) = (e^(ix) - e^(-ix))/(2i)
        Vector<Complex> coeffs = basis.Expand(f, 3);
        
        // For sin(x): c₁ = -i/2 = (0, -0.5), c₋₁ = i/2 = (0, 0.5)
        int center = 3;
        
        // c₁ ≈ -i/2
        REQUIRE(std::abs(coeffs[center + 1].real()) < 1e-5);
        REQUIRE(std::abs(coeffs[center + 1].imag() - (-0.5)) < 1e-5);
        
        // c₋₁ ≈ i/2
        REQUIRE(std::abs(coeffs[center - 1].real()) < 1e-5);
        REQUIRE(std::abs(coeffs[center - 1].imag() - 0.5) < 1e-5);
    }

    SECTION("Evaluate expansion reproduces function")
    {
        // f(x) = 2 + 3cos(x) + sin(2x)
        class TestFunction : public IRealFunction {
        public:
            Real operator()(Real x) const override {
                return 2.0 + 3.0 * std::cos(x) + std::sin(2.0 * x);
            }
        };
        
        TestFunction f;
        Vector<Complex> coeffs = basis.Expand(f, 5);
        
        for (Real x = -Constants::PI; x <= Constants::PI; x += 0.3)
        {
            Complex expanded = basis.EvaluateExpansion(coeffs, x);
            REQUIRE(std::abs(expanded.real() - f(x)) < 1e-4);
            REQUIRE(std::abs(expanded.imag()) < 1e-4);  // Should be real
        }
    }
}

TEST_CASE("ComplexFourierBasis to FourierSeries", "[fourier][basis][complex]")
{
    ComplexFourierBasis cbasis(Constants::PI);
    
    SECTION("Convert and verify")
    {
        CosFunc f(2.0);  // cos(2x)
        Vector<Complex> coeffs = cbasis.Expand(f, 5);
        FourierSeries fs = cbasis.ToFourierSeries(coeffs);
        
        // FourierSeries should reproduce cos(2x)
        for (Real x = -Constants::PI; x <= Constants::PI; x += 0.2)
        {
            REQUIRE(std::abs(fs(x) - std::cos(2.0 * x)) < 1e-4);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// Square Wave Approximation (Classic Fourier Series Test)
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("FourierBasis square wave approximation", "[fourier][basis][expansion]")
{
    FourierBasis basis(Constants::PI);
    SquareWaveFunc sw;
    
    SECTION("Coefficient structure")
    {
        Vector<Real> coeffs = basis.Expand(sw, 10);
        
        // Square wave: only odd sine harmonics are nonzero
        // bₙ = 4/(nπ) for n = 1, 3, 5, ...
        
        // c₀ (constant) should be ~0 (symmetric around 0)
        // Relax tolerance due to numerical integration of discontinuous function
        REQUIRE(std::abs(coeffs[0]) < 1e-5);
        
        // Odd cosines should be ~0 (odd function)
        // Relax tolerance due to discontinuous function integration
        for (int k = 1; k <= 10; k++)
        {
            int cos_idx = 2 * k - 1;  // Cosine indices
            REQUIRE(std::abs(coeffs[cos_idx]) < 1e-5);
        }
        
        // Check first few sine coefficients
        // b₁ = 4/π ≈ 1.273
        Real b1 = coeffs[2];  // Index 2 is sin(x)
        REQUIRE(std::abs(b1 - 4.0/Constants::PI) < 0.01);
        
        // b₂ should be ~0 (even harmonic)
        Real b2 = coeffs[4];
        REQUIRE(std::abs(b2) < 1e-6);
        
        // b₃ = 4/(3π) ≈ 0.424
        Real b3 = coeffs[6];
        REQUIRE(std::abs(b3 - 4.0/(3.0 * Constants::PI)) < 0.01);
    }

    SECTION("Convergence at discontinuity")
    {
        // At discontinuity (x=0), Fourier series converges to midpoint (0)
        for (int N : {5, 10, 20})
        {
            Vector<Real> coeffs = basis.Expand(sw, N);
            Real at_zero = basis.EvaluateExpansion(coeffs, 0.0);
            REQUIRE(std::abs(at_zero) < 0.1);  // Should be close to 0
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// IOrthogonalBasis Interface Tests
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("IOrthogonalBasis interface", "[fourier][basis][orthogonal]")
{
    SECTION("FourierBasis implements interface")
    {
        // Verify FourierBasis can be used polymorphically
        FourierBasis concrete(Constants::PI);
        OrthogonalBasis& basis = concrete;
        
        // Access through interface
        REQUIRE(std::abs(basis.DomainMin() - (-Constants::PI)) < 1e-10);
        REQUIRE(std::abs(basis.DomainMax() - Constants::PI) < 1e-10);
        REQUIRE(std::abs(basis.WeightFunction(0.0) - 1.0) < 1e-10);
        REQUIRE(std::abs(basis.Evaluate(0, 0.0) - 1.0) < 1e-10);
        REQUIRE(std::abs(basis.Normalization(0) - 2.0 * Constants::PI) < 1e-10);
    }

    SECTION("ComputeCoefficient through interface")
    {
        FourierBasis concrete(Constants::PI);
        OrthogonalBasis& basis = concrete;
        
        SinFunc f(1.0);
        Real c2 = basis.ComputeCoefficient(f, 2);  // Should be 1.0
        REQUIRE(std::abs(c2 - 1.0) < 1e-6);
    }
}

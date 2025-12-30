///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        FourierBasis.h                                                      ///
///  Description: Fourier basis functions as orthogonal function system               ///
///               Sine/cosine basis for function decomposition                        ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_FOURIER_BASIS_H
#define MML_FOURIER_BASIS_H

#include "MMLBase.h"

#include "interfaces/IFunction.h"
#include "core/OrthogonalBasis.h"

#include "base/Vector.h"
#include "core/Integration/Integration1D.h"

#include "FourierSeries.h"

#include <cmath>
#include <stdexcept>

namespace MML
{
    ///////////////////////////////////////////////////////////////////////////////////////////
    // FourierBasis - Fourier trigonometric basis on [-L, L]
    //
    // Basis functions (indexed from 0):
    //   φ₀(x) = 1                     (constant)
    //   φ₁(x) = cos(πx/L)             (first cosine)
    //   φ₂(x) = sin(πx/L)             (first sine)
    //   φ₃(x) = cos(2πx/L)            (second cosine)
    //   φ₄(x) = sin(2πx/L)            (second sine)
    //   ...
    //
    // Indexing scheme:
    //   n = 0:       constant
    //   n = 2k-1:    cos(kπx/L)  for k = 1, 2, 3, ...
    //   n = 2k:      sin(kπx/L)  for k = 1, 2, 3, ...
    //
    // Orthogonality: ∫[-L,L] φₘ(x)φₙ(x) dx = 0 for m ≠ n
    // Normalization: 
    //   ||φ₀||² = 2L
    //   ||φₙ||² = L for n > 0
    ///////////////////////////////////////////////////////////////////////////////////////////
    class FourierBasis : public OrthogonalBasis
    {
    private:
        Real _L;  // Half-period (domain is [-L, L])

    public:
        explicit FourierBasis(Real L = Constants::PI) : _L(L)
        {
            if (L <= 0.0)
                throw std::invalid_argument("FourierBasis: L must be positive");
        }

        Real HalfPeriod() const { return _L; }
        Real Period() const { return 2.0 * _L; }

        // Domain bounds
        Real DomainMin() const override { return -_L; }
        Real DomainMax() const override { return _L; }

        // Weight function (uniform for Fourier)
        Real WeightFunction(Real /*x*/) const override { return 1.0; }

        // Normalization: ||φₙ||² = ∫[-L,L] φₙ²(x) dx
        Real Normalization(int n) const override
        {
            if (n < 0)
                throw std::invalid_argument("FourierBasis::Normalization: n must be non-negative");
            
            if (n == 0)
                return 2.0 * _L;  // ∫[-L,L] 1² dx = 2L
            else
                return _L;        // ∫[-L,L] cos²(kπx/L) dx = L (also for sin²)
        }

        // Evaluate basis function φₙ at x
        Real Evaluate(int n, Real x) const override
        {
            if (n < 0)
                throw std::invalid_argument("FourierBasis::Evaluate: n must be non-negative");
            
            if (n == 0)
                return 1.0;
            
            int k = (n + 1) / 2;  // Harmonic number
            Real arg = k * Constants::PI * x / _L;
            
            if (n % 2 == 1)  // Odd index → cosine
                return std::cos(arg);
            else             // Even index → sine
                return std::sin(arg);
        }

        // Get harmonic number for basis index
        static int HarmonicNumber(int n)
        {
            if (n <= 0) return 0;
            return (n + 1) / 2;
        }

        // Is cosine basis function?
        static bool IsCosine(int n)
        {
            return (n == 0) || (n % 2 == 1);
        }

        // Is sine basis function?
        static bool IsSine(int n)
        {
            return (n > 0) && (n % 2 == 0);
        }

        // Convert to traditional (aₙ, bₙ) indexing
        // n=0 → (a, 0), n=1 → (a, 1), n=2 → (b, 1), n=3 → (a, 2), n=4 → (b, 2), ...
        static std::pair<char, int> ToTraditionalIndex(int n)
        {
            if (n == 0) return {'a', 0};
            int k = (n + 1) / 2;
            return (n % 2 == 1) ? std::make_pair('a', k) : std::make_pair('b', k);
        }

        // Convert from traditional indexing to basis index
        static int FromTraditionalIndex(char type, int k)
        {
            if (type == 'a' && k == 0) return 0;
            if (type == 'a') return 2 * k - 1;
            if (type == 'b') return 2 * k;
            throw std::invalid_argument("FromTraditionalIndex: type must be 'a' or 'b'");
        }

        ///////////////////////////////////////////////////////////////////////////////////////////
        // Expansion Methods
        ///////////////////////////////////////////////////////////////////////////////////////////

        // Expand function in Fourier basis up to N terms (2N+1 basis functions)
        Vector<Real> Expand(const IRealFunction& f, int N, Real eps = 1e-10) const
        {
            int num_terms = 2 * N + 1;  // φ₀, φ₁, φ₂, ..., φ₂ₙ
            Vector<Real> coefficients(num_terms);
            
            for (int n = 0; n < num_terms; n++)
            {
                coefficients[n] = ComputeCoefficient(f, n, eps);
            }
            
            return coefficients;
        }

        // Evaluate expansion at point x
        Real EvaluateExpansion(const Vector<Real>& coefficients, Real x) const
        {
            Real sum = 0.0;
            for (int n = 0; n < static_cast<int>(coefficients.size()); n++)
            {
                sum += coefficients[n] * Evaluate(n, x);
            }
            return sum;
        }

        // Convert expansion to FourierSeries object
        FourierSeries ToFourierSeries(const Vector<Real>& coefficients) const
        {
            int num_terms = static_cast<int>(coefficients.size());
            int N = (num_terms - 1) / 2;  // Number of harmonics
            
            Vector<Real> a(N + 1);  // a₀, a₁, ..., aₙ
            Vector<Real> b(N);      // b₁, b₂, ..., bₙ
            
            // a₀ = 2 * c₀ (because FourierSeries stores a₀, not a₀/2)
            a[0] = 2.0 * coefficients[0];
            
            for (int k = 1; k <= N; k++)
            {
                int cos_idx = 2 * k - 1;
                int sin_idx = 2 * k;
                
                if (cos_idx < num_terms)
                    a[k] = coefficients[cos_idx];
                else
                    a[k] = 0.0;
                
                if (sin_idx < num_terms)
                    b[k - 1] = coefficients[sin_idx];
                else
                    b[k - 1] = 0.0;
            }
            
            return FourierSeries(a, b, _L);
        }

        // Create expansion from FourierSeries
        Vector<Real> FromFourierSeries(const FourierSeries& fs) const
        {
            if (std::abs(fs.HalfPeriod() - _L) > 1e-10)
                throw std::invalid_argument("FromFourierSeries: period mismatch");
            
            int N = fs.NumTerms();
            int num_terms = 2 * N + 1;
            Vector<Real> coefficients(num_terms);
            
            const auto& a = fs.CosineCoefficients();
            const auto& b = fs.SineCoefficients();
            
            // c₀ = a₀/2 (FourierSeries stores full a₀)
            coefficients[0] = a[0] / 2.0;
            
            for (int k = 1; k <= N; k++)
            {
                coefficients[2 * k - 1] = a[k];  // Cosine coefficient
                coefficients[2 * k] = b[k - 1]; // Sine coefficient
            }
            
            return coefficients;
        }
    };

    ///////////////////////////////////////////////////////////////////////////////////////////
    // ComplexFourierBasis - Complex exponential basis on [-L, L]
    //
    // Basis functions: φₙ(x) = e^(inπx/L) for n = ..., -2, -1, 0, 1, 2, ...
    //
    // Orthogonality: ∫[-L,L] φₘ*(x)φₙ(x) dx = 2L·δₘₙ
    //
    // This basis uses complex inner product: ⟨f,g⟩ = ∫ f*(x)g(x) dx
    ///////////////////////////////////////////////////////////////////////////////////////////
    class ComplexFourierBasis
    {
    private:
        Real _L;  // Half-period

    public:
        explicit ComplexFourierBasis(Real L = Constants::PI) : _L(L)
        {
            if (L <= 0.0)
                throw std::invalid_argument("ComplexFourierBasis: L must be positive");
        }

        Real HalfPeriod() const { return _L; }
        Real DomainMin() const { return -_L; }
        Real DomainMax() const { return _L; }

        // Evaluate basis function e^(inπx/L)
        Complex Evaluate(int n, Real x) const
        {
            Real arg = n * Constants::PI * x / _L;
            return Complex(std::cos(arg), std::sin(arg));
        }

        // Normalization (all basis functions have same norm)
        Real Normalization() const { return 2.0 * _L; }

        // Expand function in complex exponential basis
        // Returns coefficients c₋ₙ, ..., c₋₁, c₀, c₁, ..., cₙ
        Vector<Complex> Expand(const IRealFunction& f, int N, Real eps = 1e-10) const
        {
            int num_terms = 2 * N + 1;
            Vector<Complex> coefficients(num_terms);
            
            Real norm = Normalization();
            
            for (int n = -N; n <= N; n++)
            {
                // cₙ = (1/2L) ∫[-L,L] f(x)·e^(-inπx/L) dx
                class RealPart : public IRealFunction {
                    const IRealFunction& _f;
                    Real _L;
                    int _n;
                public:
                    RealPart(const IRealFunction& f, Real L, int n) : _f(f), _L(L), _n(n) {}
                    Real operator()(Real x) const override {
                        return _f(x) * std::cos(-_n * Constants::PI * x / _L);
                    }
                };
                
                class ImagPart : public IRealFunction {
                    const IRealFunction& _f;
                    Real _L;
                    int _n;
                public:
                    ImagPart(const IRealFunction& f, Real L, int n) : _f(f), _L(L), _n(n) {}
                    Real operator()(Real x) const override {
                        return _f(x) * std::sin(-_n * Constants::PI * x / _L);
                    }
                };
                
                RealPart real_integrand(f, _L, n);
                ImagPart imag_integrand(f, _L, n);
                
                Real real_part = IntegrateTrap(real_integrand, -_L, _L, eps).value / norm;
                Real imag_part = IntegrateTrap(imag_integrand, -_L, _L, eps).value / norm;
                
                coefficients[n + N] = Complex(real_part, imag_part);
            }
            
            return coefficients;
        }

        // Evaluate expansion at point x
        Complex EvaluateExpansion(const Vector<Complex>& coefficients, Real x) const
        {
            int N = (static_cast<int>(coefficients.size()) - 1) / 2;
            Complex sum(0.0, 0.0);
            
            for (int n = -N; n <= N; n++)
            {
                sum += coefficients[n + N] * Evaluate(n, x);
            }
            
            return sum;
        }

        // Convert to real form FourierSeries
        // Using: aₙ = 2·Re(cₙ), bₙ = -2·Im(cₙ), a₀ = 2·Re(c₀)
        FourierSeries ToFourierSeries(const Vector<Complex>& coefficients) const
        {
            return FourierSeries::FromComplexCoefficients(coefficients, _L);
        }
    };

} // namespace MML

#endif // MML_FOURIER_BASIS_H
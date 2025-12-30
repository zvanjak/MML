///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        FourierSeries.h                                                     ///
///  Description: Fourier series expansion and coefficient computation                ///
///               Periodic function approximation via trigonometric series            ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_FOURIER_SERIES_H
#define MML_FOURIER_SERIES_H

#include "MMLBase.h"

#include "interfaces/IFunction.h"
#include "base/Vector.h"
#include "core/Integration/Integration1D.h"

#include <cmath>

namespace MML
{
    ///////////////////////////////////////////////////////////////////////////////////////////
    // FourierSeries - Continuous Fourier series representation
    //
    // Represents a periodic function as:
    //   f(x) = a₀/2 + Σ[aₙ·cos(nπx/L) + bₙ·sin(nπx/L)]  for n = 1, 2, ..., N
    //
    // where:
    //   - L is the half-period (function is 2L-periodic)
    //   - aₙ are cosine coefficients
    //   - bₙ are sine coefficients
    //
    // Key features:
    //   - Approximates periodic functions with spectral convergence (for smooth functions)
    //   - Closed-form derivatives and integrals
    //   - Parseval's theorem for energy computation
    //   - Conversion to/from complex exponential form
    //
    // Usage:
    //   auto sin_func = [](Real x) { return std::sin(x); };
    //   FourierSeries fs(sin_func, Constants::PI, 10);  // Approximate sin(x) on [-π, π]
    //   Real value = fs(0.5);  // Evaluate at x=0.5
    ///////////////////////////////////////////////////////////////////////////////////////////
    class FourierSeries : public IRealFunction
    {
    private:
        Vector<Real> _a;    // Cosine coefficients [a₀, a₁, a₂, ..., aₙ]
        Vector<Real> _b;    // Sine coefficients [b₁, b₂, ..., bₙ]  (b₀ = 0 by definition)
        Real         _L;    // Half-period (function has period 2L)
        int          _N;    // Number of terms

    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // Constructors
        ///////////////////////////////////////////////////////////////////////////////////////////

        // Default constructor - zero series
        FourierSeries() : _L(Constants::PI), _N(0) {}

        // Construct from coefficients
        FourierSeries(const Vector<Real>& a, const Vector<Real>& b, Real L)
            : _a(a), _b(b), _L(L)
        {
            if (_a.size() == 0)
                throw std::invalid_argument("FourierSeries: coefficient vector _a must have at least a₀");
            if (_b.size() != _a.size() - 1)
                throw std::invalid_argument("FourierSeries: _b must have size = _a.size() - 1");
            
            _N = static_cast<int>(_a.size()) - 1;  // N terms (excluding a₀)
        }

        // Construct from function by computing Fourier coefficients
        // Uses adaptive quadrature for coefficient computation
        FourierSeries(const IRealFunction& func, Real L, int N, Real eps = 1e-10)
            : _L(L), _N(N)
        {
            if (L <= 0.0)
                throw std::invalid_argument("FourierSeries: L must be positive");
            if (N < 0)
                throw std::invalid_argument("FourierSeries: N must be non-negative");

            _a.Resize(N + 1);  // a₀, a₁, ..., aₙ
            _b.Resize(N);      // b₁, ..., bₙ

            // Compute a₀ = (1/L) ∫[-L,L] f(x) dx
            Real sum_a0 = IntegrateTrap(func, -L, L, eps).value;
            _a[0] = sum_a0 / L;

            // Compute aₙ = (1/L) ∫[-L,L] f(x)·cos(nπx/L) dx
            // and bₙ = (1/L) ∫[-L,L] f(x)·sin(nπx/L) dx
            for (int n = 1; n <= N; n++)
            {
                // Cosine coefficient using helper class
                class CosIntegrand : public IRealFunction {
                    const IRealFunction& _f;
                    Real _L;
                    int _n;
                public:
                    CosIntegrand(const IRealFunction& f, Real L, int n) : _f(f), _L(L), _n(n) {}
                    Real operator()(Real x) const override {
                        return _f(x) * std::cos(_n * Constants::PI * x / _L);
                    }
                };
                CosIntegrand cos_integrand(func, L, n);
                _a[n] = IntegrateTrap(cos_integrand, -L, L, eps).value / L;

                // Sine coefficient using helper class
                class SinIntegrand : public IRealFunction {
                    const IRealFunction& _f;
                    Real _L;
                    int _n;
                public:
                    SinIntegrand(const IRealFunction& f, Real L, int n) : _f(f), _L(L), _n(n) {}
                    Real operator()(Real x) const override {
                        return _f(x) * std::sin(_n * Constants::PI * x / _L);
                    }
                };
                SinIntegrand sin_integrand(func, L, n);
                _b[n - 1] = IntegrateTrap(sin_integrand, -L, L, eps).value / L;
            }
        }

        ///////////////////////////////////////////////////////////////////////////////////////////
        // Evaluation
        ///////////////////////////////////////////////////////////////////////////////////////////

        // Evaluate Fourier series at x
        Real operator()(Real x) const override
        {
            if (_N == 0)
                return _a[0] / 2.0;

            Real sum = _a[0] / 2.0;
            
            for (int n = 1; n <= _N; n++)
            {
                Real arg = n * Constants::PI * x / _L;
                sum += _a[n] * std::cos(arg) + _b[n - 1] * std::sin(arg);
            }
            
            return sum;
        }

        ///////////////////////////////////////////////////////////////////////////////////////////
        // Calculus Operations (closed-form for Fourier series!)
        ///////////////////////////////////////////////////////////////////////////////////////////

        // Derivative: d/dx[aₙcos(nπx/L) + bₙsin(nπx/L)]
        //           = -aₙ(nπ/L)sin(nπx/L) + bₙ(nπ/L)cos(nπx/L)
        // Maps to: a'ₙ = bₙ(nπ/L), b'ₙ = -aₙ(nπ/L)
        FourierSeries Derivative() const
        {
            if (_N == 0)
            {
                // Derivative of constant is zero
                Vector<Real> new_a(1, 0.0);
                Vector<Real> new_b(0);
                return FourierSeries(new_a, new_b, _L);
            }

            Vector<Real> new_a(_N + 1);
            Vector<Real> new_b(_N);

            new_a[0] = 0.0;  // Derivative of constant term is zero

            for (int n = 1; n <= _N; n++)
            {
                Real factor = n * Constants::PI / _L;
                new_a[n] = _b[n - 1] * factor;       // a'ₙ = bₙ·(nπ/L)
                new_b[n - 1] = -_a[n] * factor;      // b'ₙ = -aₙ·(nπ/L)
            }

            return FourierSeries(new_a, new_b, _L);
        }

        // Integral: ∫[aₙcos(nπx/L) + bₙsin(nπx/L)]dx
        //         = aₙ(L/nπ)sin(nπx/L) - bₙ(L/nπ)cos(nπx/L)
        // Maps to: A'ₙ = -bₙ(L/nπ), B'ₙ = aₙ(L/nπ)
        // Note: a₀ term integrates to (a₀/2)x, but we drop it (indefinite integral constant)
        FourierSeries Integral() const
        {
            if (_N == 0)
            {
                // Integral of constant is a linear function (not in Fourier series space)
                // Return zero series (or throw exception)
                Vector<Real> new_a(1, 0.0);
                Vector<Real> new_b(0);
                return FourierSeries(new_a, new_b, _L);
            }

            Vector<Real> new_a(_N + 1);
            Vector<Real> new_b(_N);

            new_a[0] = 0.0;  // Integration constant (arbitrary)

            for (int n = 1; n <= _N; n++)
            {
                Real factor = _L / (n * Constants::PI);
                new_a[n] = -_b[n - 1] * factor;      // Aₙ = -bₙ·(L/nπ)
                new_b[n - 1] = _a[n] * factor;       // Bₙ = aₙ·(L/nπ)
            }

            return FourierSeries(new_a, new_b, _L);
        }

        ///////////////////////////////////////////////////////////////////////////////////////////
        // Coefficient Access
        ///////////////////////////////////////////////////////////////////////////////////////////

        const Vector<Real>& CosineCoefficients() const { return _a; }
        const Vector<Real>& SineCoefficients() const { return _b; }
        Real HalfPeriod() const { return _L; }
        Real Period() const { return 2.0 * _L; }
        int NumTerms() const { return _N; }

        ///////////////////////////////////////////////////////////////////////////////////////////
        // Complex Exponential Form
        ///////////////////////////////////////////////////////////////////////////////////////////

        // Convert to complex exponential form: f(x) = Σ cₙ·e^(inπx/L)
        // Relationship: cₙ = (aₙ - i·bₙ)/2 for n > 0
        //               c₀ = a₀/2
        //               c₋ₙ = (aₙ + i·bₙ)/2 = c̄ₙ (complex conjugate)
        Vector<Complex> ComplexCoefficients() const
        {
            Vector<Complex> c(2 * _N + 1);
            
            c[_N] = Complex(_a[0] / 2.0, 0.0);  // c₀
            
            for (int n = 1; n <= _N; n++)
            {
                Complex cn(_a[n] / 2.0, -_b[n - 1] / 2.0);  // cₙ = (aₙ - i·bₙ)/2
                c[_N + n] = cn;                              // Positive frequencies
                c[_N - n] = std::conj(cn);                   // Negative frequencies (conjugate)
            }
            
            return c;
        }

        // Construct from complex coefficients
        // Input: c = [c₋ₙ, ..., c₋₁, c₀, c₁, ..., cₙ] (length 2N+1)
        static FourierSeries FromComplexCoefficients(const Vector<Complex>& c, Real L)
        {
            int len = static_cast<int>(c.size());
            if (len % 2 == 0)
                throw std::invalid_argument("FromComplexCoefficients: c must have odd length");
            
            int N = (len - 1) / 2;
            Vector<Real> a(N + 1);
            Vector<Real> b(N);
            
            a[0] = 2.0 * c[N].real();  // a₀ = 2·Re(c₀)
            
            for (int n = 1; n <= N; n++)
            {
                Complex cn = c[N + n];
                a[n] = 2.0 * cn.real();     // aₙ = 2·Re(cₙ)
                b[n - 1] = -2.0 * cn.imag(); // bₙ = -2·Im(cₙ)
            }
            
            return FourierSeries(a, b, L);
        }

        ///////////////////////////////////////////////////////////////////////////////////////////
        // Energy and Norms (Parseval's Theorem)
        ///////////////////////////////////////////////////////////////////////////////////////////

        // Parseval's theorem: ∫[-L,L] |f(x)|² dx = L·(a₀²/2 + Σ(aₙ² + bₙ²))
        Real Energy() const
        {
            Real energy = (_a[0] * _a[0]) / 2.0;  // a₀²/2
            
            for (int n = 1; n <= _N; n++)
            {
                energy += _a[n] * _a[n] + _b[n - 1] * _b[n - 1];
            }
            
            return _L * energy;
        }

        // L² norm: ||f||₂ = √(Energy/2L)
        Real L2Norm() const
        {
            return std::sqrt(Energy() / (2.0 * _L));
        }

        ///////////////////////////////////////////////////////////////////////////////////////////
        // Series Operations
        ///////////////////////////////////////////////////////////////////////////////////////////

        // Truncate to M terms (M < N)
        FourierSeries Truncate(int M) const
        {
            if (M < 0)
                throw std::invalid_argument("Truncate: M must be non-negative");
            if (M > _N)
                return *this;
            
            Vector<Real> new_a(M + 1);
            Vector<Real> new_b(M);
            
            for (int i = 0; i <= M; i++)
                new_a[i] = _a[i];
            
            for (int i = 0; i < M; i++)
                new_b[i] = _b[i];
            
            return FourierSeries(new_a, new_b, _L);
        }

        // Addition
        FourierSeries operator+(const FourierSeries& other) const
        {
            if (std::abs(_L - other._L) > 1e-10)
                throw std::invalid_argument("FourierSeries::operator+: periods must match");
            
            int max_N = std::max(_N, other._N);
            Vector<Real> new_a(max_N + 1, 0.0);
            Vector<Real> new_b(max_N, 0.0);
            
            for (int i = 0; i <= _N; i++)
                new_a[i] = _a[i];
            for (int i = 0; i < _N; i++)
                new_b[i] = _b[i];
            
            for (int i = 0; i <= other._N; i++)
                new_a[i] += other._a[i];
            for (int i = 0; i < other._N; i++)
                new_b[i] += other._b[i];
            
            return FourierSeries(new_a, new_b, _L);
        }

        // Subtraction
        FourierSeries operator-(const FourierSeries& other) const
        {
            if (std::abs(_L - other._L) > 1e-10)
                throw std::invalid_argument("FourierSeries::operator-: periods must match");
            
            int max_N = std::max(_N, other._N);
            Vector<Real> new_a(max_N + 1, 0.0);
            Vector<Real> new_b(max_N, 0.0);
            
            for (int i = 0; i <= _N; i++)
                new_a[i] = _a[i];
            for (int i = 0; i < _N; i++)
                new_b[i] = _b[i];
            
            for (int i = 0; i <= other._N; i++)
                new_a[i] -= other._a[i];
            for (int i = 0; i < other._N; i++)
                new_b[i] -= other._b[i];
            
            return FourierSeries(new_a, new_b, _L);
        }

        // Scalar multiplication
        FourierSeries operator*(Real scalar) const
        {
            Vector<Real> new_a(_N + 1);
            Vector<Real> new_b(_N);
            
            for (int i = 0; i <= _N; i++)
                new_a[i] = _a[i] * scalar;
            for (int i = 0; i < _N; i++)
                new_b[i] = _b[i] * scalar;
            
            return FourierSeries(new_a, new_b, _L);
        }

        friend FourierSeries operator*(Real scalar, const FourierSeries& fs)
        {
            return fs * scalar;
        }

        ///////////////////////////////////////////////////////////////////////////////////////////
        // Convergence Analysis
        ///////////////////////////////////////////////////////////////////////////////////////////

        // Get magnitude of nth coefficient (for convergence analysis)
        Real CoefficientMagnitude(int n) const
        {
            if (n == 0)
                return std::abs(_a[0]);
            if (n > _N)
                return 0.0;
            
            return std::sqrt(_a[n] * _a[n] + _b[n - 1] * _b[n - 1]);
        }

        // Maximum coefficient magnitude (indicator of convergence)
        Real MaxCoefficientMagnitude() const
        {
            Real max_mag = std::abs(_a[0]);
            
            for (int n = 1; n <= _N; n++)
            {
                Real mag = std::sqrt(_a[n] * _a[n] + _b[n - 1] * _b[n - 1]);
                if (mag > max_mag)
                    max_mag = mag;
            }
            
            return max_mag;
        }
    };

} // namespace MML

#endif // MML_FOURIER_SERIES_H

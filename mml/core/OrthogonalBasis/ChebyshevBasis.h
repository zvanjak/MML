///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        ChebyshevBasis.h                                                    ///
///  Description: Chebyshev polynomial basis functions (first and second kind)        ///
///               Orthogonal on [-1, 1] with weight w(x) = 1/√(1-x²)                 ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_CHEBYSHEV_BASIS_H
#define MML_CHEBYSHEV_BASIS_H

#include "MMLBase.h"

#include "algorithms/ChebyshevApproximation.h"
#include "interfaces/IFunction.h"
#include "core/OrthogonalBasis.h"

#include <cmath>
#include <stdexcept>

namespace MML
{
    ///////////////////////////////////////////////////////////////////////////////////////////
    // ChebyshevBasis - Chebyshev polynomials of the first kind Tₙ(x)
    //
    // Basis functions: Tₙ(x) - Chebyshev polynomials of the first kind
    //
    // Orthogonality: ∫₋₁¹ Tₘ(x)Tₙ(x)/√(1-x²) dx = 0 for m ≠ n
    // Normalization: ∫₋₁¹ Tₙ²(x)/√(1-x²) dx = { π    if n = 0
    //                                           { π/2  if n > 0
    //
    // Domain: [-1, 1]
    // Weight function: w(x) = 1/√(1-x²)
    //
    // Recurrence relation:
    //   T₀(x) = 1
    //   T₁(x) = x
    //   Tₙ₊₁(x) = 2xTₙ(x) - Tₙ₋₁(x)
    //
    // Special property: Tₙ(cos θ) = cos(nθ)
    //
    // Applications:
    //   - Polynomial approximation (minimax property)
    //   - Spectral methods (Chebyshev collocation)
    //   - Signal processing (DCT, discrete cosine transform)
    //   - Numerical analysis (Clenshaw-Curtis quadrature)
    //
    // Properties:
    //   - |Tₙ(x)| ≤ 1 for x ∈ [-1, 1]
    //   - Extrema at xₖ = cos(kπ/n), k = 0, 1, ..., n
    //   - Zeros at xₖ = cos((2k-1)π/(2n)), k = 1, 2, ..., n
    ///////////////////////////////////////////////////////////////////////////////////////////
    class ChebyshevBasis : public OrthogonalBasis
    {
    public:
        ChebyshevBasis() = default;

        // Evaluate Chebyshev polynomial Tₙ(x)
        Real Evaluate(int n, Real x) const override
        {
            if (n < 0)
                throw std::invalid_argument("ChebyshevBasis::Evaluate: n must be non-negative");
            
            return ChebyshevT(n, x);
        }

        // Weight function w(x) = 1/√(1-x²)
        Real WeightFunction(Real x) const override 
        { 
            if (std::abs(x) >= 1.0)
            {
                // At boundaries, weight is undefined but we can use a large value
                // For practical computation, return a finite value
                return 1e10;  // Singular at x = ±1
            }
            return 1.0 / std::sqrt(1.0 - x * x); 
        }

        // Normalization: ||Tₙ||² = { π    if n = 0
        //                           { π/2  if n > 0
        Real Normalization(int n) const override
        {
            if (n < 0)
                throw std::invalid_argument("ChebyshevBasis::Normalization: n must be non-negative");
            
            if (n == 0)
                return Constants::PI;
            else
                return Constants::PI / 2.0;
        }

        // Domain bounds
        Real DomainMin() const override { return -1.0; }
        Real DomainMax() const override { return 1.0; }

        // Recurrence coefficients for three-term recurrence:
        // Tₙ₊₁(x) = (aₙx + bₙ)Tₙ(x) + cₙTₙ₋₁(x)
        // 
        // Standard form: Tₙ₊₁ = 2xTₙ - Tₙ₋₁
        // So: aₙ = 2, bₙ = 0, cₙ = -1
        void RecurrenceCoefficients(int n, Real& a, Real& b, Real& c) const
        {
            if (n < 0)
                throw std::invalid_argument("ChebyshevBasis::RecurrenceCoefficients: n must be non-negative");
            
            a = 2.0;
            b = 0.0;
            c = -1.0;
        }

        // Chebyshev extrema: xₖ = cos(kπ/n) for k = 0, 1, ..., n
        std::vector<Real> GetExtrema(int n) const
        {
            if (n < 0)
                throw std::invalid_argument("ChebyshevBasis::GetExtrema: n must be non-negative");
            
            std::vector<Real> extrema(n + 1);
            for (int k = 0; k <= n; k++)
                extrema[k] = std::cos(k * Constants::PI / n);
            return extrema;
        }

        // Chebyshev zeros (roots): xₖ = cos((2k-1)π/(2n)) for k = 1, 2, ..., n
        std::vector<Real> GetZeros(int n) const
        {
            if (n <= 0)
                throw std::invalid_argument("ChebyshevBasis::GetZeros: n must be positive");
            
            std::vector<Real> zeros(n);
            for (int k = 1; k <= n; k++)
                zeros[k-1] = std::cos((2.0*k - 1.0) * Constants::PI / (2.0 * n));
            return zeros;
        }
    };

    ///////////////////////////////////////////////////////////////////////////////////////////
    // ChebyshevBasisSecondKind - Chebyshev polynomials of the second kind Uₙ(x)
    //
    // Basis functions: Uₙ(x) - Chebyshev polynomials of the second kind
    //
    // Orthogonality: ∫₋₁¹ Uₘ(x)Uₙ(x)√(1-x²) dx = 0 for m ≠ n
    // Normalization: ∫₋₁¹ Uₙ²(x)√(1-x²) dx = π/2
    //
    // Domain: [-1, 1]
    // Weight function: w(x) = √(1-x²)
    //
    // Recurrence relation:
    //   U₀(x) = 1
    //   U₁(x) = 2x
    //   Uₙ₊₁(x) = 2xUₙ(x) - Uₙ₋₁(x)
    //
    // Special property: Uₙ(cos θ) = sin((n+1)θ)/sin(θ)
    //
    // Relation to Tₙ: Uₙ(x) = dTₙ₊₁(x)/dx / (n+1)
    //                 (1-x²)Uₙ(x) = xTₙ₊₁(x) - Tₙ₊₂(x)
    ///////////////////////////////////////////////////////////////////////////////////////////
    class ChebyshevBasisSecondKind
    {
    public:
        ChebyshevBasisSecondKind() = default;

        // Evaluate Chebyshev polynomial Uₙ(x) of the second kind
        Real Evaluate(int n, Real x) const
        {
            if (n < 0)
                throw std::invalid_argument("ChebyshevBasisSecondKind::Evaluate: n must be non-negative");
            
            return ChebyshevU(n, x);
        }

        // Weight function w(x) = √(1-x²)
        Real WeightFunction(Real x) const 
        { 
            if (std::abs(x) >= 1.0)
                return 0.0;  // Zero at boundaries
            return std::sqrt(1.0 - x * x); 
        }

        // Normalization: ||Uₙ||² = π/2 for all n ≥ 0
        Real Normalization(int n) const
        {
            if (n < 0)
                throw std::invalid_argument("ChebyshevBasisSecondKind::Normalization: n must be non-negative");
            
            return Constants::PI / 2.0;
        }

        // Domain bounds
        Real DomainMin() const { return -1.0; }
        Real DomainMax() const { return 1.0; }

        // Recurrence coefficients (same as first kind)
        void RecurrenceCoefficients(int n, Real& a, Real& b, Real& c) const
        {
            if (n < 0)
                throw std::invalid_argument("ChebyshevBasisSecondKind::RecurrenceCoefficients: n must be non-negative");
            
            a = 2.0;
            b = 0.0;
            c = -1.0;
        }

        // Zeros of Uₙ: xₖ = cos(kπ/(n+1)) for k = 1, 2, ..., n
        std::vector<Real> GetZeros(int n) const
        {
            if (n <= 0)
                throw std::invalid_argument("ChebyshevBasisSecondKind::GetZeros: n must be positive");
            
            std::vector<Real> zeros(n);
            for (int k = 1; k <= n; k++)
                zeros[k-1] = std::cos(k * Constants::PI / (n + 1.0));
            return zeros;
        }
    };

} // namespace MML

#endif // MML_CHEBYSHEV_BASIS_H

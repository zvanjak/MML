///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        OrthogonalBasis.h                                                   ///
///  Description: Orthogonal function bases (Legendre, Hermite, Laguerre)             ///
///               Gram-Schmidt orthogonalization and inner products                   ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_ORTHOGONAL_BASIS_H
#define MML_ORTHOGONAL_BASIS_H

#include <mml/MMLBase.h>

#include <mml/interfaces/IFunction.h>

#include <mml/core/Integration/Integration1D.h>

#include <cmath>
#include <stdexcept>

namespace MML
{
    ///////////////////////////////////////////////////////////////////////////////////////////
    // OrthogonalBasis - Interface for orthogonal function bases
    //
    // Defines the contract for orthogonal basis systems like:
    //   - Fourier basis (sin/cos or complex exponentials)
    //   - Legendre polynomials
    //   - Chebyshev polynomials
    //   - Hermite polynomials
    //   - Laguerre polynomials
    //
    // Key properties:
    //   - Basis functions φₙ(x) satisfy orthogonality: ⟨φₘ, φₙ⟩ = 0 for m ≠ n
    //   - Inner product defined with weight function: ⟨f,g⟩ = ∫ f(x)g(x)w(x) dx
    //   - Each basis function has a normalization constant
    ///////////////////////////////////////////////////////////////////////////////////////////
    class OrthogonalBasis
    {
    protected:
        Real CoefficientNormalization(int n, Real eps) const
        {
            if (n < 0 || !std::isfinite(eps) || eps <= 0)
                throw ArgumentError("OrthogonalBasis::ComputeCoefficient: require n >= 0 and finite positive tolerance");
            const Real norm = Normalization(n);
            if (!std::isfinite(norm) || norm <= 0)
                throw NumericalMethodError("OrthogonalBasis::ComputeCoefficient: normalization must be finite and positive");
            return norm;
        }

    public:
        virtual ~OrthogonalBasis() = default;

        // Evaluate the nth basis function at x
        virtual Real Evaluate(int n, Real x) const = 0;

        // Weight function for inner product
        virtual Real WeightFunction(Real x) const = 0;

        // Normalization constant for nth basis function: ||φₙ||²
        virtual Real Normalization(int n) const = 0;

        // Domain bounds
        virtual Real DomainMin() const = 0;
        virtual Real DomainMax() const = 0;

        // Compute expansion coefficient for function f: cₙ = ⟨f, φₙ⟩ / ||φₙ||²
        virtual Real ComputeCoefficient(const IRealFunction& f, int n, Real eps = Precision::OrthogonalBasisCoefficientTolerance) const
        {
            const Real norm = CoefficientNormalization(n, eps);
            const Real a = DomainMin(), b = DomainMax();
            if (std::isnan(a) || std::isnan(b) || a >= b)
                throw ArgumentError("OrthogonalBasis::ComputeCoefficient: invalid basis interval");
            auto contribution = [&](Real x, Real jacobian) {
                IntegrationDetail::RequireFinite(x, "Basis coefficient coordinate");
                return IntegrationDetail::Product({f(x), Evaluate(n, x), WeightFunction(x), jacobian}, norm);
            };
            if (std::isfinite(a) && std::isfinite(b))
                return Integration::IntegrateGKAdaptive([&](Real x) { return contribution(x, 1); }, a, b, eps, eps);
            // Rational maps cover the complete tails; open GK nodes never sample infinity.
            auto transformed = [&](Real t) {
                const Real distance = (1 - t) / t;
                const Real jacobian = (1 / t) / t;
                if (std::isfinite(a))
                    return contribution(a + distance, jacobian);
                if (std::isfinite(b))
                    return contribution(b - distance, jacobian);
                return contribution(distance, jacobian) + contribution(-distance, jacobian);
            };
            return Integration::IntegrateGKAdaptive(transformed, 0, 1, eps, eps);
        }
    };
  }

  #endif 
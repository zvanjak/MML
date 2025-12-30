///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        OrthogonalBasis.h                                                   ///
///  Description: Orthogonal function bases (Legendre, Hermite, Laguerre)             ///
///               Gram-Schmidt orthogonalization and inner products                   ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_ORTHOGONAL_BASIS_H
#define MML_ORTHOGONAL_BASIS_H

#include "MMLBase.h"

#include "interfaces/IFunction.h"

#include "core/Integration/Integration1D.h"

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
        virtual Real ComputeCoefficient(const IRealFunction& f, int n, Real eps = 1e-10) const
        {
            // Helper class for integrand f(x) * φₙ(x) * w(x)
            class Integrand : public IRealFunction {
                const OrthogonalBasis& _basis;
                const IRealFunction& _f;
                int _n;
            public:
                Integrand(const OrthogonalBasis& basis, const IRealFunction& f, int n)
                    : _basis(basis), _f(f), _n(n) {}
                Real operator()(Real x) const override {
                    return _f(x) * _basis.Evaluate(_n, x) * _basis.WeightFunction(x);
                }
            };

            Integrand integrand(*this, f, n);
            Real inner_product = IntegrateTrap(integrand, DomainMin(), DomainMax(), eps).value;
            return inner_product / Normalization(n);
        }
    };
  }

  #endif 
///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        LegendreBasis.h                                                     ///
///  Description: Legendre polynomial basis functions                                 ///
///               Orthogonal on [-1, 1] with weight w(x) = 1                          ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_LEGENDRE_BASIS_H
#define MML_LEGENDRE_BASIS_H

#include "MMLBase.h"

#include "base/StandardFunctions.h"
#include "interfaces/IFunction.h"
#include "core/OrthogonalBasis.h"

#include <cmath>
#include <stdexcept>

namespace MML
{
    ///////////////////////////////////////////////////////////////////////////////////////////
    // LegendreBasis - Legendre polynomial orthogonal basis
    //
    // Basis functions: Pₙ(x) - Legendre polynomials
    //
    // Orthogonality: ∫₋₁¹ Pₘ(x)Pₙ(x) dx = 0 for m ≠ n
    // Normalization: ∫₋₁¹ Pₙ²(x) dx = 2/(2n+1)
    //
    // Domain: [-1, 1]
    // Weight function: w(x) = 1 (uniform)
    //
    // Recurrence relation (Bonnet's formula):
    //   P₀(x) = 1
    //   P₁(x) = x
    //   (n+1)Pₙ₊₁(x) = (2n+1)xPₙ(x) - nPₙ₋₁(x)
    //
    // Physical applications:
    //   - Solutions to Laplace equation in spherical coordinates
    //   - Angular momentum eigenstates in quantum mechanics
    //   - Multipole expansion in electrostatics
    //   - Scattering theory
    ///////////////////////////////////////////////////////////////////////////////////////////
    class LegendreBasis : public OrthogonalBasis
    {
    public:
        LegendreBasis() = default;

        // Evaluate Legendre polynomial Pₙ(x)
        Real Evaluate(int n, Real x) const override
        {
            if (n < 0)
                throw std::invalid_argument("LegendreBasis::Evaluate: n must be non-negative");
            
            return Functions::Legendre(static_cast<unsigned int>(n), x);
        }

        // Weight function w(x) = 1 (uniform)
        Real WeightFunction(Real /*x*/) const override 
        { 
            return 1.0; 
        }

        // Normalization: ||Pₙ||² = ∫₋₁¹ Pₙ²(x) dx = 2/(2n+1)
        Real Normalization(int n) const override
        {
            if (n < 0)
                throw std::invalid_argument("LegendreBasis::Normalization: n must be non-negative");
            
            return 2.0 / (2.0 * n + 1.0);
        }

        // Domain bounds
        Real DomainMin() const override { return -1.0; }
        Real DomainMax() const override { return 1.0; }

        // Recurrence coefficients for three-term recurrence:
        // Pₙ₊₁(x) = (aₙx + bₙ)Pₙ(x) + cₙPₙ₋₁(x)
        // 
        // Standard form: (n+1)Pₙ₊₁ = (2n+1)xPₙ - nPₙ₋₁
        // So: aₙ = (2n+1)/(n+1), bₙ = 0, cₙ = -n/(n+1)
        void RecurrenceCoefficients(int n, Real& a, Real& b, Real& c) const
        {
            if (n < 0)
                throw std::invalid_argument("LegendreBasis::RecurrenceCoefficients: n must be non-negative");
            
            Real n_real = static_cast<Real>(n);
            a = (2.0 * n_real + 1.0) / (n_real + 1.0);
            b = 0.0;
            c = -n_real / (n_real + 1.0);
        }
    };

    ///////////////////////////////////////////////////////////////////////////////////////////
    // Associated Legendre functions Pₗᵐ(x) for spherical harmonics
    //
    // Used in:
    //   - Spherical harmonics Yₗᵐ(θ,φ) = √[(2l+1)(l-m)!/(4π(l+m)!)] Pₗᵐ(cos θ) e^(imφ)
    //   - Solutions to Laplace equation in spherical coordinates
    //   - Quantum angular momentum states
    ///////////////////////////////////////////////////////////////////////////////////////////
    class AssociatedLegendreBasis
    {
    private:
        int _m;  // Order (fixed for this basis)

    public:
        explicit AssociatedLegendreBasis(int m) : _m(m)
        {
            if (m < 0)
                throw std::invalid_argument("AssociatedLegendreBasis: m must be non-negative");
        }

        int Order() const { return _m; }

        // Evaluate associated Legendre function Pₗᵐ(x) for given degree l
        Real Evaluate(int l, Real x) const
        {
            if (l < _m)
                throw std::invalid_argument("AssociatedLegendreBasis::Evaluate: l must be >= m");
            
            return Functions::SphLegendre(static_cast<unsigned int>(l), static_cast<unsigned int>(_m), x);
        }

        // Domain bounds
        Real DomainMin() const { return -1.0; }
        Real DomainMax() const { return 1.0; }
    };

} // namespace MML

#endif // MML_LEGENDRE_BASIS_H

///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        LaguerreBasis.h                                                     ///
///  Description: Laguerre polynomial basis functions                                 ///
///               Orthogonal on [0, ∞) with weight w(x) = e^(-x)                     ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_LAGUERRE_BASIS_H
#define MML_LAGUERRE_BASIS_H

#include <mml/MMLBase.h>

#include <mml/base/StandardFunctions.h>
#include <mml/interfaces/IFunction.h>
#include <mml/core/OrthogonalBasis.h>

#include <cmath>
#include <limits>
#include <stdexcept>

namespace MML
{
    ///////////////////////////////////////////////////////////////////////////////////////////
    // LaguerreBasis - Laguerre polynomial orthogonal basis
    //
    // Basis functions: Lₙ(x) - Laguerre polynomials
    //
    // Orthogonality: ∫₀^∞ Lₘ(x)Lₙ(x)e^(-x) dx = 0 for m ≠ n
    // Normalization: ∫₀^∞ Lₙ²(x)e^(-x) dx = 1
    //
    // Domain: [0, ∞)
    // Weight function: w(x) = e^(-x)
    //
    // Recurrence relation:
    //   L₀(x) = 1
    //   L₁(x) = 1 - x
    //   (n+1)Lₙ₊₁(x) = (2n + 1 - x)Lₙ(x) - nLₙ₋₁(x)
    //
    // Physical applications:
    //   - Hydrogen atom wavefunctions (radial part)
    //   - Quantum harmonic oscillator in 3D
    //   - Heat conduction in semi-infinite domains
    //   - Time-dependent perturbation theory
    //
    // Properties:
    //   - Lₙ(0) = 1 for all n
    //   - Lₙ solves: xy'' + (1-x)y' + ny = 0  (Laguerre differential equation)
    //   - d/dx Lₙ(x) = -Lₙ₋₁^(1)(x) = -(Lₙ₋₁(x) + Lₙ₋₂(x) + ... + L₀(x))
    ///////////////////////////////////////////////////////////////////////////////////////////
    class LaguerreBasis : public OrthogonalBasis
    {
    public:
        LaguerreBasis() = default;

        // Evaluate Laguerre polynomial Lₙ(x)
        Real Evaluate(int n, Real x) const override
        {
            if (n < 0)
                throw ArgumentError("LaguerreBasis::Evaluate: n must be non-negative");
            
            return Functions::Laguerre(static_cast<unsigned int>(n), x);
        }

        // Weight function w(x) = e^(-x)
        Real WeightFunction(Real x) const override 
        { 
            return std::exp(-x); 
        }

        // Normalization: ||Lₙ||² = ∫₀^∞ Lₙ²(x)e^(-x) dx = 1
        Real Normalization(int n) const override
        {
            if (n < 0)
                throw ArgumentError("LaguerreBasis::Normalization: n must be non-negative");
            
            return 1.0;
        }

        // Domain bounds
        Real DomainMin() const override { return 0.0; }
        Real DomainMax() const override 
        { 
            return std::numeric_limits<Real>::infinity(); 
        }

        // Recurrence coefficients for three-term recurrence:
        // (n+1)Lₙ₊₁(x) = (2n + 1 - x)Lₙ(x) - nLₙ₋₁(x)
        // 
        // Divide by (n+1):
        // Lₙ₊₁(x) = ((2n+1)/(n+1) - x/(n+1))Lₙ(x) - (n/(n+1))Lₙ₋₁(x)
        //
        // So in form Lₙ₊₁ = (aₙx + bₙ)Lₙ + cₙLₙ₋₁:
        // aₙ = -1/(n+1), bₙ = (2n+1)/(n+1), cₙ = -n/(n+1)
        void RecurrenceCoefficients(int n, Real& a, Real& b, Real& c) const
        {
            if (n < 0)
                throw ArgumentError("LaguerreBasis::RecurrenceCoefficients: n must be non-negative");
            
            Real n_real = static_cast<Real>(n);
            Real n_plus_1 = n_real + MML::Real(1.0);
            
            a = -MML::Real(1.0) / n_plus_1;
            b = (MML::Real(2.0) * n_real + MML::Real(1.0)) / n_plus_1;
            c = -n_real / n_plus_1;
        }
    };

    ///////////////////////////////////////////////////////////////////////////////////////////
    // Associated Laguerre polynomials Lₙ^(α)(x)
    //
    // Generalized Laguerre polynomials with parameter α
    //
    // Orthogonality: ∫₀^∞ Lₘ^(α)(x)Lₙ^(α)(x)x^α e^(-x) dx = 0 for m ≠ n
    // Normalization: ∫₀^∞ [Lₙ^(α)(x)]²x^α e^(-x) dx = Γ(n+α+1)/n!
    //
    // Weight function: w(x) = x^α e^(-x)
    //
    // Physical applications:
    //   - Hydrogen atom (α = 2l + 1, where l is angular momentum quantum number)
    //   - 3D quantum harmonic oscillator
    //
    // Relation to simple Laguerre: Lₙ^(0)(x) = Lₙ(x)
    ///////////////////////////////////////////////////////////////////////////////////////////
    class AssociatedLaguerreBasis
    {
    private:
        Real _alpha;  // Order parameter (α > -1)

    public:
        explicit AssociatedLaguerreBasis(Real alpha) : _alpha(alpha)
        {
            if (!std::isfinite(alpha) || alpha <= -1)
                throw ArgumentError("AssociatedLaguerreBasis: alpha must be finite and > -1");
        }

        Real Alpha() const { return _alpha; }

        // Evaluate associated Laguerre polynomial Lₙ^(α)(x)
        Real Evaluate(int n, Real x) const
        {
            if (n < 0)
                throw ArgumentError("AssociatedLaguerreBasis::Evaluate: n must be non-negative");
            if (!std::isfinite(x))
                throw DomainError("AssociatedLaguerreBasis::Evaluate: x must be finite");
            
            // Use recurrence relation:
            // (n+1)Lₙ₊₁^(α) = (2n + α + 1 - x)Lₙ^(α) - (n + α)Lₙ₋₁^(α)
            
            if (n == 0)
                return 1.0;
            if (n == 1)
                return MML::Real(1.0) + _alpha - x;
            
            Real Ln_minus_1 = 1.0;
            Real Ln = MML::Real(1.0) + _alpha - x;
            Real Ln_plus_1 = 0.0;
            
            for (int k = 1; k < n; k++)
            {
                Real k_real = static_cast<Real>(k);
                Ln_plus_1 = ((Real(2)*k_real + _alpha + Real(1) - x) * Ln
                             - (k_real + _alpha) * Ln_minus_1) / (k_real + MML::Real(1.0));
                Ln_minus_1 = Ln;
                Ln = Ln_plus_1;
            }
            
            return Ln_plus_1;
        }

        // Weight function w(x) = x^α e^(-x)
        Real WeightFunction(Real x) const 
        { 
            if (!std::isfinite(x) || x < 0)
                throw DomainError("AssociatedLaguerreBasis::WeightFunction: x must be finite and nonnegative");
            if (x == 0)
                return _alpha < 0 ? Constants::PosInf : _alpha == 0 ? Real{1} : Real{0};
            const Real value = std::exp(_alpha * std::log(x) - x);
            IntegrationDetail::RequireFinite(value, "Associated Laguerre weight");
            return value;
        }

        // Normalization: Γ(n+α+1)/n!
        Real Normalization(int n) const
        {
            if (n < 0)
                throw ArgumentError("AssociatedLaguerreBasis::Normalization: n must be non-negative");
            
            Real result = std::tgamma(_alpha + 1);
            IntegrationDetail::RequireFinite(result, "Associated Laguerre normalization");
            for (int k = 0; k < n; ++k)
            {
                result *= 1 + _alpha / (Real(k) + 1);
                IntegrationDetail::RequireFinite(result, "Associated Laguerre normalization");
            }
            return result;
        }

        // Domain bounds
        Real DomainMin() const { return 0.0; }
        Real DomainMax() const { return std::numeric_limits<Real>::infinity(); }

        // Recurrence coefficients
        void RecurrenceCoefficients(int n, Real& a, Real& b, Real& c) const
        {
            if (n < 0)
                throw ArgumentError("AssociatedLaguerreBasis::RecurrenceCoefficients: n must be non-negative");
            
            Real n_real = static_cast<Real>(n);
            Real n_plus_1 = n_real + MML::Real(1.0);
            
            a = -MML::Real(1.0) / n_plus_1;
            b = (MML::Real(2.0) * n_real + _alpha + MML::Real(1.0)) / n_plus_1;
            c = -(n_real + _alpha) / n_plus_1;
        }
    };

} // namespace MML

#endif // MML_LAGUERRE_BASIS_H

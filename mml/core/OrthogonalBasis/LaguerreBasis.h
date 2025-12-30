///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        LaguerreBasis.h                                                     ///
///  Description: Laguerre polynomial basis functions                                 ///
///               Orthogonal on [0, ∞) with weight w(x) = e^(-x)                     ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_LAGUERRE_BASIS_H
#define MML_LAGUERRE_BASIS_H

#include "MMLBase.h"

#include "base/StandardFunctions.h"
#include "interfaces/IFunction.h"
#include "core/OrthogonalBasis.h"

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
                throw std::invalid_argument("LaguerreBasis::Evaluate: n must be non-negative");
            
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
                throw std::invalid_argument("LaguerreBasis::Normalization: n must be non-negative");
            
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
                throw std::invalid_argument("LaguerreBasis::RecurrenceCoefficients: n must be non-negative");
            
            Real n_real = static_cast<Real>(n);
            Real n_plus_1 = n_real + 1.0;
            
            a = -1.0 / n_plus_1;
            b = (2.0 * n_real + 1.0) / n_plus_1;
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
            if (alpha <= -1.0)
                throw std::invalid_argument("AssociatedLaguerreBasis: alpha must be > -1");
        }

        Real Alpha() const { return _alpha; }

        // Evaluate associated Laguerre polynomial Lₙ^(α)(x)
        Real Evaluate(int n, Real x) const
        {
            if (n < 0)
                throw std::invalid_argument("AssociatedLaguerreBasis::Evaluate: n must be non-negative");
            
            // Use recurrence relation:
            // (n+1)Lₙ₊₁^(α) = (2n + α + 1 - x)Lₙ^(α) - (n + α)Lₙ₋₁^(α)
            
            if (n == 0)
                return 1.0;
            if (n == 1)
                return 1.0 + _alpha - x;
            
            Real Ln_minus_1 = 1.0;
            Real Ln = 1.0 + _alpha - x;
            Real Ln_plus_1 = 0.0;
            
            for (int k = 1; k < n; k++)
            {
                Real k_real = static_cast<Real>(k);
                Ln_plus_1 = ((2.0*k_real + _alpha + 1.0 - x) * Ln 
                             - (k_real + _alpha) * Ln_minus_1) / (k_real + 1.0);
                Ln_minus_1 = Ln;
                Ln = Ln_plus_1;
            }
            
            return Ln_plus_1;
        }

        // Weight function w(x) = x^α e^(-x)
        Real WeightFunction(Real x) const 
        { 
            return std::pow(x, _alpha) * std::exp(-x); 
        }

        // Normalization: Γ(n+α+1)/n!
        Real Normalization(int n) const
        {
            if (n < 0)
                throw std::invalid_argument("AssociatedLaguerreBasis::Normalization: n must be non-negative");
            
            // Γ(n+α+1)/n! = (n+α)(n+α-1)...(α+1)Γ(α+1)
            // For simplicity, use: (n+α)!/n! * Γ(α+1)/α!
            
            Real gamma_alpha_plus_1 = std::tgamma(_alpha + 1.0);
            Real n_factorial = 1.0;
            for (int k = 1; k <= n; k++)
                n_factorial *= k;
            
            Real result = gamma_alpha_plus_1;
            for (int k = 1; k <= n; k++)
                result *= (_alpha + k);
            result /= n_factorial;
            
            return result;
        }

        // Domain bounds
        Real DomainMin() const { return 0.0; }
        Real DomainMax() const { return std::numeric_limits<Real>::infinity(); }

        // Recurrence coefficients
        void RecurrenceCoefficients(int n, Real& a, Real& b, Real& c) const
        {
            if (n < 0)
                throw std::invalid_argument("AssociatedLaguerreBasis::RecurrenceCoefficients: n must be non-negative");
            
            Real n_real = static_cast<Real>(n);
            Real n_plus_1 = n_real + 1.0;
            
            a = -1.0 / n_plus_1;
            b = (2.0 * n_real + _alpha + 1.0) / n_plus_1;
            c = -(n_real + _alpha) / n_plus_1;
        }
    };

} // namespace MML

#endif // MML_LAGUERRE_BASIS_H

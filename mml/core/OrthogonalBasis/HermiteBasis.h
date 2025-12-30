///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        HermiteBasis.h                                                      ///
///  Description: Hermite polynomial basis functions                                  ///
///               Orthogonal on (-∞, ∞) with weight w(x) = e^(-x²)                    ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_HERMITE_BASIS_H
#define MML_HERMITE_BASIS_H

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
    // HermiteBasis - Hermite polynomial orthogonal basis (physicist's convention)
    //
    // Basis functions: Hₙ(x) - Hermite polynomials
    //
    // Orthogonality: ∫₋∞^∞ Hₘ(x)Hₙ(x)e^(-x²) dx = 0 for m ≠ n
    // Normalization: ∫₋∞^∞ Hₙ²(x)e^(-x²) dx = 2ⁿn!√π
    //
    // Domain: (-∞, ∞)
    // Weight function: w(x) = e^(-x²)
    //
    // Recurrence relation:
    //   H₀(x) = 1
    //   H₁(x) = 2x
    //   Hₙ₊₁(x) = 2xHₙ(x) - 2nHₙ₋₁(x)
    //
    // Physical applications:
    //   - Quantum harmonic oscillator wavefunctions
    //   - Gaussian beam modes in optics
    //   - Heat equation solutions
    //   - Signal processing (Gabor transform)
    //
    // Properties:
    //   - Hₙ(-x) = (-1)ⁿHₙ(x)  (parity)
    //   - d/dx Hₙ(x) = 2nHₙ₋₁(x)
    //   - Hₙ solves: y'' - 2xy' + 2ny = 0  (Hermite differential equation)
    ///////////////////////////////////////////////////////////////////////////////////////////
    class HermiteBasis : public OrthogonalBasis
    {
    public:
        HermiteBasis() = default;

        // Evaluate Hermite polynomial Hₙ(x) using physicist's convention
        Real Evaluate(int n, Real x) const override
        {
            if (n < 0)
                throw std::invalid_argument("HermiteBasis::Evaluate: n must be non-negative");
            
            return Functions::Hermite(static_cast<unsigned int>(n), x);
        }

        // Weight function w(x) = e^(-x²)
        Real WeightFunction(Real x) const override 
        { 
            return std::exp(-x * x); 
        }

        // Normalization: ||Hₙ||² = ∫₋∞^∞ Hₙ²(x)e^(-x²) dx = 2ⁿn!√π
        Real Normalization(int n) const override
        {
            if (n < 0)
                throw std::invalid_argument("HermiteBasis::Normalization: n must be non-negative");
            
            // 2^n * n! * sqrt(π)
            Real two_to_n = std::pow(2.0, static_cast<Real>(n));
            Real n_factorial = 1.0;
            for (int k = 1; k <= n; k++)
                n_factorial *= k;
            
            return two_to_n * n_factorial * std::sqrt(Constants::PI);
        }

        // Domain bounds (unbounded, but return practical limits)
        Real DomainMin() const override 
        { 
            return -std::numeric_limits<Real>::infinity(); 
        }
        
        Real DomainMax() const override 
        { 
            return std::numeric_limits<Real>::infinity(); 
        }

        // Recurrence coefficients for three-term recurrence:
        // Hₙ₊₁(x) = (aₙx + bₙ)Hₙ(x) + cₙHₙ₋₁(x)
        // 
        // Standard form: Hₙ₊₁ = 2xHₙ - 2nHₙ₋₁
        // So: aₙ = 2, bₙ = 0, cₙ = -2n
        void RecurrenceCoefficients(int n, Real& a, Real& b, Real& c) const
        {
            if (n < 0)
                throw std::invalid_argument("HermiteBasis::RecurrenceCoefficients: n must be non-negative");
            
            a = 2.0;
            b = 0.0;
            c = -2.0 * n;
        }

        // Quantum harmonic oscillator wavefunction ψₙ(x) = Nₙ Hₙ(x) e^(-x²/2)
        // where Nₙ = 1/√(2ⁿn!√π) is normalization constant
        Real QuantumWavefunction(int n, Real x) const
        {
            if (n < 0)
                throw std::invalid_argument("HermiteBasis::QuantumWavefunction: n must be non-negative");
            
            Real Hn = Evaluate(n, x);
            Real norm_factor = 1.0 / std::sqrt(Normalization(n));
            Real gaussian = std::exp(-0.5 * x * x);
            
            return norm_factor * Hn * gaussian;
        }
    };

    ///////////////////////////////////////////////////////////////////////////////////////////
    // Probabilist's Hermite polynomials Heₙ(x)
    //
    // Weight function: w(x) = e^(-x²/2)
    // Recurrence: Heₙ₊₁(x) = xHeₙ(x) - nHeₙ₋₁(x)
    // Normalization: ∫₋∞^∞ Heₙ²(x)e^(-x²/2) dx = n!√(2π)
    //
    // Relation to physicist's: Heₙ(x) = 2^(-n/2) Hₙ(x/√2)
    //
    // Used in probability theory and statistics
    ///////////////////////////////////////////////////////////////////////////////////////////
    class ProbabilistHermiteBasis
    {
    public:
        ProbabilistHermiteBasis() = default;

        // Evaluate probabilist's Hermite polynomial Heₙ(x)
        Real Evaluate(int n, Real x) const
        {
            if (n < 0)
                throw std::invalid_argument("ProbabilistHermiteBasis::Evaluate: n must be non-negative");
            
            // Convert: Heₙ(x) = 2^(-n/2) Hₙ(x/√2)
            Real x_scaled = x / Constants::SQRT2;
            Real Hn = Functions::Hermite(static_cast<unsigned int>(n), x_scaled);
            Real scale = std::pow(2.0, -0.5 * n);
            
            return scale * Hn;
        }

        // Weight function w(x) = e^(-x²/2)
        Real WeightFunction(Real x) const 
        { 
            return std::exp(-0.5 * x * x); 
        }

        // Normalization: ||Heₙ||² = n!√(2π)
        Real Normalization(int n) const
        {
            if (n < 0)
                throw std::invalid_argument("ProbabilistHermiteBasis::Normalization: n must be non-negative");
            
            Real n_factorial = 1.0;
            for (int k = 1; k <= n; k++)
                n_factorial *= k;
            
            return n_factorial * std::sqrt(2.0 * Constants::PI);
        }

        // Domain bounds
        Real DomainMin() const { return -std::numeric_limits<Real>::infinity(); }
        Real DomainMax() const { return std::numeric_limits<Real>::infinity(); }
    };

} // namespace MML

#endif // MML_HERMITE_BASIS_H

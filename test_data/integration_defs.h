#if !defined __MML_INTEGRATION_DEFS_H
#define __MML_INTEGRATION_DEFS_H

#include <string>
#include <cmath>
#include <limits>

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#endif

namespace MML::TestBeds
{
    ///////////////////////////////////////////////////////////////////////////////////////////
    //                         INTEGRATION TEST FUNCTIONS                                     //
    ///////////////////////////////////////////////////////////////////////////////////////////

    // ========================================================================================
    // CATEGORY 1: NORMAL CONTINUOUS FUNCTIONS (smooth, well-behaved)
    // ========================================================================================

    // 1. Polynomial: f(x) = x³ - 2x + 1
    // ∫[0,2] (x³ - 2x + 1) dx = [x⁴/4 - x² + x]₀² = (4 - 4 + 2) - 0 = 2
    static Real Int_Polynomial(Real x) { return x * x * x - 2 * x + 1; }
    static Real Int_Polynomial_antideriv(Real x) { return x * x * x * x / 4 - x * x + x; }
    const static inline Real Int_Polynomial_a = 0.0;
    const static inline Real Int_Polynomial_b = 2.0;
    const static inline Real Int_Polynomial_result = 2.0;

    // 2. Exponential: f(x) = e^(-x²/2) (Gaussian bell curve, not normalized)
    // ∫[-1,1] e^(-x²/2) dx ≈ 1.7112433... (no closed form, but known numerically)
    static Real Int_Gaussian(Real x) { return std::exp(-x * x / 2); }
    const static inline Real Int_Gaussian_a = -1.0;
    const static inline Real Int_Gaussian_b = 1.0;
    const static inline Real Int_Gaussian_result = 1.7112433961055917;  // erf(1/√2) * √(2π)

    // 3. Trigonometric: f(x) = sin(x) * cos(x) = sin(2x)/2
    // ∫[0,π] sin(x)cos(x) dx = [sin²(x)/2]₀^π = 0
    static Real Int_SinCos(Real x) { return std::sin(x) * std::cos(x); }
    static Real Int_SinCos_antideriv(Real x) { return std::sin(x) * std::sin(x) / 2; }
    const static inline Real Int_SinCos_a = 0.0;
    const static inline Real Int_SinCos_b = Constants::PI;
    const static inline Real Int_SinCos_result = 0.0;

    // 4. Rational: f(x) = 1/(1 + x²) (arctangent derivative)
    // ∫[0,1] 1/(1+x²) dx = [arctan(x)]₀¹ = π/4
    static Real Int_Arctan(Real x) { return 1.0 / (1 + x * x); }
    static Real Int_Arctan_antideriv(Real x) { return std::atan(x); }
    const static inline Real Int_Arctan_a = 0.0;
    const static inline Real Int_Arctan_b = 1.0;
    const static inline Real Int_Arctan_result = Constants::PI / 4;

    // 5. Mixed: f(x) = x * e^(-x)
    // ∫[0,∞) x*e^(-x) dx = 1 (but we use [0,5] ≈ 0.9596...)
    // ∫[0,3] x*e^(-x) dx = [-e^(-x)(x+1)]₀³ = 1 - 4e^(-3) ≈ 0.8008517...
    static Real Int_XExpNegX(Real x) { return x * std::exp(-x); }
    static Real Int_XExpNegX_antideriv(Real x) { return -std::exp(-x) * (x + 1); }
    const static inline Real Int_XExpNegX_a = 0.0;
    const static inline Real Int_XExpNegX_b = 3.0;
    const static inline Real Int_XExpNegX_result = 1.0 - 4.0 * std::exp(-3.0);  // ≈ 0.8008517265

    // ========================================================================================
    // CATEGORY 2: OSCILLATORY FUNCTIONS (require more integration points)
    // ========================================================================================

    // 6. High-frequency sine: f(x) = sin(10x)
    // ∫[0,π] sin(10x) dx = [-cos(10x)/10]₀^π = (-cos(10π) + cos(0))/10 = 0
    static Real Int_FastSine(Real x) { return std::sin(10 * x); }
    static Real Int_FastSine_antideriv(Real x) { return -std::cos(10 * x) / 10; }
    const static inline Real Int_FastSine_a = 0.0;
    const static inline Real Int_FastSine_b = Constants::PI;
    const static inline Real Int_FastSine_result = 0.0;

    // 7. Damped oscillation: f(x) = e^(-x) * sin(5x)
    // ∫[0,2π] e^(-x)*sin(5x) dx (no simple closed form)
    // Result computed via integration by parts: e^(-x)(-sin(5x) - 5cos(5x))/26
    static Real Int_DampedOsc(Real x) { return std::exp(-x) * std::sin(5 * x); }
    static Real Int_DampedOsc_antideriv(Real x) { 
        return std::exp(-x) * (-std::sin(5 * x) - 5 * std::cos(5 * x)) / 26; 
    }
    const static inline Real Int_DampedOsc_a = 0.0;
    const static inline Real Int_DampedOsc_b = 2.0 * Constants::PI;
    const static inline Real Int_DampedOsc_result = 
        (1.0 - std::exp(-2 * Constants::PI) * (-std::sin(10 * Constants::PI) - 5 * std::cos(10 * Constants::PI))) / 26
        - (-5.0 / 26);  // Simplified: 5/26 * (1 - e^(-2π))

    // 8. Bessel-like: f(x) = cos(x²) (Fresnel-type integral)
    // ∫[0,√(π/2)] cos(x²) dx = √(π/8) ≈ 0.6267...
    static Real Int_CosSq(Real x) { return std::cos(x * x); }
    const static inline Real Int_CosSq_a = 0.0;
    const static inline Real Int_CosSq_b = std::sqrt(Constants::PI / 2);
    const static inline Real Int_CosSq_result = std::sqrt(Constants::PI / 8);  // Fresnel C(1) scaled

    // 9. Product of sines: f(x) = sin(x) * sin(3x)
    // Using product formula: sin(x)sin(3x) = [cos(2x) - cos(4x)]/2
    // ∫[0,π] = [sin(2x)/4 - sin(4x)/8]₀^π = 0
    static Real Int_SinProduct(Real x) { return std::sin(x) * std::sin(3 * x); }
    const static inline Real Int_SinProduct_a = 0.0;
    const static inline Real Int_SinProduct_b = Constants::PI;
    const static inline Real Int_SinProduct_result = 0.0;

    // 10. Chirp: f(x) = sin(x²) (frequency increases with x)
    // ∫[0,√(2π)] sin(x²) dx = √(π/8) ≈ 0.6267 (Fresnel S integral)
    static Real Int_Chirp(Real x) { return std::sin(x * x); }
    const static inline Real Int_Chirp_a = 0.0;
    const static inline Real Int_Chirp_b = std::sqrt(2 * Constants::PI);
    const static inline Real Int_Chirp_result = std::sqrt(Constants::PI / 8);  // Fresnel S(√(2/π))

    // ========================================================================================
    // CATEGORY 3: IMPROPER INTEGRALS (infinite limits or singularities)
    // ========================================================================================

    // 11. Integrable singularity at origin: f(x) = 1/√x
    // ∫[0,1] x^(-1/2) dx = [2√x]₀¹ = 2
    static Real Int_SqrtSing(Real x) { 
        if (x <= 0) return std::numeric_limits<Real>::infinity();
        return 1.0 / std::sqrt(x); 
    }
    static Real Int_SqrtSing_antideriv(Real x) { return 2 * std::sqrt(x); }
    const static inline Real Int_SqrtSing_a = 0.0;
    const static inline Real Int_SqrtSing_b = 1.0;
    const static inline Real Int_SqrtSing_result = 2.0;
    const static inline bool Int_SqrtSing_improper = true;
    const static inline Real Int_SqrtSing_singularity = 0.0;

    // 12. Log singularity: f(x) = -ln(x) (integrable at x=0)
    // ∫[0,1] -ln(x) dx = [x - x*ln(x)]₀¹ = 1 (using L'Hôpital at 0)
    static Real Int_LogSing(Real x) { 
        if (x <= 0) return std::numeric_limits<Real>::infinity();
        return -std::log(x); 
    }
    const static inline Real Int_LogSing_a = 0.0;
    const static inline Real Int_LogSing_b = 1.0;
    const static inline Real Int_LogSing_result = 1.0;
    const static inline bool Int_LogSing_improper = true;
    const static inline Real Int_LogSing_singularity = 0.0;

    // 13. Semi-infinite: f(x) = e^(-x) on [0,∞)
    // ∫[0,∞) e^(-x) dx = 1
    // For finite computation, ∫[0,10] e^(-x) dx ≈ 1 - e^(-10) ≈ 0.9999546...
    static Real Int_ExpDecay(Real x) { return std::exp(-x); }
    static Real Int_ExpDecay_antideriv(Real x) { return -std::exp(-x); }
    const static inline Real Int_ExpDecay_a = 0.0;
    const static inline Real Int_ExpDecay_b_inf = std::numeric_limits<Real>::infinity();
    const static inline Real Int_ExpDecay_b_finite = 20.0;  // Practical cutoff
    const static inline Real Int_ExpDecay_result_inf = 1.0;
    const static inline Real Int_ExpDecay_result_finite = 1.0 - std::exp(-20.0);
    const static inline bool Int_ExpDecay_improper = true;

    // 14. Gaussian tail: f(x) = e^(-x²) on [0,∞)
    // ∫[0,∞) e^(-x²) dx = √π/2 ≈ 0.8862269...
    static Real Int_GaussTail(Real x) { return std::exp(-x * x); }
    const static inline Real Int_GaussTail_a = 0.0;
    const static inline Real Int_GaussTail_b_inf = std::numeric_limits<Real>::infinity();
    const static inline Real Int_GaussTail_b_finite = 6.0;  // Practical cutoff (e^(-36) ≈ 0)
    const static inline Real Int_GaussTail_result = std::sqrt(Constants::PI) / 2;
    const static inline bool Int_GaussTail_improper = true;

    // 15. Algebraic decay: f(x) = 1/(1+x²) on [0,∞)
    // ∫[0,∞) 1/(1+x²) dx = π/2
    static Real Int_AlgDecay(Real x) { return 1.0 / (1 + x * x); }
    static Real Int_AlgDecay_antideriv(Real x) { return std::atan(x); }
    const static inline Real Int_AlgDecay_a = 0.0;
    const static inline Real Int_AlgDecay_b_inf = std::numeric_limits<Real>::infinity();
    const static inline Real Int_AlgDecay_b_finite = 1000.0;  // Practical cutoff
    const static inline Real Int_AlgDecay_result = Constants::PI / 2;
    const static inline bool Int_AlgDecay_improper = true;

    // ========================================================================================
    // Helper constants for known results
    // ========================================================================================
    
    // Fresnel integrals at specific points
    const static inline Real FRESNEL_C_1 = 0.7798934003768228;  // ∫[0,1] cos(πt²/2) dt
    const static inline Real FRESNEL_S_1 = 0.4382591473903548;  // ∫[0,1] sin(πt²/2) dt

} // namespace MML::TestBeds

#endif // __MML_INTEGRATION_DEFS_H

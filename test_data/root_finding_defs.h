#if !defined __MML_ROOT_FINDING_DEFS_H
#define __MML_ROOT_FINDING_DEFS_H

#include <string>
#include <cmath>
#include <vector>

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#endif

namespace MML::TestBeds
{
    ///////////////////////////////////////////////////////////////////////////////////////////
    //                    ROOT FINDING TEST FUNCTIONS AND THEIR DERIVATIVES                   //
    ///////////////////////////////////////////////////////////////////////////////////////////

    // ========================================================================================
    // CATEGORY 1: POLYNOMIAL FUNCTIONS (exact roots, well-conditioned)
    // ========================================================================================

    // 1. Simple quadratic: f(x) = x² - 2  (roots at ±√2)
    static Real Poly_xSq_minus_2(Real x) { return x * x - 2; }
    static Real Poly_xSq_minus_2_deriv(Real x) { return 2 * x; }
    // Known roots: √2 ≈ 1.41421356237, -√2 ≈ -1.41421356237

    // 2. Cubic with known factors: f(x) = x³ - 6x² + 11x - 6 = (x-1)(x-2)(x-3)
    static Real Poly_cubic_123(Real x) { return x*x*x - 6*x*x + 11*x - 6; }
    static Real Poly_cubic_123_deriv(Real x) { return 3*x*x - 12*x + 11; }
    // Known roots: 1, 2, 3

    // 3. Quadratic with integer roots: f(x) = x² - 5x + 6 = (x-2)(x-3)
    static Real Poly_quad_23(Real x) { return x*x - 5*x + 6; }
    static Real Poly_quad_23_deriv(Real x) { return 2*x - 5; }
    // Known roots: 2, 3

    // 4. Quartic: f(x) = x⁴ - 5x² + 4 = (x-1)(x+1)(x-2)(x+2)
    static Real Poly_quartic_1122(Real x) { return x*x*x*x - 5*x*x + 4; }
    static Real Poly_quartic_1122_deriv(Real x) { return 4*x*x*x - 10*x; }
    // Known roots: -2, -1, 1, 2

    // 5. Wilkinson-style polynomial: f(x) = (x-1)(x-2)...(x-5) = x⁵ - 15x⁴ + 85x³ - 225x² + 274x - 120
    static Real Poly_wilkinson5(Real x) {
        return x*x*x*x*x - 15*x*x*x*x + 85*x*x*x - 225*x*x + 274*x - 120;
    }
    static Real Poly_wilkinson5_deriv(Real x) {
        return 5*x*x*x*x - 60*x*x*x + 255*x*x - 450*x + 274;
    }
    // Known roots: 1, 2, 3, 4, 5

    // ========================================================================================
    // CATEGORY 2: TRANSCENDENTAL EQUATIONS (classical benchmarks)
    // ========================================================================================

    // 6. Fixed point: f(x) = x - cos(x)  (root where x = cos(x))
    static Real Trans_x_minus_cosx(Real x) { return x - std::cos(x); }
    static Real Trans_x_minus_cosx_deriv(Real x) { return 1 + std::sin(x); }
    // Known root: ≈ 0.7390851332151607 (Dottie number)

    // 7. Kepler's equation: f(x) = x - e*sin(x) - M  (here e=0.5, M=0.5)
    static Real Trans_kepler(Real x) { return x - 0.5 * std::sin(x) - 0.5; }
    static Real Trans_kepler_deriv(Real x) { return 1 - 0.5 * std::cos(x); }
    // Known root: ≈ 0.7581847193 (eccentric anomaly for e=0.5, M=0.5)

    // 8. Exponential decay: f(x) = exp(-x) - x  (intersection of exp(-x) and x)
    static Real Trans_exp_minus_x(Real x) { return std::exp(-x) - x; }
    static Real Trans_exp_minus_x_deriv(Real x) { return -std::exp(-x) - 1; }
    // Known root: ≈ 0.5671432904097839 (Omega constant)

    // 9. Log equation: f(x) = ln(x) - 1  (root at x = e)
    static Real Trans_lnx_minus_1(Real x) { return std::log(x) - 1; }
    static Real Trans_lnx_minus_1_deriv(Real x) { return 1.0 / x; }
    // Known root: e ≈ 2.718281828

    // 10. Trigonometric: f(x) = sin(x) - 0.5
    static Real Trans_sinx_minus_half(Real x) { return std::sin(x) - 0.5; }
    static Real Trans_sinx_minus_half_deriv(Real x) { return std::cos(x); }
    // Known roots: π/6 ≈ 0.5236, 5π/6 ≈ 2.618, etc.

    // 11. Tangent intersection: f(x) = tan(x) - x  (near x=4.5)
    static Real Trans_tanx_minus_x(Real x) { return std::tan(x) - x; }
    static Real Trans_tanx_minus_x_deriv(Real x) { 
        Real c = std::cos(x);
        return 1.0 / (c * c) - 1; 
    }
    // Known roots: 0 (trivial), ≈4.4934, ≈7.7253, etc.

    // ========================================================================================
    // CATEGORY 3: MULTIPLE/REPEATED ROOTS (challenging for most methods)
    // ========================================================================================

    // 12. Double root: f(x) = (x-1)² = x² - 2x + 1
    static Real Multi_double_root_1(Real x) { return (x - 1) * (x - 1); }
    static Real Multi_double_root_1_deriv(Real x) { return 2 * (x - 1); }
    // Known root: 1 (multiplicity 2)

    // 13. Triple root: f(x) = (x-2)³ = x³ - 6x² + 12x - 8
    static Real Multi_triple_root_2(Real x) { return (x - 2) * (x - 2) * (x - 2); }
    static Real Multi_triple_root_2_deriv(Real x) { return 3 * (x - 2) * (x - 2); }
    // Known root: 2 (multiplicity 3)

    // 14. Mixed multiplicities: f(x) = (x-1)²(x-3) = x³ - 5x² + 7x - 3
    static Real Multi_mixed_1_3(Real x) { return (x - 1) * (x - 1) * (x - 3); }
    static Real Multi_mixed_1_3_deriv(Real x) { return 2 * (x - 1) * (x - 3) + (x - 1) * (x - 1); }
    // Known roots: 1 (multiplicity 2), 3 (simple)

    // ========================================================================================
    // CATEGORY 4: CLOSELY SPACED ROOTS (ill-conditioned)
    // ========================================================================================

    // 15. Nearly coincident: f(x) = (x - 1)(x - 1.001) = x² - 2.001x + 1.001
    static Real Close_roots_1_1001(Real x) { return (x - 1) * (x - 1.001); }
    static Real Close_roots_1_1001_deriv(Real x) { return 2*x - 2.001; }
    // Known roots: 1.0, 1.001

    // 16. Very close: f(x) = (x - 2)(x - 2.0001)
    static Real Close_roots_2_20001(Real x) { return (x - 2) * (x - 2.0001); }
    static Real Close_roots_2_20001_deriv(Real x) { return 2*x - 4.0001; }
    // Known roots: 2.0, 2.0001

    // ========================================================================================
    // CATEGORY 5: PATHOLOGICAL/CHALLENGING CASES
    // ========================================================================================

    // 17. Steep exponential: f(x) = exp(x) - 10000 (root at ln(10000) ≈ 9.21)
    static Real Patho_steep_exp(Real x) { return std::exp(x) - 10000; }
    static Real Patho_steep_exp_deriv(Real x) { return std::exp(x); }
    // Known root: ln(10000) ≈ 9.210340372

    // 18. Flat near root: f(x) = (x-1)³ - 0.001 (slow convergence)
    static Real Patho_flat_cubic(Real x) { return (x - 1) * (x - 1) * (x - 1) - 0.001; }
    static Real Patho_flat_cubic_deriv(Real x) { return 3 * (x - 1) * (x - 1); }
    // Known root: 1 + ∛0.001 ≈ 1.1

    // 19. Oscillatory decay: f(x) = sin(10x) * exp(-x)  (many roots, decay)
    static Real Patho_oscill_decay(Real x) { return std::sin(10 * x) * std::exp(-x); }
    static Real Patho_oscill_decay_deriv(Real x) { 
        return std::exp(-x) * (10 * std::cos(10 * x) - std::sin(10 * x)); 
    }
    // Known roots: kπ/10 for k = 0, 1, 2, ...

    // 20. Large derivative: f(x) = atan(1000*(x-1))  (very steep at x=1)
    static Real Patho_steep_atan(Real x) { return std::atan(1000 * (x - 1)); }
    static Real Patho_steep_atan_deriv(Real x) { 
        Real arg = 1000 * (x - 1);
        return 1000 / (1 + arg * arg); 
    }
    // Known root: 1 (but extremely steep)

    // ========================================================================================
    // CATEGORY 6: FUNCTIONS WITH SINGULARITIES NEAR ROOTS
    // ========================================================================================

    // 21. Near singularity: f(x) = 1/x - 2 (root at x=0.5, singularity at x=0)
    static Real Sing_1_over_x_minus_2(Real x) { return 1.0 / x - 2; }
    static Real Sing_1_over_x_minus_2_deriv(Real x) { return -1.0 / (x * x); }
    // Known root: 0.5

    // 22. Log singularity: f(x) = ln(x) + x - 2 (root near x≈1.56, singularity at x=0)
    static Real Sing_lnx_plus_x_minus_2(Real x) { return std::log(x) + x - 2; }
    static Real Sing_lnx_plus_x_minus_2_deriv(Real x) { return 1.0 / x + 1; }
    // Known root: ≈ 1.5571455990

    // ========================================================================================
    // CATEGORY 7: PHYSICAL/ENGINEERING APPLICATIONS
    // ========================================================================================

    // 23. Van der Waals equation (simplified): f(V) = V³ - 2V² + V - 0.1
    // (Looking for molar volume at certain P,T)
    static Real Phys_vdw(Real V) { return V*V*V - 2*V*V + V - 0.1; }
    static Real Phys_vdw_deriv(Real V) { return 3*V*V - 4*V + 1; }
    // Root approximately 0.11

    // 24. Planck radiation (simplified): f(x) = x - 5*(1 - exp(-x))
    // (Wien displacement law crossover)
    static Real Phys_planck(Real x) { return x - 5 * (1 - std::exp(-x)); }
    static Real Phys_planck_deriv(Real x) { return 1 - 5 * std::exp(-x); }
    // Known root: ≈ 4.9651

    // 25. Lambert W related: f(x) = x*exp(x) - 1 (W(1) problem)
    static Real Phys_lambert(Real x) { return x * std::exp(x) - 1; }
    static Real Phys_lambert_deriv(Real x) { return (1 + x) * std::exp(x); }
    // Known root: W(1) ≈ 0.5671432904

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                              KNOWN ROOT VALUES                                         //
    ///////////////////////////////////////////////////////////////////////////////////////////

    // High precision root values for verification
    const static inline Real ROOT_SQRT2 = static_cast<Real>(1.41421356237309504880168872420969807856967187537694);
    const static inline Real ROOT_DOTTIE = static_cast<Real>(0.73908513321516064165531208767387340401341175890076);
    const static inline Real ROOT_OMEGA = static_cast<Real>(0.56714329040978387299996866221035554975381578718651);
    const static inline Real ROOT_E = static_cast<Real>(2.71828182845904523536028747135266249775724709369995);
    const static inline Real ROOT_PI_6 = static_cast<Real>(0.52359877559829887307710723054658381403286156656252);
    const static inline Real ROOT_TAN_1 = static_cast<Real>(4.49340945790906417535833613264473925167866328789787);
    const static inline Real ROOT_LAMBERT_W1 = static_cast<Real>(0.56714329040978387299996866221035554975381578718651);
    const static inline Real ROOT_WIEN = static_cast<Real>(4.96511423174427630226686427576612621889222458284869);
    const static inline Real ROOT_KEPLER_05_05 = static_cast<Real>(0.75818471929987481918927862261787169523315279313234);
    const static inline Real ROOT_LN10000 = static_cast<Real>(9.21034037197618203544374908830428643115659825936887);

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                         PRE-DEFINED BRACKET INTERVALS                                  //
    ///////////////////////////////////////////////////////////////////////////////////////////

    // Brackets for polynomial roots
    const static inline Real bracket_sqrt2_low = 1.0;
    const static inline Real bracket_sqrt2_high = 2.0;

    const static inline Real bracket_cubic123_root1_low = 0.5;
    const static inline Real bracket_cubic123_root1_high = 1.5;
    const static inline Real bracket_cubic123_root2_low = 1.5;
    const static inline Real bracket_cubic123_root2_high = 2.5;
    const static inline Real bracket_cubic123_root3_low = 2.5;
    const static inline Real bracket_cubic123_root3_high = 3.5;

    // Brackets for transcendental roots
    const static inline Real bracket_dottie_low = 0.0;
    const static inline Real bracket_dottie_high = 1.0;

    const static inline Real bracket_omega_low = 0.0;
    const static inline Real bracket_omega_high = 1.0;

    const static inline Real bracket_e_low = 2.0;
    const static inline Real bracket_e_high = 3.0;

    const static inline Real bracket_sinpi6_low = 0.0;
    const static inline Real bracket_sinpi6_high = 1.0;

    // Brackets for challenging cases
    const static inline Real bracket_double_root_low = 0.0;
    const static inline Real bracket_double_root_high = 2.0;

    const static inline Real bracket_triple_root_low = 0.0;
    const static inline Real bracket_triple_root_high = 4.0;

} // namespace MML::TestBeds

#endif // __MML_ROOT_FINDING_DEFS_H

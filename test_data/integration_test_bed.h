#if !defined __MML_INTEGRATION_TEST_BED_H
#define __MML_INTEGRATION_TEST_BED_H

#include <string>
#include <vector>
#include <functional>
#include <cmath>
#include <limits>

#include "integration_defs.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#include "interfaces/IFunction.h"
#endif

namespace MML::TestBeds
{
    ///////////////////////////////////////////////////////////////////////////////////////////
    //                       INTEGRATION TEST DATA STRUCTURES                                 //
    ///////////////////////////////////////////////////////////////////////////////////////////

    /**
     * @brief Type of integral - affects choice of integration method
     */
    enum class IntegralType
    {
        Regular,            ///< Normal definite integral, no singularities
        Oscillatory,        ///< Highly oscillatory integrand
        SemiInfinite,       ///< One limit is ±∞
        DoublyInfinite,     ///< Both limits are ±∞
        Singular,           ///< Integrand has integrable singularity
        SingularEndpoint    ///< Singularity at integration endpoint
    };

    /**
     * @brief Represents a test case for numerical integration algorithms
     */
    struct TestIntegral
    {
        std::string name;                               ///< Descriptive name
        std::function<Real(Real)> func;                ///< The integrand f(x)
        Real a;                                         ///< Lower limit
        Real b;                                         ///< Upper limit (may be infinity)
        Real exactResult;                               ///< Analytically known result
        
        // Optional components
        std::function<Real(Real)> antiderivative;      ///< F(x) where F'(x) = f(x), if known
        bool hasAntiderivative = false;                 ///< Whether antiderivative is available
        
        // Integral characteristics
        IntegralType type = IntegralType::Regular;      ///< Classification
        bool isImproper = false;                        ///< Has infinite limit or singularity
        Real singularityLocation = 0;                   ///< Location of singularity if present
        
        // Metadata
        std::string category;                           ///< smooth, oscillatory, improper, etc.
        std::string description;                        ///< Detailed description
        int difficulty = 1;                             ///< 1=easy, 2=medium, 3=hard, 4=very hard
        Real suggestedTolerance = 1e-8;                ///< Reasonable tolerance for this integral
    };

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                            WRAPPER CLASS FOR MML INTERFACE                             //
    ///////////////////////////////////////////////////////////////////////////////////////////

    /**
     * @brief Wrapper to adapt std::function to MML's IRealFunction interface
     */
    class IntegrandWrapper : public IRealFunction
    {
    private:
        std::function<Real(Real)> _func;
    public:
        IntegrandWrapper(std::function<Real(Real)> f) : _func(f) {}
        Real operator()(Real x) const override { return _func(x); }
    };

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                    CATEGORY 1: SMOOTH CONTINUOUS FUNCTIONS                             //
    ///////////////////////////////////////////////////////////////////////////////////////////

    inline TestIntegral getPolynomialIntegral()
    {
        return TestIntegral{
            "Polynomial x³ - 2x + 1",
            Int_Polynomial,
            Int_Polynomial_a,
            Int_Polynomial_b,
            Int_Polynomial_result,
            Int_Polynomial_antideriv,
            true,
            IntegralType::Regular,
            false, 0,
            "smooth",
            "Simple polynomial integral. Exact result is 2. Tests basic quadrature accuracy.",
            1, 1e-10
        };
    }

    inline TestIntegral getGaussianIntegral()
    {
        return TestIntegral{
            "Gaussian exp(-x²/2) on [-1,1]",
            Int_Gaussian,
            Int_Gaussian_a,
            Int_Gaussian_b,
            Int_Gaussian_result,
            nullptr, false,
            IntegralType::Regular,
            false, 0,
            "smooth",
            "Gaussian bell curve on finite interval. No closed-form antiderivative. "
            "Related to error function.",
            1, 1e-8
        };
    }

    inline TestIntegral getSinCosIntegral()
    {
        return TestIntegral{
            "sin(x)cos(x) on [0,π]",
            Int_SinCos,
            Int_SinCos_a,
            Int_SinCos_b,
            Int_SinCos_result,
            Int_SinCos_antideriv,
            true,
            IntegralType::Regular,
            false, 0,
            "smooth",
            "Product of sine and cosine. Result is exactly 0 due to symmetry. "
            "Good test for cancellation accuracy.",
            1, 1e-12
        };
    }

    inline TestIntegral getArctanIntegral()
    {
        return TestIntegral{
            "1/(1+x²) on [0,1]",
            Int_Arctan,
            Int_Arctan_a,
            Int_Arctan_b,
            Int_Arctan_result,
            Int_Arctan_antideriv,
            true,
            IntegralType::Regular,
            false, 0,
            "smooth",
            "Arctangent integral. Result is π/4 ≈ 0.7854. Classic test for rational functions.",
            1, 1e-10
        };
    }

    inline TestIntegral getXExpNegXIntegral()
    {
        return TestIntegral{
            "x·exp(-x) on [0,3]",
            Int_XExpNegX,
            Int_XExpNegX_a,
            Int_XExpNegX_b,
            Int_XExpNegX_result,
            Int_XExpNegX_antideriv,
            true,
            IntegralType::Regular,
            false, 0,
            "smooth",
            "Exponential decay modulated by x. Tests integration by parts equivalent. "
            "Result ≈ 0.8009.",
            1, 1e-8
        };
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                    CATEGORY 2: OSCILLATORY FUNCTIONS                                   //
    ///////////////////////////////////////////////////////////////////////////////////////////

    inline TestIntegral getFastSineIntegral()
    {
        return TestIntegral{
            "sin(10x) on [0,π]",
            Int_FastSine,
            Int_FastSine_a,
            Int_FastSine_b,
            Int_FastSine_result,
            Int_FastSine_antideriv,
            true,
            IntegralType::Oscillatory,
            false, 0,
            "oscillatory",
            "High-frequency sine. Result is exactly 0. Tests handling of oscillatory "
            "integrands. Requires more quadrature points than smooth functions.",
            2, 1e-8
        };
    }

    inline TestIntegral getDampedOscIntegral()
    {
        return TestIntegral{
            "exp(-x)·sin(5x) on [0,2π]",
            Int_DampedOsc,
            Int_DampedOsc_a,
            Int_DampedOsc_b,
            Int_DampedOsc_result,
            Int_DampedOsc_antideriv,
            true,
            IntegralType::Oscillatory,
            false, 0,
            "oscillatory",
            "Damped oscillation. Combines exponential decay with oscillation. "
            "Result computed via integration by parts formula.",
            2, 1e-6
        };
    }

    inline TestIntegral getFresnelCosIntegral()
    {
        return TestIntegral{
            "cos(x²) Fresnel-type",
            Int_CosSq,
            Int_CosSq_a,
            Int_CosSq_b,
            Int_CosSq_result,
            nullptr, false,
            IntegralType::Oscillatory,
            false, 0,
            "oscillatory",
            "Fresnel cosine integral. No closed-form antiderivative. "
            "Frequency increases with x (chirp behavior). Result ≈ 0.627.",
            3, 1e-6
        };
    }

    inline TestIntegral getSinProductIntegral()
    {
        return TestIntegral{
            "sin(x)·sin(3x) on [0,π]",
            Int_SinProduct,
            Int_SinProduct_a,
            Int_SinProduct_b,
            Int_SinProduct_result,
            nullptr, false,
            IntegralType::Oscillatory,
            false, 0,
            "oscillatory",
            "Product of sines with different frequencies. Result is 0 due to orthogonality "
            "of Fourier modes.",
            2, 1e-10
        };
    }

    inline TestIntegral getChirpIntegral()
    {
        return TestIntegral{
            "sin(x²) Fresnel-type",
            Int_Chirp,
            Int_Chirp_a,
            Int_Chirp_b,
            Int_Chirp_result,
            nullptr, false,
            IntegralType::Oscillatory,
            false, 0,
            "oscillatory",
            "Fresnel sine integral. Chirp signal where frequency grows with x. "
            "No closed-form antiderivative. Result ≈ 0.627.",
            3, 1e-6
        };
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                    CATEGORY 3: IMPROPER INTEGRALS                                      //
    ///////////////////////////////////////////////////////////////////////////////////////////

    inline TestIntegral getSqrtSingularityIntegral()
    {
        return TestIntegral{
            "1/√x on [0,1] (singular)",
            Int_SqrtSing,
            Int_SqrtSing_a,
            Int_SqrtSing_b,
            Int_SqrtSing_result,
            Int_SqrtSing_antideriv,
            true,
            IntegralType::SingularEndpoint,
            true, Int_SqrtSing_singularity,
            "improper",
            "Square root singularity at x=0. Integrand → ∞ as x → 0, but integral "
            "converges to 2. Classic improper integral test.",
            3, 1e-6
        };
    }

    inline TestIntegral getLogSingularityIntegral()
    {
        return TestIntegral{
            "-ln(x) on [0,1] (singular)",
            Int_LogSing,
            Int_LogSing_a,
            Int_LogSing_b,
            Int_LogSing_result,
            nullptr, false,
            IntegralType::SingularEndpoint,
            true, Int_LogSing_singularity,
            "improper",
            "Logarithmic singularity at x=0. Function → ∞ as x → 0⁺, but "
            "integral converges to 1. Weaker singularity than 1/√x.",
            3, 1e-6
        };
    }

    inline TestIntegral getExpDecayInfiniteIntegral()
    {
        return TestIntegral{
            "exp(-x) on [0,∞)",
            Int_ExpDecay,
            Int_ExpDecay_a,
            Int_ExpDecay_b_finite,  // Use finite cutoff for practical computation
            Int_ExpDecay_result_finite,
            Int_ExpDecay_antideriv,
            true,
            IntegralType::SemiInfinite,
            true, 0,
            "improper",
            "Exponential decay to infinity. True result is 1. Using cutoff at x=20 "
            "gives result accurate to ~10⁻⁹.",
            2, 1e-8
        };
    }

    inline TestIntegral getGaussianTailIntegral()
    {
        return TestIntegral{
            "exp(-x²) on [0,∞)",
            Int_GaussTail,
            Int_GaussTail_a,
            Int_GaussTail_b_finite,
            Int_GaussTail_result,
            nullptr, false,
            IntegralType::SemiInfinite,
            true, 0,
            "improper",
            "Gaussian tail integral. Result is √π/2 ≈ 0.8862. Decays very fast, "
            "so finite cutoff at x=6 is sufficient for high accuracy.",
            2, 1e-8
        };
    }

    inline TestIntegral getAlgebraicDecayIntegral()
    {
        return TestIntegral{
            "1/(1+x²) on [0,∞)",
            Int_AlgDecay,
            Int_AlgDecay_a,
            Int_AlgDecay_b_finite,
            Int_AlgDecay_result,
            Int_AlgDecay_antideriv,
            true,
            IntegralType::SemiInfinite,
            true, 0,
            "improper",
            "Algebraic (slow) decay. Result is π/2. Requires larger cutoff than "
            "exponential decay due to 1/x² falloff.",
            3, 1e-4  // Lower accuracy due to slow decay
        };
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                          TEST COLLECTION GENERATORS                                    //
    ///////////////////////////////////////////////////////////////////////////////////////////

    /**
     * @brief Get all smooth/regular integration test cases
     */
    inline std::vector<TestIntegral> getSmoothIntegrals()
    {
        return {
            getPolynomialIntegral(),
            getGaussianIntegral(),
            getSinCosIntegral(),
            getArctanIntegral(),
            getXExpNegXIntegral()
        };
    }

    /**
     * @brief Get all oscillatory integration test cases
     */
    inline std::vector<TestIntegral> getOscillatoryIntegrals()
    {
        return {
            getFastSineIntegral(),
            getDampedOscIntegral(),
            getFresnelCosIntegral(),
            getSinProductIntegral(),
            getChirpIntegral()
        };
    }

    /**
     * @brief Get all improper integration test cases
     */
    inline std::vector<TestIntegral> getImproperIntegrals()
    {
        return {
            getSqrtSingularityIntegral(),
            getLogSingularityIntegral(),
            getExpDecayInfiniteIntegral(),
            getGaussianTailIntegral(),
            getAlgebraicDecayIntegral()
        };
    }

    /**
     * @brief Get all integration test cases
     */
    inline std::vector<TestIntegral> getAllIntegrals()
    {
        std::vector<TestIntegral> all;
        auto smooth = getSmoothIntegrals();
        auto oscillatory = getOscillatoryIntegrals();
        auto improper = getImproperIntegrals();
        
        all.insert(all.end(), smooth.begin(), smooth.end());
        all.insert(all.end(), oscillatory.begin(), oscillatory.end());
        all.insert(all.end(), improper.begin(), improper.end());
        
        return all;
    }

    /**
     * @brief Get test cases by difficulty level
     * @param maxDifficulty Include tests with difficulty <= this value
     */
    inline std::vector<TestIntegral> getIntegralsByDifficulty(int maxDifficulty)
    {
        std::vector<TestIntegral> result;
        auto all = getAllIntegrals();
        for (const auto& test : all)
        {
            if (test.difficulty <= maxDifficulty)
                result.push_back(test);
        }
        return result;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                          UTILITY FUNCTIONS                                             //
    ///////////////////////////////////////////////////////////////////////////////////////////

    /**
     * @brief Create MML function wrapper from test integral
     */
    inline IntegrandWrapper makeIntegrand(const TestIntegral& test)
    {
        return IntegrandWrapper(test.func);
    }

    /**
     * @brief Check if computed result matches expected within tolerance
     */
    inline bool checkResult(Real computed, Real expected, Real tolerance)
    {
        if (std::abs(expected) < 1e-14)
        {
            // For expected = 0, use absolute tolerance
            return std::abs(computed) < tolerance;
        }
        // For non-zero expected, use relative tolerance
        return std::abs(computed - expected) < tolerance * std::abs(expected);
    }

} // namespace MML::TestBeds

#endif // __MML_INTEGRATION_TEST_BED_H

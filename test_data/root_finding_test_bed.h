#if !defined __MML_ROOT_FINDING_TEST_BED_H
#define __MML_ROOT_FINDING_TEST_BED_H

#include <string>
#include <vector>
#include <functional>
#include <cmath>

#include "root_finding_defs.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#include "interfaces/IFunction.h"
#endif

namespace MML::TestBeds
{
    ///////////////////////////////////////////////////////////////////////////////////////////
    //                       TEST ROOT FUNCTION DATA STRUCTURE                                //
    ///////////////////////////////////////////////////////////////////////////////////////////

    /**
     * @brief Represents a test case for root finding algorithms
     * 
     * Contains the function, its derivative, known roots, bracketing intervals,
     * and metadata for categorization and difficulty assessment.
     */
    struct TestRootFunction
    {
        std::string name;                                       ///< Descriptive name for the test
        std::function<Real(Real)> func;                        ///< The function f(x)
        std::function<Real(Real)> derivative;                  ///< f'(x) for Newton's method
        std::vector<Real> knownRoots;                          ///< Analytically known root values
        std::vector<int> multiplicities;                       ///< Multiplicity of each root (1=simple, 2=double, etc.)
        std::vector<std::pair<Real, Real>> brackets;           ///< [low, high] bracketing intervals for each root
        
        // Metadata
        std::string category;                                  ///< polynomial, transcendental, multiple_roots, etc.
        std::string description;                               ///< Detailed description of the test case
        bool hasDerivative = true;                             ///< Whether derivative is available
        bool hasSingularity = false;                           ///< Whether function has singularities
        Real singularityLocation = 0;                          ///< Location of singularity if present
        int difficulty = 1;                                    ///< 1=easy, 2=medium, 3=hard, 4=challenging
    };

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                            WRAPPER CLASS FOR MML INTERFACE                             //
    ///////////////////////////////////////////////////////////////////////////////////////////

    /**
     * @brief Wrapper to adapt std::function to MML's IRealFunction interface
     */
    class RealFunctionWrapper : public IRealFunction
    {
    private:
        std::function<Real(Real)> _func;
    public:
        RealFunctionWrapper(std::function<Real(Real)> f) : _func(f) {}
        Real operator()(Real x) const override { return _func(x); }
    };

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                              TEST CASE GENERATORS                                       //
    ///////////////////////////////////////////////////////////////////////////////////////////

    // ========================================================================================
    // POLYNOMIAL TEST CASES
    // ========================================================================================

    inline TestRootFunction getQuadraticSqrt2Test()
    {
        return TestRootFunction{
            "x² - 2 (sqrt(2) root)",
            Poly_xSq_minus_2,
            Poly_xSq_minus_2_deriv,
            { ROOT_SQRT2, -ROOT_SQRT2 },
            { 1, 1 },
            { {bracket_sqrt2_low, bracket_sqrt2_high}, {-bracket_sqrt2_high, -bracket_sqrt2_low} },
            "polynomial",
            "Simple quadratic with irrational roots at ±√2. Classic benchmark.",
            true, false, 0, 1
        };
    }

    inline TestRootFunction getCubic123Test()
    {
        return TestRootFunction{
            "x³ - 6x² + 11x - 6 (roots 1,2,3)",
            Poly_cubic_123,
            Poly_cubic_123_deriv,
            { 1.0, 2.0, 3.0 },
            { 1, 1, 1 },
            { {bracket_cubic123_root1_low, bracket_cubic123_root1_high},
              {bracket_cubic123_root2_low, bracket_cubic123_root2_high},
              {bracket_cubic123_root3_low, bracket_cubic123_root3_high} },
            "polynomial",
            "Cubic with three well-separated integer roots. Good for testing multiple root finding.",
            true, false, 0, 1
        };
    }

    inline TestRootFunction getQuadratic23Test()
    {
        return TestRootFunction{
            "x² - 5x + 6 (roots 2,3)",
            Poly_quad_23,
            Poly_quad_23_deriv,
            { 2.0, 3.0 },
            { 1, 1 },
            { {1.5, 2.5}, {2.5, 3.5} },
            "polynomial",
            "Simple quadratic with integer roots at 2 and 3.",
            true, false, 0, 1
        };
    }

    inline TestRootFunction getQuartic1122Test()
    {
        return TestRootFunction{
            "x⁴ - 5x² + 4 (roots ±1, ±2)",
            Poly_quartic_1122,
            Poly_quartic_1122_deriv,
            { -2.0, -1.0, 1.0, 2.0 },
            { 1, 1, 1, 1 },
            { {-2.5, -1.5}, {-1.5, -0.5}, {0.5, 1.5}, {1.5, 2.5} },
            "polynomial",
            "Quartic with four symmetric integer roots.",
            true, false, 0, 2
        };
    }

    inline TestRootFunction getWilkinson5Test()
    {
        return TestRootFunction{
            "Wilkinson-5 polynomial (roots 1-5)",
            Poly_wilkinson5,
            Poly_wilkinson5_deriv,
            { 1.0, 2.0, 3.0, 4.0, 5.0 },
            { 1, 1, 1, 1, 1 },
            { {0.5, 1.5}, {1.5, 2.5}, {2.5, 3.5}, {3.5, 4.5}, {4.5, 5.5} },
            "polynomial",
            "Small Wilkinson polynomial - sensitive to perturbations. Good for testing robustness.",
            true, false, 0, 2
        };
    }

    // ========================================================================================
    // TRANSCENDENTAL TEST CASES
    // ========================================================================================

    inline TestRootFunction getDottieNumberTest()
    {
        return TestRootFunction{
            "x - cos(x) (Dottie number)",
            Trans_x_minus_cosx,
            Trans_x_minus_cosx_deriv,
            { ROOT_DOTTIE },
            { 1 },
            { {bracket_dottie_low, bracket_dottie_high} },
            "transcendental",
            "Fixed point equation x = cos(x). The Dottie number is a mathematical constant.",
            true, false, 0, 1
        };
    }

    inline TestRootFunction getKeplerEquationTest()
    {
        return TestRootFunction{
            "Kepler equation (e=0.5, M=0.5)",
            Trans_kepler,
            Trans_kepler_deriv,
            { ROOT_KEPLER_05_05 },
            { 1 },
            { {0.0, 1.5} },
            "transcendental",
            "Kepler's equation E - e*sin(E) = M from orbital mechanics. Historically important.",
            true, false, 0, 2
        };
    }

    inline TestRootFunction getOmegaConstantTest()
    {
        return TestRootFunction{
            "exp(-x) - x (Omega constant)",
            Trans_exp_minus_x,
            Trans_exp_minus_x_deriv,
            { ROOT_OMEGA },
            { 1 },
            { {bracket_omega_low, bracket_omega_high} },
            "transcendental",
            "Intersection of exp(-x) and x. Root is the Omega constant Ω ≈ 0.5671.",
            true, false, 0, 1
        };
    }

    inline TestRootFunction getEulerNumberTest()
    {
        return TestRootFunction{
            "ln(x) - 1 (root at e)",
            Trans_lnx_minus_1,
            Trans_lnx_minus_1_deriv,
            { ROOT_E },
            { 1 },
            { {bracket_e_low, bracket_e_high} },
            "transcendental",
            "Simple logarithmic equation with root at Euler's number e.",
            true, true, 0, 1  // singularity at x=0
        };
    }

    inline TestRootFunction getSinHalfTest()
    {
        return TestRootFunction{
            "sin(x) - 0.5 (root at π/6)",
            Trans_sinx_minus_half,
            Trans_sinx_minus_half_deriv,
            { ROOT_PI_6 },
            { 1 },
            { {bracket_sinpi6_low, bracket_sinpi6_high} },
            "transcendental",
            "Simple trigonometric equation. Multiple roots exist (5π/6, etc.) but bracket selects π/6.",
            true, false, 0, 1
        };
    }

    inline TestRootFunction getTangentIntersectionTest()
    {
        return TestRootFunction{
            "tan(x) - x (first non-trivial root)",
            Trans_tanx_minus_x,
            Trans_tanx_minus_x_deriv,
            { 0.0, ROOT_TAN_1 },
            { 1, 1 },
            { {-0.1, 0.1}, {4.4, 4.6} },
            "transcendental",
            "Tangent intersection. Infinitely many roots near kπ. Care needed near singularities at π/2 + kπ.",
            true, true, static_cast<Real>(Constants::PI / 2), 3  // singularity near bracket
        };
    }

    // ========================================================================================
    // MULTIPLE/REPEATED ROOT TEST CASES
    // ========================================================================================

    inline TestRootFunction getDoubleRootTest()
    {
        return TestRootFunction{
            "(x-1)² double root",
            Multi_double_root_1,
            Multi_double_root_1_deriv,
            { 1.0 },
            { 2 },
            { {bracket_double_root_low, bracket_double_root_high} },
            "multiple_roots",
            "Double root at x=1. Newton's method converges only linearly, not quadratically.",
            true, false, 0, 2
        };
    }

    inline TestRootFunction getTripleRootTest()
    {
        return TestRootFunction{
            "(x-2)³ triple root",
            Multi_triple_root_2,
            Multi_triple_root_2_deriv,
            { 2.0 },
            { 3 },
            { {bracket_triple_root_low, bracket_triple_root_high} },
            "multiple_roots",
            "Triple root at x=2. Even slower convergence for Newton's method.",
            true, false, 0, 3
        };
    }

    inline TestRootFunction getMixedMultiplicityTest()
    {
        return TestRootFunction{
            "(x-1)²(x-3) mixed multiplicity",
            Multi_mixed_1_3,
            Multi_mixed_1_3_deriv,
            { 1.0, 3.0 },
            { 2, 1 },
            { {0.0, 2.0}, {2.0, 4.0} },
            "multiple_roots",
            "Root at 1 with multiplicity 2, simple root at 3. Tests handling of mixed cases.",
            true, false, 0, 3
        };
    }

    // ========================================================================================
    // CLOSELY SPACED ROOTS TEST CASES
    // ========================================================================================

    inline TestRootFunction getCloseRoots1Test()
    {
        return TestRootFunction{
            "(x-1)(x-1.001) close roots",
            Close_roots_1_1001,
            Close_roots_1_1001_deriv,
            { 1.0, 1.001 },
            { 1, 1 },
            { {0.99, 1.0005}, {1.0005, 1.01} },
            "close_roots",
            "Two roots separated by only 0.001. Tests precision and bracketing accuracy.",
            true, false, 0, 3
        };
    }

    inline TestRootFunction getCloseRoots2Test()
    {
        return TestRootFunction{
            "(x-2)(x-2.0001) very close roots",
            Close_roots_2_20001,
            Close_roots_2_20001_deriv,
            { 2.0, 2.0001 },
            { 1, 1 },
            { {1.99, 2.00005}, {2.00005, 2.01} },
            "close_roots",
            "Two roots separated by only 0.0001. Extremely challenging for bracket-based methods.",
            true, false, 0, 4
        };
    }

    // ========================================================================================
    // PATHOLOGICAL TEST CASES
    // ========================================================================================

    inline TestRootFunction getSteepExponentialTest()
    {
        return TestRootFunction{
            "exp(x) - 10000 (steep exponential)",
            Patho_steep_exp,
            Patho_steep_exp_deriv,
            { ROOT_LN10000 },
            { 1 },
            { {8.0, 10.0} },
            "pathological",
            "Very steep exponential. Large derivative causes issues for some methods.",
            true, false, 0, 2
        };
    }

    inline TestRootFunction getFlatCubicTest()
    {
        return TestRootFunction{
            "(x-1)³ - 0.001 (flat near root)",
            Patho_flat_cubic,
            Patho_flat_cubic_deriv,
            { 1.0 + std::cbrt(0.001) },  // ≈ 1.1
            { 1 },
            { {0.5, 1.5} },
            "pathological",
            "Very flat near the root. Slow convergence for derivative-based methods.",
            true, false, 0, 3
        };
    }

    inline TestRootFunction getOscillatoryDecayTest()
    {
        return TestRootFunction{
            "sin(10x)*exp(-x) (oscillatory decay)",
            Patho_oscill_decay,
            Patho_oscill_decay_deriv,
            { 0.0, static_cast<Real>(Constants::PI / 10), 
              static_cast<Real>(Constants::PI / 5), 
              static_cast<Real>(3 * Constants::PI / 10) },
            { 1, 1, 1, 1 },
            { {-0.1, 0.1}, 
              {0.25, 0.35}, 
              {0.55, 0.7}, 
              {0.85, 1.0} },
            "pathological",
            "Oscillatory function with exponential decay. Many roots, needs careful bracketing.",
            true, false, 0, 3
        };
    }

    inline TestRootFunction getSteepAtanTest()
    {
        return TestRootFunction{
            "atan(1000(x-1)) (extremely steep)",
            Patho_steep_atan,
            Patho_steep_atan_deriv,
            { 1.0 },
            { 1 },
            { {0.0, 2.0} },
            "pathological",
            "Extremely steep at x=1 (derivative ~1000). Bisection stable, Newton may overshoot.",
            true, false, 0, 3
        };
    }

    // ========================================================================================
    // SINGULARITY TEST CASES
    // ========================================================================================

    inline TestRootFunction getSingularityNearRootTest()
    {
        return TestRootFunction{
            "1/x - 2 (singularity at 0)",
            Sing_1_over_x_minus_2,
            Sing_1_over_x_minus_2_deriv,
            { 0.5 },
            { 1 },
            { {0.1, 1.0} },
            "singularity",
            "Root at 0.5 but singularity at x=0. Bracket must exclude singularity.",
            true, true, 0.0, 2
        };
    }

    inline TestRootFunction getLogSingularityTest()
    {
        return TestRootFunction{
            "ln(x) + x - 2 (log singularity)",
            Sing_lnx_plus_x_minus_2,
            Sing_lnx_plus_x_minus_2_deriv,
            { 1.5571455989976115 },  // Approximate
            { 1 },
            { {0.5, 2.5} },
            "singularity",
            "Log singularity at x=0, root near x≈1.557. Common in optimization.",
            true, true, 0.0, 2
        };
    }

    // ========================================================================================
    // PHYSICS/ENGINEERING TEST CASES
    // ========================================================================================

    inline TestRootFunction getVanDerWaalsTest()
    {
        return TestRootFunction{
            "Van der Waals cubic",
            Phys_vdw,
            Phys_vdw_deriv,
            { 0.10874 },  // Approximate numerical root
            { 1 },
            { {0.05, 0.5} },
            "physics",
            "Simplified Van der Waals equation for molar volume. Real gas thermodynamics.",
            true, false, 0, 2
        };
    }

    inline TestRootFunction getPlanckRadiationTest()
    {
        return TestRootFunction{
            "Planck radiation crossover",
            Phys_planck,
            Phys_planck_deriv,
            { ROOT_WIEN },
            { 1 },
            { {4.0, 6.0} },
            "physics",
            "Related to Wien displacement law. Root ≈ 4.965 relates to hc/kT ratio.",
            true, false, 0, 2
        };
    }

    inline TestRootFunction getLambertWTest()
    {
        return TestRootFunction{
            "x*exp(x) - 1 (Lambert W)",
            Phys_lambert,
            Phys_lambert_deriv,
            { ROOT_LAMBERT_W1 },
            { 1 },
            { {0.0, 1.0} },
            "physics",
            "Lambert W function at W(1). Appears in delay differential equations, combinatorics.",
            true, false, 0, 2
        };
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                           CATEGORY RETRIEVAL FUNCTIONS                                 //
    ///////////////////////////////////////////////////////////////////////////////////////////

    /**
     * @brief Get all available root finding test cases
     */
    inline std::vector<TestRootFunction> getAllRootFindingTests()
    {
        return {
            // Polynomials
            getQuadraticSqrt2Test(),
            getCubic123Test(),
            getQuadratic23Test(),
            getQuartic1122Test(),
            getWilkinson5Test(),
            // Transcendental
            getDottieNumberTest(),
            getKeplerEquationTest(),
            getOmegaConstantTest(),
            getEulerNumberTest(),
            getSinHalfTest(),
            getTangentIntersectionTest(),
            // Multiple roots
            getDoubleRootTest(),
            getTripleRootTest(),
            getMixedMultiplicityTest(),
            // Close roots
            getCloseRoots1Test(),
            getCloseRoots2Test(),
            // Pathological
            getSteepExponentialTest(),
            getFlatCubicTest(),
            getOscillatoryDecayTest(),
            getSteepAtanTest(),
            // Singularities
            getSingularityNearRootTest(),
            getLogSingularityTest(),
            // Physics
            getVanDerWaalsTest(),
            getPlanckRadiationTest(),
            getLambertWTest()
        };
    }

    /**
     * @brief Get polynomial test cases only
     */
    inline std::vector<TestRootFunction> getPolynomialRootTests()
    {
        return {
            getQuadraticSqrt2Test(),
            getCubic123Test(),
            getQuadratic23Test(),
            getQuartic1122Test(),
            getWilkinson5Test()
        };
    }

    /**
     * @brief Get transcendental equation test cases
     */
    inline std::vector<TestRootFunction> getTranscendentalRootTests()
    {
        return {
            getDottieNumberTest(),
            getKeplerEquationTest(),
            getOmegaConstantTest(),
            getEulerNumberTest(),
            getSinHalfTest(),
            getTangentIntersectionTest()
        };
    }

    /**
     * @brief Get multiple root test cases
     */
    inline std::vector<TestRootFunction> getMultipleRootTests()
    {
        return {
            getDoubleRootTest(),
            getTripleRootTest(),
            getMixedMultiplicityTest()
        };
    }

    /**
     * @brief Get closely spaced root test cases
     */
    inline std::vector<TestRootFunction> getCloseRootTests()
    {
        return {
            getCloseRoots1Test(),
            getCloseRoots2Test()
        };
    }

    /**
     * @brief Get pathological/challenging test cases
     */
    inline std::vector<TestRootFunction> getPathologicalRootTests()
    {
        return {
            getSteepExponentialTest(),
            getFlatCubicTest(),
            getOscillatoryDecayTest(),
            getSteepAtanTest()
        };
    }

    /**
     * @brief Get test cases with singularities
     */
    inline std::vector<TestRootFunction> getSingularityRootTests()
    {
        return {
            getSingularityNearRootTest(),
            getLogSingularityTest()
        };
    }

    /**
     * @brief Get physics/engineering application test cases
     */
    inline std::vector<TestRootFunction> getPhysicsRootTests()
    {
        return {
            getVanDerWaalsTest(),
            getPlanckRadiationTest(),
            getLambertWTest()
        };
    }

    /**
     * @brief Get easy test cases (difficulty 1)
     */
    inline std::vector<TestRootFunction> getEasyRootTests()
    {
        std::vector<TestRootFunction> result;
        for (const auto& test : getAllRootFindingTests()) {
            if (test.difficulty == 1) result.push_back(test);
        }
        return result;
    }

    /**
     * @brief Get challenging test cases (difficulty 3+)
     */
    inline std::vector<TestRootFunction> getChallengingRootTests()
    {
        std::vector<TestRootFunction> result;
        for (const auto& test : getAllRootFindingTests()) {
            if (test.difficulty >= 3) result.push_back(test);
        }
        return result;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                              UTILITY FUNCTIONS                                          //
    ///////////////////////////////////////////////////////////////////////////////////////////

    /**
     * @brief Create a custom test case with user-defined function
     */
    inline TestRootFunction createCustomRootTest(
        const std::string& name,
        std::function<Real(Real)> func,
        std::function<Real(Real)> derivative,
        Real knownRoot,
        Real bracketLow,
        Real bracketHigh,
        const std::string& description = "Custom test case")
    {
        return TestRootFunction{
            name,
            func,
            derivative,
            { knownRoot },
            { 1 },
            { {bracketLow, bracketHigh} },
            "custom",
            description,
            true, false, 0, 2
        };
    }

    /**
     * @brief Verify a computed root against the test case
     * @param test The test case
     * @param computedRoot The root found by the algorithm
     * @param rootIndex Which root to compare against (default 0)
     * @param tolerance Acceptable error tolerance
     * @return True if the computed root is within tolerance
     */
    inline bool verifyRoot(const TestRootFunction& test, Real computedRoot, 
                           size_t rootIndex = 0, Real tolerance = 1e-10)
    {
        if (rootIndex >= test.knownRoots.size()) return false;
        return std::abs(computedRoot - test.knownRoots[rootIndex]) < tolerance;
    }

    /**
     * @brief Compute the error of a root approximation
     */
    inline Real computeRootError(const TestRootFunction& test, Real computedRoot, size_t rootIndex = 0)
    {
        if (rootIndex >= test.knownRoots.size()) return std::numeric_limits<Real>::max();
        return std::abs(computedRoot - test.knownRoots[rootIndex]);
    }

    /**
     * @brief Get the function value at a computed root (should be ~0)
     */
    inline Real evaluateFunctionAtRoot(const TestRootFunction& test, Real computedRoot)
    {
        return std::abs(test.func(computedRoot));
    }

} // namespace MML::TestBeds

#endif // __MML_ROOT_FINDING_TEST_BED_H

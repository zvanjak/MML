#if !defined __MML_OPTIMIZATION_TEST_BED_H
#define __MML_OPTIMIZATION_TEST_BED_H

#include <string>
#include <vector>
#include <functional>
#include <cmath>
#include <limits>

#include "optimization_defs.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#include "base/Vector.h"
#include "interfaces/IFunction.h"
#endif

namespace MML::TestBeds
{
    ///////////////////////////////////////////////////////////////////////////////////////////
    //                      1D OPTIMIZATION TEST DATA STRUCTURE                               //
    ///////////////////////////////////////////////////////////////////////////////////////////

    /**
     * @brief Represents a test case for 1D optimization algorithms
     * 
     * Contains the function, optional derivative, known minimum location and value,
     * search range, and metadata for categorization.
     */
    struct TestOptimization1D
    {
        std::string name;                               ///< Descriptive name for the test
        std::function<Real(Real)> func;                ///< The objective function f(x)
        std::function<Real(Real)> derivative;          ///< f'(x) for gradient-based methods
        Real trueMinimumX;                             ///< Known x location of minimum
        Real trueMinimumF;                             ///< Known f(x*) value at minimum
        Real searchLow;                                ///< Lower bound of search range
        Real searchHigh;                               ///< Upper bound of search range
        
        // Metadata
        std::string category;                          ///< unimodal, multimodal, challenging
        std::string description;                       ///< Detailed description
        bool hasDerivative = true;                     ///< Whether derivative is available
        bool isMultimodal = false;                     ///< Has multiple local minima
        int difficulty = 1;                            ///< 1=easy, 2=medium, 3=hard
    };

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                      N-D OPTIMIZATION TEST DATA STRUCTURE                              //
    ///////////////////////////////////////////////////////////////////////////////////////////

    /**
     * @brief Represents a test case for N-dimensional optimization algorithms
     * 
     * Template parameter N specifies the dimension.
     */
    template<int N>
    struct TestOptimizationND
    {
        std::string name;                                      ///< Descriptive name
        std::function<Real(const VectorN<Real, N>&)> func;    ///< Objective function
        std::function<VectorN<Real, N>(const VectorN<Real, N>&)> gradient;  ///< Gradient
        std::vector<VectorN<Real, N>> trueMinima;             ///< Known minimum locations
        Real trueMinimumF;                                     ///< Known minimum value
        VectorN<Real, N> domainLow;                           ///< Lower bounds of search domain
        VectorN<Real, N> domainHigh;                          ///< Upper bounds of search domain
        VectorN<Real, N> suggestedStart;                      ///< Suggested starting point
        
        // Metadata
        std::string category;                                  ///< convex, valley, multimodal, etc.
        std::string description;                               ///< Detailed description
        bool hasGradient = false;                              ///< Whether gradient is available
        bool isMultimodal = false;                             ///< Has multiple local minima
        int numMinima = 1;                                     ///< Number of global minima
        int difficulty = 1;                                    ///< 1=easy, 2=medium, 3=hard, 4=very hard
    };

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                        WRAPPER CLASS FOR MML INTERFACE                                 //
    ///////////////////////////////////////////////////////////////////////////////////////////

    /**
     * @brief Wrapper to adapt std::function to MML's IRealFunction interface (1D)
     */
    class OptRealFunctionWrapper : public IRealFunction
    {
    private:
        std::function<Real(Real)> _func;
    public:
        OptRealFunctionWrapper(std::function<Real(Real)> f) : _func(f) {}
        Real operator()(Real x) const override { return _func(x); }
    };

    /**
     * @brief Wrapper to adapt std::function to MML's IScalarFunction interface (ND)
     */
    template<int N>
    class OptScalarFunctionWrapper : public IScalarFunction<N>
    {
    private:
        std::function<Real(const VectorN<Real, N>&)> _func;
    public:
        OptScalarFunctionWrapper(std::function<Real(const VectorN<Real, N>&)> f) : _func(f) {}
        Real operator()(const VectorN<Real, N>& x) const override { return _func(x); }
    };

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                          1D TEST CASE GENERATORS                                        //
    ///////////////////////////////////////////////////////////////////////////////////////////

    // ========================================================================================
    // UNIMODAL 1D FUNCTIONS
    // ========================================================================================

    inline TestOptimization1D getQuadratic1DTest()
    {
        return TestOptimization1D{
            "Quadratic (x-2)² + 1",
            Opt1D_Quadratic,
            Opt1D_Quadratic_deriv,
            Opt1D_Quadratic_min_x,
            Opt1D_Quadratic_min_f,
            -5.0, 10.0,
            "unimodal",
            "Simple quadratic with minimum at x=2. Easiest possible test.",
            true, false, 1
        };
    }

    inline TestOptimization1D getQuartic1DTest()
    {
        return TestOptimization1D{
            "Quartic (x-1)⁴",
            Opt1D_Quartic,
            Opt1D_Quartic_deriv,
            Opt1D_Quartic_min_x,
            Opt1D_Quartic_min_f,
            -3.0, 5.0,
            "unimodal",
            "Quartic with very flat region near minimum. Tests convergence detection.",
            true, false, 2
        };
    }

    inline TestOptimization1D getCosh1DTest()
    {
        return TestOptimization1D{
            "Hyperbolic Cosine",
            Opt1D_Cosh,
            Opt1D_Cosh_deriv,
            Opt1D_Cosh_min_x,
            Opt1D_Cosh_min_f,
            0.0, 6.0,
            "unimodal",
            "Shifted cosh function with exponential growth away from minimum.",
            true, false, 1
        };
    }

    inline TestOptimization1D getAbsValue1DTest()
    {
        return TestOptimization1D{
            "Absolute Value |x-1.5|+2",
            Opt1D_AbsValue,
            nullptr,  // No derivative at minimum
            Opt1D_AbsValue_min_x,
            Opt1D_AbsValue_min_f,
            -2.0, 5.0,
            "unimodal",
            "V-shaped function, non-differentiable at minimum. Tests robustness.",
            false, false, 2
        };
    }

    inline TestOptimization1D getLogBarrier1DTest()
    {
        return TestOptimization1D{
            "Log Barrier x² - ln(x)",
            Opt1D_LogBarrier,
            Opt1D_LogBarrier_deriv,
            Opt1D_LogBarrier_min_x,
            Opt1D_LogBarrier_min_f,
            0.1, 5.0,
            "unimodal",
            "Common in constrained optimization. Minimum at x=1/√2.",
            true, false, 2
        };
    }

    // ========================================================================================
    // MULTIMODAL 1D FUNCTIONS
    // ========================================================================================

    inline TestOptimization1D getRastrigin1DTest()
    {
        return TestOptimization1D{
            "Rastrigin 1D",
            Opt1D_Rastrigin,
            Opt1D_Rastrigin_deriv,
            Opt1D_Rastrigin_min_x,
            Opt1D_Rastrigin_min_f,
            -5.12, 5.12,
            "multimodal",
            "Highly multimodal with many local minima. Global at x=0.",
            true, true, 3
        };
    }

    inline TestOptimization1D getAckley1DTest()
    {
        return TestOptimization1D{
            "Ackley 1D",
            Opt1D_Ackley,
            nullptr,  // Complex derivative
            Opt1D_Ackley_min_x,
            Opt1D_Ackley_min_f,
            -32.768, 32.768,
            "multimodal",
            "Nearly flat outer region with sharp global minimum at x=0.",
            false, true, 3
        };
    }

    inline TestOptimization1D getSinusoidal1DTest()
    {
        return TestOptimization1D{
            "Sinusoidal sin(x)+sin(3x)/3",
            Opt1D_Sinusoidal,
            Opt1D_Sinusoidal_deriv,
            -1.9, // Approximate global minimum
            -0.87, // Approximate value
            -Constants::PI, 2 * Constants::PI,
            "multimodal",
            "Oscillating function with multiple local minima.",
            true, true, 2
        };
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                          2D TEST CASE GENERATORS                                        //
    ///////////////////////////////////////////////////////////////////////////////////////////

    inline TestOptimizationND<2> getRosenbrock2DTest()
    {
        TestOptimizationND<2> test;
        test.name = "Rosenbrock (Banana)";
        test.func = OptND_Rosenbrock2D;
        test.gradient = OptND_Rosenbrock2D_grad;
        test.trueMinima = { OptND_Rosenbrock2D_min };
        test.trueMinimumF = OptND_Rosenbrock2D_min_f;
        test.domainLow = {Domain_Rosenbrock_low, Domain_Rosenbrock_low};
        test.domainHigh = {Domain_Rosenbrock_high, Domain_Rosenbrock_high};
        test.suggestedStart = {-1.0, 1.0};
        test.category = "valley";
        test.description = "Classic banana-shaped valley. Minimum at (1,1). Standard benchmark.";
        test.hasGradient = true;
        test.isMultimodal = false;
        test.numMinima = 1;
        test.difficulty = 2;
        return test;
    }

    inline TestOptimizationND<2> getBeale2DTest()
    {
        TestOptimizationND<2> test;
        test.name = "Beale";
        test.func = OptND_Beale;
        test.trueMinima = { OptND_Beale_min };
        test.trueMinimumF = OptND_Beale_min_f;
        test.domainLow = {Domain_Beale_low, Domain_Beale_low};
        test.domainHigh = {Domain_Beale_high, Domain_Beale_high};
        test.suggestedStart = {0.0, 0.0};
        test.category = "unimodal";
        test.description = "Unimodal with flat regions. Minimum at (3, 0.5).";
        test.hasGradient = false;
        test.isMultimodal = false;
        test.numMinima = 1;
        test.difficulty = 2;
        return test;
    }

    inline TestOptimizationND<2> getBooth2DTest()
    {
        TestOptimizationND<2> test;
        test.name = "Booth";
        test.func = OptND_Booth;
        test.trueMinima = { OptND_Booth_min };
        test.trueMinimumF = OptND_Booth_min_f;
        test.domainLow = {-10.0, -10.0};
        test.domainHigh = {10.0, 10.0};
        test.suggestedStart = {0.0, 0.0};
        test.category = "convex";
        test.description = "Simple quadratic bowl. Minimum at (1, 3). Easy benchmark.";
        test.hasGradient = false;
        test.isMultimodal = false;
        test.numMinima = 1;
        test.difficulty = 1;
        return test;
    }

    inline TestOptimizationND<2> getMatyas2DTest()
    {
        TestOptimizationND<2> test;
        test.name = "Matyas";
        test.func = OptND_Matyas;
        test.trueMinima = { OptND_Matyas_min };
        test.trueMinimumF = OptND_Matyas_min_f;
        test.domainLow = {-10.0, -10.0};
        test.domainHigh = {10.0, 10.0};
        test.suggestedStart = {5.0, 5.0};
        test.category = "convex";
        test.description = "Plate-shaped, nearly flat. Minimum at origin.";
        test.hasGradient = false;
        test.isMultimodal = false;
        test.numMinima = 1;
        test.difficulty = 1;
        return test;
    }

    inline TestOptimizationND<2> getHimmelblau2DTest()
    {
        TestOptimizationND<2> test;
        test.name = "Himmelblau";
        test.func = OptND_Himmelblau;
        test.trueMinima = { 
            OptND_Himmelblau_min1, 
            OptND_Himmelblau_min2, 
            OptND_Himmelblau_min3, 
            OptND_Himmelblau_min4 
        };
        test.trueMinimumF = OptND_Himmelblau_min_f;
        test.domainLow = {Domain_Himmelblau_low, Domain_Himmelblau_low};
        test.domainHigh = {Domain_Himmelblau_high, Domain_Himmelblau_high};
        test.suggestedStart = {0.0, 0.0};
        test.category = "multimodal";
        test.description = "Four identical global minima. Tests which minimum is found.";
        test.hasGradient = false;
        test.isMultimodal = true;
        test.numMinima = 4;
        test.difficulty = 2;
        return test;
    }

    inline TestOptimizationND<2> getGoldsteinPrice2DTest()
    {
        TestOptimizationND<2> test;
        test.name = "Goldstein-Price";
        test.func = OptND_GoldsteinPrice;
        test.trueMinima = { OptND_GoldsteinPrice_min };
        test.trueMinimumF = OptND_GoldsteinPrice_min_f;
        test.domainLow = {Domain_GoldsteinPrice_low, Domain_GoldsteinPrice_low};
        test.domainHigh = {Domain_GoldsteinPrice_high, Domain_GoldsteinPrice_high};
        test.suggestedStart = {0.0, 0.0};
        test.category = "multimodal";
        test.description = "Complex landscape with several local minima. Global at (0,-1)=3.";
        test.hasGradient = false;
        test.isMultimodal = true;
        test.numMinima = 1;
        test.difficulty = 3;
        return test;
    }

    inline TestOptimizationND<2> getThreeHumpCamel2DTest()
    {
        TestOptimizationND<2> test;
        test.name = "Three-Hump Camel";
        test.func = OptND_ThreeHumpCamel;
        test.trueMinima = { OptND_ThreeHumpCamel_min };
        test.trueMinimumF = OptND_ThreeHumpCamel_min_f;
        test.domainLow = {-5.0, -5.0};
        test.domainHigh = {5.0, 5.0};
        test.suggestedStart = {1.0, 1.0};
        test.category = "multimodal";
        test.description = "Three local minima, global at origin.";
        test.hasGradient = false;
        test.isMultimodal = true;
        test.numMinima = 1;
        test.difficulty = 2;
        return test;
    }

    inline TestOptimizationND<2> getSixHumpCamel2DTest()
    {
        TestOptimizationND<2> test;
        test.name = "Six-Hump Camel";
        test.func = OptND_SixHumpCamel;
        test.trueMinima = { OptND_SixHumpCamel_min1, OptND_SixHumpCamel_min2 };
        test.trueMinimumF = OptND_SixHumpCamel_min_f;
        test.domainLow = {-3.0, -2.0};
        test.domainHigh = {3.0, 2.0};
        test.suggestedStart = {0.0, 0.0};
        test.category = "multimodal";
        test.description = "Six local minima, two global at f≈-1.0316.";
        test.hasGradient = false;
        test.isMultimodal = true;
        test.numMinima = 2;
        test.difficulty = 2;
        return test;
    }

    inline TestOptimizationND<2> getEasom2DTest()
    {
        TestOptimizationND<2> test;
        test.name = "Easom";
        test.func = OptND_Easom;
        test.trueMinima = { OptND_Easom_min };
        test.trueMinimumF = OptND_Easom_min_f;
        test.domainLow = {-100.0, -100.0};
        test.domainHigh = {100.0, 100.0};
        test.suggestedStart = {0.0, 0.0};
        test.category = "unimodal";
        test.description = "Nearly flat with sharp global minimum at (π,π). Hard to locate.";
        test.hasGradient = false;
        test.isMultimodal = false;
        test.numMinima = 1;
        test.difficulty = 3;
        return test;
    }

    inline TestOptimizationND<2> getSphere2DTest()
    {
        TestOptimizationND<2> test;
        test.name = "Sphere 2D";
        test.func = OptND_Sphere<2>;
        test.gradient = OptND_Sphere_grad<2>;
        test.trueMinima = { {0.0, 0.0} };
        test.trueMinimumF = 0.0;
        test.domainLow = {Domain_Sphere_low, Domain_Sphere_low};
        test.domainHigh = {Domain_Sphere_high, Domain_Sphere_high};
        test.suggestedStart = {4.0, -3.0};
        test.category = "convex";
        test.description = "Simplest benchmark - spherical bowl. Minimum at origin.";
        test.hasGradient = true;
        test.isMultimodal = false;
        test.numMinima = 1;
        test.difficulty = 1;
        return test;
    }

    inline TestOptimizationND<2> getRastrigin2DTest()
    {
        TestOptimizationND<2> test;
        test.name = "Rastrigin 2D";
        test.func = OptND_Rastrigin<2>;
        test.trueMinima = { {0.0, 0.0} };
        test.trueMinimumF = 0.0;
        test.domainLow = {Domain_Rastrigin_low, Domain_Rastrigin_low};
        test.domainHigh = {Domain_Rastrigin_high, Domain_Rastrigin_high};
        test.suggestedStart = {2.0, 2.0};
        test.category = "multimodal";
        test.description = "Highly multimodal with ~100 local minima. Global at origin.";
        test.hasGradient = false;
        test.isMultimodal = true;
        test.numMinima = 1;
        test.difficulty = 3;
        return test;
    }

    inline TestOptimizationND<2> getAckley2DTest()
    {
        TestOptimizationND<2> test;
        test.name = "Ackley 2D";
        test.func = OptND_Ackley<2>;
        test.trueMinima = { {0.0, 0.0} };
        test.trueMinimumF = OptND_Ackley_min_f;
        test.domainLow = {Domain_Ackley_low, Domain_Ackley_low};
        test.domainHigh = {Domain_Ackley_high, Domain_Ackley_high};
        test.suggestedStart = {10.0, 10.0};
        test.category = "multimodal";
        test.description = "Nearly flat with many local minima. Global at origin.";
        test.hasGradient = false;
        test.isMultimodal = true;
        test.numMinima = 1;
        test.difficulty = 3;
        return test;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                     HIGHER-DIMENSIONAL TEST CASE GENERATORS                            //
    ///////////////////////////////////////////////////////////////////////////////////////////

    template<int N>
    inline TestOptimizationND<N> getSphereNDTest()
    {
        TestOptimizationND<N> test;
        test.name = "Sphere " + std::to_string(N) + "D";
        test.func = OptND_Sphere<N>;
        test.gradient = OptND_Sphere_grad<N>;
        
        VectorN<Real, N> zero, low, high, start;
        for (int i = 0; i < N; ++i) {
            zero[i] = 0.0;
            low[i] = Domain_Sphere_low;
            high[i] = Domain_Sphere_high;
            start[i] = 3.0;  // Start away from minimum
        }
        
        test.trueMinima = { zero };
        test.trueMinimumF = 0.0;
        test.domainLow = low;
        test.domainHigh = high;
        test.suggestedStart = start;
        test.category = "convex";
        test.description = "N-dimensional sphere function. Simplest possible benchmark.";
        test.hasGradient = true;
        test.isMultimodal = false;
        test.numMinima = 1;
        test.difficulty = 1;
        return test;
    }

    template<int N>
    inline TestOptimizationND<N> getRosenbrockNDTest()
    {
        TestOptimizationND<N> test;
        test.name = "Rosenbrock " + std::to_string(N) + "D";
        test.func = OptND_RosenbrockND<N>;
        test.gradient = OptND_RosenbrockND_grad<N>;
        
        VectorN<Real, N> ones, low, high, start;
        for (int i = 0; i < N; ++i) {
            ones[i] = 1.0;
            low[i] = Domain_Rosenbrock_low;
            high[i] = Domain_Rosenbrock_high;
            start[i] = (i % 2 == 0) ? -1.0 : 1.0;
        }
        
        test.trueMinima = { ones };
        test.trueMinimumF = 0.0;
        test.domainLow = low;
        test.domainHigh = high;
        test.suggestedStart = start;
        test.category = "valley";
        test.description = "N-dimensional Rosenbrock with curved valley. Minimum at (1,...,1).";
        test.hasGradient = true;
        test.isMultimodal = false;
        test.numMinima = 1;
        test.difficulty = 3;
        return test;
    }

    template<int N>
    inline TestOptimizationND<N> getRastriginNDTest()
    {
        TestOptimizationND<N> test;
        test.name = "Rastrigin " + std::to_string(N) + "D";
        test.func = OptND_Rastrigin<N>;
        
        VectorN<Real, N> zero, low, high, start;
        for (int i = 0; i < N; ++i) {
            zero[i] = 0.0;
            low[i] = Domain_Rastrigin_low;
            high[i] = Domain_Rastrigin_high;
            start[i] = 2.5;
        }
        
        test.trueMinima = { zero };
        test.trueMinimumF = 0.0;
        test.domainLow = low;
        test.domainHigh = high;
        test.suggestedStart = start;
        test.category = "multimodal";
        test.description = "Highly multimodal with ~10^N local minima. Standard global opt benchmark.";
        test.hasGradient = false;
        test.isMultimodal = true;
        test.numMinima = 1;
        test.difficulty = 4;
        return test;
    }

    template<int N>
    inline TestOptimizationND<N> getAckleyNDTest()
    {
        TestOptimizationND<N> test;
        test.name = "Ackley " + std::to_string(N) + "D";
        test.func = OptND_Ackley<N>;
        
        VectorN<Real, N> zero, low, high, start;
        for (int i = 0; i < N; ++i) {
            zero[i] = 0.0;
            low[i] = Domain_Ackley_low;
            high[i] = Domain_Ackley_high;
            start[i] = 15.0;
        }
        
        test.trueMinima = { zero };
        test.trueMinimumF = OptND_Ackley_min_f;
        test.domainLow = low;
        test.domainHigh = high;
        test.suggestedStart = start;
        test.category = "multimodal";
        test.description = "Nearly flat outer region with many local minima.";
        test.hasGradient = false;
        test.isMultimodal = true;
        test.numMinima = 1;
        test.difficulty = 3;
        return test;
    }

    inline TestOptimizationND<4> getPowell4DTest()
    {
        TestOptimizationND<4> test;
        test.name = "Powell 4D";
        test.func = OptND_Powell4D;
        test.trueMinima = { OptND_Powell4D_min };
        test.trueMinimumF = OptND_Powell4D_min_f;
        test.domainLow = {-4.0, -4.0, -4.0, -4.0};
        test.domainHigh = {5.0, 5.0, 5.0, 5.0};
        test.suggestedStart = {3.0, -1.0, 0.0, 1.0};
        test.category = "ill-conditioned";
        test.description = "Ill-conditioned with narrow valley. Minimum at origin.";
        test.hasGradient = false;
        test.isMultimodal = false;
        test.numMinima = 1;
        test.difficulty = 3;
        return test;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                           CATEGORY RETRIEVAL FUNCTIONS                                 //
    ///////////////////////////////////////////////////////////////////////////////////////////

    /**
     * @brief Get all 1D optimization test cases
     */
    inline std::vector<TestOptimization1D> getAll1DOptimizationTests()
    {
        return {
            getQuadratic1DTest(),
            getQuartic1DTest(),
            getCosh1DTest(),
            getAbsValue1DTest(),
            getLogBarrier1DTest(),
            getRastrigin1DTest(),
            getAckley1DTest(),
            getSinusoidal1DTest()
        };
    }

    /**
     * @brief Get unimodal 1D test cases
     */
    inline std::vector<TestOptimization1D> getUnimodal1DTests()
    {
        return {
            getQuadratic1DTest(),
            getQuartic1DTest(),
            getCosh1DTest(),
            getAbsValue1DTest(),
            getLogBarrier1DTest()
        };
    }

    /**
     * @brief Get multimodal 1D test cases
     */
    inline std::vector<TestOptimization1D> getMultimodal1DTests()
    {
        return {
            getRastrigin1DTest(),
            getAckley1DTest(),
            getSinusoidal1DTest()
        };
    }

    /**
     * @brief Get all 2D optimization test cases
     */
    inline std::vector<TestOptimizationND<2>> getAll2DOptimizationTests()
    {
        return {
            getSphere2DTest(),
            getRosenbrock2DTest(),
            getBeale2DTest(),
            getBooth2DTest(),
            getMatyas2DTest(),
            getHimmelblau2DTest(),
            getGoldsteinPrice2DTest(),
            getThreeHumpCamel2DTest(),
            getSixHumpCamel2DTest(),
            getEasom2DTest(),
            getRastrigin2DTest(),
            getAckley2DTest()
        };
    }

    /**
     * @brief Get convex/easy 2D test cases
     */
    inline std::vector<TestOptimizationND<2>> getConvex2DTests()
    {
        return {
            getSphere2DTest(),
            getBooth2DTest(),
            getMatyas2DTest()
        };
    }

    /**
     * @brief Get multimodal 2D test cases
     */
    inline std::vector<TestOptimizationND<2>> getMultimodal2DTests()
    {
        return {
            getHimmelblau2DTest(),
            getGoldsteinPrice2DTest(),
            getThreeHumpCamel2DTest(),
            getSixHumpCamel2DTest(),
            getRastrigin2DTest(),
            getAckley2DTest()
        };
    }

    /**
     * @brief Get classic benchmark 2D tests (the famous ones)
     */
    inline std::vector<TestOptimizationND<2>> getClassicBenchmark2DTests()
    {
        return {
            getRosenbrock2DTest(),
            getBeale2DTest(),
            getHimmelblau2DTest(),
            getGoldsteinPrice2DTest(),
            getRastrigin2DTest(),
            getAckley2DTest()
        };
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                              UTILITY FUNCTIONS                                          //
    ///////////////////////////////////////////////////////////////////////////////////////////

    /**
     * @brief Verify a 1D optimization result
     */
    inline bool verify1DMinimum(const TestOptimization1D& test, Real foundX, Real tolerance = 1e-6)
    {
        return std::abs(foundX - test.trueMinimumX) < tolerance;
    }

    /**
     * @brief Verify an ND optimization result (checks if close to any known minimum)
     */
    template<int N>
    inline bool verifyNDMinimum(const TestOptimizationND<N>& test, const VectorN<Real, N>& found, 
                                Real tolerance = 1e-6)
    {
        for (const auto& minPt : test.trueMinima) {
            Real dist = 0;
            for (int i = 0; i < N; ++i) {
                Real d = found[i] - minPt[i];
                dist += d * d;
            }
            if (std::sqrt(dist) < tolerance)
                return true;
        }
        return false;
    }

    /**
     * @brief Verify an ND optimization result (Vector<Real> overload)
     */
    template<int N>
    inline bool verifyNDMinimum(const TestOptimizationND<N>& test, const Vector<Real>& found, 
                                Real tolerance = 1e-6)
    {
        for (const auto& minPt : test.trueMinima) {
            Real dist = 0;
            for (int i = 0; i < N; ++i) {
                Real d = found[i] - minPt[i];
                dist += d * d;
            }
            if (std::sqrt(dist) < tolerance)
                return true;
        }
        return false;
    }

    /**
     * @brief Compute error from nearest known minimum
     */
    template<int N>
    inline Real computeMinimumError(const TestOptimizationND<N>& test, const VectorN<Real, N>& found)
    {
        Real minDist = std::numeric_limits<Real>::max();
        for (const auto& minPt : test.trueMinima) {
            Real dist = 0;
            for (int i = 0; i < N; ++i) {
                Real d = found[i] - minPt[i];
                dist += d * d;
            }
            minDist = std::min(minDist, std::sqrt(dist));
        }
        return minDist;
    }

    /**
     * @brief Compute error from nearest known minimum (Vector<Real> overload)
     */
    template<int N>
    inline Real computeMinimumError(const TestOptimizationND<N>& test, const Vector<Real>& found)
    {
        Real minDist = std::numeric_limits<Real>::max();
        for (const auto& minPt : test.trueMinima) {
            Real dist = 0;
            for (int i = 0; i < N; ++i) {
                Real d = found[i] - minPt[i];
                dist += d * d;
            }
            minDist = std::min(minDist, std::sqrt(dist));
        }
        return minDist;
    }

    /**
     * @brief Create a custom 1D test case
     */
    inline TestOptimization1D createCustom1DTest(
        const std::string& name,
        std::function<Real(Real)> func,
        Real trueMinX, Real trueMinF,
        Real searchLow, Real searchHigh,
        const std::string& description = "Custom test case")
    {
        return TestOptimization1D{
            name, func, nullptr,
            trueMinX, trueMinF,
            searchLow, searchHigh,
            "custom", description,
            false, false, 2
        };
    }

    /**
     * @brief Generate random starting point within domain
     */
    template<int N>
    inline VectorN<Real, N> generateRandomStart(const TestOptimizationND<N>& test, unsigned int seed = 42)
    {
        VectorN<Real, N> start;
        std::srand(seed);
        for (int i = 0; i < N; ++i) {
            Real range = test.domainHigh[i] - test.domainLow[i];
            start[i] = test.domainLow[i] + range * (std::rand() / static_cast<Real>(RAND_MAX));
        }
        return start;
    }

} // namespace MML::TestBeds

#endif // __MML_OPTIMIZATION_TEST_BED_H

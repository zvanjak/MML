#if !defined __MML_INTERPOLATION_TEST_BED_H
#define __MML_INTERPOLATION_TEST_BED_H

#include <string>
#include <cmath>
#include <vector>
#include <functional>

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#include "base/Vector.h"
#endif

#include "interpolation_defs.h"

namespace MML::TestBeds
{
    ///////////////////////////////////////////////////////////////////////////////////////////
    //                       INTERPOLATION TEST DATA STRUCTURES                               //
    ///////////////////////////////////////////////////////////////////////////////////////////

    // Test data for interpolation testing
    struct TestInterpolationData
    {
        std::string _name;              // Human-readable test name
        std::string _description;       // What this test is checking
        
        Vector<Real> _xData;            // Input x coordinates
        Vector<Real> _yData;            // Input y values (may include noise)
        
        // True function for error measurement
        Real (*_trueFunction)(Real);
        std::string _trueFuncExpr;      // Mathematical expression
        
        // Domain information
        Real _xMin;
        Real _xMax;
        
        // Test metadata
        bool _hasDiscontinuousDerivative;  // Function has derivative discontinuity
        bool _hasNoise;                     // Data includes random noise
        bool _isOscillatory;               // Function oscillates rapidly
        bool _hasRungePhenomenon;          // Function causes Runge phenomenon with uniform nodes
        int  _suggestedOrder;              // Recommended interpolation order
        Real _expectedMaxError;            // Expected max error for polynomial interp (order = suggestedOrder)

        TestInterpolationData(std::string name, std::string desc,
                              const Vector<Real>& x, const Vector<Real>& y,
                              Real (*trueFunc)(Real), std::string trueExpr,
                              Real xMin, Real xMax,
                              bool discDeriv = false, bool noise = false,
                              bool oscill = false, bool runge = false,
                              int sugOrder = 4, Real expError = 1e-6)
            : _name(name), _description(desc),
              _xData(x), _yData(y),
              _trueFunction(trueFunc), _trueFuncExpr(trueExpr),
              _xMin(xMin), _xMax(xMax),
              _hasDiscontinuousDerivative(discDeriv), _hasNoise(noise),
              _isOscillatory(oscill), _hasRungePhenomenon(runge),
              _suggestedOrder(sugOrder), _expectedMaxError(expError)
        {}
    };

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                         HELPER FUNCTIONS                                               //
    ///////////////////////////////////////////////////////////////////////////////////////////

    // Generate uniformly spaced points
    static Vector<Real> generateUniformNodes(Real xMin, Real xMax, int n)
    {
        Vector<Real> nodes(n);
        Real dx = (xMax - xMin) / (n - 1);
        for (int i = 0; i < n; i++) {
            nodes[i] = xMin + i * dx;
        }
        return nodes;
    }

    // Generate Chebyshev nodes (optimal for polynomial interpolation)
    static Vector<Real> generateChebyshevNodes(Real xMin, Real xMax, int n)
    {
        Vector<Real> nodes(n);
        for (int i = 0; i < n; i++) {
            // Chebyshev nodes of the first kind: cos((2k+1)π/(2n))
            Real theta = static_cast<Real>(Constants::PI) * (2 * i + 1) / (2 * n);
            Real chebNode = std::cos(theta);  // In [-1, 1]
            // Map to [xMin, xMax]
            nodes[i] = static_cast<Real>(0.5) * ((xMax - xMin) * chebNode + xMax + xMin);
        }
        return nodes;
    }

    // Evaluate function at nodes
    static Vector<Real> evaluateAtNodes(const Vector<Real>& nodes, Real (*func)(Real))
    {
        int n = nodes.size();
        Vector<Real> values(n);
        for (int i = 0; i < n; i++) {
            values[i] = func(nodes[i]);
        }
        return values;
    }

    // Add Gaussian noise to values
    static Vector<Real> addNoise(const Vector<Real>& values, Real noiseStdDev, unsigned int seed = 42)
    {
        int n = values.size();
        Vector<Real> noisy(n);
        // Simple deterministic pseudo-noise for reproducibility
        unsigned int state = seed;
        for (int i = 0; i < n; i++) {
            // Linear congruential generator for reproducible "random" noise
            state = state * 1103515245u + 12345u;
            Real u1 = static_cast<Real>((state >> 16) & 0x7FFF) / static_cast<Real>(32767.0);
            state = state * 1103515245u + 12345u;
            Real u2 = static_cast<Real>((state >> 16) & 0x7FFF) / static_cast<Real>(32767.0);
            // Box-Muller transform (simplified)
            Real z = std::sqrt(static_cast<Real>(-2.0) * std::log(u1 + static_cast<Real>(1e-10))) 
                   * std::cos(static_cast<Real>(2.0 * Constants::PI) * u2);
            noisy[i] = values[i] + noiseStdDev * z;
        }
        return noisy;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    //           TEST DATA GENERATORS (using pre-computed static data from _defs.h)          //
    ///////////////////////////////////////////////////////////////////////////////////////////

    // 1. RUNGE FUNCTION - Classic interpolation failure case (uses static data)
    inline TestInterpolationData getRungeTest()
    {
        return TestInterpolationData(
            "Runge Function (Uniform Nodes)",
            "Classic demonstration of Runge phenomenon: polynomial interpolation fails at boundaries",
            runge_uniform_11_x, runge_uniform_11_y,
            Runge_TrueFunc, "1/(1+25x^2)",
            -1.0, 1.0,
            false, false, false, true,  // Runge phenomenon = true
            4, 0.5  // Low order works better; high order diverges
        );
    }

    // 2. RUNGE FUNCTION WITH CHEBYSHEV NODES - Optimal for polynomials (generates dynamically)
    inline TestInterpolationData getRungeChebTest()
    {
        const int n = 11;
        Vector<Real> x = generateChebyshevNodes(-1, 1, n);
        Vector<Real> y = evaluateAtNodes(x, Runge_TrueFunc);

        return TestInterpolationData(
            "Runge Function (Chebyshev Nodes)",
            "Runge function with optimal Chebyshev nodes - eliminates Runge phenomenon",
            x, y,
            Runge_TrueFunc, "1/(1+25x^2)",
            -1.0, 1.0,
            false, false, false, false,  // No Runge phenomenon with Chebyshev nodes!
            10, 1e-3  // Can use high order safely
        );
    }

    // 3. SMOOTH POLYNOMIAL - Exact reproduction test (uses static data)
    inline TestInterpolationData getSmoothPolyTest()
    {
        return TestInterpolationData(
            "Cubic Polynomial",
            "Smooth polynomial should be exactly reproduced by polynomial interpolation",
            smoothpoly_5_x, smoothpoly_5_y,
            SmoothPoly_TrueFunc, "x^3 - 2x^2 + x - 1",
            -2.0, 2.0,
            false, false, false, false,
            4, 1e-12  // Should be machine precision
        );
    }

    // 4. GAUSSIAN - Very smooth transcendental (uses static data)
    inline TestInterpolationData getGaussianTest()
    {
        return TestInterpolationData(
            "Gaussian exp(-x^2)",
            "Very smooth bell curve - tests convergence rate of interpolation methods",
            gaussian_15_x, gaussian_15_y,
            Gaussian_TrueFunc, "exp(-x^2)",
            -3.0, 3.0,
            false, false, false, false,
            6, 1e-6
        );
    }

    // 5. SINE FUNCTION - Smooth periodic (uses static data)
    inline TestInterpolationData getSinePiTest()
    {
        return TestInterpolationData(
            "Sine(pi*x)",
            "Smooth sinusoidal - one complete cycle",
            sinepi_9_x, sinepi_9_y,
            SinePi_TrueFunc, "sin(pi*x)",
            -1.0, 1.0,
            false, false, false, false,
            6, 1e-5
        );
    }

    // 6. HIGH-FREQUENCY OSCILLATION - Needs many points (uses static data)
    inline TestInterpolationData getHighFreqOscillTest()
    {
        return TestInterpolationData(
            "High Frequency sin(5*pi*x)",
            "Rapidly oscillating function - tests interpolation with insufficient sampling",
            highfreq_21_x, highfreq_21_y,
            OscillatoryHighFreq_TrueFunc, "sin(5*pi*x)",
            -1.0, 1.0,
            false, false, true, false,  // Oscillatory = true
            4, 0.1  // Low order to avoid instability
        );
    }

    // 7. ABSOLUTE VALUE - Discontinuous derivative (uses static data)
    inline TestInterpolationData getAbsValueTest()
    {
        return TestInterpolationData(
            "Absolute Value |x|",
            "Discontinuous first derivative at x=0 - polynomial interpolation struggles",
            absval_11_x, absval_11_y,
            AbsValue_TrueFunc, "|x|",
            -2.0, 2.0,
            true, false, false, false,  // Discontinuous derivative = true
            4, 0.1  // Gibbs-like oscillations expected
        );
    }

    // 8. SQUARE ROOT - Singular derivative (uses static data)
    inline TestInterpolationData getSqrtTest()
    {
        return TestInterpolationData(
            "Square Root sqrt(x+1)",
            "Infinite derivative at x=-1 (just outside domain) affects interpolation",
            sqrt_11_x, sqrt_11_y,
            SqrtShift_TrueFunc, "sqrt(x+1)",
            0.0, 4.0,
            true, false, false, false,  // Effectively discontinuous derivative at boundary
            4, 1e-3
        );
    }

    // 9. EXPONENTIAL - Classic smooth function (uses static data)
    inline TestInterpolationData getExponentialTest()
    {
        return TestInterpolationData(
            "Exponential exp(x)",
            "Smooth exponential - tests polynomial approximation of rapidly growing function",
            exp_9_x, exp_9_y,
            Exponential_TrueFunc, "exp(x)",
            -2.0, 2.0,
            false, false, false, false,
            6, 1e-6
        );
    }

    // 10. LOGARITHM - Slower convergence near singularity (uses static data)
    inline TestInterpolationData getLogarithmTest()
    {
        return TestInterpolationData(
            "Logarithm ln(x+2)",
            "Smooth but with singularity at x=-2 - convergence slows near boundary",
            log_15_x, log_15_y,
            LogShift_TrueFunc, "ln(x+2)",
            -1.0, 3.0,
            false, false, false, false,
            6, 1e-4
        );
    }

    // 11. NOISY DATA - Robustness test (generates dynamically with noise)
    inline TestInterpolationData getNoisyGaussianTest()
    {
        const int n = 21;
        Vector<Real> x = generateUniformNodes(-3, 3, n);
        Vector<Real> yClean = evaluateAtNodes(x, Gaussian_TrueFunc);
        Vector<Real> yNoisy = addNoise(yClean, static_cast<Real>(0.02));  // 2% noise

        return TestInterpolationData(
            "Noisy Gaussian",
            "Gaussian with 2% random noise - tests robustness and smoothing",
            x, yNoisy,
            Gaussian_TrueFunc, "exp(-x^2) + noise",
            -3.0, 3.0,
            false, true, false, false,  // Has noise = true
            4, 0.1  // Noise dominates at high orders
        );
    }

    // 12. STEP-LIKE FUNCTION - Very steep transition (uses static data)
    inline TestInterpolationData getStepLikeTest()
    {
        return TestInterpolationData(
            "Step-like tanh(10x)",
            "Very steep but continuous - challenges polynomial interpolation near transition",
            steplike_21_x, steplike_21_y,
            StepLike_TrueFunc, "tanh(10x)",
            -1.0, 1.0,
            false, false, false, true,  // Runge-like behavior expected
            4, 0.2
        );
    }

    // 13. CHEBYSHEV POLYNOMIAL - Exact reproduction test (uses static data)
    inline TestInterpolationData getChebyshev5Test()
    {
        return TestInterpolationData(
            "Chebyshev T_5(x)",
            "Fifth Chebyshev polynomial at Chebyshev nodes - should be exact",
            cheb5_chebNodes_7_x, cheb5_chebNodes_7_y,
            Chebyshev5_TrueFunc, "16x^5 - 20x^3 + 5x",
            -1.0, 1.0,
            false, false, false, false,
            6, 1e-12  // Machine precision expected
        );
    }

    // 14. WITCH OF AGNESI - Milder Runge phenomenon (uses static data)
    inline TestInterpolationData getWitchTest()
    {
        return TestInterpolationData(
            "Witch of Agnesi",
            "Similar to Runge but milder - 1/(1+x^2) on wider interval",
            witch_15_x, witch_15_y,
            WitchOfAgnesi_TrueFunc, "1/(1+x^2)",
            -5.0, 5.0,
            false, false, false, true,  // Mild Runge phenomenon
            6, 0.1
        );
    }

    // 15. MULTI-FREQUENCY - Complex spectral content (uses static data)
    inline TestInterpolationData getMultiFreqTest()
    {
        return TestInterpolationData(
            "Multi-frequency Sine Series",
            "sin(x) + 0.5*sin(3x) + 0.25*sin(5x) - tests interpolation of multi-frequency signals",
            multifreq_25_x, multifreq_25_y,
            MultiFreq_TrueFunc, "sin(x) + 0.5*sin(3x) + 0.25*sin(5x)",
            -static_cast<Real>(Constants::PI), static_cast<Real>(Constants::PI),
            false, false, true, false,  // Oscillatory
            8, 1e-3
        );
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    //           CUSTOM TEST DATA GENERATORS (using arbitrary parameters)                    //
    ///////////////////////////////////////////////////////////////////////////////////////////

    // Create custom test with arbitrary function and parameters
    inline TestInterpolationData createCustomTest(
        const std::string& name,
        const std::string& desc,
        Real (*trueFunc)(Real),
        const std::string& trueExpr,
        Real xMin, Real xMax, int numPoints,
        bool useChebyshevNodes = false,
        bool discDeriv = false, bool oscill = false, bool runge = false,
        int sugOrder = 4, Real expError = 1e-6)
    {
        Vector<Real> x = useChebyshevNodes 
            ? generateChebyshevNodes(xMin, xMax, numPoints)
            : generateUniformNodes(xMin, xMax, numPoints);
        Vector<Real> y = evaluateAtNodes(x, trueFunc);

        return TestInterpolationData(
            name, desc, x, y,
            trueFunc, trueExpr,
            xMin, xMax,
            discDeriv, false, oscill, runge,
            sugOrder, expError
        );
    }

    // Create noisy test from a true function
    inline TestInterpolationData createNoisyTest(
        const std::string& name,
        const std::string& desc,
        Real (*trueFunc)(Real),
        const std::string& trueExpr,
        Real xMin, Real xMax, int numPoints,
        Real noiseStdDev,
        int sugOrder = 4, Real expError = 0.1)
    {
        Vector<Real> x = generateUniformNodes(xMin, xMax, numPoints);
        Vector<Real> yClean = evaluateAtNodes(x, trueFunc);
        Vector<Real> yNoisy = addNoise(yClean, noiseStdDev);

        return TestInterpolationData(
            name, desc, x, yNoisy,
            trueFunc, trueExpr + " + noise",
            xMin, xMax,
            false, true, false, false,
            sugOrder, expError
        );
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                         TEST DATA COLLECTION                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////

    // Get all interpolation test cases
    inline std::vector<TestInterpolationData> getAllInterpolationTests()
    {
        std::vector<TestInterpolationData> tests;
        
        tests.push_back(getRungeTest());
        tests.push_back(getRungeChebTest());
        tests.push_back(getSmoothPolyTest());
        tests.push_back(getGaussianTest());
        tests.push_back(getSinePiTest());
        tests.push_back(getHighFreqOscillTest());
        tests.push_back(getAbsValueTest());
        tests.push_back(getSqrtTest());
        tests.push_back(getExponentialTest());
        tests.push_back(getLogarithmTest());
        tests.push_back(getNoisyGaussianTest());
        tests.push_back(getStepLikeTest());
        tests.push_back(getChebyshev5Test());
        tests.push_back(getWitchTest());
        tests.push_back(getMultiFreqTest());
        
        return tests;
    }

    // Get tests by category
    inline std::vector<TestInterpolationData> getRungePhenomenonTests()
    {
        std::vector<TestInterpolationData> tests;
        tests.push_back(getRungeTest());
        tests.push_back(getRungeChebTest());
        tests.push_back(getWitchTest());
        tests.push_back(getStepLikeTest());
        return tests;
    }

    inline std::vector<TestInterpolationData> getSmoothFunctionTests()
    {
        std::vector<TestInterpolationData> tests;
        tests.push_back(getSmoothPolyTest());
        tests.push_back(getGaussianTest());
        tests.push_back(getSinePiTest());
        tests.push_back(getExponentialTest());
        tests.push_back(getLogarithmTest());
        tests.push_back(getChebyshev5Test());
        return tests;
    }

    inline std::vector<TestInterpolationData> getChallengingTests()
    {
        std::vector<TestInterpolationData> tests;
        tests.push_back(getAbsValueTest());
        tests.push_back(getSqrtTest());
        tests.push_back(getHighFreqOscillTest());
        tests.push_back(getNoisyGaussianTest());
        tests.push_back(getMultiFreqTest());
        return tests;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                    EVALUATION POINTS FOR INTERPOLATION TESTING                         //
    ///////////////////////////////////////////////////////////////////////////////////////////

    // Generate evaluation points (denser than data points)
    inline Vector<Real> generateEvaluationPoints(Real xMin, Real xMax, int n)
    {
        return generateUniformNodes(xMin, xMax, n);
    }

    // Standard evaluation grid (101 points for smooth error curves)
    inline Vector<Real> getStandardEvaluationGrid(Real xMin, Real xMax)
    {
        return generateUniformNodes(xMin, xMax, 101);
    }

    // Compute maximum interpolation error
    template<typename InterpFunc>
    static Real computeMaxError(const InterpFunc& interp, 
                                 Real (*trueFunc)(Real),
                                 const Vector<Real>& evalPoints)
    {
        Real maxError = 0;
        for (int i = 0; i < evalPoints.size(); i++) {
            Real x = evalPoints[i];
            Real interpVal = interp(x);
            Real trueVal = trueFunc(x);
            Real error = std::abs(interpVal - trueVal);
            if (error > maxError) maxError = error;
        }
        return maxError;
    }

    // Compute RMS interpolation error
    template<typename InterpFunc>
    static Real computeRMSError(const InterpFunc& interp, 
                                 Real (*trueFunc)(Real),
                                 const Vector<Real>& evalPoints)
    {
        Real sumSq = 0;
        for (int i = 0; i < evalPoints.size(); i++) {
            Real x = evalPoints[i];
            Real error = interp(x) - trueFunc(x);
            sumSq += error * error;
        }
        return std::sqrt(sumSq / evalPoints.size());
    }

} // namespace MML::TestBeds

#endif // __MML_INTERPOLATION_TEST_BED_H

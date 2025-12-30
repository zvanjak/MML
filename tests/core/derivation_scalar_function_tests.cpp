#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/VectorN.h"
#include "base/Function.h"
#include "core/Derivation/DerivationScalarFunction.h"
#include "../test_data/scalar_functions_test_bed.h"
#endif

using namespace MML;
using namespace MML::Testing;
using namespace MML::Derivation;
using namespace MML::TestBeds;
using namespace Catch::Matchers;

/**
 * COMPREHENSIVE TEST SUITE: Scalar Function Derivatives
 * 
 * Part of Epic 2 (MinimalMathLibrary-ao0) - Task REAL(2.2) (MinimalMathLibrary-1fd)
 * Created after discovering catastrophic NDer2_uw bug in parametric surfaces.
 * 
 * PURPOSE:
 * Verify ALL numerical derivative formulas for scalar functions (R^N → R):
 * - NDer1/2/4/6/8 Partial: Single partial derivatives with different accuracy orders
 * - NDer1/2/4/6/8 PartialByAll: Gradient (all partials at once)
 * 
 * TEST STRATEGY:
 * - Use TestScalarFunc1 and TestScalarFunc2 with known analytical partial derivatives
 * - Test against analytical gradient components
 * - Verify increasing accuracy with higher-order methods
 * - Focus on correctness over extensive precision testing
 */

namespace   // Anonymous namespace for test helpers
{
    ///////////////////////////////////////////////////////////////////////////
    // TEST SCALAR FUNCTIONS WITH KNOWN ANALYTICAL DERIVATIVES
    ///////////////////////////////////////////////////////////////////////////
    
    // Function 1: f(x,y,z) = cos(x) + sin(y) + exp(z)
    // ∂f/∂x = -sin(x)
    // ∂f/∂y = cos(y)
    // ∂f/∂z = exp(z)
    
    // Function 2: f(x,y,z) = sin(xy)*exp(z/(y²+1))/(1+x²)
    // Much more complex with analytical derivatives provided in test_data
    
    static const TestFunctionScalar<3>& GetSimpleFunc()
    {
        return ScalarFunctionsTestBed::getTestFunctionScalar3(0);
    }
    
    static const TestFunctionScalar<3>& GetComplexFunc()
    {
        return ScalarFunctionsTestBed::getTestFunctionScalar3(1);
    }
    
    // Helper: Get analytical gradient at a point
    static VectorN<Real, 3> GetAnalyticalGradient(const TestFunctionScalar<3>& testFunc, const VectorN<Real, 3>& point)
    {
        VectorN<Real, 3> grad;
        for (int i = 0; i < 3; i++)
        {
            grad[i] = testFunc._funcDerived(point, i);
        }
        return grad;
    }

    ///////////////////////////////////////////////////////////////////////////
    // SINGLE PARTIAL DERIVATIVE: NDer1/2/4/6/8 Partial
    ///////////////////////////////////////////////////////////////////////////
    
    TEST_CASE("ScalarFunction - NDer1Partial (Forward Difference)", "[scalar_function][nder1][partial]")
    {
        auto func = GetSimpleFunc()._func;
        VectorN<Real, 3> point{REAL(1.0), REAL(0.5), REAL(2.0)};
        
        VectorN<Real, 3> grad_analytical = GetAnalyticalGradient(GetSimpleFunc(), point);
        
        // Test each partial derivative
        for (int deriv_var = 0; deriv_var < 3; deriv_var++)
        {
            Real numerical = NDer1Partial(func, deriv_var, point);
            Real analytical = grad_analytical[deriv_var];
            
            INFO("Partial derivative w.r.t. x[" << deriv_var << "]");
            REQUIRE_THAT(numerical, WithinAbs(analytical, REAL(1e-5)));  // Forward difference less accurate
        }
    }
    
    TEST_CASE("ScalarFunction - NDer2Partial (Central Difference)", "[scalar_function][nder2][partial]")
    {
        auto func = GetSimpleFunc()._func;
        VectorN<Real, 3> point{REAL(1.0), REAL(0.5), REAL(2.0)};
        
        VectorN<Real, 3> grad_analytical = GetAnalyticalGradient(GetSimpleFunc(), point);
        
        for (int deriv_var = 0; deriv_var < 3; deriv_var++)
        {
            Real numerical = NDer2Partial(func, deriv_var, point);
            Real analytical = grad_analytical[deriv_var];
            
            INFO("Partial derivative w.r.t. x[" << deriv_var << "]");
            REQUIRE_THAT(numerical, WithinAbs(analytical, REAL(1e-8)));  // Central difference better
        }
    }
    
    TEST_CASE("ScalarFunction - NDer4Partial (4th Order)", "[scalar_function][nder4][partial]")
    {
        auto func = GetSimpleFunc()._func;
        VectorN<Real, 3> point{REAL(1.0), REAL(0.5), REAL(2.0)};
        
        VectorN<Real, 3> grad_analytical = GetAnalyticalGradient(GetSimpleFunc(), point);
        
        for (int deriv_var = 0; deriv_var < 3; deriv_var++)
        {
            Real numerical = NDer4Partial(func, deriv_var, point);
            Real analytical = grad_analytical[deriv_var];
            
            INFO("Partial derivative w.r.t. x[" << deriv_var << "]");
            REQUIRE_THAT(numerical, WithinAbs(analytical, REAL(1e-10)));  // 4th order very accurate
        }
    }
    
    TEST_CASE("ScalarFunction - NDer6Partial (6th Order)", "[scalar_function][nder6][partial]")
    {
        auto func = GetSimpleFunc()._func;
        VectorN<Real, 3> point{REAL(1.0), REAL(0.5), REAL(2.0)};
        
        VectorN<Real, 3> grad_analytical = GetAnalyticalGradient(GetSimpleFunc(), point);
        
        for (int deriv_var = 0; deriv_var < 3; deriv_var++)
        {
            Real numerical = NDer6Partial(func, deriv_var, point);
            Real analytical = grad_analytical[deriv_var];
            
            INFO("Partial derivative w.r.t. x[" << deriv_var << "]");
            REQUIRE_THAT(numerical, WithinAbs(analytical, REAL(1e-11)));  // 6th order excellent
        }
    }
    
    TEST_CASE("ScalarFunction - NDer8Partial (8th Order)", "[scalar_function][nder8][partial]")
    {
        auto func = GetSimpleFunc()._func;
        VectorN<Real, 3> point{REAL(1.0), REAL(0.5), REAL(2.0)};
        
        VectorN<Real, 3> grad_analytical = GetAnalyticalGradient(GetSimpleFunc(), point);
        
        for (int deriv_var = 0; deriv_var < 3; deriv_var++)
        {
            Real numerical = NDer8Partial(func, deriv_var, point);
            Real analytical = grad_analytical[deriv_var];
            
            INFO("Partial derivative w.r.t. x[" << deriv_var << "]");
            REQUIRE_THAT(numerical, WithinAbs(analytical, REAL(1e-12)));  // 8th order best
        }
    }

    ///////////////////////////////////////////////////////////////////////////
    // GRADIENT (ALL PARTIALS): NDer1/2/4/6/8 PartialByAll
    ///////////////////////////////////////////////////////////////////////////
    
    TEST_CASE("ScalarFunction - NDer2PartialByAll (Gradient)", "[scalar_function][nder2][gradient]")
    {
        auto func = GetSimpleFunc()._func;
        VectorN<Real, 3> point{REAL(1.0), REAL(0.5), REAL(2.0)};
        
        VectorN<Real, 3> grad_numerical = NDer2PartialByAll(func, point);
        VectorN<Real, 3> grad_analytical = GetAnalyticalGradient(GetSimpleFunc(), point);
        
        for (int i = 0; i < 3; i++)
        {
            INFO("Gradient component [" << i << "]");
            REQUIRE_THAT(grad_numerical[i], WithinAbs(grad_analytical[i], REAL(1e-8)));
        }
    }
    
    TEST_CASE("ScalarFunction - NDer4PartialByAll (Gradient)", "[scalar_function][nder4][gradient]")
    {
        auto func = GetSimpleFunc()._func;
        VectorN<Real, 3> point{REAL(1.0), REAL(0.5), REAL(2.0)};
        
        VectorN<Real, 3> grad_numerical = NDer4PartialByAll(func, point);
        VectorN<Real, 3> grad_analytical = GetAnalyticalGradient(GetSimpleFunc(), point);
        
        for (int i = 0; i < 3; i++)
        {
            INFO("Gradient component [" << i << "]");
            REQUIRE_THAT(grad_numerical[i], WithinAbs(grad_analytical[i], REAL(1e-10)));
        }
    }
    
    TEST_CASE("ScalarFunction - NDer6PartialByAll (Gradient)", "[scalar_function][nder6][gradient]")
    {
        auto func = GetSimpleFunc()._func;
        VectorN<Real, 3> point{REAL(1.0), REAL(0.5), REAL(2.0)};
        
        VectorN<Real, 3> grad_numerical = NDer6PartialByAll(func, point);
        VectorN<Real, 3> grad_analytical = GetAnalyticalGradient(GetSimpleFunc(), point);
        
        for (int i = 0; i < 3; i++)
        {
            INFO("Gradient component [" << i << "]");
            REQUIRE_THAT(grad_numerical[i], WithinAbs(grad_analytical[i], REAL(1e-11)));
        }
    }
    
    TEST_CASE("ScalarFunction - NDer8PartialByAll (Gradient)", "[scalar_function][nder8][gradient]")
    {
        auto func = GetSimpleFunc()._func;
        VectorN<Real, 3> point{REAL(1.0), REAL(0.5), REAL(2.0)};
        
        VectorN<Real, 3> grad_numerical = NDer8PartialByAll(func, point);
        VectorN<Real, 3> grad_analytical = GetAnalyticalGradient(GetSimpleFunc(), point);
        
        for (int i = 0; i < 3; i++)
        {
            INFO("Gradient component [" << i << "]");
            REQUIRE_THAT(grad_numerical[i], WithinAbs(grad_analytical[i], REAL(1e-12)));
        }
    }

    ///////////////////////////////////////////////////////////////////////////
    // COMPLEX FUNCTION TESTS: Verify on more challenging function
    ///////////////////////////////////////////////////////////////////////////
    
    TEST_CASE("ScalarFunction - Complex Function NDer4Partial", "[scalar_function][nder4][complex]")
    {
        auto func = GetComplexFunc()._func;
        VectorN<Real, 3> point{REAL(1.0), REAL(0.5), REAL(2.0)};
        
        VectorN<Real, 3> grad_analytical = GetAnalyticalGradient(GetComplexFunc(), point);
        
        // Test all partial derivatives on complex function
        for (int deriv_var = 0; deriv_var < 3; deriv_var++)
        {
            Real numerical = NDer4Partial(func, deriv_var, point);
            Real analytical = grad_analytical[deriv_var];
            
            INFO("Complex function, partial w.r.t. x[" << deriv_var << "]");
            INFO("Analytical value: " << analytical);
            INFO("Numerical value: " << numerical);
            REQUIRE_THAT(numerical, WithinAbs(analytical, REAL(1e-7)));  // Complex function, slightly relaxed
        }
    }
    
    TEST_CASE("ScalarFunction - Complex Function NDer4PartialByAll", "[scalar_function][nder4][complex][gradient]")
    {
        auto func = GetComplexFunc()._func;
        VectorN<Real, 3> point{REAL(1.0), REAL(0.5), REAL(2.0)};
        
        VectorN<Real, 3> grad_numerical = NDer4PartialByAll(func, point);
        VectorN<Real, 3> grad_analytical = GetAnalyticalGradient(GetComplexFunc(), point);
        
        for (int i = 0; i < 3; i++)
        {
            INFO("Complex function gradient component [" << i << "]");
            INFO("Analytical value: " << grad_analytical[i]);
            INFO("Numerical value: " << grad_numerical[i]);
            REQUIRE_THAT(grad_numerical[i], WithinAbs(grad_analytical[i], REAL(1e-7)));
        }
    }

    ///////////////////////////////////////////////////////////////////////////
    // ACCURACY COMPARISON: Verify higher-order methods are more accurate
    ///////////////////////////////////////////////////////////////////////////
    
    TEST_CASE("ScalarFunction - Accuracy Progression", "[scalar_function][accuracy]")
    {
        auto func = GetSimpleFunc()._func;
        VectorN<Real, 3> point{REAL(1.0), REAL(0.5), REAL(2.0)};
        
        // Test on ∂f/∂x (first component, -sin(x))
        int deriv_var = 0;
        Real analytical = GetAnalyticalGradient(GetSimpleFunc(), point)[deriv_var];
        
        Real error1 = std::abs(NDer1Partial(func, deriv_var, point) - analytical);
        Real error2 = std::abs(NDer2Partial(func, deriv_var, point) - analytical);
        Real error4 = std::abs(NDer4Partial(func, deriv_var, point) - analytical);
        Real error6 = std::abs(NDer6Partial(func, deriv_var, point) - analytical);
        Real error8 = std::abs(NDer8Partial(func, deriv_var, point) - analytical);
        
        INFO("NDer1 error: " << error1);
        INFO("NDer2 error: " << error2);
        INFO("NDer4 error: " << error4);
        INFO("NDer6 error: " << error6);
        INFO("NDer8 error: " << error8);
        
        // Verify monotonic decrease in error
        REQUIRE(error2 < error1);  // Central difference better than forward
        REQUIRE(error4 < error2);  // 4th order better than 2nd
        REQUIRE(error6 < error4);  // 6th order better than 4th
        // NDer8 may not always be better due to roundoff, just verify it's accurate
        REQUIRE(error8 < 1e-11);   // 8th order still very accurate
    }

    ///////////////////////////////////////////////////////////////////////////
    // MULTIPLE TEST POINTS: Verify correctness across domain
    ///////////////////////////////////////////////////////////////////////////
    
    TEST_CASE("ScalarFunction - Multiple Test Points", "[scalar_function][robustness]")
    {
        auto func = GetSimpleFunc()._func;
        
        // Test at various points in the domain
        std::vector<VectorN<Real, 3>> test_points = {
            {REAL(0.5), REAL(0.3), REAL(1.0)},
            {REAL(1.0), REAL(0.5), REAL(2.0)},
            {REAL(1.5), REAL(1.0), REAL(1.5)},
            {REAL(0.1), REAL(0.1), REAL(0.5)},
            {REAL(2.0), REAL(1.5), REAL(3.0)}
        };
        
        for (const auto& point : test_points)
        {
            VectorN<Real, 3> grad_numerical = NDer4PartialByAll(func, point);
            VectorN<Real, 3> grad_analytical = GetAnalyticalGradient(GetSimpleFunc(), point);
            
            INFO("Testing at point (" << point[0] << ", " << point[1] << ", " << point[2] << ")");
            
            for (int i = 0; i < 3; i++)
            {
                REQUIRE_THAT(grad_numerical[i], WithinAbs(grad_analytical[i], REAL(1e-10)));
            }
        }
    }
    
    ///////////////////////////////////////////////////////////////////////////
    // GRADIENT CONSISTENCY: Verify PartialByAll equals individual Partials
    ///////////////////////////////////////////////////////////////////////////
    
    TEST_CASE("ScalarFunction - Gradient Consistency", "[scalar_function][consistency]")
    {
        auto func = GetComplexFunc()._func;
        VectorN<Real, 3> point{REAL(1.0), REAL(0.5), REAL(2.0)};
        
        // Get gradient using PartialByAll
        VectorN<Real, 3> grad_all = NDer4PartialByAll(func, point);
        
        // Get gradient components individually
        VectorN<Real, 3> grad_individual;
        for (int i = 0; i < 3; i++)
        {
            grad_individual[i] = NDer4Partial(func, i, point);
        }
        
        // They should be identical
        for (int i = 0; i < 3; i++)
        {
            INFO("Component [" << i << "]");
            REQUIRE_THAT(grad_all[i], WithinAbs(grad_individual[i], REAL(1e-14)));
        }
    }
}

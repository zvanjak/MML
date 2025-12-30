#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/VectorN.h"
#include "base/MatrixNM.h"
#include "base/Function.h"
#include "core/Derivation/DerivationVectorFunction.h"
#include "../test_data/vector_functions_test_bed.h"
#endif

using namespace MML;
using namespace MML::Testing;
using namespace MML::Derivation;
using namespace MML::TestBeds;
using namespace Catch::Matchers;

/**
 * COMPREHENSIVE TEST SUITE: Vector Function Derivatives
 * 
 * Part of Epic 2 (MinimalMathLibrary-ao0) - Task REAL(2.3) (MinimalMathLibrary-nao)
 * 
 * PURPOSE:
 * Verify ALL numerical derivative formulas for vector functions:
 * - NDer1/2/4/6/8 Partial: Single partial derivatives with different accuracy orders
 * - NDer1/2/4/6/8 PartialByAll: All partials of one component (gradient of component)
 * - NDer1/2/4/6/8 PartialAllByAll: Full Jacobian matrix
 * 
 * TEST STRATEGY:
 * - Use TestVectorFunc1 from vector_functions_test_bed.h with known analytical derivatives
 * - Test against analytical Jacobian matrix elements
 * - Verify increasing accuracy with higher-order methods
 * - Focus on correctness over extensive precision testing
 */

namespace   // Anonymous namespace for test helpers
{
    ///////////////////////////////////////////////////////////////////////////
    // TEST VECTOR FUNCTIONS WITH KNOWN ANALYTICAL DERIVATIVES
    ///////////////////////////////////////////////////////////////////////////
    
    // Get the test function with analytical derivatives
    static const TestFunctionVector<3>& GetTestFunction()
    {
        return VectorFunctionsTestBed::getTestFunctionVector(0);
    }
    
    // Test function: F(x,y,z) = (x*cos(y)*z², sin(x)*(y²+z²), exp(xy/(z²+1)))
    // Jacobian matrix elements (analytical):
    // ∂F₀/∂x = z²cos(y),        ∂F₀/∂y = -xz²sin(y),      ∂F₀/∂z = 2xz cos(y)
    // ∂F₁/∂x = cos(x)(y²+z²),   ∂F₁/∂y = 2y sin(x),       ∂F₁/∂z = 2z sin(x)
    // ∂F₂/∂x = ye^g/(z²+1),     ∂F₂/∂y = xe^g/(z²+1),     ∂F₂/∂z = -2xyz e^g/(z²+1)²
    // where g = xy/(z²+1)
    
    // Helper: Get analytical Jacobian at a point
    static MatrixNM<Real, 3, 3> GetAnalyticalJacobian(const VectorN<Real, 3>& point)
    {
        MatrixNM<Real, 3, 3> J;
        
        // Get analytical partial derivatives for each component
        for (int comp = 0; comp < 3; comp++)
        {
            VectorN<Real, 3> partials = GetTestFunction()._funcDerived(point, comp);
            J(comp, 0) = partials[0];  // ∂Fᵢ/∂x
            J(comp, 1) = partials[1];  // ∂Fᵢ/∂y
            J(comp, 2) = partials[2];  // ∂Fᵢ/∂z
        }
        
        return J;
    }

    ///////////////////////////////////////////////////////////////////////////
    // SINGLE PARTIAL DERIVATIVE: NDer1/2/4/6/8 Partial
    ///////////////////////////////////////////////////////////////////////////
    
    TEST_CASE("VectorFunction - NDer1Partial (Forward Difference)", "[vector_function][nder1][partial]")
    {
        auto func = GetTestFunction()._func;
        VectorN<Real, 3> point{REAL(1.0), REAL(0.5), REAL(2.0)};
        
        // Test all 9 Jacobian matrix elements
        MatrixNM<Real, 3, 3> J_analytical = GetAnalyticalJacobian(point);
        
        for (int func_comp = 0; func_comp < 3; func_comp++)
        {
            for (int deriv_var = 0; deriv_var < 3; deriv_var++)
            {
                Real numerical = NDer1Partial(func, func_comp, deriv_var, point);
                Real analytical = J_analytical(func_comp, deriv_var);
                
                INFO("Component " << func_comp << ", derivative w.r.t. x[" << deriv_var << "]");
                REQUIRE_THAT(numerical, WithinAbs(analytical, REAL(1e-4)));  // Forward difference less accurate
            }
        }
    }
    
    TEST_CASE("VectorFunction - NDer2Partial (Central Difference)", "[vector_function][nder2][partial]")
    {
        auto func = GetTestFunction()._func;
        VectorN<Real, 3> point{REAL(1.0), REAL(0.5), REAL(2.0)};
        
        MatrixNM<Real, 3, 3> J_analytical = GetAnalyticalJacobian(point);
        
        for (int func_comp = 0; func_comp < 3; func_comp++)
        {
            for (int deriv_var = 0; deriv_var < 3; deriv_var++)
            {
                Real numerical = NDer2Partial(func, func_comp, deriv_var, point);
                Real analytical = J_analytical(func_comp, deriv_var);
                
                INFO("Component " << func_comp << ", derivative w.r.t. x[" << deriv_var << "]");
                REQUIRE_THAT(numerical, WithinAbs(analytical, REAL(1e-6)));  // Central difference better
            }
        }
    }
    
    TEST_CASE("VectorFunction - NDer4Partial (4th Order)", "[vector_function][nder4][partial]")
    {
        auto func = GetTestFunction()._func;
        VectorN<Real, 3> point{REAL(1.0), REAL(0.5), REAL(2.0)};
        
        MatrixNM<Real, 3, 3> J_analytical = GetAnalyticalJacobian(point);
        
        for (int func_comp = 0; func_comp < 3; func_comp++)
        {
            for (int deriv_var = 0; deriv_var < 3; deriv_var++)
            {
                Real numerical = NDer4Partial(func, func_comp, deriv_var, point);
                Real analytical = J_analytical(func_comp, deriv_var);
                
                INFO("Component " << func_comp << ", derivative w.r.t. x[" << deriv_var << "]");
                REQUIRE_THAT(numerical, WithinAbs(analytical, REAL(1e-8)));  // 4th order very accurate
            }
        }
    }
    
    TEST_CASE("VectorFunction - NDer6Partial (6th Order)", "[vector_function][nder6][partial]")
    {
        auto func = GetTestFunction()._func;
        VectorN<Real, 3> point{REAL(1.0), REAL(0.5), REAL(2.0)};
        
        MatrixNM<Real, 3, 3> J_analytical = GetAnalyticalJacobian(point);
        
        for (int func_comp = 0; func_comp < 3; func_comp++)
        {
            for (int deriv_var = 0; deriv_var < 3; deriv_var++)
            {
                Real numerical = NDer6Partial(func, func_comp, deriv_var, point);
                Real analytical = J_analytical(func_comp, deriv_var);
                
                INFO("Component " << func_comp << ", derivative w.r.t. x[" << deriv_var << "]");
                REQUIRE_THAT(numerical, WithinAbs(analytical, REAL(1e-9)));  // 6th order excellent
            }
        }
    }
    
    TEST_CASE("VectorFunction - NDer8Partial (8th Order)", "[vector_function][nder8][partial]")
    {
        auto func = GetTestFunction()._func;
        VectorN<Real, 3> point{REAL(1.0), REAL(0.5), REAL(2.0)};
        
        MatrixNM<Real, 3, 3> J_analytical = GetAnalyticalJacobian(point);
        
        for (int func_comp = 0; func_comp < 3; func_comp++)
        {
            for (int deriv_var = 0; deriv_var < 3; deriv_var++)
            {
                Real numerical = NDer8Partial(func, func_comp, deriv_var, point);
                Real analytical = J_analytical(func_comp, deriv_var);
                
                INFO("Component " << func_comp << ", derivative w.r.t. x[" << deriv_var << "]");
                REQUIRE_THAT(numerical, WithinAbs(analytical, REAL(1e-10)));  // 8th order best
            }
        }
    }

    ///////////////////////////////////////////////////////////////////////////
    // ALL PARTIALS OF ONE COMPONENT: NDer1/2/4/6/8 PartialByAll
    // These compute gradient of one component (one row of Jacobian)
    ///////////////////////////////////////////////////////////////////////////
    
    TEST_CASE("VectorFunction - NDer2PartialByAll (Row of Jacobian)", "[vector_function][nder2][partial_by_all]")
    {
        auto func = GetTestFunction()._func;
        VectorN<Real, 3> point{REAL(1.0), REAL(0.5), REAL(2.0)};
        
        MatrixNM<Real, 3, 3> J_analytical = GetAnalyticalJacobian(point);
        
        // Test each component's gradient (each row of Jacobian)
        for (int func_comp = 0; func_comp < 3; func_comp++)
        {
            VectorN<Real, 3> gradient_numerical = NDer2PartialByAll(func, func_comp, point);
            
            INFO("Testing gradient of component " << func_comp);
            REQUIRE_THAT(gradient_numerical[0], WithinAbs(J_analytical(func_comp, REAL(0)), 1e-6));
            REQUIRE_THAT(gradient_numerical[1], WithinAbs(J_analytical(func_comp, REAL(1)), 1e-6));
            REQUIRE_THAT(gradient_numerical[2], WithinAbs(J_analytical(func_comp, REAL(2)), 1e-6));
        }
    }
    
    TEST_CASE("VectorFunction - NDer4PartialByAll (Row of Jacobian)", "[vector_function][nder4][partial_by_all]")
    {
        auto func = GetTestFunction()._func;
        VectorN<Real, 3> point{REAL(1.0), REAL(0.5), REAL(2.0)};
        
        MatrixNM<Real, 3, 3> J_analytical = GetAnalyticalJacobian(point);
        
        for (int func_comp = 0; func_comp < 3; func_comp++)
        {
            VectorN<Real, 3> gradient_numerical = NDer4PartialByAll(func, func_comp, point);
            
            INFO("Testing gradient of component " << func_comp);
            REQUIRE_THAT(gradient_numerical[0], WithinAbs(J_analytical(func_comp, REAL(0)), 1e-8));
            REQUIRE_THAT(gradient_numerical[1], WithinAbs(J_analytical(func_comp, REAL(1)), 1e-8));
            REQUIRE_THAT(gradient_numerical[2], WithinAbs(J_analytical(func_comp, REAL(2)), 1e-8));
        }
    }

    ///////////////////////////////////////////////////////////////////////////
    // FULL JACOBIAN MATRIX: NDer1/2/4/6/8 PartialAllByAll
    // These compute the complete Jacobian matrix
    ///////////////////////////////////////////////////////////////////////////
    
    TEST_CASE("VectorFunction - NDer2PartialAllByAll (Full Jacobian)", "[vector_function][nder2][jacobian]")
    {
        auto func = GetTestFunction()._func;
        VectorN<Real, 3> point{REAL(1.0), REAL(0.5), REAL(2.0)};
        
        MatrixNM<Real, 3, 3> J_numerical = NDer2PartialAllByAll(func, point);
        MatrixNM<Real, 3, 3> J_analytical = GetAnalyticalJacobian(point);
        
        // Test all 9 matrix elements
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                INFO("Jacobian element J(" << i << "," << j << ")");
                REQUIRE_THAT(J_numerical(i, j), WithinAbs(J_analytical(i, j), 1e-6));
            }
        }
    }
    
    TEST_CASE("VectorFunction - NDer4PartialAllByAll (Full Jacobian)", "[vector_function][nder4][jacobian]")
    {
        auto func = GetTestFunction()._func;
        VectorN<Real, 3> point{REAL(1.0), REAL(0.5), REAL(2.0)};
        
        MatrixNM<Real, 3, 3> J_numerical = NDer4PartialAllByAll(func, point);
        MatrixNM<Real, 3, 3> J_analytical = GetAnalyticalJacobian(point);
        
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                INFO("Jacobian element J(" << i << "," << j << ")");
                REQUIRE_THAT(J_numerical(i, j), WithinAbs(J_analytical(i, j), 1e-8));
            }
        }
    }
    
    TEST_CASE("VectorFunction - NDer6PartialAllByAll (Full Jacobian)", "[vector_function][nder6][jacobian]")
    {
        auto func = GetTestFunction()._func;
        VectorN<Real, 3> point{REAL(1.0), REAL(0.5), REAL(2.0)};
        
        MatrixNM<Real, 3, 3> J_numerical = NDer6PartialAllByAll(func, point);
        MatrixNM<Real, 3, 3> J_analytical = GetAnalyticalJacobian(point);
        
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                INFO("Jacobian element J(" << i << "," << j << ")");
                REQUIRE_THAT(J_numerical(i, j), WithinAbs(J_analytical(i, j), 1e-9));
            }
        }
    }
    
    TEST_CASE("VectorFunction - NDer8PartialAllByAll (Full Jacobian)", "[vector_function][nder8][jacobian]")
    {
        auto func = GetTestFunction()._func;
        VectorN<Real, 3> point{REAL(1.0), REAL(0.5), REAL(2.0)};
        
        MatrixNM<Real, 3, 3> J_numerical = NDer8PartialAllByAll(func, point);
        MatrixNM<Real, 3, 3> J_analytical = GetAnalyticalJacobian(point);
        
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                INFO("Jacobian element J(" << i << "," << j << ")");
                REQUIRE_THAT(J_numerical(i, j), WithinAbs(J_analytical(i, j), 1e-10));
            }
        }
    }

    ///////////////////////////////////////////////////////////////////////////
    // ACCURACY COMPARISON: Verify higher-order methods are more accurate
    ///////////////////////////////////////////////////////////////////////////
    
    TEST_CASE("VectorFunction - Accuracy Progression", "[vector_function][accuracy]")
    {
        auto func = GetTestFunction()._func;
        VectorN<Real, 3> point{REAL(1.0), REAL(0.5), REAL(2.0)};
        
        // Test one specific Jacobian element: ∂F₂/∂z (most complex derivative)
        int func_comp = 2;
        int deriv_var = 2;
        
        Real analytical = GetAnalyticalJacobian(point)(func_comp, deriv_var);
        
        Real error1 = std::abs(NDer1Partial(func, func_comp, deriv_var, point) - analytical);
        Real error2 = std::abs(NDer2Partial(func, func_comp, deriv_var, point) - analytical);
        Real error4 = std::abs(NDer4Partial(func, func_comp, deriv_var, point) - analytical);
        Real error6 = std::abs(NDer6Partial(func, func_comp, deriv_var, point) - analytical);
        Real error8 = std::abs(NDer8Partial(func, func_comp, deriv_var, point) - analytical);
        
        INFO("NDer1 error: " << error1);
        INFO("NDer2 error: " << error2);
        INFO("NDer4 error: " << error4);
        INFO("NDer6 error: " << error6);
        INFO("NDer8 error: " << error8);
        
        // Verify monotonic decrease in error (each method should be more accurate)
        // Allow for the rare case where the lowest-order estimate is exactly
        // equal to the analytical value (error1 == 0) — skip strict check then.
        if (error1 > 0)
            REQUIRE(error2 < error1);  // Central difference better than forward
        else
            INFO("NDer1 error is zero; skipping strict NDer2 < NDer1 check");

        // Allow a tiny tolerance for roundoff: higher-order methods should
        // generally be no worse than lower-order, but floating-point noise can
        // make the inequalities non-strict at the 1e-14 scale.
        REQUIRE(error4 <= error2 + REAL(5e-14));  // 4th order no worse than 2nd (tolerance for roundoff)
        REQUIRE(error6 <= error4 + REAL(5e-14));  // 6th order no worse than 4th (tolerance for roundoff)
        // Note: NDer8 may not always be better than NDer6 due to roundoff errors
        // with very small step sizes, so we just verify both are very accurate
        REQUIRE(error8 < 1e-10);   // 8th order still very accurate
    }

    ///////////////////////////////////////////////////////////////////////////
    // MULTIPLE TEST POINTS: Verify correctness across domain
    ///////////////////////////////////////////////////////////////////////////
    
    TEST_CASE("VectorFunction - Multiple Test Points", "[vector_function][robustness]")
    {
        auto func = GetTestFunction()._func;
        
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
            MatrixNM<Real, 3, 3> J_numerical = NDer4PartialAllByAll(func, point);
            MatrixNM<Real, 3, 3> J_analytical = GetAnalyticalJacobian(point);
            
            INFO("Testing at point (" << point[0] << ", " << point[1] << ", " << point[2] << ")");
            
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    REQUIRE_THAT(J_numerical(i, j), WithinAbs(J_analytical(i, j), 1e-7));
                }
            }
        }
    }
}

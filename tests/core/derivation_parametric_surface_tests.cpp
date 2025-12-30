#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/VectorN.h"
#include "base/Function.h"
#include "core/Derivation/DerivationParametricSurface.h"
#include "core/Surfaces.h"
#endif

using namespace MML;
using namespace MML::Testing;
using namespace MML::Derivation;
using namespace MML::Surfaces;
using namespace Catch::Matchers;

/**
 * CRITICAL TEST SUITE: Parametric Surface Derivatives
 * 
 * This test suite was created after discovering a CATASTROPHIC bug in NDer2_uw
 * that caused all differential geometry computations to be incorrect.
 * 
 * BUG FOUND: NDer2_uw used formula (f(u+h,w+h) - 2f(u,w) + f(u-h,w-h))/h²
 * CORRECT:   NDer2_uw should use (f(u+h,w+h) - f(u+h,w-h) - f(u-h,w+h) + f(u-h,w-h))/(4h²)
 * 
 * This suite verifies ALL parametric surface derivative formulas against
 * analytical derivatives computed from known surfaces.
 */

namespace MML::Tests::Core::DerivationParametricSurfaceTests
{
    ///////////////////////////////////////////////////////////////////////////
    // TEST SURFACES WITH KNOWN ANALYTICAL DERIVATIVES
    ///////////////////////////////////////////////////////////////////////////
    
    // Surface 1: Simple plane z = x + 2y
    // r(u,w) = (u, w, u + 2w)
    // r_u = (1, 0, 1), r_w = (0, 1, 2)
    // r_uu = (0,0,0), r_uw = (0,0,0), r_ww = (0,0,0)
    class SimplePlane : public IParametricSurfaceRect<3>
    {
    public:
        Real getMinU() const { return -10; }
        Real getMaxU() const { return 10; }
        Real getMinW() const { return -10; }
        Real getMaxW() const { return 10; }
        
        VectorN<Real, 3> operator()(Real u, Real w) const {
            return VectorN<Real, 3>{u, w, u + 2*w};
        }
    };
    
    // Surface 2: Sphere r(u,w) = (R*sin(u)*cos(w), R*sin(u)*sin(w), R*cos(u))
    // Known derivatives:
    // r_u = (R*cos(u)*cos(w), R*cos(u)*sin(w), -R*sin(u))
    // r_w = (-R*sin(u)*sin(w), R*sin(u)*cos(w), 0)
    // r_uu = (-R*sin(u)*cos(w), -R*sin(u)*sin(w), -R*cos(u))
    // r_ww = (-R*sin(u)*cos(w), -R*sin(u)*sin(w), 0)
    // r_uw = (-R*cos(u)*sin(w), R*cos(u)*cos(w), 0)
    
    // Surface 3: Paraboloid z = x² + y²
    // r(u,w) = (u, w, u² + w²)
    // r_u = (1, 0, 2u), r_w = (0, 1, 2w)
    // r_uu = (0, 0, 2), r_ww = (0, 0, 2), r_uw = (0, 0, 0)
    class Paraboloid : public IParametricSurfaceRect<3>
    {
    public:
        Real getMinU() const { return -5; }
        Real getMaxU() const { return 5; }
        Real getMinW() const { return -5; }
        Real getMaxW() const { return 5; }
        
        VectorN<Real, 3> operator()(Real u, Real w) const {
            return VectorN<Real, 3>{u, w, u*u + w*w};
        }
    };
    
    // Surface 4: Hyperbolic Paraboloid z = x² - y² (saddle)
    // r(u,w) = (u, w, u² - w²)
    // r_u = (1, 0, 2u), r_w = (0, 1, -2w)
    // r_uu = (0, 0, 2), r_ww = (0, 0, -2), r_uw = (0, 0, 0)
    class HyperbolicParaboloid : public IParametricSurfaceRect<3>
    {
    public:
        Real getMinU() const { return -5; }
        Real getMaxU() const { return 5; }
        Real getMinW() const { return -5; }
        Real getMaxW() const { return 5; }
        
        VectorN<Real, 3> operator()(Real u, Real w) const {
            return VectorN<Real, 3>{u, w, u*u - w*w};
        }
    };
    
    // Surface 5: Mixed product surface z = xy
    // r(u,w) = (u, w, uw)
    // r_u = (1, 0, w), r_w = (0, 1, u)
    // r_uu = (0, 0, 0), r_ww = (0, 0, 0), r_uw = (0, 0, 1) ← CRITICAL TEST!
    class MixedProductSurface : public IParametricSurfaceRect<3>
    {
    public:
        Real getMinU() const { return -5; }
        Real getMaxU() const { return 5; }
        Real getMinW() const { return -5; }
        Real getMaxW() const { return 5; }
        
        VectorN<Real, 3> operator()(Real u, Real w) const {
            return VectorN<Real, 3>{u, w, u*w};
        }
    };

    ///////////////////////////////////////////////////////////////////////////
    // FIRST-ORDER DERIVATIVES: NDer1_u, NDer1_w
    ///////////////////////////////////////////////////////////////////////////
    
    TEST_CASE("SimplePlane - First Partial Derivatives", "[parametric_surface][nder1]")
    {
        SimplePlane plane;
        Real u = REAL(1.0), w = REAL(2.0);
        
        // Analytical: r_u = (1, 0, 1)
        VectorN<Real, 3> r_u = NDer1_u(plane, u, w);
        REQUIRE_THAT(r_u[0], WithinAbs(REAL(1.0), REAL(1e-9)));
        REQUIRE_THAT(r_u[1], WithinAbs(REAL(0.0), REAL(1e-9)));
        REQUIRE_THAT(r_u[2], WithinAbs(REAL(1.0), REAL(1e-9)));
        
        // Analytical: r_w = (0, 1, 2)
        VectorN<Real, 3> r_w = NDer1_w(plane, u, w);
        REQUIRE_THAT(r_w[0], WithinAbs(REAL(0.0), REAL(1e-9)));
        REQUIRE_THAT(r_w[1], WithinAbs(REAL(1.0), REAL(1e-9)));
        REQUIRE_THAT(r_w[2], WithinAbs(REAL(2.0), REAL(1e-9)));
    }
    
    TEST_CASE("Sphere - First Partial Derivatives at Multiple Points", "[parametric_surface][nder1]")
    {
        Real R = REAL(2.0);
        Sphere sphere(R);
        
        // Test at u=π/4, w=π/3
        Real u = Constants::PI / REAL(4.0);
        Real w = Constants::PI / REAL(3.0);
        
        VectorN<Real, 3> r_u_num = NDer1_u(sphere, u, w);
        VectorN<Real, 3> r_w_num = NDer1_w(sphere, u, w);
        
        // Analytical r_u = (R*cos(u)*cos(w), R*cos(u)*sin(w), -R*sin(u))
        Real r_u_x = R * cos(u) * cos(w);  // 2 * cos(π/4) * cos(π/3) = 2 * REAL(0.707) * REAL(0.5) ≈ REAL(0.707)
        Real r_u_y = R * cos(u) * sin(w);  // 2 * cos(π/4) * sin(π/3) ≈ REAL(1.225)
        Real r_u_z = -R * sin(u);          // -2 * sin(π/4) ≈ -REAL(1.414)
        
        REQUIRE_THAT(r_u_num[0], WithinAbs(r_u_x, REAL(1e-6)));
        REQUIRE_THAT(r_u_num[1], WithinAbs(r_u_y, REAL(1e-6)));
        REQUIRE_THAT(r_u_num[2], WithinAbs(r_u_z, REAL(1e-6)));
        
        // Analytical r_w = (-R*sin(u)*sin(w), R*sin(u)*cos(w), 0)
        Real r_w_x = -R * sin(u) * sin(w);  // -2 * sin(π/4) * sin(π/3) ≈ -REAL(1.225)
        Real r_w_y = R * sin(u) * cos(w);   // 2 * sin(π/4) * cos(π/3) ≈ REAL(0.707)
        Real r_w_z = REAL(0.0);
        
        REQUIRE_THAT(r_w_num[0], WithinAbs(r_w_x, REAL(1e-6)));
        REQUIRE_THAT(r_w_num[1], WithinAbs(r_w_y, REAL(1e-6)));
        REQUIRE_THAT(r_w_num[2], WithinAbs(r_w_z, REAL(1e-9)));
    }
    
    TEST_CASE("Paraboloid - First Partial Derivatives", "[parametric_surface][nder1]")
    {
        Paraboloid paraboloid;
        Real u = REAL(1.5), w = REAL(2.5);
        
        // Analytical: r_u = (1, 0, 2u)
        VectorN<Real, 3> r_u = NDer1_u(paraboloid, u, w);
        REQUIRE_THAT(r_u[0], WithinAbs(REAL(1.0), REAL(1e-9)));
        REQUIRE_THAT(r_u[1], WithinAbs(REAL(0.0), REAL(1e-9)));
        REQUIRE_THAT(r_u[2], WithinAbs(REAL(2.0) * u, 1e-6));
        
        // Analytical: r_w = (0, 1, 2w)
        VectorN<Real, 3> r_w = NDer1_w(paraboloid, u, w);
        REQUIRE_THAT(r_w[0], WithinAbs(REAL(0.0), REAL(1e-6)));  // First derivative numerical precision
        REQUIRE_THAT(r_w[1], WithinAbs(REAL(1.0), REAL(1e-9)));
        REQUIRE_THAT(r_w[2], WithinAbs(REAL(2.0) * w, 1e-6));
    }

    ///////////////////////////////////////////////////////////////////////////
    // CENTRAL DIFFERENCE FIRST DERIVATIVES: NDer2_u, NDer2_w
    // NOTE: NDer2_u means "2-point central difference", NOT second derivative!
    ///////////////////////////////////////////////////////////////////////////
    
    TEST_CASE("SimplePlane - Central Difference First Derivatives", "[parametric_surface][nder2]")
    {
        SimplePlane plane;
        Real u = REAL(1.0), w = REAL(2.0);
        
        // NDer2_u computes (f(u+h) - f(u-h))/(2h) = first derivative using 2-point stencil
        VectorN<Real, 3> r_u_2 = NDer2_u(plane, u, w);
        VectorN<Real, 3> r_w_2 = NDer2_w(plane, u, w);
        
        // For plane z = x + 2y: r_u = (1, 0, 1), r_w = (0, 1, 2)
        REQUIRE_THAT(r_u_2[0], WithinAbs(REAL(1.0), REAL(1e-9)));
        REQUIRE_THAT(r_u_2[1], WithinAbs(REAL(0.0), REAL(1e-9)));
        REQUIRE_THAT(r_u_2[2], WithinAbs(REAL(1.0), REAL(1e-9)));
        
        REQUIRE_THAT(r_w_2[0], WithinAbs(REAL(0.0), REAL(1e-9)));
        REQUIRE_THAT(r_w_2[1], WithinAbs(REAL(1.0), REAL(1e-9)));
        REQUIRE_THAT(r_w_2[2], WithinAbs(REAL(2.0), REAL(1e-9)));
    }
    
    TEST_CASE("Sphere - Central Difference First Derivatives", "[parametric_surface][nder2]")
    {
        Real R = REAL(2.0);
        Sphere sphere(R);
        Real u = Constants::PI / REAL(4.0);
        Real w = Constants::PI / REAL(3.0);
        
        VectorN<Real, 3> r_u_2 = NDer2_u(sphere, u, w);
        VectorN<Real, 3> r_w_2 = NDer2_w(sphere, u, w);
        
        // NDer2 computes central difference: (f(x+h)-f(x-h))/(2h)
        // This approximates first derivative, not second!
        // Analytical first derivative r_u at this point:
        Real r_u_x = R * cos(u) * cos(w);
        Real r_u_y = R * cos(u) * sin(w);
        Real r_u_z = -R * sin(u);
        
        REQUIRE_THAT(r_u_2[0], WithinAbs(r_u_x, REAL(1e-6)));
        REQUIRE_THAT(r_u_2[1], WithinAbs(r_u_y, REAL(1e-6)));
        REQUIRE_THAT(r_u_2[2], WithinAbs(r_u_z, REAL(1e-6)));
    }

    ///////////////////////////////////////////////////////////////////////////
    // PURE SECOND DERIVATIVES: NDer2_uu, NDer2_ww  
    ///////////////////////////////////////////////////////////////////////////
    
    TEST_CASE("Paraboloid - Pure Second Derivatives (d2r/du2)", "[parametric_surface][nder2_uu]")
    {
        Paraboloid paraboloid;
        Real u = REAL(1.5), w = REAL(2.5);
        
        // Analytical: r_uu = (0, 0, 2)
        VectorN<Real, 3> r_uu = NDer2_uu(paraboloid, u, w);
        REQUIRE_THAT(r_uu[0], WithinAbs(REAL(0.0), REAL(1e-5)));  // Second derivative numerical precision
        REQUIRE_THAT(r_uu[1], WithinAbs(REAL(0.0), REAL(1e-5)));
        REQUIRE_THAT(r_uu[2], WithinAbs(REAL(2.0), REAL(2e-5)));
        
        // Analytical: r_ww = (0, 0, 2)
        VectorN<Real, 3> r_ww = NDer2_ww(paraboloid, u, w);
        REQUIRE_THAT(r_ww[0], WithinAbs(REAL(0.0), REAL(1e-5)));  // Second derivative numerical precision
        REQUIRE_THAT(r_ww[1], WithinAbs(REAL(0.0), REAL(1e-5)));
        REQUIRE_THAT(r_ww[2], WithinAbs(REAL(2.0), REAL(2e-5)));
    }
    
    TEST_CASE("HyperbolicParaboloid - Pure Second Derivatives", "[parametric_surface][nder2_uu]")
    {
        HyperbolicParaboloid saddle;
        Real u = REAL(1.0), w = REAL(2.0);
        
        // Analytical: r_uu = (0, 0, 2)
        VectorN<Real, 3> r_uu = NDer2_uu(saddle, u, w);
        REQUIRE_THAT(r_uu[0], WithinAbs(REAL(0.0), REAL(1e-5)));  // Second derivative numerical precision
        REQUIRE_THAT(r_uu[1], WithinAbs(REAL(0.0), REAL(1e-5)));
        REQUIRE_THAT(r_uu[2], WithinAbs(REAL(2.0), REAL(2e-5)));
        
        // Analytical: r_ww = (0, 0, -2)
        VectorN<Real, 3> r_ww = NDer2_ww(saddle, u, w);
        REQUIRE_THAT(r_ww[0], WithinAbs(REAL(0.0), REAL(1e-5)));
        REQUIRE_THAT(r_ww[1], WithinAbs(REAL(0.0), REAL(1e-5)));
        REQUIRE_THAT(r_ww[2], WithinAbs(-REAL(2.0), REAL(2e-5)));
    }
    
    TEST_CASE("Sphere - Pure Second Derivatives", "[parametric_surface][nder2_uu]")
    {
        Real R = REAL(2.0);
        Sphere sphere(R);
        Real u = Constants::PI / REAL(4.0);
        Real w = Constants::PI / REAL(3.0);
        
        // Analytical: r_uu = (-R*sin(u)*cos(w), -R*sin(u)*sin(w), -R*cos(u))
        VectorN<Real, 3> r_uu = NDer2_uu(sphere, u, w);
        Real r_uu_x = -R * sin(u) * cos(w);
        Real r_uu_y = -R * sin(u) * sin(w);
        Real r_uu_z = -R * cos(u);
        
        REQUIRE_THAT(r_uu[0], WithinAbs(r_uu_x, REAL(1e-5)));  // Sphere with trig functions, second derivative
        REQUIRE_THAT(r_uu[1], WithinAbs(r_uu_y, REAL(1e-5)));
        REQUIRE_THAT(r_uu[2], WithinAbs(r_uu_z, REAL(1e-5)));
        
        // Analytical: r_ww = (-R*sin(u)*cos(w), -R*sin(u)*sin(w), 0)
        VectorN<Real, 3> r_ww = NDer2_ww(sphere, u, w);
        Real r_ww_x = -R * sin(u) * cos(w);
        Real r_ww_y = -R * sin(u) * sin(w);
        Real r_ww_z = REAL(0.0);
        
        REQUIRE_THAT(r_ww[0], WithinAbs(r_ww_x, REAL(1e-5)));
        REQUIRE_THAT(r_ww[1], WithinAbs(r_ww_y, REAL(1e-5)));
        REQUIRE_THAT(r_ww[2], WithinAbs(r_ww_z, REAL(1e-9)));
    }

    ///////////////////////////////////////////////////////////////////////////
    // MIXED SECOND DERIVATIVE: NDer2_uw (THE CRITICAL BUG LOCATION!)
    ///////////////////////////////////////////////////////////////////////////
    
    TEST_CASE("SimplePlane - Mixed Derivative (Should Be Zero)", "[parametric_surface][nder2_uw][critical]")
    {
        SimplePlane plane;
        Real u = REAL(1.0), w = REAL(2.0);
        
        // For plane z = x + 2y, all second derivatives are zero
        // Analytical: r_uw = (0, 0, 0)
        VectorN<Real, 3> r_uw = NDer2_uw(plane, u, w);
        
        REQUIRE_THAT(r_uw[0], WithinAbs(REAL(0.0), REAL(1e-9)));
        REQUIRE_THAT(r_uw[1], WithinAbs(REAL(0.0), REAL(1e-9)));
        REQUIRE_THAT(r_uw[2], WithinAbs(REAL(0.0), REAL(1e-9)));
    }
    
    TEST_CASE("MixedProductSurface - Mixed Derivative (Critical Test!)", "[parametric_surface][nder2_uw][critical]")
    {
        MixedProductSurface surface;
        
        // Test at multiple points - this surface has r_uw = (0, 0, 1) everywhere!
        Real test_points[][2] = {
            {REAL(0.0), REAL(0.0)}, {REAL(1.0), REAL(1.0)}, {REAL(2.0), REAL(3.0)}, {-REAL(1.5), REAL(2.5)}, {REAL(3.0), -REAL(2.0)}
        };
        
        for (auto& point : test_points) {
            Real u = point[0];
            Real w = point[1];
            
            // Analytical: r(u,w) = (u, w, uw)
            // r_u = (1, 0, w), r_w = (0, 1, u)
            // r_uw = ∂²r/∂u∂w = (0, 0, 1) ← THIS IS THE CRITICAL TEST!
            VectorN<Real, 3> r_uw = NDer2_uw(surface, u, w);
            
            INFO("Testing at u=" << u << ", w=" << w);
            REQUIRE_THAT(r_uw[0], WithinAbs(REAL(0.0), REAL(1e-6)));  // Mixed partial numerical precision
            REQUIRE_THAT(r_uw[1], WithinAbs(REAL(0.0), REAL(1e-6)));
            REQUIRE_THAT(r_uw[2], WithinAbs(REAL(1.0), REAL(1e-5)));  // THE CRITICAL COMPONENT!
        }
    }
    
    TEST_CASE("Paraboloid - Mixed Derivative (Should Be Zero)", "[parametric_surface][nder2_uw][critical]")
    {
        Paraboloid paraboloid;
        Real u = REAL(1.5), w = REAL(2.5);
        
        // Analytical: r(u,w) = (u, w, u² + w²)
        // r_uw = (0, 0, 0)
        VectorN<Real, 3> r_uw = NDer2_uw(paraboloid, u, w);
        
        REQUIRE_THAT(r_uw[0], WithinAbs(REAL(0.0), REAL(1e-9)));
        REQUIRE_THAT(r_uw[1], WithinAbs(REAL(0.0), REAL(1e-9)));
        REQUIRE_THAT(r_uw[2], WithinAbs(REAL(0.0), REAL(1e-9)));
    }
    
    TEST_CASE("HyperbolicParaboloid - Mixed Derivative (Should Be Zero)", "[parametric_surface][nder2_uw][critical]")
    {
        HyperbolicParaboloid saddle;
        Real u = REAL(1.0), w = REAL(2.0);
        
        // Analytical: r(u,w) = (u, w, u² - w²)
        // r_uw = (0, 0, 0)
        VectorN<Real, 3> r_uw = NDer2_uw(saddle, u, w);
        
        REQUIRE_THAT(r_uw[0], WithinAbs(REAL(0.0), REAL(1e-5)));  // Mixed partial numerical precision
        REQUIRE_THAT(r_uw[1], WithinAbs(REAL(0.0), REAL(1e-5)));
        REQUIRE_THAT(r_uw[2], WithinAbs(REAL(0.0), REAL(1e-5)));
    }
    
    TEST_CASE("Sphere - Mixed Derivative", "[parametric_surface][nder2_uw][critical]")
    {
        Real R = REAL(2.0);
        Sphere sphere(R);
        Real u = Constants::PI / REAL(4.0);
        Real w = Constants::PI / REAL(3.0);
        
        // Analytical: r_uw = (-R*cos(u)*sin(w), R*cos(u)*cos(w), 0)
        VectorN<Real, 3> r_uw = NDer2_uw(sphere, u, w);
        
        Real r_uw_x = -R * cos(u) * sin(w);
        Real r_uw_y = R * cos(u) * cos(w);
        Real r_uw_z = REAL(0.0);
        
        REQUIRE_THAT(r_uw[0], WithinAbs(r_uw_x, REAL(1e-6)));
        REQUIRE_THAT(r_uw[1], WithinAbs(r_uw_y, REAL(1e-6)));
        REQUIRE_THAT(r_uw[2], WithinAbs(r_uw_z, REAL(1e-9)));
    }

    ///////////////////////////////////////////////////////////////////////////
    // REGRESSION TESTS: Verify Bug Fix
    ///////////////////////////////////////////////////////////////////////////
    
    TEST_CASE("Regression: NDer2_uw formula is correct (not the old buggy one)", "[parametric_surface][nder2_uw][regression]")
    {
        // This test specifically verifies the bug fix
        // Old formula: (f(u+h,w+h) - 2f(u,w) + f(u-h,w-h))/h² - WRONG!
        // New formula: (f(u+h,w+h) - f(u+h,w-h) - f(u-h,w+h) + f(u-h,w-h))/(4h²) - CORRECT!
        
        MixedProductSurface surface;
        Real u = REAL(2.0), w = REAL(3.0);
        
        // With the old buggy formula, this would give approximately:
        // (f(REAL(2.001),REAL(3.001)) - 2f(2,3) + f(REAL(1.999),REAL(2.999)))/h²
        // = ((REAL(2.001), REAL(3.001), REAL(6.005001)) - 2(2,3,6) + (REAL(1.999), REAL(2.999), REAL(5.994001)))/h²
        // ≈ (0, 0, -REAL(0.998))/h² which would be huge and wrong!
        
        // With the correct formula, we get:
        VectorN<Real, 3> r_uw = NDer2_uw(surface, u, w);
        
        // The z-component should be exactly REAL(1.0), not some huge negative number
        REQUIRE_THAT(r_uw[2], WithinAbs(REAL(1.0), REAL(1e-5)));  // Mixed partial numerical precision
        REQUIRE(r_uw[2] > REAL(0.0));  // Must be positive!
        REQUIRE(std::abs(r_uw[2]) < REAL(2.0));  // Must be reasonable magnitude
    }
    
    TEST_CASE("Integration Test: Surfaces Second Fundamental Form Uses Correct NDer2_uw", "[parametric_surface][integration]")
    {
        // This verifies that the fix propagates to actual surface geometry calculations
        Cylinder cylinder(REAL(1.0), REAL(2.0));
        Real u = Constants::PI / REAL(4.0);
        Real w = REAL(1.0);
        
        Real L, M, N;
        cylinder.GetSecondNormalFormCoefficients(u, w, L, M, N);
        
        // For cylinder: M should be exactly 0 (no twist in surface)
        // Before the fix, M was -REAL(1.5) (completely wrong!)
        // After the fix, M should be REAL(0.0)
        REQUIRE_THAT(M, WithinAbs(REAL(0.0), REAL(1e-9)));
    }
}

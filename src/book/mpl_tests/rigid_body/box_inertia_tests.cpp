///////////////////////////////////////////////////////////////////////////////////////////
// Box Inertia Tensor Tests - Rigid Body Dynamics
///////////////////////////////////////////////////////////////////////////////////////////
//
// Verifies that RigidBodyBox computes correct inertia tensors.
// Formula: I = m/3 * diag(b²+c², a²+c², a²+b²) for box with half-extents (a, b, c)
//
// Tests cover:
// - Unit cube (symmetric case)
// - Arbitrary rectangular parallelepiped
// - Inertia tensor inverse correctness
// - World-frame inertia tensor transformation
//
///////////////////////////////////////////////////////////////////////////////////////////

#include <catch2/catch_all.hpp>
#include "TestPrecision.h"
#include "TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#include "base/Vector/VectorTypes.h"
#endif

#include "RigidBody/RigidBodies.h"

using namespace MML;
using namespace MPL;
using namespace MML::Testing;

using Catch::Matchers::WithinAbs;

namespace MPL::Tests::RigidBody::BoxInertiaTests
{
    //===================================================================================
    // S4.1: Box Inertia Tensor Tests
    //===================================================================================

    TEST_CASE("Box::Inertia_tensor_unit_cube", "[rigid_body][box][inertia]")
    {
        TEST_PRECISION_INFO();

        // Unit cube: mass=1, half-extents = (0.5, 0.5, 0.5)
        // Full dimensions: 1m x 1m x 1m
        // I = m/3 * diag(b²+c², a²+c², a²+b²)
        // With a=b=c=0.5: I = 1/3 * (0.25 + 0.25) = 1/6 ≈ 0.1667 for all diagonal
        RigidBodyBox cube(1.0, 0.5, 0.5, 0.5);
        
        Real expectedI = (1.0 / 3.0) * (0.25 + 0.25);  // = 1/6
        
        auto I = cube.InertiaTensorBody();
        
        // All diagonal elements should be equal for cube
        REQUIRE_THAT(I(0, 0), WithinAbs(expectedI, 1e-12));
        REQUIRE_THAT(I(1, 1), WithinAbs(expectedI, 1e-12));
        REQUIRE_THAT(I(2, 2), WithinAbs(expectedI, 1e-12));
        
        // Off-diagonal should be zero (principal axes aligned)
        REQUIRE_THAT(I(0, 1), WithinAbs(0.0, 1e-12));
        REQUIRE_THAT(I(0, 2), WithinAbs(0.0, 1e-12));
        REQUIRE_THAT(I(1, 2), WithinAbs(0.0, 1e-12));
    }

    TEST_CASE("Box::Inertia_tensor_arbitrary_box", "[rigid_body][box][inertia]")
    {
        TEST_PRECISION_INFO();

        // Arbitrary box: mass=5, half-extents = (1.0, 0.5, 0.3)
        Real m = 5.0;
        Real a = 1.0, b = 0.5, c = 0.3;
        
        RigidBodyBox box(m, a, b, c);
        
        // I_xx = m/3 * (b² + c²)
        // I_yy = m/3 * (a² + c²)
        // I_zz = m/3 * (a² + b²)
        Real coeff = m / 3.0;
        Real I_xx = coeff * (b*b + c*c);
        Real I_yy = coeff * (a*a + c*c);
        Real I_zz = coeff * (a*a + b*b);
        
        auto I = box.InertiaTensorBody();
        
        REQUIRE_THAT(I(0, 0), WithinAbs(I_xx, 1e-12));
        REQUIRE_THAT(I(1, 1), WithinAbs(I_yy, 1e-12));
        REQUIRE_THAT(I(2, 2), WithinAbs(I_zz, 1e-12));
        
        // Verify specific values
        // I_xx = 5/3 * (0.25 + 0.09) = 5/3 * 0.34 ≈ 0.5667
        // I_yy = 5/3 * (1.00 + 0.09) = 5/3 * 1.09 ≈ 1.8167
        // I_zz = 5/3 * (1.00 + 0.25) = 5/3 * 1.25 ≈ 2.0833
        REQUIRE_THAT(I(0, 0), WithinAbs(5.0/3.0 * 0.34, 1e-10));
        REQUIRE_THAT(I(1, 1), WithinAbs(5.0/3.0 * 1.09, 1e-10));
        REQUIRE_THAT(I(2, 2), WithinAbs(5.0/3.0 * 1.25, 1e-10));
    }

    TEST_CASE("Box::Inertia_tensor_inverse", "[rigid_body][box][inertia]")
    {
        TEST_PRECISION_INFO();

        RigidBodyBox box(5.0, 0.8, 0.6, 0.4);
        
        // Verify that I * I_inv = Identity
        auto I = box.InertiaTensorBody();
        auto Iinv = box.InertiaTensorBodyInv();
        
        // For diagonal matrix, inverse is just reciprocal of diagonal
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                Real expected = (i == j) ? 1.0 : 0.0;
                Real product = 0.0;
                for (int k = 0; k < 3; ++k)
                {
                    product += I(i, k) * Iinv(k, j);
                }
                REQUIRE_THAT(product, WithinAbs(expected, 1e-10));
            }
        }
    }

    TEST_CASE("Box::World_inertia_tensor_rotated", "[rigid_body][box][inertia]")
    {
        TEST_PRECISION_INFO();

        // Create box with some rotation
        RigidBodyState state;
        state.position = Vec3Cart(0, 0, 0);
        state.velocity = Vec3Cart(0, 0, 0);
        // 45° rotation around Z axis
        Real angle = Constants::PI / 4.0;
        state.orientation = Quaternion(std::cos(angle/2), 0, 0, std::sin(angle/2));
        state.angularVel = Vec3Cart(0, 0, 0);
        
        RigidBodyBox box(2.0, 1.0, 0.5, 0.3, state);
        
        // World inertia tensor should be symmetric
        auto Iworld = box.GetWorldInertiaTensor();
        
        REQUIRE_THAT(Iworld(0, 1), WithinAbs(Iworld(1, 0), 1e-12));
        REQUIRE_THAT(Iworld(0, 2), WithinAbs(Iworld(2, 0), 1e-12));
        REQUIRE_THAT(Iworld(1, 2), WithinAbs(Iworld(2, 1), 1e-12));
        
        // For 45° rotation around Z, the world inertia should mix x and y components
        // But Z component should be unchanged
        auto Ibody = box.InertiaTensorBody();
        REQUIRE_THAT(Iworld(2, 2), WithinAbs(Ibody(2, 2), 1e-12));
    }

    TEST_CASE("Box::Inertia_tensor_thin_rod", "[rigid_body][box][inertia]")
    {
        TEST_PRECISION_INFO();

        // Thin rod approximation: very long in X, tiny in Y and Z
        // This should approach the thin rod formula: I_yy = I_zz ≈ mL²/12
        Real m = 1.0;
        Real L = 2.0;  // Total length
        Real halfL = L / 2.0;
        Real epsilon = 0.001;  // Very small cross-section
        
        RigidBodyBox rod(m, halfL, epsilon, epsilon);
        
        auto I = rod.InertiaTensorBody();
        
        // For thin rod: I_yy = I_zz = mL²/12 = 1.0 * 4.0 / 12 = 1/3
        Real I_rod_expected = m * L * L / 12.0;
        
        // I_xx should be very small (just from tiny cross-section)
        REQUIRE(I(0, 0) < I_rod_expected * 0.001);
        
        // I_yy and I_zz should approximate thin rod formula
        // Our formula: I_yy = m/3 * (a² + c²) = m/3 * (halfL² + eps²) ≈ m*L²/12
        REQUIRE_THAT(I(1, 1), WithinAbs(I_rod_expected, 0.01));  // Allow small error
        REQUIRE_THAT(I(2, 2), WithinAbs(I_rod_expected, 0.01));
    }

    TEST_CASE("Box::Inertia_tensor_flat_plate", "[rigid_body][box][inertia]")
    {
        TEST_PRECISION_INFO();

        // Flat plate in XY plane: very thin in Z
        Real m = 1.0;
        Real a = 1.0, b = 0.5;  // Half-dimensions in X and Y
        Real epsilon = 0.001;   // Very thin in Z
        
        RigidBodyBox plate(m, a, b, epsilon);
        
        auto I = plate.InertiaTensorBody();
        
        // For thin plate: I_zz ≈ m/3 * (a² + b²) = m/12 * (4a² + 4b²) for full dimensions
        // Our formula uses half-extents directly
        Real I_zz_expected = m / 3.0 * (a*a + b*b);
        
        REQUIRE_THAT(I(2, 2), WithinAbs(I_zz_expected, 1e-10));
        
        // I_xx ≈ m/3 * b² (since c ≈ 0)
        REQUIRE_THAT(I(0, 0), WithinAbs(m / 3.0 * b*b, 0.001));
        
        // I_yy ≈ m/3 * a²
        REQUIRE_THAT(I(1, 1), WithinAbs(m / 3.0 * a*a, 0.001));
    }

} // namespace MPL::Tests::RigidBody::BoxInertiaTests

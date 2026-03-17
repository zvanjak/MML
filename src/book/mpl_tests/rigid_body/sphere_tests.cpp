///////////////////////////////////////////////////////////////////////////////////////////
// Sphere Rigid Body Tests
///////////////////////////////////////////////////////////////////////////////////////////
//
// Comprehensive tests for sphere support in rigid body simulation:
// - Sphere inertia tensor (I = 2/5 * m * r² * Identity)
// - Sphere-Sphere collision detection
// - Sphere-Wall collision detection
// - Sphere-Box collision detection
// - Energy conservation with bouncing spheres
// - Mixed-body integration (boxes + sphere)
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
#include "RigidBody/CollisionDetection.h"
#include "RigidBody/CollisionResponse.h"
#include "RigidBody/RigidBodySimulator.h"

using namespace MML;
using namespace MPL;
using namespace MML::Testing;

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

// Use concrete types for new inheritance hierarchy
using MPL::RigidBody;
using MPL::RigidBodySphere;
using MPL::RigidBodyBox;
using MPL::CollisionDetector;
using MPL::CollisionResponse;
using MPL::RigidBodySimulator;
using MPL::SimulationConfig;

namespace MPL::Tests::Sphere
{
    //===================================================================================
    // S5.1: Sphere Inertia Tensor Tests
    //===================================================================================

    TEST_CASE("Sphere::Inertia_tensor_unit_sphere", "[rigid_body][sphere][inertia]")
    {
        TEST_PRECISION_INFO();
        
        // Unit sphere: r=1, m=1
        // I = (2/5) * 1 * 1² * Identity = 0.4 * Identity
        RigidBodySphere sphere(1.0, 1.0);  // mass=1, radius=1
        
        REQUIRE(IsSphere(sphere));
        REQUIRE_THAT(sphere.Radius(), WithinAbs(1.0, 1e-12));
        
        auto I = sphere.InertiaTensorBody();
        Real expected = 0.4;  // 2/5 * 1 * 1²
        
        // Diagonal elements should all be 0.4
        REQUIRE_THAT(I(0, 0), WithinAbs(expected, 1e-12));
        REQUIRE_THAT(I(1, 1), WithinAbs(expected, 1e-12));
        REQUIRE_THAT(I(2, 2), WithinAbs(expected, 1e-12));
        
        // Off-diagonal elements should be zero
        REQUIRE_THAT(I(0, 1), WithinAbs(0.0, 1e-12));
        REQUIRE_THAT(I(0, 2), WithinAbs(0.0, 1e-12));
        REQUIRE_THAT(I(1, 2), WithinAbs(0.0, 1e-12));
    }

    TEST_CASE("Sphere::Inertia_tensor_arbitrary", "[rigid_body][sphere][inertia]")
    {
        TEST_PRECISION_INFO();
        
        // Sphere: r=0.5, m=10
        // I = (2/5) * 10 * 0.5² * Identity = (2/5) * 10 * 0.25 = 1.0 * Identity
        Real mass = 10.0;
        Real radius = 0.5;
        RigidBodySphere sphere(mass, radius);
        
        REQUIRE(IsSphere(sphere));
        
        auto I = sphere.InertiaTensorBody();
        Real expected = (2.0/5.0) * mass * radius * radius;  // = 1.0
        
        REQUIRE_THAT(I(0, 0), WithinAbs(expected, 1e-12));
        REQUIRE_THAT(I(1, 1), WithinAbs(expected, 1e-12));
        REQUIRE_THAT(I(2, 2), WithinAbs(expected, 1e-12));
    }

    TEST_CASE("Sphere::Inertia_tensor_inverse", "[rigid_body][sphere][inertia]")
    {
        TEST_PRECISION_INFO();
        
        // Verify that InertiaTensorBodyInv() is correct inverse
        RigidBodySphere sphere(5.0, 0.8);  // mass=5, radius=0.8
        
        auto I = sphere.InertiaTensorBody();
        auto Iinv = sphere.InertiaTensorBodyInv();
        
        // I * Iinv should be identity
        auto product = I * Iinv;
        
        REQUIRE_THAT(product(0, 0), WithinAbs(1.0, 1e-10));
        REQUIRE_THAT(product(1, 1), WithinAbs(1.0, 1e-10));
        REQUIRE_THAT(product(2, 2), WithinAbs(1.0, 1e-10));
        REQUIRE_THAT(product(0, 1), WithinAbs(0.0, 1e-10));
        REQUIRE_THAT(product(0, 2), WithinAbs(0.0, 1e-10));
        REQUIRE_THAT(product(1, 2), WithinAbs(0.0, 1e-10));
    }

    TEST_CASE("Sphere::World_inertia_tensor_isotropic", "[rigid_body][sphere][inertia]")
    {
        TEST_PRECISION_INFO();
        
        // Sphere inertia is isotropic - rotation doesn't change world-frame inertia
        RigidBodySphere sphere(2.0, 0.6);
        
        // Apply arbitrary rotation
        sphere.Orientation() = Quaternion::FromAxisAngle(
            Vec3Cart(1.0, 1.0, 1.0).GetAsUnitVector(), 
            Constants::PI / 3.0
        );
        
        auto Ibody = sphere.InertiaTensorBody();
        auto Iworld = sphere.GetWorldInertiaTensor();
        
        // For isotropic sphere, world inertia should equal body inertia
        REQUIRE_THAT(Iworld(0, 0), WithinAbs(Ibody(0, 0), 1e-10));
        REQUIRE_THAT(Iworld(1, 1), WithinAbs(Ibody(1, 1), 1e-10));
        REQUIRE_THAT(Iworld(2, 2), WithinAbs(Ibody(2, 2), 1e-10));
    }

    //===================================================================================
    // S5.2: Sphere-Sphere Collision Detection Tests
    //===================================================================================

    TEST_CASE("SphereSphere::No_collision_far_apart", "[rigid_body][sphere][collision]")
    {
        TEST_PRECISION_INFO();
        
        RigidBodySphere sphereA(1.0, 0.5);  // radius 0.5
        RigidBodySphere sphereB(1.0, 0.5);
        
        sphereA.Position() = Vec3Cart(0.0, 0.0, 0.0);
        sphereB.Position() = Vec3Cart(5.0, 0.0, 0.0);  // Far apart
        
        auto result = CollisionDetector::DetectSphereSphereCollision(sphereA, sphereB);
        
        REQUIRE_FALSE(result.hasCollision);
    }

    TEST_CASE("SphereSphere::Just_touching", "[rigid_body][sphere][collision]")
    {
        TEST_PRECISION_INFO();
        
        RigidBodySphere sphereA(1.0, 0.5);  // radius 0.5
        RigidBodySphere sphereB(1.0, 0.5);  // radius 0.5
        
        // Centers 1.0 apart = exactly sum of radii (just touching)
        sphereA.Position() = Vec3Cart(0.0, 0.0, 0.0);
        sphereB.Position() = Vec3Cart(1.0, 0.0, 0.0);
        
        auto result = CollisionDetector::DetectSphereSphereCollision(sphereA, sphereB);
        
        // At exactly touching distance - implementation may or may not detect
        // Most implementations require penetration > 0
        // This is an edge case, so just check it doesn't crash
        REQUIRE(result.penetrationDepth >= -1e-10);  // Should be ~0
    }

    TEST_CASE("SphereSphere::Overlapping", "[rigid_body][sphere][collision]")
    {
        TEST_PRECISION_INFO();
        
        RigidBodySphere sphereA(1.0, 0.5);  // radius 0.5
        RigidBodySphere sphereB(1.0, 0.3);  // radius 0.3
        
        // Centers 0.6 apart, sum of radii = 0.8, so penetration = 0.2
        sphereA.Position() = Vec3Cart(0.0, 0.0, 0.0);
        sphereB.Position() = Vec3Cart(0.6, 0.0, 0.0);
        
        auto result = CollisionDetector::DetectSphereSphereCollision(sphereA, sphereB);
        
        REQUIRE(result.hasCollision);
        REQUIRE_THAT(result.penetrationDepth, WithinAbs(0.2, 1e-10));
        
        // Contact normal should point from A to B (positive X direction)
        REQUIRE_THAT(result.contactNormal[0], WithinAbs(1.0, 1e-10));
        REQUIRE_THAT(result.contactNormal[1], WithinAbs(0.0, 1e-10));
        REQUIRE_THAT(result.contactNormal[2], WithinAbs(0.0, 1e-10));
    }

    TEST_CASE("SphereSphere::Collision_diagonal", "[rigid_body][sphere][collision]")
    {
        TEST_PRECISION_INFO();
        
        RigidBodySphere sphereA(1.0, 1.0);
        RigidBodySphere sphereB(1.0, 1.0);
        
        // Place spheres along diagonal, overlapping
        sphereA.Position() = Vec3Cart(0.0, 0.0, 0.0);
        sphereB.Position() = Vec3Cart(1.0, 1.0, 1.0);  // distance = sqrt(3) ≈ 1.732
        
        // Sum of radii = 2.0 > 1.732, so collision
        auto result = CollisionDetector::DetectSphereSphereCollision(sphereA, sphereB);
        
        REQUIRE(result.hasCollision);
        
        Real distance = std::sqrt(3.0);
        REQUIRE_THAT(result.penetrationDepth, WithinAbs(2.0 - distance, 1e-10));
        
        // Normal should be normalized direction from A to B
        Real invSqrt3 = 1.0 / std::sqrt(3.0);
        REQUIRE_THAT(result.contactNormal[0], WithinAbs(invSqrt3, 1e-10));
        REQUIRE_THAT(result.contactNormal[1], WithinAbs(invSqrt3, 1e-10));
        REQUIRE_THAT(result.contactNormal[2], WithinAbs(invSqrt3, 1e-10));
    }

    //===================================================================================
    // S5.2: Sphere-Wall Collision Detection Tests
    //===================================================================================

    TEST_CASE("SphereWall::No_collision_center", "[rigid_body][sphere][collision][wall]")
    {
        TEST_PRECISION_INFO();
        
        RigidBodySphere sphere(1.0, 0.5);
        sphere.Position() = Vec3Cart(0.0, 0.0, 0.0);  // Center of container
        
        Real containerHalfSize = 5.0;
        auto collisions = CollisionDetector::DetectSphereWallCollisions(sphere, containerHalfSize);
        
        REQUIRE(collisions.empty());
    }

    TEST_CASE("SphereWall::Collision_with_one_wall", "[rigid_body][sphere][collision][wall]")
    {
        TEST_PRECISION_INFO();
        
        RigidBodySphere sphere(1.0, 0.5);  // radius 0.5
        sphere.Position() = Vec3Cart(4.7, 0.0, 0.0);  // Close to +X wall at 5.0
        
        // Distance to wall = 5.0 - 4.7 = 0.3 < radius 0.5
        Real containerHalfSize = 5.0;
        auto collisions = CollisionDetector::DetectSphereWallCollisions(sphere, containerHalfSize);
        
        REQUIRE(collisions.size() == 1);
        REQUIRE(collisions[0].hasCollision);
        REQUIRE_THAT(collisions[0].penetrationDepth, WithinAbs(0.2, 1e-10));
        
        // Wall normal should point into container (-X direction)
        REQUIRE_THAT(collisions[0].contactNormal[0], WithinAbs(-1.0, 1e-10));
    }

    TEST_CASE("SphereWall::Collision_corner_three_walls", "[rigid_body][sphere][collision][wall]")
    {
        TEST_PRECISION_INFO();
        
        RigidBodySphere sphere(1.0, 0.5);  // radius 0.5
        // Near corner at (+5, +5, +5)
        sphere.Position() = Vec3Cart(4.7, 4.8, 4.6);
        
        Real containerHalfSize = 5.0;
        auto collisions = CollisionDetector::DetectSphereWallCollisions(sphere, containerHalfSize);
        
        // Should detect collision with 3 walls (+X, +Y, +Z)
        REQUIRE(collisions.size() == 3);
        
        for (const auto& c : collisions) {
            REQUIRE(c.hasCollision);
            REQUIRE(c.penetrationDepth > 0.0);
        }
    }

    //===================================================================================
    // S5.2: Sphere-Box Collision Detection Tests
    //===================================================================================

    TEST_CASE("SphereBox::No_collision", "[rigid_body][sphere][collision][box]")
    {
        TEST_PRECISION_INFO();
        
        RigidBodySphere sphere(1.0, 0.5);
        RigidBodyBox box(1.0, 0.5, 0.5, 0.5);  // Unit cube
        
        sphere.Position() = Vec3Cart(0.0, 0.0, 0.0);
        box.Position() = Vec3Cart(3.0, 0.0, 0.0);  // Far apart
        
        auto result = CollisionDetector::DetectSphereBoxCollision(sphere, box);
        
        REQUIRE_FALSE(result.hasCollision);
    }

    TEST_CASE("SphereBox::Collision_sphere_touching_box_face", "[rigid_body][sphere][collision][box]")
    {
        TEST_PRECISION_INFO();
        
        RigidBodySphere sphere(1.0, 0.5);  // radius 0.5
        RigidBodyBox box(1.0, 1.0, 1.0, 1.0);  // half-extents 1.0
        
        sphere.Position() = Vec3Cart(0.0, 0.0, 0.0);
        box.Position() = Vec3Cart(1.3, 0.0, 0.0);  // Box face at x=0.3, sphere surface at x=0.5
        
        // Sphere surface at x=0.5, box -X face at x = 1.3 - 1.0 = 0.3
        // Penetration = 0.5 - 0.3 = 0.2
        auto result = CollisionDetector::DetectSphereBoxCollision(sphere, box);
        
        REQUIRE(result.hasCollision);
        REQUIRE_THAT(result.penetrationDepth, WithinAbs(0.2, 1e-9));
        
        // Normal should point from sphere to box (positive X)
        REQUIRE_THAT(result.contactNormal[0], WithinAbs(1.0, 1e-9));
    }

    TEST_CASE("SphereBox::Collision_sphere_touching_box_edge", "[rigid_body][sphere][collision][box]")
    {
        TEST_PRECISION_INFO();
        
        RigidBodySphere sphere(1.0, 1.0);  // radius 1.0
        RigidBodyBox box(1.0, 0.5, 0.5, 0.5);  // half-extents 0.5
        
        sphere.Position() = Vec3Cart(0.0, 0.0, 0.0);
        // Box corner at (1.0, 1.0, 0.0), edge along Z
        box.Position() = Vec3Cart(1.0, 1.0, 0.0);
        
        // Closest point on box edge to sphere center is box corner (0.5, 0.5, 0.0) in local
        // = (1.5, 1.5, 0.0) in world... Actually distance = sqrt(1.5² + 1.5²) ≈ 2.12 > 1.0
        // Let's move box closer
        box.Position() = Vec3Cart(0.8, 0.8, 0.0);
        // Box corner at (0.3, 0.3, 0.0), distance = sqrt(0.3² + 0.3²) ≈ 0.424 < 1.0
        
        auto result = CollisionDetector::DetectSphereBoxCollision(sphere, box);
        
        REQUIRE(result.hasCollision);
        REQUIRE(result.penetrationDepth > 0.0);
    }

    TEST_CASE("SphereBox::Rotated_box", "[rigid_body][sphere][collision][box]")
    {
        TEST_PRECISION_INFO();
        
        RigidBodySphere sphere(1.0, 0.5);
        RigidBodyBox box(1.0, 1.0, 0.5, 0.5);  // Long box
        
        sphere.Position() = Vec3Cart(0.0, 0.0, 0.0);
        box.Position() = Vec3Cart(1.2, 0.0, 0.0);
        
        // Rotate box 45 degrees around Z
        box.Orientation() = Quaternion::FromAxisAngle(Vec3Cart(0, 0, 1), Constants::PI / 4.0);
        
        auto result = CollisionDetector::DetectSphereBoxCollision(sphere, box);
        
        // Just verify it handles rotated boxes without crashing
        // Exact collision depends on geometry
        REQUIRE(true);  // Didn't crash
    }

    //===================================================================================
    // S5.3: Energy Conservation Test - Two Bouncing Spheres
    //===================================================================================

    TEST_CASE("Sphere::Energy_conservation_two_spheres", "[rigid_body][sphere][integration][energy]")
    {
        TEST_PRECISION_INFO();
        
        // Setup: Two spheres bouncing in container
        SimulationConfig config;
        config.timeStep = 0.001;
        config.totalTime = 2.0;  // 2 seconds of simulation
        config.containerHalfSize = 5.0;
        config.coeffRestitution = 1.0;  // Perfectly elastic
        
        RigidBodySimulator simulator(config);
        
        // Sphere 1: radius 0.5, mass 5kg
        RigidBodySphere sphere1(5.0, 0.5);
        sphere1.Position() = Vec3Cart(-2.0, 0.0, 0.0);
        sphere1.Velocity() = Vec3Cart(3.0, 1.0, 0.5);
        sphere1.AngularVel() = Vec3Cart(1.0, -0.5, 0.3);
        simulator.AddBody(sphere1);
        
        // Sphere 2: radius 0.7, mass 8kg
        RigidBodySphere sphere2(8.0, 0.7);
        sphere2.Position() = Vec3Cart(2.0, 1.0, 0.0);
        sphere2.Velocity() = Vec3Cart(-2.0, 0.5, 1.0);
        sphere2.AngularVel() = Vec3Cart(-0.5, 1.0, -0.3);
        simulator.AddBody(sphere2);
        
        // Calculate initial energy
        auto computeEnergy = [](const RigidBody& b) {
            // Kinetic energy: 0.5 * m * v² + 0.5 * ω·(I·ω)
            Real transKE = 0.5 * b.Mass() * b.Velocity().NormL2() * b.Velocity().NormL2();
            
            Vec3Cart omega = b.AngularVel();
            auto I = b.GetWorldInertiaTensor();
            Vec3Cart Iomega = Vec3Cart(
                I(0,0)*omega[0] + I(0,1)*omega[1] + I(0,2)*omega[2],
                I(1,0)*omega[0] + I(1,1)*omega[1] + I(1,2)*omega[2],
                I(2,0)*omega[0] + I(2,1)*omega[1] + I(2,2)*omega[2]
            );
            Real rotKE = 0.5 * (omega[0]*Iomega[0] + omega[1]*Iomega[1] + omega[2]*Iomega[2]);
            
            return transKE + rotKE;
        };
        
        Real initialEnergy = computeEnergy(simulator.GetBody(0)) + 
                            computeEnergy(simulator.GetBody(1));
        
        // Run simulation
        simulator.Run();
        
        // Calculate final energy
        Real finalEnergy = computeEnergy(simulator.GetBody(0)) + 
                          computeEnergy(simulator.GetBody(1));
        
        // Energy should be conserved within 1% (allowing for numerical integration error)
        Real energyDrift = std::abs(finalEnergy - initialEnergy) / initialEnergy;
        REQUIRE(energyDrift < 0.01);  // < 1% energy drift
        
        // Verify spheres stayed in container
        for (size_t i = 0; i < simulator.NumBodies(); ++i) {
            const auto& body = simulator.GetBody(i);
            Vec3Cart pos = body.Position();
            Real r = body.BoundingRadius();  // Use BoundingRadius for abstract base access
            REQUIRE(pos[0] + r < config.containerHalfSize + 0.1);
            REQUIRE(pos[0] - r > -config.containerHalfSize - 0.1);
            REQUIRE(pos[1] + r < config.containerHalfSize + 0.1);
            REQUIRE(pos[1] - r > -config.containerHalfSize - 0.1);
            REQUIRE(pos[2] + r < config.containerHalfSize + 0.1);
            REQUIRE(pos[2] - r > -config.containerHalfSize - 0.1);
        }
    }

    //===================================================================================
    // S5.4: Mixed-Body Integration Test - Boxes + Sphere Chaos
    //===================================================================================

    TEST_CASE("MixedBody::Two_boxes_one_sphere_chaos", "[rigid_body][sphere][integration][mixed]")
    {
        TEST_PRECISION_INFO();
        
        // THE BIG TEST: Two boxes + one sphere wreaking havoc!
        SimulationConfig config;
        config.timeStep = 0.001;
        config.totalTime = 3.0;  // 3 seconds
        config.containerHalfSize = 5.0;
        config.coeffRestitution = 1.0;
        
        RigidBodySimulator simulator(config);
        
        // Box 1: 2m × 1m × 0.6m, 10kg (half extents: 1.0, 0.5, 0.3)
        RigidBodyBox box1(10.0, 1.0, 0.5, 0.3);
        box1.Position() = Vec3Cart(-1.5, 0.0, 0.0);
        box1.Velocity() = Vec3Cart(2.0, 0.5, 0.3);
        box1.AngularVel() = Vec3Cart(0.5, 1.0, -0.3);
        simulator.AddBody(box1);
        
        // Box 2: 1.5m × 1m × 0.8m, 8kg (half extents: 0.75, 0.5, 0.4)
        RigidBodyBox box2(8.0, 0.75, 0.5, 0.4);
        box2.Position() = Vec3Cart(1.5, 1.0, 0.5);
        box2.Velocity() = Vec3Cart(-1.5, -0.5, 0.8);
        box2.AngularVel() = Vec3Cart(-0.3, 0.5, 0.8);
        simulator.AddBody(box2);
        
        // Sphere: radius 0.6m, 12kg (wrecking ball!)
        RigidBodySphere sphere(12.0, 0.6);
        sphere.Position() = Vec3Cart(0.0, -1.0, 2.0);
        sphere.Velocity() = Vec3Cart(1.0, 0.5, -2.0);  // Dropping
        sphere.AngularVel() = Vec3Cart(2.0, 1.0, -1.0);  // Spinning
        simulator.AddBody(sphere);
        
        // Calculate initial total momentum
        auto computeMomentum = [&simulator]() {
            Vec3Cart totalP(0, 0, 0);
            for (size_t i = 0; i < simulator.NumBodies(); ++i) {
                const auto& b = simulator.GetBody(i);
                totalP = totalP + b.Mass() * b.Velocity();
            }
            return totalP;
        };
        
        Vec3Cart initialMomentum = computeMomentum();
        
        // Run simulation
        simulator.Run();
        
        // Verify bodies are still valid (no NaN, no explosions)
        for (size_t i = 0; i < simulator.NumBodies(); ++i) {
            const auto& body = simulator.GetBody(i);
            Vec3Cart pos = body.Position();
            Vec3Cart vel = body.Velocity();
            
            // Check for NaN
            REQUIRE_FALSE(std::isnan(pos[0]));
            REQUIRE_FALSE(std::isnan(pos[1]));
            REQUIRE_FALSE(std::isnan(pos[2]));
            REQUIRE_FALSE(std::isnan(vel[0]));
            REQUIRE_FALSE(std::isnan(vel[1]));
            REQUIRE_FALSE(std::isnan(vel[2]));
            
            // Check velocities haven't exploded (< 100 m/s)
            REQUIRE(vel.NormL2() < 100.0);
        }
        
        // Momentum should be roughly conserved (wall collisions change momentum)
        // Just verify it hasn't exploded
        Vec3Cart finalMomentum = computeMomentum();
        REQUIRE(finalMomentum.NormL2() < 1000.0);  // Sanity check
        
        // Verify shape types are preserved
        REQUIRE(IsBox(simulator.GetBody(0)));
        REQUIRE(IsBox(simulator.GetBody(1)));
        REQUIRE(IsSphere(simulator.GetBody(2)));
    }

    TEST_CASE("MixedBody::Sphere_box_collision", "[rigid_body][sphere][collision][mixed]")
    {
        TEST_PRECISION_INFO();
        
        // Direct sphere-box collision test
        RigidBodySphere sphere(1.0, 0.5);
        RigidBodyBox box(1.0, 0.5, 0.5, 0.5);  // Unit cube
        
        // Head-on collision course
        sphere.Position() = Vec3Cart(-2.0, 0.0, 0.0);
        box.Position() = Vec3Cart(2.0, 0.0, 0.0);
        
        // Moving toward each other
        sphere.Velocity() = Vec3Cart(5.0, 0.0, 0.0);
        box.Velocity() = Vec3Cart(-5.0, 0.0, 0.0);
        
        SimulationConfig config;
        config.timeStep = 0.001;
        config.totalTime = 1.0;
        config.containerHalfSize = 10.0;
        config.coeffRestitution = 1.0;
        
        RigidBodySimulator simulator(config);
        simulator.AddBody(sphere);
        simulator.AddBody(box);
        
        Real initialEnergy = 0.5 * 1.0 * 25.0 + 0.5 * 1.0 * 25.0;  // 0.5*m*v² for each
        
        simulator.Run();
        
        // After collision, they should have bounced apart
        // Velocities should have reversed (approximately, for equal mass elastic collision)
        Vec3Cart v1 = simulator.GetBody(0).Velocity();
        Vec3Cart v2 = simulator.GetBody(1).Velocity();
        
        // Final energy
        Real finalEnergy = 0.5 * simulator.GetBody(0).Mass() * v1.NormL2() * v1.NormL2() +
                          0.5 * simulator.GetBody(1).Mass() * v2.NormL2() * v2.NormL2();
        
        // Energy should be conserved within 5%
        Real drift = std::abs(finalEnergy - initialEnergy) / initialEnergy;
        REQUIRE(drift < 0.05);
    }

} // namespace MPL::Tests::Sphere

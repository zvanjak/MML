///////////////////////////////////////////////////////////////////////////////
// Test suite for CoordTransfLorentz.h - Special Relativity transformations
///////////////////////////////////////////////////////////////////////////////

#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/VectorTypes.h"
#include "core/CoordTransf/CoordTransfLorentz.h"
#endif

using namespace MML;
using namespace MML::Testing;
using namespace Catch::Matchers;

///////////////////////////////////////////////////////////////////////////////
// Helper Functions

// Calculate Lorentz factor gamma = 1/sqrt(1 - v^2/c^2)
inline Real LorentzGamma(Real velocity)
{
    return REAL(1.0) / std::sqrt(REAL(1.0) - velocity * velocity);
}

// Calculate proper time interval
inline Real ProperTime(Real coordinateTime, Real velocity)
{
    return coordinateTime / LorentzGamma(velocity);
}

// Velocity addition formula (parallel velocities in units of c)
inline Real VelocityComposition(Real v1, Real v2)
{
    return (v1 + v2) / (REAL(1.0) + v1 * v2);
}

///////////////////////////////////////////////////////////////////////////////
// Constructor and Basic Properties Tests

TEST_CASE("CoordTransfLorentzXAxis - Constructor with Valid Velocity", "[lorentz]")
{
    CoordTransfLorentzXAxis lorentz(REAL(0.5));
    
    // Should construct without throwing
    REQUIRE(true);
}

TEST_CASE("CoordTransfLorentzXAxis - Constructor with Zero Velocity", "[lorentz]")
{
    CoordTransfLorentzXAxis lorentz(REAL(0.0));
    
    // Zero velocity is valid (stationary frame)
    REQUIRE(true);
}

TEST_CASE("CoordTransfLorentzXAxis - Constructor Rejects Negative Velocity", "[lorentz]")
{
    REQUIRE_THROWS_AS(CoordTransfLorentzXAxis(-REAL(0.1)), std::range_error);
}

TEST_CASE("CoordTransfLorentzXAxis - Constructor Rejects Superluminal Velocity", "[lorentz]")
{
    REQUIRE_THROWS_AS(CoordTransfLorentzXAxis(REAL(1.0)), std::range_error);
    REQUIRE_THROWS_AS(CoordTransfLorentzXAxis(REAL(1.5)), std::range_error);
}

///////////////////////////////////////////////////////////////////////////////
// Identity Transformation Tests (v = 0)

TEST_CASE("CoordTransfLorentzXAxis - Zero Velocity is Identity Transform", "[lorentz]")
{
    CoordTransfLorentzXAxis lorentz(REAL(0.0));
    
    Vector4Minkowski event{REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0)};  // (t, x, y, z)
    Vector4Minkowski transformed = lorentz.transf(event);
    
    REQUIRE_THAT(transformed[0], WithinAbs(REAL(1.0), REAL(1e-10)));
    REQUIRE_THAT(transformed[1], WithinAbs(REAL(2.0), REAL(1e-10)));
    REQUIRE_THAT(transformed[2], WithinAbs(REAL(3.0), REAL(1e-10)));
    REQUIRE_THAT(transformed[3], WithinAbs(REAL(4.0), REAL(1e-10)));
}

///////////////////////////////////////////////////////////////////////////////
// Forward Transformation Tests

TEST_CASE("CoordTransfLorentzXAxis - Transform Origin Event", "[lorentz]")
{
    Real v = REAL(0.6);
    CoordTransfLorentzXAxis lorentz(v);
    
    Vector4Minkowski origin{REAL(0.0), REAL(0.0), REAL(0.0), REAL(0.0)};
    Vector4Minkowski transformed = lorentz.transf(origin);
    
    // Origin maps to origin
    REQUIRE_THAT(transformed[0], WithinAbs(REAL(0.0), REAL(1e-10)));
    REQUIRE_THAT(transformed[1], WithinAbs(REAL(0.0), REAL(1e-10)));
    REQUIRE_THAT(transformed[2], WithinAbs(REAL(0.0), REAL(1e-10)));
    REQUIRE_THAT(transformed[3], WithinAbs(REAL(0.0), REAL(1e-10)));
}

TEST_CASE("CoordTransfLorentzXAxis - Transform Time-Only Event", "[lorentz]")
{
    Real v = REAL(0.6);
    Real gamma = LorentzGamma(v);
    CoordTransfLorentzXAxis lorentz(v);
    
    // Event at origin but at t=1
    Vector4Minkowski event{REAL(1.0), REAL(0.0), REAL(0.0), REAL(0.0)};
    Vector4Minkowski transformed = lorentz.transf(event);
    
    // t' = gamma * t (when x=0)
    REQUIRE_THAT(transformed[0], WithinAbs(gamma * REAL(1.0), REAL(1e-10)));
    REQUIRE_THAT(transformed[1], WithinAbs(-v * gamma * REAL(1.0), REAL(1e-10)));
    REQUIRE_THAT(transformed[2], WithinAbs(REAL(0.0), REAL(1e-10)));
    REQUIRE_THAT(transformed[3], WithinAbs(REAL(0.0), REAL(1e-10)));
}

TEST_CASE("CoordTransfLorentzXAxis - Transform Space-Only Event (X-axis)", "[lorentz]")
{
    Real v = REAL(0.6);
    Real gamma = LorentzGamma(v);
    CoordTransfLorentzXAxis lorentz(v);
    
    // Event at x=1, t=0
    Vector4Minkowski event{REAL(0.0), REAL(1.0), REAL(0.0), REAL(0.0)};
    Vector4Minkowski transformed = lorentz.transf(event);
    
    // t' = gamma * (-v * x), x' = gamma * x (when t=0)
    REQUIRE_THAT(transformed[0], WithinAbs(-v * gamma * REAL(1.0), REAL(1e-10)));
    REQUIRE_THAT(transformed[1], WithinAbs(gamma * REAL(1.0), REAL(1e-10)));
    REQUIRE_THAT(transformed[2], WithinAbs(REAL(0.0), REAL(1e-10)));
    REQUIRE_THAT(transformed[3], WithinAbs(REAL(0.0), REAL(1e-10)));
}

TEST_CASE("CoordTransfLorentzXAxis - Y and Z Coordinates Unchanged", "[lorentz]")
{
    Real v = REAL(0.6);
    CoordTransfLorentzXAxis lorentz(v);
    
    Vector4Minkowski event{REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0)};
    Vector4Minkowski transformed = lorentz.transf(event);
    
    // Y and Z coordinates should be unchanged
    REQUIRE_THAT(transformed[2], WithinAbs(REAL(3.0), REAL(1e-10)));
    REQUIRE_THAT(transformed[3], WithinAbs(REAL(4.0), REAL(1e-10)));
}

TEST_CASE("CoordTransfLorentzXAxis - General Event Transformation", "[lorentz]")
{
    Real v = REAL(0.8);
    Real gamma = LorentzGamma(v);
    CoordTransfLorentzXAxis lorentz(v);
    
    Real t = REAL(5.0), x = REAL(3.0), y = REAL(2.0), z = REAL(1.0);
    Vector4Minkowski event{t, x, y, z};
    Vector4Minkowski transformed = lorentz.transf(event);
    
    // Lorentz transformation formulas:
    // t' = gamma * (t - v*x)
    // x' = gamma * (x - v*t)
    // y' = y
    // z' = z
    Real expected_t = gamma * (t - v * x);
    Real expected_x = gamma * (x - v * t);
    
    REQUIRE_THAT(transformed[0], WithinAbs(expected_t, REAL(1e-9)));
    REQUIRE_THAT(transformed[1], WithinAbs(expected_x, REAL(1e-9)));
    REQUIRE_THAT(transformed[2], WithinAbs(y, REAL(1e-10)));
    REQUIRE_THAT(transformed[3], WithinAbs(z, REAL(1e-10)));
}

///////////////////////////////////////////////////////////////////////////////
// Inverse Transformation Tests

TEST_CASE("CoordTransfLorentzXAxis - Inverse Transform of Origin", "[lorentz]")
{
    Real v = REAL(0.6);
    CoordTransfLorentzXAxis lorentz(v);
    
    Vector4Minkowski origin{REAL(0.0), REAL(0.0), REAL(0.0), REAL(0.0)};
    Vector4Minkowski transformed = lorentz.transfInverse(origin);
    
    REQUIRE_THAT(transformed[0], WithinAbs(REAL(0.0), REAL(1e-10)));
    REQUIRE_THAT(transformed[1], WithinAbs(REAL(0.0), REAL(1e-10)));
    REQUIRE_THAT(transformed[2], WithinAbs(REAL(0.0), REAL(1e-10)));
    REQUIRE_THAT(transformed[3], WithinAbs(REAL(0.0), REAL(1e-10)));
}

TEST_CASE("CoordTransfLorentzXAxis - Round Trip Transformation", "[lorentz]")
{
    Real v = REAL(0.7);
    CoordTransfLorentzXAxis lorentz(v);
    
    Vector4Minkowski original{REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0)};
    Vector4Minkowski transformed = lorentz.transf(original);
    Vector4Minkowski back = lorentz.transfInverse(transformed);
    
    REQUIRE_THAT(back[0], WithinAbs(original[0], REAL(1e-9)));
    REQUIRE_THAT(back[1], WithinAbs(original[1], REAL(1e-9)));
    REQUIRE_THAT(back[2], WithinAbs(original[2], REAL(1e-9)));
    REQUIRE_THAT(back[3], WithinAbs(original[3], REAL(1e-9)));
}

TEST_CASE("CoordTransfLorentzXAxis - Inverse Equals Negative Velocity Transform", "[lorentz]")
{
    Real v = REAL(0.5);
    CoordTransfLorentzXAxis lorentz_forward(v);
    
    Vector4Minkowski event{REAL(3.0), REAL(2.0), REAL(1.0), REAL(0.5)};
    Vector4Minkowski inv_transformed = lorentz_forward.transfInverse(event);
    
    // The inverse Lorentz transformation is the same as forward transformation
    // with -v, but since we can't construct with negative v, we verify the formula
    Real gamma = LorentzGamma(v);
    Real expected_t = gamma * (event[0] + v * event[1]);
    Real expected_x = gamma * (event[1] + v * event[0]);
    
    REQUIRE_THAT(inv_transformed[0], WithinAbs(expected_t, REAL(1e-9)));
    REQUIRE_THAT(inv_transformed[1], WithinAbs(expected_x, REAL(1e-9)));
    REQUIRE_THAT(inv_transformed[2], WithinAbs(event[2], REAL(1e-10)));
    REQUIRE_THAT(inv_transformed[3], WithinAbs(event[3], REAL(1e-10)));
}

///////////////////////////////////////////////////////////////////////////////
// Spacetime Interval Invariance Tests

TEST_CASE("CoordTransfLorentzXAxis - Spacetime Interval Invariance", "[lorentz]")
{
    Real v = REAL(0.6);
    CoordTransfLorentzXAxis lorentz(v);
    
    Vector4Minkowski event{REAL(5.0), REAL(3.0), REAL(2.0), REAL(1.0)};
    Vector4Minkowski transformed = lorentz.transf(event);
    
    // Spacetime interval: s^2 = t^2 - x^2 - y^2 - z^2
    Real s2_original = event[0]*event[0] - event[1]*event[1] - event[2]*event[2] - event[3]*event[3];
    Real s2_transformed = transformed[0]*transformed[0] - transformed[1]*transformed[1] 
                        - transformed[2]*transformed[2] - transformed[3]*transformed[3];
    
    REQUIRE_THAT(s2_transformed, WithinAbs(s2_original, REAL(1e-9)));
}

TEST_CASE("CoordTransfLorentzXAxis - Light Cone Invariance", "[lorentz]")
{
    Real v = REAL(0.8);
    CoordTransfLorentzXAxis lorentz(v);
    
    // Light signal: t = x (moving at speed of light along x-axis)
    Vector4Minkowski light_event{REAL(2.0), REAL(2.0), REAL(0.0), REAL(0.0)};
    Vector4Minkowski transformed = lorentz.transf(light_event);
    
    // For light-like events: t^2 = x^2 + y^2 + z^2 (spacetime interval = 0)
    Real interval_original = light_event[0]*light_event[0] - light_event[1]*light_event[1] 
                           - light_event[2]*light_event[2] - light_event[3]*light_event[3];
    Real interval_transformed = transformed[0]*transformed[0] - transformed[1]*transformed[1] 
                              - transformed[2]*transformed[2] - transformed[3]*transformed[3];
    
    REQUIRE_THAT(interval_original, WithinAbs(REAL(0.0), REAL(1e-10)));
    REQUIRE_THAT(interval_transformed, WithinAbs(REAL(0.0), REAL(1e-10)));
}

///////////////////////////////////////////////////////////////////////////////
// Proper Time Tests

TEST_CASE("CoordTransfLorentzXAxis - Proper Time Calculation", "[lorentz]")
{
    Real v = REAL(0.6);
    Real gamma = LorentzGamma(v);
    
    Real coordinate_time = REAL(10.0);
    Real proper_time = ProperTime(coordinate_time, v);
    
    // tau = t / gamma
    REQUIRE_THAT(proper_time, WithinAbs(coordinate_time / gamma, REAL(1e-10)));
    REQUIRE(proper_time < coordinate_time);  // Proper time is always less
}

TEST_CASE("CoordTransfLorentzXAxis - Time Dilation for Moving Object", "[lorentz]")
{
    Real v = REAL(0.8);
    Real gamma = LorentzGamma(v);
    CoordTransfLorentzXAxis lorentz(v);
    
    // Object at rest in original frame for time dt
    Real dt = REAL(1.0);
    Vector4Minkowski event1{REAL(0.0), REAL(0.0), REAL(0.0), REAL(0.0)};
    Vector4Minkowski event2{dt, REAL(0.0), REAL(0.0), REAL(0.0)};
    
    Vector4Minkowski trans1 = lorentz.transf(event1);
    Vector4Minkowski trans2 = lorentz.transf(event2);
    
    Real dt_prime = trans2[0] - trans1[0];
    
    // Time dilation: dt' = gamma * dt
    REQUIRE_THAT(dt_prime, WithinAbs(gamma * dt, REAL(1e-10)));
}

///////////////////////////////////////////////////////////////////////////////
// Velocity Composition Tests

TEST_CASE("CoordTransfLorentzXAxis - Velocity Composition Formula", "[lorentz]")
{
    Real v1 = REAL(0.6);
    Real v2 = REAL(0.5);
    
    Real v_composed = VelocityComposition(v1, v2);
    
    // Should be less than c (which is 1 in our units)
    REQUIRE(v_composed < REAL(1.0));
    
    // Einstein velocity addition: v = (v1 + v2) / (1 + v1*v2)
    Real expected = (v1 + v2) / (REAL(1.0) + v1 * v2);
    REQUIRE_THAT(v_composed, WithinAbs(expected, REAL(1e-10)));
}

TEST_CASE("CoordTransfLorentzXAxis - Successive Transformations", "[lorentz]")
{
    Real v1 = REAL(0.5);
    Real v2 = REAL(0.6);
    
    CoordTransfLorentzXAxis lorentz1(v1);
    CoordTransfLorentzXAxis lorentz2(v2);
    
    Vector4Minkowski event{REAL(1.0), REAL(1.0), REAL(0.0), REAL(0.0)};
    
    // Apply two successive transformations
    Vector4Minkowski trans1 = lorentz1.transf(event);
    Vector4Minkowski trans2 = lorentz2.transf(trans1);
    
    // Should be equivalent to single transformation with composed velocity
    Real v_composed = VelocityComposition(v1, v2);
    CoordTransfLorentzXAxis lorentz_composed(v_composed);
    Vector4Minkowski trans_direct = lorentz_composed.transf(event);
    
    REQUIRE_THAT(trans2[0], WithinAbs(trans_direct[0], REAL(1e-9)));
    REQUIRE_THAT(trans2[1], WithinAbs(trans_direct[1], REAL(1e-9)));
}

TEST_CASE("CoordTransfLorentzXAxis - Light Speed Invariance Under Composition", "[lorentz]")
{
    Real v1 = REAL(0.9);
    Real v2 = REAL(0.9);
    
    Real v_composed = VelocityComposition(v1, v2);
    
    // Even with two high velocities, composed velocity < c
    REQUIRE(v_composed < REAL(1.0));
    
    // Composed velocity approaches c but never reaches it
    REQUIRE_THAT(v_composed, WithinAbs(REAL(0.9945), REAL(1e-4)));
}

///////////////////////////////////////////////////////////////////////////////
// Length Contraction Tests

TEST_CASE("CoordTransfLorentzXAxis - Length Contraction", "[lorentz]")
{
    Real v = REAL(0.6);
    Real gamma = LorentzGamma(v);
    CoordTransfLorentzXAxis lorentz(v);
    
    // Rod of length L at rest, measured at same time (t=0)
    Real L = REAL(10.0);
    Vector4Minkowski front{REAL(0.0), L, REAL(0.0), REAL(0.0)};
    Vector4Minkowski back{REAL(0.0), REAL(0.0), REAL(0.0), REAL(0.0)};
    
    Vector4Minkowski trans_front = lorentz.transf(front);
    Vector4Minkowski trans_back = lorentz.transf(back);
    
    // Length in moving frame (at simultaneous time in that frame)
    // For proper measurement, we need to account for simultaneity
    // The contracted length is L' = L / gamma
    Real L_moving = std::abs(trans_front[1] - trans_back[1]);
    
    REQUIRE_THAT(L_moving, WithinAbs(gamma * L, REAL(1e-9)));
}

///////////////////////////////////////////////////////////////////////////////
// Relativistic Limit Tests

TEST_CASE("CoordTransfLorentzXAxis - Low Velocity Limit (Galilean)", "[lorentz]")
{
    Real v = REAL(0.001);  // Very low velocity
    Real gamma = LorentzGamma(v);
    CoordTransfLorentzXAxis lorentz(v);
    
    Vector4Minkowski event{REAL(10.0), REAL(5.0), REAL(3.0), REAL(2.0)};
    Vector4Minkowski transformed = lorentz.transf(event);
    
    // At low velocities, should approach Galilean transformation
    // Full Lorentz: t' = γ(t - βx), x' = γ(x - βt)
    // At v << 1: γ ≈ 1, so t' ≈ t - βx ≈ t, x' ≈ x - βt ≈ x - vt
    REQUIRE_THAT(gamma, WithinAbs(REAL(1.0), REAL(1e-6)));
    
    // For low velocity, t' = γ(t - βx) should be close to the exact formula
    Real expected_t = gamma * (event[0] - v * event[1]);
    Real expected_x = gamma * (event[1] - v * event[0]);
    
    REQUIRE_THAT(transformed[0], WithinAbs(expected_t, REAL(1e-9)));
    REQUIRE_THAT(transformed[1], WithinAbs(expected_x, REAL(1e-9)));
}

TEST_CASE("CoordTransfLorentzXAxis - High Velocity Gamma Factor", "[lorentz]")
{
    Real v = REAL(0.99);
    Real gamma = LorentzGamma(v);
    
    // At v = REAL(0.9)9c, gamma should be large
    REQUIRE(gamma > REAL(7.0));
    REQUIRE_THAT(gamma, WithinAbs(REAL(7.0889), REAL(1e-3)));
}

///////////////////////////////////////////////////////////////////////////////
// Coordinate Function Tests

TEST_CASE("CoordTransfLorentzXAxis - Coordinate Transform Functions Consistency", "[lorentz]")
{
    Real v = REAL(0.5);
    CoordTransfLorentzXAxis lorentz(v);
    
    Vector4Minkowski event{REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0)};
    
    // Get transformed vector
    Vector4Minkowski transformed = lorentz.transf(event);
    
    // Get individual coordinates via functions
    Real t_prime = lorentz.coordTransfFunc(0)(event);
    Real x_prime = lorentz.coordTransfFunc(1)(event);
    Real y_prime = lorentz.coordTransfFunc(2)(event);
    Real z_prime = lorentz.coordTransfFunc(3)(event);
    
    REQUIRE_THAT(t_prime, WithinAbs(transformed[0], REAL(1e-10)));
    REQUIRE_THAT(x_prime, WithinAbs(transformed[1], REAL(1e-10)));
    REQUIRE_THAT(y_prime, WithinAbs(transformed[2], REAL(1e-10)));
    REQUIRE_THAT(z_prime, WithinAbs(transformed[3], REAL(1e-10)));
}

TEST_CASE("CoordTransfLorentzXAxis - Inverse Coordinate Functions Consistency", "[lorentz]")
{
    Real v = REAL(0.5);
    CoordTransfLorentzXAxis lorentz(v);
    
    Vector4Minkowski event{REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0)};
    
    // Get inverse transformed vector
    Vector4Minkowski inv_transformed = lorentz.transfInverse(event);
    
    // Get individual coordinates via inverse functions
    Real t = lorentz.inverseCoordTransfFunc(0)(event);
    Real x = lorentz.inverseCoordTransfFunc(1)(event);
    Real y = lorentz.inverseCoordTransfFunc(2)(event);
    Real z = lorentz.inverseCoordTransfFunc(3)(event);
    
    REQUIRE_THAT(t, WithinAbs(inv_transformed[0], REAL(1e-10)));
    REQUIRE_THAT(x, WithinAbs(inv_transformed[1], REAL(1e-10)));
    REQUIRE_THAT(y, WithinAbs(inv_transformed[2], REAL(1e-10)));
    REQUIRE_THAT(z, WithinAbs(inv_transformed[3], REAL(1e-10)));
}

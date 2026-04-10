/**
 * @file fields_tests.cpp
 * @brief Comprehensive tests for Fields.h
 * 
 * Tests cover:
 * - InverseRadialPotentialFieldCart (static functions)
 * - InverseRadialPotentialFieldSpher (static functions)
 * - InverseRadialPotentialFieldCyl (static functions)
 * - InverseRadialPotentialForceFieldCart (static functions)
 * - InverseRadialPotentialForceFieldSph (static functions)
 * - InverseRadialFieldCart class (IScalarFunction)
 * - InverseRadialFieldSpher class (IScalarFunction)
 * - InverseRadialForceFieldCart class (IVectorFunction)
 * - InverseRadialForceFieldSpher class (IVectorFunction)
 * 
 * These represent gravitational and Coulomb-type 1/r potentials and 1/r² forces.
 */

#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include <cmath>

#include "core/Fields.h"

using namespace MML;
using namespace MML::Fields;

namespace MML::Tests::Core::FieldsTests {

///////////////////////////////////////////////////////////////////////////////
//                   INVERSE RADIAL POTENTIAL - CARTESIAN                    //
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("InverseRadialPotentialFieldCart - Unit distance", "[Fields][Potential][Cartesian]") {
    // Point at distance 1 from origin along x-axis
    VectorN<Real, 3> x{1.0, 0.0, 0.0};
    Real result = InverseRadialPotentialFieldCart(x);
    REQUIRE(result == Catch::Approx(1.0));
}

TEST_CASE("InverseRadialPotentialFieldCart - Distance 2", "[Fields][Potential][Cartesian]") {
    VectorN<Real, 3> x{2.0, 0.0, 0.0};
    Real result = InverseRadialPotentialFieldCart(x);
    REQUIRE(result == Catch::Approx(0.5));
}

TEST_CASE("InverseRadialPotentialFieldCart - 3D point", "[Fields][Potential][Cartesian]") {
    // Point at (1, 1, 1), distance = sqrt(3)
    VectorN<Real, 3> x{1.0, 1.0, 1.0};
    Real result = InverseRadialPotentialFieldCart(x);
    REQUIRE(result == Catch::Approx(1.0 / std::sqrt(3.0)));
}

TEST_CASE("InverseRadialPotentialFieldCart - With constant", "[Fields][Potential][Cartesian]") {
    VectorN<Real, 3> x{2.0, 0.0, 0.0};
    Real constant = -10.0;  // e.g., -GM
    Real result = InverseRadialPotentialFieldCart(constant, x);
    REQUIRE(result == Catch::Approx(-5.0));  // -10 / 2
}

TEST_CASE("InverseRadialPotentialFieldCart - Symmetry", "[Fields][Potential][Cartesian]") {
    // Potential depends only on distance, not direction
    VectorN<Real, 3> x1{3.0, 0.0, 0.0};
    VectorN<Real, 3> x2{0.0, 3.0, 0.0};
    VectorN<Real, 3> x3{0.0, 0.0, 3.0};
    
    Real r1 = InverseRadialPotentialFieldCart(x1);
    Real r2 = InverseRadialPotentialFieldCart(x2);
    Real r3 = InverseRadialPotentialFieldCart(x3);
    
    REQUIRE(r1 == Catch::Approx(r2));
    REQUIRE(r2 == Catch::Approx(r3));
}

///////////////////////////////////////////////////////////////////////////////
//                   INVERSE RADIAL POTENTIAL - SPHERICAL                    //
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("InverseRadialPotentialFieldSpher - Unit radius", "[Fields][Potential][Spherical]") {
    // Spherical: x = (r, theta, phi), only r matters
    VectorN<Real, 3> x{1.0, 0.5, 0.0};  // r=1, arbitrary angles
    Real result = InverseRadialPotentialFieldSpher(x);
    REQUIRE(result == Catch::Approx(1.0));
}

TEST_CASE("InverseRadialPotentialFieldSpher - Radius 4", "[Fields][Potential][Spherical]") {
    VectorN<Real, 3> x{4.0, 1.0, 2.0};
    Real result = InverseRadialPotentialFieldSpher(x);
    REQUIRE(result == Catch::Approx(0.25));
}

TEST_CASE("InverseRadialPotentialFieldSpher - With constant", "[Fields][Potential][Spherical]") {
    VectorN<Real, 3> x{5.0, 0.0, 0.0};
    Real constant = -100.0;  // e.g., -GM
    Real result = InverseRadialPotentialFieldSpher(constant, x);
    REQUIRE(result == Catch::Approx(-20.0));  // -100 / 5
}

TEST_CASE("InverseRadialPotentialFieldSpher - Independent of angles", "[Fields][Potential][Spherical]") {
    // Same r, different angles should give same potential
    VectorN<Real, 3> x1{2.0, 0.0, 0.0};
    VectorN<Real, 3> x2{2.0, 1.57, 3.14};
    VectorN<Real, 3> x3{2.0, 0.5, 1.0};
    
    Real r1 = InverseRadialPotentialFieldSpher(x1);
    Real r2 = InverseRadialPotentialFieldSpher(x2);
    Real r3 = InverseRadialPotentialFieldSpher(x3);
    
    REQUIRE(r1 == Catch::Approx(r2));
    REQUIRE(r2 == Catch::Approx(r3));
}

///////////////////////////////////////////////////////////////////////////////
//                  INVERSE RADIAL POTENTIAL - CYLINDRICAL                   //
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("InverseRadialPotentialFieldCyl - Along z-axis", "[Fields][Potential][Cylindrical]") {
    // Cylindrical: x = (rho, phi, z)
    // Distance from origin = sqrt(rho² + z²)
    VectorN<Real, 3> x{0.0, 0.0, 1.0};  // rho=0, z=1
    Real result = InverseRadialPotentialFieldCyl(x);
    REQUIRE(result == Catch::Approx(1.0));
}

TEST_CASE("InverseRadialPotentialFieldCyl - In xy plane", "[Fields][Potential][Cylindrical]") {
    VectorN<Real, 3> x{2.0, 0.0, 0.0};  // rho=2, z=0
    Real result = InverseRadialPotentialFieldCyl(x);
    REQUIRE(result == Catch::Approx(0.5));  // 1/2
}

TEST_CASE("InverseRadialPotentialFieldCyl - General point", "[Fields][Potential][Cylindrical]") {
    VectorN<Real, 3> x{3.0, 1.0, 4.0};  // rho=3, z=4, distance=5
    Real result = InverseRadialPotentialFieldCyl(x);
    REQUIRE(result == Catch::Approx(0.2));  // 1/5
}

TEST_CASE("InverseRadialPotentialFieldCyl - With constant", "[Fields][Potential][Cylindrical]") {
    VectorN<Real, 3> x{3.0, 0.0, 4.0};  // distance = 5
    Real constant = -50.0;
    Real result = InverseRadialPotentialFieldCyl(constant, x);
    REQUIRE(result == Catch::Approx(-10.0));  // -50 / 5
}

///////////////////////////////////////////////////////////////////////////////
//                    INVERSE RADIAL FORCE - CARTESIAN                       //
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("InverseRadialPotentialForceFieldCart - Along x-axis", "[Fields][Force][Cartesian]") {
    VectorN<Real, 3> x{2.0, 0.0, 0.0};
    VectorN<Real, 3> F = InverseRadialPotentialForceFieldCart(x);
    
    // F = -r/|r|³ = -1/r² in direction of -r_hat
    // At (2,0,0): F = -(2,0,0)/8 = (-0.25, 0, 0)
    REQUIRE(F[0] == Catch::Approx(-0.25));
    REQUIRE(F[1] == Catch::Approx(0.0).margin(TOL(1e-14, 1e-5)));
    REQUIRE(F[2] == Catch::Approx(0.0).margin(TOL(1e-14, 1e-5)));
}

TEST_CASE("InverseRadialPotentialForceFieldCart - Unit distance", "[Fields][Force][Cartesian]") {
    VectorN<Real, 3> x{1.0, 0.0, 0.0};
    VectorN<Real, 3> F = InverseRadialPotentialForceFieldCart(x);
    
    // F = -r/|r|³ = -r/1 = -r
    REQUIRE(F[0] == Catch::Approx(-1.0));
    REQUIRE(F[1] == Catch::Approx(0.0).margin(TOL(1e-14, 1e-5)));
    REQUIRE(F[2] == Catch::Approx(0.0).margin(TOL(1e-14, 1e-5)));
}

TEST_CASE("InverseRadialPotentialForceFieldCart - 3D point", "[Fields][Force][Cartesian]") {
    VectorN<Real, 3> x{1.0, 1.0, 1.0};
    VectorN<Real, 3> F = InverseRadialPotentialForceFieldCart(x);
    
    // |r| = sqrt(3), |r|³ = 3*sqrt(3)
    // F = -r / |r|³ = -(1,1,1) / (3*sqrt(3))
    Real expected = -1.0 / (3.0 * std::sqrt(3.0));
    REQUIRE(F[0] == Catch::Approx(expected));
    REQUIRE(F[1] == Catch::Approx(expected));
    REQUIRE(F[2] == Catch::Approx(expected));
}

TEST_CASE("InverseRadialPotentialForceFieldCart - With constant", "[Fields][Force][Cartesian]") {
    VectorN<Real, 3> x{1.0, 0.0, 0.0};
    Real constant = 5.0;
    VectorN<Real, 3> F = InverseRadialPotentialForceFieldCart(constant, x);
    
    // F = -C * r / |r|³ = -5 * (1,0,0) / 1 = (-5, 0, 0)
    REQUIRE(F[0] == Catch::Approx(-5.0));
    REQUIRE(F[1] == Catch::Approx(0.0).margin(TOL(1e-14, 1e-5)));
    REQUIRE(F[2] == Catch::Approx(0.0).margin(TOL(1e-14, 1e-5)));
}

TEST_CASE("InverseRadialPotentialForceFieldCart - Points radially inward", "[Fields][Force][Cartesian]") {
    // For positive (attractive) potential, force points toward origin
    VectorN<Real, 3> x{3.0, 4.0, 0.0};
    VectorN<Real, 3> F = InverseRadialPotentialForceFieldCart(x);
    
    // Force should be anti-parallel to position vector
    // F · x < 0 means inward pointing
    Real dot = F[0] * x[0] + F[1] * x[1] + F[2] * x[2];
    REQUIRE(dot < 0);
}

///////////////////////////////////////////////////////////////////////////////
//                     INVERSE RADIAL FORCE - SPHERICAL                      //
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("InverseRadialPotentialForceFieldSph - Unit radius", "[Fields][Force][Spherical]") {
    VectorN<Real, 3> x{1.0, 0.5, 0.0};  // r=1
    VectorN<Real, 3> F = InverseRadialPotentialForceFieldSph(x);
    
    // F = (-1/r², 0, 0) in (e_r, e_theta, e_phi) basis
    REQUIRE(F[0] == Catch::Approx(-1.0));
    REQUIRE(F[1] == Catch::Approx(0.0));
    REQUIRE(F[2] == Catch::Approx(0.0));
}

TEST_CASE("InverseRadialPotentialForceFieldSph - Radius 2", "[Fields][Force][Spherical]") {
    VectorN<Real, 3> x{2.0, 1.0, 2.0};  // r=2
    VectorN<Real, 3> F = InverseRadialPotentialForceFieldSph(x);
    
    // F = (-1/4, 0, 0)
    REQUIRE(F[0] == Catch::Approx(-0.25));
    REQUIRE(F[1] == Catch::Approx(0.0));
    REQUIRE(F[2] == Catch::Approx(0.0));
}

TEST_CASE("InverseRadialPotentialForceFieldSph - With constant", "[Fields][Force][Spherical]") {
    VectorN<Real, 3> x{2.0, 0.0, 0.0};  // r=2
    Real constant = -4.0;  // Repulsive
    VectorN<Real, 3> F = InverseRadialPotentialForceFieldSph(constant, x);
    
    // F = (4/4, 0, 0) = (1, 0, 0) - outward for negative constant
    REQUIRE(F[0] == Catch::Approx(1.0));
    REQUIRE(F[1] == Catch::Approx(0.0));
    REQUIRE(F[2] == Catch::Approx(0.0));
}

TEST_CASE("InverseRadialPotentialForceFieldSph - Purely radial", "[Fields][Force][Spherical]") {
    // Force has no angular components regardless of position
    VectorN<Real, 3> x{3.0, 1.5, 2.0};
    VectorN<Real, 3> F = InverseRadialPotentialForceFieldSph(x);
    
    REQUIRE(F[1] == Catch::Approx(0.0));  // No theta component
    REQUIRE(F[2] == Catch::Approx(0.0));  // No phi component
}

///////////////////////////////////////////////////////////////////////////////
//                    INVERSE RADIAL FIELD CLASS - CART                      //
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("InverseRadialFieldCart - Default constructor", "[Fields][Class][Cartesian]") {
    InverseRadialFieldCart field;
    VectorN<Real, 3> x{1.0, 0.0, 0.0};
    
    // Default constant is -1.0, so result = -1 * 1/|r| = -1
    Real result = field(x);
    REQUIRE(result == Catch::Approx(-1.0));
}

TEST_CASE("InverseRadialFieldCart - Custom constant", "[Fields][Class][Cartesian]") {
    InverseRadialFieldCart field(10.0);  // Positive constant
    VectorN<Real, 3> x{2.0, 0.0, 0.0};
    
    // Result = 10 * 1/2 = 5
    Real result = field(x);
    REQUIRE(result == Catch::Approx(5.0));
}

TEST_CASE("InverseRadialFieldCart - Gravity-like", "[Fields][Class][Cartesian]") {
    // For gravity: Φ = -GM/r, so constant = -GM = -1 (normalized)
    InverseRadialFieldCart gravity(-1.0);
    VectorN<Real, 3> x{4.0, 0.0, 0.0};
    
    // Result = -1 * 1/4 = -0.25
    Real result = gravity(x);
    REQUIRE(result == Catch::Approx(-0.25));
}

TEST_CASE("InverseRadialFieldCart - 3D evaluation", "[Fields][Class][Cartesian]") {
    InverseRadialFieldCart field(1.0);
    VectorN<Real, 3> x{1.0, 2.0, 2.0};  // |r| = 3
    
    Real result = field(x);
    REQUIRE(result == Catch::Approx(1.0 / 3.0));
}

///////////////////////////////////////////////////////////////////////////////
//                   INVERSE RADIAL FIELD CLASS - SPHER                      //
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("InverseRadialFieldSpher - Default constructor", "[Fields][Class][Spherical]") {
    InverseRadialFieldSpher field;
    VectorN<Real, 3> x{1.0, 0.5, 1.0};  // r=1
    
    // Default constant is -1.0, so result = -1 * 1/r = -1
    Real result = field(x);
    REQUIRE(result == Catch::Approx(-1.0));
}

TEST_CASE("InverseRadialFieldSpher - Custom constant", "[Fields][Class][Spherical]") {
    InverseRadialFieldSpher field(8.0);
    VectorN<Real, 3> x{4.0, 0.0, 0.0};  // r=4
    
    // Result = 8 * 1/4 = 2
    Real result = field(x);
    REQUIRE(result == Catch::Approx(2.0));
}

TEST_CASE("InverseRadialFieldSpher - Coulomb-like", "[Fields][Class][Spherical]") {
    // For Coulomb: Φ = kQ/r, so constant = kQ = 1 (normalized)
    InverseRadialFieldSpher coulomb(1.0);
    VectorN<Real, 3> x{5.0, 1.0, 2.0};  // r=5
    
    Real result = coulomb(x);
    REQUIRE(result == Catch::Approx(0.2));
}

///////////////////////////////////////////////////////////////////////////////
//                  INVERSE RADIAL FORCE FIELD CLASS - CART                  //
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("InverseRadialForceFieldCart - Default constructor", "[Fields][Class][Cartesian]") {
    InverseRadialForceFieldCart field;
    VectorN<Real, 3> x{1.0, 0.0, 0.0};
    
    // Default constant is -1.0
    // F = -1 * (-r/|r|³) = r/|r|³ = (1,0,0)/1 = (1,0,0)
    VectorN<Real, 3> F = field(x);
    REQUIRE(F[0] == Catch::Approx(1.0));
    REQUIRE(F[1] == Catch::Approx(0.0).margin(TOL(1e-14, 1e-5)));
    REQUIRE(F[2] == Catch::Approx(0.0).margin(TOL(1e-14, 1e-5)));
}

TEST_CASE("InverseRadialForceFieldCart - Custom constant", "[Fields][Class][Cartesian]") {
    InverseRadialForceFieldCart field(2.0);
    VectorN<Real, 3> x{1.0, 0.0, 0.0};
    
    // F = 2 * (-r/|r|³) = -2*(1,0,0) = (-2,0,0)
    VectorN<Real, 3> F = field(x);
    REQUIRE(F[0] == Catch::Approx(-2.0));
    REQUIRE(F[1] == Catch::Approx(0.0).margin(TOL(1e-14, 1e-5)));
    REQUIRE(F[2] == Catch::Approx(0.0).margin(TOL(1e-14, 1e-5)));
}

TEST_CASE("InverseRadialForceFieldCart - Inverse square law", "[Fields][Class][Cartesian]") {
    InverseRadialForceFieldCart field(1.0);
    
    // Force magnitude should decrease as 1/r²
    VectorN<Real, 3> x1{1.0, 0.0, 0.0};
    VectorN<Real, 3> x2{2.0, 0.0, 0.0};
    
    VectorN<Real, 3> F1 = field(x1);
    VectorN<Real, 3> F2 = field(x2);
    
    Real mag1 = std::abs(F1[0]);  // Force along x only
    Real mag2 = std::abs(F2[0]);
    
    // |F1|/|F2| should be 4 (inverse square)
    REQUIRE(mag1 / mag2 == Catch::Approx(4.0));
}

///////////////////////////////////////////////////////////////////////////////
//                 INVERSE RADIAL FORCE FIELD CLASS - SPHER                  //
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("InverseRadialForceFieldSpher - Default constructor", "[Fields][Class][Spherical]") {
    InverseRadialForceFieldSpher field;
    VectorN<Real, 3> x{1.0, 0.5, 1.0};  // r=1
    
    // Default constant is -1.0
    // F = -1 * (-1/r², 0, 0) = (1/r², 0, 0) = (1, 0, 0)
    VectorN<Real, 3> F = field(x);
    REQUIRE(F[0] == Catch::Approx(1.0));
    REQUIRE(F[1] == Catch::Approx(0.0));
    REQUIRE(F[2] == Catch::Approx(0.0));
}

TEST_CASE("InverseRadialForceFieldSpher - Custom constant", "[Fields][Class][Spherical]") {
    InverseRadialForceFieldSpher field(3.0);
    VectorN<Real, 3> x{2.0, 0.0, 0.0};  // r=2
    
    // F = 3 * (-1/4, 0, 0) = (-0.75, 0, 0)
    VectorN<Real, 3> F = field(x);
    REQUIRE(F[0] == Catch::Approx(-0.75));
    REQUIRE(F[1] == Catch::Approx(0.0));
    REQUIRE(F[2] == Catch::Approx(0.0));
}

TEST_CASE("InverseRadialForceFieldSpher - Inverse square law", "[Fields][Class][Spherical]") {
    InverseRadialForceFieldSpher field(1.0);
    
    VectorN<Real, 3> x1{1.0, 0.5, 1.0};  // r=1
    VectorN<Real, 3> x2{3.0, 0.5, 1.0};  // r=3
    
    VectorN<Real, 3> F1 = field(x1);
    VectorN<Real, 3> F2 = field(x2);
    
    Real mag1 = std::abs(F1[0]);
    Real mag2 = std::abs(F2[0]);
    
    // |F1|/|F2| should be 9 (inverse square: 3² = 9)
    REQUIRE(mag1 / mag2 == Catch::Approx(9.0));
}

///////////////////////////////////////////////////////////////////////////////
//                          PHYSICAL CONSISTENCY                             //
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("Fields - Force is negative gradient of potential", "[Fields][Physics]") {
    // F = -∇Φ
    // For Φ = 1/r in Cartesian: ∇Φ = -r/|r|³
    // So F = -∇Φ = r/|r|³
    // But InverseRadialPotentialForceFieldCart uses F = -r/|r|³ (attractive)
    // This is F = ∇Φ (not -∇Φ), representing attractive force toward origin
    
    VectorN<Real, 3> x{Real(2.0), Real(0.0), Real(0.0)};
    
    // Numerical gradient of potential (Φ = 1/|r|)
    Real h = TOL(1e-6, 1e-3);
    VectorN<Real, 3> xph{Real(2.0) + h, Real(0.0), Real(0.0)};
    VectorN<Real, 3> xmh{Real(2.0) - h, Real(0.0), Real(0.0)};
    Real phi_ph = InverseRadialPotentialFieldCart(xph);
    Real phi_mh = InverseRadialPotentialFieldCart(xmh);
    Real dPhi_dx = (phi_ph - phi_mh) / (2 * h);
    
    // Force field F = -r/|r|³ = gradient of Φ (attractive force)
    VectorN<Real, 3> F = InverseRadialPotentialForceFieldCart(x);
    
    // The force field implementation gives F = ∇Φ (attractive toward origin)
    REQUIRE(F[0] == Catch::Approx(dPhi_dx).epsilon(TOL(1e-4, 1e-2)));
}

TEST_CASE("Fields - Spherical and Cartesian match at same point", "[Fields][Physics]") {
    // At (r, theta, phi) = (2, π/2, 0), Cartesian = (2, 0, 0)
    VectorN<Real, 3> cart{2.0, 0.0, 0.0};
    VectorN<Real, 3> spher{2.0, 1.5707963267948966, 0.0};  // r=2, theta=π/2
    
    InverseRadialFieldCart potCart(1.0);
    InverseRadialFieldSpher potSpher(1.0);
    
    Real phiCart = potCart(cart);
    Real phiSpher = potSpher(spher);
    
    REQUIRE(phiCart == Catch::Approx(phiSpher));
}

TEST_CASE("Fields - Force magnitude consistent between coordinates", "[Fields][Physics]") {
    // At same physical point, force magnitude should match
    VectorN<Real, 3> cart{2.0, 0.0, 0.0};
    VectorN<Real, 3> spher{2.0, 1.5707963267948966, 0.0};
    
    InverseRadialForceFieldCart forceCart(1.0);
    InverseRadialForceFieldSpher forceSpher(1.0);
    
    VectorN<Real, 3> Fcart = forceCart(cart);
    VectorN<Real, 3> Fspher = forceSpher(spher);
    
    Real magCart = Fcart.NormL2();
    Real magSpher = std::abs(Fspher[0]);  // Only radial component
    
    REQUIRE(magCart == Catch::Approx(magSpher));
}

} // namespace MML::Tests::Core::FieldsTests

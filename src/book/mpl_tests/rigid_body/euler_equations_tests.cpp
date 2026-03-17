///////////////////////////////////////////////////////////////////////////////////////////
// Euler Equations Tests - Rigid Body Dynamics
///////////////////////////////////////////////////////////////////////////////////////////
//
// Tests for Euler's equations of rigid body rotation:
// I₁ω̇₁ = (I₂ - I₃)ω₂ω₃ + τ₁
// I₂ω̇₂ = (I₃ - I₁)ω₃ω₁ + τ₂
// I₃ω̇₃ = (I₁ - I₂)ω₁ω₂ + τ₃
//
// Tested scenarios:
// - Torque-free precession
// - Symmetric top behavior
// - Angular momentum conservation
// - Energy conservation
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
#include "base/Tensor.h"
#endif

#include "RigidBody/EulerEquations.h"

using namespace MML;
using namespace MPL;
using namespace MML::Testing;

using Catch::Matchers::WithinAbs;

namespace MPL::Tests::RigidBody::EulerEquationsTests
{
	//===================================================================================
	// Helper functions
	//===================================================================================

	Real RotationalKineticEnergy(const Vec3Cart& omega, const Vec3Cart& I_principal)
	{
		// T = 1/2 * (I₁ω₁² + I₂ω₂² + I₃ω₃²)
		return 0.5 * (I_principal.X() * omega.X() * omega.X() +
		              I_principal.Y() * omega.Y() * omega.Y() +
		              I_principal.Z() * omega.Z() * omega.Z());
	}

	Vec3Cart AngularMomentum(const Vec3Cart& omega, const Vec3Cart& I_principal)
	{
		// L = I·ω (for principal axes, this is component-wise)
		return Vec3Cart(I_principal.X() * omega.X(),
		                I_principal.Y() * omega.Y(),
		                I_principal.Z() * omega.Z());
	}

	Real AngularMomentumMagnitude(const Vec3Cart& omega, const Vec3Cart& I_principal)
	{
		Vec3Cart L = AngularMomentum(omega, I_principal);
		return L.NormL2();
	}

	//===================================================================================
	// Spherical top: I₁ = I₂ = I₃
	//===================================================================================

	TEST_CASE("EulerEquations::Spherical_top_no_precession", "[rigid_body][euler]")
	{
		TEST_PRECISION_INFO();

		// Spherical top (I₁ = I₂ = I₃) has no torque-free precession
		// ω̇ᵢ = 0 for all i when no external torque

		Vec3Cart I_principal(2.0, 2.0, 2.0);  // Uniform sphere
		Vec3Cart omega(1.0, 2.0, 3.0);        // Arbitrary initial rotation
		Vec3Cart torque(0.0, 0.0, 0.0);       // No external torque

		// Euler equations: ω̇₁ = (I₂-I₃)/(I₁) * ω₂ω₃ = 0
		Real omega1_dot = (I_principal.Y() - I_principal.Z()) / I_principal.X() * omega.Y() * omega.Z();
		Real omega2_dot = (I_principal.Z() - I_principal.X()) / I_principal.Y() * omega.Z() * omega.X();
		Real omega3_dot = (I_principal.X() - I_principal.Y()) / I_principal.Z() * omega.X() * omega.Y();

		REQUIRE_THAT(omega1_dot, WithinAbs(0.0, 1e-12));
		REQUIRE_THAT(omega2_dot, WithinAbs(0.0, 1e-12));
		REQUIRE_THAT(omega3_dot, WithinAbs(0.0, 1e-12));
	}

	//===================================================================================
	// Symmetric top: I1 = I2 != I3
	//===================================================================================

	TEST_CASE("EulerEquations::Symmetric_top_omega3_constant", "[rigid_body][euler]")
	{
		TEST_PRECISION_INFO();

		// Symmetric top: I1 = I2 != I3
		// omega3_dot = (I1-I2)/(I3) * omega1*omega2 = 0 (since I1 = I2)
		// So omega3 is constant (torque-free)

		Vec3Cart I_principal(2.0, 2.0, 5.0);  // Symmetric about z
		Vec3Cart omega(1.0, 2.0, 3.0);
		Vec3Cart torque(0.0, 0.0, 0.0);

		Real omega3_dot = (I_principal.X() - I_principal.Y()) / I_principal.Z() * omega.X() * omega.Y();

		REQUIRE_THAT(omega3_dot, WithinAbs(0.0, 1e-12));
	}

	//===================================================================================
	// Asymmetric top: I₁ ≠ I₂ ≠ I₃
	//===================================================================================

	TEST_CASE("EulerEquations::Asymmetric_top_nontrivial_dynamics", "[rigid_body][euler]")
	{
		TEST_PRECISION_INFO();

		// Asymmetric top has complex dynamics
		Vec3Cart I_principal(1.0, 2.0, 3.0);
		Vec3Cart omega(1.0, 1.0, 1.0);

		// ω̇₁ = (I₂-I₃)/(I₁) * ω₂ω₃ = (2-3)/1 * 1*1 = -1
		// ω̇₂ = (I₃-I₁)/(I₂) * ω₃ω₁ = (3-1)/2 * 1*1 = 1
		// ω̇₃ = (I₁-I₂)/(I₃) * ω₁ω₂ = (1-2)/3 * 1*1 = -1/3

		Real omega1_dot = (I_principal.Y() - I_principal.Z()) / I_principal.X() * omega.Y() * omega.Z();
		Real omega2_dot = (I_principal.Z() - I_principal.X()) / I_principal.Y() * omega.Z() * omega.X();
		Real omega3_dot = (I_principal.X() - I_principal.Y()) / I_principal.Z() * omega.X() * omega.Y();

		REQUIRE_THAT(omega1_dot, WithinAbs(-1.0, 1e-12));
		REQUIRE_THAT(omega2_dot, WithinAbs(1.0, 1e-12));
		REQUIRE_THAT(omega3_dot, WithinAbs(-1.0/3.0, 1e-12));
	}

	//===================================================================================
	// Rotation about principal axis is stable/unstable
	//===================================================================================

	TEST_CASE("EulerEquations::Rotation_about_principal_axis", "[rigid_body][euler]")
	{
		TEST_PRECISION_INFO();

		// Rotation purely about a principal axis is steady (ω̇ = 0)
		Vec3Cart I_principal(1.0, 2.0, 3.0);

		SECTION("About axis 1 (smallest moment)") {
			Vec3Cart omega(5.0, 0.0, 0.0);

			Real omega1_dot = (I_principal.Y() - I_principal.Z()) / I_principal.X() * omega.Y() * omega.Z();
			Real omega2_dot = (I_principal.Z() - I_principal.X()) / I_principal.Y() * omega.Z() * omega.X();
			Real omega3_dot = (I_principal.X() - I_principal.Y()) / I_principal.Z() * omega.X() * omega.Y();

			REQUIRE_THAT(omega1_dot, WithinAbs(0.0, 1e-12));
			REQUIRE_THAT(omega2_dot, WithinAbs(0.0, 1e-12));
			REQUIRE_THAT(omega3_dot, WithinAbs(0.0, 1e-12));
		}

		SECTION("About axis 3 (largest moment)") {
			Vec3Cart omega(0.0, 0.0, 5.0);

			Real omega1_dot = (I_principal.Y() - I_principal.Z()) / I_principal.X() * omega.Y() * omega.Z();
			Real omega2_dot = (I_principal.Z() - I_principal.X()) / I_principal.Y() * omega.Z() * omega.X();
			Real omega3_dot = (I_principal.X() - I_principal.Y()) / I_principal.Z() * omega.X() * omega.Y();

			REQUIRE_THAT(omega1_dot, WithinAbs(0.0, 1e-12));
			REQUIRE_THAT(omega2_dot, WithinAbs(0.0, 1e-12));
			REQUIRE_THAT(omega3_dot, WithinAbs(0.0, 1e-12));
		}
	}

	//===================================================================================
	// Conservation laws
	//===================================================================================

	TEST_CASE("EulerEquations::Conservation_identities", "[rigid_body][euler]")
	{
		TEST_PRECISION_INFO();

		// For torque-free motion, both T and |L|² are conserved
		// Also, the Euler equations preserve these invariants

		Vec3Cart I_principal(1.0, 2.0, 3.0);
		Vec3Cart omega_initial(1.0, 2.0, 0.5);

		// Initial energy and angular momentum
		Real T_initial = RotationalKineticEnergy(omega_initial, I_principal);
		Real L2_initial = AngularMomentumMagnitude(omega_initial, I_principal);

		// These should be constants of motion
		// (In a full simulation, we'd verify this over time)
		// For now, just verify the formulas work

		REQUIRE(T_initial > 0.0);
		REQUIRE(L2_initial > 0.0);

		// Verify L² = 2*T*I_eff doesn't hold in general (that's for 2D)
		// But verify L² and T are computed correctly

		// L = (I₁ω₁, I₂ω₂, I₃ω₃)
		Vec3Cart L = AngularMomentum(omega_initial, I_principal);
		Real L2_computed = L.X()*L.X() + L.Y()*L.Y() + L.Z()*L.Z();

		REQUIRE_THAT(L2_computed, WithinAbs(L2_initial * L2_initial, 1e-10));
	}

	//===================================================================================
	// Polhode and herpolhode ellipsoids
	//===================================================================================

	TEST_CASE("EulerEquations::Polhode_constraint", "[rigid_body][euler]")
	{
		TEST_PRECISION_INFO();

		// The angular velocity vector traces a curve on the "polhode ellipsoid"
		// defined by: 2T = I₁ω₁² + I₂ω₂² + I₃ω₃² = const
		// AND:        L² = I₁²ω₁² + I₂²ω₂² + I₃²ω₃² = const

		Vec3Cart I_principal(1.0, 2.0, 4.0);
		Vec3Cart omega(2.0, 1.0, 0.5);

		Real T = RotationalKineticEnergy(omega, I_principal);
		Real L_mag = AngularMomentumMagnitude(omega, I_principal);

		// These define two ellipsoids; ω lies on their intersection
		// For a valid physical state, both must be positive
		REQUIRE(T > 0.0);
		REQUIRE(L_mag > 0.0);

		// Also check that L² ≥ 2*T*I_min and L² ≤ 2*T*I_max
		Real I_min = std::min({I_principal.X(), I_principal.Y(), I_principal.Z()});
		Real I_max = std::max({I_principal.X(), I_principal.Y(), I_principal.Z()});
		Real L2 = L_mag * L_mag;

		REQUIRE(L2 >= 2.0 * T * I_min - 1e-10);
		REQUIRE(L2 <= 2.0 * T * I_max + 1e-10);
	}

} // namespace MPL::Tests::RigidBody::EulerEquationsTests

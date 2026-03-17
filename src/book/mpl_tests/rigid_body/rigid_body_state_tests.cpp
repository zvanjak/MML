///////////////////////////////////////////////////////////////////////////////////////////
// Rigid Body State Tests
///////////////////////////////////////////////////////////////////////////////////////////
//
// Tests for RigidBodyState structure:
// - State initialization
// - Quaternion normalization
// - State vector conversion (for ODE integration)
// - Kinetic energy calculations
// - Angular momentum calculations
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
#include "base/Quaternions.h"
#endif

#include "RigidBody/RigidBodyState.h"

using namespace MML;
using namespace MPL;
using namespace MML::Testing;

using Catch::Matchers::WithinAbs;

namespace MPL::Tests::RigidBody::RigidBodyStateTests
{
	//===================================================================================
	// Default construction
	//===================================================================================

	TEST_CASE("RigidBodyState::Default_construction", "[rigid_body][state]")
	{
		TEST_PRECISION_INFO();

		RigidBodyState state;

		// Position at origin
		REQUIRE_THAT(state.position.X(), WithinAbs(0.0, 1e-12));
		REQUIRE_THAT(state.position.Y(), WithinAbs(0.0, 1e-12));
		REQUIRE_THAT(state.position.Z(), WithinAbs(0.0, 1e-12));

		// Zero velocity
		REQUIRE_THAT(state.velocity.X(), WithinAbs(0.0, 1e-12));
		REQUIRE_THAT(state.velocity.Y(), WithinAbs(0.0, 1e-12));
		REQUIRE_THAT(state.velocity.Z(), WithinAbs(0.0, 1e-12));

		// Identity orientation
		REQUIRE_THAT(state.orientation.w(), WithinAbs(1.0, 1e-12));
		REQUIRE_THAT(state.orientation.x(), WithinAbs(0.0, 1e-12));
		REQUIRE_THAT(state.orientation.y(), WithinAbs(0.0, 1e-12));
		REQUIRE_THAT(state.orientation.z(), WithinAbs(0.0, 1e-12));

		// Zero angular velocity
		REQUIRE_THAT(state.angularVel.X(), WithinAbs(0.0, 1e-12));
		REQUIRE_THAT(state.angularVel.Y(), WithinAbs(0.0, 1e-12));
		REQUIRE_THAT(state.angularVel.Z(), WithinAbs(0.0, 1e-12));
	}

	//===================================================================================
	// Quaternion safe normalization
	//===================================================================================

	TEST_CASE("RigidBodyState::SafeNormalized_normal_quaternion", "[rigid_body][state][quaternion]")
	{
		TEST_PRECISION_INFO();

		Quaternion q(1.0, 2.0, 3.0, 4.0);  // Not normalized
		Quaternion q_norm = SafeNormalized(q);

		Real norm = std::sqrt(q_norm.w()*q_norm.w() + q_norm.x()*q_norm.x() + 
		                      q_norm.y()*q_norm.y() + q_norm.z()*q_norm.z());

		REQUIRE_THAT(norm, WithinAbs(1.0, 1e-12));
	}

	TEST_CASE("RigidBodyState::SafeNormalized_near_zero", "[rigid_body][state][quaternion]")
	{
		TEST_PRECISION_INFO();

		// Near-zero quaternion should return identity
		Quaternion q(1e-14, 1e-14, 1e-14, 1e-14);
		Quaternion q_norm = SafeNormalized(q);

		REQUIRE_THAT(q_norm.w(), WithinAbs(1.0, 1e-12));
		REQUIRE_THAT(q_norm.x(), WithinAbs(0.0, 1e-12));
		REQUIRE_THAT(q_norm.y(), WithinAbs(0.0, 1e-12));
		REQUIRE_THAT(q_norm.z(), WithinAbs(0.0, 1e-12));
	}

	TEST_CASE("RigidBodyState::SafeNormalized_already_normalized", "[rigid_body][state][quaternion]")
	{
		TEST_PRECISION_INFO();

		// Create a normalized quaternion (90° rotation about z-axis)
		Real angle = Constants::PI / 2.0;
		Quaternion q(std::cos(angle/2), 0.0, 0.0, std::sin(angle/2));
		Quaternion q_norm = SafeNormalized(q);

		REQUIRE_THAT(q_norm.w(), WithinAbs(q.w(), 1e-12));
		REQUIRE_THAT(q_norm.x(), WithinAbs(q.x(), 1e-12));
		REQUIRE_THAT(q_norm.y(), WithinAbs(q.y(), 1e-12));
		REQUIRE_THAT(q_norm.z(), WithinAbs(q.z(), 1e-12));
	}

	//===================================================================================
	// In-place normalization
	//===================================================================================

	TEST_CASE("RigidBodyState::SafeNormalize_inplace", "[rigid_body][state][quaternion]")
	{
		TEST_PRECISION_INFO();

		Quaternion q(3.0, 4.0, 0.0, 0.0);
		SafeNormalize(q);

		Real norm = std::sqrt(q.w()*q.w() + q.x()*q.x() + q.y()*q.y() + q.z()*q.z());
		REQUIRE_THAT(norm, WithinAbs(1.0, 1e-12));

		// Check direction preserved (3:4 ratio)
		REQUIRE_THAT(q.w() / q.x(), WithinAbs(3.0/4.0, 1e-10));
	}

	//===================================================================================
	// Rotation quaternions
	//===================================================================================

	TEST_CASE("RigidBodyState::Quaternion_rotations", "[rigid_body][state][quaternion]")
	{
		TEST_PRECISION_INFO();

		SECTION("90° rotation about z-axis") {
			Real angle = Constants::PI / 2.0;
			Quaternion q(std::cos(angle/2), 0.0, 0.0, std::sin(angle/2));

			// Rotate vector (1, 0, 0) should give (0, 1, 0)
			Vec3Cart v(1.0, 0.0, 0.0);
			Vec3Cart v_rotated = q.Rotate(v);

			REQUIRE_THAT(v_rotated.X(), WithinAbs(0.0, 1e-10));
			REQUIRE_THAT(v_rotated.Y(), WithinAbs(1.0, 1e-10));
			REQUIRE_THAT(v_rotated.Z(), WithinAbs(0.0, 1e-10));
		}

		SECTION("180° rotation about x-axis") {
			Real angle = Constants::PI;
			Quaternion q(std::cos(angle/2), std::sin(angle/2), 0.0, 0.0);

			// Rotate vector (0, 1, 0) should give (0, -1, 0)
			Vec3Cart v(0.0, 1.0, 0.0);
			Vec3Cart v_rotated = q.Rotate(v);

			REQUIRE_THAT(v_rotated.X(), WithinAbs(0.0, 1e-10));
			REQUIRE_THAT(v_rotated.Y(), WithinAbs(-1.0, 1e-10));
			REQUIRE_THAT(v_rotated.Z(), WithinAbs(0.0, 1e-10));
		}

		SECTION("Identity rotation") {
			Quaternion q = Quaternion::Identity();

			Vec3Cart v(1.0, 2.0, 3.0);
			Vec3Cart v_rotated = q.Rotate(v);

			REQUIRE_THAT(v_rotated.X(), WithinAbs(v.X(), 1e-12));
			REQUIRE_THAT(v_rotated.Y(), WithinAbs(v.Y(), 1e-12));
			REQUIRE_THAT(v_rotated.Z(), WithinAbs(v.Z(), 1e-12));
		}
	}

	//===================================================================================
	// Quaternion composition
	//===================================================================================

	TEST_CASE("RigidBodyState::Quaternion_composition", "[rigid_body][state][quaternion]")
	{
		TEST_PRECISION_INFO();

		// Two 90° rotations about z-axis should equal one 180° rotation
		Real angle_90 = Constants::PI / 2.0;
		Real angle_180 = Constants::PI;

		Quaternion q_90(std::cos(angle_90/2), 0.0, 0.0, std::sin(angle_90/2));
		Quaternion q_180(std::cos(angle_180/2), 0.0, 0.0, std::sin(angle_180/2));

		// Compose two 90° rotations
		Quaternion q_composed = q_90 * q_90;
		q_composed = SafeNormalized(q_composed);

		// Apply to a test vector
		Vec3Cart v(1.0, 0.0, 0.0);
		Vec3Cart v_composed = q_composed.Rotate(v);
		Vec3Cart v_direct = q_180.Rotate(v);

		REQUIRE_THAT(v_composed.X(), WithinAbs(v_direct.X(), 1e-10));
		REQUIRE_THAT(v_composed.Y(), WithinAbs(v_direct.Y(), 1e-10));
		REQUIRE_THAT(v_composed.Z(), WithinAbs(v_direct.Z(), 1e-10));
	}

} // namespace MPL::Tests::RigidBody::RigidBodyStateTests

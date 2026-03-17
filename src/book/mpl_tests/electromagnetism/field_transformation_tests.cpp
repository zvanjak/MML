///////////////////////////////////////////////////////////////////////////////////////////
// Electromagnetic Field Transformation Tests
///////////////////////////////////////////////////////////////////////////////////////////
//
// Tests for how E and B fields transform under Lorentz boosts:
//
// For a boost along x-axis with velocity v:
//   E'_x = E_x                    B'_x = B_x
//   E'_y = γ(E_y - vB_z)          B'_y = γ(B_y + vE_z/c²)
//   E'_z = γ(E_z + vB_y)          B'_z = γ(B_z - vE_y/c²)
//
// With c = 1 (natural units).
//
// Also verifies that Lorentz invariants are preserved:
// - E² - c²B² (with c=1: E² - B²)
// - E·B
//
///////////////////////////////////////////////////////////////////////////////////////////

#include <catch2/catch_all.hpp>

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#include "base/Vector/VectorTypes.h"
#include "base/Tensor.h"
#endif

#include "Electromagnetism/EMTensor.h"
#include "SpecialRelativity/LorentzTransformation.h"

using namespace MML;
using namespace MPL;

using Catch::Matchers::WithinAbs;

namespace MPL::Tests::Electromagnetism::FieldTransformationTests
{
	//===================================================================================
	// Helper: Transform E and B fields under Lorentz boost along x-axis
	//===================================================================================

	struct FieldTransformResult {
		Vector3Cartesian E_prime;
		Vector3Cartesian B_prime;
	};

	FieldTransformResult TransformFieldsBoostX(const Vector3Cartesian& E, 
	                                            const Vector3Cartesian& B, 
	                                            Real v)
	{
		// Lorentz transformation of EM fields for boost along x with velocity v
		// Using c = 1 (natural units)

		if (v <= -1.0 || v >= 1.0)
			throw std::range_error("Invalid velocity for field transformation");

		Real gamma = 1.0 / std::sqrt(1.0 - v * v);

		Vector3Cartesian E_prime(
			E.X(),                           // E'_x = E_x
			gamma * (E.Y() - v * B.Z()),     // E'_y = γ(E_y - vB_z)
			gamma * (E.Z() + v * B.Y())      // E'_z = γ(E_z + vB_y)
		);

		Vector3Cartesian B_prime(
			B.X(),                           // B'_x = B_x
			gamma * (B.Y() + v * E.Z()),     // B'_y = γ(B_y + vE_z)
			gamma * (B.Z() - v * E.Y())      // B'_z = γ(B_z - vE_y)
		);

		return { E_prime, B_prime };
	}

	Real ComputeInvariant1(const Vector3Cartesian& E, const Vector3Cartesian& B)
	{
		// E² - B² (with c = 1)
		return E.NormL2() * E.NormL2() - B.NormL2() * B.NormL2();
	}

	Real ComputeInvariant2(const Vector3Cartesian& E, const Vector3Cartesian& B)
	{
		// E·B
		return E.X() * B.X() + E.Y() * B.Y() + E.Z() * B.Z();
	}

	//===================================================================================
	// Parallel field components unchanged
	//===================================================================================

	TEST_CASE("FieldTransform::Parallel_components_unchanged", "[electromagnetism][transform]")
	{
		// Field components parallel to boost direction are unchanged

		Vector3Cartesian E(5.0, 2.0, 3.0);
		Vector3Cartesian B(4.0, 1.0, 2.0);
		Real v = 0.6;

		auto result = TransformFieldsBoostX(E, B, v);

		// E_x and B_x should be unchanged
		REQUIRE_THAT(result.E_prime.X(), WithinAbs(E.X(), 1e-12));
		REQUIRE_THAT(result.B_prime.X(), WithinAbs(B.X(), 1e-12));
	}

	//===================================================================================
	// Perpendicular components mix E and B
	//===================================================================================

	TEST_CASE("FieldTransform::Perpendicular_components_mix", "[electromagnetism][transform]")
	{
		// E'_y = γ(E_y - vB_z), E'_z = γ(E_z + vB_y)
		// B'_y = γ(B_y + vE_z), B'_z = γ(B_z - vE_y)

		Vector3Cartesian E(0.0, 1.0, 0.0);  // Pure E_y
		Vector3Cartesian B(0.0, 0.0, 0.0);
		Real v = 0.6;
		Real gamma = 1.0 / std::sqrt(1.0 - v * v);  // γ = 1.25

		auto result = TransformFieldsBoostX(E, B, v);

		// E'_y = γ * E_y = 1.25 * 1.0 = 1.25
		REQUIRE_THAT(result.E_prime.Y(), WithinAbs(gamma * 1.0, 1e-12));

		// B'_z = γ * (-v * E_y) = 1.25 * (-0.6 * 1.0) = -0.75
		REQUIRE_THAT(result.B_prime.Z(), WithinAbs(gamma * (-v * 1.0), 1e-12));
	}

	//===================================================================================
	// Pure electric field generates magnetic field
	//===================================================================================

	TEST_CASE("FieldTransform::Pure_E_generates_B", "[electromagnetism][transform]")
	{
		// A pure electric field in one frame has a magnetic component in another

		Vector3Cartesian E(0.0, 0.0, 1.0);  // E_z
		Vector3Cartesian B(0.0, 0.0, 0.0);
		Real v = 0.8;
		Real gamma = 1.0 / std::sqrt(1.0 - v * v);  // γ ≈ 1.667

		auto result = TransformFieldsBoostX(E, B, v);

		// E'_z = γ * E_z
		REQUIRE_THAT(result.E_prime.Z(), WithinAbs(gamma * E.Z(), 1e-12));

		// B'_y = γ * v * E_z (moving charge creates magnetic field)
		REQUIRE_THAT(result.B_prime.Y(), WithinAbs(gamma * v * E.Z(), 1e-12));
	}

	//===================================================================================
	// Pure magnetic field generates electric field
	//===================================================================================

	TEST_CASE("FieldTransform::Pure_B_generates_E", "[electromagnetism][transform]")
	{
		// A pure magnetic field in one frame has an electric component in another

		Vector3Cartesian E(0.0, 0.0, 0.0);
		Vector3Cartesian B(0.0, 1.0, 0.0);  // B_y
		Real v = 0.5;
		Real gamma = 1.0 / std::sqrt(1.0 - v * v);

		auto result = TransformFieldsBoostX(E, B, v);

		// E'_z = γ * v * B_y
		REQUIRE_THAT(result.E_prime.Z(), WithinAbs(gamma * v * B.Y(), 1e-12));

		// B'_y = γ * B_y
		REQUIRE_THAT(result.B_prime.Y(), WithinAbs(gamma * B.Y(), 1e-12));
	}

	//===================================================================================
	// First Lorentz invariant: E² - B²
	//===================================================================================

	TEST_CASE("FieldTransform::Invariant1_E2_minus_B2", "[electromagnetism][transform][invariant]")
	{
		Vector3Cartesian E(1.0, 2.0, 3.0);
		Vector3Cartesian B(0.5, 1.0, 1.5);

		std::vector<Real> velocities = { 0.0, 0.3, 0.6, -0.5, 0.9, -0.9 };

		Real I1_original = ComputeInvariant1(E, B);

		for (Real v : velocities) {
			auto result = TransformFieldsBoostX(E, B, v);
			Real I1_transformed = ComputeInvariant1(result.E_prime, result.B_prime);

			DYNAMIC_SECTION("v = " << v) {
				REQUIRE_THAT(I1_transformed, WithinAbs(I1_original, 1e-10));
			}
		}
	}

	//===================================================================================
	// Second Lorentz invariant: E·B
	//===================================================================================

	TEST_CASE("FieldTransform::Invariant2_E_dot_B", "[electromagnetism][transform][invariant]")
	{
		Vector3Cartesian E(1.0, 2.0, 3.0);
		Vector3Cartesian B(0.5, 1.0, 1.5);

		std::vector<Real> velocities = { 0.0, 0.3, 0.6, -0.5, 0.9, -0.9 };

		Real I2_original = ComputeInvariant2(E, B);

		for (Real v : velocities) {
			auto result = TransformFieldsBoostX(E, B, v);
			Real I2_transformed = ComputeInvariant2(result.E_prime, result.B_prime);

			DYNAMIC_SECTION("v = " << v) {
				REQUIRE_THAT(I2_transformed, WithinAbs(I2_original, 1e-10));
			}
		}
	}

	//===================================================================================
	// Plane wave remains plane wave
	//===================================================================================

	TEST_CASE("FieldTransform::Plane_wave_stays_plane_wave", "[electromagnetism][transform]")
	{
		// For a plane wave: |E| = |B| and E·B = 0
		// These properties are preserved under Lorentz transformation

		Real amplitude = 3.0;
		Vector3Cartesian E(0.0, amplitude, 0.0);  // E_y
		Vector3Cartesian B(0.0, 0.0, amplitude);  // B_z (perpendicular)

		// Verify it's a plane wave configuration
		REQUIRE_THAT(E.NormL2(), WithinAbs(B.NormL2(), 1e-12));
		REQUIRE_THAT(ComputeInvariant2(E, B), WithinAbs(0.0, 1e-12));

		Real v = 0.6;
		auto result = TransformFieldsBoostX(E, B, v);

		// E² - B² should still be 0
		REQUIRE_THAT(ComputeInvariant1(result.E_prime, result.B_prime), WithinAbs(0.0, 1e-10));

		// E·B should still be 0
		REQUIRE_THAT(ComputeInvariant2(result.E_prime, result.B_prime), WithinAbs(0.0, 1e-10));

		// |E'| should still equal |B'|
		REQUIRE_THAT(result.E_prime.NormL2(), 
		             WithinAbs(result.B_prime.NormL2(), 1e-10));
	}

	//===================================================================================
	// Inverse transformation
	//===================================================================================

	TEST_CASE("FieldTransform::Forward_inverse_roundtrip", "[electromagnetism][transform]")
	{
		Vector3Cartesian E(1.0, 2.0, 3.0);
		Vector3Cartesian B(0.5, 1.5, 2.5);
		Real v = 0.7;

		// Transform forward
		auto result1 = TransformFieldsBoostX(E, B, v);

		// Transform back (velocity sign flips)
		auto result2 = TransformFieldsBoostX(result1.E_prime, result1.B_prime, -v);

		// Should recover original fields
		REQUIRE_THAT(result2.E_prime.X(), WithinAbs(E.X(), 1e-10));
		REQUIRE_THAT(result2.E_prime.Y(), WithinAbs(E.Y(), 1e-10));
		REQUIRE_THAT(result2.E_prime.Z(), WithinAbs(E.Z(), 1e-10));

		REQUIRE_THAT(result2.B_prime.X(), WithinAbs(B.X(), 1e-10));
		REQUIRE_THAT(result2.B_prime.Y(), WithinAbs(B.Y(), 1e-10));
		REQUIRE_THAT(result2.B_prime.Z(), WithinAbs(B.Z(), 1e-10));
	}

	//===================================================================================
	// Zero velocity transformation is identity
	//===================================================================================

	TEST_CASE("FieldTransform::Zero_velocity_is_identity", "[electromagnetism][transform]")
	{
		Vector3Cartesian E(1.0, 2.0, 3.0);
		Vector3Cartesian B(0.5, 1.5, 2.5);

		auto result = TransformFieldsBoostX(E, B, 0.0);

		// Should be unchanged
		REQUIRE_THAT(result.E_prime.X(), WithinAbs(E.X(), 1e-12));
		REQUIRE_THAT(result.E_prime.Y(), WithinAbs(E.Y(), 1e-12));
		REQUIRE_THAT(result.E_prime.Z(), WithinAbs(E.Z(), 1e-12));

		REQUIRE_THAT(result.B_prime.X(), WithinAbs(B.X(), 1e-12));
		REQUIRE_THAT(result.B_prime.Y(), WithinAbs(B.Y(), 1e-12));
		REQUIRE_THAT(result.B_prime.Z(), WithinAbs(B.Z(), 1e-12));
	}

	//===================================================================================
	// Ultra-relativistic limit
	//===================================================================================

	TEST_CASE("FieldTransform::Ultra_relativistic_amplification", "[electromagnetism][transform]")
	{
		// At high velocities, transverse field components are amplified by γ

		Vector3Cartesian E(0.0, 1.0, 0.0);
		Vector3Cartesian B(0.0, 0.0, 0.0);
		Real v = 0.99;
		Real gamma = 1.0 / std::sqrt(1.0 - v * v);  // γ ≈ 7.09

		auto result = TransformFieldsBoostX(E, B, v);

		// E'_y ≈ γ * E_y (amplified)
		REQUIRE_THAT(result.E_prime.Y(), WithinAbs(gamma * E.Y(), 1e-10));

		// Significant B field generated
		REQUIRE(std::abs(result.B_prime.Z()) > 5.0);  // Should be about γ*v ≈ 7
	}

	//===================================================================================
	// Invalid velocity
	//===================================================================================

	TEST_CASE("FieldTransform::Invalid_velocity_throws", "[electromagnetism][transform]")
	{
		Vector3Cartesian E(1.0, 0.0, 0.0);
		Vector3Cartesian B(0.0, 0.0, 0.0);

		REQUIRE_THROWS_AS(TransformFieldsBoostX(E, B, 1.0), std::range_error);
		REQUIRE_THROWS_AS(TransformFieldsBoostX(E, B, -1.0), std::range_error);
		REQUIRE_THROWS_AS(TransformFieldsBoostX(E, B, 1.5), std::range_error);
	}

} // namespace MPL::Tests::Electromagnetism::FieldTransformationTests

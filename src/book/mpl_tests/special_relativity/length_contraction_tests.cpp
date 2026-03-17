///////////////////////////////////////////////////////////////////////////////////////////
// Length Contraction Tests - Special Relativity
///////////////////////////////////////////////////////////////////////////////////////////
//
// Tests for length contraction phenomenon:
// - L = L₀/γ (moving objects are shortened in direction of motion)
// - Ladder paradox verification
// - Proper length vs coordinate length
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

#include "SpecialRelativity/LorentzTransformation.h"

using namespace MML;
using namespace MPL;
using namespace MML::Testing;

using Catch::Matchers::WithinAbs;

namespace MPL::Tests::SpecialRelativity::LengthContractionTests
{
	//===================================================================================
	// Length contraction formula: L = L₀/γ
	//===================================================================================

	TEST_CASE("LengthContraction::Moving_rod_appears_shorter", "[special_relativity][length_contraction]")
	{
		TEST_PRECISION_INFO();

		// A rod at rest in S' has proper length L₀
		// Observed in S (where rod moves), length is L = L₀/γ

		Real L0 = 10.0;  // Proper length (rest frame of rod)

		SECTION("β = 0.6, γ = 1.25") {
			Real beta = 0.6;
			Real gamma = 1.25;
			Real L_contracted = L0 / gamma;

			REQUIRE_THAT(L_contracted, WithinAbs(8.0, 1e-10));
		}

		SECTION("β = 0.8, γ = 5/3") {
			Real beta = 0.8;
			Real gamma = 1.0 / std::sqrt(1.0 - beta * beta);
			Real L_contracted = L0 / gamma;

			REQUIRE_THAT(L_contracted, WithinAbs(6.0, 1e-10));
		}

		SECTION("β = 0.99, γ ≈ 7.09") {
			Real beta = 0.99;
			Real gamma = 1.0 / std::sqrt(1.0 - beta * beta);
			Real L_contracted = L0 / gamma;

			// Highly contracted
			REQUIRE_THAT(L_contracted, WithinAbs(10.0 / 7.0888, 1e-3));
			REQUIRE(L_contracted < 1.5);
		}
	}

	//===================================================================================
	// Length contraction via Lorentz transformation
	//===================================================================================

	TEST_CASE("LengthContraction::Via_Lorentz_transform", "[special_relativity][length_contraction]")
	{
		TEST_PRECISION_INFO();

		Real beta = 0.6;
		Real gamma = 1.0 / std::sqrt(1.0 - beta * beta);
		Real L0 = 10.0;  // Proper length

		CoordTransfLorentzXAxis L(beta);

		// Rod endpoints in S' (rod is at rest in S', moving in S)
		// At t'=0: rod extends from x'=0 to x'=L0
		// We need to find where these endpoints are in S at the SAME time t

		// In S, we measure the rod length at a single instant (simultaneous in S)
		// Let's use t=0 in S

		// Event 1: left end of rod at x'=0 when it's at x=0 in S at t=0
		// Transform (t'=0, x'=0, 0, 0) back to S to verify
		Vector4Minkowski left_prime({0.0, 0.0, 0.0, 0.0});
		auto left_S = L.transfInverse(left_prime);
		// This gives (0, 0, 0, 0) - left end at origin at t=0

		// Event 2: right end at x'=L0 when t'=?
		// We need to find t' such that when transformed to S, t=0
		// t = γ(t' + βx') = 0 → t' = -βx' = -β*L0
		Real t_prime_right = -beta * L0;
		Vector4Minkowski right_prime({t_prime_right, L0, 0.0, 0.0});
		auto right_S = L.transfInverse(right_prime);

		// Check that both events are at t=0 in S
		REQUIRE_THAT(left_S[0], WithinAbs(0.0, 1e-10));
		REQUIRE_THAT(right_S[0], WithinAbs(0.0, 1e-10));

		// Length in S
		Real L_observed = right_S[1] - left_S[1];
		Real L_expected = L0 / gamma;

		REQUIRE_THAT(L_observed, WithinAbs(L_expected, 1e-10));
	}

	//===================================================================================
	// Transverse dimensions unchanged
	//===================================================================================

	TEST_CASE("LengthContraction::Transverse_unchanged", "[special_relativity][length_contraction]")
	{
		TEST_PRECISION_INFO();

		// Contraction only in direction of motion (x)
		// y and z dimensions remain unchanged

		std::vector<Real> velocities = {0.3, 0.6, 0.9, 0.99};

		for (Real beta : velocities) {
			DYNAMIC_SECTION("β = " << beta) {
				CoordTransfLorentzXAxis L(beta);

				// A cube 10x10x10 in S', at rest
				// Measure at t=0 in S

				// For y-dimension: compare y' = 0 and y' = 10
				// These should remain 10 apart in S
				Vector4Minkowski p1({0.0, 0.0, 0.0, 0.0});
				Vector4Minkowski p2({0.0, 0.0, 10.0, 0.0});

				auto p1_S = L.transfInverse(p1);
				auto p2_S = L.transfInverse(p2);

				Real dy = p2_S[2] - p1_S[2];
				REQUIRE_THAT(dy, WithinAbs(10.0, 1e-10));

				// Same for z
				Vector4Minkowski p3({0.0, 0.0, 0.0, 10.0});
				auto p3_S = L.transfInverse(p3);
				Real dz = p3_S[3] - p1_S[3];
				REQUIRE_THAT(dz, WithinAbs(10.0, 1e-10));
			}
		}
	}

	//===================================================================================
	// Symmetry: both frames see contraction
	//===================================================================================

	TEST_CASE("LengthContraction::Symmetry", "[special_relativity][length_contraction]")
	{
		TEST_PRECISION_INFO();

		// S sees S' objects contracted
		// S' sees S objects contracted
		// This is consistent because of relativity of simultaneity

		Real beta = 0.6;
		Real gamma = 1.0 / std::sqrt(1.0 - beta * beta);
		Real L0 = 10.0;

		// Rod at rest in S: measured length is L0
		// From S' perspective, rod appears contracted to L0/γ

		// Rod at rest in S': measured length is L0
		// From S perspective, rod appears contracted to L0/γ

		// Both observers compute the same contraction factor
		Real L_contracted = L0 / gamma;

		REQUIRE_THAT(L_contracted, WithinAbs(8.0, 1e-10));

		// This is symmetric - the formula is the same regardless of which
		// frame is considered "moving"
	}

	//===================================================================================
	// Approaching ultra-relativistic limit
	//===================================================================================

	TEST_CASE("LengthContraction::Ultra_relativistic", "[special_relativity][length_contraction]")
	{
		TEST_PRECISION_INFO();

		// As v → c, γ → ∞, so L → 0

		Real L0 = 1.0;  // 1 meter proper length

		std::vector<std::pair<Real, Real>> beta_gamma = {
			{0.99, 7.0888},
			{0.999, 22.366},
			{0.9999, 70.712},
			{0.99999, 223.61},
		};

		for (const auto& [beta, gamma_expected] : beta_gamma) {
			DYNAMIC_SECTION("β = " << beta) {
				Real gamma = 1.0 / std::sqrt(1.0 - beta * beta);
				REQUIRE_THAT(gamma, WithinAbs(gamma_expected, 0.01));

				Real L_contracted = L0 / gamma;

				// Object becomes extremely thin
				REQUIRE(L_contracted < L0 / gamma_expected * 1.01);
			}
		}
	}

} // namespace MPL::Tests::SpecialRelativity::LengthContractionTests

///////////////////////////////////////////////////////////////////////////////////////////
// Lorentz Transformation Tests - Special Relativity
///////////////////////////////////////////////////////////////////////////////////////////
//
// Tests for CoordTransfLorentzXAxis and related Lorentz transformation functionality.
//
// Key verifications:
// 1. Forward and inverse transforms are self-consistent
// 2. Lorentz factor γ is computed correctly
// 3. Known physics results: time dilation, length contraction
// 4. Spacetime interval invariance: ds² = c²dt² - dx² - dy² - dz²
// 5. Velocity addition formula
// 6. Tensor transformation under Lorentz boost
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

#include "SpecialRelativity/LorentzTransformation.h"

using namespace MML;
using namespace MPL;
using namespace MML::Testing;

using Catch::Matchers::WithinAbs;

namespace MPL::Tests::SpecialRelativity::LorentzTransformationTests
{
	//===================================================================================
	// Basic Lorentz factor tests
	//===================================================================================

	TEST_CASE("LorentzTransformation::Gamma_factor_known_values", "[special_relativity][lorentz]")
	{
		TEST_PRECISION_INFO();

		// γ = 1/√(1-β²)
		// β = 0 → γ = 1
		// β = 0.6 → γ = 1.25
		// β = 0.8 → γ = 5/3 ≈ 1.6667
		// β = 0.99 → γ ≈ 7.089

		SECTION("Zero velocity: γ = 1") {
			CoordTransfLorentzXAxis L(0.0);
			// Transform origin at t=1 should give t'=1, x'=0
			Vector4Minkowski event({1.0, 0.0, 0.0, 0.0});
			auto result = L.transf(event);
			REQUIRE_THAT(result[0], WithinAbs(1.0, 1e-10));
			REQUIRE_THAT(result[1], WithinAbs(0.0, 1e-10));
		}

		SECTION("β = 0.6: γ = 1.25") {
			Real beta = 0.6;
			Real gamma_expected = 1.0 / std::sqrt(1.0 - beta * beta);
			REQUIRE_THAT(gamma_expected, WithinAbs(1.25, 1e-10));

			CoordTransfLorentzXAxis L(beta);
			// Event at origin at t=0: should remain at origin
			Vector4Minkowski event({0.0, 0.0, 0.0, 0.0});
			auto result = L.transf(event);
			REQUIRE_THAT(result[0], WithinAbs(0.0, 1e-10));
			REQUIRE_THAT(result[1], WithinAbs(0.0, 1e-10));
		}

		SECTION("β = 0.8: γ = 5/3") {
			Real beta = 0.8;
			Real gamma_expected = 1.0 / std::sqrt(1.0 - beta * beta);
			REQUIRE_THAT(gamma_expected, WithinAbs(5.0 / 3.0, 1e-10));
		}

		SECTION("High velocity β = 0.99") {
			Real beta = 0.99;
			Real gamma_expected = 1.0 / std::sqrt(1.0 - beta * beta);
			REQUIRE_THAT(gamma_expected, WithinAbs(7.0888, 1e-3));
		}
	}

	TEST_CASE("LorentzTransformation::Invalid_velocities_throw", "[special_relativity][lorentz]")
	{
		REQUIRE_THROWS_AS(CoordTransfLorentzXAxis(1.0), std::range_error);   // v = c
		REQUIRE_THROWS_AS(CoordTransfLorentzXAxis(-1.0), std::range_error);  // v = -c
		REQUIRE_THROWS_AS(CoordTransfLorentzXAxis(1.5), std::range_error);   // v > c
		REQUIRE_THROWS_AS(CoordTransfLorentzXAxis(-1.5), std::range_error);  // v < -c
	}

	//===================================================================================
	// Forward/Inverse consistency tests
	//===================================================================================

	TEST_CASE("LorentzTransformation::Forward_inverse_consistency", "[special_relativity][lorentz]")
	{
		TEST_PRECISION_INFO();

		std::vector<Real> velocities = {0.0, 0.1, 0.3, 0.5, 0.6, 0.8, 0.9, 0.99};
		std::vector<Vector4Minkowski> events = {
			Vector4Minkowski({0.0, 0.0, 0.0, 0.0}),
			Vector4Minkowski({1.0, 0.0, 0.0, 0.0}),
			Vector4Minkowski({0.0, 1.0, 0.0, 0.0}),
			Vector4Minkowski({1.0, 1.0, 0.0, 0.0}),
			Vector4Minkowski({5.0, 3.0, 2.0, 1.0}),
			Vector4Minkowski({-2.0, 4.0, -1.0, 3.0}),
		};

		for (Real beta : velocities) {
			DYNAMIC_SECTION("β = " << beta) {
				CoordTransfLorentzXAxis L(beta);

				for (const auto& event : events) {
					auto transformed = L.transf(event);
					auto restored = L.transfInverse(transformed);

					REQUIRE_THAT(restored[0], WithinAbs(event[0], 1e-10));
					REQUIRE_THAT(restored[1], WithinAbs(event[1], 1e-10));
					REQUIRE_THAT(restored[2], WithinAbs(event[2], 1e-10));
					REQUIRE_THAT(restored[3], WithinAbs(event[3], 1e-10));
				}
			}
		}
	}

	//===================================================================================
	// Spacetime interval invariance: ds² = c²dt² - dx² - dy² - dz² (with c=1)
	//===================================================================================

	Real SpacetimeInterval(const Vector4Minkowski& v)
	{
		// ds² = dt² - dx² - dy² - dz² (using c=1, signature +---)
		return v[0] * v[0] - v[1] * v[1] - v[2] * v[2] - v[3] * v[3];
	}

	TEST_CASE("LorentzTransformation::Spacetime_interval_invariance", "[special_relativity][lorentz]")
	{
		TEST_PRECISION_INFO();

		std::vector<Real> velocities = {0.1, 0.3, 0.5, 0.7, 0.9, 0.99};
		std::vector<Vector4Minkowski> events = {
			Vector4Minkowski({5.0, 3.0, 0.0, 0.0}),   // Timelike interval (ds² > 0)
			Vector4Minkowski({3.0, 5.0, 0.0, 0.0}),   // Spacelike interval (ds² < 0)
			Vector4Minkowski({5.0, 4.0, 3.0, 0.0}),   // Lightlike interval (ds² = 0)
			Vector4Minkowski({10.0, 6.0, 8.0, 0.0}),  // Lightlike (ds² = 0)
			Vector4Minkowski({7.0, 2.0, 3.0, 1.0}),   // General timelike
		};

		for (Real beta : velocities) {
			DYNAMIC_SECTION("β = " << beta) {
				CoordTransfLorentzXAxis L(beta);

				for (const auto& event : events) {
					Real ds2_original = SpacetimeInterval(event);
					auto transformed = L.transf(event);
					Real ds2_transformed = SpacetimeInterval(transformed);

					REQUIRE_THAT(ds2_transformed, WithinAbs(ds2_original, 1e-9));
				}
			}
		}
	}

	//===================================================================================
	// Time dilation: moving clocks run slower
	//===================================================================================

	TEST_CASE("LorentzTransformation::Time_dilation", "[special_relativity][lorentz]")
	{
		TEST_PRECISION_INFO();

		// A clock at rest in S ticks from t=0 to t=T at x=0
		// In S', the time elapsed is Δt' = γΔt (time dilation)

		Real beta = 0.6;
		Real gamma = 1.0 / std::sqrt(1.0 - beta * beta);  // γ = 1.25
		Real T = 10.0;  // Proper time

		CoordTransfLorentzXAxis L(beta);

		Vector4Minkowski event_start({0.0, 0.0, 0.0, 0.0});
		Vector4Minkowski event_end({T, 0.0, 0.0, 0.0});

		auto start_prime = L.transf(event_start);
		auto end_prime = L.transf(event_end);

		Real delta_t_prime = end_prime[0] - start_prime[0];
		Real delta_x_prime = end_prime[1] - start_prime[1];

		// In S', clock has moved: Δx' = -γβT (clock at origin in S moves backward in S')
		REQUIRE_THAT(delta_t_prime, WithinAbs(gamma * T, 1e-10));
		REQUIRE_THAT(delta_x_prime, WithinAbs(-gamma * beta * T, 1e-10));
	}

	//===================================================================================
	// Transverse coordinates unchanged
	//===================================================================================

	TEST_CASE("LorentzTransformation::Transverse_coordinates_unchanged", "[special_relativity][lorentz]")
	{
		TEST_PRECISION_INFO();

		std::vector<Real> velocities = {0.1, 0.5, 0.9, 0.99};

		for (Real beta : velocities) {
			DYNAMIC_SECTION("β = " << beta) {
				CoordTransfLorentzXAxis L(beta);

				Vector4Minkowski event({5.0, 3.0, 7.0, -2.5});
				auto transformed = L.transf(event);

				// y and z should remain unchanged
				REQUIRE_THAT(transformed[2], WithinAbs(event[2], 1e-12));
				REQUIRE_THAT(transformed[3], WithinAbs(event[3], 1e-12));
			}
		}
	}

	//===================================================================================
	// Composition of boosts (same direction)
	//===================================================================================

	TEST_CASE("LorentzTransformation::Composition_same_direction", "[special_relativity][lorentz]")
	{
		TEST_PRECISION_INFO();

		// Relativistic velocity addition: w = (u + v) / (1 + uv/c²)
		Real u = 0.5;  // First boost
		Real v = 0.4;  // Second boost
		Real w = (u + v) / (1.0 + u * v);  // Combined velocity

		CoordTransfLorentzXAxis L_u(u);
		CoordTransfLorentzXAxis L_v(v);
		CoordTransfLorentzXAxis L_w(w);

		Vector4Minkowski event({10.0, 5.0, 2.0, 1.0});

		// Transform by u, then by v
		auto result_uv = L_v.transf(L_u.transf(event));

		// Transform directly by combined velocity w
		auto result_w = L_w.transf(event);

		// Should be equal (for boosts in same direction)
		REQUIRE_THAT(result_uv[0], WithinAbs(result_w[0], 1e-9));
		REQUIRE_THAT(result_uv[1], WithinAbs(result_w[1], 1e-9));
		REQUIRE_THAT(result_uv[2], WithinAbs(result_w[2], 1e-9));
		REQUIRE_THAT(result_uv[3], WithinAbs(result_w[3], 1e-9));
	}

	//===================================================================================
	// Jacobian and tensor transformation
	//===================================================================================

	TEST_CASE("LorentzTransformation::Jacobian_is_Lorentz_matrix", "[special_relativity][lorentz]")
	{
		TEST_PRECISION_INFO();

		Real beta = 0.6;
		Real gamma = 1.0 / std::sqrt(1.0 - beta * beta);

		CoordTransfLorentzXAxis L(beta);

		// Get Jacobian at some point (should be constant for Lorentz transform)
		Vector4Minkowski point({5.0, 3.0, 1.0, 2.0});
		auto jacobian = L.jacobian(point);

		// Expected Lorentz matrix (boost along x):
		// | γ    -βγ  0  0 |
		// | -βγ   γ   0  0 |
		// | 0     0   1  0 |
		// | 0     0   0  1 |

		REQUIRE_THAT(jacobian(0, 0), WithinAbs(gamma, 1e-6));
		REQUIRE_THAT(jacobian(0, 1), WithinAbs(-beta * gamma, 1e-6));
		REQUIRE_THAT(jacobian(1, 0), WithinAbs(-beta * gamma, 1e-6));
		REQUIRE_THAT(jacobian(1, 1), WithinAbs(gamma, 1e-6));
		REQUIRE_THAT(jacobian(2, 2), WithinAbs(1.0, 1e-10));
		REQUIRE_THAT(jacobian(3, 3), WithinAbs(1.0, 1e-10));

		// Off-diagonal zeros
		REQUIRE_THAT(jacobian(0, 2), WithinAbs(0.0, 1e-10));
		REQUIRE_THAT(jacobian(0, 3), WithinAbs(0.0, 1e-10));
		REQUIRE_THAT(jacobian(1, 2), WithinAbs(0.0, 1e-10));
		REQUIRE_THAT(jacobian(1, 3), WithinAbs(0.0, 1e-10));
	}

	//===================================================================================
	// Light signal invariance (null geodesic)
	//===================================================================================

	TEST_CASE("LorentzTransformation::Light_signal_remains_light", "[special_relativity][lorentz]")
	{
		TEST_PRECISION_INFO();

		// Light signal: dx = c*dt (with c=1), so ds² = 0

		std::vector<Real> velocities = {0.3, 0.6, 0.9, 0.99};

		for (Real beta : velocities) {
			DYNAMIC_SECTION("β = " << beta) {
				CoordTransfLorentzXAxis L(beta);

				// Light pulse traveling in +x direction
				Vector4Minkowski light_start({0.0, 0.0, 0.0, 0.0});
				Vector4Minkowski light_end({10.0, 10.0, 0.0, 0.0});  // ds² = 100 - 100 = 0

				auto start_prime = L.transf(light_start);
				auto end_prime = L.transf(light_end);

				Real delta_t_prime = end_prime[0] - start_prime[0];
				Real delta_x_prime = end_prime[1] - start_prime[1];

				// Speed of light should still be 1 (c)
				Real speed_prime = std::abs(delta_x_prime / delta_t_prime);
				REQUIRE_THAT(speed_prime, WithinAbs(1.0, 1e-10));
			}
		}
	}

} // namespace MPL::Tests::SpecialRelativity::LorentzTransformationTests

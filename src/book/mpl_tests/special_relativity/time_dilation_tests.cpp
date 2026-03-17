///////////////////////////////////////////////////////////////////////////////////////////
// Time Dilation Tests - Special Relativity
///////////////////////////////////////////////////////////////////////////////////////////
//
// Detailed tests for time dilation phenomenon:
// - Δt' = γΔτ (proper time to coordinate time)
// - Twin paradox scenarios
// - Muon lifetime extension
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

namespace MPL::Tests::SpecialRelativity::TimeDilationTests
{
	//===================================================================================
	// Time dilation formula: Δt = γΔτ
	//===================================================================================

	TEST_CASE("TimeDilation::Proper_time_to_coordinate_time", "[special_relativity][time_dilation]")
	{
		TEST_PRECISION_INFO();

		// Proper time τ is measured by clock at rest with the object
		// Coordinate time t is measured by "stationary" observer
		// Δt = γΔτ

		SECTION("β = 0.6, γ = 1.25") {
			Real beta = 0.6;
			Real gamma = 1.25;
			Real proper_time = 100.0;  // seconds in moving frame
			Real coordinate_time = gamma * proper_time;  // seconds in rest frame

			REQUIRE_THAT(coordinate_time, WithinAbs(125.0, 1e-10));
		}

		SECTION("β = 0.8, γ = 5/3") {
			Real beta = 0.8;
			Real gamma = 1.0 / std::sqrt(1.0 - beta * beta);
			Real proper_time = 60.0;  // 1 minute moving frame
			Real coordinate_time = gamma * proper_time;

			REQUIRE_THAT(coordinate_time, WithinAbs(100.0, 1e-10));
		}

		SECTION("β = 0.99, γ ≈ 7.09") {
			Real beta = 0.99;
			Real gamma = 1.0 / std::sqrt(1.0 - beta * beta);
			Real proper_time = 1.0;  // 1 year moving frame
			Real coordinate_time = gamma * proper_time;

			REQUIRE_THAT(coordinate_time, WithinAbs(7.0888, 1e-3));
		}
	}

	//===================================================================================
	// Muon lifetime extension
	//===================================================================================

	TEST_CASE("TimeDilation::Muon_lifetime_extension", "[special_relativity][time_dilation]")
	{
		TEST_PRECISION_INFO();

		// Muon rest lifetime: τ ≈ 2.2 μs
		// At v = 0.995c, γ ≈ 10
		// Extended lifetime in Earth frame: t ≈ 22 μs

		Real tau_rest = 2.2e-6;  // seconds (proper lifetime)
		Real beta = 0.995;
		Real gamma = 1.0 / std::sqrt(1.0 - beta * beta);

		Real lifetime_earth = gamma * tau_rest;

		// At this speed, gamma ≈ 10, so lifetime ≈ 22 μs
		REQUIRE(gamma > 9.9);
		REQUIRE(gamma < 10.1);
		REQUIRE_THAT(lifetime_earth, WithinAbs(22e-6, 1e-6));

		// Distance traveled: d = v * t_earth
		Real c = 3e8;  // m/s
		Real v = beta * c;
		Real distance_traveled = v * lifetime_earth;

		// Should be able to travel much farther than expected from rest lifetime
		Real distance_rest_would_predict = v * tau_rest;
		REQUIRE(distance_traveled > 9.5 * distance_rest_would_predict);
	}

	//===================================================================================
	// Symmetry of time dilation
	//===================================================================================

	TEST_CASE("TimeDilation::Symmetry_between_frames", "[special_relativity][time_dilation]")
	{
		TEST_PRECISION_INFO();

		// Key insight: Time dilation is symmetric
		// S sees S' clocks running slow by factor γ
		// S' sees S clocks running slow by factor γ
		// This is NOT a contradiction - it's about comparing clocks at same location

		Real beta = 0.6;
		Real gamma = 1.0 / std::sqrt(1.0 - beta * beta);

		CoordTransfLorentzXAxis L_to_prime(beta);
		CoordTransfLorentzXAxis L_to_rest(-beta);  // Inverse boost

		// Clock in S: ticks at x=0 from t=0 to t=10
		Vector4Minkowski S_start({0.0, 0.0, 0.0, 0.0});
		Vector4Minkowski S_end({10.0, 0.0, 0.0, 0.0});

		// In S' frame
		auto S_start_in_prime = L_to_prime.transf(S_start);
		auto S_end_in_prime = L_to_prime.transf(S_end);

		Real dt_prime = S_end_in_prime[0] - S_start_in_prime[0];

		// S' sees S clock running slow: dt' = γ * dt_proper
		// But the S clock has moved in S', so we need to account for that
		// Proper time in S is 10, coordinate time in S' is γ*10 = 12.5
		REQUIRE_THAT(dt_prime, WithinAbs(gamma * 10.0, 1e-10));
	}

	//===================================================================================
	// GPS satellite time correction (simplified)
	//===================================================================================

	TEST_CASE("TimeDilation::GPS_satellite_approximation", "[special_relativity][time_dilation]")
	{
		TEST_PRECISION_INFO();

		// GPS satellite velocity: v ≈ 3.87 km/s ≈ 1.29e-5 c
		// This is very non-relativistic, so γ ≈ 1 + β²/2

		Real v_satellite = 3870.0;  // m/s
		Real c = 299792458.0;       // m/s
		Real beta = v_satellite / c;

		// Exact gamma
		Real gamma_exact = 1.0 / std::sqrt(1.0 - beta * beta);

		// Low-velocity approximation: γ ≈ 1 + β²/2
		Real gamma_approx = 1.0 + 0.5 * beta * beta;

		// Should be very close (β is tiny)
		REQUIRE_THAT(gamma_exact, WithinAbs(gamma_approx, 1e-15));

		// Time dilation per day
		Real seconds_per_day = 86400.0;
		Real time_dilation_per_day = (gamma_exact - 1.0) * seconds_per_day;

		// Should be about 7 microseconds per day (special relativity contribution)
		// (Note: General relativity gives opposite and larger effect of ~45 μs/day)
		REQUIRE(time_dilation_per_day > 6e-6);
		REQUIRE(time_dilation_per_day < 8e-6);
	}

} // namespace MPL::Tests::SpecialRelativity::TimeDilationTests

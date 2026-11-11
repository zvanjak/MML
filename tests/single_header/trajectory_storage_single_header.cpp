#include <MML.h>
#include <catch2/catch_test_macros.hpp>

TEST_CASE("Single header exposes bounded trajectory storage", "[single_header][trajectory-storage]")
{
	MML::ODESystemSolution ode(0, 1, 1, 0);
	ode.fillValues(0, 0, MML::Vector<MML::Real>{3});
	ode.fillValues(99999, 1, MML::Vector<MML::Real>{4});
	REQUIRE(ode.size() == 100000);
	REQUIRE_NOTHROW(ode.fillValues(100000, 2, MML::Vector<MML::Real>{5}));
	REQUIRE(ode.size() == 100001);
	REQUIRE_THROWS_AS(ode.reserveForPoint(ode.max_size()), MML::ResourceLimitError);
	REQUIRE(ode.getXValue(99999, 0) == MML::Real(4));
	REQUIRE(MML::to_string(MML::AlgorithmStatus::ResourceLimitExceeded) == "ResourceLimitExceeded");
}

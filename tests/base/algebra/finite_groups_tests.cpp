#include <catch2/catch_test_macros.hpp>

#include <mml/base/Algebra_base.h>

#include <vector>

namespace MML::Tests::Base::AlgebraTests
{
	TEST_CASE("CyclicGroup exposes rotations and polygon permutations", "[Algebra][CyclicGroup][Base]")
	{
		Algebra::CyclicGroup group(5);

		REQUIRE(group.order() == 5);
		REQUIRE(group.generator() == group.rotation(1));
		REQUIRE(group.rotation(-1) == group.rotation(4));
		REQUIRE(group.compose(group.rotation(2), group.rotation(4)) == group.rotation(1));
		REQUIRE(group.inverse(group.rotation(2)) == group.rotation(3));
		REQUIRE(group.permutation(group.rotation(2)).images() == std::vector<int>{2, 3, 4, 0, 1});
		REQUIRE_THROWS_AS(Algebra::CyclicGroup(0), ArgumentError);
	}

	TEST_CASE("DihedralGroup exposes rotations reflections and polygon actions", "[Algebra][DihedralGroup][Base]")
	{
		Algebra::DihedralGroup group(4);
		const auto rotation = group.rotation(1);
		const auto reflection = group.reflection();

		REQUIRE(group.order() == 8);
		REQUIRE(group.inverse(rotation) == group.rotation(3));
		REQUIRE(group.inverse(group.reflection(2)) == group.reflection(2));
		REQUIRE(group.compose(group.compose(reflection, rotation), reflection) == group.inverse(rotation));
		REQUIRE(group.permutation(rotation).images() == std::vector<int>{1, 2, 3, 0});
		REQUIRE(group.permutation(reflection).images() == std::vector<int>{0, 3, 2, 1});
		REQUIRE_THROWS_AS(Algebra::DihedralGroup(1), ArgumentError);
	}
}
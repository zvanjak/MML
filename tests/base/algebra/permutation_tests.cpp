#include <catch2/catch_test_macros.hpp>

#include <mml/base/Algebra_base.h>

#include <array>
#include <vector>

namespace MML::Tests::Base::AlgebraTests
{
	using Algebra::DynamicPermutation;
	using Algebra::Permutation;
	using Algebra::PermutationParity;

	TEST_CASE("Permutation validates images", "[Algebra][Permutation][Base]")
	{
		REQUIRE_THROWS_AS(Permutation<3>({0, 0, 2}), ArgumentError);
		REQUIRE_THROWS_AS(Permutation<3>({0, 1, 3}), ArgumentError);
		REQUIRE_THROWS_AS(DynamicPermutation({0, -1, 2}), ArgumentError);
		REQUIRE_THROWS_AS(Permutation<4>::from_cycles({{0, 1}, {1, 2}}), ArgumentError);
	}

	TEST_CASE("Permutation identity inverse and composition follow function order", "[Algebra][Permutation][Base]")
	{
		const auto cycle = Permutation<3>::from_cycles({{0, 1, 2}});
		const auto swap = Permutation<3>::from_cycles({{0, 1}});
		const auto identity = Permutation<3>::identity();

		REQUIRE(cycle.compose(identity) == cycle);
		REQUIRE(identity.compose(cycle) == cycle);
		REQUIRE(cycle.compose(cycle.inverse()) == identity);
		REQUIRE(cycle.inverse().compose(cycle) == identity);
		REQUIRE(cycle.compose(swap).apply(0) == cycle.apply(swap.apply(0)));
		REQUIRE(cycle.compose(swap).compose(identity) == cycle.compose(swap.compose(identity)));
	}

	TEST_CASE("Permutation reports cycles order and parity", "[Algebra][Permutation][Base]")
	{
		const auto permutation = Permutation<6>::from_cycles({{0, 2, 4}, {1, 3}});
		const auto decomposition = permutation.cycles();

		REQUIRE(decomposition == std::vector<std::vector<int>>{{0, 2, 4}, {1, 3}});
		REQUIRE(permutation.order() == 6);
		REQUIRE(permutation.transposition_count() == 3);
		REQUIRE(permutation.parity() == PermutationParity::Odd);
		REQUIRE(permutation.sign() == -1);
		REQUIRE(Permutation<6>::identity().order() == 1);
		REQUIRE(Permutation<6>::identity().parity() == PermutationParity::Even);
	}

	TEST_CASE("Permutation applies to polygon vertex positions", "[Algebra][Permutation][Base]")
	{
		const auto quarterTurn = Permutation<4>::from_cycles({{0, 1, 2, 3}});
		const std::array<char, 4> vertices{'A', 'B', 'C', 'D'};

		REQUIRE(quarterTurn.apply(0) == 1);
		REQUIRE(quarterTurn.apply(vertices) == std::array<char, 4>{'D', 'A', 'B', 'C'});
	}

	TEST_CASE("DynamicPermutation mirrors fixed-size behavior", "[Algebra][Permutation][Base]")
	{
		const auto permutation = DynamicPermutation::from_cycles(5, {{0, 3, 1}, {2, 4}});
		const auto identity = DynamicPermutation::identity(5);

		REQUIRE(permutation.compose(permutation.inverse()) == identity);
		REQUIRE(permutation.cycles() == std::vector<std::vector<int>>{{0, 3, 1}, {2, 4}});
		REQUIRE(permutation.order() == 6);
		REQUIRE(permutation.apply(std::vector<int>{10, 20, 30, 40, 50}) ==
			std::vector<int>{20, 40, 50, 10, 30});
		REQUIRE_THROWS_AS(permutation.compose(DynamicPermutation::identity(4)), ArgumentError);
	}

	TEST_CASE("FiniteGroup stores exact structure for Core validation", "[Algebra][FiniteGroup][Base]")
	{
		Algebra::FiniteGroup<int> group(
			{0, 1, 2}, 0,
			[](const int& left, const int& right) { return (left + right) % 3; },
			[](const int& element) { return (3 - element) % 3; });

		REQUIRE(group.order() == 3);
		REQUIRE(group.contains(2));
		REQUIRE_FALSE(group.contains(3));
		REQUIRE(group.compose(1, 2) == 0);
		REQUIRE(group.inverse(2) == 1);
	}
}
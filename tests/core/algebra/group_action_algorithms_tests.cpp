#include <catch2/catch_test_macros.hpp>

#include <mml/base/Algebra_base.h>

#include <cstddef>
#include <vector>

namespace MML::Tests::Core::AlgebraTests
{
	std::vector<std::vector<int>> BinaryColorings(int size)
	{
		std::vector<std::vector<int>> colorings;
		const int count = 1 << size;
		for (int mask = 0; mask < count; mask++) {
			std::vector<int> coloring(static_cast<std::size_t>(size));
			for (int vertex = 0; vertex < size; vertex++)
				coloring[static_cast<std::size_t>(vertex)] = (mask >> vertex) & 1;
			colorings.push_back(std::move(coloring));
		}
		return colorings;
	}

	TEST_CASE("Polygon action satisfies orbit stabilizer under D5", "[Algebra][GroupAction][Core]")
	{
		Algebra::DihedralGroup group(5);
		auto action = Algebra::MakeGroupAction<Algebra::DihedralElement, int>(
			[&group](const Algebra::DihedralElement& element, const int& vertex) {
				return group.permutation(element).apply(vertex);
			});

		const auto orbit = Algebra::Orbit(group, action, 0);
		const auto stabilizer = Algebra::Stabilizer(group, action, 0);

		REQUIRE(orbit == std::vector<int>{0, 1, 2, 3, 4});
		REQUIRE(stabilizer.size() == 2);
		REQUIRE(orbit.size() * stabilizer.size() == group.order());
		REQUIRE_FALSE(Algebra::IsInvariant(group, action, 0));
	}

	TEST_CASE("Group action invariance supports object relations and predicates", "[Algebra][GroupAction][Core]")
	{
		Algebra::CyclicGroup group(4);
		auto coloringAction = Algebra::MakeGroupAction<Algebra::CyclicElement, std::vector<int>>(
			[&group](const Algebra::CyclicElement& element, const std::vector<int>& coloring) {
				return group.permutation(element).apply(coloring);
			});
		auto vertexAction = Algebra::MakeGroupAction<Algebra::CyclicElement, int>(
			[&group](const Algebra::CyclicElement& element, const int& vertex) {
				return group.permutation(element).apply(vertex);
			});

		REQUIRE(Algebra::IsInvariant(group, coloringAction, std::vector<int>{1, 1, 1, 1}));
		REQUIRE_FALSE(Algebra::IsInvariant(group, coloringAction, std::vector<int>{1, 0, 0, 0}));
		REQUIRE(Algebra::IsPredicateInvariant(group, vertexAction, std::vector<int>{0, 1, 2, 3},
			[](int vertex) { return vertex >= 0 && vertex < 4; }));
		REQUIRE_FALSE(Algebra::IsPredicateInvariant(group, vertexAction, std::vector<int>{0, 1, 2, 3},
			[](int vertex) { return vertex == 0; }));
	}

	TEST_CASE("Burnside counts binary necklaces under C4", "[Algebra][Burnside][Core]")
	{
		Algebra::CyclicGroup group(4);
		const auto colorings = BinaryColorings(4);
		auto action = Algebra::MakeGroupAction<Algebra::CyclicElement, std::vector<int>>(
			[&group](const Algebra::CyclicElement& element, const std::vector<int>& coloring) {
				return group.permutation(element).apply(coloring);
			});

		REQUIRE(Algebra::FixedPointCount(group, action, group.identity(), colorings) == 16);
		REQUIRE(Algebra::FixedPointCount(group, action, group.rotation(1), colorings) == 2);
		REQUIRE(Algebra::BurnsideCount(group, action, colorings) == 6);
		REQUIRE(Algebra::OrbitPartition(group, action, colorings).size() == 6);
	}

	TEST_CASE("Burnside counts binary bracelets under D4", "[Algebra][Burnside][Core]")
	{
		Algebra::DihedralGroup group(4);
		const auto colorings = BinaryColorings(4);
		auto action = Algebra::MakeGroupAction<Algebra::DihedralElement, std::vector<int>>(
			[&group](const Algebra::DihedralElement& element, const std::vector<int>& coloring) {
				return group.permutation(element).apply(coloring);
			});

		REQUIRE(Algebra::BurnsideCount(group, action, colorings) == 6);
		REQUIRE(Algebra::OrbitPartition(group, action, colorings).size() == 6);
	}

	TEST_CASE("Burnside rejects a non-closed finite object set", "[Algebra][Burnside][Core]")
	{
		Algebra::CyclicGroup group(3);
		auto action = Algebra::MakeGroupAction<Algebra::CyclicElement, int>(
			[&group](const Algebra::CyclicElement& element, const int& vertex) {
				return group.permutation(element).apply(vertex);
			});

		REQUIRE_THROWS_AS(Algebra::BurnsideCount(group, action, std::vector<int>{0, 1}), DomainError);
	}
}
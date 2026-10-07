#include <catch2/catch_test_macros.hpp>

#include <mml/base/Algebra_base.h>

#include <algorithm>
#include <vector>

namespace MML::Tests::Core::AlgebraTests
{
	TEST_CASE("CyclicGroup generator powers and algorithms cover C6", "[Algebra][CyclicGroup][Core]")
	{
		Algebra::CyclicGroup group(6);
		std::vector<Algebra::CyclicElement> powers;
		auto power = group.identity();
		for (std::size_t exponent = 0; exponent < group.order(); exponent++) {
			powers.push_back(power);
			power = group.compose(power, group.generator());
		}

		REQUIRE(powers == group.elements());
		REQUIRE(Algebra::element_order(group, group.generator()) == 6);
		REQUIRE(Algebra::element_order(group, group.rotation(2)) == 3);
		REQUIRE(Algebra::GeneratedSubgroup(group, {group.rotation(2)}).size() == 3);
		REQUIRE(Algebra::ConjugacyClasses(group).size() == 6);
		REQUIRE(Algebra::CheckGroupLaws(group));
	}

	TEST_CASE("DihedralGroup satisfies relations and finite group algorithms", "[Algebra][DihedralGroup][Core]")
	{
		Algebra::DihedralGroup group(4);
		const auto rotation = group.generator();
		const auto reflection = group.reflection();

		REQUIRE(Algebra::CheckGroupLaws(group));
		REQUIRE(Algebra::element_order(group, rotation) == 4);
		REQUIRE(Algebra::element_order(group, reflection) == 2);
		REQUIRE(Algebra::GeneratedSubgroup(group, {rotation}).size() == 4);
		REQUIRE(Algebra::GeneratedSubgroup(group, {rotation, reflection}).size() == 8);

		const auto classes = Algebra::ConjugacyClasses(group);
		std::vector<std::size_t> classSizes;
		for (const auto& conjugacyClass : classes)
			classSizes.push_back(conjugacyClass.size());
		std::sort(classSizes.begin(), classSizes.end());
		REQUIRE(classSizes == std::vector<std::size_t>{1, 1, 2, 2, 2});
	}

	TEST_CASE("Cayley table stores closed element indices", "[Algebra][FiniteGroup][Core]")
	{
		Algebra::CyclicGroup group(4);
		const auto table = Algebra::CayleyTable(group);

		REQUIRE(table.size() == 4);
		REQUIRE(table[1] == std::vector<std::size_t>{1, 2, 3, 0});
		for (const auto& row : table)
			for (std::size_t index : row)
				REQUIRE(index < group.order());
	}

	TEST_CASE("Cayley graph data exposes generator-labelled edges", "[Algebra][FiniteGroup][Core]")
	{
		Algebra::DihedralGroup group(3);
		const auto data = Algebra::MakeCayleyGraphData(group, {group.generator(), group.reflection()});

		REQUIRE(data.elements.size() == 6);
		REQUIRE(data.generators.size() == 2);
		REQUIRE(data.edges.size() == 12);
		REQUIRE(data.edges[0].from == 0);
		REQUIRE(data.edges[0].to == 1);
		REQUIRE(data.edges[0].generator == 0);
		REQUIRE(data.edges[1].to == 3);
		REQUIRE(data.edges[1].generator == 1);
	}
}
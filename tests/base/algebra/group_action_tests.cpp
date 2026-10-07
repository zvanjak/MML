#include <catch2/catch_test_macros.hpp>

#include <mml/base/Algebra_base.h>

namespace MML::Tests::Base::AlgebraTests
{
	TEST_CASE("GroupAction wraps a typed left action", "[Algebra][GroupAction][Base]")
	{
		Algebra::CyclicGroup group(5);
		auto action = Algebra::MakeGroupAction<Algebra::CyclicElement, int>(
			[&group](const Algebra::CyclicElement& element, const int& vertex) {
				return group.permutation(element).apply(vertex);
			});

		REQUIRE(action.apply(group.identity(), 2) == 2);
		REQUIRE(action(group.rotation(2), 4) == 1);
		REQUIRE(action(group.compose(group.rotation(2), group.rotation(3)), 4) ==
			action(group.rotation(2), action(group.rotation(3), 4)));
		REQUIRE_THROWS_AS(
			(Algebra::GroupAction<Algebra::CyclicElement, int>({})), ArgumentError);
	}
}
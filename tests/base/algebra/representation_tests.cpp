#include <catch2/catch_test_macros.hpp>

#include <mml/base/Algebra_base.h>

namespace MML::Tests::Base::AlgebraTests
{
	TEST_CASE("Representation wraps matrix evaluation and vector application", "[Algebra][Representation][Base]")
	{
		Algebra::CyclicGroup group(2);
		auto representation = Algebra::MakeRepresentation<Algebra::CyclicElement, Real, 2>(
			[](const Algebra::CyclicElement& element) {
				return element.exponent == 0
					? MatrixNM<Real, 2, 2>::Identity()
					: MatrixNM<Real, 2, 2>{{-1, 0}, {0, -1}};
			});

		REQUIRE(representation.matrix(group.identity()) == MatrixNM<Real, 2, 2>::Identity());
		REQUIRE(representation.apply(group.generator(), VectorN<Real, 2>{2, -3}) ==
			VectorN<Real, 2>{-2, 3});
		REQUIRE_THROWS_AS(
			(Algebra::Representation<Algebra::CyclicElement, Real, 2>({})), ArgumentError);
	}
}
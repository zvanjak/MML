#include <catch2/catch_test_macros.hpp>

#include <mml/base/Algebra_base.h>

namespace MML::Tests::Core::AlgebraTests
{
	TEST_CASE("Modular structures satisfy ring and field laws", "[Algebra][PrimeField][Core]")
	{
		Algebra::ModularRing<6> ring;
		Algebra::PrimeField<5> field;

		REQUIRE(Algebra::CheckRingLaws(ring));
		REQUIRE(Algebra::CheckRingLaws(field));
		REQUIRE(Algebra::CheckFieldLaws(field));
	}

	TEST_CASE("Prime field helpers find primitive roots and quadratic residues", "[Algebra][PrimeField][Core]")
	{
		using F7 = Algebra::PrimeFieldElement<7>;
		const auto generator = Algebra::PrimitiveRoot<7>();

		REQUIRE(generator == F7(3));
		REQUIRE(generator.pow(6) == F7(1));
		REQUIRE(Algebra::LegendreSymbol<7>(F7(0)) == 0);
		REQUIRE(Algebra::LegendreSymbol<7>(F7(2)) == 1);
		REQUIRE(Algebra::LegendreSymbol<7>(F7(3)) == -1);
		REQUIRE(Algebra::IsQuadraticResidue<7>(F7(4)));
		REQUIRE_FALSE(Algebra::IsQuadraticResidue<7>(F7(5)));
	}

	template<int Prime>
	void CheckMatrixInverse(const MatrixNM<Algebra::PrimeFieldElement<Prime>, 2, 2>& matrix)
	{
		using Field = Algebra::PrimeFieldElement<Prime>;
		const auto inverse = Algebra::FieldMatrixInverse(matrix);
		const auto product = matrix * inverse;

		REQUIRE(product(0, 0) == Field(1));
		REQUIRE(product(0, 1) == Field(0));
		REQUIRE(product(1, 0) == Field(0));
		REQUIRE(product(1, 1) == Field(1));
		REQUIRE(Algebra::FieldMatrixDeterminant(matrix) != Field(0));
	}

	TEST_CASE("Exact field matrices invert over GF2 GF3 and GF5", "[Algebra][PrimeField][Matrix][Core]")
	{
		CheckMatrixInverse<2>({{1, 1}, {1, 0}});
		CheckMatrixInverse<3>({{1, 2}, {2, 2}});
		CheckMatrixInverse<5>({{1, 2}, {3, 4}});

		using F5 = Algebra::PrimeFieldElement<5>;
		REQUIRE_THROWS_AS((Algebra::FieldMatrixInverse(MatrixNM<F5, 2, 2>{{1, 2}, {2, 4}})), DomainError);
	}
}
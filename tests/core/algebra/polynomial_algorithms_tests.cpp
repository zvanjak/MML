#include <catch2/catch_test_macros.hpp>

#include <mml/base/Algebra_base.h>

namespace MML::Tests::Core::AlgebraTests
{
	using F2 = Algebra::PrimeFieldElement<2>;
	using F3 = Algebra::PrimeFieldElement<3>;

	struct CoreGF4Modulus
	{
		static Algebra::Polynomial<F2> modulus() { return {1, 1, 1}; }
	};

	struct CoreGF8Modulus
	{
		static Algebra::Polynomial<F2> modulus() { return {1, 1, 0, 1}; }
	};

	struct CoreGF9Modulus
	{
		static Algebra::Polynomial<F3> modulus() { return {1, 0, 1}; }
	};

	TEST_CASE("Polynomial GCD and extended GCD satisfy Bezout identity", "[Algebra][Polynomial][Core]")
	{
		using Polynomial = Algebra::Polynomial<F2>;
		const Polynomial left{1, 0, 1, 1};
		const Polynomial right{1, 1, 1};
		const auto result = Algebra::PolynomialExtendedGcd(left, right);

		REQUIRE(Algebra::PolynomialGcd(left, right) == result.gcd);
		REQUIRE(result.left_coefficient * left + result.right_coefficient * right == result.gcd);
		REQUIRE(result.gcd == Polynomial::one());
	}

	TEST_CASE("Polynomial irreducibility recognizes extension moduli", "[Algebra][Polynomial][Core]")
	{
		REQUIRE(Algebra::IsIrreducible(CoreGF4Modulus::modulus()));
		REQUIRE(Algebra::IsIrreducible(CoreGF8Modulus::modulus()));
		REQUIRE(Algebra::IsIrreducible(CoreGF9Modulus::modulus()));
		REQUIRE_FALSE(Algebra::IsIrreducible(Algebra::Polynomial<F2>{1, 0, 1}));
	}

	TEST_CASE("Small extension fields satisfy exhaustive field laws", "[Algebra][FiniteField][Core]")
	{
		Algebra::ExtensionField<2, 2, CoreGF4Modulus> gf4;
		Algebra::ExtensionField<2, 3, CoreGF8Modulus> gf8;
		Algebra::ExtensionField<3, 2, CoreGF9Modulus> gf9;

		REQUIRE(Algebra::CheckFieldLaws(gf4));
		REQUIRE(Algebra::CheckFieldLaws(gf8));
		REQUIRE(Algebra::CheckFieldLaws(gf9));
		for (const auto& element : gf9.elements())
			if (element != gf9.zero())
				REQUIRE(gf9.multiply(element, gf9.inverse(element)) == gf9.one());
	}
}
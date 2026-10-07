#include <catch2/catch_test_macros.hpp>

#include <mml/base/Algebra_base.h>

namespace MML::Tests::Base::AlgebraTests
{
	using F2 = Algebra::PrimeFieldElement<2>;
	using F3 = Algebra::PrimeFieldElement<3>;

	struct GF4Modulus
	{
		static Algebra::Polynomial<F2> modulus() { return {1, 1, 1}; }
	};

	struct GF8Modulus
	{
		static Algebra::Polynomial<F2> modulus() { return {1, 1, 0, 1}; }
	};

	struct GF9Modulus
	{
		static Algebra::Polynomial<F3> modulus() { return {1, 0, 1}; }
	};

	TEST_CASE("Exact Polynomial canonicalizes and divides over fields", "[Algebra][Polynomial][Base]")
	{
		using Polynomial = Algebra::Polynomial<Algebra::PrimeFieldElement<5>>;
		const Polynomial dividend{4, 0, 1};
		const Polynomial divisor{1, 1};
		const auto division = dividend.divmod(divisor);

		REQUIRE(Polynomial({1, 2, 0, 0}) == Polynomial({1, 2}));
		REQUIRE(dividend == division.first * divisor + division.second);
		REQUIRE(division.second.degree() < divisor.degree());
		REQUIRE(dividend.evaluate({2}) == Algebra::PrimeFieldElement<5>(3));
		REQUIRE_THROWS_AS(dividend.divmod(Polynomial::zero()), DivisionByZeroError);
	}

	TEST_CASE("GF4 arithmetic follows x squared equals x plus one", "[Algebra][FiniteField][Base]")
	{
		using GF4 = Algebra::FiniteFieldElement<2, 2, GF4Modulus>;
		const GF4 alpha{0, 1};
		const std::array<GF4, 4> elements{GF4(0), GF4(1), alpha, alpha + GF4(1)};
		const int expectedProducts[4][4] = {
			{0, 0, 0, 0},
			{0, 1, 2, 3},
			{0, 2, 3, 1},
			{0, 3, 1, 2}
		};

		REQUIRE(alpha * alpha == alpha + GF4(1));
		REQUIRE(alpha.pow(3) == GF4(1));
		REQUIRE(alpha.inverse() == alpha + GF4(1));
		for (int row = 0; row < 4; row++)
			for (int column = 0; column < 4; column++)
				REQUIRE(elements[static_cast<std::size_t>(row)] * elements[static_cast<std::size_t>(column)] ==
					elements[static_cast<std::size_t>(expectedProducts[row][column])]);
		REQUIRE_THROWS_AS(GF4(0).inverse(), DivisionByZeroError);
	}

	TEST_CASE("GF8 and GF9 reduce by documented irreducible polynomials", "[Algebra][FiniteField][Base]")
	{
		using GF8 = Algebra::FiniteFieldElement<2, 3, GF8Modulus>;
		using GF9 = Algebra::FiniteFieldElement<3, 2, GF9Modulus>;
		const GF8 alpha8{0, 1, 0};
		const GF9 alpha9{0, 1};

		REQUIRE(alpha8.pow(3) == alpha8 + GF8(1));
		REQUIRE(alpha8.pow(7) == GF8(1));
		REQUIRE(alpha9 * alpha9 == GF9(2));
		REQUIRE(alpha9.pow(4) == GF9(1));
		REQUIRE(Algebra::ExtensionField<2, 3, GF8Modulus>().order() == 8);
		REQUIRE(Algebra::ExtensionField<3, 2, GF9Modulus>().order() == 9);
	}
}
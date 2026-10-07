#include <catch2/catch_test_macros.hpp>

#include <mml/base/Algebra_base.h>

namespace MML::Tests::Base::AlgebraTests
{
	TEST_CASE("ModInt normalizes and computes exact ring arithmetic", "[Algebra][ModInt][Base]")
	{
		using Z12 = Algebra::ModInt<12>;

		static_assert(!Z12::is_field);
		static_assert(Algebra::IsPrime<13>::value);
		static_assert(!Algebra::IsPrime<12>::value);
		REQUIRE(Z12(-1).value() == 11);
		REQUIRE(Z12(10) + Z12(5) == Z12(3));
		REQUIRE(Z12(2) - Z12(5) == Z12(9));
		REQUIRE(Z12(7) * Z12(8) == Z12(8));
		REQUIRE(Z12(5).pow(3) == Z12(5));
		REQUIRE(Z12(5).inverse() == Z12(5));
		REQUIRE_THROWS_AS(Z12(6).inverse(), DomainError);
		REQUIRE_THROWS_AS(Z12(0).inverse(), DivisionByZeroError);
	}

	TEST_CASE("PrimeFieldElement supports exact field division", "[Algebra][PrimeField][Base]")
	{
		using F7 = Algebra::PrimeFieldElement<7>;

		static_assert(F7::is_field);
		REQUIRE(F7(3).inverse() == F7(5));
		REQUIRE(F7(3) / F7(2) == F7(5));
		REQUIRE(F7(6).pow(2) == F7(1));
		for (int value = 1; value < 7; value++)
			REQUIRE(F7(value) * F7(value).inverse() == F7(1));
	}

	TEST_CASE("Modular structure adapters expose finite ring and field domains", "[Algebra][PrimeField][Base]")
	{
		Algebra::ModularRing<8> ring;
		Algebra::PrimeField<5> field;

		REQUIRE(ring.order() == 8);
		REQUIRE(ring.multiply({2}, {4}) == Algebra::ModInt<8>(0));
		REQUIRE(field.order() == 5);
		REQUIRE(field.inverse({3}) == Algebra::PrimeFieldElement<5>(2));
	}
}
///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML) Tests                            ///
///                                                                                   ///
///  File:        number_theory_tests.cpp                                             ///
///  Description: Tests for NumberTheory.h (gcd, modular arithmetic, primes, CRT)     ///
///////////////////////////////////////////////////////////////////////////////////////////

#include "../../mml/base/NumberTheory.h"

#include <catch2/catch_test_macros.hpp>

#include <cstdint>
#include <vector>

using namespace MML;
using namespace MML::NumberTheory;

namespace MML::Tests::Base::NumberTheoryTests {

TEST_CASE("NumberTheory_gcd_lcm", "[numbertheory]") {
	REQUIRE(Gcd(48, 36) == 12);
	REQUIRE(Gcd(-48, 36) == 12);
	REQUIRE(Gcd(17, 5) == 1);
	REQUIRE(Gcd(0, 0) == 0);
	REQUIRE(Gcd(0, 7) == 7);

	REQUIRE(Lcm(4, 6) == 12);
	REQUIRE(Lcm(21, 6) == 42);
	REQUIRE(Lcm(0, 5) == 0);
}

TEST_CASE("NumberTheory_extended_gcd_identity", "[numbertheory]") {
	long long x, y;
	long long g = ExtendedGcd(240, 46, x, y);
	REQUIRE(g == 2);
	REQUIRE(240 * x + 46 * y == g);

	g = ExtendedGcd(1000003, 17, x, y);
	REQUIRE(g == 1);
	REQUIRE(1000003 * x + 17 * y == 1);
}

TEST_CASE("NumberTheory_modpow", "[numbertheory]") {
	REQUIRE(ModPow(2, 10, 1000) == 24);
	REQUIRE(ModPow(3, 0, 7) == 1);
	REQUIRE(ModPow(0, 5, 7) == 0);
	// Fermat's little theorem: a^(p-1) == 1 (mod p) for prime p not dividing a
	const std::uint64_t p = 1000000007ull;
	REQUIRE(ModPow(2, p - 1, p) == 1);
	REQUIRE(ModPow(123456789, p - 1, p) == 1);
}

TEST_CASE("NumberTheory_mod_inverse", "[numbertheory]") {
	REQUIRE(ModInverse(3, 11) == 4);        // 3*4 = 12 == 1 (mod 11)
	REQUIRE((3 * ModInverse(3, 11)) % 11 == 1);

	const long long p = 1000000007ll;
	long long inv = ModInverse(2, p);
	REQUIRE((2 * inv) % p == 1);

	REQUIRE_THROWS_AS(ModInverse(2, 4), DomainError); // gcd(2,4) != 1
}

TEST_CASE("NumberTheory_is_prime", "[numbertheory]") {
	REQUIRE(IsPrime(2));
	REQUIRE(IsPrime(3));
	REQUIRE(IsPrime(97));
	REQUIRE(IsPrime(7919));
	REQUIRE(IsPrime(1000000007ull));
	REQUIRE(IsPrime(1000000009ull));

	REQUIRE_FALSE(IsPrime(0));
	REQUIRE_FALSE(IsPrime(1));
	REQUIRE_FALSE(IsPrime(4));
	REQUIRE_FALSE(IsPrime(561));   // Carmichael number
	REQUIRE_FALSE(IsPrime(1105));  // Carmichael number

	// large semiprime must be recognized as composite
	REQUIRE_FALSE(IsPrime(1000000007ull * 1000000009ull));
}

TEST_CASE("NumberTheory_factorize", "[numbertheory]") {
	auto f360 = Factorize(360);
	REQUIRE(f360 == std::vector<std::uint64_t>{ 2, 2, 2, 3, 3, 5 });

	auto pw = FactorizePowers(360);
	REQUIRE(pw.size() == 3);
	REQUIRE(pw[0] == std::pair<std::uint64_t, int>{ 2, 3 });
	REQUIRE(pw[1] == std::pair<std::uint64_t, int>{ 3, 2 });
	REQUIRE(pw[2] == std::pair<std::uint64_t, int>{ 5, 1 });

	// Pollard rho must split a large semiprime and recover the factors
	const std::uint64_t a = 1000000007ull, b = 1000000009ull;
	auto sp = Factorize(a * b);
	REQUIRE(sp == std::vector<std::uint64_t>{ a, b });

	REQUIRE(Factorize(1).empty());
	REQUIRE(Factorize(13) == std::vector<std::uint64_t>{ 13 });
}

TEST_CASE("NumberTheory_multiplicative_functions", "[numbertheory]") {
	REQUIRE(EulerTotient(1) == 1);
	REQUIRE(EulerTotient(10) == 4);
	REQUIRE(EulerTotient(36) == 12);
	REQUIRE(EulerTotient(1000000007) == 1000000006); // prime

	REQUIRE(MobiusMu(1) == 1);
	REQUIRE(MobiusMu(2) == -1);
	REQUIRE(MobiusMu(6) == 1);
	REQUIRE(MobiusMu(12) == 0);   // 4 | 12
	REQUIRE(MobiusMu(30) == -1);

	REQUIRE(DivisorCount(1) == 1);
	REQUIRE(DivisorCount(28) == 6);
	REQUIRE(DivisorSum(6) == 12);
	REQUIRE(DivisorSum(28) == 56); // perfect number
}

TEST_CASE("NumberTheory_prime_enumeration", "[numbertheory]") {
	auto primes = PrimesUpTo(30);
	REQUIRE(primes == std::vector<std::uint64_t>{ 2, 3, 5, 7, 11, 13, 17, 19, 23, 29 });
	REQUIRE(PrimesUpTo(1).empty());

	REQUIRE(NextPrime(13) == 17);
	REQUIRE(NextPrime(1000000) == 1000003);
}

TEST_CASE("NumberTheory_chinese_remainder", "[numbertheory]") {
	auto r1 = ChineseRemainder({ 2, 3, 2 }, { 3, 5, 7 });
	REQUIRE(r1.first == 23);
	REQUIRE(r1.second == 105);

	// non-coprime but consistent moduli
	auto r2 = ChineseRemainder({ 1, 4 }, { 6, 9 });
	REQUIRE(r2.second == 18);
	REQUIRE(r2.first % 6 == 1);
	REQUIRE(r2.first % 9 == 4);

	// inconsistent congruences must throw
	REQUIRE_THROWS_AS(ChineseRemainder({ 0, 1 }, { 2, 4 }), DomainError);
}

} // namespace MML::Tests::Base::NumberTheoryTests

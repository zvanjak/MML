#include "../../mml/base/Combinatorics.h"

#include <catch2/catch_test_macros.hpp>
#include <climits>
#include <cstdint>
#include <limits>
#include <vector>

namespace MML::Tests::Base::ExactIntegerRegression {

TEST_CASE("Exact integers safe baseline findings", "[exact-integer][exact-safe-red]") {
	SECTION("B01 zero-residue witnesses") {
		REQUIRE(NumberTheory::IsPrime(73));
		REQUIRE(NumberTheory::IsPrime(193));
		REQUIRE(NumberTheory::NextPrime(71) == 73);
	}
	SECTION("B03 next prime exhaustion") {
		REQUIRE_THROWS_AS(NumberTheory::NextPrime(UINT64_MAX), DomainError);
	}
	SECTION("B10 minimum numerator identity") {
		REQUIRE(Rational(LLONG_MIN) + Rational(0) == Rational(LLONG_MIN));
	}
	SECTION("B10 portable exact ordering") {
		REQUIRE(Rational(LLONG_MAX, 2) > Rational(LLONG_MAX - 2, 3));
	}
	SECTION("B11 Catalan representable result") {
		REQUIRE(Combinatorics::Catalan(35) == 3116285494907301262LL);
	}
}

TEST_CASE("Exact primality independent sieve and witness factors", "[exact-integer][numbertheory]") {
	constexpr std::size_t limit = 100000;
	std::vector<bool> prime(limit + 1, true);
	prime[0] = prime[1] = false;
	for (std::size_t p = 2; p <= limit / p; ++p)
		if (prime[p])
			for (std::size_t m = p * p; m <= limit; m += p) prime[m] = false;
	for (std::size_t n = 0; n <= limit; ++n) {
		INFO("candidate " << n);
		REQUIRE(NumberTheory::IsPrime(n) == prime[n]);
	}
	// Every prime divisor of the seven deterministic Miller-Rabin witness constants.
	for (std::uint64_t p : {2ULL, 3ULL, 5ULL, 13ULL, 19ULL, 73ULL, 193ULL,
							 407521ULL, 299210837ULL}) {
		REQUIRE(NumberTheory::IsPrime(p));
		REQUIRE(NumberTheory::Factorize(p) == std::vector<std::uint64_t>{p});
	}
	REQUIRE(NumberTheory::Factorize(0).empty());
	REQUIRE(NumberTheory::Factorize(1).empty());
	REQUIRE(NumberTheory::FactorizePowers(0).empty());
	constexpr std::uint64_t lastPrime = 18446744073709551557ULL;
	REQUIRE(NumberTheory::IsPrime(lastPrime));
	REQUIRE_FALSE(NumberTheory::IsPrime(UINT64_MAX));
	REQUIRE_FALSE(NumberTheory::IsPrime(341550071728321ULL));
	REQUIRE_FALSE(NumberTheory::IsPrime(3825123056546413051ULL));
	REQUIRE(NumberTheory::NextPrime(lastPrime - 1) == lastPrime);
	REQUIRE_THROWS_AS(NumberTheory::NextPrime(lastPrime), DomainError);
	REQUIRE_THROWS_AS(NumberTheory::NextPrime(UINT64_MAX - 1), DomainError);
	REQUIRE(NumberTheory::NextPrime(0) == 2);
	REQUIRE(NumberTheory::Factorize(UINT64_MAX) ==
		std::vector<std::uint64_t>{3, 5, 17, 257, 641, 65537, 6700417});
	// Above 2^63: rho's square-plus-c must use modular, not wrapping, addition.
	REQUIRE(NumberTheory::Factorize(18446743979220271189ULL) ==
		std::vector<std::uint64_t>{4294967279ULL, 4294967291ULL});
}

TEST_CASE("Exact signed and modular full range boundaries", "[exact-integer][numbertheory]") {
	using namespace NumberTheory;
	REQUIRE(Gcd(LLONG_MIN, 6) == 2);
	REQUIRE(Gcd(LLONG_MIN, LLONG_MAX) == 1);
	REQUIRE_THROWS_AS(Gcd(LLONG_MIN, 0), DomainError);
	REQUIRE_THROWS_AS(Gcd(LLONG_MIN, LLONG_MIN), DomainError);
	REQUIRE(Lcm(LLONG_MIN, 0) == 0);
	REQUIRE(Lcm(-LLONG_MAX, LLONG_MAX) == LLONG_MAX);
	REQUIRE_THROWS_AS(Lcm(LLONG_MIN, 1), DomainError);
	REQUIRE_THROWS_AS(Lcm(LLONG_MAX, 2), DomainError);
	long long x = 0, y = 0;
	REQUIRE(ExtendedGcd(LLONG_MIN, -1, x, y) == 1);
	REQUIRE(x == 0);
	REQUIRE(y == -1);
	REQUIRE(ExtendedGcd(-240, 46, x, y) == 2);
	REQUIRE(-240 * x + 46 * y == 2);
	REQUIRE_THROWS_AS(ExtendedGcd(LLONG_MIN, 0, x, y), DomainError);
	REQUIRE_THROWS_AS(MulMod(1, 2, 0), DomainError);
	REQUIRE_THROWS_AS(ModPow(2, 0, 0), DomainError);
	REQUIRE(MulMod(UINT64_MAX, UINT64_MAX, UINT64_MAX) == 0);
	REQUIRE(MulMod(UINT64_MAX - 1, UINT64_MAX - 1, UINT64_MAX) == 1);
	REQUIRE(MulMod(UINT64_MAX - 2, UINT64_MAX - 3, UINT64_MAX) == 6);
	// Independently computed with System.Numerics.BigInteger, not this implementation.
	REQUIRE(MulMod(UINT64_MAX, UINT64_MAX - 1, 18446744073709551557ULL) == 3306);
	REQUIRE(MulMod(1ULL << 63, 1ULL << 63, UINT64_MAX) == 4611686018427387904ULL);
	REQUIRE(MulMod(12345678901234567890ULL, 9876543210987654321ULL, UINT64_MAX) ==
		6743105841750238095ULL);
	REQUIRE(MulMod(UINT64_MAX, UINT64_MAX, 9223372036854775809ULL) == 9);
	REQUIRE(ModPow(UINT64_MAX - 1, 2, UINT64_MAX) == 1);
	REQUIRE(ModPow(0, 0, 1) == 0);
	REQUIRE(ModInverse(LLONG_MIN, LLONG_MAX) == LLONG_MAX - 1);
	REQUIRE(ModInverse(2, LLONG_MAX) == LLONG_MAX / 2 + 1);
	REQUIRE(ModInverse(LLONG_MAX - 1, LLONG_MAX) == LLONG_MAX - 1);
	REQUIRE_THROWS_AS(ModInverse(1, 0), DomainError);
	REQUIRE_THROWS_AS(ModInverse(1, -1), DomainError);
	REQUIRE(detail::MulModS(LLONG_MIN, LLONG_MIN, LLONG_MAX) == 1);
	REQUIRE_THROWS_AS(detail::MulModS(1, 2, 0), DomainError);
	REQUIRE(ChineseRemainder({LLONG_MIN}, {LLONG_MAX}) ==
		std::pair<long long, long long>{LLONG_MAX - 1, LLONG_MAX});
	REQUIRE_THROWS_AS(ChineseRemainder({LLONG_MAX - 1, 1}, {LLONG_MAX, 2}), DomainError);
}

TEST_CASE("Exact CRT prechecked growth and canonical sums", "[exact-integer][numbertheory]") {
	using namespace NumberTheory;
	REQUIRE_THROWS_AS(ChineseRemainder({0, 0}, {LLONG_MAX, 2}), DomainError);
	REQUIRE_THROWS_AS(ChineseRemainder({0}, {0}), DomainError);
	REQUIRE_THROWS_AS(ChineseRemainder({0}, {-1}), DomainError);
	REQUIRE_THROWS_AS(ChineseRemainder({0, 1}, {2, 4}), DomainError);
	// Already-canonical residues near the signed limit must not overflow normalization.
	REQUIRE(ChineseRemainder({LLONG_MAX - 2, 0}, {LLONG_MAX, 1}) ==
		std::pair<long long, long long>{LLONG_MAX - 2, LLONG_MAX});
	REQUIRE(ChineseRemainder({-1, -1}, {2147483647LL, 2147483629LL}) ==
		std::pair<long long, long long>{4611685975477714962LL, 4611685975477714963LL});
	REQUIRE(DivisorSum(1LL << 61) == (1LL << 62) - 1);
	REQUIRE_THROWS_AS(DivisorSum(3LL << 61), DomainError);
	REQUIRE_THROWS_AS(PrimesUpTo(UINT64_MAX), DomainError);
}

TEST_CASE("Exact shared preoperation signed arithmetic boundaries", "[exact-integer]") {
	using namespace ExactInteger;
	REQUIRE(Add(LLONG_MIN, 0) == LLONG_MIN);
	REQUIRE(Add(LLONG_MIN, LLONG_MAX) == -1);
	REQUIRE(Add(LLONG_MAX - 1, 1) == LLONG_MAX);
	REQUIRE_THROWS_AS(Add(LLONG_MAX, 1), DomainError);
	REQUIRE_THROWS_AS(Add(LLONG_MIN, -1), DomainError);
	REQUIRE(Subtract(LLONG_MIN, LLONG_MIN) == 0);
	REQUIRE(Subtract(-1, LLONG_MIN) == LLONG_MAX);
	REQUIRE_THROWS_AS(Subtract(0, LLONG_MIN), DomainError);
	REQUIRE_THROWS_AS(Subtract(LLONG_MIN, 1), DomainError);
	REQUIRE(Multiply(LLONG_MIN, 1) == LLONG_MIN);
	REQUIRE(Multiply(LLONG_MIN, 0) == 0);
	REQUIRE(Multiply(-LLONG_MAX, -1) == LLONG_MAX);
	REQUIRE(Multiply(LLONG_MIN / 2, 2) == LLONG_MIN);
	REQUIRE_THROWS_AS(Multiply(LLONG_MIN, -1), DomainError);
	REQUIRE_THROWS_AS(Multiply(LLONG_MIN, 2), DomainError);
	REQUIRE_THROWS_AS(Multiply(LLONG_MAX, -2), DomainError);
	REQUIRE_THROWS_AS(Multiply(-LLONG_MAX, -2), DomainError);
}

TEST_CASE("Exact Stirling entries match an independent small triangle", "[exact-integer][combinatorics]") {
	std::vector<long long> first{1}, second{1};
	for (int n = 0; n <= 20; ++n) {
		for (int k = 0; k <= n + 1; ++k) {
			INFO("n=" << n << " k=" << k);
			REQUIRE(Combinatorics::StirlingFirst(n, k) == (k <= n ? first[k] : 0));
			REQUIRE(Combinatorics::StirlingSecond(n, k) == (k <= n ? second[k] : 0));
		}
		if (n == 20) break;
		std::vector<long long> f(n + 2, 0), s(n + 2, 0);
		for (int k = 1; k <= n + 1; ++k) {
			f[k] = first[k - 1] + (k <= n ? n * first[k] : 0);
			s[k] = second[k - 1] + (k <= n ? k * second[k] : 0);
		}
		first = std::move(f); second = std::move(s);
	}
}

TEST_CASE("Exact combinatorics last valid first invalid sequences", "[exact-integer][combinatorics]") {
	using namespace Combinatorics;
	REQUIRE(Fibonacci(92) == 7540113804746346429LL);
	REQUIRE_THROWS_AS(Fibonacci(93), DomainError);
	REQUIRE(Lucas(90) == 6440026026380244498LL);
	REQUIRE_THROWS_AS(Lucas(91), DomainError);
	REQUIRE(Catalan(35) == 3116285494907301262LL);
	REQUIRE_THROWS_AS(Catalan(36), DomainError);
	REQUIRE(BellNumber(25) == 4638590332229999353LL);
	REQUIRE_THROWS_AS(BellNumber(26), DomainError);
	REQUIRE(StirlingFirst(21, 1) == 2432902008176640000LL);
	REQUIRE_THROWS_AS(StirlingFirst(22, 1), DomainError);
	REQUIRE(StirlingFirst(21, 2) == 8752948036761600000LL);
	REQUIRE_THROWS_AS(StirlingFirst(22, 2), DomainError);
	REQUIRE(StirlingSecond(64, 2) == LLONG_MAX);
	REQUIRE_THROWS_AS(StirlingSecond(65, 2), DomainError);
	REQUIRE(StirlingFirst(INT_MAX, 0) == 0);
	REQUIRE(StirlingSecond(INT_MAX, 0) == 0);
	REQUIRE_THROWS_AS(StirlingFirst(INT_MAX, INT_MAX / 2), DomainError);
	REQUIRE_THROWS_AS(StirlingSecond(INT_MAX, INT_MAX / 2), DomainError);
	REQUIRE(StirlingFirst(100, 99) == 4950);
	REQUIRE(StirlingSecond(100, 99) == 4950);
	REQUIRE(StirlingFirst(INT_MAX, INT_MAX) == 1);
	REQUIRE(StirlingSecond(INT_MAX, INT_MAX) == 1);
	REQUIRE(StirlingSecond(INT_MAX, 1) == 1);
	REQUIRE(StirlingFirst(INT_MAX, INT_MAX - 1) == 2305843005992468481LL);
	REQUIRE(StirlingSecond(INT_MAX, INT_MAX - 1) == 2305843005992468481LL);
	REQUIRE(ZigzagNumber(23) == 1015423886506852352LL);
	REQUIRE_THROWS_AS(ZigzagNumber(24), DomainError);
	REQUIRE(Tangent(12) == 1015423886506852352LL);
	REQUIRE_THROWS_AS(Tangent(13), DomainError);
	REQUIRE(PartitionCount(405) == 9147679068859117602LL);
	REQUIRE_THROWS_AS(PartitionCount(406), DomainError);
	for (int n : {INT_MAX, 1000000}) {
		REQUIRE_THROWS_AS(Fibonacci(n), DomainError);
		REQUIRE_THROWS_AS(Lucas(n), DomainError);
		REQUIRE_THROWS_AS(Catalan(n), DomainError);
		REQUIRE_THROWS_AS(BellNumber(n), DomainError);
		REQUIRE_THROWS_AS(ZigzagNumber(n), DomainError);
		REQUIRE_THROWS_AS(Tangent(n), DomainError);
		REQUIRE_THROWS_AS(PartitionCount(n), DomainError);
	}
}

TEST_CASE("Exact combinatorics checked products indices and independent partition oracle",
		  "[exact-integer][combinatorics]") {
	using namespace Combinatorics;
	REQUIRE(FallingFactorial(LLONG_MIN, 1) == LLONG_MIN);
	REQUIRE(RisingFactorial(LLONG_MIN, 1) == LLONG_MIN);
	REQUIRE_THROWS_AS(FallingFactorial(LLONG_MIN, 2), DomainError);
	REQUIRE_THROWS_AS(RisingFactorial(LLONG_MAX, 2), DomainError);
	REQUIRE_THROWS_AS(FallingFactorial(LLONG_MAX, 2), DomainError);
	REQUIRE_THROWS_AS(RisingFactorial(LLONG_MIN, 2), DomainError);
	REQUIRE(FallingFactorial(21, 22) == 0);
	REQUIRE(RisingFactorial(-21, 22) == 0);
	REQUIRE(Multinomial({INT_MAX, 1}) == 2147483648LL);
	REQUIRE(Multinomial({INT_MAX, 0}) == 1);
	REQUIRE_THROWS_AS(Multinomial({INT_MAX, INT_MAX}), DomainError);
	REQUIRE(BinomialCoefficient(66, 33) == 7219428434016265740LL);
	REQUIRE_THROWS_AS(BinomialCoefficient(67, 33), DomainError);
	REQUIRE_THROWS_AS(detail::MulCheckedLL(LLONG_MIN, -1), DomainError);
	REQUIRE_THROWS_AS(detail::MulCheckedLL(LLONG_MAX, 2), DomainError);
	REQUIRE(detail::MulCheckedLL(LLONG_MIN, 1) == LLONG_MIN);
	// Independent generating-function coin-change oracle over the small domain.
	std::vector<long long> p(101, 0);
	p[0] = 1;
	for (int part = 1; part <= 100; ++part)
		for (int n = part; n <= 100; ++n) p[n] += p[n - part];
	for (int n = 0; n <= 100; ++n) REQUIRE(PartitionCount(n) == p[n]);
	REQUIRE(PartitionCount(-1) == 0);
	REQUIRE_THROWS_AS(BellNumber(-1), DomainError);
	REQUIRE_THROWS_AS(StirlingFirst(-1, 0), DomainError);
	REQUIRE_THROWS_AS(StirlingSecond(0, -1), DomainError);
	REQUIRE(Bernoulli(INT_MAX) == Rational(0));
	REQUIRE_THROWS_AS(Bernoulli(INT_MAX - 1), DomainError);
}

TEST_CASE("Exact rational minimum domain identities reduced arithmetic and ordering",
		  "[exact-integer][rational]") {
	const Rational min(LLONG_MIN);
	REQUIRE(Rational(LLONG_MIN, 1) == min);
	REQUIRE(Rational(LLONG_MIN, 3).Num() == LLONG_MIN);
	REQUIRE(Rational(LLONG_MIN, LLONG_MIN) == Rational(1));
	REQUIRE(Rational(0, LLONG_MIN) == Rational(0));
	REQUIRE(Rational(2, LLONG_MIN) == Rational(-1, 1LL << 62));
	REQUIRE(min + 0 == min);
	REQUIRE(min - 0 == min);
	REQUIRE(min - min == Rational(0));
	REQUIRE(min * 1 == min);
	REQUIRE(min * 0 == Rational(0));
	REQUIRE(min / 1 == min);
	REQUIRE(min / min == Rational(1));
	REQUIRE(Rational(0) - min / 2 == Rational(1LL << 62));
	REQUIRE_THROWS_AS(-min, DomainError);
	REQUIRE_THROWS_AS(min / -1, DomainError);
	REQUIRE_THROWS_AS(min - 1, DomainError);
	REQUIRE_THROWS_AS(Rational(LLONG_MIN, -1), DomainError);
	REQUIRE_THROWS_AS(Rational(1, LLONG_MIN), DomainError);
	REQUIRE_THROWS_AS(Rational(1) / min, DomainError);
	REQUIRE_THROWS_AS(min / 0, DivisionByZeroError);
	REQUIRE(Rational(LLONG_MAX, 2) + Rational(LLONG_MAX, 2) == Rational(LLONG_MAX));
	REQUIRE(Rational(LLONG_MAX / 2, 2) - Rational(LLONG_MAX - 2, 3) ==
		Rational(-4611686018427387901LL, 6));
	REQUIRE_THROWS_AS(Rational(LLONG_MAX, 2) - Rational(LLONG_MAX - 2, 3), DomainError);
	REQUIRE(Rational(1, LLONG_MAX) + Rational(LLONG_MAX - 1, LLONG_MAX) == Rational(1));
	REQUIRE(Rational(1, LLONG_MAX) / Rational(2, LLONG_MAX) == Rational(1, 2));
	REQUIRE(Rational(LLONG_MIN, 3) * Rational(3, 2) == Rational(LLONG_MIN / 2));
	const std::vector<Rational> ordered{
		min, Rational(-LLONG_MAX, 2), Rational(LLONG_MIN, 3),
		Rational(-1, LLONG_MAX), Rational(0), Rational(1, LLONG_MAX),
		Rational(LLONG_MAX - 2, 3), Rational(LLONG_MAX, 2), Rational(LLONG_MAX)};
	for (std::size_t i = 0; i < ordered.size(); ++i) {
		const auto& a = ordered[i];
		REQUIRE_FALSE(a < a);
		REQUIRE(a <= a);
		for (std::size_t j = 0; j < ordered.size(); ++j) {
			REQUIRE((a < ordered[j]) == (i < j));
			REQUIRE((a < ordered[j]) == (ordered[j] > a));
		}
	}
	REQUIRE(Rational(-LLONG_MAX, 2) < Rational(LLONG_MIN, 3));
	REQUIRE(min < Rational(LLONG_MIN + 1));
	REQUIRE(Rational(LLONG_MAX, 2) > Rational(LLONG_MAX - 2, 3));
	REQUIRE(Rational(1, LLONG_MAX) < Rational(1, LLONG_MAX - 1));
	Rational unchanged(LLONG_MAX);
	REQUIRE_THROWS_AS(unchanged += Rational(1), DomainError);
	REQUIRE(unchanged == Rational(LLONG_MAX));
}

#if defined(__SIZEOF_INT128__)
TEST_CASE("Exact portable two limb arithmetic agrees with native int128", "[exact-integer][int128]") {
	using U128 = unsigned __int128;
	std::uint64_t state = 0x7349a2b8c117d053ULL;
	auto next = [&]() {
		state ^= state << 13; state ^= state >> 7; state ^= state << 17;
		return state;
	};
	for (int i = 0; i < 2000; ++i) {
		const auto a = next(), b = next(), d = next() | std::uint64_t{1};
		const auto expected = static_cast<U128>(a) * b;
		const auto actual = ExactInteger::WideUnsigned::Product(a, b);
		REQUIRE(actual.high == static_cast<std::uint64_t>(expected >> 64));
		REQUIRE(actual.low == static_cast<std::uint64_t>(expected));
		std::uint64_t remainder;
		const auto quotient = actual.Divide(d, remainder);
		REQUIRE(quotient.high == static_cast<std::uint64_t>((expected / d) >> 64));
		REQUIRE(quotient.low == static_cast<std::uint64_t>(expected / d));
		REQUIRE(remainder == static_cast<std::uint64_t>(expected % d));
		REQUIRE(NumberTheory::MulMod(a, b, d) == static_cast<std::uint64_t>(expected % d));
		const long long sa = static_cast<long long>(a), sb = static_cast<long long>(b);
		long long x, y;
		const auto g = NumberTheory::ExtendedGcd(sa, sb, x, y);
		REQUIRE(static_cast<__int128>(sa) * x + static_cast<__int128>(sb) * y == g);
		const auto da = static_cast<long long>((next() >> 1) | 1ULL);
		const auto db = static_cast<long long>((next() >> 1) | 1ULL);
		REQUIRE((Rational(sa, da) < Rational(sb, db)) ==
			(static_cast<__int128>(sa) * db < static_cast<__int128>(sb) * da));
	}
}

TEST_CASE("Exact rational arithmetic agrees with independent native wide oracle", "[exact-integer][int128]") {
	using I128 = __int128;
	std::uint64_t state = 0xb451cd739f286a01ULL;
	auto next = [&]() {
		state ^= state << 13; state ^= state >> 7; state ^= state << 17;
		return state;
	};
	auto verify = [](I128 n, I128 d, const auto& operation) {
		if (d < 0) { n = -n; d = -d; }
		I128 a = n < 0 ? -n : n, b = d;
		while (b != 0) { const auto r = a % b; a = b; b = r; }
		n /= a; d /= a;
		if (n < LLONG_MIN || n > LLONG_MAX || d > LLONG_MAX) {
			REQUIRE_THROWS_AS(operation(), DomainError);
		}
		else {
			const auto actual = operation();
			REQUIRE(actual.Num() == static_cast<long long>(n));
			REQUIRE(actual.Den() == static_cast<long long>(d));
		}
	};
	for (int i = 0; i < 1500; ++i) {
		const auto a = static_cast<long long>(next());
		const auto c = static_cast<long long>(next() | 1ULL);
		const auto b = static_cast<long long>(i % 2 ? (next() >> 1) | 1ULL : next() % 100 + 1);
		const auto d = static_cast<long long>(i % 2 ? b : next() % 100 + 1);
		const Rational left(a, b), right(c, d);
		verify(static_cast<I128>(a) * d + static_cast<I128>(c) * b,
			static_cast<I128>(b) * d, [&]() { return left + right; });
		verify(static_cast<I128>(a) * d - static_cast<I128>(c) * b,
			static_cast<I128>(b) * d, [&]() { return left - right; });
		verify(static_cast<I128>(a) * c, static_cast<I128>(b) * d,
			[&]() { return left * right; });
		verify(static_cast<I128>(a) * d, static_cast<I128>(b) * c,
			[&]() { return left / right; });
	}
}

TEST_CASE("Exact CRT full range residues satisfy independently checked congruences",
		  "[exact-integer][int128]") {
	std::uint64_t state = 0xf238914ac5076bd1ULL;
	auto next = [&]() {
		state ^= state << 13; state ^= state >> 7; state ^= state << 17;
		return state;
	};
	for (int i = 0; i < 500; ++i) {
		const auto m = static_cast<long long>(next() % 3000000000ULL + 1);
		const auto n = static_cast<long long>(next() % 3000000000ULL + 1);
		const auto a = static_cast<long long>(next()), b = static_cast<long long>(next());
		long long g = m, h = n;
		while (h != 0) { const auto r = g % h; g = h; h = r; }
		const auto ra = static_cast<long long>((static_cast<__int128>(a) % m + m) % m);
		const auto rb = static_cast<long long>((static_cast<__int128>(b) % n + n) % n);
		if ((ra - rb) % g != 0) {
			REQUIRE_THROWS_AS(NumberTheory::ChineseRemainder({a, b}, {m, n}), DomainError);
		}
		else {
			const auto result = NumberTheory::ChineseRemainder({a, b}, {m, n});
			REQUIRE(result.second == static_cast<long long>(static_cast<__int128>(m / g) * n));
			REQUIRE(result.first >= 0);
			REQUIRE(result.first < result.second);
			REQUIRE(result.first % m == ra);
			REQUIRE(result.first % n == rb);
		}
	}
}
#endif

} // namespace MML::Tests::Base::ExactIntegerRegression

#include "../TestPrecision.h"
#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/base/StandardFunctions.h>
#include <mml/base/DiracDeltaFunction.h>
#endif
#include <catch2/catch_all.hpp>
#include <cfenv>

using namespace MML;
using Catch::Matchers::WithinRel;

namespace {
	int LastFiniteFactorialArgument() {
		const Real logMax = std::log(std::numeric_limits<Real>::max());
		for (int n = 2; n <= std::numeric_limits<Real>::max_exponent; ++n) {
			if (std::lgamma(Real(n + 1)) > logMax)
				return n - 1;
		}
		throw NumericalMethodError("Factorial test could not locate the Real range boundary");
	}

	class FloatingEnvironmentGuard {
		std::fenv_t saved;
	public:
		FloatingEnvironmentGuard() { REQUIRE(std::feholdexcept(&saved) == 0); }
		~FloatingEnvironmentGuard() { CHECK(std::fesetenv(&saved) == 0); }
	};

	void CheckStepBoundary(int sharpness) {
		CAPTURE(sharpness);
		const Real height = static_cast<Real>(sharpness);
		const Real halfWidth = Real(0.5) / height;
		const DiracStep delta(sharpness);
		CHECK(delta(0) == height);
		CHECK(delta(-Real(0)) == height);
		CHECK(delta(halfWidth) == height);
		CHECK(delta(-halfWidth) == height);
		CHECK(delta(std::nextafter(halfWidth, Real(0))) == height);
		CHECK(delta(-std::nextafter(halfWidth, Real(0))) == height);
		CHECK(delta(std::nextafter(halfWidth, std::numeric_limits<Real>::infinity())) == 0);
		CHECK(delta(-std::nextafter(halfWidth, std::numeric_limits<Real>::infinity())) == 0);
	}

	void CheckCauchyValues(int sharpness) {
		const Real height = static_cast<Real>(sharpness) / Constants::PI;
		const DiracSqr delta(sharpness);
		for (Real offset : {Real(0), Real(0.5), Real(-0.5), Real(1), Real(-1), Real(2), Real(-2)}) {
			CAPTURE(sharpness, offset);
			const Real value = delta(offset / static_cast<Real>(sharpness));
			CHECK(std::isfinite(value));
			CHECK(value > 0);
			CHECK_THAT(value, WithinRel(height / (1 + offset * offset), Precision::RoundoffZeroThreshold));
		}
	}
}

TEST_CASE("Floating factorial preserves every finite multiplication result", "[integer-boundaries][Factorial][factorial-range-safe]") {
	const int lastFinite = LastFiniteFactorialArgument();
	CAPTURE(lastFinite, std::numeric_limits<Real>::max_exponent);
	REQUIRE(std::isfinite(Functions::Factorial(lastFinite)));
	Real expected = 1;
	for (int n = 0; n <= lastFinite; ++n) {
		CAPTURE(n);
		if (n >= 2)
			expected *= Real(n);
		CHECK(Functions::Factorial(n) == expected);
	}
}

TEST_CASE("Floating factorial checks overflow before multiplication in every rounding mode", "[integer-boundaries][Factorial][factorial-range-safe]") {
	const int firstOverflow = LastFiniteFactorialArgument() + 1;
	FloatingEnvironmentGuard environment;
	for (int mode : {FE_TONEAREST, FE_DOWNWARD, FE_UPWARD, FE_TOWARDZERO}) {
		CAPTURE(mode);
		REQUIRE(std::fesetround(mode) == 0);
		for (int argument : {firstOverflow, std::numeric_limits<Real>::max_exponent}) {
			CAPTURE(argument);
			REQUIRE(std::feclearexcept(FE_ALL_EXCEPT) == 0);
			volatile int input = argument;
			const Real value = Functions::Factorial(input);
			const int overflow = std::fetestexcept(FE_OVERFLOW);
			CHECK(value == Constants::PosInf);
			CHECK(overflow == 0);
		}
	}
}

TEST_CASE("Floating factorial terminates at extreme positive integer inputs", "[integer-boundaries][Factorial][factorial-extreme]") {
	for (int n : {std::numeric_limits<int>::max() / 2, std::numeric_limits<int>::max() - 1, std::numeric_limits<int>::max()}) {
		CAPTURE(n);
		CHECK(Functions::Factorial(n) == Constants::PosInf);
	}
}

TEST_CASE("Floating factorial rejects negative integer boundaries", "[integer-boundaries][Factorial][factorial-range-safe]") {
	for (int n : {std::numeric_limits<int>::min(), -100, -1}) {
		CAPTURE(n);
		CHECK_THROWS_AS(Functions::Factorial(n), DomainError);
	}
}

TEST_CASE("Exact factorial retains every representable signed 64-bit value", "[integer-boundaries][integer-safe]") {
	constexpr long long expected[] = {
		1LL, 1LL, 2LL, 6LL, 24LL, 120LL, 720LL, 5040LL, 40320LL, 362880LL,
		3628800LL, 39916800LL, 479001600LL, 6227020800LL, 87178291200LL,
		1307674368000LL, 20922789888000LL, 355687428096000LL,
		6402373705728000LL, 121645100408832000LL, 2432902008176640000LL
	};
	for (int n = 0; n <= 20; ++n) {
		CAPTURE(n);
		CHECK(Functions::FactorialInt(n) == expected[n]);
	}
	for (int n : {std::numeric_limits<int>::min(), -100, -1})
		CHECK_THROWS_AS(Functions::FactorialInt(n), DomainError);
}

TEST_CASE("Exact factorial rejects overflow before multiplication", "[integer-boundaries][integer-overflow]") {
	for (int n : {21, 22, 30, 100, std::numeric_limits<int>::max() / 2, std::numeric_limits<int>::max()}) {
		CAPTURE(n);
		CHECK_THROWS_AS(Functions::FactorialInt(n), DomainError);
	}
}

TEST_CASE("Dirac step preserves inclusive ordinary support boundaries", "[integer-boundaries][integer-safe]") {
	for (int n : {1, 10, 46340, std::numeric_limits<int>::max() / 2})
		CheckStepBoundary(n);
}

TEST_CASE("Dirac step doubles sharpness without signed integer overflow", "[integer-boundaries][integer-overflow]") {
	for (int n : {std::numeric_limits<int>::max() / 2 + 1, std::numeric_limits<int>::max() - 1, std::numeric_limits<int>::max()})
		CheckStepBoundary(n);
}

TEST_CASE("Dirac Cauchy preserves ordinary scaled values", "[integer-boundaries][integer-safe]") {
	for (int n : {1, 10, 46340})
		CheckCauchyValues(n);
}

TEST_CASE("Dirac Cauchy squares sharpness without signed integer overflow", "[integer-boundaries][integer-overflow]") {
	for (int n : {46341, 50000, std::numeric_limits<int>::max() / 2, std::numeric_limits<int>::max()})
		CheckCauchyValues(n);
}

TEST_CASE("Other Dirac kernels retain their values at maximum sharpness", "[integer-boundaries][integer-safe]") {
	const int n = std::numeric_limits<int>::max();
	const Real height = static_cast<Real>(n);
	const Real x = Real(1) / height;
	const DiracExp gaussian(n);
	const DiracSin sinc(n);
	CHECK_THAT(gaussian(0), WithinRel(height / std::sqrt(2 * Constants::PI), Precision::RoundoffZeroThreshold));
	CHECK_THAT(gaussian(x), WithinRel(height * Real(0.24197072451914334979783019293556L), Precision::RoundoffZeroThreshold));
	CHECK(gaussian(-x) == gaussian(x));
	CHECK_THAT(sinc(0), WithinRel(height / Constants::PI, Precision::RoundoffZeroThreshold));
	CHECK_THAT(sinc(x), WithinRel(height * Real(0.26784853340116378615474613099694L), Precision::RoundoffZeroThreshold));
	CHECK(sinc(-x) == sinc(x));
}

TEST_CASE("Dirac kernels still reject nonpositive sharpness", "[integer-boundaries][integer-safe]") {
	for (int n : {std::numeric_limits<int>::min(), -1, 0}) {
		CAPTURE(n);
		CHECK_THROWS_AS(DiracStep(n), ArgumentError);
		CHECK_THROWS_AS(DiracExp(n), ArgumentError);
		CHECK_THROWS_AS(DiracSqr(n), ArgumentError);
		CHECK_THROWS_AS(DiracSin(n), ArgumentError);
	}
}

#include "../TestPrecision.h"
#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/base/StandardFunctions.h>
#endif
#include <catch2/catch_all.hpp>
#include <array>
#include <complex>
#include <type_traits>

using namespace MML;

namespace {
	template<typename T>
	struct ReciprocalCase {
		const char* name;
		std::complex<T> (*function)(std::complex<T>);
		std::complex<T> (*denominator)(std::complex<T>);
		std::complex<T> expected;
		bool odd;
	};

	template<typename T>
	auto ReciprocalCases() {
		using Z = std::complex<T>;
		using namespace Functions;
		return std::array<ReciprocalCase<T>, 6>{{
			{"Sec", Sec<Z>, Cos<Z>, Z(T(0.498337030555186785213805891772169534L), T(0.591083841721045048050391692974335072L)), false},
			{"Csc", Csc<Z>, Sin<Z>, Z(T(0.621518017170428421234907805855920148L), T(-0.303931001628426450334485604509703273L)), true},
			{"Ctg", Ctg<Z>, Tan<Z>, Z(T(0.217621561854402681365134243605238074L), T(-0.868014142895924948635849208916273888L)), true},
			{"Sech", Sech<Z>, Cosh<Z>, Z(T(0.498337030555186785213805891772169534L), T(-0.591083841721045048050391692974335072L)), false},
			{"Csch", Csch<Z>, Sinh<Z>, Z(T(0.303931001628426450334485604509703273L), T(-0.621518017170428421234907805855920148L)), true},
			{"Ctgh", Ctgh<Z>, Tanh<Z>, Z(T(0.868014142895924948635849208916273888L), T(-0.217621561854402681365134243605238074L)), true}
		}};
	}

	template<typename T>
	void CheckComplex(std::complex<T> actual, std::complex<T> expected) {
		CHECK(std::isfinite(actual.real()));
		CHECK(std::isfinite(actual.imag()));
		// Keep the component type: MML's scalar comparison helpers use the build's Real.
		CHECK(std::abs(actual - expected) <= PrecisionValues<T>::RoundoffZeroThreshold * std::max(T(1), std::abs(expected)));
	}
}

TEMPLATE_TEST_CASE("Complex reciprocals match independent references and symmetry", "[complex-reciprocals]", float, double, long double) {
	using Z = std::complex<TestType>;
	const Z z(1, 1);
	for (const auto& entry : ReciprocalCases<TestType>()) {
		CAPTURE(entry.name, sizeof(TestType));
		STATIC_REQUIRE(std::is_same_v<decltype(entry.function(z)), Z>);
		CheckComplex(entry.function(z), entry.expected);
		CheckComplex(entry.function(std::conj(z)), std::conj(entry.expected));
		CheckComplex(entry.function(-z), entry.odd ? -entry.expected : entry.expected);
	}
}

TEMPLATE_TEST_CASE("Odd complex reciprocal guards use magnitude and strict machine epsilon", "[complex-reciprocals]", float, double, long double) {
	using T = TestType;
	using Z = std::complex<T>;
	const T eps = std::numeric_limits<T>::epsilon();
	for (const auto& entry : ReciprocalCases<T>()) {
		if (!entry.odd)
			continue;
		CAPTURE(entry.name, sizeof(T));
		for (Z z : {Z(0), Z(eps / 2, 0), Z(0, eps / 2), Z(eps / 2, eps / 2)}) {
			CAPTURE(z);
			CHECK_THROWS_AS(entry.function(z), DomainError);
		}
		for (Z z : {Z(eps, 0), Z(0, eps), Z(2 * eps, 0), Z(0, 2 * eps), Z(3 * eps / 4, 3 * eps / 4)}) {
			CAPTURE(z);
			const Z denominator = entry.denominator(z);
			// libm can round a denominator at the boundary to either side of epsilon.
			if (std::abs(denominator) < eps)
				CHECK_THROWS_AS(entry.function(z), DomainError);
			else
				CheckComplex(entry.function(z), Z(1) / denominator);
		}
	}
}

TEMPLATE_TEST_CASE("Complex secant guards the magnitude near its poles", "[complex-reciprocals]", float, double, long double) {
	using T = TestType;
	using Z = std::complex<T>;
	const T eps = std::numeric_limits<T>::epsilon();
	const T halfPi = std::acos(T(0));
	const Z inside(halfPi, eps / 2);
	REQUIRE(std::abs(std::cos(inside)) < eps);
	CHECK_THROWS_AS(Functions::Sec(inside), DomainError);
	CHECK_THROWS_AS(Functions::Sec(std::conj(inside)), DomainError);
	const Z outside(halfPi, 2 * eps);
	REQUIRE(std::abs(std::cos(outside)) >= eps);
	CheckComplex(Functions::Sec(outside), Z(1) / std::cos(outside));
}

TEMPLATE_TEST_CASE("Real reciprocal values and existing guards are preserved", "[complex-reciprocals]", float, double, long double) {
	using T = TestType;
	const T x = 1;
	const T eps = std::numeric_limits<T>::epsilon();
	CHECK(Functions::Sec(x) == T(1) / std::cos(x));
	CHECK(Functions::Csc(x) == T(1) / std::sin(x));
	CHECK(Functions::Ctg(x) == T(1) / std::tan(x));
	CHECK(Functions::Sech(x) == T(1) / std::cosh(x));
	CHECK(Functions::Csch(x) == T(1) / std::sinh(x));
	CHECK(Functions::Ctgh(x) == T(1) / std::tanh(x));
	CHECK_THROWS_AS(Functions::Csc(eps / 2), DomainError);
	CHECK_THROWS_AS(Functions::Ctg(eps / 2), DomainError);
	CHECK_THROWS_AS(Functions::Csch(eps / 2), DomainError);
	CHECK_THROWS_AS(Functions::Ctgh(eps / 2), DomainError);
	CHECK_NOTHROW(Functions::Csc(eps));
	CHECK_NOTHROW(Functions::Ctg(eps));
	CHECK_NOTHROW(Functions::Csch(eps));
	CHECK_NOTHROW(Functions::Ctgh(eps));
	CHECK_THROWS_AS(Functions::Sec(std::acos(T(0))), DomainError);
	CHECK(Functions::Sec(T(0)) == 1);
	CHECK(Functions::Sech(T(0)) == 1);
}

TEMPLATE_TEST_CASE("Complex hyperbolic secant retains its unguarded reciprocal behavior", "[complex-reciprocals]", float, double, long double) {
	using T = TestType;
	using Z = std::complex<T>;
	const Z nearPole(0, std::acos(T(0)));
	REQUIRE(std::abs(std::cosh(nearPole)) < std::numeric_limits<T>::epsilon());
	CheckComplex(Functions::Sech(nearPole), Z(1) / std::cosh(nearPole));
}

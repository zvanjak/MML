#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "../../TestPrecision.h"

#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/core/Derivation/ADAdapters.h>
#include <mml/core/Derivation/ADJacobians.h>
#endif

#include <cmath>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>

using namespace MML;
using namespace MML::AD;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {
	struct Identity {
		template<typename T> T operator()(T x) const { return x; }
	};
	struct ScalarExpression {
		template<typename V> auto operator()(const V& x) const { return x[0] * x[0]; }
	};
	struct VectorExpression {
		template<typename V> auto operator()(const V& x) const { return x; }
	};
	template<typename T, typename = void> struct HasDerivative : std::false_type {};
	template<typename T> struct HasDerivative<T, std::void_t<decltype(CalcDerivativeForwardAD(Identity{}, std::declval<T>()))>> : std::true_type {};
	template<typename V, typename = void> struct HasGradient : std::false_type {};
	template<typename V> struct HasGradient<V, std::void_t<decltype(CalcGradientForwardAD(ScalarExpression{}, std::declval<V>()))>> : std::true_type {};
	template<typename V, typename = void> struct HasJacobian : std::false_type {};
	template<typename V> struct HasJacobian<V, std::void_t<decltype(CalcJacobianForwardAD(VectorExpression{}, std::declval<V>()))>> : std::true_type {};
	template<typename V, typename = void> struct HasHessian : std::false_type {};
	template<typename V> struct HasHessian<V, std::void_t<decltype(CalcHessianForwardAD(ScalarExpression{}, std::declval<V>()))>> : std::true_type {};
}

TEST_CASE("Forward AD helpers reject integer deduction across scalar and vector surfaces", "[ad-domain-contract]") {
	CHECK_FALSE(HasDerivative<int>::value);
	CHECK_FALSE(HasDerivative<bool>::value);
	CHECK_FALSE(HasDerivative<unsigned long long>::value);
	CHECK_FALSE(AD::Detail::HasFloatingLeaf<Dual<int>>::value);
	CHECK_FALSE(HasGradient<std::vector<int>>::value);
	CHECK_FALSE(HasJacobian<std::vector<int>>::value);
	CHECK_FALSE(HasJacobian<Vector<int>>::value);
	CHECK_FALSE(HasHessian<std::vector<int>>::value);
	CHECK_FALSE(HasHessian<Vector<int>>::value);
	CHECK(HasDerivative<float>::value);
	CHECK(HasDerivative<double>::value);
	CHECK(HasDerivative<long double>::value);
	CHECK(HasDerivative<Dual<Real>>::value);
	static_assert(std::is_same_v<decltype(CalcDerivativeForwardAD(Identity{}, float(2))), std::pair<float, float>>);
	const auto result = CalcDerivativeForwardAD([](auto x) { using std::sin; return sin(x); }, Real(2));
	CHECK_THAT(result.first, WithinRel(std::sin(Real(2)), TOL(1e-12, 1e-4)));
	CHECK_THAT(result.second, WithinRel(std::cos(Real(2)), TOL(1e-12, 1e-4)));
	const Dual<Real> x(2, 1);
	CHECK((x + 3).value == 5);
	CHECK((3 * x).deriv == 3);
}

TEST_CASE("Raw inverse trig endpoints expose signed singularities rather than finite clamps", "[ad-domain-contract]") {
	for (Real endpoint : {Real(-1), Real(1)}) {
		for (Real seed : {Real(-1), Real(1)}) {
			const auto a = asin(Dual<Real>(endpoint, seed));
			const auto c = acos(Dual<Real>(endpoint, seed));
			CHECK(a.value == std::asin(endpoint));
			CHECK(c.value == std::acos(endpoint));
			CHECK(std::isinf(a.deriv));
			CHECK(std::isinf(c.deriv));
			CHECK(std::signbit(a.deriv) == std::signbit(seed));
			CHECK(std::signbit(c.deriv) != std::signbit(seed));
		}
		CHECK(std::isnan(asin(Dual<Real>(endpoint, 0)).deriv));
		CHECK(std::isnan(acos(Dual<Real>(endpoint, 0)).deriv));
		Tape tape;
		ADVar x(&tape, endpoint);
		asin(x).backward();
		CHECK(std::isinf(x.adjoint()));
		CHECK_FALSE(std::signbit(x.adjoint()));
		acos(x).backward();
		CHECK(std::isinf(x.adjoint()));
		CHECK(std::signbit(x.adjoint()));
	}
}

TEST_CASE("Inverse trig invalid inputs and near endpoints retain raw IEEE and interior derivatives", "[ad-domain-contract]") {
	for (Real x : {Real(-1.5), Real(1.5)}) {
		CHECK(std::isnan(asin(Dual<Real>(x, 1)).value));
		CHECK(std::isnan(asin(Dual<Real>(x, 1)).deriv));
		CHECK(std::isnan(acos(Dual<Real>(x, 1)).deriv));
	}
	for (Real x : {std::nextafter(Real(-1), Real(0)), Real(0), std::nextafter(Real(1), Real(0))}) {
		const Real expected = 1 / std::sqrt((1 - x) * (1 + x));
		CHECK_THAT(asin(Dual<Real>(x, 1)).deriv, WithinRel(expected, TOL(1e-12, 1e-4)));
		CHECK_THAT(acos(Dual<Real>(x, 1)).deriv, WithinRel(-expected, TOL(1e-12, 1e-4)));
	}
}

TEST_CASE("Checked AD fields reject singular first and nested inverse trig derivatives", "[ad-domain-contract]") {
	for (Real endpoint : {Real(-1), Real(1)}) {
		const auto a = MakeADScalarField<1>([](const auto& x) { using std::asin; return asin(x[0]); });
		const auto c = MakeADScalarField<1>([](const auto& x) { using std::acos; return acos(x[0]); });
		CHECK_THROWS_AS(a.EvaluatePartial({endpoint}, 0), DomainError);
		CHECK_THROWS_AS(a.EvaluateSecondPartial({endpoint}, 0, 0), DomainError);
		CHECK_THROWS_AS(c.EvaluatePartial({endpoint}, 0), DomainError);
		CHECK_THROWS_AS(c.EvaluateSecondPartial({endpoint}, 0, 0), DomainError);
		const auto h = CalcHessianForwardAD([](const auto& x) { using std::asin; return asin(x[0]); }, std::vector<Real>{endpoint});
		CHECK_FALSE(std::isfinite(h(0, 0)));
	}
}

TEST_CASE("Fixed zero power has exact zero forward and reverse derivatives at every finite base", "[ad-domain-contract]") {
	for (Real value : {Real(-2), Real(0), Real(2)}) {
		const Dual<Real> x(value, 1);
		for (const auto& y : {pow(x, 0), pow(x, Dual<Real>(0)), pow(value, Dual<Real>(0))}) {
			CHECK(y.value == 1);
			CHECK(y.deriv == 0);
		}
		Tape tape;
		ADVar variable(&tape, value);
		auto y = pow(variable, 0);
		CHECK(y.value() == 1);
		y.backward();
		CHECK(variable.adjoint() == 0);
		auto zero = ADVar::constant(&tape, 2) - ADVar::constant(&tape, 2);
		CHECK_NOTHROW(y = pow(variable, zero));
		y.backward();
		CHECK(variable.adjoint() == 0);
	}
}

TEST_CASE("Undefined fixed real powers are rejected consistently by all overloads", "[ad-domain-contract]") {
	for (const auto& point : {std::pair<Real, Real>{Real(-2), Real(0.5)}, {Real(-2), Real(-0.5)}, {Real(0), Real(-1)}, {Real(0), Real(-0.5)}}) {
		const auto [base, exponent] = point;
		CHECK_THROWS_AS(pow(Dual<Real>(base, 1), exponent), DomainError);
		CHECK_THROWS_AS(pow(Dual<Real>(base, 1), Dual<Real>(exponent)), DomainError);
		CHECK_THROWS_AS(pow(base, Dual<Real>(exponent)), DomainError);
		Tape tape;
		CHECK_THROWS_AS(pow(ADVar(&tape, base), exponent), DomainError);
		CHECK_THROWS_AS(pow(ADVar(base), ADVar(exponent)), DomainError);
		CHECK_THROWS_AS(pow(base, ADVar::constant(&tape, exponent)), DomainError);
	}
}

TEST_CASE("Variable exponents require positive bases including hidden nested seeds", "[ad-domain-contract]") {
	for (Real base : {Real(-2), Real(0)}) {
		for (Real exponent : {Real(0), Real(2), Real(3)}) {
			CHECK_THROWS_AS(pow(Dual<Real>(base, 1), Dual<Real>(exponent, 1)), DomainError);
			CHECK_THROWS_AS(pow(base, Dual<Real>(exponent, 1)), DomainError);
			Tape tape;
			ADVar y(&tape, exponent);
			CHECK_THROWS_AS(pow(base, y), DomainError);
			CHECK_THROWS_AS(pow(ADVar(&tape, base), y), DomainError);
			CHECK_THROWS_AS(pow(base, y + 1), DomainError);
		}
		using D = Dual<Real>;
		using D2 = Dual<D>;
		CHECK_THROWS_AS(pow(D2(D(base), D(0)), D2(D(3, 1), D(0))), DomainError);
		CHECK_THROWS_AS(pow(D2(D(base), D(0)), D2(D(3), D(0, 1))), DomainError);
	}
}

TEST_CASE("Valid fixed negative-base powers and zero-base polynomials preserve derivatives", "[ad-domain-contract]") {
	for (Real exponent : {Real(-2), Real(-1), Real(1), Real(2), Real(3)}) {
		const auto y = pow(Dual<Real>(-2, 1), exponent);
		CHECK_THAT(y.value, WithinRel(std::pow(Real(-2), exponent), TOL(1e-12, 1e-4)));
		const Real expected = exponent * std::pow(Real(-2), exponent - 1);
		CHECK_THAT(y.deriv, WithinRel(expected, TOL(1e-12, 1e-4)));
		Tape tape;
		ADVar x(&tape, -2);
		auto fixed = ADVar::constant(&tape, exponent - 1) + ADVar::constant(&tape, 1);
		auto reverse = pow(x, fixed);
		reverse.backward();
		CHECK_THAT(x.adjoint(), WithinRel(expected, TOL(1e-12, 1e-4)));
	}
	for (Real exponent : {Real(0), Real(1), Real(2), Real(3)}) {
		const auto h = CalcHessianForwardAD([=](const auto& x) { using std::pow; return pow(x[0], exponent); }, std::vector<Real>{0});
		CHECK(h(0, 0) == (exponent == 2 ? 2 : 0));
		const auto d = CalcDerivativeForwardAD([=](auto x) { using std::pow; return pow(x, exponent); }, Real(0));
		CHECK(d.second == (exponent == 1 ? 1 : 0));
	}
}

TEST_CASE("Variable power Hessians retain mixed seeds at a zero exponent and positive base", "[ad-domain-contract]") {
	const auto f = [](const auto& x) { using std::pow; return pow(x[0], x[1]); };
	const auto h = CalcHessianForwardAD(f, std::vector<Real>{2, 0});
	CHECK(h(0, 0) == 0);
	CHECK_THAT(h(0, 1), WithinRel(Real(0.5), TOL(1e-12, 1e-4)));
	CHECK(h(1, 0) == h(0, 1));
	CHECK_THAT(h(1, 1), WithinRel(std::log(Real(2)) * std::log(Real(2)), TOL(1e-12, 1e-4)));
	const auto gradient = CalcGradientReverseAD(f, std::vector<Real>{2, 0});
	CHECK(gradient[0] == 0);
	CHECK_THAT(gradient[1], WithinRel(std::log(Real(2)), TOL(1e-12, 1e-4)));
}

TEST_CASE("Recorded fixed powers follow constant ancestry without a backward sweep", "[ad-domain-contract]") {
	Tape tape;
	ADVar unrelated(&tape, 17), x(&tape, -2);
	const auto fixed = sin(ADVar::constant(&tape, 0)) + ADVar::constant(&tape, 2);
	CHECK(tape.IsConstant(fixed.index()));
	CHECK_FALSE(tape.IsConstant(x.index()));
	auto y = pow(x, fixed);
	y.backward();
	CHECK(y.value() == 4);
	CHECK(x.adjoint() == -4);
	CHECK(unrelated.adjoint() == 0);
	const auto variableZero = x - x;
	const auto size = tape.size();
	CHECK_THROWS_AS(pow(x, variableZero), DomainError);
	CHECK(tape.size() == size);
	CHECK_THROWS_AS(tape.IsConstant(-1), ArgumentError);
	CHECK_THROWS_AS(tape.IsConstant(static_cast<int>(tape.size())), ArgumentError);
}

TEST_CASE("Fixed powers at zero preserve third derivatives and singular higher jets", "[ad-domain-contract]") {
	using D = Dual<Real>;
	using D2 = Dual<D>;
	const Dual<D2> x(D2(D(0, 1), D(1, 0)), D2(D(1, 0), D(0, 0)));
	for (int exponent : {0, 1, 2, 3}) {
		const auto y = pow(x, exponent);
		CHECK(y.value.value.value == (exponent == 0 ? 1 : 0));
		CHECK(y.deriv.value.value == (exponent == 1 ? 1 : 0));
		CHECK(y.deriv.deriv.value == (exponent == 2 ? 2 : 0));
		CHECK(y.deriv.deriv.deriv == (exponent == 3 ? 6 : 0));
	}
}

TEST_CASE("Power inputs are finite while valid singular and overflowing results remain explicit", "[ad-domain-contract]") {
	for (Real bad : {std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity()}) {
		CHECK_THROWS_AS(pow(Dual<Real>(bad, 1), 0), DomainError);
		CHECK_THROWS_AS(pow(Dual<Real>(2, 1), bad), DomainError);
		CHECK_THROWS_AS(pow(ADVar(bad), 0), DomainError);
		CHECK_THROWS_AS(pow(ADVar(2), bad), DomainError);
	}
	CHECK(std::isinf(pow(Dual<Real>(0, 1), Real(0.5)).deriv));
	using D = Dual<Real>;
	const auto nested = pow(Dual<D>(D(0, 1), D(1)), Real(0.5));
	CHECK_FALSE(std::isfinite(nested.deriv.deriv));
	const auto overflow = pow(Dual<Real>(2, 1), Real(2 * std::numeric_limits<Real>::max_exponent));
	CHECK(std::isinf(overflow.value));
	CHECK(std::isinf(overflow.deriv));
}

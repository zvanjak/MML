#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "../../TestPrecision.h"

#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/core/Derivation/ADJacobians.h>
#endif

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace MML;
using namespace MML::AD;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

TEST_CASE("Forward free and compound division retain large and tiny quotient derivatives", "[ad-quotient-baseline][ad-quotient-contract]") {
	for (Real scale : {Real(1), TOL(1e200, 1e20), TOL(1e-200, 1e-20), TOL(1e300, 1e30), TOL(1e-300, 1e-30)}) {
		for (Real direction : {Real(-1), Real(1)}) {
			const Dual<Real> numerator(scale, direction), denominator(scale, 0);
			const auto ordinary = numerator / denominator;
			auto compound = numerator;
			compound /= denominator;
			CHECK(ordinary.value == 1);
			CHECK(compound.value == 1);
			CHECK_THAT(ordinary.deriv, WithinRel(direction / scale, TOL(1e-12, 1e-4)));
			CHECK_THAT(compound.deriv, WithinRel(direction / scale, TOL(1e-12, 1e-4)));
			const Dual<Real> activeDenominator(scale, direction);
			const auto inverse = scale / activeDenominator;
			CHECK(inverse.value == 1);
			CHECK_THAT(inverse.deriv, WithinRel(-direction / scale, TOL(1e-12, 1e-4)));
		}
	}
}

TEST_CASE("Reverse division retains both representable partial derivatives across scales", "[ad-quotient-baseline][ad-quotient-contract]") {
	for (Real scale : {Real(1), TOL(1e200, 1e20), TOL(1e-200, 1e-20), TOL(1e300, 1e30), TOL(1e-300, 1e-30)}) {
		Tape tape;
		ADVar numerator(&tape, scale), denominator(&tape, scale);
		auto quotient = numerator / denominator;
		quotient.backward();
		CHECK(quotient.value() == 1);
		CHECK_THAT(numerator.adjoint(), WithinRel(1 / scale, TOL(1e-12, 1e-4)));
		CHECK_THAT(denominator.adjoint(), WithinRel(-1 / scale, TOL(1e-12, 1e-4)));
	}
}

TEST_CASE("Reverse gradient accepts a constant output with the original input dimension", "[ad-constant-baseline][ad-quotient-contract]") {
	int calls = 0;
	const auto constant = [&](const std::vector<ADVar>&) { ++calls; return ADVar(7); };
	std::vector<Real> gradient;
	REQUIRE_NOTHROW(gradient = CalcGradientReverseAD(constant, std::vector<Real>{2, 3, 4}));
	CHECK(calls == 1);
	REQUIRE(gradient.size() == 3);
	for (Real value : gradient) CHECK(value == 0);
}

TEST_CASE("Reverse Jacobians accept interleaved constant and variable output rows", "[ad-constant-baseline][ad-quotient-contract]") {
	int calls = 0;
	const auto mixed = [&](const std::vector<ADVar>& x) {
		++calls;
		return std::vector<ADVar>{ADVar(7), x[0] / x[1], ADVar(0), x[0] + x[1]};
	};
	Matrix<Real> jacobian;
	REQUIRE_NOTHROW(jacobian = CalcJacobianReverseAD(mixed, Vector<Real>{3, 2}));
	REQUIRE(jacobian.rows() == 4);
	REQUIRE(jacobian.cols() == 2);
	CHECK(calls == 5);
	const Real expected[4][2] = {{0, 0}, {Real(0.5), Real(-0.75)}, {0, 0}, {1, 1}};
	for (int row = 0; row < 4; ++row)
		for (int col = 0; col < 2; ++col) CHECK(jacobian(row, col) == expected[row][col]);
}

TEST_CASE("AD division preserves representable quotients and tangents at floating range limits", "[ad-quotient-contract]") {
	for (Real scale : {std::numeric_limits<Real>::denorm_min(), std::numeric_limits<Real>::min(), std::numeric_limits<Real>::max() / 4}) {
		const Dual<Real> numerator(scale, scale), denominator(scale, 0);
		const auto ordinary = numerator / denominator;
		auto compound = numerator;
		compound /= denominator;
		CHECK(ordinary.value == 1);
		CHECK(ordinary.deriv == 1);
		CHECK(compound.value == ordinary.value);
		CHECK(compound.deriv == ordinary.deriv);
		const auto reciprocal = scale / Dual<Real>(scale, scale);
		CHECK(reciprocal.value == 1);
		CHECK(reciprocal.deriv == -1);
	}
}

TEST_CASE("Forward scalar division and self-aliasing share quotient arithmetic", "[ad-quotient-contract]") {
	for (Real scale : {TOL(1e200, 1e20), TOL(1e-200, 1e-20), std::numeric_limits<Real>::denorm_min()}) {
		auto x = Dual<Real>(scale, scale);
		const auto freeScalar = x / scale;
		auto compoundScalar = x;
		compoundScalar /= scale;
		CHECK(freeScalar.value == 1);
		CHECK(freeScalar.deriv == 1);
		CHECK(compoundScalar.value == freeScalar.value);
		CHECK(compoundScalar.deriv == freeScalar.deriv);
		const auto self = x / x;
		x /= x;
		CHECK(self.value == 1);
		CHECK(self.deriv == 0);
		CHECK(x.value == 1);
		CHECK(x.deriv == 0);
	}
}

TEST_CASE("Forward quotient cancellation retains the small determinant of nearby values", "[ad-quotient-contract]") {
	for (Real scale : {Real(1), TOL(1e300, 1e30), TOL(1e-300, 1e-30)}) {
		const Real nearby = std::nextafter(scale, std::numeric_limits<Real>::infinity());
		const Real expected = (scale - nearby) / scale;
		const Dual<Real> a(nearby, scale), b(scale, scale);
		const auto quotient = a / b;
		auto compound = a;
		compound /= b;
		REQUIRE(expected != 0);
		CHECK_THAT(quotient.deriv, WithinRel(expected, TOL(1e-12, 1e-4)));
		CHECK(compound.deriv == quotient.deriv);
	}

}

TEST_CASE("Quotient tangents preserve numerator and denominator signs and exact zero numerators", "[ad-quotient-contract]") {
	for (Real scale : {TOL(1e200, 1e20), TOL(1e-200, 1e-20)})
		for (Real aSign : {Real(-1), Real(1)})
			for (Real bSign : {Real(-1), Real(1)}) {
				const Real a = aSign * scale, b = bSign * scale;
				const auto numerator = Dual<Real>(a, 1) / Dual<Real>(b, 0);
				const auto denominator = a / Dual<Real>(b, 1);
				CHECK(numerator.value == aSign / bSign);
				CHECK_THAT(numerator.deriv, WithinRel(bSign / scale, TOL(1e-12, 1e-4)));
				CHECK_THAT(denominator.deriv, WithinRel(-aSign / scale, TOL(1e-12, 1e-4)));
				const auto zero = Real(0) / Dual<Real>(b, 1);
				CHECK(zero.value == 0);
				CHECK(zero.deriv == 0);
			}
}

TEST_CASE("Reverse quotient propagation retains finite chained adjoints without product overflow or underflow", "[ad-quotient-contract]") {
	struct Case { Real a, b, weight, da, db; };
	const std::vector<Case> cases{
		{TOL(1e100, 1e10), TOL(1e200, 1e20), TOL(1e300, 1e30), TOL(1e100, 1e10), Real(-1)},
		{TOL(1e-100, 1e-10), TOL(1e-200, 1e-20), TOL(1e-300, 1e-30), TOL(1e-100, 1e-10), Real(-1)},
		{std::numeric_limits<Real>::denorm_min(), std::numeric_limits<Real>::denorm_min(),
			std::numeric_limits<Real>::denorm_min(), Real(1), Real(-1)}
	};
	for (const auto& test : cases) {
		Tape tape;
		ADVar a(&tape, test.a), b(&tape, test.b);
		const auto output = (a / b) * test.weight;
		CHECK(std::isfinite(output.value()));
		output.backward();
		CHECK_THAT(a.adjoint(), WithinRel(test.da, TOL(1e-12, 1e-4)));
		CHECK_THAT(b.adjoint(), WithinRel(test.db, TOL(1e-12, 1e-4)));
	}
}

TEST_CASE("Reverse free compound scalar and self division preserve graph derivatives", "[ad-quotient-contract]") {
	for (Real scale : {TOL(1e200, 1e20), TOL(1e-200, 1e-20)}) {
		Tape tape;
		ADVar a(&tape, scale), b(&tape, scale);
		auto quotient = a;
		quotient /= b;
		quotient.backward();
		CHECK(quotient.value() == 1);
		CHECK_THAT(a.adjoint(), WithinRel(1 / scale, TOL(1e-12, 1e-4)));
		CHECK_THAT(b.adjoint(), WithinRel(-1 / scale, TOL(1e-12, 1e-4)));
		auto scalar = a / scale;
		scalar.backward();
		CHECK_THAT(a.adjoint(), WithinRel(1 / scale, TOL(1e-12, 1e-4)));
		const auto inverse = scale / b;
		inverse.backward();
		CHECK_THAT(b.adjoint(), WithinRel(-1 / scale, TOL(1e-12, 1e-4)));
	}
	for (Real scale : {TOL(1e200, 1e20), -TOL(1e200, 1e20), std::numeric_limits<Real>::denorm_min()}) {
		Tape tape;
		ADVar x(&tape, scale);
		auto self = x;
		self /= self;
		self.backward();
		CHECK(self.value() == 1);
		CHECK(x.adjoint() == 0);
	}
}

TEST_CASE("Forward and reverse high-level Jacobians keep finite small quotient entries", "[ad-quotient-contract]") {
	const auto forward = [](const std::vector<Dual<Real>>& x) {
		return std::vector<Dual<Real>>{x[0] / x[1], Dual<Real>(7)};
	};
	const auto reverse = [](const std::vector<ADVar>& x) {
		return std::vector<ADVar>{x[0] / x[1], ADVar(7)};
	};
	for (Real scale : {TOL(1e200, 1e20), TOL(1e-200, 1e-20)}) {
		const std::vector<Real> point{scale, scale};
		const auto f = CalcJacobianForwardAD(forward, point);
		const auto r = CalcJacobianReverseAD(reverse, point);
		REQUIRE(f.rows() == 2);
		REQUIRE(f.cols() == 2);
		REQUIRE(r.rows() == 2);
		REQUIRE(r.cols() == 2);
		CHECK_THAT(f(0, 0), WithinRel(1 / scale, TOL(1e-12, 1e-4)));
		CHECK_THAT(f(0, 1), WithinRel(-1 / scale, TOL(1e-12, 1e-4)));
		CHECK_THAT(r(0, 0), WithinRel(1 / scale, TOL(1e-12, 1e-4)));
		CHECK_THAT(r(0, 1), WithinRel(-1 / scale, TOL(1e-12, 1e-4)));
		CHECK(f(1, 0) == 0);
		CHECK(f(1, 1) == 0);
		CHECK(r(1, 0) == 0);
		CHECK(r(1, 1) == 0);
	}
}

TEST_CASE("Quotients preserve nested-dual derivatives and rational Hessians", "[ad-quotient-contract]") {
	using Dual2 = Dual<Dual<Real>>;
	for (Real scale : {TOL(1e200, 1e20), TOL(1e-200, 1e-20), std::numeric_limits<Real>::denorm_min()}) {
		const Dual2 numerator(Dual<Real>(scale, 0), Dual<Real>(0, 0));
		const Dual2 denominator(Dual<Real>(scale, scale), Dual<Real>(scale, 0));
		const auto quotient = numerator / denominator;
		auto compound = numerator;
		compound /= denominator;
		CHECK(quotient.value.value == 1);
		CHECK(quotient.value.deriv == -1);
		CHECK(quotient.deriv.value == -1);
		CHECK_THAT(quotient.deriv.deriv, WithinRel(Real(2), Precision::EqualityTolerance));
		CHECK(compound.value.value == quotient.value.value);
		CHECK(compound.value.deriv == quotient.value.deriv);
		CHECK(compound.deriv.value == quotient.deriv.value);
		CHECK(compound.deriv.deriv == quotient.deriv.deriv);
	}

	const auto freeQuotient = [](const std::vector<Dual2>& x) { return x[0] / x[1]; };
	const auto compoundQuotient = [](const std::vector<Dual2>& x) { auto q = x[0]; q /= x[1]; return q; };
	const auto hessian = CalcHessianForwardAD(freeQuotient, std::vector<Real>{3, 2});
	const auto compoundHessian = CalcHessianForwardAD(compoundQuotient, std::vector<Real>{3, 2});
	const Real expected[2][2] = {{0, Real(-0.25)}, {Real(-0.25), Real(0.75)}};
	for (int row = 0; row < 2; ++row)
		for (int col = 0; col < 2; ++col) {
			CHECK_THAT(hessian(row, col), WithinAbs(expected[row][col], TOL(1e-12, 1e-4)));
			CHECK(compoundHessian(row, col) == hessian(row, col));
		}
}

TEST_CASE("Nested quotient tangents retain large products and nearby-value cancellation", "[ad-quotient-contract]") {
	using Dual2 = Dual<Dual<Real>>;
	const Dual2 a(Dual<Real>(TOL(1e300, 1e30), 0), Dual<Real>(0, 0));
	const Dual2 b(Dual<Real>(TOL(1e200, 1e20), 0), Dual<Real>(TOL(1e300, 1e30), 0));
	const auto quotient = a / b;
	CHECK_THAT(quotient.value.value, WithinRel(TOL(1e100, 1e10), TOL(1e-12, 1e-4)));
	CHECK(quotient.value.deriv == 0);
	CHECK_THAT(quotient.deriv.value, WithinRel(-TOL(1e200, 1e20), TOL(1e-12, 1e-4)));
	CHECK(quotient.deriv.deriv == 0);
	const Real below = std::nextafter(Real(2), Real(0));
	const Dual2 nearbyNumerator(Dual<Real>(2, 0), Dual<Real>(1, 0));
	const Dual2 nearbyDenominator(Dual<Real>(below, 0), Dual<Real>(1, 0));
	const auto cancellation = nearbyNumerator / nearbyDenominator;
	const Real expected = ((below - 2) / below) / below;
	REQUIRE(expected != 0);
	CHECK_THAT(cancellation.deriv.value, WithinRel(expected, TOL(1e-12, 1e-4)));
	CHECK(cancellation.deriv.deriv == 0);
	auto self = b;
	self /= self;
	CHECK(self.value.value == 1);
	CHECK(self.value.deriv == 0);
	CHECK(self.deriv.value == 0);
	CHECK(self.deriv.deriv == 0);
}

TEST_CASE("Constant helpers ignore previous callback adjoints and preserve recorded constants", "[ad-quotient-contract]") {
	const auto constantAfterWork = [](std::vector<ADVar>& x) {
		(x[0] * x[0]).backward();
		return ADVar(7);
	};
	const auto gradient = CalcGradientReverseAD(constantAfterWork, std::vector<Real>{3});
	REQUIRE(gradient.size() == 1);
	CHECK(gradient[0] == 0);
	const auto rows = [](std::vector<ADVar>& x) {
		(x[0] * x[0]).backward();
		return std::vector<ADVar>{ADVar(7), ADVar::constant(x[0].tape(), 5), x[0] * 0};
	};
	const auto jacobian = CalcJacobianReverseAD(rows, std::vector<Real>{3});
	REQUIRE(jacobian.rows() == 3);
	REQUIRE(jacobian.cols() == 1);
	for (int row = 0; row < 3; ++row) CHECK(jacobian(row, 0) == 0);
	const auto empty = CalcGradientReverseAD([](const auto&) { return ADVar(7); }, std::vector<Real>{});
	CHECK(empty.empty());
}

TEST_CASE("Reverse constant helpers preserve callback exceptions and reject nonfinite constants", "[ad-quotient-contract]") {
	const std::vector<Real> point{2, 3};
	for (Real invalid : {std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity()}) {
		const auto scalar = [&](const auto&) { return ADVar(invalid); };
		const auto vector = [&](const auto&) { return std::vector<ADVar>{ADVar(invalid)}; };
		CHECK_THROWS_AS(CalcGradientReverseAD(scalar, point), NumericalMethodError);
		CHECK_THROWS_AS(CalcJacobianReverseAD(vector, point), NumericalMethodError);
	}
	const auto throwingScalar = [](const auto&) -> ADVar { throw std::logic_error("callback witness"); };
	const auto throwingVector = [](const auto&) -> std::vector<ADVar> { throw std::logic_error("callback witness"); };
	CHECK_THROWS_AS(CalcGradientReverseAD(throwingScalar, point), std::logic_error);
	CHECK_THROWS_AS(CalcJacobianReverseAD(throwingVector, point), std::logic_error);
	ADVar constant(7);
	CHECK_THROWS_AS(constant.backward(), ArgumentError);
	Tape tape;
	ADVar x(&tape, 1);
	CHECK_THROWS_AS(x / 0, DivisionByZeroError);
}

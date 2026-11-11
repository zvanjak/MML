#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "../../TestPrecision.h"

#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/core/Derivation/Jacobians.h>
#endif

#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <vector>

using namespace MML;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {
	Real LargeJacobianMagnitude() {
		return std::scalbn(Real(1.5), std::numeric_limits<Real>::max_exponent - 1);
	}
	Real SmallJacobianMagnitude() {
		return std::scalbn(Real(1), std::numeric_limits<Real>::min_exponent - 2);
	}
	class LargeIdentity : public IVectorFunction<1> {
	public:
		mutable int calls = 0;
		VectorN<Real, 1> operator()(const VectorN<Real, 1>& x) const override { ++calls; return x; }
	};
	class LargeRectangle : public IVectorFunctionNM<2, 3> {
	public:
		mutable int calls = 0;
		VectorN<Real, 3> operator()(const VectorN<Real, 2>& x) const override {
			++calls;
			return {x[0] / 2 + x[1] / 4, x[0] / 4 - x[1] / 2, -x[0] / 4 + x[1] / 4};
		}
	};
}

TEST_CASE("ODE and map Jacobians balance central steps and retain representable separation", "[jacobian-contract][float-final]") {
	const Vector<Real> point{Real(1.25), Real(-0.75)};
	for (bool ode : {false, true}) {
		for (Real step : {Real(0), Precision::CbrtEps / 2}) {
			CAPTURE(ode, step);
			std::vector<Vector<Real>> samples;
			const auto map = [&](const Vector<Real>& x) {
				samples.push_back(x);
				return Vector<Real>{std::sin(x[0]), std::cos(x[1])};
			};
			Matrix<Real> jacobian;
			if (ode) {
				const auto rhs = [&](Real t, const Vector<Real>& x, Vector<Real>& f) {
					CHECK(t == Real(0.5));
					f = map(x);
				};
				Derivation::CalcJacobianODE(rhs, Real(0.5), point, jacobian, step);
			} else {
				Derivation::CalcJacobianMap(map, point, jacobian, step);
			}
			REQUIRE(samples.size() == 4);
			REQUIRE(jacobian.rows() == 2);
			REQUIRE(jacobian.cols() == 2);
			for (int column = 0; column < 2; ++column) {
				const Real h = (step == 0 ? Precision::CbrtEps : step) * std::max(Real(1), std::abs(point[column]));
				CHECK(samples[2 * column][column] == point[column] + h);
				CHECK(samples[2 * column + 1][column] == point[column] - h);
				CHECK(samples[2 * column][1 - column] == point[1 - column]);
				CHECK(samples[2 * column + 1][1 - column] == point[1 - column]);
			}
			CHECK_THAT(jacobian(0, 0), WithinAbs(std::cos(point[0]), TOL(1e-9, Precision::RoundoffZeroThreshold)));
			CHECK_THAT(jacobian(1, 1), WithinAbs(-std::sin(point[1]), TOL(1e-9, Precision::RoundoffZeroThreshold)));
			CHECK(jacobian(0, 1) == 0);
			CHECK(jacobian(1, 0) == 0);
		}
		for (Real scale : {LargeJacobianMagnitude(), Real(1), SmallJacobianMagnitude()}) {
			const Vector<Real> x{scale};
			const auto identity = [](const Vector<Real>& v) { return v; };
			Matrix<Real> jacobian;
			if (ode)
				Derivation::CalcJacobianODE([&](Real, const Vector<Real>& v, Vector<Real>& f) { f = identity(v); }, 0, x, jacobian);
			else
				Derivation::CalcJacobianMap(identity, x, jacobian);
			CHECK(jacobian(0, 0) == 1);
		}
	}
}

TEST_CASE("Public full Jacobians preserve a large finite identity derivative", "[jacobian-contract-baseline][jacobian-contract]") {
	for (Real scale : {-LargeJacobianMagnitude(), LargeJacobianMagnitude(), Real(1), SmallJacobianMagnitude()}) {
		CAPTURE(scale);
		REQUIRE(std::isfinite(scale));
		REQUIRE(scale != Real(0));
		LargeIdentity function;
		int count = 5;
		const auto fixed = Derivation::CalcJacobian(function, VectorN<Real, 1>{scale}, Real(0), &count);
		CHECK_THAT(fixed(0, 0), WithinRel(Real(1), TOL(1e-9, 1e-3)));
		CHECK(function.calls == 4);
		CHECK(count == 9);
		int calls = 0;
		const auto dynamic = [&](const Vector<Real>& x) { ++calls; return x; };
		const Vector<Real> point{scale};
		const auto square = Derivation::CalcJacobian(dynamic, point);
		CHECK_THAT(square(0, 0), WithinRel(Real(1), TOL(1e-9, 1e-3)));
		CHECK(calls == 4);
		Matrix<Real> inPlace;
		Derivation::CalcJacobianInPlace(dynamic, point, inPlace);
		CHECK_THAT(inPlace(0, 0), WithinRel(Real(1), TOL(1e-9, 1e-3)));
		CHECK(calls == 8);
	}
}

TEST_CASE("Rectangular Jacobians keep coupled entries shapes and callback counts across scales", "[jacobian-contract-baseline][jacobian-contract]") {
	for (Real scale : {LargeJacobianMagnitude(), Real(1), SmallJacobianMagnitude()}) {
		CAPTURE(scale);
		LargeRectangle function;
		const Real step = scale / 1024;
		REQUIRE(std::isfinite(scale));
		REQUIRE(std::isfinite(step));
		REQUIRE(step > Real(0));
		REQUIRE(scale + step > scale);
		REQUIRE(scale - step < scale);
		const auto fixed = Derivation::CalcJacobian(function, VectorN<Real, 2>{scale, scale / 2}, step);
		CHECK(function.calls == 8);
		int calls = 0;
		const auto dynamic = [&](const Vector<Real>& x) {
			++calls;
			return Vector<Real>{x[0] / 2 + x[1] / 4, x[0] / 4 - x[1] / 2, -x[0] / 4 + x[1] / 4};
		};
		const auto matrix = Derivation::CalcJacobian(dynamic, Vector<Real>{scale, scale / 2}, 3, step);
		CHECK(matrix.rows() == 3);
		CHECK(matrix.cols() == 2);
		CHECK(calls == 8);
		const Real expected[3][2] = {{Real(0.5), Real(0.25)}, {Real(0.25), Real(-0.5)}, {Real(-0.25), Real(0.25)}};
		for (int row = 0; row < 3; ++row)
			for (int col = 0; col < 2; ++col) {
				CHECK_THAT(fixed(row, col), WithinAbs(expected[row][col], TOL(1e-9, 1e-3)));
				CHECK_THAT(matrix(row, col), WithinAbs(expected[row][col], TOL(1e-9, 1e-3)));
			}
	}
}

TEST_CASE("Jacobian difference quotients survive overflowing differences and denominators", "[jacobian-contract-baseline][jacobian-contract]") {
	const Real magnitude = LargeJacobianMagnitude(), step = magnitude / 4;
	REQUIRE(std::isfinite(magnitude));
	REQUIRE(std::isfinite(step));
	REQUIRE(step > Real(0));
	REQUIRE(std::isinf(magnitude - (-magnitude)));
	REQUIRE(std::isinf(12 * step));
	REQUIRE(std::isinf(2 * magnitude));
	int calls = 0;
	const auto function = [&](const Vector<Real>& x) { ++calls; return Vector<Real>{std::copysign(magnitude, x[0])}; };
	const Vector<Real> origin{0};
	CHECK_THAT(Derivation::NDer4PartialDyn(function, 0, 0, origin, step), WithinRel(Real(14) / 3, TOL(1e-10, 1e-3)));
	CHECK(calls == 4);
	CHECK_THAT(Derivation::NDer2PartialDyn(function, 0, 0, origin, magnitude), WithinRel(Real(1), TOL(1e-10, 1e-3)));
	CHECK(calls == 6);
	CHECK_THAT(Derivation::CalcJacobianOrder2(function, origin, magnitude)(0, 0), WithinRel(Real(1), TOL(1e-10, 1e-3)));
	CHECK(calls == 8);
	Matrix<Real> map;
	Derivation::CalcJacobianMap(function, origin, map, magnitude);
	CHECK_THAT(map(0, 0), WithinRel(Real(1), TOL(1e-10, 1e-3)));
	CHECK(calls == 10);
	const auto ode = [&](Real, const Vector<Real>& x, Vector<Real>& rhs) { rhs = function(x); };
	CHECK_THAT(Derivation::CalcJacobianODE(ode, 0, origin, magnitude)(0, 0), WithinRel(Real(1), TOL(1e-10, 1e-3)));
	CHECK(calls == 12);
}

TEST_CASE("Fourth-order Jacobians do not erase small differences behind equal huge outer samples", "[jacobian-contract-baseline][jacobian-contract]") {
	const Real step = SmallJacobianMagnitude(), magnitude = LargeJacobianMagnitude();
	REQUIRE(step > Real(0));
	REQUIRE(step < std::numeric_limits<Real>::min());
	REQUIRE(std::isfinite(magnitude));
	int magnitudeExponent;
	std::frexp(magnitude, &magnitudeExponent);
	REQUIRE(std::scalbn(step, -magnitudeExponent) == Real(0));
	int calls = 0;
	const auto function = [&](const Vector<Real>& x) {
		++calls;
		return Vector<Real>{std::abs(x[0]) == 2 * step ? magnitude : std::copysign(step, x[0])};
	};
	REQUIRE(function(Vector<Real>{2 * step})[0] == magnitude);
	REQUIRE(function(Vector<Real>{-2 * step})[0] == magnitude);
	REQUIRE(function(Vector<Real>{step})[0] == step);
	REQUIRE(function(Vector<Real>{-step})[0] == -step);
	calls = 0;
	CHECK_THAT(Derivation::NDer4PartialDyn(function, 0, 0, Vector<Real>{0}, step),
		WithinRel(Real(4) / 3, TOL(1e-10, 1e-3)));
	CHECK(calls == 4);
}

TEST_CASE("Jacobian stencils fail before callbacks on collapsed coordinates", "[jacobian-contract-baseline][jacobian-contract]") {
	int calls = 0;
	const auto function = [&](const Vector<Real>& x) { ++calls; return x; };
	const Vector<Real> point{Real(1e20)};
	CHECK_THROWS_AS(Derivation::CalcJacobian(function, point, Real(0.1)), NumericalMethodError);
	CHECK(calls == 0);
	calls = 0;
	CHECK_THROWS_AS(Derivation::NDer4PartialDyn(function, 0, 0, point, Real(0.1)), NumericalMethodError);
	CHECK(calls == 0);
	calls = 0;
	CHECK_THROWS_AS(Derivation::NDer2PartialDyn(function, 0, 0, point, Real(0.1)), NumericalMethodError);
	CHECK(calls == 0);
}

TEST_CASE("Public Jacobians reject nonfinite final entries rather than returning matrices", "[jacobian-contract-baseline][jacobian-contract]") {
	const Real magnitude = LargeJacobianMagnitude();
	REQUIRE(std::isfinite(magnitude));
	int calls = 0;
	const auto function = [&](const Vector<Real>& x) { ++calls; return Vector<Real>{std::copysign(magnitude, x[0])}; };
	CHECK_THROWS_AS(Derivation::CalcJacobian(function, Vector<Real>{0}, Real(0.01)), NumericalMethodError);
	CHECK(calls == 4);
	CHECK_THROWS_AS(Derivation::NDer2PartialDyn(function, 0, 0, Vector<Real>{0}, Real(0.01)), NumericalMethodError);
	CHECK(calls == 6);
}

TEST_CASE("Jacobian families reject invalid steps coordinates and collapsed outer stencils", "[jacobian-contract]") {
	int calls = 0;
	const auto function = [&](const Vector<Real>& x) { ++calls; return x; };
	const auto ode = [&](Real, const Vector<Real>& x, Vector<Real>& rhs) { ++calls; rhs = x; };
	const Vector<Real> point{1};
	Matrix<Real> matrix;
	for (Real h : {Real(-1), std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity()}) {
		LargeIdentity fixed;
		LargeRectangle rectangle;
		CHECK_THROWS_AS(Derivation::CalcJacobian(fixed, VectorN<Real, 1>{1}, h), ArgumentError);
		CHECK_THROWS_AS(Derivation::CalcJacobian(rectangle, VectorN<Real, 2>{1, 2}, h), ArgumentError);
		CHECK_THROWS_AS(Derivation::CalcJacobian(function, point, h), ArgumentError);
		CHECK_THROWS_AS(Derivation::NDer4PartialDyn(function, 0, 0, point, h), ArgumentError);
		CHECK_THROWS_AS(Derivation::NDer2PartialDyn(function, 0, 0, point, h), ArgumentError);
		CHECK_THROWS_AS(Derivation::CalcJacobianOrder2(function, point, h), ArgumentError);
		CHECK_THROWS_AS(Derivation::CalcJacobianODE(ode, 0, point, h), ArgumentError);
		CHECK_THROWS_AS(Derivation::CalcJacobianMap(function, point, matrix, h), ArgumentError);
		CHECK(fixed.calls == 0);
		CHECK(rectangle.calls == 0);
	}
	LargeIdentity fixed;
	CHECK_THROWS_AS(Derivation::CalcJacobian(fixed, VectorN<Real, 1>{Real(1e20)}, Real(0.1)), NumericalMethodError);
	CHECK(fixed.calls == 0);
	CHECK_THROWS_AS(Derivation::CalcJacobian(function, Vector<Real>{std::numeric_limits<Real>::max()}), NumericalMethodError);
	CHECK_THROWS_AS(Derivation::CalcJacobian(function, Vector<Real>{std::numeric_limits<Real>::infinity()}), ArgumentError);
	CHECK_THROWS_AS(Derivation::CalcJacobianODE(ode, 0, point, std::numeric_limits<Real>::denorm_min()), NumericalMethodError);
	CHECK_THROWS_AS(Derivation::CalcJacobianMap(function, point, matrix, std::numeric_limits<Real>::denorm_min()), NumericalMethodError);
	CHECK(calls == 0);
}

TEST_CASE("Jacobian callback dimensions and nonfinite values fail explicitly", "[jacobian-contract]") {
	const Vector<Real> point{1};
	Matrix<Real> matrix;
	for (bool nonfinite : {false, true}) {
		const auto function = [&](const Vector<Real>&) {
			return nonfinite ? Vector<Real>{std::numeric_limits<Real>::quiet_NaN()} : Vector<Real>{1, 2};
		};
		const auto ode = [&](Real, const Vector<Real>& x, Vector<Real>& rhs) { rhs = function(x); };
		if (nonfinite) {
			CHECK_THROWS_AS(Derivation::CalcJacobian(function, point), NumericalMethodError);
			CHECK_THROWS_AS(Derivation::CalcJacobianOrder2(function, point), NumericalMethodError);
			CHECK_THROWS_AS(Derivation::NDer4PartialDyn(function, 0, 0, point), NumericalMethodError);
			CHECK_THROWS_AS(Derivation::NDer2PartialDyn(function, 0, 0, point), NumericalMethodError);
			CHECK_THROWS_AS(Derivation::CalcJacobianODE(ode, 0, point), NumericalMethodError);
			CHECK_THROWS_AS(Derivation::CalcJacobianMap(function, point, matrix), NumericalMethodError);
		}
		else {
			CHECK_THROWS_AS(Derivation::CalcJacobian(function, point), VectorDimensionError);
			CHECK_THROWS_AS(Derivation::CalcJacobianOrder2(function, point), VectorDimensionError);
			CHECK_THROWS_AS(Derivation::NDer4PartialDyn(function, 2, 0, point), VectorDimensionError);
			CHECK_THROWS_AS(Derivation::NDer2PartialDyn(function, 2, 0, point), VectorDimensionError);
			CHECK_THROWS_AS(Derivation::CalcJacobianODE(ode, 0, point), VectorDimensionError);
			CHECK_THROWS_AS(Derivation::CalcJacobianMap(function, point, matrix), VectorDimensionError);
		}
	}
	int calls = 0;
	const auto changing = [&](const Vector<Real>&) { return ++calls == 1 ? Vector<Real>{1, 2} : Vector<Real>{1}; };
	CHECK_THROWS_AS(Derivation::NDer4PartialDyn(changing, 0, 0, point), VectorDimensionError);
	const auto identity = [](const Vector<Real>& x) { return x; };
	CHECK_THROWS_AS(Derivation::NDer4PartialDyn(identity, -1, 0, point), ArgumentError);
	CHECK_THROWS_AS(Derivation::NDer2PartialDyn(identity, 0, 1, point), ArgumentError);
	CHECK_THROWS_AS(Derivation::CalcJacobianOrder2(identity, Vector<Real>{}), ArgumentError);
}

TEST_CASE("Jacobian callback exceptions and attempted-evaluation counters are preserved", "[jacobian-contract]") {
	int count = 7;
	const auto throwing = [](const Vector<Real>&) -> Vector<Real> { throw std::logic_error("callback witness"); };
	CHECK_THROWS_AS(Derivation::CalcJacobian(throwing, Vector<Real>{1}, Real(0), &count), std::logic_error);
	CHECK(count == 8);
	CHECK_THROWS_AS(Derivation::NDer4PartialDyn(throwing, 0, 0, Vector<Real>{1}), std::logic_error);
	CHECK_THROWS_AS(Derivation::NDer2PartialDyn(throwing, 0, 0, Vector<Real>{1}), std::logic_error);
}

TEST_CASE("Jacobian component steps and relative ODE map steps keep legacy sampling and costs", "[jacobian-contract]") {
	const Vector<Real> point{2, 3};
	std::vector<Vector<Real>> samples;
	const auto identity = [&](const Vector<Real>& x) { samples.push_back(x); return x; };
	const Real automatic = std::pow(std::numeric_limits<Real>::epsilon(), Real(1) / 5);
	Derivation::NDer4PartialDyn(identity, 0, 0, point);
	REQUIRE(samples.size() == 4);
	CHECK(samples[0][0] == point[0] + 2 * automatic);
	CHECK(samples[1][0] == point[0] + automatic);
	CHECK(samples[2][0] == point[0] - automatic);
	CHECK(samples[3][0] == point[0] - 2 * automatic);
	samples.clear();
	const auto ode = [&](Real, const Vector<Real>& x, Vector<Real>& rhs) { rhs = identity(x); };
	Derivation::CalcJacobianODE(ode, 0, point, Real(0.125));
	REQUIRE(samples.size() == 4);
	CHECK(samples[0][0] == Real(2.25));
	CHECK(samples[2][1] == Real(3.375));
	samples.clear();
	Matrix<Real> map;
	Derivation::CalcJacobianMap(identity, point, map, Real(0.125));
	REQUIRE(samples.size() == 4);
	CHECK(samples[0][0] == Real(2.25));
	CHECK(samples[2][1] == Real(3.375));
	samples.clear();
	const auto order2 = Derivation::CalcJacobianOrder2(identity, point);
	CHECK(samples.size() == 8);
	CHECK_THAT(order2(0, 0), WithinAbs(Real(1), TOL(1e-9, 1e-3)));
	CHECK_THAT(order2(1, 1), WithinAbs(Real(1), TOL(1e-9, 1e-3)));
}

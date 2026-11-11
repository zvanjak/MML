#include <catch2/catch_all.hpp>
#include "../../TestPrecision.h"

#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/core/Integration/PathIntegration.h>
#include <mml/core/Integration/SurfaceIntegration.h>
#endif

#include <array>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

using namespace MML;

namespace MML::Tests::GeometricIntegrationContracts
{
	struct Surface : IParametricSurfaceRect<3> {
		VectorN<Real, 3> a{1, 0, 0}, b{0, 1, 0};
		Real lo = 0, hi = 1;
		mutable int calls = 0;
		bool poison = false;
		Real getMinU() const override { return lo; }
		Real getMaxU() const override { return hi; }
		Real getMinW() const override { return lo; }
		Real getMaxW() const override { return hi; }
		VectorN<Real, 3> operator()(Real u, Real v) const override {
			++calls;
			if (!(lo <= u && u <= hi && lo <= v && v <= hi))
				throw std::out_of_range("Surface stencil left its declared domain");
			if (poison) return {std::numeric_limits<Real>::quiet_NaN(), 0, 0};
			return a * u + b * v;
		}
	};

	struct CountingField : IVectorFunction<3> {
		VectorN<Real, 3> value{0, 0, 1};
		mutable int calls = 0;
		bool fail = false;
		VectorN<Real, 3> operator()(const VectorN<Real, 3>&) const override {
			++calls;
			if (fail) throw std::runtime_error("field callback");
			return value;
		}
	};

	struct Line : IParametricCurve<3>, ICurveDerivativeProvider<3> {
		mutable int calls = 0, derivative_calls = 0;
		Real getMinT() const override { return 0; }
		Real getMaxT() const override { return 1; }
		VectorN<Real, 3> operator()(Real t) const override { ++calls; return {t, 0, 0}; }
		VectorN<Real, 3> EvaluateDerivative(Real, int order) const override {
			++derivative_calls;
			return order == 1 ? VectorN<Real, 3>{1, 0, 0} : VectorN<Real, 3>{0, 0, 0};
		}
	};

	struct Density : IRealFunction {
		mutable int calls = 0;
		Real operator()(Real t) const override { ++calls; return -1 - t; }
	};
	struct Scalar : IScalarFunction<3> {
		mutable int calls = 0;
		Real operator()(const VectorN<Real, 3>& p) const override { ++calls; return -1 - p[0]; }
	};

	void CheckFixed(const IntegrationResult& result, int calls) {
		CHECK(result.status == AlgorithmStatus::Completed);
		CHECK(result.IsUsable());
		CHECK_FALSE(result.converged);
		CHECK_FALSE(result.error_estimated);
		CHECK_FALSE(HasEstimate(result.error_estimate));
		CHECK(result.iterations == calls);
		CHECK(result.function_evaluations == calls);
		CHECK(std::isfinite(result.value));
	}

	void CheckFiniteInputs(const Surface& surface, const CountingField& field) {
		REQUIRE(std::isfinite(surface.lo));
		REQUIRE(std::isfinite(surface.hi));
		REQUIRE(surface.lo < surface.hi);
		for (int i = 0; i < 3; ++i) {
			REQUIRE(std::isfinite(surface.a[i]));
			REQUIRE(std::isfinite(surface.b[i]));
			REQUIRE(std::isfinite(field.value[i]));
		}
	}
}

using namespace MML::Tests::GeometricIntegrationContracts;

TEST_CASE("Small Jacobians retain their physical flux", "[geometric-integration]") {
	Surface surface;
	surface.a = {Real(1e-8), 0, 0};
	surface.b = {0, Real(1e-8), 0};
	CountingField field;
	field.value = {0, 0, Real(1e16)};
	CheckFiniteInputs(surface, field);
	const auto result = SurfaceIntegration::SurfaceIntegralDetailed(field, surface, 0, 1, 0, 1, 3, 4);
	CHECK_THAT(result.value, Catch::Matchers::WithinAbs(Real(1), TOL(1e-8, 32 * Precision::Eps)));
	CHECK(field.calls == 12);
	CHECK(surface.calls == 60);
	CheckFixed(result, field.calls);
}

TEST_CASE("Surface stencils stay in narrow declared parameter domains", "[geometric-integration]") {
	const int exponent = std::min(100, (std::numeric_limits<Real>::max_exponent - 8) / 2);
	for (const Real width : {Real(1e-7), std::ldexp(Real(1), -exponent)}) {
		CAPTURE(width);
		Surface surface;
		surface.hi = width;
		CountingField field;
		field.value = {0, 0, (1 / width) / width};
		CheckFiniteInputs(surface, field);
		REQUIRE(field.value[2] > Real(0));
		const auto result = SurfaceIntegration::SurfaceIntegralDetailed(field, surface, 0, width, 0, width, 2, 3);
		CHECK_THAT(result.value, Catch::Matchers::WithinAbs(Real(1), Real(1e-8)));
		CHECK(field.calls == 6);
		CHECK(surface.calls == 30);
		CheckFixed(result, field.calls);
	}
}

TEST_CASE("Weighted surface flux avoids cross product range loss", "[geometric-integration]") {
	const int tangentExponent = std::numeric_limits<Real>::max_exponent == 1024
		? 600 : 5 * std::numeric_limits<Real>::max_exponent / 8;
	const int domainExponent = std::numeric_limits<Real>::max_exponent == 1024
		? 100 : 5 * std::numeric_limits<Real>::max_exponent / 32;
	for (int sign : {-1, 1}) {
		CAPTURE(sign, tangentExponent, domainExponent);
		Surface surface;
		const Real scale = std::ldexp(Real(1), tangentExponent * sign);
		surface.a = {scale, 0, 0};
		surface.b = {0, scale, 0};
		surface.hi = std::ldexp(Real(1), -domainExponent * sign);
		CountingField field;
		field.value = {0, 0, std::ldexp(Real(1), -2 * (tangentExponent - domainExponent) * sign)};
		CheckFiniteInputs(surface, field);
		REQUIRE(scale > Real(0));
		REQUIRE(field.value[2] > Real(0));
		if (sign < 0)
			REQUIRE(scale * scale == Real(0));
		else
			REQUIRE(std::isinf(scale * scale));
		const auto result = SurfaceIntegration::SurfaceIntegralDetailed(
			field, surface, 0, surface.hi, 0, surface.hi, 2, 2);
		CHECK_THAT(result.value, Catch::Matchers::WithinAbs(Real(1), Real(1e-8)));
		CHECK(field.calls == 4);
		CHECK(surface.calls == 20);
		CheckFixed(result, field.calls);
	}
}

TEST_CASE("General oriented flux obeys independent determinant scaling", "[geometric-integration]") {
	const int exponent = std::min(400, (std::numeric_limits<Real>::max_exponent - 8) / 2);
	for (int p : {-exponent, 0, exponent}) for (int winding : {-1, 1}) {
		CAPTURE(p, winding);
		Surface surface;
		surface.a = VectorN<Real, 3>{1, 2, -1} * std::ldexp(Real(winding), p);
		surface.b = VectorN<Real, 3>{-2, 1, 3} * std::ldexp(Real(1), p);
		CountingField field;
		field.value = VectorN<Real, 3>{2, -3, 5} * std::ldexp(Real(1), -2 * p);
		CheckFiniteInputs(surface, field);
		REQUIRE(surface.a[0] != Real(0));
		REQUIRE(surface.b[1] != Real(0));
		REQUIRE(field.value[0] > Real(0));
		const auto result = SurfaceIntegration::SurfaceIntegralDetailed(field, surface, 0, 1, 0, 1, 2, 3);
		CHECK_THAT(result.value, Catch::Matchers::WithinAbs(Real(42 * winding), TOL(1e-7, 32 * Precision::Eps * 42)));
		CHECK(field.calls == 6);
		CHECK(surface.calls == 30);
		CheckFixed(result, field.calls);
	}
}

TEST_CASE("Anisotropic flux remains invariant on observable coordinate axes", "[geometric-integration]") {
	const int exponent = std::min(400, (std::numeric_limits<Real>::max_exponent - 8) / 2);
	for (int p : {-exponent, 0, exponent}) for (int q : {-exponent, 0, exponent}) for (int axis = 0; axis < 3; ++axis) {
		CAPTURE(p, q, axis);
		Surface surface;
		surface.a = {0, 0, 0};
		surface.b = {0, 0, 0};
		surface.a[axis] = std::ldexp(Real(1), p);
		surface.b[(axis + 1) % 3] = std::ldexp(Real(1), q);
		CountingField field;
		field.value = {0, 0, 0};
		field.value[(axis + 2) % 3] = std::ldexp(Real(1), -p - q);
		CheckFiniteInputs(surface, field);
		REQUIRE(surface.a[axis] > Real(0));
		REQUIRE(surface.b[(axis + 1) % 3] > Real(0));
		REQUIRE(field.value[(axis + 2) % 3] > Real(0));
		const auto result = SurfaceIntegration::SurfaceIntegralDetailed(field, surface, 0, 1, 0, 1, 2, 3);
		CHECK_THAT(result.value, Catch::Matchers::WithinAbs(Real(1), Real(1e-8)));
		CHECK(field.calls == 6);
		CHECK(surface.calls == 30);
		CheckFixed(result, field.calls);
	}
}

TEST_CASE("Representable endpoint midpoints use in-domain one-sided stencils", "[geometric-integration]") {
	Surface surface;
	surface.lo = 1;
	surface.hi = std::nextafter(std::nextafter(Real(1), Real(2)), Real(2));
	const Real width = surface.hi - surface.lo;
	CountingField field;
	field.value = {0, 0, (1 / width) / width};
	const auto result = SurfaceIntegration::SurfaceIntegralDetailed(
		field, surface, surface.lo, surface.hi, surface.lo, surface.hi, 2, 2);
	CHECK_THAT(result.value, Catch::Matchers::WithinAbs(Real(1), Real(1e-8)));
	CHECK(surface.calls == 20);
	CheckFixed(result, field.calls);
}

TEST_CASE("Finite surface limits need not have a representable full width", "[geometric-integration]") {
	Surface surface;
	surface.lo = -std::numeric_limits<Real>::max();
	surface.hi = std::numeric_limits<Real>::max();
	const int tangentExponent = std::numeric_limits<Real>::max_exponent == 1024
		? 600 : 5 * std::numeric_limits<Real>::max_exponent / 8;
	const Real scale = std::ldexp(Real(1), -tangentExponent);
	surface.a = {scale, 0, 0};
	surface.b = {0, scale, 0};
	CountingField field;
	field.value = {0, 0, std::ldexp(Real(1), 2 * (tangentExponent - std::numeric_limits<Real>::max_exponent))};
	CheckFiniteInputs(surface, field);
	REQUIRE(scale > Real(0));
	REQUIRE(field.value[2] > Real(0));
	REQUIRE(std::isinf(surface.hi - surface.lo));
	const long double normalizedLimit = std::ldexp(static_cast<long double>(surface.hi),
		-std::numeric_limits<Real>::max_exponent);
	const Real expected = TOL(4, static_cast<Real>(4 * normalizedLimit * normalizedLimit));
	const auto result = SurfaceIntegration::SurfaceIntegralDetailed(
		field, surface, surface.lo, surface.hi, surface.lo, surface.hi, 1, 1);
	CHECK_THAT(result.value, Catch::Matchers::WithinAbs(expected, Real(1e-8)));
	CHECK(field.calls == 1);
	CHECK(surface.calls == 5);
	CheckFixed(result, field.calls);
}

TEST_CASE("Triangle mesh flux uses oriented area without a normal cutoff", "[geometric-integration]") {
	const int extreme = std::min(400, (std::numeric_limits<Real>::max_exponent - 8) / 2);
	for (int exponent : {-extreme, -40, 0, extreme}) for (int winding : {-1, 1}) {
		CAPTURE(exponent, winding);
		const Real scale = std::ldexp(Real(1), exponent);
		REQUIRE(std::isfinite(scale));
		REQUIRE(scale > Real(0));
		const Triangle3D triangle(Pnt3Cart(0, 0, 0), Pnt3Cart(winding * scale, 0, 0), Pnt3Cart(0, scale, 0));
		CountingField field;
		field.value = {0, 0, std::ldexp(Real(1), -2 * exponent)};
		REQUIRE(std::isfinite(field.value[2]));
		REQUIRE(field.value[2] > Real(0));
		const auto result = SurfaceIntegration::SurfaceIntegralDetailed(field, triangle, Real(1e-8), 0);
		CHECK_THAT(result.value, Catch::Matchers::WithinAbs(Real(winding) / 2, Real(1e-10)));
		CHECK(result.converged);
		CHECK(result.error_estimated);
		CHECK(result.function_evaluations == field.calls);
		CHECK(field.calls == 5);
	}
}

TEST_CASE("Exactly degenerate surfaces retain finite zero flux and actual field counts", "[geometric-integration]") {
	Surface surface;
	surface.b = surface.a;
	CountingField field;
	const auto fixed = SurfaceIntegration::SurfaceIntegralDetailed(field, surface, 0, 1, 0, 1, 2, 3);
	CHECK(fixed.value == 0);
	CHECK(field.calls == 6);
	CheckFixed(fixed, field.calls);
	field.calls = 0;
	const Triangle3D triangle(Pnt3Cart(0, 0, 0), Pnt3Cart(1, 0, 0), Pnt3Cart(2, 0, 0));
	const auto mesh = SurfaceIntegration::SurfaceIntegralDetailed(field, triangle, Real(1e-8), 0);
	CHECK(mesh.value == 0);
	CHECK(mesh.converged);
	CHECK(mesh.function_evaluations == field.calls);
	CHECK(field.calls == 5);
}

TEST_CASE("Parametric surface invalid controls fail before callbacks", "[geometric-integration]") {
	Surface surface;
	CountingField field;
	field.fail = true;
	const Real nan = std::numeric_limits<Real>::quiet_NaN();
	const Real inf = std::numeric_limits<Real>::infinity();
	for (Real invalid : {nan, inf, -inf}) {
		CHECK_THROWS_AS(SurfaceIntegration::SurfaceIntegralDetailed(field, surface, invalid, 1, 0, 1), NumericInputError);
		CHECK_THROWS_AS(SurfaceIntegration::SurfaceIntegralDetailed(field, surface, 0, 1, 0, invalid), NumericInputError);
	}
	CHECK_THROWS_AS(SurfaceIntegration::SurfaceIntegralDetailed(field, surface, 1, 0, 0, 1), ArgumentError);
	CHECK_THROWS_AS(SurfaceIntegration::SurfaceIntegralDetailed(field, surface, 0, 1, 0, 1, 0, 2), ArgumentError);
	CHECK_THROWS_AS(SurfaceIntegration::SurfaceIntegralDetailed(field, surface, 0, 1, 0, 1,
		std::numeric_limits<int>::max(), 2), ResourceLimitError);
	CHECK(surface.calls == 0);
	CHECK(field.calls == 0);
}

TEST_CASE("Surface integration validates declared bounds before sampling", "[geometric-integration]") {
	Surface surface;
	CountingField field;
	CHECK_THROWS_AS(SurfaceIntegration::SurfaceIntegralDetailed(field, surface, -1, 1, 0, 1), ArgumentError);
	CHECK_THROWS_AS(SurfaceIntegration::SurfaceIntegralDetailed(field, surface, 0, 1, 0, 2), ArgumentError);
	surface.lo = std::numeric_limits<Real>::quiet_NaN();
	CHECK_THROWS_AS(SurfaceIntegration::SurfaceIntegralDetailed(field, surface, 0, 1, 0, 1), ArgumentError);
	CHECK(surface.calls == 0);
	CHECK(field.calls == 0);
}

TEST_CASE("Nonfinite surface data cannot produce Completed results", "[geometric-integration]") {
	Surface surface;
	CountingField field;
	surface.poison = true;
	CHECK_THROWS_AS(SurfaceIntegration::SurfaceIntegralDetailed(field, surface, 0, 1, 0, 1, 1, 1), NumericalMethodError);
	surface.poison = false;
	field.value[2] = std::numeric_limits<Real>::infinity();
	CHECK_THROWS_AS(SurfaceIntegration::SurfaceIntegralDetailed(field, surface, 0, 1, 0, 1, 1, 1), NumericalMethodError);
	field.fail = true;
	CHECK_THROWS_AS(SurfaceIntegration::SurfaceIntegralDetailed(field, surface, 0, 1, 0, 1, 1, 1), std::runtime_error);
}

TEST_CASE("Insufficient distinct surface stencil nodes fail explicitly", "[geometric-integration]") {
	Surface surface;
	surface.lo = 1;
	surface.hi = std::nextafter(Real(1), Real(2));
	CountingField field;
	CHECK_THROWS_AS(SurfaceIntegration::SurfaceIntegralDetailed(
		field, surface, surface.lo, surface.hi, surface.lo, surface.hi, 1, 1), NumericalMethodError);
	CHECK(field.calls == 0);
}

TEST_CASE("Unsigned path measures normalize limits without changing scalar signs", "[geometric-integration]") {
	Line line;
	Density density;
	Scalar scalar;
	CountingField field;
	field.value = {1, 0, 0};
	for (IntegrationMethod method : {TRAP, SIMPSON, ROMBERG, GAUSS10, GAUSS10KRONROD21}) {
		line.calls = line.derivative_calls = field.calls = scalar.calls = density.calls = 0;
		const auto length = PathIntegration::ParametricCurveLengthDetailed(line, 1, 0, method);
		CHECK(length.function_evaluations == line.derivative_calls);
		CHECK(line.calls == 0);
		line.derivative_calls = 0;
		const auto mass = PathIntegration::ParametricCurveMassDetailed(line, density, 1, 0, method);
		CHECK(mass.function_evaluations == line.derivative_calls);
		CHECK(mass.function_evaluations == density.calls);
		line.calls = line.derivative_calls = 0;
		const auto integral = PathIntegration::LineIntegralDetailed(scalar, line, 1, 0, method);
		CHECK(integral.function_evaluations == scalar.calls);
		CHECK(integral.function_evaluations == line.calls);
		line.calls = line.derivative_calls = 0;
		const auto work = PathIntegration::LineIntegralDetailed(field, line, 1, 0, method);
		CHECK(work.function_evaluations == field.calls);
		CHECK(work.function_evaluations == line.calls);
		CHECK(work.function_evaluations == line.derivative_calls);
		CHECK_THAT(length.value, Catch::Matchers::WithinAbs(Real(1), Real(1e-8)));
		CHECK_THAT(mass.value, Catch::Matchers::WithinAbs(Real(-1.5), Real(1e-8)));
		CHECK_THAT(integral.value, Catch::Matchers::WithinAbs(Real(-1.5), Real(1e-8)));
		CHECK_THAT(work.value, Catch::Matchers::WithinAbs(Real(-1), Real(1e-8)));
		CHECK(length.IsUsable());
		CHECK(mass.IsUsable());
		CHECK(integral.IsUsable());
		CHECK(work.IsUsable());
	}
	CHECK_THAT(PathIntegration::ParametricCurveLength(line, 1, 0), Catch::Matchers::WithinAbs(Real(1), Real(1e-8)));
	CHECK_THAT(PathIntegration::ParametricCurveMass(line, density, 1, 0), Catch::Matchers::WithinAbs(Real(-1.5), Real(1e-8)));
	CHECK_THAT(PathIntegration::LineIntegral(scalar, line, 1, 0), Catch::Matchers::WithinAbs(Real(-1.5), Real(1e-8)));
	CHECK_THAT(PathIntegration::LineIntegralDetailed(scalar, line, 1, 0).value,
		Catch::Matchers::WithinAbs(Real(-1.5), Real(1e-8)));
	const Real nan = std::numeric_limits<Real>::quiet_NaN();
	CHECK_THROWS_AS(PathIntegration::ParametricCurveLengthDetailed(line, nan, 1), NumericInputError);
	CHECK_THROWS_AS(PathIntegration::LineIntegralDetailed(scalar, line, 1, nan), NumericInputError);
}

TEST_CASE("Zero path intervals retain quadrature callback accounting", "[geometric-integration]") {
	Line line;
	Density density;
	Scalar scalar;
	for (IntegrationMethod method : {TRAP, SIMPSON, ROMBERG, GAUSS10, GAUSS10KRONROD21}) {
		line.calls = line.derivative_calls = density.calls = scalar.calls = 0;
		const auto length = PathIntegration::ParametricCurveLengthDetailed(line, Real(.5), Real(.5), method);
		CHECK(length.function_evaluations == line.derivative_calls);
		CHECK(line.calls == 0);
		line.derivative_calls = 0;
		const auto mass = PathIntegration::ParametricCurveMassDetailed(line, density, Real(.5), Real(.5), method);
		CHECK(mass.function_evaluations == density.calls);
		CHECK(mass.function_evaluations == line.derivative_calls);
		line.calls = line.derivative_calls = 0;
		const auto integral = PathIntegration::LineIntegralDetailed(scalar, line, Real(.5), Real(.5), method);
		CHECK(integral.function_evaluations == scalar.calls);
		CHECK(integral.function_evaluations == line.calls);
		CHECK(integral.function_evaluations == line.derivative_calls);
		for (const auto& result : {length, mass, integral}) {
			CHECK(result.value == 0);
			CHECK(result.IsUsable());
		}
	}
}

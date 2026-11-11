#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "../../TestPrecision.h"

#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/base/Function.h>
#include <mml/core/Derivation.h>
#endif

#include <functional>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace MML;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

TEST_CASE("Complex step quotient retains finite extremes and rejects true overflow", "[step-contract][precision]") {
	const Real maximum = std::numeric_limits<Real>::max();
	CHECK(Derivation::ComplexStep([](Complex z) { return z; }, Real(0), Real(0.5)) == Real(1));
	CHECK(Derivation::ComplexStep([=](Complex) { return Complex(0, maximum / 2); }, Real(0), Real(0.5)) == maximum);
	CHECK_THROWS_AS(Derivation::ComplexStep([=](Complex) { return Complex(0, maximum); }, Real(0), Real(0.5)),
		NumericalMethodError);
	CHECK(Derivation::ComplexStep([](Complex z) { return z; }, Real(0), std::numeric_limits<Real>::denorm_min()) == Real(1));
}

namespace {
	template<typename Invoke>
	void CheckInvalidSteps(Invoke invoke, int& calls) {
		for (Real step : {Real(0), Real(-1), std::numeric_limits<Real>::infinity(), std::numeric_limits<Real>::quiet_NaN()})
			for (bool estimate : {false, true}) {
				calls = 0;
				Real error = 17;
				CHECK_THROWS_AS(invoke(step, estimate ? &error : nullptr), NumericInputError);
				CHECK(calls == 0);
				CHECK(error == 17);
			}
	}

	class StepCurve : public IParametricCurve<2> {
	public:
		mutable int calls = 0;
		VectorN<Real, 2> operator()(Real t) const override { ++calls; return {t, 2 * t}; }
		Real getMinT() const override { return -10; }
		Real getMaxT() const override { return 10; }
	};
	class StepSurface : public IParametricSurfaceRect<2> {
	public:
		mutable int calls = 0;
		VectorN<Real, 2> operator()(Real u, Real w) const override { ++calls; return {u + w, u - w}; }
		Real getMinU() const override { return -10; }
		Real getMaxU() const override { return 10; }
		Real getMinW() const override { return -10; }
		Real getMaxW() const override { return 10; }
	};
	class StepTensor2 : public ITensorField2<2> {
	public:
		mutable int calls = 0;
		StepTensor2() : ITensorField2<2>(0, 2) {}
		Real component(int, int, const VectorN<Real, 2>& x) const override { ++calls; return x[0]; }
		Tensor2<2> operator()(const VectorN<Real, 2>&) const override { return Tensor2<2>(0, 2); }
	};
	class StepTensor3 : public ITensorField3<2> {
	public:
		mutable int calls = 0;
		Real component(int, int, int, const VectorN<Real, 2>& x) const override { ++calls; return x[0]; }
		Tensor3<2> operator()(const VectorN<Real, 2>&) const override { return Tensor3<2>(0, 3); }
	};
	class StepTensor4 : public ITensorField4<2> {
	public:
		mutable int calls = 0;
		Real component(int, int, int, int, const VectorN<Real, 2>& x) const override { ++calls; return x[0]; }
		Tensor4<2> operator()(const VectorN<Real, 2>&) const override { return Tensor4<2>(0, 4); }
	};
	class StepComplex : public IComplexFunction {
	public:
		mutable int calls = 0;
		Complex operator()(Complex z) const override { ++calls; return z * z; }
	};
}

TEST_CASE("Scalar and vector partials reject invalid explicit steps for every stencil order", "[step-contract-baseline][step-contract]") {
	int calls = 0;
	ScalarFunctionFromStdFunc<2> scalar([&](const VectorN<Real, 2>& x) { ++calls; return x[0] + x[1]; });
	VectorFunctionFromStdFunc<2> vector([&](const VectorN<Real, 2>& x) { ++calls; return x; });
	const VectorN<Real, 2> point{1, 2};
	const std::vector<std::function<void(Real, Real*)>> derivatives = {
		[&](Real h, Real* e) { Derivation::NDer1Partial(scalar, 0, point, h, e); },
		[&](Real h, Real* e) { Derivation::NDer2Partial(scalar, 0, point, h, e); },
		[&](Real h, Real* e) { Derivation::NDer4Partial(scalar, 0, point, h, e); },
		[&](Real h, Real* e) { Derivation::NDer6Partial(scalar, 0, point, h, e); },
		[&](Real h, Real* e) { Derivation::NDer8Partial(scalar, 0, point, h, e); },
		[&](Real h, Real* e) { Derivation::NDer1Partial(vector, 0, 0, point, h, e); },
		[&](Real h, Real* e) { Derivation::NDer2Partial(vector, 0, 0, point, h, e); },
		[&](Real h, Real* e) { Derivation::NDer4Partial(vector, 0, 0, point, h, e); },
		[&](Real h, Real* e) { Derivation::NDer6Partial(vector, 0, 0, point, h, e); },
		[&](Real h, Real* e) { Derivation::NDer8Partial(vector, 0, 0, point, h, e); }
	};
	for (const auto& derivative : derivatives) CheckInvalidSteps(derivative, calls);
}

TEST_CASE("Curve and surface first partials reject invalid steps before callbacks", "[step-contract-baseline][step-contract]") {
	StepCurve curve;
	const std::vector<std::function<void(Real, Real*)>> curves = {
		[&](Real h, Real* e) { Derivation::NDer1(curve, 1, h, e); },
		[&](Real h, Real* e) { Derivation::NDer2(curve, 1, h, e); },
		[&](Real h, Real* e) { Derivation::NDer4(curve, 1, h, e); },
		[&](Real h, Real* e) { Derivation::NDer6(curve, 1, h, e); },
		[&](Real h, Real* e) { Derivation::NDer8(curve, 1, h, e); }
	};
	for (const auto& derivative : curves) CheckInvalidSteps(derivative, curve.calls);
	StepSurface surface;
	const std::vector<std::function<void(Real, Real*)>> surfaces = {
		[&](Real h, Real* e) { Derivation::NDer1_u(surface, 1, 2, h, e); },
		[&](Real h, Real* e) { Derivation::NDer1_w(surface, 1, 2, h, e); },
		[&](Real h, Real* e) { Derivation::NDer2_u(surface, 1, 2, h, e); },
		[&](Real h, Real* e) { Derivation::NDer2_w(surface, 1, 2, h, e); }
	};
	for (const auto& derivative : surfaces) CheckInvalidSteps(derivative, surface.calls);
}

TEST_CASE("Tensor component steps are validated at ranks two three and four", "[step-contract-baseline][step-contract]") {
	const VectorN<Real, 2> point{1, 2};
	StepTensor2 rank2;
	StepTensor3 rank3;
	StepTensor4 rank4;
	CheckInvalidSteps([&](Real h, Real* e) { Derivation::NDer1Partial(rank2, 0, 0, 0, point, h, e); }, rank2.calls);
	CheckInvalidSteps([&](Real h, Real* e) { Derivation::NDer2Partial(rank2, 0, 0, 0, point, h, e); }, rank2.calls);
	CheckInvalidSteps([&](Real h, Real* e) { Derivation::NDer4Partial(rank2, 0, 0, 0, point, h, e); }, rank2.calls);
	CheckInvalidSteps([&](Real h, Real* e) { Derivation::NDer1Partial(rank3, 0, 0, 0, 0, point, h, e); }, rank3.calls);
	CheckInvalidSteps([&](Real h, Real* e) { Derivation::NDer2Partial(rank3, 0, 0, 0, 0, point, h, e); }, rank3.calls);
	CheckInvalidSteps([&](Real h, Real* e) { Derivation::NDer4Partial(rank3, 0, 0, 0, 0, point, h, e); }, rank3.calls);
	CheckInvalidSteps([&](Real h, Real* e) { Derivation::NDer1Partial(rank4, 0, 0, 0, 0, 0, point, h, e); }, rank4.calls);
	CheckInvalidSteps([&](Real h, Real* e) { Derivation::NDer2Partial(rank4, 0, 0, 0, 0, 0, point, h, e); }, rank4.calls);
	CheckInvalidSteps([&](Real h, Real* e) { Derivation::NDer4Partial(rank4, 0, 0, 0, 0, 0, point, h, e); }, rank4.calls);
}

TEST_CASE("Complex-step methods reject invalid explicit steps before callbacks", "[step-contract-baseline][step-contract]") {
	int calls = 0;
	const auto function = [&](Complex z) { ++calls; return z * z; };
	CheckInvalidSteps([&](Real h, Real*) { Derivation::ComplexStep(function, 1, h); }, calls);
	CheckInvalidSteps([&](Real h, Real*) { Derivation::ComplexStep2(function, 1, h); }, calls);
}

TEST_CASE("Aggregate partials validate steps without changing error outputs", "[step-contract]") {
	int calls = 0;
	ScalarFunctionFromStdFunc<2> scalar([&](const VectorN<Real, 2>& x) { ++calls; return x[0] + x[1]; });
	VectorFunctionFromStdFunc<2> vector([&](const VectorN<Real, 2>& x) { ++calls; return x; });
	const VectorN<Real, 2> point{1, 2};
	VectorN<Real, 2> errors{17, 17};
	MatrixNM<Real, 2, 2> matrixErrors;
	for (int i = 0; i < 2; ++i)
		for (int j = 0; j < 2; ++j) matrixErrors(i, j) = 17;
	const std::vector<std::function<void(Real, Real*)>> aggregates = {
		[&](Real h, Real* e) { Derivation::NDer1PartialByAll(scalar, point, h, e ? &errors : nullptr); },
		[&](Real h, Real* e) { Derivation::NDer2PartialByAll(scalar, point, h, e ? &errors : nullptr); },
		[&](Real h, Real* e) { Derivation::NDer4PartialByAll(scalar, point, h, e ? &errors : nullptr); },
		[&](Real h, Real* e) { Derivation::NDer6PartialByAll(scalar, point, h, e ? &errors : nullptr); },
		[&](Real h, Real* e) { Derivation::NDer8PartialByAll(scalar, point, h, e ? &errors : nullptr); },
		[&](Real h, Real* e) { Derivation::NDer1PartialByAll(vector, 0, point, h, e ? &errors : nullptr); },
		[&](Real h, Real* e) { Derivation::NDer2PartialByAll(vector, 0, point, h, e ? &errors : nullptr); },
		[&](Real h, Real* e) { Derivation::NDer4PartialByAll(vector, 0, point, h, e ? &errors : nullptr); },
		[&](Real h, Real* e) { Derivation::NDer6PartialByAll(vector, 0, point, h, e ? &errors : nullptr); },
		[&](Real h, Real* e) { Derivation::NDer8PartialByAll(vector, 0, point, h, e ? &errors : nullptr); },
		[&](Real h, Real* e) { Derivation::NDer1PartialAllByAll(vector, point, h, e ? &matrixErrors : nullptr); },
		[&](Real h, Real* e) { Derivation::NDer2PartialAllByAll(vector, point, h, e ? &matrixErrors : nullptr); },
		[&](Real h, Real* e) { Derivation::NDer4PartialAllByAll(vector, point, h, e ? &matrixErrors : nullptr); },
		[&](Real h, Real* e) { Derivation::NDer6PartialAllByAll(vector, point, h, e ? &matrixErrors : nullptr); },
		[&](Real h, Real* e) { Derivation::NDer8PartialAllByAll(vector, point, h, e ? &matrixErrors : nullptr); }
	};
	for (const auto& derivative : aggregates) {
		CheckInvalidSteps(derivative, calls);
		CHECK(errors[0] == 17);
		CHECK(errors[1] == 17);
		for (int i = 0; i < 2; ++i)
			for (int j = 0; j < 2; ++j) CHECK(matrixErrors(i, j) == 17);
	}
}

TEST_CASE("Real complex and shifted one-sided wrappers retain the same invalid-step policy", "[step-contract]") {
	int calls = 0;
	RealFunctionFromStdFunc function([&](Real x) { ++calls; return x * x; });
	const std::vector<std::function<void(Real, Real*)>> real = {
		[&](Real h, Real* e) { Derivation::NDer1(function, 1, h, e); },
		[&](Real h, Real* e) { Derivation::NDer2(function, 1, h, e); },
		[&](Real h, Real* e) { Derivation::NDer4(function, 1, h, e); },
		[&](Real h, Real* e) { Derivation::NDer6(function, 1, h, e); },
		[&](Real h, Real* e) { Derivation::NDer8(function, 1, h, e); },
		[&](Real h, Real* e) { Derivation::NDer1Left(function, 1, h, e); },
		[&](Real h, Real* e) { Derivation::NDer1Right(function, 1, h, e); },
		[&](Real h, Real* e) { Derivation::NDer2Left(function, 1, h, e); },
		[&](Real h, Real* e) { Derivation::NDer2Right(function, 1, h, e); },
		[&](Real h, Real* e) { Derivation::NDer4Left(function, 1, h, e); },
		[&](Real h, Real* e) { Derivation::NDer4Right(function, 1, h, e); },
		[&](Real h, Real* e) { Derivation::NDer6Left(function, 1, h, e); },
		[&](Real h, Real* e) { Derivation::NDer6Right(function, 1, h, e); },
		[&](Real h, Real* e) { Derivation::NDer8Left(function, 1, h, e); },
		[&](Real h, Real* e) { Derivation::NDer8Right(function, 1, h, e); }
	};
	for (const auto& derivative : real) CheckInvalidSteps(derivative, calls);
	StepComplex complex;
	CheckInvalidSteps([&](Real h, Real* e) { Derivation::NDer1Complex(complex, Complex(1, 2), h, e); }, complex.calls);
	CheckInvalidSteps([&](Real h, Real* e) { Derivation::NDer2Complex(complex, Complex(1, 2), h, e); }, complex.calls);
	CheckInvalidSteps([&](Real h, Real* e) { Derivation::NDer4Complex(complex, Complex(1, 2), h, e); }, complex.calls);
	CheckInvalidSteps([&](Real h, Real* e) { Derivation::NDer6Complex(complex, Complex(1, 2), h, e); }, complex.calls);
	DerivativeConfig config;
	config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
	CHECK(Derivation::NDer4Detailed(function, 1, Real(0), config).status == AlgorithmStatus::InvalidInput);
	CHECK(Derivation::NDer4ComplexDetailed(complex, Complex(1, 2), Real(0), config).status == AlgorithmStatus::InvalidInput);
	CHECK(calls == 0);
	CHECK(complex.calls == 0);
	config.step = 0;
	CHECK(Derivation::NDer4Detailed(function, 1, config).IsSuccess());
	CHECK(Derivation::NDer4ComplexDetailed(complex, Complex(1, 2), config).IsSuccess());
}

TEST_CASE("Complex-step checks samples squared steps and final division while retaining tiny first steps", "[step-contract]") {
	int calls = 0;
	const auto square = [&](Complex z) { ++calls; return z * z; };
	for (Real x : {std::numeric_limits<Real>::infinity(), std::numeric_limits<Real>::quiet_NaN()}) {
		CHECK_THROWS_AS(Derivation::ComplexStep(square, x), NumericInputError);
		CHECK_THROWS_AS(Derivation::ComplexStep2(square, x), NumericInputError);
	}
	for (Real h : {std::numeric_limits<Real>::denorm_min(), std::numeric_limits<Real>::max()})
		CHECK_THROWS_AS(Derivation::ComplexStep2(square, 1, h), NumericalMethodError);
	CHECK(calls == 0);
	const auto identity = [](Complex z) { return z; };
	CHECK(Derivation::ComplexStep(identity, 1, std::numeric_limits<Real>::denorm_min()) == 1);
	CHECK_THAT(Derivation::ComplexStep(square, 1), WithinRel(Real(2), TOL(1e-12, 1e-4)));
	CHECK_THAT(Derivation::ComplexStep2(square, 0), WithinRel(Real(2), TOL(1e-12, 1e-4)));
	CHECK(calls == 3);
	for (Complex invalid : {Complex(std::numeric_limits<Real>::infinity(), 1), Complex(1, std::numeric_limits<Real>::quiet_NaN())}) {
		const auto nonfinite = [&](Complex) { return invalid; };
		CHECK_THROWS_AS(Derivation::ComplexStep(nonfinite, 1), NumericalMethodError);
		CHECK_THROWS_AS(Derivation::ComplexStep2(nonfinite, 1), NumericalMethodError);
	}
	const auto overflowing = [](Complex) { return Complex(0, std::numeric_limits<Real>::max()); };
	CHECK_THROWS_AS(Derivation::ComplexStep(overflowing, 1, std::numeric_limits<Real>::min()), NumericalMethodError);
	const auto throwing = [](Complex) -> Complex { throw std::logic_error("callback witness"); };
	CHECK_THROWS_AS(Derivation::ComplexStep(throwing, 1), std::logic_error);
	CHECK_THROWS_AS(Derivation::ComplexStep2(throwing, 1), std::logic_error);
}

TEST_CASE("Legacy one-sided functions differentiate shifted centers rather than the named endpoint", "[step-contract]") {
	RealFunction square([](Real x) { return x * x; });
	const Real x = 1, h = Real(0.125);
	const std::vector<std::pair<int, std::function<Real(Real*)>>> shifted = {
		{-3, [&](Real* e) { return Derivation::NDer2Left(square, x, h, e); }},
		{ 3, [&](Real* e) { return Derivation::NDer2Right(square, x, h, e); }},
		{-4, [&](Real* e) { return Derivation::NDer4Left(square, x, h, e); }},
		{ 4, [&](Real* e) { return Derivation::NDer4Right(square, x, h, e); }},
		{-5, [&](Real* e) { return Derivation::NDer6Left(square, x, h, e); }},
		{ 5, [&](Real* e) { return Derivation::NDer6Right(square, x, h, e); }},
		{-6, [&](Real* e) { return Derivation::NDer8Left(square, x, h, e); }},
		{ 6, [&](Real* e) { return Derivation::NDer8Right(square, x, h, e); }}
	};
	for (const auto& [shift, derivative] : shifted)
		for (bool estimate : {false, true}) {
			Real error = 17;
			CHECK_THAT(derivative(estimate ? &error : nullptr), WithinAbs(2 * (x + shift * h), TOL(1e-10, 1e-3)));
			CHECK(error >= 0);
		}
	CHECK_THAT(Derivation::NDer4Left(square, 1, Real(0.01)), WithinAbs(Real(1.92), TOL(1e-10, 1e-3)));
	CHECK_THAT(Derivation::NDer1Left(square, x, h), WithinAbs(2 * x - h, TOL(1e-10, 1e-3)));
	CHECK_THAT(Derivation::NDer1Right(square, x, h), WithinAbs(2 * x + h, TOL(1e-10, 1e-3)));
}

TEST_CASE("First-order Right error estimation can cross the nominal domain boundary", "[step-contract]") {
	std::vector<Real> samples;
	RealFunctionFromStdFunc bounded([&](Real x) {
		samples.push_back(x);
		if (x < 0) throw DomainError("negative domain");
		return std::sqrt(x);
	});
	CHECK(std::isfinite(Derivation::NDer1Right(bounded, 0, Real(0.125))));
	REQUIRE(samples.size() == 2);
	CHECK(samples[0] == Real(0.125));
	CHECK(samples[1] == 0);
	samples.clear();
	Real error = 17;
	CHECK_THROWS_AS(Derivation::NDer1Right(bounded, 0, Real(0.125), &error), DomainError);
	REQUIRE(samples.size() == 3);
	CHECK(samples[2] == Real(-0.125));
	CHECK(error == 17);
}

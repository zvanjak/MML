#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "../../TestPrecision.h"

#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/core/Derivation/DerivationRealFunction.h>
#include <mml/core/Derivation/DerivationTensorField.h>
#endif

#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace MML;
using Catch::Matchers::WithinAbs;

namespace {
	class RefinementCallbackError : public std::runtime_error {
	public:
		RefinementCallbackError() : std::runtime_error("Richardson callback failure") {}
	};

	class RefinementFunction : public IRealFunction {
	public:
		std::function<Real(Real)> value = [](Real x) { return x; };
		mutable std::vector<Real> arguments;
		Real operator()(Real x) const override {
			arguments.push_back(x);
			return value(x);
		}
	};

	Real Refine(int order, const IRealFunction& function, Real x, Real h,
		Real* error = nullptr, int depth = 10, Real shrink = Real(1.4)) {
		switch (order) {
		case 1: return Derivation::NDerRichardson(function, x, h, error, depth, shrink);
		case 2: return Derivation::NSecDerRichardson(function, x, h, error, depth, shrink);
		default: return Derivation::NThirdDerRichardson(function, x, h, error, depth, shrink);
		}
	}

	class DefaultTensor3 : public ITensorField3<2> {
	public:
		mutable int calls = 0;
		Real component(int i, int j, int k, const VectorN<Real, 2>& x) const override {
			++calls;
			return x[0] + 2 * x[1] + i + 2 * j + 3 * k;
		}
		Tensor3<2> operator()(const VectorN<Real, 2>&) const override { return Tensor3<2>(0, 3); }
	};

	class DefaultTensor4 : public ITensorField4<2> {
	public:
		mutable int calls = 0;
		Real component(int i, int j, int k, int l, const VectorN<Real, 2>& x) const override {
			++calls;
			return x[0] + 2 * x[1] + i + 2 * j + 3 * k + 4 * l;
		}
		Tensor4<2> operator()(const VectorN<Real, 2>&) const override { return Tensor4<2>(0, 4); }
	};
}

TEST_CASE("Richardson rejects invalid step shrink and center before evaluating callbacks",
	"[derivative-refinement-baseline][derivative-refinement]") {
	for (int order : {1, 2, 3}) {
		INFO("derivative order " << order);
		RefinementFunction function;
		for (Real step : {Real(0), Real(-1), std::numeric_limits<Real>::quiet_NaN(),
			std::numeric_limits<Real>::infinity()}) {
			Real error = -1;
			CHECK_THROWS_AS(Refine(order, function, 0, step, &error), NumericInputError);
			CHECK(error == Real(-1));
			CHECK(function.arguments.empty());
			function.arguments.clear();
		}
		for (Real shrink : {Real(0), Real(1), Real(-2), std::numeric_limits<Real>::quiet_NaN(),
			std::numeric_limits<Real>::infinity(), std::numeric_limits<Real>::max()}) {
			CHECK_THROWS_AS(Refine(order, function, 0, Real(0.125), nullptr, 3, shrink), NumericInputError);
			CHECK(function.arguments.empty());
			function.arguments.clear();
		}
		for (Real center : {std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity()}) {
			CHECK_THROWS_AS(Refine(order, function, center, Real(0.125)), NumericInputError);
			CHECK(function.arguments.empty());
			function.arguments.clear();
		}
	}
}

TEST_CASE("Richardson rejects nonpositive tableau budgets without allocation or callbacks", "[derivative-refinement]") {
	for (int order : {1, 2, 3}) {
		RefinementFunction function;
		for (int depth : {0, -1, std::numeric_limits<int>::min()}) {
			CHECK_THROWS_AS(Refine(order, function, 0, Real(0.125), nullptr, depth), NumericInputError);
			CHECK(function.arguments.empty());
		}
	}
}

TEST_CASE("Richardson fails explicitly on initially collapsed or overflowing stencil coordinates",
	"[derivative-refinement-baseline][derivative-refinement]") {
	for (int order : {1, 2, 3}) {
		for (Real sign : {Real(-1), Real(1)}) {
			RefinementFunction function;
			Real error = -1;
			CHECK_THROWS_AS(Refine(order, function, sign * Real(1e20), Real(0.1), &error), NumericalMethodError);
			CHECK(error == Real(-1));
			CHECK(function.arguments.empty());
			function.arguments.clear();
			CHECK_THROWS_AS(Refine(order, function, sign * std::numeric_limits<Real>::max(),
				std::numeric_limits<Real>::max() / 4), NumericalMethodError);
			CHECK(function.arguments.empty());
		}
	}
}

TEST_CASE("Richardson validates a refined stencil before evaluating its collapsed samples",
	"[derivative-refinement-baseline][derivative-refinement]") {
	for (int order : {1, 2, 3}) {
		RefinementFunction function;
		const Real center = std::ldexp(Real(1), std::numeric_limits<Real>::digits - 1);
		const Real gap = std::nextafter(center, std::numeric_limits<Real>::infinity()) - center;
		Real error = -1;
		CHECK_THROWS_AS(Refine(order, function, center, gap, &error, 3, 4), NumericalMethodError);
		CHECK(error == Real(-1));
		CHECK(function.arguments.size() == static_cast<size_t>(order == 1 ? 2 : order == 2 ? 3 : 4));
		for (Real argument : function.arguments) CHECK(std::isfinite(argument));
	}
}

TEST_CASE("Richardson rejects nonfinite callback data and arithmetic rather than returning a result",
	"[derivative-refinement-baseline][derivative-refinement]") {
	for (int order : {1, 2, 3}) {
		for (Real invalid : {std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity(),
			-std::numeric_limits<Real>::infinity()}) {
			RefinementFunction function;
			function.value = [=](Real) { return invalid; };
			Real error = -1;
			CHECK_THROWS_AS(Refine(order, function, 0, Real(0.125), &error), NumericalMethodError);
			CHECK(error == Real(-1));
			CHECK(function.arguments.size() == 1);
		}
		RefinementFunction overflowing;
		overflowing.value = [](Real x) { return std::copysign(std::numeric_limits<Real>::max(), x); };
		CHECK_THROWS_AS(Refine(order, overflowing, 0, Real(0.125)), NumericalMethodError);
	}
}

TEST_CASE("Richardson one-level estimates retain the no-error-information sentinel", "[derivative-refinement]") {
	RefinementFunction function;
	for (int order : {1, 2, 3}) {
		function.arguments.clear();
		function.value = [=](Real x) {
			return order == 1 ? x : order == 2 ? x * x : x * x * x;
		};
		Real error = 0;
		const Real derivative = Refine(order, function, 0, Real(0.125), &error, 1);
		CHECK(derivative == (order == 1 ? Real(1) : order == 2 ? Real(2) : Real(6)));
		CHECK(error == std::numeric_limits<Real>::max());
		CHECK(function.arguments.size() == static_cast<size_t>(order == 1 ? 2 : order == 2 ? 3 : 4));
	}
}

TEST_CASE("Richardson analytic derivatives remain accurate with finite nonnegative estimates", "[derivative-refinement]") {
	for (Real x : {Real(-1), Real(0), Real(1)}) {
		RefinementFunction exponential;
		exponential.value = [](Real t) { return std::exp(t); };
		for (int order : {1, 2, 3}) {
			Real error = -1;
			CHECK_THAT(Refine(order, exponential, x, Real(0.25), &error),
				WithinAbs(std::exp(x), TOL(1e-7, 1e-3)));
			CHECK(std::isfinite(error));
			CHECK(error >= 0);
			CHECK(error < std::numeric_limits<Real>::max());
		}
	}
	const Real scale = TOL(1, 32);
	const auto polynomial = [scale](Real x) { return x * x * x * x * x / scale; };
	CHECK_THAT(Derivation::NDerRichardson(polynomial, Real(2)), WithinAbs(Real(80) / scale, TOL(1e-7, 1e-3)));
	CHECK_THAT(Derivation::NSecDerRichardson(polynomial, Real(2)), WithinAbs(Real(160) / scale, TOL(1e-7, 1e-3)));
	CHECK_THAT(Derivation::NThirdDerRichardson(polynomial, Real(2)), WithinAbs(Real(240) / scale, TOL(1e-7, 1e-3)));
}

TEST_CASE("Richardson validates every callback sample including later refinement levels", "[derivative-refinement]") {
	for (int order : {1, 2, 3}) {
		const int initialSamples = order == 1 ? 2 : order == 2 ? 3 : 4;
		const int nextSamples = order == 3 ? 4 : 2;
		for (int failAt = 1; failAt <= initialSamples + nextSamples; ++failAt) {
			for (Real invalid : {std::numeric_limits<Real>::quiet_NaN(),
				std::numeric_limits<Real>::infinity(), -std::numeric_limits<Real>::infinity()}) {
				RefinementFunction function;
				int calls = 0;
				function.value = [&](Real x) { return ++calls == failAt ? invalid : std::exp(x); };
				Real error = -1;
				CHECK_THROWS_AS(Refine(order, function, 0, Real(0.25), &error), NumericalMethodError);
				CHECK(calls == failAt);
				CHECK(error == Real(-1));
			}
		}
	}
}

TEST_CASE("Richardson rejects step underflow and third-derivative outer sample collapse", "[derivative-refinement]") {
	for (int order : {1, 2, 3}) {
		RefinementFunction function;
		Real error = -1;
		CHECK_THROWS_AS(Refine(order, function, 0, std::numeric_limits<Real>::denorm_min(), &error, 2, 4),
			NumericalMethodError);
		CHECK(error == Real(-1));
		CHECK(function.arguments.size() == static_cast<size_t>(order == 1 ? 2 : order == 2 ? 3 : 4));
	}
	RefinementFunction function;
	const Real center = std::ldexp(Real(1), std::numeric_limits<Real>::digits - 1);
	const Real gap = std::nextafter(center, std::numeric_limits<Real>::infinity()) - center;
	CHECK_THROWS_AS(Refine(3, function, center, gap * Real(0.6)), NumericalMethodError);
	CHECK(function.arguments.empty());
}

TEST_CASE("Richardson rejects extrapolation overflow from finite base estimates", "[derivative-refinement]") {
	RefinementFunction function;
	const Real large = std::numeric_limits<Real>::max() / 4;
	function.value = [=](Real x) {
		return std::copysign(large, std::abs(x) == Real(1) ? x : -x);
	};
	Real error = -1;
	CHECK_THROWS_AS(Refine(1, function, 0, 1, &error, 2, std::nextafter(Real(1), Real(2))), NumericalMethodError);
	CHECK(function.arguments.size() == 4);
	CHECK(error == Real(-1));
}

TEST_CASE("Richardson preserves callback exception types and leaves error outputs untouched", "[derivative-refinement]") {
	for (int order : {1, 2, 3}) {
		for (int failAt : {1, 4}) {
			RefinementFunction function;
			int calls = 0;
			function.value = [&](Real x) {
				if (++calls == failAt) throw RefinementCallbackError();
				return std::exp(x);
			};
			Real error = -1;
			CHECK_THROWS_AS(Refine(order, function, 0, Real(0.25), &error), RefinementCallbackError);
			CHECK(calls == failAt);
			CHECK(error == Real(-1));
		}
	}
}

TEST_CASE("Callable Richardson overloads and adaptive aliases share refinement safety", "[derivative-refinement]") {
	const auto identity = [](Real x) { return x; };
	Real error = -1;
	CHECK_THROWS_AS(Derivation::NDerRichardson(identity, Real(1e20), Real(0.1), &error), NumericalMethodError);
	CHECK_THROWS_AS(Derivation::NSecDerRichardson(identity, Real(1e20), Real(0.1), &error), NumericalMethodError);
	CHECK_THROWS_AS(Derivation::NThirdDerRichardson(identity, Real(1e20), Real(0.1), &error), NumericalMethodError);
	CHECK_THROWS_AS(Derivation::NDerAdaptive(identity, Real(1e20), &error), NumericalMethodError);
	CHECK_THROWS_AS(Derivation::NSecDerAdaptive(identity, Real(1e20), &error), NumericalMethodError);
	CHECK_THROWS_AS(Derivation::NThirdDerAdaptive(identity, Real(1e20), &error), NumericalMethodError);
	CHECK(error == Real(-1));
	const auto exponential = [](Real x) { return std::exp(x); };
	CHECK_THAT(Derivation::NDerAdaptive(exponential, Real(0), &error), WithinAbs(Real(1), TOL(1e-7, 1e-3)));
	CHECK_THAT(Derivation::NSecDerAdaptive(exponential, Real(0), &error), WithinAbs(Real(1), TOL(1e-7, 1e-3)));
	CHECK_THAT(Derivation::NThirdDerAdaptive(exponential, Real(0), &error), WithinAbs(Real(1), TOL(1e-7, 1e-3)));
}

TEST_CASE("Default rank-three and rank-four tensor partial overloads instantiate and match explicit steps",
	"[derivative-refinement][tensor-default-overloads]") {
	DefaultTensor3 rank3;
	DefaultTensor4 rank4;
	const VectorN<Real, 2> point{Real(1.25), Real(-2)};
	for (int axis : {0, 1}) {
		auto check = [&](auto& field, auto implicit, auto explicitStep, Real defaultStep, int samples, int errorSamples) {
			for (bool estimate : {false, true}) {
				Real error = -1, controlError = -1;
				field.calls = 0;
				const Real value = implicit(estimate ? &error : nullptr);
				CHECK(field.calls == (estimate ? errorSamples : samples));
				const Real control = explicitStep(Derivation::ScaleStep(defaultStep, point[axis]),
					estimate ? &controlError : nullptr);
				CHECK(value == control);
				CHECK_THAT(value, WithinAbs(Real(axis + 1), TOL(1e-7, 1e-3)));
				CHECK(error == controlError);
				if (estimate) CHECK(error >= 0);
			}
		};
		check(rank3, [&](Real* error) { return Derivation::NDer1Partial(rank3, 1, 0, 1, axis, point, error); },
			[&](Real h, Real* error) { return Derivation::NDer1Partial(rank3, 1, 0, 1, axis, point, h, error); },
			Derivation::NDer1_h, 2, 3);
		check(rank3, [&](Real* error) { return Derivation::NDer2Partial(rank3, 1, 0, 1, axis, point, error); },
			[&](Real h, Real* error) { return Derivation::NDer2Partial(rank3, 1, 0, 1, axis, point, h, error); },
			Derivation::NDer2_h, 2, 4);
		check(rank3, [&](Real* error) { return Derivation::NDer4Partial(rank3, 1, 0, 1, axis, point, error); },
			[&](Real h, Real* error) { return Derivation::NDer4Partial(rank3, 1, 0, 1, axis, point, h, error); },
			Derivation::NDer4_h, 4, 6);
		check(rank4, [&](Real* error) { return Derivation::NDer1Partial(rank4, 1, 0, 1, 0, axis, point, error); },
			[&](Real h, Real* error) { return Derivation::NDer1Partial(rank4, 1, 0, 1, 0, axis, point, h, error); },
			Derivation::NDer1_h, 2, 3);
		check(rank4, [&](Real* error) { return Derivation::NDer2Partial(rank4, 1, 0, 1, 0, axis, point, error); },
			[&](Real h, Real* error) { return Derivation::NDer2Partial(rank4, 1, 0, 1, 0, axis, point, h, error); },
			Derivation::NDer2_h, 2, 4);
		check(rank4, [&](Real* error) { return Derivation::NDer4Partial(rank4, 1, 0, 1, 0, axis, point, error); },
			[&](Real h, Real* error) { return Derivation::NDer4Partial(rank4, 1, 0, 1, 0, axis, point, h, error); },
			Derivation::NDer4_h, 4, 6);
		CHECK_THAT(Derivation::NDer1Partial(rank3, 1, 0, 1, axis, point), WithinAbs(Real(axis + 1), TOL(1e-7, 1e-3)));
		CHECK_THAT(Derivation::NDer2Partial(rank3, 1, 0, 1, axis, point), WithinAbs(Real(axis + 1), TOL(1e-7, 1e-3)));
		CHECK_THAT(Derivation::NDer4Partial(rank3, 1, 0, 1, axis, point), WithinAbs(Real(axis + 1), TOL(1e-7, 1e-3)));
		CHECK_THAT(Derivation::NDer1Partial(rank4, 1, 0, 1, 0, axis, point), WithinAbs(Real(axis + 1), TOL(1e-7, 1e-3)));
		CHECK_THAT(Derivation::NDer2Partial(rank4, 1, 0, 1, 0, axis, point), WithinAbs(Real(axis + 1), TOL(1e-7, 1e-3)));
		CHECK_THAT(Derivation::NDer4Partial(rank4, 1, 0, 1, 0, axis, point), WithinAbs(Real(axis + 1), TOL(1e-7, 1e-3)));
	}
}

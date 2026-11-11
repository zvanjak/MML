#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "../../TestPrecision.h"

#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/core/Derivation/DerivationRealFunction.h>
#include <mml/core/Derivation/DerivationScalarFunction.h>
#include <mml/core/Derivation/DerivationParametricCurve.h>
#include <mml/core/Derivation/DerivationParametricSurface.h>
#endif

#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

using namespace MML;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {
	struct ErrorWitness {
		int derivative, accuracy, degree;
		Real leading;
		int valueSamples, errorSamples;
	};
	constexpr ErrorWitness witnesses[] = {
		{2, 2, 4, 2, 3, 5}, {2, 4, 6, -8, 5, 7},
		{3, 2, 5, 30, 4, 6}, {3, 4, 7, -294, 6, 8}
	};
	using DirectDerivative = Real (*)(const IRealFunction&, Real, Real, Real*);
	const DirectDerivative directDerivatives[] = {
		Derivation::NSecDer2, Derivation::NSecDer4, Derivation::NThirdDer2, Derivation::NThirdDer4
	};

	class ErrorPolynomial : public IRealFunction {
	public:
		int degree = 4;
		mutable int calls = 0;
		Real operator()(Real x) const override { ++calls; return std::pow(x, MML::Real(degree)); }
	};

	class ErrorCurve : public IParametricCurve<3> {
	public:
		int degree = 4;
		mutable int calls = 0;
		Real getMinT() const override { return -4; }
		Real getMaxT() const override { return 4; }
		VectorN<Real, 3> operator()(Real t) const override {
			++calls;
			const Real value = std::pow(t, MML::Real(degree));
			return {value, 2 * value, -2 * value};
		}
	};

	class ErrorScalar : public IScalarFunction<2> {
	public:
		std::function<Real(const VectorN<Real, 2>&)> value;
		mutable std::vector<VectorN<Real, 2>> arguments;
		Real operator()(const VectorN<Real, 2>& point) const override {
			arguments.push_back(point);
			return value(point);
		}
	};

	class ErrorSurface : public IParametricSurfaceRect<3> {
	public:
		std::function<Real(Real, Real)> value;
		mutable std::vector<std::pair<Real, Real>> arguments;
		Real getMinU() const override { return -4; }
		Real getMaxU() const override { return 4; }
		Real getMinW() const override { return -4; }
		Real getMaxW() const override { return 4; }
		VectorN<Real, 3> operator()(Real u, Real w) const override {
			arguments.emplace_back(u, w);
			const Real result = value(u, w);
			return {result, 2 * result, -2 * result};
		}
	};

	DerivativeResult<Real> RealEstimate(int method, const IRealFunction& function, Real h, bool estimate = true) {
		DerivativeConfig config;
		config.estimate_error = estimate;
		switch (method) {
		case 0: return Derivation::NSecDer2Detailed(function, 0, h, config);
		case 1: return Derivation::NSecDer4Detailed(function, 0, h, config);
		case 2: return Derivation::NThirdDer2Detailed(function, 0, h, config);
		default: return Derivation::NThirdDer4Detailed(function, 0, h, config);
		}
	}

	VectorN<Real, 3> CurveEstimate(int method, const IParametricCurve<3>& curve, Real h, Real* error) {
		switch (method) {
		case 0: return Derivation::NSecDer2(curve, 0, h, error);
		case 1: return Derivation::NSecDer4(curve, 0, h, error);
		case 2: return Derivation::NThirdDer2(curve, 0, h, error);
		default: return Derivation::NThirdDer4(curve, 0, h, error);
		}
	}
}

TEST_CASE("Real higher derivative estimates have the correct leading coefficient and step power", "[higher-error-baseline][higher-error]") {
	for (int method = 0; method < 4; ++method) {
		const auto witness = witnesses[method];
		ErrorPolynomial function;
		function.degree = witness.degree;
		Real previousError = 0, previousStep = 0;
		for (Real h : {Real(0.25), Real(0.125), Real(0.0625), Real(0.1)}) {
			function.calls = 0;
			const auto result = RealEstimate(method, function, h);
			const Real expected = witness.leading * std::pow(h, MML::Real(witness.accuracy));
			CHECK_THAT(result.value, WithinRel(expected, TOL(1e-10, 1e-3)));
			CHECK_THAT(result.error, WithinRel(std::abs(expected), TOL(1e-10, 1e-3)));
			CHECK(result.function_evaluations == witness.errorSamples);
			CHECK(function.calls == witness.errorSamples);
			if (previousError > 0)
				CHECK_THAT(previousError / result.error, WithinRel(std::pow(previousStep / h, Real(witness.accuracy)), TOL(1e-10, 1e-3)));
			previousError = result.error;
			previousStep = h;
			function.calls = 0;
			const auto valueOnly = RealEstimate(method, function, h, false);
			CHECK(valueOnly.value == result.value);
			CHECK_FALSE(HasEstimate(valueOnly.error));
			CHECK(function.calls == witness.valueSamples);
			CHECK(valueOnly.function_evaluations == witness.valueSamples);
			Real directError = -1;
			CHECK(directDerivatives[method](function, 0, h, &directError) == result.value);
			CHECK(directError == result.error);
		}
	}
}

TEST_CASE("Curve higher derivative errors match the vector norm of polynomial truncation", "[higher-error-baseline][higher-error]") {
	for (int method = 0; method < 4; ++method) {
		const auto witness = witnesses[method];
		ErrorCurve curve;
		curve.degree = witness.degree;
		for (Real h : {Real(0.25), Real(0.125), Real(0.1)}) {
			curve.calls = 0;
			Real error = -1;
			const auto result = CurveEstimate(method, curve, h, &error);
			const Real expected = witness.leading * std::pow(h, MML::Real(witness.accuracy));
			CHECK_THAT(result[0], WithinRel(expected, TOL(1e-10, 1e-3)));
			CHECK_THAT(error, WithinRel(3 * std::abs(expected), TOL(1e-10, 1e-3)));
			CHECK(curve.calls == witness.errorSamples);
			curve.calls = 0;
			CHECK(CurveEstimate(method, curve, h, nullptr) == result);
			CHECK(curve.calls == witness.valueSamples);
		}
	}
}

TEST_CASE("Pure scalar second partial estimates share correct truncation scaling on either axis", "[higher-error-baseline][higher-error]") {
	for (int axis : {0, 1}) {
		for (int accuracy : {2, 4}) {
			ErrorScalar function;
			function.value = [=](const auto& point) { return std::pow(point[axis], Real(accuracy + 2)); };
			for (Real h : {Real(0.25), Real(0.125), Real(0.1)}) {
				function.arguments.clear();
				Real error = -1;
				const auto evaluate = [&](Real* out) {
					return accuracy == 2 ? Derivation::NSecDer2Partial(function, axis, axis, VectorN<Real, 2>{0, 0}, h, out)
						: Derivation::NSecDer4Partial(function, axis, axis, VectorN<Real, 2>{0, 0}, h, out);
				};
				const Real value = evaluate(&error);
				CHECK_THAT(error, WithinRel(std::abs(value), TOL(1e-10, 1e-3)));
				CHECK(function.arguments.size() == static_cast<size_t>(accuracy == 2 ? 5 : 7));
				function.arguments.clear();
				CHECK(evaluate(nullptr) == value);
				CHECK(function.arguments.size() == static_cast<size_t>(accuracy == 2 ? 3 : 5));
			}
		}
	}
}

TEST_CASE("Surface pure second errors estimate truncation instead of half the derivative", "[higher-error-baseline][higher-error]") {
	for (int axis : {0, 1}) {
		ErrorSurface surface;
		surface.value = [=](Real u, Real w) { return std::pow(axis == 0 ? u : w, Real(4)); };
		for (Real h : {Real(0.25), Real(0.125), Real(0.1)}) {
			surface.arguments.clear();
			Real error = -1;
			const auto evaluate = [&](Real* out) {
				return axis == 0 ? Derivation::NDer2_uu(surface, 0, 0, h, out) : Derivation::NDer2_ww(surface, 0, 0, h, out);
			};
			const auto value = evaluate(&error);
			CHECK_THAT(error, WithinRel(value.NormL2(), TOL(1e-10, 1e-3)));
			CHECK(surface.arguments.size() == 5);
			surface.arguments.clear();
			CHECK(evaluate(nullptr) == value);
			CHECK(surface.arguments.size() == 3);
		}
	}
}

TEST_CASE("Mixed scalar estimates include truncation with four optional in-region samples", "[higher-error-baseline][higher-error]") {
	for (int axis : {0, 1}) {
		for (int accuracy : {2, 4}) {
			ErrorScalar function;
			function.value = [=](const auto& point) { return std::pow(point[axis], Real(accuracy + 1)) * point[1 - axis]; };
			const VectorN<Real, 2> point = axis == 0 ? VectorN<Real, 2>{0, 1} : VectorN<Real, 2>{1, 0};
			for (Real h : {Real(0.25), Real(0.125), Real(0.1)}) {
				function.arguments.clear();
				Real error = -1;
				const auto evaluate = [&](Real* out) {
					return accuracy == 2 ? Derivation::NSecDer2Partial(function, 0, 1, point, h, out)
						: Derivation::NSecDer4Partial(function, 0, 1, point, h, out);
				};
				const Real value = evaluate(&error);
				CHECK_THAT(error, WithinRel(std::abs(value), TOL(1e-10, 1e-3)));
				const int count = accuracy == 2 ? 4 : 8;
				REQUIRE(function.arguments.size() == static_cast<size_t>(count + 4));
				const Real reach = accuracy == 2 ? h : 2 * h;
				for (const auto& sample : function.arguments)
					for (int coordinate : {0, 1}) {
						CHECK(sample[coordinate] >= point[coordinate] - reach);
						CHECK(sample[coordinate] <= point[coordinate] + reach);
					}
				function.arguments.clear();
				CHECK(evaluate(nullptr) == value);
				CHECK(function.arguments.size() == static_cast<size_t>(count));
			}
		}
	}
}

TEST_CASE("Anisotropic surface mixed estimates include truncation without widening the stencil", "[higher-error-baseline][higher-error]") {
	ErrorSurface surface;
	surface.value = [](Real u, Real w) { return u * u * u * w + u * w * w * w; };
	for (Real hu : {Real(0.25), Real(0.125)}) {
		const Real hw = hu / 2;
		Real error = -1;
		surface.arguments.clear();
		const auto value = Derivation::NDer2_uw(surface, 0, 0, hu, hw, &error);
		CHECK_THAT(error, WithinRel(3 * (hu * hu + hw * hw), TOL(1e-10, 1e-3)));
		REQUIRE(surface.arguments.size() == 8);
		for (const auto& [u, w] : surface.arguments) {
			CHECK(std::abs(u) <= hu);
			CHECK(std::abs(w) <= hw);
		}
		surface.arguments.clear();
		CHECK(Derivation::NDer2_uw(surface, 0, 0, hu, hw, nullptr) == value);
		CHECK(surface.arguments.size() == 4);
	}
}

TEST_CASE("Higher truncation stencils annihilate all lower-degree polynomial moments", "[higher-error]") {
	const Real h = TOL(0.125, 0.25);
	for (int method = 0; method < 4; ++method) {
		const auto witness = witnesses[method];
		for (int degree = 0; degree < witness.degree; ++degree) {
			ErrorPolynomial function;
			function.degree = degree;
			ErrorCurve curve;
			curve.degree = degree;
			const Real exact = degree == witness.derivative ? Real(witness.derivative == 2 ? 2 : 6) : Real(0);
			const auto real = RealEstimate(method, function, h);
			CHECK_THAT(real.value, WithinAbs(exact, TOL(1e-10, 1e-3)));
			CHECK(real.error >= 0);
			CHECK(real.error < TOL(1e-8, 1e-3));
			Real vectorError = -1;
			const auto vector = CurveEstimate(method, curve, h, &vectorError);
			CHECK_THAT(vector[0], WithinAbs(exact, TOL(1e-10, 1e-3)));
			CHECK(vectorError >= 0);
			CHECK(vectorError < TOL(1e-8, 1e-3));
		}
	}
}

TEST_CASE("Surface quadratic second derivatives have roundoff-only residual estimates", "[higher-error]") {
	ErrorSurface surface;
	surface.value = [](Real u, Real w) { return u * u + w * w; };
	Real errorU = -1, errorW = -1;
	CHECK(Derivation::NDer2_uu(surface, 0, 0, Real(0.125), &errorU)[0] == Real(2));
	CHECK(Derivation::NDer2_ww(surface, 0, 0, Real(0.125), &errorW)[0] == Real(2));
	CHECK(errorU >= 0);
	CHECK(errorW >= 0);
	CHECK(errorU < TOL(1e-10, 1e-3));
	CHECK(errorW < TOL(1e-10, 1e-3));
}

TEST_CASE("Mixed estimates include roundoff from both grids and fourth-order outer samples", "[higher-error]") {
	ErrorScalar function;
	function.value = [](const auto&) { return Real(1); };
	const Real h = Real(0.125);
	Real error = -1;
	CHECK(Derivation::NSecDer2Partial(function, 0, 1, VectorN<Real, 2>{0, 0}, h, &error) == Real(0));
	CHECK_THAT(error, WithinRel(Real(20) / 3 * Constants::Eps / (h * h), TOL(1e-10, 1e-3)));
	CHECK(Derivation::NSecDer4Partial(function, 0, 1, VectorN<Real, 2>{0, 0}, h, &error) == Real(0));
	CHECK_THAT(error, WithinRel(Real(68) / 9 * Constants::Eps / (h * h), TOL(1e-10, 1e-3)));
	ErrorSurface surface;
	surface.value = [](Real, Real) { return Real(1); };
	CHECK(Derivation::NDer2_uw(surface, 0, 0, h, h / 2, &error)[0] == Real(0));
	CHECK_THAT(error, WithinRel(Real(20) * Constants::Eps / (h * h / 2), TOL(1e-10, 1e-3)));
}

TEST_CASE("Mixed error refinement fails before callbacks if half-step coordinates collapse", "[higher-error]") {
	const Real center = std::ldexp(Real(1), std::numeric_limits<Real>::digits - 1);
	const Real h = std::nextafter(center, std::numeric_limits<Real>::infinity()) - center;
	for (int accuracy : {2, 4}) {
		ErrorScalar function;
		function.value = [](const auto&) { return Real(0); };
		Real error = -1;
		if (accuracy == 2)
			CHECK_THROWS_AS(Derivation::NSecDer2Partial(function, 0, 1, VectorN<Real, 2>{center, 0}, h, &error), NumericalMethodError);
		else
			CHECK_THROWS_AS(Derivation::NSecDer4Partial(function, 0, 1, VectorN<Real, 2>{0, center}, h, &error), NumericalMethodError);
		CHECK(function.arguments.empty());
		CHECK(error == Real(-1));
	}
	ErrorSurface surface;
	surface.value = [](Real, Real) { return Real(0); };
	for (int axis : {0, 1}) {
		Real error = -1;
		CHECK_THROWS_AS(Derivation::NDer2_uw(surface, axis == 0 ? center : Real(0),
			axis == 1 ? center : Real(0), h, h, &error), NumericalMethodError);
		CHECK(surface.arguments.empty());
		CHECK(error == Real(-1));
	}
}

TEST_CASE("Optional mixed error callbacks propagate failures without writing an error result", "[higher-error]") {
	for (int accuracy : {2, 4}) {
		const int baseCalls = accuracy == 2 ? 4 : 8;
		for (int extra = 1; extra <= 4; ++extra) {
			for (bool nonfinite : {false, true}) {
				ErrorScalar function;
				function.value = [&](const auto& point) {
					if (function.arguments.size() == static_cast<size_t>(baseCalls + extra)) {
						if (nonfinite) return std::numeric_limits<Real>::quiet_NaN();
						throw std::domain_error("mixed error callback failure");
					}
					return point[0] * point[1];
				};
				Real error = -1;
				auto evaluate = [&]() {
					return accuracy == 2 ? Derivation::NSecDer2Partial(function, 0, 1, VectorN<Real, 2>{0, 0}, Real(0.125), &error)
						: Derivation::NSecDer4Partial(function, 0, 1, VectorN<Real, 2>{0, 0}, Real(0.125), &error);
				};
				if (nonfinite) CHECK_THROWS_AS(evaluate(), NumericalMethodError);
				else CHECK_THROWS_AS(evaluate(), std::domain_error);
				CHECK(error == Real(-1));
			}
		}
	}
	for (bool nonfinite : {false, true}) {
		ErrorSurface surface;
		surface.value = [&](Real u, Real w) {
			if (surface.arguments.size() == 5) {
				if (nonfinite) return std::numeric_limits<Real>::infinity();
				throw std::domain_error("surface mixed error callback failure");
			}
			return u * w;
		};
		Real error = -1;
		if (nonfinite) CHECK_THROWS_AS(Derivation::NDer2_uw(surface, 0, 0, Real(0.125), &error), NumericalMethodError);
		else CHECK_THROWS_AS(Derivation::NDer2_uw(surface, 0, 0, Real(0.125), &error), std::domain_error);
		CHECK(error == Real(-1));
	}
}

TEST_CASE("Touched higher derivative entry points reject invalid explicit steps before callbacks", "[higher-error]") {
	for (Real invalid : {Real(0), Real(-1), std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity()}) {
		for (bool estimate : {false, true}) {
			Real error = -1;
			Real* output = estimate ? &error : nullptr;
			for (int method = 0; method < 4; ++method) {
				ErrorCurve curve;
				CHECK_THROWS_AS(CurveEstimate(method, curve, invalid, output), NumericInputError);
				CHECK(curve.calls == 0);
			}
			ErrorScalar function;
			function.value = [](const auto&) { return Real(0); };
			for (int other : {0, 1}) {
				CHECK_THROWS_AS(Derivation::NSecDer2Partial(function, 0, other, VectorN<Real, 2>{0, 0}, invalid, output), NumericInputError);
				CHECK_THROWS_AS(Derivation::NSecDer4Partial(function, 0, other, VectorN<Real, 2>{0, 0}, invalid, output), NumericInputError);
			}
			CHECK(function.arguments.empty());
			ErrorSurface surface;
			surface.value = [](Real, Real) { return Real(0); };
			CHECK_THROWS_AS(Derivation::NDer2_uu(surface, 0, 0, invalid, output), NumericInputError);
			CHECK_THROWS_AS(Derivation::NDer2_ww(surface, 0, 0, invalid, output), NumericInputError);
			CHECK_THROWS_AS(Derivation::NDer2_uw(surface, 0, 0, invalid, Real(0.125), output), NumericInputError);
			CHECK_THROWS_AS(Derivation::NDer2_uw(surface, 0, 0, Real(0.125), invalid, output), NumericInputError);
			CHECK(surface.arguments.empty());
			CHECK(error == Real(-1));
		}
	}
}

TEST_CASE("Mixed estimates preserve axis symmetry and default-step overloads", "[higher-error]") {
	ErrorScalar function;
	function.value = [](const auto& point) { return point[0] * point[1]; };
	const VectorN<Real, 2> point{Real(0.25), Real(-0.5)};
	for (int accuracy : {2, 4}) {
		Real forwardError = -1, reverseError = -1;
		const Real forward = accuracy == 2 ? Derivation::NSecDer2Partial(function, 0, 1, point, &forwardError)
			: Derivation::NSecDer4Partial(function, 0, 1, point, &forwardError);
		const Real reverse = accuracy == 2 ? Derivation::NSecDer2Partial(function, 1, 0, point, &reverseError)
			: Derivation::NSecDer4Partial(function, 1, 0, point, &reverseError);
		CHECK_THAT(forward, WithinAbs(Real(1), TOL(1e-5, 1e-2)));
		CHECK_THAT(reverse, WithinAbs(forward, TOL(1e-5, 1e-2)));
		CHECK(std::isfinite(forwardError));
		CHECK(std::isfinite(reverseError));
		CHECK(forwardError >= 0);
		CHECK(reverseError >= 0);
	}
}

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "../../TestPrecision.h"

#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/base/ChebyshevApproximation.h>
#include <mml/base/InterpolatedFunction.h>
#endif

#include <array>
#include <cmath>
#include <functional>
#include <limits>

using namespace MML;

namespace {
	void CheckNumericalValue(Real actual, Real expected) {
		CHECK(std::isfinite(actual));
		CHECK_THAT(actual, Catch::Matchers::WithinAbs(
			expected, Precision::EqualityTolerance * std::max(Real(1), std::abs(expected))));
	}

	template<typename Interpolator>
	void CheckBarycentricScale(const Interpolator& interpolator, Real scale) {
		for (int i = 0; i <= 24; ++i) {
			const Real u = Real(i) / 8;
			CheckNumericalValue(interpolator(u * scale), 1 + u * u);
		}
		Real maximum = 0;
		for (int i = 0; i < 4; ++i) {
			CHECK(std::isfinite(interpolator.getWeight(i)));
			CHECK(interpolator.getWeight(i) != 0);
			maximum = std::max(maximum, std::abs(interpolator.getWeight(i)));
		}
		CheckNumericalValue(maximum, 1);
	}
}

TEST_CASE("Two-point splines preserve natural clamped and mixed boundary conditions",
	"[interpolation-numerics][interpolation][spline]") {
	for (bool descending : {false, true}) {
		CAPTURE(descending);
		const Vector<Real> x = descending ? Vector<Real>{2, 0} : Vector<Real>{0, 2};
		const Vector<Real> y = descending ? Vector<Real>{8, 0} : Vector<Real>{0, 8};
		const Real natural = SplineInterpRealFunc::NaturalBoundaryDerivative;
		const SplineInterpRealFunc linear(x, y);
		const SplineInterpRealFunc clamped(x, y, descending ? Real(12) : Real(0), descending ? Real(0) : Real(12));
		const SplineInterpRealFunc naturalAtZero(x, y, descending ? 12 : natural, descending ? natural : 12);
		const SplineInterpRealFunc naturalAtTwo(x, y, descending ? natural : 0, descending ? 0 : natural);
		SplineInterpRealFunc reinitialized(x, y);
		reinitialized.initSecDerivs(&x[0], &y[0], descending ? Real(12) : Real(0), descending ? Real(0) : Real(12));
		for (int i = 0; i <= 8; ++i) {
			const Real query = Real(i) / 4;
			CAPTURE(query);
			CheckNumericalValue(linear(query), 4 * query);
			CheckNumericalValue(linear.Derivative(query), 4);
			CheckNumericalValue(linear.SecondDerivative(query), 0);
			for (const auto* cubic : std::array<const SplineInterpRealFunc*, 3>{&clamped, &naturalAtZero, &reinitialized}) {
				CheckNumericalValue((*cubic)(query), query * query * query);
				CheckNumericalValue(cubic->Derivative(query), 3 * query * query);
				CheckNumericalValue(cubic->SecondDerivative(query), 6 * query);
			}
			CheckNumericalValue(naturalAtTwo(query), 3 * query * query - query * query * query / 2);
			CheckNumericalValue(naturalAtTwo.Derivative(query), 6 * query - 3 * query * query / 2);
			CheckNumericalValue(naturalAtTwo.SecondDerivative(query), 6 - 3 * query);
		}
		CheckNumericalValue(linear.Integrate(0, 2), 8);
		CheckNumericalValue(clamped.Integrate(0, 2), 4);
		CheckNumericalValue(naturalAtZero.Integrate(0, 2), 4);
		CheckNumericalValue(naturalAtTwo.Integrate(0, 2), 6);
	}
}

TEST_CASE("Chebyshev antiderivative retains the highest coefficient", "[interpolation-numerics][chebyshev]") {
	ChebyshevApproximation linear(Vector<Real>{0, 1}, -1, 1);
	auto integral = linear.Integral();
	CHECK(integral.NumCoefficients() == 3);
	for (int i = 0; i <= 16; ++i) {
		const Real x = -1 + Real(i) / 8;
		CheckNumericalValue(integral(x), (x * x - 1) / 2);
	}
}

TEST_CASE("Chebyshev minimal degrees integrate monomials on mapped domains",
	"[interpolation-numerics][chebyshev][unsafe-baseline]") {
	const std::array<Vector<Real>, 7> powers{
		Vector<Real>{2}, Vector<Real>{0, 1}, Vector<Real>{1, 0, Real(1) / 2},
		Vector<Real>{0, Real(3) / 4, 0, Real(1) / 4},
		Vector<Real>{Real(3) / 4, 0, Real(1) / 2, 0, Real(1) / 8},
		Vector<Real>{0, Real(5) / 8, 0, Real(5) / 16, 0, Real(1) / 16},
		Vector<Real>{Real(5) / 8, 0, Real(15) / 32, 0, Real(3) / 16, 0, Real(1) / 32}
	};
	for (const auto& domain : {std::array<Real, 2>{-1, 1}, std::array<Real, 2>{2, 6}}) {
		const Real halfWidth = (domain[1] - domain[0]) / 2;
		for (int power = 0; power < int(powers.size()); ++power) {
			CAPTURE(power, domain[0], domain[1]);
			ChebyshevApproximation function(powers[power], domain[0], domain[1]);
			auto integral = function.Integral();
			auto derivative = function.Derivative();
			CHECK(integral.NumCoefficients() == function.NumTerms() + 1);
			for (int i = 0; i <= 16; ++i) {
				const Real u = -1 + Real(i) / 8;
				const Real x = domain[0] + halfWidth * (u + 1);
				CheckNumericalValue(integral(x),
					halfWidth * (std::pow(u, MML::Real(power + 1)) - std::pow(Real(-1), MML::Real(power + 1))) / (power + 1));
				CheckNumericalValue(derivative(x), power == 0 ? 0 : power * std::pow(u, MML::Real(power - 1)) / halfWidth);
			}
		}
	}
	ChebyshevApproximation empty;
	CheckNumericalValue(empty.Integral()(0), 0);
	CheckNumericalValue(empty.Derivative()(0), 0);
}

TEST_CASE("Chebyshev calculus uses only the active series", "[interpolation-numerics][chebyshev]") {
	ChebyshevApproximation function(Vector<Real>{6, 2, 4, 8}, 2, 6);
	function.SetNumTerms(2);
	auto derivative = function.Derivative();
	auto integral = function.Integral();
	for (int i = 0; i <= 16; ++i) {
		const Real x = 2 + Real(i) / 4;
		CheckNumericalValue(function(x), x - 1);
		CheckNumericalValue(derivative(x), 1);
		CheckNumericalValue(integral(x), (x * x / 2 - x));
	}
	function.SetNumTerms(1);
	CheckNumericalValue(function.Derivative()(3), 0);
	CheckNumericalValue(function.Integral()(3), 3);
}

TEST_CASE("Chebyshev fixed-term success is not a convergence claim", "[interpolation-numerics][chebyshev]") {
	auto result = ApproximateChebyshevDetailed([](Real x) { return x * x; }, -1, 1, 1, 9, 1);
	CHECK(result.IsSuccess());
	CHECK_FALSE(result.converged);
	CHECK(result.max_error_estimate > Real(1) / 2);
	CHECK(result.samples_used == 10);
}

TEST_CASE("Chebyshev diagnostics reject invalid sampling and nonfinite data",
	"[interpolation-numerics][chebyshev]") {
	int calls = 0;
	std::function<Real(Real)> function = [&](Real x) { ++calls; return x; };
	ChebyshevApproximation approximation(Vector<Real>{0, 1}, -1, 1);
	for (int count : {0, -1}) {
		CHECK_THROWS_AS(approximation.MaxError(function, count), ArgumentError);
		CHECK(calls == 0);
	}
	CHECK_THROWS_AS(ApproximateChebyshevDetailed(function, -1, 1, 4, 1, 1), ArgumentError);
	CHECK_THROWS_AS(ApproximateChebyshevDetailed(function, -1, 1, 4, 9, 0), ArgumentError);
	CHECK(calls == 0);
	for (Real bad : {std::numeric_limits<Real>::infinity(), std::numeric_limits<Real>::quiet_NaN()}) {
		std::function<Real(Real)> invalid = [=](Real) { return bad; };
		CHECK_THROWS_AS(ChebyshevApproximation(invalid, -1, 1, 4), NumericInputError);
		CHECK_THROWS_AS(approximation.MaxError(invalid, 4), NumericInputError);
		CHECK_THROWS_AS(ApproximateChebyshevDetailed(invalid, -1, 1, 4), NumericInputError);
		ChebyshevApproximationConfig config;
		config.tolerance = bad;
		CHECK_THROWS_AS(ApproximateChebyshevAdaptive(function, -1, 1, config), NumericInputError);
	}
	std::function<Real(Real)> invalidValidation = [](Real x) {
		return x == -1 ? std::numeric_limits<Real>::quiet_NaN() : x;
	};
	CHECK_THROWS_AS(ApproximateChebyshevDetailed(invalidValidation, -1, 1, 4), NumericInputError);
}

TEST_CASE("Chebyshev diagnostics retain representable large RMS errors", "[interpolation-numerics][chebyshev]") {
	const Real scale = std::sqrt(std::numeric_limits<Real>::max()) * 16;
	auto result = ApproximateChebyshevDetailed([=](Real x) { return scale * x; }, -1, 1, 1, 3, 1);
	CheckNumericalValue(result.max_error_estimate / scale, 1);
	CheckNumericalValue(result.rms_error_estimate / scale, std::sqrt(Real(2) / 3));
}

TEST_CASE("Chebyshev finite construction and diagnostics cover both callable interfaces",
	"[interpolation-numerics][chebyshev]") {
	struct NonfiniteFunction : IRealFunction {
		Real value;
		explicit NonfiniteFunction(Real bad) : value(bad) {}
		Real operator()(Real) const override { return value; }
	};
	for (Real bad : {std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity()}) {
		NonfiniteFunction function(bad);
		CHECK_THROWS_AS(ChebyshevApproximation(function, -1, 1, 4), NumericInputError);
		ChebyshevApproximation linear(Vector<Real>{0, 1}, -1, 1);
		CHECK_THROWS_AS(linear.MaxError(function, 4), NumericInputError);
		CHECK_THROWS_AS(ChebyshevApproximation(Vector<Real>{0, bad}, -1, 1), NumericInputError);
		CHECK_THROWS_AS(ChebyshevApproximation(Vector<Real>{0, 1}, bad, 1), NumericInputError);
	}
	const Real largest = std::numeric_limits<Real>::max();
	ChebyshevApproximation mapped(Vector<Real>{0, 1}, largest / 2, largest);
	CheckNumericalValue(mapped(largest / 2), -1);
	CheckNumericalValue(mapped(largest), 1);
	CHECK_THROWS_AS(ChebyshevApproximation(Vector<Real>{0, 1}, -largest, largest), ArgumentError);
	ChebyshevApproximation largeConstant([=](Real) { return largest / 8; }, -1, 1, 8);
	CheckNumericalValue(largeConstant(0) / (largest / 8), 1);
}

TEST_CASE("Rational exact nodes take precedence over proximity", "[interpolation-numerics][rational]") {
	for (Real spacing : {std::numeric_limits<Real>::epsilon() / 2, Real(1)}) {
		const Real offset = spacing == 1 ? Real(1) / std::numeric_limits<Real>::epsilon() : 0;
		Vector<Real> x{offset, offset + spacing};
		Vector<Real> y{2, 7};
		RationalInterpRealFunc ascending(x, y, 2);
		RationalInterpRealFunc descending(Vector<Real>{x[1], x[0]}, Vector<Real>{7, 2}, 2);
		CHECK(ascending(x[0]) == 2);
		CHECK(ascending(x[1]) == 7);
		CHECK(descending(x[0]) == 2);
		CHECK(descending(x[1]) == 7);
	}
	const Real spacing = std::numeric_limits<Real>::epsilon() / 2;
	RationalInterpRealFunc compressed(Vector<Real>{0, spacing}, Vector<Real>{1, Real(1) / 2}, 2);
	CheckNumericalValue(compressed(spacing / 2), Real(2) / 3);
}

TEST_CASE("Rational constant corrections are removable but genuine poles still fail",
	"[interpolation-numerics][rational][float-final]") {
	for (Real value : {Real(0), Real(1), Real(-3)}) {
		for (bool descending : {false, true}) {
			Vector<Real> x = descending ? Vector<Real>{3, 2, 1, 0} : Vector<Real>{0, 1, 2, 3};
			RationalInterpRealFunc function(x, Vector<Real>{value, value, value, value}, 4);
			for (Real query : {Real(-1), Real(1) / 2, Real(5) / 4, Real(4)}) {
				CheckNumericalValue(function(query), value);
				CHECK(function.getLastErrorEst() == 0);
			}
		}
	}
	for (bool descending : {false, true}) {
		RationalInterpRealFunc pole(descending ? Vector<Real>{1, 0} : Vector<Real>{0, 1},
			descending ? Vector<Real>{-1, 1} : Vector<Real>{1, -1}, 2);
		CHECK_THROWS_AS(pole(Real(1) / 2), RealFuncInterpRuntimeError);
	}
}

TEST_CASE("Rational low-degree curves do not manufacture float poles", "[interpolation-numerics][rational][float-final]") {
	for (bool descending : {false, true}) {
		for (Real scale : {Real(1) / 8, Real(1), Real(8)}) {
			Vector<Real> nodes{0, 1, 2, 3}, values{scale, scale / 2, scale / 3, scale / 4};
			if (descending) {
				std::reverse(nodes.begin(), nodes.end());
				std::reverse(values.begin(), values.end());
			}
			RationalInterpRealFunc reciprocal(nodes, values, 4);
			for (int i = 0; i <= 24; ++i) {
				const Real x = Real(i) / 8;
				CAPTURE(descending, scale, x);
				CheckNumericalValue(reciprocal(x), scale / (1 + x));
				CHECK(std::isfinite(reciprocal.getLastErrorEst()));
			}
		}
	}
}

TEST_CASE("Spline integrals follow physical bounds in either node order", "[interpolation-numerics][spline]") {
	for (bool descending : {false, true}) {
		Vector<Real> x = descending ? Vector<Real>{3, 2, Real(1) / 2, 0} : Vector<Real>{0, Real(1) / 2, 2, 3};
		Vector<Real> y(x.size());
		for (int i = 0; i < x.size(); ++i) y[i] = x[i] * x[i] * x[i] - 2 * x[i] + 1;
		auto derivative = [](Real x) { return 3 * x * x - 2; };
		auto primitive = [](Real x) { return x * x * x * x / 4 - x * x + x; };
		SplineInterpRealFunc spline(x, y, derivative(x[0]), derivative(x[x.size() - 1]));
		for (Real a : {Real(0), Real(1) / 4, Real(1), Real(2)}) {
			for (Real b : {Real(1) / 2, Real(3) / 2, Real(3)}) {
				CheckNumericalValue(spline.Integrate(a, b), primitive(b) - primitive(a));
				CheckNumericalValue(spline.Integrate(b, a), primitive(a) - primitive(b));
			}
		}
		CheckNumericalValue(spline.Integrate(1, 1), 0);
		CHECK_THROWS_AS(spline.Integrate(-1, 1), DomainError);
		CHECK_THROWS_AS(spline.Integrate(1, 4), DomainError);
		CHECK_THROWS_AS(spline.Integrate(std::numeric_limits<Real>::quiet_NaN(), 1), DomainError);
	}
}

TEST_CASE("Barycentric interpolation preserves subnormal affine values", "[interpolation-numerics][barycentric]") {
	const Real tiny = std::numeric_limits<Real>::denorm_min();
	for (Real direction : {Real(-1), Real(1)}) {
		BarycentricPolynomialInterp polynomial(Vector<Real>{0, 1}, Vector<Real>{0, 1});
		BarycentricRationalInterp rational(Vector<Real>{0, 1}, Vector<Real>{0, 1}, 1);
		const Real query = direction * tiny;
		CHECK(polynomial(query) == query);
		CHECK(rational(query) == query);
		CHECK(polynomial.EvaluateDetailed(query).value == query);
		CHECK(rational.EvaluateDetailed(query).value == query);
	}
}

TEST_CASE("Barycentric relative weights and values are stable across node scales",
	"[interpolation-numerics][barycentric]") {
	for (Real scale : {Real(1), std::sqrt(std::numeric_limits<Real>::min()),
		std::sqrt(std::numeric_limits<Real>::max()) / 16}) {
		Vector<Real> x{0, scale, 2 * scale, 3 * scale};
		Vector<Real> y{1, 2, 5, 10};
		BarycentricPolynomialInterp polynomial(x, y);
		BarycentricRationalInterp rational(x, y, 3);
		CheckBarycentricScale(polynomial, scale);
		CheckBarycentricScale(rational, scale);
		BarycentricRationalInterp descending(Vector<Real>{x[3], x[2], x[1], x[0]},
			Vector<Real>{10, 5, 2, 1}, 3);
		CheckBarycentricScale(descending, scale);
	}
}

TEST_CASE("Floater Hormann blending reproduces lower-degree polynomials at every order",
	"[interpolation-numerics][barycentric]") {
	for (Real scale : {Real(1), std::sqrt(std::numeric_limits<Real>::min()),
		std::sqrt(std::numeric_limits<Real>::max()) / 16}) {
		for (bool descending : {false, true}) {
			for (int order = 0; order <= 3; ++order) {
				Vector<Real> x(5), y(5);
				for (int i = 0; i < 5; ++i) {
					const Real u = Real(descending ? 4 - i : i);
					x[i] = u * scale;
					y[i] = 1 + std::pow(u, MML::Real(order));
				}
				BarycentricRationalInterp function(x, y, order);
				for (int i = 0; i <= 32; ++i) {
					const Real u = Real(i) / 8;
					CheckNumericalValue(function(u * scale), 1 + std::pow(u, MML::Real(order)));
				}
				if (order == 2) {
					const std::array<Real, 5> weights{Real(1) / 4, -Real(3) / 4, 1, -Real(3) / 4, Real(1) / 4};
					for (int i = 0; i < 5; ++i) CheckNumericalValue(function.getWeight(i), weights[i]);
				}
			}
		}
	}
}

TEST_CASE("Barycentric extreme spans ordinates and failure policies stay explicit",
	"[interpolation-numerics][barycentric]") {
	const Real largest = std::numeric_limits<Real>::max();
	BarycentricPolynomialInterp polynomial(Vector<Real>{-largest, largest}, Vector<Real>{-1, 1});
	BarycentricRationalInterp rational(Vector<Real>{-largest, largest}, Vector<Real>{-1, 1}, 1);
	for (Real fraction : {Real(-1), -Real(1) / 2, Real(0), Real(1) / 2, Real(1)}) {
		CheckNumericalValue(polynomial(fraction * largest), fraction);
		CheckNumericalValue(rational(fraction * largest), fraction);
	}
	BarycentricPolynomialInterp distant(Vector<Real>{-largest, -largest * (Real(3) / 4)}, Vector<Real>{-1, -Real(3) / 4});
	CheckNumericalValue(distant(largest), 1);
	BarycentricPolynomialInterp constant(Vector<Real>{-1, 0, 1}, Vector<Real>{largest, largest, largest});
	CHECK(constant(Real(1) / 4) == largest);
	const Real tiny = std::numeric_limits<Real>::denorm_min();
	BarycentricPolynomialInterp minimal(Vector<Real>{0, tiny}, Vector<Real>{0, 1});
	CheckNumericalValue(minimal(-tiny), -1);
	CheckNumericalValue(minimal(2 * tiny), 2);
	CHECK_THROWS_AS(BarycentricPolynomialInterp(
		Vector<Real>{0, std::numeric_limits<Real>::quiet_NaN()}, Vector<Real>{0, 1}), RealFuncInterpInitError);
	auto checkFailure = [](const auto& interpolator) {
		const Real bad = std::numeric_limits<Real>::quiet_NaN();
		CHECK_THROWS_AS(interpolator(bad), RealFuncInterpRuntimeError);
		InterpolationConfig config;
		config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
		auto result = interpolator.EvaluateDetailed(bad, config);
		CHECK_FALSE(result.IsSuccess());
		CHECK_FALSE(result.error_message.empty());
		CHECK(result.function_evaluations == 1);
		config.exception_policy = EvaluationExceptionPolicy::Propagate;
		CHECK_THROWS_AS(interpolator.EvaluateDetailed(bad, config), RealFuncInterpRuntimeError);
	};
	checkFailure(polynomial);
	checkFailure(rational);
}

TEST_CASE("Closed linear curves traverse the closing chord with arc-length timing",
	"[interpolation-numerics][parametric]") {
	Matrix<Real> points(3, 2, {0, 0, 3, 0, 3, 4});
	LinInterpParametricCurve<2> closed(points, true);
	const std::array<Real, 5> times{0, Real(3) / 12, Real(7) / 12, Real(19) / 24, 1};
	const std::array<std::array<Real, 2>, 5> expected{{{0, 0}, {3, 0}, {3, 4}, {Real(3) / 2, 2}, {0, 0}}};
	for (int i = 0; i < int(times.size()); ++i) {
		for (Real shift : {Real(-2), Real(0), Real(3)}) {
			auto actual = closed(times[i] + shift);
			for (int j = 0; j < 2; ++j) CheckNumericalValue(actual[j], expected[i][j]);
		}
	}
	const Real before = std::nextafter(Real(1), Real(0));
	for (int j = 0; j < 2; ++j) CheckNumericalValue(closed(before)[j], 0);
	LinInterpParametricCurve<2> open(points);
	CheckNumericalValue(open(-1)[0], 0);
	CheckNumericalValue(open(2)[0], 3);
	CheckNumericalValue(open(2)[1], 4);
	Matrix<Real> repeated(3, 2, {2, 3, 2, 3, 2, 3});
	LinInterpParametricCurve<2> constant(repeated, true);
	CheckNumericalValue(constant(Real(1) / 2)[0], 2);
	CheckNumericalValue(constant(Real(1) / 2)[1], 3);
	Matrix<Real> duplicated(4, 2, {0, 0, 3, 0, 3, 4, 0, 0});
	LinInterpParametricCurve<2> alreadyClosed(duplicated, true);
	for (Real t : times)
		for (int j = 0; j < 2; ++j) CheckNumericalValue(alreadyClosed(t)[j], closed(t)[j]);
}

TEST_CASE("Closed spline custom ranges normalize before wrapping and preserve range rejection",
	"[interpolation-numerics][parametric]") {
	Matrix<Real> points(4, 2, {0, 0, 1, 0, 1, 1, 0, 1});
	SplineInterpParametricCurve<2> normalized(points, true);
	for (const auto& domain : {std::array<Real, 2>{2, 4}, std::array<Real, 2>{-3, 5},
		std::array<Real, 2>{Real(1) / 4, Real(3) / 4}}) {
		SplineInterpParametricCurve<2> mapped(domain[0], domain[1], points, true);
		for (int i = 0; i <= 32; ++i) {
			const Real u = Real(i) / 32;
			auto expected = normalized(u);
			auto actual = mapped(domain[0] + (domain[1] - domain[0]) * u);
			for (int j = 0; j < 2; ++j) CheckNumericalValue(actual[j], expected[j]);
		}
		for (int i = 0; i <= 4; ++i) {
			const Real t = domain[0] + (domain[1] - domain[0]) * Real(i) / 4;
			for (int j = 0; j < 2; ++j) CheckNumericalValue(mapped(t)[j], points(i % 4, j));
		}
		for (int j = 0; j < 2; ++j)
			CheckNumericalValue(mapped(std::nextafter(domain[1], domain[0]))[j], points(0, j));
		CHECK_THROWS_AS(mapped(std::nextafter(domain[0], -std::numeric_limits<Real>::infinity())), DomainError);
		CHECK_THROWS_AS(mapped(std::nextafter(domain[1], std::numeric_limits<Real>::infinity())), DomainError);
	}
}

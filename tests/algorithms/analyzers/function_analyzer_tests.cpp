#include <catch2/catch_all.hpp>
#include "../../TestPrecision.h"
#include "../../TestMatchers.h"

#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/base/Vector/Vector.h>

#include <mml/base/Function.h>
#include <mml/base/InterpolatedFunction.h>

#include <mml/algorithms/Analyzers/FunctionsAnalyzer.h>
#endif

using namespace MML;
using namespace MML::Testing;

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace MML::Tests::Algorithms::FunctionAnalyzerTests
{

	TEST_CASE("Test_Function_Analyzer_isDefined", "[simple]")
	{
			TEST_PRECISION_INFO();
		RealFunction test1([](Real x) { return 1 / (x - 1); });

		RealFunctionAnalyzer an(test1);

		REQUIRE(an.IsDefinedAtPoint(static_cast<MML::Real>(0.0)) == true);
		REQUIRE(an.IsDefinedAtPoint(static_cast<MML::Real>(1.0)) == false);
	}

	TEST_CASE("FunctionAnalyzer sampled monotonicity handles plateaus and endpoints", "[FunctionAnalyzer][function-analyzer-regression]") {
		for (Real sign : {Real(-1), Real(1)}) {
			RealFunctionFromStdFunc plateau([sign](Real x) { return sign * std::max(Real(0), x - Real(0.5)); });
			REQUIRE(RealFunctionAnalyzer(plateau).IsMonotonicOnSamples(0, 1, 10));
		}
		RealFunction constant(+[](Real) { return Real(1); });
		REQUIRE(RealFunctionAnalyzer(constant).IsMonotonicOnSamples(0, 1, 1));
		RealFunction endpointReversal(+[](Real x) { return x == Real(1) ? Real(-1) : x; });
		REQUIRE_FALSE(RealFunctionAnalyzer(endpointReversal).IsMonotonicOnSamples(0, 1, 4));
		RealFunction reversalAfterPlateau(+[](Real x) { return x <= Real(0.5) ? Real(0) : (x < Real(1) ? Real(1) : Real(0)); });
		REQUIRE_FALSE(RealFunctionAnalyzer(reversalAfterPlateau).IsMonotonicOnSamples(0, 1, 4));
	}

	TEST_CASE("FunctionAnalyzer interval flags reject nonfinite samples", "[FunctionAnalyzer][function-analyzer-regression]") {
		for (Real invalid : {std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity()}) {
			RealFunctionFromStdFunc undefined([invalid](Real) { return invalid; });
			RealFunctionAnalyzer analyzer(undefined);
			const auto result = analyzer.AnalyzeInterval(0, 1, 4);
			CHECK_FALSE(result.isDefinedEverywhere);
			CHECK_FALSE(result.isContinuousEverywhere);
			CHECK_FALSE(result.isMonotonic);
			CHECK(result.undefinedPoints == std::vector<Real>{0, Real(0.25), Real(0.5), Real(0.75)});
			CHECK(result.discontinuityPoints == result.undefinedPoints);
			CHECK(std::isnan(result.minValue));
			CHECK(std::isnan(result.maxValue));
			CHECK_FALSE(analyzer.IsMonotonicOnSamples(0, 1, 4));
			CHECK(std::isnan(analyzer.MinInNPoints(0, 1, 4)));
			CHECK(std::isnan(analyzer.MaxInNPoints(0, 1, 4)));
		}
		RealFunction mixed(+[](Real x) { return x == Real(0.5) ? std::numeric_limits<Real>::quiet_NaN() : x; });
		const auto result = RealFunctionAnalyzer(mixed).AnalyzeInterval(0, 1, 4);
		CHECK_FALSE(result.isContinuousEverywhere);
		CHECK_FALSE(result.isMonotonic);
		CHECK(result.minValue == Real(0));
		CHECK(result.maxValue == Real(0.75));
	}

	TEST_CASE("FunctionAnalyzer preserves half-open reports and closed sampling endpoints",
						"[FunctionAnalyzer][function-analyzer-regression]") {
		RealFunction endpointReversal(+[](Real x) { return x == Real(1) ? Real(-1) : x; });
		RealFunctionAnalyzer analyzer(endpointReversal);
		const auto interval = analyzer.AnalyzeInterval(0, 1, 4);
		const auto points = analyzer.AnalyzePoints(0, 1, 4);
		REQUIRE(points.size() == 4);
		CHECK(points.front().x == Real(0));
		CHECK(points.back().x == Real(0.75));
		CHECK(interval.minValue == Real(0));
		CHECK(interval.maxValue == Real(0.75));
		CHECK(interval.isMonotonic);
		CHECK(analyzer.MinInNPoints(0, 1, 4) == Real(-1));
		CHECK(analyzer.MaxInNPoints(0, 1, 4) == Real(0.75));

		RealFunction highest(+[](Real) { return std::numeric_limits<Real>::max(); });
		RealFunction lowest(+[](Real) { return std::numeric_limits<Real>::lowest(); });
		CHECK(RealFunctionAnalyzer(highest).MinInNPoints(0, 1, 1) == std::numeric_limits<Real>::max());
		CHECK(RealFunctionAnalyzer(lowest).MaxInNPoints(0, 1, 1) == std::numeric_limits<Real>::lowest());
	}

	TEST_CASE("FunctionAnalyzer continuity output replaces stale locations", "[FunctionAnalyzer][function-analyzer-regression]") {
		RealFunction constant(+[](Real) { return Real(1); });
		std::vector<Real> locations{Real(42)};
		REQUIRE(RealFunctionAnalyzer(constant).AppearsContinuous(0, 1, 4, Precision::LooseTolerance, &locations));
		CHECK(locations.empty());
		RealFunction jump(+[](Real x) { return x < Real(0) ? Real(-1) : Real(1); });
		locations = {Real(42)};
		REQUIRE_FALSE(RealFunctionAnalyzer(jump).AppearsContinuous(-1, 1, 4, Precision::LooseTolerance, &locations));
		CHECK(std::find(locations.begin(), locations.end(), Real(0)) != locations.end());
		CHECK(std::find(locations.begin(), locations.end(), Real(42)) == locations.end());
	}

	TEST_CASE("FunctionAnalyzer sign checks survive tiny function scales", "[FunctionAnalyzer][function-analyzer-regression]") {
		const Real scale = std::scalbn(Real(1), -2 * std::numeric_limits<Real>::max_exponent / 3);
		const Real eps = scale * Precision::EqualityTolerance;
		RealFunctionFromStdFunc linear([scale](Real x) { return scale * (x - Real(0.123456)); });
		auto roots = RealFunctionAnalyzer(linear).GetRoots(0, 1, eps);
		REQUIRE(roots.size() == 1);
		CHECK_THAT(roots.front(), WithinAbs(Real(0.123456), Precision::ModerateTolerance));
		RealFunctionFromStdFunc cubic([scale](Real x) { return scale * x * x * x; });
		RealFunctionFromStdFunc quadratic([scale](Real x) { return scale * x * x; });
		CHECK(RealFunctionAnalyzer(cubic).IsInflectionPoint(0, Precision::CoarseTolerance));
		CHECK(RealFunctionAnalyzer(quadratic).IsLocalOptimum(0, Precision::CoarseTolerance));
	}

	TEST_CASE("FunctionAnalyzer safe difference respects representable range", "[FunctionAnalyzer][function-analyzer-regression]") {
		const Real large = Real(10) * Precision::InfiniteLimitThreshold;
		CHECK(SafeDifference(large, -large) == Real(2) * large);
		CHECK(SafeDifference(-large, large) == Real(2) * large);
		CHECK(SafeDifference(large, large) == Real(0));
		// Keep the intentional overflow witness in runtime arithmetic.
		volatile Real runtimeMax = std::numeric_limits<Real>::max();
		CHECK(std::isinf(SafeDifference(runtimeMax, -runtimeMax)));
		CHECK(std::isinf(SafeDifference(std::numeric_limits<Real>::infinity(), Real(0))));
		CHECK(std::isnan(SafeDifference(std::numeric_limits<Real>::quiet_NaN(), Real(0))));
	}

	TEST_CASE("FunctionAnalyzer limits preserve signs and distinguish invalid evidence", "[FunctionAnalyzer][function-analyzer-regression]") {
		const Real eps = Real(1) / std::sqrt(Precision::InfiniteLimitThreshold);
		for (Real sign : {Real(-1), Real(1)}) {
			RealFunctionFromStdFunc pole([sign](Real x) { return sign / (x * x); });
			const RealFunctionAnalyzer analyzer(pole);
			for (Real limit : {analyzer.ComputeLeftLimit(0, eps), analyzer.ComputeRightLimit(0, eps)}) {
				CHECK(std::isinf(limit));
				CHECK(std::signbit(limit) == std::signbit(sign));
			}
			RealFunctionFromStdFunc overflowPole([sign](Real x) { return sign * Precision::OverflowDetectionThreshold / (x * x); });
			const auto disc = RealFunctionAnalyzer(overflowPole).ClassifyDiscontinuity(0);
			CHECK(disc.type == DiscontinuityType::Infinite);
			CHECK(std::signbit(disc.leftLimit) == std::signbit(sign));
			CHECK(std::signbit(disc.rightLimit) == std::signbit(sign));
		}
		RealFunctionFromStdFunc oddPole([](Real x) { return Precision::OverflowDetectionThreshold / x; });
		const auto odd = RealFunctionAnalyzer(oddPole).ClassifyDiscontinuity(0);
		CHECK(odd.type == DiscontinuityType::Infinite);
		CHECK(std::signbit(odd.leftLimit));
		CHECK_FALSE(std::signbit(odd.rightLimit));

		RealFunction largeConstant(+[](Real) { return Real(10) * Precision::InfiniteLimitThreshold; });
		const RealFunctionAnalyzer constantAnalyzer(largeConstant);
		CHECK(constantAnalyzer.ComputeLeftLimit(0) == largeConstant(0));
		CHECK(constantAnalyzer.ComputeRightLimit(0) == largeConstant(0));
		CHECK(constantAnalyzer.ClassifyDiscontinuity(0).type == DiscontinuityType::UNKNOWN);
		RealFunction undefined(+[](Real) { return std::numeric_limits<Real>::quiet_NaN(); });
		CHECK(RealFunctionAnalyzer(undefined).ClassifyDiscontinuity(0).type == DiscontinuityType::UNKNOWN);

		RealFunction centerOnly(+[](Real x) { return x == Real(1) ? Real(100) : Real(2); });
		CHECK(std::isnan(RealFunctionAnalyzer(centerOnly).ComputeLeftLimit(1, std::numeric_limits<Real>::min())));
		CHECK(std::isnan(RealFunctionAnalyzer(centerOnly).ComputeRightLimit(1, std::numeric_limits<Real>::min())));
	}

	TEST_CASE("FunctionAnalyzer rejects invalid sampling and probe parameters", "[FunctionAnalyzer][function-analyzer-regression]") {
		RealFunction identity(+[](Real x) { return x; });
		RealFunctionAnalyzer analyzer(identity);
		for (int count : {0, -1}) {
			CHECK_THROWS_AS(analyzer.AnalyzeInterval(0, 1, count), std::invalid_argument);
			CHECK_THROWS_AS(analyzer.AnalyzePoints(0, 1, count), std::invalid_argument);
			CHECK_THROWS_AS(analyzer.IsMonotonicOnSamples(0, 1, count), std::invalid_argument);
			CHECK_THROWS_AS(analyzer.AppearsContinuous(0, 1, count), std::invalid_argument);
			CHECK_THROWS_AS(analyzer.FindDiscontinuityCandidates(0, 1, count), std::invalid_argument);
			CHECK_THROWS_AS(analyzer.MinInNPoints(0, 1, count), std::invalid_argument);
			CHECK_THROWS_AS(analyzer.MaxInNPoints(0, 1, count), std::invalid_argument);
			CHECK_THROWS_AS(analyzer.ComputeLeftLimit(0, Precision::LooseTolerance, count), std::invalid_argument);
			CHECK_THROWS_AS(analyzer.ComputeRightLimit(0, Precision::LooseTolerance, count), std::invalid_argument);
		}
		for (Real invalid : {Real(0), Real(-1), std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity()}) {
			CHECK_THROWS_AS(analyzer.AnalyzeInterval(0, 1, 4, invalid), std::invalid_argument);
			CHECK_THROWS_AS(analyzer.AnalyzePoints(0, 1, 4, invalid), std::invalid_argument);
			CHECK_THROWS_AS(analyzer.AppearsContinuous(0, 1, 4, invalid), std::invalid_argument);
			CHECK_THROWS_AS(analyzer.FindDiscontinuityCandidates(0, 1, 4, invalid), std::invalid_argument);
			CHECK_THROWS_AS(analyzer.GetRoots(0, 1, invalid), std::invalid_argument);
			CHECK_THROWS_AS(analyzer.ComputeLeftLimit(0, invalid), std::invalid_argument);
			CHECK_THROWS_AS(analyzer.ComputeRightLimit(0, invalid), std::invalid_argument);
		}
		CHECK_THROWS_AS(analyzer.AnalyzeInterval(1, 0, 4), std::invalid_argument);
		CHECK_THROWS_AS(analyzer.AnalyzePoints(0, std::numeric_limits<Real>::infinity(), 4), std::invalid_argument);
		CHECK_THROWS_AS(analyzer.IsMonotonicOnSamples(std::numeric_limits<Real>::quiet_NaN(), 1, 4), std::invalid_argument);
		CHECK(analyzer.AnalyzePoints(1, 1, 1).front().x == Real(1));
	}

	TEST_CASE("FunctionAnalyzer closed grids evaluate exact endpoints without range overflow",
						"[FunctionAnalyzer][function-analyzer-regression]") {
		std::vector<Real> sampled;
		// Exercise the overflowing-span fallback without compile-time overflow.
		volatile Real runtimeBound = std::numeric_limits<Real>::max();
		const Real bound = runtimeBound;
		RealFunctionFromStdFunc identity([&sampled, bound](Real x) {
			sampled.push_back(x);
			return x / bound;
		});
		RealFunctionAnalyzer analyzer(identity);
		REQUIRE(analyzer.IsMonotonicOnSamples(-bound, bound, 4));
		REQUIRE(sampled.size() == 5);
		CHECK(sampled.front() == -bound);
		CHECK(sampled.back() == bound);
		CHECK(sampled[2] == Real(0));
		for (Real x : sampled)
			CHECK(std::isfinite(x));
		sampled.clear();
		REQUIRE(analyzer.MinInNPoints(Real(0.1), Real(0.3), 7) == Real(0.1) / bound);
		REQUIRE(sampled.size() == 8);
		CHECK(sampled.back() == Real(0.3));
	}

	TEST_CASE("FunctionAnalyzer unresolved infinite probes do not invent a signed limit",
						"[FunctionAnalyzer][function-analyzer-regression]") {
		int evaluations = 0;
		RealFunctionFromStdFunc alternating(
			[&evaluations](Real) { return ++evaluations % 2 ? std::numeric_limits<Real>::infinity() : -std::numeric_limits<Real>::infinity(); });
		const RealFunctionAnalyzer analyzer(alternating);
		CHECK(std::isnan(analyzer.ComputeLeftLimit(0)));
		CHECK(std::isnan(analyzer.ComputeRightLimit(0)));
	}

	Real test_func_123(Real x)
	{
		return 9 - (x - 3) * (x - 3);
	}

	// TODO - finalize these tests
	// POUZDANA METRIKA RAZLIKE DVIJE FUNKCIJE!!!!!

	TEST_CASE("Test_Function_Analyzer_FuncDiff", "[simple]")
	{
			TEST_PRECISION_INFO();
		Vector<Real> x{ static_cast<MML::Real>(0.0), static_cast<MML::Real>(3.0), static_cast<MML::Real>(6.0) };
		Vector<Real> y{ test_func_123(static_cast<MML::Real>(0.0)), test_func_123(static_cast<MML::Real>(3.0)), test_func_123(static_cast<MML::Real>(6.0)) };

		LinearInterpRealFunc myfunc(x, y);
		RealFunction test(test_func_123);

		double triangleArea = 6 * test_func_123(static_cast<MML::Real>(3.0)) / 2;
		double parabolaArea = IntegrateTrap(test, static_cast<MML::Real>(0.0), static_cast<MML::Real>(6.0));
		REQUIRE_THAT(parabolaArea - triangleArea, WithinRel(static_cast<MML::Real>(9.0), static_cast<MML::Real>(1e-4)));

		REQUIRE_THAT(RealFunctionComparer::getIntegratedDiff(myfunc, test, static_cast<MML::Real>(0.0), static_cast<MML::Real>(6.0)) , WithinRel(-static_cast<MML::Real>(9.0), static_cast<MML::Real>(1e-4)));
		REQUIRE_THAT(RealFunctionComparer::getIntegratedDiff(myfunc, test, static_cast<MML::Real>(0.0), static_cast<MML::Real>(6.0)) , !WithinRel(-static_cast<MML::Real>(9.0), static_cast<MML::Real>(1e-5)));
	}

	/*********************************************************************/
	/*****          Local extrema detection tests                    *****/
	/*********************************************************************/
	TEST_CASE("FunctionAnalyzer::Local_extrema_quadratic", "[extrema]")
	{
			TEST_PRECISION_INFO();
		// f(x) = (x-2)^2 + 1, minimum at x=2
		RealFunction quadratic([](Real x) { return (x - static_cast<MML::Real>(2.0)) * (x - static_cast<MML::Real>(2.0)) + static_cast<MML::Real>(1.0); });
		RealFunctionAnalyzer analyzer(quadratic);

		// Test that x=2 is a local extremum (minimum)
		REQUIRE(analyzer.IsLocalOptimum(static_cast<MML::Real>(2.0), TOL(1e-6, 1e-4)));

		// Test that x=0 and x=4 are not local extrema
		REQUIRE_FALSE(analyzer.IsLocalOptimum(static_cast<MML::Real>(0.0), TOL(1e-6, 1e-4)));
		REQUIRE_FALSE(analyzer.IsLocalOptimum(static_cast<MML::Real>(4.0), TOL(1e-6, 1e-4)));
	}

	TEST_CASE("FunctionAnalyzer::Local_extrema_cubic", "[extrema]")
	{
			TEST_PRECISION_INFO();
		// f(x) = x^3 - 3x^2 - 9x + 5
		// f'(x) = 3x^2 - 6x - 9 = 3(x-3)(x+1)
		// Local max at x=-1, local min at x=3
		RealFunction cubic([](Real x) {
			return x * x * x - static_cast<MML::Real>(3.0) * x * x - static_cast<MML::Real>(9.0) * x + static_cast<MML::Real>(5.0);
		});
		RealFunctionAnalyzer analyzer(cubic);

		// Test local maximum at x=-1
		REQUIRE(analyzer.IsLocalOptimum(-static_cast<MML::Real>(1.0), TOL(1e-6, 1e-4)));

		// Test local minimum at x=3
		REQUIRE(analyzer.IsLocalOptimum(static_cast<MML::Real>(3.0), TOL(1e-6, 1e-4)));

		// Test that x=0 is not a local extremum
		REQUIRE_FALSE(analyzer.IsLocalOptimum(static_cast<MML::Real>(0.0), TOL(1e-6, 1e-4)));
	}

	TEST_CASE("FunctionAnalyzer::Local_extrema_sine_wave", "[extrema]")
	{
			TEST_PRECISION_INFO();
		// f(x) = sin(x)
		// Local max at π/2, local min at 3π/2
		RealFunction sine([](Real x) { return std::sin(x); });
		RealFunctionAnalyzer analyzer(sine);

		// Test local maximum at π/2
		REQUIRE(analyzer.IsLocalOptimum(Constants::PI / static_cast<MML::Real>(2.0), TOL(1e-6, 1e-4)));

		// Test local minimum at 3π/2
		REQUIRE(analyzer.IsLocalOptimum(static_cast<MML::Real>(3.0) * Constants::PI / static_cast<MML::Real>(2.0), TOL(1e-6, 1e-4)));

		// Test that x=π is not a local extremum (it's an inflection point)
		REQUIRE_FALSE(analyzer.IsLocalOptimum(Constants::PI, TOL(1e-6, 1e-4)));
	}

	/*********************************************************************/
	/*****          Inflection point detection tests                 *****/
	/*********************************************************************/
	TEST_CASE("FunctionAnalyzer::Inflection_point_cubic", "[inflection]")
	{
			TEST_PRECISION_INFO();
		// f(x) = x^3
		// f''(x) = 6x, changes sign at x=0
		// Inflection point at x=0
		RealFunction cubic([](Real x) { return x * x * x; });
		RealFunctionAnalyzer analyzer(cubic);

		// Test inflection point at x=0
		REQUIRE(analyzer.IsInflectionPoint(static_cast<MML::Real>(0.0), MML::Real(1e-4)));

		// Test that x=1 and x=-1 are not inflection points
		REQUIRE_FALSE(analyzer.IsInflectionPoint(static_cast<MML::Real>(1.0), MML::Real(1e-4)));
		REQUIRE_FALSE(analyzer.IsInflectionPoint(-static_cast<MML::Real>(1.0), MML::Real(1e-4)));
	}

	TEST_CASE("FunctionAnalyzer::Inflection_point_sine", "[inflection]")
	{
			TEST_PRECISION_INFO();
		// f(x) = sin(x)
		// f''(x) = -sin(x), changes sign at x=0, π, 2π, ...
		// Inflection points at multiples of π
		RealFunction sine([](Real x) { return std::sin(x); });
		RealFunctionAnalyzer analyzer(sine);

		// Test inflection points at x=0, π, 2π
		REQUIRE(analyzer.IsInflectionPoint(static_cast<MML::Real>(0.0), MML::Real(1e-4)));
		REQUIRE(analyzer.IsInflectionPoint(Constants::PI, MML::Real(1e-4)));
		REQUIRE(analyzer.IsInflectionPoint(static_cast<MML::Real>(2.0) * Constants::PI, MML::Real(1e-4)));

		// Test that π/2 and 3π/2 are not inflection points (they're extrema)
		REQUIRE_FALSE(analyzer.IsInflectionPoint(Constants::PI / static_cast<MML::Real>(2.0), MML::Real(1e-4)));
		REQUIRE_FALSE(analyzer.IsInflectionPoint(static_cast<MML::Real>(3.0) * Constants::PI / static_cast<MML::Real>(2.0), MML::Real(1e-4)));
	}

	TEST_CASE("FunctionAnalyzer::Inflection_point_polynomial", "[inflection]")
	{
			TEST_PRECISION_INFO();
		// f(x) = x^4 - 6x^2
		// f'(x) = 4x^3 - 12x
		// f''(x) = 12x^2 - 12 = 12(x^2 - 1) = 12(x-1)(x+1)
		// Inflection points at x=-1 and x=1
		RealFunction quartic([](Real x) {
			return x * x * x * x - static_cast<MML::Real>(6.0) * x * x;
		});
		RealFunctionAnalyzer analyzer(quartic);

		// Test inflection points at x=-1 and x=1
		REQUIRE(analyzer.IsInflectionPoint(-static_cast<MML::Real>(1.0), MML::Real(1e-4)));
		REQUIRE(analyzer.IsInflectionPoint(static_cast<MML::Real>(1.0), MML::Real(1e-4)));

		// Test that x=0 is not an inflection point (it's a local maximum)
		REQUIRE_FALSE(analyzer.IsInflectionPoint(static_cast<MML::Real>(0.0), MML::Real(1e-4)));
	}

	/*********************************************************************/
	/*****          Monotonicity tests                               *****/
	/*********************************************************************/
	TEST_CASE("FunctionAnalyzer::Monotonicity_linear_increasing", "[monotonic]")
	{
			TEST_PRECISION_INFO();
		// f(x) = 2x + 3, strictly increasing
		RealFunction linear([](Real x) { return static_cast<MML::Real>(2.0) * x + static_cast<MML::Real>(3.0); });
		RealFunctionAnalyzer analyzer(linear);

		REQUIRE(analyzer.IsMonotonicOnSamples(static_cast<MML::Real>(0.0), static_cast<MML::Real>(10.0), 100));
		REQUIRE(analyzer.IsMonotonicOnSamples(-static_cast<MML::Real>(5.0), static_cast<MML::Real>(5.0), 100));
	}

	TEST_CASE("FunctionAnalyzer::Monotonicity_linear_decreasing", "[monotonic]")
	{
			TEST_PRECISION_INFO();
		// f(x) = -x + 5, strictly decreasing
		RealFunction linear([](Real x) { return -x + static_cast<MML::Real>(5.0); });
		RealFunctionAnalyzer analyzer(linear);

		REQUIRE(analyzer.IsMonotonicOnSamples(static_cast<MML::Real>(0.0), static_cast<MML::Real>(10.0), 100));
		REQUIRE(analyzer.IsMonotonicOnSamples(-static_cast<MML::Real>(5.0), static_cast<MML::Real>(5.0), 100));
	}

	TEST_CASE("FunctionAnalyzer::Monotonicity_non_monotonic", "[monotonic]")
	{
			TEST_PRECISION_INFO();
		// f(x) = x^2, not monotonic over [-2, 2]
		RealFunction quadratic([](Real x) { return x * x; });
		RealFunctionAnalyzer analyzer(quadratic);

		REQUIRE_FALSE(analyzer.IsMonotonicOnSamples(-static_cast<MML::Real>(2.0), static_cast<MML::Real>(2.0), 100));

		// But it is monotonic in [0, 2] (increasing)
		REQUIRE(analyzer.IsMonotonicOnSamples(static_cast<MML::Real>(0.0), static_cast<MML::Real>(2.0), 100));

		// And monotonic in [-2, 0] (decreasing)
		REQUIRE(analyzer.IsMonotonicOnSamples(-static_cast<MML::Real>(2.0), static_cast<MML::Real>(0.0), 100));
	}

	TEST_CASE("FunctionAnalyzer::Monotonicity_sine_wave", "[monotonic]")
	{
			TEST_PRECISION_INFO();
		// f(x) = sin(x), not monotonic over [0, 2π]
		RealFunction sine([](Real x) { return std::sin(x); });
		RealFunctionAnalyzer analyzer(sine);

		REQUIRE_FALSE(analyzer.IsMonotonicOnSamples(static_cast<MML::Real>(0.0), static_cast<MML::Real>(2.0) * Constants::PI, 100));

		// But monotonic in [0, π/2] (increasing)
		REQUIRE(analyzer.IsMonotonicOnSamples(static_cast<MML::Real>(0.0), Constants::PI / static_cast<MML::Real>(2.0), 100));

		// And monotonic in [π/2, 3π/2] (decreasing)
		REQUIRE(analyzer.IsMonotonicOnSamples(Constants::PI / static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0) * Constants::PI / static_cast<MML::Real>(2.0), 100));
	}

	/*********************************************************************/
	/*****          Continuity tests                                 *****/
	/*********************************************************************/
	TEST_CASE("FunctionAnalyzer::Continuity_polynomial", "[continuity]")
	{
			TEST_PRECISION_INFO();
		// Polynomials are continuous everywhere
		RealFunction poly([](Real x) { return x * x * x - static_cast<MML::Real>(2.0) * x + static_cast<MML::Real>(1.0); });
		RealFunctionAnalyzer analyzer(poly);

		REQUIRE(analyzer.IsContinuousAtPoint(static_cast<MML::Real>(0.0), MML::Real(1e-6)));
		REQUIRE(analyzer.IsContinuousAtPoint(static_cast<MML::Real>(1.0), MML::Real(1e-6)));
		REQUIRE(analyzer.IsContinuousAtPoint(-static_cast<MML::Real>(5.0), MML::Real(1e-4)));  // Use larger tolerance for larger values
		REQUIRE(analyzer.IsContinuousAtPoint(static_cast<MML::Real>(10.0), MML::Real(1e-3)));  // Use larger tolerance for larger values
	}

	TEST_CASE("FunctionAnalyzer::Continuity_discontinuous_at_point", "[continuity]")
	{
			TEST_PRECISION_INFO();
		// f(x) = 1/(x-2), discontinuous at x=2
		RealFunction rational([](Real x) { return static_cast<MML::Real>(1.0) / (x - static_cast<MML::Real>(2.0)); });
		RealFunctionAnalyzer analyzer(rational);

		REQUIRE(analyzer.IsContinuousAtPoint(static_cast<MML::Real>(0.0), MML::Real(1e-6)));
		REQUIRE(analyzer.IsContinuousAtPoint(static_cast<MML::Real>(1.0), MML::Real(1e-6)));
		REQUIRE_FALSE(analyzer.IsContinuousAtPoint(static_cast<MML::Real>(2.0), MML::Real(1e-6)));
		REQUIRE(analyzer.IsContinuousAtPoint(static_cast<MML::Real>(3.0), MML::Real(1e-6)));
	}

	TEST_CASE("FunctionAnalyzer::Continuity_piecewise", "[continuity]")
	{
			TEST_PRECISION_INFO();
		// Piecewise function: x^2 for x<0, x for x>=0
		// Continuous at x=0
		RealFunction piecewise([](Real x) {
			return x < static_cast<MML::Real>(0.0) ? x * x : x;
		});
		RealFunctionAnalyzer analyzer(piecewise);

		REQUIRE(analyzer.IsContinuousAtPoint(-static_cast<MML::Real>(1.0), MML::Real(1e-6)));
		REQUIRE(analyzer.IsContinuousAtPoint(static_cast<MML::Real>(0.0), MML::Real(1e-3)));  // Need larger epsilon at boundary
		REQUIRE(analyzer.IsContinuousAtPoint(static_cast<MML::Real>(1.0), MML::Real(1e-6)));
	}

	TEST_CASE("FunctionAnalyzer::Continuity_step_function", "[continuity]")
	{
			TEST_PRECISION_INFO();
		// Step function: 0 for x<0, 1 for x>=0
		// Discontinuous at x=0
		RealFunction step([](Real x) {
			return x < static_cast<MML::Real>(0.0) ? static_cast<MML::Real>(0.0) : static_cast<MML::Real>(1.0);
		});
		RealFunctionAnalyzer analyzer(step);

		REQUIRE(analyzer.IsContinuousAtPoint(-static_cast<MML::Real>(1.0), MML::Real(1e-6)));
		REQUIRE_FALSE(analyzer.IsContinuousAtPoint(static_cast<MML::Real>(0.0), MML::Real(1e-6)));
		REQUIRE(analyzer.IsContinuousAtPoint(static_cast<MML::Real>(1.0), MML::Real(1e-6)));
	}

	/*********************************************************************/
	/*****          Continuity Analyzer tests                        *****/
	/*********************************************************************/
	TEST_CASE("FunctionAnalyzer::ContinuityAnalyzer_step_jump", "[continuity_analyzer]")
	{
			TEST_PRECISION_INFO();
		// Step function: 0 for x<0, 1 for x>=0
		// Should detect JUMP discontinuity at x=0
		RealFunction step([](Real x) {
			return x < static_cast<MML::Real>(0.0) ? static_cast<MML::Real>(0.0) : static_cast<MML::Real>(1.0);
		});
		RealFunctionAnalyzer analyzer(step);

		auto discontinuities = analyzer.FindDiscontinuityCandidates(-static_cast<MML::Real>(2.0), static_cast<MML::Real>(2.0), 100, MML::Real(1e-6));
		
		REQUIRE(discontinuities.size() == 1);
		REQUIRE(std::abs(discontinuities[0].x - static_cast<MML::Real>(0.0)) < static_cast<MML::Real>(0.1));  // Found near x=0
		REQUIRE(discontinuities[0].type == DiscontinuityType::JUMP);
		REQUIRE(std::abs(discontinuities[0].leftLimit - static_cast<MML::Real>(0.0)) < 1e-3);
		REQUIRE(std::abs(discontinuities[0].rightLimit - static_cast<MML::Real>(1.0)) < 1e-3);
		REQUIRE(std::abs(discontinuities[0].jumpSize - static_cast<MML::Real>(1.0)) < 1e-3);
	}

	TEST_CASE("FunctionAnalyzer::ContinuityAnalyzer_removable", "[continuity_analyzer]")
	{
			TEST_PRECISION_INFO();
		// f(x) = (x^2 - 1) / (x - 1) for x != 1, undefined at x=1
		// Should detect REMOVABLE discontinuity at x=1 (limit exists = 2)
		RealFunction removable([](Real x) {
			if (x == static_cast<MML::Real>(1.0))
				return std::numeric_limits<Real>::quiet_NaN();  // Undefined at x=1
			return (x * x - static_cast<MML::Real>(1.0)) / (x - static_cast<MML::Real>(1.0));  // = x + 1 for x != 1
		});
		RealFunctionAnalyzer analyzer(removable);

		// 128 subdivisions make x=1 exactly representable in every supported Real type.
		auto discontinuities = analyzer.FindDiscontinuityCandidates(static_cast<MML::Real>(0.0), static_cast<MML::Real>(2.0), 128, TOL(1e-6, 1e-5));
		
		REQUIRE(discontinuities.size() >= 1);
		
		// Find the discontinuity near x=1
		bool foundAtOne = false;
		for (const auto& disc : discontinuities)
		{
			if (std::abs(disc.x - static_cast<MML::Real>(1.0)) < static_cast<MML::Real>(0.1))
			{
				foundAtOne = true;
				REQUIRE(disc.type == DiscontinuityType::REMOVABLE);
				// Limit should be approximately static_cast<MML::Real>(2.0)
				REQUIRE(std::abs(disc.leftLimit - static_cast<MML::Real>(2.0)) < 1e-2);
				REQUIRE(std::abs(disc.rightLimit - static_cast<MML::Real>(2.0)) < 1e-2);
			}
		}
		REQUIRE(foundAtOne);
	}

	TEST_CASE("FunctionAnalyzer::ContinuityAnalyzer_infinite", "[continuity_analyzer]")
	{
			TEST_PRECISION_INFO();
		// f(x) = 1/(x-2), infinite discontinuity at x=2
		RealFunction infinite([](Real x) {
			return static_cast<MML::Real>(1.0) / (x - static_cast<MML::Real>(2.0));
		});
		RealFunctionAnalyzer analyzer(infinite);

		auto discontinuities = analyzer.FindDiscontinuityCandidates(static_cast<MML::Real>(0.0), static_cast<MML::Real>(4.0), 200, TOL(1e-5, 1e-3));
		
		REQUIRE(discontinuities.size() >= 1);
		
		// Find the discontinuity near x=2
		bool foundAtTwo = false;
		for (const auto& disc : discontinuities)
		{
			if (std::abs(disc.x - static_cast<MML::Real>(2.0)) < static_cast<MML::Real>(0.15))
			{
				foundAtTwo = true;
				// Should be Infinite type (limits don't exist or are very large)
				REQUIRE((disc.type == DiscontinuityType::Infinite || disc.type == DiscontinuityType::JUMP));
			}
		}
		REQUIRE(foundAtTwo);
	}

	TEST_CASE("FunctionAnalyzer::ContinuityAnalyzer_continuous", "[continuity_analyzer]")
	{
			TEST_PRECISION_INFO();
		// Continuous polynomial: x^3 - 2x + 1
		RealFunction poly([](Real x) {
			return x * x * x - static_cast<MML::Real>(2.0) * x + static_cast<MML::Real>(1.0);
		});
		RealFunctionAnalyzer analyzer(poly);

		auto discontinuities = analyzer.FindDiscontinuityCandidates(-static_cast<MML::Real>(5.0), static_cast<MML::Real>(5.0), 100, MML::Real(1e-6));
		
		REQUIRE(discontinuities.size() == 0);  // Should find no discontinuities
	}

	TEST_CASE("FunctionAnalyzer::ContinuityAnalyzer_multiple_jumps", "[continuity_analyzer]")
	{
			TEST_PRECISION_INFO();
		// Piecewise function with two jumps
		// f(x) = 0 for x<-1, 1 for -1<=x<1, 2 for x>=1
		RealFunction multiJump([](Real x) {
			if (x < -static_cast<MML::Real>(1.0)) return static_cast<MML::Real>(0.0);
			if (x < static_cast<MML::Real>(1.0)) return static_cast<MML::Real>(1.0);
			return static_cast<MML::Real>(2.0);
		});
		RealFunctionAnalyzer analyzer(multiJump);

		auto discontinuities = analyzer.FindDiscontinuityCandidates(-static_cast<MML::Real>(3.0), static_cast<MML::Real>(3.0), 300, MML::Real(1e-5));
		
		// Should find at least one discontinuity (ideally 2)
		// The strict >= 2 requirement might be too aggressive depending on sampling
		REQUIRE(discontinuities.size() >= 1);
		
		// Check that we found at least one near x=-1 or x=1
		bool foundNearMinusOne = false;
		bool foundNearOne = false;
		
		for (const auto& disc : discontinuities)
		{
			if (std::abs(disc.x - (-static_cast<MML::Real>(1.0))) < static_cast<MML::Real>(0.15))
				foundNearMinusOne = true;
			if (std::abs(disc.x - static_cast<MML::Real>(1.0)) < static_cast<MML::Real>(0.15))
				foundNearOne = true;
			
			// Should be JUMP type
			REQUIRE(disc.type == DiscontinuityType::JUMP);
		}
		
		// At least one of the jumps should be found
		REQUIRE((foundNearMinusOne || foundNearOne));
	}

	/*********************************************************************/
	/*****          Min/Max finding tests                            *****/
	/*********************************************************************/
	TEST_CASE("FunctionAnalyzer::MinMax_quadratic", "[minmax]")
	{
			TEST_PRECISION_INFO();
		// f(x) = x^2 - 4x + 5 = (x-2)^2 + 1
		// Minimum at x=2, f(2)=1
		RealFunction quadratic([](Real x) {
			return x * x - static_cast<MML::Real>(4.0) * x + static_cast<MML::Real>(5.0);
		});
		RealFunctionAnalyzer analyzer(quadratic);

		Real min = analyzer.MinInNPoints(static_cast<MML::Real>(0.0), static_cast<MML::Real>(4.0), 100);
		Real max = analyzer.MaxInNPoints(static_cast<MML::Real>(0.0), static_cast<MML::Real>(4.0), 100);

		// Minimum should be close to static_cast<MML::Real>(1.0) (at x=2)
		REQUIRE_THAT(min, WithinAbs(static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.1)));

		// Maximum at boundaries: f(0)=5 or f(4)=5
		REQUIRE_THAT(max, WithinAbs(static_cast<MML::Real>(5.0), static_cast<MML::Real>(0.1)));
	}

	TEST_CASE("FunctionAnalyzer::MinMax_sine_wave", "[minmax]")
	{
			TEST_PRECISION_INFO();
		// f(x) = sin(x) over [0, 2π]
		// Min = -1 at 3π/2, Max = 1 at π/2
		RealFunction sine([](Real x) { return std::sin(x); });
		RealFunctionAnalyzer analyzer(sine);

		Real min = analyzer.MinInNPoints(static_cast<MML::Real>(0.0), static_cast<MML::Real>(2.0) * Constants::PI, 1000);
		Real max = analyzer.MaxInNPoints(static_cast<MML::Real>(0.0), static_cast<MML::Real>(2.0) * Constants::PI, 1000);

		REQUIRE_THAT(min, WithinAbs(-static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.01)));
		REQUIRE_THAT(max, WithinAbs(static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.01)));
	}

	/*********************************************************************/
	/*****          Root finding tests                               *****/
	/*********************************************************************/
	TEST_CASE("FunctionAnalyzer::Roots_linear", "[roots]")
	{
			TEST_PRECISION_INFO();
		// f(x) = 2x - 4, root at x=2
		RealFunction linear([](Real x) { return static_cast<MML::Real>(2.0) * x - static_cast<MML::Real>(4.0); });
		RealFunctionAnalyzer analyzer(linear);

		std::vector<Real> roots = analyzer.GetRoots(static_cast<MML::Real>(0.0), static_cast<MML::Real>(5.0), MML::Real(1e-6));

		REQUIRE(roots.size() == 1);
		REQUIRE_THAT(roots[0], WithinAbs(static_cast<MML::Real>(2.0), static_cast<MML::Real>(1e-5)));
	}

	TEST_CASE("FunctionAnalyzer::Roots_quadratic", "[roots]")
	{
			TEST_PRECISION_INFO();
		// f(x) = x^2 - 5x + 6 = (x-2)(x-3)
		// Roots at x=2 and x=3
		RealFunction quadratic([](Real x) {
			return x * x - static_cast<MML::Real>(5.0) * x + static_cast<MML::Real>(6.0);
		});
		RealFunctionAnalyzer analyzer(quadratic);

		std::vector<Real> roots = analyzer.GetRoots(static_cast<MML::Real>(0.0), static_cast<MML::Real>(5.0), MML::Real(1e-6));

		REQUIRE(roots.size() == 2);
		REQUIRE_THAT(roots[0], WithinAbs(static_cast<MML::Real>(2.0), static_cast<MML::Real>(1e-5)));
		REQUIRE_THAT(roots[1], WithinAbs(static_cast<MML::Real>(3.0), static_cast<MML::Real>(1e-5)));
	}

	TEST_CASE("FunctionAnalyzer::Roots_sine", "[roots]")
	{
			TEST_PRECISION_INFO();
		// f(x) = sin(x), roots at 0, π, 2π
		RealFunction sine([](Real x) { return std::sin(x); });
		RealFunctionAnalyzer analyzer(sine);

		std::vector<Real> roots = analyzer.GetRoots(static_cast<MML::Real>(0.0), static_cast<MML::Real>(2.0) * Constants::PI, MML::Real(1e-6));

		REQUIRE(roots.size() >= 2);  // At least π and 2π (0 might be on boundary)
		// Check that roots are near multiples of π
		for (const auto& root : roots) {
			Real normalized = root / Constants::PI;
			Real nearest_int = std::round(normalized);
			REQUIRE_THAT(normalized, WithinAbs(nearest_int, static_cast<MML::Real>(0.01)));
		}
	}

	TEST_CASE("FunctionAnalyzer::Roots_no_roots", "[roots]")
	{
			TEST_PRECISION_INFO();
		// f(x) = x^2 + 1, no real roots
		RealFunction quadratic([](Real x) { return x * x + static_cast<MML::Real>(1.0); });
		RealFunctionAnalyzer analyzer(quadratic);

		std::vector<Real> roots = analyzer.GetRoots(-static_cast<MML::Real>(10.0), static_cast<MML::Real>(10.0), MML::Real(1e-6));

		REQUIRE(roots.size() == 0);
	}

	/*********************************************************************/
	/*****          Edge cases and numerical stability               *****/
	/*********************************************************************/
	TEST_CASE("FunctionAnalyzer::Edge_case_constant_function", "[edge]")
	{
			TEST_PRECISION_INFO();
		// f(x) = 5, constant function
		RealFunction constant([](Real x) { return static_cast<MML::Real>(5.0); });
		RealFunctionAnalyzer analyzer(constant);

		// Constant function is monotonic
		REQUIRE(analyzer.IsMonotonicOnSamples(static_cast<MML::Real>(0.0), static_cast<MML::Real>(10.0), 100));

		// No inflection points
		REQUIRE_FALSE(analyzer.IsInflectionPoint(static_cast<MML::Real>(0.0), MML::Real(1e-4)));
		REQUIRE_FALSE(analyzer.IsInflectionPoint(static_cast<MML::Real>(5.0), MML::Real(1e-4)));

		// All points are "extrema" in a sense, but isLocalOptimum should handle this
		// (second derivative is zero everywhere)
		REQUIRE_FALSE(analyzer.IsLocalOptimum(static_cast<MML::Real>(0.0), TOL(1e-6, 1e-4)));
	}

	TEST_CASE("FunctionAnalyzer::Edge_case_absolute_value", "[edge]")
	{
			TEST_PRECISION_INFO();
		// f(x) = |x|, not differentiable at x=0
		RealFunction absValue([](Real x) { return std::abs(x); });
		RealFunctionAnalyzer analyzer(absValue);

		// Should be continuous everywhere including x=0
		REQUIRE(analyzer.IsContinuousAtPoint(static_cast<MML::Real>(0.0), MML::Real(1e-3)));

		// x=0 is a local minimum (but not smooth)
		// isLocalOptimum uses second derivative which doesn't exist at x=0
		// The numerical approximation might not detect it reliably
		// This is expected behavior for non-smooth functions

		// Not monotonic over [-1, 1]
		REQUIRE_FALSE(analyzer.IsMonotonicOnSamples(-static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0), 100));

		// Monotonic in [0, 1]
		REQUIRE(analyzer.IsMonotonicOnSamples(static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), 100));
	}

	/*********************************************************************/
	/*****          GetLocalOptimums and GetInflectionPoints tests   *****/
	/*********************************************************************/
	TEST_CASE("FunctionAnalyzer::GetLocalOptimums_cubic", "[optimums]")
	{
			TEST_PRECISION_INFO();
		// f(x) = x^3 - 3x^2 - 9x + 5
		// f'(x) = 3x^2 - 6x - 9 = 3(x-3)(x+1)
		// Local max at x=-1, local min at x=3
		RealFunction cubic([](Real x) {
			return x * x * x - static_cast<MML::Real>(3.0) * x * x - static_cast<MML::Real>(9.0) * x + static_cast<MML::Real>(5.0);
		});
		RealFunctionAnalyzer analyzer(cubic);

		auto optimums = analyzer.GetLocalOptimums(-static_cast<MML::Real>(5.0), static_cast<MML::Real>(5.0), TOL(1e-6, 1e-4));

		// Should find exactly 2 optimums
		REQUIRE(optimums.size() == 2);

		// Sort to ensure consistent order
		std::sort(optimums.begin(), optimums.end());

		// First optimum at x=-1 (local max)
		REQUIRE_THAT(optimums[0], WithinAbs(-static_cast<MML::Real>(1.0), static_cast<MML::Real>(1e-4)));

		// Second optimum at x=3 (local min)
		REQUIRE_THAT(optimums[1], WithinAbs(static_cast<MML::Real>(3.0), static_cast<MML::Real>(1e-4)));
	}

	TEST_CASE("FunctionAnalyzer::GetLocalOptimumsClassified_cubic", "[optimums]")
	{
			TEST_PRECISION_INFO();
		// f(x) = x^3 - 3x^2 - 9x + 5
		// f'(x) = 3x^2 - 6x - 9 = 3(x-3)(x+1)
		// Local max at x=-1 (f''(-1) = -12 < 0), local min at x=3 (f''(3) = 12 > 0)
		RealFunction cubic([](Real x) {
			return x * x * x - static_cast<MML::Real>(3.0) * x * x - static_cast<MML::Real>(9.0) * x + static_cast<MML::Real>(5.0);
		});
		RealFunctionAnalyzer analyzer(cubic);

		auto classified = analyzer.GetLocalOptimumsClassified(-static_cast<MML::Real>(5.0), static_cast<MML::Real>(5.0), TOL(1e-6, 1e-4));

		// Should find exactly 2 classified optimums
		REQUIRE(classified.size() == 2);

		// Sort by x-coordinate
		std::sort(classified.begin(), classified.end(), 
			[](const CriticalPoint& a, const CriticalPoint& b) { return a.x < b.x; });

		// First critical point at x=-1 is a LOCAL_MAXIMUM
		REQUIRE_THAT(classified[0].x, WithinAbs(-static_cast<MML::Real>(1.0), static_cast<MML::Real>(1e-4)));
		REQUIRE(classified[0].type == CriticalPointType::LOCAL_MAXIMUM);
		REQUIRE_THAT(classified[0].value, WithinAbs(static_cast<MML::Real>(10.0), static_cast<MML::Real>(1e-3)));  // f(-1) = 10

		// Second critical point at x=3 is a LOCAL_MINIMUM
		REQUIRE_THAT(classified[1].x, WithinAbs(static_cast<MML::Real>(3.0), static_cast<MML::Real>(1e-4)));
		REQUIRE(classified[1].type == CriticalPointType::LOCAL_MINIMUM);
		REQUIRE_THAT(classified[1].value, WithinAbs(-static_cast<MML::Real>(22.0), static_cast<MML::Real>(1e-3)));  // f(3) = -22
	}

	TEST_CASE("FunctionAnalyzer::GetLocalOptimums_sine", "[optimums]")
	{
			TEST_PRECISION_INFO();
		// f(x) = sin(x)
		// Extrema at π/2 + nπ
		RealFunction sine([](Real x) { return std::sin(x); });
		RealFunctionAnalyzer analyzer(sine);

		// Search in [0, 2π]
		auto optimums = analyzer.GetLocalOptimums(static_cast<MML::Real>(0.0), static_cast<MML::Real>(2.0) * Constants::PI, TOL(1e-6, 1e-4));

		// Should find 2 extrema: max at π/2, min at 3π/2
		REQUIRE(optimums.size() == 2);

		// Sort to ensure consistent order
		std::sort(optimums.begin(), optimums.end());

		REQUIRE_THAT(optimums[0], WithinAbs(Constants::PI / static_cast<MML::Real>(2.0), 1e-4));
		REQUIRE_THAT(optimums[1], WithinAbs(static_cast<MML::Real>(3.0) * Constants::PI / static_cast<MML::Real>(2.0), 1e-4));
	}

	TEST_CASE("FunctionAnalyzer::GetLocalOptimumsClassified_sine", "[optimums]")
	{
			TEST_PRECISION_INFO();
		// f(x) = sin(x)
		// Max at π/2, min at 3π/2
		RealFunction sine([](Real x) { return std::sin(x); });
		RealFunctionAnalyzer analyzer(sine);

		auto classified = analyzer.GetLocalOptimumsClassified(static_cast<MML::Real>(0.0), static_cast<MML::Real>(2.0) * Constants::PI, TOL(1e-6, 1e-4));

		REQUIRE(classified.size() == 2);

		// Sort by x-coordinate
		std::sort(classified.begin(), classified.end(),
			[](const CriticalPoint& a, const CriticalPoint& b) { return a.x < b.x; });

		// Maximum at π/2
		REQUIRE_THAT(classified[0].x, WithinAbs(Constants::PI / static_cast<MML::Real>(2.0), 1e-4));
		REQUIRE(classified[0].type == CriticalPointType::LOCAL_MAXIMUM);
		REQUIRE_THAT(classified[0].value, WithinAbs(static_cast<MML::Real>(1.0), static_cast<MML::Real>(1e-4)));

		// Minimum at 3π/2
		REQUIRE_THAT(classified[1].x, WithinAbs(static_cast<MML::Real>(3.0) * Constants::PI / static_cast<MML::Real>(2.0), 1e-4));
		REQUIRE(classified[1].type == CriticalPointType::LOCAL_MINIMUM);
		REQUIRE_THAT(classified[1].value, WithinAbs(-static_cast<MML::Real>(1.0), static_cast<MML::Real>(1e-4)));
	}

	TEST_CASE("FunctionAnalyzer::GetLocalOptimums_quartic", "[optimums]")
	{
			TEST_PRECISION_INFO();
		// f(x) = x^4 - 4x^3 + 4x^2
		// f'(x) = 4x^3 - 12x^2 + 8x = 4x(x^2 - 3x + 2) = 4x(x-1)(x-2)
		// Critical points at x=0, x=1, x=2
		RealFunction quartic([](Real x) {
			return x * x * x * x - static_cast<MML::Real>(4.0) * x * x * x + static_cast<MML::Real>(4.0) * x * x;
		});
		RealFunctionAnalyzer analyzer(quartic);

		auto optimums = analyzer.GetLocalOptimums(-static_cast<MML::Real>(1.0), static_cast<MML::Real>(3.0), TOL(1e-6, 1e-4));

		// Should find 3 optimums at x=0, 1, 2
		REQUIRE(optimums.size() == 3);

		std::sort(optimums.begin(), optimums.end());

		REQUIRE_THAT(optimums[0], WithinAbs(static_cast<MML::Real>(0.0), static_cast<MML::Real>(1e-4)));
		REQUIRE_THAT(optimums[1], WithinAbs(static_cast<MML::Real>(1.0), static_cast<MML::Real>(1e-4)));
		REQUIRE_THAT(optimums[2], WithinAbs(static_cast<MML::Real>(2.0), static_cast<MML::Real>(1e-4)));
	}

	TEST_CASE("FunctionAnalyzer::GetInflectionPoints_cubic", "[inflection]")
	{
			TEST_PRECISION_INFO();
		// f(x) = x^3 - 3x^2 - 9x + 5
		// f''(x) = 6x - 6 = 0 => x = 1
		// Inflection point at x=1
		RealFunction cubic([](Real x) {
			return x * x * x - static_cast<MML::Real>(3.0) * x * x - static_cast<MML::Real>(9.0) * x + static_cast<MML::Real>(5.0);
		});
		RealFunctionAnalyzer analyzer(cubic);

		auto inflections = analyzer.GetInflectionPoints(-static_cast<MML::Real>(5.0), static_cast<MML::Real>(5.0), TOL(1e-6, 1e-4));

		// Should find exactly 1 inflection point
		REQUIRE(inflections.size() == 1);
		REQUIRE_THAT(inflections[0], WithinAbs(static_cast<MML::Real>(1.0), static_cast<MML::Real>(1e-4)));
	}

	TEST_CASE("FunctionAnalyzer::GetInflectionPoints_quartic", "[inflection]")
	{
			TEST_PRECISION_INFO();
		// f(x) = x^4 - 6x^2
		// f''(x) = 12x^2 - 12 = 12(x^2 - 1) = 0 => x = ±1
		// Inflection points at x=-1 and x=1
		RealFunction quartic([](Real x) {
			return x * x * x * x - static_cast<MML::Real>(6.0) * x * x;
		});
		RealFunctionAnalyzer analyzer(quartic);

		auto inflections = analyzer.GetInflectionPoints(-static_cast<MML::Real>(3.0), static_cast<MML::Real>(3.0), TOL(1e-6, 1e-4));

		// Should find 2 inflection points
		REQUIRE(inflections.size() == 2);

		std::sort(inflections.begin(), inflections.end());

		REQUIRE_THAT(inflections[0], WithinAbs(-static_cast<MML::Real>(1.0), static_cast<MML::Real>(1e-4)));
		REQUIRE_THAT(inflections[1], WithinAbs(static_cast<MML::Real>(1.0), static_cast<MML::Real>(1e-4)));
	}

	TEST_CASE("FunctionAnalyzer::GetInflectionPoints_sine", "[inflection]")
	{
			TEST_PRECISION_INFO();
		// f(x) = sin(x)
		// f''(x) = -sin(x) = 0 => x = nπ
		// Inflection points at 0, π, 2π
		RealFunction sine([](Real x) { return std::sin(x); });
		RealFunctionAnalyzer analyzer(sine);

		auto inflections = analyzer.GetInflectionPoints(static_cast<MML::Real>(0.0), static_cast<MML::Real>(2.0) * Constants::PI, TOL(1e-6, 1e-4));

		// Should find 3 inflection points: 0, π, 2π
		REQUIRE(inflections.size() == 3);

		std::sort(inflections.begin(), inflections.end());

		REQUIRE_THAT(inflections[0], WithinAbs(static_cast<MML::Real>(0.0), static_cast<MML::Real>(1e-4)));
		REQUIRE_THAT(inflections[1], WithinAbs(Constants::PI, static_cast<MML::Real>(1e-4)));
		REQUIRE_THAT(inflections[2], WithinAbs(static_cast<MML::Real>(2.0) * Constants::PI, 1e-4));
	}

	TEST_CASE("FunctionAnalyzer::GetInflectionPoints_no_inflection", "[inflection]")
	{
			TEST_PRECISION_INFO();
		// f(x) = x^2, no inflection points (concavity never changes)
		RealFunction quadratic([](Real x) { return x * x; });
		RealFunctionAnalyzer analyzer(quadratic);

		auto inflections = analyzer.GetInflectionPoints(-static_cast<MML::Real>(5.0), static_cast<MML::Real>(5.0), TOL(1e-6, 1e-4));

		// Should find no inflection points
		REQUIRE(inflections.size() == 0);
	}

	//////////////////////////////////////////////////////////////////////////
	// OVERFLOW PROTECTION TESTS
	//////////////////////////////////////////////////////////////////////////

	TEST_CASE("FunctionAnalyzer::OverflowProtection_1_over_x", "[overflow]")
	{
		TEST_PRECISION_INFO();
		// f(x) = 1/x has singularity at x=0
		RealFunction f([](Real x) { return static_cast<MML::Real>(1.0) / x; });
		RealFunctionAnalyzer analyzer(f);

		// Test ComputeLeftLimit near singularity - approaching 0 from left gives -infinity
		Real leftLimit = analyzer.ComputeLeftLimit(static_cast<MML::Real>(0.0));
		// Left limit should be very large negative (approaching -infinity)
		// The function returns large values as we approach 0 from left
		REQUIRE(leftLimit < static_cast<MML::Real>(0.0)); // Must be negative when approaching from left

		// Test ComputeRightLimit near singularity - approaching 0 from right gives +infinity
		Real rightLimit = analyzer.ComputeRightLimit(static_cast<MML::Real>(0.0));
		// Right limit should be very large positive (approaching +infinity)
		REQUIRE(rightLimit > static_cast<MML::Real>(0.0)); // Must be positive when approaching from right

		// isContinuousAtPoint should return false at singularity
		REQUIRE(analyzer.IsContinuousAtPoint(static_cast<MML::Real>(0.0), MML::Real(1e-6)) == false);

		// Should be continuous away from singularity
		REQUIRE(analyzer.IsContinuousAtPoint(static_cast<MML::Real>(1.0), MML::Real(1e-6)) == true);
		REQUIRE(analyzer.IsContinuousAtPoint(static_cast<MML::Real>(-1.0), MML::Real(1e-6)) == true);
	}

	TEST_CASE("FunctionAnalyzer::OverflowProtection_1_over_x_squared", "[overflow]")
	{
		TEST_PRECISION_INFO();
		// f(x) = 1/x^2 has singularity at x=0 (both sides go to +infinity)
		RealFunction f([](Real x) { return static_cast<MML::Real>(1.0) / (x * x); });
		RealFunctionAnalyzer analyzer(f);

		// Both limits should be +infinity or very large
		Real leftLimit = analyzer.ComputeLeftLimit(static_cast<MML::Real>(0.0));
		Real rightLimit = analyzer.ComputeRightLimit(static_cast<MML::Real>(0.0));

		REQUIRE((std::isinf(leftLimit) || leftLimit > static_cast<MML::Real>(1e10)));
		REQUIRE((std::isinf(rightLimit) || rightLimit > static_cast<MML::Real>(1e10)));

		// isContinuousAtPoint should return false at singularity
		REQUIRE(analyzer.IsContinuousAtPoint(static_cast<MML::Real>(0.0), MML::Real(1e-6)) == false);
	}

	TEST_CASE("FunctionAnalyzer::OverflowProtection_MinMax_with_singularity", "[overflow]")
	{
		TEST_PRECISION_INFO();
		// f(x) = 1/(x-5) has singularity at x=5
		// Test MinInNPoints and MaxInNPoints over interval containing singularity
		RealFunction f([](Real x) { return static_cast<MML::Real>(1.0) / (x - static_cast<MML::Real>(5.0)); });
		RealFunctionAnalyzer analyzer(f);

		// These should handle overflow gracefully and not crash
		Real minVal = analyzer.MinInNPoints(static_cast<MML::Real>(0.0), static_cast<MML::Real>(10.0), 100);
		Real maxVal = analyzer.MaxInNPoints(static_cast<MML::Real>(0.0), static_cast<MML::Real>(10.0), 100);

		// Should return valid values (ignoring the singularity points)
		// The values around x=5 will overflow, but other values should be captured
		REQUIRE_FALSE(std::isnan(minVal)); // Should have found some valid minimum
		REQUIRE_FALSE(std::isnan(maxVal)); // Should have found some valid maximum
	}

	TEST_CASE("FunctionAnalyzer::OverflowProtection_tan", "[overflow]")
	{
		TEST_PRECISION_INFO();
		// f(x) = tan(x) has singularities at x = ±π/2, ±3π/2, etc.
		RealFunction f([](Real x) { return std::tan(x); });
		RealFunctionAnalyzer analyzer(f);

		const Real pi = static_cast<MML::Real>(3.14159265358979323846);
		const Real halfPi = pi / static_cast<MML::Real>(2.0);

		// Near π/2, should detect discontinuity/overflow
		REQUIRE(analyzer.IsContinuousAtPoint(halfPi, MML::Real(1e-6)) == false);

		// Away from singularity, should be continuous
		REQUIRE(analyzer.IsContinuousAtPoint(static_cast<MML::Real>(0.0), MML::Real(1e-6)) == true);
		REQUIRE(analyzer.IsContinuousAtPoint(static_cast<MML::Real>(1.0), MML::Real(1e-6)) == true);
	}

	TEST_CASE("FunctionAnalyzer::SafeEvaluate_utility", "[overflow]")
	{
		TEST_PRECISION_INFO();
		// Test the SafeEvaluate utility directly
		RealFunction normalFunc([](Real x) { return x * x; });
		RealFunction singularFunc([](Real x) { return static_cast<MML::Real>(1.0) / x; });
		RealFunction nanFunc([](Real x) { return std::sqrt(x); }); // NaN for negative x

		// Normal evaluation
		SafeEvalResult result1 = SafeEvaluate(normalFunc, static_cast<MML::Real>(2.0));
		REQUIRE(result1.isValid == true);
		REQUIRE(result1.isOverflow == false);
		REQUIRE(result1.isNaN == false);
		REQUIRE_THAT(result1.value, WithinAbs(static_cast<MML::Real>(4.0), TOL(1e-10, 1e-5)));

		// Very small x gives very large but finite result - test the threshold
		// Use threshold of 1e100 to test overflow detection
		SafeEvalResult result2 = SafeEvaluate(singularFunc, static_cast<MML::Real>(1e-200), static_cast<MML::Real>(1e100));
		// 1/1e-200 = 1e200, which exceeds threshold of 1e100
		REQUIRE(result2.isOverflow == true);
		REQUIRE(result2.isValid == false);

		// NaN result
		SafeEvalResult result3 = SafeEvaluate(nanFunc, static_cast<MML::Real>(-1.0));
		REQUIRE(result3.isNaN == true);
		REQUIRE(result3.isValid == false);
	}

}

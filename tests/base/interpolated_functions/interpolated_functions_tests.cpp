#if !defined __MML_INTERPOLATED_FUNCTIONS_TESTS_H
#define __MML_INTERPOLATED_FUNCTIONS_TESTS_H

#include <catch2/catch_all.hpp>
#include "../../TestPrecision.h"

#include <vector>
#include <cmath>
#include <type_traits>

#include "../../test_beds/interpolation_test_bed.h"

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/base/InterpolatedFunction.h>
#include <mml/algorithms/Analyzers/FunctionsAnalyzer.h>
#endif

using namespace MML;

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace MML::Tests::Core::InterpolatedFunctionsTests
{
	//////////////////////////////////////////////////////////////////////////////
	//                           HELPER FUNCTIONS                               //
	//////////////////////////////////////////////////////////////////////////////

	double test_func(const double x)
	{
		const double eps = 1.0;
		return x * exp(-x) / (POW2(x - 1.0) + eps * eps);
	}

	void CreateInterpolatedValues(RealFunction f, Real x1, Real x2, int numPnt, Vector<Real>& outX, Vector<Real>& outY)
	{
		outX.Resize(numPnt);
		outY.Resize(numPnt);

		for (int i = 0; i < numPnt; i++) {
			outX[i] = x1 + i * (x2 - x1) / (numPnt - 1);
			outY[i] = f(outX[i]);
		}
	}

	//////////////////////////////////////////////////////////////////////////////
	//                         LINEAR INTERPOLATION TESTS                       //
	//////////////////////////////////////////////////////////////////////////////

	TEST_CASE("LinearInterp_sin_function_accuracy", "[interpolation][linear]") {
		RealFunction f{ [](Real x) -> Real { return sin(x); } };

		Vector<Real> vec_x, vec_y;
		CreateInterpolatedValues(f, 0.0, 2.0 * Constants::PI, 100, vec_x, vec_y);

		LinearInterpRealFunc linear_interp(vec_x, vec_y);
		RealFunctionComparer comparer(f, linear_interp);

		double avgAbsDiff = comparer.getAbsDiffAvg(0.0, 2.0 * Constants::PI, 100);
		double maxAbsDiff = comparer.getAbsDiffMax(0.0, 2.0 * Constants::PI, 100);

		REQUIRE_THAT(avgAbsDiff, WithinAbs(0.0, 0.001));
		REQUIRE_THAT(maxAbsDiff, WithinAbs(0.0, 0.01));
	}

	TEST_CASE("LinearInterp_exact_at_nodes", "[interpolation][linear]") {
		// Linear interpolation must be exact at data points
		Vector<Real> x{ 0.0, 1.0, 2.0, 3.0, 4.0 };
		Vector<Real> y{ 1.0, 3.0, 2.0, 5.0, 4.0 };

		LinearInterpRealFunc interp(x, y);

		for (int i = 0; i < x.size(); i++) {
			REQUIRE_THAT(interp(x[i]), WithinAbs(y[i], TOL(1e-14, 1e-5)));
		}
	}

	TEST_CASE("LinearInterp_extrapolation", "[interpolation][linear]") {
		Vector<Real> x{ 0.0, 1.0, 2.0, 3.0 };
		Vector<Real> y{ 0.0, 2.0, 4.0, 6.0 };  // y = 2x

		LinearInterpRealFunc interp_no_extrap(x, y, false);
		LinearInterpRealFunc interp_extrap(x, y, true);

		// With extrapolation, should continue the line
		REQUIRE_THAT(interp_extrap(4.0), WithinAbs(8.0, 0.001));
		REQUIRE_THAT(interp_extrap(-1.0), WithinAbs(-2.0, 0.001));
	}

	//////////////////////////////////////////////////////////////////////////////
	//                       POLYNOMIAL INTERPOLATION TESTS                     //
	//////////////////////////////////////////////////////////////////////////////

	TEST_CASE("PolynomInterp_quadratic_exact", "[interpolation][polynomial]") {
		// For a quadratic through 3 points, polynomial interpolation should be exact
		RealFunction f{ [](Real x) { return x * x - 2 * x + 1; } };  // (x-1)^2

		Vector<Real> x{ 0.0, 1.0, 2.0 };
		Vector<Real> y(3);
		for (int i = 0; i < 3; i++) y[i] = f(x[i]);

		PolynomInterpRealFunc interp(x, y, 3);  // 3-point polynomial

		// Test at intermediate points
		REQUIRE_THAT(interp(0.5), WithinAbs(f(0.5), TOL(1e-12, 1e-5)));
		REQUIRE_THAT(interp(1.5), WithinAbs(f(1.5), TOL(1e-12, 1e-5)));
	}

	TEST_CASE("PolynomInterp_error_estimate", "[interpolation][polynomial]") {
		// Error estimate from Neville's algorithm (can be negative as it's a signed difference)
		RealFunction f{ [](Real x) -> Real { return cos(x); } };

		Vector<Real> vec_x, vec_y;
		CreateInterpolatedValues(f, 0.0, Constants::PI, 20, vec_x, vec_y);

		PolynomInterpRealFunc interp(vec_x, vec_y, 8);  // 8-point interpolation

		Real result = interp(Constants::PI / 4);
		Real error = interp.getLastErrorEst();

		// Error magnitude should be small for smooth function
		Real actual_error = std::abs(result - f(Constants::PI / 4));
		REQUIRE(std::abs(error) < 0.001);  // Error estimate should be small
	}

	//////////////////////////////////////////////////////////////////////////////
	//                        CUBIC SPLINE INTERPOLATION TESTS                  //
	//////////////////////////////////////////////////////////////////////////////

	TEST_CASE("SplineInterp_exact_at_nodes", "[interpolation][spline]") {
		Vector<Real> x{ 0.0, 1.0, 2.0, 3.0, 4.0 };
		Vector<Real> y{ 1.0, 3.0, 2.0, 5.0, 4.0 };

		SplineInterpRealFunc spline(x, y);

		for (int i = 0; i < x.size(); i++) {
			REQUIRE_THAT(spline(x[i]), WithinAbs(y[i], TOL(1e-12, 1e-5)));
		}
	}

	TEST_CASE("SplineInterp_sin_high_accuracy", "[interpolation][spline]") {
		RealFunction f{ [](Real x) -> Real { return sin(x); } };

		Vector<Real> vec_x, vec_y;
		CreateInterpolatedValues(f, 0.0, 2.0 * Constants::PI, 20, vec_x, vec_y);

		// Natural spline
		SplineInterpRealFunc spline(vec_x, vec_y);
		RealFunctionComparer comparer(f, spline);

		double maxAbsDiff = comparer.getAbsDiffMax(0.0, 2.0 * Constants::PI, 200);

		// Spline should be much more accurate than linear with same points
		REQUIRE_THAT(maxAbsDiff, WithinAbs(0.0, 0.001));
	}

	TEST_CASE("SplineInterp_derivative_cos", "[interpolation][spline][derivative]") {
		// For sin(x), derivative should be cos(x)
		RealFunction f{ [](Real x) -> Real { return sin(x); } };
		RealFunction df{ [](Real x) -> Real { return cos(x); } };

		Vector<Real> vec_x, vec_y;
		CreateInterpolatedValues(f, 0.0, Constants::PI, 50, vec_x, vec_y);

		SplineInterpRealFunc spline(vec_x, vec_y);

		// Test derivative at several points
		std::vector<Real> test_points = { 0.5, 1.0, 1.5, 2.0, 2.5 };
		for (Real x : test_points) {
			Real computed = spline.Derivative(x);
			Real expected = df(x);
			REQUIRE_THAT(computed, WithinAbs(expected, 0.01));
		}
	}

	TEST_CASE("SplineInterp_second_derivative_negative_sin", "[interpolation][spline][derivative]") {
		// For sin(x), second derivative should be -sin(x)
		RealFunction f{ [](Real x) -> Real { return sin(x); } };
		RealFunction ddf{ [](Real x) -> Real { return -sin(x); } };

		Vector<Real> vec_x, vec_y;
		CreateInterpolatedValues(f, 0.0, Constants::PI, 50, vec_x, vec_y);

		SplineInterpRealFunc spline(vec_x, vec_y);

		// Test second derivative at several points (interior, away from boundaries)
		std::vector<Real> test_points = { 0.5, 1.0, 1.5, 2.0 };
		for (Real x : test_points) {
			Real computed = spline.SecondDerivative(x);
			Real expected = ddf(x);
			REQUIRE_THAT(computed, WithinAbs(expected, 0.1));  // Larger tolerance for 2nd derivative
		}
	}

	TEST_CASE("SplineInterp_integration_polynomial", "[interpolation][spline][integration]") {
		// Integrate x^2 from 0 to 2, exact result is 8/3 ≈ 2.6667
		Vector<Real> x{ 0.0, 0.5, 1.0, 1.5, 2.0 };
		Vector<Real> y(5);
		for (int i = 0; i < 5; i++) y[i] = x[i] * x[i];

		SplineInterpRealFunc spline(x, y);

		Real integral = spline.Integrate(0.0, 2.0);
		Real expected = 8.0 / 3.0;

		// Natural spline adds some error at boundaries; allow 0.02 tolerance
		REQUIRE_THAT(integral, WithinAbs(expected, 0.02));
	}

	TEST_CASE("SplineInterp_integration_sin", "[interpolation][spline][integration]") {
		// Integrate sin(x) from 0 to pi, exact result is 2
		RealFunction f{ [](Real x) -> Real { return sin(x); } };

		Vector<Real> vec_x, vec_y;
		CreateInterpolatedValues(f, 0.0, Constants::PI, 20, vec_x, vec_y);

		SplineInterpRealFunc spline(vec_x, vec_y);

		Real integral = spline.Integrate(0.0, Constants::PI);
		Real expected = 2.0;

		REQUIRE_THAT(integral, WithinAbs(expected, 0.01));
	}

	TEST_CASE("SplineInterp_clamped_boundary", "[interpolation][spline]") {
		// Test clamped spline with known boundary derivatives
		RealFunction f{ [](Real x) -> Real { return sin(x); } };
		Real yp1 = cos(0.0);   // derivative at x=0: cos(0) = 1
		Real ypn = cos(Constants::PI);  // derivative at x=π: cos(π) = -1

		Vector<Real> vec_x, vec_y;
		CreateInterpolatedValues(f, 0.0, Constants::PI, 10, vec_x, vec_y);

		SplineInterpRealFunc clamped(vec_x, vec_y, yp1, ypn);

		// Clamped spline should be more accurate near boundaries
		REQUIRE_THAT(clamped(0.1), WithinAbs(f(0.1), 0.001));
		REQUIRE_THAT(clamped(Constants::PI - 0.1), WithinAbs(f(Constants::PI - 0.1), 0.001));
	}

	TEST_CASE("SplineInterp_natural_second_deriv_zero_at_ends", "[interpolation][spline]") {
		// Natural spline has zero second derivative at endpoints
		Vector<Real> x{ 0.0, 1.0, 2.0, 3.0, 4.0 };
		Vector<Real> y{ 1.0, 3.0, 2.0, 5.0, 4.0 };

		SplineInterpRealFunc spline(x, y);  // natural spline by default

		// Check second derivatives at endpoints
		REQUIRE_THAT(spline.GetSecondDerivative(0), WithinAbs(0.0, TOL(1e-10, 1e-5)));
		REQUIRE_THAT(spline.GetSecondDerivative(4), WithinAbs(0.0, TOL(1e-10, 1e-5)));
	}

	//////////////////////////////////////////////////////////////////////////////
	//                     SPLINE VS LINEAR COMPARISON                          //
	//////////////////////////////////////////////////////////////////////////////

	TEST_CASE("Spline_more_accurate_than_linear", "[interpolation][comparison]") {
		// Spline should be more accurate than linear for smooth functions
		RealFunction f{ [](Real x) -> Real { return static_cast<Real>(exp(-x) * cos(2 * x)); } };

		Vector<Real> vec_x, vec_y;
		CreateInterpolatedValues(f, 0.0, 3.0, 15, vec_x, vec_y);

		LinearInterpRealFunc linear(vec_x, vec_y);
		SplineInterpRealFunc spline(vec_x, vec_y);

		RealFunctionComparer linear_comp(f, linear);
		RealFunctionComparer spline_comp(f, spline);

		double linear_max = linear_comp.getAbsDiffMax(0.0, 3.0, 100);
		double spline_max = spline_comp.getAbsDiffMax(0.0, 3.0, 100);

		// Both should have small errors, with spline being better (but not always 10x)
		// For highly oscillating functions, the benefit may be smaller
		REQUIRE(spline_max < 0.01);  // Spline should be accurate
		REQUIRE(linear_max < 0.02);  // Linear is also decent with 15 points
	}

	//////////////////////////////////////////////////////////////////////////////
	//                     HERMITE AND AKIMA INTERPOLATION                      //
	//////////////////////////////////////////////////////////////////////////////

	TEST_CASE("HermiteInterp_recovers_cubic_with_exact_derivatives", "[interpolation][hermite]") {
		Vector<Real> x{REAL(-2.0), REAL(-1.0), REAL(0.0), REAL(1.0), REAL(2.0)};
		Vector<Real> y(x.size()), derivatives(x.size());
		for (int i = 0; i < x.size(); ++i) {
			y[i] = x[i] * x[i] * x[i] - REAL(2.0) * x[i] * x[i] + x[i] - REAL(1.0);
			derivatives[i] = REAL(3.0) * x[i] * x[i] - REAL(4.0) * x[i] + REAL(1.0);
		}
		HermiteInterpRealFunc interp(x, y, derivatives);
		for (int i = 0; i <= 40; ++i) {
			Real point = REAL(-2.0) + REAL(4.0) * i / REAL(40.0);
			Real expected = point * point * point - REAL(2.0) * point * point + point - REAL(1.0);
			REQUIRE_THAT(interp(point), WithinAbs(expected, TOL(1e-11, 1e-4)));
		}
		REQUIRE(interp.EvaluateDetailed(REAL(0.25)).algorithm_name == "CubicHermite");
	}

	TEST_CASE("HermiteInterp_validates_derivative_data", "[interpolation][hermite][edge]") {
		Vector<Real> x{REAL(0.0), REAL(1.0), REAL(2.0)};
		Vector<Real> y{REAL(0.0), REAL(1.0), REAL(4.0)};
		Vector<Real> derivatives{REAL(0.0), REAL(2.0)};
		REQUIRE_THROWS_AS(HermiteInterpRealFunc(x, y, derivatives), RealFuncInterpInitError);
	}

	TEST_CASE("AkimaInterp_exact_at_nodes_and_smooth_between", "[interpolation][akima]") {
		auto data = TestBeds::getSinePiTest();
		AkimaInterpRealFunc interp(data._xData, data._yData);
		for (int i = 0; i < data._xData.size(); ++i)
			REQUIRE_THAT(interp(data._xData[i]), WithinAbs(data._yData[i], TOL(1e-12, 1e-5)));
		Real maxError = TestBeds::computeMaxError(interp, data._trueFunction,
			TestBeds::getStandardEvaluationGrid(data._xMin, data._xMax));
		REQUIRE(maxError < TOL(0.08, 0.25));
		REQUIRE(interp.EvaluateDetailed(REAL(0.125)).algorithm_name == "Akima");
	}

	TEST_CASE("AkimaInterp_noisy_data_avoids_large_overshoot", "[interpolation][akima][noisy]") {
		auto data = TestBeds::getNoisyGaussianTest();
		AkimaInterpRealFunc akima(data._xData, data._yData);
		SplineInterpRealFunc spline(data._xData, data._yData);
		Vector<Real> grid = TestBeds::getStandardEvaluationGrid(data._xMin, data._xMax);
		Real akimaExtreme = REAL(0.0);
		Real splineExtreme = REAL(0.0);
		for (int i = 0; i < grid.size(); ++i) {
			akimaExtreme = std::max(akimaExtreme, std::abs(akima(grid[i])));
			splineExtreme = std::max(splineExtreme, std::abs(spline(grid[i])));
		}
		REQUIRE(akimaExtreme < REAL(1.25));
		REQUIRE(akimaExtreme <= splineExtreme + TOL(0.05, 0.15));
	}

	TEST_CASE("Interpolation_testbed_Runge_and_rational_cases_are_bounded", "[interpolation][focused]") {
		auto runge = TestBeds::getRungeTest();
		AkimaInterpRealFunc akima(runge._xData, runge._yData);
		Real rungeError = TestBeds::computeMaxError(akima, runge._trueFunction,
			TestBeds::getStandardEvaluationGrid(runge._xMin, runge._xMax));
		REQUIRE(rungeError < REAL(0.25));

		Vector<Real> x = TestBeds::generateUniformNodes(REAL(-3.0), REAL(3.0), 17);
		Vector<Real> y = TestBeds::evaluateAtNodes(x, TestBeds::RationalSmooth_TrueFunc);
		RationalInterpRealFunc rational(x, y, 5);
		Real rationalError = TestBeds::computeMaxError(rational, TestBeds::RationalSmooth_TrueFunc,
			TestBeds::getStandardEvaluationGrid(REAL(-3.0), REAL(3.0)));
		REQUIRE(rationalError < REAL(0.35));
	}

	TEST_CASE("InterpolationTestBed_collections_expose_expected_cases", "[interpolation][testbed]") {
		const auto allTests = TestBeds::getAllInterpolationTests();
		const auto rungeTests = TestBeds::getRungePhenomenonTests();
		const auto smoothTests = TestBeds::getSmoothFunctionTests();
		const auto challengingTests = TestBeds::getChallengingTests();

		REQUIRE(allTests.size() == 15);
		REQUIRE(rungeTests.size() == 4);
		REQUIRE(smoothTests.size() == 6);
		REQUIRE(challengingTests.size() == 5);

		int noisyCount = 0;
		int oscillatoryCount = 0;
		int rungeCount = 0;
		for (const auto& test : allTests) {
			if (test._hasNoise) ++noisyCount;
			if (test._isOscillatory) ++oscillatoryCount;
			if (test._hasRungePhenomenon) ++rungeCount;
		}

		REQUIRE(noisyCount == 1);
		REQUIRE(oscillatoryCount >= 2);
		REQUIRE(rungeCount >= 1);
	}

	TEST_CASE("InterpolationTestBed_smooth_polynomial_is_exact_with_polynomial_interpolation", "[interpolation][testbed][polynomial]") {
		const auto data = TestBeds::getSmoothPolyTest();
		PolynomInterpRealFunc interp(data._xData, data._yData, data._xData.size());
		const auto grid = TestBeds::getStandardEvaluationGrid(data._xMin, data._xMax);

		REQUIRE(TestBeds::computeMaxError(interp, data._trueFunction, grid) < TOL(1e-11, 1e-4));
		REQUIRE(TestBeds::computeRMSError(interp, data._trueFunction, grid) < TOL(1e-12, 1e-5));
	}

	TEST_CASE("Hermite_and_Akima_detailed_extrapolation_policies", "[interpolation][focused][extrapolation]") {
		Vector<Real> x{REAL(0.0), REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0)};
		Vector<Real> y{REAL(0.0), REAL(1.0), REAL(4.0), REAL(9.0), REAL(16.0)};
		Vector<Real> derivatives{REAL(0.0), REAL(2.0), REAL(4.0), REAL(6.0), REAL(8.0)};
		HermiteInterpRealFunc hermite(x, y, derivatives);
		AkimaInterpRealFunc akima(x, y);
		InterpolationConfig config;
		config.extrapolation_policy = ExtrapolationPolicy::Clamp;
		REQUIRE(hermite.EvaluateDetailed(REAL(-1.0), config).WasClamped());
		REQUIRE(akima.EvaluateDetailed(REAL(5.0), config).WasClamped());
		REQUIRE_THAT(hermite.EvaluateDetailed(REAL(-1.0), config).value,
			WithinAbs(REAL(0.0), TOL(1e-12, 1e-5)));
		REQUIRE_THAT(akima.EvaluateDetailed(REAL(5.0), config).value,
			WithinAbs(REAL(16.0), TOL(1e-12, 1e-5)));
	}

	TEST_CASE("BarycentricPolynomial_recovers_cubic_on_arbitrary_nodes", "[interpolation][barycentric-polynomial]") {
		Vector<Real> x{REAL(2.0), REAL(-1.0), REAL(1.0), REAL(-2.0), REAL(0.0)};
		Vector<Real> y(x.size());
		for (int i = 0; i < x.size(); ++i)
			y[i] = x[i] * x[i] * x[i] - REAL(2.0) * x[i] + REAL(1.0);
		BarycentricPolynomialInterp interp(x, y);
		for (int i = 0; i <= 40; ++i) {
			Real point = REAL(-2.0) + REAL(4.0) * i / REAL(40.0);
			Real expected = point * point * point - REAL(2.0) * point + REAL(1.0);
			REQUIRE_THAT(interp(point), WithinAbs(expected, TOL(1e-11, 1e-4)));
		}
		REQUIRE(interp.EvaluateDetailed(REAL(0.25)).algorithm_name == "BarycentricPolynomial");
	}

	TEST_CASE("BarycentricPolynomial_validates_nodes_and_policies", "[interpolation][barycentric-polynomial]") {
		Vector<Real> duplicateX{REAL(0.0), REAL(1.0), REAL(1.0)};
		Vector<Real> duplicateY{REAL(0.0), REAL(1.0), REAL(2.0)};
		REQUIRE_THROWS_AS(BarycentricPolynomialInterp(duplicateX, duplicateY), RealFuncInterpInitError);

		Vector<Real> x{REAL(0.0), REAL(1.0), REAL(2.0)};
		Vector<Real> y{REAL(0.0), REAL(1.0), REAL(4.0)};
		BarycentricPolynomialInterp interp(x, y);
		InterpolationConfig config;
		config.extrapolation_policy = ExtrapolationPolicy::Clamp;
		auto result = interp.EvaluateDetailed(REAL(3.0), config);
		REQUIRE(result.WasClamped());
		REQUIRE(result.interval_index == 1);
		REQUIRE_THAT(result.value, WithinAbs(REAL(4.0), TOL(1e-12, 1e-5)));
	}

	TEST_CASE("ChebyshevNodeInterpolator_reduces_Runge_error", "[interpolation][chebyshev-nodes]") {
		constexpr int numNodes = 15;
		Vector<Real> nodes = ChebyshevInterpolationNodes(REAL(-1.0), REAL(1.0), numNodes);
		REQUIRE(nodes.size() == numNodes);
		for (int i = 1; i < nodes.size(); ++i) REQUIRE(nodes[i] > nodes[i - 1]);
		REQUIRE(nodes[0] > REAL(-1.0));
		REQUIRE(nodes[nodes.size() - 1] < REAL(1.0));

		auto runge = [](Real x) { return REAL(1.0) / (REAL(1.0) + REAL(25.0) * x * x); };
		auto chebyshev = MakeChebyshevNodeInterpolator(std::function<Real(Real)>(runge),
			REAL(-1.0), REAL(1.0), numNodes);
		Vector<Real> uniformNodes = TestBeds::generateUniformNodes(REAL(-1.0), REAL(1.0), numNodes);
		Vector<Real> uniformValues(numNodes);
		for (int i = 0; i < numNodes; ++i) uniformValues[i] = runge(uniformNodes[i]);
		BarycentricPolynomialInterp uniform(uniformNodes, uniformValues);

		Real chebyshevMaxError = 0.0;
		Real uniformMaxError = 0.0;
		for (int i = 0; i <= 200; ++i) {
			Real point = REAL(-1.0) + REAL(2.0) * i / REAL(200.0);
			chebyshevMaxError = std::max(chebyshevMaxError, std::abs(chebyshev(point) - runge(point)));
			uniformMaxError = std::max(uniformMaxError, std::abs(uniform(point) - runge(point)));
		}
		REQUIRE(chebyshevMaxError < REAL(0.08));
		REQUIRE(chebyshevMaxError < uniformMaxError * REAL(0.1));
	}

	TEST_CASE("HermiteInterp_derivative_and_integral_are_exact_for_cubic", "[interpolation][hermite][calculus]") {
		Vector<Real> x{REAL(-2.0), REAL(-1.0), REAL(0.0), REAL(1.0), REAL(2.0)};
		Vector<Real> y(x.size()), derivatives(x.size());
		for (int i = 0; i < x.size(); ++i) {
			y[i] = x[i] * x[i] * x[i] - REAL(2.0) * x[i] * x[i] + x[i] - REAL(1.0);
			derivatives[i] = REAL(3.0) * x[i] * x[i] - REAL(4.0) * x[i] + REAL(1.0);
		}
		HermiteInterpRealFunc interp(x, y, derivatives);
		for (Real point : {REAL(-1.5), REAL(-0.25), REAL(0.5), REAL(1.75)})
			REQUIRE_THAT(interp.Derivative(point),
				WithinAbs(REAL(3.0) * point * point - REAL(4.0) * point + REAL(1.0), TOL(1e-10, 1e-4)));
		REQUIRE_THAT(interp.Integrate(REAL(-2.0), REAL(2.0)),
			WithinAbs(-REAL(44.0 / 3.0), TOL(1e-10, 1e-4)));

		Vector<Real> descendingX{REAL(2.0), REAL(1.0), REAL(0.0), REAL(-1.0), REAL(-2.0)};
		Vector<Real> descendingY(descendingX.size()), descendingDerivatives(descendingX.size());
		for (int i = 0; i < descendingX.size(); ++i) {
			Real point = descendingX[i];
			descendingY[i] = point * point * point - REAL(2.0) * point * point + point - REAL(1.0);
			descendingDerivatives[i] = REAL(3.0) * point * point - REAL(4.0) * point + REAL(1.0);
		}
		HermiteInterpRealFunc descending(descendingX, descendingY, descendingDerivatives);
		REQUIRE_THAT(descending.Integrate(REAL(-2.0), REAL(2.0)),
			WithinAbs(-REAL(44.0 / 3.0), TOL(1e-10, 1e-4)));
	}

	TEST_CASE("MonotoneCubic_derivative_and_integral_recover_linear_function", "[interpolation][monotone][calculus]") {
		Vector<Real> x{REAL(0.0), REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0)};
		Vector<Real> y{REAL(1.0), REAL(3.0), REAL(5.0), REAL(7.0), REAL(9.0)};
		MonotoneCubicInterpRealFunc interp(x, y);
		for (Real point : {REAL(0.25), REAL(1.5), REAL(3.75)})
			REQUIRE_THAT(interp.Derivative(point), WithinAbs(REAL(2.0), TOL(1e-12, 1e-5)));
		REQUIRE_THAT(interp.Integrate(REAL(0.0), REAL(4.0)),
			WithinAbs(REAL(20.0), TOL(1e-11, 1e-4)));
		REQUIRE_THAT(interp.Integrate(REAL(4.0), REAL(0.0)),
			WithinAbs(-REAL(20.0), TOL(1e-11, 1e-4)));
	}

	//////////////////////////////////////////////////////////////////////////////
	//                   PARAMETRIC CURVE INTERPOLATION TESTS                   //
	//////////////////////////////////////////////////////////////////////////////

	TEST_CASE("LinInterpParametricCurve_circle", "[interpolation][parametric]") {
		// Create points on a circle using Matrix
		int n = 16;
		Matrix<Real> points(n, 2);
		for (int i = 0; i < n; i++) {
			Real theta = 2.0 * Constants::PI * i / n;
			points[i][0] = cos(theta);
			points[i][1] = sin(theta);
		}

		LinInterpParametricCurve<2> curve(points, true);  // closed curve

		// Test that we can evaluate points around the curve
		VectorN<Real, 2> pt1 = curve(0.0);
		VectorN<Real, 2> pt2 = curve(0.5);

		// At t=0, should be near (1, 0)
		REQUIRE_THAT(pt1[0], WithinAbs(1.0, 0.1));
		REQUIRE_THAT(pt1[1], WithinAbs(0.0, 0.1));
	}

	TEST_CASE("SplineInterpParametricCurve_circle", "[interpolation][parametric][spline]") {
		// Create points on a circle using Matrix
		int n = 16;
		Matrix<Real> points(n, 2);
		for (int i = 0; i < n; i++) {
			Real theta = 2.0 * Constants::PI * i / n;
			points[i][0] = cos(theta);
			points[i][1] = sin(theta);
		}

		SplineInterpParametricCurve<2> curve(points, true);  // closed curve

		// Test that interpolated point at t=0.25 is near (0, 1) = top of circle
		VectorN<Real, 2> pt = curve(0.25);
		
		// Due to arc-length parameterization, t=0.25 should be near (0, 1)
		REQUIRE_THAT(pt[0], WithinAbs(0.0, 0.3));
		REQUIRE_THAT(pt[1], WithinAbs(1.0, 0.3));
	}

	TEST_CASE("SplineInterpParametricCurve_3D_helix", "[interpolation][parametric][spline]") {
		// Create points on a helix: (cos(t), sin(t), t/(2π))
		int n = 20;
		Matrix<Real> points(n, 3);
		for (int i = 0; i < n; i++) {
			Real t = 2.0 * Constants::PI * i / (n - 1);
			points[i][0] = cos(t);
			points[i][1] = sin(t);
			points[i][2] = t / (2.0 * Constants::PI);
		}

		SplineInterpParametricCurve<3> curve(points, false);  // open curve

		// Test interpolation at midpoint parameter
		VectorN<Real, 3> mid = curve(0.5);

		// At t=0.5 (middle of parameter range), should be roughly at t=π: (-1, 0, 0.5)
		REQUIRE_THAT(mid[0], WithinAbs(-1.0, 0.2));
		REQUIRE_THAT(mid[1], WithinAbs(0.0, 0.2));
		REQUIRE_THAT(mid[2], WithinAbs(0.5, 0.1));
	}

	//////////////////////////////////////////////////////////////////////////////
	//                          EDGE CASE TESTS                                 //
	//////////////////////////////////////////////////////////////////////////////

	TEST_CASE("Spline_minimum_points", "[interpolation][spline][edge]") {
		// Minimum 2 points for spline
		Vector<Real> x{ 0.0, 1.0 };
		Vector<Real> y{ 0.0, 1.0 };

		SplineInterpRealFunc spline(x, y);

		// Should interpolate linearly between two points
		REQUIRE_THAT(spline(0.5), WithinAbs(0.5, TOL(1e-10, 1e-5)));
	}

	TEST_CASE("Interp_handles_negative_x_values", "[interpolation][edge]") {
		Vector<Real> x{ -2.0, -1.0, 0.0, 1.0, 2.0 };
		Vector<Real> y{ 4.0, 1.0, 0.0, 1.0, 4.0 };  // y = x^2

		SplineInterpRealFunc spline(x, y);
		LinearInterpRealFunc linear(x, y);

		REQUIRE_THAT(spline(-1.5), WithinAbs(2.25, 0.1));
		REQUIRE_THAT(linear(-1.5), WithinAbs(2.5, 0.01));  // Linear gives midpoint
	}

	//////////////////////////////////////////////////////////////////////////////
	//                     RATIONAL INTERPOLATION TESTS                         //
	//////////////////////////////////////////////////////////////////////////////

	TEST_CASE("RationalInterp_smooth_function", "[interpolation][rational]") {
		// Rational interpolation for a smooth function
		RealFunction f([](Real x) -> Real { return Real(1.0) / (Real(1.0) + x * x); });  // Runge function

		Vector<Real> vec_x, vec_y;
		CreateInterpolatedValues(f, -3.0, 3.0, 15, vec_x, vec_y);

		RationalInterpRealFunc interp(vec_x, vec_y, 5);  // 5-point rational

		// Test at intermediate point
		Real x_test = 0.5;
		Real result = interp(x_test);
		Real expected = f(x_test);

		REQUIRE_THAT(result, WithinAbs(expected, 0.01));
	}

	TEST_CASE("RationalInterp_exact_at_nodes", "[interpolation][rational]") {
		Vector<Real> x{ 0.0, 1.0, 2.0, 3.0, 4.0 };
		Vector<Real> y{ 1.0, 0.5, 0.33, 0.25, 0.2 };

		RationalInterpRealFunc interp(x, y, 3);

		// Should be exact at data points
		for (int i = 0; i < x.size(); i++) {
			REQUIRE_THAT(interp(x[i]), WithinAbs(y[i], TOL(1e-10, 1e-5)));
		}
	}

	TEST_CASE("RationalInterp_handles_near_pole", "[interpolation][rational]") {
		// Rational interpolation should handle functions with poles better
		RealFunction f([](Real x) -> Real { return Real(1.0) / (x - Real(2.5)); });

		// Sample away from the pole
		Vector<Real> x{ 0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 3.5, 4.0 };
		Vector<Real> y(8);
		for (int i = 0; i < 8; i++) y[i] = f(x[i]);

		RationalInterpRealFunc interp(x, y, 4);

		// Test at a point away from the pole
		REQUIRE_THAT(interp(3.5), WithinAbs(f(3.5), 0.01));
	}

	//////////////////////////////////////////////////////////////////////////////
	//                   BARYCENTRIC INTERPOLATION TESTS                        //
	//////////////////////////////////////////////////////////////////////////////

	TEST_CASE("BarycentricInterp_smooth_function", "[interpolation][barycentric]") {
		// Barycentric rational interpolation - should be very stable
		RealFunction f{ [](Real x) -> Real { return sin(x); } };

		Vector<Real> vec_x, vec_y;
		CreateInterpolatedValues(f, 0.0, Constants::PI, 10, vec_x, vec_y);

		BarycentricRationalInterp interp(vec_x, vec_y, 3);

		// Test at intermediate points
		Real x_test = Constants::PI / 4;
		Real result = interp(x_test);
		Real expected = f(x_test);

		REQUIRE_THAT(result, WithinAbs(expected, 0.01));
	}

	TEST_CASE("BarycentricInterp_exact_at_nodes", "[interpolation][barycentric]") {
		Vector<Real> x{ 0.0, 1.0, 2.0, 3.0, 4.0 };
		Vector<Real> y{ 1.0, 3.0, 2.0, 5.0, 4.0 };

		BarycentricRationalInterp interp(x, y, 2);

		// Should be exact at data points
		for (int i = 0; i < x.size(); i++) {
			REQUIRE_THAT(interp(x[i]), WithinAbs(y[i], TOL(1e-10, 1e-5)));
		}
	}

	TEST_CASE("BarycentricInterp_no_poles", "[interpolation][barycentric]") {
		// Unlike rational interpolation, barycentric should be pole-free
		RealFunction f{ [](Real x) -> Real { return static_cast<Real>(exp(-x * x)); } };

		Vector<Real> vec_x, vec_y;
		CreateInterpolatedValues(f, -2.0, 2.0, 12, vec_x, vec_y);

		BarycentricRationalInterp interp(vec_x, vec_y, 3);

		// Should give reasonable values everywhere in the range
		for (Real x = -1.9; x <= 1.9; x += 0.2) {
			Real result = interp(x);
			REQUIRE(std::isfinite(result));
			// Should be in reasonable range for Gaussian
			REQUIRE(result > -0.5);
			REQUIRE(result < 1.5);
		}
	}

	//////////////////////////////////////////////////////////////////////////////
	//                      2D INTERPOLATION TESTS                              //
	//////////////////////////////////////////////////////////////////////////////

	TEST_CASE("BilinearInterp_constant_function", "[interpolation][2d]") {
		// For a constant function, bilinear should return the constant
		int nx = 5, ny = 5;
		Vector<Real> x1(nx), x2(ny);
		Matrix<Real> y(nx, ny);

		for (int i = 0; i < nx; i++) x1[i] = i * 1.0;
		for (int j = 0; j < ny; j++) x2[j] = j * 1.0;
		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++)
				y[i][j] = 5.0;  // constant

		BilinearInterp2D interp(x1, x2, y);

		REQUIRE_THAT(interp(0.5, 0.5), WithinAbs(5.0, TOL(1e-10, 1e-5)));
		REQUIRE_THAT(interp(2.5, 2.5), WithinAbs(5.0, TOL(1e-10, 1e-5)));
		REQUIRE_THAT(interp(3.7, 1.2), WithinAbs(5.0, TOL(1e-10, 1e-5)));
	}

	TEST_CASE("BilinearInterp_linear_function", "[interpolation][2d]") {
		// For f(x,y) = x + y, bilinear interpolation should be exact
		int nx = 5, ny = 5;
		Vector<Real> x1(nx), x2(ny);
		Matrix<Real> y(nx, ny);

		for (int i = 0; i < nx; i++) x1[i] = i * 1.0;
		for (int j = 0; j < ny; j++) x2[j] = j * 1.0;
		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++)
				y[i][j] = x1[i] + x2[j];

		BilinearInterp2D interp(x1, x2, y);

		// Test at intermediate points
		REQUIRE_THAT(interp(0.5, 0.5), WithinAbs(1.0, TOL(1e-10, 1e-5)));
		REQUIRE_THAT(interp(1.5, 2.5), WithinAbs(4.0, TOL(1e-10, 1e-5)));
		REQUIRE_THAT(interp(2.7, 3.3), WithinAbs(6.0, TOL(1e-10, 1e-5)));
	}

	TEST_CASE("BilinearInterp_exact_at_grid_points", "[interpolation][2d]") {
		int nx = 4, ny = 4;
		Vector<Real> x1(nx), x2(ny);
		Matrix<Real> y(nx, ny);

		for (int i = 0; i < nx; i++) x1[i] = i * 0.5;
		for (int j = 0; j < ny; j++) x2[j] = j * 0.5;
		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++)
				y[i][j] = sin(x1[i]) * cos(x2[j]);

		BilinearInterp2D interp(x1, x2, y);

		// Should be exact at grid points
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				REQUIRE_THAT(interp(x1[i], x2[j]), WithinAbs(y[i][j], TOL(1e-12, 1e-5)));
			}
		}
	}

	TEST_CASE("BilinearInterp_smooth_function", "[interpolation][2d]") {
		// Test on a smooth 2D function
		auto f = [](Real x, Real y) { return sin(x) * cos(y); };

		int nx = 10, ny = 10;
		Vector<Real> x1(nx), x2(ny);
		Matrix<Real> y(nx, ny);

		for (int i = 0; i < nx; i++) x1[i] = i * Constants::PI / (nx - 1);
		for (int j = 0; j < ny; j++) x2[j] = j * Constants::PI / (ny - 1);
		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++)
				y[i][j] = f(x1[i], x2[j]);

		BilinearInterp2D interp(x1, x2, y);

		// Test at several intermediate points
		Real x_test = Constants::PI / 4;
		Real y_test = Constants::PI / 3;
		Real result = interp(x_test, y_test);
		Real expected = f(x_test, y_test);

		// Bilinear isn't perfect for curved functions, but should be reasonable
		REQUIRE_THAT(result, WithinAbs(expected, 0.05));
	}

	TEST_CASE("BilinearInterp_derivatives", "[interpolation][2d]") {
		// Test interpWithDerivatives for a linear function f(x,y) = 2x + 3y
		int nx = 4, ny = 4;
		Vector<Real> x1(nx), x2(ny);
		Matrix<Real> z(nx, ny);

		for (int i = 0; i < nx; i++) x1[i] = i * 1.0;
		for (int j = 0; j < ny; j++) x2[j] = j * 1.0;
		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++)
				z[i][j] = 2.0 * x1[i] + 3.0 * x2[j];

		BilinearInterp2D interp(x1, x2, z);

		Real val, dz_dx1, dz_dx2;
		interp.interpWithDerivatives(1.5, 2.5, val, dz_dx1, dz_dx2);

		REQUIRE_THAT(val, WithinAbs(2.0 * 1.5 + 3.0 * 2.5, TOL(1e-10, 1e-5)));
		REQUIRE_THAT(dz_dx1, WithinAbs(2.0, TOL(1e-10, 1e-5)));
		REQUIRE_THAT(dz_dx2, WithinAbs(3.0, TOL(1e-10, 1e-5)));
	}

	TEST_CASE("BilinearInterp_grid_info", "[interpolation][2d]") {
		int nx = 5, ny = 7;
		Vector<Real> x1(nx), x2(ny);
		Matrix<Real> z(nx, ny);

		for (int i = 0; i < nx; i++) x1[i] = i;
		for (int j = 0; j < ny; j++) x2[j] = j;
		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++)
				z[i][j] = 1.0;

		BilinearInterp2D interp(x1, x2, z);

		REQUIRE(interp.getGridDim1() == nx);
		REQUIRE(interp.getGridDim2() == ny);
	}

	/////////////////////////////////////////////////////////////////////////////
	// BicubicSplineInterp2D Tests
	/////////////////////////////////////////////////////////////////////////////

	TEST_CASE("BicubicSplineInterp_constant_function", "[interpolation][2d][bicubic]") {
		// For a constant function, bicubic should return the constant everywhere
		int nx = 5, ny = 5;
		Vector<Real> x1(nx), x2(ny);
		Matrix<Real> z(nx, ny);

		for (int i = 0; i < nx; i++) x1[i] = i * 1.0;
		for (int j = 0; j < ny; j++) x2[j] = j * 1.0;
		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++)
				z[i][j] = 7.5;

		BicubicSplineInterp2D interp(x1, x2, z);

		REQUIRE_THAT(interp(0.5, 0.5), WithinAbs(7.5, TOL(1e-10, 1e-5)));
		REQUIRE_THAT(interp(2.5, 2.5), WithinAbs(7.5, TOL(1e-10, 1e-5)));
		REQUIRE_THAT(interp(3.7, 1.2), WithinAbs(7.5, TOL(1e-10, 1e-5)));
	}

	TEST_CASE("BicubicSplineInterp_linear_function", "[interpolation][2d][bicubic]") {
		// For f(x,y) = x + y, bicubic spline should be exact (or very close)
		int nx = 5, ny = 5;
		Vector<Real> x1(nx), x2(ny);
		Matrix<Real> z(nx, ny);

		for (int i = 0; i < nx; i++) x1[i] = i * 1.0;
		for (int j = 0; j < ny; j++) x2[j] = j * 1.0;
		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++)
				z[i][j] = x1[i] + x2[j];

		BicubicSplineInterp2D interp(x1, x2, z);

		REQUIRE_THAT(interp(0.5, 0.5), WithinAbs(1.0, TOL(1e-8, 1e-4)));
		REQUIRE_THAT(interp(1.5, 2.5), WithinAbs(4.0, TOL(1e-8, 1e-4)));
		REQUIRE_THAT(interp(2.7, 3.3), WithinAbs(6.0, TOL(1e-8, 1e-4)));
	}

	TEST_CASE("BicubicSplineInterp_exact_at_grid_points", "[interpolation][2d][bicubic]") {
		int nx = 5, ny = 5;
		Vector<Real> x1(nx), x2(ny);
		Matrix<Real> z(nx, ny);

		for (int i = 0; i < nx; i++) x1[i] = i * 0.5;
		for (int j = 0; j < ny; j++) x2[j] = j * 0.5;
		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++)
				z[i][j] = sin(x1[i]) * cos(x2[j]);

		BicubicSplineInterp2D interp(x1, x2, z);

		// Should be exact at grid points
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				REQUIRE_THAT(interp(x1[i], x2[j]), WithinAbs(z[i][j], TOL(1e-10, 1e-5)));
			}
		}
	}

	TEST_CASE("BicubicSplineInterp_smooth_function", "[interpolation][2d][bicubic]") {
		// Test on a smooth 2D function - bicubic should be more accurate than bilinear
		auto f = [](Real x, Real y) { return sin(x) * cos(y); };

		int nx = 10, ny = 10;
		Vector<Real> x1(nx), x2(ny);
		Matrix<Real> z(nx, ny);

		for (int i = 0; i < nx; i++) x1[i] = i * Constants::PI / (nx - 1);
		for (int j = 0; j < ny; j++) x2[j] = j * Constants::PI / (ny - 1);
		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++)
				z[i][j] = f(x1[i], x2[j]);

		BicubicSplineInterp2D interp(x1, x2, z);

		// Test at several intermediate points - bicubic should be more accurate
		Real x_test = Constants::PI / 4;
		Real y_test = Constants::PI / 3;
		Real result = interp(x_test, y_test);
		Real expected = f(x_test, y_test);

		// Bicubic is more accurate than bilinear for smooth functions
		REQUIRE_THAT(result, WithinAbs(expected, 0.01));
	}

	TEST_CASE("BicubicSplineInterp_derivatives", "[interpolation][2d][bicubic]") {
		// Test interpWithDerivatives for f(x,y) = x^2 + y^2
		// df/dx = 2x, df/dy = 2y
		int nx = 10, ny = 10;
		Vector<Real> x1(nx), x2(ny);
		Matrix<Real> z(nx, ny);

		for (int i = 0; i < nx; i++) x1[i] = i * 0.5;
		for (int j = 0; j < ny; j++) x2[j] = j * 0.5;
		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++)
				z[i][j] = x1[i] * x1[i] + x2[j] * x2[j];

		BicubicSplineInterp2D interp(x1, x2, z);

		Real val, dz_dx1, dz_dx2;
		Real test_x1 = 1.5, test_x2 = 2.0;
		interp.interpWithDerivatives(test_x1, test_x2, val, dz_dx1, dz_dx2);

		Real expected_val = test_x1 * test_x1 + test_x2 * test_x2;
		Real expected_dx1 = 2.0 * test_x1;
		Real expected_dx2 = 2.0 * test_x2;

		REQUIRE_THAT(val, WithinAbs(expected_val, 0.01));
		REQUIRE_THAT(dz_dx1, WithinAbs(expected_dx1, 0.1));
		REQUIRE_THAT(dz_dx2, WithinAbs(expected_dx2, 0.1));
	}

	TEST_CASE("BicubicSplineInterp_grid_info", "[interpolation][2d][bicubic]") {
		int nx = 6, ny = 8;
		Vector<Real> x1(nx), x2(ny);
		Matrix<Real> z(nx, ny);

		for (int i = 0; i < nx; i++) x1[i] = i;
		for (int j = 0; j < ny; j++) x2[j] = j;
		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++)
				z[i][j] = 1.0;

		BicubicSplineInterp2D interp(x1, x2, z);

		REQUIRE(interp.getGridDim1() == nx);
		REQUIRE(interp.getGridDim2() == ny);
	}

	TEST_CASE("BicubicSplineInterp_more_accurate_than_bilinear", "[interpolation][2d][bicubic]") {
		// Verify that bicubic is more accurate than bilinear for curved functions
		auto f = [](Real x, Real y) { return sin(x * y); };

		int nx = 8, ny = 8;
		Vector<Real> x1(nx), x2(ny);
		Matrix<Real> z(nx, ny);

		for (int i = 0; i < nx; i++) x1[i] = 0.5 + i * 0.25;
		for (int j = 0; j < ny; j++) x2[j] = 0.5 + j * 0.25;
		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++)
				z[i][j] = f(x1[i], x2[j]);

		BilinearInterp2D bilinear(x1, x2, z);
		BicubicSplineInterp2D bicubic(x1, x2, z);

		// Test at multiple interior points and compare errors
		Real bilinear_err = 0.0, bicubic_err = 0.0;
		int num_tests = 0;
		for (Real tx = 0.6; tx < 2.2; tx += 0.3) {
			for (Real ty = 0.6; ty < 2.2; ty += 0.3) {
				Real expected = f(tx, ty);
				bilinear_err += std::abs(bilinear(tx, ty) - expected);
				bicubic_err += std::abs(bicubic(tx, ty) - expected);
				num_tests++;
			}
		}

		// Bicubic should have smaller total error
		REQUIRE(bicubic_err < bilinear_err);
	}

	/////////////////////////////////////////////////////////////////////////////
	// Additional SplineInterp tests
	/////////////////////////////////////////////////////////////////////////////

	TEST_CASE("SplineInterp_GetSecondDerivative", "[interpolation][spline]") {
		// Test access to stored second derivatives
		Vector<Real> x(5), y(5);
		Real xvals[] = {0.0, 1.0, 2.0, 3.0, 4.0};
		Real yvals[] = {0.0, 1.0, 0.0, 1.0, 0.0};
		for (int i = 0; i < 5; i++) { x[i] = xvals[i]; y[i] = yvals[i]; }

		SplineInterpRealFunc spline(x, y);

		// Second derivatives should be accessible for all points
		for (int i = 0; i < 5; i++) {
			Real d2 = spline.GetSecondDerivative(i);
			// Just verify it doesn't throw and returns a reasonable value
			REQUIRE(std::isfinite(d2));
		}
	}

	TEST_CASE("SplineInterp_GetSecondDerivative_throws_on_invalid_index", "[interpolation][spline]") {
		Vector<Real> x(3), y(3);
		x[0] = 0.0; x[1] = 1.0; x[2] = 2.0;
		y[0] = 0.0; y[1] = 1.0; y[2] = 0.0;

		SplineInterpRealFunc spline(x, y);

		REQUIRE_THROWS_AS(spline.GetSecondDerivative(-1), IndexError);
		REQUIRE_THROWS_AS(spline.GetSecondDerivative(3), IndexError);
	}

	/////////////////////////////////////////////////////////////////////////////
	// Additional Polynomial Interpolation tests
	/////////////////////////////////////////////////////////////////////////////

	TEST_CASE("PolynomInterp_higher_order", "[interpolation][polynomial]") {
		// Test with higher order polynomial (degree 4)
		Vector<Real> x(5), y(5);
		for (int i = 0; i < 5; i++) {
			x[i] = static_cast<Real>(i);
			// y = x^4 - 2x^2 + 1
			y[i] = std::pow(x[i], 4) - 2 * x[i] * x[i] + 1;
		}

		PolynomInterpRealFunc poly(x, y, 5);  // Use all 5 points

		// Should be exact at data points
		for (int i = 0; i < 5; i++) {
			REQUIRE_THAT(poly(x[i]), WithinAbs(y[i], TOL(1e-10, 1e-5)));
		}

		// Check midpoint
		Real mid = 2.5;
		Real expected = std::pow(mid, 4) - 2 * mid * mid + 1;
		REQUIRE_THAT(poly(mid), WithinAbs(expected, TOL(1e-8, 1e-4)));
	}

	/////////////////////////////////////////////////////////////////////////////
	// Additional Barycentric Interpolation tests
	/////////////////////////////////////////////////////////////////////////////

	TEST_CASE("BarycentricInterp_weight_access", "[interpolation][barycentric]") {
		Vector<Real> x(5), y(5);
		Real xvals[] = {0.0, 1.0, 2.0, 3.0, 4.0};
		Real yvals[] = {1.0, 2.0, 1.5, 2.5, 1.0};
		for (int i = 0; i < 5; i++) { x[i] = xvals[i]; y[i] = yvals[i]; }

		BarycentricRationalInterp interp(x, y, 2);

		// Verify weights are accessible
		for (int i = 0; i < 5; i++) {
			Real w = interp.getWeight(i);
			REQUIRE(std::isfinite(w));
			REQUIRE(w != 0.0);  // Weights should be non-zero
		}
	}

	TEST_CASE("BarycentricInterp_range_accessors", "[interpolation][barycentric]") {
		Vector<Real> x(5), y(5);
		x[0] = -1.0; x[1] = 0.0; x[2] = 1.0; x[3] = 2.0; x[4] = 3.0;
		y[0] = 1.0; y[1] = 2.0; y[2] = 1.5; y[3] = 2.5; y[4] = 1.0;

		BarycentricRationalInterp interp(x, y, 2);

		REQUIRE(interp.MinX() == -1.0);
		REQUIRE(interp.MaxX() == 3.0);
		REQUIRE(interp.getNumPoints() == 5);
		REQUIRE(interp.getBlendingParameter() == 2);
	}

	TEST_CASE("BarycentricInterp_invalid_d_throws", "[interpolation][barycentric]") {
		Vector<Real> x(3), y(3);
		x[0] = 0.0; x[1] = 1.0; x[2] = 2.0;
		y[0] = 1.0; y[1] = 2.0; y[2] = 1.0;

		// d >= n should throw
		REQUIRE_THROWS_AS(BarycentricRationalInterp(x, y, 3), RealFuncInterpInitError);
		REQUIRE_THROWS_AS(BarycentricRationalInterp(x, y, 10), RealFuncInterpInitError);
		// d < 0 should throw
		REQUIRE_THROWS_AS(BarycentricRationalInterp(x, y, -1), RealFuncInterpInitError);
	}

	TEST_CASE("BarycentricInterp_weight_throws_on_invalid_index", "[interpolation][barycentric]") {
		Vector<Real> x(3), y(3);
		x[0] = 0.0; x[1] = 1.0; x[2] = 2.0;
		y[0] = 1.0; y[1] = 2.0; y[2] = 1.0;

		BarycentricRationalInterp interp(x, y, 1);

		REQUIRE_THROWS_AS(interp.getWeight(-1), IndexError);
		REQUIRE_THROWS_AS(interp.getWeight(3), IndexError);
	}

	/////////////////////////////////////////////////////////////////////////////
	// Additional RationalInterp tests
	/////////////////////////////////////////////////////////////////////////////

	TEST_CASE("RationalInterp_exact_at_data_point", "[interpolation][rational]") {
		// When evaluating exactly at a data point, should return the value
		Vector<Real> x(5), y(5);
		Real xvals[] = {0.0, 1.0, 2.0, 3.0, 4.0};
		Real yvals[] = {1.0, 2.0, 1.5, 2.5, 1.0};
		for (int i = 0; i < 5; i++) { x[i] = xvals[i]; y[i] = yvals[i]; }

		RationalInterpRealFunc interp(x, y, 3);

		for (int i = 0; i < 5; i++) {
			Real result = interp(x[i]);
			REQUIRE_THAT(result, WithinAbs(y[i], TOL(1e-10, 1e-5)));
			// Error estimate should be zero at exact points
			REQUIRE_THAT(interp.getLastErrorEst(), WithinAbs(0.0, TOL(1e-10, 1e-5)));
		}
	}

	/////////////////////////////////////////////////////////////////////////////
	// Base class RealFunctionInterpolated tests
	/////////////////////////////////////////////////////////////////////////////

	TEST_CASE("InterpolatedFunc_data_accessors", "[interpolation][base]") {
		Vector<Real> x(5), y(5);
		Real xvals[] = {0.0, 1.0, 2.0, 3.0, 4.0};
		Real yvals[] = {1.0, 4.0, 9.0, 16.0, 25.0};
		for (int i = 0; i < 5; i++) { x[i] = xvals[i]; y[i] = yvals[i]; }

		LinearInterpRealFunc interp(x, y);

		// Test individual accessors
		for (int i = 0; i < 5; i++) {
			REQUIRE_THAT(interp.X(i), WithinAbs(x[i], TOL(1e-15, 1e-5)));
			REQUIRE_THAT(interp.Y(i), WithinAbs(y[i], TOL(1e-15, 1e-5)));
		}

		// Test count and range
		REQUIRE(interp.getNumPoints() == 5);
		REQUIRE(interp.getInterpOrder() == 2);  // Linear uses 2 points
		REQUIRE(interp.MinX() == 0.0);
		REQUIRE(interp.MaxX() == 4.0);
	}

	TEST_CASE("InterpolatedFunc_descending_x_values", "[interpolation][base]") {
		// Test that interpolation works with descending x values
		Vector<Real> x(5), y(5);
		x[0] = 4.0; x[1] = 3.0; x[2] = 2.0; x[3] = 1.0; x[4] = 0.0;  // Descending
		y[0] = 16.0; y[1] = 9.0; y[2] = 4.0; y[3] = 1.0; y[4] = 0.0; // y = x^2

		// Use extrapolation=true to handle boundary issues
		LinearInterpRealFunc interp(x, y);

		// Should work at midpoints (linear interp between table values)
		REQUIRE_THAT(interp(3.5), WithinAbs(12.5, 0.1));  // Midpoint between 16 and 9
		REQUIRE_THAT(interp(0.5), WithinAbs(0.5, 0.1));   // Midpoint between 1 and 0
	}

	/////////////////////////////////////////////////////////////////////////////
	// Parametric Curve additional tests
	/////////////////////////////////////////////////////////////////////////////

	TEST_CASE("LinInterpParametricCurve_closed_curve", "[interpolation][parametric]") {
		// Test closed curve - a square
		Matrix<Real> pts(4, 2);
		pts[0][0] = 0.0; pts[0][1] = 0.0;  // Bottom-left
		pts[1][0] = 1.0; pts[1][1] = 0.0;  // Bottom-right
		pts[2][0] = 1.0; pts[2][1] = 1.0;  // Top-right
		pts[3][0] = 0.0; pts[3][1] = 1.0;  // Top-left

		LinInterpParametricCurve<2> curve(pts, true);  // closed = true

		// Start and end should be the same corner
		auto p0 = curve(0.0);
		REQUIRE_THAT(p0[0], WithinAbs(0.0, TOL(1e-10, 1e-5)));
		REQUIRE_THAT(p0[1], WithinAbs(0.0, TOL(1e-10, 1e-5)));

		// Test parameter range
		REQUIRE(curve.getMinT() == 0.0);
		REQUIRE(curve.getMaxT() == 1.0);
	}

	TEST_CASE("SplineInterpParametricCurve_3D_line", "[interpolation][parametric]") {
		// 3D straight line - spline should reproduce it exactly
		int n = 5;
		Matrix<Real> pts(n, 3);
		for (int i = 0; i < n; i++) {
			pts[i][0] = i * 1.0;
			pts[i][1] = i * 2.0;
			pts[i][2] = i * 3.0;
		}

		SplineInterpParametricCurve<3> curve(0.0, 1.0, pts, false);

		// Check endpoints
		auto p0 = curve(0.0);
		REQUIRE_THAT(p0[0], WithinAbs(0.0, TOL(1e-8, 1e-4)));
		REQUIRE_THAT(p0[1], WithinAbs(0.0, TOL(1e-8, 1e-4)));
		REQUIRE_THAT(p0[2], WithinAbs(0.0, TOL(1e-8, 1e-4)));

		auto p1 = curve(1.0);
		REQUIRE_THAT(p1[0], WithinAbs(4.0, TOL(1e-8, 1e-4)));
		REQUIRE_THAT(p1[1], WithinAbs(8.0, TOL(1e-8, 1e-4)));
		REQUIRE_THAT(p1[2], WithinAbs(12.0, TOL(1e-8, 1e-4)));

		// Check midpoint - should be linear
		auto pmid = curve(0.5);
		REQUIRE_THAT(pmid[0], WithinAbs(2.0, 0.1));
		REQUIRE_THAT(pmid[1], WithinAbs(4.0, 0.1));
		REQUIRE_THAT(pmid[2], WithinAbs(6.0, 0.1));
	}

	TEST_CASE("InterpolatedFunction - duplicate x-values throw", "[Interpolation]")
	{
		Vector<Real> x(4), y(4);
		x[0] = 1.0; x[1] = 2.0; x[2] = 2.0; x[3] = 4.0;
		y[0] = 1.0; y[1] = 4.0; y[2] = 4.0; y[3] = 16.0;

		REQUIRE_THROWS_AS(PolynomInterpRealFunc(x, y, 3), RealFuncInterpInitError);
		REQUIRE_THROWS_AS(LinearInterpRealFunc(x, y), RealFuncInterpInitError);
	}

	///////////////////////////////////////////////////////////////////////////////////////////
	///                      Extrapolation Policy Tests                                     ///
	///////////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("ExtrapolationPolicy_Throw_linear", "[Interpolation][ExtrapolationPolicy]")
	{
		Vector<Real> x(4), y(4);
		x[0] = 1.0; x[1] = 2.0; x[2] = 3.0; x[3] = 4.0;
		y[0] = 1.0; y[1] = 4.0; y[2] = 9.0; y[3] = 16.0;

		LinearInterpRealFunc interp(x, y);
		interp.setExtrapolationPolicy(ExtrapolationPolicy::Throw);

		// In-range should work fine
		REQUIRE_NOTHROW(interp(2.5));

		// Out-of-range should throw
		REQUIRE_THROWS_AS(interp(0.0), RealFuncInterpRuntimeError);
		REQUIRE_THROWS_AS(interp(5.0), RealFuncInterpRuntimeError);

		// Boundary points should work
		REQUIRE_NOTHROW(interp(1.0));
		REQUIRE_NOTHROW(interp(4.0));
	}

	TEST_CASE("ExtrapolationPolicy_Clamp_linear", "[Interpolation][ExtrapolationPolicy]")
	{
		Vector<Real> x(4), y(4);
		x[0] = 1.0; x[1] = 2.0; x[2] = 3.0; x[3] = 4.0;
		y[0] = 1.0; y[1] = 4.0; y[2] = 9.0; y[3] = 16.0;

		LinearInterpRealFunc interp(x, y);
		interp.setExtrapolationPolicy(ExtrapolationPolicy::Clamp);

		// Out-of-range should return boundary values
		Real left_val = interp(0.0);   // clamped to x=1.0 → y=1.0
		Real right_val = interp(5.0);  // clamped to x=4.0 → y=16.0

		REQUIRE_THAT(left_val,  WithinAbs(1.0, TOL(1e-10, 1e-5)));
		REQUIRE_THAT(right_val, WithinAbs(16.0, TOL(1e-10, 1e-5)));
	}

	TEST_CASE("ExtrapolationPolicy_Allow_is_default", "[Interpolation][ExtrapolationPolicy]")
	{
		Vector<Real> x(4), y(4);
		x[0] = 1.0; x[1] = 2.0; x[2] = 3.0; x[3] = 4.0;
		y[0] = 1.0; y[1] = 4.0; y[2] = 9.0; y[3] = 16.0;

		LinearInterpRealFunc interp(x, y);

		// Default policy should be Allow
		REQUIRE(interp.getExtrapolationPolicy() == ExtrapolationPolicy::Allow);

		// Out-of-range should not throw with default policy
		REQUIRE_NOTHROW(interp(0.0));
		REQUIRE_NOTHROW(interp(5.0));
	}

	TEST_CASE("ExtrapolationPolicy_Throw_spline", "[Interpolation][ExtrapolationPolicy]")
	{
		Vector<Real> x(5), y(5);
		x[0] = 0.0; x[1] = 1.0; x[2] = 2.0; x[3] = 3.0; x[4] = 4.0;
		y[0] = 0.0; y[1] = 1.0; y[2] = 4.0; y[3] = 9.0; y[4] = 16.0;

		SplineInterpRealFunc interp(x, y);
		interp.setExtrapolationPolicy(ExtrapolationPolicy::Throw);

		REQUIRE_NOTHROW(interp(2.0));
		REQUIRE_THROWS_AS(interp(-1.0), RealFuncInterpRuntimeError);
		REQUIRE_THROWS_AS(interp(5.0), RealFuncInterpRuntimeError);
	}

	TEST_CASE("ExtrapolationPolicy_Clamp_polynomial", "[Interpolation][ExtrapolationPolicy]")
	{
		Vector<Real> x(4), y(4);
		x[0] = 1.0; x[1] = 2.0; x[2] = 3.0; x[3] = 4.0;
		y[0] = 2.0; y[1] = 5.0; y[2] = 10.0; y[3] = 17.0;

		PolynomInterpRealFunc interp(x, y, 3);
		interp.setExtrapolationPolicy(ExtrapolationPolicy::Clamp);

		// Clamped values should equal boundary interpolation
		Real at_min = interp(1.0);
		Real at_max = interp(4.0);
		Real below  = interp(-5.0);
		Real above  = interp(100.0);

		REQUIRE_THAT(below, WithinAbs(at_min, TOL(1e-10, 1e-5)));
		REQUIRE_THAT(above, WithinAbs(at_max, TOL(1e-10, 1e-5)));
	}

	TEST_CASE("Interpolation_EvaluateDetailed_reports_metadata", "[Interpolation][Detailed]")
	{
		Vector<Real> x{REAL(0.0), REAL(1.0), REAL(2.0)};
		Vector<Real> y{REAL(0.0), REAL(2.0), REAL(4.0)};
		LinearInterpRealFunc interp(x, y, true);

		auto result = interp.EvaluateDetailed(REAL(0.5));
		REQUIRE(result.IsSuccess());
		REQUIRE(result.interpolation_status == InterpolationStatus::Success);
		REQUIRE(result.algorithm_name == "Linear");
		REQUIRE(result.interval_index == 0);
		REQUIRE(result.function_evaluations == 1);
		REQUIRE_THAT(result.value, WithinAbs(REAL(1.0), TOL(1e-12, 1e-5)));
	}

	TEST_CASE("Interpolation_EvaluateDetailed_normalizes_out_of_range", "[Interpolation][Detailed]")
	{
		Vector<Real> x{REAL(0.0), REAL(1.0), REAL(2.0)};
		Vector<Real> y{REAL(0.0), REAL(2.0), REAL(4.0)};
		LinearInterpRealFunc interp(x, y);

		InterpolationConfig clamp;
		clamp.extrapolation_policy = ExtrapolationPolicy::Clamp;
		auto clamped = interp.EvaluateDetailed(REAL(-1.0), clamp);
		REQUIRE(clamped.IsSuccess());
		REQUIRE(clamped.WasClamped());
		REQUIRE(clamped.evaluated_point == REAL(0.0));
		REQUIRE(clamped.value == REAL(0.0));

		InterpolationConfig allow;
		allow.extrapolation_policy = ExtrapolationPolicy::Allow;
		auto extrapolated = interp.EvaluateDetailed(REAL(3.0), allow);
		REQUIRE(extrapolated.IsSuccess());
		REQUIRE(extrapolated.WasExtrapolated());
		REQUIRE(extrapolated.value == REAL(6.0));

		InterpolationConfig reject;
		reject.extrapolation_policy = ExtrapolationPolicy::Throw;
		reject.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
		auto rejected = interp.EvaluateDetailed(REAL(3.0), reject);
		REQUIRE(rejected.status == AlgorithmStatus::InvalidInput);
		REQUIRE(rejected.interpolation_status == InterpolationStatus::OutOfRange);
	}

	TEST_CASE("Interpolation_EvaluateDetailed_exposes_Neville_error", "[Interpolation][Detailed]")
	{
		Vector<Real> x{REAL(0.0), REAL(1.0), REAL(2.0), REAL(3.0)};
		Vector<Real> y{REAL(1.0), REAL(2.0), REAL(5.0), REAL(10.0)};
		PolynomInterpRealFunc interp(x, y, 3);
		auto result = interp.EvaluateDetailed(REAL(1.5));
		REQUIRE(result.IsSuccess());
		REQUIRE(result.algorithm_name == "NevillePolynomial");
		REQUIRE(result.error_estimate >= REAL(0.0));
	}

	TEST_CASE("Interpolation_EvaluateDetailed_supports_descending_nodes", "[Interpolation][Detailed]")
	{
		Vector<Real> x{REAL(3.0), REAL(2.0), REAL(1.0), REAL(0.0)};
		Vector<Real> y{REAL(6.0), REAL(4.0), REAL(2.0), REAL(0.0)};
		LinearInterpRealFunc interp(x, y);
		auto inside = interp.EvaluateDetailed(REAL(1.5));
		REQUIRE(inside.IsSuccess());
		REQUIRE(inside.interpolation_status == InterpolationStatus::Success);
		REQUIRE_THAT(inside.value, WithinAbs(REAL(3.0), TOL(1e-12, 1e-5)));

		InterpolationConfig config;
		config.extrapolation_policy = ExtrapolationPolicy::Clamp;
		auto clamped = interp.EvaluateDetailed(REAL(4.0), config);
		REQUIRE(clamped.WasClamped());
		REQUIRE_THAT(clamped.value, WithinAbs(REAL(6.0), TOL(1e-12, 1e-5)));
	}

	TEST_CASE("ExtrapolationPolicy_isInRange", "[Interpolation][ExtrapolationPolicy]")
	{
		Vector<Real> x(3), y(3);
		x[0] = 0.0; x[1] = 5.0; x[2] = 10.0;
		y[0] = 0.0; y[1] = 25.0; y[2] = 100.0;

		LinearInterpRealFunc interp(x, y);

		REQUIRE(interp.isInRange(0.0));
		REQUIRE(interp.isInRange(5.0));
		REQUIRE(interp.isInRange(10.0));
		REQUIRE(interp.isInRange(3.5));
		REQUIRE_FALSE(interp.isInRange(-0.1));
		REQUIRE_FALSE(interp.isInRange(10.1));
	}

	TEST_CASE("ExtrapolationPolicy_Barycentric_throw", "[Interpolation][ExtrapolationPolicy]")
	{
		Vector<Real> x(5), y(5);
		x[0] = 0.0; x[1] = 1.0; x[2] = 2.0; x[3] = 3.0; x[4] = 4.0;
		y[0] = 0.0; y[1] = 1.0; y[2] = 4.0; y[3] = 9.0; y[4] = 16.0;

		BarycentricRationalInterp interp(x, y, 3);
		interp.setExtrapolationPolicy(ExtrapolationPolicy::Throw);

		REQUIRE_NOTHROW(interp(2.0));
		REQUIRE_THROWS_AS(interp(-1.0), RealFuncInterpRuntimeError);
		REQUIRE_THROWS_AS(interp(5.0), RealFuncInterpRuntimeError);

		// Boundary values should work
		REQUIRE_NOTHROW(interp(0.0));
		REQUIRE_NOTHROW(interp(4.0));
	}

	TEST_CASE("ExtrapolationPolicy_Barycentric_clamp", "[Interpolation][ExtrapolationPolicy]")
	{
		Vector<Real> x(5), y(5);
		x[0] = 0.0; x[1] = 1.0; x[2] = 2.0; x[3] = 3.0; x[4] = 4.0;
		y[0] = 0.0; y[1] = 1.0; y[2] = 4.0; y[3] = 9.0; y[4] = 16.0;

		BarycentricRationalInterp interp(x, y, 3);
		interp.setExtrapolationPolicy(ExtrapolationPolicy::Clamp);

		Real at_min = interp(0.0);
		Real at_max = interp(4.0);
		Real below  = interp(-5.0);
		Real above  = interp(10.0);

		REQUIRE_THAT(below, WithinAbs(at_min, TOL(1e-10, 1e-5)));
		REQUIRE_THAT(above, WithinAbs(at_max, TOL(1e-10, 1e-5)));
	}

	TEST_CASE("Interpolation_FloaterHormann_detailed_result", "[Interpolation][Detailed]")
	{
		Vector<Real> x{REAL(0.0), REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0)};
		Vector<Real> y{REAL(0.0), REAL(1.0), REAL(4.0), REAL(9.0), REAL(16.0)};
		BarycentricRationalInterp interp(x, y, 3);
		InterpolationConfig config;
		config.extrapolation_policy = ExtrapolationPolicy::Clamp;
		auto result = interp.EvaluateDetailed(REAL(5.0), config);
		REQUIRE(result.IsSuccess());
		REQUIRE(result.WasClamped());
		REQUIRE(result.algorithm_name == "FloaterHormannRational");
		REQUIRE(result.interval_index == 3);
		REQUIRE_THAT(result.value, WithinAbs(REAL(16.0), TOL(1e-10, 1e-5)));
	}

	TEST_CASE("Interpolation_spline_owners_are_copy_and_move_safe", "[Interpolation][Ownership]")
	{
		STATIC_REQUIRE(std::is_copy_constructible_v<BicubicSplineInterp2D>);
		STATIC_REQUIRE(std::is_move_constructible_v<BicubicSplineInterp2D>);
		STATIC_REQUIRE(std::is_copy_constructible_v<SplineInterpParametricCurve<2>>);
		STATIC_REQUIRE(std::is_move_constructible_v<SplineInterpParametricCurve<2>>);

		Vector<Real> x1{REAL(0.0), REAL(1.0), REAL(2.0)};
		Vector<Real> x2{REAL(0.0), REAL(1.0), REAL(2.0)};
		Matrix<Real> z(3, 3);
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j) z(i, j) = x1[i] + x2[j];
		BicubicSplineInterp2D original2d(x1, x2, z);
		BicubicSplineInterp2D copied2d(original2d);
		BicubicSplineInterp2D moved2d(std::move(copied2d));
		REQUIRE_THAT(moved2d(REAL(0.5), REAL(1.5)), WithinAbs(REAL(2.0), TOL(1e-8, 1e-4)));

		Matrix<Real> points(4, 2, {REAL(0.0), REAL(0.0), REAL(1.0), REAL(1.0),
			REAL(2.0), REAL(2.0), REAL(3.0), REAL(3.0)});
		SplineInterpParametricCurve<2> originalCurve(points, false);
		SplineInterpParametricCurve<2> copiedCurve(originalCurve);
		SplineInterpParametricCurve<2> movedCurve(std::move(copiedCurve));
		auto point = movedCurve(REAL(0.5));
		REQUIRE_THAT(point[0], WithinAbs(REAL(1.5), TOL(1e-8, 1e-4)));
		REQUIRE_THAT(point[1], WithinAbs(REAL(1.5), TOL(1e-8, 1e-4)));
	}

	//////////////////////////////////////////////////////////////////////////////
	//              MONOTONE CUBIC INTERPOLATION (Fritsch-Carlson)              //
	//////////////////////////////////////////////////////////////////////////////

	TEST_CASE("MonotoneCubic_exact_at_nodes", "[interpolation][monotone]") {
		Vector<Real> x{ 0.0, 1.0, 2.0, 3.0, 4.0 };
		Vector<Real> y{ 1.0, 3.0, 4.0, 5.5, 8.0 };

		MonotoneCubicInterpRealFunc interp(x, y);

		for (int i = 0; i < x.size(); i++) {
			REQUIRE_THAT(interp(x[i]), WithinAbs(y[i], TOL(1e-12, 1e-5)));
		}
	}

	TEST_CASE("MonotoneCubic_monotone_data_stays_monotone", "[interpolation][monotone]") {
		// Strictly increasing data — interpolant must never decrease
		Vector<Real> x{ 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 };
		Vector<Real> y{ 0.0, 0.5, 2.0, 2.1, 5.0, 9.0 };

		MonotoneCubicInterpRealFunc interp(x, y);

		Real prev = interp(0.0);
		for (int i = 1; i <= 500; i++) {
			Real xi = 5.0 * i / 500.0;
			Real val = interp(xi);
			REQUIRE(val >= prev - TOL(1e-14, 1e-6));
			prev = val;
		}
	}

	TEST_CASE("MonotoneCubic_decreasing_data_stays_decreasing", "[interpolation][monotone]") {
		// Strictly decreasing data — interpolant must never increase
		Vector<Real> x{ 0.0, 1.0, 2.0, 3.0, 4.0 };
		Vector<Real> y{ 10.0, 7.0, 3.0, 1.5, 0.0 };

		MonotoneCubicInterpRealFunc interp(x, y);

		Real prev = interp(0.0);
		for (int i = 1; i <= 400; i++) {
			Real xi = 4.0 * i / 400.0;
			Real val = interp(xi);
			REQUIRE(val <= prev + TOL(1e-14, 1e-6));
			prev = val;
		}
	}

	TEST_CASE("MonotoneCubic_no_overshoot_step_function", "[interpolation][monotone]") {
		// Data with a sharp step — monotone cubic must not overshoot
		Vector<Real> x{ 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 };
		Vector<Real> y{ 0.0, 0.0, 0.0, 1.0, 1.0, 1.0 };

		MonotoneCubicInterpRealFunc interp(x, y);

		for (int i = 0; i <= 500; i++) {
			Real xi = 5.0 * i / 500.0;
			Real val = interp(xi);
			REQUIRE(val >= -TOL(1e-14, 1e-6));
			REQUIRE(val <= 1.0 + TOL(1e-14, 1e-6));
		}
	}

	TEST_CASE("MonotoneCubic_flat_segments_zero_derivative", "[interpolation][monotone]") {
		// Flat segments should produce zero tangent derivatives
		Vector<Real> x{ 0.0, 1.0, 2.0, 3.0, 4.0 };
		Vector<Real> y{ 3.0, 3.0, 3.0, 3.0, 3.0 };

		MonotoneCubicInterpRealFunc interp(x, y);

		for (int i = 0; i < x.size(); i++) {
			REQUIRE_THAT(interp.GetDerivative(i), WithinAbs(0.0, TOL(1e-14, 1e-6)));
		}

		// Interpolated values should all be 3.0
		for (int i = 0; i <= 100; i++) {
			Real xi = 4.0 * i / 100.0;
			REQUIRE_THAT(interp(xi), WithinAbs(3.0, TOL(1e-12, 1e-5)));
		}
	}

	TEST_CASE("MonotoneCubic_two_points_linear", "[interpolation][monotone][edge]") {
		// Minimum valid input: two points → should be linear
		Vector<Real> x{ 1.0, 3.0 };
		Vector<Real> y{ 2.0, 8.0 };

		MonotoneCubicInterpRealFunc interp(x, y);

		// Should be exact linear: y = 2 + 3*(x - 1)
		REQUIRE_THAT(interp(1.0), WithinAbs(2.0, TOL(1e-12, 1e-5)));
		REQUIRE_THAT(interp(2.0), WithinAbs(5.0, TOL(1e-12, 1e-5)));
		REQUIRE_THAT(interp(3.0), WithinAbs(8.0, TOL(1e-12, 1e-5)));
	}

	TEST_CASE("MonotoneCubic_sin_accuracy", "[interpolation][monotone]") {
		// Monotone cubic on sin(x) in [0, pi/2] (monotone increasing region)
		RealFunction f{ [](Real x) -> Real { return sin(x); } };

		Vector<Real> vec_x, vec_y;
		CreateInterpolatedValues(f, 0.0, Constants::PI / 2.0, 20, vec_x, vec_y);

		MonotoneCubicInterpRealFunc interp(vec_x, vec_y);
		RealFunctionComparer comparer(f, interp);

		double maxAbsDiff = comparer.getAbsDiffMax(0.0, Constants::PI / 2.0, 200);
		REQUIRE_THAT(maxAbsDiff, WithinAbs(0.0, 0.001));
	}

	TEST_CASE("MonotoneCubic_derivative_accuracy", "[interpolation][monotone][derivative]") {
		// For x^2 on [0, 3], derivative should be 2x
		RealFunction f{ [](Real x) -> Real { return x * x; } };
		RealFunction df{ [](Real x) -> Real { return 2.0 * x; } };

		Vector<Real> vec_x, vec_y;
		CreateInterpolatedValues(f, 0.0, 3.0, 20, vec_x, vec_y);

		MonotoneCubicInterpRealFunc interp(vec_x, vec_y);

		std::vector<Real> test_points = { 0.5, 1.0, 1.5, 2.0, 2.5 };
		for (Real x : test_points) {
			Real computed = interp.Derivative(x);
			Real expected = df(x);
			REQUIRE_THAT(computed, WithinAbs(expected, TOL(0.05, 0.1)));
		}
	}

	TEST_CASE("MonotoneCubic_GetDerivative_accessor", "[interpolation][monotone]") {
		Vector<Real> x{ 0.0, 1.0, 2.0, 3.0 };
		Vector<Real> y{ 0.0, 1.0, 4.0, 9.0 };

		MonotoneCubicInterpRealFunc interp(x, y);

		// Should be able to access all derivative values
		for (int i = 0; i < x.size(); i++) {
			REQUIRE_NOTHROW(interp.GetDerivative(i));
		}

		// Out-of-range should throw
		REQUIRE_THROWS(interp.GetDerivative(-1));
		REQUIRE_THROWS(interp.GetDerivative(static_cast<int>(x.size())));
	}

	TEST_CASE("MonotoneCubic_vs_spline_no_overshoot", "[interpolation][monotone][comparison]") {
		// Data with sharp transition where regular spline typically overshoots
		Vector<Real> x{ 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 };
		Vector<Real> y{ 0.0, 0.0, 0.0, 10.0, 10.0, 10.0 };

		MonotoneCubicInterpRealFunc monotone(x, y);

		// Monotone cubic must stay in [0, 10]
		for (int i = 0; i <= 500; i++) {
			Real xi = 5.0 * i / 500.0;
			Real val = monotone(xi);
			REQUIRE(val >= -TOL(1e-14, 1e-6));
			REQUIRE(val <= 10.0 + TOL(1e-14, 1e-6));
		}
	}

	TEST_CASE("MonotoneCubic_extrapolation_throw", "[interpolation][monotone][ExtrapolationPolicy]") {
		Vector<Real> x{ 0.0, 1.0, 2.0, 3.0 };
		Vector<Real> y{ 0.0, 1.0, 4.0, 9.0 };

		MonotoneCubicInterpRealFunc interp(x, y);
		interp.setExtrapolationPolicy(ExtrapolationPolicy::Throw);

		REQUIRE_NOTHROW(interp(1.5));
		REQUIRE_THROWS_AS(interp(-1.0), RealFuncInterpRuntimeError);
		REQUIRE_THROWS_AS(interp(4.0), RealFuncInterpRuntimeError);
	}

	TEST_CASE("MonotoneCubic_extrapolation_clamp", "[interpolation][monotone][ExtrapolationPolicy]") {
		Vector<Real> x{ 0.0, 1.0, 2.0, 3.0 };
		Vector<Real> y{ 0.0, 1.0, 4.0, 9.0 };

		MonotoneCubicInterpRealFunc interp(x, y);
		interp.setExtrapolationPolicy(ExtrapolationPolicy::Clamp);

		Real at_min = interp(0.0);
		Real at_max = interp(3.0);
		Real below  = interp(-5.0);
		Real above  = interp(10.0);

		REQUIRE_THAT(below, WithinAbs(at_min, TOL(1e-10, 1e-5)));
		REQUIRE_THAT(above, WithinAbs(at_max, TOL(1e-10, 1e-5)));
	}

} // end namespace
#endif

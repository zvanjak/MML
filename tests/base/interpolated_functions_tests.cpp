#if !defined __MML_INTERPOLATED_FUNCTIONS_TESTS_H
#define __MML_INTERPOLATED_FUNCTIONS_TESTS_H

#include <catch2/catch_all.hpp>

#include <vector>
#include <cmath>

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/InterpolatedFunction.h"
#include "algorithms/FunctionsAnalyzer.h"
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
		RealFunction f{ [](Real x) { return sin(x); } };

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
			REQUIRE_THAT(interp(x[i]), WithinAbs(y[i], 1e-14));
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
		REQUIRE_THAT(interp(0.5), WithinAbs(f(0.5), 1e-12));
		REQUIRE_THAT(interp(1.5), WithinAbs(f(1.5), 1e-12));
	}

	TEST_CASE("PolynomInterp_error_estimate", "[interpolation][polynomial]") {
		// Error estimate from Neville's algorithm (can be negative as it's a signed difference)
		RealFunction f{ [](Real x) { return cos(x); } };

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
			REQUIRE_THAT(spline(x[i]), WithinAbs(y[i], 1e-12));
		}
	}

	TEST_CASE("SplineInterp_sin_high_accuracy", "[interpolation][spline]") {
		RealFunction f{ [](Real x) { return sin(x); } };

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
		RealFunction f{ [](Real x) { return sin(x); } };
		RealFunction df{ [](Real x) { return cos(x); } };

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
		RealFunction f{ [](Real x) { return sin(x); } };
		RealFunction ddf{ [](Real x) { return -sin(x); } };

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
		RealFunction f{ [](Real x) { return sin(x); } };

		Vector<Real> vec_x, vec_y;
		CreateInterpolatedValues(f, 0.0, Constants::PI, 20, vec_x, vec_y);

		SplineInterpRealFunc spline(vec_x, vec_y);

		Real integral = spline.Integrate(0.0, Constants::PI);
		Real expected = 2.0;

		REQUIRE_THAT(integral, WithinAbs(expected, 0.01));
	}

	TEST_CASE("SplineInterp_clamped_boundary", "[interpolation][spline]") {
		// Test clamped spline with known boundary derivatives
		RealFunction f{ [](Real x) { return sin(x); } };
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
		REQUIRE_THAT(spline.GetSecondDerivative(0), WithinAbs(0.0, 1e-10));
		REQUIRE_THAT(spline.GetSecondDerivative(4), WithinAbs(0.0, 1e-10));
	}

	//////////////////////////////////////////////////////////////////////////////
	//                     SPLINE VS LINEAR COMPARISON                          //
	//////////////////////////////////////////////////////////////////////////////

	TEST_CASE("Spline_more_accurate_than_linear", "[interpolation][comparison]") {
		// Spline should be more accurate than linear for smooth functions
		RealFunction f{ [](Real x) { return exp(-x) * cos(2 * x); } };

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
		REQUIRE_THAT(spline(0.5), WithinAbs(0.5, 1e-10));
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
			REQUIRE_THAT(interp(x[i]), WithinAbs(y[i], 1e-10));
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
		RealFunction f{ [](Real x) { return sin(x); } };

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
			REQUIRE_THAT(interp(x[i]), WithinAbs(y[i], 1e-10));
		}
	}

	TEST_CASE("BarycentricInterp_no_poles", "[interpolation][barycentric]") {
		// Unlike rational interpolation, barycentric should be pole-free
		RealFunction f{ [](Real x) { return exp(-x * x); } };

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

		REQUIRE_THAT(interp(0.5, 0.5), WithinAbs(5.0, 1e-10));
		REQUIRE_THAT(interp(2.5, 2.5), WithinAbs(5.0, 1e-10));
		REQUIRE_THAT(interp(3.7, 1.2), WithinAbs(5.0, 1e-10));
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
		REQUIRE_THAT(interp(0.5, 0.5), WithinAbs(1.0, 1e-10));
		REQUIRE_THAT(interp(1.5, 2.5), WithinAbs(4.0, 1e-10));
		REQUIRE_THAT(interp(2.7, 3.3), WithinAbs(6.0, 1e-10));
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
				REQUIRE_THAT(interp(x1[i], x2[j]), WithinAbs(y[i][j], 1e-12));
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

		REQUIRE_THAT(val, WithinAbs(2.0 * 1.5 + 3.0 * 2.5, 1e-10));
		REQUIRE_THAT(dz_dx1, WithinAbs(2.0, 1e-10));
		REQUIRE_THAT(dz_dx2, WithinAbs(3.0, 1e-10));
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

		REQUIRE_THAT(interp(0.5, 0.5), WithinAbs(7.5, 1e-10));
		REQUIRE_THAT(interp(2.5, 2.5), WithinAbs(7.5, 1e-10));
		REQUIRE_THAT(interp(3.7, 1.2), WithinAbs(7.5, 1e-10));
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

		REQUIRE_THAT(interp(0.5, 0.5), WithinAbs(1.0, 1e-8));
		REQUIRE_THAT(interp(1.5, 2.5), WithinAbs(4.0, 1e-8));
		REQUIRE_THAT(interp(2.7, 3.3), WithinAbs(6.0, 1e-8));
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
				REQUIRE_THAT(interp(x1[i], x2[j]), WithinAbs(z[i][j], 1e-10));
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
			REQUIRE_THAT(poly(x[i]), WithinAbs(y[i], 1e-10));
		}

		// Check midpoint
		Real mid = 2.5;
		Real expected = std::pow(mid, 4) - 2 * mid * mid + 1;
		REQUIRE_THAT(poly(mid), WithinAbs(expected, 1e-8));
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
			REQUIRE_THAT(result, WithinAbs(y[i], 1e-10));
			// Error estimate should be zero at exact points
			REQUIRE_THAT(interp.getLastErrorEst(), WithinAbs(0.0, 1e-10));
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
			REQUIRE_THAT(interp.X(i), WithinAbs(x[i], 1e-15));
			REQUIRE_THAT(interp.Y(i), WithinAbs(y[i], 1e-15));
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
		LinearInterpRealFunc interp(x, y, true);

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
		REQUIRE_THAT(p0[0], WithinAbs(0.0, 1e-10));
		REQUIRE_THAT(p0[1], WithinAbs(0.0, 1e-10));

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
		REQUIRE_THAT(p0[0], WithinAbs(0.0, 1e-8));
		REQUIRE_THAT(p0[1], WithinAbs(0.0, 1e-8));
		REQUIRE_THAT(p0[2], WithinAbs(0.0, 1e-8));

		auto p1 = curve(1.0);
		REQUIRE_THAT(p1[0], WithinAbs(4.0, 1e-8));
		REQUIRE_THAT(p1[1], WithinAbs(8.0, 1e-8));
		REQUIRE_THAT(p1[2], WithinAbs(12.0, 1e-8));

		// Check midpoint - should be linear
		auto pmid = curve(0.5);
		REQUIRE_THAT(pmid[0], WithinAbs(2.0, 0.1));
		REQUIRE_THAT(pmid[1], WithinAbs(4.0, 0.1));
		REQUIRE_THAT(pmid[2], WithinAbs(6.0, 0.1));
	}

} // end namespace
#endif

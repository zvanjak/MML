#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/ComplexFunction.h"
#include "core/Derivation.h"
#include "algorithms/RootFinding.h"
#include "core/ComplexAnalysis.h"
#endif

using namespace MML;
using namespace MML::Testing;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace MML::Tests::Core::ComplexAnalysisTests
{
	/*********************************************************************/
	/*****     IComplexFunction / ComplexFunction wrappers           *****/
	/*********************************************************************/

	TEST_CASE("ComplexFunction_pointer_wrapper", "[complex][interface]")
	{
		// f(z) = z^2
		ComplexFunction f([](Complex z) -> Complex { return z * z; });
		Complex z(3.0, 4.0);
		Complex result = f(z);
		// (3+4i)^2 = 9 + 24i - 16 = -7 + 24i
		REQUIRE_THAT(result.real(), WithinAbs(-7.0, TOL(1e-12, 1e-5)));
		REQUIRE_THAT(result.imag(), WithinAbs(24.0, TOL(1e-12, 1e-5)));
	}

	TEST_CASE("ComplexFunctionFromStdFunc_lambda_wrapper", "[complex][interface]")
	{
		ComplexFunctionFromStdFunc f([](Complex z) -> Complex { return z * z * z; });
		Complex z(1.0, 1.0);
		Complex result = f(z);
		// (1+i)^3 = (1+i)(1+2i-1) = (1+i)(2i) = 2i + 2i^2 = -2 + 2i
		REQUIRE_THAT(result.real(), WithinAbs(-2.0, TOL(1e-12, 1e-5)));
		REQUIRE_THAT(result.imag(), WithinAbs(2.0, TOL(1e-12, 1e-5)));
	}

	TEST_CASE("RealToComplexFunction_wrapper", "[complex][interface]")
	{
		RealToComplexFunctionFromStdFunc gamma([](Real t) -> Complex {
			return Complex(std::cos(t), std::sin(t));
		});
		Complex result = gamma(Constants::PI / 2.0);
		REQUIRE_THAT(result.real(), WithinAbs(0.0, TOL(1e-12, 1e-5)));
		REQUIRE_THAT(result.imag(), WithinAbs(1.0, TOL(1e-12, 1e-5)));
	}

	/*********************************************************************/
	/*****            Complex Derivatives                            *****/
	/*********************************************************************/

	TEST_CASE("NDer1Complex_z_squared", "[complex][derivation]")
	{
		// f(z) = z^2, f'(z) = 2z
		ComplexFunctionFromStdFunc f([](Complex z) -> Complex { return z * z; });
		Complex z(2.0, 1.0);
		Complex expected(4.0, 2.0);  // 2*(2+i) = 4+2i

		Complex der = Derivation::NDer1Complex(f, z);
		REQUIRE_THAT(der.real(), WithinAbs(expected.real(), TOL(1e-6, 5e-3)));
		REQUIRE_THAT(der.imag(), WithinAbs(expected.imag(), TOL(1e-6, 5e-3)));
	}

	TEST_CASE("NDer2Complex_z_squared", "[complex][derivation]")
	{
		ComplexFunctionFromStdFunc f([](Complex z) -> Complex { return z * z; });
		Complex z(2.0, 1.0);
		Complex expected(4.0, 2.0);

		Complex der = Derivation::NDer2Complex(f, z);
		REQUIRE_THAT(der.real(), WithinAbs(expected.real(), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(der.imag(), WithinAbs(expected.imag(), TOL(1e-10, 1e-5)));
	}

	TEST_CASE("NDer4Complex_z_squared", "[complex][derivation]")
	{
		ComplexFunctionFromStdFunc f([](Complex z) -> Complex { return z * z; });
		Complex z(2.0, 1.0);
		Complex expected(4.0, 2.0);

		Complex der = Derivation::NDer4Complex(f, z);
		REQUIRE_THAT(der.real(), WithinAbs(expected.real(), TOL(1e-12, 1e-5)));
		REQUIRE_THAT(der.imag(), WithinAbs(expected.imag(), TOL(1e-12, 1e-5)));
	}

	TEST_CASE("NDer6Complex_z_squared", "[complex][derivation]")
	{
		ComplexFunctionFromStdFunc f([](Complex z) -> Complex { return z * z; });
		Complex z(2.0, 1.0);
		Complex expected(4.0, 2.0);

		Complex der = Derivation::NDer6Complex(f, z);
		REQUIRE_THAT(der.real(), WithinAbs(expected.real(), TOL(1e-12, 1e-5)));
		REQUIRE_THAT(der.imag(), WithinAbs(expected.imag(), TOL(1e-12, 1e-5)));
	}

	TEST_CASE("NDer4Complex_exp", "[complex][derivation]")
	{
		// f(z) = e^z, f'(z) = e^z
		ComplexFunctionFromStdFunc f([](Complex z) -> Complex { return std::exp(z); });
		Complex z(1.0, Constants::PI / 4.0);
		Complex expected = std::exp(z);

		Complex der = Derivation::NDer4Complex(f, z);
		REQUIRE_THAT(der.real(), WithinAbs(expected.real(), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(der.imag(), WithinAbs(expected.imag(), TOL(1e-10, 1e-5)));
	}

	TEST_CASE("NDer2Complex_sin", "[complex][derivation]")
	{
		// f(z) = sin(z), f'(z) = cos(z)
		ComplexFunctionFromStdFunc f([](Complex z) -> Complex { return std::sin(z); });
		Complex z(1.0, 0.5);
		Complex expected = std::cos(z);

		Complex der = Derivation::NDer2Complex(f, z);
		REQUIRE_THAT(der.real(), WithinAbs(expected.real(), TOL(1e-9, 1e-4)));
		REQUIRE_THAT(der.imag(), WithinAbs(expected.imag(), TOL(1e-9, 1e-4)));
	}

	TEST_CASE("NDer4ComplexDetailed_returns_error", "[complex][derivation]")
	{
		ComplexFunctionFromStdFunc f([](Complex z) -> Complex { return z * z; });
		Complex z(1.0, 1.0);

		DerivativeConfig config;
		config.estimate_error = true;
		auto result = Derivation::NDer4ComplexDetailed(f, z, config);

		REQUIRE(result.status == AlgorithmStatus::Success);
		REQUIRE(result.error >= 0.0);
		// Error should be very small for a polynomial
		REQUIRE(result.error < TOL(1e-8, 1e-4));
		// Value should be 2z = 2+2i
		REQUIRE_THAT(result.value.real(), WithinAbs(2.0, TOL(1e-10, 1e-5)));
		REQUIRE_THAT(result.value.imag(), WithinAbs(2.0, TOL(1e-10, 1e-5)));
	}

	/*********************************************************************/
	/*****         Complex Root Finding - Newton                     *****/
	/*********************************************************************/

	TEST_CASE("NewtonComplex_z_squared_plus_1", "[complex][rootfinding]")
	{
		// f(z) = z^2 + 1, roots at z = ±i
		ComplexFunctionFromStdFunc f([](Complex z) -> Complex { return z * z + Complex(1.0, 0.0); });

		// Starting near +i
		auto result = RootFinding::FindRootNewtonComplex(f, Complex(0.1, 1.2));
		REQUIRE(result.converged);
		REQUIRE_THAT(result.root.real(), WithinAbs(0.0, TOL(1e-8, 1e-4)));
		REQUIRE_THAT(result.root.imag(), WithinAbs(1.0, TOL(1e-8, 1e-4)));
	}

	TEST_CASE("NewtonComplex_z_squared_plus_1_minus_i", "[complex][rootfinding]")
	{
		ComplexFunctionFromStdFunc f([](Complex z) -> Complex { return z * z + Complex(1.0, 0.0); });

		// Starting near -i
		auto result = RootFinding::FindRootNewtonComplex(f, Complex(-0.1, -0.8));
		REQUIRE(result.converged);
		REQUIRE_THAT(result.root.real(), WithinAbs(0.0, TOL(1e-8, 1e-4)));
		REQUIRE_THAT(result.root.imag(), WithinAbs(-1.0, TOL(1e-8, 1e-4)));
	}

	TEST_CASE("NewtonComplex_cubic", "[complex][rootfinding]")
	{
		// f(z) = z^3 - 1, roots at z = 1, e^(2πi/3), e^(4πi/3)
		ComplexFunctionFromStdFunc f([](Complex z) -> Complex { return z * z * z - Complex(1.0, 0.0); });

		auto result = RootFinding::FindRootNewtonComplex(f, Complex(1.1, 0.1));
		REQUIRE(result.converged);
		// Should converge to real root z = 1
		Real abs_fz = std::abs(result.function_value);
		REQUIRE(abs_fz < TOL(1e-10, 1e-5));
	}

	/*********************************************************************/
	/*****         Complex Root Finding - Muller                     *****/
	/*********************************************************************/

	TEST_CASE("Muller_z_squared_plus_1", "[complex][rootfinding]")
	{
		// f(z) = z^2 + 1, roots at z = ±i
		ComplexFunctionFromStdFunc f([](Complex z) -> Complex { return z * z + Complex(1.0, 0.0); });

		// Muller can find complex roots from real starting points
		auto result = RootFinding::FindRootMuller(f, Complex(0.5, 0.5));
		REQUIRE(result.converged);
		Real abs_fz = std::abs(result.function_value);
		REQUIRE(abs_fz < TOL(1e-8, 1e-4));
		// Root should be +i or -i
		REQUIRE_THAT(std::abs(result.root), WithinAbs(1.0, TOL(1e-8, 1e-4)));
	}

	TEST_CASE("Muller_three_point_cubic", "[complex][rootfinding]")
	{
		// f(z) = z^3 - 1
		ComplexFunctionFromStdFunc f([](Complex z) -> Complex { return z * z * z - Complex(1.0, 0.0); });

		auto result = RootFinding::FindRootMuller(f, Complex(0.5, 0.5), Complex(0.8, 0.2), Complex(1.1, -0.1));
		REQUIRE(result.converged);
		Real abs_fz = std::abs(result.function_value);
		REQUIRE(abs_fz < TOL(1e-8, 1e-4));
	}

	/*********************************************************************/
	/*****            Contour Classes                                *****/
	/*********************************************************************/

	TEST_CASE("CircleContour_unit_circle", "[complex][contour]")
	{
		ComplexAnalysis::CircleContour unit_circle(Complex(0.0, 0.0), 1.0);

		// At t=0: γ(0) = 1
		Complex z0 = unit_circle(0.0);
		REQUIRE_THAT(z0.real(), WithinAbs(1.0, TOL(1e-14, 1e-5)));
		REQUIRE_THAT(z0.imag(), WithinAbs(0.0, TOL(1e-14, 1e-5)));

		// At t=π/2: γ(π/2) = i
		Complex z1 = unit_circle(Constants::PI / 2.0);
		REQUIRE_THAT(z1.real(), WithinAbs(0.0, TOL(1e-14, 1e-5)));
		REQUIRE_THAT(z1.imag(), WithinAbs(1.0, TOL(1e-14, 1e-5)));

		// Derivative at t=0: γ'(0) = i
		Complex dz0 = unit_circle.derivative(0.0);
		REQUIRE_THAT(dz0.real(), WithinAbs(0.0, TOL(1e-14, 1e-5)));
		REQUIRE_THAT(dz0.imag(), WithinAbs(1.0, TOL(1e-14, 1e-5)));
	}

	TEST_CASE("LineSegmentContour_basic", "[complex][contour]")
	{
		Complex z1(0.0, 0.0);
		Complex z2(2.0, 2.0);
		ComplexAnalysis::LineSegmentContour seg(z1, z2);

		// At t=0: start
		Complex s = seg(0.0);
		REQUIRE_THAT(s.real(), WithinAbs(0.0, TOL(1e-14, 1e-5)));
		REQUIRE_THAT(s.imag(), WithinAbs(0.0, TOL(1e-14, 1e-5)));

		// At t=0.5: midpoint
		Complex m = seg(0.5);
		REQUIRE_THAT(m.real(), WithinAbs(1.0, TOL(1e-14, 1e-5)));
		REQUIRE_THAT(m.imag(), WithinAbs(1.0, TOL(1e-14, 1e-5)));

		// At t=1: end
		Complex e = seg(1.0);
		REQUIRE_THAT(e.real(), WithinAbs(2.0, TOL(1e-14, 1e-5)));
		REQUIRE_THAT(e.imag(), WithinAbs(2.0, TOL(1e-14, 1e-5)));

		// Derivative is constant: z2 - z1 = 2+2i
		Complex d = seg.derivative(0.5);
		REQUIRE_THAT(d.real(), WithinAbs(2.0, TOL(1e-14, 1e-5)));
		REQUIRE_THAT(d.imag(), WithinAbs(2.0, TOL(1e-14, 1e-5)));
	}

	/*********************************************************************/
	/*****          Contour Integration                              *****/
	/*********************************************************************/

	TEST_CASE("ContourIntegral_z_around_origin", "[complex][integration]")
	{
		// ∮ z dz around unit circle = 0 (z is analytic)
		ComplexFunctionFromStdFunc f([](Complex z) -> Complex { return z; });
		ComplexAnalysis::CircleContour unit_circle(Complex(0.0, 0.0), 1.0);

		auto result = ComplexAnalysis::ContourIntegral(f, unit_circle);
		REQUIRE_THAT(std::abs(result.value), WithinAbs(0.0, TOL(1e-8, 1e-4)));
	}

	TEST_CASE("ContourIntegral_one_over_z_around_origin", "[complex][integration]")
	{
		// ∮ (1/z) dz around unit circle = 2πi
		ComplexFunctionFromStdFunc f([](Complex z) -> Complex { return REAL(1.0) / z; });
		ComplexAnalysis::CircleContour unit_circle(Complex(0.0, 0.0), 1.0);

		auto result = ComplexAnalysis::ContourIntegral(f, unit_circle);
		Complex expected(0.0, 2.0 * Constants::PI);

		REQUIRE_THAT(result.value.real(), WithinAbs(expected.real(), TOL(1e-8, 1e-4)));
		REQUIRE_THAT(result.value.imag(), WithinAbs(expected.imag(), TOL(1e-8, 1e-4)));
	}

	TEST_CASE("ContourIntegral_z_squared_analytic", "[complex][integration]")
	{
		// ∮ z^2 dz = 0 (analytic function, Cauchy theorem)
		ComplexFunctionFromStdFunc f([](Complex z) -> Complex { return z * z; });
		ComplexAnalysis::CircleContour circle(Complex(1.0, 0.0), 2.0);

		auto result = ComplexAnalysis::ContourIntegral(f, circle);
		REQUIRE_THAT(std::abs(result.value), WithinAbs(0.0, TOL(1e-7, 1e-3)));
	}

	TEST_CASE("ContourIntegral_line_segment", "[complex][integration]")
	{
		// ∫ z dz from 0 to 1+i along a line segment
		// Antiderivative: z^2/2, so integral = (1+i)^2/2 - 0 = (2i)/2 = i
		ComplexFunctionFromStdFunc f([](Complex z) -> Complex { return z; });
		ComplexAnalysis::LineSegmentContour seg(Complex(0.0, 0.0), Complex(1.0, 1.0));

		auto result = ComplexAnalysis::ContourIntegral(f, seg);
		REQUIRE_THAT(result.value.real(), WithinAbs(0.0, TOL(1e-10, 1e-5)));
		REQUIRE_THAT(result.value.imag(), WithinAbs(1.0, TOL(1e-10, 1e-5)));
	}

	/*********************************************************************/
	/*****             Winding Number                                *****/
	/*********************************************************************/

	TEST_CASE("WindingNumber_inside", "[complex][winding]")
	{
		// Winding number of unit circle around origin = 1
		ComplexAnalysis::CircleContour unit_circle(Complex(0.0, 0.0), 1.0);
		int n = ComplexAnalysis::WindingNumber(unit_circle, Complex(0.0, 0.0));
		REQUIRE(n == 1);
	}

	TEST_CASE("WindingNumber_outside", "[complex][winding]")
	{
		// Winding number of unit circle around point outside = 0
		ComplexAnalysis::CircleContour unit_circle(Complex(0.0, 0.0), 1.0);
		int n = ComplexAnalysis::WindingNumber(unit_circle, Complex(5.0, 0.0));
		REQUIRE(n == 0);
	}

	TEST_CASE("WindingNumber_off_center_inside", "[complex][winding]")
	{
		// Circle centered at (2,0) with radius 3, point at origin is inside
		ComplexAnalysis::CircleContour circle(Complex(2.0, 0.0), 3.0);
		int n = ComplexAnalysis::WindingNumber(circle, Complex(0.0, 0.0));
		REQUIRE(n == 1);
	}

	/*********************************************************************/
	/*****          Cauchy Integral Formula                           *****/
	/*********************************************************************/

	TEST_CASE("CauchyIntegralFormula_exp", "[complex][cauchy]")
	{
		// f(z) = e^z, evaluate at z0 = 0 using Cauchy integral
		// Should give f(0) = 1
		ComplexFunctionFromStdFunc f([](Complex z) -> Complex { return std::exp(z); });
		ComplexAnalysis::CircleContour contour(Complex(0.0, 0.0), 1.0);

		Complex result = ComplexAnalysis::CauchyIntegralFormula(f, contour, Complex(0.0, 0.0));
		REQUIRE_THAT(result.real(), WithinAbs(1.0, 1e-6));
		REQUIRE_THAT(result.imag(), WithinAbs(0.0, 1e-6));
	}

	TEST_CASE("CauchyIntegralFormula_z_squared", "[complex][cauchy]")
	{
		// f(z) = z^2, evaluate at z0 = 1+i
		// Should give (1+i)^2 = 2i
		ComplexFunctionFromStdFunc f([](Complex z) -> Complex { return z * z; });
		ComplexAnalysis::CircleContour contour(Complex(1.0, 1.0), 2.0);

		Complex result = ComplexAnalysis::CauchyIntegralFormula(f, contour, Complex(1.0, 1.0));
		REQUIRE_THAT(result.real(), WithinAbs(0.0, 1e-6));
		REQUIRE_THAT(result.imag(), WithinAbs(2.0, 1e-6));
	}

	TEST_CASE("CauchyDerivative_exp_first", "[complex][cauchy]")
	{
		// f(z) = e^z, f'(0) = 1
		ComplexFunctionFromStdFunc f([](Complex z) -> Complex { return std::exp(z); });
		ComplexAnalysis::CircleContour contour(Complex(0.0, 0.0), 1.0);

		Complex result = ComplexAnalysis::CauchyDerivative(f, contour, Complex(0.0, 0.0), 1);
		REQUIRE_THAT(result.real(), WithinAbs(1.0, 1e-5));
		REQUIRE_THAT(result.imag(), WithinAbs(0.0, 1e-5));
	}

	TEST_CASE("CauchyDerivative_exp_second", "[complex][cauchy]")
	{
		// f(z) = e^z, f''(z) = e^z, f''(0) = 1
		// Uses exp to avoid Simpson aliasing with polynomial integrands
		ComplexFunctionFromStdFunc f([](Complex z) -> Complex { return std::exp(z); });
		ComplexAnalysis::CircleContour contour(Complex(0.0, 0.0), 1.0);

		Complex result = ComplexAnalysis::CauchyDerivative(f, contour, Complex(0.0, 0.0), 2);
		REQUIRE_THAT(result.real(), WithinAbs(1.0, 1e-4));
		REQUIRE_THAT(result.imag(), WithinAbs(0.0, 1e-4));
	}

	/*********************************************************************/
	/*****            Residue Computation                            *****/
	/*********************************************************************/

	TEST_CASE("Residue_one_over_z", "[complex][residue]")
	{
		// f(z) = 1/z has Res(f, 0) = 1
		ComplexFunctionFromStdFunc f([](Complex z) -> Complex { return REAL(1.0) / z; });

		Complex res = ComplexAnalysis::Residue(f, Complex(0.0, 0.0), 0.5);
		REQUIRE_THAT(res.real(), WithinAbs(1.0, 1e-6));
		REQUIRE_THAT(res.imag(), WithinAbs(0.0, 1e-6));
	}

	TEST_CASE("Residue_one_over_z_squared_plus_1_at_i", "[complex][residue]")
	{
		// f(z) = 1/(z^2 + 1) = 1/((z-i)(z+i))
		// Res(f, i) = 1/(i + i) = 1/(2i) = -i/2
		ComplexFunctionFromStdFunc f([](Complex z) -> Complex { return REAL(1.0) / (z * z + REAL(1.0)); });

		Complex res = ComplexAnalysis::Residue(f, Complex(0.0, 1.0), 0.5);
		REQUIRE_THAT(res.real(), WithinAbs(0.0, 1e-5));
		REQUIRE_THAT(res.imag(), WithinAbs(-0.5, 1e-5));
	}

	TEST_CASE("ResidueSimplePole_one_over_z", "[complex][residue]")
	{
		// f(z) = 1/z, Res(f, 0) = lim_{z→0} z·(1/z) = 1
		ComplexFunctionFromStdFunc f([](Complex z) -> Complex { return REAL(1.0) / z; });

		Complex res = ComplexAnalysis::ResidueSimplePole(f, Complex(0.0, 0.0));
		REQUIRE_THAT(res.real(), WithinAbs(1.0, 1e-5));
		REQUIRE_THAT(res.imag(), WithinAbs(0.0, 1e-5));
	}

	/*********************************************************************/
	/*****           Argument Principle / Count Zeros                *****/
	/*********************************************************************/

	TEST_CASE("ArgumentPrinciple_z_squared_plus_1", "[complex][argument]")
	{
		// f(z) = z^2 + 1, roots at ±i
		// A circle of radius 2 around origin contains both roots → N = 2
		ComplexFunctionFromStdFunc f([](Complex z) -> Complex { return z * z + REAL(1.0); });
		ComplexAnalysis::CircleContour contour(Complex(0.0, 0.0), 2.0);

		int n = ComplexAnalysis::ArgumentPrinciple(f, contour);
		REQUIRE(n == 2);
	}

	TEST_CASE("CountZeros_z_cubed", "[complex][argument]")
	{
		// f(z) = z^3 has a triple root at 0
		// Circle radius 1 around origin → N = 3
		ComplexFunctionFromStdFunc f([](Complex z) -> Complex { return z * z * z; });
		ComplexAnalysis::CircleContour contour(Complex(0.0, 0.0), 1.0);

		int n = ComplexAnalysis::CountZeros(f, contour);
		REQUIRE(n == 3);
	}

	TEST_CASE("CountZeros_polynomial_partial", "[complex][argument]")
	{
		// f(z) = z^2 + 1, roots at ±i
		// Circle of radius 0.5 centered at i should contain exactly 1 root
		ComplexFunctionFromStdFunc f([](Complex z) -> Complex { return z * z + REAL(1.0); });
		ComplexAnalysis::CircleContour contour(Complex(0.0, 1.0), 0.5);

		int n = ComplexAnalysis::CountZeros(f, contour);
		REQUIRE(n == 1);
	}

	TEST_CASE("CountZeros_no_zeros_inside", "[complex][argument]")
	{
		// f(z) = z^2 + 1, roots at ±i
		// Circle of radius 0.5 centered at (5, 0) contains no roots
		ComplexFunctionFromStdFunc f([](Complex z) -> Complex { return z * z + REAL(1.0); });
		ComplexAnalysis::CircleContour contour(Complex(5.0, 0.0), 0.5);

		int n = ComplexAnalysis::CountZeros(f, contour);
		REQUIRE(n == 0);
	}

} // namespace MML::Tests::Core::ComplexAnalysisTests

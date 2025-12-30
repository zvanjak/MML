#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#include "core/Integration/PathIntegration.h"
#include "core/Curves.h"
#include "base/Geometry3D.h"
#endif

using namespace MML;
using namespace MML::Testing;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace MML::Tests::Algorithms::PathIntegrationTests
{
	/*****************************************************************************
	 *                     ARC LENGTH TESTS
	 * Test ParametricCurveLength for curves with known analytical formulas
	 *****************************************************************************/

	TEST_CASE("ArcLength::circle_unit_radius", "[path_integration][arc_length]")
	{
			TEST_PRECISION_INFO();
		// Circle: r(t) = (cos t, sin t, 0), t ∈ [0, 2π]
		// |r'(t)| = 1, so L = 2π
		ParametricCurve<3> circle([](Real t) {
			return VectorN<Real, 3>{std::cos(t), std::sin(t), REAL(0.0)};
		});

		Real length = PathIntegration::ParametricCurveLength(circle, REAL(0.0), REAL(2.0) * Constants::PI);
		Real expected = REAL(2.0) * Constants::PI;

		REQUIRE_THAT(length, WithinRel(expected, REAL(1e-4)));
	}

	TEST_CASE("ArcLength::circle_radius_R", "[path_integration][arc_length]")
	{
			TEST_PRECISION_INFO();
		// Circle radius R=3: r(t) = (3cos t, 3sin t, 0)
		// |r'(t)| = 3, so L = 6π
		Real R = REAL(3.0);
		ParametricCurveFromStdFunc<3> circle([R](Real t) {
			return VectorN<Real, 3>{R * std::cos(t), R * std::sin(t), REAL(0.0)};
		});

		Real length = PathIntegration::ParametricCurveLength(circle, REAL(0.0), REAL(2.0) * Constants::PI);
		Real expected = REAL(2.0) * Constants::PI * R;

		REQUIRE_THAT(length, WithinRel(expected, REAL(1e-4)));
	}

	TEST_CASE("ArcLength::semicircle", "[path_integration][arc_length]")
	{
			TEST_PRECISION_INFO();
		// Semicircle: L = π
		ParametricCurve<3> semicircle([](Real t) {
			return VectorN<Real, 3>{std::cos(t), std::sin(t), REAL(0.0)};
		});

		Real length = PathIntegration::ParametricCurveLength(semicircle, REAL(0.0), Constants::PI);
		Real expected = Constants::PI;

		REQUIRE_THAT(length, WithinRel(expected, REAL(1e-4)));
	}

	TEST_CASE("ArcLength::helix_standard", "[path_integration][arc_length]")
	{
			TEST_PRECISION_INFO();
		// Helix: r(t) = (cos t, sin t, t), t ∈ [0, 2π]
		// |r'(t)| = √(sin²t + cos²t + 1) = √2
		// L = 2π√2
		ParametricCurve<3> helix([](Real t) {
			return VectorN<Real, 3>{std::cos(t), std::sin(t), t};
		});

		Real length = PathIntegration::ParametricCurveLength(helix, REAL(0.0), REAL(2.0) * Constants::PI);
		Real expected = REAL(2.0) * Constants::PI * std::sqrt(REAL(2.0));

		REQUIRE_THAT(length, WithinRel(expected, REAL(1e-4)));
	}

	TEST_CASE("ArcLength::helix_general", "[path_integration][arc_length]")
	{
			TEST_PRECISION_INFO();
		// General helix: r(t) = (a cos t, a sin t, b t)
		// |r'(t)| = √(a² + b²)
		// L = 2π√(a² + b²)
		Real a = REAL(2.0), b = REAL(3.0);
		ParametricCurveFromStdFunc<3> helix([a, b](Real t) {
			return VectorN<Real, 3>{a * std::cos(t), a * std::sin(t), b * t};
		});

		Real length = PathIntegration::ParametricCurveLength(helix, REAL(0.0), REAL(2.0) * Constants::PI);
		Real expected = REAL(2.0) * Constants::PI * std::sqrt(a * a + b * b);

		REQUIRE_THAT(length, WithinRel(expected, REAL(1e-4)));
	}

	TEST_CASE("ArcLength::straight_line", "[path_integration][arc_length]")
	{
			TEST_PRECISION_INFO();
		// Line from (0,0,0) to (3,4,0): length = 5
		ParametricCurve<3> line([](Real t) {
			return VectorN<Real, 3>{REAL(3.0) * t, REAL(4.0) * t, REAL(0.0)};
		});

		Real length = PathIntegration::ParametricCurveLength(line, REAL(0.0), REAL(1.0));
		Real expected = REAL(5.0);

		REQUIRE_THAT(length, WithinRel(expected, REAL(1e-6)));
	}

	TEST_CASE("ArcLength::parabola", "[path_integration][arc_length]")
	{
			TEST_PRECISION_INFO();
		// Parabola: r(t) = (t, t², 0), t ∈ [0, 1]
		// |r'(t)| = √(1 + 4t²)
		// L = ∫[0,1] √(1 + 4t²) dt = (√5 + asinh(2)/2) / 2 ≈ REAL(1.4789)
		ParametricCurve<3> parabola([](Real t) {
			return VectorN<Real, 3>{t, t * t, REAL(0.0)};
		});

		Real length = PathIntegration::ParametricCurveLength(parabola, REAL(0.0), REAL(1.0));
		// Analytical: (√5 + asinh(2)/2) / 2
		Real expected = (std::sqrt(REAL(5.0)) + std::asinh(REAL(2.0)) / REAL(2.0)) / REAL(2.0);
		// Approx REAL(1.4789)

		REQUIRE_THAT(length, WithinRel(expected, REAL(1e-3)));
	}

	TEST_CASE("ArcLength::2D_circle", "[path_integration][arc_length][2D]")
	{
			TEST_PRECISION_INFO();
		// 2D circle
		ParametricCurve<2> circle2D([](Real t) {
			return VectorN<Real, 2>{std::cos(t), std::sin(t)};
		});

		Real length = PathIntegration::ParametricCurveLength(circle2D, REAL(0.0), REAL(2.0) * Constants::PI);
		Real expected = REAL(2.0) * Constants::PI;

		REQUIRE_THAT(length, WithinRel(expected, REAL(1e-4)));
	}

	/*****************************************************************************
	 *                     CURVE MASS TESTS
	 * Test ParametricCurveMass with various density functions
	 *****************************************************************************/

	TEST_CASE("CurveMass::uniform_density", "[path_integration][curve_mass]")
	{
			TEST_PRECISION_INFO();
		// Uniform density ρ=1 should equal arc length
		ParametricCurve<3> circle([](Real t) {
			return VectorN<Real, 3>{std::cos(t), std::sin(t), REAL(0.0)};
		});

		RealFunction unit_density([](Real) { return REAL(1.0); });

		Real mass = PathIntegration::ParametricCurveMass(circle, unit_density, REAL(0.0), REAL(2.0) * Constants::PI);
		Real expected = REAL(2.0) * Constants::PI;  // Same as arc length

		REQUIRE_THAT(mass, WithinRel(expected, REAL(1e-4)));
	}

	TEST_CASE("CurveMass::constant_density", "[path_integration][curve_mass]")
	{
			TEST_PRECISION_INFO();
		// Constant density ρ=5 on circle: M = 5 * 2π = 10π
		ParametricCurve<3> circle([](Real t) {
			return VectorN<Real, 3>{std::cos(t), std::sin(t), REAL(0.0)};
		});

		RealFunction const_density([](Real) { return REAL(5.0); });

		Real mass = PathIntegration::ParametricCurveMass(circle, const_density, REAL(0.0), REAL(2.0) * Constants::PI);
		Real expected = REAL(5.0) * REAL(2.0) * Constants::PI;

		REQUIRE_THAT(mass, WithinRel(expected, REAL(1e-4)));
	}

	TEST_CASE("CurveMass::helix_exponential_density", "[path_integration][curve_mass]")
	{
			TEST_PRECISION_INFO();
		// Helix with exponential density: ρ(t) = e^(-t)
		// |r'(t)| = √2 for helix r(t) = (cos t, sin t, t)
		// M = √2 ∫[0,2π] e^(-t) dt = √2 * (1 - e^(-2π))
		ParametricCurve<3> helix([](Real t) {
			return VectorN<Real, 3>{std::cos(t), std::sin(t), t};
		});

		RealFunction exp_density([](Real t) { return std::exp(-t); });

		Real mass = PathIntegration::ParametricCurveMass(helix, exp_density, REAL(0.0), REAL(2.0) * Constants::PI);
		Real expected = std::sqrt(REAL(2.0)) * (REAL(1.0) - std::exp(-REAL(2.0) * Constants::PI));

		REQUIRE_THAT(mass, WithinRel(expected, REAL(1e-3)));
	}

	TEST_CASE("CurveMass::line_linear_density", "[path_integration][curve_mass]")
	{
			TEST_PRECISION_INFO();
		// Line from (0,0,0) to (1,0,0) with ρ(t) = t
		// |r'(t)| = 1
		// M = ∫[0,1] t dt = 1/2
		ParametricCurve<3> line([](Real t) {
			return VectorN<Real, 3>{t, REAL(0.0), REAL(0.0)};
		});

		RealFunction linear_density([](Real t) { return t; });

		Real mass = PathIntegration::ParametricCurveMass(line, linear_density, REAL(0.0), REAL(1.0));
		Real expected = REAL(0.5);

		REQUIRE_THAT(mass, WithinAbs(expected, REAL(1e-6)));
	}

	TEST_CASE("CurveMass::semicircle_sinusoidal", "[path_integration][curve_mass]")
	{
			TEST_PRECISION_INFO();
		// Semicircle with ρ(t) = sin(t) for t ∈ [0, π]
		// |r'(t)| = 1
		// M = ∫[0,π] sin(t) dt = 2
		ParametricCurve<3> semicircle([](Real t) {
			return VectorN<Real, 3>{std::cos(t), std::sin(t), REAL(0.0)};
		});

		RealFunction sin_density([](Real t) { return std::sin(t); });

		Real mass = PathIntegration::ParametricCurveMass(semicircle, sin_density, REAL(0.0), Constants::PI);
		Real expected = REAL(2.0);

		REQUIRE_THAT(mass, WithinRel(expected, REAL(1e-4)));
	}

	/*****************************************************************************
	 *                 SCALAR LINE INTEGRAL TESTS
	 * Test LineIntegral for scalar fields: ∫_C f ds
	 *****************************************************************************/

	TEST_CASE("ScalarLineIntegral::constant_field", "[path_integration][scalar_integral]")
	{
			TEST_PRECISION_INFO();
		// Constant field f=1 gives arc length
		ScalarFunction<3> unit_field([](const VectorN<Real, 3>&) { return REAL(1.0); });

		ParametricCurve<3> circle([](Real t) {
			return VectorN<Real, 3>{std::cos(t), std::sin(t), REAL(0.0)};
		});

		Real integral = PathIntegration::LineIntegral(unit_field, circle, REAL(0.0), REAL(2.0) * Constants::PI, 1e-5);
		Real expected = REAL(2.0) * Constants::PI;

		REQUIRE_THAT(integral, WithinRel(expected, REAL(1e-3)));
	}

	TEST_CASE("ScalarLineIntegral::linear_field_quarter_circle", "[path_integration][scalar_integral]")
	{
			TEST_PRECISION_INFO();
		// f(x,y,z) = x + y on quarter circle (radius 3)
		// r(t) = (3cos t, 3sin t, 0), t ∈ [0, π/2]
		// f(r(t)) = 3cos t + 3sin t
		// |r'(t)| = 3
		// ∫ = ∫[0,π/2] (3cos t + 3sin t) * 3 dt = 9[sin t - cos t]₀^(π/2) = 18
		Real R = REAL(3.0);
		ScalarFunction<3> linear_field([](const VectorN<Real, 3>& r) {
			return r[0] + r[1];  // x + y
		});

		ParametricCurveFromStdFunc<3> quarter_circle([R](Real t) {
			return VectorN<Real, 3>{R * std::cos(t), R * std::sin(t), REAL(0.0)};
		});

		Real integral = PathIntegration::LineIntegral(linear_field, quarter_circle, REAL(0.0), Constants::PI / REAL(2.0), 1e-5);
		Real expected = REAL(18.0);

		REQUIRE_THAT(integral, WithinRel(expected, REAL(1e-3)));
	}

	TEST_CASE("ScalarLineIntegral::radial_field_circle", "[path_integration][scalar_integral]")
	{
			TEST_PRECISION_INFO();
		// f(x,y,z) = √(x² + y²) = |r| on unit circle
		// On unit circle, f = 1 always
		// So integral = arc length = 2π
		ScalarFunction<3> radial_field([](const VectorN<Real, 3>& r) {
			return std::sqrt(r[0] * r[0] + r[1] * r[1]);
		});

		ParametricCurve<3> circle([](Real t) {
			return VectorN<Real, 3>{std::cos(t), std::sin(t), REAL(0.0)};
		});

		Real integral = PathIntegration::LineIntegral(radial_field, circle, REAL(0.0), REAL(2.0) * Constants::PI, 1e-5);
		Real expected = REAL(2.0) * Constants::PI;

		REQUIRE_THAT(integral, WithinRel(expected, REAL(1e-3)));
	}

	TEST_CASE("ScalarLineIntegral::height_on_helix", "[path_integration][scalar_integral]")
	{
			TEST_PRECISION_INFO();
		// f(x,y,z) = z on helix r(t) = (cos t, sin t, t)
		// f(r(t)) = t
		// |r'(t)| = √2
		// ∫[0,2π] t * √2 dt = √2 * [t²/2]₀^(2π) = √2 * 2π² = 2√2 π²
		ScalarFunction<3> height_field([](const VectorN<Real, 3>& r) {
			return r[2];  // z coordinate
		});

		ParametricCurve<3> helix([](Real t) {
			return VectorN<Real, 3>{std::cos(t), std::sin(t), t};
		});

		Real integral = PathIntegration::LineIntegral(height_field, helix, REAL(0.0), REAL(2.0) * Constants::PI, 1e-5);
		Real expected = std::sqrt(REAL(2.0)) * REAL(2.0) * Constants::PI * Constants::PI;

		REQUIRE_THAT(integral, WithinRel(expected, REAL(1e-3)));
	}

	/*****************************************************************************
	 *                 VECTOR LINE INTEGRAL (WORK) TESTS
	 * Test LineIntegral for vector fields: ∫_C F·dr
	 *****************************************************************************/

	TEST_CASE("VectorLineIntegral::constant_tangent_field", "[path_integration][vector_integral]")
	{
			TEST_PRECISION_INFO();
		// F = (1, 0, 0) along x-axis from 0 to L
		// r(t) = (t, 0, 0), t ∈ [0, L]
		// r'(t) = (1, 0, 0)
		// F·r' = 1
		// ∫[0,L] 1 dt = L
		Real L = REAL(5.0);
		VectorFunction<3> const_field([](const VectorN<Real, 3>&) {
			return VectorN<Real, 3>{REAL(1.0), REAL(0.0), REAL(0.0)};
		});

		ParametricCurve<3> x_line([](Real t) {
			return VectorN<Real, 3>{t, REAL(0.0), REAL(0.0)};
		});

		Real work = PathIntegration::LineIntegral(const_field, x_line, REAL(0.0), L, 1e-5);
		Real expected = L;

		REQUIRE_THAT(work, WithinRel(expected, REAL(1e-6)));
	}

	TEST_CASE("VectorLineIntegral::rotation_field_circle", "[path_integration][vector_integral]")
	{
			TEST_PRECISION_INFO();
		// Rotation field: F = (-y, x, 0)
		// Circle: r(t) = (R cos t, R sin t, 0)
		// r'(t) = (-R sin t, R cos t, 0)
		// F(r(t)) = (-R sin t, R cos t, 0)
		// F·r' = R²sin²t + R²cos²t = R²
		// ∫[0,2π] R² dt = 2πR²
		Real R = REAL(2.0);
		VectorFunction<3> rotation_field([](const VectorN<Real, 3>& r) {
			return VectorN<Real, 3>{-r[1], r[0], REAL(0.0)};
		});

		ParametricCurveFromStdFunc<3> circle([R](Real t) {
			return VectorN<Real, 3>{R * std::cos(t), R * std::sin(t), REAL(0.0)};
		});

		Real circulation = PathIntegration::LineIntegral(rotation_field, circle, REAL(0.0), REAL(2.0) * Constants::PI, 1e-5);
		Real expected = REAL(2.0) * Constants::PI * R * R;

		REQUIRE_THAT(circulation, WithinRel(expected, REAL(1e-3)));
	}

	TEST_CASE("VectorLineIntegral::conservative_field_path_independence", "[path_integration][vector_integral][conservative]")
	{
			TEST_PRECISION_INFO();
		// Conservative field: F = ∇φ where φ = x² + y² + z²
		// F = (2x, 2y, 2z)
		// Work = φ(end) - φ(start)
		VectorFunction<3> grad_field([](const VectorN<Real, 3>& r) {
			return VectorN<Real, 3>{REAL(2.0) * r[0], REAL(2.0) * r[1], REAL(2.0) * r[2]};
		});

		// Path 1: Straight line from (1,0,0) to (0,1,1)
		ParametricCurve<3> path1([](Real t) {
			return VectorN<Real, 3>{REAL(1.0) - t, t, t};
		});

		// Path 2: Two segments
		// (1,0,0) -> (0,0,0) -> (0,1,1)
		ParametricCurve<3> path2a([](Real t) {
			return VectorN<Real, 3>{REAL(1.0) - t, REAL(0.0), REAL(0.0)};
		});
		ParametricCurve<3> path2b([](Real t) {
			return VectorN<Real, 3>{REAL(0.0), t, t};
		});

		Real work1 = PathIntegration::LineIntegral(grad_field, path1, REAL(0.0), REAL(1.0), 1e-5);
		Real work2a = PathIntegration::LineIntegral(grad_field, path2a, REAL(0.0), REAL(1.0), 1e-5);
		Real work2b = PathIntegration::LineIntegral(grad_field, path2b, REAL(0.0), REAL(1.0), 1e-5);
		Real work2 = work2a + work2b;

		// Both should equal φ(end) - φ(start) = (0+1+1) - (1+0+0) = 1
		Real expected = REAL(1.0);

		REQUIRE_THAT(work1, WithinRel(expected, REAL(1e-3)));
		REQUIRE_THAT(work2, WithinRel(expected, REAL(1e-3)));
		REQUIRE_THAT(work1, WithinRel(work2, REAL(1e-3)));  // Path independence
	}

	TEST_CASE("VectorLineIntegral::closed_loop_conservative", "[path_integration][vector_integral][conservative]")
	{
			TEST_PRECISION_INFO();
		// Conservative field: closed loop integral = 0
		// F = ∇φ = (2x, 2y, 0)
		VectorFunction<3> grad_field([](const VectorN<Real, 3>& r) {
			return VectorN<Real, 3>{REAL(2.0) * r[0], REAL(2.0) * r[1], REAL(0.0)};
		});

		// Closed circle
		ParametricCurve<3> circle([](Real t) {
			return VectorN<Real, 3>{std::cos(t), std::sin(t), REAL(0.0)};
		});

		Real work = PathIntegration::LineIntegral(grad_field, circle, REAL(0.0), REAL(2.0) * Constants::PI, 1e-5);

		REQUIRE_THAT(work, WithinAbs(REAL(0.0), REAL(1e-3)));
	}

	TEST_CASE("VectorLineIntegral::radial_field_line", "[path_integration][vector_integral]")
	{
			TEST_PRECISION_INFO();
		// Radial field: F = r (points outward)
		// Line from (1,0,0) to (2,0,0): r(t) = (1+t, 0, 0), t ∈ [0,1]
		// r'(t) = (1, 0, 0)
		// F(r(t)) = (1+t, 0, 0)
		// F·r' = 1+t
		// ∫[0,1] (1+t) dt = [t + t²/2]₀¹ = REAL(1.5)
		VectorFunction<3> radial_field([](const VectorN<Real, 3>& r) {
			return r;  // F = r
		});

		ParametricCurve<3> line([](Real t) {
			return VectorN<Real, 3>{REAL(1.0) + t, REAL(0.0), REAL(0.0)};
		});

		Real work = PathIntegration::LineIntegral(radial_field, line, REAL(0.0), REAL(1.0), 1e-5);
		Real expected = REAL(1.5);

		REQUIRE_THAT(work, WithinRel(expected, REAL(1e-4)));
	}

	TEST_CASE("VectorLineIntegral::inverse_square_field", "[path_integration][vector_integral][physics]")
	{
			TEST_PRECISION_INFO();
		// Gravitational-like: F = -k * r / |r|³ = -k * r̂ / |r|² (points toward origin)
		// This is conservative with φ = -k/|r|
		// Work = φ(end) - φ(start) = -k/2 - (-k/1) = -5 + 10 = 5? 
		// BUT: ∫F·dr = ∫(-k/r²)dr = k/r evaluated at r=2 minus at r=1 = k/2 - k = -k/2
		// The work is NEGATIVE because we're moving against the (inward) force
		Real k = REAL(10.0);
		VectorFunctionFromStdFunc<3> inv_sq_field([k](const VectorN<Real, 3>& r) {
			Real r_mag = r.NormL2();
			Real factor = -k / (r_mag * r_mag * r_mag);
			return VectorN<Real, 3>{factor * r[0], factor * r[1], factor * r[2]};
		});

		// Line from (1,0,0) to (2,0,0)
		ParametricCurve<3> line([](Real t) {
			return VectorN<Real, 3>{REAL(1.0) + t, REAL(0.0), REAL(0.0)};
		});

		Real work = PathIntegration::LineIntegral(inv_sq_field, line, REAL(0.0), REAL(1.0), 1e-5);
		// Work = ∫(-k/r²)dr from r=1 to r=2 = k[1/r]₁² = k(1/2 - 1) = -k/2 = -5
		Real expected = -k / REAL(2.0);  // = -5

		REQUIRE_THAT(work, WithinRel(expected, REAL(1e-3)));
	}

	/*****************************************************************************
	 *                     EDGE CASES AND SPECIAL TESTS
	 *****************************************************************************/

	TEST_CASE("PathIntegration::zero_length_curve", "[path_integration][edge_cases]")
	{
			TEST_PRECISION_INFO();
		// Zero interval should give zero
		ParametricCurve<3> circle([](Real t) {
			return VectorN<Real, 3>{std::cos(t), std::sin(t), REAL(0.0)};
		});

		Real length = PathIntegration::ParametricCurveLength(circle, REAL(0.0), REAL(0.0));

		REQUIRE_THAT(length, WithinAbs(REAL(0.0), REAL(1e-10)));
	}

	TEST_CASE("PathIntegration::reversed_orientation", "[path_integration][orientation]")
	{
			TEST_PRECISION_INFO();
		// Vector integral with reversed path should negate
		VectorFunction<3> field([](const VectorN<Real, 3>& r) {
			return VectorN<Real, 3>{r[0], r[1], REAL(0.0)};
		});

		// Forward: (0,0,0) to (1,1,0)
		ParametricCurve<3> forward([](Real t) {
			return VectorN<Real, 3>{t, t, REAL(0.0)};
		});

		// Backward: (1,1,0) to (0,0,0)
		ParametricCurve<3> backward([](Real t) {
			return VectorN<Real, 3>{REAL(1.0) - t, REAL(1.0) - t, REAL(0.0)};
		});

		Real work_fwd = PathIntegration::LineIntegral(field, forward, REAL(0.0), REAL(1.0), 1e-5);
		Real work_bwd = PathIntegration::LineIntegral(field, backward, REAL(0.0), REAL(1.0), 1e-5);

		REQUIRE_THAT(work_fwd, WithinRel(-work_bwd, REAL(1e-3)));
	}

	TEST_CASE("PathIntegration::scalar_integral_symmetry", "[path_integration][symmetry]")
	{
			TEST_PRECISION_INFO();
		// Scalar integral is path-direction independent
		ScalarFunction<3> field([](const VectorN<Real, 3>& r) {
			return r[0] * r[0] + r[1] * r[1];
		});

		ParametricCurve<3> forward([](Real t) {
			return VectorN<Real, 3>{t, t, REAL(0.0)};
		});

		ParametricCurve<3> backward([](Real t) {
			return VectorN<Real, 3>{REAL(1.0) - t, REAL(1.0) - t, REAL(0.0)};
		});

		Real int_fwd = PathIntegration::LineIntegral(field, forward, REAL(0.0), REAL(1.0), 1e-5);
		Real int_bwd = PathIntegration::LineIntegral(field, backward, REAL(0.0), REAL(1.0), 1e-5);

		REQUIRE_THAT(int_fwd, WithinRel(int_bwd, REAL(1e-3)));  // Same value
	}

	/*****************************************************************************
	 *                     PHYSICS APPLICATION TESTS
	 *****************************************************************************/

	TEST_CASE("Physics::gravitational_work", "[path_integration][physics]")
	{
			TEST_PRECISION_INFO();
		// Work done moving in uniform gravitational field
		// F = (0, 0, -mg), moving from z=0 to z=h
		// Work = -mgh (negative = work done against gravity)
		Real m = REAL(1.0), g = REAL(9.81), h = REAL(10.0);

		VectorFunctionFromStdFunc<3> gravity([m, g](const VectorN<Real, 3>&) {
			return VectorN<Real, 3>{REAL(0.0), REAL(0.0), -m * g};
		});

		// Vertical line from (0,0,0) to (0,0,h)
		ParametricCurveFromStdFunc<3> vertical([h](Real t) {
			return VectorN<Real, 3>{REAL(0.0), REAL(0.0), h * t};
		});

		Real work = PathIntegration::LineIntegral(gravity, vertical, REAL(0.0), REAL(1.0), 1e-6);
		Real expected = -m * g * h;

		REQUIRE_THAT(work, WithinRel(expected, REAL(1e-5)));
	}

	TEST_CASE("Physics::electric_field_work", "[path_integration][physics]")
	{
			TEST_PRECISION_INFO();
		// Uniform electric field E = (E₀, 0, 0)
		// Work on charge q moving from x=0 to x=d: W = qE₀d
		Real E0 = REAL(100.0), q = 1.6e-19, d = REAL(0.01);

		VectorFunctionFromStdFunc<3> E_field([E0, q](const VectorN<Real, 3>&) {
			return VectorN<Real, 3>{q * E0, REAL(0.0), REAL(0.0)};  // Force = qE
		});

		ParametricCurveFromStdFunc<3> path([d](Real t) {
			return VectorN<Real, 3>{d * t, REAL(0.0), REAL(0.0)};
		});

		Real work = PathIntegration::LineIntegral(E_field, path, REAL(0.0), REAL(1.0), 1e-6);
		Real expected = q * E0 * d;

		REQUIRE_THAT(work, WithinRel(expected, REAL(1e-5)));
	}

}  // namespace

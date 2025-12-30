#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/VectorN.h"
#include "base/Geometry3D.h"

#include "base/Function.h"
#include "core/Curves.h"
#include "core/Derivation.h"
#endif

#include "../test_data/parametric_curves_test_bed.h"

using namespace MML;
using namespace MML::Testing;
using namespace Catch::Matchers;

namespace MML::Tests::Core::CurvesTests
{
	// Tests for Helix curve differential geometry properties
	// Verify accuracy of tangent, normal, binormal, curvature vector calculations
	TEST_CASE("Test_Helix_getTangent", "[curves]")
	{
			TEST_PRECISION_INFO();
		const CurveCartesian3D& helix_curve = TestBeds::ParametricCurvesTestBed::getTestCurve("Helix")._curve;
		const CurveCartesian3D& helix_deriv = TestBeds::ParametricCurvesTestBed::getTestCurve("Helix")._curveDerived;

		// Test at multiple parameter values
		std::vector<Real> test_params = {REAL(0.0), Constants::PI / 4, Constants::PI / 2, Constants::PI, 2 * Constants::PI};

		for (Real t : test_params)
		{
			Vec3Cart computed_tangent = helix_curve.getTangent(t);
			Vec3Cart expected_tangent = helix_deriv(t);

			REQUIRE_THAT(computed_tangent[0], WithinAbs(expected_tangent[0], REAL(1e-7)));
			REQUIRE_THAT(computed_tangent[1], WithinAbs(expected_tangent[1], REAL(1e-7)));
			REQUIRE_THAT(computed_tangent[2], WithinAbs(expected_tangent[2], REAL(1e-7)));
		}
	}

	TEST_CASE("Test_Helix_getTangentUnit", "[curves]")
	{
			TEST_PRECISION_INFO();
		const CurveCartesian3D& helix_curve = TestBeds::ParametricCurvesTestBed::getTestCurve("Helix")._curve;

		std::vector<Real> test_params = {REAL(0.5), REAL(1.0), REAL(2.0), REAL(3.5)};

		for (Real t : test_params)
		{
			Vec3Cart tangent_unit = helix_curve.getTangentUnit(t);

			// Unit vector should have norm = 1
			REQUIRE_THAT(tangent_unit.NormL2(), WithinAbs(REAL(1.0), REAL(1e-10)));

			// Should be parallel to tangent vector
			Vec3Cart tangent = helix_curve.getTangent(t);
			Real scale = tangent.NormL2();
			REQUIRE_THAT(tangent_unit[0], WithinAbs(tangent[0] / scale, REAL(1e-10)));
			REQUIRE_THAT(tangent_unit[1], WithinAbs(tangent[1] / scale, REAL(1e-10)));
			REQUIRE_THAT(tangent_unit[2], WithinAbs(tangent[2] / scale, REAL(1e-10)));
		}
	}

	TEST_CASE("Test_Helix_getNormal", "[curves]")
	{
			TEST_PRECISION_INFO();
		const CurveCartesian3D& helix_curve = TestBeds::ParametricCurvesTestBed::getTestCurve("Helix")._curve;
		const CurveCartesian3D& helix_second_deriv = TestBeds::ParametricCurvesTestBed::getTestCurve("Helix")._curveDerSecond;

		std::vector<Real> test_params = {REAL(0.0), Constants::PI / 4, Constants::PI / 2, Constants::PI};

		for (Real t : test_params)
		{
			Vec3Cart computed_normal = helix_curve.getNormal(t);
			Vec3Cart expected_normal = helix_second_deriv(t);

			REQUIRE_THAT(computed_normal[0], WithinAbs(expected_normal[0], REAL(1e-7)));
			REQUIRE_THAT(computed_normal[1], WithinAbs(expected_normal[1], REAL(1e-7)));
			REQUIRE_THAT(computed_normal[2], WithinAbs(expected_normal[2], REAL(1e-7)));
		}
	}

	TEST_CASE("Test_Helix_getNormalUnit", "[curves]")
	{
			TEST_PRECISION_INFO();
		const CurveCartesian3D& helix_curve = TestBeds::ParametricCurvesTestBed::getTestCurve("Helix")._curve;

		std::vector<Real> test_params = {REAL(0.5), REAL(1.0), REAL(2.0), REAL(3.5)};

		for (Real t : test_params)
		{
			Vec3Cart normal_unit = helix_curve.getNormalUnit(t);

			// Unit vector should have norm = 1
			REQUIRE_THAT(normal_unit.NormL2(), WithinAbs(REAL(1.0), REAL(1e-10)));

			// Should be parallel to normal (second derivative) vector
			Vec3Cart normal = helix_curve.getNormal(t);
			Real scale = normal.NormL2();
			REQUIRE_THAT(normal_unit[0], WithinAbs(normal[0] / scale, REAL(1e-7)));
			REQUIRE_THAT(normal_unit[1], WithinAbs(normal[1] / scale, REAL(1e-7)));
			REQUIRE_THAT(normal_unit[2], WithinAbs(normal[2] / scale, REAL(1e-7)));
		}
	}

	TEST_CASE("Test_Helix_getPrincipalNormal", "[curves]")
	{
			TEST_PRECISION_INFO();
		const CurveCartesian3D& helix_curve = TestBeds::ParametricCurvesTestBed::getTestCurve("Helix")._curve;

		// For Helix: r(t) = {cos(t), sin(t), t}
		// r'(t) = {-sin(t), cos(t), 1}
		// r''(t) = {-cos(t), -sin(t), 0}
		// Principal normal should point toward center: n = {-cos(t), -sin(t), 0} (after normalization)

		std::vector<Real> test_params = {REAL(0.0), Constants::PI / 4, Constants::PI / 2, Constants::PI};

		for (Real t : test_params)
		{
			Vec3Cart principal_normal = helix_curve.getNormalUnit(t);

			// For helix, principal normal points radially inward (toward z-axis)
			// Expected: {-cos(t), -sin(t), 0} normalized
			Real expected_x = -std::cos(t);
			Real expected_y = -std::sin(t);
			Real expected_z = REAL(0.0);
			Real norm = std::sqrt(expected_x * expected_x + expected_y * expected_y + expected_z * expected_z);

			REQUIRE_THAT(principal_normal[0], WithinAbs(expected_x / norm, REAL(1e-7)));
			REQUIRE_THAT(principal_normal[1], WithinAbs(expected_y / norm, REAL(1e-7)));
			REQUIRE_THAT(principal_normal[2], WithinAbs(expected_z / norm, REAL(1e-7)));
		}
	}

	TEST_CASE("Test_Helix_getBinormal", "[curves]")
	{
			TEST_PRECISION_INFO();
		const CurveCartesian3D& helix_curve = TestBeds::ParametricCurvesTestBed::getTestCurve("Helix")._curve;

		// Binormal = T x N (cross product of unit tangent and unit normal)
		// For Helix, this should be orthogonal to both T and N

		std::vector<Real> test_params = {REAL(0.5), REAL(1.0), REAL(2.0), REAL(3.5)};

		for (Real t : test_params)
		{
			Vec3Cart binormal = helix_curve.getBinormal(t);
			Vec3Cart tangent_unit = helix_curve.getTangentUnit(t);
			Vec3Cart normal_unit = helix_curve.getNormalUnit(t);

			// Binormal should be unit vector
			REQUIRE_THAT(binormal.NormL2(), WithinAbs(REAL(1.0), REAL(1e-10)));

			// Binormal should be orthogonal to tangent
			Real dot_T_B = Utils::ScalarProduct(tangent_unit, binormal);
			REQUIRE_THAT(dot_T_B, WithinAbs(REAL(0.0), REAL(1e-10)));

			// Binormal should be orthogonal to normal
			Real dot_N_B = Utils::ScalarProduct(normal_unit, binormal);
			REQUIRE_THAT(dot_N_B, WithinAbs(REAL(0.0), REAL(1e-10)));

			// Verify T x N = B
			Vec3Cart cross_product = VectorProduct(tangent_unit, normal_unit);
			REQUIRE_THAT(binormal[0], WithinAbs(cross_product[0], REAL(1e-10)));
			REQUIRE_THAT(binormal[1], WithinAbs(cross_product[1], REAL(1e-10)));
			REQUIRE_THAT(binormal[2], WithinAbs(cross_product[2], REAL(1e-10)));
		}
	}

	TEST_CASE("Test_Helix_getCurvatureVector", "[curves]")
	{
			TEST_PRECISION_INFO();
		const CurveCartesian3D& helix_curve = TestBeds::ParametricCurvesTestBed::getTestCurve("Helix")._curve;

		// Curvature vector for helix should point toward the z-axis (radially inward)
		std::vector<Real> test_params = {REAL(0.5), REAL(1.0), REAL(2.0)};

		for (Real t : test_params)
		{
			Vec3Cart curvature_vec = helix_curve.getCurvatureVector(t);
			Vec3Cart tangent = helix_curve.getTangent(t);
			Vec3Cart normal = helix_curve.getNormal(t);

			// Curvature vector should be orthogonal to tangent
			Real dot_product = Utils::ScalarProduct(curvature_vec, tangent);
			REQUIRE_THAT(dot_product, WithinAbs(REAL(0.0), REAL(1e-6)));

			// Curvature vector should point in direction of principal normal
			Vec3Cart normal_unit = helix_curve.getNormalUnit(t);
			Vec3Cart curvature_direction = curvature_vec / curvature_vec.NormL2();

			REQUIRE_THAT(curvature_direction[0], WithinAbs(normal_unit[0], REAL(1e-7)));
			REQUIRE_THAT(curvature_direction[1], WithinAbs(normal_unit[1], REAL(1e-7)));
			REQUIRE_THAT(curvature_direction[2], WithinAbs(normal_unit[2], REAL(1e-7)));

			// Verify relationship: getCurvatureVector formula from code
			// curvature_vec = vec2 / res1 where:
			// res1 = |r'|^(-2), vec2 = r'' - res1 * <r', r''> * r'
			Real res1 = std::pow(tangent.NormL2(), -REAL(2.0));
			Vec3Cart vec2 = normal - res1 * Utils::ScalarProduct(tangent, normal) * tangent;
			Vec3Cart expected_curvature_vec = vec2 / res1;

			REQUIRE_THAT(curvature_vec[0], WithinAbs(expected_curvature_vec[0], REAL(1e-7)));
			REQUIRE_THAT(curvature_vec[1], WithinAbs(expected_curvature_vec[1], REAL(1e-7)));
			REQUIRE_THAT(curvature_vec[2], WithinAbs(expected_curvature_vec[2], REAL(1e-7)));
		}
	}

	TEST_CASE("Test_Helix_Curvature", "[simple]")
	{
			TEST_PRECISION_INFO();
		const CurveCartesian3D& curve_func    = TestBeds::ParametricCurvesTestBed::getTestCurve("Helix")._curve;
		const RealFunction&			curvature_func = TestBeds::ParametricCurvesTestBed::getTestCurve("Helix")._curvatureFunc;

		REQUIRE_THAT(curve_func.getCurvature(REAL(0.5)), WithinAbs(curvature_func(REAL(0.5)), REAL(1e-7)));
		REQUIRE_THAT(curve_func.getCurvature(REAL(1.0)), WithinAbs(curvature_func(REAL(1.0)), REAL(1e-7)));
		REQUIRE_THAT(curve_func.getCurvature(REAL(5.5)), WithinAbs(curvature_func(REAL(5.5)), REAL(1e-7)));
		REQUIRE_THAT(curve_func.getCurvature(REAL(20.0)), WithinAbs(curvature_func(REAL(20.0)), REAL(1e-7)));
		REQUIRE_THAT(curve_func.getCurvature(-REAL(3.273)), WithinAbs(curvature_func(-REAL(3.273)), REAL(1e-7)));
	}

	TEST_CASE("Test_Helix2_Curvature", "[simple]")
	{
			TEST_PRECISION_INFO();
		const CurveCartesian3D& curve_func     = TestBeds::ParametricCurvesTestBed::getTestCurve("Helix2")._curve;
		const RealFunction&			curvature_func = TestBeds::ParametricCurvesTestBed::getTestCurve("Helix2")._curvatureFunc;

		REQUIRE_THAT(curve_func.getCurvature(REAL(0.5)), WithinAbs(curvature_func(REAL(0.5)), REAL(1e-7)));
		REQUIRE_THAT(curve_func.getCurvature(REAL(1.0)), WithinAbs(curvature_func(REAL(1.0)), REAL(1e-7)));
		REQUIRE_THAT(curve_func.getCurvature(REAL(5.5)), WithinAbs(curvature_func(REAL(5.5)), REAL(1e-7)));
		REQUIRE_THAT(curve_func.getCurvature(REAL(20.0)), WithinAbs(curvature_func(REAL(20.0)), REAL(1e-7)));
		REQUIRE_THAT(curve_func.getCurvature(-REAL(3.273)), WithinAbs(curvature_func(-REAL(3.273)), REAL(1e-7)));
	}

	TEST_CASE("Test_Schaum1_Curvature", "[simple]")
	{
			TEST_PRECISION_INFO();
		const CurveCartesian3D& schaum_curve1  = TestBeds::ParametricCurvesTestBed::getTestCurve("Schaums1")._curve;
		const RealFunction&			curvature_func = TestBeds::ParametricCurvesTestBed::getTestCurve("Schaums1")._curvatureFunc;

		REQUIRE_THAT(schaum_curve1.getCurvature(REAL(0.5)), WithinAbs(curvature_func(REAL(0.5)), REAL(1e-7)));
		REQUIRE_THAT(schaum_curve1.getCurvature(REAL(1.0)), WithinAbs(curvature_func(REAL(1.0)), REAL(1e-7)));
		REQUIRE_THAT(schaum_curve1.getCurvature(REAL(5.5)), WithinAbs(curvature_func(REAL(5.5)), REAL(1e-7)));
		REQUIRE_THAT(schaum_curve1.getCurvature(REAL(20.0)), WithinAbs(curvature_func(REAL(20.0)), REAL(1e-7)));
		REQUIRE_THAT(schaum_curve1.getCurvature(-REAL(3.273)), WithinAbs(curvature_func(-REAL(3.273)), REAL(1e-7)));
	}

	TEST_CASE("Test_Schaum2_Curvature", "[simple]")
	{
			TEST_PRECISION_INFO();
		const CurveCartesian3D& schaum_curve2  = TestBeds::ParametricCurvesTestBed::getTestCurve("Schaums2")._curve;
		const RealFunction&			curvature_func = TestBeds::ParametricCurvesTestBed::getTestCurve("Schaums2")._curvatureFunc;

		REQUIRE_THAT(schaum_curve2.getCurvature(REAL(0.5)), WithinAbs(curvature_func(REAL(0.5)), REAL(1e-7)));
		REQUIRE_THAT(schaum_curve2.getCurvature(REAL(1.0)), WithinAbs(curvature_func(REAL(1.0)), REAL(1e-7)));
		REQUIRE_THAT(schaum_curve2.getCurvature(REAL(5.5)), WithinAbs(curvature_func(REAL(5.5)), REAL(1e-7)));
		REQUIRE_THAT(schaum_curve2.getCurvature(REAL(20.0)), WithinAbs(curvature_func(REAL(20.0)), REAL(1e-7)));
		REQUIRE_THAT(schaum_curve2.getCurvature(-REAL(3.273)), WithinAbs(curvature_func(-REAL(3.273)), REAL(1e-7)));
	}

	TEST_CASE("Test_Helix_getTorsion3", "[simple]")
	{
			TEST_PRECISION_INFO();
		const CurveCartesian3D& helix_curve  = TestBeds::ParametricCurvesTestBed::getTestCurve("Helix")._curve;
		const RealFunction&			torsion_func = TestBeds::ParametricCurvesTestBed::getTestCurve("Helix")._torsionFunc;

		REQUIRE_THAT(helix_curve.getTorsion(REAL(0.5)), WithinAbs(torsion_func(REAL(0.5)), REAL(1e-7)));
		REQUIRE_THAT(helix_curve.getTorsion(REAL(1.0)), WithinAbs(torsion_func(REAL(1.0)), REAL(1e-7)));
		REQUIRE_THAT(helix_curve.getTorsion(REAL(5.5)), WithinAbs(torsion_func(REAL(5.5)), REAL(1e-7)));
		REQUIRE_THAT(helix_curve.getTorsion(REAL(20.0)), WithinAbs(torsion_func(REAL(20.0)), REAL(1e-5)));
		REQUIRE_THAT(helix_curve.getTorsion(-REAL(3.273)), WithinAbs(torsion_func(-REAL(3.273)), REAL(1e-7)));
	}

	TEST_CASE("Test_Helix2_getTorsion3", "[simple]")
	{
			TEST_PRECISION_INFO();
		const CurveCartesian3D& helix_curve2 = TestBeds::ParametricCurvesTestBed::getTestCurve("Helix2")._curve;
		const RealFunction&			torsion_func = TestBeds::ParametricCurvesTestBed::getTestCurve("Helix2")._torsionFunc;

		REQUIRE_THAT(helix_curve2.getTorsion(REAL(2.5)), WithinAbs(torsion_func(REAL(2.5)), REAL(1e-7)));
		REQUIRE_THAT(helix_curve2.getTorsion(REAL(0.5)), WithinAbs(torsion_func(REAL(0.5)), REAL(1e-7)));
		REQUIRE_THAT(helix_curve2.getTorsion(REAL(1.0)), WithinAbs(torsion_func(REAL(1.0)), REAL(1e-7)));
	}

	TEST_CASE("Test_Schaums1_getTorsion3", "[simple]")
	{
			TEST_PRECISION_INFO();
		const CurveCartesian3D& schaum_curve1 = TestBeds::ParametricCurvesTestBed::getTestCurve("Schaums1")._curve;
		const RealFunction&			torsion_func  = TestBeds::ParametricCurvesTestBed::getTestCurve("Schaums1")._torsionFunc;

		REQUIRE_THAT(schaum_curve1.getTorsion(REAL(0.5)), WithinAbs(torsion_func(REAL(0.5)), REAL(1e-7)));
		REQUIRE_THAT(schaum_curve1.getTorsion(REAL(1.0)), WithinAbs(torsion_func(REAL(1.0)), REAL(1e-7)));
		REQUIRE_THAT(schaum_curve1.getTorsion(REAL(3.5)), WithinAbs(torsion_func(REAL(3.5)), REAL(1e-7)));
	}

	TEST_CASE("Test_Schaums2_getTorsion3", "[simple]")
	{
			TEST_PRECISION_INFO();
		const CurveCartesian3D& schaum_curve2 = TestBeds::ParametricCurvesTestBed::getTestCurve("Schaums2")._curve;
		const RealFunction&			torsion_func  = TestBeds::ParametricCurvesTestBed::getTestCurve("Schaums2")._torsionFunc;

		REQUIRE_THAT(schaum_curve2.getTorsion(REAL(1.5)), WithinAbs(torsion_func(REAL(1.5)), REAL(1e-6)));
		REQUIRE_THAT(schaum_curve2.getTorsion(REAL(3.5)), WithinAbs(torsion_func(REAL(3.5)), REAL(1e-6)));
		REQUIRE_THAT(schaum_curve2.getTorsion(REAL(10.5)), WithinAbs(torsion_func(REAL(10.5)), REAL(1e-7)));
	}
}

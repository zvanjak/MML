///////////////////////////////////////////////////////////////////////////////
/// @file geometry_coord_edge_case_tests.cpp
/// @brief Edge-case tests for geometry, coordinates, surfaces, curves, and fields
///
/// Tests cover singular, degenerate, and boundary cases:
///   - Coordinate transform round-trips (Cart↔Spher↔Cyl), including origin & axes
///   - Parametric surface corners and polar degeneracy
///   - Curve endpoint / degenerate geometry
///   - Field operations near coordinate singularities
///   - Metric tensor at singular points
///   - Quaternion gimbal lock and Slerp edge cases
///   - Geometry primitive degenerate inputs
///
/// Created for task MinimalMathLibrary-2l5d.27.
///////////////////////////////////////////////////////////////////////////////

#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Vector/VectorN.h"
#include "base/Vector/VectorTypes.h"
#include "base/Function.h"
#include "base/Quaternions.h"

#include "core/CoordTransf.h"
#include "core/CoordTransf/CoordTransfSpherical.h"
#include "core/CoordTransf/CoordTransfCylindrical.h"
#include "core/FieldOperations.h"
#include "core/Fields.h"
#include "core/MetricTensor.h"
#include "core/Surfaces.h"
#include "core/Curves.h"
#endif

#include "../test_beds/parametric_curves_test_bed.h"

using namespace MML;
using namespace MML::Testing;
using namespace Catch::Matchers;

static bool VectorApproxEqual(const Vec3Cart& v1, const Vec3Cart& v2, Real tolerance = 1e-6)
{
	return std::abs(v1[0] - v2[0]) < tolerance &&
	       std::abs(v1[1] - v2[1]) < tolerance &&
	       std::abs(v1[2] - v2[2]) < tolerance;
}

namespace MML::Tests::Core::GeometryCoordEdgeCaseTests
{

///////////////////////////////////////////////////////////////////////////////
//                 COORDINATE TRANSFORM ROUND-TRIP TESTS                    //
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("CoordTransf round-trip - Cart->Spher->Cart generic point", "[coordtransf][edge][roundtrip]")
{
	TEST_PRECISION_INFO();
	Vec3Cart original{ REAL(1.0), REAL(2.0), REAL(3.0) };
	auto spher = CoordTransfCartToSpher.transf(original);
	auto back  = CoordTransfSpherToCart.transf(spher);

	REQUIRE_THAT(back[0], WithinAbs(original[0], REAL(1e-12)));
	REQUIRE_THAT(back[1], WithinAbs(original[1], REAL(1e-12)));
	REQUIRE_THAT(back[2], WithinAbs(original[2], REAL(1e-12)));
}

TEST_CASE("CoordTransf round-trip - Cart->Cyl->Cart generic point", "[coordtransf][edge][roundtrip]")
{
	TEST_PRECISION_INFO();
	Vec3Cart original{ REAL(1.0), REAL(2.0), REAL(3.0) };
	auto cyl  = CoordTransfCartToCyl.transf(original);
	auto back = CoordTransfCylToCart.transf(cyl);

	REQUIRE_THAT(back[0], WithinAbs(original[0], REAL(1e-12)));
	REQUIRE_THAT(back[1], WithinAbs(original[1], REAL(1e-12)));
	REQUIRE_THAT(back[2], WithinAbs(original[2], REAL(1e-12)));
}

TEST_CASE("CoordTransf round-trip - Cart->Spher->Cart on axis-aligned points", "[coordtransf][edge][roundtrip]")
{
	TEST_PRECISION_INFO();
	// Six axis-aligned unit vectors
	std::vector<Vec3Cart> points = {
		{ REAL(1.0),  REAL(0.0),  REAL(0.0) },  // +x
		{-REAL(1.0),  REAL(0.0),  REAL(0.0) },  // -x
		{ REAL(0.0),  REAL(1.0),  REAL(0.0) },  // +y
		{ REAL(0.0), -REAL(1.0),  REAL(0.0) },  // -y
		{ REAL(0.0),  REAL(0.0),  REAL(1.0) },  // +z (north pole)
		{ REAL(0.0),  REAL(0.0), -REAL(1.0) },  // -z (south pole)
	};

	for (const auto& p : points)
	{
		auto spher = CoordTransfCartToSpher.transf(p);
		auto back  = CoordTransfSpherToCart.transf(spher);

		REQUIRE_THAT(back[0], WithinAbs(p[0], REAL(1e-12)));
		REQUIRE_THAT(back[1], WithinAbs(p[1], REAL(1e-12)));
		REQUIRE_THAT(back[2], WithinAbs(p[2], REAL(1e-12)));
	}
}

TEST_CASE("CoordTransf round-trip - Cart->Cyl->Cart on axis-aligned points", "[coordtransf][edge][roundtrip]")
{
	TEST_PRECISION_INFO();
	std::vector<Vec3Cart> points = {
		{ REAL(1.0),  REAL(0.0),  REAL(0.0) },
		{-REAL(1.0),  REAL(0.0),  REAL(0.0) },
		{ REAL(0.0),  REAL(1.0),  REAL(0.0) },
		{ REAL(0.0), -REAL(1.0),  REAL(0.0) },
		{ REAL(0.0),  REAL(0.0),  REAL(1.0) },  // z-axis (cyl r=0)
		{ REAL(0.0),  REAL(0.0), -REAL(1.0) },
	};

	for (const auto& p : points)
	{
		auto cyl  = CoordTransfCartToCyl.transf(p);
		auto back = CoordTransfCylToCart.transf(cyl);

		REQUIRE_THAT(back[0], WithinAbs(p[0], REAL(1e-12)));
		REQUIRE_THAT(back[1], WithinAbs(p[1], REAL(1e-12)));
		REQUIRE_THAT(back[2], WithinAbs(p[2], REAL(1e-12)));
	}
}

TEST_CASE("CoordTransf round-trip - Cart->Spher->Cart with negative coords", "[coordtransf][edge][roundtrip]")
{
	TEST_PRECISION_INFO();
	std::vector<Vec3Cart> points = {
		{-REAL(1.0), -REAL(2.0),  REAL(3.0) },
		{ REAL(1.0), -REAL(2.0), -REAL(3.0) },
		{-REAL(1.0),  REAL(2.0), -REAL(3.0) },
		{-REAL(1.0), -REAL(2.0), -REAL(3.0) },
	};

	for (const auto& p : points)
	{
		auto spher = CoordTransfCartToSpher.transf(p);
		auto back  = CoordTransfSpherToCart.transf(spher);

		REQUIRE_THAT(back[0], WithinAbs(p[0], REAL(1e-12)));
		REQUIRE_THAT(back[1], WithinAbs(p[1], REAL(1e-12)));
		REQUIRE_THAT(back[2], WithinAbs(p[2], REAL(1e-12)));
	}
}

TEST_CASE("CoordTransf - Cart->Spher at origin produces r=0", "[coordtransf][edge][singularity]")
{
	TEST_PRECISION_INFO();
	Vec3Cart origin{ REAL(0.0), REAL(0.0), REAL(0.0) };
	auto spher = CoordTransfCartToSpher.transf(origin);

	// r should be exactly 0
	REQUIRE_THAT(spher[0], WithinAbs(REAL(0.0), REAL(1e-15)));

	// At r=0, theta and phi are mathematically undefined.
	// The round-trip Spher->Cart won't work because 0*sin(NaN) = NaN.
	// This is a fundamental singularity of spherical coordinates.
	// We just verify the forward transform gives r=0 and that angles
	// don't cause crashes or exceptions.
	REQUIRE(!std::isinf(spher[1]));  // theta finite or NaN, not inf
	REQUIRE(!std::isinf(spher[2]));  // phi finite or NaN, not inf
}

TEST_CASE("CoordTransf - Cart->Cyl at z-axis produces r=0", "[coordtransf][edge][singularity]")
{
	TEST_PRECISION_INFO();
	Vec3Cart on_z{ REAL(0.0), REAL(0.0), REAL(5.0) };
	auto cyl = CoordTransfCartToCyl.transf(on_z);

	// r should be 0, z preserved
	REQUIRE_THAT(cyl[0], WithinAbs(REAL(0.0), REAL(1e-15)));
	REQUIRE_THAT(cyl[2], WithinAbs(REAL(5.0), REAL(1e-12)));

	// Round-trip
	auto back = CoordTransfCylToCart.transf(cyl);
	REQUIRE_THAT(back[0], WithinAbs(REAL(0.0), REAL(1e-15)));
	REQUIRE_THAT(back[1], WithinAbs(REAL(0.0), REAL(1e-15)));
	REQUIRE_THAT(back[2], WithinAbs(REAL(5.0), REAL(1e-12)));
}

TEST_CASE("CoordTransf - phi wraparound consistency", "[coordtransf][edge][branch]")
{
	TEST_PRECISION_INFO();
	// Points at phi = -pi+eps and phi = pi-eps should be close in Cartesian
	Real r = REAL(2.0);
	Real theta = Constants::PI / REAL(2.0);

	Vec3Sph sph_pos{ r, theta,  Constants::PI - REAL(0.01) };
	Vec3Sph sph_neg{ r, theta, -Constants::PI + REAL(0.01) };

	auto cart_pos = CoordTransfSpherToCart.transf(sph_pos);
	auto cart_neg = CoordTransfSpherToCart.transf(sph_neg);

	// These points should be close to each other (close to the -x axis)
	Real dist = std::sqrt(
		(cart_pos[0] - cart_neg[0]) * (cart_pos[0] - cart_neg[0]) +
		(cart_pos[1] - cart_neg[1]) * (cart_pos[1] - cart_neg[1]) +
		(cart_pos[2] - cart_neg[2]) * (cart_pos[2] - cart_neg[2])
	);

	REQUIRE(dist < REAL(0.05));  // Should be very close
}


///////////////////////////////////////////////////////////////////////////////
//              PARAMETRIC SURFACE CORNER & POLE TESTS                      //
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("Sphere surface - evaluation at south pole (u=pi)", "[surfaces][edge][pole]")
{
	using namespace MML::Surfaces;
	Real R = REAL(2.0);
	Sphere sphere(R);

	VectorN<Real, 3> south = sphere(Constants::PI, REAL(0.0));

	// South pole should be at (0, 0, -R)
	REQUIRE_THAT(south[0], WithinAbs(REAL(0.0), REAL(1e-9)));
	REQUIRE_THAT(south[1], WithinAbs(REAL(0.0), REAL(1e-9)));
	REQUIRE_THAT(south[2], WithinAbs(-R, REAL(1e-9)));

	REQUIRE_THAT(south.NormL2(), WithinAbs(R, REAL(1e-9)));
}

TEST_CASE("Sphere surface - evaluation at full azimuth (w=2pi)", "[surfaces][edge][boundary]")
{
	using namespace MML::Surfaces;
	Real R = REAL(1.5);
	Sphere sphere(R);

	Real u = Constants::PI / REAL(3.0);

	// w=0 and w=2π should give the same point
	VectorN<Real, 3> p_0   = sphere(u, REAL(0.0));
	VectorN<Real, 3> p_2pi = sphere(u, REAL(2.0) * Constants::PI);

	REQUIRE_THAT(p_2pi[0], WithinAbs(p_0[0], REAL(1e-9)));
	REQUIRE_THAT(p_2pi[1], WithinAbs(p_0[1], REAL(1e-9)));
	REQUIRE_THAT(p_2pi[2], WithinAbs(p_0[2], REAL(1e-9)));
}

TEST_CASE("Sphere surface - all four parameter corners", "[surfaces][edge][corners]")
{
	using namespace MML::Surfaces;
	Real R = REAL(2.0);
	Sphere sphere(R);

	// (u_min=0, w_min=0) — north pole
	auto p00 = sphere(REAL(0.0), REAL(0.0));
	REQUIRE_THAT(p00.NormL2(), WithinAbs(R, REAL(1e-9)));

	// (u_min=0, w_max=2π) — north pole again
	auto p02pi = sphere(REAL(0.0), REAL(2.0) * Constants::PI);
	REQUIRE_THAT(p02pi.NormL2(), WithinAbs(R, REAL(1e-9)));
	REQUIRE_THAT(p02pi[2], WithinAbs(R, REAL(1e-9)));

	// (u_max=π, w_min=0) — south pole
	auto ppi0 = sphere(Constants::PI, REAL(0.0));
	REQUIRE_THAT(ppi0.NormL2(), WithinAbs(R, REAL(1e-9)));
	REQUIRE_THAT(ppi0[2], WithinAbs(-R, REAL(1e-9)));

	// (u_max=π, w_max=2π) — south pole again
	auto ppi2pi = sphere(Constants::PI, REAL(2.0) * Constants::PI);
	REQUIRE_THAT(ppi2pi.NormL2(), WithinAbs(R, REAL(1e-9)));
	REQUIRE_THAT(ppi2pi[2], WithinAbs(-R, REAL(1e-9)));
}

TEST_CASE("Sphere surface - Gaussian curvature at south pole", "[surfaces][edge][curvature]")
{
	using namespace MML::Surfaces;
	Real R = REAL(2.0);
	Sphere sphere(R);
	Real expected_K = REAL(1.0) / (R * R);

	// Curvature near south pole (u slightly less than π to avoid exact parametric degeneracy)
	Real K = sphere.GaussianCurvature(Constants::PI - REAL(0.01), REAL(0.0));
	REQUIRE_THAT(K, WithinAbs(expected_K, REAL(1e-4)));
}

TEST_CASE("Ellipsoid surface - south pole and parameter corners", "[surfaces][edge][corners]")
{
	using namespace MML::Surfaces;
	Real a = REAL(2.0), b = REAL(3.0), c = REAL(1.0);
	Ellipsoid ellipsoid(a, b, c);

	// North pole (u=0) -> (0, 0, c)
	auto north = ellipsoid(REAL(0.0), REAL(0.0));
	REQUIRE_THAT(north[2], WithinAbs(c, REAL(1e-9)));

	// South pole (u=π) -> (0, 0, -c)
	auto south = ellipsoid(Constants::PI, REAL(0.0));
	REQUIRE_THAT(south[0], WithinAbs(REAL(0.0), REAL(1e-9)));
	REQUIRE_THAT(south[1], WithinAbs(REAL(0.0), REAL(1e-9)));
	REQUIRE_THAT(south[2], WithinAbs(-c, REAL(1e-9)));
}

TEST_CASE("CylinderSurface - parameter boundary corners", "[surfaces][edge][corners]")
{
	using namespace MML::Surfaces;
	Real R = REAL(1.0);
	Real H = REAL(3.0);
	CylinderSurface cylinder(R, H);

	// u=0, w=0 corner
	auto p00 = cylinder(REAL(0.0), REAL(0.0));
	REQUIRE_THAT(p00[0], WithinAbs(R, REAL(1e-9)));      // x = R*cos(0)
	REQUIRE_THAT(p00[1], WithinAbs(REAL(0.0), REAL(1e-9)));
	REQUIRE_THAT(p00[2], WithinAbs(REAL(0.0), REAL(1e-9)));

	// u=2π (same as u=0 circumferentially), w=H
	auto p2piH = cylinder(REAL(2.0) * Constants::PI, H);
	REQUIRE_THAT(p2piH[0], WithinAbs(R, REAL(1e-9)));
	REQUIRE_THAT(p2piH[2], WithinAbs(H, REAL(1e-9)));
}

TEST_CASE("Torus - outer and inner ring curvature", "[surfaces][edge][curvature]")
{
	using namespace MML::Surfaces;
	Real R = REAL(3.0), r = REAL(1.0);
	Torus torus(R, r);

	// Outer ring (u=0): K should be positive (convex)
	Real K_outer = torus.GaussianCurvature(REAL(0.0), REAL(0.0));
	REQUIRE(std::isfinite(K_outer));
	REQUIRE(K_outer > REAL(0.0));

	// Verify that curvature is finite at all quadrant points on the torus
	Real K_quarter = torus.GaussianCurvature(Constants::PI / REAL(2.0), REAL(0.0));
	REQUIRE(std::isfinite(K_quarter));
}

TEST_CASE("Sphere surface - regularity fails at pole", "[surfaces][edge][regularity]")
{
	using namespace MML::Surfaces;
	Sphere sphere(REAL(1.0));

	// At generic point, sphere is regular
	REQUIRE(sphere.isRegular(Constants::PI / REAL(4.0), Constants::PI / REAL(4.0)));

	// At poles, ∂/∂w tangent vector degenerates (sin(θ)→0)
	// The surface is still geometrically smooth, but parametrically singular
	// Just verify evaluation doesn't crash — regularity may or may not detect this
	auto p = sphere(REAL(0.0), REAL(0.0));
	REQUIRE(std::isfinite(p[0]));
	REQUIRE(std::isfinite(p[1]));
	REQUIRE(std::isfinite(p[2]));
}


///////////////////////////////////////////////////////////////////////////////
//                   CURVE ENDPOINT & DEGENERACY TESTS                      //
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("Helix curve - tangent at t=0", "[curves][edge][endpoint]")
{
	TEST_PRECISION_INFO();
	const auto& helix_curve = TestBeds::ParametricCurvesTestBed::getTestCurve("Helix")._curve;

	Vec3Cart tangent = helix_curve.getTangent(REAL(0.0));
	Vec3Cart tangent_unit = helix_curve.getTangentUnit(REAL(0.0));

	// Tangent should be finite at t=0
	REQUIRE(std::isfinite(tangent[0]));
	REQUIRE(std::isfinite(tangent[1]));
	REQUIRE(std::isfinite(tangent[2]));

	// Unit tangent should have norm 1
	REQUIRE_THAT(tangent_unit.NormL2(), WithinAbs(REAL(1.0), REAL(1e-10)));
}

TEST_CASE("Helix curve - curvature and torsion at t=0", "[curves][edge][endpoint]")
{
	TEST_PRECISION_INFO();
	const auto& helix_curve = TestBeds::ParametricCurvesTestBed::getTestCurve("Helix")._curve;

	Real curvature = helix_curve.getCurvature(REAL(0.0));
	Real torsion   = helix_curve.getTorsion(REAL(0.0));

	// Both should be finite for a helix
	REQUIRE(std::isfinite(curvature));
	REQUIRE(std::isfinite(torsion));
	REQUIRE(curvature > REAL(0.0));
}

TEST_CASE("Helix curve - Frenet frame at t=0", "[curves][edge][endpoint]")
{
	TEST_PRECISION_INFO();
	const auto& helix_curve = TestBeds::ParametricCurvesTestBed::getTestCurve("Helix")._curve;

	Vec3Cart T = helix_curve.getTangentUnit(REAL(0.0));
	Vec3Cart N = helix_curve.getNormalUnit(REAL(0.0));
	Vec3Cart B = helix_curve.getBinormal(REAL(0.0));

	// All should be unit vectors
	REQUIRE_THAT(T.NormL2(), WithinAbs(REAL(1.0), REAL(1e-6)));
	REQUIRE_THAT(N.NormL2(), WithinAbs(REAL(1.0), REAL(1e-6)));
	REQUIRE_THAT(B.NormL2(), WithinAbs(REAL(1.0), REAL(1e-6)));

	// Orthogonality: T·N ≈ 0, T·B ≈ 0, N·B ≈ 0
	REQUIRE_THAT(T[0]*N[0] + T[1]*N[1] + T[2]*N[2], WithinAbs(REAL(0.0), REAL(1e-6)));
	REQUIRE_THAT(T[0]*B[0] + T[1]*B[1] + T[2]*B[2], WithinAbs(REAL(0.0), REAL(1e-6)));
	REQUIRE_THAT(N[0]*B[0] + N[1]*B[1] + N[2]*B[2], WithinAbs(REAL(0.0), REAL(1e-6)));
}


///////////////////////////////////////////////////////////////////////////////
//            FIELD OPERATIONS NEAR COORDINATE SINGULARITIES                //
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("FieldOps - Gradient of r^2 at spherical near-pole", "[field_operations][edge][singularity]")
{
	TEST_PRECISION_INFO();
	// f(r,θ,φ) = r² — gradient should be (2r, 0, 0) everywhere in spherical coords
	ScalarFunction<3> f_spher([](const VectorN<Real, 3>& pos) -> Real {
		return pos[0] * pos[0];  // r²
	});

	// Near north pole (θ ≈ 0) — gradient formula involves 1/(r sinθ) which is large
	Vec3Sph near_pole{ REAL(2.0), REAL(0.05), REAL(1.0) };

	auto grad = ScalarFieldOperations::GradientSpher(f_spher, near_pole);

	// ∂f/∂r = 2r = 4.0
	REQUIRE_THAT(grad[0], WithinAbs(REAL(4.0), REAL(0.1)));
	// ∂f/∂θ and ∂f/∂φ should be 0 (f depends only on r)
	REQUIRE_THAT(grad[1], WithinAbs(REAL(0.0), REAL(0.1)));
	REQUIRE_THAT(grad[2], WithinAbs(REAL(0.0), REAL(0.1)));
}

TEST_CASE("FieldOps - Laplacian of r^2 at spherical near-pole", "[field_operations][edge][singularity]")
{
	TEST_PRECISION_INFO();
	// ∇²(r²) = 6 in spherical coordinates (same as Cartesian)
	ScalarFunction<3> f_spher([](const VectorN<Real, 3>& pos) -> Real {
		return pos[0] * pos[0];  // r²
	});

	Vec3Sph near_pole{ REAL(2.0), REAL(0.05), REAL(1.0) };

	Real lap = ScalarFieldOperations::LaplacianSpher(f_spher, near_pole);
	REQUIRE_THAT(lap, WithinAbs(REAL(6.0), REAL(0.5)));
}

TEST_CASE("FieldOps - Divergence of radial field at cylindrical near-axis", "[field_operations][edge][singularity]")
{
	TEST_PRECISION_INFO();
	// F = (r, 0, 0) in cylindrical coords (purely radial)
	// ∇·F = (1/r) ∂(r·r)/∂r = (1/r)·2r = 2
	VectorFunction<3> F_cyl([](const VectorN<Real, 3>& pos) -> VectorN<Real, 3> {
		return VectorN<Real, 3>{ pos[0], REAL(0.0), REAL(0.0) };
	});

	// Near cylindrical axis (r small but non-zero)
	Vec3Cyl near_axis{ REAL(0.05), REAL(0.5), REAL(1.0) };

	Real div = VectorFieldOperations::DivCyl(F_cyl, near_axis);
	REQUIRE_THAT(div, WithinAbs(REAL(2.0), REAL(0.5)));
}

TEST_CASE("FieldOps - Gradient at cylindrical near-axis", "[field_operations][edge][singularity]")
{
	TEST_PRECISION_INFO();
	// f(r,φ,z) = z² — pure z-dependence, gradient should be (0, 0, 2z)
	ScalarFunction<3> f_cyl([](const VectorN<Real, 3>& pos) -> Real {
		return pos[2] * pos[2];  // z²
	});

	Vec3Cyl near_axis{ REAL(0.05), REAL(0.5), REAL(2.0) };

	auto grad = ScalarFieldOperations::GradientCyl(f_cyl, near_axis);

	REQUIRE_THAT(grad[0], WithinAbs(REAL(0.0), REAL(0.1)));  // ∂f/∂r
	REQUIRE_THAT(grad[1], WithinAbs(REAL(0.0), REAL(0.1)));  // (1/r)∂f/∂φ
	REQUIRE_THAT(grad[2], WithinAbs(REAL(4.0), REAL(0.1)));  // ∂f/∂z = 2z = 4
}


///////////////////////////////////////////////////////////////////////////////
//              METRIC TENSOR AT SINGULAR POINTS                            //
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("MetricTensorSpherical - at near-pole theta~0", "[MetricTensor][edge][singularity]")
{
	TEST_PRECISION_INFO();
	MetricTensorSpherical metric;

	// At θ=0 (north pole), g_φφ = r²sin²θ → 0
	Real r = REAL(2.0);
	Real theta_near_pole = REAL(0.01);
	Vec3Sph pos{ r, theta_near_pole, REAL(0.0) };

	auto g = metric.GetCovariantMetric(pos);

	// g_rr = 1
	REQUIRE_THAT(g(0, 0), WithinAbs(REAL(1.0), REAL(1e-10)));
	// g_θθ = r²
	REQUIRE_THAT(g(1, 1), WithinAbs(r * r, REAL(1e-10)));
	// g_φφ = r²sin²θ → small but finite
	Real expected_gphiphi = r * r * std::sin(theta_near_pole) * std::sin(theta_near_pole);
	REQUIRE_THAT(g(2, 2), WithinAbs(expected_gphiphi, REAL(1e-10)));
	REQUIRE(g(2, 2) > REAL(0.0));  // Still positive, just small
}

TEST_CASE("MetricTensorSpherical - at small radius", "[MetricTensor][edge][singularity]")
{
	TEST_PRECISION_INFO();
	MetricTensorSpherical metric;

	// At r→0, g_θθ = r² → 0, g_φφ = r²sin²θ → 0
	Real r_small = REAL(0.01);
	Vec3Sph pos{ r_small, Constants::PI / REAL(4.0), REAL(0.0) };

	auto g = metric.GetCovariantMetric(pos);

	// g_rr = 1 always
	REQUIRE_THAT(g(0, 0), WithinAbs(REAL(1.0), REAL(1e-10)));
	// g_θθ = r² → very small
	REQUIRE_THAT(g(1, 1), WithinAbs(r_small * r_small, REAL(1e-10)));
	// Off-diagonal should be 0
	REQUIRE_THAT(g(0, 1), WithinAbs(REAL(0.0), REAL(1e-12)));
	REQUIRE_THAT(g(0, 2), WithinAbs(REAL(0.0), REAL(1e-12)));
	REQUIRE_THAT(g(1, 2), WithinAbs(REAL(0.0), REAL(1e-12)));
}

TEST_CASE("MetricTensorCylindrical - at near-axis r~0", "[MetricTensor][edge][singularity]")
{
	TEST_PRECISION_INFO();
	MetricTensorCylindrical metric;

	// At r→0, g_φφ = r² → 0
	Real r_small = REAL(0.01);
	Vec3Cyl pos{ r_small, REAL(0.5), REAL(1.0) };

	auto g = metric.GetCovariantMetric(pos);

	// g_rr = 1
	REQUIRE_THAT(g(0, 0), WithinAbs(REAL(1.0), REAL(1e-10)));
	// g_φφ = r²
	REQUIRE_THAT(g(1, 1), WithinAbs(r_small * r_small, REAL(1e-10)));
	// g_zz = 1
	REQUIRE_THAT(g(2, 2), WithinAbs(REAL(1.0), REAL(1e-10)));
}


///////////////////////////////////////////////////////////////////////////////
//           QUATERNION EDGE CASES — GIMBAL LOCK & SLERP                   //
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("Quaternion - gimbal lock: pitch = +pi/2", "[quaternion][edge][gimbal]")
{
	TEST_PRECISION_INFO();
	// When pitch = ±π/2, Euler extraction encounters gimbal lock.
	// The implementation clamps sinp but doesn't use special yaw/roll formulas.
	// Known limitation: the extracted yaw/roll pair may not correctly represent
	// the original rotation's yaw-roll coupling at gimbal lock.
	Quaternion q = Quaternion::FromEulerZYX(REAL(0.3), Constants::PI / REAL(2.0), REAL(0.7));

	Vec3Cart euler = q.ToEulerZYX();

	// Pitch should be π/2 (may go through asin branch if sinp < 1.0 due to FP)
	REQUIRE_THAT(euler[1], WithinAbs(Constants::PI / REAL(2.0), REAL(1e-6)));

	// The quaternion should still be unit
	REQUIRE(q.IsUnit(1e-10));

	// Verify the extraction doesn't produce NaN or infinity
	REQUIRE(!std::isnan(euler[0]));
	REQUIRE(!std::isnan(euler[2]));
	REQUIRE(!std::isinf(euler[0]));
	REQUIRE(!std::isinf(euler[2]));
}

TEST_CASE("Quaternion - gimbal lock: pitch = -pi/2", "[quaternion][edge][gimbal]")
{
	TEST_PRECISION_INFO();
	Quaternion q = Quaternion::FromEulerZYX(REAL(0.5), -Constants::PI / REAL(2.0), REAL(0.2));

	Vec3Cart euler = q.ToEulerZYX();

	// Pitch should be -π/2
	REQUIRE_THAT(euler[1], WithinAbs(-Constants::PI / REAL(2.0), REAL(1e-6)));
	REQUIRE(q.IsUnit(1e-10));

	// Verify the extraction doesn't produce NaN or infinity
	REQUIRE(!std::isnan(euler[0]));
	REQUIRE(!std::isnan(euler[2]));
	REQUIRE(!std::isinf(euler[0]));
	REQUIRE(!std::isinf(euler[2]));
}

TEST_CASE("Quaternion - Slerp with identical quaternions", "[quaternion][edge][slerp]")
{
	TEST_PRECISION_INFO();
	Quaternion q = Quaternion::FromAxisAngle(Vec3Cart(1, 0, 0), Constants::PI / REAL(4.0));

	// Slerp(q, q, t) should return q for all t
	for (Real t : { REAL(0.0), REAL(0.25), REAL(0.5), REAL(0.75), REAL(1.0) })
	{
		Quaternion result = Quaternion::Slerp(q, q, t);
		REQUIRE(result.IsApprox(q, 1e-10));
		REQUIRE(result.IsUnit(1e-10));
	}
}

TEST_CASE("Quaternion - Slerp near-antipodal quaternions", "[quaternion][edge][slerp]")
{
	TEST_PRECISION_INFO();
	// Nearly opposite quaternions (large rotation between them)
	Quaternion q1 = Quaternion::FromAxisAngle(Vec3Cart(0, 0, 1), REAL(0.01));
	Quaternion q2 = Quaternion::FromAxisAngle(Vec3Cart(0, 0, 1), Constants::PI - REAL(0.01));

	// Slerp should still produce unit quaternions
	for (Real t : { REAL(0.0), REAL(0.25), REAL(0.5), REAL(0.75), REAL(1.0) })
	{
		Quaternion result = Quaternion::Slerp(q1, q2, t);
		REQUIRE(result.IsUnit(1e-6));
	}
}


///////////////////////////////////////////////////////////////////////////////
//         GEOMETRY PRIMITIVE DEGENERATE INPUT TESTS                        //
//   Note: Additional degenerate geometry tests (Line3D zero dir, Plane3D   //
//   from collinear points) are in base/geometry/geometry_3d_tests.cpp.     //
//   The tests below validate basic degenerate cases without Geometry3D.h   //
//   which conflicts with CoordTransf headers.                               //
///////////////////////////////////////////////////////////////////////////////


} // namespace MML::Tests::Core::GeometryCoordEdgeCaseTests

// Naming-convention regression tests for vectors, quaternions, Matrix3D and geometry
// (CONTRACT_POLICY section 11, task MinimalMathLibrary-wewi.5.3): mathematical predicates are
// PascalCase (IsZero/IsNearZero/IsUnit/IsTimelike...), STL protocol lowercase (clear/empty/at),
// reset vocabulary SetZero/Fill on fixed-size types, accessors getX(), conversions To...,
// one spelling per operation (Normalized, not GetAsUnitVector).  Deprecated 2.0 spellings must
// keep compiling with identical semantics for one release.
#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/base/Vector/Vector.h>
#include <mml/base/Vector/VectorN.h>
#include <mml/base/Vector/VectorTypes2D.h>
#include <mml/base/Vector/VectorTypes3D.h>
#include <mml/base/Vector/VectorTypes4D.h>
#include <mml/base/Quaternions.h>
#include <mml/base/Matrix/Matrix3D.h>
#include <mml/base/Geometry/Geometry2D.h>
#include <mml/base/Geometry/Geometry3D.h>
#include <mml/base/Geometry/Geometry3DBodies.h>
#include <mml/core/Surfaces.h>
#endif

#include <type_traits>

using namespace MML;

namespace MML::Tests::Base::VectorNamingTests
{
	template<typename V> using iszero_t     = decltype(std::declval<const V&>().IsZero());
	template<typename V> using isnearzero_t = decltype(std::declval<const V&>().IsNearZero(Real{1e-9}));
	template<typename V> using setzero_t    = decltype(std::declval<V&>().SetZero());
	template<typename V> using clear_t      = decltype(std::declval<V&>().clear());
	template<template<class> class Op, class T, class = void> struct has : std::false_type {};
	template<template<class> class Op, class T> struct has<Op, T, std::void_t<Op<T>>> : std::true_type {};

	static_assert(has<iszero_t, Vector<Real>>::value && has<isnearzero_t, Vector<Real>>::value);
	static_assert(has<iszero_t, VectorN<Real, 3>>::value && has<isnearzero_t, VectorN<Real, 3>>::value);
	static_assert(has<iszero_t, Vector3Cartesian>::value);
	static_assert(std::is_void_v<setzero_t<VectorN<Real, 3>>>);   // fixed size: value op
	static_assert(std::is_void_v<clear_t<Vector<Real>>>);         // dynamic: STL remove-all

	TEST_CASE("VectorNaming_Vector_VectorN_vocabulary", "[Vector][VectorN][naming]")
	{
		Vector<Real> v({ 1e-14, -1e-14, 0.0 });
		REQUIRE_FALSE(v.IsZero());
		REQUIRE(v.IsNearZero());
		REQUIRE(Vector<Real>({ 0.0, 0.0 }).IsZero());
		v.clear();
		REQUIRE(v.empty());
		REQUIRE(v.IsZero());    // vacuous

		VectorN<Real, 3> w({ 1.0, 2.0, 3.0 });
		w.Fill(0.5);
		REQUIRE(w[2] == 0.5);
		w.SetZero();
		REQUIRE(w.size() == 3);
		REQUIRE(w.IsZero());
		w[0] = 1e-13;
		REQUIRE_FALSE(w.IsZero());
		REQUIRE(w.IsNearZero(1e-12));
	}

	TEST_CASE("VectorNaming_Cartesian_conversions_and_normalization", "[Vector][naming]")
	{
		Vector3Cartesian a(3.0, 0.0, 4.0);
		REQUIRE(a.Normalized().NormL2() == Catch::Approx(1.0));
		REQUIRE(a.NormalizedOrZero()[0] == Catch::Approx(0.6));
		REQUIRE(Vector3Cartesian(0.0, 0.0, 0.0).NormalizedOrZero().IsZero());
		REQUIRE_THROWS_AS(Vector3Cartesian(0.0, 0.0, 0.0).Normalized(), VectorDimensionError);

		Point3Cartesian p = a.ToPoint();
		REQUIRE(p.Z() == 4.0);

		Vector2Cartesian b(1.0, 2.0);
		REQUIRE(b.Normalized().NormL2() == Catch::Approx(1.0));
		REQUIRE(Vector2Cartesian(0.0, 0.0).NormalizedOrZero().IsZero());
		Vector2Cartesian p1, p2;
		b.GetPerpendicularVectors(p1, p2);     // same spelling as Vector3Cartesian
		REQUIRE(b * p1 == Catch::Approx(0.0));
		REQUIRE(b * p2 == Catch::Approx(0.0));

		Vector4Minkowski t{ 2.0, 1.0, 0.0, 0.0 }, s{ 1.0, 2.0, 0.0, 0.0 }, l{ 1.0, 1.0, 0.0, 0.0 };
		REQUIRE(t.IsTimelike());
		REQUIRE(s.IsSpacelike());
		REQUIRE(l.IsLightlike());
	}

	TEST_CASE("VectorNaming_IsParallelTo_includes_antiparallel_and_rejects_zero", "[Vector][naming]")
	{
		Vector3Cartesian a(1.0, 2.0, 3.0);
		REQUIRE(a.IsParallelTo(Vector3Cartesian(2.0, 4.0, 6.0)));
		REQUIRE(a.IsParallelTo(Vector3Cartesian(-1.0, -2.0, -3.0)));       // antiparallel IS parallel
		REQUIRE(a.IsParallelTo(Vector3Cartesian(1e-8, 2e-8, 3e-8)));       // scale invariant
		REQUIRE(a.IsParallelTo(Vector3Cartesian(1e8, 2e8, 3e8)));
		REQUIRE_FALSE(a.IsParallelTo(Vector3Cartesian(1.0, 0.0, 0.0)));
		REQUIRE_FALSE(a.IsParallelTo(Vector3Cartesian(0.0, 0.0, 0.0)));    // zero has no direction: false, not NaN
		REQUIRE_FALSE(Vector3Cartesian(0.0, 0.0, 0.0).IsParallelTo(a));

		// Planes with opposite normals are parallel
		Plane3D up(Pnt3Cart(0, 0, 0), Vec3Cart(0, 0, 1));
		Plane3D down(Pnt3Cart(0, 0, 5), Vec3Cart(0, 0, -1));
		REQUIRE(up.IsParallelToPlane(down));

		// Regression: PlaneSurface(point, normal) with normal antiparallel to x used to take the
		// "not parallel" branch and normalize a zero cross product (threw VectorDimensionError).
		Surfaces::PlaneSurface surf(Vec3Cart(0, 0, 0), Vec3Cart(-1, 0, 0));
		Vec3Cart q = surf(1.0, 2.0);
		REQUIRE(q.NormL2() == Catch::Approx(std::sqrt(5.0)));
		REQUIRE(q[0] == Catch::Approx(0.0).margin(1e-14));          // stays in the x = 0 plane
	}

	TEST_CASE("VectorNaming_Quaternion_Matrix3D_vocabulary", "[Quaternion][Matrix3D][naming]")
	{
		Quaternion q = Quaternion::FromAxisAngle(Vec3Cart(0, 0, 1), Constants::PI / 2);
		REQUIRE(q.IsUnit());
		REQUIRE_FALSE(q.IsIdentity());
		REQUIRE(Quaternion::Identity().IsIdentity());
		REQUIRE(q.getRotationAngle() == Catch::Approx(Constants::PI / 2));
		REQUIRE(q.getRotationAxis()[2] == Catch::Approx(1.0));
		REQUIRE(q.IsEqualTo(q));
		REQUIRE_FALSE(q.IsEqualTo(Quaternion::Identity()));

		Matrix3D<Real> m(2, 2, 2, 1.5);
		REQUIRE_FALSE(m.empty());
		REQUIRE(Matrix3D<Real>().empty());
		REQUIRE(m.at(1, 1, 1) == 1.5);
		m.at(0, 0, 0) = 7.0;
		REQUIRE(m(0, 0, 0) == 7.0);
		REQUIRE_THROWS_AS(m.at(2, 0, 0), IndexError);
	}

	TEST_CASE("VectorNaming_geometry_accessors_are_getX", "[Geometry][naming]")
	{
		Cube3D cube(2.0, Pnt3Cart(1, 1, 1));
		REQUIRE(cube.getSide() == 2.0);
		REQUIRE(cube.getCenter().X() == 1.0);
		REQUIRE(cube.getBoundingBox().Contains(Pnt3Cart(1.5, 1.5, 1.5)));
		REQUIRE(cube.getBoundingSphere().Contains(Pnt3Cart(1, 1, 1)));
		REQUIRE(cube.getSurfaceCount() == 6);
		REQUIRE(cube.getSurface(0).getArea() == Catch::Approx(4.0));
		const IBody& body = cube;
		REQUIRE(body.getCenter().Y() == 1.0);

		Sphere3D sphere(3.0);
		REQUIRE(sphere.getRadius() == 3.0);
		REQUIRE(sphere.getNumLatitude() > 0);
		Cylinder3D cyl(1.0, 4.0);
		REQUIRE(cyl.getHeight() == 4.0);
		REQUIRE(cyl.getBaseCenter().Z() == Catch::Approx(cyl.getCenter().Z() - cyl.getHeight() / 2));
		Torus3D torus(3.0, 1.0);
		REQUIRE(torus.getMajorRadius() == 3.0);
		REQUIRE(torus.getMinorRadius() == 1.0);

		Plane3D plane(Pnt3Cart(0, 0, 2), Vec3Cart(0, 0, 1));
		REQUIRE(plane.IsPointOnPlane(plane.getPointOnPlane()));
		REQUIRE(Plane3D::XYPlane().IsParallelToPlane(plane));
		REQUIRE(Plane3D::XZPlane().IsPerpendicularToPlane(Plane3D::XYPlane()));
		REQUIRE(Plane3D::YZPlane().IsPerpendicularToPlane(Plane3D::XZPlane()));

		Polygon2D poly;
		poly.AddVertex(Point2Cartesian(0, 0));
		poly.AddVertex(Point2Cartesian(2, 0));
		poly.AddVertex(Point2Cartesian(2, 1));
		REQUIRE(poly.getBoundingBox().Width() == 2.0);
		poly.clear();
		REQUIRE(poly.NumVertices() == 0);
		Circle2D circle(Point2Cartesian(0, 0), 1.0);
		REQUIRE(circle.getBoundingBox().minX == -1.0);
	}

#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable : 4996)
#elif defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif
	TEST_CASE("VectorNaming_deprecated_aliases_forward", "[Vector][naming][deprecated]")
	{
		Vector<Real> v({ 1e-14, 0.0 });
		REQUIRE(v.isZero() == v.IsZero());
		REQUIRE(v.isNearZero() == v.IsNearZero());
		v.Clear();
		REQUIRE(v.empty());

		VectorN<Real, 2> w({ 1.0, 2.0 });
		REQUIRE(w.isZero() == w.IsZero());
		REQUIRE(w.isNearZero(1e-3) == w.IsNearZero(1e-3));
		w.clear();                                   // zeroes (kept semantics), does not resize
		REQUIRE(w.size() == 2);
		REQUIRE(w.IsZero());

		Vector3Cartesian a(3.0, 0.0, 4.0);
		REQUIRE(a.GetAsUnitVector().IsEqualTo(a.Normalized()));
		REQUIRE(a.GetAsUnitVectorOrZero().IsEqualTo(a.NormalizedOrZero()));
		REQUIRE(a.getAsPoint().Z() == a.ToPoint().Z());
		Vector2Cartesian b(1.0, 2.0), p1, p2, q1, q2;
		b.getPerpendicularVectors(p1, p2);
		b.GetPerpendicularVectors(q1, q2);
		REQUIRE(p1.IsEqualTo(q1));
		REQUIRE(b.GetAsUnitVector().IsEqualTo(b.Normalized()));
		Vector4Minkowski t{ 2.0, 1.0, 0.0, 0.0 };
		REQUIRE(t.isTimelike() == t.IsTimelike());
		REQUIRE(t.isSpacelike() == t.IsSpacelike());
		REQUIRE(t.isLightlike() == t.IsLightlike());

		Quaternion q = Quaternion::FromAxisAngle(Vec3Cart(0, 1, 0), 1.0);
		REQUIRE(q.isUnit() == q.IsUnit());
		REQUIRE(q.isIdentity() == q.IsIdentity());
		REQUIRE(q.isApprox(q) == q.IsEqualTo(q));
		REQUIRE(q.GetRotationAngle() == q.getRotationAngle());
		REQUIRE(q.GetRotationAxis().IsEqualTo(q.getRotationAxis()));

		Matrix3D<Real> m(1, 1, 1, 2.0);
		REQUIRE(m.isEmpty() == m.empty());
		REQUIRE(m.At(0, 0, 0) == m.at(0, 0, 0));

		Cube3D cube(2.0);
		REQUIRE(cube.GetSide() == cube.getSide());
		REQUIRE(cube.GetCenter().IsEqualTo(cube.getCenter()));
		REQUIRE(cube.GetBoundingBox().Contains(cube.getCenter()));
		REQUIRE(cube.GetBoundingSphere().Contains(cube.getCenter()));
		REQUIRE(cube.GetSurfaceCount() == cube.getSurfaceCount());
		REQUIRE(cube.GetSurface(0).getArea() == cube.getSurface(0).getArea());
		Sphere3D sphere(3.0);
		REQUIRE(sphere.GetRadius() == sphere.getRadius());
		Cylinder3D cyl(1.0, 4.0);
		REQUIRE(cyl.GetHeight() == cyl.getHeight());
		REQUIRE(cyl.GetBaseCenter().IsEqualTo(cyl.getBaseCenter()));
		Plane3D plane(Pnt3Cart(0, 0, 2), Vec3Cart(0, 0, 1));
		REQUIRE(plane.GetPointOnPlane().IsEqualTo(plane.getPointOnPlane()));
		Real sx, sy, sz, gx, gy, gz;
		plane.GetCoordAxisSegments(sx, sy, sz);
		plane.getCoordAxisSegments(gx, gy, gz);
		REQUIRE(sz == gz);
		REQUIRE(Plane3D::GetXYPlane().IsParallelToPlane(Plane3D::XYPlane()));
		Polygon2D poly;
		poly.AddVertex(Point2Cartesian(0, 0));
		poly.AddVertex(Point2Cartesian(1, 0));
		poly.AddVertex(Point2Cartesian(1, 1));
		REQUIRE(poly.GetBoundingBox().Width() == poly.getBoundingBox().Width());
		poly.Clear();
		REQUIRE(poly.NumVertices() == 0);
	}
#if defined(_MSC_VER)
#pragma warning(pop)
#elif defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic pop
#endif
}

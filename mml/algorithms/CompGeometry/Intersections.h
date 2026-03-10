///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        comp_geometry/Intersections.h                                       ///
///  Description: Intersection algorithms for segments, lines, and rays             ///
///               Delegates to base Geometry classes for core algorithms             ///
///               Provides convenience overloads for point-based APIs                ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_COMP_GEOMETRY_INTERSECTIONS_H
#define MML_COMP_GEOMETRY_INTERSECTIONS_H

#include "mml/algorithms/CompGeometry/CompGeometryBase.h"
#include "mml/MMLBase.h"
#include "mml/base/Geometry/Geometry2D.h"
#include "mml/base/Geometry/Geometry3D.h"

#include <cmath>

namespace MML {
namespace CompGeometry {

	/// @brief Intersection algorithms for geometric primitives
	/// @details Provides static methods for intersection tests.
	///          Core algorithms are implemented in the geometry classes (Line2D, Line3D, etc.)
	///          This class provides convenience overloads for point-based APIs and
	///          additional algorithms like ray-triangle intersection.
	class Intersections {
	private:
		static constexpr Real EPSILON = Constants::GEOMETRY_EPSILON;

	public:
		// ========================================================================
		// Segment-Segment Intersection (2D)
		// ========================================================================
		// Delegates to SegmentLine2D::Intersection()
		// ========================================================================

		/// @brief Intersect two segments given by endpoints
		/// @param p1, q1 Endpoints of first segment
		/// @param p2, q2 Endpoints of second segment
		/// @return SegmentIntersection with type, point, or overlap endpoints
		static SegmentIntersection IntersectSegments(
			const Point2Cartesian& p1, const Point2Cartesian& q1,
			const Point2Cartesian& p2, const Point2Cartesian& q2)
		{
			SegmentLine2D seg1(p1, q1);
			SegmentLine2D seg2(p2, q2);
			return seg1.Intersection(seg2);
		}

		/// @brief Intersect two SegmentLine2D objects
		static SegmentIntersection IntersectSegments(const SegmentLine2D& seg1, const SegmentLine2D& seg2) {
			return seg1.Intersection(seg2);
		}

		// ========================================================================
		// Line-Line Intersection (2D)
		// ========================================================================
		// Delegates to Line2D::Intersection()
		// ========================================================================

		/// @brief Intersect two 2D lines given by point and direction
		static LineIntersection2D IntersectLines2D(
			const Point2Cartesian& p1, const Vector2Cartesian& d1,
			const Point2Cartesian& p2, const Vector2Cartesian& d2)
		{
			Line2D line1(p1, d1);
			Line2D line2(p2, d2);
			return line1.Intersection(line2);
		}

		/// @brief Intersect two Line2D objects
		static LineIntersection2D IntersectLines2D(const Line2D& line1, const Line2D& line2) {
			return line1.Intersection(line2);
		}

		// ========================================================================
		// Line-Line Intersection (3D)
		// ========================================================================
		// Delegates to Line3D::Intersection()
		// ========================================================================

		/// @brief Intersect two 3D lines given by point and direction
		static LineIntersection3D IntersectLines3D(
			const Point3Cartesian& p1, const Vector3Cartesian& d1,
			const Point3Cartesian& p2, const Vector3Cartesian& d2)
		{
			Line3D line1(p1, d1);
			Line3D line2(p2, d2);
			return line1.Intersection(line2);
		}

		/// @brief Intersect two Line3D objects
		static LineIntersection3D IntersectLines3D(const Line3D& line1, const Line3D& line2) {
			return line1.Intersection(line2);
		}

		// ========================================================================
		// Ray-Triangle Intersection (Moller-Trumbore Algorithm)
		// ========================================================================
		// Fast ray-triangle intersection using barycentric coordinates.
		// Complexity: O(1)
		// Reference: Moller & Trumbore, "Fast, Minimum Storage Ray-Triangle Intersection"
		// ========================================================================

		/// @brief Intersect ray with triangle (Moller-Trumbore algorithm)
		/// @param rayOrigin Origin of the ray
		/// @param rayDir Direction of the ray (not necessarily unit)
		/// @param v0, v1, v2 Triangle vertices
		/// @param cullBackface If true, backface hits are rejected
		/// @return RayTriangleHit with hit status, t parameter, barycentric coords
		static RayTriangleHit IntersectRayTriangle(
			const Point3Cartesian& rayOrigin,
			const Vector3Cartesian& rayDir,
			const Point3Cartesian& v0,
			const Point3Cartesian& v1,
			const Point3Cartesian& v2,
			bool cullBackface = false)
		{
			RayTriangleHit result;

			Vector3Cartesian edge1(v0, v1);
			Vector3Cartesian edge2(v0, v2);
			Vector3Cartesian pvec = VectorProduct(rayDir, edge2);
			Real det = ScalarProduct(edge1, pvec);

			if (cullBackface) {
				if (det < EPSILON) return result;
			} else {
				if (std::abs(det) < EPSILON) return result;
			}

			Real invDet = 1.0 / det;
			Vector3Cartesian tvec(v0, rayOrigin);

			Real u = ScalarProduct(tvec, pvec) * invDet;
			if (u < 0.0 || u > 1.0) return result;

			Vector3Cartesian qvec = VectorProduct(tvec, edge1);
			Real v = ScalarProduct(rayDir, qvec) * invDet;
			if (v < 0.0 || u + v > 1.0) return result;

			Real t = ScalarProduct(edge2, qvec) * invDet;
			if (t < 0.0) return result;

			result.hit = true;
			result.t = t;
			result.u = u;
			result.v = v;
			result.point = rayOrigin + t * rayDir;
			return result;
		}

		/// @brief Intersect ray with Triangle3D
		static RayTriangleHit IntersectRayTriangle(
			const Point3Cartesian& rayOrigin,
			const Vector3Cartesian& rayDir,
			const Triangle3D& triangle,
			bool cullBackface = false)
		{
			return IntersectRayTriangle(rayOrigin, rayDir, 
				triangle.Pnt1(), triangle.Pnt2(), triangle.Pnt3(), cullBackface);
		}

		/// @brief Intersect Line3D (as ray) with triangle vertices
		static RayTriangleHit IntersectRayTriangle(
			const Line3D& ray,
			const Point3Cartesian& v0,
			const Point3Cartesian& v1,
			const Point3Cartesian& v2,
			bool cullBackface = false)
		{
			return IntersectRayTriangle(ray.StartPoint(), ray.Direction(), v0, v1, v2, cullBackface);
		}

		/// @brief Intersect Line3D (as ray) with Triangle3D
		static RayTriangleHit IntersectRayTriangle(
			const Line3D& ray,
			const Triangle3D& triangle,
			bool cullBackface = false)
		{
			return IntersectRayTriangle(ray.StartPoint(), ray.Direction(), 
				triangle.Pnt1(), triangle.Pnt2(), triangle.Pnt3(), cullBackface);
		}
	};

} // namespace CompGeometry
} // namespace MML

#endif // MML_COMP_GEOMETRY_INTERSECTIONS_H

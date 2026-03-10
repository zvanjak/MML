///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        comp_geometry/CompGeometryBase.h                                    ///
///  Description: Type definitions for computational geometry algorithms             ///
///               Enums and result structs used across multiple algorithms           ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                        ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_COMP_GEOMETRY_BASE_H
#define MML_COMP_GEOMETRY_BASE_H

#include "mml/MMLBase.h"
#include "mml/base/Geometry/Geometry2D.h"
#include "mml/base/Geometry/Geometry3D.h"

#include <tuple>

namespace MML {

	// ============================================================================
	// NOTE: Basic intersection types (LineIntersection2D, LineIntersection3D,
	// SegmentIntersection) are now defined in base/Geometry/Geometry2D.h and
	// base/Geometry/Geometry3D.h as they are fundamental to the geometry classes.
	// ============================================================================

	// ============================================================================
	// Ray-Triangle Intersection Result
	// ============================================================================

	/// @brief Result of ray-triangle intersection test
	/// @details Contains hit status, distance along ray, barycentric coordinates,
	///          and the 3D intersection point if the ray hits the triangle.
	struct RayTriangleHit {
		bool hit = false;          ///< True if ray hits the triangle
		Real t = 0.0;              ///< Distance along ray to hit point (ray_origin + t * direction)
		Real u = 0.0;              ///< Barycentric coordinate u (weight for v1)
		Real v = 0.0;              ///< Barycentric coordinate v (weight for v2)
		Point3Cartesian point;     ///< Intersection point in 3D space

		/// @brief Check if ray hit the triangle
		bool Hit() const { return hit; }

		/// @brief Get barycentric coordinate for v0 (1 - u - v)
		Real BarycentricW() const { return 1.0 - u - v; }

		/// @brief Get barycentric coordinates as a tuple (u, v, w) where w = 1-u-v
		/// @details Point = w*v0 + u*v1 + v*v2
		std::tuple<Real, Real, Real> GetBarycentricCoords() const {
			return std::make_tuple(u, v, 1.0 - u - v);
		}
	};

	// ============================================================================
	// Closest Pair Result
	// ============================================================================

	/// Result of closest pair of points algorithm
	struct ClosestPairResult {
		Point2Cartesian p1, p2;
		Real distance;

		ClosestPairResult()
			: distance(std::numeric_limits<Real>::max()) {}
		ClosestPairResult(const Point2Cartesian& a, const Point2Cartesian& b, Real d)
			: p1(a), p2(b), distance(d) {}
	};

} // namespace MML

#endif // MML_COMP_GEOMETRY_BASE_H

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

#include <mml/MMLBase.h>
#include <mml/base/Geometry/Geometry2D.h>
#include <mml/base/Geometry/Geometry3D.h>
#include <mml/base/Geometry/GeometryBase/GeometryPredicates.h>

#include <tuple>

namespace MML {

	namespace detail {
		inline GeometryBinary GeometryDistance2D(const Point2Cartesian& a, const Point2Cartesian& b) {
			const std::array<GeometryBinary, 2> difference{
				GeometryBinary::Difference(a.X(), b.X()), GeometryBinary::Difference(a.Y(), b.Y())};
			return GeometryBinary::SumProducts(difference, difference).Sqrt();
		}

		inline Vector2Cartesian GeometryUnitDifference2D(const Point2Cartesian& a, const Point2Cartesian& b) {
			const auto x = GeometryBinary::Difference(b.X(), a.X()), y = GeometryBinary::Difference(b.Y(), a.Y());
			if (x.fraction == 0 && y.fraction == 0)
				throw GeometryError("Geometry direction - coincident points");
			const int exponent = x.fraction == 0 ? y.exponent : y.fraction == 0 ? x.exponent : std::max(x.exponent, y.exponent);
			const Real dx = std::ldexp(x.fraction, x.exponent - exponent);
			const Real dy = std::ldexp(y.fraction, y.exponent - exponent);
			const Real length = std::hypot(dx, dy);
			return {dx / length, dy / length};
		}

		inline bool GeometryCircumcenter2D(const Point2Cartesian& a, const Point2Cartesian& b,
			const Point2Cartesian& c, Point2Cartesian& center) {
			using B = GeometryBinary;
			auto determinant = GeometryExact2D::Determinant(a, b, c);
			if (determinant.fraction == 0)
				return false;
			++determinant.exponent;
			const auto ax = B::Difference(b.X(), a.X()), ay = B::Difference(b.Y(), a.Y());
			const auto bx = B::Difference(c.X(), a.X()), by = B::Difference(c.Y(), a.Y());
			const auto aSquared = B::SumProducts<2>({ax, ay}, {ax, ay});
			const auto bSquared = B::SumProducts<2>({bx, by}, {bx, by});
			// Local circumcenter equations; separate binary exponents avoid squared-length overflow.
			const auto x = B::Quotient(B::DifferenceOfProducts(by, aSquared, ay, bSquared), determinant);
			const auto y = B::Quotient(B::DifferenceOfProducts(ax, bSquared, bx, aSquared), determinant);
			center = {B::Add(B(a.X()), x).CheckedValue(), B::Add(B(a.Y()), y).CheckedValue()};
			return true;
		}
	}

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
		Real BarycentricW() const { return MML::Real(1.0) - u - v; }

		/// @brief Get barycentric coordinates as a tuple (u, v, w) where w = 1-u-v
		/// @details Point = w*v0 + u*v1 + v*v2
		std::tuple<Real, Real, Real> GetBarycentricCoords() const {
			return std::make_tuple(u, v, static_cast<MML::Real>(1.0) - u - v);
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

///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        comp_geometry/ConvexHull.h                                          ///
///  Description: 2D Convex Hull algorithm (Andrew's Monotone Chain)                 ///
///               Complexity: O(n log n)                                              ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_COMP_GEOMETRY_CONVEX_HULL_H
#define MML_COMP_GEOMETRY_CONVEX_HULL_H

#include "mml/MMLBase.h"
#include "mml/base/Geometry/Geometry2D.h"

#include <algorithm>
#include <vector>

namespace MML {
namespace CompGeometry {

	/// @brief 2D Convex Hull using Andrew's Monotone Chain Algorithm
	/// @details Computes the convex hull of a set of 2D points.
	///          Complexity: O(n log n) time, O(n) space
	class ConvexHull2D {
	private:
		// Use centralized geometry epsilon from Constants
		static constexpr Real EPSILON = Constants::GEOMETRY_EPSILON;

		/// Cross product of vectors OA and OB (2D cross product = z-component of 3D cross)
		static Real Cross(const Point2Cartesian& O, const Point2Cartesian& A, const Point2Cartesian& B) {
			return (A.X() - O.X()) * (B.Y() - O.Y()) - (A.Y() - O.Y()) * (B.X() - O.X());
		}

	public:
		// ========================================================================
		// Convex Hull - Andrew's Monotone Chain Algorithm
		// ========================================================================
		// Computes convex hull of a set of 2D points in O(n log n) time.
		// Returns vertices in counter-clockwise order.
		// ========================================================================

		/// Compute convex hull of a set of points
		/// @param points Input points (will be sorted internally)
		/// @return Polygon2D representing the convex hull in CCW order
		static Polygon2D Compute(std::vector<Point2Cartesian> points) {
			int n = static_cast<int>(points.size());
			if (n < 3) {
				// Degenerate cases
				if (n == 0)
					return Polygon2D();
				if (n == 1)
					return Polygon2D({points[0]});
				if (n == 2)
					return Polygon2D({points[0], points[1]});
			}

			// Sort points lexicographically (by x, then by y)
			std::sort(points.begin(), points.end(), [](const Point2Cartesian& a, const Point2Cartesian& b) {
				return (a.X() < b.X()) || (a.X() == b.X() && a.Y() < b.Y());
			});

			// Remove duplicates
			auto last = std::unique(points.begin(), points.end(), [](const Point2Cartesian& a, const Point2Cartesian& b) {
				return std::abs(a.X() - b.X()) < EPSILON && std::abs(a.Y() - b.Y()) < EPSILON;
			});
			points.erase(last, points.end());
			n = static_cast<int>(points.size());

			if (n < 3) {
				return Polygon2D(points);
			}

			std::vector<Point2Cartesian> hull(2 * n);
			int k = 0;

			// Build lower hull
			for (int i = 0; i < n; ++i) {
				while (k >= 2 && Cross(hull[k - 2], hull[k - 1], points[i]) <= EPSILON)
					--k;
				hull[k++] = points[i];
			}

			// Build upper hull
			int lowerSize = k + 1;
			for (int i = n - 2; i >= 0; --i) {
				while (k >= lowerSize && Cross(hull[k - 2], hull[k - 1], points[i]) <= EPSILON)
					--k;
				hull[k++] = points[i];
			}

			// Remove last point (it's the same as first)
			hull.resize(k - 1);

			return Polygon2D(hull);
		}

		/// Convenience method for Polygon2D input
		static Polygon2D Compute(const Polygon2D& polygon) { 
			return Compute(polygon.Vertices()); 
		}
	};

} // namespace CompGeometry
} // namespace MML

#endif // MML_COMP_GEOMETRY_CONVEX_HULL_H

#ifndef MML_POLYGON_OPS_H
#define MML_POLYGON_OPS_H

/////////////////////////////////////////////////////////////////////////////////////////
///  @file      PolygonOps.h
///  @brief     Polygon Clipping and Boolean Operations
///  @details   Part of the MML Computational Geometry module
///  @warning   Boolean operations are APPROXIMATE for non-convex polygons!
/////////////////////////////////////////////////////////////////////////////////////////

#include "mml/MMLBase.h"
#include "mml/base/Geometry/Geometry2D.h"
#include "mml/algorithms/CompGeometry/CompGeometryBase.h"
#include "mml/algorithms/CompGeometry/ConvexHull.h"
#include "mml/algorithms/CompGeometry/Intersections.h"

#include <algorithm>
#include <vector>

namespace MML::CompGeometry {

//////////////////////////////////////////////////////////////////////////////////////////
/// @brief Polygon clipping and boolean operations
/// @details Provides Sutherland-Hodgman clipping (exact for convex clip regions)
///          and approximate boolean operations for general polygons.
/// @warning Boolean ops are APPROXIMATE for non-convex polygons! Use Clipper2/CGAL
///          for exact operations on arbitrary polygons.
//////////////////////////////////////////////////////////////////////////////////////////
class PolygonOps {
	// Use centralized geometry epsilon from Constants
	static constexpr Real EPSILON = Constants::GEOMETRY_EPSILON;

public:
	// ========================================================================
	// Sutherland-Hodgman Polygon Clipping
	// ========================================================================
	// Clips a subject polygon against a CONVEX clipping polygon.
	// Time complexity: O(n * m) where n = subject vertices, m = clip edges
	// ========================================================================

	/// @brief Clip a polygon against a convex clipping region
	/// @param subject The polygon to clip
	/// @param clipRegion The convex clipping polygon (must be convex!)
	/// @return The clipped polygon, or empty polygon if fully clipped
	static Polygon2D ClipPolygon(const Polygon2D& subject, const Polygon2D& clipRegion) {
		if (subject.NumVertices() < 3 || clipRegion.NumVertices() < 3)
			return Polygon2D();

		std::vector<Point2Cartesian> output = subject.Vertices();

		int clipN = clipRegion.NumVertices();
		for (int i = 0; i < clipN; i++) {
			if (output.empty())
				break;

			const Point2Cartesian& edgeStart = clipRegion[i];
			const Point2Cartesian& edgeEnd = clipRegion[(i + 1) % clipN];

			output = ClipAgainstEdge(output, edgeStart, edgeEnd);
		}

		return output.size() >= 3 ? Polygon2D(output) : Polygon2D();
	}

	/// @brief Clip polygon to axis-aligned rectangle
	/// @param subject The polygon to clip
	/// @param minX Left boundary
	/// @param minY Bottom boundary
	/// @param maxX Right boundary
	/// @param maxY Top boundary
	/// @return The clipped polygon
	static Polygon2D ClipToRect(const Polygon2D& subject, Real minX, Real minY, Real maxX, Real maxY) {
		// CCW order: bottom-left -> bottom-right -> top-right -> top-left
		Polygon2D rect({
			Point2Cartesian(minX, minY), // bottom-left
			Point2Cartesian(maxX, minY), // bottom-right
			Point2Cartesian(maxX, maxY), // top-right
			Point2Cartesian(minX, maxY)	 // top-left
		});
		return ClipPolygon(subject, rect);
	}

	// ========================================================================
	// Boolean Polygon Operations - APPROXIMATE (Convex Hull Fallback)
	// ========================================================================
	// @warning APPROXIMATE IMPLEMENTATIONS - READ CAREFULLY!
	//
	// These operations are EXACT only for CONVEX polygons. For non-convex
	// polygons, they use CONVEX HULL APPROXIMATION which may:
	//   - Include regions that should be excluded (intersection)
	//   - Exclude regions that should be included (difference)
	//   - Produce incorrect topology (union)
	//
	// For robust exact boolean operations on general (non-convex) polygons,
	// use a dedicated computational geometry library such as:
	//   - CGAL (Computational Geometry Algorithms Library)
	//   - Clipper/Clipper2 (Angus Johnson's polygon clipping library)
	//   - boost::geometry
	//
	// Current implementation:
	//   - Convex + Convex: EXACT (Sutherland-Hodgman clipping)
	//   - Any + Convex clip: EXACT for intersection (Sutherland-Hodgman)
	//   - Non-convex cases: APPROXIMATE (convex hull of relevant points)
	//
	// Time complexity: O((n + k) log n) where k = intersections
	// ========================================================================

	/// @brief Compute approximate intersection of two polygons (A ∩ B)
	/// @warning APPROXIMATE for non-convex polygons! Uses convex hull fallback.
	/// @details For convex polygons or convex clip region B, this is EXACT using
	///          Sutherland-Hodgman algorithm. For non-convex cases, returns the
	///          convex hull of intersection points as an APPROXIMATION.
	/// @param polyA First polygon (subject)
	/// @param polyB Second polygon (clip region)
	/// @return Vector of result polygons (usually 0 or 1 for approximate case)
	/// @note For exact boolean ops on non-convex polygons, use Clipper2 or CGAL.
	static std::vector<Polygon2D> PolygonIntersectionApprox(const Polygon2D& polyA, const Polygon2D& polyB) {
		std::vector<Polygon2D> result;

		// Simple case: use Sutherland-Hodgman for convex clip region
		if (polyB.IsConvex()) {
			Polygon2D clipped = ClipPolygon(polyA, polyB);
			if (clipped.NumVertices() >= 3)
				result.push_back(clipped);
			return result;
		}

		// Fallback: if both are convex, result is convex
		if (polyA.IsConvex() && polyB.IsConvex()) {
			Polygon2D clipped = ClipPolygon(polyA, polyB);
			if (clipped.NumVertices() >= 3)
				result.push_back(clipped);
			return result;
		}

		// For general non-convex case, use point-in-polygon tests
		// Build intersection by collecting vertices of A inside B and vice versa
		std::vector<Point2Cartesian> intersectionPts;

		// Add vertices of A that are inside B
		for (int i = 0; i < polyA.NumVertices(); i++) {
			if (polyB.Contains(polyA[i]))
				intersectionPts.push_back(polyA[i]);
		}

		// Add vertices of B that are inside A
		for (int i = 0; i < polyB.NumVertices(); i++) {
			if (polyA.Contains(polyB[i]))
				intersectionPts.push_back(polyB[i]);
		}

		// Add intersection points
		for (int i = 0; i < polyA.NumVertices(); i++) {
			for (int j = 0; j < polyB.NumVertices(); j++) {
				auto isect = Intersections::IntersectSegments(
					polyA[i], polyA[(i + 1) % polyA.NumVertices()], 
					polyB[j], polyB[(j + 1) % polyB.NumVertices()]);
				if (isect.IsPointIntersection())
					intersectionPts.push_back(isect.point);
			}
		}

		if (intersectionPts.size() >= 3) {
			// WARNING: APPROXIMATE - Form convex hull of intersection points
			// This is NOT the true intersection for non-convex polygons!
			Polygon2D hull = ConvexHull2D::Compute(intersectionPts);
			if (hull.NumVertices() >= 3)
				result.push_back(hull);
		}

		return result;
	}

	/// @brief Compute approximate union of two polygons (A ∪ B)
	/// @warning APPROXIMATE! Returns convex hull for non-convex polygons.
	/// @details For two convex polygons, the union's convex hull is correct.
	///          For non-convex polygons, this returns the convex hull which
	///          INCLUDES regions that should be excluded (concavities).
	/// @param polyA First polygon
	/// @param polyB Second polygon
	/// @return Vector containing single approximate result polygon
	/// @note For exact boolean ops on non-convex polygons, use Clipper2 or CGAL.
	static std::vector<Polygon2D> PolygonUnionApprox(const Polygon2D& polyA, const Polygon2D& polyB) {
		std::vector<Polygon2D> result;

		// Collect all vertices
		std::vector<Point2Cartesian> allPts;
		for (int i = 0; i < polyA.NumVertices(); i++)
			allPts.push_back(polyA[i]);
		for (int i = 0; i < polyB.NumVertices(); i++)
			allPts.push_back(polyB[i]);

		// For convex polygons, union is convex hull (this is exact for convex!)
		if (polyA.IsConvex() && polyB.IsConvex()) {
			Polygon2D hull = ConvexHull2D::Compute(allPts);
			if (hull.NumVertices() >= 3)
				result.push_back(hull);
			return result;
		}

		// WARNING: APPROXIMATE - For non-convex case, return convex hull
		// This INCLUDES concave regions that should be excluded!
		// Full non-convex union requires complex boundary tracing (Weiler-Atherton)
		Polygon2D hull = ConvexHull2D::Compute(allPts);
		if (hull.NumVertices() >= 3)
			result.push_back(hull);

		return result;
	}

	/// @brief Compute approximate difference of two polygons (A - B)
	/// @warning APPROXIMATE! Uses convex hull of exterior points.
	/// @details This finds points of A outside B and intersection points,
	///          then returns their convex hull. This is NOT the true difference
	///          for non-convex polygons and may produce incorrect results.
	/// @param polyA Polygon to subtract from
	/// @param polyB Polygon to subtract
	/// @return Vector of result polygons (may be empty, or approximate)
	/// @note For exact boolean ops on non-convex polygons, use Clipper2 or CGAL.
	static std::vector<Polygon2D> PolygonDifferenceApprox(const Polygon2D& polyA, const Polygon2D& polyB) {
		std::vector<Polygon2D> result;

		// Check if polygons don't overlap
		auto intersection = PolygonIntersectionApprox(polyA, polyB);
		if (intersection.empty()) {
			// No overlap - return A unchanged
			result.push_back(polyA);
			return result;
		}

		// Collect vertices of A that are outside B
		std::vector<Point2Cartesian> outsidePts;
		for (int i = 0; i < polyA.NumVertices(); i++) {
			if (!polyB.Contains(polyA[i]))
				outsidePts.push_back(polyA[i]);
		}

		// Add intersection points on boundary of A
		for (int i = 0; i < polyA.NumVertices(); i++) {
			for (int j = 0; j < polyB.NumVertices(); j++) {
				auto isect = Intersections::IntersectSegments(
					polyA[i], polyA[(i + 1) % polyA.NumVertices()], 
					polyB[j], polyB[(j + 1) % polyB.NumVertices()]);
				if (isect.IsPointIntersection())
					outsidePts.push_back(isect.point);
			}
		}

		if (outsidePts.size() >= 3) {
			// WARNING: APPROXIMATE - Form polygon from outside points using convex hull
			// Proper difference may result in multiple non-convex polygons!
			Polygon2D hull = ConvexHull2D::Compute(outsidePts);
			if (hull.NumVertices() >= 3)
				result.push_back(hull);
		}

		return result;
	}

private:
	// ========================================================================
	// Helper functions
	// ========================================================================

	// Cross product helper for orientation tests
	static Real Cross(const Point2Cartesian& o, const Point2Cartesian& a, const Point2Cartesian& b) {
		return (a.X() - o.X()) * (b.Y() - o.Y()) - (a.Y() - o.Y()) * (b.X() - o.X());
	}

	// Intersect segment (p1->p2) with infinite line through (l1, l2)
	// Returns intersection point and whether it exists
	static std::pair<Point2Cartesian, bool> IntersectSegmentWithLine(const Point2Cartesian& p1, const Point2Cartesian& p2,
																	 const Point2Cartesian& l1, const Point2Cartesian& l2) {
		// Direction vectors
		Real dx = p2.X() - p1.X();
		Real dy = p2.Y() - p1.Y();
		Real lx = l2.X() - l1.X();
		Real ly = l2.Y() - l1.Y();

		Real denom = dx * ly - dy * lx;
		if (std::abs(denom) < EPSILON) {
			// Parallel - no intersection
			return {Point2Cartesian(), false};
		}

		// Parameter t for segment p1->p2
		Real t = ((l1.X() - p1.X()) * ly - (l1.Y() - p1.Y()) * lx) / denom;

		// Only accept if within segment bounds (with small tolerance)
		if (t < -EPSILON || t > 1.0 + EPSILON) {
			return {Point2Cartesian(), false};
		}

		// Clamp to [0, 1]
		t = std::max(Real(0), std::min(Real(1), t));

		Point2Cartesian intersection(p1.X() + t * dx, p1.Y() + t * dy);
		return {intersection, true};
	}

	// Clip polygon against a single edge (half-plane)
	static std::vector<Point2Cartesian> ClipAgainstEdge(const std::vector<Point2Cartesian>& polygon, 
														const Point2Cartesian& edgeStart,
														const Point2Cartesian& edgeEnd) {
		std::vector<Point2Cartesian> output;
		int n = static_cast<int>(polygon.size());
		if (n == 0)
			return output;

		// Edge direction for inside test (left side is inside for CCW clip polygon)
		auto isInside = [&](const Point2Cartesian& p) { return Cross(edgeStart, edgeEnd, p) >= -EPSILON; };

		for (int i = 0; i < n; i++) {
			const Point2Cartesian& current = polygon[i];
			const Point2Cartesian& next = polygon[(i + 1) % n];

			bool currentInside = isInside(current);
			bool nextInside = isInside(next);

			if (currentInside) {
				output.push_back(current);
				if (!nextInside) {
					// Exiting: add intersection with clip edge (infinite line)
					auto [point, found] = IntersectSegmentWithLine(current, next, edgeStart, edgeEnd);
					if (found)
						output.push_back(point);
				}
			} else if (nextInside) {
				// Entering: add intersection with clip edge (infinite line)
				auto [point, found] = IntersectSegmentWithLine(current, next, edgeStart, edgeEnd);
				if (found)
					output.push_back(point);
			}
		}

		return output;
	}
};

} // namespace MML::CompGeometry

#endif // MML_POLYGON_OPS_H

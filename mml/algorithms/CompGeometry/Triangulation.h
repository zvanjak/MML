#ifndef MML_TRIANGULATION_H
#define MML_TRIANGULATION_H

/////////////////////////////////////////////////////////////////////////////////////////
///  @file      Triangulation.h
///  @brief     2D Triangulation Algorithms: Ear Clipping and Delaunay
///  @details   Part of the MML Computational Geometry module
/////////////////////////////////////////////////////////////////////////////////////////

#include "mml/MMLBase.h"
#include "mml/base/Geometry/Geometry2D.h"
#include "mml/algorithms/CompGeometry/CompGeometryBase.h"

#include <algorithm>
#include <array>
#include <map>
#include <random>
#include <vector>

namespace MML::CompGeometry {

//////////////////////////////////////////////////////////////////////////////////////////
/// @brief Delaunay triangulation result structure
/// @details Contains the triangulated mesh with adjacency information
//////////////////////////////////////////////////////////////////////////////////////////
struct DelaunayTriangulation {
	std::vector<Point2Cartesian> points;               ///< Input points
	std::vector<std::array<int, 3>> triangles;         ///< Triangle vertex indices (CCW order)
	std::vector<std::array<int, 3>> neighbors;         ///< Adjacent triangle indices (-1 if no neighbor)

	/// @brief Get number of triangles in the triangulation
	int NumTriangles() const { return static_cast<int>(triangles.size()); }

	/// @brief Get number of points
	int NumPoints() const { return static_cast<int>(points.size()); }

	/// @brief Get the triangle vertices as Point2Cartesian
	std::array<Point2Cartesian, 3> GetTrianglePoints(int triIndex) const {
		const auto& tri = triangles[triIndex];
		return {points[tri[0]], points[tri[1]], points[tri[2]]};
	}

	/// @brief Get triangle as Triangle2D object
	Triangle2D GetTriangle(int triIndex) const {
		const auto& tri = triangles[triIndex];
		return Triangle2D(points[tri[0]], points[tri[1]], points[tri[2]]);
	}

	/// @brief Find the triangle containing a given point
	/// @param p The query point
	/// @return Index of containing triangle, or -1 if outside convex hull
	int LocatePoint(const Point2Cartesian& p) const {
		for (int i = 0; i < NumTriangles(); i++) {
			const auto& tri = triangles[i];
			const Point2Cartesian& a = points[tri[0]];
			const Point2Cartesian& b = points[tri[1]];
			const Point2Cartesian& c = points[tri[2]];

			// Check if point is inside triangle using barycentric coordinates
			Real d1 = Sign(p, a, b);
			Real d2 = Sign(p, b, c);
			Real d3 = Sign(p, c, a);

			bool hasNeg = (d1 < 0) || (d2 < 0) || (d3 < 0);
			bool hasPos = (d1 > 0) || (d2 > 0) || (d3 > 0);

			if (!(hasNeg && hasPos))
				return i;
		}
		return -1;
	}

	/// @brief Interpolate values at the input points to a query point
	/// @param p Query point
	/// @param values Values at each input point (same order as points)
	/// @param defaultVal Value to return if p is outside triangulation
	/// @return Linearly interpolated value using barycentric coordinates
	Real Interpolate(const Point2Cartesian& p, const std::vector<Real>& values, 
					 Real defaultVal = 0.0) const {
		constexpr Real EPSILON = 1e-10;
		if (values.size() != points.size())
			return defaultVal;

		int triIdx = LocatePoint(p);
		if (triIdx < 0)
			return defaultVal;

		const auto& tri = triangles[triIdx];
		const Point2Cartesian& a = points[tri[0]];
		const Point2Cartesian& b = points[tri[1]];
		const Point2Cartesian& c = points[tri[2]];

		// Compute barycentric coordinates
		Real denom = (b.Y() - c.Y()) * (a.X() - c.X()) + (c.X() - b.X()) * (a.Y() - c.Y());
		if (std::abs(denom) < EPSILON)
			return defaultVal;

		Real w1 = ((b.Y() - c.Y()) * (p.X() - c.X()) + (c.X() - b.X()) * (p.Y() - c.Y())) / denom;
		Real w2 = ((c.Y() - a.Y()) * (p.X() - c.X()) + (a.X() - c.X()) * (p.Y() - c.Y())) / denom;
		Real w3 = 1.0 - w1 - w2;

		return w1 * values[tri[0]] + w2 * values[tri[1]] + w3 * values[tri[2]];
	}

	/// @brief Get convex hull edges (boundary of triangulation)
	/// @return Vector of vertex index pairs representing boundary edges
	std::vector<std::pair<int, int>> GetConvexHullEdges() const {
		std::vector<std::pair<int, int>> hull;
		for (size_t i = 0; i < triangles.size(); i++) {
			for (int j = 0; j < 3; j++) {
				if (neighbors[i][j] == -1) {
					// This edge is on the boundary
					int v1 = triangles[i][j];
					int v2 = triangles[i][(j + 1) % 3];
					hull.push_back({v1, v2});
				}
			}
		}
		return hull;
	}

private:
	static Real Sign(const Point2Cartesian& p1, const Point2Cartesian& p2, 
					 const Point2Cartesian& p3) {
		return (p1.X() - p3.X()) * (p2.Y() - p3.Y()) - 
			   (p2.X() - p3.X()) * (p1.Y() - p3.Y());
	}
};

//////////////////////////////////////////////////////////////////////////////////////////
/// @brief 2D Triangulation algorithms
/// @details Provides ear clipping for simple polygons and Delaunay triangulation
//////////////////////////////////////////////////////////////////////////////////////////
class Triangulation {
	// Use centralized geometry epsilon from Constants
	static constexpr Real EPSILON = Constants::GEOMETRY_EPSILON;

public:
	// ========================================================================
	// Ear Clipping Triangulation
	// ========================================================================
	// Triangulates a simple polygon using the ear clipping method.
	// Works for both convex and non-convex simple polygons.
	// Time complexity: O(n²)
	// ========================================================================

	/// @brief Triangulate a simple polygon using ear clipping
	/// @param polygon The input polygon (must be simple - no self-intersections)
	/// @return Vector of triangles forming the triangulation
	static std::vector<Triangle2D> EarClipTriangulation(const Polygon2D& polygon) {
		std::vector<Triangle2D> triangles;
		int n = polygon.NumVertices();

		if (n < 3)
			return triangles;

		if (n == 3) {
			triangles.push_back(Triangle2D(polygon[0], polygon[1], polygon[2]));
			return triangles;
		}

		// Copy vertices to working list
		std::vector<Point2Cartesian> vertices = polygon.Vertices();

		// Ensure CCW orientation (ear clipping assumes CCW)
		Real signedArea = 0;
		for (int i = 0; i < n; i++) {
			int j = (i + 1) % n;
			signedArea += vertices[i].X() * vertices[j].Y();
			signedArea -= vertices[i].Y() * vertices[j].X();
		}
		if (signedArea < 0) {
			std::reverse(vertices.begin(), vertices.end());
		}

		// Clip ears until only a triangle remains
		while (vertices.size() > 3) {
			n = static_cast<int>(vertices.size());
			bool earFound = false;

			for (int i = 0; i < n; i++) {
				if (IsEar(vertices, i, n)) {
					// Found an ear - add triangle and remove vertex
					int prev = (i - 1 + n) % n;
					int next = (i + 1) % n;

					triangles.push_back(Triangle2D(vertices[prev], vertices[i], vertices[next]));
					vertices.erase(vertices.begin() + i);
					earFound = true;
					break;
				}
			}

			if (!earFound) {
				// Polygon may be self-intersecting or degenerate
				// Fall back to fan triangulation for remaining vertices
				for (size_t i = 1; i < vertices.size() - 1; i++) {
					triangles.push_back(Triangle2D(vertices[0], vertices[i], vertices[i + 1]));
				}
				break;
			}
		}

		// Add final triangle
		if (vertices.size() == 3) {
			triangles.push_back(Triangle2D(vertices[0], vertices[1], vertices[2]));
		}

		return triangles;
	}

	// ========================================================================
	// Delaunay Triangulation (Bowyer-Watson Algorithm)
	// ========================================================================
	// Computes the Delaunay triangulation of a point set.
	// Time complexity: O(n log n) expected, O(n²) worst case
	// ========================================================================

	/// @brief Compute Delaunay triangulation using Bowyer-Watson algorithm
	/// @param points Input point set (at least 3 non-collinear points)
	/// @return DelaunayTriangulation structure with triangles and adjacency
	static DelaunayTriangulation ComputeDelaunay(std::vector<Point2Cartesian> points) {
		DelaunayTriangulation result;

		int n = static_cast<int>(points.size());
		if (n < 3) {
			result.points = std::move(points);
			return result;
		}

		// Find bounding box
		Real minX = points[0].X(), maxX = points[0].X();
		Real minY = points[0].Y(), maxY = points[0].Y();
		for (const auto& p : points) {
			minX = std::min(minX, p.X());
			maxX = std::max(maxX, p.X());
			minY = std::min(minY, p.Y());
			maxY = std::max(maxY, p.Y());
		}

		// Create super-triangle that contains all points
		// Make it large enough to ensure all points fit well inside
		Real dx = maxX - minX;
		Real dy = maxY - minY;
		Real deltaMax = std::max(dx, dy);
		
		// Ensure deltaMax is large enough for numerical stability
		deltaMax = std::max<Real>(deltaMax, std::max<Real>(std::abs(minX), std::abs(minY)) * Real(1e-10));
		deltaMax = std::max<Real>(deltaMax, std::max<Real>(std::abs(maxX), std::abs(maxY)) * Real(1e-10));
		deltaMax = std::max<Real>(deltaMax, Real(1e-10));  // Absolute minimum
		
		Real midX = (minX + maxX) / Real(2.0);
		Real midY = (minY + maxY) / Real(2.0);

		// Super-triangle vertices in CCW order (very large to contain all points)
		Point2Cartesian st1(midX - 20 * deltaMax, midY - deltaMax);
		Point2Cartesian st2(midX + 20 * deltaMax, midY - deltaMax);
		Point2Cartesian st3(midX, midY + 20 * deltaMax);

		// Working point list: original points then super-triangle vertices
		std::vector<Point2Cartesian> workPoints = points;
		int st1Idx = static_cast<int>(workPoints.size());
		int st2Idx = st1Idx + 1;
		int st3Idx = st1Idx + 2;
		workPoints.push_back(st1);
		workPoints.push_back(st2);
		workPoints.push_back(st3);

		// Triangle representation: vertex indices and validity flag
		struct Tri {
			std::array<int, 3> v;
			bool valid = true;
		};
		std::vector<Tri> tris;
		tris.push_back({{st1Idx, st2Idx, st3Idx}, true});

		// Incrementally add each point
		for (int i = 0; i < n; i++) {
			const Point2Cartesian& p = workPoints[i];

			// Find all triangles whose circumcircle contains the new point
			std::vector<int> badTriangles;
			for (size_t t = 0; t < tris.size(); t++) {
				if (!tris[t].valid) continue;

				const Point2Cartesian& a = workPoints[tris[t].v[0]];
				const Point2Cartesian& b = workPoints[tris[t].v[1]];
				const Point2Cartesian& c = workPoints[tris[t].v[2]];

				if (InCircumcircle(p, a, b, c)) {
					badTriangles.push_back(static_cast<int>(t));
				}
			}

			// Find boundary polygon: edges that appear exactly once
			std::map<EdgeKey, int> edgeCount;
			std::map<EdgeKey, std::pair<int, int>> edgeDir;  // Original direction
			
			for (int t : badTriangles) {
				for (int j = 0; j < 3; j++) {
					int v1 = tris[t].v[j];
					int v2 = tris[t].v[(j + 1) % 3];
					EdgeKey key(v1, v2);
					edgeCount[key]++;
					edgeDir[key] = {v1, v2};
				}
			}

			// Collect boundary edges (edges appearing exactly once)
			std::vector<std::pair<int, int>> polygon;
			for (const auto& [key, count] : edgeCount) {
				if (count == 1) {
					polygon.push_back(edgeDir[key]);
				}
			}

			// Remove bad triangles
			for (int t : badTriangles) {
				tris[t].valid = false;
			}

			// Create new triangles from the new point to each boundary edge
			for (const auto& [v1, v2] : polygon) {
				tris.push_back({{i, v1, v2}, true});
			}
		}

		// Extract final triangulation (exclude triangles using super-triangle vertices)
		result.points = points;  // Only original points
		
		for (const auto& tri : tris) {
			if (!tri.valid) continue;

			// Skip if any vertex is a super-triangle vertex
			bool usesSuperTri = false;
			for (int j = 0; j < 3; j++) {
				if (tri.v[j] >= n) {
					usesSuperTri = true;
					break;
				}
			}

			if (!usesSuperTri) {
				result.triangles.push_back(tri.v);
			}
		}

		// Build neighbor information
		BuildNeighborInfo(result);

		return result;
	}

	/// @brief Compute Delaunay triangulation and return as vector of Triangle2D
	/// @param points Input point set
	/// @return Vector of triangles
	static std::vector<Triangle2D> DelaunayTriangulate(const std::vector<Point2Cartesian>& points) {
		auto dt = ComputeDelaunay(std::vector<Point2Cartesian>(points));
		std::vector<Triangle2D> result;
		result.reserve(dt.triangles.size());
		for (int i = 0; i < dt.NumTriangles(); i++) {
			result.push_back(dt.GetTriangle(i));
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

	// Check if vertex at index i is an "ear" in the polygon
	static bool IsEar(const std::vector<Point2Cartesian>& vertices, int i, int n) {
		int prev = (i - 1 + n) % n;
		int next = (i + 1) % n;

		const Point2Cartesian& a = vertices[prev];
		const Point2Cartesian& b = vertices[i];
		const Point2Cartesian& c = vertices[next];

		// Check if the triangle is convex (CCW orientation means convex for CCW polygon)
		if (Cross(a, b, c) <= EPSILON)
			return false; // Reflex vertex, not an ear

		// Check that no other vertex is inside this triangle
		for (int j = 0; j < n; j++) {
			if (j == prev || j == i || j == next)
				continue;

			if (IsPointInTriangle(vertices[j], a, b, c))
				return false;
		}

		return true;
	}

	// Point-in-triangle test using barycentric coordinates
	static bool IsPointInTriangle(const Point2Cartesian& p, const Point2Cartesian& a, 
								  const Point2Cartesian& b, const Point2Cartesian& c) {
		Real d1 = Cross(a, b, p);
		Real d2 = Cross(b, c, p);
		Real d3 = Cross(c, a, p);

		bool hasNeg = (d1 < -EPSILON) || (d2 < -EPSILON) || (d3 < -EPSILON);
		bool hasPos = (d1 > EPSILON) || (d2 > EPSILON) || (d3 > EPSILON);

		return !(hasNeg && hasPos);
	}

	// Check if point is inside circumcircle of triangle (Delaunay property)
	static bool InCircumcircle(const Point2Cartesian& p,
							   const Point2Cartesian& a,
							   const Point2Cartesian& b,
							   const Point2Cartesian& c) {
		// For a CCW-oriented triangle abc, point p is inside circumcircle if:
		// | ax-px  ay-py  (ax-px)²+(ay-py)² |
		// | bx-px  by-py  (bx-px)²+(by-py)² | > 0
		// | cx-px  cy-py  (cx-px)²+(cy-py)² |
		
		Real ax = a.X() - p.X();
		Real ay = a.Y() - p.Y();
		Real bx = b.X() - p.X();
		Real by = b.Y() - p.Y();
		Real cx = c.X() - p.X();
		Real cy = c.Y() - p.Y();
		
		Real ap = ax * ax + ay * ay;
		Real bp = bx * bx + by * by;
		Real cp = cx * cx + cy * cy;

		Real det = ax * (by * cp - bp * cy) -
				   ay * (bx * cp - bp * cx) +
				   ap * (bx * cy - by * cx);

		return det > EPSILON;
	}

	// Unique edge representation for Bowyer-Watson
	struct EdgeKey {
		int v1, v2;
		EdgeKey(int a, int b) : v1(std::min(a, b)), v2(std::max(a, b)) {}
		bool operator==(const EdgeKey& other) const { 
			return v1 == other.v1 && v2 == other.v2; 
		}
		bool operator<(const EdgeKey& other) const {
			return v1 < other.v1 || (v1 == other.v1 && v2 < other.v2);
		}
	};

	// Build adjacency information for triangulation
	static void BuildNeighborInfo(DelaunayTriangulation& dt) {
		int numTri = static_cast<int>(dt.triangles.size());
		dt.neighbors.resize(numTri);

		// Initialize all neighbors to -1
		for (auto& n : dt.neighbors) {
			n = {-1, -1, -1};
		}

		// Build edge-to-triangle map
		std::map<EdgeKey, std::vector<std::pair<int, int>>> edgeToTri;  // edge -> [(triIdx, edgeIdx)]
		for (int i = 0; i < numTri; i++) {
			for (int j = 0; j < 3; j++) {
				int v1 = dt.triangles[i][j];
				int v2 = dt.triangles[i][(j + 1) % 3];
				EdgeKey key(v1, v2);
				edgeToTri[key].push_back({i, j});
			}
		}

		// Set neighbors
		for (const auto& [edge, triList] : edgeToTri) {
			if (triList.size() == 2) {
				int t1 = triList[0].first;
				int e1 = triList[0].second;
				int t2 = triList[1].first;
				int e2 = triList[1].second;
				dt.neighbors[t1][e1] = t2;
				dt.neighbors[t2][e2] = t1;
			}
		}
	}
};

} // namespace MML::CompGeometry

#endif // MML_TRIANGULATION_H

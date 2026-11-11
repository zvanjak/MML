#ifndef MML_TRIANGULATION_H
#define MML_TRIANGULATION_H

/////////////////////////////////////////////////////////////////////////////////////////
///  @file      Triangulation.h
///  @brief     2D Triangulation Algorithms: Ear Clipping and Delaunay
///  @details   Part of the MML Computational Geometry module
/////////////////////////////////////////////////////////////////////////////////////////

#include <mml/MMLBase.h>
#include <mml/base/AlgorithmTypes.h>
#include <mml/base/Geometry/Geometry2D.h>
#include <mml/algorithms/CompGeometry/CompGeometryBase.h>
#include <mml/algorithms/CompGeometry/RobustPredicates.h>

#include <algorithm>
#include <array>
#include <map>
#include <limits>
#include <numeric>
#include <string>
#include <utility>
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
		MML::detail::ValidateGeometryCoordinates({p.X(), p.Y()}, "DelaunayTriangulation::LocatePoint");
		for (int i = 0; i < NumTriangles(); i++) {
			const auto& tri = triangles[i];
			const Point2Cartesian& a = points[tri[0]];
			const Point2Cartesian& b = points[tri[1]];
			const Point2Cartesian& c = points[tri[2]];

			if (RobustPredicates::Orientation2D(a, b, c) == 0) continue;
			int d1 = RobustPredicates::Orientation2D(p, a, b);
			int d2 = RobustPredicates::Orientation2D(p, b, c);
			int d3 = RobustPredicates::Orientation2D(p, c, a);

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
		if (values.size() != points.size())
			throw ArgumentError("DelaunayTriangulation::Interpolate - one value per input point required");
		for (Real value : values)
			MML::detail::ValidateGeometryCoordinates({value}, "DelaunayTriangulation::Interpolate");

		int triIdx = LocatePoint(p);
		if (triIdx < 0)
			return defaultVal;

		const auto& tri = triangles[triIdx];
		const Point2Cartesian& a = points[tri[0]];
		const Point2Cartesian& b = points[tri[1]];
		const Point2Cartesian& c = points[tri[2]];

		if (values[tri[0]] == values[tri[1]] && values[tri[1]] == values[tri[2]])
			return values[tri[0]];
		using B = MML::detail::GeometryBinary;
		using Exact = MML::detail::GeometryExact2D;
		const std::array<B, 3> areas{Exact::Determinant(p, b, c), Exact::Determinant(a, p, c), Exact::Determinant(a, b, p)};
		const std::array<B, 3> samples{B(values[tri[0]]), B(values[tri[1]]), B(values[tri[2]])};
		const auto numerator = B::SumProducts(areas, samples);
		return B::Quotient(numerator, Exact::Determinant(a, b, c)).CheckedValue();
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

};

//////////////////////////////////////////////////////////////////////////////////////////
/// @brief 2D Triangulation algorithms
/// @details Provides ear clipping for simple polygons and Delaunay triangulation
//////////////////////////////////////////////////////////////////////////////////////////
class Triangulation {
public:
	// ========================================================================
	// Ear Clipping Triangulation
	// ========================================================================
	// Triangulates a simple polygon using the ear clipping method.
	// Works for both convex and non-convex simple polygons.
	// Time complexity: O(n²)
	// ========================================================================

	/// @brief Result of ear-clipping triangulation with diagnostics
	struct EarClipResult {
		std::vector<Triangle2D> triangles;      ///< Triangles produced
		bool used_fan_fallback = false;         ///< True if no ear existed at some stage (input not a simple polygon)
		int unclipped_vertices = 0;             ///< Vertices remaining when the fallback engaged (0 if none)
		AlgorithmStatus status = AlgorithmStatus::Success;
		std::string error_message;
		bool IsSuccess() const noexcept { return status == AlgorithmStatus::Success; }
	};

	/// @brief Triangulate a simple polygon using ear clipping, with diagnostics
	/// @param polygon The input polygon (must be simple - no self-intersections)
	/// @details For a simple polygon the two-ears theorem guarantees an ear at every stage.
	///          If none is found the input is self-intersecting or degenerate; the remaining
	///          vertices are fan-triangulated so the output is still a complete set of triangles,
	///          but the result is flagged (status InvalidInput, used_fan_fallback = true) because
	///          a fan of a non-simple polygon is not a triangulation of it. Until analysis v4
	///          (ALGORITHMS A10) this fallback was silent.
	static EarClipResult EarClipTriangulationDetailed(const Polygon2D& polygon) {
		EarClipResult result;
		std::vector<Triangle2D>& triangles = result.triangles;
		int n = polygon.NumVertices();

		if (n < 3)
			return result;

		// Copy vertices to working list
		std::vector<Point2Cartesian> vertices = polygon.Vertices();
		for (const auto& p : vertices)
			MML::detail::ValidateGeometryCoordinates({p.X(), p.Y()}, "EarClipTriangulation");
		// At a lexicographic extreme of a simple polygon, the local turn
		// has the polygon's winding, without global shoelace cancellation.
		const auto extreme = std::min_element(vertices.begin(), vertices.end(), [](const auto& a, const auto& b) {
			return a.X() < b.X() || (a.X() == b.X() && a.Y() < b.Y());
		});
		const int pivot = static_cast<int>(extreme - vertices.begin());
		if (RobustPredicates::Orientation2D(vertices[(pivot + n - 1) % n], vertices[pivot], vertices[(pivot + 1) % n]) < 0) {
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
				// Not a simple polygon (self-intersecting or degenerate): fan the rest and say so
				result.used_fan_fallback = true;
				result.unclipped_vertices = static_cast<int>(vertices.size());
				result.status = AlgorithmStatus::InvalidInput;
				result.error_message = "EarClipTriangulation: no ear found with " + std::to_string(vertices.size()) +
					" vertices remaining - polygon is not simple (self-intersecting or degenerate); remaining vertices fan-triangulated";
				for (size_t i = 1; i < vertices.size() - 1; i++) {
					triangles.push_back(Triangle2D(vertices[0], vertices[i], vertices[i + 1]));
				}
				return result;
			}
		}

		// Add final triangle
		if (vertices.size() == 3) {
			triangles.push_back(Triangle2D(vertices[0], vertices[1], vertices[2]));
		}

		return result;
	}

	/// @brief Triangulate a simple polygon using ear clipping
	/// @param polygon The input polygon (must be simple - no self-intersections)
	/// @return Vector of triangles forming the triangulation
	/// @throws GeometryError if the polygon is not simple (no ear could be found). Use
	///         EarClipTriangulationDetailed to receive the fan-triangulated fallback with a status instead.
	static std::vector<Triangle2D> EarClipTriangulation(const Polygon2D& polygon) {
		EarClipResult result = EarClipTriangulationDetailed(polygon);
		if (result.used_fan_fallback)
			throw GeometryError(result.error_message);
		return std::move(result.triangles);
	}

	// ========================================================================
	// Delaunay Triangulation (real-point hull insertion and Lawson flips)
	// ========================================================================
	// Computes the Delaunay triangulation of a point set.
	// Sorting O(n log n); hull scans and legalization can be quadratic.
	// ========================================================================

	/// @brief Compute Delaunay triangulation without finite auxiliary vertices
	/// @param points Finite input points, preserved in their original order.
	/// Exact duplicates use their first input index; collinear sets have no faces.
	/// @return DelaunayTriangulation structure with triangles and adjacency
	static DelaunayTriangulation ComputeDelaunay(std::vector<Point2Cartesian> points) {
		DelaunayTriangulation result;

		if (points.size() > static_cast<std::size_t>(std::numeric_limits<int>::max() / 2))
			throw ArgumentError("ComputeDelaunay - mesh too large for integer indices");
		for (const auto& p : points)
			MML::detail::ValidateGeometryCoordinates({p.X(), p.Y()}, "ComputeDelaunay");
		result.points = std::move(points);
		const auto& input = result.points;
		std::vector<int> order(input.size());
		std::iota(order.begin(), order.end(), 0);
		std::sort(order.begin(), order.end(), [&](int a, int b) {
			if (input[a].X() != input[b].X()) return input[a].X() < input[b].X();
			if (input[a].Y() != input[b].Y()) return input[a].Y() < input[b].Y();
			return a < b;
		});
		order.erase(std::unique(order.begin(), order.end(), [&](int a, int b) {
			return input[a].X() == input[b].X() && input[a].Y() == input[b].Y();
		}), order.end());
		if (order.size() < 3) return result;
		std::size_t seed = 2;
		while (seed < order.size() && RobustPredicates::Orientation2D(input[order[0]], input[order[1]], input[order[seed]]) == 0)
			++seed;
		if (seed == order.size()) return result;

		std::map<EdgeKey, std::pair<int, int>> boundary;
		const auto addFace = [&](int a, int b, int c) {
			if (RobustPredicates::Orientation2D(input[a], input[b], input[c]) < 0) std::swap(a, b);
			const std::array<int, 3> face{a, b, c};
			result.triangles.push_back(face);
			for (int e = 0; e < 3; ++e) {
				const int first = face[e], second = face[(e + 1) % 3];
				const EdgeKey key(first, second);
				if (boundary.erase(key) == 0) boundary.emplace(key, std::make_pair(first, second));
			}
		};
		// Retain every vertex of the initial collinear chain, not just its endpoints.
		for (std::size_t i = 1; i < seed; ++i)
			addFace(order[i - 1], order[i], order[seed]);
		for (std::size_t i = seed + 1; i < order.size(); ++i) {
			// Lexicographic insertion puts each new point outside the current hull.
			std::vector<std::pair<int, int>> visible;
			for (const auto& entry : boundary) {
				const auto [a, b] = entry.second;
				if (RobustPredicates::Orientation2D(input[a], input[b], input[order[i]]) < 0)
					visible.push_back({a, b});
			}
			if (visible.empty()) throw GeometryError("ComputeDelaunay - no visible hull edge for a distinct point");
			for (const auto& [a, b] : visible) addFace(a, order[i], b);
		}
		LegalizeEdges(result);
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

	// Check if vertex at index i is an "ear" in the polygon
	static bool IsEar(const std::vector<Point2Cartesian>& vertices, int i, int n) {
		int prev = (i - 1 + n) % n;
		int next = (i + 1) % n;

		const Point2Cartesian& a = vertices[prev];
		const Point2Cartesian& b = vertices[i];
		const Point2Cartesian& c = vertices[next];

		// Check if the triangle is convex (CCW orientation means convex for CCW polygon)
		if (RobustPredicates::Orientation2D(a, b, c) <= 0)
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
		int d1 = RobustPredicates::Orientation2D(a, b, p);
		int d2 = RobustPredicates::Orientation2D(b, c, p);
		int d3 = RobustPredicates::Orientation2D(c, a, p);

		bool hasNeg = (d1 < 0) || (d2 < 0) || (d3 < 0);
		bool hasPos = (d1 > 0) || (d2 > 0) || (d3 > 0);

		return !(hasNeg && hasPos);
	}

	// Undirected edge key; triangle-local edges retain their CCW direction.
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

	using EdgeMap = std::map<EdgeKey, std::vector<std::pair<int, int>>>;

	static EdgeMap BuildEdgeMap(const DelaunayTriangulation& dt) {
		EdgeMap edges;
		for (int t = 0; t < dt.NumTriangles(); ++t)
			for (int e = 0; e < 3; ++e)
				edges[EdgeKey(dt.triangles[t][e], dt.triangles[t][(e + 1) % 3])].push_back({t, e});
		for (const auto& entry : edges)
			if (entry.second.size() > 2) throw GeometryError("ComputeDelaunay - nonmanifold edge");
		return edges;
	}

	static void LegalizeEdges(DelaunayTriangulation& dt) {
		auto edges = BuildEdgeMap(dt);
		std::vector<EdgeKey> pending;
		for (const auto& entry : edges)
			if (entry.second.size() == 2) pending.push_back(entry.first);
		const auto detach = [&](int t) {
			for (int e = 0; e < 3; ++e) {
				const EdgeKey key(dt.triangles[t][e], dt.triangles[t][(e + 1) % 3]);
				auto& owners = edges.at(key);
				owners.erase(std::remove_if(owners.begin(), owners.end(), [t](const auto& owner) { return owner.first == t; }), owners.end());
				if (owners.empty()) edges.erase(key);
			}
		};
		// Lawson (1977): flip strictly illegal, convex interior edges; exact
		// cocircular ties stay unchanged, so legalization cannot oscillate.
		while (!pending.empty()) {
			const auto key = pending.back();
			pending.pop_back();
			const auto found = edges.find(key);
			if (found == edges.end() || found->second.size() != 2) continue;
			const auto [first, e] = found->second[0];
			const auto [second, f] = found->second[1];
			const auto triangle = dt.triangles[first];
			const int a = triangle[e], b = triangle[(e + 1) % 3], c = triangle[(e + 2) % 3];
			const int d = dt.triangles[second][(f + 2) % 3];
			const auto& p = dt.points;
			if (RobustPredicates::Orientation2D(p[c], p[d], p[b]) <= 0 ||
				RobustPredicates::Orientation2D(p[d], p[c], p[a]) <= 0) continue;
			if (RobustPredicates::InCircle2D(p[a], p[b], p[c], p[d]) <= 0) continue;
			detach(first);
			detach(second);
			dt.triangles[first] = {c, d, b};
			dt.triangles[second] = {d, c, a};
			for (int t : {first, second})
				for (int j = 0; j < 3; ++j) {
					const EdgeKey changed(dt.triangles[t][j], dt.triangles[t][(j + 1) % 3]);
					edges[changed].push_back({t, j});
					pending.push_back(changed);
				}
		}
	}

	// Build adjacency information for triangulation
	static void BuildNeighborInfo(DelaunayTriangulation& dt) {
		int numTri = static_cast<int>(dt.triangles.size());
		dt.neighbors.resize(numTri);

		// Initialize all neighbors to -1
		for (auto& n : dt.neighbors) {
			n = {-1, -1, -1};
		}

		// Build edge-to-triangle map
		const auto edgeToTri = BuildEdgeMap(dt);

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

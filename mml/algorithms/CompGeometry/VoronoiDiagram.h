#ifndef MML_VORONOI_DIAGRAM_H
#define MML_VORONOI_DIAGRAM_H

/////////////////////////////////////////////////////////////////////////////////////////
///  @file      VoronoiDiagram.h
///  @brief     Voronoi Diagram computation (dual of Delaunay triangulation)
///  @details   Part of the MML Computational Geometry module
/////////////////////////////////////////////////////////////////////////////////////////

#include <mml/MMLBase.h>
#include <mml/base/Geometry/Geometry2D.h>
#include <mml/base/Vector/VectorN.h>
#include <mml/algorithms/CompGeometry/CompGeometryBase.h>
#include <mml/algorithms/CompGeometry/Triangulation.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <set>
#include <vector>

namespace MML::CompGeometry {

//////////////////////////////////////////////////////////////////////////////////////////
/// @brief Voronoi edge connecting two Voronoi vertices (or extending to infinity)
//////////////////////////////////////////////////////////////////////////////////////////
struct VoronoiEdge {
	int v1;                      ///< First vertex index (-1 if extends to infinity)
	int v2;                      ///< Second vertex index (-1 if extends to infinity)
	int site1;                   ///< First adjacent site
	int site2;                   ///< Second adjacent site
	Vector2Cartesian direction;  ///< Direction for infinite edges
	
	bool IsInfinite() const { return v1 < 0 || v2 < 0; }
};

//////////////////////////////////////////////////////////////////////////////////////////
/// @brief Voronoi diagram structure - the dual of Delaunay triangulation
/// 
/// Each Voronoi cell contains all points closer to its site than to any other site.
/// - Voronoi vertices are circumcenters of Delaunay triangles
/// - Voronoi edges connect circumcenters of adjacent Delaunay triangles
/// - Convex hull sites have unbounded (infinite) cells
//////////////////////////////////////////////////////////////////////////////////////////
struct VoronoiDiagram {
	std::vector<Point2Cartesian> sites;     ///< Original input points
	std::vector<Point2Cartesian> vertices;  ///< Voronoi vertices (circumcenters)
	std::vector<VoronoiEdge> edges;         ///< Voronoi edges
	std::vector<std::vector<int>> cells;    ///< Vertex indices for each cell (CCW order)
	std::vector<bool> cellIsInfinite;       ///< True if cell extends to infinity

	int NumSites() const { return static_cast<int>(sites.size()); }
	int NumVertices() const { return static_cast<int>(vertices.size()); }
	int NumEdges() const { return static_cast<int>(edges.size()); }

	/// @brief Find the site (input point) nearest to a query point
	/// Uses linear search with scaled distances; rounded ties keep the first site.
	int FindNearestSite(const Point2Cartesian& p) const {
		MML::detail::ValidateGeometryCoordinates({p.X(), p.Y()}, "FindNearestSite");
		if (sites.empty()) return -1;
		
		int nearest = 0;
		auto minDistance = MML::detail::GeometryDistance2D(sites[0], p);
		
		for (int i = 1; i < NumSites(); i++) {
			const auto distance = MML::detail::GeometryDistance2D(sites[i], p);
			const bool closer = distance.fraction == 0 ? minDistance.fraction != 0
				: minDistance.fraction != 0 && (distance.exponent < minDistance.exponent
					|| (distance.exponent == minDistance.exponent && distance.fraction < minDistance.fraction));
			if (closer) {
				minDistance = distance;
				nearest = i;
			}
		}
		return nearest;
	}

	/// @brief Get the Voronoi cell as a polygon (for finite cells only)
	/// Returns empty polygon if cell is infinite or invalid
	Polygon2D GetCell(int siteIndex) const {
		if (siteIndex < 0 || siteIndex >= NumSites()) return Polygon2D();
		if (cellIsInfinite[siteIndex]) return Polygon2D();
		
		const auto& cellVertices = cells[siteIndex];
		if (cellVertices.size() < 3) return Polygon2D();
		
		std::vector<Point2Cartesian> pts;
		pts.reserve(cellVertices.size());
		for (int vi : cellVertices) {
			if (vi >= 0 && vi < NumVertices()) {
				pts.push_back(vertices[vi]);
			}
		}
		
		if (pts.size() < 3) return Polygon2D();
		return Polygon2D(pts);
	}

	/// @brief Get cell vertices in CCW order (may include -1 for infinite directions)
	const std::vector<int>& GetCellVertexIndices(int siteIndex) const {
		static const std::vector<int> empty;
		if (siteIndex < 0 || siteIndex >= NumSites()) return empty;
		return cells[siteIndex];
	}

	/// @brief Check if a cell is infinite (unbounded)
	bool IsCellInfinite(int siteIndex) const {
		if (siteIndex < 0 || siteIndex >= NumSites()) return true;
		return cellIsInfinite[siteIndex];
	}

	/// @brief Intersect a site's closed nearest-site halfplanes with a finite bounding box.
	/// Invalid site indices and zero-area intersections return an empty polygon.
	Polygon2D GetClippedCell(int siteIndex, Real minX, Real minY, Real maxX, Real maxY) const {
		MML::detail::ValidateGeometryCoordinates({minX, minY, maxX, maxY}, "GetClippedCell");
		if (minX > maxX || minY > maxY)
			throw GeometryError("GetClippedCell - bounding box endpoints must be ordered");
		if (siteIndex < 0 || siteIndex >= NumSites()) return Polygon2D();
		for (const auto& point : sites)
			MML::detail::ValidateGeometryCoordinates({point.X(), point.Y()}, "GetClippedCell sites");
		if (minX == maxX || minY == maxY) return Polygon2D();

		// Starting at the requested rectangle avoids an overflowing artificial outer box.
		std::vector<Point2Cartesian> points{
			{minX, minY}, {maxX, minY}, {maxX, maxY}, {minX, maxY}
		};
		const auto& site = sites[siteIndex];
		using B = MML::detail::GeometryBinary;
		for (int j = 0; j < NumSites(); j++) {
			if (j == siteIndex) continue;
			const auto& other = sites[j];
			const std::array<B, 2> normal{B::Difference(other.X(), site.X()), B::Difference(other.Y(), site.Y())};
			if (normal[0].fraction == 0 && normal[1].fraction == 0) continue;
			auto halfSquared = B::SumProducts(normal, normal);
			--halfSquared.exponent;
			points = ClipAgainstBisector(points, site, normal, halfSquared);
			if (points.size() < 3) return Polygon2D();
		}
		for (std::size_t i = 2; i < points.size(); ++i)
			if (MML::detail::GeometryOrientation2D::Sign(points[0], points[1], points[i]) != 0)
				return Polygon2D(points);
		return Polygon2D();
	}

private:
	static MML::detail::GeometryBinary BisectorDistance(const Point2Cartesian& point,
		const Point2Cartesian& site, const std::array<MML::detail::GeometryBinary, 2>& normal,
		const MML::detail::GeometryBinary& halfSquared) {
		using B = MML::detail::GeometryBinary;
		// |other-site|^2/2 - (point-site).(other-site): no rounded global midpoint.
		const std::array<B, 3> a{halfSquared, B::Difference(site.X(), point.X()), B::Difference(site.Y(), point.Y())};
		const std::array<B, 3> b{B(1), normal[0], normal[1]};
		const auto distance = B::SumProducts(a, b);
		B magnitude;
		for (std::size_t i = 0; i < a.size(); ++i) {
			auto term = B::Product(a[i], b[i]);
			term.fraction = std::abs(term.fraction);
			magnitude = B::Add(magnitude, term);
		}
		if (distance.fraction != 0 &&
			std::abs(B::Quotient(distance, magnitude).Value()) <= Precision::RoundoffZeroThreshold)
			return {};
		return distance;
	}

	static std::vector<Point2Cartesian> ClipAgainstBisector(const std::vector<Point2Cartesian>& points,
		const Point2Cartesian& site, const std::array<MML::detail::GeometryBinary, 2>& normal,
		const MML::detail::GeometryBinary& halfSquared) {
		using B = MML::detail::GeometryBinary;
		std::vector<Point2Cartesian> result;
		const auto append = [&](const Point2Cartesian& point) {
			if (result.empty() || point.X() != result.back().X() || point.Y() != result.back().Y())
				result.push_back(point);
		};
		for (std::size_t i = 0; i < points.size(); ++i) {
			const auto& first = points[i];
			const auto& second = points[(i + 1) % points.size()];
			auto d1 = BisectorDistance(first, site, normal, halfSquared);
			auto d2 = BisectorDistance(second, site, normal, halfSquared);
			const bool inside1 = d1.fraction >= 0, inside2 = d2.fraction >= 0;
			if (inside1) append(first);
			if (inside1 != inside2) {
				d1.fraction = std::abs(d1.fraction);
				d2.fraction = std::abs(d2.fraction);
				const auto denominator = B::Add(d1, d2);
				const auto interpolate = [&](Real a, Real b) {
					// The interpolation weight may underflow even when the vertex is representable.
					const auto numerator = B::SumProducts<2>({B(a), B(b)}, {d2, d1});
					const Real value = B::Quotient(numerator, denominator).Value();
					return std::clamp(value, std::min(a, b), std::max(a, b));
				};
				append({interpolate(first.X(), second.X()), interpolate(first.Y(), second.Y())});
			}
		}
		if (result.size() > 1 && result.front().X() == result.back().X() && result.front().Y() == result.back().Y())
			result.pop_back();
		return result;
	}
};

//////////////////////////////////////////////////////////////////////////////////////////
/// @brief Voronoi diagram computation
/// @details Constructs Voronoi diagram as the dual of Delaunay triangulation
//////////////////////////////////////////////////////////////////////////////////////////
class Voronoi {
public:
	/// @brief Compute Voronoi diagram from a set of points
	/// 
	/// The Voronoi diagram is constructed as the dual of the Delaunay triangulation:
	/// - Each Voronoi vertex is the circumcenter of a Delaunay triangle
	/// - Each Voronoi edge connects circumcenters of adjacent Delaunay triangles
	/// - Sites on the convex hull have infinite (unbounded) cells
	/// 
	/// @param points Input points (sites)
	/// @return VoronoiDiagram structure
	static VoronoiDiagram ComputeVoronoi(std::vector<Point2Cartesian> points) {
		if (points.size() > static_cast<std::size_t>(std::numeric_limits<int>::max() / 2))
			throw ArgumentError("ComputeVoronoi - too many points for integer indices");
		for (const auto& point : points)
			MML::detail::ValidateGeometryCoordinates({point.X(), point.Y()}, "ComputeVoronoi");
		VoronoiDiagram vd;
		vd.sites = points;
		
		int n = static_cast<int>(points.size());
		std::vector<int> order(n), unique, representative(n);
		std::iota(order.begin(), order.end(), 0);
		std::sort(order.begin(), order.end(), [&](int a, int b) {
			if (points[a].X() != points[b].X()) return points[a].X() < points[b].X();
			if (points[a].Y() != points[b].Y()) return points[a].Y() < points[b].Y();
			return a < b;
		});
		for (int index : order) {
			if (unique.empty() || points[index].X() != points[unique.back()].X() || points[index].Y() != points[unique.back()].Y())
				unique.push_back(index);
			representative[index] = unique.back();
		}
		bool collinear = true;
		for (std::size_t i = 2; i < unique.size(); ++i)
			if (RobustPredicates::Orientation2D(points[unique[0]], points[unique[1]], points[unique[i]]) != 0) {
				collinear = false;
				break;
			}
		vd.cells.resize(n);
		vd.cellIsInfinite.resize(n, collinear);
		if (collinear) {
			for (std::size_t i = 1; i < unique.size(); ++i) {
				const int first = unique[i - 1], second = unique[i];
				const auto direction = MML::detail::GeometryUnitDifference2D(points[first], points[second]);
				vd.edges.push_back({-1, -1, first, second, Vector2Cartesian(-direction.Y(), direction.X())});
			}
			return vd;
		}
		
		// Compute Delaunay triangulation
		auto dt = Triangulation::ComputeDelaunay(std::move(points));
		
		if (dt.NumTriangles() == 0) {
			throw GeometryError("ComputeVoronoi - incomplete Delaunay topology for noncollinear sites");
		}
		
		// Compute circumcenters (Voronoi vertices)
		vd.vertices.reserve(dt.NumTriangles());
		for (int i = 0; i < dt.NumTriangles(); i++) {
			const auto& tri = dt.triangles[i];
			const Point2Cartesian& a = dt.points[tri[0]];
			const Point2Cartesian& b = dt.points[tri[1]];
			const Point2Cartesian& c = dt.points[tri[2]];
			
			Point2Cartesian circumcenter = ComputeCircumcenter(a, b, c);
			vd.vertices.push_back(circumcenter);
		}
		
		// Initialize cells
		vd.cells.resize(n);
		vd.cellIsInfinite.resize(n, false);
		
		// Build map from site to triangles containing it
		std::vector<std::vector<int>> siteToTriangles(n);
		for (int t = 0; t < dt.NumTriangles(); t++) {
			for (int v : dt.triangles[t]) {
				if (v >= 0 && v < n) {
					siteToTriangles[v].push_back(t);
				}
			}
		}
		
		// Build Voronoi edges from Delaunay edges
		std::set<std::pair<int, int>> processedEdges;
		for (int t = 0; t < dt.NumTriangles(); t++) {
			for (int e = 0; e < 3; e++) {
				int v1 = dt.triangles[t][e];
				int v2 = dt.triangles[t][(e + 1) % 3];
				
				// Skip super-triangle vertices
				if (v1 >= n || v2 >= n) continue;
				
				auto edgeKey = std::make_pair(std::min(v1, v2), std::max(v1, v2));
				if (processedEdges.count(edgeKey)) continue;
				processedEdges.insert(edgeKey);
				
				int neighbor = dt.neighbors[t][e];
				
				VoronoiEdge vedge;
				vedge.site1 = v1;
				vedge.site2 = v2;
				vedge.v1 = t;  // Current triangle's circumcenter
				
				if (neighbor >= 0) {
					// Interior edge: connects two circumcenters
					vedge.v2 = neighbor;
				} else {
					// Boundary edge: extends to infinity
					vedge.v2 = -1;
					vd.cellIsInfinite[v1] = true;
					vd.cellIsInfinite[v2] = true;
					
					// Direction perpendicular to the Delaunay edge, pointing outward
					const Point2Cartesian& p1 = dt.points[v1];
					const Point2Cartesian& p2 = dt.points[v2];
					
					// The third vertex of the triangle
					int v3 = dt.triangles[t][(e + 2) % 3];
					const Point2Cartesian& p3 = dt.points[v3];
					
					// Direction perpendicular to edge
					const auto direction = MML::detail::GeometryUnitDifference2D(p1, p2);
					Vector2Cartesian perp(-direction.Y(), direction.X());
					if (RobustPredicates::Orientation2D(p1, p2, p3) > 0) {
						perp = Vector2Cartesian(-perp.X(), -perp.Y());
					}
					vedge.direction = perp;
				}
				
				vd.edges.push_back(vedge);
			}
		}
		
		// Build cell vertex lists (circumcenters around each site in CCW order)
		for (int s = 0; s < n; s++) {
			const auto& triangles = siteToTriangles[s];
			if (triangles.empty()) {
				vd.cellIsInfinite[s] = true;
				continue;
			}
			
			// Order triangles around the site
			std::vector<int> ordered;
			std::vector<bool> used(triangles.size(), false);
			
			// Start with first triangle
			ordered.push_back(triangles[0]);
			used[0] = true;
			
			// Try to build a chain
			for (size_t iter = 0; iter < triangles.size() - 1; iter++) {
				int lastTri = ordered.back();
				bool found = false;
				
				for (size_t i = 0; i < triangles.size(); i++) {
					if (used[i]) continue;
					int candTri = triangles[i];
					
					// Check if triangles share an edge (are neighbors)
					for (int edge = 0; edge < 3; edge++) {
						if (dt.neighbors[lastTri][edge] == candTri) {
							ordered.push_back(candTri);
							used[i] = true;
							found = true;
							break;
						}
					}
					if (found) break;
				}
				
				if (!found) {
					// Gap in the chain - cell is infinite
					vd.cellIsInfinite[s] = true;
					break;
				}
			}
			
			// Check if chain closes (first and last triangles are neighbors)
			if (ordered.size() == triangles.size() && ordered.size() >= 2) {
				int first = ordered.front();
				int last = ordered.back();
				bool closes = false;
				for (int edge = 0; edge < 3; edge++) {
					if (dt.neighbors[last][edge] == first) {
						closes = true;
						break;
					}
				}
				if (!closes) {
					vd.cellIsInfinite[s] = true;
				}
			}
			
			vd.cells[s] = ordered;
		}
		for (int s = 0; s < n; ++s) {
			if (representative[s] != s) {
				vd.cells[s] = vd.cells[representative[s]];
				vd.cellIsInfinite[s] = vd.cellIsInfinite[representative[s]];
			}
		}
		
		return vd;
	}

private:
	/// @brief Compute circumcenter of a triangle
	static Point2Cartesian ComputeCircumcenter(const Point2Cartesian& a, 
	                                           const Point2Cartesian& b, 
	                                           const Point2Cartesian& c) {
		Point2Cartesian center;
		if (!MML::detail::GeometryCircumcenter2D(a, b, c, center))
			throw GeometryError("ComputeVoronoi - degenerate Delaunay face has no circumcenter");
		return center;
	}
};

} // namespace MML::CompGeometry

#endif // MML_VORONOI_DIAGRAM_H

#ifndef MML_VORONOI_DIAGRAM_H
#define MML_VORONOI_DIAGRAM_H

/////////////////////////////////////////////////////////////////////////////////////////
///  @file      VoronoiDiagram.h
///  @brief     Voronoi Diagram computation (dual of Delaunay triangulation)
///  @details   Part of the MML Computational Geometry module
/////////////////////////////////////////////////////////////////////////////////////////

#include "mml/MMLBase.h"
#include "mml/base/Geometry/Geometry2D.h"
#include "mml/base/Vector/VectorN.h"
#include "mml/algorithms/CompGeometry/CompGeometryBase.h"
#include "mml/algorithms/CompGeometry/Triangulation.h"

#include <algorithm>
#include <cmath>
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
	/// Uses simple linear search - O(n)
	int FindNearestSite(const Point2Cartesian& p) const {
		if (sites.empty()) return -1;
		
		int nearest = 0;
		Real minDist2 = (sites[0].X() - p.X()) * (sites[0].X() - p.X()) + 
		                (sites[0].Y() - p.Y()) * (sites[0].Y() - p.Y());
		
		for (int i = 1; i < NumSites(); i++) {
			Real dist2 = (sites[i].X() - p.X()) * (sites[i].X() - p.X()) + 
			             (sites[i].Y() - p.Y()) * (sites[i].Y() - p.Y());
			if (dist2 < minDist2) {
				minDist2 = dist2;
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

	/// @brief Get a clipped version of cell for visualization
	/// Clips infinite cells to the specified bounding box
	Polygon2D GetClippedCell(int siteIndex, Real minX, Real minY, Real maxX, Real maxY) const {
		if (siteIndex < 0 || siteIndex >= NumSites()) return Polygon2D();
		
		// For finite cells, just return the cell (optionally clip to bounds)
		if (!cellIsInfinite[siteIndex]) {
			Polygon2D cell = GetCell(siteIndex);
			if (cell.NumVertices() > 0) {
				return ClipToRect(cell, minX, minY, maxX, maxY);
			}
			return Polygon2D();
		}
		
		// For infinite cells, we need to construct the cell with clipped boundaries
		Real expand = 2.0 * std::max(maxX - minX, maxY - minY);
		Real bigMinX = minX - expand;
		Real bigMaxX = maxX + expand;
		Real bigMinY = minY - expand;
		Real bigMaxY = maxY + expand;
		
		// Start with a large box centered on the site
		const Point2Cartesian& site = sites[siteIndex];
		std::vector<Point2Cartesian> cellPts = {
			Point2Cartesian(bigMinX, bigMinY),
			Point2Cartesian(bigMaxX, bigMinY),
			Point2Cartesian(bigMaxX, bigMaxY),
			Point2Cartesian(bigMinX, bigMaxY)
		};
		Polygon2D cell(cellPts);
		
		// Clip against perpendicular bisectors with all other sites
		for (int j = 0; j < NumSites(); j++) {
			if (j == siteIndex) continue;
			
			const Point2Cartesian& other = sites[j];
			// Midpoint between sites
			Real midX = (site.X() + other.X()) / 2.0;
			Real midY = (site.Y() + other.Y()) / 2.0;
			
			// Normal pointing from other toward site (we keep the side containing site)
			Real nx = site.X() - other.X();
			Real ny = site.Y() - other.Y();
			Real len = std::sqrt(nx * nx + ny * ny);
			if (len < 1e-10) continue;
			nx /= len;
			ny /= len;
			
			// Clip the cell to keep only points where (p - mid) · n >= 0
			std::vector<Point2Cartesian> clipped;
			int n = cell.NumVertices();
			for (int k = 0; k < n; k++) {
				const Point2Cartesian& p1 = cell[k];
				const Point2Cartesian& p2 = cell[(k + 1) % n];
				
				Real d1 = (p1.X() - midX) * nx + (p1.Y() - midY) * ny;
				Real d2 = (p2.X() - midX) * nx + (p2.Y() - midY) * ny;
				
				if (d1 >= -1e-10) {
					clipped.push_back(p1);
				}
				
				// Edge crosses the bisector line
				if ((d1 > 1e-10 && d2 < -1e-10) || (d1 < -1e-10 && d2 > 1e-10)) {
					Real t = d1 / (d1 - d2);
					Real ix = p1.X() + t * (p2.X() - p1.X());
					Real iy = p1.Y() + t * (p2.Y() - p1.Y());
					clipped.push_back(Point2Cartesian(ix, iy));
				}
			}
			
			if (clipped.size() < 3) {
				return Polygon2D();  // Cell was clipped away
			}
			cell = Polygon2D(clipped);
		}
		
		// Final clip to the requested bounding box
		return ClipToRect(cell, minX, minY, maxX, maxY);
	}

private:
	// Simple rectangle clipping helper
	static Polygon2D ClipToRect(const Polygon2D& poly, Real minX, Real minY, Real maxX, Real maxY) {
		std::vector<Point2Cartesian> pts = poly.Vertices();
		
		// Clip against left edge
		pts = ClipAgainstLine(pts, minX, 0, 1, 0);
		// Clip against right edge
		pts = ClipAgainstLine(pts, maxX, 0, -1, 0);
		// Clip against bottom edge
		pts = ClipAgainstLine(pts, 0, minY, 0, 1);
		// Clip against top edge
		pts = ClipAgainstLine(pts, 0, maxY, 0, -1);
		
		return pts.size() >= 3 ? Polygon2D(pts) : Polygon2D();
	}
	
	static std::vector<Point2Cartesian> ClipAgainstLine(const std::vector<Point2Cartesian>& pts,
														Real px, Real py, Real nx, Real ny) {
		std::vector<Point2Cartesian> result;
		int n = static_cast<int>(pts.size());
		if (n == 0) return result;
		
		for (int i = 0; i < n; i++) {
			const Point2Cartesian& p1 = pts[i];
			const Point2Cartesian& p2 = pts[(i + 1) % n];
			
			Real d1 = (p1.X() - px) * nx + (p1.Y() - py) * ny;
			Real d2 = (p2.X() - px) * nx + (p2.Y() - py) * ny;
			
			if (d1 >= -1e-10) {
				result.push_back(p1);
			}
			
			if ((d1 > 1e-10 && d2 < -1e-10) || (d1 < -1e-10 && d2 > 1e-10)) {
				Real t = d1 / (d1 - d2);
				Real ix = p1.X() + t * (p2.X() - p1.X());
				Real iy = p1.Y() + t * (p2.Y() - p1.Y());
				result.push_back(Point2Cartesian(ix, iy));
			}
		}
		return result;
	}
};

//////////////////////////////////////////////////////////////////////////////////////////
/// @brief Voronoi diagram computation
/// @details Constructs Voronoi diagram as the dual of Delaunay triangulation
//////////////////////////////////////////////////////////////////////////////////////////
class Voronoi {
	// Use centralized geometry epsilon from Constants
	static constexpr Real EPSILON = Constants::GEOMETRY_EPSILON;

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
		VoronoiDiagram vd;
		vd.sites = points;
		
		int n = static_cast<int>(points.size());
		if (n < 2) {
			// Degenerate cases
			vd.cells.resize(n);
			vd.cellIsInfinite.resize(n, true);
			return vd;
		}
		
		if (n == 2) {
			// Two points: single perpendicular bisector, both cells infinite
			vd.cells.resize(2);
			vd.cellIsInfinite = {true, true};
			
			// Add the edge between the two infinite cells
			VoronoiEdge edge;
			edge.v1 = -1;
			edge.v2 = -1;
			edge.site1 = 0;
			edge.site2 = 1;
			// Direction perpendicular to line between sites
			Real dx = points[1].X() - points[0].X();
			Real dy = points[1].Y() - points[0].Y();
			edge.direction = Vector2Cartesian(-dy, dx);
			vd.edges.push_back(edge);
			return vd;
		}
		
		// Compute Delaunay triangulation
		auto dt = Triangulation::ComputeDelaunay(std::move(points));
		
		if (dt.NumTriangles() == 0) {
			// All points collinear - each cell is infinite
			vd.cells.resize(n);
			vd.cellIsInfinite.resize(n, true);
			return vd;
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
					Real dx = p2.X() - p1.X();
					Real dy = p2.Y() - p1.Y();
					
					// The third vertex of the triangle
					int v3 = dt.triangles[t][(e + 2) % 3];
					const Point2Cartesian& p3 = dt.points[v3];
					
					// Direction perpendicular to edge
					Vector2Cartesian perp(-dy, dx);
					Real len = std::sqrt(perp.X() * perp.X() + perp.Y() * perp.Y());
					if (len > 1e-10) {
						perp = Vector2Cartesian(perp.X() / len, perp.Y() / len);
					}
					
					// Check if perp points away from the third vertex
					Real midX = (p1.X() + p2.X()) / 2.0;
					Real midY = (p1.Y() + p2.Y()) / 2.0;
					Real toThird = (p3.X() - midX) * perp.X() + (p3.Y() - midY) * perp.Y();
					if (toThird > 0) {
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
		
		return vd;
	}

private:
	/// @brief Compute circumcenter of a triangle
	static Point2Cartesian ComputeCircumcenter(const Point2Cartesian& a, 
	                                           const Point2Cartesian& b, 
	                                           const Point2Cartesian& c) {
		Real ax = a.X(), ay = a.Y();
		Real bx = b.X(), by = b.Y();
		Real cx = c.X(), cy = c.Y();
		
		Real d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by));
		if (std::abs(d) < 1e-14) {
			// Degenerate triangle - return centroid
			return Point2Cartesian((ax + bx + cx) / 3.0, (ay + by + cy) / 3.0);
		}
		
		Real ux = ((ax*ax + ay*ay) * (by - cy) + (bx*bx + by*by) * (cy - ay) + 
		           (cx*cx + cy*cy) * (ay - by)) / d;
		Real uy = ((ax*ax + ay*ay) * (cx - bx) + (bx*bx + by*by) * (ax - cx) + 
		           (cx*cx + cy*cy) * (bx - ax)) / d;
		
		return Point2Cartesian(ux, uy);
	}
};

} // namespace MML::CompGeometry

#endif // MML_VORONOI_DIAGRAM_H

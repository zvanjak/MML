///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Geometry2DPolygon.h                                                 ///
///  Description: 2D Polygon class with computational geometry operations            ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                        ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_GEOMETRY_2D_POLYGON_H
#define MML_GEOMETRY_2D_POLYGON_H

#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <vector>

#include "mml/MMLBase.h"
#include "mml/base/Vector/VectorTypes.h"
#include "mml/base/Geometry/Geometry.h"
#include "mml/base/Geometry/Geometry2DCore/Geometry2DLines.h"
#include "mml/base/Geometry/Geometry2DCore/Geometry2DTriangle.h"

namespace MML
{
	/// @brief 2D polygon defined by an ordered sequence of vertices
	/// @details Professional computational geometry polygon class supporting
	///          simple and convex polygons. Provides comprehensive operations
	///          including area (shoelace formula), centroid, bounding box,
	///          point containment (ray casting and winding number), convexity test,
	///          simplicity test, and basic triangularization.
	class Polygon2D
	{
	private:
		std::vector<Point2Cartesian> _vertices;

		/// @brief Helper for segment intersection test (proper intersection only)
		static bool SegmentsIntersect(const Point2Cartesian& p1, const Point2Cartesian& p2,
		                              const Point2Cartesian& p3, const Point2Cartesian& p4)
		{
			auto sign = [](Real val) { return (val > 0) - (val < 0); };
			
			// Compute cross products for orientation tests
			Real d1 = (p3.X() - p1.X()) * (p2.Y() - p1.Y()) - (p3.Y() - p1.Y()) * (p2.X() - p1.X());
			Real d2 = (p4.X() - p1.X()) * (p2.Y() - p1.Y()) - (p4.Y() - p1.Y()) * (p2.X() - p1.X());
			Real d3 = (p1.X() - p3.X()) * (p4.Y() - p3.Y()) - (p1.Y() - p3.Y()) * (p4.X() - p3.X());
			Real d4 = (p2.X() - p3.X()) * (p4.Y() - p3.Y()) - (p2.Y() - p3.Y()) * (p4.X() - p3.X());
			
			// Check if segments straddle each other
			if (sign(d1) * sign(d2) < 0 && sign(d3) * sign(d4) < 0)
				return true;
			
			return false;
		}

	public:
		/// @brief Default constructor - creates empty polygon
		Polygon2D() {}
		
		/// @brief Constructs polygon from vector of vertices
		/// @param vertices The vertices in order
		Polygon2D(const std::vector<Point2Cartesian>& vertices) : _vertices(vertices) {}
		
		/// @brief Constructs polygon from initializer list
		/// @param list The vertices in order
		Polygon2D(std::initializer_list<Point2Cartesian> list) : _vertices(list) {}

		/// @brief Gets the number of vertices
		int NumVertices() const { return static_cast<int>(_vertices.size()); }
		
		/// @brief Gets const reference to vertices vector
		const std::vector<Point2Cartesian>& Vertices() const { return _vertices; }
		/// @brief Gets mutable reference to vertices vector
		std::vector<Point2Cartesian>& Vertices() { return _vertices; }
		
		/// @brief Index operator for vertex access (const)
		const Point2Cartesian& operator[](int i) const { return _vertices[i]; }
		/// @brief Index operator for vertex access (mutable)
		Point2Cartesian& operator[](int i) { return _vertices[i]; }
		
		/// @brief Gets vertex by index with bounds checking
		/// @throws VectorAccessBoundsError if index out of range
		const Point2Cartesian& Vertex(int i) const { 
			if (i < 0 || i >= NumVertices())
				throw VectorAccessBoundsError("Polygon2D::Vertex - index out of range", i, NumVertices());
			return _vertices[i]; 
		}
		
		/// @brief Adds a vertex to the polygon
		void AddVertex(const Point2Cartesian& p) { _vertices.push_back(p); }
		
		/// @brief Removes all vertices
		void Clear() { _vertices.clear(); }

		/// @brief Gets edge as a segment
		/// @param i Edge index (0 to NumVertices-1)
		/// @throws VectorAccessBoundsError if index out of range
		SegmentLine2D Edge(int i) const {
			if (i < 0 || i >= NumVertices())
				throw VectorAccessBoundsError("Polygon2D::Edge - index out of range", i, NumVertices());
			return SegmentLine2D(_vertices[i], _vertices[(i + 1) % NumVertices()]);
		}

		/// @brief Calculates the perimeter
		Real Perimeter() const {
			if (NumVertices() < 2) return 0.0;
			
			Real perimeter = 0.0;
			int n = NumVertices();
			for (int i = 0; i < n; i++) {
				perimeter += _vertices[i].Dist(_vertices[(i + 1) % n]);
			}
			return perimeter;
		}

		/// @brief Calculates the signed area using shoelace formula
		/// @return Positive for CCW, negative for CW orientation
		Real SignedArea() const {
			if (NumVertices() < 3) return 0.0;
			
			Real area = 0.0;
			int n = NumVertices();
			for (int i = 0; i < n; i++) {
				area += _vertices[i].X() * _vertices[(i + 1) % n].Y();
				area -= _vertices[i].Y() * _vertices[(i + 1) % n].X();
			}
			return area / 2.0;
		}

		/// @brief Calculates the unsigned area
		Real Area() const {
			return std::abs(SignedArea());
		}

		/// @brief Calculates the centroid (geometric center)
		/// @throws GeometryError if polygon is empty or degenerate
		Point2Cartesian Centroid() const {
			if (NumVertices() == 0)
				throw GeometryError("Polygon2D::Centroid - empty polygon");
			
			if (NumVertices() == 1)
				return _vertices[0];
			
			if (NumVertices() == 2)
				return Point2Cartesian((_vertices[0].X() + _vertices[1].X()) / 2.0,
				                      (_vertices[0].Y() + _vertices[1].Y()) / 2.0);
			
			// For polygons with 3+ vertices, use area-weighted formula
			Real cx = 0.0, cy = 0.0;
			Real signedArea = 0.0;
			int n = NumVertices();
			
			for (int i = 0; i < n; i++) {
				int j = (i + 1) % n;
				Real cross = _vertices[i].X() * _vertices[j].Y() - _vertices[j].X() * _vertices[i].Y();
				signedArea += cross;
				cx += (_vertices[i].X() + _vertices[j].X()) * cross;
				cy += (_vertices[i].Y() + _vertices[j].Y()) * cross;
			}
			
			signedArea /= 2.0;
			if (std::abs(signedArea) < Constants::GEOMETRY_EPSILON)
				throw GeometryError("Polygon2D::Centroid - degenerate polygon (zero area)");
			
			cx /= (6.0 * signedArea);
			cy /= (6.0 * signedArea);
			
			return Point2Cartesian(cx, cy);
		}

		/// @brief Axis-aligned bounding box structure
		struct BoundingBox {
			Real minX, maxX, minY, maxY;
			
			/// @brief Gets the width of the bounding box
			Real Width() const { return maxX - minX; }
			/// @brief Gets the height of the bounding box
			Real Height() const { return maxY - minY; }
			/// @brief Gets the center of the bounding box
			Point2Cartesian Center() const { return Point2Cartesian((minX + maxX) / 2.0, (minY + maxY) / 2.0); }
		};
		
		/// @brief Computes the axis-aligned bounding box
		/// @throws GeometryError if polygon is empty
		BoundingBox GetBoundingBox() const {
			if (NumVertices() == 0)
				throw GeometryError("Polygon2D::GetBoundingBox - empty polygon");
			
			BoundingBox box;
			box.minX = box.maxX = _vertices[0].X();
			box.minY = box.maxY = _vertices[0].Y();
			
			for (int i = 1; i < NumVertices(); i++) {
				box.minX = std::min(box.minX, _vertices[i].X());
				box.maxX = std::max(box.maxX, _vertices[i].X());
				box.minY = std::min(box.minY, _vertices[i].Y());
				box.maxY = std::max(box.maxY, _vertices[i].Y());
			}
			
			return box;
		}

		/// @brief Tests if vertices are ordered counter-clockwise
		bool IsCounterClockwise() const {
			return SignedArea() > 0.0;
		}
		
		/// @brief Tests if vertices are ordered clockwise
		bool IsClockwise() const {
			return SignedArea() < 0.0;
		}

		/// @brief Reverses the vertex order
		void Reverse() {
			std::reverse(_vertices.begin(), _vertices.end());
		}

		/// @brief Tests if the polygon is simple (no self-intersections)
		bool IsSimple() const {
			int n = NumVertices();
			if (n < 3) return false;
			
			// Check all non-adjacent edge pairs for intersections
			for (int i = 0; i < n; i++) {
				Point2Cartesian p1 = _vertices[i];
				Point2Cartesian p2 = _vertices[(i + 1) % n];
				
				// Check against non-adjacent edges
				for (int j = i + 2; j < n; j++) {
					// Skip if edges are adjacent (share a vertex)
					if (j == (i + n - 1) % n) continue;
					
					Point2Cartesian p3 = _vertices[j];
					Point2Cartesian p4 = _vertices[(j + 1) % n];
					
					if (SegmentsIntersect(p1, p2, p3, p4))
						return false;
				}
			}
			
			return true;
		}

		/// @brief Tests if the polygon is convex
		bool IsConvex() const {
			int n = NumVertices();
			if (n < 3) return false;
			
			bool hasPositive = false;
			bool hasNegative = false;
			
			for (int i = 0; i < n; i++) {
				Point2Cartesian p0 = _vertices[i];
				Point2Cartesian p1 = _vertices[(i + 1) % n];
				Point2Cartesian p2 = _vertices[(i + 2) % n];
				
				// Compute cross product of edges (p0->p1) and (p1->p2)
				Vector2Cartesian v1(p0, p1);
				Vector2Cartesian v2(p1, p2);
				Real cross = v1.X() * v2.Y() - v1.Y() * v2.X();
				
				if (cross > Constants::GEOMETRY_EPSILON) hasPositive = true;
				if (cross < -Constants::GEOMETRY_EPSILON) hasNegative = true;
				
				// If we have both signs, polygon is not convex
				if (hasPositive && hasNegative)
					return false;
			}
			
			return true;
		}

		/// @brief Tests if a point is inside the polygon using ray-casting
		/// @details Works for simple polygons (both convex and non-convex)
		/// @param point The point to test
		/// @return True if point is inside the polygon
		bool Contains(const Point2Cartesian& point) const {
			int n = NumVertices();
			if (n < 3) return false;
			
			// Ray casting: count intersections of horizontal ray from point to the right
			int intersections = 0;
			
			for (int i = 0; i < n; i++) {
				Point2Cartesian p1 = _vertices[i];
				Point2Cartesian p2 = _vertices[(i + 1) % n];
				
				// Check if ray intersects edge
				if ((p1.Y() > point.Y()) != (p2.Y() > point.Y())) {
					// Compute x-coordinate of intersection
					Real xIntersect = p1.X() + (point.Y() - p1.Y()) * (p2.X() - p1.X()) / (p2.Y() - p1.Y());
					
					if (point.X() < xIntersect)
						intersections++;
				}
			}
			
			// Odd number of intersections = inside
			return (intersections % 2) == 1;
		}

		/// @brief Computes the winding number for a point
		/// @details More robust than ray-casting for complex polygons
		/// @param point The point to test
		/// @return Winding number (non-zero means inside)
		int WindingNumber(const Point2Cartesian& point) const {
			int n = NumVertices();
			if (n < 3) return 0;
			
			int wn = 0; // winding number counter
			
			for (int i = 0; i < n; i++) {
				Point2Cartesian p1 = _vertices[i];
				Point2Cartesian p2 = _vertices[(i + 1) % n];
				
				if (p1.Y() <= point.Y()) {
					if (p2.Y() > point.Y()) {
						// Upward crossing
						Real cross = (p2.X() - p1.X()) * (point.Y() - p1.Y()) - (point.X() - p1.X()) * (p2.Y() - p1.Y());
						if (cross > 0) {
							wn++; // Point is left of edge
						}
					}
				} else {
					if (p2.Y() <= point.Y()) {
						// Downward crossing
						Real cross = (p2.X() - p1.X()) * (point.Y() - p1.Y()) - (point.X() - p1.X()) * (p2.Y() - p1.Y());
						if (cross < 0) {
							wn--; // Point is right of edge
						}
					}
				}
			}
			
			return wn;
		}

		/// @brief Decomposes polygon into triangles
		/// @details Uses fan triangulation (works correctly for convex polygons)
		/// @return Vector of triangles covering the polygon
		std::vector<Triangle2D> Triangularization() const {
			std::vector<Triangle2D> triangles;
			int n = NumVertices();
			
			if (n < 3) return triangles;
			if (n == 3) {
				triangles.push_back(Triangle2D(_vertices[0], _vertices[1], _vertices[2]));
				return triangles;
			}
			
			// Simple fan triangulation from first vertex (works for convex polygons)
			// For non-convex, this is a placeholder - full ear clipping would be more complex
			for (int i = 1; i < n - 1; i++) {
				triangles.push_back(Triangle2D(_vertices[0], _vertices[i], _vertices[i + 1]));
			}
			
			return triangles;
		}
	};

} // namespace MML

#endif // MML_GEOMETRY_2D_POLYGON_H

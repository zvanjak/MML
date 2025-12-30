///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Geometry2D.h                                                        ///
///  Description: 2D geometry classes (Triangle2D, Polygon2D, Circle, Ellipse)        ///
///               Intersection, area, perimeter calculations                          ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_GEOMETRY_2D_H
#define MML_GEOMETRY_2D_H

#include "MMLBase.h"

#include "base/VectorTypes.h"
#include "base/Geometry.h"

namespace MML
{
	// Geometry epsilon for floating-point comparisons
	namespace {
		constexpr Real GEOMETRY_EPSILON = 1e-10;
		
		inline bool AreEqual(Real a, Real b, Real eps = GEOMETRY_EPSILON) {
			return std::abs(a - b) < eps;
		}
	}
	class Line2D
	{
	private:
		Point2Cartesian _point;
		Vector2Cartesian _direction;					// unit vector in line direction

	public:
		// Constructors
		Line2D(const Point2Cartesian& pnt, const Vector2Cartesian& dir)
		{
			_point = pnt;
			_direction = dir.GetAsUnitVector();
		}

		Line2D(const Point2Cartesian& a, const Point2Cartesian& b)
		{
			Vector2Cartesian dir(a, b);
			_point = a;
			_direction = dir.GetAsUnitVector();
		}

		// Accessors
		Point2Cartesian   StartPoint() const { return _point; }
		Point2Cartesian&  StartPoint() { return _point; }

		Vector2Cartesian  Direction() const { return _direction; }
		Vector2Cartesian& Direction() { return _direction; }

		// Point on line at parameter t
		Point2Cartesian operator()(Real t) const
		{
			Vector2Cartesian dist = t * _direction;
			Point2Cartesian ret = _point + dist;
			return ret;
		}

		// Distance from point to line
		Real DistanceToPoint(const Point2Cartesian& p) const
		{
			Vector2Cartesian v(_point, p);
			// Distance = |v × direction| (2D cross product magnitude)
			Real cross = v.X() * _direction.Y() - v.Y() * _direction.X();
			return std::abs(cross);
		}

		// Closest point on line to given point
		Point2Cartesian ClosestPoint(const Point2Cartesian& p) const
		{
			Vector2Cartesian v(_point, p);
			// Project v onto direction
			Real t = v.X() * _direction.X() + v.Y() * _direction.Y();
			return (*this)(t);
		}

		// Create perpendicular line through given point
		Line2D Perpendicular(const Point2Cartesian& throughPoint) const
		{
			// Perpendicular direction in 2D: rotate 90 degrees
			Vector2Cartesian perpDir(-_direction.Y(), _direction.X());
			return Line2D(throughPoint, perpDir);
		}

		// Check if point lies on line (within epsilon tolerance)
		bool Contains(const Point2Cartesian& p, Real epsilon = GEOMETRY_EPSILON) const
		{
			return DistanceToPoint(p) < epsilon;
		}

		// Intersection with another line (returns false if parallel)
		bool Intersects(const Line2D& other, Point2Cartesian* intersection = nullptr) const
		{
			// Lines: P1 + t*D1 and P2 + s*D2
			// Solve: P1 + t*D1 = P2 + s*D2
			Real det = _direction.X() * other._direction.Y() - _direction.Y() * other._direction.X();
			
			if (std::abs(det) < GEOMETRY_EPSILON)
				return false; // Parallel or coincident
			
			if (intersection) {
				Vector2Cartesian dp(other._point, _point);
				Real t = (dp.X() * other._direction.Y() - dp.Y() * other._direction.X()) / det;
				*intersection = (*this)(t);
			}
			
			return true;
		}
	};

	class SegmentLine2D
	{
	private:
		Point2Cartesian _point1;
		Point2Cartesian _point2;

	public:
		// Constructors
		SegmentLine2D(const Point2Cartesian& pnt1, const Point2Cartesian& pnt2) : _point1(pnt1), _point2(pnt2)
		{ }

		SegmentLine2D(const Point2Cartesian& pnt1, const Vector2Cartesian& direction, Real t) : _point1(pnt1)
		{
			_point2 = pnt1 + direction * t;
		}

		// Accessors
		Point2Cartesian  StartPoint() const { return _point1; }
		Point2Cartesian& StartPoint() { return _point1; }

		Point2Cartesian  EndPoint()  const { return _point2; }
		Point2Cartesian& EndPoint() { return _point2; }

		// Point on segment at parameter t ∈ [0,1]
		Point2Cartesian PointOnSegment(Real t) const
		{
			if (t < 0.0 || t > 1.0)
				throw GeometryError("SegmentLine2D::PointOnSegment t must be in [0,1]");

			Vector2Cartesian dist = t * Direction();
			Point2Cartesian ret = _point1 + dist;
			return ret;
		}

		// Basic properties
		Real                Length()    const { return _point1.Dist(_point2); }
		Vector2Cartesian    Direction() const { return Vector2Cartesian(_point1, _point2); }
		
		// Midpoint of segment
		Point2Cartesian Midpoint() const
		{
			return Point2Cartesian((_point1.X() + _point2.X()) / 2.0,
			                      (_point1.Y() + _point2.Y()) / 2.0);
		}

		// Distance from point to segment
		Real DistanceToPoint(const Point2Cartesian& p) const
		{
			Vector2Cartesian v(_point1, _point2);
			Vector2Cartesian w(_point1, p);
			
			Real lenSq = v.NormL2() * v.NormL2();
			if (lenSq < GEOMETRY_EPSILON)
				return _point1.Dist(p); // Degenerate segment
			
			// Project w onto v
			Real t = (w.X() * v.X() + w.Y() * v.Y()) / lenSq;
			t = std::clamp(t, REAL(0.0), REAL(1.0)); // Clamp to segment
			
			Point2Cartesian projection = _point1 + v * t;
			return p.Dist(projection);
		}

		// Closest point on segment to given point
		Point2Cartesian ClosestPoint(const Point2Cartesian& p) const
		{
			Vector2Cartesian v(_point1, _point2);
			Vector2Cartesian w(_point1, p);
			
			Real lenSq = v.NormL2() * v.NormL2();
			if (lenSq < GEOMETRY_EPSILON)
				return _point1; // Degenerate segment
			
			Real t = (w.X() * v.X() + w.Y() * v.Y()) / lenSq;
			t = std::clamp(t, REAL(0.0), REAL(1.0));
			
			return _point1 + v * t;
		}

		// Segment-segment intersection test
		bool Intersects(const SegmentLine2D& other, Point2Cartesian* intersection = nullptr) const
		{
			Vector2Cartesian r = Direction();
			Vector2Cartesian s = other.Direction();
			Vector2Cartesian pq(_point1, other._point1);
			
			Real rxs = r.X() * s.Y() - r.Y() * s.X();
			Real pqxr = pq.X() * r.Y() - pq.Y() * r.X();
			
			// Check if parallel
			if (std::abs(rxs) < GEOMETRY_EPSILON) {
				return false; // Parallel or collinear (collinear overlap not handled)
			}
			
			// Compute parameters
			Real t = (pq.X() * s.Y() - pq.Y() * s.X()) / rxs;
			Real u = pqxr / rxs;
			
			// Check if intersection is within both segments
			if (t >= 0.0 && t <= 1.0 && u >= 0.0 && u <= 1.0) {
				if (intersection) {
					*intersection = _point1 + r * t;
				}
				return true;
			}
			
			return false;
		}
	};

	// Triangle2D - Professional 2D triangle class with cached side lengths
	class Triangle2D
	{
	private:
		Point2Cartesian _pnt1, _pnt2, _pnt3;
		mutable Real _a, _b, _c; // Cached side lengths
		mutable bool _sidesComputed;

		void ComputeSides() const
		{
			if (!_sidesComputed) {
				_a = _pnt1.Dist(_pnt2);
				_b = _pnt2.Dist(_pnt3);
				_c = _pnt3.Dist(_pnt1);
				_sidesComputed = true;
			}
		}

	public:
		// Constructor
		Triangle2D(const Point2Cartesian& pnt1, const Point2Cartesian& pnt2, const Point2Cartesian& pnt3) 
			: _pnt1(pnt1), _pnt2(pnt2), _pnt3(pnt3), _a(0), _b(0), _c(0), _sidesComputed(false)
		{ }

		// Vertex accessors
		Point2Cartesian  Pnt1() const { return _pnt1; }
		Point2Cartesian& Pnt1() { _sidesComputed = false; return _pnt1; }
		Point2Cartesian  Pnt2() const { return _pnt2; }
		Point2Cartesian& Pnt2() { _sidesComputed = false; return _pnt2; }
		Point2Cartesian  Pnt3() const { return _pnt3; }
		Point2Cartesian& Pnt3() { _sidesComputed = false; return _pnt3; }

		// Side lengths (cached for efficiency)
		Real A() const { ComputeSides(); return _a; }
		Real B() const { ComputeSides(); return _b; }
		Real C() const { ComputeSides(); return _c; }

		// Perimeter
		Real Perimeter() const { return A() + B() + C(); }

		// Area using Heron's formula
		Real Area() const
		{
			Real a = A(), b = B(), c = C();
			Real s = (a + b + c) / 2.0;
			return std::sqrt(s * (s - a) * (s - b) * (s - c));
		}

		// Signed area (positive = CCW, negative = CW)
		Real SignedArea() const
		{
			Real area = (_pnt2.X() - _pnt1.X()) * (_pnt3.Y() - _pnt1.Y()) -
			           (_pnt3.X() - _pnt1.X()) * (_pnt2.Y() - _pnt1.Y());
			return area / 2.0;
		}

		// Centroid (geometric center)
		Point2Cartesian Centroid() const
		{
			return Point2Cartesian((_pnt1.X() + _pnt2.X() + _pnt3.X()) / 3.0,
			                      (_pnt1.Y() + _pnt2.Y() + _pnt3.Y()) / 3.0);
		}

		// Circumcenter (center of circumscribed circle)
		Point2Cartesian Circumcenter() const
		{
			Real D = 2.0 * (_pnt1.X() * (_pnt2.Y() - _pnt3.Y()) +
			               _pnt2.X() * (_pnt3.Y() - _pnt1.Y()) +
			               _pnt3.X() * (_pnt1.Y() - _pnt2.Y()));
			
			if (std::abs(D) < GEOMETRY_EPSILON)
				throw GeometryError("Triangle2D::Circumcenter - degenerate triangle");
			
			Real p1Sq = _pnt1.X() * _pnt1.X() + _pnt1.Y() * _pnt1.Y();
			Real p2Sq = _pnt2.X() * _pnt2.X() + _pnt2.Y() * _pnt2.Y();
			Real p3Sq = _pnt3.X() * _pnt3.X() + _pnt3.Y() * _pnt3.Y();
			
			Real ux = (p1Sq * (_pnt2.Y() - _pnt3.Y()) +
			          p2Sq * (_pnt3.Y() - _pnt1.Y()) +
			          p3Sq * (_pnt1.Y() - _pnt2.Y())) / D;
			
			Real uy = (p1Sq * (_pnt3.X() - _pnt2.X()) +
			          p2Sq * (_pnt1.X() - _pnt3.X()) +
			          p3Sq * (_pnt2.X() - _pnt1.X())) / D;
			
			return Point2Cartesian(ux, uy);
		}

		// Circumradius (radius of circumscribed circle)
		Real CircumRadius() const
		{
			Real a = A(), b = B(), c = C();
			Real area = Area();
			
			if (area < GEOMETRY_EPSILON)
				throw GeometryError("Triangle2D::CircumRadius - degenerate triangle");
			
			return (a * b * c) / (4.0 * area);
		}

		// Incenter (center of inscribed circle)
		Point2Cartesian Incenter() const
		{
			Real a = A(), b = B(), c = C();
			Real perimeter = a + b + c;
			
			if (perimeter < GEOMETRY_EPSILON)
				throw GeometryError("Triangle2D::Incenter - degenerate triangle");
			
			Real x = (a * _pnt3.X() + b * _pnt1.X() + c * _pnt2.X()) / perimeter;
			Real y = (a * _pnt3.Y() + b * _pnt1.Y() + c * _pnt2.Y()) / perimeter;
			
			return Point2Cartesian(x, y);
		}

		// Inradius (radius of inscribed circle)
		Real InRadius() const
		{
			Real area = Area();
			Real s = Perimeter() / 2.0;
			
			if (s < GEOMETRY_EPSILON)
				throw GeometryError("Triangle2D::InRadius - degenerate triangle");
			
			return area / s;
		}

		// Angles in radians
		Real AngleA() const // Angle at vertex 3 (pnt3) - opposite to side a
		{
			Real a = A(), b = B(), c = C();
			// Law of cosines: a² = b² + c² - 2bc·cos(A)
			Real cosA = (b * b + c * c - a * a) / (2.0 * b * c);
			return std::acos(std::clamp(cosA, REAL(-1.0), REAL(1.0)));
		}

		Real AngleB() const // Angle at vertex 1 (pnt1) - opposite to side b
		{
			Real a = A(), b = B(), c = C();
			Real cosB = (a * a + c * c - b * b) / (2.0 * a * c);
			return std::acos(std::clamp(cosB, REAL(-1.0), REAL(1.0)));
		}

		Real AngleC() const // Angle at vertex 2 (pnt2) - opposite to side c
		{
			Real a = A(), b = B(), c = C();
			Real cosC = (a * a + b * b - c * c) / (2.0 * a * b);
			return std::acos(std::clamp(cosC, REAL(-1.0), REAL(1.0)));
		}

		// Point containment using barycentric coordinates
		bool Contains(const Point2Cartesian& p, Real epsilon = GEOMETRY_EPSILON) const
		{
			// Compute barycentric coordinates
			Real denom = (_pnt2.Y() - _pnt3.Y()) * (_pnt1.X() - _pnt3.X()) +
			            (_pnt3.X() - _pnt2.X()) * (_pnt1.Y() - _pnt3.Y());
			
			if (std::abs(denom) < epsilon)
				return false; // Degenerate triangle
			
			Real a = ((_pnt2.Y() - _pnt3.Y()) * (p.X() - _pnt3.X()) +
			         (_pnt3.X() - _pnt2.X()) * (p.Y() - _pnt3.Y())) / denom;
			
			Real b = ((_pnt3.Y() - _pnt1.Y()) * (p.X() - _pnt3.X()) +
			         (_pnt1.X() - _pnt3.X()) * (p.Y() - _pnt3.Y())) / denom;
			
			Real c = 1.0 - a - b;
			
			// Point is inside if all barycentric coordinates are non-negative
			return (a >= -epsilon) && (b >= -epsilon) && (c >= -epsilon);
		}

		// Triangle classification with epsilon tolerance
		bool IsRight(Real epsilon = GEOMETRY_EPSILON) const
		{
			Real a = A(), b = B(), c = C();
			Real aSq = a * a, bSq = b * b, cSq = c * c;
			
			return AreEqual(aSq + bSq, cSq, epsilon) ||
			       AreEqual(aSq + cSq, bSq, epsilon) ||
			       AreEqual(bSq + cSq, aSq, epsilon);
		}

		bool IsIsosceles(Real epsilon = GEOMETRY_EPSILON) const
		{
			Real a = A(), b = B(), c = C();
			return AreEqual(a, b, epsilon) || AreEqual(a, c, epsilon) || AreEqual(b, c, epsilon);
		}

		bool IsEquilateral(Real epsilon = GEOMETRY_EPSILON) const
		{
			Real a = A(), b = B(), c = C();
			return AreEqual(a, b, epsilon) && AreEqual(a, c, epsilon);
		}

		// Orientation
		bool IsCounterClockwise() const
		{
			return SignedArea() > 0.0;
		}

		bool IsClockwise() const
		{
			return SignedArea() < 0.0;
		}
	};

	// Polygon2D - Professional computational geometry polygon class
	// Represents a 2D polygon as an ordered sequence of vertices
	class Polygon2D
	{
	private:
		std::vector<Point2Cartesian> _vertices;

		// Helper: Check if two line segments intersect (proper intersection, not endpoint touching)
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
		// Constructors
		Polygon2D() {}
		
		Polygon2D(const std::vector<Point2Cartesian>& vertices) : _vertices(vertices) {}
		
		Polygon2D(std::initializer_list<Point2Cartesian> list) : _vertices(list) {}

		// Vertex access
		int NumVertices() const { return static_cast<int>(_vertices.size()); }
		
		const std::vector<Point2Cartesian>& Vertices() const { return _vertices; }
		std::vector<Point2Cartesian>& Vertices() { return _vertices; }
		
		const Point2Cartesian& operator[](int i) const { return _vertices[i]; }
		Point2Cartesian& operator[](int i) { return _vertices[i]; }
		
		const Point2Cartesian& Vertex(int i) const { 
			if (i < 0 || i >= NumVertices())
				throw VectorAccessBoundsError("Polygon2D::Vertex - index out of range", i, NumVertices());
			return _vertices[i]; 
		}
		
		void AddVertex(const Point2Cartesian& p) { _vertices.push_back(p); }
		
		void Clear() { _vertices.clear(); }

		// Edge access - returns directed edge as segment
		SegmentLine2D Edge(int i) const {
			if (i < 0 || i >= NumVertices())
				throw VectorAccessBoundsError("Polygon2D::Edge - index out of range", i, NumVertices());
			return SegmentLine2D(_vertices[i], _vertices[(i + 1) % NumVertices()]);
		}

		// Basic properties
		Real Perimeter() const {
			if (NumVertices() < 2) return 0.0;
			
			Real perimeter = 0.0;
			int n = NumVertices();
			for (int i = 0; i < n; i++) {
				perimeter += _vertices[i].Dist(_vertices[(i + 1) % n]);
			}
			return perimeter;
		}

		// Signed area using shoelace formula (positive = counterclockwise, negative = clockwise)
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

		// Unsigned area
		Real Area() const {
			return std::abs(SignedArea());
		}

		// Centroid (geometric center)
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
			if (std::abs(signedArea) < GEOMETRY_EPSILON)
				throw GeometryError("Polygon2D::Centroid - degenerate polygon (zero area)");
			
			cx /= (6.0 * signedArea);
			cy /= (6.0 * signedArea);
			
			return Point2Cartesian(cx, cy);
		}

		// Bounding box
		struct BoundingBox {
			Real minX, maxX, minY, maxY;
			
			Real Width() const { return maxX - minX; }
			Real Height() const { return maxY - minY; }
			Point2Cartesian Center() const { return Point2Cartesian((minX + maxX) / 2.0, (minY + maxY) / 2.0); }
		};
		
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

		// Orientation: counterclockwise (CCW) or clockwise (CW)
		bool IsCounterClockwise() const {
			return SignedArea() > 0.0;
		}
		
		bool IsClockwise() const {
			return SignedArea() < 0.0;
		}

		// Reverse vertex order
		void Reverse() {
			std::reverse(_vertices.begin(), _vertices.end());
		}

		// Simple polygon check - no self-intersections (excluding consecutive edges)
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

		// Convex polygon check
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
				
				if (cross > GEOMETRY_EPSILON) hasPositive = true;
				if (cross < -GEOMETRY_EPSILON) hasNegative = true;
				
				// If we have both signs, polygon is not convex
				if (hasPositive && hasNegative)
					return false;
			}
			
			return true;
		}

		// Point-in-polygon test using ray-casting algorithm
		// Works for simple polygons (both convex and non-convex)
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

		// Winding number algorithm for point-in-polygon (more robust than ray-casting)
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

		// Legacy compatibility
		std::vector<Point2Cartesian> Points() const { return _vertices; }
		std::vector<Point2Cartesian>& Points() { return _vertices; }
		bool IsInside(const Point2Cartesian& pnt) const { return Contains(pnt); }

		// Triangularization (ear clipping algorithm for simple polygons)
		// Note: This is a basic implementation - for production use, consider more robust algorithms
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
}
#endif

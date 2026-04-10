///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Geometry2DLines.h                                                   ///
///  Description: 2D line and segment classes with intersection support              ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                        ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_GEOMETRY_2D_LINES_H
#define MML_GEOMETRY_2D_LINES_H

#include <algorithm>
#include <cmath>
#include <vector>

#include "mml/MMLBase.h"
#include "mml/base/Vector/VectorTypes.h"
#include "mml/base/Geometry/Geometry.h"

namespace MML
{
	/// @brief Helper function for floating-point comparisons with tolerance
	/// @param a First value
	/// @param b Second value
	/// @param eps Tolerance (default: GEOMETRY_EPSILON)
	/// @return True if |a - b| < eps
	inline bool AreEqual(Real a, Real b, Real eps = Constants::GEOMETRY_EPSILON) {
		return std::abs(a - b) < eps;
	}

	// ============================================================================
	// 2D Line Intersection Types
	// ============================================================================

	/// @brief Classification of 2D line-line intersection
	enum class LineIntersectionType2D {
		Point,		///< Lines intersect at a single point
		Parallel,	///< Lines are parallel (no intersection)
		Coincident	///< Lines are the same (infinite intersections)
	};

	/// @brief Result of 2D line-line intersection
	struct LineIntersection2D {
		LineIntersectionType2D type = LineIntersectionType2D::Parallel;
		Point2Cartesian point;	///< Intersection point (valid only for Point type)
		Real t1 = 0.0;			///< Parameter on first line: point = p1 + t1 * d1
		Real t2 = 0.0;			///< Parameter on second line: point = p2 + t2 * d2

		/// @brief Check if lines intersect at a single point
		bool HasIntersection() const { return type == LineIntersectionType2D::Point; }
		/// @brief Check if intersection is a single point
		bool IsPoint() const { return type == LineIntersectionType2D::Point; }
		/// @brief Check if lines are parallel (no intersection)
		bool IsParallel() const { return type == LineIntersectionType2D::Parallel; }
		/// @brief Check if lines are coincident (infinite intersections)
		bool IsCoincident() const { return type == LineIntersectionType2D::Coincident; }
	};

	// ============================================================================
	// 2D Segment Intersection Types
	// ============================================================================

	/// @brief Classification of segment-segment intersection
	enum class SegmentIntersectionType {
		None,	   ///< No intersection
		Point,	   ///< Single point intersection
		Collinear, ///< Segments are collinear (may or may not overlap)
		Overlap	   ///< Collinear and overlapping
	};

	/// @brief Result of segment-segment intersection
	struct SegmentIntersection {
		SegmentIntersectionType type = SegmentIntersectionType::None;
		Point2Cartesian point;		  ///< Intersection point (for Point type)
		Point2Cartesian overlapStart; ///< Start of overlap (for Overlap type)
		Point2Cartesian overlapEnd;	  ///< End of overlap (for Overlap type)

		/// @brief Check if segments intersect (point or overlap)
		bool HasIntersection() const { 
			return type == SegmentIntersectionType::Point || type == SegmentIntersectionType::Overlap; 
		}
		/// @brief Check if intersection is a single point
		bool IsPointIntersection() const { return type == SegmentIntersectionType::Point; }
		/// @brief Check if segments overlap
		bool IsOverlap() const { return type == SegmentIntersectionType::Overlap; }
	};

	/// @brief 2D infinite line defined by a point and direction vector
	/// @details Represents an infinite line in 2D space using parametric form: P(t) = point + t * direction.
	///          The direction is automatically normalized to a unit vector.
	///          Provides methods for point distance, closest point, perpendicular lines, and line-line intersection.
	class Line2D
	{
	private:
		Point2Cartesian _point;
		Vector2Cartesian _direction;					// unit vector in line direction

	public:
		/// @brief Constructs a line from a point and direction vector
		/// @param pnt A point on the line
		/// @param dir Direction vector (will be normalized)
		Line2D(const Point2Cartesian& pnt, const Vector2Cartesian& dir)
		{
			_point = pnt;
			_direction = dir.GetAsUnitVector();
		}

		/// @brief Constructs a line through two points
		/// @param a First point on the line
		/// @param b Second point on the line
		Line2D(const Point2Cartesian& a, const Point2Cartesian& b)
		{
			Vector2Cartesian dir(a, b);
			_point = a;
			_direction = dir.GetAsUnitVector();
		}

		/// @brief Gets the start point (const version)
		Point2Cartesian   StartPoint() const { return _point; }
		/// @brief Gets the start point (mutable reference)
		Point2Cartesian&  StartPoint() { return _point; }

		/// @brief Gets the unit direction vector (const version)
		Vector2Cartesian  Direction() const { return _direction; }
		/// @brief Gets the unit direction vector (mutable reference)
		Vector2Cartesian& Direction() { return _direction; }

		/// @brief Gets point on line at parameter t
		/// @param t Parameter value (signed distance along direction)
		/// @return Point at StartPoint + t * Direction
		Point2Cartesian operator()(Real t) const
		{
			Vector2Cartesian dist = t * _direction;
			Point2Cartesian ret = _point + dist;
			return ret;
		}

		/// @brief Calculates perpendicular distance from point to line
		/// @param p The point
		/// @return Shortest distance from p to the line
		Real DistanceToPoint(const Point2Cartesian& p) const
		{
			Vector2Cartesian v(_point, p);
			// Distance = |v × direction| (2D cross product magnitude)
			Real cross = v.X() * _direction.Y() - v.Y() * _direction.X();
			return std::abs(cross);
		}

		/// @brief Finds the closest point on the line to a given point
		/// @param p The point
		/// @return Point on line closest to p
		Point2Cartesian ClosestPoint(const Point2Cartesian& p) const
		{
			Vector2Cartesian v(_point, p);
			// Project v onto direction
			Real t = v.X() * _direction.X() + v.Y() * _direction.Y();
			return (*this)(t);
		}

		/// @brief Creates a perpendicular line through a given point
		/// @param throughPoint The point through which the perpendicular passes
		/// @return Line perpendicular to this one, passing through throughPoint
		Line2D Perpendicular(const Point2Cartesian& throughPoint) const
		{
			// Perpendicular direction in 2D: rotate 90 degrees
			Vector2Cartesian perpDir(-_direction.Y(), _direction.X());
			return Line2D(throughPoint, perpDir);
		}

		/// @brief Checks if a point lies on the line
		/// @param p The point to test
		/// @param epsilon Tolerance for the test
		/// @return True if point is within epsilon of the line
		bool Contains(const Point2Cartesian& p, Real epsilon = PrecisionValues<Real>::GeometryEpsilon) const
		{
			return DistanceToPoint(p) < epsilon;
		}

		/// @brief Computes full intersection information with another line
		/// @param other The other line
		/// @return LineIntersection2D with type, point, and parameters
		/// @details Handles all cases: intersection point, parallel, coincident
		LineIntersection2D Intersection(const Line2D& other) const
		{
			LineIntersection2D result;

			Real cross = _direction.X() * other._direction.Y() - _direction.Y() * other._direction.X();
			Real dx = other._point.X() - _point.X();
			Real dy = other._point.Y() - _point.Y();

			if (std::abs(cross) < PrecisionValues<Real>::GeometryEpsilon) {
				// Lines are parallel - check if coincident
				Real crossP = dx * _direction.Y() - dy * _direction.X();
				if (std::abs(crossP) < PrecisionValues<Real>::GeometryEpsilon)
					result.type = LineIntersectionType2D::Coincident;
				else
					result.type = LineIntersectionType2D::Parallel;
				return result;
			}

			Real t = (dx * other._direction.Y() - dy * other._direction.X()) / cross;
			result.type = LineIntersectionType2D::Point;
			result.point = (*this)(t);
			result.t1 = t;
			result.t2 = (dx * _direction.Y() - dy * _direction.X()) / cross;
			return result;
		}

		/// @brief Tests for intersection with another line (convenience method)
		/// @param other The other line
		/// @param intersection Optional output parameter for intersection point
		/// @return True if lines intersect at a single point, false if parallel/coincident
		bool Intersects(const Line2D& other, Point2Cartesian* intersection = nullptr) const
		{
			auto result = Intersection(other);
			if (intersection && result.HasIntersection())
				*intersection = result.point;
			return result.HasIntersection();
		}
	};

	/// @brief 2D line segment defined by two endpoints
	/// @details Represents a finite line segment in 2D space between two points.
	///          Provides methods for point-on-segment access, distance calculations,
	///          closest point computation, and segment-segment intersection tests.
	class SegmentLine2D
	{
	private:
		Point2Cartesian _point1;
		Point2Cartesian _point2;

	public:
		/// @brief Constructs a segment from two endpoints
		/// @param pnt1 Start point
		/// @param pnt2 End point
		SegmentLine2D(const Point2Cartesian& pnt1, const Point2Cartesian& pnt2) : _point1(pnt1), _point2(pnt2)
		{ }

		/// @brief Constructs a segment from a start point, direction, and length
		/// @param pnt1 Start point
		/// @param direction Direction vector (not necessarily unit)
		/// @param t Length parameter (endpoint = pnt1 + direction * t)
		SegmentLine2D(const Point2Cartesian& pnt1, const Vector2Cartesian& direction, Real t) : _point1(pnt1)
		{
			_point2 = pnt1 + direction * t;
		}

		/// @brief Gets the start point (const version)
		Point2Cartesian  StartPoint() const { return _point1; }
		/// @brief Gets the start point (mutable reference)
		Point2Cartesian& StartPoint() { return _point1; }

		/// @brief Gets the end point (const version)
		Point2Cartesian  EndPoint()  const { return _point2; }
		/// @brief Gets the end point (mutable reference)
		Point2Cartesian& EndPoint() { return _point2; }

		/// @brief Gets a point on the segment at parameter t ∈ [0,1]
		/// @param t Parameter value (0 = start, 1 = end)
		/// @return Point on segment at the specified parameter
		/// @throws GeometryError if t is outside [0,1]
		Point2Cartesian PointOnSegment(Real t) const
		{
			if (t < 0.0 || t > 1.0)
				throw GeometryError("SegmentLine2D::PointOnSegment t must be in [0,1]");

			Vector2Cartesian dist = t * Direction();
			Point2Cartesian ret = _point1 + dist;
			return ret;
		}

		/// @brief Gets the length of the segment
		Real                Length()    const { return _point1.Dist(_point2); }
		/// @brief Gets the direction vector from start to end (not normalized)
		Vector2Cartesian    Direction() const { return Vector2Cartesian(_point1, _point2); }
		
		/// @brief Gets the midpoint of the segment
		Point2Cartesian Midpoint() const
		{
			return Point2Cartesian((_point1.X() + _point2.X()) / 2.0,
			                      (_point1.Y() + _point2.Y()) / 2.0);
		}

		/// @brief Calculates the distance from a point to the segment
		/// @param p The point
		/// @return Shortest distance from p to any point on the segment
		Real DistanceToPoint(const Point2Cartesian& p) const
		{
			Vector2Cartesian v(_point1, _point2);
			Vector2Cartesian w(_point1, p);
			
			Real lenSq = v.NormL2() * v.NormL2();
			if (lenSq < Constants::GEOMETRY_EPSILON)
				return _point1.Dist(p); // Degenerate segment
			
			// Project w onto v
			Real t = (w.X() * v.X() + w.Y() * v.Y()) / lenSq;
			t = std::clamp(t, REAL(0.0), REAL(1.0)); // Clamp to segment
			
			Point2Cartesian projection = _point1 + v * t;
			return p.Dist(projection);
		}

		/// @brief Finds the closest point on the segment to a given point
		/// @param p The point
		/// @return Point on segment closest to p
		Point2Cartesian ClosestPoint(const Point2Cartesian& p) const
		{
			Vector2Cartesian v(_point1, _point2);
			Vector2Cartesian w(_point1, p);
			
			Real lenSq = v.NormL2() * v.NormL2();
			if (lenSq < Constants::GEOMETRY_EPSILON)
				return _point1; // Degenerate segment
			
			Real t = (w.X() * v.X() + w.Y() * v.Y()) / lenSq;
			t = std::clamp(t, REAL(0.0), REAL(1.0));
			
			return _point1 + v * t;
		}

	private:
		/// @brief Cross product of vectors OA and OB (internal helper)
		static Real Cross(const Point2Cartesian& O, const Point2Cartesian& A, const Point2Cartesian& B) {
			return (A.X() - O.X()) * (B.Y() - O.Y()) - (A.Y() - O.Y()) * (B.X() - O.X());
		}

		/// @brief Check if point q lies on segment pr (assuming collinearity)
		static bool OnSegment(const Point2Cartesian& p, const Point2Cartesian& q, const Point2Cartesian& r) {
			return q.X() <= std::max(p.X(), r.X()) + Constants::GEOMETRY_EPSILON && 
			       q.X() >= std::min(p.X(), r.X()) - Constants::GEOMETRY_EPSILON &&
			       q.Y() <= std::max(p.Y(), r.Y()) + Constants::GEOMETRY_EPSILON && 
			       q.Y() >= std::min(p.Y(), r.Y()) - Constants::GEOMETRY_EPSILON;
		}

	public:
		/// @brief Computes full intersection information with another segment
		/// @param other The other segment
		/// @return SegmentIntersection with type, point, or overlap endpoints
		/// @details Handles all cases: crossing, endpoint touching, collinear, overlapping
		SegmentIntersection Intersection(const SegmentLine2D& other) const
		{
			SegmentIntersection result;
			constexpr Real EPSILON = Constants::GEOMETRY_EPSILON;

			const Point2Cartesian& p1 = _point1;
			const Point2Cartesian& q1 = _point2;
			const Point2Cartesian& p2 = other._point1;
			const Point2Cartesian& q2 = other._point2;

			Real d1 = Cross(p2, q2, p1);
			Real d2 = Cross(p2, q2, q1);
			Real d3 = Cross(p1, q1, p2);
			Real d4 = Cross(p1, q1, q2);

			// General case: segments cross
			if (((d1 > EPSILON && d2 < -EPSILON) || (d1 < -EPSILON && d2 > EPSILON)) &&
				((d3 > EPSILON && d4 < -EPSILON) || (d3 < -EPSILON && d4 > EPSILON))) {
				Real dx1 = q1.X() - p1.X();
				Real dy1 = q1.Y() - p1.Y();
				Real dx2 = q2.X() - p2.X();
				Real dy2 = q2.Y() - p2.Y();
				Real dx12 = p2.X() - p1.X();
				Real dy12 = p2.Y() - p1.Y();

				Real denom = dx1 * dy2 - dy1 * dx2;
				Real t = (dx12 * dy2 - dy12 * dx2) / denom;

				result.type = SegmentIntersectionType::Point;
				result.point = Point2Cartesian(p1.X() + t * dx1, p1.Y() + t * dy1);
				return result;
			}

			// Check for collinearity
			bool collinear = (std::abs(d1) <= EPSILON && std::abs(d2) <= EPSILON);

			if (collinear) {
				bool p2_on_seg1 = OnSegment(p1, p2, q1);
				bool q2_on_seg1 = OnSegment(p1, q2, q1);
				bool p1_on_seg2 = OnSegment(p2, p1, q2);
				bool q1_on_seg2 = OnSegment(p2, q1, q2);

				if (p2_on_seg1 || q2_on_seg1 || p1_on_seg2 || q1_on_seg2) {
					std::vector<Point2Cartesian> candidates;
					if (p2_on_seg1) candidates.push_back(p2);
					if (q2_on_seg1) candidates.push_back(q2);
					if (p1_on_seg2) candidates.push_back(p1);
					if (q1_on_seg2) candidates.push_back(q1);

					if (candidates.size() >= 2) {
						std::sort(candidates.begin(), candidates.end(), 
							[&p1, &q1, EPSILON](const Point2Cartesian& a, const Point2Cartesian& b) {
								Real ta = (std::abs(q1.X() - p1.X()) > EPSILON) 
									? (a.X() - p1.X()) / (q1.X() - p1.X())
									: (a.Y() - p1.Y()) / (q1.Y() - p1.Y());
								Real tb = (std::abs(q1.X() - p1.X()) > EPSILON) 
									? (b.X() - p1.X()) / (q1.X() - p1.X())
									: (b.Y() - p1.Y()) / (q1.Y() - p1.Y());
								return ta < tb;
							});

						result.overlapStart = candidates.front();
						result.overlapEnd = candidates.back();

						if (result.overlapStart.Dist(result.overlapEnd) < EPSILON) {
							result.type = SegmentIntersectionType::Point;
							result.point = result.overlapStart;
						} else {
							result.type = SegmentIntersectionType::Overlap;
						}
						return result;
					}
				}
				result.type = SegmentIntersectionType::Collinear;
				return result;
			}

			// Check for endpoint touches
			if (std::abs(d1) <= EPSILON && OnSegment(p2, p1, q2)) {
				result.type = SegmentIntersectionType::Point;
				result.point = p1;
				return result;
			}
			if (std::abs(d2) <= EPSILON && OnSegment(p2, q1, q2)) {
				result.type = SegmentIntersectionType::Point;
				result.point = q1;
				return result;
			}
			if (std::abs(d3) <= EPSILON && OnSegment(p1, p2, q1)) {
				result.type = SegmentIntersectionType::Point;
				result.point = p2;
				return result;
			}
			if (std::abs(d4) <= EPSILON && OnSegment(p1, q2, q1)) {
				result.type = SegmentIntersectionType::Point;
				result.point = q2;
				return result;
			}

			return result; // No intersection
		}

		/// @brief Tests for intersection with another segment (convenience method)
		/// @param other The other segment
		/// @param intersection Optional output parameter for intersection point
		/// @return True if segments intersect (point or overlap)
		bool Intersects(const SegmentLine2D& other, Point2Cartesian* intersection = nullptr) const
		{
			auto result = Intersection(other);
			if (intersection && result.IsPointIntersection())
				*intersection = result.point;
			return result.HasIntersection();
		}
	};

} // namespace MML

#endif // MML_GEOMETRY_2D_LINES_H

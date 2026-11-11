///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Geometry3DLines.h                                                   ///
///  Description: 3D line primitives and intersection types                           ///
///               Line3D (infinite line), SegmentLine3D (finite segment)              ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_GEOMETRY_3D_LINES_H
#define MML_GEOMETRY_3D_LINES_H

#include <mml/MMLBase.h>
#include <mml/base/Vector/VectorN.h>
#include <mml/base/Vector/VectorTypes3D.h>
#include <mml/base/Geometry/Geometry.h>
#include <mml/base/BaseUtils.h>

namespace MML
{

	// ============================================================================
	// 3D Line Intersection Types
	// ============================================================================

	/// @brief Classification of 3D line-line intersection
	enum class LineIntersectionType3D {
		Point,		///< Lines intersect at a single point
		Parallel,	///< Lines are parallel (no intersection)
		Coincident,	///< Lines are the same (infinite intersections)
		Skew		///< Lines are skew (not parallel but don't intersect)
	};

	/// @brief Result of 3D line-line intersection
	struct LineIntersection3D {
		LineIntersectionType3D type = LineIntersectionType3D::Parallel;
		Point3Cartesian point1;	///< Closest point on first line (also intersection point if intersecting)
		Point3Cartesian point2;	///< Closest point on second line
		Real t1 = 0.0;			///< Parameter on first line
		Real t2 = 0.0;			///< Parameter on second line
		Real distance = 0.0;	///< Distance between closest points (0 if intersecting)

		/// @brief Check if lines intersect at a single point
		bool HasIntersection() const { return type == LineIntersectionType3D::Point; }
		/// @brief Check if intersection is a single point
		bool IsIntersecting() const { return type == LineIntersectionType3D::Point; }
		/// @brief Check if lines are skew (closest approach but no intersection)
		bool AreSkew() const { return type == LineIntersectionType3D::Skew; }
		/// @brief Check if lines are skew
		bool IsSkew() const { return type == LineIntersectionType3D::Skew; }
		/// @brief Check if lines are parallel
		bool IsParallel() const { return type == LineIntersectionType3D::Parallel; }
		/// @brief Check if lines are coincident
		bool IsCoincident() const { return type == LineIntersectionType3D::Coincident; }
		
		/// @brief Get intersection point (alias for point1, valid when IsIntersecting)
		const Point3Cartesian& Point() const { return point1; }
	};

	// ============================================================================
	// 3D Line Primitives
	// ============================================================================

	/// @brief Infinite line in 3D space
	/// @details Represented by a point on the line and a unit direction vector.
	///          Provides operations for point distance, line-line distance,
	///          intersection detection, and perpendicular line construction.
	class Line3D
	{
	private:
		Pnt3Cart  _point;
		Vec3Cart _direction;

		Real DistanceGeometry(const Line3D& line, Vec3Cart& normal, Real& normalLength) const
		{
			const Vec3Cart delta(_point, line._point);
			const Real length1 = std::hypot(_direction.X(), _direction.Y(), _direction.Z());
			const Real length2 = std::hypot(line._direction.X(), line._direction.Y(), line._direction.Z());
			if (!std::isfinite(delta.X()) || !std::isfinite(delta.Y()) || !std::isfinite(delta.Z()) ||
				!std::isfinite(length1) || !std::isfinite(length2) || length1 == 0 || length2 == 0)
				throw GeometryError("Line3D distance - finite points and nonzero finite directions required");

			Real distance;
			if (IsParallel(line)) {
				normalLength = 0;
				const Vec3Cart perpendicular = VectorProduct(delta, _direction / length1);
				distance = std::hypot(perpendicular.X(), perpendicular.Y(), perpendicular.Z());
			} else {
				normal = VectorProduct(_direction, line._direction);
				normalLength = std::hypot(normal.X(), normal.Y(), normal.Z());
				if (!std::isfinite(normalLength) || normalLength == 0)
					throw GeometryError("Line3D distance - nonparallel direction cross product is not representable");
				normal = normal / normalLength;
				distance = Abs(ScalarProduct(normal, delta));
			}
			if (!std::isfinite(distance))
				throw GeometryError("Line3D distance - separation is not representable");
			return distance;
		}

	public:
		/// @brief Default constructor
		Line3D() {}
		
		/// @brief Constructs a line from point and direction
		/// @param pnt A point on the line
		/// @param dir Direction vector (will be normalized)
		/// @throws GeometryError if direction is null vector
		Line3D(const Pnt3Cart& pnt, const Vec3Cart dir)
		{
			// check for null vector as direction
			if (dir.X() == 0.0 && dir.Y() == 0.0 && dir.Z() == 0.0)
				throw GeometryError("Line3D ctor - null vector as direction");

			_point = pnt;
			_direction = dir.Normalized();
		}
		
		/// @brief Constructs a line through two points
		/// @param a First point
		/// @param b Second point
		/// @throws GeometryError if points are identical
		Line3D(const Pnt3Cart& a, const Pnt3Cart& b)
		{
			// check for same points
			if (a == b)
				throw GeometryError("Line3D ctor - same points");

			Vec3Cart dir(a, b);
			_point = a;
			_direction = dir.Normalized();
		}

		/// @brief Gets the reference point on the line (const)
		Pnt3Cart   StartPoint() const { return _point; }
		/// @brief Gets the reference point on the line (mutable)
		Pnt3Cart& StartPoint() { return _point; }

		/// @brief Gets the unit direction vector (const)
		Vec3Cart  Direction() const { return _direction; }
		/// @brief Gets the unit direction vector (mutable)
		Vec3Cart& Direction() { return _direction; }

		/// @brief Gets a point on the line at parameter t
		/// @param t Parameter value (point = startPoint + t * direction)
		Pnt3Cart operator()(Real t) const { return _point + t * _direction; }

		/// @brief Tests if two lines are equal (same geometric line, not just parallel)
		/// @param b The other line
		/// @param eps Tolerance for comparison
		/// @details Two lines are equal if they share the same points (each start point lies
		///          on the other line) AND their directions are either equal or anti-parallel
		///          (opposite directions represent the same geometric line).
		bool AreEqual(const Line3D& b, Real eps = Precision::Line3DAreEqualTolerance) const
		{
			if (!IsPointOnLine(b.StartPoint()) || !b.IsPointOnLine(StartPoint()))
				return false;
			
			// Check if directions are equal OR anti-parallel (same line, opposite direction)
			Vec3Cart negDir = -b.Direction();
			return Direction().IsEqualTo(b.Direction(), eps) || Direction().IsEqualTo(negDir, eps);
		}
		
		/// @brief Equality operator
		bool operator==(const Line3D& b) const
		{
			return AreEqual(b, Precision::Line3DAreEqualTolerance);
		}

		/// @brief Tests if this line is perpendicular to another line
		bool IsPerpendicular(const Line3D& b, Real eps = Precision::Line3DIsPerpendicularTolerance) const
		{
			// |cos(angle)| < eps - without Abs, any negative dot (e.g. anti-parallel) passed
			return Abs(ScalarProduct(Direction(), b.Direction())) < eps;
		}
		
		/// @brief Tests if this line is parallel to another line
		/// @details Two lines are parallel if their directions are either equal or anti-parallel
		bool IsParallel(const Line3D& b, Real eps = Precision::Line3DIsParallelTolerance) const
		{
			Vec3Cart negDir = -b.Direction();
			return Direction().IsEqualTo(b.Direction(), eps) || Direction().IsEqualTo(negDir, eps);
		}

		/// @brief Tests if a point lies on this line
		/// @param pnt The point to test
		/// @param eps Tolerance for the test
		bool IsPointOnLine(const Pnt3Cart& pnt, Real eps = Precision::Line3DIsPointOnLineTolerance) const
		{
			if (pnt == StartPoint())
				return true;	// point is start point of line

			// check if point is on line by checking if the vector from start point to point is parallel OR anti-parallel to direction vector
			// A line extends infinitely in both directions, so we use cross product which is zero for both parallel and anti-parallel vectors
			Vec3Cart vecFromStartToPnt(_point, pnt);
			Vec3Cart crossProd = VectorProduct(vecFromStartToPnt, _direction);
			return crossProd.NormL2() < eps;
		}

		/// @brief Calculates the distance from a point to this line
		/// @param pnt The point
		/// @return Perpendicular distance from point to line
		Real Dist(const Pnt3Cart& pnt) const
		{
			// Bronshtein 3.394
			const Real a = pnt.X();
			const Real b = pnt.Y();
			const Real c = pnt.Z();

			const Real x1 = StartPoint().X();
			const Real y1 = StartPoint().Y();
			const Real z1 = StartPoint().Z();

			const Real l = Direction().X();
			const Real m = Direction().Y();
			const Real n = Direction().Z();

			Real numer = POW2((a - x1) * m - (b - y1) * l) + POW2((b - y1) * n - (c - z1) * m) + POW2((c - z1) * l - (a - x1) * n);
			Real denom = l * l + m * m + n * n;

			return sqrt(numer / denom);
		}
	
		/// @brief Finds the nearest point on this line to a given point
		/// @param pnt The external point
		/// @return Point on line closest to pnt
		Pnt3Cart NearestPointOnLine(const Pnt3Cart& pnt) const
		{
			// https://math.stackexchange.com/questions/1521128/given-a-line-and-a-point-in-3d-how-to-find-the-closest-point-on-the-line         
			Vec3Cart line_dir = this->Direction();
			Vec3Cart rad_vec_AP(StartPoint(), pnt);

			Real t = rad_vec_AP * line_dir / POW2(line_dir.NormL2());

			return StartPoint() + t * line_dir;
		}

		/// @brief Calculates the distance between two lines
		/// @param line The other line
		/// @return Shortest distance, using the IsParallel tolerance for the parallel approximation
		Real Dist(const Line3D& line) const
		{
			Vec3Cart normal;
			Real normalLength;
			return DistanceGeometry(line, normal, normalLength);
		}
		
		/// @brief Calculates distance between two lines with nearest points
		/// @param line The other line
		/// @param out_dist Output: shortest distance between lines
		/// @param out_line1_pnt Output: nearest point on this line
		/// @param out_line2_pnt Output: nearest point on other line
		/// @return false if IsParallel is true (nearest-point outputs are left unchanged)
		bool Dist(const Line3D& line, Real& out_dist, Pnt3Cart& out_line1_pnt, Pnt3Cart& out_line2_pnt) const
		{
			if (IsParallel(line)) {
				out_dist = Dist(line);
				return false;
			}
			const auto result = Intersection(line);
			out_dist = result.distance;
			out_line1_pnt = result.point1;
			out_line2_pnt = result.point2;
			return true;
		}
		
		/// @brief Computes full intersection information with another line
		/// @param line The other line
		/// @return LineIntersection3D with type, closest points, parameters, and distance
		/// @details Parallel classification uses IsParallel. Its representative pair fixes
		///          the second reference point and projects it onto the first line.
		/// @throws GeometryError if inputs or required results are not representable
		LineIntersection3D Intersection(const Line3D& line) const
		{
			LineIntersection3D result;
			Vec3Cart normal;
			Real normalLength;
			result.distance = DistanceGeometry(line, normal, normalLength);
			const Vec3Cart delta(_point, line._point);
			if (normalLength == 0) {
				const Real length = std::hypot(_direction.X(), _direction.Y(), _direction.Z());
				result.t1 = ScalarProduct(_direction / length, delta) / length;
				result.type = result.distance < Precision::Line3DIsPointOnLineTolerance
					? LineIntersectionType3D::Coincident : LineIntersectionType3D::Parallel;
			} else {
				// Skew-line cross-product solve avoids cancellation in a*c-b*b.
				// https://en.wikipedia.org/wiki/Skew_lines#Nearest_points
				result.t1 = ScalarProduct(VectorProduct(line._direction, normal), delta) / normalLength;
				result.t2 = ScalarProduct(VectorProduct(_direction, normal), delta) / normalLength;
				result.type = result.distance < Precision::Line3DIntersectionTolerance
					? LineIntersectionType3D::Point : LineIntersectionType3D::Skew;
			}
			if (!std::isfinite(result.t1) || !std::isfinite(result.t2))
				throw GeometryError("Line3D intersection - closest-point parameters are not representable");
			result.point1 = (*this)(result.t1);
			result.point2 = line(result.t2);
			if (!std::isfinite(result.point1.X()) || !std::isfinite(result.point1.Y()) || !std::isfinite(result.point1.Z()) ||
				!std::isfinite(result.point2.X()) || !std::isfinite(result.point2.Y()) || !std::isfinite(result.point2.Z()))
				throw GeometryError("Line3D intersection - closest points are not representable");
			return result;
		}

		/// @brief Finds the intersection point of two lines (convenience method)
		/// @param line The other line
		/// @param out_inter_pnt Output: intersection point if exists
		/// @return true if lines intersect, false if parallel or skew
		bool Intersects(const Line3D& line, Pnt3Cart& out_inter_pnt) const
		{
			auto result = Intersection(line);
			if (result.HasIntersection()) {
				out_inter_pnt = result.point1;
				return true;
			}
			if (result.IsCoincident()) {
				out_inter_pnt = StartPoint();
				return true;
			}
			return false;
		}

		/// @brief Creates a perpendicular line through a given point
		/// @param pnt The point (must not be on this line)
		/// @return Line perpendicular to this line, passing through pnt
		/// @throws GeometryError if point is on the line
		Line3D PerpendicularLineThroughPoint(const Pnt3Cart& pnt) const
		{
			if (IsPointOnLine(pnt))
				throw GeometryError("Line3D::PerpendicularLineThroughPoint - point is on the line");

			return Line3D(pnt, Vec3Cart(pnt, NearestPointOnLine(pnt)).Normalized());
		}
	};

	/// @brief Line segment in 3D space
	/// @details Represents a finite line segment defined by two endpoints.
	///          Provides length, direction, and point-to-segment distance calculations.
	class SegmentLine3D
	{
	private:
		Pnt3Cart _point1;
		Pnt3Cart _point2;

	public:
		/// @brief Constructs a segment from two points
		/// @param pnt1 Start point
		/// @param pnt2 End point
		SegmentLine3D(Pnt3Cart pnt1, Pnt3Cart pnt2) : _point1(pnt1), _point2(pnt2)
		{
		}
		
		/// @brief Constructs a segment from point, direction, and length
		/// @param pnt1 Start point
		/// @param direction Direction vector
		/// @param t Length along direction
		SegmentLine3D(Pnt3Cart pnt1, Vec3Cart direction, Real t)
		{
			_point1 = pnt1;
			_point2 = pnt1 + t * direction;
		}

		/// @brief Gets the start point (const)
		Pnt3Cart   StartPoint() const { return _point1; }
		/// @brief Gets the start point (mutable)
		Pnt3Cart& StartPoint() { return _point1; }

		/// @brief Gets the end point (const)
		Pnt3Cart   EndPoint() const { return _point2; }
		/// @brief Gets the end point (mutable)
		Pnt3Cart& EndPoint() { return _point2; }

		/// @brief Gets a point on the segment at parameter t
		/// @param t Parameter in [0, 1]
		/// @throws GeometryError if t is outside [0, 1]
		Pnt3Cart		operator()(Real t) const
		{
			if (t < 0.0 || t > 1.0)
				throw GeometryError("SegmentLine3D::PointOnSegment - t not in [0,1]");

			return _point1 + t * Direction();
		}

		/// @brief Gets the length of the segment
		Real              Length()    const { return _point1.Dist(_point2); }
		/// @brief Gets the direction vector from start to end
		Vec3Cart  Direction() const { return Vec3Cart(_point1, _point2); }

		/// @brief Calculates the distance from a point to this segment
		/// @param pnt The point
		/// @return Shortest distance from point to segment
		Real Dist(const Pnt3Cart& pnt) const
		{
			Vec3Cart segDir = Direction();
			Real segLength = Length();
			
			if (segLength == 0.0)
				return _point1.Dist(pnt);  // degenerate segment
			
			// Project point onto line containing segment
			Vec3Cart vecToPnt(_point1, pnt);
			Real t = (vecToPnt * segDir) / (segLength * segLength);
			
			// Clamp t to [0, 1] to stay on segment
			if (t < 0.0)
				return _point1.Dist(pnt);  // closest to start point
			else if (t > 1.0)
				return _point2.Dist(pnt);  // closest to end point
			else
			{
				// Closest point is on the segment
				Pnt3Cart closestPoint = _point1 + t * segDir;
				return closestPoint.Dist(pnt);
			}
		}
	};

} // namespace MML

#endif // MML_GEOMETRY_3D_LINES_H

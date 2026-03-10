///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Geometry2DCircleBox.h                                               ///
///  Description: 2D Circle and Box (AABB) classes                                    ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                        ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_GEOMETRY_2D_CIRCLE_BOX_H
#define MML_GEOMETRY_2D_CIRCLE_BOX_H

#include <algorithm>
#include <cmath>
#include <vector>

#include "mml/MMLBase.h"
#include "mml/base/Vector/VectorTypes.h"
#include "mml/base/Geometry/Geometry.h"
#include "mml/base/Geometry/Geometry2DCore/Geometry2DLines.h"
#include "mml/base/Geometry/Geometry2DCore/Geometry2DTriangle.h"
#include "mml/base/Geometry/Geometry2DCore/Geometry2DPolygon.h"

namespace MML
{
	/// @brief Circle positioned in 2D space with center and radius
	/// @details Provides comprehensive circle operations including point containment,
	///          closest point, tangent vectors, arc/sector calculations, and
	///          intersection tests with other circles and lines.
	class Circle2D
	{
	private:
		Point2Cartesian _center;
		Real _radius;

	public:
		/// @brief Default constructor - unit circle at origin
		Circle2D() : _center(0, 0), _radius(1) {}
		/// @brief Constructs a circle from center point and radius
		/// @param center Center point
		/// @param radius Radius (must be non-negative)
		/// @throws GeometryError if radius is negative
		Circle2D(const Point2Cartesian& center, Real radius) : _center(center), _radius(radius)
		{
			if (radius < 0)
				throw GeometryError("Circle2D - radius must be non-negative");
		}
		/// @brief Constructs a circle from center coordinates and radius
		/// @param centerX Center x-coordinate
		/// @param centerY Center y-coordinate
		/// @param radius Radius (must be non-negative)
		/// @throws GeometryError if radius is negative
		Circle2D(Real centerX, Real centerY, Real radius) : _center(centerX, centerY), _radius(radius)
		{
			if (radius < 0)
				throw GeometryError("Circle2D - radius must be non-negative");
		}

		/// @brief Creates a circle from two diameter endpoints
		static Circle2D FromDiameter(const Point2Cartesian& p1, const Point2Cartesian& p2)
		{
			Point2Cartesian center((p1.X() + p2.X()) / 2, (p1.Y() + p2.Y()) / 2);
			Real radius = p1.Dist(p2) / 2;
			return Circle2D(center, radius);
		}

		/// @brief Creates the circumcircle of three points
		static Circle2D FromThreePoints(const Point2Cartesian& p1, const Point2Cartesian& p2, const Point2Cartesian& p3)
		{
			// Circumcircle of triangle
			Triangle2D tri(p1, p2, p3);
			return Circle2D(tri.Circumcenter(), tri.CircumRadius());
		}

		/// @brief Creates a unit circle centered at origin
		static Circle2D UnitCircle() { return Circle2D(Point2Cartesian(0, 0), 1); }

		/// @brief Gets the center (const version)
		Point2Cartesian  Center() const { return _center; }
		/// @brief Gets the center (mutable reference)
		Point2Cartesian& Center() { return _center; }
		/// @brief Gets the radius (const version)
		Real  Radius() const { return _radius; }
		/// @brief Gets the radius (mutable reference)
		Real& Radius() { return _radius; }
		/// @brief Gets the diameter
		Real  Diameter() const { return 2 * _radius; }

		/// @brief Calculates the area
		Real Area() const { return Constants::PI * _radius * _radius; }
		/// @brief Calculates the circumference
		Real Circumference() const { return 2 * Constants::PI * _radius; }
		/// @brief Alias for Circumference()
		Real Perimeter() const { return Circumference(); }

		/// @brief Gets the point on the circle at angle theta
		/// @param theta Angle from positive x-axis (in radians)
		Point2Cartesian PointAt(Real theta) const
		{
			return Point2Cartesian(_center.X() + _radius * std::cos(theta),
			                       _center.Y() + _radius * std::sin(theta));
		}

		/// @brief Calculates arc length for given central angle
		/// @param theta Central angle in radians
		Real ArcLength(Real theta) const { return _radius * std::abs(theta); }

		/// @brief Calculates sector area for given central angle
		/// @param theta Central angle in radians
		Real SectorArea(Real theta) const { return 0.5 * _radius * _radius * std::abs(theta); }

		/// @brief Calculates chord length for given central angle
		/// @param theta Central angle in radians
		Real ChordLength(Real theta) const { return 2 * _radius * std::sin(std::abs(theta) / 2); }

		/// @brief Gets distance from center to a point
		Real DistanceFromCenter(const Point2Cartesian& p) const { return _center.Dist(p); }

		/// @brief Gets signed distance from circle boundary to point
		/// @return Negative if inside, positive if outside
		Real SignedDistance(const Point2Cartesian& p) const { return DistanceFromCenter(p) - _radius; }
		/// @brief Gets absolute distance from circle boundary to point
		Real DistanceToPoint(const Point2Cartesian& p) const { return std::abs(SignedDistance(p)); }

		/// @brief Tests if a point is inside or on the circle
		/// @param p The point to test
		/// @param epsilon Tolerance for boundary cases
		bool Contains(const Point2Cartesian& p, Real epsilon = Constants::GEOMETRY_EPSILON) const
		{
			return DistanceFromCenter(p) <= _radius + epsilon;
		}

		/// @brief Tests if a point is on the circle boundary
		/// @param p The point to test
		/// @param epsilon Tolerance for the test
		bool IsOnBoundary(const Point2Cartesian& p, Real epsilon = Constants::GEOMETRY_EPSILON) const
		{
			return std::abs(DistanceFromCenter(p) - _radius) < epsilon;
		}

		/// @brief Finds the closest point on the circle to a given point
		/// @param p The point
		/// @return Point on circle closest to p
		Point2Cartesian ClosestPoint(const Point2Cartesian& p) const
		{
			Real dist = DistanceFromCenter(p);
			if (dist < Constants::GEOMETRY_EPSILON)
			{
				// Point at center - return arbitrary point on circle
				return PointAt(0);
			}
			// Scale vector from center to point to radius length
			Real scale = _radius / dist;
			return Point2Cartesian(_center.X() + scale * (p.X() - _center.X()),
			                       _center.Y() + scale * (p.Y() - _center.Y()));
		}

		/// @brief Gets the tangent vector at angle theta
		/// @param theta Angle in radians
		/// @return Unit tangent vector (perpendicular to radius)
		Vector2Cartesian TangentAt(Real theta) const
		{
			return Vector2Cartesian(-std::sin(theta), std::cos(theta));
		}

		/// @brief Tests for intersection with another circle
		bool Intersects(const Circle2D& other) const
		{
			Real d = _center.Dist(other._center);
			return d <= _radius + other._radius && d >= std::abs(_radius - other._radius);
		}

		/// @brief Tests for intersection with a line
		bool Intersects(const Line2D& line) const
		{
			return line.DistanceToPoint(_center) <= _radius;
		}

		/// @brief Gets the axis-aligned bounding box
		Polygon2D::BoundingBox GetBoundingBox() const
		{
			Polygon2D::BoundingBox box;
			box.minX = _center.X() - _radius;
			box.maxX = _center.X() + _radius;
			box.minY = _center.Y() - _radius;
			box.maxY = _center.Y() + _radius;
			return box;
		}

		/// @brief Creates a translated copy of the circle
		/// @param v Translation vector
		Circle2D Translated(const Vector2Cartesian& v) const
		{
			return Circle2D(_center + v, _radius);
		}

		/// @brief Creates a scaled copy of the circle (from center)
		/// @param factor Scale factor (must be positive)
		/// @throws GeometryError if factor is not positive
		Circle2D Scaled(Real factor) const
		{
			if (factor <= 0)
				throw GeometryError("Circle2D::Scaled - factor must be positive");
			return Circle2D(_center, _radius * factor);
		}
	};

	/// @brief 2D axis-aligned bounding box (AABB)
	/// @details Represents a rectangular region aligned with the coordinate axes.
	///          Provides comprehensive operations including point/box containment,
	///          distance calculations, intersection/union, and transformation methods.
	class Box2D
	{
	private:
		Point2Cartesian _min;  // Bottom-left corner
		Point2Cartesian _max;  // Top-right corner

	public:
		/// @brief Default constructor - unit box from (0,0) to (1,1)
		Box2D() : _min(0, 0), _max(1, 1) {}
		
		/// @brief Constructs a box from corner points
		/// @param min Bottom-left corner
		/// @param max Top-right corner
		Box2D(const Point2Cartesian& min, const Point2Cartesian& max) : _min(min), _max(max)
		{
			Normalize();
		}
		
		/// @brief Constructs a box from coordinates
		/// @param minX Minimum x-coordinate
		/// @param minY Minimum y-coordinate
		/// @param maxX Maximum x-coordinate
		/// @param maxY Maximum y-coordinate
		Box2D(Real minX, Real minY, Real maxX, Real maxY) : _min(minX, minY), _max(maxX, maxY)
		{
			Normalize();
		}

		/// @brief Creates a box from center and dimensions
		static Box2D FromCenterAndSize(const Point2Cartesian& center, Real width, Real height)
		{
			Real halfW = width / 2;
			Real halfH = height / 2;
			return Box2D(center.X() - halfW, center.Y() - halfH,
			            center.X() + halfW, center.Y() + halfH);
		}

		/// @brief Creates a box from center and half-extents
		static Box2D FromCenterAndHalfExtents(const Point2Cartesian& center, Real halfWidth, Real halfHeight)
		{
			return Box2D(center.X() - halfWidth, center.Y() - halfHeight,
			            center.X() + halfWidth, center.Y() + halfHeight);
		}

		/// @brief Creates a bounding box for a set of points
		/// @throws GeometryError if points vector is empty
		static Box2D FromPoints(const std::vector<Point2Cartesian>& points)
		{
			if (points.empty())
				throw GeometryError("Box2D::FromPoints - empty point set");
			
			Real minX = points[0].X(), maxX = points[0].X();
			Real minY = points[0].Y(), maxY = points[0].Y();
			
			for (size_t i = 1; i < points.size(); ++i)
			{
				minX = std::min(minX, points[i].X());
				maxX = std::max(maxX, points[i].X());
				minY = std::min(minY, points[i].Y());
				maxY = std::max(maxY, points[i].Y());
			}
			
			return Box2D(minX, minY, maxX, maxY);
		}

		/// @brief Creates a unit box from (0,0) to (1,1)
		static Box2D UnitBox() { return Box2D(0, 0, 1, 1); }

		/// @brief Ensures min <= max for all coordinates
		void Normalize()
		{
			if (_min.X() > _max.X()) std::swap(_min.X(), _max.X());
			if (_min.Y() > _max.Y()) std::swap(_min.Y(), _max.Y());
		}

		/// @brief Gets the minimum corner (const version)
		Point2Cartesian  Min() const { return _min; }
		/// @brief Gets the minimum corner (mutable reference)
		Point2Cartesian& Min() { return _min; }
		/// @brief Gets the maximum corner (const version)
		Point2Cartesian  Max() const { return _max; }
		/// @brief Gets the maximum corner (mutable reference)
		Point2Cartesian& Max() { return _max; }

		/// @brief Gets the minimum x-coordinate
		Real MinX() const { return _min.X(); }
		/// @brief Gets the minimum y-coordinate
		Real MinY() const { return _min.Y(); }
		/// @brief Gets the maximum x-coordinate
		Real MaxX() const { return _max.X(); }
		/// @brief Gets the maximum y-coordinate
		Real MaxY() const { return _max.Y(); }

		/// @brief Gets the width
		Real Width() const { return _max.X() - _min.X(); }
		/// @brief Gets the height
		Real Height() const { return _max.Y() - _min.Y(); }
		/// @brief Gets the size as a vector
		Vector2Cartesian Size() const { return Vector2Cartesian(Width(), Height()); }
		/// @brief Gets the half-extents
		Vector2Cartesian HalfExtents() const { return Vector2Cartesian(Width() / 2, Height() / 2); }
		/// @brief Gets the diagonal length
		Real Diagonal() const { return std::sqrt(Width() * Width() + Height() * Height()); }
		/// @brief Gets the aspect ratio (width/height)
		Real AspectRatio() const { return Width() / Height(); }

		/// @brief Calculates the area
		Real Area() const { return Width() * Height(); }
		/// @brief Calculates the perimeter
		Real Perimeter() const { return 2 * (Width() + Height()); }

		/// @brief Gets the center point
		Point2Cartesian Center() const
		{
			return Point2Cartesian((_min.X() + _max.X()) / 2, (_min.Y() + _max.Y()) / 2);
		}

		/// @brief Gets the bottom-left corner
		Point2Cartesian BottomLeft() const { return _min; }
		/// @brief Gets the bottom-right corner
		Point2Cartesian BottomRight() const { return Point2Cartesian(_max.X(), _min.Y()); }
		/// @brief Gets the top-left corner
		Point2Cartesian TopLeft() const { return Point2Cartesian(_min.X(), _max.Y()); }
		/// @brief Gets the top-right corner
		Point2Cartesian TopRight() const { return _max; }

		/// @brief Gets all four corners as a vector
		std::vector<Point2Cartesian> Corners() const
		{
			return { BottomLeft(), BottomRight(), TopRight(), TopLeft() };
		}

		/// @brief Converts to a polygon
		Polygon2D ToPolygon() const
		{
			return Polygon2D({ BottomLeft(), BottomRight(), TopRight(), TopLeft() });
		}

		/// @brief Tests if a point is inside the box
		/// @param p The point to test
		/// @param epsilon Tolerance for boundary cases
		bool Contains(const Point2Cartesian& p, Real epsilon = Constants::GEOMETRY_EPSILON) const
		{
			return p.X() >= _min.X() - epsilon && p.X() <= _max.X() + epsilon &&
			       p.Y() >= _min.Y() - epsilon && p.Y() <= _max.Y() + epsilon;
		}

		/// @brief Tests if another box is fully contained within this box
		/// @param other The box to test
		/// @param epsilon Tolerance for boundary cases
		bool Contains(const Box2D& other, Real epsilon = Constants::GEOMETRY_EPSILON) const
		{
			return other._min.X() >= _min.X() - epsilon && other._max.X() <= _max.X() + epsilon &&
			       other._min.Y() >= _min.Y() - epsilon && other._max.Y() <= _max.Y() + epsilon;
		}

		/// @brief Calculates distance from point to box
		/// @param p The point
		/// @return Distance (0 if inside)
		Real DistanceToPoint(const Point2Cartesian& p) const
		{
			Real dx = std::max({ _min.X() - p.X(), REAL(0.0), p.X() - _max.X() });
			Real dy = std::max({ _min.Y() - p.Y(), REAL(0.0), p.Y() - _max.Y() });
			return std::sqrt(dx * dx + dy * dy);
		}

		/// @brief Finds the closest point on the box boundary to a given point
		/// @param p The point
		/// @return Point on box closest to p
		Point2Cartesian ClosestPoint(const Point2Cartesian& p) const
		{
			return Point2Cartesian(std::clamp(p.X(), _min.X(), _max.X()),
			                       std::clamp(p.Y(), _min.Y(), _max.Y()));
		}

		/// @brief Tests for intersection with another box
		bool Intersects(const Box2D& other) const
		{
			return _min.X() <= other._max.X() && _max.X() >= other._min.X() &&
			       _min.Y() <= other._max.Y() && _max.Y() >= other._min.Y();
		}

		/// @brief Computes the intersection of two boxes
		/// @return Intersection box, or empty box (0,0,0,0) if no intersection
		Box2D Intersection(const Box2D& other) const
		{
			Real minX = std::max(_min.X(), other._min.X());
			Real minY = std::max(_min.Y(), other._min.Y());
			Real maxX = std::min(_max.X(), other._max.X());
			Real maxY = std::min(_max.Y(), other._max.Y());
			
			if (minX > maxX || minY > maxY)
				return Box2D(0, 0, 0, 0);  // Empty box
			
			return Box2D(minX, minY, maxX, maxY);
		}

		/// @brief Computes the union of two boxes
		/// @return Smallest box containing both boxes
		Box2D Union(const Box2D& other) const
		{
			return Box2D(std::min(_min.X(), other._min.X()),
			            std::min(_min.Y(), other._min.Y()),
			            std::max(_max.X(), other._max.X()),
			            std::max(_max.Y(), other._max.Y()));
		}

		/// @brief Expands the box to include a point
		/// @param p The point to include
		void ExpandToInclude(const Point2Cartesian& p)
		{
			_min.X() = std::min(_min.X(), p.X());
			_min.Y() = std::min(_min.Y(), p.Y());
			_max.X() = std::max(_max.X(), p.X());
			_max.Y() = std::max(_max.Y(), p.Y());
		}

		/// @brief Creates a box expanded by a uniform margin
		/// @param margin Expansion amount in all directions
		Box2D Expanded(Real margin) const
		{
			return Box2D(_min.X() - margin, _min.Y() - margin,
			            _max.X() + margin, _max.Y() + margin);
		}

		/// @brief Tests if the box is degenerate (zero area)
		/// @param epsilon Tolerance for the test
		bool IsDegenerate(Real epsilon = Constants::GEOMETRY_EPSILON) const
		{
			return Width() < epsilon || Height() < epsilon;
		}

		/// @brief Creates a translated copy of the box
		/// @param v Translation vector
		Box2D Translated(const Vector2Cartesian& v) const
		{
			return Box2D(_min + v, _max + v);
		}

		/// @brief Creates a scaled copy of the box (from center)
		/// @param factor Scale factor
		Box2D Scaled(Real factor) const
		{
			Point2Cartesian c = Center();
			Vector2Cartesian half = HalfExtents();
			return Box2D::FromCenterAndHalfExtents(c, half.X() * factor, half.Y() * factor);
		}

		/// @brief Tests for intersection with a circle
		bool Intersects(const Circle2D& circle) const
		{
			Point2Cartesian closest = ClosestPoint(circle.Center());
			return circle.Center().Dist(closest) <= circle.Radius();
		}
	};

} // namespace MML

#endif // MML_GEOMETRY_2D_CIRCLE_BOX_H

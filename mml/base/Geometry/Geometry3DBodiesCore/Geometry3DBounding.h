///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Geometry3DBounding.h                                                ///
///  Description: 3D Bounding volumes (BoundingSphere3D, Box3D AABB)                  ///
///               Fast spatial queries and collision detection                        ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                        ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

/// @file Geometry3DBounding.h
/// @brief 3D bounding volume classes for spatial queries.
/// - BoundingSphere3D: Rotation-invariant bounding sphere
/// - Box3D: Axis-Aligned Bounding Box (AABB)
/// @see Geometry3DBodies.h for the aggregate header

#if !defined MML_GEOMETRY_3D_BOUNDING_H
#define MML_GEOMETRY_3D_BOUNDING_H

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "mml/MMLBase.h"
#include "mml/MMLExceptions.h"
#include "mml/base/Vector/VectorTypes.h"
#include "mml/base/Geometry/Geometry3D.h"

namespace MML {
	// Forward declaration
	class Box3D;

	/// @brief Bounding sphere for fast spatial queries and collision detection.
	/// A bounding sphere is defined by a center point and radius. While less tight-fitting
	/// than an AABB, sphere-sphere intersection tests are extremely fast (single distance
	/// comparison), making them ideal for broad-phase collision detection.
	/// @note Bounding spheres are rotation-invariant, unlike AABBs.

	class BoundingSphere3D {
	private:
		Pnt3Cart _center; ///< Center point of the sphere
		Real _radius;	  ///< Radius (must be non-negative)

	public:
		/// @name Constructors
		/// @{

		/// /** @brief Default constructor: unit sphere at origin. */

		BoundingSphere3D()
			: _center(0, 0, 0)
			, _radius(0) {}

		/// @brief Construct from center and radius.
		/// @param center Center point
		/// @param radius Sphere radius (must be ≥ 0)
		/// @throws ArgumentError if radius is negative

		BoundingSphere3D(const Pnt3Cart& center, Real radius)
			: _center(center)
			, _radius(radius) {
			if (radius < 0.0)
				throw ArgumentError("BoundingSphere3D: radius must be non-negative");
		}
		/// @}

		/// @name Accessors
		/// @{
		Pnt3Cart Center() const { return _center; }
		Real Radius() const { return _radius; }
		/// @}

		/// /** @brief Calculate sphere volume: V = 4/3 π r³ */

		Real Volume() const { return (4.0 / 3.0) * Constants::PI * _radius * _radius * _radius; }

		/// @brief Test if a point lies inside or on the sphere.
		/// @param pnt Point to test
		/// @return true if distance from center ≤ radius

		bool Contains(const Pnt3Cart& pnt) const { return _center.Dist(pnt) <= _radius; }

		/// @brief Test intersection with another bounding sphere.
		/// @param other The other sphere
		/// @return true if spheres overlap or touch

		bool Intersects(const BoundingSphere3D& other) const {
			Real centerDist = _center.Dist(other._center);
			return centerDist <= (_radius + other._radius);
		}

		/// @brief Expand sphere to include a point.
		/// @param pnt Point to include
		/// @note Only increases radius; center remains fixed.

		void ExpandToInclude(const Pnt3Cart& pnt) {
			Real dist = _center.Dist(pnt);
			if (dist > _radius)
				_radius = dist;
		}

		std::string ToString() const {
			std::ostringstream oss;
			oss << "BoundingSphere3D[Center=(" << _center.X() << "," << _center.Y() << "," << _center.Z() << ")"
				<< ", Radius=" << _radius << "]";
			return oss.str();
		}
	};

	/// @brief Axis-Aligned Bounding Box (AABB) for 3D spatial queries.
	/// An AABB is defined by its minimum and maximum corner points, with edges
	/// parallel to the coordinate axes. AABBs provide:
	/// - Fast intersection tests (separating axis theorem)
	/// - Tight bounds for axis-aligned geometry
	/// - Efficient spatial partitioning (octrees, BVH)
	/// ## Factory Methods
	/// - FromCenterAndSize() - From center point and dimensions
	/// - FromCenterAndHalfExtents() - From center and half-sizes
	/// - FromPoints() - Minimal box containing point set
	/// - UnitBox(), CenteredUnitCube() - Standard unit boxes
	/// @note Coordinates are automatically normalized (min ≤ max) on construction.

	class Box3D {
	private:
		Pnt3Cart _min; ///< Minimum corner (smallest x, y, z)
		Pnt3Cart _max; ///< Maximum corner (largest x, y, z)

	public:
		/// @name Constructors
		/// @{

		/// /** @brief Default constructor: unit cube [0,1]³. */

		Box3D()
			: _min(0, 0, 0)
			, _max(1, 1, 1) {}

		/// @brief Construct from min/max corner points.
		/// @param min Minimum corner point
		/// @param max Maximum corner point
		/// @note Automatically normalizes if min > max in any dimension.

		Box3D(const Pnt3Cart& min, const Pnt3Cart& max)
			: _min(min)
			, _max(max) {
			Normalize();
		}

		/// @brief Construct from explicit coordinates.
		/// @param minX,minY,minZ Minimum corner coordinates
		/// @param maxX,maxY,maxZ Maximum corner coordinates

		Box3D(Real minX, Real minY, Real minZ, Real maxX, Real maxY, Real maxZ)
			: _min(minX, minY, minZ)
			, _max(maxX, maxY, maxZ) {
			Normalize();
		}
		/// @}

		/// @name Factory Methods
		/// @{

		/// @brief Create box from center point and full dimensions.
		/// @param center Center point of the box
		/// @param width Size in X direction
		/// @param height Size in Y direction
		/// @param depth Size in Z direction

		static Box3D FromCenterAndSize(const Pnt3Cart& center, Real width, Real height, Real depth) {
			Real halfW = width / 2;
			Real halfH = height / 2;
			Real halfD = depth / 2;
			return Box3D(center.X() - halfW, center.Y() - halfH, center.Z() - halfD, center.X() + halfW, center.Y() + halfH,
						 center.Z() + halfD);
		}

		/// /** @brief Create box from center and half-extents (scalars). */

		static Box3D FromCenterAndHalfExtents(const Pnt3Cart& center, Real halfWidth, Real halfHeight, Real halfDepth) {
			return Box3D(center.X() - halfWidth, center.Y() - halfHeight, center.Z() - halfDepth, center.X() + halfWidth,
						 center.Y() + halfHeight, center.Z() + halfDepth);
		}

		/// /** @brief Create box from center and half-extents vector. */

		static Box3D FromCenterAndHalfExtents(const Pnt3Cart& center, const Vec3Cart& halfExtents) {
			return Box3D(center.X() - halfExtents.X(), center.Y() - halfExtents.Y(), center.Z() - halfExtents.Z(),
						 center.X() + halfExtents.X(), center.Y() + halfExtents.Y(), center.Z() + halfExtents.Z());
		}

		/// @brief Create minimal AABB containing a set of points.
		/// @param points Vector of points to bound
		/// @throws GeometryError if points is empty

		static Box3D FromPoints(const std::vector<Pnt3Cart>& points) {
			if (points.empty())
				throw GeometryError("Box3D::FromPoints - empty point set");

			Real minX = points[0].X(), maxX = points[0].X();
			Real minY = points[0].Y(), maxY = points[0].Y();
			Real minZ = points[0].Z(), maxZ = points[0].Z();

			for (size_t i = 1; i < points.size(); ++i) {
				minX = std::min(minX, points[i].X());
				maxX = std::max(maxX, points[i].X());
				minY = std::min(minY, points[i].Y());
				maxY = std::max(maxY, points[i].Y());
				minZ = std::min(minZ, points[i].Z());
				maxZ = std::max(maxZ, points[i].Z());
			}

			return Box3D(minX, minY, minZ, maxX, maxY, maxZ);
		}

		/// /** @brief Unit box [0,1]³. */

		static Box3D UnitBox() { return Box3D(0, 0, 0, 1, 1, 1); }

		/// /** @brief Unit cube [0,1]³ (alias for UnitBox). */

		static Box3D UnitCube() { return Box3D(0, 0, 0, 1, 1, 1); }

		/// /** @brief Centered unit cube [-0.5, 0.5]³. */

		static Box3D CenteredUnitCube() { return Box3D(-0.5, -0.5, -0.5, 0.5, 0.5, 0.5); }
		/// @}

		/// /** @brief Ensure min ≤ max in all dimensions (swaps if needed). */

		void Normalize() {
			if (_min.X() > _max.X())
				std::swap(_min.X(), _max.X());
			if (_min.Y() > _max.Y())
				std::swap(_min.Y(), _max.Y());
			if (_min.Z() > _max.Z())
				std::swap(_min.Z(), _max.Z());
		}

		/// @name Accessors
		/// @{
		Pnt3Cart Min() const { return _min; }
		Pnt3Cart& Min() { return _min; }
		Pnt3Cart Max() const { return _max; }
		Pnt3Cart& Max() { return _max; }

		Real MinX() const { return _min.X(); }
		Real MinY() const { return _min.Y(); }
		Real MinZ() const { return _min.Z(); }
		Real MaxX() const { return _max.X(); }
		Real MaxY() const { return _max.Y(); }
		Real MaxZ() const { return _max.Z(); }
		/// @}

		/// @name Dimensions
		/// @{
		Real Width() const { return _max.X() - _min.X(); }	///< Size in X direction
		Real Height() const { return _max.Y() - _min.Y(); } ///< Size in Y direction
		Real Depth() const { return _max.Z() - _min.Z(); }	///< Size in Z direction
		Vec3Cart Size() const { return Vec3Cart(Width(), Height(), Depth()); }
		Vec3Cart HalfExtents() const { return Vec3Cart(Width() / 2, Height() / 2, Depth() / 2); }
		Real Diagonal() const { return std::sqrt(Width() * Width() + Height() * Height() + Depth() * Depth()); }
		/// @}

		/// @name Geometric Properties
		/// @{

		/// /** @brief Calculate box volume: V = W × H × D */

		Real Volume() const { return Width() * Height() * Depth(); }
		/// /** @brief Calculate surface area: A = 2(WH + HD + DW) */

		Real SurfaceArea() const { return 2 * (Width() * Height() + Height() * Depth() + Depth() * Width()); }

		/// /** @brief Get the center point of the box. */

		Pnt3Cart Center() const { return Pnt3Cart((_min.X() + _max.X()) / 2, (_min.Y() + _max.Y()) / 2, (_min.Z() + _max.Z()) / 2); }
		/// @}

		/// @name Corner Access
		/// @{
		/// /** @brief Get corner at (minX, minY, minZ). */

		Pnt3Cart Corner000() const { return _min; }
		Pnt3Cart Corner100() const { return Pnt3Cart(_max.X(), _min.Y(), _min.Z()); }
		Pnt3Cart Corner010() const { return Pnt3Cart(_min.X(), _max.Y(), _min.Z()); }
		Pnt3Cart Corner110() const { return Pnt3Cart(_max.X(), _max.Y(), _min.Z()); }
		Pnt3Cart Corner001() const { return Pnt3Cart(_min.X(), _min.Y(), _max.Z()); }
		Pnt3Cart Corner101() const { return Pnt3Cart(_max.X(), _min.Y(), _max.Z()); }
		Pnt3Cart Corner011() const { return Pnt3Cart(_min.X(), _max.Y(), _max.Z()); }
		/// /** @brief Get corner at (maxX, maxY, maxZ). */

		Pnt3Cart Corner111() const { return _max; }

		/// /** @brief Get all 8 corners as a vector. */

		std::vector<Pnt3Cart> Corners() const {
			return {Corner000(), Corner100(), Corner010(), Corner110(), Corner001(), Corner101(), Corner011(), Corner111()};
		}
		/// @}

		/// @name Spatial Queries
		/// @{

		/// @brief Test if a point is inside or on the box boundary.
		/// @param p Point to test
		/// @param epsilon Tolerance for boundary tests

		bool Contains(const Pnt3Cart& p, Real epsilon = Constants::GEOMETRY_EPSILON) const {
			return p.X() >= _min.X() - epsilon && p.X() <= _max.X() + epsilon && p.Y() >= _min.Y() - epsilon &&
				   p.Y() <= _max.Y() + epsilon && p.Z() >= _min.Z() - epsilon && p.Z() <= _max.Z() + epsilon;
		}

		/// @brief Test if another box is fully contained within this box.
		/// @param other Box to test
		/// @param epsilon Tolerance for boundary tests

		bool Contains(const Box3D& other, Real epsilon = Constants::GEOMETRY_EPSILON) const {
			return other._min.X() >= _min.X() - epsilon && other._max.X() <= _max.X() + epsilon && other._min.Y() >= _min.Y() - epsilon &&
				   other._max.Y() <= _max.Y() + epsilon && other._min.Z() >= _min.Z() - epsilon && other._max.Z() <= _max.Z() + epsilon;
		}

		/// @brief Calculate distance from a point to the box surface.
		/// @param p Point to measure from
		/// @return 0 if inside, positive distance otherwise

		Real DistanceToPoint(const Pnt3Cart& p) const {
			Real dx = std::max({_min.X() - p.X(), REAL(0.0), p.X() - _max.X()});
			Real dy = std::max({_min.Y() - p.Y(), REAL(0.0), p.Y() - _max.Y()});
			Real dz = std::max({_min.Z() - p.Z(), REAL(0.0), p.Z() - _max.Z()});
			return std::sqrt(dx * dx + dy * dy + dz * dz);
		}

		/// @brief Find the closest point on the box to a given point.
		/// @param p External point
		/// @return Closest point on or inside the box

		Pnt3Cart ClosestPoint(const Pnt3Cart& p) const {
			return Pnt3Cart(std::clamp(p.X(), _min.X(), _max.X()), std::clamp(p.Y(), _min.Y(), _max.Y()),
							std::clamp(p.Z(), _min.Z(), _max.Z()));
		}

		/// @brief Test if two boxes overlap.
		/// @param other Box to test against
		/// @return true if boxes intersect (including touching)

		bool Intersects(const Box3D& other) const {
			return _min.X() <= other._max.X() && _max.X() >= other._min.X() && _min.Y() <= other._max.Y() && _max.Y() >= other._min.Y() &&
				   _min.Z() <= other._max.Z() && _max.Z() >= other._min.Z();
		}
		/// @}

		/// @name Set Operations
		/// @{

		/// @brief Compute intersection of two boxes.
		/// @param other Box to intersect with
		/// @return Intersection box, or empty box (0,0,0,0,0,0) if disjoint

		Box3D Intersection(const Box3D& other) const {
			Real minX = std::max(_min.X(), other._min.X());
			Real minY = std::max(_min.Y(), other._min.Y());
			Real minZ = std::max(_min.Z(), other._min.Z());
			Real maxX = std::min(_max.X(), other._max.X());
			Real maxY = std::min(_max.Y(), other._max.Y());
			Real maxZ = std::min(_max.Z(), other._max.Z());

			if (minX > maxX || minY > maxY || minZ > maxZ)
				return Box3D(0, 0, 0, 0, 0, 0); // Empty box

			return Box3D(minX, minY, minZ, maxX, maxY, maxZ);
		}

		/// @brief Compute union (smallest enclosing box) of two boxes.
		/// @param other Box to merge with
		/// @return Minimal AABB containing both boxes

		Box3D Union(const Box3D& other) const {
			return Box3D(std::min(_min.X(), other._min.X()), std::min(_min.Y(), other._min.Y()), std::min(_min.Z(), other._min.Z()),
						 std::max(_max.X(), other._max.X()), std::max(_max.Y(), other._max.Y()), std::max(_max.Z(), other._max.Z()));
		}

		/// @brief Expand box in-place to include a point.
		/// @param p Point to include

		void ExpandToInclude(const Pnt3Cart& p) {
			_min.X() = std::min(_min.X(), p.X());
			_min.Y() = std::min(_min.Y(), p.Y());
			_min.Z() = std::min(_min.Z(), p.Z());
			_max.X() = std::max(_max.X(), p.X());
			_max.Y() = std::max(_max.Y(), p.Y());
			_max.Z() = std::max(_max.Z(), p.Z());
		}

		/// @brief Create expanded copy with uniform margin.
		/// @param margin Distance to expand in all directions

		Box3D Expanded(Real margin) const {
			return Box3D(_min.X() - margin, _min.Y() - margin, _min.Z() - margin, _max.X() + margin, _max.Y() + margin, _max.Z() + margin);
		}
		/// @}

		/// @name Transformations & Utilities
		/// @{

		/// @brief Check if box is degenerate (zero volume in any dimension).
		/// @param epsilon Tolerance for near-zero dimension

		bool IsDegenerate(Real epsilon = Constants::GEOMETRY_EPSILON) const {
			return Width() < epsilon || Height() < epsilon || Depth() < epsilon;
		}

		/// @brief Create translated copy.
		/// @param v Translation vector

		Box3D Translated(const Vec3Cart& v) const { return Box3D(_min + v, _max + v); }

		/// @brief Create scaled copy (about center).
		/// @param factor Scale factor

		Box3D Scaled(Real factor) const {
			Pnt3Cart c = Center();
			Vec3Cart half = HalfExtents();
			return Box3D::FromCenterAndHalfExtents(c, half.X() * factor, half.Y() * factor, half.Z() * factor);
		}
		/// @}

		std::string ToString() const {
			std::ostringstream oss;
			oss << "Box3D[Min=(" << _min.X() << "," << _min.Y() << "," << _min.Z() << ")"
				<< ", Max=(" << _max.X() << "," << _max.Y() << "," << _max.Z() << ")"
				<< ", Size=(" << Width() << "×" << Height() << "×" << Depth() << ")]";
			return oss.str();
		}
	};

} // namespace MML

#endif

///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Geometry3DPlane.h                                                   ///
///  Description: Plane3D class for 3D plane representation and operations            ///
///               Point/line/plane tests, distances, and intersections               ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_GEOMETRY_3D_PLANE_H
#define MML_GEOMETRY_3D_PLANE_H

#include <mml/MMLBase.h>
#include <mml/base/Vector/VectorN.h>
#include <mml/base/Vector/VectorTypes3D.h>
#include <mml/base/Geometry/Geometry.h>
#include <mml/base/BaseUtils.h>
#include <mml/base/Geometry/Geometry3D/Geometry3DLines.h>

namespace MML
{

	/// @brief Plane in 3D space
	/// @details Represented in general form Ax + By + Cz + D = 0 where (A,B,C) is the unit normal.
	///          Provides multiple construction methods, point/line/plane operations,
	///          distance calculations, and intersection computations.
	class Plane3D
	{
	private:
		Real _coefA, _coefB, _coefC, _coefD;

	public:
		/// @brief Constructs a plane from a point and normal vector
		/// @param a A point on the plane
		/// @param normal The plane normal (will be normalized)
		/// @throws GeometryError if normal is null vector
		Plane3D(const Pnt3Cart& a, const Vec3Cart& normal)
		{
			// same effective threshold GetAsUnitVector uses, so degenerate normals fail here with GeometryError
			if (normal.IsNearZero(Precision::RoundoffZeroThreshold))
				throw GeometryError("Plane3D ctor - normal is null vector");

			Vec3Cart unitNormal = normal.Normalized();

			_coefA = unitNormal.X();
			_coefB = unitNormal.Y();
			_coefC = unitNormal.Z();
			_coefD = -(a.X() * unitNormal.X() + a.Y() * unitNormal.Y() + a.Z() * unitNormal.Z());
		}
		
		/// @brief Constructs a plane from three points
		/// @param a First point
		/// @param b Second point
		/// @param c Third point
		Plane3D(const Pnt3Cart& a, const Pnt3Cart& b, const Pnt3Cart& c)
			: Plane3D(a, VectorProduct(Vec3Cart(a, b), Vec3Cart(a, c)))
		{
		}
		
		/// @brief Constructs a plane in Hesse normal form
		/// @param alpha Direction cosine angle for x
		/// @param beta Direction cosine angle for y
		/// @param gamma Direction cosine angle for z
		/// @param d Distance from origin
		Plane3D(Real alpha, Real beta, Real gamma, Real d)
		{
			_coefA = cos(alpha);
			_coefB = cos(beta);
			_coefC = cos(gamma);
			_coefD = -d;
		}
		
		/// @brief Constructs a plane from intercepts on coordinate axes
		/// @param seg_x X-axis intercept
		/// @param seg_y Y-axis intercept
		/// @param seg_z Z-axis intercept
		/// @throws GeometryError if any intercept is zero
		Plane3D(Real seg_x, Real seg_y, Real seg_z)
		{
			if (seg_x == 0 || seg_y == 0 || seg_z == 0)
				throw GeometryError("Plane3D ctor - zero segment");

			_coefA = 1 / seg_x;
			_coefB = 1 / seg_y;
			_coefC = 1 / seg_z;
			_coefD = -1;
		}

		/// @brief Gets the XY plane (z = 0)
		static Plane3D XYPlane() { return Plane3D(Pnt3Cart(0, 0, 0), Vec3Cart(0, 0, 1)); }
		/// @deprecated Use XYPlane() (CONTRACT_POLICY §11: factories are PascalCase nouns).
		[[deprecated("use XYPlane()")]]
		static Plane3D GetXYPlane()  { return XYPlane(); }
		/// @brief Gets the XZ plane (y = 0)
		static Plane3D XZPlane() { return Plane3D(Pnt3Cart(0, 0, 0), Vec3Cart(0, 1, 0)); }
		/// @deprecated Use XZPlane() (CONTRACT_POLICY §11: factories are PascalCase nouns).
		[[deprecated("use XZPlane()")]]
		static Plane3D GetXZPlane()  { return XZPlane(); }
		/// @brief Gets the YZ plane (x = 0)
		static Plane3D YZPlane() { return Plane3D(Pnt3Cart(0, 0, 0), Vec3Cart(1, 0, 0)); }
		/// @deprecated Use YZPlane() (CONTRACT_POLICY §11: factories are PascalCase nouns).
		[[deprecated("use YZPlane()")]]
		static Plane3D GetYZPlane()  { return YZPlane(); }

		/// @brief Gets coefficient A (const)
		Real  A() const { return _coefA; }
		/// @brief Gets coefficient A (mutable)
		Real& A() { return _coefA; }
		/// @brief Gets coefficient B (const)
		Real  B() const { return _coefB; }
		/// @brief Gets coefficient B (mutable)
		Real& B() { return _coefB; }
		/// @brief Gets coefficient C (const)
		Real  C() const { return _coefC; }
		/// @brief Gets coefficient C (mutable)
		Real& C() { return _coefC; }
		/// @brief Gets coefficient D (const)
		Real  D() const { return _coefD; }
		/// @brief Gets coefficient D (mutable)
		Real& D() { return _coefD; }

		/// @brief Gets the unit normal vector
		Vec3Cart	Normal() const { return Vec3Cart(_coefA, _coefB, _coefC); }
		
		/// @brief Gets an arbitrary point on the plane
		Pnt3Cart	getPointOnPlane() const {
			if (_coefA != 0.0)
				return Pnt3Cart(-_coefD / _coefA, 0, 0);
			else  if (_coefB != 0.0)
				return Pnt3Cart(0, -_coefD / _coefB, 0);
			else
				return Pnt3Cart(0, 0, -_coefD / _coefC);
		}
		/// @deprecated Use getPointOnPlane() (CONTRACT_POLICY §11: accessors are getX()).
		[[deprecated("use getPointOnPlane()")]]
		Pnt3Cart	GetPointOnPlane() const { return getPointOnPlane(); }

		/// @brief Gets the intercepts on coordinate axes
		/// @param outseg_x Output: X-axis intercept (PosInf if parallel)
		/// @param outseg_y Output: Y-axis intercept (PosInf if parallel)
		/// @param outseg_z Output: Z-axis intercept (PosInf if parallel)
		void getCoordAxisSegments(Real& outseg_x, Real& outseg_y, Real& outseg_z) const
		{
			if (_coefA != 0.0)
				outseg_x = -_coefD / _coefA;
			else
				outseg_x = Constants::PosInf;

			if (_coefB != 0.0)
				outseg_y = -_coefD / _coefB;
			else
				outseg_y = Constants::PosInf;

			if (_coefC != 0.0)
				outseg_z = -_coefD / _coefC;
			else
				outseg_z = Constants::PosInf;
		}
		/// @deprecated Use getCoordAxisSegments().
		[[deprecated("use getCoordAxisSegments()")]]
		void GetCoordAxisSegments(Real& outseg_x, Real& outseg_y, Real& outseg_z) const { getCoordAxisSegments(outseg_x, outseg_y, outseg_z); }

		/// @brief Tests if a point lies on this plane
		/// @param pnt The point to test
		/// @param defEps Tolerance for the test
		bool IsPointOnPlane(const Pnt3Cart& pnt, Real defEps = Precision::Plane3DIsPointOnPlaneTolerance) const
		{
			return std::abs(_coefA * pnt.X() + _coefB * pnt.Y() + _coefC * pnt.Z() + _coefD) < defEps;
		}
		
		/// @brief Calculates the distance from a point to this plane
		Real DistToPoint(const Pnt3Cart& pnt) const
		{
			Real a = _coefA * pnt.X() + _coefB * pnt.Y() + _coefC * pnt.Z() + _coefD;
			Real b = sqrt(_coefA * _coefA + _coefB * _coefB + _coefC * _coefC);

			return std::abs(a / b);
		}

		/// @brief Projects a point onto this plane
		/// @param pnt The point to project
		/// @return The projected point on the plane
		/// @throws GeometryError if plane normal is zero
		Pnt3Cart ProjectionToPlane(const Pnt3Cart& pnt) const
		{
			// Plane: n·x + D = 0
			// Projection: p_proj = pnt - n * (n·pnt + D) / |n|^2
			Vec3Cart n = Normal();
			
			Real n_norm2 = n.NormL2() * n.NormL2();
			if (n_norm2 == 0.0)
				throw GeometryError("ProjectionToPlane: plane normal is zero vector");

			Real dist = (n * pnt + D()) / n_norm2;
			
			return pnt - n * dist;
		}

		/// @brief Tests if a line lies entirely on this plane
		bool IsLineOnPlane(const Line3D& line) const
		{
			// get two points
			const Pnt3Cart pnt1 = line.StartPoint();
			const Pnt3Cart pnt2 = line(1.0);

			if (IsPointOnPlane(pnt1) && IsPointOnPlane(pnt2))
				return true;
			else
				return false;
		}
		
		/// @brief Calculates the angle between a line and this plane
		/// @return Angle in radians
		Real AngleToLine(const Line3D& line) const
		{
			// angle between line and normal to plane
			return Constants::PI / 2.0 - line.Direction().AngleToVector(this->Normal());
		}
		
		/// @brief Finds the intersection point of a line with this plane
		/// @param line The line
		/// @param out_inter_pnt Output: intersection point if exists
		/// @return false if line is parallel to plane
		bool IntersectionWithLine(const Line3D& line, Pnt3Cart& out_inter_pnt, Real eps = Precision::Line3DIsParallelTolerance) const
		{
			Vec3Cart line_dir = line.Direction();
			Vec3Cart plane_normal = Normal();

			// Check if the line is parallel to the plane
			Real dot = ScalarProduct(line_dir, plane_normal);
			if (std::abs(dot) < eps) {
				return false;
			}

			// Calculate the distance between the point on the line and the plane
			double t = -(plane_normal * line.StartPoint() + D()) / dot;
			Pnt3Cart inter_pnt = line.StartPoint() + line_dir * t;
			out_inter_pnt = inter_pnt;
			
			return true;
		}

		/// @brief Tests if this plane is parallel to another plane
		bool IsParallelToPlane(const Plane3D& plane) const
		{
			Vec3Cart norm1(_coefA, _coefB, _coefC);
			Vec3Cart norm2(plane._coefA, plane._coefB, plane._coefC);

			return norm1.IsParallelTo(norm2);
		}
		
		/// @brief Tests if this plane is perpendicular to another plane
		bool IsPerpendicularToPlane(const Plane3D& plane) const
		{
			Vec3Cart norm1(_coefA, _coefB, _coefC);
			Vec3Cart norm2(plane._coefA, plane._coefB, plane._coefC);

			return norm1.IsPerpendicularTo(norm2);
		}
		
		/// @brief Calculates the angle between two planes
		/// @return Angle in radians
		Real AngleToPlane(const Plane3D& plane) const
		{
			// angle between normals of two planes
			return this->Normal().AngleToVector(plane.Normal());
		}
		
		/// @brief Calculates the distance between two planes
		/// @return Distance (0 if planes intersect)
		Real DistToPlane(const Plane3D& plane) const
		{
			if (IsParallelToPlane(plane))
			{
				// distance between two parallel planes is the distance between any point on one plane and the other plane
				return DistToPoint(plane.getPointOnPlane());
			}
			else
				return 0.0;
		}

		/// @brief Finds the intersection line of two planes
		/// @param plane The other plane
		/// @param out_inter_line Output: intersection line if exists
		/// @return false if planes are parallel
		bool IntersectionWithPlane(const Plane3D& plane, Line3D& out_inter_line) const
		{
			Vec3Cart n1 = Normal();
			Vec3Cart n2 = plane.Normal();
			Vec3Cart dir = VectorProduct(n1, n2);

			if (dir.NormL2() < Precision::Vec3CartIsParallelTolerance) {
				// Planes are parallel
				return false;
			}

			// Solve for a point on the intersection line:
			// x = ( (d2 n1 - d1 n2) x (n1 x n2) ) / |n1 x n2|^2
			double d1 = -D();
			double d2 = -plane.D();

			Vec3Cart n1xn2 = dir;
			Vec3Cart p = VectorProduct((d2 * n1 - d1 * n2), n1xn2) / POW2(n1xn2.NormL2());

			out_inter_line = Line3D(Pnt3Cart(p.X(), p.Y(), p.Z()), dir.Normalized());
			return true;
		}
	};

} // namespace MML

#endif // MML_GEOMETRY_3D_PLANE_H

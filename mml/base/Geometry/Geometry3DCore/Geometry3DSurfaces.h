///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Geometry3DSurfaces.h                                                ///
///  Description: 3D surface primitives: Triangle3D, TriangleSurface3D, RectSurface3D ///
///               Parametric surface implementations for integration and analysis     ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_GEOMETRY_3D_SURFACES_H
#define MML_GEOMETRY_3D_SURFACES_H

#include "mml/MMLBase.h"
#include "mml/interfaces/IFunction.h"
#include "mml/base/Vector/VectorN.h"
#include "mml/base/Vector/VectorTypes.h"
#include "mml/base/Geometry/Geometry.h"
#include "mml/base/BaseUtils.h"
#include "mml/base/Geometry/Geometry3DCore/Geometry3DLines.h"
#include "mml/base/Geometry/Geometry3DCore/Geometry3DPlane.h"

namespace MML
{

	/// @brief Triangle in 3D space
	/// @details Represents a triangle defined by three vertices in 3D space.
	///          Provides side lengths, area calculation, triangle classification,
	///          centroid, and point containment testing using barycentric coordinates.
	class Triangle3D
	{
	protected:
		Pnt3Cart _pnt1, _pnt2, _pnt3;
	public:
		/// @brief Default constructor
		Triangle3D() {}
		
		/// @brief Constructs a triangle from three vertices
		Triangle3D(Pnt3Cart pnt1, Pnt3Cart pnt2, Pnt3Cart pnt3)
			: _pnt1(pnt1), _pnt2(pnt2), _pnt3(pnt3)
		{
		}

		/// @brief Gets side length A (from pnt1 to pnt2)
		Real A() const { return _pnt1.Dist(_pnt2); }
		/// @brief Gets side length B (from pnt2 to pnt3)
		Real B() const { return _pnt2.Dist(_pnt3); }
		/// @brief Gets side length C (from pnt3 to pnt1)
		Real C() const { return _pnt3.Dist(_pnt1); }

		/// @brief Gets vertex 1 (const)
		Pnt3Cart	Pnt1() const { return _pnt1; }
		/// @brief Gets vertex 1 (mutable)
		Pnt3Cart& Pnt1()			 { return _pnt1; }
		/// @brief Gets vertex 2 (const)
		Pnt3Cart	Pnt2() const { return _pnt2; }
		/// @brief Gets vertex 2 (mutable)
		Pnt3Cart& Pnt2()			 { return _pnt2; }
		/// @brief Gets vertex 3 (const)
		Pnt3Cart	Pnt3() const { return _pnt3; }
		/// @brief Gets vertex 3 (mutable)
		Pnt3Cart& Pnt3()			 { return _pnt3; }

		/// @brief Calculates the area using Heron's formula
		Real Area() const
		{
			Real s = (A() + B() + C()) / 2.0;
			return sqrt(s * (s - A()) * (s - B()) * (s - C()));
		}

		/// @brief Tests if this is a right triangle
		bool IsRight(Real eps = Defaults::Triangle3DIsRightTolerance) const
		{
			// Check if a^2 + b^2 = c^2 (Pythagorean theorem) with tolerance
			Real a = A(), b = B(), c = C();
			return (std::abs(a*a + b*b - c*c) < eps) || 
						(std::abs(a*a + c*c - b*b) < eps) || 
						(std::abs(b*b + c*c - a*a) < eps);
		}
		
		/// @brief Tests if this is an isosceles triangle (two equal sides)
		bool IsIsosceles(Real eps = Defaults::Triangle3DIsIsoscelesTolerance) const
		{
			Real a = A(), b = B(), c = C();
			return (std::abs(a - b) < eps) || (std::abs(a - c) < eps) || (std::abs(b - c) < eps);
		}
		
		/// @brief Tests if this is an equilateral triangle (all sides equal)
		bool IsEquilateral(Real eps = Defaults::Triangle3DIsEquilateralTolerance) const
		{
			Real a = A(), b = B(), c = C();
			return (std::abs(a - b) < eps) && (std::abs(a - c) < eps);
		}

		/// @brief Gets the centroid (center of mass)
		Pnt3Cart Centroid() const
		{
			return (_pnt1 + _pnt2 + _pnt3) / 3.0;
		}

		/// @brief Gets the plane containing this triangle
		Plane3D getDefinedPlane() const
		{
			return Plane3D(_pnt1, _pnt2, _pnt3);
		}

		/// @brief Tests if a point is inside the triangle
		/// @param pnt The point to test
		/// @param eps Tolerance for boundary cases
		/// @details Uses barycentric coordinates for containment testing
		bool IsPointInside(const Pnt3Cart& pnt, Real eps = Defaults::Triangle3DIsPointInsideTolerance) const
		{
			// check if point is inside triangle using barycentric coordinates
			Vec3Cart v0 = Vec3Cart(_pnt1, _pnt2);
			Vec3Cart v1 = Vec3Cart(_pnt1, _pnt3);
			Vec3Cart v2 = Vec3Cart(_pnt1, pnt);
			
			Real d00 = ScalarProduct(v0, v0);
			Real d01 = ScalarProduct(v0, v1);
			Real d11 = ScalarProduct(v1, v1);
			Real d20 = ScalarProduct(v2, v0);
			Real d21 = ScalarProduct(v2, v1);
			Real denom = d00 * d11 - d01 * d01;
			
			if (denom == 0.0)
				return false; // degenerate triangle
			
			Real v = (d11 * d20 - d01 * d21) / denom;
			
			if (v < 0.0 || v > 1.0)
				return false;
			
			Real w = (d00 * d21 - d01 * d20) / denom;
			
			if (w < 0.0 || w > 1.0)
				return false;
			
			return (v + w <= 1.0 + eps); // allow for small numerical errors
		}
	};

	/// @brief Triangular parametric surface in 3D
	/// @details Extends Triangle3D to implement IParametricSurface interface.
	///          Establishes a local coordinate system for parametric mapping.
	///          The longest side is used as the base for the local x-axis.
	class TriangleSurface3D : public Triangle3D, IParametricSurface<3>
	{
	public:
		Real _minX, _maxX, _minY, _maxY;
		Pnt3Cart _origin;
		Pnt3Cart _center;
		Vec3Cart _localX, _localY;
		Vec3Cart _normal;
		Real _pnt3XCoord;

		/// @brief Constructs a triangular surface from three vertices
		/// @details Automatically reorders vertices so the longest side becomes the base.
		/// @param pnt1 First vertex
		/// @param pnt2 Second vertex
		/// @param pnt3 Third vertex
		TriangleSurface3D(Pnt3Cart pnt1, Pnt3Cart pnt2, Pnt3Cart pnt3)
		{
			// CHECK THAT 1-2 side is THE LONGEST ONE!!!
			Real a = pnt1.Dist(pnt2);
			Real b = pnt2.Dist(pnt3);
			Real c = pnt3.Dist(pnt1);

			if (a >= b && a >= c)
			{
				; // nice, we are all good
			}
			else if (b >= a && b >= c)
			{
				// rotate points one place, so that 'b' side (pnt2-pnt3) is at the beginning
				Pnt3Cart tmp = pnt1;
				pnt1 = pnt2;
				pnt2 = pnt3;
				pnt3 = tmp;
			}
			else
			{
				// rotate points two places (c is longest, so pnt3-pnt1 becomes new pnt1-pnt2)
				Pnt3Cart tmp = pnt1;
				pnt1 = pnt3;
				pnt3 = pnt2;
				pnt2 = tmp;
			}

			_pnt1 = pnt1;		// initializing points in base Triangle3D class
			_pnt2 = pnt2;
			_pnt3 = pnt3;

			// calculate min and max
			// calculate center
			_center = (pnt1 + pnt2 + pnt3) / 3.0;

			// calculate local coordinate system
			// we will take x-axis to be from pnt1 to pnt2
			_localX = Vec3Cart(pnt1, pnt2).GetAsUnitVector();

			// we will calculate y-axis as following
			// calculate perpendicular vector to x-axis, that goes through pnt3
			Line3D	 lineX(pnt1, pnt2);
			Pnt3Cart pntOnLine = lineX.NearestPointOnLine(pnt3);

			_localY = Vec3Cart(pntOnLine, pnt3).GetAsUnitVector();

			_minX = 0.0;
			_maxX = pnt1.Dist(pnt2);
			_minY = 0.0;
			_maxY = pntOnLine.Dist(pnt3);
			_pnt3XCoord = pntOnLine.Dist(pnt1);

			_normal = VectorProduct(_localX, _localY).GetAsUnitVector();

			_origin = pnt1;
		}

		/// @brief Gets minimum u parameter
		virtual Real getMinU() const { return _minX; }
		/// @brief Gets maximum u parameter
		virtual Real getMaxU() const { return _maxX; }
		/// @brief Gets minimum w parameter for given u
		virtual Real getMinW(Real u) const { return _minY; }
		/// @brief Gets maximum w parameter for given u (varies with u due to triangle shape)
		virtual Real getMaxW(Real u) const
		{
			// this depends on value of u (which is localX)
			if (u < _pnt3XCoord)
				return _minY + (u / _pnt3XCoord) * (_maxY - _minY);
			else
				return _minY + ((_maxX - u) / (_maxX - _pnt3XCoord)) * (_maxY - _minY);
		}

		/// @brief Evaluates the surface at parameters (u, w)
		/// @param u First parameter (along local x-axis)
		/// @param w Second parameter (along local y-axis)
		/// @return 3D point on the surface
		VectorN<Real, 3> operator()(Real u, Real w) const
		{
			Pnt3Cart ret = _origin + u * _localX + w * _localY;
			return VectorN<Real, 3>({ ret.X(), ret.Y(), ret.Z() });
		}
	};

	/// @brief Rectangular parametric surface in 3D
	/// @details Represents a planar rectangular region in 3D space.
	///          Implements IParametricSurfaceRect interface for integration and analysis.
	///          All four corners must be coplanar.
	class RectSurface3D : public IParametricSurfaceRect<3>
	{
	public:
		// all constructors are expected to set properly these values
		Pnt3Cart _pnt1, _pnt2, _pnt3, _pnt4;
		Real _minX, _maxX, _minY, _maxY;
		Pnt3Cart _center;
		Vec3Cart _localX, _localY;

		/// @brief Default constructor
		RectSurface3D() {}
		
		/// @brief Constructs a rectangular surface from four coplanar corners
		/// @param pnt1 First corner
		/// @param pnt2 Second corner (along x-direction from pnt1)
		/// @param pnt3 Third corner (diagonal from pnt1)
		/// @param pnt4 Fourth corner (along y-direction from pnt1)
		/// @throws GeometryError if points are not coplanar
		RectSurface3D(Pnt3Cart pnt1, Pnt3Cart pnt2, Pnt3Cart pnt3, Pnt3Cart pnt4)
			: _pnt1(pnt1), _pnt2(pnt2), _pnt3(pnt3), _pnt4(pnt4)
		{
			// Validate coplanarity: all 4 points must lie in the same plane
			Plane3D plane(pnt1, pnt2, pnt3);
			if (!plane.IsPointOnPlane(pnt4))
				throw GeometryError("RectSurface3D ctor - points are not coplanar");

			// podesiti min i max
			Real lenX = _pnt1.Dist(_pnt2);
			Real lenY = _pnt1.Dist(_pnt4);
			_minX = -lenX / 2;
			_maxX = lenX / 2;
			_minY = -lenY / 2;
			_maxY = lenY / 2;

			// calculate center
			_center = (_pnt1 + _pnt2 + _pnt3 + _pnt4) / 4.0;

			// calculate local coordinate system
			_localX = Vec3Cart(_pnt1, _pnt2).GetAsUnitVector();
			_localY = Vec3Cart(_pnt1, _pnt4).GetAsUnitVector();
		}

		/// @brief Gets minimum u parameter
		virtual Real getMinU() const { return _minX; }
		/// @brief Gets maximum u parameter
		virtual Real getMaxU() const { return _maxX; }
		/// @brief Gets minimum w parameter
		virtual Real getMinW() const { return _minY; }
		/// @brief Gets maximum w parameter
		virtual Real getMaxW() const { return _maxY; }

		/// @brief Gets the surface normal vector
		Vec3Cart getNormal() const {
			return VectorProduct(Vec3Cart(_pnt1, _pnt2), Vec3Cart(_pnt1, _pnt4)).GetAsUnitVector();
		}
		
		/// @brief Gets the center point
		Pnt3Cart getCenter() const { return _center; }

		/// @brief Gets the area of the rectangle
		Real getArea() const {
			return VectorProduct(Vec3Cart(_pnt1, _pnt2), Vec3Cart(_pnt1, _pnt4)).NormL2();
		}
		
		/// @brief Evaluates the surface at parameters (u, w)
		/// @param u First parameter (along local x-axis)
		/// @param w Second parameter (along local y-axis)
		/// @return 3D point on the surface
		VectorN<Real, 3> operator()(Real u, Real w) const {
			Pnt3Cart ret = _center + u * _localX + w * _localY;
			return VectorN<Real, 3>({ ret.X(), ret.Y(), ret.Z() });
		}
	};

} // namespace MML

#endif // MML_GEOMETRY_3D_SURFACES_H

///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Geometry2DTriangle.h                                                ///
///  Description: 2D Triangle class with comprehensive geometric operations          ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                        ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_GEOMETRY_2D_TRIANGLE_H
#define MML_GEOMETRY_2D_TRIANGLE_H

#include <algorithm>
#include <cmath>

#include <mml/MMLBase.h>
#include <mml/base/Geometry/Geometry.h>
#include <mml/base/Geometry/Geometry2D/Geometry2DLines.h>

namespace MML
{
	/// @brief 2D triangle defined by three vertices
	/// @details Side lengths are computed from the current vertices without cached state.
	///          Provides comprehensive geometric computations including area (signed and unsigned),
	///          perimeter, centroid, circumcenter, incenter, circumradius, inradius, angles,
	///          point containment tests, and triangle classification (right, isosceles, equilateral).
	class Triangle2D
	{
	private:
		Point2Cartesian _pnt1, _pnt2, _pnt3;

	public:
		/// @brief Constructs a triangle from three vertices
		/// @param pnt1 First vertex
		/// @param pnt2 Second vertex
		/// @param pnt3 Third vertex
		Triangle2D(const Point2Cartesian& pnt1, const Point2Cartesian& pnt2, const Point2Cartesian& pnt3) 
			: _pnt1(pnt1), _pnt2(pnt2), _pnt3(pnt3)
		{ }

		/// @brief Gets first vertex (const version)
		Point2Cartesian  Pnt1() const { return _pnt1; }
		/// @brief Gets first vertex (mutable reference)
		Point2Cartesian& Pnt1() { return _pnt1; }
		/// @brief Gets second vertex (const version)
		Point2Cartesian  Pnt2() const { return _pnt2; }
		/// @brief Gets second vertex (mutable reference)
		Point2Cartesian& Pnt2() { return _pnt2; }
		/// @brief Gets third vertex (const version)
		Point2Cartesian  Pnt3() const { return _pnt3; }
		/// @brief Gets third vertex (mutable reference)
		Point2Cartesian& Pnt3() { return _pnt3; }

		/// @brief Gets length of side a (P1 to P2)
		Real A() const { return _pnt1.Dist(_pnt2); }
		/// @brief Gets length of side b (P2 to P3)
		Real B() const { return _pnt2.Dist(_pnt3); }
		/// @brief Gets length of side c (P3 to P1)
		Real C() const { return _pnt3.Dist(_pnt1); }

		/// @brief Calculates the perimeter
		Real Perimeter() const { return A() + B() + C(); }

		/// @brief Calculates the unsigned area from the coordinate determinant
		Real Area() const
		{
			return std::abs(SignedArea());
		}

		/// @brief Calculates the signed area (positive = CCW, negative = CW)
		Real SignedArea() const
		{
			Real area = (_pnt2.X() - _pnt1.X()) * (_pnt3.Y() - _pnt1.Y()) -
			           (_pnt3.X() - _pnt1.X()) * (_pnt2.Y() - _pnt1.Y());
			return area / MML::Real(2.0);
		}

		/// @brief Calculates the centroid (geometric center)
		/// @return Point at (P1 + P2 + P3) / 3
		Point2Cartesian Centroid() const
		{
			return Point2Cartesian((_pnt1.X() + _pnt2.X() + _pnt3.X()) / MML::Real(3.0),
			                      (_pnt1.Y() + _pnt2.Y() + _pnt3.Y()) / MML::Real(3.0));
		}

		/// @brief Calculates the circumcenter (center of circumscribed circle)
		/// @throws GeometryError if triangle is degenerate
		Point2Cartesian Circumcenter() const
		{
			Real D = MML::Real(2.0) * (_pnt1.X() * (_pnt2.Y() - _pnt3.Y()) +
			               _pnt2.X() * (_pnt3.Y() - _pnt1.Y()) +
			               _pnt3.X() * (_pnt1.Y() - _pnt2.Y()));
			
			if (std::abs(D) < Precision::GeometryEpsilon)
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

		/// @brief Calculates the circumradius (radius of circumscribed circle)
		/// @throws GeometryError if triangle is degenerate
		Real CircumRadius() const
		{
			Real a = A(), b = B(), c = C();
			Real area = Area();
			
			if (area < Precision::GeometryEpsilon)
				throw GeometryError("Triangle2D::CircumRadius - degenerate triangle");
			
			return (a * b * c) / (MML::Real(4.0) * area);
		}

		/// @brief Calculates the incenter (center of inscribed circle)
		/// @throws GeometryError if triangle is degenerate
		Point2Cartesian Incenter() const
		{
			Real a = A(), b = B(), c = C();
			Real perimeter = a + b + c;
			
			if (perimeter < Precision::GeometryEpsilon)
				throw GeometryError("Triangle2D::Incenter - degenerate triangle");
			
			Real x = (a * _pnt3.X() + b * _pnt1.X() + c * _pnt2.X()) / perimeter;
			Real y = (a * _pnt3.Y() + b * _pnt1.Y() + c * _pnt2.Y()) / perimeter;
			
			return Point2Cartesian(x, y);
		}

		/// @brief Calculates the inradius (radius of inscribed circle)
		/// @throws GeometryError if triangle is degenerate
		Real InRadius() const
		{
			Real area = Area();
			Real s = Perimeter() / MML::Real(2.0);
			
			if (s < Precision::GeometryEpsilon)
				throw GeometryError("Triangle2D::InRadius - degenerate triangle");
			
			return area / s;
		}

		/// @brief Calculates angle A at vertex 3 (opposite to side a) in radians
		Real AngleA() const
		{
			Real a = A(), b = B(), c = C();
			// Law of cosines: a² = b² + c² - 2bc·cos(A)
			Real cosA = (b * b + c * c - a * a) / (MML::Real(2.0) * b * c);
			return std::acos(std::clamp(cosA, static_cast<MML::Real>(-1.0), static_cast<MML::Real>(1.0)));
		}

		/// @brief Calculates angle B at vertex 1 (opposite to side b) in radians
		Real AngleB() const
		{
			Real a = A(), b = B(), c = C();
			Real cosB = (a * a + c * c - b * b) / (MML::Real(2.0) * a * c);
			return std::acos(std::clamp(cosB, static_cast<MML::Real>(-1.0), static_cast<MML::Real>(1.0)));
		}

		/// @brief Calculates angle C at vertex 2 (opposite to side c) in radians
		Real AngleC() const
		{
			Real a = A(), b = B(), c = C();
			Real cosC = (a * a + b * b - c * c) / (MML::Real(2.0) * a * b);
			return std::acos(std::clamp(cosC, static_cast<MML::Real>(-1.0), static_cast<MML::Real>(1.0)));
		}

		/// @brief Tests if a point is inside the triangle using barycentric coordinates
		/// @param p The point to test
		/// @param epsilon Tolerance for boundary cases
		/// @return True if point is inside or on the boundary
		bool Contains(const Point2Cartesian& p, Real epsilon = Precision::GeometryEpsilon) const
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
			
			Real c = MML::Real(1.0) - a - b;
			
			// Point is inside if all barycentric coordinates are non-negative
			return (a >= -epsilon) && (b >= -epsilon) && (c >= -epsilon);
		}

		/// @brief Tests if the triangle is a right triangle
		/// @param epsilon Tolerance for the test
		bool IsRight(Real epsilon = Precision::GeometryEpsilon) const
		{
			Real a = A(), b = B(), c = C();
			Real aSq = a * a, bSq = b * b, cSq = c * c;
			
			return AreEqual(aSq + bSq, cSq, epsilon) ||
			       AreEqual(aSq + cSq, bSq, epsilon) ||
			       AreEqual(bSq + cSq, aSq, epsilon);
		}

		/// @brief Tests if the triangle is isosceles (at least two equal sides)
		/// @param epsilon Tolerance for the test
		bool IsIsosceles(Real epsilon = Precision::GeometryEpsilon) const
		{
			Real a = A(), b = B(), c = C();
			return AreEqual(a, b, epsilon) || AreEqual(a, c, epsilon) || AreEqual(b, c, epsilon);
		}

		/// @brief Tests if the triangle is equilateral (all sides equal)
		/// @param epsilon Tolerance for the test
		bool IsEquilateral(Real epsilon = Precision::GeometryEpsilon) const
		{
			Real a = A(), b = B(), c = C();
			return AreEqual(a, b, epsilon) && AreEqual(a, c, epsilon);
		}

		/// @brief Tests if vertices are ordered counter-clockwise
		bool IsCounterClockwise() const
		{
			return SignedArea() > 0.0;
		}

		/// @brief Tests if vertices are ordered clockwise
		bool IsClockwise() const
		{
			return SignedArea() < 0.0;
		}
	};

} // namespace MML

#endif // MML_GEOMETRY_2D_TRIANGLE_H

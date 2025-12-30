///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Geometry.h                                                          ///
///  Description: Basic geometric primitives (Point2D, Point3D, Line, Plane)          ///
///               Coordinate representations and geometric operations                 ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_GEOMETRY_H
#define MML_GEOMETRY_H

#include "MMLBase.h"


namespace MML
{
	class Point2Cartesian
	{
	private:
		Real _x, _y;

	public:
		Real  X() const { return _x; }
		Real& X()				{ return _x; }
		Real  Y() const { return _y; }
		Real& Y()				{ return _y; }

		Point2Cartesian() : _x(0), _y(0) {}
		Point2Cartesian(Real x, Real y) : _x(x), _y(y) {}

		Real Dist(const Point2Cartesian& b) const { return sqrt(POW2(b.X() - X()) + POW2(b.Y() - Y())); }

		bool	operator==(const Point2Cartesian& b) const { return (X() == b.X()) && (Y() == b.Y()); }
		bool	operator!=(const Point2Cartesian& b) const { return (X() != b.X()) || (Y() != b.Y()); }
		bool  IsEqual(const Point2Cartesian& b, Real eps = Defaults::Pnt2CartIsEqualTolerance) const { return Dist(b) < eps; }

		Point2Cartesian  operator+(const Point2Cartesian& b) const { return Point2Cartesian(X() + b.X(), Y() + b.Y()); }
		Point2Cartesian& operator+=(const Point2Cartesian& b)	{	_x += b.X(); _y += b.Y(); return *this;	}

		Point2Cartesian  operator-(const Point2Cartesian& b) const { return Point2Cartesian(X() - b.X(), Y() - b.Y()); }
		Point2Cartesian& operator-=(const Point2Cartesian& b) { _x -= b.X(); _y -= b.Y(); return *this; }

		Point2Cartesian  operator*(Real b) const { return Point2Cartesian(X() * b, Y() * b); }
		Point2Cartesian& operator*=(Real b) { _x *= b; _y *= b; return *this; }

		Point2Cartesian	 operator/(Real b) const { return Point2Cartesian(X() / b, Y() / b); }
		Point2Cartesian& operator/=(Real b) { _x /= b; _y /= b; return *this; }

		friend Point2Cartesian operator*(Real a, const Point2Cartesian& b) { return Point2Cartesian(a * b.X(), a * b.Y()); }
	};	

	class Point2Polar
	{
	private:
		Real _r, _phi;

	public:
		Real  R() const		{ return _r; }
		Real& R()					{ return _r; }
		Real  Phi() const { return _phi; }
		Real& Phi()				{ return _phi; }

		Point2Polar() : _r(0), _phi(0) {}
		Point2Polar(Real r, Real phi) : _r(r), _phi(phi) {}
		Point2Polar(const Point2Cartesian& pnt) 
		{ 
			_r = sqrt(POW2(pnt.X()) + POW2(pnt.Y()));
			_phi = atan2(pnt.Y(), pnt.X());
		}

		bool	operator==(const Point2Polar& b) const { return (R() == b.R()) && (Phi() == b.Phi()); }
		bool	operator!=(const Point2Polar& b) const { return (R() != b.R()) || (Phi() != b.Phi()); }
		bool  IsEqual(const Point2Polar& b, Real eps = Defaults::Pnt2PolarIsEqualTolerance) const { return Dist(b) < eps; }

		Point2Cartesian TransfToCart() const {  return Point2Cartesian(R() * cos(Phi()), R() * sin(Phi())); }

		Real Dist(const Point2Polar& b) const { return sqrt(R() * R() + b.R() * b.R() - 2 * R() * b.R() * cos(b.Phi() - Phi())); }
	};

	class Point3Cartesian
	{
	private:
		Real _x, _y, _z;

	public:
		Real  X() const { return _x; }
		Real& X()				{ return _x; }
		Real  Y() const { return _y; }
		Real& Y()				{ return _y; }
		Real  Z() const { return _z; }
		Real& Z()				{ return _z; }

		Point3Cartesian() : _x(0), _y(0), _z(0) {}
		Point3Cartesian(Real x, Real y, Real z) : _x(x), _y(y), _z(z) {}

		Real Dist(const Point3Cartesian& b) const
		{ 
			return sqrt(POW2(b.X() - X()) + POW2(b.Y() - Y()) + POW2(b.Z() - Z())); 
		}

		bool	operator==(const Point3Cartesian& b) const {return (X() == b.X()) && (Y() == b.Y()) && (Z() == b.Z()); }
		bool	operator!=(const Point3Cartesian& b) const {return (X() != b.X()) || (Y() != b.Y()) || (Z() != b.Z()); }
		bool  IsEqual(const Point3Cartesian& b, Real eps = Defaults::Pnt3CartIsEqualTolerance) const { return Dist(b) < eps; }

		Point3Cartesian operator+(const Point3Cartesian& b) const { return Point3Cartesian(X() + b.X(), Y() + b.Y(), Z() + b.Z()); }
		Point3Cartesian& operator+=(const Point3Cartesian& b) { _x += b.X(); _y += b.Y(); _z += b.Z(); return *this; }

		Point3Cartesian operator-(const Point3Cartesian& b) const { return Point3Cartesian(X() - b.X(), Y() - b.Y(), Z() - b.Z()); }
		Point3Cartesian& operator-=(const Point3Cartesian& b) { _x -= b.X(); _y -= b.Y(); _z -= b.Z(); return *this; }

		Point3Cartesian operator*(Real b) const { return Point3Cartesian(X() * b, Y() * b, Z() * b); }
		Point3Cartesian& operator*=(Real b) { _x *= b; _y *= b; _z *= b; return *this; }

		Point3Cartesian operator/(Real b) const { return Point3Cartesian(X() / b, Y() / b, Z() / b); }
		Point3Cartesian& operator/=(Real b) { _x /= b; _y /= b; _z /= b; return *this; }

		friend Point3Cartesian operator*(Real a, const Point3Cartesian& b) { return Point3Cartesian(a * b.X(), a * b.Y(), a * b.Z()); }
	};

	class Point3Spherical
	{
	private:
		Real _r, _theta, _phi;

	public:
		Real  R() const			{ return _r; }
		Real& R()						{ return _r; }
		Real  Theta() const { return _theta; }
		Real& Theta()				{ return _theta; }
		Real  Phi() const		{ return _phi; }
		Real& Phi()					{ return _phi; }

		Point3Spherical() : _r(0), _theta(0), _phi(0) {}
		Point3Spherical(Real r, Real theta, Real phi) : _r(r), _theta(theta), _phi(phi) {}
		Point3Spherical(const Point3Cartesian& pnt)
		{
			_r = sqrt(POW2(pnt.X()) + POW2(pnt.Y()) + POW2(pnt.Z()));
			_theta = atan2(sqrt(POW2(pnt.X()) + POW2(pnt.Y())), pnt.Z());
			_phi = atan2(pnt.Y(), pnt.X());
		}
		
		bool	operator==(const Point3Spherical& b) const { return (R() == b.R()) && (Theta() == b.Theta()) && (Phi() == b.Phi()); }
		bool	operator!=(const Point3Spherical& b) const { return (R() != b.R()) || (Theta() != b.Theta()) || (Phi() != b.Phi()); }
		bool  IsEqual(const Point3Spherical& b, Real eps = Defaults::Pnt3SphIsEqualTolerance) const { return Dist(b) < eps; }

		Point3Cartesian TransfToCart() const
		{
			return Point3Cartesian(R() * sin(Theta()) * cos(Phi()), R() * sin(Theta()) * sin(Phi()), R() * cos(Theta()));
		}
		
		Real Dist(const Point3Spherical& b) const
		{
			return sqrt(R() * R() + b.R() * b.R() - 2 * R() * b.R() * cos(b.Theta() - Theta()) * cos(b.Phi() - Phi()));
		}
	};

	class Point3Cylindrical
	{
	private:
		Real _r, _phi, _z;

	public:
		Real  R() const		{ return _r; }
		Real& R()					{ return _r; }
		Real  Phi() const { return _phi; }
		Real& Phi()				{ return _phi; }
		Real  Z() const		{ return _z; }
		Real& Z()					{ return _z; }

		Point3Cylindrical() : _r(0), _phi(0), _z(0) {}
		Point3Cylindrical(Real r, Real phi, Real z) : _r(r), _phi(phi), _z(z) {}
		Point3Cylindrical(const Point3Cartesian& pnt)
		{
			_r = sqrt(POW2(pnt.X()) + POW2(pnt.Y()));
			_phi = atan2(pnt.Y(), pnt.X());
			_z = pnt.Z();
		}
		
		bool	operator==(const Point3Cylindrical& b) const { return (R() == b.R()) && (Phi() == b.Phi()) && (Z() == b.Z()); }
		bool	operator!=(const Point3Cylindrical& b) const { return (R() != b.R()) || (Phi() != b.Phi()) || (Z() != b.Z()); }
		bool  IsEqual(const Point3Cylindrical& b, Real eps = Defaults::Pnt3CylIsEqualTolerance) const { return Dist(b) < eps; }

		Point3Cartesian TransfToCart() const
		{
			return Point3Cartesian(R() * cos(Phi()), R() * sin(Phi()), Z());
		}
		
		Real Dist(const Point3Cylindrical& b) const
		{
			return sqrt(R() * R() + b.R() * b.R() - 2 * R() * b.R() * cos(b.Phi() - Phi()) + POW2(b.Z() - Z()));
		}
	};

	class Triangle
	{
	private:
		Real _a, _b, _c;

	public:
		Real  A() const { return _a; }
		Real& A()				{ return _a; }
		Real  B() const { return _b; }
		Real& B()				{ return _b; }
		Real  C() const { return _c; }
		Real& C()				{ return _c; }

		Triangle() : _a(0.0), _b(0.0), _c(0.0){}
		Triangle(Real a, Real b, Real c) : _a(a), _b(b), _c(c) {}

		Real Area() const
		{
			Real s = (A() + B() + C()) / 2.0;
			return sqrt(s * (s - A()) * (s - B()) * (s - C()));
		}
		
		bool IsRight() const
		{
			return (hypot(A(), B()) == C() || hypot(A(), C()) == B() || hypot(B(), C()) == A());
		}
		bool IsIsosceles() const		// two sides are the same length
		{
			return (A() == B()) || (A() == C()) || (B() == C());
		}
		bool IsEquilateral() const	// all sides are the same length
		{
			return (A() == B()) && (A() == C());
		}
	};

	typedef Point2Cartesian		Pnt2Cart;
	typedef Point2Polar				Pnt2Pol;
	typedef Point3Cartesian		Pnt3Cart;
	typedef Point3Spherical		Pnt3Sph;
	typedef Point3Cylindrical Pnt3Cyl;
}
#endif

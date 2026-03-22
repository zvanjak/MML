///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        GeometryPoints.h                                                    ///
///  Description: Point classes in 2D and 3D coordinate systems                       ///
///               (Cartesian, Polar, Spherical, Cylindrical)                          ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_GEOMETRY_POINTS_H
#define MML_GEOMETRY_POINTS_H

#include <cmath>

#include "mml/MMLBase.h"


namespace MML {
	/// @brief 2D point in Cartesian coordinates (x, y)
	/// @details Represents a point in 2D Euclidean space using Cartesian coordinates.
	///          Supports basic arithmetic operations (addition, subtraction, scalar multiplication/division),
	///          distance calculations, and equality comparisons with configurable tolerance.
	class Point2Cartesian {
	private:
		Real _x, _y;

	public:
		/// @brief Gets the x-coordinate (const version)
		Real X() const { return _x; }
		/// @brief Gets the x-coordinate (mutable reference)
		Real& X() { return _x; }
		/// @brief Gets the y-coordinate (const version)
		Real Y() const { return _y; }
		/// @brief Gets the y-coordinate (mutable reference)
		Real& Y() { return _y; }

		/// @brief Default constructor - initializes to origin (0, 0)
		Point2Cartesian()
			: _x(0)
			, _y(0) {}
		/// @brief Constructs a point from Cartesian coordinates
		/// @param x The x-coordinate
		/// @param y The y-coordinate
		Point2Cartesian(Real x, Real y)
			: _x(x)
			, _y(y) {}

		/// @brief Computes Euclidean distance to another point
		/// @param b The other point
		/// @return Distance between this point and b
		Real Dist(const Point2Cartesian& b) const { return sqrt(POW2(b.X() - X()) + POW2(b.Y() - Y())); }

		/// @brief Exact equality comparison
		bool operator==(const Point2Cartesian& b) const { return (X() == b.X()) && (Y() == b.Y()); }
		/// @brief Exact inequality comparison
		bool operator!=(const Point2Cartesian& b) const { return (X() != b.X()) || (Y() != b.Y()); }
		/// @brief Approximate equality comparison with tolerance
		/// @param b The other point
		/// @param eps Tolerance for comparison (default: Defaults::Pnt2CartIsEqualTolerance)
		/// @return True if distance between points is less than eps
		bool IsEqualTo(const Point2Cartesian& b, Real eps = Defaults::Pnt2CartIsEqualTolerance) const { return Dist(b) < eps; }

		/// @brief Vector addition of two points
		Point2Cartesian operator+(const Point2Cartesian& b) const { return Point2Cartesian(X() + b.X(), Y() + b.Y()); }
		/// @brief In-place vector addition
		Point2Cartesian& operator+=(const Point2Cartesian& b) {
			_x += b.X();
			_y += b.Y();
			return *this;
		}

		/// @brief Vector subtraction of two points
		Point2Cartesian operator-(const Point2Cartesian& b) const { return Point2Cartesian(X() - b.X(), Y() - b.Y()); }
		/// @brief In-place vector subtraction
		Point2Cartesian& operator-=(const Point2Cartesian& b) {
			_x -= b.X();
			_y -= b.Y();
			return *this;
		}

		/// @brief Scalar multiplication (point * scalar)
		Point2Cartesian operator*(Real b) const { return Point2Cartesian(X() * b, Y() * b); }
		/// @brief In-place scalar multiplication
		Point2Cartesian& operator*=(Real b) {
			_x *= b;
			_y *= b;
			return *this;
		}

		/// @brief Scalar division
		Point2Cartesian operator/(Real b) const { return Point2Cartesian(X() / b, Y() / b); }
		/// @brief In-place scalar division
		Point2Cartesian& operator/=(Real b) {
			_x /= b;
			_y /= b;
			return *this;
		}

		/// @brief Scalar multiplication (scalar * point)
		friend Point2Cartesian operator*(Real a, const Point2Cartesian& b) { return Point2Cartesian(a * b.X(), a * b.Y()); }
	};

	/// @brief 2D point in polar coordinates (r, φ)
	/// @details Represents a point in 2D space using polar coordinates where r is the radial
	///          distance from origin and φ (phi) is the angle from the positive x-axis in radians.
	///          Supports conversion to/from Cartesian coordinates and distance calculations.
	class Point2Polar {
	private:
		Real _r, _phi;

	public:
		/// @brief Gets the radial coordinate (const version)
		Real R() const { return _r; }
		/// @brief Gets the radial coordinate (mutable reference)
		Real& R() { return _r; }
		/// @brief Gets the angular coordinate φ in radians (const version)
		Real Phi() const { return _phi; }
		/// @brief Gets the angular coordinate φ in radians (mutable reference)
		Real& Phi() { return _phi; }

		/// @brief Default constructor - initializes to origin (0, 0)
		Point2Polar()
			: _r(0)
			, _phi(0) {}
		/// @brief Constructs a point from polar coordinates
		/// @param r Radial distance from origin
		/// @param phi Angle from positive x-axis (in radians)
		Point2Polar(Real r, Real phi)
			: _r(r)
			, _phi(phi) {}
		/// @brief Constructs a polar point from a Cartesian point
		/// @param pnt The Cartesian point to convert
		Point2Polar(const Point2Cartesian& pnt) {
			_r = sqrt(POW2(pnt.X()) + POW2(pnt.Y()));
			_phi = atan2(pnt.Y(), pnt.X());
		}

		/// @brief Exact equality comparison
		bool operator==(const Point2Polar& b) const { return (R() == b.R()) && (Phi() == b.Phi()); }
		/// @brief Exact inequality comparison
		bool operator!=(const Point2Polar& b) const { return (R() != b.R()) || (Phi() != b.Phi()); }
		/// @brief Approximate equality comparison with tolerance
		/// @param b The other point
		/// @param eps Tolerance for comparison (default: Defaults::Pnt2PolarIsEqualTolerance)
		/// @return True if distance between points is less than eps
		bool IsEqualTo(const Point2Polar& b, Real eps = Defaults::Pnt2PolarIsEqualTolerance) const { return Dist(b) < eps; }

		/// @brief Converts this polar point to Cartesian coordinates
		/// @return Equivalent Point2Cartesian (x = r*cos(φ), y = r*sin(φ))
		Point2Cartesian TransfToCart() const { return Point2Cartesian(R() * cos(Phi()), R() * sin(Phi())); }

		/// @brief Computes distance to another polar point using law of cosines
		/// @param b The other point
		/// @return Euclidean distance between the two points
		Real Dist(const Point2Polar& b) const { return sqrt(R() * R() + b.R() * b.R() - 2 * R() * b.R() * cos(b.Phi() - Phi())); }
	};

	/// @brief 3D point in Cartesian coordinates (x, y, z)
	/// @details Represents a point in 3D Euclidean space using Cartesian coordinates.
	///          Supports basic arithmetic operations (addition, subtraction, scalar multiplication/division),
	///          distance calculations, and equality comparisons with configurable tolerance.
	class Point3Cartesian {
	private:
		Real _x, _y, _z;

	public:
		/// @brief Gets the x-coordinate (const version)
		Real X() const { return _x; }
		/// @brief Gets the x-coordinate (mutable reference)
		Real& X() { return _x; }
		/// @brief Gets the y-coordinate (const version)
		Real Y() const { return _y; }
		/// @brief Gets the y-coordinate (mutable reference)
		Real& Y() { return _y; }
		/// @brief Gets the z-coordinate (const version)
		Real Z() const { return _z; }
		/// @brief Gets the z-coordinate (mutable reference)
		Real& Z() { return _z; }

		/// @brief Default constructor - initializes to origin (0, 0, 0)
		Point3Cartesian()
			: _x(0)
			, _y(0)
			, _z(0) {}
		/// @brief Constructs a point from Cartesian coordinates
		/// @param x The x-coordinate
		/// @param y The y-coordinate
		/// @param z The z-coordinate
		Point3Cartesian(Real x, Real y, Real z)
			: _x(x)
			, _y(y)
			, _z(z) {}

		/// @brief Computes Euclidean distance to another point
		/// @param b The other point
		/// @return Distance between this point and b
		Real Dist(const Point3Cartesian& b) const { return sqrt(POW2(b.X() - X()) + POW2(b.Y() - Y()) + POW2(b.Z() - Z())); }

		/// @brief Exact equality comparison
		bool operator==(const Point3Cartesian& b) const { return (X() == b.X()) && (Y() == b.Y()) && (Z() == b.Z()); }
		/// @brief Exact inequality comparison
		bool operator!=(const Point3Cartesian& b) const { return (X() != b.X()) || (Y() != b.Y()) || (Z() != b.Z()); }
		/// @brief Approximate equality comparison with tolerance
		/// @param b The other point
		/// @param eps Tolerance for comparison (default: Defaults::Pnt3CartIsEqualTolerance)
		/// @return True if distance between points is less than eps
		bool IsEqualTo(const Point3Cartesian& b, Real eps = Defaults::Pnt3CartIsEqualTolerance) const { return Dist(b) < eps; }

		/// @brief Vector addition of two points
		Point3Cartesian operator+(const Point3Cartesian& b) const { return Point3Cartesian(X() + b.X(), Y() + b.Y(), Z() + b.Z()); }
		/// @brief In-place vector addition
		Point3Cartesian& operator+=(const Point3Cartesian& b) {
			_x += b.X();
			_y += b.Y();
			_z += b.Z();
			return *this;
		}

		/// @brief Vector subtraction of two points
		Point3Cartesian operator-(const Point3Cartesian& b) const { return Point3Cartesian(X() - b.X(), Y() - b.Y(), Z() - b.Z()); }
		/// @brief In-place vector subtraction
		Point3Cartesian& operator-=(const Point3Cartesian& b) {
			_x -= b.X();
			_y -= b.Y();
			_z -= b.Z();
			return *this;
		}

		/// @brief Scalar multiplication (point * scalar)
		Point3Cartesian operator*(Real b) const { return Point3Cartesian(X() * b, Y() * b, Z() * b); }
		/// @brief In-place scalar multiplication
		Point3Cartesian& operator*=(Real b) {
			_x *= b;
			_y *= b;
			_z *= b;
			return *this;
		}

		/// @brief Scalar division
		Point3Cartesian operator/(Real b) const { return Point3Cartesian(X() / b, Y() / b, Z() / b); }
		/// @brief In-place scalar division
		Point3Cartesian& operator/=(Real b) {
			_x /= b;
			_y /= b;
			_z /= b;
			return *this;
		}

		/// @brief Scalar multiplication (scalar * point)
		friend Point3Cartesian operator*(Real a, const Point3Cartesian& b) { return Point3Cartesian(a * b.X(), a * b.Y(), a * b.Z()); }
	};

	/// @brief 3D point in spherical coordinates (r, θ, φ)
	/// @details Represents a point in 3D space using spherical coordinates where:
	///          - r: radial distance from origin
	///          - θ (theta): polar angle from positive z-axis (0 to π)
	///          - φ (phi): azimuthal angle from positive x-axis in xy-plane (0 to 2π)
	///          Uses physics convention (ISO 31-11). Supports conversion to Cartesian coordinates.
	class Point3Spherical {
	private:
		Real _r, _theta, _phi;

	public:
		/// @brief Gets the radial coordinate (const version)
		Real R() const { return _r; }
		/// @brief Gets the radial coordinate (mutable reference)
		Real& R() { return _r; }
		/// @brief Gets the polar angle θ in radians (const version)
		Real Theta() const { return _theta; }
		/// @brief Gets the polar angle θ in radians (mutable reference)
		Real& Theta() { return _theta; }
		/// @brief Gets the azimuthal angle φ in radians (const version)
		Real Phi() const { return _phi; }
		/// @brief Gets the azimuthal angle φ in radians (mutable reference)
		Real& Phi() { return _phi; }

		/// @brief Default constructor - initializes to origin
		Point3Spherical()
			: _r(0)
			, _theta(0)
			, _phi(0) {}
		/// @brief Constructs a point from spherical coordinates
		/// @param r Radial distance from origin
		/// @param theta Polar angle from z-axis (in radians, 0 to π)
		/// @param phi Azimuthal angle in xy-plane (in radians, 0 to 2π)
		Point3Spherical(Real r, Real theta, Real phi)
			: _r(r)
			, _theta(theta)
			, _phi(phi) {}
		/// @brief Constructs a spherical point from a Cartesian point
		/// @param pnt The Cartesian point to convert
		Point3Spherical(const Point3Cartesian& pnt) {
			_r = sqrt(POW2(pnt.X()) + POW2(pnt.Y()) + POW2(pnt.Z()));
			_theta = atan2(sqrt(POW2(pnt.X()) + POW2(pnt.Y())), pnt.Z());
			_phi = atan2(pnt.Y(), pnt.X());
		}

		/// @brief Exact equality comparison
		bool operator==(const Point3Spherical& b) const { return (R() == b.R()) && (Theta() == b.Theta()) && (Phi() == b.Phi()); }
		/// @brief Exact inequality comparison
		bool operator!=(const Point3Spherical& b) const { return (R() != b.R()) || (Theta() != b.Theta()) || (Phi() != b.Phi()); }
		/// @brief Approximate equality comparison with tolerance
		/// @param b The other point
		/// @param eps Tolerance for comparison (default: Defaults::Pnt3SphIsEqualTolerance)
		/// @return True if distance between points is less than eps
		bool IsEqualTo(const Point3Spherical& b, Real eps = Defaults::Pnt3SphIsEqualTolerance) const { return Dist(b) < eps; }

		/// @brief Converts this spherical point to Cartesian coordinates
		/// @return Equivalent Point3Cartesian using spherical-to-Cartesian transformation
		Point3Cartesian TransfToCart() const {
			return Point3Cartesian(R() * sin(Theta()) * cos(Phi()), R() * sin(Theta()) * sin(Phi()), R() * cos(Theta()));
		}

		/// @brief Computes distance to another spherical point
		/// @param b The other point
		/// @return Euclidean distance between the two points (using spherical law of cosines)
		Real Dist(const Point3Spherical& b) const {
			Real cosGamma = cos(Theta()) * cos(b.Theta()) + sin(Theta()) * sin(b.Theta()) * cos(b.Phi() - Phi());
			return sqrt(R() * R() + b.R() * b.R() - 2 * R() * b.R() * cosGamma);
		}
	};

	/// @brief 3D point in cylindrical coordinates (r, φ, z)
	/// @details Represents a point in 3D space using cylindrical coordinates where:
	///          - r: radial distance from z-axis
	///          - φ (phi): azimuthal angle from positive x-axis in xy-plane (in radians)
	///          - z: height along z-axis
	///          Supports conversion to Cartesian coordinates and distance calculations.
	class Point3Cylindrical {
	private:
		Real _r, _phi, _z;

	public:
		/// @brief Gets the radial coordinate (const version)
		Real R() const { return _r; }
		/// @brief Gets the radial coordinate (mutable reference)
		Real& R() { return _r; }
		/// @brief Gets the azimuthal angle φ in radians (const version)
		Real Phi() const { return _phi; }
		/// @brief Gets the azimuthal angle φ in radians (mutable reference)
		Real& Phi() { return _phi; }
		/// @brief Gets the z-coordinate (const version)
		Real Z() const { return _z; }
		/// @brief Gets the z-coordinate (mutable reference)
		Real& Z() { return _z; }

		/// @brief Default constructor - initializes to origin
		Point3Cylindrical()
			: _r(0)
			, _phi(0)
			, _z(0) {}
		/// @brief Constructs a point from cylindrical coordinates
		/// @param r Radial distance from z-axis
		/// @param phi Azimuthal angle from x-axis (in radians)
		/// @param z Height along z-axis
		Point3Cylindrical(Real r, Real phi, Real z)
			: _r(r)
			, _phi(phi)
			, _z(z) {}
		/// @brief Constructs a cylindrical point from a Cartesian point
		/// @param pnt The Cartesian point to convert
		Point3Cylindrical(const Point3Cartesian& pnt) {
			_r = sqrt(POW2(pnt.X()) + POW2(pnt.Y()));
			_phi = atan2(pnt.Y(), pnt.X());
			_z = pnt.Z();
		}

		/// @brief Exact equality comparison
		bool operator==(const Point3Cylindrical& b) const { return (R() == b.R()) && (Phi() == b.Phi()) && (Z() == b.Z()); }
		/// @brief Exact inequality comparison
		bool operator!=(const Point3Cylindrical& b) const { return (R() != b.R()) || (Phi() != b.Phi()) || (Z() != b.Z()); }
		/// @brief Approximate equality comparison with tolerance
		/// @param b The other point
		/// @param eps Tolerance for comparison (default: Defaults::Pnt3CylIsEqualTolerance)
		/// @return True if distance between points is less than eps
		bool IsEqualTo(const Point3Cylindrical& b, Real eps = Defaults::Pnt3CylIsEqualTolerance) const { return Dist(b) < eps; }

		/// @brief Converts this cylindrical point to Cartesian coordinates
		/// @return Equivalent Point3Cartesian using cylindrical-to-Cartesian transformation
		Point3Cartesian TransfToCart() const { return Point3Cartesian(R() * cos(Phi()), R() * sin(Phi()), Z()); }

		/// @brief Computes distance to another cylindrical point
		/// @param b The other point
		/// @return Euclidean distance between the two points
		Real Dist(const Point3Cylindrical& b) const {
			return sqrt(R() * R() + b.R() * b.R() - 2 * R() * b.R() * cos(b.Phi() - Phi()) + POW2(b.Z() - Z()));
		}
	};

	//=============================================================================
	// TYPE ALIASES FOR POINTS
	//=============================================================================

	/// @brief Alias for Point2Cartesian
	typedef Point2Cartesian Pnt2Cart;
	/// @brief Alias for Point2Polar
	typedef Point2Polar Pnt2Pol;
	/// @brief Alias for Point3Cartesian
	typedef Point3Cartesian Pnt3Cart;
	/// @brief Alias for Point3Spherical
	typedef Point3Spherical Pnt3Sph;
	/// @brief Alias for Point3Cylindrical
	typedef Point3Cylindrical Pnt3Cyl;

} // namespace MML
#endif

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

#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <limits>

#include <mml/MMLBase.h>


namespace MML {
	namespace detail {
		inline void ValidateGeometryCoordinates(std::initializer_list<Real> coordinates, const char* context) {
			for (Real coordinate : coordinates)
				if (!std::isfinite(coordinate))
					throw GeometryError(std::string(context) + " - finite coordinates required");
		}

		inline Real CheckedGeometryDistance(Real distance, const char* context, bool nonzero = false) {
			if (!std::isfinite(distance) || (distance == 0 && nonzero))
				throw GeometryError(std::string(context) + " - distance is not representable");
			return distance;
		}

		inline Real GeometryAngleDifferenceError(Real first, Real second, Real difference) {
			// Knuth/Shewchuk TwoDiff: retain the low part of a finite angular subtraction.
			const Real virtualSecond = first - difference;
			return (first - (difference + virtualSecond)) + (virtualSecond - second);
		}

		inline Real GeometrySineChord(Real first, Real second) {
			const Real difference = first - second;
			if (std::abs(difference) >= 3 * Constants::PI / 2 && std::abs(difference) <= 5 * Constants::PI / 2) {
				const Real residual = (difference - std::copysign(2 * Constants::PI, difference))
					+ GeometryAngleDifferenceError(first, second, difference);
				return -GeometrySineChord(residual, 0);
			}
			// Retain subnormal angular differences before multiplying by a large radius.
			const Real underflow = std::numeric_limits<Real>::min();
			if (std::abs(difference) < underflow) return difference;
			const Real half = std::isfinite(difference) ? difference / 2 : first / 2 - second / 2;
			return 2 * std::sin(half);
		}

		inline Real GeometryCosineChord(Real first, Real second) {
			const Real difference = first - second;
			if (std::abs(difference) >= Constants::PI / 2 && std::abs(difference) <= 3 * Constants::PI / 2) {
				const Real residual = (difference - std::copysign(Constants::PI, difference))
					+ GeometryAngleDifferenceError(first, second, difference);
				return (difference < 0 ? 1 : -1) * GeometrySineChord(residual, 0);
			}
			const Real half = std::isfinite(difference) ? difference / 2 : first / 2 - second / 2;
			return 2 * std::cos(half);
		}

		inline Real GeometryAngularComponent(Real firstRadius, Real secondRadius, Real chord) {
			return (std::sqrt(std::max(firstRadius, secondRadius)) * chord)
				* std::sqrt(std::min(firstRadius, secondRadius));
		}

		inline Real GeometryRadialDistance(Real radius1, Real angle1, Real radius2, Real angle2,
			Real height, const char* context) {
			ValidateGeometryCoordinates({radius1, angle1, radius2, angle2}, context);
			const Real first = std::abs(radius1), second = std::abs(radius2);
			const Real chord = std::signbit(radius1) != std::signbit(radius2)
				? GeometryCosineChord(angle1, angle2) : GeometrySineChord(angle1, angle2);
			const Real radial = first - second;
			const Real angular = GeometryAngularComponent(first, second, chord);
			return CheckedGeometryDistance(std::hypot(radial, angular, height), context,
				height != 0 || radial != 0 || (first != 0 && second != 0 && chord != 0));
		}
	}

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
		Real Dist(const Point2Cartesian& b) const {
			detail::ValidateGeometryCoordinates({X(), Y(), b.X(), b.Y()}, "Point2Cartesian::Dist");
			return detail::CheckedGeometryDistance(std::hypot(b.X() - X(), b.Y() - Y()), "Point2Cartesian::Dist");
		}

		/// @brief Exact equality comparison
		bool operator==(const Point2Cartesian& b) const { return (X() == b.X()) && (Y() == b.Y()); }
		/// @brief Exact inequality comparison
		bool operator!=(const Point2Cartesian& b) const { return (X() != b.X()) || (Y() != b.Y()); }
		/// @brief Approximate equality comparison with tolerance
		/// @param b The other point
		/// @param eps Tolerance for comparison (default: Precision::Pnt2CartIsEqualTolerance)
		/// @return True if distance between points is less than eps
		bool IsEqualTo(const Point2Cartesian& b, Real eps = Precision::Pnt2CartIsEqualTolerance) const { return Dist(b) < eps; }

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
			_r = Point2Cartesian().Dist(pnt);
			_phi = atan2(pnt.Y(), pnt.X());
		}

		/// @brief Exact equality comparison
		bool operator==(const Point2Polar& b) const { return (R() == b.R()) && (Phi() == b.Phi()); }
		/// @brief Exact inequality comparison
		bool operator!=(const Point2Polar& b) const { return (R() != b.R()) || (Phi() != b.Phi()); }
		/// @brief Approximate equality comparison with tolerance
		/// @param b The other point
		/// @param eps Tolerance for comparison (default: Precision::Pnt2PolarIsEqualTolerance)
		/// @return True if distance between points is less than eps
		bool IsEqualTo(const Point2Polar& b, Real eps = Precision::Pnt2PolarIsEqualTolerance) const { return Dist(b) < eps; }

		/// @brief Converts this polar point to Cartesian coordinates
		/// @return Equivalent Point2Cartesian (x = r*cos(φ), y = r*sin(φ))
		Point2Cartesian TransfToCart() const { return Point2Cartesian(R() * cos(Phi()), R() * sin(Phi())); }

		/// @brief Computes distance using radial differences and a scaled angular chord
		/// @param b The other point
		/// @return Euclidean distance between the two points
		Real Dist(const Point2Polar& b) const {
			return detail::GeometryRadialDistance(R(), Phi(), b.R(), b.Phi(), 0, "Point2Polar::Dist");
		}
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
		Real Dist(const Point3Cartesian& b) const {
			detail::ValidateGeometryCoordinates({X(), Y(), Z(), b.X(), b.Y(), b.Z()}, "Point3Cartesian::Dist");
			return detail::CheckedGeometryDistance(std::hypot(b.X() - X(), b.Y() - Y(), b.Z() - Z()),
				"Point3Cartesian::Dist");
		}

		/// @brief Exact equality comparison
		bool operator==(const Point3Cartesian& b) const { return (X() == b.X()) && (Y() == b.Y()) && (Z() == b.Z()); }
		/// @brief Exact inequality comparison
		bool operator!=(const Point3Cartesian& b) const { return (X() != b.X()) || (Y() != b.Y()) || (Z() != b.Z()); }
		/// @brief Approximate equality comparison with tolerance
		/// @param b The other point
		/// @param eps Tolerance for comparison (default: Precision::Pnt3CartIsEqualTolerance)
		/// @return True if distance between points is less than eps
		bool IsEqualTo(const Point3Cartesian& b, Real eps = Precision::Pnt3CartIsEqualTolerance) const { return Dist(b) < eps; }

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
			_r = Point3Cartesian().Dist(pnt);
			_theta = atan2(std::hypot(pnt.X(), pnt.Y()), pnt.Z());
			_phi = atan2(pnt.Y(), pnt.X());
		}

		/// @brief Exact equality comparison
		bool operator==(const Point3Spherical& b) const { return (R() == b.R()) && (Theta() == b.Theta()) && (Phi() == b.Phi()); }
		/// @brief Exact inequality comparison
		bool operator!=(const Point3Spherical& b) const { return (R() != b.R()) || (Theta() != b.Theta()) || (Phi() != b.Phi()); }
		/// @brief Approximate equality comparison with tolerance
		/// @param b The other point
		/// @param eps Tolerance for comparison (default: Precision::Pnt3SphIsEqualTolerance)
		/// @return True if distance between points is less than eps
		bool IsEqualTo(const Point3Spherical& b, Real eps = Precision::Pnt3SphIsEqualTolerance) const { return Dist(b) < eps; }

		/// @brief Converts this spherical point to Cartesian coordinates
		/// @return Equivalent Point3Cartesian using spherical-to-Cartesian transformation
		Point3Cartesian TransfToCart() const {
			return Point3Cartesian(R() * sin(Theta()) * cos(Phi()), R() * sin(Theta()) * sin(Phi()), R() * cos(Theta()));
		}

		/// @brief Computes distance to another spherical point
		/// @param b The other point
		/// @return Euclidean distance from radial and half-angle chord components
		Real Dist(const Point3Spherical& b) const {
			detail::ValidateGeometryCoordinates({R(), Theta(), Phi(), b.R(), b.Theta(), b.Phi()}, "Point3Spherical::Dist");
			const Real first = std::abs(R()), second = std::abs(b.R());
			const Real sine1 = std::abs(Theta()) == Constants::PI ? 0 : std::sin(Theta());
			const Real sine2 = std::abs(b.Theta()) == Constants::PI ? 0 : std::sin(b.Theta());
			const Real theta1 = sine1 < 0 ? -Theta() : Theta();
			const Real theta2 = sine2 < 0 ? -b.Theta() : b.Theta();
			const bool oppositeRadii = std::signbit(R()) != std::signbit(b.R());
			Real meridianChord;
			if (!oppositeRadii)
				meridianChord = detail::GeometrySineChord(theta1, theta2);
			else if (theta1 >= 0 && theta1 <= Constants::PI && theta2 >= 0 && theta2 <= Constants::PI)
				meridianChord = detail::GeometrySineChord(Constants::PI - std::max(theta1, theta2), std::min(theta1, theta2));
			else
				meridianChord = detail::GeometryCosineChord(theta1, -theta2);
			const bool oppositeMeridians = (sine1 < 0) != (sine2 < 0);
			const Real azimuthChord = oppositeMeridians != oppositeRadii
				? detail::GeometryCosineChord(Phi(), b.Phi()) : detail::GeometrySineChord(Phi(), b.Phi());
			const Real radial = first - second;
			const Real meridian = detail::GeometryAngularComponent(first, second, meridianChord);
			const Real azimuth = detail::GeometryAngularComponent(first * std::abs(sine1), second * std::abs(sine2), azimuthChord);
			return detail::CheckedGeometryDistance(std::hypot(radial, meridian, azimuth), "Point3Spherical::Dist",
				radial != 0 || (first != 0 && second != 0 &&
					(meridianChord != 0 || (sine1 != 0 && sine2 != 0 && azimuthChord != 0))));
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
			detail::ValidateGeometryCoordinates({pnt.X(), pnt.Y(), pnt.Z()}, "Point3Cylindrical conversion");
			_r = Point2Cartesian().Dist(Point2Cartesian(pnt.X(), pnt.Y()));
			_phi = atan2(pnt.Y(), pnt.X());
			_z = pnt.Z();
		}

		/// @brief Exact equality comparison
		bool operator==(const Point3Cylindrical& b) const { return (R() == b.R()) && (Phi() == b.Phi()) && (Z() == b.Z()); }
		/// @brief Exact inequality comparison
		bool operator!=(const Point3Cylindrical& b) const { return (R() != b.R()) || (Phi() != b.Phi()) || (Z() != b.Z()); }
		/// @brief Approximate equality comparison with tolerance
		/// @param b The other point
		/// @param eps Tolerance for comparison (default: Precision::Pnt3CylIsEqualTolerance)
		/// @return True if distance between points is less than eps
		bool IsEqualTo(const Point3Cylindrical& b, Real eps = Precision::Pnt3CylIsEqualTolerance) const { return Dist(b) < eps; }

		/// @brief Converts this cylindrical point to Cartesian coordinates
		/// @return Equivalent Point3Cartesian using cylindrical-to-Cartesian transformation
		Point3Cartesian TransfToCart() const { return Point3Cartesian(R() * cos(Phi()), R() * sin(Phi()), Z()); }

		/// @brief Computes distance to another cylindrical point
		/// @param b The other point
		/// @return Euclidean distance between the two points
		Real Dist(const Point3Cylindrical& b) const {
			detail::ValidateGeometryCoordinates({Z(), b.Z()}, "Point3Cylindrical::Dist");
			return detail::GeometryRadialDistance(R(), Phi(), b.R(), b.Phi(), b.Z() - Z(), "Point3Cylindrical::Dist");
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

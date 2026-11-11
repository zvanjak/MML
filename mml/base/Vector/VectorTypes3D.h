///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        VectorTypes3D.h                                                     ///
///  Description: 3D Cartesian, spherical, and cylindrical vector types              ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_VECTOR_TYPES_3D_H
#define MML_VECTOR_TYPES_3D_H

#include <mml/MMLBase.h>
#include <mml/MMLExceptions.h>

#include <mml/base/Vector/VectorN.h>
#include <mml/base/Geometry/GeometryBase/GeometryPoints.h>
#include <mml/base/BaseUtils.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>

namespace MML
{
	/// @brief 3D Cartesian vector with X, Y, Z components.
	class Vector3Cartesian : public VectorN<Real, 3>
	{
		bool TryUnitDirection(Vector3Cartesian& unit) const
		{
			if (!std::isfinite(X()) || !std::isfinite(Y()) || !std::isfinite(Z()))
				return false;
			const Real scale = std::max({std::abs(X()), std::abs(Y()), std::abs(Z())});
			if (scale == 0)
				return false;
			unit = Vector3Cartesian(X() / scale, Y() / scale, Z() / scale);
			unit /= std::hypot(unit.X(), unit.Y(), unit.Z());
			return true;
		}

	public:
		/// @brief Returns X component (const).
		Real  X() const { return val(0); }
		/// @brief Returns X component (non-const).
		Real& X()				{ return val(0); }
		/// @brief Returns Y component (const).
		Real  Y() const { return val(1); }
		/// @brief Returns Y component (non-const).
		Real& Y()				{ return val(1); }
		/// @brief Returns Z component (const).
		Real  Z() const { return val(2); }
		/// @brief Returns Z component (non-const).
		Real& Z()				{ return val(2); }

		/// @brief Default constructor (zero vector).
		Vector3Cartesian() : VectorN<Real, 3>{ 0.0, 0.0, 0.0 } {}
		/// @brief Constructs from VectorN<Real,3>.
		/// @param b Vector to copy
		Vector3Cartesian(const VectorN<Real, 3>& b) : VectorN<Real, 3>{ b } {}
		/// @brief Constructs from X, Y, Z components.
		/// @param x X component
		/// @param y Y component
		/// @param z Z component
		Vector3Cartesian(Real x, Real y, Real z) : VectorN<Real, 3>{ x, y, z } {}
		/// @brief Constructs from initializer list.
		/// @param list Initializer list of values
		Vector3Cartesian(std::initializer_list<Real> list) : VectorN<Real, 3>(list) { }
		/// @brief Constructs vector from two points (b - a).
		/// @param a Start point
		/// @param b End point
		Vector3Cartesian(const Point3Cartesian& a, const Point3Cartesian& b)
		{
			val(0) = b.X() - a.X();
			val(1) = b.Y() - a.Y();
			val(2) = b.Z() - a.Z();
		}
		/// @brief Constructs from Point3Cartesian (converts point to position vector).
		/// @param a Point to convert
		Vector3Cartesian(const Point3Cartesian& a)
		{
			val(0) = a.X();
			val(1) = a.Y();
			val(2) = a.Z();
		}

		/// @brief Converts vector to Point3Cartesian.
		Point3Cartesian ToPoint() const
		{
			return Point3Cartesian(val(0), val(1), val(2));
		}
		/// @deprecated Use ToPoint() (CONTRACT_POLICY §11: conversions are PascalCase To… nouns).
		[[deprecated("use ToPoint()")]]
		Point3Cartesian getAsPoint() const { return ToPoint(); }
		
		/// @brief Returns normalized vector (unit length).
		/// @throws VectorDimensionError if the vector norm is near zero
		Vector3Cartesian Normalized() const
		{
			Real norm = NormL2();
			if (norm < Precision::RoundoffZeroThreshold)
				throw VectorDimensionError("Vector3Cartesian::Normalized - cannot normalize near-zero vector", 3, 0);
			return Vector3Cartesian{ (*this) / norm };
		}
		/// @brief Returns normalized vector, or the zero vector when norm is near zero.
		Vector3Cartesian NormalizedOrZero() const
		{
			Real norm = NormL2();
			if (norm < Precision::RoundoffZeroThreshold)
				return Vector3Cartesian{ static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0) };
			return Vector3Cartesian{ (*this) / norm };
		}
		/// @deprecated Use Normalized() — identical semantics (CONTRACT_POLICY §11: one spelling per operation).
		[[deprecated("use Normalized()")]]
		Vector3Cartesian GetAsUnitVector() const { return Normalized(); }
		/// @deprecated Use NormalizedOrZero().
		[[deprecated("use NormalizedOrZero()")]]
		Vector3Cartesian GetAsUnitVectorOrZero() const { return NormalizedOrZero(); }
		/// @brief Returns unit vector at given position (for Cartesian, position doesn't matter).
		/// @param pos Position (unused)
		Vector3Cartesian GetAsUnitVectorAtPos(const Vector3Cartesian& pos) const
		{
			return Normalized();
		}

		/// @brief Scalar product (dot product) via operator*.
		/// @param b Vector to multiply
		Real operator*(const Vector3Cartesian& b) const
		{
			return X()*b.X() + Y()*b.Y() + Z()*b.Z();
		}

		/// @brief Unary minus (negation).
		Vector3Cartesian operator-() const
		{
			Vector3Cartesian ret;
			for (int i = 0; i < 3; i++)
				ret[i] = -val(i);
			return ret;
		}

		/// @brief Vector addition.
		/// @param b Vector to add
		Vector3Cartesian operator+(const Vector3Cartesian& b) const
		{
			Vector3Cartesian ret;
			for (int i = 0; i < 3; i++)
				ret[i] = val(i) + b[i];
			return ret;
		}
		/// @brief Vector subtraction.
		/// @param b Vector to subtract
		Vector3Cartesian operator-(const Vector3Cartesian& b) const
		{
			Vector3Cartesian ret;
			for (int i = 0; i < 3; i++)
				ret[i] = val(i) - b[i];
			return ret;
		}

		/// @brief Scalar multiplication (vector * scalar).
		/// @param b Scalar multiplier
		Vector3Cartesian operator*(Real b) const
		{
			Vector3Cartesian ret;
			for (int i = 0; i < 3; i++)
				ret[i] = val(i) * b;
			return ret;
		}
		/// @brief Scalar division (vector / scalar).
		/// @param b Scalar divisor
		Vector3Cartesian operator/(Real b) const
		{
			Vector3Cartesian ret;
			for (int i = 0; i < 3; i++)
				ret[i] = val(i) / b;
			return ret;
		}
		/// @brief Scalar multiplication (scalar * vector).
		/// @param a Scalar multiplier
		/// @param b Vector
		friend Vector3Cartesian operator*(Real a, const Vector3Cartesian& b)
		{
			Vector3Cartesian ret;
			for (int i = 0; i < 3; i++)
				ret[i] = a * b[i];
			return ret;
		}

		/// @brief Exact equality comparison.
		/// @param b Vector to compare
		bool operator==(const Vector3Cartesian& b) const
		{
			return (X() == b.X()) && (Y() == b.Y()) && (Z() == b.Z());
		}
		/// @brief Inequality comparison.
		/// @param b Vector to compare
		bool operator!=(const Vector3Cartesian& b) const
		{
			return (X() != b.X()) || (Y() != b.Y()) || (Z() != b.Z());
		}
		/// @brief Checks equality within tolerance.
		/// @param b Vector to compare
		/// @param absEps Absolute tolerance
		bool IsEqualTo(const Vector3Cartesian& b, Real absEps = Precision::Vec3CartIsEqualTolerance) const
		{
			if (!std::isfinite(absEps) || absEps < 0)
				return false;
			return (std::abs(X() - b.X()) < absEps) && (std::abs(Y() - b.Y()) < absEps) && (std::abs(Z() - b.Z()) < absEps);
		}

		/// @brief Point + vector operation.
		friend Point3Cartesian operator+(const Point3Cartesian& a, const Vector3Cartesian& b) { return Point3Cartesian(a.X() + b[0], a.Y() + b[1], a.Z() + b[2]); }
		/// @brief Point - vector operation.
		friend Point3Cartesian operator-(const Point3Cartesian& a, const Vector3Cartesian& b) { return Point3Cartesian(a.X() - b[0], a.Y() - b[1], a.Z() - b[2]); }

		/// @brief Checks if vectors are parallel (same OR opposite direction).
		/// @details Criterion: ‖a × b‖ ≤ eps · ‖a‖ · ‖b‖, i.e. |sin θ| ≤ eps, scale-invariant and
		/// antiparallel-inclusive (parallel lines/planes have parallel direction vectors regardless
		/// of orientation). The zero vector has no direction and is parallel to nothing (returns false).
		/// @param b Vector to compare
		/// @param eps Tolerance on |sin θ|
		bool IsParallelTo(const Vector3Cartesian& b, Real eps = Precision::Vec3CartIsParallelTolerance) const
		{
			Vector3Cartesian first, second;
			if (!std::isfinite(eps) || eps < 0 || !TryUnitDirection(first) || !b.TryUnitDirection(second))
				return false;
			return VectorProduct(first, second).NormL2() <= eps;
		}
		/// @brief Checks if vectors are perpendicular.
		/// @param b Vector to compare
		/// @param eps Dimensionless tolerance on the absolute cosine of the angle
		bool IsPerpendicularTo(const Vector3Cartesian& b, Real eps = Precision::OrthogonalityTolerance) const
		{
			Vector3Cartesian first, second;
			if (!std::isfinite(eps) || eps < 0 || !TryUnitDirection(first) || !b.TryUnitDirection(second))
				return false;
			return std::abs(first.ScalarProduct(second)) <= eps;
		}
		/// @brief Computes two vectors perpendicular to this vector.
		/// @param v1 First perpendicular vector (output)
		/// @param v2 Second perpendicular vector (output)
		/// @return true if successful, false if zero vector
		bool GetPerpendicularVectors(Vector3Cartesian& v1, Vector3Cartesian& v2) const
		{
			Vector3Cartesian direction;
			if (!TryUnitDirection(direction))
				return false;
			
			// find a vector that is not parallel to this vector
			Vector3Cartesian not_parallel;
			if (std::abs(X()) <= std::abs(Y()) && std::abs(X()) <= std::abs(Z()))
				not_parallel = Vector3Cartesian(1.0, 0.0, 0.0);
			else if (std::abs(Y()) <= std::abs(X()) && std::abs(Y()) <= std::abs(Z()))
				not_parallel = Vector3Cartesian(0.0, 1.0, 0.0);
			else
				not_parallel = Vector3Cartesian(0.0, 0.0, 1.0);
			
			v1 = VectorProduct(direction, not_parallel).Normalized();
			v2 = VectorProduct(direction, v1).Normalized();
			
			return true;
		}

		/// @brief Computes angle between this vector and another.
		/// @param b Other vector
		/// @return Angle in radians [0, p]
		Real AngleToVector(const Vector3Cartesian& b)
		{
			if (!std::isfinite(X()) || !std::isfinite(Y()) || !std::isfinite(Z()) ||
				!std::isfinite(b.X()) || !std::isfinite(b.Y()) || !std::isfinite(b.Z()))
				throw NumericInputError("AngleToVector requires finite vectors");
			Vector3Cartesian first, second;
			if (!TryUnitDirection(first) || !b.TryUnitDirection(second))
				throw DivisionByZeroError("AngleToVector: cannot compute angle with zero vector");
			return std::atan2(VectorProduct(first, second).NormL2(), first.ScalarProduct(second));
		}

		/// @brief Scalar product (member method).
		/// @param b Vector to multiply
		Real ScalarProduct(const Vector3Cartesian& b) const
		{
			return (*this) * b;
		}

		/// @brief Scalar product (friend function).
		/// @param a First vector
		/// @param b Second vector
		friend Real ScalarProduct(const Vector3Cartesian& a, const Vector3Cartesian& b)
		{
			return a * b;
		}

		/// @brief Vector (cross) product.
		/// @param a First vector
		/// @param b Second vector
		friend Vector3Cartesian VectorProduct(const Vector3Cartesian& a, const Vector3Cartesian& b)
		{
			Vector3Cartesian ret;

			ret.X() = a.Y() * b.Z() - a.Z() * b.Y();
			ret.Y() = a.Z() * b.X() - a.X() * b.Z();
			ret.Z() = a.X() * b.Y() - a.Y() * b.X();

			return ret;
		}
	};

	/// @brief 3D Spherical vector with R (radius), Theta (polar angle), Phi (azimuthal angle) components.
	class Vector3Spherical : public VectorN<Real, 3>
	{
	public:
		/// @brief Returns R component (const).
		Real  R()     const { return val(0); }
		/// @brief Returns R component (non-const).
		Real& R()						{ return val(0); }
		/// @brief Returns Theta component (const).
		Real  Theta() const { return val(1); }
		/// @brief Returns Theta component (non-const).
		Real& Theta()				{ return val(1); }
		/// @brief Returns Phi component (const).
		Real  Phi()   const { return val(2); }
		/// @brief Returns Phi component (non-const).
		Real& Phi()					{ return val(2); }

		/// @brief Default constructor (zero vector).
		Vector3Spherical() : VectorN<Real, 3>{ 0.0, 0.0, 0.0 } {}
		/// @brief Constructs from VectorN<Real,3>.
		/// @param b Vector to copy
		Vector3Spherical(const VectorN<Real, 3>& b) : VectorN<Real, 3>{ b[0], b[1], b[2] } {}
		/// @brief Constructs from R, Theta, Phi components.
		/// @param r Radius
		/// @param theta Polar angle
		/// @param phi Azimuthal angle
		Vector3Spherical(Real r, Real theta, Real phi) : VectorN<Real, 3>{ r, theta, phi } {}
		/// @brief Constructs from initializer list.
		/// @param list Initializer list of values
		Vector3Spherical(std::initializer_list<Real> list) : VectorN<Real, 3>(list) { }

		/// @brief Unary minus (negation, flips direction in spherical coordinates).
		Vector3Spherical operator-() const
		{
			return Vector3Spherical(R(), Constants::PI - Theta(), Phi() + Constants::PI);
		}

		/// @brief Vector addition (converts to/from Cartesian).
		/// @param b Vector to add
		Vector3Spherical operator+(const Vector3Spherical& b) const
		{
			// Convert both to Cartesian
			Real x1, y1, z1, x2, y2, z2;
			Utils::SphericalToCartesian(R(), Theta(), Phi(), x1, y1, z1);
			Utils::SphericalToCartesian(b.R(), b.Theta(), b.Phi(), x2, y2, z2);
			
			// Add and convert back
			Real r, theta, phi;
			Utils::CartesianToSpherical(x1 + x2, y1 + y2, z1 + z2, r, theta, phi);
			return Vector3Spherical(r, theta, phi);
		}
		/// @brief Vector subtraction (converts to/from Cartesian).
		/// @param b Vector to subtract
		Vector3Spherical operator-(const Vector3Spherical& b) const
		{
			// Convert both to Cartesian
			Real x1, y1, z1, x2, y2, z2;
			Utils::SphericalToCartesian(R(), Theta(), Phi(), x1, y1, z1);
			Utils::SphericalToCartesian(b.R(), b.Theta(), b.Phi(), x2, y2, z2);
			
			// Subtract and convert back
			Real r, theta, phi;
			Utils::CartesianToSpherical(x1 - x2, y1 - y2, z1 - z2, r, theta, phi);
			return Vector3Spherical(r, theta, phi);
		}
		Vector3Spherical& operator+=(const Vector3Spherical& b) { *this = *this + b; return *this; }
		Vector3Spherical& operator-=(const Vector3Spherical& b) { *this = *this - b; return *this; }
		Vector3Spherical& operator*=(Real b) { *this = *this * b; return *this; }
		Vector3Spherical& operator/=(Real b) { *this = *this / b; return *this; }

		/// @brief Scalar multiplication (scales radius).
		/// @param b Scalar multiplier
		Vector3Spherical operator*(Real b) const
		{
			return Vector3Spherical(R() * b, Theta(), Phi());
		}
		/// @brief Scalar division (divides radius).
		/// @param b Scalar divisor
		Vector3Spherical operator/(Real b) const
		{
			if (b == 0.0)
				throw DivisionByZeroError("Division by zero in Vector3Spherical division.");
			return Vector3Spherical(R() / b, Theta(), Phi());
		}
		/// @brief Scalar multiplication (scalar * vector).
		/// @param a Scalar multiplier
		/// @param b Vector
		friend Vector3Spherical operator*(Real a, const Vector3Spherical& b)
		{
			return Vector3Spherical(b.R() * a, b.Theta(), b.Phi());
		}

		/// @brief Exact equality comparison.
		/// @param b Vector to compare
		bool operator==(const Vector3Spherical& b) const
		{
			return (R() == b.R()) && (Theta() == b.Theta()) && (Phi() == b.Phi());
		}
		/// @brief Inequality comparison.
		/// @param b Vector to compare
		bool operator!=(const Vector3Spherical& b) const
		{
			return (R() != b.R()) || (Theta() != b.Theta()) || (Phi() != b.Phi());
		}
		/// @brief Checks equality within tolerance with angle wrap-around awareness.
		/// @param b Vector to compare
		/// @param absEps Absolute tolerance
		/// @note For spherical vectors, theta ∈ [0, π] and phi ∈ [-π, π].
		///       Phi comparison uses AnglesAreEqual() to handle wrap-around at ±π.
		bool IsEqualTo(const Vector3Spherical& b, Real absEps = Precision::Vec3SphIsEqualTolerance) const
		{
			if (!std::isfinite(absEps) || absEps < 0)
				return false;
			// R and theta are straightforward comparisons
			// Phi needs angle-aware comparison for wrap-around at ±π
			return (std::abs(R() - b.R()) < absEps) && 
			       (std::abs(Theta() - b.Theta()) < absEps) && 
			       AnglesAreEqual(Phi(), b.Phi(), absEps);
		}
		
		/// @brief Returns unit vector at given position (global spherical result).
		/// @param pos Position in spherical coordinates
		Vector3Spherical GetAsUnitVectorAtPosGlobalSpher(const Vector3Spherical& pos) const
		{
			if (!std::isfinite(R()) || !std::isfinite(Theta()) || !std::isfinite(Phi()) ||
				!std::isfinite(pos.Theta()) || !std::isfinite(pos.Phi()))
				throw NumericInputError("Spherical unit directions require finite components and position angles");
			if (R() == 0)
				return Vector3Spherical(0, 0, 0);
			// Convert this vector to Cartesian (as a displacement from the origin)
			Real r = std::copysign(Real{1}, R());
			Real theta = Theta();
			Real phi = Phi();

			// Spherical to Cartesian (displacement from origin)
			Real x = r * std::sin(theta) * std::cos(phi);
			Real y = r * std::sin(theta) * std::sin(phi);
			Real z = r * std::cos(theta);

			// Now, express this vector in the local spherical basis at 'pos'
			// The local basis at 'pos' consists of:
			//   e_r     = [sinθ cosφ, sinθ sinφ, cosθ]
			//   e_theta = [cosθ cosφ, cosθ sinφ, -sinθ]
			//   e_phi   = [-sinφ, cosφ, 0]
			// at the position (pos.Theta(), pos.Phi())

			Real st = std::sin(pos.Theta());
			Real ct = std::cos(pos.Theta());
			Real sp = std::sin(pos.Phi());
			Real cp = std::cos(pos.Phi());

			// Local basis vectors at pos (expressed in Cartesian coordinates)
			Real e_r[3] = { st * cp, st * sp, ct };
			Real e_theta[3] = { ct * cp, ct * sp, -st };
			Real e_phi[3] = { -sp,   cp,    0 };

			// Project (x, y, z) onto the local basis at pos
			Real v_r = x * e_r[0] + y * e_r[1] + z * e_r[2];
			Real v_theta = x * e_theta[0] + y * e_theta[1] + z * e_theta[2];
			Real v_phi = x * e_phi[0] + y * e_phi[1] + z * e_phi[2];

			// Normalize
			Real norm = std::hypot(v_r, v_theta, v_phi);
			if (norm == 0.0)
				return Vector3Spherical(0.0, 0.0, 0.0);

			// Clamp to [-1, 1] to prevent NaN from floating-point rounding
			Real cos_theta = std::max(Real{-1.0}, std::min(Real{1.0}, v_r / norm));
			return Vector3Spherical(1.0, std::acos(cos_theta), std::atan2(v_phi, v_theta));
		}
		/// @brief Returns unit vector at given position (local spherical basis result).
		/// @param pos Position in spherical coordinates
		Vector3Spherical GetAsUnitVectorAtPos(const Vector3Spherical& pos) const
		{
			if (!std::isfinite(R()) || !std::isfinite(Theta()) || !std::isfinite(Phi()) ||
				!std::isfinite(pos.Theta()) || !std::isfinite(pos.Phi()))
				throw NumericInputError("Spherical unit directions require finite components and position angles");
			if (R() == 0)
				return Vector3Spherical(0, 0, 0);
			// Step 1: Convert *this to Cartesian (as a displacement from the origin)
			Real r = std::copysign(Real{1}, R());
			Real theta = Theta();
			Real phi = Phi();

			Real x = r * std::sin(theta) * std::cos(phi);
			Real y = r * std::sin(theta) * std::sin(phi);
			Real z = r * std::cos(theta);

			// Step 2: Compute local spherical basis at pos
			Real st = std::sin(pos.Theta());
			Real ct = std::cos(pos.Theta());
			Real sp = std::sin(pos.Phi());
			Real cp = std::cos(pos.Phi());

			// Local basis vectors at pos
			Real e_r[3] = { st * cp, st * sp, ct };
			Real e_theta[3] = { ct * cp, ct * sp, -st };
			Real e_phi[3] = { -sp,   cp,    0 };

			// Step 3: Project (x, y, z) onto the local basis at pos
			Real v_r = x * e_r[0] + y * e_r[1] + z * e_r[2];
			Real v_theta = x * e_theta[0] + y * e_theta[1] + z * e_theta[2];
			Real v_phi = x * e_phi[0] + y * e_phi[1] + z * e_phi[2];

			// Step 4: Normalize
			Real norm = std::hypot(v_r, v_theta, v_phi);
			if (norm == 0.0)
				return Vector3Spherical(0.0, 0.0, 0.0);

			return Vector3Spherical(v_r / norm, v_theta / norm, v_phi / norm);
		}

		/// @brief Prints vector with angles in degrees.
		/// @param stream Output stream
		/// @param width Field width
		/// @param precision Decimal precision
		std::ostream& PrintDeg(std::ostream& stream, int width, int precision) const
		{
			stream << "[ ";
			stream << std::fixed << std::setw(width) << std::setprecision(precision);
			stream << R();
			stream << ", " << Theta() * 180.0 / Constants::PI;
			stream << ", " << Phi() * 180.0 / Constants::PI << " ]" << std::endl;

			return stream;
		}
	};

	/// @brief 3D Cylindrical vector with R (radial), Phi (azimuthal), Z components.
	class Vector3Cylindrical : public VectorN<Real, 3>
	{
	public:
		/// @brief Returns R component (const).
		Real  R()   const { return val(0); }
		/// @brief Returns R component (non-const).
		Real& R()					{ return val(0); }
		/// @brief Returns Phi component (const).
		Real  Phi() const { return val(1); }
		/// @brief Returns Phi component (non-const).
		Real& Phi()				{ return val(1); }
		/// @brief Returns Z component (const).
		Real  Z()   const { return val(2); }
		/// @brief Returns Z component (non-const).
		Real& Z()					{ return val(2); }

		/// @brief Default constructor (zero vector).
		Vector3Cylindrical() : VectorN<Real, 3>{ 0.0, 0.0, 0.0 } {}
		/// @brief Constructs from VectorN<Real,3>.
		/// @param b Vector to copy
		Vector3Cylindrical(const VectorN<Real, 3>& b) : VectorN<Real, 3>{ b[0], b[1], b[2] } {}
		/// @brief Constructs from R, Phi, Z components.
		/// @param r Radial distance
		/// @param phi Azimuthal angle
		/// @param z Vertical coordinate
		Vector3Cylindrical(Real r, Real phi, Real z) : VectorN<Real, 3>{ r, phi, z } {}
		/// @brief Constructs from initializer list.
		/// @param list Initializer list of values
		Vector3Cylindrical(std::initializer_list<Real> list) : VectorN<Real, 3>(list) { }

		/// @brief Unary minus (negation, adjusts angle and flips Z).
		Vector3Cylindrical operator-() const
		{
			return Vector3Cylindrical(R(), Phi() + Constants::PI, -Z());
		}

		/// @brief Vector addition (converts to/from Cartesian).
		/// @param b Vector to add
		Vector3Cylindrical operator+(const Vector3Cylindrical& b) const
		{
			// Convert both to Cartesian
			Real x1, y1, z1, x2, y2, z2;
			Utils::CylindricalToCartesian(R(), Phi(), Z(), x1, y1, z1);
			Utils::CylindricalToCartesian(b.R(), b.Phi(), b.Z(), x2, y2, z2);
			
			// Add and convert back
			Real r, phi, z;
			Utils::CartesianToCylindrical(x1 + x2, y1 + y2, z1 + z2, r, phi, z);
			return Vector3Cylindrical(r, phi, z);
		}
		/// @brief Vector subtraction (converts to/from Cartesian).
		/// @param b Vector to subtract
		Vector3Cylindrical operator-(const Vector3Cylindrical& b) const
		{
			// Convert both to Cartesian
			Real x1, y1, z1, x2, y2, z2;
			Utils::CylindricalToCartesian(R(), Phi(), Z(), x1, y1, z1);
			Utils::CylindricalToCartesian(b.R(), b.Phi(), b.Z(), x2, y2, z2);
			
			// Subtract and convert back
			Real r, phi, z;
			Utils::CartesianToCylindrical(x1 - x2, y1 - y2, z1 - z2, r, phi, z);
			return Vector3Cylindrical(r, phi, z);
		}
		Vector3Cylindrical& operator+=(const Vector3Cylindrical& b) { *this = *this + b; return *this; }
		Vector3Cylindrical& operator-=(const Vector3Cylindrical& b) { *this = *this - b; return *this; }
		Vector3Cylindrical& operator*=(Real b) { *this = *this * b; return *this; }
		Vector3Cylindrical& operator/=(Real b) { *this = *this / b; return *this; }

		/// @brief Scalar multiplication (scales R and Z).
		/// @param b Scalar multiplier
		Vector3Cylindrical operator*(Real b) const
		{
			return Vector3Cylindrical(R() * b, Phi(), Z() * b);
		}
		/// @brief Scalar division (divides R and Z).
		/// @param b Scalar divisor
		Vector3Cylindrical operator/(Real b) const
		{
			if (b == 0.0)
				throw DivisionByZeroError("Division by zero in Vector3Cylindrical division.");
			return Vector3Cylindrical(R() / b, Phi(), Z() / b);
		}
		/// @brief Scalar multiplication (scalar * vector).
		/// @param a Scalar multiplier
		/// @param b Vector
		friend Vector3Cylindrical operator*(Real a, const Vector3Cylindrical& b)
		{
			return Vector3Cylindrical(b.R() * a, b.Phi(), b.Z() * a);
		}

		/// @brief Returns unit vector at given position (local cylindrical basis result).
		/// @param pos Position in cylindrical coordinates
		Vector3Cylindrical GetAsUnitVectorAtPos(const Vector3Cylindrical& pos) const
		{
			if (!std::isfinite(R()) || !std::isfinite(Phi()) || !std::isfinite(Z()) || !std::isfinite(pos.Phi()))
				throw NumericInputError("Cylindrical unit directions require finite components and position angles");
			const Real scale = std::max(std::abs(R()), std::abs(Z()));
			if (scale == 0)
				return Vector3Cylindrical(0, 0, 0);
			// Step 1: Convert *this to Cartesian
			Real x = (R() / scale) * std::cos(Phi());
			Real y = (R() / scale) * std::sin(Phi());
			Real z = Z() / scale;

			// Step 2: Normalize in Cartesian
			Real norm = std::hypot(x, y, z);
			if (norm == static_cast<MML::Real>(0.0))
				return Vector3Cylindrical(static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0));

			x /= norm;
			y /= norm;
			z /= norm;

			// Step 3: Compute local cylindrical basis at pos and project
			Real cp = std::cos(pos.Phi());
			Real sp = std::sin(pos.Phi());

			// ê_R   = (cos(φ_pos), sin(φ_pos), 0)
			// ê_Phi = (-sin(φ_pos), cos(φ_pos), 0)
			// ê_Z   = (0, 0, 1)
			Real v_R = x * cp + y * sp;
			Real v_Phi = -x * sp + y * cp;
			Real v_Z = z;

			return Vector3Cylindrical(v_R, v_Phi, v_Z);
		}
	};

	typedef Vector3Cartesian   Vec3Cart;
	typedef Vector3Spherical   Vec3Sph;
	typedef Vector3Cylindrical Vec3Cyl;
}

#endif

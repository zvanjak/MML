///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        VectorTypes.h                                                       ///
///  Description: Type aliases for 2D, 3D, 4D vectors (Vec2, Vec3, Vec4)              ///
///               Complex vector types and related utilities                          ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_VECTOR_TYPES_H
#define MML_VECTOR_TYPES_H

#include "MMLBase.h"
#include "MMLExceptions.h"

#include "base/Vector/VectorN.h"
#include "base/Geometry/Geometry.h"
#include "base/BaseUtils.h"

// Standard headers - include what we use
#include <algorithm>
#include <iomanip>
#include <iostream>

namespace MML
{
	/// @brief 2D Cartesian vector with X, Y components.
	class Vector2Cartesian : public VectorN<Real, 2>
	{
	public:
		/// @brief Default constructor (zero vector).
		Vector2Cartesian() {}
		/// @brief Constructs from X, Y components.
		/// @param x X component
		/// @param y Y component
		Vector2Cartesian(Real x, Real y)
		{
			_val[0] = x;
			_val[1] = y;
		}
		/// @brief Constructs from VectorN<Real,2>.
		/// @param b Vector to copy
		Vector2Cartesian(const VectorN<Real, 2>& b) : VectorN<Real, 2>{ b[0], b[1] } {}
		/// @brief Constructs vector from two points (b - a).
		/// @param a Start point
		/// @param b End point
		Vector2Cartesian(const Point2Cartesian& a, const Point2Cartesian& b)
		{
			_val[0] = b.X() - a.X();
			_val[1] = b.Y() - a.Y();
		}
		/// @brief Constructs from initializer list.
		/// @param list Initializer list of values
		Vector2Cartesian(std::initializer_list<Real> list) : VectorN<Real, 2>(list) {}

		/// @brief Returns X component (const).
		Real  X() const { return _val[0]; }
		/// @brief Returns X component (non-const).
		Real& X()				{ return _val[0]; }
		/// @brief Returns Y component (const).
		Real  Y() const { return _val[1]; }
		/// @brief Returns Y component (non-const).
		Real& Y()				{ return _val[1]; }
		
		/// @brief Unary minus (negation).
		Vector2Cartesian operator-() const
		{
			Vector2Cartesian ret;
			for (int i = 0; i < 2; i++)
				ret._val[i] = -_val[i];
			return ret;
		}
		
		/// @brief Vector addition.
		/// @param b Vector to add
		Vector2Cartesian operator+(const Vector2Cartesian& b) const
		{
			Vector2Cartesian ret;
			for (int i = 0; i < 2; i++)
				ret._val[i] = _val[i] + b[i];
			return ret;
		}
		/// @brief Vector subtraction.
		/// @param b Vector to subtract
		Vector2Cartesian operator-(const Vector2Cartesian& b) const
		{
			Vector2Cartesian ret;
			for (int i = 0; i < 2; i++)
				ret._val[i] = _val[i] - b[i];
			return ret;
		}

		/// @brief Scalar multiplication (vector * scalar).
		/// @param b Scalar multiplier
		Vector2Cartesian operator*(Real b) const
		{
			Vector2Cartesian ret;
			for (int i = 0; i < 2; i++)
				ret._val[i] = _val[i] * b;
			return ret;
		}
		/// @brief Scalar division (vector / scalar).
		/// @param b Scalar divisor
		Vector2Cartesian operator/(Real b) const
		{
			Vector2Cartesian ret;
			for (int i = 0; i < 2; i++)
				ret._val[i] = _val[i] / b;
			return ret;
		}
		/// @brief Scalar multiplication (scalar * vector).
		/// @param a Scalar multiplier
		/// @param b Vector
		friend Vector2Cartesian operator*(Real a, const Vector2Cartesian& b)
		{
			Vector2Cartesian ret;
			for (int i = 0; i < 2; i++)
				ret._val[i] = a * b[i];
			return ret;
		}

		/// @brief Exact equality comparison.
		/// @param b Vector to compare
		bool operator==(const Vector2Cartesian& b) const
		{
			return (X() == b.X()) && (Y() == b.Y());
		}
		/// @brief Inequality comparison.
		/// @param b Vector to compare
		bool operator!=(const Vector2Cartesian& b) const
		{
			return (X() != b.X()) || (Y() != b.Y());
		}
		/// @brief Checks equality within tolerance.
		/// @param b Vector to compare
		/// @param absEps Absolute tolerance
		bool IsEqualTo(const Vector2Cartesian& b, Real absEps = Defaults::Vec2CartIsEqualTolerance) const
		{
			return (std::abs(X() - b.X()) < absEps) && (std::abs(Y() - b.Y()) < absEps);
		}

		/// @brief Returns normalized vector (unit length).
		Vector2Cartesian Normalized() const
		{
			Vector2Cartesian ret;
			Real norm = NormL2();
			if (norm > 0.0)
			{
				for (int i = 0; i < 2; i++)
					ret._val[i] = _val[i] / norm;
			}
			else
			{
				ret._val[0] = 0.0;
				ret._val[1] = 0.0;
			}
			return ret;
		}
		/// @brief Returns unit vector in this direction.
		Vector2Cartesian GetAsUnitVector() const
		{
			VectorN<Real, 2> res = (*this) / NormL2();

			return Vector2Cartesian(res[0], res[1]);
		}
		/// @brief Returns unit vector at given position (for Cartesian, position doesn't matter).
		/// @param pos Position (unused in Cartesian coordinates)
		Vector2Cartesian GetAsUnitVectorAtPos(const Vector2Cartesian& pos) const
		{
			return Vector2Cartesian{ (*this) / NormL2() };
		}
		
		/// @brief Computes two perpendicular vectors.
		/// @param v1 First perpendicular vector (output)
		/// @param v2 Second perpendicular vector (output)
		void getPerpendicularVectors(Vector2Cartesian& v1, Vector2Cartesian& v2) const
		{
			v1 = Vector2Cartesian(-Y(), X());
			v2 = Vector2Cartesian(Y(), -X());
		}

		/// @brief Scalar product (dot product).
		/// @param b Vector to multiply
		Real operator*(const Vector2Cartesian& b) const
		{
			return X() * b.X() + Y() * b.Y();
		}

		/// @brief Scalar product (friend function).
		/// @param a First vector
		/// @param b Second vector
		friend Real ScalarProduct(const Vector2Cartesian& a, const Vector2Cartesian& b)
		{
			return a * b;
		}

		/// @brief Point + vector operation.
		friend Point2Cartesian operator+(const Point2Cartesian& a, const Vector2Cartesian& b) { return Point2Cartesian(a.X() + b[0], a.Y() + b[1]); }
		/// @brief Point - vector operation.
		friend Point2Cartesian operator-(const Point2Cartesian& a, const Vector2Cartesian& b) { return Point2Cartesian(a.X() - b[0], a.Y() - b[1]); }
	};

	/// @brief 2D Polar vector with R (radius), Phi (angle) components.
	class Vector2Polar : public VectorN<Real, 2>
	{
	public:
		/// @brief Returns R component (const).
		Real  R() const		{ return _val[0]; }
		/// @brief Returns R component (non-const).
		Real& R()					{ return _val[0]; }
		/// @brief Returns Phi component (const).
		Real  Phi() const { return _val[1]; }
		/// @brief Returns Phi component (non-const).
		Real& Phi()				{ return _val[1]; }

		/// @brief Default constructor (zero vector).
		Vector2Polar() {}
		/// @brief Constructs from R, Phi components.
		/// @param r Radius
		/// @param phi Angle
		Vector2Polar(Real r, Real phi)
		{
			_val[0] = r;
			_val[1] = phi;
		}
		/// @brief Constructs from VectorN<Real,2>.
		/// @param b Vector to copy
		Vector2Polar(const VectorN<Real, 2>& b) : VectorN<Real, 2>{ b[0], b[1] } {}

		/// @brief Unary minus (negation, adds p to angle).
		Vector2Polar operator-() const
		{
			// Negating a polar vector means keeping the radius and adding pi to the angle
			return Vector2Polar(R(), Phi() + Constants::PI);
		}
		
		/// @brief Vector addition (converts to Cartesian, adds, converts back).
		/// @param b Vector to add
		Vector2Polar operator+(const Vector2Polar& b) const
		{
			Real r1 = R();
			Real phi1 = Phi();
			Real r2 = b.R();
			Real phi2 = b.Phi();

			Real delta = phi2 - phi1;
			Real r	 = std::sqrt(r1 * r1 + r2 * r2 + 2 * r1 * r2 * std::cos(delta));
			Real phi = phi1 + std::atan2(r2 * std::sin(delta), r1 + r2 * std::cos(delta));

			return Vector2Polar(r, phi);
		}
		/// @brief Vector subtraction (converts to Cartesian, subtracts, converts back).
		/// @param b Vector to subtract
		Vector2Polar operator-(const Vector2Polar& b) const
		{
			Real r1 = R();
			Real phi1 = Phi();
			Real r2 = b.R();
			Real phi2 = b.Phi();

			Real delta = phi2 - phi1;
			Real r	 = std::sqrt(r1 * r1 + r2 * r2 - 2 * r1 * r2 * std::cos(delta));
			Real phi = phi1 + std::atan2(-r2 * std::sin(delta), r1 - r2 * std::cos(delta));

			return Vector2Polar(r, phi);
		}
		/// @brief Scalar multiplication (scales radius).
		/// @param b Scalar multiplier
		Vector2Polar operator*(Real b) const
		{
			return Vector2Polar(R() * b, Phi());
		}
		/// @brief Scalar division (divides radius).
		/// @param b Scalar divisor
		Vector2Polar operator/(Real b) const
		{
			return Vector2Polar(R() / b, Phi());
		}

		/// @brief Scalar multiplication (scalar * vector).
		/// @param a Scalar multiplier
		/// @param b Vector
		friend Vector2Polar operator*(Real a, const Vector2Polar& b)
		{
			return Vector2Polar(b.R() * a, b.Phi());
		}

		/// @brief Returns unit vector at given position.
		/// @param pos Position (unused)
		Vector2Polar GetAsUnitVectorAtPos(const Vector2Polar& /*pos*/) const
		{
			// Returns a unit vector in the direction of this vector (ignores pos)
			return Vector2Polar(1.0, Phi());
		}
	};

	/// @brief 3D Cartesian vector with X, Y, Z components.
	class Vector3Cartesian : public VectorN<Real, 3>
	{
	public:
		/// @brief Returns X component (const).
		Real  X() const { return _val[0]; }
		/// @brief Returns X component (non-const).
		Real& X()				{ return _val[0]; }
		/// @brief Returns Y component (const).
		Real  Y() const { return _val[1]; }
		/// @brief Returns Y component (non-const).
		Real& Y()				{ return _val[1]; }
		/// @brief Returns Z component (const).
		Real  Z() const { return _val[2]; }
		/// @brief Returns Z component (non-const).
		Real& Z()				{ return _val[2]; }

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
			_val[0] = b.X() - a.X();
			_val[1] = b.Y() - a.Y();
			_val[2] = b.Z() - a.Z();
		}
		/// @brief Constructs from Point3Cartesian (converts point to position vector).
		/// @param a Point to convert
		Vector3Cartesian(const Point3Cartesian& a)
		{
			_val[0] = a.X();
			_val[1] = a.Y();
			_val[2] = a.Z();
		}

		/// @brief Converts vector to Point3Cartesian.
		Point3Cartesian getAsPoint()
		{
			return Point3Cartesian(_val[0], _val[1], _val[2]);
		}
		
		/// @brief Returns normalized vector (unit length).
		Vector3Cartesian Normalized() const
		{
			Vector3Cartesian ret;
			Real norm = NormL2();
			if (norm > 0.0)
			{
				for (int i = 0; i < 3; i++)
					ret._val[i] = _val[i] / norm;
			}
			else
			{
				ret._val[0] = 0.0;
				ret._val[1] = 0.0;
				ret._val[2] = 0.0;
			}
			return ret;
		}
		Vector3Cartesian GetAsUnitVector() const
		{
			if(NormL2() == 0.0)
			{
				return Vector3Cartesian{ REAL(0.0), REAL(0.0), REAL(0.0) };
			}
			return Vector3Cartesian{ (*this) / NormL2() };
		}
		/// @brief Returns unit vector at given position (for Cartesian, position doesn't matter).
		/// @param pos Position (unused)
		Vector3Cartesian GetAsUnitVectorAtPos(const Vector3Cartesian& pos) const
		{
			return GetAsUnitVector();
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
				ret._val[i] = -_val[i];
			return ret;
		}

		/// @brief Vector addition.
		/// @param b Vector to add
		Vector3Cartesian operator+(const Vector3Cartesian& b) const
		{
			Vector3Cartesian ret;
			for (int i = 0; i < 3; i++)
				ret._val[i] = _val[i] + b[i];
			return ret;
		}
		/// @brief Vector subtraction.
		/// @param b Vector to subtract
		Vector3Cartesian operator-(const Vector3Cartesian& b) const
		{
			Vector3Cartesian ret;
			for (int i = 0; i < 3; i++)
				ret._val[i] = _val[i] - b[i];
			return ret;
		}

		/// @brief Scalar multiplication (vector * scalar).
		/// @param b Scalar multiplier
		Vector3Cartesian operator*(Real b) const
		{
			Vector3Cartesian ret;
			for (int i = 0; i < 3; i++)
				ret._val[i] = _val[i] * b;
			return ret;
		}
		/// @brief Scalar division (vector / scalar).
		/// @param b Scalar divisor
		Vector3Cartesian operator/(Real b) const
		{
			Vector3Cartesian ret;
			for (int i = 0; i < 3; i++)
				ret._val[i] = _val[i] / b;
			return ret;
		}
		/// @brief Scalar multiplication (scalar * vector).
		/// @param a Scalar multiplier
		/// @param b Vector
		friend Vector3Cartesian operator*(Real a, const Vector3Cartesian& b)
		{
			Vector3Cartesian ret;
			for (int i = 0; i < 3; i++)
				ret._val[i] = a * b[i];
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
		bool IsEqualTo(const Vector3Cartesian& b, Real absEps = Defaults::Vec3CartIsEqualTolerance) const
		{
			return (std::abs(X() - b.X()) < absEps) && (std::abs(Y() - b.Y()) < absEps) && (std::abs(Z() - b.Z()) < absEps);
		}

		/// @brief Point + vector operation.
		friend Point3Cartesian operator+(const Point3Cartesian& a, const Vector3Cartesian& b) { return Point3Cartesian(a.X() + b[0], a.Y() + b[1], a.Z() + b[2]); }
		/// @brief Point - vector operation.
		friend Point3Cartesian operator-(const Point3Cartesian& a, const Vector3Cartesian& b) { return Point3Cartesian(a.X() - b[0], a.Y() - b[1], a.Z() - b[2]); }

		/// @brief Checks if vectors are parallel.
		/// @param b Vector to compare
		/// @param eps Tolerance
		bool IsParallelTo(const Vector3Cartesian& b, Real eps = Defaults::Vec3CartIsParallelTolerance) const
		{
			Real norm1 = NormL2();
			Real norm2 = b.NormL2();

			return std::abs(X() / norm1 - b.X() / norm2) < eps &&
				std::abs(Y() / norm1 - b.Y() / norm2) < eps &&
				std::abs(Z() / norm1 - b.Z() / norm2) < eps;
		}
		/// @brief Checks if vectors are perpendicular.
		/// @param b Vector to compare
		/// @param eps Tolerance
		bool IsPerpendicularTo(const Vector3Cartesian& b, Real eps = Defaults::OrthogonalityTolerance) const
		{
			if (std::abs(this->ScalarProduct(b)) < eps)
				return true;
			else
				return false;
		}
		/// @brief Computes two vectors perpendicular to this vector.
		/// @param v1 First perpendicular vector (output)
		/// @param v2 Second perpendicular vector (output)
		/// @return true if successful, false if zero vector
		bool GetPerpendicularVectors(Vector3Cartesian& v1, Vector3Cartesian& v2) const
		{
			if (isZero())
				return false;
			
			// find a vector that is not parallel to this vector
			Vector3Cartesian not_parallel;
			if (std::abs(X()) <= std::abs(Y()) && std::abs(X()) <= std::abs(Z()))
				not_parallel = Vector3Cartesian(1.0, 0.0, 0.0);
			else if (std::abs(Y()) <= std::abs(X()) && std::abs(Y()) <= std::abs(Z()))
				not_parallel = Vector3Cartesian(0.0, 1.0, 0.0);
			else
				not_parallel = Vector3Cartesian(0.0, 0.0, 1.0);
			
			v1 = VectorProduct(*this, not_parallel).GetAsUnitVector();
			v2 = VectorProduct(*this, v1).GetAsUnitVector();
			
			return true;
		}

		/// @brief Computes angle between this vector and another.
		/// @param b Other vector
		/// @return Angle in radians [0, p]
		Real AngleToVector(const Vector3Cartesian& b)
		{
			Real normA = NormL2(), normB = b.NormL2();
			if (normA < std::numeric_limits<Real>::epsilon() ||
			    normB < std::numeric_limits<Real>::epsilon())
				throw DivisionByZeroError("AngleToVector: cannot compute angle with zero vector");
			Real cos_phi = this->ScalarProduct(b) / (normA * normB);
			// Clamp to [-1, 1] to prevent NaN from floating-point rounding errors
			cos_phi = std::max(Real{-1.0}, std::min(Real{1.0}, cos_phi));
			return std::acos(cos_phi);
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
		Real  R()     const { return _val[0]; }
		/// @brief Returns R component (non-const).
		Real& R()						{ return _val[0]; }
		/// @brief Returns Theta component (const).
		Real  Theta() const { return _val[1]; }
		/// @brief Returns Theta component (non-const).
		Real& Theta()				{ return _val[1]; }
		/// @brief Returns Phi component (const).
		Real  Phi()   const { return _val[2]; }
		/// @brief Returns Phi component (non-const).
		Real& Phi()					{ return _val[2]; }

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
		/// @note For spherical vectors, theta ? [0, p] and phi ? [-p, p].
		///       Phi comparison uses AnglesAreEqual() to handle wrap-around at �p.
		bool IsEqualTo(const Vector3Spherical& b, Real absEps = Defaults::Vec3SphIsEqualTolerance) const
		{
			// R and theta are straightforward comparisons
			// Phi needs angle-aware comparison for wrap-around at �p
			return (std::abs(R() - b.R()) < absEps) && 
			       (std::abs(Theta() - b.Theta()) < absEps) && 
			       AnglesAreEqual(Phi(), b.Phi(), absEps);
		}
		
		/// @brief Returns unit vector at given position (global spherical result).
		/// @param pos Position in spherical coordinates
		Vector3Spherical GetAsUnitVectorAtPosGlobalSpher(const Vector3Spherical& pos) const
		{
			// Convert this vector to Cartesian (as a displacement from the origin)
			Real r = R();
			Real theta = Theta();
			Real phi = Phi();

			// Spherical to Cartesian (displacement from origin)
			Real x = r * std::sin(theta) * std::cos(phi);
			Real y = r * std::sin(theta) * std::sin(phi);
			Real z = r * std::cos(theta);

			// Now, express this vector in the local spherical basis at 'pos'
			// The local basis at 'pos' consists of:
			//   e_r     = [sin? cosf, sin? sinf, cos?]
			//   e_theta = [cos? cosf, cos? sinf, -sin?]
			//   e_phi   = [-sinf, cosf, 0]
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
			Real norm = std::sqrt(v_r * v_r + v_theta * v_theta + v_phi * v_phi);
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
			// Step 1: Convert *this to Cartesian (as a displacement from the origin)
			Real r = R();
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
			Real norm = std::sqrt(v_r * v_r + v_theta * v_theta + v_phi * v_phi);
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
		Real  R()   const { return _val[0]; }
		/// @brief Returns R component (non-const).
		Real& R()					{ return _val[0]; }
		/// @brief Returns Phi component (const).
		Real  Phi() const { return _val[1]; }
		/// @brief Returns Phi component (non-const).
		Real& Phi()				{ return _val[1]; }
		/// @brief Returns Z component (const).
		Real  Z()   const { return _val[2]; }
		/// @brief Returns Z component (non-const).
		Real& Z()					{ return _val[2]; }

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
			return Vector3Cylindrical(R(), Phi() + Constants::PI, Z());
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
			// Step 1: Convert *this to Cartesian
			Real x = R() * std::cos(Phi());
			Real y = R() * std::sin(Phi());
			Real z = Z();

			// Step 2: Normalize in Cartesian
			Real norm = std::sqrt(x * x + y * y + z * z);
			if (norm == REAL(0.0))
				return Vector3Cylindrical(REAL(0.0), REAL(0.0), REAL(0.0));

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

	/// @brief 4D Minkowski spacetime vector with T (time), X, Y, Z components.
	class Vector4Minkowski : public VectorN<Real, 4>
	{
	public:
		/// @brief Returns T component (const).
		Real  T() const { return _val[0]; }
		/// @brief Returns T component (non-const).
		Real& T()				{ return _val[0]; }
		/// @brief Returns X component (const).
		Real  X() const { return _val[1]; }
		/// @brief Returns X component (non-const).
		Real& X()				{ return _val[1]; }
		/// @brief Returns Y component (const).
		Real  Y() const { return _val[2]; }
		/// @brief Returns Y component (non-const).
		Real& Y()				{ return _val[2]; }
		/// @brief Returns Z component (const).
		Real  Z() const { return _val[3]; }
		/// @brief Returns Z component (non-const).
		Real& Z()				{ return _val[3]; }

		/// @brief Default constructor (zero 4-vector).
		Vector4Minkowski() : VectorN<Real, 4>{ 0.0, 0.0, 0.0, 0.0 } {}
		/// @brief Constructs from initializer list.
		/// @param list Initializer list of values
		Vector4Minkowski(std::initializer_list<Real> list) : VectorN<Real, 4>(list) { }

		/// @brief Minkowski scalar product (T*T - X*X - Y*Y - Z*Z).
		/// @param a First 4-vector
		/// @param b Second 4-vector
		friend Real ScalarProduct(const Vector4Minkowski& a, const Vector4Minkowski& b)
		{
			return a.T() * b.T() - a.X() * b.X() - a.Y() * b.Y() - a.Z() * b.Z();
		}

		/// @brief Computes Minkowski norm (spacetime interval).
		Real Norm() const
		{
			if(T() == 0.0 && X() == 0.0 && Y() == 0.0 && Z() == 0.0)
				return 0.0;

			// Minkowski norm: sqrt(T^2 - X^2 - Y^2 - Z^2)
			if (isTimelike())
				return sqrt(T() * T() - X() * X() - Y() * Y() - Z() * Z());
			else if (isSpacelike())
				return -sqrt(X() * X() + Y() * Y() + Z() * Z() - T() * T());
			else 
				return 0.0; // Lightlike vectors have zero Minkowski norm
		}

		/// @brief Computes Minkowski distance to another 4-vector.
		/// @param b Other 4-vector
		Real Distance(const Vector4Minkowski& b) const
		{
			return sqrt((T() - b.T()) * (T() - b.T()) - (X() - b.X()) * (X() - b.X()) -
								(Y() - b.Y()) * (Y() - b.Y()) - (Z() - b.Z()) * (Z() - b.Z()));
		}
		
		/// @brief Checks if 4-vector is timelike (inside light cone).
		bool isTimelike() const
		{
			return (-T() * T() + X() * X() + Y() * Y() + Z() * Z()) < 0;
		}
		/// @brief Checks if 4-vector is spacelike (outside light cone).
		bool isSpacelike() const
		{
			return (-T() * T() + X() * X() + Y() * Y() + Z() * Z()) > 0;
		}
		/// @brief Checks if 4-vector is lightlike (on light cone).
		bool isLightlike() const
		{
			return (T() * T() - X() * X() - Y() * Y() - Z() * Z()) == 0;
		}
	};

	typedef Vector2Cartesian    Vec2Cart;
	typedef Vector2Polar				Vec2Pol;
	typedef Vector3Cartesian    Vec3Cart;
	typedef Vector3Spherical    Vec3Sph;
	typedef Vector3Cylindrical  Vec3Cyl;
	typedef Vector4Minkowski    Vec4Mink;
}
#endif

///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        VectorTypes2D.h                                                     ///
///  Description: 2D Cartesian and polar vector types and aliases                     ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_VECTOR_TYPES_2D_H
#define MML_VECTOR_TYPES_2D_H

#include <mml/MMLBase.h>
#include <mml/MMLExceptions.h>

#include <mml/base/Vector/VectorN.h>
#include <mml/base/Geometry/GeometryBase/GeometryPoints.h>

#include <cmath>
#include <limits>

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
			val(0) = x;
			val(1) = y;
		}
		/// @brief Constructs from VectorN<Real,2>.
		/// @param b Vector to copy
		Vector2Cartesian(const VectorN<Real, 2>& b) : VectorN<Real, 2>{ b[0], b[1] } {}
		/// @brief Constructs vector from two points (b - a).
		/// @param a Start point
		/// @param b End point
		Vector2Cartesian(const Point2Cartesian& a, const Point2Cartesian& b)
		{
			val(0) = b.X() - a.X();
			val(1) = b.Y() - a.Y();
		}
		/// @brief Constructs from initializer list.
		/// @param list Initializer list of values
		Vector2Cartesian(std::initializer_list<Real> list) : VectorN<Real, 2>(list) {}

		/// @brief Returns X component (const).
		Real  X() const { return val(0); }
		/// @brief Returns X component (non-const).
		Real& X()				{ return val(0); }
		/// @brief Returns Y component (const).
		Real  Y() const { return val(1); }
		/// @brief Returns Y component (non-const).
		Real& Y()				{ return val(1); }
		
		/// @brief Unary minus (negation).
		Vector2Cartesian operator-() const
		{
			Vector2Cartesian ret;
			for (int i = 0; i < 2; i++)
				ret[i] = -val(i);
			return ret;
		}
		
		/// @brief Vector addition.
		/// @param b Vector to add
		Vector2Cartesian operator+(const Vector2Cartesian& b) const
		{
			Vector2Cartesian ret;
			for (int i = 0; i < 2; i++)
				ret[i] = val(i) + b[i];
			return ret;
		}
		/// @brief Vector subtraction.
		/// @param b Vector to subtract
		Vector2Cartesian operator-(const Vector2Cartesian& b) const
		{
			Vector2Cartesian ret;
			for (int i = 0; i < 2; i++)
				ret[i] = val(i) - b[i];
			return ret;
		}

		/// @brief Scalar multiplication (vector * scalar).
		/// @param b Scalar multiplier
		Vector2Cartesian operator*(Real b) const
		{
			Vector2Cartesian ret;
			for (int i = 0; i < 2; i++)
				ret[i] = val(i) * b;
			return ret;
		}
		/// @brief Scalar division (vector / scalar).
		/// @param b Scalar divisor
		Vector2Cartesian operator/(Real b) const
		{
			Vector2Cartesian ret;
			for (int i = 0; i < 2; i++)
				ret[i] = val(i) / b;
			return ret;
		}
		/// @brief Scalar multiplication (scalar * vector).
		/// @param a Scalar multiplier
		/// @param b Vector
		friend Vector2Cartesian operator*(Real a, const Vector2Cartesian& b)
		{
			Vector2Cartesian ret;
			for (int i = 0; i < 2; i++)
				ret[i] = a * b[i];
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
		bool IsEqualTo(const Vector2Cartesian& b, Real absEps = Precision::Vec2CartIsEqualTolerance) const
		{
			if (!std::isfinite(absEps) || absEps < 0)
				return false;
			return (std::abs(X() - b.X()) < absEps) && (std::abs(Y() - b.Y()) < absEps);
		}

		/// @brief Returns normalized vector (unit length).
		/// @throws VectorDimensionError if the vector norm is near zero
		Vector2Cartesian Normalized() const
		{
			Real norm = NormL2();
			if (norm < Precision::RoundoffZeroThreshold)
				throw VectorDimensionError("Vector2Cartesian::Normalized - cannot normalize near-zero vector", 2, 0);
			return Vector2Cartesian{ (*this) / norm };
		}
		/// @brief Returns normalized vector, or the zero vector when norm is near zero.
		Vector2Cartesian NormalizedOrZero() const
		{
			Real norm = NormL2();
			if (norm < Precision::RoundoffZeroThreshold)
				return Vector2Cartesian{ static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0) };
			return Vector2Cartesian{ (*this) / norm };
		}
		/// @deprecated Use Normalized() — identical semantics (CONTRACT_POLICY §11: one spelling per operation).
		[[deprecated("use Normalized()")]]
		Vector2Cartesian GetAsUnitVector() const { return Normalized(); }
		/// @deprecated Use NormalizedOrZero().
		[[deprecated("use NormalizedOrZero()")]]
		Vector2Cartesian GetAsUnitVectorOrZero() const { return NormalizedOrZero(); }
		/// @brief Returns unit vector at given position (for Cartesian, position doesn't matter).
		/// @param pos Position (unused in Cartesian coordinates)
		Vector2Cartesian GetAsUnitVectorAtPos(const Vector2Cartesian& pos) const
		{
			return Normalized();
		}
		
		/// @brief Computes two perpendicular vectors.
		/// @param v1 First perpendicular vector (output)
		/// @param v2 Second perpendicular vector (output)
		void GetPerpendicularVectors(Vector2Cartesian& v1, Vector2Cartesian& v2) const
		{
			v1 = Vector2Cartesian(-Y(), X());
			v2 = Vector2Cartesian(Y(), -X());
		}
		/// @deprecated Use GetPerpendicularVectors() (same spelling as Vector3Cartesian).
		[[deprecated("use GetPerpendicularVectors()")]]
		void getPerpendicularVectors(Vector2Cartesian& v1, Vector2Cartesian& v2) const { GetPerpendicularVectors(v1, v2); }

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
		Vector2Polar Combine(const Vector2Polar& b, Real sign) const
		{
			if (!std::isfinite(R()) || !std::isfinite(Phi()) || !std::isfinite(b.R()) || !std::isfinite(b.Phi()))
				throw NumericInputError("Vector2Polar arithmetic requires finite radii and angles");
			const Real scale = std::max(std::abs(R()), std::abs(b.R()));
			if (scale == 0)
				return Vector2Polar(0, Phi());
			Real delta = b.Phi() - Phi();
			if (!std::isfinite(delta) || std::abs(delta) > Constants::PI)
				delta = std::atan2(std::sin(b.Phi()) * std::cos(Phi()) - std::cos(b.Phi()) * std::sin(Phi()),
					std::cos(b.Phi()) * std::cos(Phi()) + std::sin(b.Phi()) * std::sin(Phi()));
			const Real half = delta / Real{2};
			const Real first = R() / scale;
			const Real second = sign * (b.R() / scale);
			// A midpoint frame avoids the 1-cos(delta) cancellation of near-equal vectors.
			const Real x = (first + second) * std::cos(half);
			const Real y = (second - first) * std::sin(half);
			const Real angle = std::atan2(std::sin(Phi()), std::cos(Phi())) + half + std::atan2(y, x);
			return Vector2Polar(scale * std::hypot(x, y), std::atan2(std::sin(angle), std::cos(angle)));
		}

	public:
		/// @brief Returns R component (const).
		Real  R() const		{ return val(0); }
		/// @brief Returns R component (non-const).
		Real& R()					{ return val(0); }
		/// @brief Returns Phi component (const).
		Real  Phi() const { return val(1); }
		/// @brief Returns Phi component (non-const).
		Real& Phi()				{ return val(1); }

		/// @brief Default constructor (zero vector).
		Vector2Polar() {}
		/// @brief Constructs from R, Phi components.
		/// @param r Radius
		/// @param phi Angle
		Vector2Polar(Real r, Real phi)
		{
			val(0) = r;
			val(1) = phi;
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
			return Combine(b, Real{1});
		}
		/// @brief Vector subtraction (converts to Cartesian, subtracts, converts back).
		/// @param b Vector to subtract
		Vector2Polar operator-(const Vector2Polar& b) const
		{
			return Combine(b, Real{-1});
		}
		Vector2Polar& operator+=(const Vector2Polar& b) { *this = *this + b; return *this; }
		Vector2Polar& operator-=(const Vector2Polar& b) { *this = *this - b; return *this; }
		Vector2Polar& operator*=(Real b) { *this = *this * b; return *this; }
		Vector2Polar& operator/=(Real b) { *this = *this / b; return *this; }
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

	typedef Vector2Cartesian Vec2Cart;
	typedef Vector2Polar     Vec2Pol;
}

#endif

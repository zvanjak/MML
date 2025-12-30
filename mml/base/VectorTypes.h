///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        VectorTypes.h                                                       ///
///  Description: Type aliases for 2D, 3D, 4D vectors (Vec2, Vec3, Vec4)              ///
///               Complex vector types and related utilities                          ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_VECTOR_TYPES_H
#define MML_VECTOR_TYPES_H

#include "MMLBase.h"

#include "base/VectorN.h"
#include "base/Geometry.h"

namespace MML
{
	class Vector2Cartesian : public VectorN<Real, 2>
	{
	public:
		Vector2Cartesian() {}
		Vector2Cartesian(Real x, Real y)
		{
			_val[0] = x;
			_val[1] = y;
		}
		Vector2Cartesian(const VectorN<Real, 2>& b) : VectorN<Real, 2>{ b[0], b[1] } {}
		Vector2Cartesian(const Point2Cartesian& a, const Point2Cartesian& b)
		{
			_val[0] = b.X() - a.X();
			_val[1] = b.Y() - a.Y();
		}
		Vector2Cartesian(std::initializer_list<Real> list) : VectorN<Real, 2>(list) {}

		Real  X() const { return _val[0]; }
		Real& X()				{ return _val[0]; }
		Real  Y() const { return _val[1]; }
		Real& Y()				{ return _val[1]; }
		
		// unary minus operator
		Vector2Cartesian operator-() const
		{
			Vector2Cartesian ret;
			for (int i = 0; i < 2; i++)
				ret._val[i] = -_val[i];
			return ret;
		}
		
		// arithmetic operations
		Vector2Cartesian operator+(const Vector2Cartesian& b) const
		{
			Vector2Cartesian ret;
			for (int i = 0; i < 2; i++)
				ret._val[i] = _val[i] + b[i];
			return ret;
		}
		Vector2Cartesian operator-(const Vector2Cartesian& b) const
		{
			Vector2Cartesian ret;
			for (int i = 0; i < 2; i++)
				ret._val[i] = _val[i] - b[i];
			return ret;
		}

		Vector2Cartesian operator*(Real b) const
		{
			Vector2Cartesian ret;
			for (int i = 0; i < 2; i++)
				ret._val[i] = _val[i] * b;
			return ret;
		}
		Vector2Cartesian operator/(Real b) const
		{
			Vector2Cartesian ret;
			for (int i = 0; i < 2; i++)
				ret._val[i] = _val[i] / b;
			return ret;
		}
		friend Vector2Cartesian operator*(Real a, const Vector2Cartesian& b)
		{
			Vector2Cartesian ret;
			for (int i = 0; i < 2; i++)
				ret._val[i] = a * b[i];
			return ret;
		}

		// equality operators
		bool operator==(const Vector2Cartesian& b) const
		{
			return (X() == b.X()) && (Y() == b.Y());
		}
		bool operator!=(const Vector2Cartesian& b) const
		{
			return (X() != b.X()) || (Y() != b.Y());
		}
		bool IsEqual(const Vector2Cartesian& b, Real absEps = Defaults::Vec2CartIsEqualTolerance) const
		{
			return (std::abs(X() - b.X()) < absEps) && (std::abs(Y() - b.Y()) < absEps);
		}

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
		Vector2Cartesian GetAsUnitVector() const
		{
			VectorN<Real, 2> res = (*this) / NormL2();

			return Vector2Cartesian(res[0], res[1]);
		}
		Vector2Cartesian GetAsUnitVectorAtPos(const Vector2Cartesian& pos) const
		{
			return Vector2Cartesian{ (*this) / NormL2() };
		}
		
		// returns two vectors that are perpendicular to this vector
		void getPerpendicularVectors(Vector2Cartesian& v1, Vector2Cartesian& v2) const
		{
			v1 = Vector2Cartesian(-Y(), X());
			v2 = Vector2Cartesian(Y(), -X());
		}

		// For Cartesian vector, we will enable operator* to represent standard scalar product
		Real operator*(const Vector2Cartesian& b) const
		{
			return X() * b.X() + Y() * b.Y();
		}

		friend Real ScalarProduct(const Vector2Cartesian& a, const Vector2Cartesian& b)
		{
			return a * b;
		}

		friend Point2Cartesian operator+(const Point2Cartesian& a, const Vector2Cartesian& b) { return Point2Cartesian(a.X() + b[0], a.Y() + b[1]); }
		friend Point2Cartesian operator-(const Point2Cartesian& a, const Vector2Cartesian& b) { return Point2Cartesian(a.X() - b[0], a.Y() - b[1]); }
	};

	class Vector2Polar : public VectorN<Real, 2>
	{
	public:
		Real  R() const		{ return _val[0]; }
		Real& R()					{ return _val[0]; }
		Real  Phi() const { return _val[1]; }
		Real& Phi()				{ return _val[1]; }

		Vector2Polar() {}
		Vector2Polar(Real r, Real phi)
		{
			_val[0] = r;
			_val[1] = phi;
		}
		Vector2Polar(const VectorN<Real, 2>& b) : VectorN<Real, 2>{ b[0], b[1] } {}

		// unary minus operator
		Vector2Polar operator-() const
		{
			// Negating a polar vector means keeping the radius and adding pi to the angle
			return Vector2Polar(R(), Phi() + Constants::PI);
		}
		
		// arithmetic operations
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
		Vector2Polar operator*(Real b) const
		{
			return Vector2Polar(R() * b, Phi());
		}
		Vector2Polar operator/(Real b) const
		{
			return Vector2Polar(R() / b, Phi());
		}

		friend Vector2Polar operator*(Real a, const Vector2Polar& b)
		{
			return Vector2Polar(b.R() * a, b.Phi());
		}

		Vector2Polar GetAsUnitVectorAtPos(const Vector2Polar& /*pos*/) const
		{
			// Returns a unit vector in the direction of this vector (ignores pos)
			return Vector2Polar(1.0, Phi());
		}
	};

	class Vector3Cartesian : public VectorN<Real, 3>
	{
	public:
		Real  X() const { return _val[0]; }
		Real& X()				{ return _val[0]; }
		Real  Y() const { return _val[1]; }
		Real& Y()				{ return _val[1]; }
		Real  Z() const { return _val[2]; }
		Real& Z()				{ return _val[2]; }

		Vector3Cartesian() : VectorN<Real, 3>{ 0.0, 0.0, 0.0 } {}
		Vector3Cartesian(const VectorN<Real, 3>& b) : VectorN<Real, 3>{ b } {}
		Vector3Cartesian(Real x, Real y, Real z) : VectorN<Real, 3>{ x, y, z } {}
		Vector3Cartesian(std::initializer_list<Real> list) : VectorN<Real, 3>(list) { }
		Vector3Cartesian(const Point3Cartesian& a, const Point3Cartesian& b)
		{
			_val[0] = b.X() - a.X();
			_val[1] = b.Y() - a.Y();
			_val[2] = b.Z() - a.Z();
		}
		// we'll provide a constructor that takes a Point3Cartesian, 
		// so it can be used to convert a Cartesian point to a vector
		Vector3Cartesian(const Point3Cartesian& a)
		{
			_val[0] = a.X();
			_val[1] = a.Y();
			_val[2] = a.Z();
		}

		Point3Cartesian getAsPoint()
		{
			return Point3Cartesian(_val[0], _val[1], _val[2]);
		}
		
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
		Vector3Cartesian GetAsUnitVectorAtPos(const Vector3Cartesian& pos) const
		{
			return GetAsUnitVector();
		}

		// For Cartesian vector, we will enable operator* to represent standard scalar product
		Real operator*(const Vector3Cartesian& b) const
		{
			return X()*b.X() + Y()*b.Y() + Z()*b.Z();
		}

		// unary minus operator
		Vector3Cartesian operator-() const
		{
			Vector3Cartesian ret;
			for (int i = 0; i < 3; i++)
				ret._val[i] = -_val[i];
			return ret;
		}

		// arithmetic operations
		Vector3Cartesian operator+(const Vector3Cartesian& b) const
		{
			Vector3Cartesian ret;
			for (int i = 0; i < 3; i++)
				ret._val[i] = _val[i] + b[i];
			return ret;
		}
		Vector3Cartesian operator-(const Vector3Cartesian& b) const
		{
			Vector3Cartesian ret;
			for (int i = 0; i < 3; i++)
				ret._val[i] = _val[i] - b[i];
			return ret;
		}

		Vector3Cartesian operator*(Real b) const
		{
			Vector3Cartesian ret;
			for (int i = 0; i < 3; i++)
				ret._val[i] = _val[i] * b;
			return ret;
		}
		Vector3Cartesian operator/(Real b) const
		{
			Vector3Cartesian ret;
			for (int i = 0; i < 3; i++)
				ret._val[i] = _val[i] / b;
			return ret;
		}
		friend Vector3Cartesian operator*(Real a, const Vector3Cartesian& b)
		{
			Vector3Cartesian ret;
			for (int i = 0; i < 3; i++)
				ret._val[i] = a * b[i];
			return ret;
		}

		// equality operators
		bool operator==(const Vector3Cartesian& b) const
		{
			return (X() == b.X()) && (Y() == b.Y()) && (Z() == b.Z());
		}
		bool operator!=(const Vector3Cartesian& b) const
		{
			return (X() != b.X()) || (Y() != b.Y()) || (Z() != b.Z());
		}
		bool IsEqual(const Vector3Cartesian& b, Real absEps = Defaults::Vec3CartIsEqualTolerance) const
		{
			return (std::abs(X() - b.X()) < absEps) && (std::abs(Y() - b.Y()) < absEps) && (std::abs(Z() - b.Z()) < absEps);
		}

		friend Point3Cartesian operator+(const Point3Cartesian& a, const Vector3Cartesian& b) { return Point3Cartesian(a.X() + b[0], a.Y() + b[1], a.Z() + b[2]); }
		friend Point3Cartesian operator-(const Point3Cartesian& a, const Vector3Cartesian& b) { return Point3Cartesian(a.X() - b[0], a.Y() - b[1], a.Z() - b[2]); }

		bool IsParallelTo(const Vector3Cartesian& b, Real eps = Defaults::Vec3CartIsParallelTolerance) const
		{
			Real norm1 = NormL2();
			Real norm2 = b.NormL2();

			return std::abs(X() / norm1 - b.X() / norm2) < eps &&
				std::abs(Y() / norm1 - b.Y() / norm2) < eps &&
				std::abs(Z() / norm1 - b.Z() / norm2) < eps;
		}
		bool IsPerpendicularTo(const Vector3Cartesian& b, Real eps = 1e-15) const
		{
			if (std::abs(this->ScalarProduct(b)) < eps)
				return true;
			else
				return false;
		}
		bool GetPerpendicularVectors(Vector3Cartesian& v1, Vector3Cartesian& v2) const
		{
			if (IsNullVec())
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

		Real AngleToVector(const Vector3Cartesian& b)
		{
			Real cos_phi = this->ScalarProduct(b) / (NormL2() * b.NormL2());

			return acos(cos_phi);
		}

		Real ScalarProduct(const Vector3Cartesian& b) const
		{
			return (*this) * b;
		}

		friend Real ScalarProduct(const Vector3Cartesian& a, const Vector3Cartesian& b)
		{
			return a * b;
		}

		friend Vector3Cartesian VectorProduct(const Vector3Cartesian& a, const Vector3Cartesian& b)
		{
			Vector3Cartesian ret;

			ret.X() = a.Y() * b.Z() - a.Z() * b.Y();
			ret.Y() = a.Z() * b.X() - a.X() * b.Z();
			ret.Z() = a.X() * b.Y() - a.Y() * b.X();

			return ret;
		}
	};

	class Vector3Spherical : public VectorN<Real, 3>
	{
	public:
		Real  R()     const { return _val[0]; }
		Real& R()						{ return _val[0]; }
		Real  Theta() const { return _val[1]; }
		Real& Theta()				{ return _val[1]; }
		Real  Phi()   const { return _val[2]; }
		Real& Phi()					{ return _val[2]; }

		Vector3Spherical() : VectorN<Real, 3>{ 0.0, 0.0, 0.0 } {}
		Vector3Spherical(const VectorN<Real, 3>& b) : VectorN<Real, 3>{ b[0], b[1], b[2] } {}
		Vector3Spherical(Real r, Real theta, Real phi) : VectorN<Real, 3>{ r, theta, phi } {}
		Vector3Spherical(std::initializer_list<Real> list) : VectorN<Real, 3>(list) { }

		// unary minus operator
		Vector3Spherical operator-() const
		{
			return Vector3Spherical(R(), Constants::PI - Theta(), Phi() + Constants::PI);
		}

		// arithmetic operations
		Vector3Spherical operator+(const Vector3Spherical& b) const
		{
			Real r1 = R();
			Real theta1 = Theta();
			Real phi1 = Phi();
			Real r2 = b.R();
			Real theta2 = b.Theta();
			Real phi2 = b.Phi();

			// Convert spherical to Cartesian coordinates
			Real x1 = r1 * std::sin(theta1) * std::cos(phi1);
			Real y1 = r1 * std::sin(theta1) * std::sin(phi1);
			Real z1 = r1 * std::cos(theta1);
			Real x2 = r2 * std::sin(theta2) * std::cos(phi2);
			Real y2 = r2 * std::sin(theta2) * std::sin(phi2);
			Real z2 = r2 * std::cos(theta2);
			
			// Add Cartesian coordinates
			Real x = x1 + x2;
			Real y = y1 + y2;
			Real z = z1 + z2;
			
			// Convert back to spherical coordinates
			Real r = std::sqrt(x * x + y * y + z * z);
			if (r == 0.0)
				return Vector3Spherical(0.0, 0.0, 0.0);
			Real theta = std::acos(z / r);
			Real phi = std::atan2(y, x);
			
			return Vector3Spherical(r, theta, phi);
		}
		Vector3Spherical operator-(const Vector3Spherical& b) const
		{
			Real r1 = R();
			Real theta1 = Theta();
			Real phi1 = Phi();
			Real r2 = b.R();
			Real theta2 = b.Theta();
			Real phi2 = b.Phi();
			// Convert spherical to Cartesian coordinates
			Real x1 = r1 * std::sin(theta1) * std::cos(phi1);
			Real y1 = r1 * std::sin(theta1) * std::sin(phi1);
			Real z1 = r1 * std::cos(theta1);
			Real x2 = r2 * std::sin(theta2) * std::cos(phi2);
			Real y2 = r2 * std::sin(theta2) * std::sin(phi2);
			Real z2 = r2 * std::cos(theta2);
			// Subtract Cartesian coordinates
			Real x = x1 - x2;
			Real y = y1 - y2;
			Real z = z1 - z2;
			// Convert back to spherical coordinates
			Real r = std::sqrt(x * x + y * y + z * z);
			if (r == 0.0)
				return Vector3Spherical(0.0, 0.0, 0.0);
			Real theta = std::acos(z / r);
			Real phi = std::atan2(y, x);
			return Vector3Spherical(r, theta, phi);
		}

		Vector3Spherical operator*(Real b) const
		{
			return Vector3Spherical(R() * b, Theta(), Phi());
		}
		Vector3Spherical operator/(Real b) const
		{
			if (b == 0.0)
				throw std::runtime_error("Division by zero in Vector3Spherical division.");
			return Vector3Spherical(R() / b, Theta(), Phi());
		}
		friend Vector3Spherical operator*(Real a, const Vector3Spherical& b)
		{
			return Vector3Spherical(b.R() * a, b.Theta(), b.Phi());
		}

		// equality operators
		bool operator==(const Vector3Spherical& b) const
		{
			return (R() == b.R()) && (Theta() == b.Theta()) && (Phi() == b.Phi());
		}
		bool operator!=(const Vector3Spherical& b) const
		{
			return (R() != b.R()) || (Theta() != b.Theta()) || (Phi() != b.Phi());
		}
		bool IsEqual(const Vector3Spherical& b, Real absEps = Defaults::Vec3SphIsEqualTolerance) const
		{
			return (std::abs(R() - b.R()) < absEps) && (std::abs(Theta() - b.Theta()) < absEps) && (std::abs(Phi() - b.Phi()) < absEps);
		}
		
		// for spherical vector with, with components defined in local spherical coordinate system
		// at given positon 'pos', returns the same vector but with unit length
		// and result is in spherical coordinates
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
			Real norm = std::sqrt(v_r * v_r + v_theta * v_theta + v_phi * v_phi);
			if (norm == 0.0)
				return Vector3Spherical(0.0, 0.0, 0.0);

			return Vector3Spherical(1.0, std::acos(v_r / norm), std::atan2(v_phi, v_theta));
		}
		// for spherical vector with, with components defined in local spherical coordinate system
		// at given positon 'pos', returns the same vector but with unit length
		// and result is in local spherical basis at 'pos'
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

	class Vector3Cylindrical : public VectorN<Real, 3>
	{
	public:
		Real  R()   const { return _val[0]; }
		Real& R()					{ return _val[0]; }
		Real  Phi() const { return _val[1]; }
		Real& Phi()				{ return _val[1]; }
		Real  Z()   const { return _val[2]; }
		Real& Z()					{ return _val[2]; }

		Vector3Cylindrical() : VectorN<Real, 3>{ 0.0, 0.0, 0.0 } {}
		Vector3Cylindrical(const VectorN<Real, 3>& b) : VectorN<Real, 3>{ b[0], b[1], b[2] } {}
		Vector3Cylindrical(Real r, Real phi, Real z) : VectorN<Real, 3>{ r, phi, z } {}
		Vector3Cylindrical(std::initializer_list<Real> list) : VectorN<Real, 3>(list) { }

		// unary minus operator
		Vector3Cylindrical operator-() const
		{
			return Vector3Cylindrical(R(), Phi() + Constants::PI, Z());
		}

		// arithmetic operations
		Vector3Cylindrical operator+(const Vector3Cylindrical& b) const
		{
			Real r1 = R();
			Real phi1 = Phi();
			Real z1 = Z();
			Real r2 = b.R();
			Real phi2 = b.Phi();
			Real z2 = b.Z();
			// Convert cylindrical to Cartesian coordinates
			Real x1 = r1 * std::cos(phi1);
			Real y1 = r1 * std::sin(phi1);
			Real x2 = r2 * std::cos(phi2);
			Real y2 = r2 * std::sin(phi2);
			// Add Cartesian coordinates
			Real x = x1 + x2;
			Real y = y1 + y2;
			Real z = z1 + z2;
			// Convert back to cylindrical coordinates
			Real r = std::sqrt(x * x + y * y);
			if (r == 0.0)
				return Vector3Cylindrical(0.0, 0.0, z);
			Real phi = std::atan2(y, x);
			return Vector3Cylindrical(r, phi, z);
		}
		Vector3Cylindrical operator-(const Vector3Cylindrical& b) const
		{
			Real r1 = R();
			Real phi1 = Phi();
			Real z1 = Z();
			Real r2 = b.R();
			Real phi2 = b.Phi();
			Real z2 = b.Z();
			// Convert cylindrical to Cartesian coordinates
			Real x1 = r1 * std::cos(phi1);
			Real y1 = r1 * std::sin(phi1);
			Real x2 = r2 * std::cos(phi2);
			Real y2 = r2 * std::sin(phi2);
			// Subtract Cartesian coordinates
			Real x = x1 - x2;
			Real y = y1 - y2;
			Real z = z1 - z2;
			// Convert back to cylindrical coordinates
			Real r = std::sqrt(x * x + y * y);
			if (r == 0.0)
				return Vector3Cylindrical(0.0, 0.0, z);
			Real phi = std::atan2(y, x);
			return Vector3Cylindrical(r, phi, z);
		}

		Vector3Cylindrical operator*(Real b) const
		{
			return Vector3Cylindrical(R() * b, Phi(), Z() * b);
		}
		Vector3Cylindrical operator/(Real b) const
		{
			if (b == 0.0)
				throw std::runtime_error("Division by zero in Vector3Cylindrical division.");
			return Vector3Cylindrical(R() / b, Phi(), Z() / b);
		}
		friend Vector3Cylindrical operator*(Real a, const Vector3Cylindrical& b)
		{
			return Vector3Cylindrical(b.R() * a, b.Phi(), b.Z() * a);
		}

		Vector3Cylindrical GetAsUnitVectorAtPos(const Vector3Cylindrical& pos) const
		{
			return Vector3Cylindrical{ R(), Phi() / pos.R(), Z() };
		}
	};

	class Vector4Minkowski : public VectorN<Real, 4>
	{
	public:
		Real  T() const { return _val[0]; }
		Real& T()				{ return _val[0]; }
		Real  X() const { return _val[1]; }
		Real& X()				{ return _val[1]; }
		Real  Y() const { return _val[2]; }
		Real& Y()				{ return _val[2]; }
		Real  Z() const { return _val[3]; }
		Real& Z()				{ return _val[3]; }

		Vector4Minkowski() : VectorN<Real, 4>{ 0.0, 0.0, 0.0, 0.0 } {}
		Vector4Minkowski(std::initializer_list<Real> list) : VectorN<Real, 4>(list) { }

		friend Real ScalarProduct(const Vector4Minkowski& a, const Vector4Minkowski& b)
		{
			return a.T() * b.T() - a.X() * b.X() - a.Y() * b.Y() - a.Z() * b.Z();
		}

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

		Real Distance(const Vector4Minkowski& b) const
		{
			return sqrt((T() - b.T()) * (T() - b.T()) - (X() - b.X()) * (X() - b.X()) -
									(Y() - b.Y()) * (Y() - b.Y()) - (Z() - b.Z()) * (Z() - b.Z()));
		}
		
		bool isTimelike() const
		{
			return (-T() * T() + X() * X() + Y() * Y() + Z() * Z()) < 0;
		}
		bool isSpacelike() const
		{
			return (-T() * T() + X() * X() + Y() * Y() + Z() * Z()) > 0;
		}
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

///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Quaternions.h                                                       ///
///  Description: Quaternion class for 3D rotations and orientation                   ///
///               SLERP interpolation, conversion to/from Euler angles                ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_QUATERNIONS_H
#define MML_QUATERNIONS_H

#include "MMLBase.h"
#include "base/Vector/Vector.h"
#include "base/Vector/VectorTypes.h"
#include "base/Matrix/MatrixNM.h"

// Standard headers - include what we use
#include <iostream>

namespace MML
{
	///////////////////////////////////////////////////////////////////////////////
	// QUATERNION CONVENTIONS USED IN THIS FILE
	//
	//   Algebra:       Hamilton convention (ij = k, ji = -k)
	//                  NOT JPL/aerospace convention (ij = -k)
	//
	//   Storage order:  [w, x, y, z]  where w = scalar part, (x,y,z) = vector part
	//                  Identity quaternion = [1, 0, 0, 0]
	//
	//   Rotation:      Active (alibi) rotation via  v' = q * [0,v] * q⁻¹
	//                  A quaternion q = [cos(α/2), sin(α/2)·û] rotates by angle α
	//                  about axis û using the right-hand rule.
	//
	//   Angle units:   All angles in RADIANS (axis-angle, Euler angles, SLERP)
	//
	//   Euler angles:  Tait-Bryan (three different axes):
	//                    FromEulerXYZ(roll, pitch, yaw)  → R = Rx·Ry·Rz
	//                    FromEulerXZY(x, z, y)           → R = Rx·Rz·Ry
	//                    FromEulerYXZ(y, x, z)           → R = Ry·Rx·Rz
	//                    FromEulerYZX(y, z, x)           → R = Ry·Rz·Rx
	//                    FromEulerZXY(z, x, y)           → R = Rz·Rx·Ry
	//                    FromEulerZYX(yaw, pitch, roll)  → R = Rz·Ry·Rx
	//                  Proper Euler (first and last axis same):
	//                    FromEulerZXZ(alpha, beta, gamma) → R = Rz·Rx·Rz
	//                    FromEulerZYZ(alpha, beta, gamma) → R = Rz·Ry·Rz
	//                  Note: parameter ORDER matches the name, not the
	//                  multiplication order.
	//
	//   Rotation matrix: ToRotationMatrix() returns R such that R*v rotates
	//                  the vector (column-vector convention, active rotation).
	//
	//   SLERP:         Takes shortest path — negates q2 when dot(q1,q2) < 0.
	//
	// See also: CoordTransf.h, CoordTransfSpherical.h
	///////////////////////////////////////////////////////////////////////////////

	/// @brief Quaternion class for 3D rotations and spatial orientation
	/// @details Implements quaternions as q = w + xi + yj + zk with Hamilton multiplication.
	///          Storage order: [w, x, y, z] where w is scalar, (x,y,z) is vector part.
	///          Unit quaternions represent rotations via q*v*q^(-1) formula.
	///          Supports axis-angle, Euler angles, matrix conversions, and SLERP interpolation.
	class Quaternion
	{
	private:
		Real _data[4];  // [w, x, y, z] storage

	public:
		/// @brief Default constructor - creates identity quaternion [1, 0, 0, 0] (no rotation)
		Quaternion() : _data{1, 0, 0, 0} {}

		/// @brief Construct from components q = w + xi + yj + zk
		Quaternion(Real w, Real x, Real y, Real z) : _data{w, x, y, z} {}

		/// @brief Construct from scalar and vector parts q = w + v
		Quaternion(Real w, const Vec3Cart& vec) 
			: _data{w, vec[0], vec[1], vec[2]} {}

		/// @brief Construct pure imaginary quaternion from vector q = 0 + v
		explicit Quaternion(const Vec3Cart& vec)
			: _data{0, vec[0], vec[1], vec[2]} {}

		/// @brief Copy constructor
		Quaternion(const Quaternion& q) 
			: _data{q._data[0], q._data[1], q._data[2], q._data[3]} {}

		/// @brief Create identity quaternion (no rotation)
		static Quaternion Identity()
		{
			return Quaternion(1, 0, 0, 0);
		}

		/// @brief Create rotation quaternion from axis-angle
		/// @param axis Unit vector defining rotation axis
		/// @param angle Rotation angle in radians (right-hand rule)
		static Quaternion FromAxisAngle(const Vec3Cart& axis, Real angle)
		{
			Real halfAngle = angle / 2.0;
			Real sinHalf = std::sin(halfAngle);
			Real cosHalf = std::cos(halfAngle);
			
			return Quaternion(cosHalf, 
											 axis[0] * sinHalf,
											 axis[1] * sinHalf,
											 axis[2] * sinHalf);
		}

		/// @brief Create rotation quaternion from Euler angles (ZYX convention)
		/// @param yaw Rotation around Z-axis (radians)
		/// @param pitch Rotation around Y-axis (radians)
		/// @param roll Rotation around X-axis (radians)
		/// @details Applied as R = Rz(yaw) * Ry(pitch) * Rx(roll)
		static Quaternion FromEulerZYX(Real yaw, Real pitch, Real roll)
		{
			Real cy = std::cos(yaw * 0.5);
			Real sy = std::sin(yaw * 0.5);
			Real cp = std::cos(pitch * 0.5);
			Real sp = std::sin(pitch * 0.5);
			Real cr = std::cos(roll * 0.5);
			Real sr = std::sin(roll * 0.5);

			return Quaternion(
				cr * cp * cy + sr * sp * sy,
				sr * cp * cy - cr * sp * sy,
				cr * sp * cy + sr * cp * sy,
				cr * cp * sy - sr * sp * cy
			);
		}

		/// @brief Create rotation quaternion from Euler angles (XYZ convention)
		/// @param roll Rotation around X-axis (radians)
		/// @param pitch Rotation around Y-axis (radians)
		/// @param yaw Rotation around Z-axis (radians)
		/// @details Applied as R = Rx(roll) * Ry(pitch) * Rz(yaw)
		static Quaternion FromEulerXYZ(Real roll, Real pitch, Real yaw)
		{
			Real cr = std::cos(roll * 0.5);
			Real sr = std::sin(roll * 0.5);
			Real cp = std::cos(pitch * 0.5);
			Real sp = std::sin(pitch * 0.5);
			Real cy = std::cos(yaw * 0.5);
			Real sy = std::sin(yaw * 0.5);

			return Quaternion(
				cr * cp * cy - sr * sp * sy,
				sr * cp * cy + cr * sp * sy,
				cr * sp * cy - sr * cp * sy,
				cr * cp * sy + sr * sp * cy
			);
		}

		/// @brief Create rotation quaternion from Euler angles (XZY convention)
		/// @param xAngle Rotation around X-axis (radians)
		/// @param zAngle Rotation around Z-axis (radians)
		/// @param yAngle Rotation around Y-axis (radians)
		/// @details Applied as R = Rx(xAngle) * Rz(zAngle) * Ry(yAngle)
		static Quaternion FromEulerXZY(Real xAngle, Real zAngle, Real yAngle)
		{
			Real c1 = std::cos(xAngle * 0.5), s1 = std::sin(xAngle * 0.5);
			Real c2 = std::cos(zAngle * 0.5), s2 = std::sin(zAngle * 0.5);
			Real c3 = std::cos(yAngle * 0.5), s3 = std::sin(yAngle * 0.5);

			return Quaternion(
				c1*c2*c3 + s1*s2*s3,
				s1*c2*c3 - c1*s2*s3,
				c1*c2*s3 - s1*s2*c3,
				s1*c2*s3 + c1*s2*c3
			);
		}

		/// @brief Create rotation quaternion from Euler angles (YXZ convention)
		/// @param yAngle Rotation around Y-axis (radians)
		/// @param xAngle Rotation around X-axis (radians)
		/// @param zAngle Rotation around Z-axis (radians)
		/// @details Applied as R = Ry(yAngle) * Rx(xAngle) * Rz(zAngle)
		static Quaternion FromEulerYXZ(Real yAngle, Real xAngle, Real zAngle)
		{
			Real c1 = std::cos(yAngle * 0.5), s1 = std::sin(yAngle * 0.5);
			Real c2 = std::cos(xAngle * 0.5), s2 = std::sin(xAngle * 0.5);
			Real c3 = std::cos(zAngle * 0.5), s3 = std::sin(zAngle * 0.5);

			return Quaternion(
				c1*c2*c3 + s1*s2*s3,
				c1*s2*c3 + s1*c2*s3,
				s1*c2*c3 - c1*s2*s3,
				c1*c2*s3 - s1*s2*c3
			);
		}

		/// @brief Create rotation quaternion from Euler angles (YZX convention)
		/// @param yAngle Rotation around Y-axis (radians)
		/// @param zAngle Rotation around Z-axis (radians)
		/// @param xAngle Rotation around X-axis (radians)
		/// @details Applied as R = Ry(yAngle) * Rz(zAngle) * Rx(xAngle)
		static Quaternion FromEulerYZX(Real yAngle, Real zAngle, Real xAngle)
		{
			Real c1 = std::cos(yAngle * 0.5), s1 = std::sin(yAngle * 0.5);
			Real c2 = std::cos(zAngle * 0.5), s2 = std::sin(zAngle * 0.5);
			Real c3 = std::cos(xAngle * 0.5), s3 = std::sin(xAngle * 0.5);

			return Quaternion(
				c1*c2*c3 - s1*s2*s3,
				c1*c2*s3 + s1*s2*c3,
				s1*c2*c3 + c1*s2*s3,
				c1*s2*c3 - s1*c2*s3
			);
		}

		/// @brief Create rotation quaternion from Euler angles (ZXY convention)
		/// @param zAngle Rotation around Z-axis (radians)
		/// @param xAngle Rotation around X-axis (radians)
		/// @param yAngle Rotation around Y-axis (radians)
		/// @details Applied as R = Rz(zAngle) * Rx(xAngle) * Ry(yAngle)
		static Quaternion FromEulerZXY(Real zAngle, Real xAngle, Real yAngle)
		{
			Real c1 = std::cos(zAngle * 0.5), s1 = std::sin(zAngle * 0.5);
			Real c2 = std::cos(xAngle * 0.5), s2 = std::sin(xAngle * 0.5);
			Real c3 = std::cos(yAngle * 0.5), s3 = std::sin(yAngle * 0.5);

			return Quaternion(
				c1*c2*c3 - s1*s2*s3,
				c1*s2*c3 - s1*c2*s3,
				c1*c2*s3 + s1*s2*c3,
				c1*s2*s3 + s1*c2*c3
			);
		}

		/// @brief Create rotation quaternion from Euler angles (ZXZ proper Euler convention)
		/// @param alpha First rotation around Z-axis (radians)
		/// @param beta Rotation around X-axis (radians)
		/// @param gamma Second rotation around Z-axis (radians)
		/// @details Applied as R = Rz(alpha) * Rx(beta) * Rz(gamma). Classical physics convention.
		static Quaternion FromEulerZXZ(Real alpha, Real beta, Real gamma)
		{
			Real cb = std::cos(beta * 0.5),  sb = std::sin(beta * 0.5);
			Real cag = std::cos((alpha + gamma) * 0.5), sag = std::sin((alpha + gamma) * 0.5);
			Real cdg = std::cos((alpha - gamma) * 0.5), sdg = std::sin((alpha - gamma) * 0.5);

			return Quaternion(cb * cag, sb * cdg, sb * sdg, cb * sag);
		}

		/// @brief Create rotation quaternion from Euler angles (ZYZ proper Euler convention)
		/// @param alpha First rotation around Z-axis (radians)
		/// @param beta Rotation around Y-axis (radians)
		/// @param gamma Second rotation around Z-axis (radians)
		/// @details Applied as R = Rz(alpha) * Ry(beta) * Rz(gamma).
		static Quaternion FromEulerZYZ(Real alpha, Real beta, Real gamma)
		{
			Real cb = std::cos(beta * 0.5),  sb = std::sin(beta * 0.5);
			Real cag = std::cos((alpha + gamma) * 0.5), sag = std::sin((alpha + gamma) * 0.5);
			Real cdg = std::cos((alpha - gamma) * 0.5), sdg = std::sin((alpha - gamma) * 0.5);

			return Quaternion(cb * cag, -sb * sdg, sb * cdg, cb * sag);
		}

		/// @brief Get scalar component w (const)
		Real w() const { return _data[0]; }
		/// @brief Get x component (const)
		Real x() const { return _data[1]; }
		/// @brief Get y component (const)
		Real y() const { return _data[2]; }
		/// @brief Get z component (const)
		Real z() const { return _data[3]; }

		/// @brief Get scalar component w (mutable)
		Real& w() { return _data[0]; }
		/// @brief Get x component (mutable)
		Real& x() { return _data[1]; }
		/// @brief Get y component (mutable)
		Real& y() { return _data[2]; }
		/// @brief Get z component (mutable)
		Real& z() { return _data[3]; }

		/// @brief Access by index [0]=w, [1]=x, [2]=y, [3]=z (const)
		Real operator[](int i) const { return _data[i]; }
		/// @brief Access by index [0]=w, [1]=x, [2]=y, [3]=z (mutable)
		Real& operator[](int i) { return _data[i]; }

		/// @brief Get scalar part
		Real Scalar() const { return _data[0]; }

		/// @brief Get vector part as Vec3Cart
		Vec3Cart Vector() const 
		{ 
			return Vec3Cart(_data[1], _data[2], _data[3]); 
		}

		/// @brief Quaternion addition
		Quaternion operator+(const Quaternion& q) const
		{
			return Quaternion(_data[0] + q._data[0],
											 _data[1] + q._data[1],
											 _data[2] + q._data[2],
											 _data[3] + q._data[3]);
		}

		// Subtraction
		Quaternion operator-(const Quaternion& q) const
		{
			return Quaternion(_data[0] - q._data[0],
											 _data[1] - q._data[1],
											 _data[2] - q._data[2],
											 _data[3] - q._data[3]);
		}

		// Negation
		Quaternion operator-() const
		{
			return Quaternion(-_data[0], -_data[1], -_data[2], -_data[3]);
		}

		// Scalar multiplication
		Quaternion operator*(Real scalar) const
		{
			return Quaternion(_data[0] * scalar,
											 _data[1] * scalar,
											 _data[2] * scalar,
											 _data[3] * scalar);
		}

		// Scalar division
		Quaternion operator/(Real scalar) const
		{
			Real inv = 1.0 / scalar;
			return (*this) * inv;
		}

		/// @brief Quaternion multiplication (Hamilton product)
		/// @note Non-commutative: p*q != q*p in general
		Quaternion operator*(const Quaternion& q) const
		{
			return Quaternion(
				_data[0] * q._data[0] - _data[1] * q._data[1] - _data[2] * q._data[2] - _data[3] * q._data[3],
				_data[0] * q._data[1] + _data[1] * q._data[0] + _data[2] * q._data[3] - _data[3] * q._data[2],
				_data[0] * q._data[2] - _data[1] * q._data[3] + _data[2] * q._data[0] + _data[3] * q._data[1],
				_data[0] * q._data[3] + _data[1] * q._data[2] - _data[2] * q._data[1] + _data[3] * q._data[0]
			);
		}

		/// @brief In-place addition
		Quaternion& operator+=(const Quaternion& q)
		{
			_data[0] += q._data[0];
			_data[1] += q._data[1];
			_data[2] += q._data[2];
			_data[3] += q._data[3];
			return *this;
		}

		/// @brief In-place subtraction
		Quaternion& operator-=(const Quaternion& q)
		{
			_data[0] -= q._data[0];
			_data[1] -= q._data[1];
			_data[2] -= q._data[2];
			_data[3] -= q._data[3];
			return *this;
		}

		/// @brief In-place scalar multiplication
		Quaternion& operator*=(Real scalar)
		{
			_data[0] *= scalar;
			_data[1] *= scalar;
			_data[2] *= scalar;
			_data[3] *= scalar;
			return *this;
		}

		/// @brief In-place quaternion multiplication
		Quaternion& operator*=(const Quaternion& q)
		{
			*this = (*this) * q;
			return *this;
		}

		/// @brief Conjugate q* = w - xi - yj - zk
		/// @note For unit quaternions: q* = q^(-1)
		Quaternion Conjugate() const
		{
			return Quaternion(_data[0], -_data[1], -_data[2], -_data[3]);
		}

		/// @brief Squared norm ||q||² = w² + x² + y² + z²
		Real NormSquared() const
		{
			return _data[0] * _data[0] + 
						 _data[1] * _data[1] + 
						 _data[2] * _data[2] + 
						 _data[3] * _data[3];
		}

		/// @brief Norm (magnitude) ||q|| = sqrt(w² + x² + y² + z²)
		Real Norm() const
		{
			return std::sqrt(NormSquared());
		}

		/// @brief Inverse q^(-1) = q* / ||q||²
		/// @note For unit quaternions: q^(-1) = q*
		Quaternion Inverse() const
		{
			Real normSq = NormSquared();
			if (normSq < PrecisionValues<Real>::QuaternionZeroThreshold)
				throw QuaternionError("Cannot invert zero quaternion");
			
			Real invNormSq = 1.0 / normSq;
			return Quaternion(_data[0] * invNormSq,
											 -_data[1] * invNormSq,
											 -_data[2] * invNormSq,
											 -_data[3] * invNormSq);
		}

		/// @brief Normalize quaternion to unit length
		void Normalize()
		{
			Real norm = Norm();
			if (norm < PrecisionValues<Real>::QuaternionZeroThreshold)
				throw QuaternionError("Cannot normalize zero quaternion");
			
			Real invNorm = 1.0 / norm;
			_data[0] *= invNorm;
			_data[1] *= invNorm;
			_data[2] *= invNorm;
			_data[3] *= invNorm;
		}

		/// @brief Return normalized copy
		Quaternion Normalized() const
		{
			Quaternion result(*this);
			result.Normalize();
			return result;
		}

		/// @brief Check if quaternion is unit (within tolerance)
		bool isUnit(Real tolerance = PrecisionValues<Real>::DefaultTolerance) const
		{
			return std::abs(NormSquared() - 1.0) < tolerance;
		}

		/// @brief Check if quaternion is identity
		bool isIdentity(Real tolerance = PrecisionValues<Real>::DefaultTolerance) const
		{
			return std::abs(_data[0] - 1.0) < tolerance &&
						 std::abs(_data[1]) < tolerance &&
						 std::abs(_data[2]) < tolerance &&
						 std::abs(_data[3]) < tolerance;
		}

		/// @brief Dot product q1 · q2 = w1*w2 + x1*x2 + y1*y2 + z1*z2
		Real Dot(const Quaternion& q) const
		{
			return _data[0] * q._data[0] + 
						 _data[1] * q._data[1] + 
						 _data[2] * q._data[2] + 
						 _data[3] * q._data[3];
		}

		/// @brief Rotate a 3D vector using this quaternion
		/// @param v The vector to rotate
		/// @return Rotated vector
		/// @note Assumes this is a unit rotation quaternion. Formula: v' = q * [0, v] * q^(-1)
		Vec3Cart Rotate(const Vec3Cart& v) const
		{
			// Optimized version avoiding quaternion construction
			// v' = v + 2w(u × v) + 2u × (u × v)
			// where q = [w, u]
			
			Vec3Cart u(_data[1], _data[2], _data[3]);
			Real w = _data[0];
			
			Vec3Cart uCrossV = VectorProduct(u, v);
			Vec3Cart uCrossUCrossV = VectorProduct(u, uCrossV);
			
			return v + uCrossV * (2.0 * w) + uCrossUCrossV * 2.0;
		}

		/// @brief Get rotation axis (for non-identity rotation quaternions)
		/// @return Unit vector along rotation axis
		Vec3Cart GetRotationAxis() const
		{
			if (isIdentity())
				return Vec3Cart(0, 0, 1);
			
			Vec3Cart axis(_data[1], _data[2], _data[3]);
			Real vecNorm = axis.NormL2();
			
			if (vecNorm < PrecisionValues<Real>::QuaternionZeroThreshold)
				return Vec3Cart(0, 0, 1);
			
			return axis / vecNorm;
		}

		/// @brief Get rotation angle in radians
		/// @return Rotation angle θ = 2*acos(w)
		Real GetRotationAngle() const
		{
			Real w = _data[0];
			
			// Clamp to [-1, 1] for numerical stability
			if (w > 1.0) w = 1.0;
			if (w < -1.0) w = -1.0;
			
			return 2.0 * std::acos(w);
		}

		/// @brief Get axis-angle representation
		/// @param axis Output rotation axis
		/// @param angle Output rotation angle in radians
		void ToAxisAngle(Vec3Cart& axis, Real& angle) const
		{
			angle = GetRotationAngle();
			axis = GetRotationAxis();
		}

		/// @brief Convert to Euler angles (ZYX convention)
		/// @return Vector [yaw, pitch, roll] in radians
		Vec3Cart ToEulerZYX() const
		{
			Real w = _data[0], x = _data[1], y = _data[2], z = _data[3];
			
			// Roll (x-axis rotation)
			Real sinr_cosp = 2.0 * (w * x + y * z);
			Real cosr_cosp = 1.0 - 2.0 * (x * x + y * y);
			Real roll = std::atan2(sinr_cosp, cosr_cosp);
			
			// Pitch (y-axis rotation)
			Real sinp = 2.0 * (w * y - z * x);
			Real pitch;
			if (std::abs(sinp) >= 1.0)
				pitch = std::copysign(Constants::PI / 2.0, sinp); // Gimbal lock
			else
				pitch = std::asin(sinp);
			
			// Yaw (z-axis rotation)
			Real siny_cosp = 2.0 * (w * z + x * y);
			Real cosy_cosp = 1.0 - 2.0 * (y * y + z * z);
			Real yaw = std::atan2(siny_cosp, cosy_cosp);
			
			return Vec3Cart(yaw, pitch, roll);
		}

		/// @brief Convert to Euler angles (XYZ convention)
		/// @return Vector [roll, pitch, yaw] in radians
		/// @details Extracts angles for R = Rx(roll) * Ry(pitch) * Rz(yaw)
		Vec3Cart ToEulerXYZ() const
		{
			Real w = _data[0], x = _data[1], y = _data[2], z = _data[3];

			// Pitch (y-axis rotation)
			Real sinp = 2.0 * (x * z + w * y);
			Real pitch;
			if (std::abs(sinp) >= 1.0)
				pitch = std::copysign(Constants::PI / 2.0, sinp); // Gimbal lock
			else
				pitch = std::asin(sinp);

			// Roll (x-axis rotation)
			Real sinr = -2.0 * (y * z - w * x);
			Real cosr = 1.0 - 2.0 * (x * x + y * y);
			Real roll = std::atan2(sinr, cosr);

			// Yaw (z-axis rotation)
			Real siny = -2.0 * (x * y - w * z);
			Real cosy = 1.0 - 2.0 * (y * y + z * z);
			Real yaw = std::atan2(siny, cosy);

			return Vec3Cart(roll, pitch, yaw);
		}

		/// @brief Convert to Euler angles (ZXZ proper Euler convention)
		/// @return Vector [alpha, beta, gamma] in radians
		/// @details Extracts angles for R = Rz(alpha) * Rx(beta) * Rz(gamma).
		///          beta in [0, pi]. Gimbal lock when beta ≈ 0 or beta ≈ pi.
		Vec3Cart ToEulerZXZ() const
		{
			Real w = _data[0], x = _data[1], y = _data[2], z = _data[3];

			// beta from R22 = cos(beta) = 1 - 2(x² + y²)
			Real cosBeta = 1.0 - 2.0 * (x*x + y*y);
			cosBeta = std::clamp(cosBeta, static_cast<Real>(-1.0), static_cast<Real>(1.0));
			Real beta = std::acos(cosBeta);

			Real alpha, gamma;
			if (std::abs(std::sin(beta)) < 1e-10) {
				// Gimbal lock: beta ≈ 0 or π, assign all rotation to alpha
				alpha = std::atan2(2.0 * (x*y + w*z), 1.0 - 2.0 * (y*y + z*z));
				gamma = 0.0;
			} else {
				// alpha = atan2(R02, -R12) = atan2(2(xz+wy), 2(wx-yz))
				alpha = std::atan2(2.0 * (x*z + w*y), 2.0 * (w*x - y*z));
				// gamma = atan2(R20, R21) = atan2(2(xz-wy), 2(yz+wx))
				gamma = std::atan2(2.0 * (x*z - w*y), 2.0 * (y*z + w*x));
			}

			return Vec3Cart(alpha, beta, gamma);
		}

		/// @brief Convert to Euler angles (ZYZ proper Euler convention)
		/// @return Vector [alpha, beta, gamma] in radians
		/// @details Extracts angles for R = Rz(alpha) * Ry(beta) * Rz(gamma).
		///          beta in [0, pi]. Gimbal lock when beta ≈ 0 or beta ≈ pi.
		Vec3Cart ToEulerZYZ() const
		{
			Real w = _data[0], x = _data[1], y = _data[2], z = _data[3];

			// beta from R22 = cos(beta) = 1 - 2(x² + y²)
			Real cosBeta = 1.0 - 2.0 * (x*x + y*y);
			cosBeta = std::clamp(cosBeta, static_cast<Real>(-1.0), static_cast<Real>(1.0));
			Real beta = std::acos(cosBeta);

			Real alpha, gamma;
			if (std::abs(std::sin(beta)) < 1e-10) {
				// Gimbal lock: beta ≈ 0 or π, assign all rotation to alpha
				alpha = std::atan2(2.0 * (x*y + w*z), 1.0 - 2.0 * (y*y + z*z));
				gamma = 0.0;
			} else {
				// alpha = atan2(R12, R02) = atan2(2(yz-wx), 2(xz+wy))
				alpha = std::atan2(2.0 * (y*z - w*x), 2.0 * (x*z + w*y));
				// gamma = atan2(R21, -R20) = atan2(2(yz+wx), 2(wy-xz))
				gamma = std::atan2(2.0 * (y*z + w*x), 2.0 * (w*y - x*z));
			}

			return Vec3Cart(alpha, beta, gamma);
		}

		/// @brief Convert to 3×3 rotation matrix
		/// @return Orthogonal matrix R such that R*v rotates vector v
		MatrixNM<Real, 3, 3> ToRotationMatrix() const
		{
			Real w = _data[0], x = _data[1], y = _data[2], z = _data[3];
			
			MatrixNM<Real, 3, 3> mat;
			
			// First row
			mat[0][0] = 1.0 - 2.0 * (y*y + z*z);
			mat[0][1] = 2.0 * (x*y - w*z);
			mat[0][2] = 2.0 * (x*z + w*y);
			
			// Second row
			mat[1][0] = 2.0 * (x*y + w*z);
			mat[1][1] = 1.0 - 2.0 * (x*x + z*z);
			mat[1][2] = 2.0 * (y*z - w*x);
			
			// Third row
			mat[2][0] = 2.0 * (x*z - w*y);
			mat[2][1] = 2.0 * (y*z + w*x);
			mat[2][2] = 1.0 - 2.0 * (x*x + y*y);
			
			return mat;
		}

		/// @brief Create quaternion from 3×3 rotation matrix
		/// @param mat Orthogonal rotation matrix
		/// @return Quaternion representation
		/// @details Uses Shepperd's method for numerical stability
		static Quaternion FromRotationMatrix(const MatrixNM<Real, 3, 3>& mat)
		{
			Real trace = mat[0][0] + mat[1][1] + mat[2][2];
			
			if (trace > 0.0)
			{
				// w is the largest component
				Real s = std::sqrt(trace + 1.0) * 2.0;  // s = 4*w
				Real w = 0.25 * s;
				Real x = (mat[2][1] - mat[1][2]) / s;
				Real y = (mat[0][2] - mat[2][0]) / s;
				Real z = (mat[1][0] - mat[0][1]) / s;
				return Quaternion(w, x, y, z);
			}
			else if (mat[0][0] > mat[1][1] && mat[0][0] > mat[2][2])
			{
				// x is the largest component
				Real s = std::sqrt(1.0 + mat[0][0] - mat[1][1] - mat[2][2]) * 2.0;  // s = 4*x
				Real w = (mat[2][1] - mat[1][2]) / s;
				Real x = 0.25 * s;
				Real y = (mat[0][1] + mat[1][0]) / s;
				Real z = (mat[0][2] + mat[2][0]) / s;
				return Quaternion(w, x, y, z);
			}
			else if (mat[1][1] > mat[2][2])
			{
				// y is the largest component
				Real s = std::sqrt(1.0 + mat[1][1] - mat[0][0] - mat[2][2]) * 2.0;  // s = 4*y
				Real w = (mat[0][2] - mat[2][0]) / s;
				Real x = (mat[0][1] + mat[1][0]) / s;
				Real y = 0.25 * s;
				Real z = (mat[1][2] + mat[2][1]) / s;
				return Quaternion(w, x, y, z);
			}
			else
			{
				// z is the largest component
				Real s = std::sqrt(1.0 + mat[2][2] - mat[0][0] - mat[1][1]) * 2.0;  // s = 4*z
				Real w = (mat[1][0] - mat[0][1]) / s;
				Real x = (mat[0][2] + mat[2][0]) / s;
				Real y = (mat[1][2] + mat[2][1]) / s;
				Real z = 0.25 * s;
				return Quaternion(w, x, y, z);
			}
		}

		/// @brief Linear interpolation (faster but not constant angular velocity)
		/// @param q1 Start quaternion
		/// @param q2 End quaternion
		/// @param t Interpolation parameter [0,1]
		/// @return Interpolated quaternion (should be normalized for rotations)
		static Quaternion Lerp(const Quaternion& q1, const Quaternion& q2, Real t)
		{
			return q1 * (1.0 - t) + q2 * t;
		}

		/// @brief Spherical linear interpolation (constant angular velocity)
		/// @param q1 Start quaternion
		/// @param q2 End quaternion
		/// @param t Interpolation parameter [0,1]
		/// @return Interpolated quaternion along shortest arc on 4D unit sphere
		static Quaternion Slerp(const Quaternion& q1, const Quaternion& q2, Real t)
		{
			Quaternion q2_adjusted = q2;
			
			// Compute dot product
			Real dot = q1.Dot(q2);
			
			// If dot < 0, negate q2 to take shortest path
			if (dot < 0.0)
			{
				q2_adjusted = -q2_adjusted;
				dot = -dot;
			}
			
			// If quaternions are very close, use linear interpolation
			if (dot > 0.9995)
			{
				return Lerp(q1, q2_adjusted, t).Normalized();
			}
			
			// Perform slerp
			Real theta = std::acos(dot);
			Real sinTheta = std::sin(theta);

			// Guard against near-zero sinTheta for numerical stability
			if (std::abs(sinTheta) < Real{1e-6})
			{
				return Lerp(q1, q2_adjusted, t).Normalized();
			}
			
			Real w1 = std::sin((1.0 - t) * theta) / sinTheta;
			Real w2 = std::sin(t * theta) / sinTheta;
			
			return q1 * w1 + q2_adjusted * w2;
		}

		/// @brief Exact equality test
		bool operator==(const Quaternion& q) const
		{
			return _data[0] == q._data[0] &&
						 _data[1] == q._data[1] &&
						 _data[2] == q._data[2] &&
						 _data[3] == q._data[3];
		}

		/// @brief Inequality test
		bool operator!=(const Quaternion& q) const
		{
			return !(*this == q);
		}

		/// @brief Check approximate equality within tolerance
		bool isApprox(const Quaternion& q, Real tolerance = PrecisionValues<Real>::DefaultTolerance) const
		{
			return std::abs(_data[0] - q._data[0]) < tolerance &&
						 std::abs(_data[1] - q._data[1]) < tolerance &&
						 std::abs(_data[2] - q._data[2]) < tolerance &&
						 std::abs(_data[3] - q._data[3]) < tolerance;
		}

		/// @brief Stream output operator
		friend std::ostream& operator<<(std::ostream& os, const Quaternion& q)
		{
			os << "[" << q._data[0] << ", " 
				 << q._data[1] << "i, " 
				 << q._data[2] << "j, " 
				 << q._data[3] << "k]";
			return os;
		}

		/// @brief Print to output stream
		void Print(std::ostream& os = std::cout) const
		{
			os << *this;
		}
	};

	/// @brief Scalar * Quaternion operator
	inline Quaternion operator*(Real scalar, const Quaternion& q)
	{
		return q * scalar;
	}

} // namespace MML

#endif // MML_QUATERNIONS_H
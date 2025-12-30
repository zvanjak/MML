///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Quaternions.h                                                       ///
///  Description: Quaternion class for 3D rotations and orientation                   ///
///               SLERP interpolation, conversion to/from Euler angles                ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_QUATERNIONS_H
#define MML_QUATERNIONS_H

#include "MMLBase.h"
#include "base/Vector.h"
#include "base/VectorTypes.h"
#include "base/MatrixNM.h"

namespace MML
{
	/*****************************************************************************
	 * QUATERNION CONVENTIONS AND MATHEMATICAL BACKGROUND
	 * ====================================================
	 * 
	 * REPRESENTATION:
	 * A quaternion is represented as: q = w + xi + yj + zk
	 * where:
	 *   - w is the scalar (real) part
	 *   - (x, y, z) is the vector (imaginary) part
	 *   - i, j, k are the imaginary units satisfying:
	 *     i² = j² = k² = ijk = -1
	 *     ij = k,  jk = i,  ki = j
	 *     ji = -k, kj = -i, ik = -j
	 * 
	 * STORAGE ORDER:
	 * We store quaternions as [w, x, y, z] where:
	 *   _data[0] = w  (scalar/real part)
	 *   _data[1] = x  (i component)
	 *   _data[2] = y  (j component)
	 *   _data[3] = z  (k component)
	 * 
	 * MULTIPLICATION CONVENTION (Hamilton convention):
	 * For quaternions p = (w1, x1, y1, z1) and q = (w2, x2, y2, z2):
	 * 
	 * p * q = (w1*w2 - x1*x2 - y1*y2 - z1*z2,
	 *          w1*x2 + x1*w2 + y1*z2 - z1*y2,
	 *          w1*y2 - x1*z2 + y1*w2 + z1*x2,
	 *          w1*z2 + x1*y2 - y1*x2 + z1*w2)
	 * 
	 * ROTATION CONVENTION:
	 * To rotate a vector v by angle θ around unit axis u:
	 *   1. Create rotation quaternion: q = [cos(θ/2), sin(θ/2)*u]
	 *   2. Extend vector to quaternion: v_quat = [0, v]
	 *   3. Apply rotation: v' = q * v_quat * q^(-1)
	 * 
	 * This follows RIGHT-HANDED rotation (counter-clockwise when looking
	 * along the axis toward the origin).
	 * 
	 * NORMALIZATION:
	 * Unit quaternions (||q|| = 1) represent rotations. We provide:
	 *   - IsUnit(): Check if quaternion is normalized
	 *   - Normalize(): Normalize to unit length
	 *   - Normalized(): Return normalized copy
	 * 
	 * AXIS-ANGLE CONVENTION:
	 * When constructing from axis-angle:
	 *   - Axis must be a UNIT vector
	 *   - Angle is in RADIANS
	 *   - Positive angle = right-hand rule rotation
	 * 
	 * EULER ANGLES CONVENTION:
	 * We support multiple Euler angle conventions:
	 *   - ZYX (Yaw-Pitch-Roll): Common in aerospace
	 *   - XYZ: Alternative rotation order
	 *   - Angles are in RADIANS
	 *   - Applied in INTRINSIC order (each rotation in rotating frame)
	 * 
	 * INTERPOLATION:
	 * - Slerp (Spherical Linear Interpolation): Shortest path on 4D unit sphere
	 *   Used for smooth rotation interpolation between orientations
	 * - Lerp (Linear Interpolation): Faster but not constant angular velocity
	 *   Must be followed by normalization for rotations
	 * 
	 *****************************************************************************/

	class Quaternion
	{
	private:
		Real _data[4];  // [w, x, y, z] storage

	public:
		/********************************************************************
		 * CONSTRUCTORS
		 ********************************************************************/

		// Default: Identity quaternion [1, 0, 0, 0] (no rotation)
		Quaternion() : _data{1, 0, 0, 0} {}

		// From components: q = w + xi + yj + zk
		Quaternion(Real w, Real x, Real y, Real z) : _data{w, x, y, z} {}

		// From scalar and vector parts: q = w + v
		Quaternion(Real w, const Vec3Cart& vec) 
			: _data{w, vec[0], vec[1], vec[2]} {}

		// Pure imaginary quaternion from vector: q = 0 + v
		explicit Quaternion(const Vec3Cart& vec)
			: _data{0, vec[0], vec[1], vec[2]} {}

		// Copy constructor
		Quaternion(const Quaternion& q) 
			: _data{q._data[0], q._data[1], q._data[2], q._data[3]} {}

		/********************************************************************
		 * STATIC FACTORY METHODS FOR ROTATIONS
		 ********************************************************************/

		// Create identity quaternion (no rotation)
		static Quaternion Identity()
		{
			return Quaternion(1, 0, 0, 0);
		}

		// Create rotation quaternion from axis-angle
		// axis: UNIT vector defining rotation axis
		// angle: rotation angle in RADIANS (right-hand rule)
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

		// Create rotation quaternion from Euler angles (ZYX convention - Yaw, Pitch, Roll)
		// yaw: rotation around Z-axis (radians)
		// pitch: rotation around Y-axis (radians)
		// roll: rotation around X-axis (radians)
		// Applied as: R = Rz(yaw) * Ry(pitch) * Rx(roll)
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

		// Create rotation quaternion from Euler angles (XYZ convention)
		// roll: rotation around X-axis (radians)
		// pitch: rotation around Y-axis (radians)
		// yaw: rotation around Z-axis (radians)
		// Applied as: R = Rx(roll) * Ry(pitch) * Rz(yaw)
		static Quaternion FromEulerXYZ(Real roll, Real pitch, Real yaw)
		{
			Real cr = std::cos(roll * 0.5);
			Real sr = std::sin(roll * 0.5);
			Real cp = std::cos(pitch * 0.5);
			Real sp = std::sin(pitch * 0.5);
			Real cy = std::cos(yaw * 0.5);
			Real sy = std::sin(yaw * 0.5);

			return Quaternion(
				cr * cp * cy + sr * sp * sy,
				sr * cp * cy - cr * sp * sy,
				cr * sp * cy + sr * cp * sy,
				cr * cp * sy - sr * sp * cy
			);
		}

		/********************************************************************
		 * ACCESSORS
		 ********************************************************************/

		Real w() const { return _data[0]; }
		Real x() const { return _data[1]; }
		Real y() const { return _data[2]; }
		Real z() const { return _data[3]; }

		Real& w() { return _data[0]; }
		Real& x() { return _data[1]; }
		Real& y() { return _data[2]; }
		Real& z() { return _data[3]; }

		// Access by index: [0]=w, [1]=x, [2]=y, [3]=z
		Real operator[](int i) const { return _data[i]; }
		Real& operator[](int i) { return _data[i]; }

		// Get scalar part
		Real Scalar() const { return _data[0]; }

		// Get vector part as Vec3Cart
		Vec3Cart Vector() const 
		{ 
			return Vec3Cart(_data[1], _data[2], _data[3]); 
		}

		/********************************************************************
		 * BASIC OPERATIONS
		 ********************************************************************/

		// Addition
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

		// Quaternion multiplication (Hamilton product)
		// NOTE: Non-commutative! p*q != q*p in general
		Quaternion operator*(const Quaternion& q) const
		{
			return Quaternion(
				_data[0] * q._data[0] - _data[1] * q._data[1] - _data[2] * q._data[2] - _data[3] * q._data[3],
				_data[0] * q._data[1] + _data[1] * q._data[0] + _data[2] * q._data[3] - _data[3] * q._data[2],
				_data[0] * q._data[2] - _data[1] * q._data[3] + _data[2] * q._data[0] + _data[3] * q._data[1],
				_data[0] * q._data[3] + _data[1] * q._data[2] - _data[2] * q._data[1] + _data[3] * q._data[0]
			);
		}

		// In-place operations
		Quaternion& operator+=(const Quaternion& q)
		{
			_data[0] += q._data[0];
			_data[1] += q._data[1];
			_data[2] += q._data[2];
			_data[3] += q._data[3];
			return *this;
		}

		Quaternion& operator-=(const Quaternion& q)
		{
			_data[0] -= q._data[0];
			_data[1] -= q._data[1];
			_data[2] -= q._data[2];
			_data[3] -= q._data[3];
			return *this;
		}

		Quaternion& operator*=(Real scalar)
		{
			_data[0] *= scalar;
			_data[1] *= scalar;
			_data[2] *= scalar;
			_data[3] *= scalar;
			return *this;
		}

		Quaternion& operator*=(const Quaternion& q)
		{
			*this = (*this) * q;
			return *this;
		}

		/********************************************************************
		 * QUATERNION-SPECIFIC OPERATIONS
		 ********************************************************************/

		// Conjugate: q* = w - xi - yj - zk
		// For unit quaternions: q* = q^(-1)
		Quaternion Conjugate() const
		{
			return Quaternion(_data[0], -_data[1], -_data[2], -_data[3]);
		}

		// Squared norm: ||q||² = w² + x² + y² + z²
		Real NormSquared() const
		{
			return _data[0] * _data[0] + 
						 _data[1] * _data[1] + 
						 _data[2] * _data[2] + 
						 _data[3] * _data[3];
		}

		// Norm (magnitude): ||q|| = sqrt(w² + x² + y² + z²)
		Real Norm() const
		{
			return std::sqrt(NormSquared());
		}

		// Inverse: q^(-1) = q* / ||q||²
		// For unit quaternions: q^(-1) = q*
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

		// Normalize quaternion to unit length
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

		// Return normalized copy
		Quaternion Normalized() const
		{
			Quaternion result(*this);
			result.Normalize();
			return result;
		}

		// Check if quaternion is unit (within tolerance)
		bool IsUnit(Real tolerance = PrecisionValues<Real>::DefaultTolerance) const
		{
			return std::abs(NormSquared() - 1.0) < tolerance;
		}

		// Check if quaternion is identity
		bool IsIdentity(Real tolerance = PrecisionValues<Real>::DefaultTolerance) const
		{
			return std::abs(_data[0] - 1.0) < tolerance &&
						 std::abs(_data[1]) < tolerance &&
						 std::abs(_data[2]) < tolerance &&
						 std::abs(_data[3]) < tolerance;
		}

		// Dot product: q1 · q2 = w1*w2 + x1*x2 + y1*y2 + z1*z2
		Real Dot(const Quaternion& q) const
		{
			return _data[0] * q._data[0] + 
						 _data[1] * q._data[1] + 
						 _data[2] * q._data[2] + 
						 _data[3] * q._data[3];
		}

		/********************************************************************
		 * ROTATION OPERATIONS
		 ********************************************************************/

		// Rotate a 3D vector using this quaternion
		// Assumes this quaternion is a unit rotation quaternion
		// Formula: v' = q * [0, v] * q^(-1)
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

		// Get rotation axis (for non-identity rotation quaternions)
		// Returns unit vector along rotation axis
		Vec3Cart GetRotationAxis() const
		{
			if (IsIdentity())
				return Vec3Cart(0, 0, 1);  // Arbitrary axis for zero rotation
			
			Vec3Cart axis(_data[1], _data[2], _data[3]);
			Real vecNorm = axis.NormL2();
			
			if (vecNorm < PrecisionValues<Real>::QuaternionZeroThreshold)
				return Vec3Cart(0, 0, 1);  // Arbitrary axis
			
			return axis / vecNorm;
		}

		// Get rotation angle (in radians)
		Real GetRotationAngle() const
		{
			// For unit quaternion q = [cos(θ/2), sin(θ/2)*axis]
			// angle θ = 2 * acos(w)
			Real w = _data[0];
			
			// Clamp to [-1, 1] for numerical stability
			if (w > 1.0) w = 1.0;
			if (w < -1.0) w = -1.0;
			
			return 2.0 * std::acos(w);
		}

		// Get axis-angle representation
		void ToAxisAngle(Vec3Cart& axis, Real& angle) const
		{
			angle = GetRotationAngle();
			axis = GetRotationAxis();
		}

		// Convert to Euler angles (ZYX convention: Yaw, Pitch, Roll)
		// Returns [yaw, pitch, roll] in radians
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

		// Convert to 3×3 rotation matrix
		// Returns orthogonal matrix R such that R*v rotates vector v
		// Compatible with CoordTransf3D transformation matrices
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

		// Create quaternion from 3×3 rotation matrix
		// Matrix must be orthogonal (rotation matrix)
		// Uses Shepperd's method for numerical stability
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

		/********************************************************************
		 * INTERPOLATION
		 ********************************************************************/

		// Linear interpolation (faster but not constant angular velocity)
		// t ∈ [0, 1]: 0 returns this, 1 returns q
		// Result should be normalized for rotations
		static Quaternion Lerp(const Quaternion& q1, const Quaternion& q2, Real t)
		{
			return q1 * (1.0 - t) + q2 * t;
		}

		// Spherical linear interpolation (constant angular velocity)
		// t ∈ [0, 1]: 0 returns this, 1 returns q
		// Interpolates along shortest arc on 4D unit sphere
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
			
			Real w1 = std::sin((1.0 - t) * theta) / sinTheta;
			Real w2 = std::sin(t * theta) / sinTheta;
			
			return q1 * w1 + q2_adjusted * w2;
		}

		/********************************************************************
		 * COMPARISON
		 ********************************************************************/

		bool operator==(const Quaternion& q) const
		{
			return _data[0] == q._data[0] &&
						 _data[1] == q._data[1] &&
						 _data[2] == q._data[2] &&
						 _data[3] == q._data[3];
		}

		bool operator!=(const Quaternion& q) const
		{
			return !(*this == q);
		}

		// Check approximate equality
		bool IsApprox(const Quaternion& q, Real tolerance = PrecisionValues<Real>::DefaultTolerance) const
		{
			return std::abs(_data[0] - q._data[0]) < tolerance &&
						 std::abs(_data[1] - q._data[1]) < tolerance &&
						 std::abs(_data[2] - q._data[2]) < tolerance &&
						 std::abs(_data[3] - q._data[3]) < tolerance;
		}

		/********************************************************************
		 * OUTPUT
		 ********************************************************************/

		friend std::ostream& operator<<(std::ostream& os, const Quaternion& q)
		{
			os << "[" << q._data[0] << ", " 
				 << q._data[1] << "i, " 
				 << q._data[2] << "j, " 
				 << q._data[3] << "k]";
			return os;
		}

		void Print(std::ostream& os = std::cout) const
		{
			os << *this;
		}
	};

	// Scalar * Quaternion
	inline Quaternion operator*(Real scalar, const Quaternion& q)
	{
		return q * scalar;
	}

} // namespace MML

#endif // MML_QUATERNIONS_H
///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        CoordTransf3D.h                                                     ///
///  Description: 3D coordinate transformations (rotation, scaling, translation)      ///
///               Euler angles, quaternion rotations, rigid body transformations      ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_COORD_TRANSF_3D_H
#define MML_COORD_TRANSF_3D_H

#include "MMLBase.h"

#include "base/Quaternions.h"
#include "core/CoordTransf.h"
#include "core/LinAlgEqSolvers.h"

namespace MML
{
	/// @brief 3D Cartesian rotation around X-axis
	/// 
	/// Rotates 3D Cartesian coordinates around the X-axis by a fixed angle.
	/// Rotation matrix: [1 0 0; 0 cos(?) -sin(?); 0 sin(?) cos(?)]
	/// X-coordinate remains unchanged, Y and Z coordinates rotate in YZ-plane.
	class CoordTransfCart3DRotationXAxis : 
		public CoordTransfWithInverse<Vec3Cart, Vec3Cart, 3>
	{
	private:
		Real    _angle;
		MatrixNM<Real, 3, 3>  _transf;
		MatrixNM<Real, 3, 3>  _inverse;

		const ScalarFunctionFromStdFunc<3> _f1, _f2, _f3;
		const ScalarFunctionFromStdFunc<3> _fInverse1, _fInverse2, _fInverse3;

	public:
		/// @brief Constructor for X-axis rotation
		/// @param inAngle Rotation angle in radians (counterclockwise looking down positive X-axis)
		CoordTransfCart3DRotationXAxis(Real inAngle) : _angle(inAngle),
			_f1([this](const VectorN<Real, 3>& q) { return func1(q); }),
			_f2([this](const VectorN<Real, 3>& q) { return func2(q); }),
			_f3([this](const VectorN<Real, 3>& q) { return func3(q); }),
			_fInverse1([this](const VectorN<Real, 3>& q) { return funcInverse1(q); }),
			_fInverse2([this](const VectorN<Real, 3>& q) { return funcInverse2(q); }),
			_fInverse3([this](const VectorN<Real, 3>& q) { return funcInverse3(q); })
		{
			_transf[0][0] = 1.0;
			_transf[1][1] = cos(_angle);
			_transf[1][2] = -sin(_angle);
			_transf[2][1] = sin(_angle);
			_transf[2][2] = cos(_angle);

			_inverse[0][0] = 1.0;
			_inverse[1][1] = cos(_angle);
			_inverse[1][2] = sin(_angle);
			_inverse[2][1] = -sin(_angle);
			_inverse[2][2] = cos(_angle);
		}

		/// @brief X component unchanged by X-axis rotation
		/// @param q Input Cartesian coordinates
		/// @return x coordinate (unchanged)
		Real func1(const VectorN<Real, 3>& q) const { return q[0]; }
		
		/// @brief Rotated Y component: y' = y�cos(?) - z�sin(?)
		/// @param q Input Cartesian coordinates
		/// @return y component after X-axis rotation
		Real func2(const VectorN<Real, 3>& q) const { return q[1] * cos(_angle) - q[2] * sin(_angle); }
		
		/// @brief Rotated Z component: z' = y�sin(?) + z�cos(?)
		/// @param q Input Cartesian coordinates
		/// @return z component after X-axis rotation
		Real func3(const VectorN<Real, 3>& q) const { return q[1] * sin(_angle) + q[2] * cos(_angle); }

		Real funcInverse1(const VectorN<Real, 3>& q) const { return q[0]; }
		Real funcInverse2(const VectorN<Real, 3>& q) const { return q[1] * cos(_angle) + q[2] * sin(_angle); }
		Real funcInverse3(const VectorN<Real, 3>& q) const { return -q[1] * sin(_angle) + q[2] * cos(_angle); }

		Vec3Cart    transf(const Vec3Cart& q) const { return Vec3Cart{ func1(q), func2(q), func3(q) }; }
		Vec3Cart    transfInverse(const Vec3Cart& q) const { return Vec3Cart{ funcInverse1(q), funcInverse2(q), funcInverse3(q) }; }

		const IScalarFunction<3>& coordTransfFunc(int i) const
		{
			if (i == 0) return _f1;
			else if (i == 1) return _f2;
			else return _f3;
		}
		const IScalarFunction<3>& inverseCoordTransfFunc(int i)const
		{
			if (i == 0) return _fInverse1;
			else if (i == 1) return _fInverse2;
			else return _fInverse3;
		}
	};
	/// @brief 3D Cartesian rotation around Y-axis
	/// 
	/// Rotates 3D Cartesian coordinates around the Y-axis by a fixed angle.
	/// Rotation matrix: [cos(?) 0 sin(?); 0 1 0; -sin(?) 0 cos(?)]
	/// Y-coordinate remains unchanged, X and Z coordinates rotate in XZ-plane.
	class CoordTransfCart3DRotationYAxis : public CoordTransfWithInverse<Vec3Cart, Vec3Cart, 3>
	{
	private:
		Real    _angle;
		MatrixNM<Real, 3, 3>  _transf;
		MatrixNM<Real, 3, 3>  _inverse;

		const ScalarFunctionFromStdFunc<3> _f1, _f2, _f3;
		const ScalarFunctionFromStdFunc<3> _fInverse1, _fInverse2, _fInverse3;

	public:
		CoordTransfCart3DRotationYAxis(Real inAngle) : _angle(inAngle),
			_f1([this](const VectorN<Real, 3>& q) { return func1(q); }),
			_f2([this](const VectorN<Real, 3>& q) { return func2(q); }),
			_f3([this](const VectorN<Real, 3>& q) { return func3(q); }),
			_fInverse1([this](const VectorN<Real, 3>& q) { return funcInverse1(q); }),
			_fInverse2([this](const VectorN<Real, 3>& q) { return funcInverse2(q); }),
			_fInverse3([this](const VectorN<Real, 3>& q) { return funcInverse3(q); }) 
		{
			_transf[0][0] = cos(_angle);
			_transf[0][2] = sin(_angle);
			_transf[1][1] = 1.0;
			_transf[2][0] = -sin(_angle);
			_transf[2][2] = cos(_angle);

			_inverse[0][0] = cos(_angle);
			_inverse[0][2] = -sin(_angle);
			_inverse[1][1] = 1.0;
			_inverse[2][0] = sin(_angle);
			_inverse[2][2] = cos(_angle);
		}

		Real func1(const VectorN<Real, 3>& q) const { return q[0] * cos(_angle) + q[2] * sin(_angle); }
		Real func2(const VectorN<Real, 3>& q) const { return q[1]; }
		Real func3(const VectorN<Real, 3>& q) const { return -q[0] * sin(_angle) + q[2] * cos(_angle); }

		Real funcInverse1(const VectorN<Real, 3>& q) const { return q[0] * cos(_angle) - q[2] * sin(_angle); }
		Real funcInverse2(const VectorN<Real, 3>& q) const { return q[1]; }
		Real funcInverse3(const VectorN<Real, 3>& q) const { return q[0] * sin(_angle) + q[2] * cos(_angle); }

		Vec3Cart    transf(const Vec3Cart& q) const { return Vec3Cart{ func1(q), func2(q), func3(q) }; }
		Vec3Cart    transfInverse(const Vec3Cart& q) const { return Vec3Cart{ funcInverse1(q), funcInverse2(q), funcInverse3(q) }; }

		const IScalarFunction<3>& coordTransfFunc(int i) const
		{
			if (i == 0) return _f1;
			else if (i == 1) return _f2;
			else return _f3;
		}
		const IScalarFunction<3>& inverseCoordTransfFunc(int i) const
		{
			if (i == 0) return _fInverse1;
			else if (i == 1) return _fInverse2;
			else return _fInverse3;
		}
	};
	/// @brief 3D Cartesian rotation around Z-axis
	/// 
	/// Rotates 3D Cartesian coordinates around the Z-axis by a fixed angle.
	/// Rotation matrix: [cos(?) -sin(?) 0; sin(?) cos(?) 0; 0 0 1]
	/// Z-coordinate remains unchanged, X and Y coordinates rotate in XY-plane.
	class CoordTransfCart3DRotationZAxis : public CoordTransfWithInverse<Vec3Cart, Vec3Cart, 3>
	{
	private:
		Real    _angle;
		MatrixNM<Real, 3, 3>  _transf;
		MatrixNM<Real, 3, 3>  _inverse;

		const ScalarFunctionFromStdFunc<3> _f1, _f2, _f3;
		const ScalarFunctionFromStdFunc<3> _fInverse1, _fInverse2, _fInverse3;

	public:
		CoordTransfCart3DRotationZAxis(Real inAngle) : _angle(inAngle),
			_f1([this](const VectorN<Real, 3>& q) { return func1(q); }),
			_f2([this](const VectorN<Real, 3>& q) { return func2(q); }),
			_f3([this](const VectorN<Real, 3>& q) { return func3(q); }),
			_fInverse1([this](const VectorN<Real, 3>& q) { return funcInverse1(q); }),
			_fInverse2([this](const VectorN<Real, 3>& q) { return funcInverse2(q); }),
			_fInverse3([this](const VectorN<Real, 3>& q) { return funcInverse3(q); }) 
		{
			_transf[0][0] = cos(_angle);
			_transf[0][1] = -sin(_angle);
			_transf[1][0] = sin(_angle);
			_transf[1][1] = cos(_angle);
			_transf[2][2] = 1.0;

			_inverse[0][0] = cos(_angle);
			_inverse[0][1] = sin(_angle);
			_inverse[1][0] = -sin(_angle);
			_inverse[1][1] = cos(_angle);
			_inverse[2][2] = 1.0;
		}

		Real func1(const VectorN<Real, 3>& q) const { return q[0] * cos(_angle) - q[1] * sin(_angle); }
		Real func2(const VectorN<Real, 3>& q) const { return q[0] * sin(_angle) + q[1] * cos(_angle); }
		Real func3(const VectorN<Real, 3>& q) const { return q[2]; }

		Real funcInverse1(const VectorN<Real, 3>& q) const { return q[0] * cos(_angle) + q[1] * sin(_angle); }
		Real funcInverse2(const VectorN<Real, 3>& q) const { return -q[0] * sin(_angle) + q[1] * cos(_angle); }
		Real funcInverse3(const VectorN<Real, 3>& q) const { return q[2]; }

		Vec3Cart    transf(const Vec3Cart& q) const { return Vec3Cart{ func1(q), func2(q), func3(q) }; }
		Vec3Cart    transfInverse(const Vec3Cart& q) const { return Vec3Cart{ funcInverse1(q), funcInverse2(q), funcInverse3(q) }; }

		const IScalarFunction<3>& coordTransfFunc(int i) const
		{
			if (i == 0) return _f1;
			else if (i == 1) return _f2;
			else return _f3;
		}
		const IScalarFunction<3>& inverseCoordTransfFunc(int i) const
		{
			if (i == 0) return _fInverse1;
			else if (i == 1) return _fInverse2;
			else return _fInverse3;
		}
	};

	/// @brief General 3D rotation using Euler angles (ZXZ convention)
	/// 
	/// Implements 3D rotation using Euler angles with ZXZ convention (physics standard).
	/// Three successive rotations: a around Z, � around X', ? around Z''.
	/// This convention avoids gimbal lock in many physical applications.
	class CoordTransfCart3DRotationEuler : public CoordTransfWithInverse<Vec3Cart, Vec3Cart, 3>
	{
	private:
		Real _alpha; ///< First rotation angle around Z axis (radians)
		Real _beta;  ///< Second rotation angle around X' axis (radians)
		Real _gamma; ///< Third rotation angle around Z'' axis (radians)

		MatrixNM<Real, 3, 3> _transf;
		MatrixNM<Real, 3, 3> _inverse;

		const ScalarFunctionFromStdFunc<3> _f1, _f2, _f3;
		const ScalarFunctionFromStdFunc<3> _fInverse1, _fInverse2, _fInverse3;

		/// @brief Build ZXZ Euler angle rotation matrix
		/// @param alpha First angle (Z-axis rotation)
		/// @param beta Second angle (X'-axis rotation)
		/// @param gamma Third angle (Z''-axis rotation)
		/// @return 3x3 rotation matrix
		static MatrixNM<Real, 3, 3> buildEulerZXZMatrix(Real alpha, Real beta, Real gamma)
		{
			Real ca = std::cos(alpha), sa = std::sin(alpha);
			Real cb = std::cos(beta), sb = std::sin(beta);
			Real cg = std::cos(gamma), sg = std::sin(gamma);

			MatrixNM<Real, 3, 3> m;
			m(0, 0) = ca * cb * cg - sa * sg;
			m(0, 1) = -ca * cb * sg - sa * cg;
			m(0, 2) = ca * sb;
			m(1, 0) = sa * cb * cg + ca * sg;
			m(1, 1) = -sa * cb * sg + ca * cg;
			m(1, 2) = sa * sb;
			m(2, 0) = -sb * cg;
			m(2, 1) = sb * sg;
			m(2, 2) = cb;
			return m;
		}

	public:
		/// @brief Constructor using ZXZ Euler angles (physics convention)
		/// @param alpha First rotation around Z-axis (radians)
		/// @param beta Second rotation around X'-axis (radians)
		/// @param gamma Third rotation around Z''-axis (radians)
		CoordTransfCart3DRotationEuler(Real alpha, Real beta, Real gamma)
			: _alpha(alpha), _beta(beta), _gamma(gamma),
				_transf(buildEulerZXZMatrix(alpha, beta, gamma)),
				_inverse(buildEulerZXZMatrix(-gamma, -beta, -alpha)), // Inverse: reverse order and sign
				_f1([this](const VectorN<Real, 3>& q) { return func1(q); }),
				_f2([this](const VectorN<Real, 3>& q) { return func2(q); }),
				_f3([this](const VectorN<Real, 3>& q) { return func3(q); }),
				_fInverse1([this](const VectorN<Real, 3>& q) { return funcInverse1(q); }),
				_fInverse2([this](const VectorN<Real, 3>& q) { return funcInverse2(q); }),
				_fInverse3([this](const VectorN<Real, 3>& q) { return funcInverse3(q); })
		{	}

		Real func1(const VectorN<Real, 3>& q) const { return (_transf * q)[0]; }
		Real func2(const VectorN<Real, 3>& q) const { return (_transf * q)[1]; }
		Real func3(const VectorN<Real, 3>& q) const { return (_transf * q)[2]; }

		Real funcInverse1(const VectorN<Real, 3>& q) const { return (_inverse * q)[0]; }
		Real funcInverse2(const VectorN<Real, 3>& q) const { return (_inverse * q)[1]; }
		Real funcInverse3(const VectorN<Real, 3>& q) const { return (_inverse * q)[2]; }

		Vec3Cart transf(const Vec3Cart& q) const override
		{
			return Vec3Cart{ func1(q), func2(q), func3(q) };
		}
		Vec3Cart transfInverse(const Vec3Cart& q) const override
		{
			return Vec3Cart{ funcInverse1(q), funcInverse2(q), funcInverse3(q) };
		}

		const IScalarFunction<3>& coordTransfFunc(int i) const override
		{
			if (i == 0) return _f1;
			else if (i == 1) return _f2;
			else return _f3;
		}
		const IScalarFunction<3>& inverseCoordTransfFunc(int i) const override
		{
			if (i == 0) return _fInverse1;
			else if (i == 1) return _fInverse2;
			else return _fInverse3;
		}
	};

	/// @brief Transformation to orthogonal Cartesian coordinate system
	/// 
	/// Transforms from original Cartesian system to another orthogonal Cartesian
	/// system defined by three mutually orthogonal basis vectors.
	/// Verifies orthogonality and provides right/left-handed orientation check.
	class CoordTransf3DCartOrthogonal : public CoordTransfWithInverse<Vec3Cart, Vec3Cart, 3>
	{
	private:
		Vec3Cart _base[3];

		MatrixNM<Real, 3, 3> _transf;
		MatrixNM<Real, 3, 3> _transfInverse;

		const ScalarFunctionFromStdFunc<3> _f1, _f2, _f3;
		const ScalarFunctionFromStdFunc<3> _fInverse1, _fInverse2, _fInverse3;

		Real func1(const VectorN<Real, 3>& q) const { return (_transf * q)[0]; }
		Real func2(const VectorN<Real, 3>& q) const { return (_transf * q)[1]; }
		Real func3(const VectorN<Real, 3>& q) const { return (_transf * q)[2]; }

		Real funcInverse1(const VectorN<Real, 3>& q) const { return (_transfInverse * q)[0]; }
		Real funcInverse2(const VectorN<Real, 3>& q) const { return (_transfInverse * q)[1]; }
		Real funcInverse3(const VectorN<Real, 3>& q) const { return (_transfInverse * q)[2]; }

	public:
		/// @brief Constructor for orthogonal transformation
		/// @param b1 First basis vector (must be orthogonal to b2, b3)
		/// @param b2 Second basis vector (must be orthogonal to b1, b3)
		/// @param b3 Third basis vector (must be orthogonal to b1, b2)
		/// @throws GeometryError if basis vectors are not mutually orthogonal
		CoordTransf3DCartOrthogonal(const VectorN<Real, 3>& b1, const VectorN<Real, 3>& b2, const VectorN<Real, 3>& b3) :
			_f1([this](const VectorN<Real, 3>& q) { return func1(q); }),
			_f2([this](const VectorN<Real, 3>& q) { return func2(q); }),
			_f3([this](const VectorN<Real, 3>& q) { return func3(q); }),
			_fInverse1([this](const VectorN<Real, 3>& q) { return funcInverse1(q); }),
			_fInverse2([this](const VectorN<Real, 3>& q) { return funcInverse2(q); }),
			_fInverse3([this](const VectorN<Real, 3>& q) { return funcInverse3(q); }) 
		{
			// Convert to Vec3Cart for easier operations
			Vec3Cart base1(b1);
			Vec3Cart base2(b2);
			Vec3Cart base3(b3);

			// Check orthogonality (dot products should be zero)
			constexpr Real eps = PrecisionValues<Real>::OrthogonalityTolerance;
			if (std::abs(base1 * base2) > eps ||
					std::abs(base1 * base3) > eps ||
					std::abs(base2 * base3) > eps)
			{
				throw GeometryError("CoordTransf3DCartOrthogonal: Provided base vectors are not mutually orthogonal.");
			}

			_base[0] = b1;
			_base[1] = b2;
			_base[2] = b3;

			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++)
				{
					_transf(i, j) = _base[i][j];
					_transfInverse(i, j) = _base[j][i];
				}
			}
		}

		MatrixNM<Real, 3, 3> getTransfMatrix() { return _transf; }
		MatrixNM<Real, 3, 3> getInvTransfMatrix() { return _transfInverse; }

		Vec3Cart    transf(const Vec3Cart& q) const { return Vec3Cart{ func1(q), func2(q), func3(q) }; }
		Vec3Cart    transfInverse(const Vec3Cart& q) const { return Vec3Cart{ funcInverse1(q), funcInverse2(q), funcInverse3(q) }; }

		const IScalarFunction<3>& coordTransfFunc(int i) const
		{
			if (i == 0) return _f1;
			else if (i == 1) return _f2;
			else return _f3;
		}
		const IScalarFunction<3>& inverseCoordTransfFunc(int i) const
		{
			if (i == 0) return _fInverse1;
			else if (i == 1) return _fInverse2;
			else return _fInverse3;
		}

		bool IsRightHanded()
		{
			Vec3Cart cross = VectorProduct(_base[0], _base[1]);
			if (cross.ScalarProduct(_base[2]) > 0.0)
				return true;
			else
				return false;
		}
	};

	/// @brief General 3D Cartesian transformation from matrix
	/// 
	/// General 3D Cartesian coordinate transformation specified by an arbitrary 3x3 matrix.
	/// Handles any linear transformation including rotation, scaling, shearing, and combinations.
	/// The inverse transformation uses the matrix inverse (computed via LU decomposition).
	/// @note For orthogonal matrices (rotations), the inverse equals the transpose, but this
	///       class handles the general non-orthogonal case correctly.
	class CoordTransf3DCartGeneral : public CoordTransfWithInverse<Vec3Cart, Vec3Cart, 3>
	{
	private:
		Vec3Cart _base[3];

		MatrixNM<Real, 3, 3> _transf;
		MatrixNM<Real, 3, 3> _transfInverse;

		const ScalarFunctionFromStdFunc<3> _f1, _f2, _f3;
		const ScalarFunctionFromStdFunc<3> _fInverse1, _fInverse2, _fInverse3;

		Real func1(const VectorN<Real, 3>& q) const { return (_transf * q)[0]; }
		Real func2(const VectorN<Real, 3>& q) const { return (_transf * q)[1]; }
		Real func3(const VectorN<Real, 3>& q) const { return (_transf * q)[2]; }

		Real funcInverse1(const VectorN<Real, 3>& q) const { return (_transfInverse * q)[0]; }
		Real funcInverse2(const VectorN<Real, 3>& q) const { return (_transfInverse * q)[1]; }
		Real funcInverse3(const VectorN<Real, 3>& q) const { return (_transfInverse * q)[2]; }

	public:
		/// @brief Constructor from transformation matrix
		/// @param transfMat 3x3 transformation matrix (must be non-singular)
		/// @throws SingularMatrixError if matrix is singular (determinant ≈ 0)
		CoordTransf3DCartGeneral(const MatrixNM<Real, 3, 3>& transfMat) :
			_f1([this](const VectorN<Real, 3>& q) { return func1(q); }),
			_f2([this](const VectorN<Real, 3>& q) { return func2(q); }),
			_f3([this](const VectorN<Real, 3>& q) { return func3(q); }),
			_fInverse1([this](const VectorN<Real, 3>& q) { return funcInverse1(q); }),
			_fInverse2([this](const VectorN<Real, 3>& q) { return funcInverse2(q); }),
			_fInverse3([this](const VectorN<Real, 3>& q) { return funcInverse3(q); })
		{
			_transf = transfMat;

			// Store basis vectors (rows of transformation matrix)
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++)
				{
					_base[i][j] = _transf(i, j);
				}
			}

			// Compute proper matrix inverse (throws SingularMatrixError if singular)
			_transfInverse = _transf.GetInverse();
		}

		MatrixNM<Real, 3, 3> getTransfMatrix() { return _transf; }
		MatrixNM<Real, 3, 3> getInvTransfMatrix() { return _transfInverse; }

		Vec3Cart    transf(const Vec3Cart& q) const { return Vec3Cart{ func1(q), func2(q), func3(q) }; }
		Vec3Cart    transfInverse(const Vec3Cart& q) const { return Vec3Cart{ funcInverse1(q), funcInverse2(q), funcInverse3(q) }; }

		const IScalarFunction<3>& coordTransfFunc(int i) const
		{
			if (i == 0) return _f1;
			else if (i == 1) return _f2;
			else return _f3;
		}
		const IScalarFunction<3>& inverseCoordTransfFunc(int i) const
		{
			if (i == 0) return _fInverse1;
			else if (i == 1) return _fInverse2;
			else return _fInverse3;
		}

		bool IsRightHanded()
		{
			Vec3Cart cross = VectorProduct(_base[0], _base[1]);
			if (cross.ScalarProduct(_base[2]) > 0.0)
				return true;
			else
				return false;
		}
	};

	/// @brief Cartesian to oblique 3D coordinate transformation
	/// 
	/// Transforms from Cartesian to oblique (non-orthogonal) coordinate system.
	/// Uses dual basis for transformation: q? = dual_i � r
	/// Handles arbitrary non-orthogonal basis with automatic dual basis computation.
	class CoordTransfCartesianToOblique3D : public CoordTransfWithInverse<Vec3Cart, Vec3Cart, 3>
	{
	private:
		Vec3Cart _base[3];                                ///< Oblique basis vectors (covariant)
		Vec3Cart _dual[3];                                ///< Dual basis vectors (contravariant)
		MatrixNM<Real, 3, 3> _baseMat;                    ///< Matrix with basis vectors as columns
		MatrixNM<Real, 3, 3> _baseMatInv;                 ///< Inverse of basis matrix

		const ScalarFunctionFromStdFunc<3> _f1, _f2, _f3;
		const ScalarFunctionFromStdFunc<3> _fInverse1, _fInverse2, _fInverse3;

		Real func1(const VectorN<Real, 3>& q) const { return Utils::ScalarProduct<3>(q, static_cast<const VectorN<Real, 3>&>(_dual[0])); }
		Real func2(const VectorN<Real, 3>& q) const { return Utils::ScalarProduct<3>(q, static_cast<const VectorN<Real, 3>&>(_dual[1])); }
		Real func3(const VectorN<Real, 3>& q) const { return Utils::ScalarProduct<3>(q, static_cast<const VectorN<Real, 3>&>(_dual[2])); }

		Real funcInverse1(const VectorN<Real, 3>& q) const { return (_baseMat * q)[0]; }
		Real funcInverse2(const VectorN<Real, 3>& q) const { return (_baseMat * q)[1]; }
		Real funcInverse3(const VectorN<Real, 3>& q) const { return (_baseMat * q)[2]; }

	public:
		/// @brief Constructor for oblique coordinate system
		/// @param b1 First basis vector (need not be orthogonal)
		/// @param b2 Second basis vector (need not be orthogonal)
		/// @param b3 Third basis vector (need not be orthogonal)
		/// @throws SingularMatrixError if basis vectors are linearly dependent
		CoordTransfCartesianToOblique3D(const VectorN<Real, 3>& b1, const VectorN<Real, 3>& b2, const VectorN<Real, 3>& b3) :
			_f1([this](const VectorN<Real, 3>& q) { return func1(q); }),
			_f2([this](const VectorN<Real, 3>& q) { return func2(q); }),
			_f3([this](const VectorN<Real, 3>& q) { return func3(q); }),
			_fInverse1([this](const VectorN<Real, 3>& q) { return funcInverse1(q); }),
			_fInverse2([this](const VectorN<Real, 3>& q) { return funcInverse2(q); }),
			_fInverse3([this](const VectorN<Real, 3>& q) { return funcInverse3(q); })
		{
			_base[0] = b1; _base[1] = b2; _base[2] = b3;

			for (int i = 0; i < 3; ++i) {
				_baseMat(0, i) = _base[i][0];
				_baseMat(1, i) = _base[i][1];
				_baseMat(2, i) = _base[i][2];
			}
			// Check for non-degeneracy
			Real V = ScalarProduct(_base[0], VectorProduct(_base[1], _base[2]));
			if (std::abs(V) < PrecisionValues<Real>::LinearDependenceTolerance)
				throw SingularMatrixError("Oblique base vectors are linearly dependent!", V);

		// Compute dual basis
		_dual[0] = (1 / V) * VectorProduct(_base[1], _base[2]);
		_dual[1] = (1 / V) * VectorProduct(_base[2], _base[0]);
		_dual[2] = (1 / V) * VectorProduct(_base[0], _base[1]);

		// Compute inverse of base matrix for alternative transformation method
		_baseMatInv = _baseMat.GetInverse();
	}		/// @brief Get oblique basis vector
		/// @param i Index (0, 1, or 2)
		/// @return Oblique basis vector i
		Vec3Cart    Base(int i) { return _base[i]; }
		
		/// @brief Get dual basis vector (contravariant)
		/// @param i Index (0, 1, or 2)
		/// @return Dual basis vector i
		Vec3Cart    Dual(int i) { return _dual[i]; }

		/// @brief Transform Cartesian to oblique coordinates
		/// @param cartesian Cartesian position vector
		/// @return Oblique coordinates q? = dual_i � cartesian
		Vec3Cart transf(const Vec3Cart& cartesian) const {
			// q^i = dual_i � cartesian
			return Vec3Cart{
					ScalarProduct(cartesian, _dual[0]),
					ScalarProduct(cartesian, _dual[1]),
					ScalarProduct(cartesian, _dual[2])
			};
			// Or, if you have _baseMatInv: return _baseMatInv * cartesian;
		}

		/// @brief Transform oblique to Cartesian coordinates
		/// @param oblique Oblique coordinates q?
		/// @return Cartesian position r = q? � base_i
		Vec3Cart transfInverse(const Vec3Cart& oblique) const {
			// r = q^i * base_i
			return _baseMat * oblique;
		}
		
		const IScalarFunction<3>& coordTransfFunc(int i) const
		{
			if (i == 0) return _f1;
			else if (i == 1) return _f2;
			else return _f3;
		}
		const IScalarFunction<3>& inverseCoordTransfFunc(int i) const
		{
			if (i == 0) return _fInverse1;
			else if (i == 1) return _fInverse2;
			else return _fInverse3;
		}

		bool IsRightHanded()
		{
			Vec3Cart cross = VectorProduct(_base[0], _base[1]);
			if (ScalarProduct(cross, _base[2]) > 0.0)
				return true;
			else
				return false;
		}
	};

	/*****************************************************************************
	 * QUATERNION-BASED 3D ROTATION COORDINATE TRANSFORMATION
	 * 
	 * This class provides 3D coordinate transformations using quaternions,
	 * offering several advantages over matrix-based rotations:
	 * 
	 * ADVANTAGES:
	 * - Compact representation: 4 values instead of 9 matrix elements
	 * - No gimbal lock issues (unlike Euler angles)
	 * - Efficient composition: quaternion multiplication is faster than matrix
	 * - Smooth interpolation: SLERP provides natural rotation interpolation
	 * - Numerical stability: quaternions are easier to renormalize
	 * 
	 * USAGE:
	 * 1. Construct from quaternion:
	 *    Quaternion q = Quaternion::FromAxisAngle(Vec3Cart(0,0,1), PI/4);
	 *    CoordTransfCart3DRotationQuaternion rot(q);
	 * 
	 * 2. Construct from axis-angle:
	 *    CoordTransfCart3DRotationQuaternion rot(Vec3Cart(1,0,0), PI/2);
	 * 
	 * 3. Construct from rotation matrix:
	 *    MatrixNM<Real,3,3> mat = ...;
	 *    CoordTransfCart3DRotationQuaternion rot(mat);
	 * 
	 * 4. Apply transformation:
	 *    Vec3Cart v_rotated = rot.transf(v);
	 *    Vec3Cart v_original = rot.transfInverse(v_rotated);
	 * 
	 * INTEGRATION WITH OTHER ROTATION CLASSES:
	 * This class seamlessly integrates with existing CoordTransf3D rotation
	 * classes (RotationXAxis, RotationYAxis, etc.) through shared matrix
	 * representation, allowing mixing of quaternion and matrix rotations.
	 *****************************************************************************/
	class CoordTransfCart3DRotationQuaternion : 
		public CoordTransfWithInverse<Vec3Cart, Vec3Cart, 3>
	{
	private:
		Quaternion _quat;
		MatrixNM<Real, 3, 3> _transf;
		MatrixNM<Real, 3, 3> _inverse;

		const ScalarFunctionFromStdFunc<3> _f1, _f2, _f3;
		const ScalarFunctionFromStdFunc<3> _fInverse1, _fInverse2, _fInverse3;

		// Initialize transformation matrices from quaternion
		void initializeMatrices()
		{
			// Ensure quaternion is normalized for proper rotation
			if (!_quat.IsUnit(PrecisionValues<Real>::DefaultTolerance))
				_quat.Normalize();

			// Convert quaternion to rotation matrix
			_transf = _quat.ToRotationMatrix();

			// Inverse of rotation quaternion is its conjugate
			// For matrices: inverse of orthogonal matrix is its transpose
			Quaternion quat_inv = _quat.Conjugate();
			_inverse = quat_inv.ToRotationMatrix();
		}

	public:
		// Construct from quaternion
		explicit CoordTransfCart3DRotationQuaternion(const Quaternion& q) 
			: _quat(q),
				_f1([this](const VectorN<Real, 3>& v) { return func1(v); }),
				_f2([this](const VectorN<Real, 3>& v) { return func2(v); }),
				_f3([this](const VectorN<Real, 3>& v) { return func3(v); }),
				_fInverse1([this](const VectorN<Real, 3>& v) { return funcInverse1(v); }),
				_fInverse2([this](const VectorN<Real, 3>& v) { return funcInverse2(v); }),
				_fInverse3([this](const VectorN<Real, 3>& v) { return funcInverse3(v); })
		{
			initializeMatrices();
		}

		// Construct from axis-angle representation
		CoordTransfCart3DRotationQuaternion(const Vec3Cart& axis, Real angle)
			: _quat(Quaternion::FromAxisAngle(axis, angle)),
				_f1([this](const VectorN<Real, 3>& v) { return func1(v); }),
				_f2([this](const VectorN<Real, 3>& v) { return func2(v); }),
				_f3([this](const VectorN<Real, 3>& v) { return func3(v); }),
				_fInverse1([this](const VectorN<Real, 3>& v) { return funcInverse1(v); }),
				_fInverse2([this](const VectorN<Real, 3>& v) { return funcInverse2(v); }),
				_fInverse3([this](const VectorN<Real, 3>& v) { return funcInverse3(v); })
		{
			initializeMatrices();
		}

		// Construct from rotation matrix
		explicit CoordTransfCart3DRotationQuaternion(const MatrixNM<Real, 3, 3>& rotMatrix)
			: _quat(Quaternion::FromRotationMatrix(rotMatrix)),
				_f1([this](const VectorN<Real, 3>& v) { return func1(v); }),
				_f2([this](const VectorN<Real, 3>& v) { return func2(v); }),
				_f3([this](const VectorN<Real, 3>& v) { return func3(v); }),
				_fInverse1([this](const VectorN<Real, 3>& v) { return funcInverse1(v); }),
				_fInverse2([this](const VectorN<Real, 3>& v) { return funcInverse2(v); }),
				_fInverse3([this](const VectorN<Real, 3>& v) { return funcInverse3(v); })
		{
			initializeMatrices();
		}

		// Construct from Euler angles (ZYX convention: yaw, pitch, roll)
		static CoordTransfCart3DRotationQuaternion FromEulerZYX(Real yaw, Real pitch, Real roll)
		{
			return CoordTransfCart3DRotationQuaternion(
				Quaternion::FromEulerZYX(yaw, pitch, roll)
			);
		}

		// Construct from Euler angles (XYZ convention)
		static CoordTransfCart3DRotationQuaternion FromEulerXYZ(Real roll, Real pitch, Real yaw)
		{
			return CoordTransfCart3DRotationQuaternion(
				Quaternion::FromEulerXYZ(roll, pitch, yaw)
			);
		}

		// Transformation functions using matrix representation
		Real func1(const VectorN<Real, 3>& v) const 
		{ 
			return _transf[0][0] * v[0] + _transf[0][1] * v[1] + _transf[0][2] * v[2]; 
		}
		
		Real func2(const VectorN<Real, 3>& v) const 
		{ 
			return _transf[1][0] * v[0] + _transf[1][1] * v[1] + _transf[1][2] * v[2]; 
		}
		
		Real func3(const VectorN<Real, 3>& v) const 
		{ 
			return _transf[2][0] * v[0] + _transf[2][1] * v[1] + _transf[2][2] * v[2]; 
		}

		// Inverse transformation functions
		Real funcInverse1(const VectorN<Real, 3>& v) const 
		{ 
			return _inverse[0][0] * v[0] + _inverse[0][1] * v[1] + _inverse[0][2] * v[2]; 
		}
		
		Real funcInverse2(const VectorN<Real, 3>& v) const 
		{ 
			return _inverse[1][0] * v[0] + _inverse[1][1] * v[1] + _inverse[1][2] * v[2]; 
		}
		
		Real funcInverse3(const VectorN<Real, 3>& v) const 
		{ 
			return _inverse[2][0] * v[0] + _inverse[2][1] * v[1] + _inverse[2][2] * v[2]; 
		}

		// Apply transformation (using quaternion directly for efficiency)
		Vec3Cart transf(const Vec3Cart& v) const 
		{ 
			return _quat.Rotate(v);
		}

		// Apply inverse transformation
		Vec3Cart transfInverse(const Vec3Cart& v) const 
		{ 
			return _quat.Conjugate().Rotate(v);
		}

		// Interface implementation
		const IScalarFunction<3>& coordTransfFunc(int i) const
		{
			if (i == 0) return _f1;
			else if (i == 1) return _f2;
			else return _f3;
		}

		const IScalarFunction<3>& inverseCoordTransfFunc(int i) const
		{
			if (i == 0) return _fInverse1;
			else if (i == 1) return _fInverse2;
			else return _fInverse3;
		}

		// Accessors
		const Quaternion& GetQuaternion() const { return _quat; }
		const MatrixNM<Real, 3, 3>& GetTransformationMatrix() const { return _transf; }
		const MatrixNM<Real, 3, 3>& GetInverseMatrix() const { return _inverse; }

		// Get rotation axis and angle
		Vec3Cart GetRotationAxis() const { return _quat.GetRotationAxis(); }
		Real GetRotationAngle() const { return _quat.GetRotationAngle(); }

		// Composition: compose this rotation with another
		// Returns a new transformation representing this rotation followed by other
		CoordTransfCart3DRotationQuaternion Compose(
			const CoordTransfCart3DRotationQuaternion& other) const
		{
			return CoordTransfCart3DRotationQuaternion(other._quat * _quat);
		}

		// Interpolation between two rotations
		// t=0 returns this, t=1 returns other
		// Uses SLERP for smooth interpolation
		CoordTransfCart3DRotationQuaternion Interpolate(
			const CoordTransfCart3DRotationQuaternion& other, Real t) const
		{
			Quaternion interpolated = Quaternion::Slerp(_quat, other._quat, t);
			return CoordTransfCart3DRotationQuaternion(interpolated);
		}
	};
}

#endif

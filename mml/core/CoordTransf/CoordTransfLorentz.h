///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        CoordTransfLorentz.h                                                ///
///  Description: Lorentz transformations for special relativity                      ///
///               Boost transformations and Minkowski spacetime                       ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_COORD_TRANSF_LORENTZ_H
#define MML_COORD_TRANSF_LORENTZ_H

#include "MMLBase.h"

#include "interfaces/ICoordTransf.h"

#include "base/Vector/Vector.h"
#include "base/Vector/VectorTypes.h"
#include "base/Matrix/Matrix.h"
#include "base/Matrix/MatrixNM.h"
#include "base/Tensor.h"
#include "mml/base/Geometry/Geometry.h"

#include "core/CoordTransf.h"

namespace MML {
	/// @brief Lorentz boost transformation along X-axis in Minkowski spacetime
	/// @details Implements special relativistic coordinate transformation for inertial frames
	///          moving with constant velocity β (in units of c) along x-axis.
	///          Uses Lorentz factor γ = 1/√(1-β²), transforms (ct, x, y, z) coordinates.
	class CoordTransfLorentzXAxis : public CoordTransfWithInverse<Vector4Minkowski, Vector4Minkowski, 4> {
	private:
		/// @brief Velocity in units of c (speed of light), must be in [0, 1)
		Real _velocity;
		MatrixNM<Real, 4, 4> _transf;
		MatrixNM<Real, 4, 4> _inverse;

		/// @brief Scalar functions for forward transformation components
		const ScalarFunctionFromStdFunc<4> _f1, _f2, _f3, _f4;
		/// @brief Scalar functions for inverse transformation components
		const ScalarFunctionFromStdFunc<4> _fInv1, _fInv2, _fInv3, _fInv4;

	public:
		/// @brief Construct Lorentz boost along X-axis
		/// @param inVelocity Velocity in units of c, must be in [0, 1)
		/// @throws std::range_error if velocity is outside valid range
		CoordTransfLorentzXAxis(Real inVelocity)
				: _velocity(inVelocity)
				, _f1([this](const VectorN<Real, 4>& q) { return func1(q); })
				, _f2([this](const VectorN<Real, 4>& q) { return func2(q); })
				, _f3([this](const VectorN<Real, 4>& q) { return func3(q); })
				, _f4([this](const VectorN<Real, 4>& q) { return func4(q); })
				, _fInv1([this](const VectorN<Real, 4>& q) { return funcInverse1(q); })
				, _fInv2([this](const VectorN<Real, 4>& q) { return funcInverse2(q); })
				, _fInv3([this](const VectorN<Real, 4>& q) { return funcInverse3(q); })
				, _fInv4([this](const VectorN<Real, 4>& q) { return funcInverse4(q); }) {
			if (_velocity < 0.0 || _velocity >= 1.0)
				throw std::range_error("Invalid velocity for Lorentz transformation: must be in [0, 1).");

			Real gamma = REAL(1.0) / sqrt(REAL(1.0) - _velocity * _velocity);
			Real beta = _velocity;

			_transf[0][0] = gamma;
			_transf[0][1] = -beta * gamma;
			_transf[1][0] = -beta * gamma;
			_transf[1][1] = gamma;
			_transf[2][2] = 1.0;
			_transf[3][3] = 1.0;

			_inverse[0][0] = gamma;
			_inverse[0][1] = beta * gamma;
			_inverse[1][0] = beta * gamma;
			_inverse[1][1] = gamma;
			_inverse[2][2] = 1.0;
			_inverse[3][3] = 1.0;
		}

		/// @brief Forward transformation component functions (apply _transf matrix)
		Real func1(const VectorN<Real, 4>& q) const { return (_transf * q)[0]; }
		Real func2(const VectorN<Real, 4>& q) const { return (_transf * q)[1]; }
		Real func3(const VectorN<Real, 4>& q) const { return (_transf * q)[2]; }
		Real func4(const VectorN<Real, 4>& q) const { return (_transf * q)[3]; }

		/// @brief Inverse transformation component functions (apply _inverse matrix)
		Real funcInverse1(const VectorN<Real, 4>& q) const { return (_inverse * q)[0]; }
		Real funcInverse2(const VectorN<Real, 4>& q) const { return (_inverse * q)[1]; }
		Real funcInverse3(const VectorN<Real, 4>& q) const { return (_inverse * q)[2]; }
		Real funcInverse4(const VectorN<Real, 4>& q) const { return (_inverse * q)[3]; }

		/// @brief Apply forward Lorentz transformation
		/// @param q Minkowski 4-vector (ct, x, y, z)
		/// @return Transformed 4-vector in boosted frame
		Vector4Minkowski transf(const Vector4Minkowski& q) const { return Vector4Minkowski{func1(q), func2(q), func3(q), func4(q)}; }
		/// @brief Apply inverse Lorentz transformation (boost in opposite direction)
		/// @param q Minkowski 4-vector in boosted frame
		/// @return Transformed 4-vector in original frame
		Vector4Minkowski transfInverse(const Vector4Minkowski& q) const {
			return Vector4Minkowski{funcInverse1(q), funcInverse2(q), funcInverse3(q), funcInverse4(q)};
		}

		/// @brief Get i-th component of forward transformation as scalar function
		const IScalarFunction<4>& coordTransfFunc(int i) const {
			if (i == 0)
				return _f1;
			else if (i == 1)
				return _f2;
			else if (i == 2)
				return _f3;
			else
				return _f4;
		}
		/// @brief Get i-th component of inverse transformation as scalar function
		const IScalarFunction<4>& inverseCoordTransfFunc(int i) const {
			if (i == 0)
				return _fInv1;
			else if (i == 1)
				return _fInv2;
			else if (i == 2)
				return _fInv3;
			else
				return _fInv4;
		}
	};
} // namespace MML

#endif
///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        CoordTransfLorentz.h                                                ///
///  Description: Lorentz transformations for special relativity                      ///
///               Boost transformations and Minkowski spacetime                       ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_COORD_TRANSF_LORENTZ_H
#define MML_COORD_TRANSF_LORENTZ_H

#include "MMLBase.h"

#include "interfaces/ICoordTransf.h"

#include "base/Vector.h"
#include "base/VectorTypes.h"
#include "base/Matrix.h"
#include "base/MatrixNM.h"
#include "base/Tensor.h"
#include "base/Geometry.h"

#include "core/CoordTransf.h"

namespace MML
{
	class CoordTransfLorentzXAxis : public CoordTransfWithInverse<Vector4Minkowski, Vector4Minkowski, 4>
	{
	private:
		Real    _velocity;			// expressed in units of c (speed of light)
		MatrixNM<Real, 4, 4>  _transf;
		MatrixNM<Real, 4, 4>  _inverse;

		const ScalarFunctionFromStdFunc<4> _f1, _f2, _f3, _f4;
		const ScalarFunctionFromStdFunc<4> _fInv1, _fInv2, _fInv3, _fInv4;

	public:
		CoordTransfLorentzXAxis(Real inVelocity) : _velocity(inVelocity),
			_f1(std::function<Real(const VectorN<Real, 4>&)> { std::bind(&CoordTransfLorentzXAxis::func1, this, std::placeholders::_1) }),
			_f2(std::function<Real(const VectorN<Real, 4>&)> { std::bind(&CoordTransfLorentzXAxis::func2, this, std::placeholders::_1) }),
			_f3(std::function<Real(const VectorN<Real, 4>&)> { std::bind(&CoordTransfLorentzXAxis::func3, this, std::placeholders::_1) }),
			_f4(std::function<Real(const VectorN<Real, 4>&)> { std::bind(&CoordTransfLorentzXAxis::func4, this, std::placeholders::_1) }),
			_fInv1(std::function<Real(const VectorN<Real, 4>&)> { std::bind(&CoordTransfLorentzXAxis::funcInverse1, this, std::placeholders::_1) }),
			_fInv2(std::function<Real(const VectorN<Real, 4>&)> { std::bind(&CoordTransfLorentzXAxis::funcInverse2, this, std::placeholders::_1) }),
			_fInv3(std::function<Real(const VectorN<Real, 4>&)> { std::bind(&CoordTransfLorentzXAxis::funcInverse3, this, std::placeholders::_1) }),
			_fInv4(std::function<Real(const VectorN<Real, 4>&)> { std::bind(&CoordTransfLorentzXAxis::funcInverse4, this, std::placeholders::_1) })
		{
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

		Real func1(const VectorN<Real, 4>& q) const { return (_transf * q)[0]; }
		Real func2(const VectorN<Real, 4>& q) const { return (_transf * q)[1]; }
		Real func3(const VectorN<Real, 4>& q) const { return (_transf * q)[2]; }
		Real func4(const VectorN<Real, 4>& q) const { return (_transf * q)[3]; }

		Real funcInverse1(const VectorN<Real, 4>& q) const { return (_inverse * q)[0]; }
		Real funcInverse2(const VectorN<Real, 4>& q) const { return (_inverse * q)[1]; }
		Real funcInverse3(const VectorN<Real, 4>& q) const { return (_inverse * q)[2]; }
		Real funcInverse4(const VectorN<Real, 4>& q) const { return (_inverse * q)[3]; }

		Vector4Minkowski    transf(const Vector4Minkowski& q) const { return Vector4Minkowski{ func1(q), func2(q), func3(q), func4(q) }; }
		Vector4Minkowski    transfInverse(const Vector4Minkowski& q) const { return Vector4Minkowski{ funcInverse1(q), funcInverse2(q), funcInverse3(q), funcInverse4(q) }; }

		const IScalarFunction<4>& coordTransfFunc(int i) const
		{
			if (i == 0) return _f1;
			else if (i == 1) return _f2;
			else if (i == 2) return _f3;
			else return _f4;
		}
		const IScalarFunction<4>& inverseCoordTransfFunc(int i) const
		{
			if (i == 0) return _fInv1;
			else if (i == 1) return _fInv2;
			else if (i == 2) return _fInv3;
			else return _fInv4;
		}
	};
}

#endif 
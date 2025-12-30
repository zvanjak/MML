///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        CoordTransf2D.h                                                     ///
///  Description: 2D coordinate transformations (Polar, rotation, scaling)            ///
///               Translation, reflection, and composite transformations              ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_COORD_TRANSF_2D_H
#define MML_COORD_TRANSF_2D_H

#include "MMLBase.h"

#include "core/CoordTransf.h"


namespace MML
{
	class CoordTransfPolarToCartesian2D : public CoordTransfWithInverse<Vector2Polar, Vector2Cartesian, 2>
	{
		// q[0] = r     - radial distance
		// q[1] = phi   - polar angle
	public:
		static Real func1(const VectorN<Real, 2>& q) { return q[0] * cos(q[1]); }
		static Real func2(const VectorN<Real, 2>& q) { return q[0] * sin(q[1]); }

		// q[0] = x
		// q[1] = y
		static Real funcInverse1(const VectorN<Real, 2>& q) { return sqrt(q[0] * q[0] + q[1] * q[1]); }
		static Real funcInverse2(const VectorN<Real, 2>& q) { return atan2(q[1], q[0]); }

		inline static ScalarFunction<2> _func[2] = { ScalarFunction<2>{func1},
																								 ScalarFunction<2>{func2}
		};

		inline static ScalarFunction<2> _funcInverse[2] = { ScalarFunction<2>{funcInverse1},
																												ScalarFunction<2>{funcInverse2}
		};

		Vector2Cartesian     transf(const Vector2Polar& q) const { return Vector2Cartesian{ func1(q), func2(q) }; }
		Vector2Polar         transfInverse(const Vector2Cartesian& q) const { return Vector2Polar{ funcInverse1(q), funcInverse2(q) }; }

		const IScalarFunction<2>& coordTransfFunc(int i) const { return _func[i]; }
		const IScalarFunction<2>& inverseCoordTransfFunc(int i) const { return _funcInverse[i]; }
	};

	class CoordTransfCart2DRotation : public CoordTransfWithInverse<Vector2Cartesian, Vector2Cartesian, 2>
	{
	private:
		Real    _angle;
		MatrixNM<Real, 2, 2>  _transf;
		MatrixNM<Real, 2, 2>  _inverse;

		const ScalarFunctionFromStdFunc<2> _f1;
		const ScalarFunctionFromStdFunc<2> _f2;

		const ScalarFunctionFromStdFunc<2> _fInverse1;
		const ScalarFunctionFromStdFunc<2> _fInverse2;

	public:
		CoordTransfCart2DRotation(Real inAngle) : 
			_angle(inAngle),
			_f1([this](const VectorN<Real, 2>& q) { return func1(q); }),
			_f2([this](const VectorN<Real, 2>& q) { return func2(q); }),
			_fInverse1([this](const VectorN<Real, 2>& q) { return funcInverse1(q); }),
			_fInverse2([this](const VectorN<Real, 2>& q) { return funcInverse2(q); })
			//_angle(inAngle),
			//_f1(std::function<Real(const VectorN<Real, 2>&)> { std::bind(&CoordTransfCart2DRotation::func1, this, std::placeholders::_1) }),
			//_f2(std::function<Real(const VectorN<Real, 2>&)> { std::bind(&CoordTransfCart2DRotation::func2, this, std::placeholders::_1) }),
			//_fInverse1(std::function<Real(const VectorN<Real, 2>&)> { std::bind(&CoordTransfCart2DRotation::funcInverse1, this, std::placeholders::_1) }),
			//_fInverse2(std::function<Real(const VectorN<Real, 2>&)> { std::bind(&CoordTransfCart2DRotation::funcInverse2, this, std::placeholders::_1) })
		{
			_transf[0][0] = cos(_angle);
			_transf[0][1] = -sin(_angle);
			_transf[1][0] = sin(_angle);
			_transf[1][1] = cos(_angle);

			_inverse[0][0] = cos(_angle);
			_inverse[0][1] = sin(_angle);
			_inverse[1][0] = -sin(_angle);
			_inverse[1][1] = cos(_angle);
		}

		Real func1(const VectorN<Real, 2>& q) const { return _transf[0][0] * q[0] + _transf[0][1] * q[1]; }
		Real func2(const VectorN<Real, 2>& q) const { return (_transf * q)[1]; }

		Real funcInverse1(const VectorN<Real, 2>& q) const { return (_inverse * q)[0]; }
		Real funcInverse2(const VectorN<Real, 2>& q) const { return (_inverse * q)[1]; }

		Vector2Cartesian    transf(const Vector2Cartesian& q) const { return Vector2Cartesian{ func1(q), func2(q) }; }
		Vector2Cartesian    transfInverse(const Vector2Cartesian& q) const { return Vector2Cartesian{ funcInverse1(q), funcInverse2(q) }; }

		const IScalarFunction<2>& coordTransfFunc(int i) const
		{
			if (i == 0) return _f1;
			else return _f2;
		}
		const IScalarFunction<2>& inverseCoordTransfFunc(int i) const
		{
			if (i == 0) return _fInverse1;
			else return _fInverse2;
		}
	};

	// oblique transformation in 2D, where Y axis is not at angle 90 degrees to X axis
	// but is at some angle _angle to X axis
	class CoordTransfObliqueToCartesian2D : public CoordTransfWithInverse<Vector2Cartesian, Vector2Cartesian, 2>
	{
	private:
		Real    _angle;

		const ScalarFunctionFromStdFunc<2> _f1;
		const ScalarFunctionFromStdFunc<2> _f2;

		const ScalarFunctionFromStdFunc<2> _fInverse1;
		const ScalarFunctionFromStdFunc<2> _fInverse2;

	public:
		CoordTransfObliqueToCartesian2D(Real inAngle)
			: _angle(inAngle),
			_f1([this](const VectorN<Real, 2>& q) { return func1(q); }),
			_f2([this](const VectorN<Real, 2>& q) { return func2(q); }),
			_fInverse1([this](const VectorN<Real, 2>& q) { return funcInverse1(q); }),
			_fInverse2([this](const VectorN<Real, 2>& q) { return funcInverse2(q); })
		{
		}

		// Forward: oblique (q) → Cartesian (x, y)
		Real func1(const VectorN<Real, 2>& q) const { return q[0] + q[1] * std::cos(_angle); }
		Real func2(const VectorN<Real, 2>& q) const { return q[1] * std::sin(_angle); }

		// Inverse: Cartesian (x, y) → oblique (q)
		Real funcInverse2(const VectorN<Real, 2>& q) const { return q[1] / std::sin(_angle); }
		Real funcInverse1(const VectorN<Real, 2>& q) const { return q[0] - funcInverse2(q) * std::cos(_angle); }

		Vector2Cartesian transf(const Vector2Cartesian& q) const override { return Vector2Cartesian{ func1(q), func2(q) }; }
		Vector2Cartesian transfInverse(const Vector2Cartesian& q) const override { return Vector2Cartesian{ funcInverse1(q), funcInverse2(q) }; }

		const IScalarFunction<2>& coordTransfFunc(int i) const override
		{
			return (i == 0) ? _f1 : _f2;
		}
		const IScalarFunction<2>& inverseCoordTransfFunc(int i) const override
		{
			return (i == 0) ? _fInverse1 : _fInverse2;
		}
	};


	class CoordTransfCartesianToOblique2D : public CoordTransfWithInverse<Vector2Cartesian, Vector2Cartesian, 2>
	{
	private:
		Real    _angle;

		const ScalarFunctionFromStdFunc<2> _f1;
		const ScalarFunctionFromStdFunc<2> _f2;

		const ScalarFunctionFromStdFunc<2> _fInverse1;
		const ScalarFunctionFromStdFunc<2> _fInverse2;

	public:
		CoordTransfCartesianToOblique2D(Real inAngle)
			: _angle(inAngle),
			_f1([this](const VectorN<Real, 2>& q) { return func1(q); }),
			_f2([this](const VectorN<Real, 2>& q) { return func2(q); }),
			_fInverse1([this](const VectorN<Real, 2>& q) { return funcInverse1(q); }),
			_fInverse2([this](const VectorN<Real, 2>& q) { return funcInverse2(q); })
		{	}

		// Forward: Cartesian (x, y) → oblique (q)
		Real func2(const VectorN<Real, 2>& q) const { return q[1] / std::sin(_angle); }
		Real func1(const VectorN<Real, 2>& q) const { return q[0] - func2(q) * std::cos(_angle); }

		// Inverse: oblique (q) → Cartesian (x, y)
		Real funcInverse1(const VectorN<Real, 2>& q) const { return q[0] + q[1] * std::cos(_angle); }
		Real funcInverse2(const VectorN<Real, 2>& q) const { return q[1] * std::sin(_angle); }

		Vector2Cartesian transf(const Vector2Cartesian& q) const override { return Vector2Cartesian{ func1(q), func2(q) }; }
		Vector2Cartesian transfInverse(const Vector2Cartesian& q) const override { return Vector2Cartesian{ funcInverse1(q), funcInverse2(q) }; }

		const IScalarFunction<2>& coordTransfFunc(int i) const override
		{
			return (i == 0) ? _f1 : _f2;
		}
		const IScalarFunction<2>& inverseCoordTransfFunc(int i) const override
		{
			return (i == 0) ? _fInverse1 : _fInverse2;
		}
	};

}

#endif

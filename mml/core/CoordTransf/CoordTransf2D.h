///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        CoordTransf2D.h                                                     ///
///  Description: 2D coordinate transformations (Polar, rotation, scaling)            ///
///               Translation, reflection, and composite transformations              ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_COORD_TRANSF_2D_H
#define MML_COORD_TRANSF_2D_H

#include "MMLBase.h"

#include "core/CoordTransf.h"


namespace MML
{
	/// @brief 2D polar to Cartesian coordinate transformation
	/// 
	/// Transforms between polar coordinates (r, f) and Cartesian coordinates (x, y).
	/// Forward: x = r�cos(f), y = r�sin(f)
	/// Inverse: r = v(x�+y�), f = atan2(y, x)
	class CoordTransfPolarToCartesian2D : public CoordTransfWithInverse<Vector2Polar, Vector2Cartesian, 2>
	{
		// q[0] = r     - radial distance
		// q[1] = phi   - polar angle
	public:
		/// @brief Compute x coordinate from polar: x = r�cos(f)
		/// @param q Polar coordinates (r, f)
		/// @return x coordinate
		static Real func1(const VectorN<Real, 2>& q) { return q[0] * cos(q[1]); }
		
		/// @brief Compute y coordinate from polar: y = r�sin(f)
		/// @param q Polar coordinates (r, f)
		/// @return y coordinate
		static Real func2(const VectorN<Real, 2>& q) { return q[0] * sin(q[1]); }

		// q[0] = x
		// q[1] = y
		/// @brief Compute radial distance: r = v(x�+y�)
		/// @param q Cartesian coordinates (x, y)
		/// @return Radial distance r
		static Real funcInverse1(const VectorN<Real, 2>& q) { return sqrt(q[0] * q[0] + q[1] * q[1]); }
		
		/// @brief Compute polar angle: f = atan2(y, x)
		/// @param q Cartesian coordinates (x, y)
		/// @return Polar angle f in radians
		static Real funcInverse2(const VectorN<Real, 2>& q) { return atan2(q[1], q[0]); }

		inline static ScalarFunction<2> _func[2] = { ScalarFunction<2>{func1},
																								 ScalarFunction<2>{func2}
		};

		inline static ScalarFunction<2> _funcInverse[2] = { ScalarFunction<2>{funcInverse1},
																												ScalarFunction<2>{funcInverse2}
		};

		/// @brief Transform polar coordinates to Cartesian
		/// @param q Polar coordinates (r, f)
		/// @return Cartesian coordinates (x, y)
		Vector2Cartesian     transf(const Vector2Polar& q) const { return Vector2Cartesian{ func1(q), func2(q) }; }
		
		/// @brief Transform Cartesian coordinates to polar
		/// @param q Cartesian coordinates (x, y)
		/// @return Polar coordinates (r, f)
		Vector2Polar         transfInverse(const Vector2Cartesian& q) const { return Vector2Polar{ funcInverse1(q), funcInverse2(q) }; }

		/// @brief Get forward transformation function component
		/// @param i Component index (0 or 1)
		/// @return Scalar function for component i
		const IScalarFunction<2>& coordTransfFunc(int i) const { return _func[i]; }
		
		/// @brief Get inverse transformation function component
		/// @param i Component index (0 or 1)
		/// @return Inverse scalar function for component i
		const IScalarFunction<2>& inverseCoordTransfFunc(int i) const { return _funcInverse[i]; }
	};

	/// @brief 2D Cartesian rotation transformation
	/// 
	/// Rotates Cartesian coordinates by a fixed angle around the origin.
	/// Uses rotation matrix: [cos(?) -sin(?); sin(?) cos(?)]
	/// Inverse transformation is rotation by -? (transpose of rotation matrix).
	class CoordTransfCart2DRotation : public CoordTransfWithInverse<Vector2Cartesian, Vector2Cartesian, 2>
	{
	private:
		Real    _angle;                                   ///< Rotation angle in radians
		MatrixNM<Real, 2, 2>  _transf;                    ///< Forward rotation matrix
		MatrixNM<Real, 2, 2>  _inverse;                   ///< Inverse rotation matrix (transpose)

		const ScalarFunctionFromStdFunc<2> _f1;
		const ScalarFunctionFromStdFunc<2> _f2;

		const ScalarFunctionFromStdFunc<2> _fInverse1;
		const ScalarFunctionFromStdFunc<2> _fInverse2;

	public:
		/// @brief Constructor for 2D rotation transformation
		/// @param inAngle Rotation angle in radians (counterclockwise positive)
		CoordTransfCart2DRotation(Real inAngle) : 
			_angle(inAngle),
			_f1([this](const VectorN<Real, 2>& q) { return func1(q); }),
			_f2([this](const VectorN<Real, 2>& q) { return func2(q); }),
			_fInverse1([this](const VectorN<Real, 2>& q) { return funcInverse1(q); }),
			_fInverse2([this](const VectorN<Real, 2>& q) { return funcInverse2(q); })
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

		/// @brief Apply rotation: x' = cos(?)x - sin(?)y
		/// @param q Input Cartesian coordinates
		/// @return x component after rotation
		Real func1(const VectorN<Real, 2>& q) const { return _transf[0][0] * q[0] + _transf[0][1] * q[1]; }
		
		/// @brief Apply rotation: y' = sin(?)x + cos(?)y
		/// @param q Input Cartesian coordinates
		/// @return y component after rotation
		Real func2(const VectorN<Real, 2>& q) const { return (_transf * q)[1]; }

		/// @brief Apply inverse rotation (rotation by -?)
		/// @param q Rotated Cartesian coordinates
		/// @return x component after inverse rotation
		Real funcInverse1(const VectorN<Real, 2>& q) const { return (_inverse * q)[0]; }
		
		/// @brief Apply inverse rotation (rotation by -?)
		/// @param q Rotated Cartesian coordinates
		/// @return y component after inverse rotation
		Real funcInverse2(const VectorN<Real, 2>& q) const { return (_inverse * q)[1]; }

		/// @brief Apply rotation transformation
		/// @param q Input Cartesian coordinates
		/// @return Rotated Cartesian coordinates
		Vector2Cartesian    transf(const Vector2Cartesian& q) const { return Vector2Cartesian{ func1(q), func2(q) }; }
		
		/// @brief Apply inverse rotation transformation
		/// @param q Rotated Cartesian coordinates
		/// @return Original Cartesian coordinates
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

	/// @brief Oblique to Cartesian 2D coordinate transformation
	/// 
	/// Transforms from oblique coordinates (where Y-axis is at angle ? to X-axis)
	/// to standard Cartesian coordinates (orthogonal axes).
	/// Forward: x = q1 + q2�cos(?), y = q2�sin(?)
	/// Useful for non-orthogonal coordinate systems.
	class CoordTransfObliqueToCartesian2D : public CoordTransfWithInverse<Vector2Cartesian, Vector2Cartesian, 2>
	{
	private:
		Real    _angle;

		const ScalarFunctionFromStdFunc<2> _f1;
		const ScalarFunctionFromStdFunc<2> _f2;

		const ScalarFunctionFromStdFunc<2> _fInverse1;
		const ScalarFunctionFromStdFunc<2> _fInverse2;

	public:
		/// @brief Constructor for oblique to Cartesian transformation
		/// @param inAngle Angle between oblique Y-axis and X-axis (in radians)
		CoordTransfObliqueToCartesian2D(Real inAngle)
			: _angle(inAngle),
			_f1([this](const VectorN<Real, 2>& q) { return func1(q); }),
			_f2([this](const VectorN<Real, 2>& q) { return func2(q); }),
			_fInverse1([this](const VectorN<Real, 2>& q) { return funcInverse1(q); }),
			_fInverse2([this](const VectorN<Real, 2>& q) { return funcInverse2(q); })
		{
			if (std::abs(std::sin(inAngle)) < Constants::Eps)
				throw MML::ArgumentError("CoordTransfObliqueToCartesian2D: angle must not be a multiple of pi (sin(angle) ~ 0)");
		}

		/// @brief Transform oblique x-coordinate: x = q1 + q2�cos(?)
		/// @param q Oblique coordinates (q1, q2)
		/// @return Cartesian x coordinate
		Real func1(const VectorN<Real, 2>& q) const { return q[0] + q[1] * std::cos(_angle); }
		
		/// @brief Transform oblique y-coordinate: y = q2�sin(?)
		/// @param q Oblique coordinates (q1, q2)
		/// @return Cartesian y coordinate
		Real func2(const VectorN<Real, 2>& q) const { return q[1] * std::sin(_angle); }

		/// @brief Inverse transform y: q2 = y/sin(?)
		/// @param q Cartesian coordinates (x, y)
		/// @return Oblique q2 coordinate
		Real funcInverse2(const VectorN<Real, 2>& q) const { return q[1] / std::sin(_angle); }
		
		/// @brief Inverse transform x: q1 = x - q2�cos(?)
		/// @param q Cartesian coordinates (x, y)
		/// @return Oblique q1 coordinate
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


	/// @brief Cartesian to oblique 2D coordinate transformation
	/// 
	/// Inverse of CoordTransfObliqueToCartesian2D.
	/// Transforms from standard Cartesian coordinates to oblique coordinates
	/// where Y-axis makes angle ? with X-axis.
	class CoordTransfCartesianToOblique2D : public CoordTransfWithInverse<Vector2Cartesian, Vector2Cartesian, 2>
	{
	private:
		Real    _angle;

		const ScalarFunctionFromStdFunc<2> _f1;
		const ScalarFunctionFromStdFunc<2> _f2;

		const ScalarFunctionFromStdFunc<2> _fInverse1;
		const ScalarFunctionFromStdFunc<2> _fInverse2;

	public:
		/// @brief Constructor for Cartesian to oblique transformation
		/// @param inAngle Angle between oblique Y-axis and X-axis (in radians)
		CoordTransfCartesianToOblique2D(Real inAngle)
			: _angle(inAngle),
			_f1([this](const VectorN<Real, 2>& q) { return func1(q); }),
			_f2([this](const VectorN<Real, 2>& q) { return func2(q); }),
			_fInverse1([this](const VectorN<Real, 2>& q) { return funcInverse1(q); }),
			_fInverse2([this](const VectorN<Real, 2>& q) { return funcInverse2(q); })
		{
			if (std::abs(std::sin(inAngle)) < Constants::Eps)
				throw MML::ArgumentError("CoordTransfCartesianToOblique2D: angle must not be a multiple of pi (sin(angle) ~ 0)");
		}

		/// @brief Transform y to oblique: q2 = y/sin(?)
		/// @param q Cartesian coordinates (x, y)
		/// @return Oblique q2 coordinate
		Real func2(const VectorN<Real, 2>& q) const { return q[1] / std::sin(_angle); }
		
		/// @brief Transform x to oblique: q1 = x - q2�cos(?)
		/// @param q Cartesian coordinates (x, y)
		/// @return Oblique q1 coordinate
		Real func1(const VectorN<Real, 2>& q) const { return q[0] - func2(q) * std::cos(_angle); }

		/// @brief Inverse transform to Cartesian: x = q1 + q2�cos(?)
		/// @param q Oblique coordinates (q1, q2)
		/// @return Cartesian x coordinate
		Real funcInverse1(const VectorN<Real, 2>& q) const { return q[0] + q[1] * std::cos(_angle); }
		
		/// @brief Inverse transform to Cartesian: y = q2�sin(?)
		/// @param q Oblique coordinates (q1, q2)
		/// @return Cartesian y coordinate
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

///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Algebra/LieGroups/SO2.h                                             ///
///  Description: Special orthogonal group SO(2)                                      ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_ALGEBRA_SO2_H
#define MML_ALGEBRA_SO2_H

#include <mml/MMLBase.h>
#include <mml/base/Matrix/MatrixNM.h>
#include <mml/base/Vector/VectorN.h>

#include <cmath>

namespace MML
{
	namespace Algebra
	{
		/// @brief Active right-handed planar rotation represented by a normalized angle.
		class SO2
		{
			Real _angle = REAL(0.0);

		public:
			SO2() = default;
			explicit SO2(Real angle) : _angle(normalizeAngle(angle)) { }

			static SO2 Identity() { return SO2(); }
			static SO2 Exp(Real angularDisplacement) { return SO2(angularDisplacement); }

			Real angle() const noexcept { return _angle; }
			Real log() const noexcept { return _angle; }

			MatrixNM<Real, 2, 2> matrix() const
			{
				const Real cosine = std::cos(_angle);
				const Real sine = std::sin(_angle);
				return {{cosine, -sine}, {sine, cosine}};
			}

			SO2 compose(const SO2& right) const { return SO2(_angle + right._angle); }
			SO2 inverse() const { return SO2(-_angle); }
			VectorN<Real, 2> apply(const VectorN<Real, 2>& vector) const { return matrix() * vector; }

			Real geodesic_distance(const SO2& other) const
			{
				return std::abs(normalizeAngle(other._angle - _angle));
			}

			static SO2 Interpolate(const SO2& from, const SO2& to, Real parameter)
			{
				const Real displacement = normalizeAngle(to._angle - from._angle);
				return SO2(from._angle + parameter * displacement);
			}
		};
	}
}

#endif // MML_ALGEBRA_SO2_H
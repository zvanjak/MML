///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Algebra/LieGroups/SE2.h                                             ///
///  Description: Planar rigid-motion group SE(2)                                     ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_ALGEBRA_SE2_H
#define MML_ALGEBRA_SE2_H

#include <mml/MMLExceptions.h>
#include <mml/base/Algebra/LieGroups/SO2.h>
#include <mml/base/DifferentialGeometry/Point.h>

#include <cmath>

namespace MML::Algebra
{
	class SE2
	{
		SO2 _rotation;
		VectorN<Real, 2> _translation;

	public:
		SE2() = default;
		SE2(const SO2& rotation, const VectorN<Real, 2>& translation)
			: _rotation(rotation), _translation(translation) { }

		static SE2 Identity() { return SE2(); }
		const SO2& rotation() const noexcept { return _rotation; }
		const VectorN<Real, 2>& translation() const noexcept { return _translation; }

		VectorN<Real, 2> apply_point(const VectorN<Real, 2>& point) const
		{
			return _rotation.apply(point) + _translation;
		}
		VectorN<Real, 2> apply_vector(const VectorN<Real, 2>& vector) const { return _rotation.apply(vector); }

		template<class FromFrame, class ToFrame>
		Point<Real, 2, ToFrame> apply_point(const Point<Real, 2, FromFrame>& point) const
		{
			return Point<Real, 2, ToFrame>(apply_point(point.coordinates()));
		}
		template<class FromFrame, class ToFrame>
		TangentVector<2, ToFrame> apply_vector(const TangentVector<2, FromFrame>& vector) const
		{
			return TangentVector<2, ToFrame>(apply_vector(vector.components()));
		}

		SE2 compose(const SE2& right) const
		{
			return SE2(_rotation.compose(right._rotation), _rotation.apply(right._translation) + _translation);
		}
		SE2 inverse() const
		{
			const SO2 inverseRotation = _rotation.inverse();
			return SE2(inverseRotation, inverseRotation.apply(-_translation));
		}

		MatrixNM<Real, 3, 3> homogeneous_matrix() const
		{
			MatrixNM<Real, 3, 3> result = MatrixNM<Real, 3, 3>::Identity();
			const auto rotationMatrix = _rotation.matrix();
			for (int row = 0; row < 2; row++) {
				for (int column = 0; column < 2; column++)
					result(row, column) = rotationMatrix(row, column);
				result(row, 2) = _translation[row];
			}
			return result;
		}

		static SE2 FromHomogeneousMatrix(const MatrixNM<Real, 3, 3>& matrix,
			Real tolerance = PrecisionValues<Real>::OrthogonalityTolerance)
		{
			if (std::abs(matrix(2, 0)) > tolerance || std::abs(matrix(2, 1)) > tolerance ||
				std::abs(matrix(2, 2) - REAL(1.0)) > tolerance)
				throw DomainError("SE2 homogeneous matrix must have bottom row [0,0,1]");
			const SO2 rotation(std::atan2(matrix(1, 0), matrix(0, 0)));
			if (!rotation.matrix().IsEqualTo({{matrix(0, 0), matrix(0, 1)},
				{matrix(1, 0), matrix(1, 1)}}, tolerance))
				throw DomainError("SE2 homogeneous matrix has an invalid rotation block");
			return SE2(rotation, {matrix(0, 2), matrix(1, 2)});
		}
	};
}

#endif // MML_ALGEBRA_SE2_H
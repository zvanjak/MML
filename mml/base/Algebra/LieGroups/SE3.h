///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Algebra/LieGroups/SE3.h                                             ///
///  Description: Spatial rigid-motion group SE(3)                                    ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_ALGEBRA_SE3_H
#define MML_ALGEBRA_SE3_H

#include <mml/MMLExceptions.h>
#include <mml/base/Algebra/LieGroups/SO3.h>
#include <mml/base/DifferentialGeometry/Point.h>

#include <cmath>

namespace MML::Algebra
{
	class SE3
	{
		SO3 _rotation;
		VectorN<Real, 3> _translation;

	public:
		SE3() = default;
		SE3(const SO3& rotation, const VectorN<Real, 3>& translation)
			: _rotation(rotation), _translation(translation) { }

		static SE3 Identity() { return SE3(); }
		const SO3& rotation() const noexcept { return _rotation; }
		const VectorN<Real, 3>& translation() const noexcept { return _translation; }

		VectorN<Real, 3> apply_point(const VectorN<Real, 3>& point) const
		{
			return _rotation.apply(point) + _translation;
		}
		VectorN<Real, 3> apply_vector(const VectorN<Real, 3>& vector) const { return _rotation.apply(vector); }

		template<class FromFrame, class ToFrame>
		Point<Real, 3, ToFrame> apply_point(const Point<Real, 3, FromFrame>& point) const
		{
			return Point<Real, 3, ToFrame>(apply_point(point.coordinates()));
		}
		template<class FromFrame, class ToFrame>
		TangentVector<3, ToFrame> apply_vector(const TangentVector<3, FromFrame>& vector) const
		{
			return TangentVector<3, ToFrame>(apply_vector(vector.components()));
		}

		SE3 compose(const SE3& right) const
		{
			return SE3(_rotation.compose(right._rotation), _rotation.apply(right._translation) + _translation);
		}
		SE3 inverse() const
		{
			const SO3 inverseRotation = _rotation.inverse();
			return SE3(inverseRotation, inverseRotation.apply(-_translation));
		}

		MatrixNM<Real, 4, 4> homogeneous_matrix() const
		{
			MatrixNM<Real, 4, 4> result = MatrixNM<Real, 4, 4>::Identity();
			const auto rotationMatrix = _rotation.matrix();
			for (int row = 0; row < 3; row++) {
				for (int column = 0; column < 3; column++)
					result(row, column) = rotationMatrix(row, column);
				result(row, 3) = _translation[row];
			}
			return result;
		}

		static SE3 FromHomogeneousMatrix(const MatrixNM<Real, 4, 4>& matrix,
			Real tolerance = PrecisionValues<Real>::OrthogonalityTolerance)
		{
			for (int column = 0; column < 3; column++)
				if (std::abs(matrix(3, column)) > tolerance)
					throw DomainError("SE3 homogeneous matrix must have bottom row [0,0,0,1]");
			if (std::abs(matrix(3, 3) - REAL(1.0)) > tolerance)
				throw DomainError("SE3 homogeneous matrix must have bottom row [0,0,0,1]");
			MatrixNM<Real, 3, 3> rotationMatrix;
			for (int row = 0; row < 3; row++)
				for (int column = 0; column < 3; column++)
					rotationMatrix(row, column) = matrix(row, column);
			return SE3(SO3::FromMatrix(rotationMatrix, tolerance),
				{matrix(0, 3), matrix(1, 3), matrix(2, 3)});
		}
	};
}

#endif // MML_ALGEBRA_SE3_H
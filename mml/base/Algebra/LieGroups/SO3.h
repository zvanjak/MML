///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Algebra/LieGroups/SO3.h                                             ///
///  Description: Quaternion-backed special orthogonal group SO(3)                    ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_ALGEBRA_SO3_H
#define MML_ALGEBRA_SO3_H

#include <mml/MMLBase.h>
#include <mml/MMLExceptions.h>
#include <mml/base/Quaternions.h>
#include <mml/base/Matrix/MatrixNM.h>
#include <mml/base/Vector/VectorN.h>

#include <algorithm>
#include <cmath>

namespace MML
{
	namespace Algebra
	{
		/// @brief Active right-handed 3D rotation with a canonical unit quaternion.
		class SO3
		{
			Quaternion _quaternion;

			static Quaternion Canonicalize(Quaternion quaternion)
			{
				quaternion.Normalize();
				const bool negativeHemisphere = quaternion.w() < REAL(0.0);
				const bool ambiguousNegative = std::abs(quaternion.w()) <= PrecisionValues<Real>::DefaultTolerance &&
					(quaternion.x() < REAL(0.0) ||
					 (quaternion.x() == REAL(0.0) && quaternion.y() < REAL(0.0)) ||
					 (quaternion.x() == REAL(0.0) && quaternion.y() == REAL(0.0) && quaternion.z() < REAL(0.0)));
				return negativeHemisphere || ambiguousNegative ? -quaternion : quaternion;
			}

			explicit SO3(const Quaternion& quaternion) : _quaternion(Canonicalize(quaternion)) { }

			static Real Determinant(const MatrixNM<Real, 3, 3>& matrix)
			{
				return matrix(0, 0) * (matrix(1, 1) * matrix(2, 2) - matrix(1, 2) * matrix(2, 1))
					- matrix(0, 1) * (matrix(1, 0) * matrix(2, 2) - matrix(1, 2) * matrix(2, 0))
					+ matrix(0, 2) * (matrix(1, 0) * matrix(2, 1) - matrix(1, 1) * matrix(2, 0));
			}

			static Real Dot(const VectorN<Real, 3>& left, const VectorN<Real, 3>& right)
			{
				return left[0] * right[0] + left[1] * right[1] + left[2] * right[2];
			}

			static VectorN<Real, 3> Cross(const VectorN<Real, 3>& left, const VectorN<Real, 3>& right)
			{
				return {
					left[1] * right[2] - left[2] * right[1],
					left[2] * right[0] - left[0] * right[2],
					left[0] * right[1] - left[1] * right[0]
				};
			}

			static Real Norm(const VectorN<Real, 3>& vector)
			{
				return std::sqrt(Dot(vector, vector));
			}

		public:
			SO3() = default;

			static SO3 Identity() { return SO3(); }

			static bool IsRotationMatrix(const MatrixNM<Real, 3, 3>& matrix,
				Real tolerance = PrecisionValues<Real>::OrthogonalityTolerance)
			{
				const auto orthogonality = matrix.transpose() * matrix;
				return orthogonality.IsEqualTo(MatrixNM<Real, 3, 3>::Identity(), tolerance) &&
					std::abs(Determinant(matrix) - REAL(1.0)) <= tolerance;
			}

			static SO3 FromMatrix(const MatrixNM<Real, 3, 3>& matrix,
				Real tolerance = PrecisionValues<Real>::OrthogonalityTolerance)
			{
				if (!IsRotationMatrix(matrix, tolerance))
					throw DomainError("SO3::FromMatrix requires an orthogonal matrix with determinant one");
				return SO3(Quaternion::FromRotationMatrix(matrix));
			}

			static SO3 FromAxisAngle(const VectorN<Real, 3>& axis, Real angle)
			{
				const Real norm = Norm(axis);
				if (norm <= PrecisionValues<Real>::NumericalZeroThreshold)
					throw ArgumentError("SO3 axis must be nonzero");
				const VectorN<Real, 3> unitAxis = axis / norm;
				return SO3(Quaternion::FromAxisAngle(
					Vec3Cart(unitAxis[0], unitAxis[1], unitAxis[2]), angle));
			}

			static MatrixNM<Real, 3, 3> Hat(const VectorN<Real, 3>& omega)
			{
				return {{0, -omega[2], omega[1]},
					{omega[2], 0, -omega[0]},
					{-omega[1], omega[0], 0}};
			}

			static VectorN<Real, 3> Vee(const MatrixNM<Real, 3, 3>& skew,
				Real tolerance = PrecisionValues<Real>::DefaultTolerance)
			{
				if (!skew.IsEqualTo(-skew.transpose(), tolerance))
					throw DomainError("SO3::Vee requires a skew-symmetric matrix");
				return {skew(2, 1), skew(0, 2), skew(1, 0)};
			}

			static SO3 Exp(const VectorN<Real, 3>& rotationVector)
			{
				const Real angle = Norm(rotationVector);
				const Real angleSquared = angle * angle;
				const Real halfAngle = angle / REAL(2.0);
				const Real vectorScale = angle <= REAL(1e-8)
					? REAL(0.5) - angleSquared / REAL(48.0)
					: std::sin(halfAngle) / angle;
				return SO3(Quaternion(std::cos(halfAngle),
					rotationVector[0] * vectorScale,
					rotationVector[1] * vectorScale,
					rotationVector[2] * vectorScale));
			}

			MatrixNM<Real, 3, 3> matrix() const { return _quaternion.ToRotationMatrix(); }

			VectorN<Real, 3> log() const
			{
				const Real vectorNorm = std::sqrt(_quaternion.x() * _quaternion.x() +
					_quaternion.y() * _quaternion.y() + _quaternion.z() * _quaternion.z());
				if (vectorNorm <= REAL(1e-12))
					return {REAL(2.0) * _quaternion.x(), REAL(2.0) * _quaternion.y(), REAL(2.0) * _quaternion.z()};
				const Real angle = REAL(2.0) * std::atan2(vectorNorm,
					std::clamp(_quaternion.w(), REAL(-1.0), REAL(1.0)));
				const Real scale = angle / vectorNorm;
				return {_quaternion.x() * scale, _quaternion.y() * scale, _quaternion.z() * scale};
			}

			SO3 compose(const SO3& right) const { return SO3(_quaternion * right._quaternion); }
			SO3 inverse() const { return SO3(_quaternion.Conjugate()); }
			VectorN<Real, 3> apply(const VectorN<Real, 3>& vector) const { return matrix() * vector; }

			Real geodesic_distance(const SO3& other) const
			{
				return Norm(inverse().compose(other).log());
			}

			static SO3 Slerp(const SO3& from, const SO3& to, Real parameter)
			{
				return SO3(Quaternion::Slerp(from._quaternion, to._quaternion, parameter));
			}

			static SO3 Project(const MatrixNM<Real, 3, 3>& matrix)
			{
				VectorN<Real, 3> first{matrix(0, 0), matrix(1, 0), matrix(2, 0)};
				VectorN<Real, 3> second{matrix(0, 1), matrix(1, 1), matrix(2, 1)};
				const VectorN<Real, 3> originalThird{matrix(0, 2), matrix(1, 2), matrix(2, 2)};
				const Real firstNorm = Norm(first);
				if (firstNorm <= PrecisionValues<Real>::NumericalZeroThreshold)
					throw DomainError("SO3::Project requires a full-rank matrix");
				first = first / firstNorm;
				second = second - first * Dot(first, second);
				const Real secondNorm = Norm(second);
				if (secondNorm <= PrecisionValues<Real>::NumericalZeroThreshold)
					throw DomainError("SO3::Project requires linearly independent columns");
				second = second / secondNorm;
				VectorN<Real, 3> third = Cross(first, second);
				if (Dot(third, originalThird) < REAL(0.0)) {
					second = -second;
					third = -third;
				}

				MatrixNM<Real, 3, 3> projected;
				for (int row = 0; row < 3; row++) {
					projected(row, 0) = first[row];
					projected(row, 1) = second[row];
					projected(row, 2) = third[row];
				}
				return FromMatrix(projected, REAL(100.0) * PrecisionValues<Real>::OrthogonalityTolerance);
			}

			const Quaternion& quaternion() const noexcept { return _quaternion; }
		};
	}
}

#endif // MML_ALGEBRA_SO3_H
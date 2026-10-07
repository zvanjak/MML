///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        VectorSpaces/Basis.h                                                ///
///  Description: Fixed-size vector-space basis and basis-change utilities            ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_VECTOR_SPACES_BASIS_H
#define MML_VECTOR_SPACES_BASIS_H

#include <mml/MMLBase.h>

#include <mml/base/Matrix/MatrixNM.h>
#include <mml/base/Vector/VectorN.h>
#include <mml/base/Algebra/Algorithms/FieldMatrixAlgorithms.h>

#include <array>
#include <initializer_list>
#include <type_traits>

namespace MML::VectorSpaces
{
	template<class Scalar, int N>
	class Basis
	{
		MatrixNM<Scalar, N, N> _matrix;
		MatrixNM<Scalar, N, N> _inverse;
		Scalar _determinant;

		void validate(Real tolerance)
		{
			_determinant = Algebra::FieldMatrixDeterminant(_matrix);
			if (std::abs(_determinant) <= tolerance)
				throw SingularMatrixError("Basis vectors must be linearly independent");
			_inverse = _matrix.GetInverse();
		}

	public:
		using scalar_type = Scalar;
		static constexpr int Dimension = N;

		Basis() : _matrix(MatrixNM<Scalar, N, N>::Identity()), _inverse(MatrixNM<Scalar, N, N>::Identity()), _determinant(Scalar{ 1 }) { }

		explicit Basis(const MatrixNM<Scalar, N, N>& columns, Real tolerance = Defaults::MatrixIsEqualTolerance)
			: _matrix(columns), _inverse(), _determinant()
		{
			validate(tolerance);
		}

		explicit Basis(const std::array<VectorN<Scalar, N>, N>& vectors, Real tolerance = Defaults::MatrixIsEqualTolerance)
			: _matrix(), _inverse(), _determinant()
		{
			for (int col = 0; col < N; col++)
				for (int row = 0; row < N; row++)
					_matrix(row, col) = vectors[col][row];
			validate(tolerance);
		}

		Basis(std::initializer_list<VectorN<Scalar, N>> vectors, Real tolerance = Defaults::MatrixIsEqualTolerance)
			: _matrix(), _inverse(), _determinant()
		{
			if (static_cast<int>(vectors.size()) != N)
				throw VectorDimensionError("Basis initializer list must contain exactly N vectors", N, static_cast<int>(vectors.size()));

			int col = 0;
			for (const auto& vector : vectors) {
				for (int row = 0; row < N; row++)
					_matrix(row, col) = vector[row];
				col++;
			}
			validate(tolerance);
		}

		static Basis Standard()
		{
			return Basis();
		}

		const MatrixNM<Scalar, N, N>& matrix() const noexcept { return _matrix; }
		const MatrixNM<Scalar, N, N>& inverseMatrix() const noexcept { return _inverse; }
		Scalar determinant() const noexcept { return _determinant; }

		int orientationSign(Real tolerance = Defaults::MatrixIsEqualTolerance) const requires MMLArithmetic<Scalar>
		{
			if (_determinant > Scalar(tolerance))
				return 1;
			if (_determinant < -Scalar(tolerance))
				return -1;
			return 0;
		}

		VectorN<Scalar, N> coordinatesInStandard(const VectorN<Scalar, N>& coordsInBasis) const
		{
			return _matrix * coordsInBasis;
		}

		VectorN<Scalar, N> coordinatesFromStandard(const VectorN<Scalar, N>& vector) const
		{
			return _inverse * vector;
		}
	};

	template<class Scalar, int N>
	VectorN<Scalar, N> ChangeBasis(const VectorN<Scalar, N>& coords,
		const Basis<Scalar, N>& fromBasis,
		const Basis<Scalar, N>& toBasis)
	{
		return toBasis.coordinatesFromStandard(fromBasis.coordinatesInStandard(coords));
	}
}

#endif // MML_VECTOR_SPACES_BASIS_H
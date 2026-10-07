///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Properties.h                                                        ///
///  Description: Matrix shape and algebraic property predicates                     ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_MATRIX_ALG_PROPERTIES_H
#define MML_MATRIX_ALG_PROPERTIES_H

#include <mml/algorithms/MatrixAnalysisTypes.h>

#include <algorithm>
#include <cmath>
#include <complex>

namespace MML::MatrixAlg
{
	namespace Detail
	{
		template<MMLScalar Scalar>
		constexpr MatrixComparisonTolerance<MatrixMagnitude<Scalar>> DefaultDiagonalTolerance()
		{
			using Magnitude = MatrixMagnitude<Scalar>;
			constexpr Magnitude tolerance = PrecisionValues<Magnitude>::IsMatrixDiagonalTolerance;
			return {tolerance, tolerance};
		}

		template<MMLScalar Scalar>
		constexpr MatrixComparisonTolerance<MatrixMagnitude<Scalar>> DefaultSymmetryTolerance()
		{
			using Magnitude = MatrixMagnitude<Scalar>;
			constexpr Magnitude tolerance = PrecisionValues<Magnitude>::IsMatrixSymmetricTolerance;
			return {tolerance, tolerance};
		}

		template<MMLScalar Scalar>
		MatrixMagnitude<Scalar> EntryScale(const Matrix<Scalar>& matrix)
		{
			using Magnitude = MatrixMagnitude<Scalar>;
			Magnitude scale{1};
			for (int row = 0; row < matrix.rows(); ++row)
				for (int col = 0; col < matrix.cols(); ++col)
					scale = std::max(scale, static_cast<Magnitude>(std::abs(matrix(row, col))));
			return scale;
		}

		template<MMLReal Magnitude>
		void ValidateTolerance(const MatrixComparisonTolerance<Magnitude>& tolerance)
		{
			if (tolerance.absolute < Magnitude{} || tolerance.relative < Magnitude{})
				throw DomainError("Matrix comparison tolerances must be nonnegative");
		}

		template<MMLReal Magnitude>
		bool IsApproximatelyZero(Magnitude value, Magnitude scale,
								 const MatrixComparisonTolerance<Magnitude>& tolerance)
		{
			return value <= tolerance.absolute + tolerance.relative * scale;
		}

		template<MMLScalar Scalar>
		Scalar Conjugate(const Scalar& value)
		{
			if constexpr (MMLComplex<Scalar>)
				return std::conj(value);
			else
				return value;
		}
	}

	template<MMLScalar Scalar>
	int Rows(const Matrix<Scalar>& matrix) noexcept { return matrix.rows(); }

	template<MMLScalar Scalar>
	int Cols(const Matrix<Scalar>& matrix) noexcept { return matrix.cols(); }

	template<MMLScalar Scalar>
	bool IsSquare(const Matrix<Scalar>& matrix) noexcept { return matrix.rows() == matrix.cols(); }

	template<MMLScalar Scalar>
	bool IsTall(const Matrix<Scalar>& matrix) noexcept { return matrix.rows() > matrix.cols(); }

	template<MMLScalar Scalar>
	bool IsWide(const Matrix<Scalar>& matrix) noexcept { return matrix.rows() < matrix.cols(); }

	template<MMLScalar Scalar>
	bool IsUpperHessenberg(const Matrix<Scalar>& matrix,
						  MatrixComparisonTolerance<MatrixMagnitude<Scalar>> tolerance = Detail::DefaultDiagonalTolerance<Scalar>())
	{
		Detail::ValidateTolerance(tolerance);
		const auto scale = Detail::EntryScale(matrix);
		for (int row = 0; row < matrix.rows(); ++row)
			for (int col = 0; col < matrix.cols() && col + 1 < row; ++col)
				if (!Detail::IsApproximatelyZero(static_cast<MatrixMagnitude<Scalar>>(std::abs(matrix(row, col))), scale, tolerance))
					return false;
		return true;
	}

	template<MMLScalar Scalar>
	bool IsUpperTriangular(const Matrix<Scalar>& matrix,
						   MatrixComparisonTolerance<MatrixMagnitude<Scalar>> tolerance = Detail::DefaultDiagonalTolerance<Scalar>())
	{
		Detail::ValidateTolerance(tolerance);
		const auto scale = Detail::EntryScale(matrix);
		for (int row = 1; row < matrix.rows(); ++row)
			for (int col = 0; col < matrix.cols() && col < row; ++col)
				if (!Detail::IsApproximatelyZero(static_cast<MatrixMagnitude<Scalar>>(std::abs(matrix(row, col))), scale, tolerance))
					return false;
		return true;
	}

	template<MMLScalar Scalar>
	bool IsLowerTriangular(const Matrix<Scalar>& matrix,
						   MatrixComparisonTolerance<MatrixMagnitude<Scalar>> tolerance = Detail::DefaultDiagonalTolerance<Scalar>())
	{
		Detail::ValidateTolerance(tolerance);
		const auto scale = Detail::EntryScale(matrix);
		for (int row = 0; row < matrix.rows(); ++row)
			for (int col = row + 1; col < matrix.cols(); ++col)
				if (!Detail::IsApproximatelyZero(static_cast<MatrixMagnitude<Scalar>>(std::abs(matrix(row, col))), scale, tolerance))
					return false;
		return true;
	}

	template<MMLScalar Scalar>
	bool IsDiagonal(const Matrix<Scalar>& matrix,
					MatrixComparisonTolerance<MatrixMagnitude<Scalar>> tolerance = Detail::DefaultDiagonalTolerance<Scalar>())
	{
		Detail::ValidateTolerance(tolerance);
		const auto scale = Detail::EntryScale(matrix);
		for (int row = 0; row < matrix.rows(); ++row)
			for (int col = 0; col < matrix.cols(); ++col)
				if (row != col && !Detail::IsApproximatelyZero(static_cast<MatrixMagnitude<Scalar>>(std::abs(matrix(row, col))), scale, tolerance))
					return false;
		return true;
	}

	template<MMLScalar Scalar>
	bool IsSymmetric(const Matrix<Scalar>& matrix,
					 MatrixComparisonTolerance<MatrixMagnitude<Scalar>> tolerance = Detail::DefaultSymmetryTolerance<Scalar>())
	{
		Detail::ValidateTolerance(tolerance);
		if (!IsSquare(matrix))
			return false;
		const auto scale = Detail::EntryScale(matrix);
		for (int row = 0; row < matrix.rows(); ++row)
			for (int col = row + 1; col < matrix.cols(); ++col)
				if (!Detail::IsApproximatelyZero(static_cast<MatrixMagnitude<Scalar>>(std::abs(matrix(row, col) - matrix(col, row))), scale, tolerance))
					return false;
		return true;
	}

	template<MMLScalar Scalar>
	bool IsSkewSymmetric(const Matrix<Scalar>& matrix,
						 MatrixComparisonTolerance<MatrixMagnitude<Scalar>> tolerance = Detail::DefaultSymmetryTolerance<Scalar>())
	{
		Detail::ValidateTolerance(tolerance);
		if (!IsSquare(matrix))
			return false;
		const auto scale = Detail::EntryScale(matrix);
		for (int row = 0; row < matrix.rows(); ++row)
			for (int col = row; col < matrix.cols(); ++col)
				if (!Detail::IsApproximatelyZero(static_cast<MatrixMagnitude<Scalar>>(std::abs(matrix(row, col) + matrix(col, row))), scale, tolerance))
					return false;
		return true;
	}

	template<MMLScalar Scalar>
	bool IsHermitian(const Matrix<Scalar>& matrix,
					 MatrixComparisonTolerance<MatrixMagnitude<Scalar>> tolerance = Detail::DefaultSymmetryTolerance<Scalar>())
	{
		Detail::ValidateTolerance(tolerance);
		if (!IsSquare(matrix))
			return false;
		const auto scale = Detail::EntryScale(matrix);
		for (int row = 0; row < matrix.rows(); ++row)
			for (int col = row; col < matrix.cols(); ++col)
				if (!Detail::IsApproximatelyZero(static_cast<MatrixMagnitude<Scalar>>(std::abs(matrix(row, col) - Detail::Conjugate(matrix(col, row)))), scale, tolerance))
					return false;
		return true;
	}

	template<MMLScalar Scalar>
	bool IsSkewHermitian(const Matrix<Scalar>& matrix,
						 MatrixComparisonTolerance<MatrixMagnitude<Scalar>> tolerance = Detail::DefaultSymmetryTolerance<Scalar>())
	{
		Detail::ValidateTolerance(tolerance);
		if (!IsSquare(matrix))
			return false;
		const auto scale = Detail::EntryScale(matrix);
		for (int row = 0; row < matrix.rows(); ++row)
			for (int col = row; col < matrix.cols(); ++col)
				if (!Detail::IsApproximatelyZero(static_cast<MatrixMagnitude<Scalar>>(std::abs(matrix(row, col) + Detail::Conjugate(matrix(col, row)))), scale, tolerance))
					return false;
		return true;
	}

	template<MMLScalar Scalar>
	bool IsDiagonallyDominant(const Matrix<Scalar>& matrix,
							MatrixComparisonTolerance<MatrixMagnitude<Scalar>> tolerance = Detail::DefaultDiagonalTolerance<Scalar>())
	{
		Detail::ValidateTolerance(tolerance);
		if (!IsSquare(matrix) || matrix.rows() == 0)
			return false;
		const auto scale = Detail::EntryScale(matrix);
		const auto margin = tolerance.absolute + tolerance.relative * scale;
		for (int row = 0; row < matrix.rows(); ++row)
		{
			MatrixMagnitude<Scalar> offDiagonal{};
			for (int col = 0; col < matrix.cols(); ++col)
				if (col != row)
					offDiagonal += static_cast<MatrixMagnitude<Scalar>>(std::abs(matrix(row, col)));
			if (static_cast<MatrixMagnitude<Scalar>>(std::abs(matrix(row, row))) <= offDiagonal + margin)
				return false;
		}
		return true;
	}

	template<MMLScalar Scalar>
	bool IsOrthogonal(const Matrix<Scalar>& matrix,
					  MatrixComparisonTolerance<MatrixMagnitude<Scalar>> tolerance = Detail::DefaultSymmetryTolerance<Scalar>())
	{
		Detail::ValidateTolerance(tolerance);
		if (!IsSquare(matrix) || matrix.rows() == 0)
			return false;
		const int size = matrix.rows();
		const auto scale = Detail::EntryScale(matrix);
		for (int row = 0; row < size; ++row)
			for (int col = 0; col < size; ++col)
			{
				Scalar product{};
				for (int index = 0; index < size; ++index)
					product += matrix(index, row) * matrix(index, col);
				const Scalar expected = row == col ? Scalar{1} : Scalar{};
				if (!Detail::IsApproximatelyZero(static_cast<MatrixMagnitude<Scalar>>(std::abs(product - expected)), scale, tolerance))
					return false;
			}
		return true;
	}

	template<MMLScalar Scalar>
	bool IsUnitary(const Matrix<Scalar>& matrix,
					MatrixComparisonTolerance<MatrixMagnitude<Scalar>> tolerance = Detail::DefaultSymmetryTolerance<Scalar>())
	{
		Detail::ValidateTolerance(tolerance);
		if (!IsSquare(matrix) || matrix.rows() == 0)
			return false;
		const int size = matrix.rows();
		const auto scale = Detail::EntryScale(matrix);
		for (int row = 0; row < size; ++row)
			for (int col = 0; col < size; ++col)
			{
				Scalar product{};
				for (int index = 0; index < size; ++index)
					product += Detail::Conjugate(matrix(index, row)) * matrix(index, col);
				const Scalar expected = row == col ? Scalar{1} : Scalar{};
				if (!Detail::IsApproximatelyZero(static_cast<MatrixMagnitude<Scalar>>(std::abs(product - expected)), scale, tolerance))
					return false;
			}
		return true;
	}

	template<MMLScalar Scalar>
	bool IsZero(const Matrix<Scalar>& matrix,
					MatrixComparisonTolerance<MatrixMagnitude<Scalar>> tolerance = Detail::DefaultDiagonalTolerance<Scalar>())
	{
		Detail::ValidateTolerance(tolerance);
		const auto scale = Detail::EntryScale(matrix);
		for (int row = 0; row < matrix.rows(); ++row)
			for (int col = 0; col < matrix.cols(); ++col)
				if (!Detail::IsApproximatelyZero(static_cast<MatrixMagnitude<Scalar>>(std::abs(matrix(row, col))), scale, tolerance))
					return false;
		return true;
	}

	template<MMLScalar Scalar>
	bool IsNilpotent(const Matrix<Scalar>& matrix,
						 MatrixComparisonTolerance<MatrixMagnitude<Scalar>> tolerance = Detail::DefaultDiagonalTolerance<Scalar>())
	{
		if (!IsSquare(matrix) || matrix.rows() == 0)
			throw MatrixDimensionError("IsNilpotent - matrix must be non-empty and square", matrix.rows(), matrix.cols(), -1, -1);
		Matrix<Scalar> power = matrix;
		for (int exponent = 1; exponent <= matrix.rows(); ++exponent)
		{
			if (IsZero(power, tolerance))
				return true;
			if (exponent < matrix.rows())
				power = power * matrix;
		}
		return false;
	}

	template<MMLScalar Scalar>
	bool IsUnipotent(const Matrix<Scalar>& matrix,
						 MatrixComparisonTolerance<MatrixMagnitude<Scalar>> tolerance = Detail::DefaultDiagonalTolerance<Scalar>())
	{
		if (!IsSquare(matrix) || matrix.rows() == 0)
			throw MatrixDimensionError("IsUnipotent - matrix must be non-empty and square", matrix.rows(), matrix.cols(), -1, -1);
		return IsNilpotent(matrix - Matrix<Scalar>::Identity(matrix.rows()), tolerance);
	}
}

#endif // MML_MATRIX_ALG_PROPERTIES_H
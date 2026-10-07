///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Algebra.h                                                           ///
///  Description: Direct matrix algebra operations                                    ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_MATRIX_ALG_ALGEBRA_H
#define MML_MATRIX_ALG_ALGEBRA_H

#include <mml/algorithms/MatrixAlg/Measures.h>
#include <mml/core/LinAlgEqSolvers.h>

#include <algorithm>
#include <cmath>

namespace MML::MatrixAlg
{
	template<MMLScalar Scalar>
	Scalar Determinant(const Matrix<Scalar>& matrix)
	{
		if (!IsSquare(matrix) || matrix.rows() == 0)
			throw MatrixDimensionError("Determinant - matrix must be non-empty and square", matrix.rows(), matrix.cols(), -1, -1);
		try
		{
			LUSolver<Scalar> solver(matrix);
			return solver.det();
		}
		catch (const SingularMatrixError&)
		{
			return Scalar{};
		}
	}

	template<MMLScalar Scalar>
	Matrix<Scalar> Inverse(const Matrix<Scalar>& matrix)
	{
		if (!IsSquare(matrix) || matrix.rows() == 0)
			throw MatrixDimensionError("Inverse - matrix must be non-empty and square", matrix.rows(), matrix.cols(), -1, -1);
		LUSolver<Scalar> solver(matrix);
		Matrix<Scalar> inverse;
		solver.inverse(inverse);
		return inverse;
	}

	template<MMLScalar Scalar>
	int RankGaussian(const Matrix<Scalar>& matrix,
					 MatrixMagnitude<Scalar> tolerance = static_cast<MatrixMagnitude<Scalar>>(Defaults::RankAlgEPS))
	{
		using Magnitude = MatrixMagnitude<Scalar>;
		if (matrix.rows() == 0 || matrix.cols() == 0)
			throw MatrixDimensionError("RankGaussian - matrix must be non-empty", matrix.rows(), matrix.cols(), -1, -1);
		if (tolerance < Magnitude{})
			throw DomainError("RankGaussian tolerance must be nonnegative");

		Matrix<Scalar> echelon = matrix;
		int pivotRow = 0;
		for (int col = 0; col < echelon.cols() && pivotRow < echelon.rows(); ++col)
		{
			int bestRow = pivotRow;
			Magnitude bestMagnitude{};
			for (int row = pivotRow; row < echelon.rows(); ++row)
			{
				const auto magnitude = static_cast<Magnitude>(std::abs(echelon(row, col)));
				if (magnitude > bestMagnitude)
				{
					bestMagnitude = magnitude;
					bestRow = row;
				}
			}
			if (bestMagnitude <= tolerance)
				continue;
			if (bestRow != pivotRow)
				for (int index = 0; index < echelon.cols(); ++index)
					std::swap(echelon(pivotRow, index), echelon(bestRow, index));

			for (int row = pivotRow + 1; row < echelon.rows(); ++row)
			{
				const Scalar factor = echelon(row, col) / echelon(pivotRow, col);
				for (int index = col; index < echelon.cols(); ++index)
					echelon(row, index) -= factor * echelon(pivotRow, index);
			}
			++pivotRow;
		}
		return pivotRow;
	}

	template<MMLScalar Scalar>
	FaddeevLeVerrierResult<Scalar> FaddeevLeVerrier(const Matrix<Scalar>& matrix)
	{
		if (!IsSquare(matrix) || matrix.rows() == 0)
			throw MatrixDimensionError("FaddeevLeVerrier - matrix must be non-empty and square", matrix.rows(), matrix.cols(), -1, -1);

		const int size = matrix.rows();
		FaddeevLeVerrierResult<Scalar> result;
		result.characteristicCoefficients.Resize(size + 1);
		result.characteristicCoefficients[0] = Scalar{1};

		Matrix<Scalar> recurrence = Matrix<Scalar>::Identity(size);
		for (int coefficientIndex = 1; coefficientIndex <= size; ++coefficientIndex)
		{
			const Matrix<Scalar> product = matrix * recurrence;
			const Scalar coefficient = -Trace(product) / static_cast<Scalar>(coefficientIndex);
			result.characteristicCoefficients[coefficientIndex] = coefficient;

			if (coefficientIndex == size)
			{
				result.determinant = (size % 2 == 0) ? coefficient : -coefficient;
				if (Determinant(matrix) != Scalar{})
				{
					Matrix<Scalar> inverse(size, size);
					for (int row = 0; row < size; ++row)
						for (int col = 0; col < size; ++col)
							inverse(row, col) = -recurrence(row, col) / coefficient;
					result.inverse = std::move(inverse);
				}
			}
			else
			{
				recurrence = product;
				for (int index = 0; index < size; ++index)
					recurrence(index, index) += coefficient;
			}
		}
		return result;
	}
}

#endif // MML_MATRIX_ALG_ALGEBRA_H
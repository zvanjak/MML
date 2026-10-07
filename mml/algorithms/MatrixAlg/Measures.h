///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Measures.h                                                          ///
///  Description: Matrix scalar measures and norms                                    ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_MATRIX_ALG_MEASURES_H
#define MML_MATRIX_ALG_MEASURES_H

#include <mml/algorithms/MatrixAlg/Properties.h>

#include <algorithm>
#include <cmath>

namespace MML::MatrixAlg
{
	template<MMLScalar Scalar>
	MatrixMagnitude<Scalar> Sparsity(const Matrix<Scalar>& matrix,
									 MatrixMagnitude<Scalar> threshold = PrecisionValues<MatrixMagnitude<Scalar>>::EigenSolverZeroThreshold)
	{
		if (matrix.rows() == 0 || matrix.cols() == 0)
			throw MatrixDimensionError("Sparsity - matrix must be non-empty", matrix.rows(), matrix.cols(), -1, -1);
		if (threshold < MatrixMagnitude<Scalar>{})
			throw DomainError("Sparsity threshold must be nonnegative");
		int zeros = 0;
		for (int row = 0; row < matrix.rows(); ++row)
			for (int col = 0; col < matrix.cols(); ++col)
				if (std::abs(matrix(row, col)) <= threshold)
					++zeros;
		return static_cast<MatrixMagnitude<Scalar>>(zeros) /
			   static_cast<MatrixMagnitude<Scalar>>(matrix.rows() * matrix.cols());
	}

	template<MMLScalar Scalar>
	Scalar Trace(const Matrix<Scalar>& matrix)
	{
		if (!IsSquare(matrix) || matrix.rows() == 0)
			throw MatrixDimensionError("Trace - matrix must be non-empty and square", matrix.rows(), matrix.cols(), -1, -1);
		Scalar trace{};
		for (int index = 0; index < matrix.rows(); ++index)
			trace += matrix(index, index);
		return trace;
	}

	template<MMLScalar Scalar>
	MatrixMagnitude<Scalar> FrobeniusNorm(const Matrix<Scalar>& matrix)
	{
		if (matrix.rows() == 0 || matrix.cols() == 0)
			throw MatrixDimensionError("FrobeniusNorm - matrix must be non-empty", matrix.rows(), matrix.cols(), -1, -1);
		MatrixMagnitude<Scalar> norm{};
		for (int row = 0; row < matrix.rows(); ++row)
			for (int col = 0; col < matrix.cols(); ++col)
				norm = std::hypot(norm, static_cast<MatrixMagnitude<Scalar>>(std::abs(matrix(row, col))));
		return norm;
	}

	template<MMLScalar Scalar>
	MatrixMagnitude<Scalar> OneNorm(const Matrix<Scalar>& matrix)
	{
		if (matrix.rows() == 0 || matrix.cols() == 0)
			throw MatrixDimensionError("OneNorm - matrix must be non-empty", matrix.rows(), matrix.cols(), -1, -1);
		MatrixMagnitude<Scalar> maximum{};
		for (int col = 0; col < matrix.cols(); ++col)
		{
			MatrixMagnitude<Scalar> sum{};
			for (int row = 0; row < matrix.rows(); ++row)
				sum += static_cast<MatrixMagnitude<Scalar>>(std::abs(matrix(row, col)));
			maximum = std::max(maximum, sum);
		}
		return maximum;
	}

	template<MMLScalar Scalar>
	MatrixMagnitude<Scalar> InfinityNorm(const Matrix<Scalar>& matrix)
	{
		if (matrix.rows() == 0 || matrix.cols() == 0)
			throw MatrixDimensionError("InfinityNorm - matrix must be non-empty", matrix.rows(), matrix.cols(), -1, -1);
		MatrixMagnitude<Scalar> maximum{};
		for (int row = 0; row < matrix.rows(); ++row)
		{
			MatrixMagnitude<Scalar> sum{};
			for (int col = 0; col < matrix.cols(); ++col)
				sum += static_cast<MatrixMagnitude<Scalar>>(std::abs(matrix(row, col)));
			maximum = std::max(maximum, sum);
		}
		return maximum;
	}
}

#endif // MML_MATRIX_ALG_MEASURES_H
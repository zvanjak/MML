///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Decompositions.h                                                    ///
///  Description: Matrix decompositions and SVD construction helpers                  ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_MATRIX_ALG_DECOMPOSITIONS_H
#define MML_MATRIX_ALG_DECOMPOSITIONS_H

#include <mml/algorithms/MatrixAlg/Properties.h>
#include <mml/core/LinAlgEqSolvers.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>

namespace MML::MatrixAlg
{
	template<MMLScalar Scalar>
	LUDecomposition<Scalar> LUDecompose(const Matrix<Scalar>& matrix)
	{
		if (!IsSquare(matrix) || matrix.rows() == 0)
			throw MatrixDimensionError("LUDecompose - matrix must be non-empty and square", matrix.rows(), matrix.cols(), -1, -1);

		LUSolver<Scalar> solver(matrix);
		const Matrix<Scalar>& combined = solver.GetDecomposition();
		const int size = matrix.rows();
		LUDecomposition<Scalar> result;
		result.L.Resize(size, size);
		result.U.Resize(size, size);
		for (int row = 0; row < size; ++row)
			for (int col = 0; col < size; ++col)
			{
				if (row > col)
					result.L(row, col) = combined(row, col);
				else
				{
					result.L(row, col) = row == col ? Scalar{1} : Scalar{};
					result.U(row, col) = combined(row, col);
				}
			}

		result.permutation.resize(size);
		for (int row = 0; row < size; ++row)
			result.permutation[row] = row;
		const auto& pivots = solver.GetPivotIndices();
		for (int step = 0; step < size; ++step)
			std::swap(result.permutation[step], result.permutation[pivots[step]]);
		result.determinant = solver.det();
		return result;
	}

	template<MMLScalar Scalar>
	QRDecomposition<Scalar> QRDecompose(const Matrix<Scalar>& matrix)
	{
		if (matrix.rows() == 0 || matrix.cols() == 0)
			throw MatrixDimensionError("QRDecompose - matrix must be non-empty", matrix.rows(), matrix.cols(), -1, -1);
		if (matrix.rows() < matrix.cols())
			throw MatrixDimensionError("QRDecompose - economy QR requires rows >= columns", matrix.rows(), matrix.cols(), -1, -1);
		QRSolver<Scalar> solver(matrix);
		return {solver.GetQ(), solver.GetR()};
	}

	template<MMLScalar Scalar>
	CholeskyDecomposition<Scalar> CholeskyDecompose(
			const Matrix<Scalar>& matrix,
			MatrixMagnitude<Scalar> tolerance = PrecisionValues<MatrixMagnitude<Scalar>>::IsMatrixSymmetricTolerance)
	{
		CholeskySolver<Scalar> solver(matrix, tolerance);
		return {solver.L()};
	}

	namespace Detail
	{
		template<MMLReal Scalar>
		Matrix<Scalar> CompleteOrthonormalBasis(const Matrix<Scalar>& economy)
		{
			const int dimension = economy.rows();
			if (economy.cols() > dimension)
				throw MatrixDimensionError("CompleteOrthonormalBasis - too many supplied columns", economy.rows(), economy.cols(), -1, -1);

			Matrix<Scalar> full(dimension, dimension);
			for (int row = 0; row < dimension; ++row)
				for (int col = 0; col < economy.cols(); ++col)
					full(row, col) = economy(row, col);

			const Scalar completionTolerance = std::numeric_limits<Scalar>::epsilon() * static_cast<Scalar>(dimension * 16);
			int completedColumns = economy.cols();
			for (int basisIndex = 0; basisIndex < dimension && completedColumns < dimension; ++basisIndex)
			{
				Vector<Scalar> candidate(dimension);
				candidate[basisIndex] = Scalar{1};
				for (int pass = 0; pass < 2; ++pass)
					for (int col = 0; col < completedColumns; ++col)
					{
						Scalar projection{};
						for (int row = 0; row < dimension; ++row)
							projection += full(row, col) * candidate[row];
						for (int row = 0; row < dimension; ++row)
							candidate[row] -= projection * full(row, col);
					}

				const Scalar norm = candidate.NormL2();
				if (norm <= completionTolerance)
					continue;
				for (int row = 0; row < dimension; ++row)
					full(row, completedColumns) = candidate[row] / norm;
				++completedColumns;
			}

			if (completedColumns != dimension)
				throw MatrixNumericalError("CompleteOrthonormalBasis - failed to complete basis");
			return full;
		}

		template<MMLScalar Scalar>
		MatrixMagnitude<Scalar> ResolveSVDThreshold(const Vector<MatrixMagnitude<Scalar>>& singularValues,
											 int rows, int cols,
											 std::optional<MatrixMagnitude<Scalar>> requested)
		{
			using Magnitude = MatrixMagnitude<Scalar>;
			if (requested.has_value())
			{
				if (!std::isfinite(*requested) || *requested < Magnitude{})
					throw DomainError("SVD threshold must be finite and nonnegative");
				return *requested;
			}
			return static_cast<Magnitude>(std::max(rows, cols)) *
				   std::numeric_limits<Magnitude>::epsilon() * singularValues[0];
		}
	}

	template<MMLReal Scalar>
	SVDDecomposition<Scalar> SVDDecompose(const Matrix<Scalar>& matrix,
									  std::optional<MatrixMagnitude<Scalar>> threshold = std::nullopt)
	{
		if (matrix.rows() == 0 || matrix.cols() == 0)
			throw MatrixDimensionError("SVDDecompose - matrix must be non-empty", matrix.rows(), matrix.cols(), -1, -1);

		SVDDecomposition<Scalar> result;
		const int minDimension = std::min(matrix.rows(), matrix.cols());
		if (matrix.rows() >= matrix.cols())
		{
			SVDecompositionSolver<Scalar> solver(matrix);
			result.U = Detail::CompleteOrthonormalBasis(solver.getU());
			result.V = solver.getV();
			const auto values = solver.getW();
			result.singularValues.Resize(minDimension);
			for (int index = 0; index < minDimension; ++index)
				result.singularValues[index] = values[index];
		}
		else
		{
			SVDecompositionSolver<Scalar> solver(matrix.transpose());
			result.U = solver.getV();
			result.V = Detail::CompleteOrthonormalBasis(solver.getU());
			const auto values = solver.getW();
			result.singularValues.Resize(minDimension);
			for (int index = 0; index < minDimension; ++index)
				result.singularValues[index] = values[index];
		}

		result.threshold = Detail::ResolveSVDThreshold<Scalar>(result.singularValues, matrix.rows(), matrix.cols(), threshold);
		for (int index = 0; index < minDimension; ++index)
			if (result.singularValues[index] > result.threshold)
				++result.rank;
		return result;
	}

	inline SVDDecomposition<Complex> SVDDecompose(
			const Matrix<Complex>& matrix,
			std::optional<Real> threshold = std::nullopt)
	{
		const auto decomposition = ComplexSVDecompositionSolver::Decompose(matrix);
		SVDDecomposition<Complex> result;
		result.U = decomposition.U;
		result.V = decomposition.V;
		result.singularValues = decomposition.singularValues;
		result.threshold = Detail::ResolveSVDThreshold<Complex>(
			result.singularValues, matrix.rows(), matrix.cols(), threshold);
		for (int index = 0; index < result.singularValues.size(); ++index)
			if (result.singularValues[index] > result.threshold)
				++result.rank;
		return result;
	}
}

#endif // MML_MATRIX_ALG_DECOMPOSITIONS_H
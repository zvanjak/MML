///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Subspaces.h                                                         ///
///  Description: SVD-derived subspaces, conditioning, and stability                  ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_MATRIX_ALG_SUBSPACES_H
#define MML_MATRIX_ALG_SUBSPACES_H

#include <mml/algorithms/MatrixAlg/Decompositions.h>
#include <mml/algorithms/MatrixAlg/Measures.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>

namespace MML::MatrixAlg
{
	namespace Detail
	{
		template<MMLScalar Scalar>
		Matrix<Scalar> PseudoInverseFromSVD(const SVDDecomposition<Scalar>& svd, int rows, int cols)
		{
			Matrix<Scalar> result(cols, rows);
			for (int row = 0; row < result.rows(); ++row)
				for (int col = 0; col < result.cols(); ++col)
				{
					Scalar sum{};
					for (int index = 0; index < svd.singularValues.size(); ++index)
						if (svd.singularValues[index] > svd.threshold)
							sum += svd.V(row, index) * Conjugate(svd.U(col, index)) / svd.singularValues[index];
					result(row, col) = sum;
				}
			return result;
		}

		template<MMLScalar Scalar>
		MatrixAlg::FundamentalSubspaces<Scalar> FundamentalSubspacesFromSVD(
				const SVDDecomposition<Scalar>& svd, int rows, int cols)
		{
			MatrixAlg::FundamentalSubspaces<Scalar> result;
			result.rank = svd.rank;
			result.threshold = svd.threshold;
			result.columnSpace.Resize(rows, svd.rank);
			result.rowSpace.Resize(cols, svd.rank);
			result.nullSpace.Resize(cols, cols - svd.rank);
			result.leftNullSpace.Resize(rows, rows - svd.rank);

			for (int row = 0; row < rows; ++row)
			{
				for (int col = 0; col < svd.rank; ++col)
					result.columnSpace(row, col) = svd.U(row, col);
				for (int col = svd.rank; col < rows; ++col)
					result.leftNullSpace(row, col - svd.rank) = svd.U(row, col);
			}
			for (int row = 0; row < cols; ++row)
			{
				for (int col = 0; col < svd.rank; ++col)
					result.rowSpace(row, col) = svd.V(row, col);
				for (int col = svd.rank; col < cols; ++col)
					result.nullSpace(row, col - svd.rank) = svd.V(row, col);
			}
			return result;
		}

		template<MMLScalar Scalar>
		MatrixMagnitude<Scalar> ConditionNumberFromSVD(const SVDDecomposition<Scalar>& svd)
		{
			using Magnitude = MatrixMagnitude<Scalar>;
			if (svd.rank < svd.singularValues.size())
				return std::numeric_limits<Magnitude>::infinity();
			return svd.singularValues[0] / svd.singularValues[svd.singularValues.size() - 1];
		}

		template<MMLReal Magnitude>
		MatrixStability StabilityFromConditionNumber(Magnitude condition)
		{
			if (!std::isfinite(condition))
				return MatrixStability::Singular;
			if (condition < static_cast<Magnitude>(1e4))
				return MatrixStability::WellConditioned;
			if (condition < static_cast<Magnitude>(1e8))
				return MatrixStability::ModeratelyConditioned;
			return MatrixStability::IllConditioned;
		}

		template<MMLReal Magnitude>
		std::optional<int> DigitsLostFromConditionNumber(Magnitude condition)
		{
			if (!std::isfinite(condition))
				return std::nullopt;
			return static_cast<int>(std::floor(std::log10(std::max(condition, Magnitude{1}))));
		}
	}

	template<MMLScalar Scalar>
	Vector<MatrixMagnitude<Scalar>> SingularValues(const Matrix<Scalar>& matrix)
	{
		return SVDDecompose(matrix).singularValues;
	}

	template<MMLScalar Scalar>
	int Rank(const Matrix<Scalar>& matrix,
			 std::optional<MatrixMagnitude<Scalar>> threshold = std::nullopt)
	{
		return SVDDecompose(matrix, threshold).rank;
	}

	template<MMLScalar Scalar>
	int Nullity(const Matrix<Scalar>& matrix,
				std::optional<MatrixMagnitude<Scalar>> threshold = std::nullopt)
	{
		return matrix.cols() - SVDDecompose(matrix, threshold).rank;
	}

	template<MMLScalar Scalar>
	MatrixMagnitude<Scalar> ConditionNumber(const Matrix<Scalar>& matrix,
										std::optional<MatrixMagnitude<Scalar>> threshold = std::nullopt)
	{
		const auto svd = SVDDecompose(matrix, threshold);
		return Detail::ConditionNumberFromSVD(svd);
	}

	template<MMLScalar Scalar>
	Matrix<Scalar> PseudoInverse(const Matrix<Scalar>& matrix,
								 std::optional<MatrixMagnitude<Scalar>> threshold = std::nullopt)
	{
		const auto svd = SVDDecompose(matrix, threshold);
		return Detail::PseudoInverseFromSVD(svd, matrix.rows(), matrix.cols());
	}

	template<MMLScalar Scalar>
	MatrixAlg::FundamentalSubspaces<Scalar> FundamentalSubspacesOf(const Matrix<Scalar>& matrix,
														 std::optional<MatrixMagnitude<Scalar>> threshold = std::nullopt)
	{
		const auto svd = SVDDecompose(matrix, threshold);
		return Detail::FundamentalSubspacesFromSVD(svd, matrix.rows(), matrix.cols());
	}

	template<MMLScalar Scalar>
	Matrix<Scalar> NullSpace(const Matrix<Scalar>& matrix,
							 std::optional<MatrixMagnitude<Scalar>> threshold = std::nullopt)
	{
		return FundamentalSubspacesOf(matrix, threshold).nullSpace;
	}

	template<MMLScalar Scalar>
	Matrix<Scalar> ColumnSpace(const Matrix<Scalar>& matrix,
							   std::optional<MatrixMagnitude<Scalar>> threshold = std::nullopt)
	{
		return FundamentalSubspacesOf(matrix, threshold).columnSpace;
	}

	template<MMLScalar Scalar>
	Matrix<Scalar> RowSpace(const Matrix<Scalar>& matrix,
							std::optional<MatrixMagnitude<Scalar>> threshold = std::nullopt)
	{
		return FundamentalSubspacesOf(matrix, threshold).rowSpace;
	}

	template<MMLScalar Scalar>
	Matrix<Scalar> LeftNullSpace(const Matrix<Scalar>& matrix,
								 std::optional<MatrixMagnitude<Scalar>> threshold = std::nullopt)
	{
		return FundamentalSubspacesOf(matrix, threshold).leftNullSpace;
	}

	template<MMLScalar Scalar>
	MatrixMagnitude<Scalar> ConditionNumber1(const Matrix<Scalar>& matrix,
										 std::optional<MatrixMagnitude<Scalar>> threshold = std::nullopt)
	{
		using Magnitude = MatrixMagnitude<Scalar>;
		const auto svd = SVDDecompose(matrix, threshold);
		if (svd.rank < svd.singularValues.size())
			return std::numeric_limits<Magnitude>::infinity();
		return OneNorm(matrix) * OneNorm(Detail::PseudoInverseFromSVD(svd, matrix.rows(), matrix.cols()));
	}

	template<MMLScalar Scalar>
	MatrixMagnitude<Scalar> ConditionNumberInfinity(const Matrix<Scalar>& matrix,
												std::optional<MatrixMagnitude<Scalar>> threshold = std::nullopt)
	{
		using Magnitude = MatrixMagnitude<Scalar>;
		const auto svd = SVDDecompose(matrix, threshold);
		if (svd.rank < svd.singularValues.size())
			return std::numeric_limits<Magnitude>::infinity();
		return InfinityNorm(matrix) * InfinityNorm(Detail::PseudoInverseFromSVD(svd, matrix.rows(), matrix.cols()));
	}

	template<MMLScalar Scalar>
	MatrixStability AssessStability(const Matrix<Scalar>& matrix,
									 std::optional<MatrixMagnitude<Scalar>> threshold = std::nullopt)
	{
		return Detail::StabilityFromConditionNumber(ConditionNumber(matrix, threshold));
	}

	template<MMLScalar Scalar>
	std::optional<int> ExpectedDigitsLost(const Matrix<Scalar>& matrix,
										  std::optional<MatrixMagnitude<Scalar>> threshold = std::nullopt)
	{
		return Detail::DigitsLostFromConditionNumber(ConditionNumber(matrix, threshold));
	}
}

#endif // MML_MATRIX_ALG_SUBSPACES_H
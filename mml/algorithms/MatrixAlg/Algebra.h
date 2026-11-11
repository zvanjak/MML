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

#include <mml/algorithms/MatrixAlg/MatrixInversion.h>
#include <mml/algorithms/MatrixAlg/Measures.h>
#include <mml/core/LinAlgEqSolvers.h>

#include <algorithm>
#include <cmath>
#include <limits>

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
	int RankGaussian(const Matrix<Scalar>& matrix,
					 MatrixMagnitude<Scalar> tolerance = static_cast<MatrixMagnitude<Scalar>>(Precision::RankAlgEPS))
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

		using Magnitude = MatrixMagnitude<Scalar>;
		const auto requireFinite = [](const Scalar& value) {
			if (!std::isfinite(std::real(value)) || !std::isfinite(std::imag(value)))
				throw MatrixNumericalError("FaddeevLeVerrier - non-finite input, recurrence, or inverse arithmetic");
		};
		const auto scaleBinary = [](const Scalar& value, long long exponent) -> Scalar {
			if (value == Scalar{})
				return Scalar{};
			if (exponent > std::numeric_limits<int>::max() || exponent < std::numeric_limits<int>::min()) // tol-ok: integer scalbn exponent range, not a numerical tolerance.
				throw MatrixNumericalError("FaddeevLeVerrier - coefficient exponent is not representable");
			const auto scaleComponent = [exponent](Magnitude component) {
				const Magnitude scaled = std::scalbn(component, static_cast<int>(exponent));
				if (!std::isfinite(scaled) || (component != Magnitude{} && scaled == Magnitude{}))
					throw MatrixNumericalError("FaddeevLeVerrier - scaling overflow or underflow");
				return scaled;
			};
			if constexpr (MMLComplex<Scalar>)
				return Scalar{scaleComponent(value.real()), scaleComponent(value.imag())};
			else
				return scaleComponent(value);
		};

		const int size = matrix.rows();
		Magnitude largest{};
		for (int row = 0; row < size; ++row)
			for (int col = 0; col < size; ++col) {
				requireFinite(matrix(row, col));
				largest = std::max(largest, std::abs(std::real(matrix(row, col))));
				largest = std::max(largest, std::abs(std::imag(matrix(row, col))));
			}
		int scaleExponent = 0;
		if (largest != Magnitude{})
			std::frexp(largest, &scaleExponent);
		Matrix<Scalar> normalized(size, size);
		for (int row = 0; row < size; ++row)
			for (int col = 0; col < size; ++col) {
				normalized(row, col) = scaleBinary(matrix(row, col), -scaleExponent);
				if (scaleBinary(normalized(row, col), scaleExponent) != matrix(row, col))
					throw MatrixNumericalError("FaddeevLeVerrier - normalization loses input precision");
			}

		FaddeevLeVerrierResult<Scalar> result;
		result.characteristicCoefficients.resize(size + 1);
		result.characteristicCoefficients[0] = Scalar{1};

		// Newton identities on A / 2^e; the kth coefficient scales back by 2^(k*e).
		Matrix<Scalar> recurrence = Matrix<Scalar>::Identity(size);
		for (int coefficientIndex = 1; coefficientIndex <= size; ++coefficientIndex)
		{
			const Matrix<Scalar> product = normalized * recurrence;
			for (int row = 0; row < size; ++row)
				for (int col = 0; col < size; ++col)
					requireFinite(product(row, col));
			Scalar coefficient{};
			for (int index = 0; index < size; ++index)
				coefficient -= product(index, index) / ScalarDetail::ScalarFromInteger<Scalar>(coefficientIndex);
			requireFinite(coefficient);
			result.characteristicCoefficients[coefficientIndex] =
				scaleBinary(coefficient, static_cast<long long>(coefficientIndex) * scaleExponent);

			if (coefficientIndex == size)
			{
				const Scalar constant = result.characteristicCoefficients[coefficientIndex];
				result.determinant = (size % 2 == 0) ? constant : -constant;
			}
			else
			{
				recurrence = product;
				for (int index = 0; index < size; ++index) {
					recurrence(index, index) += coefficient;
					requireFinite(recurrence(index, index));
				}
			}
		}

		// Polynomial cancellation is not a reliable inverse formula, even when LU succeeds.
		Matrix<Scalar> inverse;
		try {
			inverse = Inverse(matrix);
		}
		catch (const SingularMatrixError&) {
			return result;
		}
		for (int row = 0; row < size; ++row)
			for (int col = 0; col < size; ++col)
				requireFinite(inverse(row, col));
		const Matrix<Scalar> leftProduct = matrix * inverse;
		const Matrix<Scalar> rightProduct = inverse * matrix;
		const Magnitude tolerance = PrecisionValues<Magnitude>::MatrixIsEqualTolerance * static_cast<Magnitude>(size);
		for (int row = 0; row < size; ++row)
			for (int col = 0; col < size; ++col) {
				const Scalar expected = (row == col) ? Scalar{1} : Scalar{};
				requireFinite(leftProduct(row, col));
				requireFinite(rightProduct(row, col));
				if (std::abs(leftProduct(row, col) - expected) > tolerance ||
					std::abs(rightProduct(row, col) - expected) > tolerance)
					throw MatrixNumericalError("FaddeevLeVerrier - LU inverse failed two-sided identity validation");
			}
		result.inverse = std::move(inverse);
		return result;
	}
}

#endif // MML_MATRIX_ALG_ALGEBRA_H
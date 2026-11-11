///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        LinAlgComplexSVD.h                                                  ///
///  Description: One-sided Jacobi SVD for dense complex matrices                     ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_LINEAR_ALG_COMPLEX_SVD_H
#define MML_LINEAR_ALG_COMPLEX_SVD_H

#include <mml/MMLBase.h>
#include <mml/base/Matrix/Matrix.h>
#include <mml/base/Vector/Vector.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <utility>

namespace MML
{
	/// Dense complex SVD computed directly from A by one-sided Jacobi rotations.
	/// Returns full unitary factors A = U diag(sigma) V* for any nonempty m x n matrix.
	class ComplexSVDecompositionSolver
	{
	public:
		struct Result
		{
			Matrix<Complex> U;
			Vector<Real> singularValues;
			Matrix<Complex> V;
			int sweeps = 0;
			Real maxCorrelation = static_cast<MML::Real>(0.0);
		};

		static Result Decompose(
				const Matrix<Complex>& matrix,
				Real tolerance = PrecisionValues<Real>::EigenSolverConvergenceTolerance,
				int maxSweeps = 100)
		{
			if (matrix.rows() == 0 || matrix.cols() == 0)
				throw MatrixDimensionError("ComplexSVDecompositionSolver - matrix must be non-empty",
					matrix.rows(), matrix.cols(), -1, -1);
			if (!std::isfinite(tolerance) || tolerance < static_cast<MML::Real>(0.0))
				throw DomainError("ComplexSVDecompositionSolver - tolerance must be finite and nonnegative");
			if (maxSweeps <= 0)
				throw DomainError("ComplexSVDecompositionSolver - maxSweeps must be positive");

			if (matrix.rows() >= matrix.cols())
				return DecomposeTall(matrix, tolerance, maxSweeps);

			Result adjointResult = DecomposeTall(Adjoint(matrix), tolerance, maxSweeps);
			Result result;
			result.U = std::move(adjointResult.V);
			result.V = std::move(adjointResult.U);
			result.singularValues = std::move(adjointResult.singularValues);
			result.sweeps = adjointResult.sweeps;
			result.maxCorrelation = adjointResult.maxCorrelation;
			return result;
		}

	private:
		static Matrix<Complex> Adjoint(const Matrix<Complex>& matrix)
		{
			Matrix<Complex> result(matrix.cols(), matrix.rows());
			for (int row = 0; row < matrix.rows(); ++row)
				for (int col = 0; col < matrix.cols(); ++col)
					result(col, row) = std::conj(matrix(row, col));
			return result;
		}

		static Complex NormalizedColumnInnerProduct(const Matrix<Complex>& matrix,
			int left, int right, Real leftNorm, Real rightNorm)
		{
			Complex result{};
			for (int row = 0; row < matrix.rows(); ++row)
				result += std::conj(matrix(row, left) / leftNorm) * (matrix(row, right) / rightNorm);
			return result;
		}

		static Real ColumnNorm(const Matrix<Complex>& matrix, int col)
		{
			Real result = static_cast<MML::Real>(0.0);
			for (int row = 0; row < matrix.rows(); ++row)
				result = std::hypot(result, std::abs(matrix(row, col)));
			return result;
		}

		static void RotateColumns(Matrix<Complex>& matrix, int p, int q, Real cosine, Real sine, Complex phase)
		{
			for (int row = 0; row < matrix.rows(); ++row) {
				const Complex left = matrix(row, p);
				const Complex right = matrix(row, q);
				matrix(row, p) = cosine * left - sine * std::conj(phase) * right;
				matrix(row, q) = sine * phase * left + cosine * right;
			}
		}

		static Complex SubtractCancellation(const Complex& left, const Complex& right)
		{
			const auto component = [](Real a, Real b) {
				const Real difference = a - b;
				return std::abs(difference) <= Precision::Eps *
					(std::abs(a) + std::abs(b)) ? Real(0) : difference;
			};
			return {component(left.real(), right.real()), component(left.imag(), right.imag())};
		}

		static void RotateWorkingColumns(Matrix<Complex>& matrix, int p, int q,
			Real cosine, Real sine, Real scaledSine, Complex phase, Real normP, Real normQ)
		{
			// Form the cancellation in unit-column coordinates before restoring the
			// small norm. Discard only componentwise cancellation roundoff, not a small
			// column norm: otherwise exact dependencies chase underflow for hundreds
			// of sweeps in extended-range types, while independent tiny entries must survive.
			const bool smallerLeft = normP <= normQ;
			for (int row = 0; row < matrix.rows(); ++row) {
				const Complex left = matrix(row, p), right = matrix(row, q);
				if (smallerLeft) {
					matrix(row, p) = normP * SubtractCancellation(cosine * (left / normP),
						scaledSine * std::conj(phase) * (right / normQ));
					matrix(row, q) = sine * phase * left + cosine * right;
				} else {
					matrix(row, p) = cosine * left - sine * std::conj(phase) * right;
					matrix(row, q) = normQ * SubtractCancellation(scaledSine * phase * (left / normP),
						-cosine * (right / normQ));
				}
			}
		}

		static Matrix<Complex> CompleteUnitaryBasis(const Matrix<Complex>& leading, int leadingColumns)
		{
			const int dimension = leading.rows();
			Matrix<Complex> result(dimension, dimension);
			for (int row = 0; row < dimension; ++row)
				for (int col = 0; col < leadingColumns; ++col)
					result(row, col) = leading(row, col);

			const Real completionTolerance = Precision::SVDCompletionTolerancePerDimension * dimension;
			int completed = leadingColumns;
			for (int basisIndex = 0; basisIndex < dimension && completed < dimension; ++basisIndex) {
				Vector<Complex> candidate(dimension);
				candidate[basisIndex] = Complex{static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0)};
				for (int pass = 0; pass < 2; ++pass)
					for (int col = 0; col < completed; ++col) {
						Complex projection{};
						for (int row = 0; row < dimension; ++row)
							projection += std::conj(result(row, col)) * candidate[row];
						for (int row = 0; row < dimension; ++row)
							candidate[row] -= projection * result(row, col);
					}

				const Real norm = candidate.NormL2();
				if (norm <= completionTolerance)
					continue;
				for (int row = 0; row < dimension; ++row)
					result(row, completed) = candidate[row] / norm;
				++completed;
			}

			if (completed != dimension)
				throw MatrixNumericalError("ComplexSVDecompositionSolver - failed to complete unitary basis");
			return result;
		}

		static Result DecomposeTall(const Matrix<Complex>& matrix, Real tolerance, int maxSweeps)
		{
			const int rows = matrix.rows();
			const int cols = matrix.cols();
			Real scale = static_cast<MML::Real>(0.0);
			for (int row = 0; row < rows; ++row)
				for (int col = 0; col < cols; ++col) {
					const Real magnitude = std::abs(matrix(row, col));
					if (!std::isfinite(magnitude))
						throw MatrixNumericalError("ComplexSVDecompositionSolver - non-finite input");
					scale = std::max(scale, magnitude);
				}

			Matrix<Complex> working(rows, cols);
			Matrix<Complex> right = Matrix<Complex>::Identity(cols);
			if (scale == static_cast<MML::Real>(0.0)) {
				Result zero;
				zero.U = Matrix<Complex>::Identity(rows);
				zero.V = std::move(right);
				zero.singularValues.resize(cols);
				return zero;
			}
			for (int row = 0; row < rows; ++row)
				for (int col = 0; col < cols; ++col)
					working(row, col) = matrix(row, col) / scale;

			const Real convergenceTolerance = std::max(tolerance, static_cast<MML::Real>(32.0) * std::numeric_limits<Real>::epsilon());
			Real maxCorrelation = std::numeric_limits<Real>::infinity();
			int sweep = 0;
			for (; sweep < maxSweeps; ++sweep) {
				maxCorrelation = static_cast<MML::Real>(0.0);
				bool rotated = false;
				for (int p = 0; p < cols - 1; ++p)
					for (int q = p + 1; q < cols; ++q) {
						const Real normP = ColumnNorm(working, p);
						const Real normQ = ColumnNorm(working, q);
						if (normP == static_cast<MML::Real>(0.0) || normQ == static_cast<MML::Real>(0.0))
							continue;
						const Complex normalizedGamma = NormalizedColumnInnerProduct(working, p, q, normP, normQ);
						const Real correlation = std::abs(normalizedGamma);
						maxCorrelation = std::max(maxCorrelation, correlation);
						if (correlation <= convergenceTolerance)
							continue;

						const Real pairScale = std::max(normP, normQ);
						const Real scaledP = normP / pairScale, scaledQ = normQ / pairScale;
						const Real delta = scaledQ * scaledQ - scaledP * scaledP;
						const Real twiceGamma = Real(2) * correlation * scaledP * scaledQ;
						// Equivalent to the usual zeta formula without forming a huge zeta.
						const Real tangent = delta == Real(0) ? Real(1) :
							std::copysign(twiceGamma / (std::abs(delta) + std::hypot(delta, twiceGamma)), delta);
						const Real cosine = static_cast<MML::Real>(1.0) / std::hypot(static_cast<MML::Real>(1.0), tangent);
						const Real sine = tangent * cosine;
						// Compute sine*(largerNorm/smallerNorm) without first rounding
						// a subnormal sine; it is the coefficient used for cancellation.
						const Real scaledSine = delta == Real(0) ? sine :
							std::copysign((Real(2) * correlation /
								(std::abs(delta) + std::hypot(delta, twiceGamma))) * cosine, delta);
						const Complex phase = normalizedGamma / correlation;
						RotateWorkingColumns(working, p, q, cosine, sine, scaledSine, phase, normP, normQ);
						RotateColumns(right, p, q, cosine, sine, phase);
						rotated = true;
					}
				if (!rotated || maxCorrelation <= convergenceTolerance)
					break;
			}
			if (sweep == maxSweeps)
				throw ConvergenceError("ComplexSVDecompositionSolver - Jacobi sweeps did not converge", maxSweeps, maxCorrelation);

			Vector<Real> singularValues(cols);
			for (int col = 0; col < cols; ++col)
				singularValues[col] = ColumnNorm(working, col) * scale;
			for (int left = 0; left < cols - 1; ++left) {
				int largest = left;
				for (int rightIndex = left + 1; rightIndex < cols; ++rightIndex)
					if (singularValues[rightIndex] > singularValues[largest])
						largest = rightIndex;
				if (largest != left) {
					std::swap(singularValues[left], singularValues[largest]);
					for (int row = 0; row < rows; ++row)
						std::swap(working(row, left), working(row, largest));
					for (int row = 0; row < cols; ++row)
						std::swap(right(row, left), right(row, largest));
				}
			}

			Matrix<Complex> economyLeft(rows, cols);
			int normalizedColumns = 0;
			for (int col = 0; col < cols; ++col) {
				const Real scaledNorm = ColumnNorm(working, col);
				if (scaledNorm == Real(0))
					break;
				for (int row = 0; row < rows; ++row)
					economyLeft(row, col) = working(row, col) / scaledNorm;
				++normalizedColumns;
			}

			Result result;
			result.U = CompleteUnitaryBasis(economyLeft, normalizedColumns);
			result.V = std::move(right);
			result.singularValues = std::move(singularValues);
			result.sweeps = sweep;
			result.maxCorrelation = maxCorrelation;
			return result;
		}
	};
}

#endif // MML_LINEAR_ALG_COMPLEX_SVD_H

///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        ComplexEigenSolver.h                                                ///
///  Description: Shifted-QR eigensolver for dense general complex matrices           ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_COMPLEX_EIGEN_SOLVER_H
#define MML_COMPLEX_EIGEN_SOLVER_H

#include <mml/MMLBase.h>
#include <mml/base/Matrix/Matrix.h>
#include <mml/base/Vector/Vector.h>
#include <mml/base/AlgorithmTypes.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <string>

namespace MML
{
	/// Dense complex eigensolver based on shifted unitary QR iteration.
	/// The iteration computes a complex Schur form and obtains right eigenvectors
	/// by triangular back-substitution followed by the accumulated unitary transform.
	class ComplexEigenSolver
	{
	public:
		struct Result
		{
			Vector<Complex> eigenvalues;
			Matrix<Complex> eigenvectors;
			bool converged = false;
			int iterations = 0;
			Real maxResidual = REAL(0.0);
			AlgorithmStatus status = AlgorithmStatus::AlgorithmSpecificFailure;
			std::string algorithmName = "ComplexQR";
			std::string errorMessage;
		};

		static Result Solve(
				const Matrix<Complex>& matrix,
				Real tolerance = PrecisionValues<Real>::EigenSolverConvergenceTolerance,
				int maxIterations = 2000)
		{
			ValidateInput(matrix, tolerance, maxIterations);
			const int size = matrix.rows();
			Result result;
			result.eigenvalues.Resize(size);
			result.eigenvectors.Resize(size, size);

			Real scale = REAL(0.0);
			for (int row = 0; row < size; ++row)
				for (int col = 0; col < size; ++col)
					scale = std::max(scale, static_cast<Real>(std::abs(matrix(row, col))));
			if (scale == REAL(0.0)) {
				result.eigenvectors = Matrix<Complex>::Identity(size);
				result.converged = true;
				result.status = AlgorithmStatus::Success;
				return result;
			}

			Matrix<Complex> schur(size, size);
			for (int row = 0; row < size; ++row)
				for (int col = 0; col < size; ++col)
					schur(row, col) = matrix(row, col) / scale;
			Matrix<Complex> schurVectors = Matrix<Complex>::Identity(size);

			int activeSize = size;
			while (activeSize > 1 && result.iterations < maxIterations) {
				if (CanDeflate(schur, activeSize, tolerance)) {
					for (int col = 0; col < activeSize - 1; ++col)
						schur(activeSize - 1, col) = Complex{};
					--activeSize;
					continue;
				}

				const Complex shift = WilkinsonShift(schur, activeSize);
				Matrix<Complex> shifted(activeSize, activeSize);
				for (int row = 0; row < activeSize; ++row)
					for (int col = 0; col < activeSize; ++col)
						shifted(row, col) = schur(row, col) - (row == col ? shift : Complex{});

				Matrix<Complex> unitary;
				Matrix<Complex> upper;
				HouseholderQR(shifted, unitary, upper);
				const Matrix<Complex> next = upper * unitary;
				for (int row = 0; row < activeSize; ++row)
					for (int col = 0; col < activeSize; ++col)
						schur(row, col) = next(row, col) + (row == col ? shift : Complex{});

				if (activeSize < size) {
					Matrix<Complex> transformed(activeSize, size - activeSize);
					for (int row = 0; row < activeSize; ++row)
						for (int col = activeSize; col < size; ++col) {
							Complex value{};
							for (int index = 0; index < activeSize; ++index)
								value += std::conj(unitary(index, row)) * schur(index, col);
							transformed(row, col - activeSize) = value;
						}
					for (int row = 0; row < activeSize; ++row)
						for (int col = activeSize; col < size; ++col)
							schur(row, col) = transformed(row, col - activeSize);
				}

				Matrix<Complex> accumulated(size, activeSize);
				for (int row = 0; row < size; ++row)
					for (int col = 0; col < activeSize; ++col) {
						Complex value{};
						for (int index = 0; index < activeSize; ++index)
							value += schurVectors(row, index) * unitary(index, col);
						accumulated(row, col) = value;
					}
				for (int row = 0; row < size; ++row)
					for (int col = 0; col < activeSize; ++col)
						schurVectors(row, col) = accumulated(row, col);
				++result.iterations;
			}

			while (activeSize > 1 && CanDeflate(schur, activeSize, tolerance)) {
				for (int col = 0; col < activeSize - 1; ++col)
					schur(activeSize - 1, col) = Complex{};
				--activeSize;
			}
			result.converged = activeSize <= 1;
			result.status = result.converged ? AlgorithmStatus::Success : AlgorithmStatus::MaxIterationsExceeded;
			if (!result.converged)
				result.errorMessage = "Complex QR iteration did not converge within " +
					std::to_string(maxIterations) + " iterations";

			for (int index = 0; index < size; ++index)
				result.eigenvalues[index] = schur(index, index) * scale;
			result.eigenvectors = ComputeEigenvectors(schur, schurVectors, tolerance);
			result.maxResidual = MaximumResidual(matrix, result.eigenvalues, result.eigenvectors);
			return result;
		}

	private:
		static void ValidateInput(const Matrix<Complex>& matrix, Real tolerance, int maxIterations)
		{
			if (matrix.rows() == 0 || matrix.rows() != matrix.cols())
				throw MatrixDimensionError("ComplexEigenSolver::Solve - matrix must be non-empty and square",
					matrix.rows(), matrix.cols(), -1, -1);
			if (!std::isfinite(tolerance) || tolerance < REAL(0.0))
				throw DomainError("ComplexEigenSolver::Solve - tolerance must be finite and nonnegative");
			if (maxIterations < 0)
				throw DomainError("ComplexEigenSolver::Solve - maxIterations must be nonnegative");
			for (int row = 0; row < matrix.rows(); ++row)
				for (int col = 0; col < matrix.cols(); ++col)
					if (!std::isfinite(matrix(row, col).real()) || !std::isfinite(matrix(row, col).imag()))
						throw MatrixNumericalError("ComplexEigenSolver::Solve - non-finite input");
		}

		static bool CanDeflate(const Matrix<Complex>& matrix, int activeSize, Real tolerance)
		{
			Real lowerNorm = REAL(0.0);
			for (int col = 0; col < activeSize - 1; ++col)
				lowerNorm = std::hypot(lowerNorm, static_cast<Real>(std::abs(matrix(activeSize - 1, col))));
			const Real scale = std::max(REAL(1.0),
				static_cast<Real>(std::abs(matrix(activeSize - 1, activeSize - 1))));
			return lowerNorm <= tolerance * scale;
		}

		static Complex WilkinsonShift(const Matrix<Complex>& matrix, int activeSize)
		{
			const Complex a = matrix(activeSize - 2, activeSize - 2);
			const Complex b = matrix(activeSize - 2, activeSize - 1);
			const Complex c = matrix(activeSize - 1, activeSize - 2);
			const Complex d = matrix(activeSize - 1, activeSize - 1);
			const Complex midpoint = (a + d) / REAL(2.0);
			const Complex discriminant = std::sqrt((a - d) * (a - d) / REAL(4.0) + b * c);
			const Complex first = midpoint + discriminant;
			const Complex second = midpoint - discriminant;
			return std::abs(first - d) <= std::abs(second - d) ? first : second;
		}

		static void HouseholderQR(
				const Matrix<Complex>& matrix, Matrix<Complex>& unitary, Matrix<Complex>& upper)
		{
			const int size = matrix.rows();
			upper = matrix;
			unitary = Matrix<Complex>::Identity(size);
			for (int col = 0; col < size - 1; ++col) {
				Real norm = REAL(0.0);
				for (int row = col; row < size; ++row)
					norm = std::hypot(norm, static_cast<Real>(std::abs(upper(row, col))));
				if (norm <= std::numeric_limits<Real>::min())
					continue;

				Vector<Complex> reflector(size - col);
				for (int row = col; row < size; ++row)
					reflector[row - col] = upper(row, col);
				const Complex phase = std::abs(reflector[0]) > REAL(0.0)
					? reflector[0] / std::abs(reflector[0]) : Complex{REAL(1.0), REAL(0.0)};
				reflector[0] += phase * norm;
				const Real reflectorNorm = reflector.NormL2();
				if (reflectorNorm <= std::numeric_limits<Real>::min())
					continue;
				for (int index = 0; index < reflector.size(); ++index)
					reflector[index] /= reflectorNorm;

				for (int targetCol = col; targetCol < size; ++targetCol) {
					Complex projection{};
					for (int row = col; row < size; ++row)
						projection += std::conj(reflector[row - col]) * upper(row, targetCol);
					for (int row = col; row < size; ++row)
						upper(row, targetCol) -= REAL(2.0) * reflector[row - col] * projection;
				}
				for (int row = col + 1; row < size; ++row)
					upper(row, col) = Complex{};

				for (int row = 0; row < size; ++row) {
					Complex projection{};
					for (int index = col; index < size; ++index)
						projection += unitary(row, index) * reflector[index - col];
					for (int index = col; index < size; ++index)
						unitary(row, index) -= REAL(2.0) * projection * std::conj(reflector[index - col]);
				}
			}
		}

		static Matrix<Complex> ComputeEigenvectors(
				const Matrix<Complex>& schur, const Matrix<Complex>& schurVectors, Real tolerance)
		{
			const int size = schur.rows();
			Matrix<Complex> result(size, size);
			const Real safeDenominator = std::max(tolerance,
				std::numeric_limits<Real>::epsilon() * static_cast<Real>(size * 16));
			for (int eigenIndex = 0; eigenIndex < size; ++eigenIndex) {
				const Complex eigenvalue = schur(eigenIndex, eigenIndex);
				Vector<Complex> schurVector(size);
				schurVector[eigenIndex] = Complex{REAL(1.0), REAL(0.0)};
				for (int row = eigenIndex - 1; row >= 0; --row) {
					Complex sum{};
					for (int col = row + 1; col <= eigenIndex; ++col)
						sum += schur(row, col) * schurVector[col];
					Complex denominator = schur(row, row) - eigenvalue;
					if (std::abs(denominator) < safeDenominator) {
						const Complex phase = std::abs(denominator) > REAL(0.0)
							? denominator / std::abs(denominator) : Complex{REAL(1.0), REAL(0.0)};
						denominator = phase * safeDenominator;
					}
					schurVector[row] = -sum / denominator;
				}

				Vector<Complex> eigenvector(size);
				for (int row = 0; row < size; ++row)
					for (int col = 0; col < size; ++col)
						eigenvector[row] += schurVectors(row, col) * schurVector[col];
				const Real norm = eigenvector.NormL2();
				if (norm > REAL(0.0))
					for (int row = 0; row < size; ++row)
						result(row, eigenIndex) = eigenvector[row] / norm;
			}
			return result;
		}

		static Real MaximumResidual(
				const Matrix<Complex>& matrix, const Vector<Complex>& eigenvalues,
				const Matrix<Complex>& eigenvectors)
		{
			Real matrixNorm = REAL(0.0);
			for (int row = 0; row < matrix.rows(); ++row)
				for (int col = 0; col < matrix.cols(); ++col)
					matrixNorm = std::hypot(matrixNorm, static_cast<Real>(std::abs(matrix(row, col))));
			Real maximum = REAL(0.0);
			for (int index = 0; index < matrix.rows(); ++index) {
				const Vector<Complex> vector = eigenvectors.VectorFromColumn(index);
				const Real residual = (matrix * vector - eigenvalues[index] * vector).NormL2() /
					std::max(REAL(1.0), matrixNorm * vector.NormL2());
				maximum = std::max(maximum, residual);
			}
			return maximum;
		}
	};
}

#endif // MML_COMPLEX_EIGEN_SOLVER_H

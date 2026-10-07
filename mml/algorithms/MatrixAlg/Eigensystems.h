///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Eigensystems.h                                                      ///
///  Description: Matrix eigensystem operations and Hessenberg facade                 ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_MATRIX_ALG_EIGENSYSTEMS_H
#define MML_MATRIX_ALG_EIGENSYSTEMS_H

#include <mml/algorithms/MatrixAlg/Properties.h>
#include <mml/algorithms/Eigen/ComplexEigenSolver.h>
#include <mml/algorithms/Eigen/EigenSolver.h>
#include <mml/algorithms/Eigen/HessenbergReduction.h>
#include <mml/algorithms/Eigen/HermitianMatEigenSolverJacobi.h>
#include <mml/algorithms/Eigen/SymmMatEigenSolverJacobi.h>
#include <mml/base/Matrix/MatrixSym.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <string>

namespace MML::MatrixAlg
{
	namespace Detail
	{
		template<MMLScalar Scalar>
		MatrixMagnitude<Scalar> MaximumEigenResidual(
				const Matrix<Scalar>& matrix,
				const Vector<MatrixComplexScalar<Scalar>>& eigenvalues,
				const Matrix<MatrixComplexScalar<Scalar>>& eigenvectors)
		{
			using Magnitude = MatrixMagnitude<Scalar>;
			using ComplexScalar = MatrixComplexScalar<Scalar>;
			Magnitude matrixNorm{};
			for (int row = 0; row < matrix.rows(); ++row)
				for (int col = 0; col < matrix.cols(); ++col)
					matrixNorm = std::hypot(matrixNorm, static_cast<Magnitude>(std::abs(matrix(row, col))));

			Magnitude maximum{};
			for (int eigenIndex = 0; eigenIndex < matrix.rows(); ++eigenIndex)
			{
				Magnitude residual{};
				Magnitude vectorNorm{};
				for (int row = 0; row < matrix.rows(); ++row)
				{
					ComplexScalar product{};
					for (int col = 0; col < matrix.cols(); ++col)
						product += static_cast<ComplexScalar>(matrix(row, col)) * eigenvectors(col, eigenIndex);
					residual = std::hypot(residual, static_cast<Magnitude>(
						std::abs(product - eigenvalues[eigenIndex] * eigenvectors(row, eigenIndex))));
					vectorNorm = std::hypot(vectorNorm,
						static_cast<Magnitude>(std::abs(eigenvectors(row, eigenIndex))));
				}
				maximum = std::max(maximum, residual / std::max(Magnitude{1}, matrixNorm * vectorNorm));
			}
			return maximum;
		}
	}

	inline SelfAdjointEigensystemResult<Real> SymmetricEigensystem(
			const Matrix<Real>& matrix,
			Real tolerance = PrecisionValues<Real>::EigenSolverConvergenceTolerance,
			int maxIterations = 100)
	{
		if (matrix.rows() == 0 || matrix.rows() != matrix.cols())
			throw MatrixDimensionError("SymmetricEigensystem - matrix must be non-empty and square",
				matrix.rows(), matrix.cols(), -1, -1);
		if (!std::isfinite(tolerance) || tolerance < REAL(0.0))
			throw DomainError("SymmetricEigensystem - tolerance must be finite and nonnegative");
		if (maxIterations < 0)
			throw DomainError("SymmetricEigensystem - maxIterations must be nonnegative");
		if (!IsSymmetric(matrix, MatrixComparisonTolerance<Real>{tolerance, tolerance}))
			throw MatrixDimensionError("SymmetricEigensystem - matrix must be symmetric",
				matrix.rows(), matrix.cols(), -1, -1);
		MatrixSym<Real> symmetric(matrix.rows());
		for (int row = 0; row < matrix.rows(); ++row)
			for (int col = row; col < matrix.cols(); ++col)
				symmetric(row, col) = REAL(0.5) * (matrix(row, col) + matrix(col, row));
		const auto decomposition = SymmMatEigenSolverJacobi::Solve(symmetric, tolerance, maxIterations);
		SelfAdjointEigensystemResult<Real> result;
		result.eigenvalues = decomposition.eigenvalues;
		result.eigenvectors = decomposition.eigenvectors;
		result.converged = decomposition.converged;
		result.iterations = decomposition.iterations;
		result.algorithmName = "SymmetricJacobi";
		result.status = result.converged ? AlgorithmStatus::Success : AlgorithmStatus::MaxIterationsExceeded;
		if (!result.converged)
			result.errorMessage = "Symmetric Jacobi iteration did not converge within " +
				std::to_string(maxIterations) + " sweeps";
		Vector<Complex> eigenvalues(result.eigenvalues.size());
		Matrix<Complex> eigenvectors(result.eigenvectors.rows(), result.eigenvectors.cols());
		for (int index = 0; index < eigenvalues.size(); ++index)
			eigenvalues[index] = Complex{result.eigenvalues[index], REAL(0.0)};
		for (int row = 0; row < eigenvectors.rows(); ++row)
			for (int col = 0; col < eigenvectors.cols(); ++col)
				eigenvectors(row, col) = Complex{result.eigenvectors(row, col), REAL(0.0)};
		result.maxResidual = Detail::MaximumEigenResidual(matrix, eigenvalues, eigenvectors);
		return result;
	}

	inline SelfAdjointEigensystemResult<Complex> HermitianEigensystem(
			const Matrix<Complex>& matrix,
			Real tolerance = PrecisionValues<Real>::EigenSolverConvergenceTolerance,
			int maxIterations = 100)
	{
		if (matrix.rows() == 0 || matrix.rows() != matrix.cols())
			throw MatrixDimensionError("HermitianEigensystem - matrix must be non-empty and square",
				matrix.rows(), matrix.cols(), -1, -1);
		if (!std::isfinite(tolerance) || tolerance < REAL(0.0))
			throw DomainError("HermitianEigensystem - tolerance must be finite and nonnegative");
		if (maxIterations < 0)
			throw DomainError("HermitianEigensystem - maxIterations must be nonnegative");
		if (!IsHermitian(matrix, MatrixComparisonTolerance<Real>{tolerance, tolerance}))
			throw MatrixDimensionError("HermitianEigensystem - matrix must be Hermitian",
				matrix.rows(), matrix.cols(), -1, -1);
		Matrix<Complex> hermitian(matrix.rows(), matrix.cols());
		for (int row = 0; row < matrix.rows(); ++row)
		{
			hermitian(row, row) = Complex{matrix(row, row).real(), REAL(0.0)};
			for (int col = row + 1; col < matrix.cols(); ++col)
			{
				const Complex value = REAL(0.5) * (matrix(row, col) + std::conj(matrix(col, row)));
				hermitian(row, col) = value;
				hermitian(col, row) = std::conj(value);
			}
		}
		const auto decomposition = HermitianMatEigenSolverJacobi::Solve(hermitian, tolerance, maxIterations);
		SelfAdjointEigensystemResult<Complex> result;
		result.eigenvalues = decomposition.eigenvalues;
		result.eigenvectors = decomposition.eigenvectors;
		result.converged = decomposition.converged;
		result.iterations = decomposition.iterations;
		result.algorithmName = decomposition.algorithm_name;
		result.status = result.converged ? AlgorithmStatus::Success : AlgorithmStatus::MaxIterationsExceeded;
		if (!result.converged)
			result.errorMessage = "Hermitian Jacobi iteration did not converge within " +
				std::to_string(maxIterations) + " sweeps";
		Vector<Complex> eigenvalues(result.eigenvalues.size());
		for (int index = 0; index < eigenvalues.size(); ++index)
			eigenvalues[index] = Complex{result.eigenvalues[index], REAL(0.0)};
		result.maxResidual = Detail::MaximumEigenResidual(matrix, eigenvalues, result.eigenvectors);
		return result;
	}

	inline Vector<Real> SymmetricEigenvalues(
			const Matrix<Real>& matrix,
			Real tolerance = PrecisionValues<Real>::EigenSolverConvergenceTolerance,
			int maxIterations = 100)
	{
		return SymmetricEigensystem(matrix, tolerance, maxIterations).eigenvalues;
	}

	inline Vector<Real> HermitianEigenvalues(
			const Matrix<Complex>& matrix,
			Real tolerance = PrecisionValues<Real>::EigenSolverConvergenceTolerance,
			int maxIterations = 100)
	{
		return HermitianEigensystem(matrix, tolerance, maxIterations).eigenvalues;
	}

	template<MMLScalar Scalar>
	EigensystemResult<Scalar> Eigensystem(
			const Matrix<Scalar>& matrix,
			MatrixMagnitude<Scalar> tolerance = PrecisionValues<MatrixMagnitude<Scalar>>::EigenSolverConvergenceTolerance,
			int maxIterations = 1000)
	{
		using ComplexScalar = MatrixComplexScalar<Scalar>;
		if (matrix.rows() == 0 || matrix.rows() != matrix.cols())
			throw MatrixDimensionError("Eigensystem - matrix must be non-empty and square",
				matrix.rows(), matrix.cols(), -1, -1);
		if (!std::isfinite(tolerance) || tolerance < MatrixMagnitude<Scalar>{})
			throw DomainError("Eigensystem - tolerance must be finite and nonnegative");
		if (maxIterations < 0)
			throw DomainError("Eigensystem - maxIterations must be nonnegative");

		EigensystemResult<Scalar> result;
		if constexpr (MMLComplex<Scalar>)
		{
			if (IsHermitian(matrix, MatrixComparisonTolerance<MatrixMagnitude<Scalar>>{tolerance, tolerance}))
			{
				const auto selfAdjoint = HermitianEigensystem(matrix, tolerance, maxIterations);
				result.eigenvalues.Resize(selfAdjoint.eigenvalues.size());
				for (int index = 0; index < result.eigenvalues.size(); ++index)
					result.eigenvalues[index] = ComplexScalar{selfAdjoint.eigenvalues[index], REAL(0.0)};
				result.eigenvectors = selfAdjoint.eigenvectors;
				result.converged = selfAdjoint.converged;
				result.iterations = selfAdjoint.iterations;
				result.maxResidual = selfAdjoint.maxResidual;
				result.status = selfAdjoint.status;
				result.algorithmName = selfAdjoint.algorithmName;
				result.errorMessage = selfAdjoint.errorMessage;
				return result;
			}

			const auto decomposition = ComplexEigenSolver::Solve(matrix, tolerance, maxIterations);
			result.eigenvalues = decomposition.eigenvalues;
			result.eigenvectors = decomposition.eigenvectors;
			result.converged = decomposition.converged;
			result.iterations = decomposition.iterations;
			result.maxResidual = decomposition.maxResidual;
			result.status = decomposition.status;
			result.algorithmName = decomposition.algorithmName;
			result.errorMessage = decomposition.errorMessage;
			return result;
		}
		else
		{
			if (IsSymmetric(matrix, MatrixComparisonTolerance<MatrixMagnitude<Scalar>>{tolerance, tolerance}))
			{
				const auto selfAdjoint = SymmetricEigensystem(matrix, tolerance, maxIterations);
				result.eigenvalues.Resize(selfAdjoint.eigenvalues.size());
				for (int index = 0; index < result.eigenvalues.size(); ++index)
					result.eigenvalues[index] = ComplexScalar{selfAdjoint.eigenvalues[index], REAL(0.0)};
				result.eigenvectors.Resize(matrix.rows(), matrix.cols());
				for (int row = 0; row < matrix.rows(); ++row)
					for (int col = 0; col < matrix.cols(); ++col)
						result.eigenvectors(row, col) = ComplexScalar{selfAdjoint.eigenvectors(row, col), REAL(0.0)};
				result.converged = selfAdjoint.converged;
				result.iterations = selfAdjoint.iterations;
				result.maxResidual = selfAdjoint.maxResidual;
				result.status = selfAdjoint.status;
				result.algorithmName = selfAdjoint.algorithmName;
				result.errorMessage = selfAdjoint.errorMessage;
				return result;
			}

			const auto legacy = EigenSolver::Solve(matrix, tolerance, maxIterations);
			result.eigenvalues.Resize(matrix.rows());
			result.eigenvectors.Resize(matrix.rows(), matrix.cols());
			for (int index = 0; index < matrix.rows(); ++index)
				result.eigenvalues[index] = ComplexScalar{legacy.eigenvalues[index].real, legacy.eigenvalues[index].imag};
			for (int index = 0; index < matrix.rows(); ++index)
			{
				if (legacy.isComplexPair[index] && legacy.eigenvalues[index].imag > tolerance && index + 1 < matrix.rows())
				{
					for (int row = 0; row < matrix.rows(); ++row)
					{
						const ComplexScalar value{legacy.eigenvectors(row, index), legacy.eigenvectors(row, index + 1)};
						result.eigenvectors(row, index) = value;
						result.eigenvectors(row, index + 1) = std::conj(value);
					}
					++index;
				}
				else if (!legacy.isComplexPair[index])
				{
					for (int row = 0; row < matrix.rows(); ++row)
						result.eigenvectors(row, index) = ComplexScalar{legacy.eigenvectors(row, index), REAL(0.0)};
				}
			}
			result.converged = legacy.converged;
			result.iterations = legacy.iterations;
			result.status = result.converged ? AlgorithmStatus::Success : AlgorithmStatus::MaxIterationsExceeded;
			result.algorithmName = "GeneralRealQR";
			if (!result.converged)
				result.errorMessage = "General real QR iteration did not converge within " +
					std::to_string(maxIterations) + " iterations";
			result.maxResidual = Detail::MaximumEigenResidual(matrix, result.eigenvalues, result.eigenvectors);
			return result;
		}
	}

	template<MMLScalar Scalar>
	Vector<MatrixComplexScalar<Scalar>> Eigenvalues(
			const Matrix<Scalar>& matrix,
			MatrixMagnitude<Scalar> tolerance = PrecisionValues<MatrixMagnitude<Scalar>>::EigenSolverConvergenceTolerance,
			int maxIterations = 1000)
	{
		return Eigensystem(matrix, tolerance, maxIterations).eigenvalues;
	}

	template<MMLScalar Scalar>
	MatrixMagnitude<Scalar> SpectralRadius(
			const Matrix<Scalar>& matrix,
			MatrixMagnitude<Scalar> tolerance = PrecisionValues<MatrixMagnitude<Scalar>>::EigenSolverConvergenceTolerance,
			int maxIterations = 1000)
	{
		MatrixMagnitude<Scalar> maximum{};
		for (const auto& eigenvalue : Eigensystem(matrix, tolerance, maxIterations).eigenvalues)
			maximum = std::max(maximum, static_cast<MatrixMagnitude<Scalar>>(std::abs(eigenvalue)));
		return maximum;
	}

	/// @brief Result of Hessenberg reduction: H = Q^T * A * Q
	using HessenbergResult = MML::HessenbergResult;

	/// Reduce matrix A to upper Hessenberg form using Householder reflections.
	/// Delegates to the canonical public Eigen implementation.
	inline HessenbergResult ReduceToHessenberg(const Matrix<Real>& A)
	{
		return MML::ReduceToHessenberg(A);
	}
}

#endif // MML_MATRIX_ALG_EIGENSYSTEMS_H
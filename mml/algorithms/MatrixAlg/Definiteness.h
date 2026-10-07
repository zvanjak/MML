///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Definiteness.h                                                      ///
///  Description: Matrix definiteness classification and predicates                  ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_MATRIX_ALG_DEFINITENESS_H
#define MML_MATRIX_ALG_DEFINITENESS_H

#include <mml/algorithms/MatrixAlg/Properties.h>
#include <mml/algorithms/Eigen/HermitianMatEigenSolverJacobi.h>
#include <mml/algorithms/Eigen/SymmMatEigenSolverJacobi.h>
#include <mml/base/Matrix/MatrixSym.h>

#include <cmath>
#include <complex>

namespace MML::MatrixAlg
{
	namespace Detail
	{
		template<MMLReal Magnitude>
		Definiteness ClassifyEigenvalueSigns(const Vector<Magnitude>& eigenvalues, Magnitude tolerance)
		{
			int numPositive = 0;
			int numNegative = 0;
			int numZero = 0;

			for (int index = 0; index < eigenvalues.size(); ++index)
			{
				if (eigenvalues[index] > tolerance)
					++numPositive;
				else if (eigenvalues[index] < -tolerance)
					++numNegative;
				else
					++numZero;
			}

			if (numZero == eigenvalues.size())
				return Definiteness::ZeroSemidefinite;
			if (numNegative == 0 && numZero == 0)
				return Definiteness::PositiveDefinite;
			if (numNegative == 0)
				return Definiteness::PositiveSemidefinite;
			if (numPositive == 0 && numZero == 0)
				return Definiteness::NegativeDefinite;
			if (numPositive == 0)
				return Definiteness::NegativeSemidefinite;
			return Definiteness::Indefinite;
		}
	}

	/// Classify the definiteness of a symmetric matrix.
	inline Definiteness ClassifyDefiniteness(
			const MatrixSym<Real>& A,
			Real tol = PrecisionValues<Real>::EigenSolverConvergenceTolerance)
	{
		if (A.rows() == 0)
			throw MatrixDimensionError("ClassifyDefiniteness - matrix must be non-empty", 0, 0, -1, -1);
		if (!std::isfinite(tol) || tol < REAL(0.0))
			throw DomainError("ClassifyDefiniteness tolerance must be finite and nonnegative");

		auto result = SymmMatEigenSolverJacobi::Solve(A, tol);
		if (!result.converged)
			throw MatrixNumericalError("ClassifyDefiniteness - symmetric eigensolver did not converge");

		return Detail::ClassifyEigenvalueSigns(result.eigenvalues, tol);
	}

	/// Classify definiteness of a symmetric or Hermitian matrix using general storage.
	template<MMLScalar Scalar>
	Definiteness ClassifyDefiniteness(
			const Matrix<Scalar>& A,
			MatrixMagnitude<Scalar> tol = PrecisionValues<MatrixMagnitude<Scalar>>::EigenSolverConvergenceTolerance)
	{
		using Magnitude = MatrixMagnitude<Scalar>;
		int n = A.rows();
		if (A.cols() != n || n == 0)
			throw MatrixDimensionError("ClassifyDefiniteness - matrix must be non-empty and square", n, A.cols(), -1, -1);
		if (!std::isfinite(tol) || tol < Magnitude{})
			throw DomainError("ClassifyDefiniteness tolerance must be finite and nonnegative");

		if constexpr (MMLComplex<Scalar>)
		{
			if (!IsHermitian(A, MatrixComparisonTolerance<Magnitude>{tol, tol}))
				throw MatrixDimensionError("ClassifyDefiniteness - matrix must be Hermitian", n, A.cols(), -1, -1);
			Matrix<Complex> hermitian(n, n);
			for (int row = 0; row < n; ++row)
			{
				hermitian(row, row) = Complex{A(row, row).real(), REAL(0.0)};
				for (int col = row + 1; col < n; ++col)
				{
					const Complex value = REAL(0.5) * (A(row, col) + std::conj(A(col, row)));
					hermitian(row, col) = value;
					hermitian(col, row) = std::conj(value);
				}
			}
			const auto result = HermitianMatEigenSolverJacobi::Solve(hermitian, tol);
			if (!result.converged)
				throw MatrixNumericalError("ClassifyDefiniteness - Hermitian eigensolver did not converge");
			return Detail::ClassifyEigenvalueSigns(result.eigenvalues, tol);
		}
		else
		{
			if (!IsSymmetric(A, MatrixComparisonTolerance<Magnitude>{tol, tol}))
				throw MatrixDimensionError("ClassifyDefiniteness - matrix must be symmetric", n, A.cols(), -1, -1);
			MatrixSym<Real> symmetric(n);
			for (int row = 0; row < n; ++row)
				for (int col = row; col < n; ++col)
					symmetric(row, col) = REAL(0.5) * (A(row, col) + A(col, row));
			return ClassifyDefiniteness(symmetric, tol);
		}
	}

	inline bool IsPositiveDefinite(const MatrixSym<Real>& A,
			Real tol = PrecisionValues<Real>::EigenSolverConvergenceTolerance)
	{
		return ClassifyDefiniteness(A, tol) == Definiteness::PositiveDefinite;
	}

	template<MMLScalar Scalar>
	bool IsPositiveDefinite(const Matrix<Scalar>& A,
			MatrixMagnitude<Scalar> tol = PrecisionValues<MatrixMagnitude<Scalar>>::EigenSolverConvergenceTolerance)
	{
		return ClassifyDefiniteness(A, tol) == Definiteness::PositiveDefinite;
	}

	inline bool IsPositiveSemiDefinite(const MatrixSym<Real>& A,
			Real tol = PrecisionValues<Real>::EigenSolverConvergenceTolerance)
	{
		auto def = ClassifyDefiniteness(A, tol);
		return def == Definiteness::PositiveDefinite || def == Definiteness::PositiveSemidefinite ||
			   def == Definiteness::ZeroSemidefinite;
	}

	template<MMLScalar Scalar>
	bool IsPositiveSemiDefinite(const Matrix<Scalar>& A,
			MatrixMagnitude<Scalar> tol = PrecisionValues<MatrixMagnitude<Scalar>>::EigenSolverConvergenceTolerance)
	{
		auto def = ClassifyDefiniteness(A, tol);
		return def == Definiteness::PositiveDefinite || def == Definiteness::PositiveSemidefinite ||
			   def == Definiteness::ZeroSemidefinite;
	}

	inline bool IsNegativeDefinite(const MatrixSym<Real>& A,
			Real tol = PrecisionValues<Real>::EigenSolverConvergenceTolerance)
	{
		return ClassifyDefiniteness(A, tol) == Definiteness::NegativeDefinite;
	}

	template<MMLScalar Scalar>
	bool IsNegativeDefinite(const Matrix<Scalar>& A,
			MatrixMagnitude<Scalar> tol = PrecisionValues<MatrixMagnitude<Scalar>>::EigenSolverConvergenceTolerance)
	{
		return ClassifyDefiniteness(A, tol) == Definiteness::NegativeDefinite;
	}

	inline bool IsNegativeSemiDefinite(const MatrixSym<Real>& A,
			Real tol = PrecisionValues<Real>::EigenSolverConvergenceTolerance)
	{
		auto def = ClassifyDefiniteness(A, tol);
		return def == Definiteness::NegativeDefinite || def == Definiteness::NegativeSemidefinite ||
			   def == Definiteness::ZeroSemidefinite;
	}

	template<MMLScalar Scalar>
	bool IsNegativeSemiDefinite(const Matrix<Scalar>& A,
			MatrixMagnitude<Scalar> tol = PrecisionValues<MatrixMagnitude<Scalar>>::EigenSolverConvergenceTolerance)
	{
		auto def = ClassifyDefiniteness(A, tol);
		return def == Definiteness::NegativeDefinite || def == Definiteness::NegativeSemidefinite ||
			   def == Definiteness::ZeroSemidefinite;
	}

	inline bool IsIndefinite(const MatrixSym<Real>& A,
			Real tol = PrecisionValues<Real>::EigenSolverConvergenceTolerance)
	{
		return ClassifyDefiniteness(A, tol) == Definiteness::Indefinite;
	}

	template<MMLScalar Scalar>
	bool IsIndefinite(const Matrix<Scalar>& A,
			MatrixMagnitude<Scalar> tol = PrecisionValues<MatrixMagnitude<Scalar>>::EigenSolverConvergenceTolerance)
	{
		return ClassifyDefiniteness(A, tol) == Definiteness::Indefinite;
	}
}

#endif // MML_MATRIX_ALG_DEFINITENESS_H
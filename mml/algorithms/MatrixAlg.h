///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        MatrixAlg.h                                                         ///
///  Description: Matrix algebra algorithms requiring eigensolvers                    ///
///               Hessenberg reduction, definiteness classification                   ///
///               This is the "top of pyramid" - uses solvers from /algorithms        ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_MATRIX_ALG_H
#define MML_MATRIX_ALG_H

#include "base/Vector/Vector.h"
#include "base/Matrix/Matrix.h"
#include "base/Matrix/MatrixSym.h"
#include "core/MatrixUtils.h"
#include "algorithms/EigenSystemSolvers.h"

#include <cmath>

namespace MML {
	/// @namespace MatrixAlg
	/// @brief Collection of matrix algebra algorithms requiring eigensolvers
	///
	/// This namespace provides high-level matrix operations that depend on
	/// eigenvalue solvers from /algorithms. For operations that only need
	/// LU or SVD decomposition, use Utils namespace in MatrixUtils.h.
	///
	/// Contains:
	/// - Hessenberg reduction
	/// - Matrix definiteness classification
	///
	/// For SVD-based utilities (rank, condition number, null space, etc.),
	/// see Utils namespace in MatrixUtils.h.
	namespace MatrixAlg {
		// =========================================================================
		// HESSENBERG REDUCTION - Delegating to EigenSolverHelpers
		// =========================================================================

		/// @struct HessenbergResult
		/// @brief Result of Hessenberg reduction: H = Q^T * A * Q
		/// @note This is a type alias for EigenSolverHelpers::HessenbergResult
		using HessenbergResult = EigenSolverHelpers::HessenbergResult;

		/// Reduce matrix A to upper Hessenberg form using Householder reflections.
		///
		/// MATHEMATICAL DEFINITION:
		/// Find orthogonal Q such that H = Q^T * A * Q is upper Hessenberg,
		/// meaning H(i,j) = 0 for all i > j + 1.
		///
		/// The transformation preserves eigenvalues (similarity transformation).
		///
		/// @param A Input matrix (n x n)
		/// @return HessenbergResult containing H and Q where H = Q^T * A * Q
		///
		/// ALGORITHM (Householder):
		/// For k = 0, 1, ..., n-3:
		/// 1. Let x = A(k+1:n, k) be the subdiagonal column
		/// 2. Compute Householder vector v to zero x(1:end)
		/// 3. Apply P = I - 2*v*v^T from left: A = P*A
		/// 4. Apply P from right: A = A*P (similarity)
		/// 5. Accumulate Q = Q * P
		///
		/// Complexity: O(10n³/3) flops
		///
		/// REFERENCES:
		/// - Golub & Van Loan, "Matrix Computations", 4th ed., Section 7.4
		///
		/// @note Implementation: Delegates to EigenSolverHelpers::ReduceToHessenberg()
		///       which is the canonical implementation used by all eigen solvers.
		inline HessenbergResult ReduceToHessenberg(const Matrix<Real>& A) {
			return EigenSolverHelpers::ReduceToHessenberg(A);
		}

		// =========================================================================
		// MATRIX SYMMETRY CHECKS - Delegating to Utils
		// =========================================================================

		/// Check if matrix is symmetric within tolerance: A[i][j] ≈ A[j][i]
		/// @see Utils::IsSymmetric for implementation
		inline bool IsSymmetric(const Matrix<Real>& A, Real tol = 1e-10) { return Utils::IsSymmetric(A, tol); }

		/// Check if matrix is skew-symmetric: A[i][j] ≈ -A[j][i]
		/// @see Utils::IsSkewSymmetric for implementation
		inline bool IsSkewSymmetric(const Matrix<Real>& A, Real tol = 1e-10) { return Utils::IsSkewSymmetric(A, tol); }

		// =========================================================================
		// MATRIX DEFINITENESS CHECKS
		// =========================================================================

		/// @enum Definiteness
		/// @brief Classification of matrix definiteness
		enum class Definiteness {
			PositiveDefinite,	  ///< All eigenvalues > 0
			PositiveSemiDefinite, ///< All eigenvalues >= 0, at least one = 0
			NegativeDefinite,	  ///< All eigenvalues < 0
			NegativeSemiDefinite, ///< All eigenvalues <= 0, at least one = 0
			Indefinite			  ///< Has both positive and negative eigenvalues
		};

		/// Classify the definiteness of a symmetric matrix.
		///
		/// MATHEMATICAL DEFINITION:
		/// - Positive definite: x^T A x > 0 for all x ≠ 0 (all eigenvalues > 0)
		/// - Positive semi-definite: x^T A x >= 0 for all x (eigenvalues >= 0)
		/// - Negative definite: x^T A x < 0 for all x ≠ 0 (all eigenvalues < 0)
		/// - Negative semi-definite: x^T A x <= 0 for all x (eigenvalues <= 0)
		/// - Indefinite: Has both positive and negative eigenvalues
		///
		/// @param A Symmetric matrix (only upper triangle used)
		/// @param tol Tolerance for zero eigenvalue detection
		/// @return Definiteness classification
		///
		/// @note Matrix must be symmetric. For non-symmetric matrices,
		/// consider using (A + A^T)/2.
		inline Definiteness ClassifyDefiniteness(const MatrixSym<Real>& A, Real tol = 1e-10) {
			// Compute eigenvalues using Jacobi method
			auto result = SymmMatEigenSolverJacobi::Solve(A, tol);

			int n = A.rows();
			int numPositive = 0;
			int numNegative = 0;
			int numZero = 0;

			for (int i = 0; i < n; i++) {
				if (result.eigenvalues[i] > tol)
					numPositive++;
				else if (result.eigenvalues[i] < -tol)
					numNegative++;
				else
					numZero++;
			}

			if (numNegative == 0 && numZero == 0)
				return Definiteness::PositiveDefinite;
			else if (numNegative == 0 && numPositive > 0)
				return Definiteness::PositiveSemiDefinite;
			else if (numPositive == 0 && numZero == 0)
				return Definiteness::NegativeDefinite;
			else if (numPositive == 0 && numNegative > 0)
				return Definiteness::NegativeSemiDefinite;
			else
				return Definiteness::Indefinite;
		}

		/// Classify definiteness of a general matrix (uses symmetric part)
		inline Definiteness ClassifyDefiniteness(const Matrix<Real>& A, Real tol = 1e-10) {
			int n = A.rows();
			if (A.cols() != n)
				throw MatrixDimensionError("ClassifyDefiniteness - matrix must be square", n, A.cols(), -1, -1);

			// Convert to symmetric matrix (A + A^T) / 2
			MatrixSym<Real> symA(n);
			for (int i = 0; i < n; i++)
				for (int j = i; j < n; j++)
					symA(i, j) = 0.5 * (A(i, j) + A(j, i));

			return ClassifyDefiniteness(symA, tol);
		}

		/// Check if symmetric matrix is positive definite.
		///
		/// A matrix is positive definite if x^T A x > 0 for all x ≠ 0,
		/// which is equivalent to all eigenvalues being positive.
		///
		/// @param A Symmetric matrix
		/// @param tol Tolerance for eigenvalue comparison
		/// @return true if A is positive definite
		inline bool IsPositiveDefinite(const MatrixSym<Real>& A, Real tol = 1e-10) {
			return ClassifyDefiniteness(A, tol) == Definiteness::PositiveDefinite;
		}

		/// Check if general matrix is positive definite (uses symmetric part)
		inline bool IsPositiveDefinite(const Matrix<Real>& A, Real tol = 1e-10) {
			return ClassifyDefiniteness(A, tol) == Definiteness::PositiveDefinite;
		}

		/// Check if symmetric matrix is positive semi-definite.
		///
		/// A matrix is positive semi-definite if x^T A x >= 0 for all x,
		/// which is equivalent to all eigenvalues being non-negative.
		///
		/// @param A Symmetric matrix
		/// @param tol Tolerance for eigenvalue comparison
		/// @return true if A is positive semi-definite (includes positive definite)
		inline bool IsPositiveSemiDefinite(const MatrixSym<Real>& A, Real tol = 1e-10) {
			auto def = ClassifyDefiniteness(A, tol);
			return def == Definiteness::PositiveDefinite || def == Definiteness::PositiveSemiDefinite;
		}

		inline bool IsPositiveSemiDefinite(const Matrix<Real>& A, Real tol = 1e-10) {
			auto def = ClassifyDefiniteness(A, tol);
			return def == Definiteness::PositiveDefinite || def == Definiteness::PositiveSemiDefinite;
		}

		/// Check if symmetric matrix is negative definite.
		///
		/// A matrix is negative definite if x^T A x < 0 for all x ≠ 0,
		/// which is equivalent to all eigenvalues being negative.
		///
		/// @param A Symmetric matrix
		/// @param tol Tolerance for eigenvalue comparison
		/// @return true if A is negative definite
		inline bool IsNegativeDefinite(const MatrixSym<Real>& A, Real tol = 1e-10) {
			return ClassifyDefiniteness(A, tol) == Definiteness::NegativeDefinite;
		}

		inline bool IsNegativeDefinite(const Matrix<Real>& A, Real tol = 1e-10) {
			return ClassifyDefiniteness(A, tol) == Definiteness::NegativeDefinite;
		}

		/// Check if symmetric matrix is negative semi-definite.
		///
		/// @param A Symmetric matrix
		/// @param tol Tolerance for eigenvalue comparison
		/// @return true if A is negative semi-definite (includes negative definite)
		inline bool IsNegativeSemiDefinite(const MatrixSym<Real>& A, Real tol = 1e-10) {
			auto def = ClassifyDefiniteness(A, tol);
			return def == Definiteness::NegativeDefinite || def == Definiteness::NegativeSemiDefinite;
		}

		inline bool IsNegativeSemiDefinite(const Matrix<Real>& A, Real tol = 1e-10) {
			auto def = ClassifyDefiniteness(A, tol);
			return def == Definiteness::NegativeDefinite || def == Definiteness::NegativeSemiDefinite;
		}

		/// Check if matrix is indefinite (has both positive and negative eigenvalues).
		///
		/// @param A Symmetric matrix
		/// @param tol Tolerance for eigenvalue comparison
		/// @return true if A is indefinite
		inline bool IsIndefinite(const MatrixSym<Real>& A, Real tol = 1e-10) {
			return ClassifyDefiniteness(A, tol) == Definiteness::Indefinite;
		}

		inline bool IsIndefinite(const Matrix<Real>& A, Real tol = 1e-10) {
			return ClassifyDefiniteness(A, tol) == Definiteness::Indefinite;
		}

	} // namespace MatrixAlg

} // namespace MML

#endif // MML_MATRIX_ALG_H
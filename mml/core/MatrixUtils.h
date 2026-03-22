///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        MatrixUtils.h                                                       ///
///  Description: Matrix utility functions - property checks, norms, special matrices ///
///               Uses Core solvers for det, rank, ... calculations                   ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_MATRIXUTILS_H
#define MML_MATRIXUTILS_H

#include "MMLBase.h"

#include "base/Vector/Vector.h"
#include "base/Vector/VectorN.h"
#include "base/Matrix/Matrix.h"
#include "base/Matrix/MatrixNM.h"
#include "base/Polynom.h"

#include "core/LinAlgEqSolvers.h"

#include <limits>

namespace MML
{
	namespace Utils
	{
		// =========================================================================
		// MATRIX PROPERTY CHECKS (Simple O(n²) operations)
		// =========================================================================

		/**
		 * @brief Check if matrix H is upper Hessenberg within tolerance.
		 * 
		 * Upper Hessenberg: H(i,j) = 0 for all i > j + 1
		 * (all elements below the first subdiagonal are zero)
		 * 
		 * @param H Matrix to check
		 * @param tol Tolerance for zero comparison
		 * @return true if H is upper Hessenberg
		 */
		template<class Type>
		bool IsUpperHessenberg(const Matrix<Type>& H, Type tol = 1e-10)
		{
			int n = H.rows();
			for (int i = 2; i < n; i++)
				for (int j = 0; j < i - 1; j++)
					if (std::abs(H(i, j)) > tol)
						return false;
			return true;
		}

		/**
		 * @brief Check if matrix is upper triangular within tolerance.
		 * 
		 * @param A Matrix to check
		 * @param tol Tolerance for zero comparison
		 * @return true if A is upper triangular
		 */
		template<class Type>
		bool IsUpperTriangular(const Matrix<Type>& A, Type tol = 1e-10)
		{
			int n = A.rows();
			for (int i = 1; i < n; i++)
				for (int j = 0; j < i; j++)
					if (std::abs(A(i, j)) > tol)
						return false;
			return true;
		}

		/**
		 * @brief Check if matrix is lower triangular within tolerance.
		 * 
		 * @param A Matrix to check
		 * @param tol Tolerance for zero comparison
		 * @return true if A is lower triangular
		 */
		template<class Type>
		bool IsLowerTriangular(const Matrix<Type>& A, Type tol = 1e-10)
		{
			int n = A.rows();
			for (int i = 0; i < n; i++)
				for (int j = i + 1; j < n; j++)
					if (std::abs(A(i, j)) > tol)
						return false;
			return true;
		}

		/**
		 * @brief Check if matrix is diagonal within tolerance.
		 * 
		 * @param A Matrix to check
		 * @param tol Tolerance for zero comparison
		 * @return true if A is diagonal
		 */
		template<class Type>
		bool IsDiagonal(const Matrix<Type>& A, Type tol = 1e-10)
		{
			int n = A.rows();
			for (int i = 0; i < n; i++)
				for (int j = 0; j < n; j++)
					if (i != j && std::abs(A(i, j)) > tol)
						return false;
			return true;
		}

		/**
		 * @brief Check if matrix is symmetric within tolerance: A[i][j] ≈ A[j][i]
		 * 
		 * @param A Matrix to check
		 * @param tol Tolerance for comparison
		 * @return true if A is symmetric
		 */
		template<class Type>
		bool IsSymmetric(const Matrix<Type>& A, Type tol = 1e-10)
		{
			int n = A.rows();
			if (A.cols() != n) return false;  // Must be square
			
			for (int i = 0; i < n; i++)
				for (int j = i + 1; j < n; j++)
					if (std::abs(A(i, j) - A(j, i)) > tol)
						return false;
			return true;
		}

		/**
		 * @brief Check if matrix is skew-symmetric: A[i][j] ≈ -A[j][i]
		 * 
		 * @param A Matrix to check
		 * @param tol Tolerance for comparison
		 * @return true if A is skew-symmetric
		 */
		template<class Type>
		bool IsSkewSymmetric(const Matrix<Type>& A, Type tol = 1e-10)
		{
			int n = A.rows();
			if (A.cols() != n) return false;  // Must be square
			
			for (int i = 0; i < n; i++)
			{
				// Diagonal must be zero for skew-symmetric
				if (std::abs(A(i, i)) > tol)
					return false;
				for (int j = i + 1; j < n; j++)
					if (std::abs(A(i, j) + A(j, i)) > tol)
						return false;
			}
			return true;
		}

		// =========================================================================
		// MATRIX SCALAR QUANTITIES (Simple O(n²) operations)
		// =========================================================================

		/**
		 * @brief Compute trace (sum of diagonal elements)
		 * 
		 * @param A Input matrix
		 * @return Sum of diagonal elements
		 */
		template<class Type>
		Type Trace(const Matrix<Type>& A)
		{
			Type sum = 0;
			int n = std::min(A.rows(), A.cols());
			for (int i = 0; i < n; i++)
				sum += A(i, i);
			return sum;
		}

		/**
		 * @brief Compute Frobenius norm: ||A||_F = sqrt(sum of all a_ij²)
		 * 
		 * @param A Input matrix
		 * @return Frobenius norm
		 */
		template<class Type>
		Real FrobeniusNorm(const Matrix<Type>& A)
		{
			Real sum = 0.0;
			for (int i = 0; i < A.rows(); i++)
				for (int j = 0; j < A.cols(); j++)
					sum += std::abs(A(i, j)) * std::abs(A(i, j));
			return std::sqrt(sum);
		}

		/**
		 * @brief Compute infinity norm (maximum absolute row sum): ||A||_∞
		 * 
		 * @param A Input matrix
		 * @return Infinity norm
		 */
		template<class Type>
		Real InfinityNorm(const Matrix<Type>& A)
		{
			Real maxRowSum = 0.0;
			for (int i = 0; i < A.rows(); i++)
			{
				Real rowSum = 0.0;
				for (int j = 0; j < A.cols(); j++)
					rowSum += std::abs(A(i, j));
				maxRowSum = std::max(maxRowSum, rowSum);
			}
			return maxRowSum;
		}

		/**
		 * @brief Compute 1-norm (maximum absolute column sum): ||A||_1
		 * 
		 * @param A Input matrix
		 * @return 1-norm
		 */
		template<class Type>
		Real OneNorm(const Matrix<Type>& A)
		{
			Real maxColSum = 0.0;
			for (int j = 0; j < A.cols(); j++)
			{
				Real colSum = 0.0;
				for (int i = 0; i < A.rows(); i++)
					colSum += std::abs(A(i, j));
				maxColSum = std::max(maxColSum, colSum);
			}
			return maxColSum;
		}

		/**
		 * @brief Compute maximum absolute difference between two matrices
		 * 
		 * @param A First matrix
		 * @param B Second matrix
		 * @return Maximum |A(i,j) - B(i,j)|
		 */
		template<class Type>
		Real MaxAbsDiff(const Matrix<Type>& A, const Matrix<Type>& B)
		{
			Real maxDiff = 0.0;
			int rows = std::min(A.rows(), B.rows());
			int cols = std::min(A.cols(), B.cols());
			for (int i = 0; i < rows; i++)
				for (int j = 0; j < cols; j++)
					maxDiff = std::max(maxDiff, static_cast<Real>(std::abs(A(i, j) - B(i, j))));
			return maxDiff;
		}

		// =========================================================================
		// MATRIX TRANSFORMATIONS (Simple operations)
		// =========================================================================

		/**
		 * @brief Compute Q^T * A * Q (similarity transformation)
		 * 
		 * @param Q Transformation matrix
		 * @param A Matrix to transform
		 * @return Q^T * A * Q
		 */
		template<class Type>
		Matrix<Type> SimilarityTransform(const Matrix<Type>& Q, const Matrix<Type>& A)
		{
			int n = A.rows();
			Matrix<Type> temp(n, n);
			Matrix<Type> result(n, n);
			
			// temp = A * Q
			for (int i = 0; i < n; i++)
				for (int j = 0; j < n; j++)
				{
					temp(i, j) = 0;
					for (int k = 0; k < n; k++)
						temp(i, j) += A(i, k) * Q(k, j);
				}
			
			// result = Q^T * temp
			for (int i = 0; i < n; i++)
				for (int j = 0; j < n; j++)
				{
					result(i, j) = 0;
					for (int k = 0; k < n; k++)
						result(i, j) += Q(k, i) * temp(k, j);
				}
			
			return result;
		}

		// =========================================================================
		// SPECIAL MATRIX PROPERTIES (No solver dependencies)
		// =========================================================================

		/**
		 * @brief Check if matrix is nilpotent (A^n = 0 for some n)
		 */
		template<class Type>
		bool IsNilpotent(const Matrix<Type>& A)
		{
			int N = A.rows();
			if (N != A.cols())
				throw MatrixDimensionError("IsNilpotent - must be square matrix", N, A.cols(), -1, -1);
			
			Matrix<Type> M = A;
			for (int k = 1; k < N; k++)
			{
				M = A * M;
				if (M.isZero())
					return true;
			}
			return false;
		}

		/**
		 * @brief Check if matrix is unipotent ((A-I)^n = 0 for some n)
		 */
		template<class Type>
		bool IsUnipotent(const Matrix<Type>& A)
		{
			int N = A.rows();
			if (N != A.cols())
				throw MatrixDimensionError("IsUnipotent - must be square matrix", N, A.cols(), -1, -1);

			Matrix<Type> B = A - Matrix<Type>::Identity(N);
			Matrix<Type> M = B;
			for (int k = 1; k < N; k++)
			{
				if (M.isZero())
					return true;
				M = B * M;
			}
			return M.isZero();
		}

		/**
		 * @brief Faddeev-Leverrier algorithm for characteristic polynomial, determinant and inverse matrix
		 * 
		 * Computes the characteristic polynomial coefficients, determinant, and inverse of a square matrix.
		 * The inverse is only valid if the matrix is invertible (det != 0).
		 * 
		 * Characteristic polynomial: det(λI - A) = λ^n + c_{n-1}λ^{n-1} + ... + c_1λ + c_0
		 * where c_0 = (-1)^n * det(A)
		 */
		inline void FaddeevAlg(const Matrix<Real>& A, PolynomRealFunc& outCharPoly, Real& outDet, Matrix<Real>& outInv)
		{
			int N = A.rows();

			auto identMat = Matrix<Real>::Identity(N);
			Matrix<Real> nullMat(N, N);

			// Faddeev-Leverrier recurrence:
			// B_0 = I
			// c_{n-k} = -1/k * Tr(A * B_{k-1})  for k=1,...,n
			// B_k = A*B_{k-1} + c_{n-k}*I       for k=1,...,n
			// 
			// At the end: A^(-1) = -B_{n-1} / c_0 (when A is invertible)
			
			Vector<Matrix<Real>> B(N + 1, nullMat);
			outCharPoly.SetDegree(N);
			
			B[0] = identMat;
			outCharPoly[N] = 1.0;  // Leading coefficient is always 1
			
			for (int k = 1; k <= N; k++)
			{
				outCharPoly[N - k] = -1.0 / k * Trace(A * B[k - 1]);
				B[k] = A * B[k - 1] + outCharPoly[N - k] * identMat;
			}

			outDet = std::pow(-1.0, N) * outCharPoly[0];
			// A^(-1) = -B_{n-1} / c_0
			outInv = -B[N - 1] / outCharPoly[0];
		}

		// =========================================================================
		// RANK CALCULATION (Gaussian elimination - no external solver)
		// =========================================================================

		/**
		 * @brief Calculate the rank of a matrix using Gaussian elimination
		 * 
		 * Note: For numerical stability with ill-conditioned matrices,
		 * consider using MatrixAlg::Rank() which uses SVD.
		 */
		template<class Type>
		int RankGaussian(const Matrix<Type>& A, Real EPS = Defaults::RankAlgEPS)
		{
			int rows = A.rows();
			int cols = A.cols();
			Matrix<Type> mat = A;			// Make a copy to preserve the original

			int rank = 0;
			std::vector<bool> row_selected(rows, false);

			for (int col = 0; col < cols; ++col)
			{
				// Partial pivoting: find row with largest absolute value in this column
				int pivot_row = -1;
				Real best = 0;
				for (int row = 0; row < rows; ++row)
				{
					Real val = std::abs(mat(row, col));
					if (!row_selected[row] && val > EPS && val > best)
					{
						best = val;
						pivot_row = row;
					}
				}
				if (pivot_row == -1)
					continue;

				++rank;
				row_selected[pivot_row] = true;

				// Eliminate below
				for (int row = 0; row < rows; ++row)
				{
					if (row != pivot_row)
					{
						Type factor = mat(row, col) / mat(pivot_row, col);
						for (int k = col; k < cols; ++k)
							mat(row, k) -= factor * mat(pivot_row, k);
					}
				}
			}
			return rank;
		}

		// =========================================================================
		// DETERMINANT (uses LU decomposition)
		// =========================================================================

		/**
		 * @brief Calculate the determinant of a matrix using LU decomposition
		 * 
		 * Uses LU decomposition for efficient O(n³) computation.
		 * For singular matrices, returns 0 instead of throwing.
		 * 
		 * @param A Input matrix (must be square)
		 * @return Determinant value
		 */
		template<class Type>
		Type Det(const Matrix<Type>& A)
		{
			try
			{
				LUSolver<Type> lu(A);
				return lu.det();
			}
			catch(const SingularMatrixError&)
			{
				return Type(0);
			}
		}

		/**
		 * @brief Calculate the determinant of a fixed-size MatrixNM
		 * 
		 * Uses closed-form expressions for dimensions 1, 2, and 3 for maximum
		 * efficiency. For larger dimensions, converts to dynamic Matrix and uses
		 * LU decomposition.
		 * 
		 * @tparam Type Element type
		 * @tparam N    Matrix dimension (N×N)
		 * @param A     Input square MatrixNM
		 * @return      Determinant value
		 */
		template<class Type, int N>
		Type Det(const MatrixNM<Type, N, N>& A)
		{
			if constexpr (N == 1)
			{
				return A(0, 0);
			}
			else if constexpr (N == 2)
			{
				// |a b| = ad - bc
				// |c d|
				return A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0);
			}
			else if constexpr (N == 3)
			{
				// Sarrus' rule (cofactor expansion along first row)
				return A(0, 0) * (A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1))
				     - A(0, 1) * (A(1, 0) * A(2, 2) - A(1, 2) * A(2, 0))
				     + A(0, 2) * (A(1, 0) * A(2, 1) - A(1, 1) * A(2, 0));
			}
			else
			{
				// Convert to dynamic Matrix and use LU-based Det
				Matrix<Type> M(N, N);
				for (int i = 0; i < N; i++)
					for (int j = 0; j < N; j++)
						M(i, j) = A(i, j);
				return Det(M);
			}
		}

		// =========================================================================
		// SVD-BASED LINEAR ALGEBRA UTILITIES
		// =========================================================================

		/**
		 * @struct SVDResult
		 * @brief Result of SVD decomposition with linear algebra utilities
		 */
		struct SVDResult
		{
			Vector<Real> singularValues;  ///< Singular values in descending order
			Matrix<Real> U;               ///< Left singular vectors (m x n)
			Matrix<Real> V;               ///< Right singular vectors (n x n)
			int rank;                     ///< Numerical rank
			Real conditionNumber;         ///< Condition number (σ_max / σ_min)
			
			SVDResult() : rank(0), conditionNumber(0.0) {}
		};

		/**
		 * Compute SVD decomposition and derive linear algebra quantities.
		 * 
		 * Computes A = U * diag(σ) * V^T where:
		 * - U: m×n matrix of left singular vectors
		 * - σ: singular values in descending order
		 * - V: n×n matrix of right singular vectors
		 * 
		 * @param A Input matrix (m x n)
		 * @param tol Tolerance for rank determination
		 * @return SVDResult with singular values, U, V, rank, and condition number
		 */
		inline SVDResult ComputeSVD(const Matrix<Real>& A, Real tol = -1.0)
		{
			SVDecompositionSolver<Real> svd(A);
			SVDResult result;
			
			result.singularValues = svd.getW();
			result.U = svd.getU();
			result.V = svd.getV();
			result.rank = svd.Rank(tol);
			
			// Condition number = σ_max / σ_min
			int n = result.singularValues.size();
			if (n > 0 && result.singularValues[n-1] > 0)
				result.conditionNumber = result.singularValues[0] / result.singularValues[n-1];
			else
				result.conditionNumber = std::numeric_limits<Real>::infinity();
			
			return result;
		}

		/**
		 * Compute singular values of a matrix.
		 * 
		 * @param A Input matrix
		 * @return Vector of singular values in descending order
		 */
		inline Vector<Real> SingularValues(const Matrix<Real>& A)
		{
			SVDecompositionSolver<Real> svd(A);
			return svd.getW();
		}

		/**
		 * Compute numerical rank of a matrix using SVD.
		 * 
		 * The rank is the number of singular values above the threshold.
		 * This is more numerically stable than Gaussian elimination for
		 * ill-conditioned matrices.
		 * 
		 * @param A Input matrix
		 * @param tol Threshold for zero singular values (-1 for automatic)
		 * @return Numerical rank
		 */
		inline int RankSVD(const Matrix<Real>& A, Real tol = -1.0)
		{
			SVDecompositionSolver<Real> svd(A);
			return svd.Rank(tol);
		}

		/**
		 * Compute nullity (dimension of null space) of a matrix.
		 * 
		 * Nullity = n - rank, where n is the number of columns.
		 * 
		 * @param A Input matrix
		 * @param tol Threshold for zero singular values
		 * @return Nullity (dimension of kernel)
		 */
		inline int Nullity(const Matrix<Real>& A, Real tol = -1.0)
		{
			SVDecompositionSolver<Real> svd(A);
			return svd.Nullity(tol);
		}

		/**
		 * Compute condition number of a matrix.
		 * 
		 * The condition number κ(A) = σ_max / σ_min measures how sensitive
		 * the solution of Ax = b is to perturbations in A and b.
		 * 
		 * - κ ≈ 1: Well-conditioned
		 * - κ >> 1: Ill-conditioned
		 * - κ = ∞: Singular matrix
		 * 
		 * Rule of thumb: You lose log10(κ) decimal digits of precision
		 * when solving linear systems.
		 * 
		 * @param A Input matrix
		 * @return Condition number (infinity if singular)
		 */
		inline Real ConditionNumber(const Matrix<Real>& A)
		{
			SVDecompositionSolver<Real> svd(A);
			Vector<Real> w = svd.getW();
			int n = w.size();
			
			if (n == 0 || w[n-1] <= 0.0)
				return std::numeric_limits<Real>::infinity();
			
			return w[0] / w[n-1];
		}

		/**
		 * Compute orthonormal basis for the null space (kernel) of A.
		 * 
		 * The null space is {x : Ax = 0}. The columns of the returned matrix
		 * form an orthonormal basis for this space.
		 * 
		 * @param A Input matrix (m x n)
		 * @param tol Threshold for zero singular values
		 * @return Matrix whose columns are orthonormal basis of null(A)
		 *         Returns empty matrix (n x 0) if null space is trivial
		 */
		inline Matrix<Real> NullSpace(const Matrix<Real>& A, Real tol = -1.0)
		{
			SVDecompositionSolver<Real> svd(A);
			return svd.Nullspace(tol);
		}

		/**
		 * Alias for NullSpace - kernel is the same as null space
		 */
		inline Matrix<Real> Kernel(const Matrix<Real>& A, Real tol = -1.0)
		{
			return NullSpace(A, tol);
		}

		/**
		 * Compute orthonormal basis for the column space (range) of A.
		 * 
		 * The column space is {Ax : x ∈ R^n} = span of columns of A.
		 * The columns of the returned matrix form an orthonormal basis.
		 * 
		 * @param A Input matrix (m x n)
		 * @param tol Threshold for zero singular values
		 * @return Matrix whose columns are orthonormal basis of col(A)
		 */
		inline Matrix<Real> ColumnSpace(const Matrix<Real>& A, Real tol = -1.0)
		{
			SVDecompositionSolver<Real> svd(A);
			return svd.Range(tol);
		}

		/**
		 * Alias for ColumnSpace - range is the same as column space
		 */
		inline Matrix<Real> Range(const Matrix<Real>& A, Real tol = -1.0)
		{
			return ColumnSpace(A, tol);
		}

		/**
		 * Compute the row space of A.
		 * 
		 * The row space is the column space of A^T.
		 * 
		 * @param A Input matrix
		 * @param tol Threshold for zero singular values
		 * @return Matrix whose columns are orthonormal basis of row(A)
		 */
		inline Matrix<Real> RowSpace(const Matrix<Real>& A, Real tol = -1.0)
		{
			// Row space of A = Column space of A^T
			return ColumnSpace(A.transpose(), tol);
		}

		/**
		 * Compute the left null space of A.
		 * 
		 * The left null space is {y : y^T A = 0} = null(A^T).
		 * 
		 * @param A Input matrix
		 * @param tol Threshold for zero singular values
		 * @return Matrix whose columns are orthonormal basis of null(A^T)
		 */
		inline Matrix<Real> LeftNullSpace(const Matrix<Real>& A, Real tol = -1.0)
		{
			// Left null space = Null space of A^T
			return NullSpace(A.transpose(), tol);
		}

		/**
		 * @struct FundamentalSubspaces
		 * @brief The four fundamental subspaces of a matrix
		 * 
		 * For an m×n matrix A:
		 * - Column space: dim = r (rank)
		 * - Row space: dim = r
		 * - Null space: dim = n - r
		 * - Left null space: dim = m - r
		 * 
		 * Key relationships:
		 * - Column space ⊥ Left null space (in R^m)
		 * - Row space ⊥ Null space (in R^n)
		 */
		struct FundamentalSubspaces
		{
			Matrix<Real> columnSpace;    ///< Orthonormal basis for col(A)
			Matrix<Real> rowSpace;       ///< Orthonormal basis for row(A)
			Matrix<Real> nullSpace;      ///< Orthonormal basis for null(A)
			Matrix<Real> leftNullSpace;  ///< Orthonormal basis for null(A^T)
			int rank;                    ///< Rank of the matrix
		};

		/**
		 * Compute all four fundamental subspaces of a matrix.
		 * 
		 * @param A Input matrix (m x n)
		 * @param tol Threshold for zero singular values
		 * @return FundamentalSubspaces containing bases for all four spaces
		 */
		inline FundamentalSubspaces ComputeFundamentalSubspaces(const Matrix<Real>& A, Real tol = -1.0)
		{
			FundamentalSubspaces result;
			
			SVDecompositionSolver<Real> svd(A);
			result.columnSpace = svd.Range(tol);
			result.nullSpace = svd.Nullspace(tol);
			result.rank = svd.Rank(tol);
			
			// For row space and left null space, use SVD of A^T
			SVDecompositionSolver<Real> svdT(A.transpose());
			result.rowSpace = svdT.Range(tol);
			result.leftNullSpace = svdT.Nullspace(tol);
			
			return result;
		}

		/**
		 * Compute the pseudoinverse (Moore-Penrose inverse) of a matrix.
		 * 
		 * For A = U * diag(σ) * V^T, the pseudoinverse is:
		 * A⁺ = V * diag(1/σ) * U^T
		 * 
		 * Properties:
		 * - A * A⁺ * A = A
		 * - A⁺ * A * A⁺ = A⁺
		 * - For full-rank square A: A⁺ = A⁻¹
		 * - x = A⁺ * b is the least-squares solution to Ax = b
		 * 
		 * @param A Input matrix (m x n)
		 * @param tol Threshold for zero singular values
		 * @return Pseudoinverse A⁺ (n x m)
		 */
		inline Matrix<Real> PseudoInverse(const Matrix<Real>& A, Real tol = -1.0)
		{
			int m = A.rows();
			int n = A.cols();
			
			SVDecompositionSolver<Real> svd(A);
			Vector<Real> w = svd.getW();
			Matrix<Real> U = svd.getU();
			Matrix<Real> V = svd.getV();
			
			// Compute threshold
			Real eps = std::numeric_limits<Real>::epsilon();
			Real threshold = (tol >= 0.0) ? tol : 0.5 * std::sqrt(m + n + 1.0) * w[0] * eps;
			
			// Compute A⁺ = V * diag(1/σ) * U^T
			Matrix<Real> result(n, m);
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < m; j++)
				{
					Real sum = 0.0;
					for (int k = 0; k < n; k++)
					{
						if (w[k] > threshold)
							sum += V(i, k) * U(j, k) / w[k];
					}
					result(i, j) = sum;
				}
			}
			
			return result;
		}
	}
} // namespace MML

#endif




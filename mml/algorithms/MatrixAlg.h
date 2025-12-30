///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        MatrixAlg.h                                                         ///
///  Description: Matrix algebra algorithms and utilities                             ///
///               Decompositions, norms, definiteness, rank, condition number         ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_MATRIX_ALG_H
#define MML_MATRIX_ALG_H

#include "base/Vector.h"
#include "base/Matrix.h"
#include "base/MatrixSym.h"
#include "core/LinAlgEqSolvers.h"
#include "algorithms/EigenSystemSolvers.h"

#include <cmath>
#include <algorithm>
#include <limits>

namespace MML
{
  /**
   * @namespace MatrixAlg
   * @brief Collection of matrix algebra algorithms and utilities
   */
  namespace MatrixAlg
  {
    // =========================================================================
    // MATRIX PROPERTY CHECKS
    // =========================================================================

    /**
     * Check if matrix H is upper Hessenberg within tolerance.
     * 
     * Upper Hessenberg: H(i,j) = 0 for all i > j + 1
     * (all elements below the first subdiagonal are zero)
     * 
     * @param H Matrix to check
     * @param tol Tolerance for zero comparison
     * @return true if H is upper Hessenberg
     */
    inline bool IsUpperHessenberg(const Matrix<Real>& H, Real tol = 1e-10)
    {
      int n = H.RowNum();
      for (int i = 2; i < n; i++)
        for (int j = 0; j < i - 1; j++)
          if (std::abs(H(i, j)) > tol)
            return false;
      return true;
    }

    /**
     * Check if matrix Q is orthogonal: Q^T * Q = I
     * 
     * @param Q Matrix to check
     * @param tol Tolerance for identity comparison
     * @return true if Q is orthogonal
     */
    inline bool IsOrthogonal(const Matrix<Real>& Q, Real tol = 1e-10)
    {
      int n = Q.RowNum();
      for (int i = 0; i < n; i++)
      {
        for (int j = 0; j < n; j++)
        {
          Real dot = 0.0;
          for (int k = 0; k < n; k++)
            dot += Q(k, i) * Q(k, j);
          
          Real expected = (i == j) ? 1.0 : 0.0;
          if (std::abs(dot - expected) > tol)
            return false;
        }
      }
      return true;
    }

    /**
     * Check if matrix is upper triangular within tolerance.
     * 
     * @param A Matrix to check
     * @param tol Tolerance for zero comparison
     * @return true if A is upper triangular
     */
    inline bool IsUpperTriangular(const Matrix<Real>& A, Real tol = 1e-10)
    {
      int n = A.RowNum();
      for (int i = 1; i < n; i++)
        for (int j = 0; j < i; j++)
          if (std::abs(A(i, j)) > tol)
            return false;
      return true;
    }

    /**
     * Check if matrix is lower triangular within tolerance.
     * 
     * @param A Matrix to check
     * @param tol Tolerance for zero comparison
     * @return true if A is lower triangular
     */
    inline bool IsLowerTriangular(const Matrix<Real>& A, Real tol = 1e-10)
    {
      int n = A.RowNum();
      for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
          if (std::abs(A(i, j)) > tol)
            return false;
      return true;
    }

    /**
     * Check if matrix is diagonal within tolerance.
     * 
     * @param A Matrix to check
     * @param tol Tolerance for zero comparison
     * @return true if A is diagonal
     */
    inline bool IsDiagonal(const Matrix<Real>& A, Real tol = 1e-10)
    {
      int n = A.RowNum();
      for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
          if (i != j && std::abs(A(i, j)) > tol)
            return false;
      return true;
    }

    // =========================================================================
    // MATRIX SCALAR QUANTITIES
    // =========================================================================

    /**
     * Compute trace (sum of diagonal elements)
     * 
     * @param A Input matrix
     * @return Sum of diagonal elements
     */
    inline Real Trace(const Matrix<Real>& A)
    {
      Real sum = 0.0;
      int n = std::min(A.RowNum(), A.ColNum());
      for (int i = 0; i < n; i++)
        sum += A(i, i);
      return sum;
    }

    /**
     * Compute Frobenius norm: ||A||_F = sqrt(sum of all a_ij^2)
     * 
     * @param A Input matrix
     * @return Frobenius norm
     */
    inline Real FrobeniusNorm(const Matrix<Real>& A)
    {
      Real sum = 0.0;
      for (int i = 0; i < A.RowNum(); i++)
        for (int j = 0; j < A.ColNum(); j++)
          sum += A(i, j) * A(i, j);
      return std::sqrt(sum);
    }

    /**
     * Compute infinity norm (maximum absolute row sum): ||A||_∞
     * 
     * @param A Input matrix
     * @return Infinity norm
     */
    inline Real InfinityNorm(const Matrix<Real>& A)
    {
      Real maxRowSum = 0.0;
      for (int i = 0; i < A.RowNum(); i++)
      {
        Real rowSum = 0.0;
        for (int j = 0; j < A.ColNum(); j++)
          rowSum += std::abs(A(i, j));
        maxRowSum = std::max(maxRowSum, rowSum);
      }
      return maxRowSum;
    }

    /**
     * Compute 1-norm (maximum absolute column sum): ||A||_1
     * 
     * @param A Input matrix
     * @return 1-norm
     */
    inline Real OneNorm(const Matrix<Real>& A)
    {
      Real maxColSum = 0.0;
      for (int j = 0; j < A.ColNum(); j++)
      {
        Real colSum = 0.0;
        for (int i = 0; i < A.RowNum(); i++)
          colSum += std::abs(A(i, j));
        maxColSum = std::max(maxColSum, colSum);
      }
      return maxColSum;
    }

    /**
     * Compute maximum absolute difference between two matrices
     * 
     * @param A First matrix
     * @param B Second matrix
     * @return Maximum |A(i,j) - B(i,j)|
     */
    inline Real MaxAbsDiff(const Matrix<Real>& A, const Matrix<Real>& B)
    {
      Real maxDiff = 0.0;
      int n = A.RowNum();
      for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
          maxDiff = std::max(maxDiff, std::abs(A(i, j) - B(i, j)));
      return maxDiff;
    }

    // =========================================================================
    // MATRIX TRANSFORMATIONS
    // =========================================================================

    /**
     * Compute Q^T * A * Q (similarity transformation)
     * 
     * @param Q Transformation matrix
     * @param A Matrix to transform
     * @return Q^T * A * Q
     */
    inline Matrix<Real> SimilarityTransform(const Matrix<Real>& Q, const Matrix<Real>& A)
    {
      int n = A.RowNum();
      Matrix<Real> temp(n, n);
      Matrix<Real> result(n, n);
      
      // temp = A * Q
      for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
        {
          temp(i, j) = 0.0;
          for (int k = 0; k < n; k++)
            temp(i, j) += A(i, k) * Q(k, j);
        }
      
      // result = Q^T * temp
      for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
        {
          result(i, j) = 0.0;
          for (int k = 0; k < n; k++)
            result(i, j) += Q(k, i) * temp(k, j);
        }
      
      return result;
    }

    // =========================================================================
    // HESSENBERG REDUCTION
    // =========================================================================

    /**
     * @struct HessenbergResult
     * @brief Result of Hessenberg reduction: H = Q^T * A * Q
     */
    struct HessenbergResult
    {
      Matrix<Real> H;    // Upper Hessenberg matrix
      Matrix<Real> Q;    // Orthogonal transformation matrix
      
      HessenbergResult(int n) : H(n, n), Q(n, n) {}
    };

    /**
     * Reduce matrix A to upper Hessenberg form using Householder reflections.
     * 
     * MATHEMATICAL DEFINITION:
     * Find orthogonal Q such that H = Q^T * A * Q is upper Hessenberg,
     * meaning H(i,j) = 0 for all i > j + 1.
     * 
     * The transformation preserves eigenvalues (similarity transformation).
     * 
     * @param A Input matrix (n x n)
     * @return HessenbergResult containing H and Q where H = Q^T * A * Q
     * 
     * ALGORITHM (Householder):
     * For k = 0, 1, ..., n-3:
     *   1. Let x = A(k+1:n, k) be the subdiagonal column
     *   2. Compute Householder vector v to zero x(1:end)
     *   3. Apply P = I - 2*v*v^T from left: A = P*A
     *   4. Apply P from right: A = A*P (similarity)
     *   5. Accumulate Q = Q * P
     * 
     * Complexity: O(10n³/3) flops
     * 
     * REFERENCES:
     * - Golub & Van Loan, "Matrix Computations", 4th ed., Section 7.4
     */
    inline HessenbergResult ReduceToHessenberg(const Matrix<Real>& A)
    {
      int n = A.RowNum();
      HessenbergResult result(n);
      
      if (n <= 2)
      {
        result.H = A;
        result.Q = Matrix<Real>::GetUnitMatrix(n);
        return result;
      }
      
      // Copy A to H (we'll transform H in place)
      result.H = A;
      // Initialize Q as identity
      result.Q = Matrix<Real>::GetUnitMatrix(n);
      
      // Temporary storage for Householder vector
      Vector<Real> v(n);
      
      // For each column k, zero out elements below the subdiagonal
      for (int k = 0; k < n - 2; k++)
      {
        // Compute the norm of the column below the diagonal
        Real sigma = 0.0;
        for (int i = k + 1; i < n; i++)
          sigma += result.H(i, k) * result.H(i, k);
        sigma = std::sqrt(sigma);
        
        if (sigma < 1e-30)
          continue;  // Column already zero below subdiagonal
        
        // Choose sign to avoid cancellation: sign opposite to H(k+1, k)
        if (result.H(k + 1, k) > 0.0)
          sigma = -sigma;
        
        // Build Householder vector v = [0, ..., 0, H(k+1,k) - sigma, H(k+2,k), ..., H(n-1,k)]
        // But we only need elements k+1 to n-1
        Real h_k1_k_minus_sigma = result.H(k + 1, k) - sigma;
        
        // Compute ||v||^2
        Real vNormSq = h_k1_k_minus_sigma * h_k1_k_minus_sigma;
        for (int i = k + 2; i < n; i++)
          vNormSq += result.H(i, k) * result.H(i, k);
        
        if (vNormSq < 1e-30)
          continue;
        
        // beta = 2 / ||v||^2
        Real beta = 2.0 / vNormSq;
        
        // Store v in temporary array
        v[k + 1] = h_k1_k_minus_sigma;
        for (int i = k + 2; i < n; i++)
          v[i] = result.H(i, k);
        
        // Apply Householder from left: H = (I - beta*v*v^T) * H
        // H(i,j) -= beta * v(i) * sum_m(v(m) * H(m,j))
        for (int j = k; j < n; j++)
        {
          Real dot = 0.0;
          for (int i = k + 1; i < n; i++)
            dot += v[i] * result.H(i, j);
          
          for (int i = k + 1; i < n; i++)
            result.H(i, j) -= beta * v[i] * dot;
        }
        
        // Apply Householder from right: H = H * (I - beta*v*v^T)
        // H(i,j) -= beta * H(i,m) * v(m) * v(j)
        for (int i = 0; i < n; i++)
        {
          Real dot = 0.0;
          for (int j = k + 1; j < n; j++)
            dot += result.H(i, j) * v[j];
          
          for (int j = k + 1; j < n; j++)
            result.H(i, j) -= beta * dot * v[j];
        }
        
        // Accumulate Q: Q = Q * (I - beta*v*v^T)
        // Q(i,j) -= beta * Q(i,m) * v(m) * v(j)
        for (int i = 0; i < n; i++)
        {
          Real dot = 0.0;
          for (int j = k + 1; j < n; j++)
            dot += result.Q(i, j) * v[j];
          
          for (int j = k + 1; j < n; j++)
            result.Q(i, j) -= beta * dot * v[j];
        }
        
        // After applying transformation, H(k+1, k) should be sigma
        // and H(k+2:n-1, k) should be 0
        // Clean up the zeros explicitly
        result.H(k + 1, k) = sigma;
        for (int i = k + 2; i < n; i++)
          result.H(i, k) = 0.0;
      }
      
      return result;
    }

    // =========================================================================
    // MATRIX SYMMETRY CHECKS
    // =========================================================================

    /**
     * Check if matrix is symmetric within tolerance: A[i][j] ≈ A[j][i]
     * 
     * @param A Matrix to check
     * @param tol Tolerance for comparison
     * @return true if A is symmetric
     */
    inline bool IsSymmetric(const Matrix<Real>& A, Real tol = 1e-10)
    {
      int n = A.RowNum();
      if (A.ColNum() != n) return false;  // Must be square
      
      for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
          if (std::abs(A(i, j) - A(j, i)) > tol)
            return false;
      return true;
    }

    /**
     * Check if matrix is skew-symmetric: A[i][j] ≈ -A[j][i]
     * 
     * @param A Matrix to check
     * @param tol Tolerance for comparison
     * @return true if A is skew-symmetric
     */
    inline bool IsSkewSymmetric(const Matrix<Real>& A, Real tol = 1e-10)
    {
      int n = A.RowNum();
      if (A.ColNum() != n) return false;  // Must be square
      
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
    // MATRIX DEFINITENESS CHECKS
    // =========================================================================
    
    /**
     * @enum Definiteness
     * @brief Classification of matrix definiteness
     */
    enum class Definiteness
    {
      PositiveDefinite,      ///< All eigenvalues > 0
      PositiveSemiDefinite,  ///< All eigenvalues >= 0, at least one = 0
      NegativeDefinite,      ///< All eigenvalues < 0
      NegativeSemiDefinite,  ///< All eigenvalues <= 0, at least one = 0
      Indefinite            ///< Has both positive and negative eigenvalues
    };

    /**
     * Classify the definiteness of a symmetric matrix.
     * 
     * MATHEMATICAL DEFINITION:
     * - Positive definite: x^T A x > 0 for all x ≠ 0 (all eigenvalues > 0)
     * - Positive semi-definite: x^T A x >= 0 for all x (eigenvalues >= 0)
     * - Negative definite: x^T A x < 0 for all x ≠ 0 (all eigenvalues < 0)
     * - Negative semi-definite: x^T A x <= 0 for all x (eigenvalues <= 0)
     * - Indefinite: Has both positive and negative eigenvalues
     * 
     * @param A Symmetric matrix (only upper triangle used)
     * @param tol Tolerance for zero eigenvalue detection
     * @return Definiteness classification
     * 
     * @note Matrix must be symmetric. For non-symmetric matrices, 
     *       consider using (A + A^T)/2.
     */
    inline Definiteness ClassifyDefiniteness(const MatrixSym<Real>& A, Real tol = 1e-10)
    {
      // Compute eigenvalues using Jacobi method
      auto result = SymmMatEigenSolverJacobi::Solve(A, tol);
      
      int n = A.RowNum();
      int numPositive = 0;
      int numNegative = 0;
      int numZero = 0;
      
      for (int i = 0; i < n; i++)
      {
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

    /**
     * Classify definiteness of a general matrix (uses symmetric part)
     */
    inline Definiteness ClassifyDefiniteness(const Matrix<Real>& A, Real tol = 1e-10)
    {
      int n = A.RowNum();
      if (A.ColNum() != n)
        throw MatrixDimensionError("ClassifyDefiniteness - matrix must be square", n, A.ColNum(), -1, -1);
      
      // Convert to symmetric matrix (A + A^T) / 2
      MatrixSym<Real> symA(n);
      for (int i = 0; i < n; i++)
        for (int j = i; j < n; j++)
          symA(i, j) = 0.5 * (A(i, j) + A(j, i));
      
      return ClassifyDefiniteness(symA, tol);
    }

    /**
     * Check if symmetric matrix is positive definite.
     * 
     * A matrix is positive definite if x^T A x > 0 for all x ≠ 0,
     * which is equivalent to all eigenvalues being positive.
     * 
     * @param A Symmetric matrix
     * @param tol Tolerance for eigenvalue comparison
     * @return true if A is positive definite
     */
    inline bool IsPositiveDefinite(const MatrixSym<Real>& A, Real tol = 1e-10)
    {
      return ClassifyDefiniteness(A, tol) == Definiteness::PositiveDefinite;
    }

    /**
     * Check if general matrix is positive definite (uses symmetric part)
     */
    inline bool IsPositiveDefinite(const Matrix<Real>& A, Real tol = 1e-10)
    {
      return ClassifyDefiniteness(A, tol) == Definiteness::PositiveDefinite;
    }

    /**
     * Check if symmetric matrix is positive semi-definite.
     * 
     * A matrix is positive semi-definite if x^T A x >= 0 for all x,
     * which is equivalent to all eigenvalues being non-negative.
     * 
     * @param A Symmetric matrix
     * @param tol Tolerance for eigenvalue comparison
     * @return true if A is positive semi-definite (includes positive definite)
     */
    inline bool IsPositiveSemiDefinite(const MatrixSym<Real>& A, Real tol = 1e-10)
    {
      auto def = ClassifyDefiniteness(A, tol);
      return def == Definiteness::PositiveDefinite || def == Definiteness::PositiveSemiDefinite;
    }

    inline bool IsPositiveSemiDefinite(const Matrix<Real>& A, Real tol = 1e-10)
    {
      auto def = ClassifyDefiniteness(A, tol);
      return def == Definiteness::PositiveDefinite || def == Definiteness::PositiveSemiDefinite;
    }

    /**
     * Check if symmetric matrix is negative definite.
     * 
     * A matrix is negative definite if x^T A x < 0 for all x ≠ 0,
     * which is equivalent to all eigenvalues being negative.
     * 
     * @param A Symmetric matrix
     * @param tol Tolerance for eigenvalue comparison
     * @return true if A is negative definite
     */
    inline bool IsNegativeDefinite(const MatrixSym<Real>& A, Real tol = 1e-10)
    {
      return ClassifyDefiniteness(A, tol) == Definiteness::NegativeDefinite;
    }

    inline bool IsNegativeDefinite(const Matrix<Real>& A, Real tol = 1e-10)
    {
      return ClassifyDefiniteness(A, tol) == Definiteness::NegativeDefinite;
    }

    /**
     * Check if symmetric matrix is negative semi-definite.
     * 
     * @param A Symmetric matrix
     * @param tol Tolerance for eigenvalue comparison
     * @return true if A is negative semi-definite (includes negative definite)
     */
    inline bool IsNegativeSemiDefinite(const MatrixSym<Real>& A, Real tol = 1e-10)
    {
      auto def = ClassifyDefiniteness(A, tol);
      return def == Definiteness::NegativeDefinite || def == Definiteness::NegativeSemiDefinite;
    }

    inline bool IsNegativeSemiDefinite(const Matrix<Real>& A, Real tol = 1e-10)
    {
      auto def = ClassifyDefiniteness(A, tol);
      return def == Definiteness::NegativeDefinite || def == Definiteness::NegativeSemiDefinite;
    }

    /**
     * Check if matrix is indefinite (has both positive and negative eigenvalues).
     * 
     * @param A Symmetric matrix
     * @param tol Tolerance for eigenvalue comparison
     * @return true if A is indefinite
     */
    inline bool IsIndefinite(const MatrixSym<Real>& A, Real tol = 1e-10)
    {
      return ClassifyDefiniteness(A, tol) == Definiteness::Indefinite;
    }

    inline bool IsIndefinite(const Matrix<Real>& A, Real tol = 1e-10)
    {
      return ClassifyDefiniteness(A, tol) == Definiteness::Indefinite;
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
      SVDecompositionSolver svd(A);
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
      SVDecompositionSolver svd(A);
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
    inline int Rank(const Matrix<Real>& A, Real tol = -1.0)
    {
      SVDecompositionSolver svd(A);
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
      SVDecompositionSolver svd(A);
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
      SVDecompositionSolver svd(A);
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
      SVDecompositionSolver svd(A);
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
      SVDecompositionSolver svd(A);
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
      return ColumnSpace(A.GetTranspose(), tol);
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
      return NullSpace(A.GetTranspose(), tol);
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
      
      SVDecompositionSolver svd(A);
      result.columnSpace = svd.Range(tol);
      result.nullSpace = svd.Nullspace(tol);
      result.rank = svd.Rank(tol);
      
      // For row space and left null space, use SVD of A^T
      SVDecompositionSolver svdT(A.GetTranspose());
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
      int m = A.RowNum();
      int n = A.ColNum();
      
      SVDecompositionSolver svd(A);
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

  } // namespace MatrixAlg

} // namespace MML

#endif // MML_MATRIX_ALG_H
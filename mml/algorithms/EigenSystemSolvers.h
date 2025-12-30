///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        EigenSystemSolvers.h                                                ///
///  Description: Eigenvalue/eigenvector solvers (Symmetric, General, Power Method)   ///
///               QR algorithm, Jacobi iteration, inverse iteration                   ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined( MML_EIGENSYSTEM_SOLVERS_H )
#define MML_EIGENSYSTEM_SOLVERS_H

#include "MMLBase.h"
#include "base/Vector.h"
#include "base/Matrix.h"
#include "base/MatrixSym.h"
#include "algorithms/EigenSolverHelpers.h"  // For building blocks used by EigenSolver

namespace MML
{
	/***************************************************************************************************
	 * JACOBI EIGENSOLVER FOR SYMMETRIC MATRICES
	 * 
	 * Classical iterative method using Givens rotations to diagonalize symmetric matrices.
	 * 
	 * ALGORITHM:
	 * 1. Find largest off-diagonal element a_pq
	 * 2. Compute Givens rotation to zero out a_pq
	 * 3. Apply rotation: A' = J^T * A * J
	 * 4. Accumulate rotations to build eigenvector matrix
	 * 5. Repeat until off-diagonal norm < tolerance
	 * 
	 * COMPLEXITY: O(n³) per sweep, typically 5-10 sweeps needed
	 * 
	 * PROS:
	 * - Simple, robust algorithm
	 * - Automatically computes all eigenvalues and eigenvectors
	 * - Excellent for small to medium matrices (n < 100)
	 * - Numerically very stable
	 * 
	 * CONS:
	 * - Slower than QR method for large matrices
	 * - O(n³) per sweep vs O(n³) total for QR
	 * 
	 * REFERENCES:
	 * - Golub & Van Loan, "Matrix Computations", Section 8.5
	 * - Numerical Recipes, Section 11.1
	 ***************************************************************************************************/
	class SymmMatEigenSolverJacobi
  {
	public:
		// Result structure for eigenvalue decomposition
		struct Result
		{
			Vector<Real> eigenvalues;      // Eigenvalues in ascending order
			Matrix<Real> eigenvectors;     // Column i = eigenvector for eigenvalue i
			int iterations;                // Number of sweeps performed
			Real residual;                 // Final off-diagonal norm
			bool converged;                // True if converged within tolerance

			Result(int n) : eigenvalues(n), eigenvectors(n, n), iterations(0), residual(0.0), converged(false) {}
		};

		/**
		 * Solve symmetric eigenvalue problem: A*v = λ*v
		 * 
		 * @param A        Symmetric matrix (only upper triangle is used)
		 * @param tol      Convergence tolerance (default: 1e-10)
		 * @param maxIter  Maximum number of sweeps (default: 100)
		 * @return Result containing eigenvalues, eigenvectors, and convergence info
		 * 
		 * POSTCONDITIONS:
		 * - eigenvalues[i] <= eigenvalues[i+1] (ascending order)
		 * - eigenvectors.Column(i) corresponds to eigenvalues[i]
		 * - Eigenvectors are orthonormal: V^T * V = I
		 * - A * V = V * diag(λ)
		 */
		static Result Solve(const MatrixSym<Real>& A, Real tol = 1e-10, int maxIter = 100)
		{
			int n = A.RowNum();
			Result result(n);

			// Initialize: Copy A to working matrix, V to identity
			Matrix<Real> D(n, n);
			for (int i = 0; i < n; i++)
				for (int j = 0; j < n; j++)
					D[i][j] = (i <= j) ? A(i, j) : A(j, i);

			Matrix<Real> V = Matrix<Real>::GetUnitMatrix(n);

			// Jacobi iterations
			int sweep = 0;
			Real offDiagNorm = OffDiagonalNorm(D, n);

			while (offDiagNorm > tol && sweep < maxIter)
			{
				// One sweep: process all off-diagonal elements
				for (int p = 0; p < n - 1; p++)
				{
					for (int q = p + 1; q < n; q++)
					{
					if (std::abs(D[p][q]) < PrecisionValues<Real>::EigenSolverZeroThreshold)  // Skip if already zero
							continue;
						// Compute Givens rotation parameters
						Real theta = ComputeRotationAngle(D, p, q);
						Real c = std::cos(theta);
						Real s = std::sin(theta);

						// Apply rotation
						ApplyRotation(D, V, p, q, c, s, n);
					}
				}

				sweep++;
				offDiagNorm = OffDiagonalNorm(D, n);
			}

			// Extract eigenvalues from diagonal
			for (int i = 0; i < n; i++)
				result.eigenvalues[i] = D[i][i];

			// Copy eigenvectors
			result.eigenvectors = V;

			// Sort eigenvalues and eigenvectors
			SortEigenvalues(result.eigenvalues, result.eigenvectors);

			// Set convergence info
			result.iterations = sweep;
			result.residual = offDiagNorm;
			result.converged = (offDiagNorm <= tol);

			return result;
		}

		/**
		 * Solve for a regular Matrix (extracts symmetric part)
		 * Note: Only use if you're certain the matrix is symmetric!
		 */
		static Result Solve(const Matrix<Real>& A, Real tol = 1e-10, int maxIter = 100)
		{
			// Convert to symmetric matrix (average A and A^T)
			int n = A.RowNum();
			MatrixSym<Real> symA(n);
			for (int i = 0; i < n; i++)
				for (int j = i; j < n; j++)
					symA(i, j) = 0.5 * (A[i][j] + A[j][i]);
			
			return Solve(symA, tol, maxIter);
		}

	private:
		// Compute rotation angle to zero out A[p][q]
		static Real ComputeRotationAngle(const Matrix<Real>& A, int p, int q)
		{
			Real theta;
			
			if (std::abs(A[p][q]) < PrecisionValues<Real>::EigenSolverZeroThreshold)
			{
				theta = 0.0;
			}
			else
			{
				Real tau = (A[q][q] - A[p][p]) / (2.0 * A[p][q]);
				Real t;
				
				// Choose sign to avoid loss of significance
				if (tau >= 0)
					t = 1.0 / (tau + std::sqrt(1.0 + tau * tau));
				else
					t = -1.0 / (-tau + std::sqrt(1.0 + tau * tau));
				
				theta = std::atan(t);
			}

			return theta;
		}
		
		// Apply Givens rotation J(p,q,θ) to A: A' = J^T * A * J
		static void ApplyRotation(Matrix<Real>& A, Matrix<Real>& V, int p, int q, Real c, Real s, int n)
		{
			// Update A
			for (int i = 0; i < n; i++)
			{
				if (i != p && i != q)
				{
					Real Aip = A[i][p];
					Real Aiq = A[i][q];
					A[i][p] = c * Aip - s * Aiq;
					A[p][i] = A[i][p];
					A[i][q] = s * Aip + c * Aiq;
					A[q][i] = A[i][q];
				}
			}

			// Update diagonal elements
			Real App = A[p][p];
			Real Aqq = A[q][q];
			Real Apq = A[p][q];

			A[p][p] = c * c * App - 2.0 * s * c * Apq + s * s * Aqq;
			A[q][q] = s * s * App + 2.0 * s * c * Apq + c * c * Aqq;
			A[p][q] = 0.0;
			A[q][p] = 0.0;

			// Accumulate rotation in eigenvector matrix V
			for (int i = 0; i < n; i++)
			{
				Real Vip = V[i][p];
				Real Viq = V[i][q];
				V[i][p] = c * Vip - s * Viq;
				V[i][q] = s * Vip + c * Viq;
			}
		}
		
		// Compute off-diagonal norm: sqrt(sum of squares of off-diagonal elements)
		static Real OffDiagonalNorm(const Matrix<Real>& A, int n)
		{
			Real sum = 0.0;
			for (int i = 0; i < n - 1; i++)
				for (int j = i + 1; j < n; j++)
					sum += A[i][j] * A[i][j];
			
			return std::sqrt(sum);
		}
		
		// Sort eigenvalues and eigenvectors in ascending order
		static void SortEigenvalues(Vector<Real>& eigenvalues, Matrix<Real>& eigenvectors)
		{
			int n = eigenvalues.size();
			
			// Simple selection sort (eigenvalues should be small n)
			for (int i = 0; i < n - 1; i++)
			{
				int minIdx = i;
				for (int j = i + 1; j < n; j++)
				{
					if (eigenvalues[j] < eigenvalues[minIdx])
						minIdx = j;
				}
				
				if (minIdx != i)
				{
					// Swap eigenvalues
					std::swap(eigenvalues[i], eigenvalues[minIdx]);
					
					// Swap corresponding eigenvectors (columns)
					for (int row = 0; row < n; row++)
						std::swap(eigenvectors[row][i], eigenvectors[row][minIdx]);
				}
			}
		}
  };

  /**
   * @class SymmMatEigenSolverQR
   * @brief QR algorithm for symmetric matrix eigenvalue decomposition
   *
   * Algorithm: Tridiagonal reduction + Implicit QR iteration with Wilkinson shift
   * 1. Reduce symmetric matrix to tridiagonal form using Householder reflections
   * 2. Apply implicit QR algorithm with Wilkinson shift to tridiagonal matrix
   * 3. Accumulate transformations to obtain eigenvectors
   *
   * Complexity: O(n³) for reduction + O(n²) per QR iteration (typically O(n) iterations)
   *             Total: O(n³) - more efficient than Jacobi for large matrices
   *
   * Pros: - Faster convergence than Jacobi (especially for larger matrices)
   *       - Cubic convergence rate with Wilkinson shift
   *       - Industry-standard algorithm
   * Cons: - More complex implementation
   *       - Requires careful handling of deflation
   * 
   * REFERENCES:
   * - Golub & Van Loan, "Matrix Computations", 4th ed., Sections 8.3-8.4
   * - Numerical Recipes, 3rd ed., Section 11.3
   * - Wilkinson, "The Algebraic Eigenvalue Problem"
   */
  class SymmMatEigenSolverQR
  {
  public:
    struct Result
    {
      Vector<Real>  eigenvalues;     // Sorted eigenvalues (ascending)
      Matrix<Real>  eigenvectors;    // Corresponding eigenvectors (columns)
      int           iterations;      // Number of QR iterations performed
      bool          converged;       // True if algorithm converged
      Real          residual;        // Final off-diagonal norm

      Result(int n) : eigenvalues(n), eigenvectors(n, n), iterations(0), converged(false), residual(0.0) {}
    };

    // Main interface: Solve symmetric eigenvalue problem
    static Result Solve(const MatrixSym<Real>& A, Real tol = 1e-10, int maxIter = 1000)
    {
      int n = A.RowNum();
      Result result(n);

      if (n == 0)
      {
        result.converged = true;
        return result;
      }

      if (n == 1)
      {
        result.eigenvalues[0] = A(0, 0);
        result.eigenvectors(0, 0) = 1.0;
        result.converged = true;
        return result;
      }

      // Step 1: Copy to full matrix for working storage
      Matrix<Real> T(n, n);
      for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
          T(i, j) = (i <= j) ? A(i, j) : A(j, i);

      // Step 2: Initialize eigenvector matrix to identity
      Matrix<Real> Q = Matrix<Real>::GetUnitMatrix(n);

      // Step 3: Reduce to tridiagonal form using Householder reflections
      // Store diagonal in diag, subdiagonal in subdiag
      Vector<Real> diag(n);
      Vector<Real> subdiag(n);
      
      TridiagonalReduction(T, Q, n);

      // Extract tridiagonal elements
      for (int i = 0; i < n; i++)
        diag[i] = T(i, i);
      for (int i = 0; i < n - 1; i++)
        subdiag[i] = T(i + 1, i);
      subdiag[n - 1] = 0.0;

      // Step 4: Apply implicit QR algorithm with Wilkinson shift
      ImplicitQRAlgorithm(diag, subdiag, Q, n, tol, maxIter,
                          result.iterations, result.converged);

      // Step 5: Compute final residual (off-diagonal norm)
      result.residual = 0.0;
      for (int i = 0; i < n - 1; i++)
        result.residual += subdiag[i] * subdiag[i];
      result.residual = std::sqrt(result.residual);

      // Step 6: Copy results
      result.eigenvalues = diag;
      result.eigenvectors = Q;

      // Step 7: Sort eigenvalues and eigenvectors in ascending order
      SortEigenpairs(result.eigenvalues, result.eigenvectors, n);

      return result;
    }
    
    // Overload for regular Matrix (will symmetrize)
    static Result Solve(const Matrix<Real>& A, Real tol = 1e-10, int maxIter = 1000)
    {
      int n = A.RowNum();
      MatrixSym<Real> symA(n);
      for (int i = 0; i < n; i++)
        for (int j = i; j < n; j++)
          symA(i, j) = 0.5 * (A[i][j] + A[j][i]);
      
      return Solve(symA, tol, maxIter);
    }

  private:
    /**
     * Reduce symmetric matrix to tridiagonal form using Householder reflections.
     * 
     * Algorithm: For k = 0, 1, ..., n-3:
     *   1. Compute Householder vector v to zero out T(k+2:n, k)
     *   2. Apply transformation: T = H * T * H where H = I - 2*v*v^T
     *   3. Accumulate: Q = Q * H
     * 
     * Result: T becomes tridiagonal, Q accumulates the orthogonal transformations
     * 
     * Based on Golub & Van Loan, "Matrix Computations", Algorithm 8.3.1
     */
    static void TridiagonalReduction(Matrix<Real>& T, Matrix<Real>& Q, int n)
    {
      for (int k = 0; k < n - 2; k++)
      {
        // Compute ||x||^2 where x = T(k+1:n-1, k)
        Real sigma = 0.0;
        for (int i = k + 1; i < n; i++)
          sigma += T(i, k) * T(i, k);

        if (sigma < PrecisionValues<Real>::DivisionSafetyThreshold)  // Column already zeroed
          continue;

        Real xNorm = std::sqrt(sigma);
        
        // alpha = -sign(T[k+1,k]) * ||x||
        // Choose sign to maximize |u[k+1]| for numerical stability
        Real alpha = (T(k + 1, k) >= 0.0) ? -xNorm : xNorm;
        
        // Householder vector: u = x - alpha*e1, then v = u/||u||
        Real u_k1 = T(k + 1, k) - alpha;
        
        // ||u||^2 = sigma + alpha^2 - 2*alpha*T[k+1,k]
        Real uNormSq = sigma + alpha * alpha - 2.0 * alpha * T(k + 1, k);
        Real uNorm = std::sqrt(uNormSq);
        
        if (uNorm < PrecisionValues<Real>::DivisionSafetyThreshold)
          continue;

        // Build normalized Householder vector v
        Vector<Real> v(n, 0.0);
        v[k + 1] = u_k1 / uNorm;
        for (int i = k + 2; i < n; i++)
          v[i] = T(i, k) / uNorm;

        // Apply H = I - 2*v*v^T to T from left and right
        // For symmetric matrix: T' = T - v*w^T - w*v^T
        // where w = 2*(T*v - (v^T*T*v)*v)
        
        // Compute p = T * v (restricted to rows/cols k+1 to n-1)
        Vector<Real> p(n, 0.0);
        for (int i = k + 1; i < n; i++)
          for (int j = k + 1; j < n; j++)
            p[i] += T(i, j) * v[j];

        // K = v^T * p
        Real K = 0.0;
        for (int i = k + 1; i < n; i++)
          K += v[i] * p[i];

        // w = 2*(p - K*v)
        Vector<Real> w(n, 0.0);
        for (int i = k + 1; i < n; i++)
          w[i] = 2.0 * (p[i] - K * v[i]);

        // Update T: T' = T - v*w^T - w*v^T
        for (int i = k + 1; i < n; i++)
          for (int j = k + 1; j < n; j++)
            T(i, j) -= v[i] * w[j] + w[i] * v[j];

        // Set subdiagonal explicitly (result of Householder on column k)
        T(k + 1, k) = alpha;
        T(k, k + 1) = alpha;
        
        // Zero elements below subdiagonal
        for (int i = k + 2; i < n; i++)
        {
          T(i, k) = 0.0;
          T(k, i) = 0.0;
        }

        // Accumulate transformation: Q = Q * H = Q * (I - 2*v*v^T)
        // Q_new[:, j] = Q[:, j] - 2*(Q*v)*v[j]
        for (int i = 0; i < n; i++)
        {
          Real dot = 0.0;
          for (int m = k + 1; m < n; m++)
            dot += Q(i, m) * v[m];
          
          for (int j = k + 1; j < n; j++)
            Q(i, j) -= 2.0 * dot * v[j];
        }
      }
    }

    /**
     * Implicit QR algorithm with Wilkinson shift for tridiagonal matrices.
     * 
     * This is the workhorse of the eigenvalue computation. Uses:
     * - Wilkinson shift for cubic convergence
     * - Implicit QR step via bulge chasing with Givens rotations
     * - Deflation when subdiagonal elements become negligible
     */
    static void ImplicitQRAlgorithm(Vector<Real>& d, Vector<Real>& e, Matrix<Real>& Q,
                                     int n, Real tol, int maxIter,
                                     int& iterations, bool& converged)
    {
      iterations = 0;
      converged = false;

      // Work on unreduced portion [lo, hi]
      int lo = 0;
      int hi = n - 1;

      while (hi > 0 && iterations < maxIter)
      {
        // Deflate: Find the largest unreduced block
        // Check if e[hi-1] is negligible
        Real threshold = tol * (std::abs(d[hi - 1]) + std::abs(d[hi]));
        if (threshold < tol) threshold = tol;
        
        if (std::abs(e[hi - 1]) <= threshold)
        {
          e[hi - 1] = 0.0;
          hi--;
          continue;
        }

        // Find lo: start of the unreduced block ending at hi
        lo = hi - 1;
        while (lo > 0)
        {
          threshold = tol * (std::abs(d[lo - 1]) + std::abs(d[lo]));
          if (threshold < tol) threshold = tol;
          
          if (std::abs(e[lo - 1]) <= threshold)
          {
            e[lo - 1] = 0.0;
            break;
          }
          lo--;
        }

        // Now work on the unreduced block [lo, hi]
        
        // Compute Wilkinson shift (eigenvalue of trailing 2x2 closer to d[hi])
        Real shift = WilkinsonShift(d[hi - 1], e[hi - 1], d[hi]);

        // Implicit QR step with shift (bulge chasing)
        ImplicitQRStep(d, e, Q, lo, hi, shift, n);

        iterations++;
      }

      converged = (hi <= 0) || IsConverged(e, n, tol);
    }

    /**
     * Compute Wilkinson shift: eigenvalue of 2x2 trailing submatrix closer to d[hi]
     * 
     * For matrix [[a, b], [b, c]], eigenvalues are:
     *   λ = (a+c)/2 ± sqrt(((a-c)/2)^2 + b^2)
     * 
     * We want the one closer to c (= d[hi])
     */
    static Real WilkinsonShift(Real a, Real b, Real c)
    {
      Real d = (a - c) * 0.5;
      Real b_sq = b * b;
      
      if (std::abs(d) < PrecisionValues<Real>::DivisionSafetyThreshold)
        return c - std::abs(b);

      Real sign_d = (d >= 0.0) ? 1.0 : -1.0;
      Real r = std::sqrt(d * d + b_sq);
      
      // Shift = c - b^2 / (d + sign(d)*sqrt(d^2+b^2))
      return c - b_sq / (d + sign_d * r);
    }

    /**
     * Implicit QR step with shift using Givens rotations (bulge chasing).
     * 
     * Creates a bulge at position (1,0) by introducing shift, then chases
     * it down the diagonal using Givens rotations.
     * 
     * For Givens G = [c s; -s c], the similarity transformation G*T*G^T
     * on the 2x2 block [d_k e_k; e_k d_{k+1}] gives:
     *   d'_k     = c²*d_k + 2*c*s*e_k + s²*d_{k+1}
     *   d'_{k+1} = s²*d_k - 2*c*s*e_k + c²*d_{k+1}
     *   e'_k     = c*s*(d_{k+1} - d_k) + (c² - s²)*e_k
     */
    static void ImplicitQRStep(Vector<Real>& d, Vector<Real>& e, Matrix<Real>& Q,
                                int lo, int hi, Real shift, int n)
    {
      // Initial bulge: the first rotation is applied to [d[lo]-shift; e[lo]]
      Real x = d[lo] - shift;
      Real z = e[lo];

      for (int k = lo; k < hi; k++)
      {
        // Compute Givens rotation to eliminate z
        // G = [c s; -s c] such that G * [x; z] = [r; 0]
        Real c, s;
        ComputeGivens(x, z, c, s);

        // Update the previous off-diagonal (if not first iteration)
        if (k > lo)
        {
          // The element e[k-1] becomes sqrt(x² + z²) after rotation
          e[k - 1] = std::sqrt(x * x + z * z);
        }

        // Save current values
        Real d_k = d[k];
        Real d_k1 = d[k + 1];
        Real e_k = e[k];

        // Apply similarity transformation G * [d_k e_k; e_k d_{k+1}] * G^T
        // d'_k     = c²*d_k + 2*c*s*e_k + s²*d_{k+1}
        // d'_{k+1} = s²*d_k - 2*c*s*e_k + c²*d_{k+1}
        // e'_k     = c*s*(d_{k+1} - d_k) + (c² - s²)*e_k
        d[k] = c * c * d_k + 2.0 * c * s * e_k + s * s * d_k1;
        d[k + 1] = s * s * d_k - 2.0 * c * s * e_k + c * c * d_k1;
        e[k] = c * s * (d_k1 - d_k) + (c * c - s * s) * e_k;

        // Bulge chasing: the rotation creates a fill-in at position (k+2, k)
        if (k < hi - 1)
        {
          // The bulge element is s * e[k+1] (from applying G to column k+1)
          x = e[k];
          z = s * e[k + 1];
          e[k + 1] = c * e[k + 1];
        }

        // Accumulate rotation in eigenvector matrix
        // We're computing Q * G^T (to get eigenvectors of original matrix)
        // Q_new[:, k] = c*Q[:, k] + s*Q[:, k+1]
        // Q_new[:, k+1] = -s*Q[:, k] + c*Q[:, k+1]
        for (int i = 0; i < n; i++)
        {
          Real qik = Q(i, k);
          Real qik1 = Q(i, k + 1);
          Q(i, k) = c * qik + s * qik1;
          Q(i, k + 1) = -s * qik + c * qik1;
        }
      }
    }

    /**
     * Compute Givens rotation coefficients c and s such that:
     *   [c  s] [a]   [r]
     *   [-s c] [b] = [0]
     * 
     * where r = sqrt(a^2 + b^2) (or -sqrt if needed for sign consistency)
     * 
     * Standard formulas: c = a/r, s = b/r
     * This zeros out b: c*a + s*b = (a^2 + b^2)/r = r
     *                  -s*a + c*b = (-ab + ab)/r = 0
     */
    static void ComputeGivens(Real a, Real b, Real& c, Real& s)
    {
      if (std::abs(b) < PrecisionValues<Real>::DivisionSafetyThreshold)
      {
        c = (a >= 0.0) ? 1.0 : -1.0;
        s = 0.0;
      }
      else if (std::abs(a) < PrecisionValues<Real>::DivisionSafetyThreshold)
      {
        c = 0.0;
        s = (b >= 0.0) ? 1.0 : -1.0;
      }
      else if (std::abs(b) > std::abs(a))
      {
        // Use t = a/b to avoid overflow
        Real t = a / b;
        Real u = std::sqrt(1.0 + t * t);
        if (b < 0.0) u = -u;
        s = 1.0 / u;
        c = s * t;
      }
      else
      {
        // Use t = b/a to avoid overflow
        Real t = b / a;
        Real u = std::sqrt(1.0 + t * t);
        if (a < 0.0) u = -u;
        c = 1.0 / u;
        s = c * t;
      }
    }

    /**
     * Check if the tridiagonal matrix has converged (all subdiagonals negligible)
     */
    static bool IsConverged(const Vector<Real>& e, int n, Real tol)
    {
      for (int i = 0; i < n - 1; i++)
      {
        if (std::abs(e[i]) > tol)
          return false;
      }
      return true;
    }

    /**
     * Sort eigenvalues and corresponding eigenvectors in ascending order
     */
    static void SortEigenpairs(Vector<Real>& eigenvalues, Matrix<Real>& eigenvectors, int n)
    {
      // Selection sort (n is typically small)
      for (int i = 0; i < n - 1; i++)
      {
        int minIdx = i;
        for (int j = i + 1; j < n; j++)
        {
          if (eigenvalues[j] < eigenvalues[minIdx])
            minIdx = j;
        }
        
        if (minIdx != i)
        {
          // Swap eigenvalues
          std::swap(eigenvalues[i], eigenvalues[minIdx]);
          
          // Swap corresponding eigenvectors (columns)
          for (int row = 0; row < n; row++)
            std::swap(eigenvectors(row, i), eigenvectors(row, minIdx));
        }
      }
    }
  };

  /***************************************************************************************************
   * GENERAL (NON-SYMMETRIC) EIGENSOLVER
   * 
   * Complete eigensolver for general (non-symmetric) matrices using QR iteration
   * with Francis double-shift and Wilkinson shift strategies.
   * 
   * ALGORITHM:
   * 1. Reduce matrix to upper Hessenberg form H = Q^T * A * Q
   * 2. Apply QR iteration with shifts until H converges to quasi-upper-triangular
   * 3. Extract eigenvalues from diagonal and 2x2 blocks
   * 4. Compute eigenvectors via back-substitution
   * 
   * FEATURES:
   * - Handles both real and complex eigenvalues
   * - Complex eigenvalues are stored as conjugate pairs
   * - Uses Wilkinson shift for real eigenvalues (cubic convergence)
   * - Uses Francis double-shift for complex eigenvalue pairs
   * - Deflation to reduce problem size as eigenvalues converge
   * - Exceptional shifts to prevent stagnation
   * 
   * COMPLEXITY: O(n³) - dominated by Hessenberg reduction and QR iterations
   * 
   * OUTPUT FORMAT:
   * - eigenvalues: Array of {real, imag} pairs
   * - eigenvectors: Matrix where column i corresponds to eigenvalue i
   * - For complex eigenvalues λ = a±bi:
   *   - Stored as consecutive pairs with same real part
   *   - Columns i, i+1 represent real and imaginary parts of eigenvector
   * 
   * REFERENCES:
   * - Golub & Van Loan, "Matrix Computations", 4th ed., Chapter 7
   * - Francis, "The QR Transformation—A Unitary Analogue to the LR Transformation"
   * - Wilkinson, "The Algebraic Eigenvalue Problem"
   ***************************************************************************************************/
  class EigenSolver
  {
  public:
    /**
     * @struct ComplexEigenvalue
     * @brief Represents a complex eigenvalue λ = real + imag*i
     */
    struct ComplexEigenvalue
    {
      Real real;
      Real imag;
      
      ComplexEigenvalue(Real r = 0.0, Real i = 0.0) : real(r), imag(i) {}
      
      bool isComplex(Real tol = 1e-12) const { return std::abs(imag) > tol; }
      Real magnitude() const { return std::sqrt(real * real + imag * imag); }
    };

    /**
     * @struct Result
     * @brief Complete eigensolution result
     */
    struct Result
    {
      std::vector<ComplexEigenvalue> eigenvalues;
      Matrix<Real> eigenvectors;   // Column i is eigenvector for eigenvalue i
      std::vector<bool> isComplexPair;  // True if columns i,i+1 form complex pair
      bool converged;
      int iterations;
      Real maxResidual;
      
      Result(int n) : eigenvectors(n, n), converged(false), iterations(0), maxResidual(0) {}
    };

    /**
     * Solve the complete eigenvalue problem for a general matrix.
     * 
     * @param A Input matrix (n x n)
     * @param tol Convergence tolerance (default: 1e-10)
     * @param maxIter Maximum total QR iterations (default: 1000)
     * @return Result with all eigenvalues and eigenvectors
     * 
     * POSTCONDITIONS:
     * - eigenvalues.size() == n
     * - eigenvectors has n columns
     * - For real eigenvalue λ: A*v ≈ λ*v
     * - For complex pair λ±μi with eigenvector (vr ± i*vi):
     *   - eigenvalues[k] = {λ, μ}, eigenvalues[k+1] = {λ, -μ}
     *   - eigenvectors.Column(k) = vr, Column(k+1) = vi
     */
    static Result Solve(const Matrix<Real>& A, Real tol = 1e-10, int maxIter = 1000)
    {
      int n = A.RowNum();
      Result result(n);
      
      if (n == 0)
      {
        result.converged = true;
        return result;
      }
      
      if (n == 1)
      {
        result.eigenvalues.push_back(ComplexEigenvalue(A(0, 0), 0.0));
        result.eigenvectors(0, 0) = 1.0;
        result.isComplexPair.push_back(false);
        result.converged = true;
        return result;
      }
      
      // Step 1: Reduce to Hessenberg form
      auto hess = EigenSolverHelpers::ReduceToHessenberg(A);
      Matrix<Real> H = hess.H;
      Matrix<Real> Q = hess.Q;  // Accumulated transformation
      
      // Step 2: QR iteration with deflation
      int iLo = 0;        // Start of active region
      int iHi = n - 1;    // End of active region
      int totalIter = 0;
      
      while (iHi > iLo && totalIter < maxIter)
      {
        // Check for deflation at bottom of active region
        auto defl = EigenSolverHelpers::CheckDeflation(H, iLo, iHi, tol);
        
        if (defl.canDeflate)
        {
          // Apply deflation
          EigenSolverHelpers::ApplyDeflation(H, defl.deflationIndex);
          
          // Shrink active region
          if (defl.deflationIndex == iHi)
          {
            iHi--;
          }
          else if (defl.deflationIndex == iHi - 1 && defl.blockSize == 2)
          {
            iHi -= 2;  // 2x2 block deflated
          }
          else
          {
            // Middle deflation - just continue
            iHi = defl.deflationIndex - 1;
          }
          continue;
        }
        
        // Determine if we should use single or double shift
        // Check if bottom 2x2 has complex eigenvalues
        bool useDoubleShift = false;
        if (iHi >= iLo + 1)
        {
          Real a = H(iHi - 1, iHi - 1);
          Real b = H(iHi - 1, iHi);
          Real c = H(iHi, iHi - 1);
          Real d = H(iHi, iHi);
          auto eig = EigenSolverHelpers::Eigenvalues2x2(a, b, c, d);
          useDoubleShift = eig.isComplex;
          
          // Special case: if we have only a 2x2 block with complex eigenvalues,
          // it's already in final form - we're done with this block
          if (useDoubleShift && iHi == iLo + 1)
          {
            // 2x2 block with complex eigenvalues is irreducible - accept it
            break;
          }
        }
        
        if (useDoubleShift && iHi >= iLo + 2)
        {
          // Apply Francis double-shift
          auto step = EigenSolverHelpers::FrancisDoubleShift(H, iLo, iHi, true);
          H = step.H;
          
          // Update accumulated Q
          Matrix<Real> newQ(n, n);
          for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
            {
              newQ(i, j) = 0.0;
              for (int k = 0; k < n; k++)
                newQ(i, j) += Q(i, k) * step.Q(k, j);
            }
          Q = newQ;
          totalIter += 2;  // Double shift counts as 2 iterations
        }
        else
        {
          // Apply single Wilkinson shift
          auto step = EigenSolverHelpers::SingleQRStep(H, true);
          H = step.H;
          
          // Update accumulated Q
          Matrix<Real> newQ(n, n);
          for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
            {
              newQ(i, j) = 0.0;
              for (int k = 0; k < n; k++)
                newQ(i, j) += Q(i, k) * step.Q(k, j);
            }
          Q = newQ;
          totalIter++;
        }
        
        // Apply exceptional shift if stuck (every 30 iterations without progress)
        if (totalIter % 30 == 29 && iHi >= 2)
        {
          // Exceptional shift to break cycles
          Real exceptionalShift = std::abs(H(iHi, iHi - 1)) + std::abs(H(iHi - 1, iHi - 2));
          auto step = EigenSolverHelpers::SingleQRStep(H, true, &exceptionalShift);
          H = step.H;
          
          Matrix<Real> newQ(n, n);
          for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
            {
              newQ(i, j) = 0.0;
              for (int k = 0; k < n; k++)
                newQ(i, j) += Q(i, k) * step.Q(k, j);
            }
          Q = newQ;
          totalIter++;
        }
      }
      
      // Step 3: Extract eigenvalues from quasi-triangular form
      auto eigResult = EigenSolverHelpers::ExtractEigenvalues(H, tol);
      
      // Convert to our ComplexEigenvalue format
      for (const auto& e : eigResult.eigenvalues)
        result.eigenvalues.push_back(ComplexEigenvalue(e.real, e.imag));
      
      // Step 4: Compute eigenvectors
      auto vecResult = EigenSolverHelpers::ComputeEigenvectorsFromSchur(H, Q, tol);
      result.eigenvectors = vecResult.vectors;
      result.isComplexPair = vecResult.isComplexPair;
      
      // Step 5: Compute residuals and check convergence
      result.iterations = totalIter;
      result.maxResidual = 0.0;
      result.converged = (totalIter < maxIter);
      
      // Verify eigenvector quality
      for (size_t i = 0; i < result.eigenvalues.size(); i++)
      {
        if (!result.isComplexPair[i] || (i > 0 && result.isComplexPair[i-1]))
        {
          Vector<Real> v(n);
          for (int j = 0; j < n; j++)
            v[j] = result.eigenvectors(j, static_cast<int>(i));
          
          Real residual = EigenSolverHelpers::EigenvectorResidual(
            A, v, result.eigenvalues[i].real, result.eigenvalues[i].imag);
          result.maxResidual = std::max(result.maxResidual, residual);
        }
      }
      
      return result;
    }
  };

} // namespace MML

#endif // MML_EIGENSYSTEM_SOLVERS_H

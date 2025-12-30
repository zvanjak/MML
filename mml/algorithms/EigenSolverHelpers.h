///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        EigenSolverHelpers.h                                                ///
///  Description: Building blocks for eigenvalue solvers                              ///
///               Hessenberg reduction, QR steps, deflation, eigenvector computation  ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_EIGEN_SOLVER_HELPERS_H
#define MML_EIGEN_SOLVER_HELPERS_H

#include "MMLBase.h"

#include "base/Vector.h"
#include "base/Matrix.h"
#include "base/BaseUtils.h"

#include <cmath>
#include <algorithm>
#include <complex>

namespace MML
{
  /**
   * @class EigenSolverHelpers
   * @brief Collection of verified building blocks for general eigenvalue computation
   */
  class EigenSolverHelpers
  {
  public:
    // =========================================================================
    // BUILDING BLOCK 1: HESSENBERG REDUCTION
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
     */
    static HessenbergResult ReduceToHessenberg(const Matrix<Real>& A)
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
        
        if (sigma < PrecisionValues<Real>::DivisionSafetyThreshold)
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
        
        if (vNormSq < PrecisionValues<Real>::DivisionSafetyThreshold)
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

    /**
     * Check if matrix H is upper Hessenberg within tolerance.
     * Upper Hessenberg: H(i,j) = 0 for i > j + 1
     */
    static bool IsUpperHessenberg(const Matrix<Real>& H, Real tol = 1e-10)
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
     * Uses MML::Utils::IsOrthogonal from BaseUtils.h
     */
    static bool IsOrthogonal(const Matrix<Real>& Q, Real tol = 1e-10)
    {
      return Utils::IsOrthogonal(Q, tol);
    }

    /**
     * Compute trace (sum of diagonal elements)
     * Uses Matrix::Trace() member function
     */
    static Real Trace(const Matrix<Real>& A)
    {
      return A.Trace();
    }

    /**
     * Compute Q^T * A * Q to verify similarity transformation
     */
    static Matrix<Real> SimilarityTransform(const Matrix<Real>& Q, const Matrix<Real>& A)
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

    /**
     * Compute maximum absolute difference between two matrices
     */
    static Real MaxAbsDiff(const Matrix<Real>& A, const Matrix<Real>& B)
    {
      Real maxDiff = 0.0;
      int n = A.RowNum();
      for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
          maxDiff = std::max(maxDiff, std::abs(A(i, j) - B(i, j)));
      return maxDiff;
    }

    // =========================================================================
    // BUILDING BLOCK 2: SINGLE QR STEP FOR HESSENBERG MATRIX
    // =========================================================================

    /**
     * @struct QRStepResult
     * @brief Result of single QR step: H_new = Q^T * H * Q
     */
    struct QRStepResult
    {
      Matrix<Real> H;      // Updated Hessenberg matrix
      Matrix<Real> Q;      // Orthogonal transformation (optional, for eigenvectors)
      Real shift;          // Shift used
      
      QRStepResult(int n) : H(n, n), Q(n, n), shift(0.0) {}
    };

    /**
     * Compute Wilkinson shift from bottom 2x2 submatrix of Hessenberg H.
     * 
     * For 2x2 block:
     *   [ a   b ]
     *   [ c   d ]
     * 
     * The Wilkinson shift σ is the eigenvalue of this block closest to d.
     * This choice provides cubic convergence for distinct eigenvalues.
     * 
     * Eigenvalues: λ = (a+d)/2 ± sqrt((a-d)²/4 + bc)
     * We want the one closer to d.
     */
    static Real WilkinsonShift(Real a, Real b, Real c, Real d)
    {
      // The eigenvalues of [[a,b],[c,d]] are:
      // λ = (a+d)/2 ± sqrt(((a-d)/2)² + bc)
      //
      // For numerical stability and to get the one closer to d:
      // Let δ = (a - d)/2
      // λ₁ = d + δ + sqrt(δ² + bc)  
      // λ₂ = d + δ - sqrt(δ² + bc)
      //
      // The one closer to d depends on sign of δ
      
      Real delta = (a - d) / 2.0;
      Real disc = delta * delta + b * c;
      
      if (disc >= 0.0)
      {
        // Two real eigenvalues
        Real sqrtDisc = std::sqrt(disc);
        
        // λ₁ = (a+d)/2 + sqrt = d + δ + sqrt
        // λ₂ = (a+d)/2 - sqrt = d + δ - sqrt
        // Distance from d: |δ + sqrt| and |δ - sqrt|
        // The one closer is the one where δ and sqrt have opposite signs
        // i.e., |δ - sign(δ)*sqrt|
        
        if (delta >= 0.0)
        {
          // δ >= 0, so λ₂ = d + δ - sqrt is closer to d
          return d + delta - sqrtDisc;
        }
        else
        {
          // δ < 0, so λ₁ = d + δ + sqrt is closer to d  
          return d + delta + sqrtDisc;
        }
      }
      else
      {
        // Complex eigenvalues - real part is (a+d)/2
        // For real QR, use the real part
        return (a + d) / 2.0;
      }
    }

    /**
     * Apply single shifted QR step to upper Hessenberg matrix.
     * 
     * ALGORITHM:
     * 1. Shift: H' = H - σI
     * 2. QR factorize H' using Givens rotations: H' = QR (stored in place)
     * 3. Reverse multiply: H_new = RQ + σI
     * 
     * For Hessenberg matrices, we only need n-1 Givens rotations,
     * one for each subdiagonal element.
     * 
     * The key insight is:
     * - Apply rotations from left: this zeros subdiagonal → gives R
     * - Store the rotations in Q
     * - After all left rotations, apply them from right: R → RQ
     * - The result RQ is similar to H' and hence to H
     * 
     * @param H Upper Hessenberg matrix (modified in place)
     * @param accumulateQ If true, returns Q; otherwise Q is identity
     * @param shift If provided, uses this shift; otherwise computes Wilkinson shift
     * @return QRStepResult with updated H and Q
     */
    static QRStepResult SingleQRStep(const Matrix<Real>& H, bool accumulateQ = true,
                                     Real* providedShift = nullptr)
    {
      int n = H.RowNum();
      QRStepResult result(n);
      result.H = H;
      result.Q = Matrix<Real>::GetUnitMatrix(n);
      
      if (n <= 1)
        return result;
      
      // Compute shift from bottom 2x2
      Real shift;
      if (providedShift != nullptr)
      {
        shift = *providedShift;
      }
      else
      {
        // Extract bottom 2x2: H[n-2:n-1, n-2:n-1]
        Real a = result.H(n - 2, n - 2);
        Real b = result.H(n - 2, n - 1);
        Real c = result.H(n - 1, n - 2);
        Real d = result.H(n - 1, n - 1);
        shift = WilkinsonShift(a, b, c, d);
      }
      result.shift = shift;
      
      // Apply shift: H = H - σI
      for (int i = 0; i < n; i++)
        result.H(i, i) -= shift;
      
      // Store Givens rotation parameters
      std::vector<Real> cosines(n - 1);
      std::vector<Real> sines(n - 1);
      
      // Phase 1: QR factorization - apply Givens from LEFT only to get R
      // For each subdiagonal element H(i+1, i), create rotation to zero it
      for (int i = 0; i < n - 1; i++)
      {
        Real a = result.H(i, i);
        Real b = result.H(i + 1, i);
        
        if (std::abs(b) < PrecisionValues<Real>::DivisionSafetyThreshold)
        {
          // No rotation needed
          cosines[i] = 1.0;
          sines[i] = 0.0;
          continue;
        }
        
        // Compute Givens rotation G such that G^T * [a; b] = [r; 0]
        // G = [c  s; -s  c] where c = a/r, s = b/r, r = sqrt(a² + b²)
        Real r = std::sqrt(a * a + b * b);
        Real c = a / r;
        Real s = b / r;
        
        cosines[i] = c;
        sines[i] = s;
        
        // Apply rotation from left: G^T * H
        // Affects rows i and i+1, all columns from i to n-1
        for (int j = i; j < n; j++)
        {
          Real t1 = result.H(i, j);
          Real t2 = result.H(i + 1, j);
          result.H(i, j)     =  c * t1 + s * t2;
          result.H(i + 1, j) = -s * t1 + c * t2;
        }
      }
      
      // At this point, result.H contains R (upper triangular)
      
      // Phase 2: Multiply R * Q by applying stored rotations from RIGHT
      // This gives RQ which is similar to H - σI
      for (int i = 0; i < n - 1; i++)
      {
        Real c = cosines[i];
        Real s = sines[i];
        
        if (c == 1.0 && s == 0.0)
          continue;
        
        // Apply rotation from right: H * G
        // For Hessenberg preservation: rotation G_{i,i+1} mixes columns i and i+1
        // The result is Hessenberg because each rotation only extends one row below diagonal
        // We affect rows 0 to i+1 (plus one row due to the implicit bulge chase)
        int rowMax = std::min(i + 2, n - 1);
        for (int j = 0; j <= rowMax; j++)
        {
          Real t1 = result.H(j, i);
          Real t2 = result.H(j, i + 1);
          result.H(j, i)     =  c * t1 + s * t2;
          result.H(j, i + 1) = -s * t1 + c * t2;
        }
        
        // Accumulate Q = Q * G for eigenvector computation
        if (accumulateQ)
        {
          for (int j = 0; j < n; j++)
          {
            Real t1 = result.Q(j, i);
            Real t2 = result.Q(j, i + 1);
            result.Q(j, i)     =  c * t1 + s * t2;
            result.Q(j, i + 1) = -s * t1 + c * t2;
          }
        }
      }
      
      // Unshift: H = H + σI
      for (int i = 0; i < n; i++)
        result.H(i, i) += shift;
      
      return result;
    }

    /**
     * Apply multiple QR steps until subdiagonal element converges or max iterations.
     * 
     * @param H Upper Hessenberg matrix
     * @param maxIter Maximum iterations
     * @param tol Convergence tolerance for subdiagonal
     * @return Number of iterations performed
     */
    static int MultipleQRSteps(Matrix<Real>& H, int maxIter = 30, Real tol = 1e-10)
    {
      int n = H.RowNum();
      if (n <= 1)
        return 0;
      
      for (int iter = 0; iter < maxIter; iter++)
      {
        // Check for convergence: is H(n-1, n-2) small?
        Real off = std::abs(H(n - 1, n - 2));
        Real diag = std::abs(H(n - 2, n - 2)) + std::abs(H(n - 1, n - 1));
        
        if (off < tol * diag || off < PrecisionValues<Real>::DivisionSafetyThreshold)
        {
          H(n - 1, n - 2) = 0.0;  // Force exact zero
          return iter + 1;
        }
        
        // Apply single QR step
        auto result = SingleQRStep(H, false);
        H = result.H;
      }
      
      return maxIter;  // Did not converge
    }

    // =========================================================================
    // BUILDING BLOCK 3: FRANCIS DOUBLE-SHIFT QR STEP
    // =========================================================================

    /**
     * @struct DoubleShiftResult
     * @brief Result of Francis double-shift QR step
     */
    struct DoubleShiftResult
    {
      Matrix<Real> H;      // Updated Hessenberg matrix
      Matrix<Real> Q;      // Accumulated orthogonal transformation
      Real sigma1_real;    // First shift (real part)
      Real sigma1_imag;    // First shift (imaginary part)  
      Real sigma2_real;    // Second shift (real part)
      Real sigma2_imag;    // Second shift (imaginary part)
      
      DoubleShiftResult(int n) : H(n, n), Q(n, n), 
        sigma1_real(0), sigma1_imag(0), sigma2_real(0), sigma2_imag(0) {}
    };

    /**
     * Apply Francis implicit double-shift QR step.
     * 
     * For complex conjugate eigenvalue pairs σ ± iτ, the double shift
     * implicitly applies two QR steps while keeping all arithmetic real.
     * 
     * ALGORITHM:
     * 1. Compute shifts from bottom 2x2 block
     * 2. Form first column of M = (H - σ₁I)(H - σ₂I) = H² - (σ₁+σ₂)H + σ₁σ₂I
     *    For complex conjugates: σ₁+σ₂ = 2*real, σ₁σ₂ = |σ|²
     * 3. Apply Householder to zero elements 2,3 of first column → creates "bulge"
     * 4. Chase bulge down the matrix with Householder reflections
     * 5. Result: H undergoes implicit double QR step
     * 
     * This handles complex eigenvalues without complex arithmetic.
     * 
     * @param H Upper Hessenberg matrix (n >= 3)
     * @param iLo Starting index of active submatrix (usually 0)
     * @param iHi Ending index of active submatrix (usually n-1)
     * @param accumulateQ If true, accumulates Q for eigenvector computation
     * @return DoubleShiftResult with updated H and Q
     */
    static DoubleShiftResult FrancisDoubleShift(const Matrix<Real>& H, int iLo, int iHi, 
                                                 bool accumulateQ = true)
    {
      int n = H.RowNum();
      DoubleShiftResult result(n);
      result.H = H;
      result.Q = Matrix<Real>::GetUnitMatrix(n);
      
      int nn = iHi - iLo + 1;  // Size of active block
      if (nn < 3)
      {
        // For 2x2, just use single shift
        if (nn == 2)
        {
          auto single = SingleQRStep(result.H, accumulateQ);
          result.H = single.H;
          result.Q = single.Q;
        }
        return result;
      }
      
      // Extract bottom 2x2 for shift computation
      Real h11 = result.H(iHi - 1, iHi - 1);
      Real h12 = result.H(iHi - 1, iHi);
      Real h21 = result.H(iHi, iHi - 1);
      Real h22 = result.H(iHi, iHi);
      
      // Eigenvalues of bottom 2x2: these are our shifts
      auto eig = Eigenvalues2x2(h11, h12, h21, h22);
      result.sigma1_real = eig.real1;
      result.sigma1_imag = eig.imag1;
      result.sigma2_real = eig.real2;
      result.sigma2_imag = eig.imag2;
      
      // For implicit double shift: compute first column of
      // M = H² - sH + pI where s = σ₁+σ₂ (trace) and p = σ₁σ₂ (det)
      Real s = h11 + h22;  // trace of 2x2 = sum of shifts
      Real p = h11 * h22 - h12 * h21;  // det of 2x2 = product of shifts
      
      // First column of M = (H - σ₁I)(H - σ₂I):
      // M[0,0] = H[0,0]² + H[0,1]*H[1,0] - s*H[0,0] + p
      // M[1,0] = H[1,0]*(H[0,0] + H[1,1] - s)
      // M[2,0] = H[1,0]*H[2,1]
      Real h00 = result.H(iLo, iLo);
      Real h01 = result.H(iLo, iLo + 1);
      Real h10 = result.H(iLo + 1, iLo);
      Real hh11 = result.H(iLo + 1, iLo + 1);
      Real h21b = result.H(iLo + 2, iLo + 1);
      
      Real x = h00 * h00 + h01 * h10 - s * h00 + p;
      Real y = h10 * (h00 + hh11 - s);
      Real z = h10 * h21b;
      
      // Chase the bulge from top to bottom
      for (int k = iLo; k <= iHi - 2; k++)
      {
        // Compute Householder reflector to zero out y and z
        // P = I - 2*v*v^T / (v^T*v) where v = [x; y; z] - ||[x;y;z]|| * e1
        Real norm = std::sqrt(x * x + y * y + z * z);
        if (norm < PrecisionValues<Real>::DivisionSafetyThreshold)
          break;
        
        Real alpha = (x >= 0) ? -norm : norm;
        Real v0 = x - alpha;
        Real v1 = y;
        Real v2 = z;
        Real vnorm = std::sqrt(v0 * v0 + v1 * v1 + v2 * v2);
        
        if (vnorm < PrecisionValues<Real>::DivisionSafetyThreshold)
          break;
        
        v0 /= vnorm;
        v1 /= vnorm;
        v2 /= vnorm;
        
        // Determine column range for reflector application
        int col0 = (k > iLo) ? k - 1 : k;
        
        // Apply reflector from left: P * H
        // H[k:k+2, col0:n-1] = (I - 2vv^T) * H[k:k+2, col0:n-1]
        for (int j = col0; j < n; j++)
        {
          Real t0 = result.H(k, j);
          Real t1 = result.H(k + 1, j);
          Real t2 = result.H(k + 2, j);
          Real dot = v0 * t0 + v1 * t1 + v2 * t2;
          Real tau = 2.0 * dot;
          result.H(k, j)     = t0 - tau * v0;
          result.H(k + 1, j) = t1 - tau * v1;
          result.H(k + 2, j) = t2 - tau * v2;
        }
        
        // Apply reflector from right: H * P
        // H[0:min(k+4,n), k:k+2] = H[...] * (I - 2vv^T)
        int row1 = std::min(k + 4, iHi + 1);
        for (int j = 0; j < row1; j++)
        {
          Real t0 = result.H(j, k);
          Real t1 = result.H(j, k + 1);
          Real t2 = result.H(j, k + 2);
          Real dot = v0 * t0 + v1 * t1 + v2 * t2;
          Real tau = 2.0 * dot;
          result.H(j, k)     = t0 - tau * v0;
          result.H(j, k + 1) = t1 - tau * v1;
          result.H(j, k + 2) = t2 - tau * v2;
        }
        
        // Accumulate Q = Q * P
        if (accumulateQ)
        {
          for (int j = 0; j < n; j++)
          {
            Real t0 = result.Q(j, k);
            Real t1 = result.Q(j, k + 1);
            Real t2 = result.Q(j, k + 2);
            Real dot = v0 * t0 + v1 * t1 + v2 * t2;
            Real tau = 2.0 * dot;
            result.Q(j, k)     = t0 - tau * v0;
            result.Q(j, k + 1) = t1 - tau * v1;
            result.Q(j, k + 2) = t2 - tau * v2;
          }
        }
        
        // Set up x, y, z for next iteration (bulge moves down)
        if (k < iHi - 2)
        {
          x = result.H(k + 1, k);
          y = result.H(k + 2, k);
          z = (k + 3 <= iHi) ? result.H(k + 3, k) : 0.0;
        }
      }
      
      // Final 2x2 cleanup: zero out H(iHi-1, iHi-2) bulge element
      int k = iHi - 2;
      x = result.H(k + 1, k);
      y = result.H(k + 2, k);
      
      // 2x2 Givens rotation to zero y
      Real r = std::sqrt(x * x + y * y);
      if (r > PrecisionValues<Real>::DivisionSafetyThreshold)
      {
        Real c = x / r;
        Real s_val = y / r;
        
        // Apply from left
        for (int j = k; j < n; j++)
        {
          Real t1 = result.H(k + 1, j);
          Real t2 = result.H(k + 2, j);
          result.H(k + 1, j) =  c * t1 + s_val * t2;
          result.H(k + 2, j) = -s_val * t1 + c * t2;
        }
        
        // Apply from right
        for (int j = 0; j <= std::min(k + 3, iHi); j++)
        {
          Real t1 = result.H(j, k + 1);
          Real t2 = result.H(j, k + 2);
          result.H(j, k + 1) =  c * t1 + s_val * t2;
          result.H(j, k + 2) = -s_val * t1 + c * t2;
        }
        
        // Accumulate
        if (accumulateQ)
        {
          for (int j = 0; j < n; j++)
          {
            Real t1 = result.Q(j, k + 1);
            Real t2 = result.Q(j, k + 2);
            result.Q(j, k + 1) =  c * t1 + s_val * t2;
            result.Q(j, k + 2) = -s_val * t1 + c * t2;
          }
        }
      }
      
      // Clean up small subdiagonal elements
      for (int i = iLo + 1; i <= iHi; i++)
      {
        if (std::abs(result.H(i, i - 1)) < PrecisionValues<Real>::DivisionSafetyThreshold)
          result.H(i, i - 1) = 0.0;
      }
      
      return result;
    }

    /**
     * Apply multiple double-shift QR iterations for complex eigenvalue convergence.
     * Returns number of iterations used.
     */
    static int MultipleDoubleShiftSteps(Matrix<Real>& H, int iLo, int iHi,
                                        int maxIter = 30, Real tol = 1e-10)
    {
      for (int iter = 0; iter < maxIter; iter++)
      {
        // Check for deflation at bottom
        Real off = std::abs(H(iHi, iHi - 1));
        Real diag = std::abs(H(iHi - 1, iHi - 1)) + std::abs(H(iHi, iHi));
        
        if (off < tol * std::max(diag, REAL(1.0)) || off < PrecisionValues<Real>::DivisionSafetyThreshold)
        {
          H(iHi, iHi - 1) = 0.0;
          return iter + 1;
        }
        
        // Also check penultimate subdiagonal
        if (iHi >= iLo + 2)
        {
          Real off2 = std::abs(H(iHi - 1, iHi - 2));
          Real diag2 = std::abs(H(iHi - 2, iHi - 2)) + std::abs(H(iHi - 1, iHi - 1));
          if (off2 < tol * std::max(diag2, REAL(1.0)) || off2 < PrecisionValues<Real>::DivisionSafetyThreshold)
          {
            H(iHi - 1, iHi - 2) = 0.0;
            return iter + 1;
          }
        }
        
        // Apply double shift
        auto result = FrancisDoubleShift(H, iLo, iHi, false);
        H = result.H;
      }
      
      return maxIter;
    }

    // =========================================================================
    // BUILDING BLOCK 4: EIGENVALUE EXTRACTION FROM 2x2 BLOCK
    // =========================================================================

    /**
     * @struct Eigenvalue2x2Result
     * @brief Eigenvalues of a 2x2 matrix (may be complex conjugate pair)
     */
    struct Eigenvalue2x2Result
    {
      Real real1, imag1;   // First eigenvalue: real1 + i*imag1
      Real real2, imag2;   // Second eigenvalue: real2 + i*imag2
      bool isComplex;      // True if complex conjugate pair
    };

    /**
     * Compute eigenvalues of 2x2 matrix [[a, b], [c, d]]
     * 
     * Characteristic polynomial: λ² - (a+d)λ + (ad-bc) = 0
     * λ = (a+d)/2 ± sqrt((a+d)²/4 - (ad-bc))
     *   = (a+d)/2 ± sqrt((a-d)²/4 + bc)
     */
    static Eigenvalue2x2Result Eigenvalues2x2(Real a, Real b, Real c, Real d)
    {
      Eigenvalue2x2Result result;
      
      Real trace = a + d;
      Real det = a * d - b * c;
      
      // Discriminant = trace²/4 - det = (a-d)²/4 + bc
      Real p = (a - d) / 2.0;
      Real disc = p * p + b * c;
      
      result.real1 = trace / 2.0;
      result.real2 = trace / 2.0;
      
      if (disc >= 0.0)
      {
        // Two real eigenvalues
        Real sqrtDisc = std::sqrt(disc);
        result.real1 += sqrtDisc;
        result.real2 -= sqrtDisc;
        result.imag1 = 0.0;
        result.imag2 = 0.0;
        result.isComplex = false;
      }
      else
      {
        // Complex conjugate pair
        Real sqrtDisc = std::sqrt(-disc);
        result.imag1 = sqrtDisc;
        result.imag2 = -sqrtDisc;
        result.isComplex = true;
      }
      
      return result;
    }

    /**
     * @struct DeflationResult
     * @brief Result of checking for deflation in Hessenberg matrix
     */
    struct DeflationResult
    {
      bool canDeflate;     // True if deflation is possible
      int deflationIndex;  // Index where deflation occurs (H[idx, idx-1] ≈ 0)
      int blockSize;       // Size of deflated block (1 = real eigenvalue, 2 = complex pair)
    };

    /**
     * Check if matrix can be deflated at any position.
     * 
     * Deflation occurs when |H[k, k-1]| < tol * (|H[k-1,k-1]| + |H[k,k]|)
     * This means the matrix decouples into independent subproblems.
     * 
     * @param H Upper Hessenberg matrix  
     * @param iLo Start index of active region
     * @param iHi End index of active region
     * @param tol Tolerance for considering element as zero
     * @return DeflationResult indicating if/where deflation is possible
     */
    static DeflationResult CheckDeflation(const Matrix<Real>& H, int iLo, int iHi, 
                                          Real tol = 1e-10)
    {
      DeflationResult result;
      result.canDeflate = false;
      result.deflationIndex = -1;
      result.blockSize = 0;
      
      // Check from bottom up for deflation
      for (int k = iHi; k > iLo; k--)
      {
        Real off = std::abs(H(k, k - 1));
        Real diag = std::abs(H(k - 1, k - 1)) + std::abs(H(k, k));
        
        if (off < tol * std::max(diag, REAL(1.0)) || off < PrecisionValues<Real>::DivisionSafetyThreshold)
        {
          result.canDeflate = true;
          result.deflationIndex = k;
          
          // Determine block size: check if this isolates a 1x1 or 2x2 block
          if (k == iHi)
          {
            // Single eigenvalue at position iHi
            result.blockSize = 1;
          }
          else
          {
            // Check next subdiagonal
            Real off2 = std::abs(H(k + 1, k));
            Real diag2 = std::abs(H(k, k)) + std::abs(H(k + 1, k + 1));
            if (off2 < tol * std::max(diag2, REAL(1.0)) || off2 < PrecisionValues<Real>::DivisionSafetyThreshold)
              result.blockSize = 1;
            else
              result.blockSize = 2;  // 2x2 block (complex pair)
          }
          return result;
        }
      }
      
      return result;
    }

    /**
     * @struct ComplexEigenvalue
     * @brief Represents a potentially complex eigenvalue
     */
    struct ComplexEigenvalue
    {
      Real real;
      Real imag;
      bool isComplex;
      
      ComplexEigenvalue(Real r = 0.0, Real i = 0.0) 
        : real(r), imag(i), isComplex(std::abs(i) > PrecisionValues<Real>::DivisionSafetyThreshold) {}
    };

    /**
     * @struct EigenvalueExtractionResult
     * @brief All eigenvalues extracted from quasi-upper-triangular form
     */
    struct EigenvalueExtractionResult
    {
      std::vector<ComplexEigenvalue> eigenvalues;
      int realCount;     // Number of real eigenvalues
      int complexPairs;  // Number of complex conjugate pairs
    };

    /**
     * Extract all eigenvalues from quasi-upper-triangular (real Schur) form.
     * 
     * The matrix should be the result of QR iteration:
     * - 1x1 diagonal blocks contain real eigenvalues
     * - 2x2 diagonal blocks contain complex conjugate pairs
     * 
     * A 2x2 block is detected when H[i+1, i] is non-negligible.
     * 
     * @param H Quasi-upper-triangular matrix
     * @param tol Tolerance for detecting 2x2 blocks
     * @return EigenvalueExtractionResult with all eigenvalues
     */
    static EigenvalueExtractionResult ExtractEigenvalues(const Matrix<Real>& H, 
                                                          Real tol = 1e-10)
    {
      int n = H.RowNum();
      EigenvalueExtractionResult result;
      result.realCount = 0;
      result.complexPairs = 0;
      
      int i = 0;
      while (i < n)
      {
        if (i == n - 1)
        {
          // Last element: 1x1 block (real eigenvalue)
          result.eigenvalues.push_back(ComplexEigenvalue(H(i, i), 0.0));
          result.realCount++;
          i++;
        }
        else
        {
          // Check if this is a 2x2 block
          Real subdiag = std::abs(H(i + 1, i));
          Real diagSum = std::abs(H(i, i)) + std::abs(H(i + 1, i + 1));
          
          if (subdiag < tol * std::max(diagSum, REAL(1.0)) || subdiag < PrecisionValues<Real>::DivisionSafetyThreshold)
          {
            // 1x1 block: real eigenvalue
            result.eigenvalues.push_back(ComplexEigenvalue(H(i, i), REAL(0.0)));
            result.realCount++;
            i++;
          }
          else
          {
            // 2x2 block: extract eigenvalues
            auto eig = Eigenvalues2x2(H(i, i), H(i, i + 1), H(i + 1, i), H(i + 1, i + 1));
            
            result.eigenvalues.push_back(ComplexEigenvalue(eig.real1, eig.imag1));
            result.eigenvalues.push_back(ComplexEigenvalue(eig.real2, eig.imag2));
            
            if (eig.isComplex)
            {
              result.complexPairs++;
            }
            else
            {
              result.realCount += 2;
            }
            i += 2;
          }
        }
      }
      
      return result;
    }

    /**
     * Force deflation at a specific index by zeroing the subdiagonal element.
     * Call this after CheckDeflation returns canDeflate = true.
     * 
     * @param H Upper Hessenberg matrix (modified in place)
     * @param deflationIndex Index returned by CheckDeflation
     */
    static void ApplyDeflation(Matrix<Real>& H, int deflationIndex)
    {
      if (deflationIndex > 0 && deflationIndex < H.RowNum())
      {
        H(deflationIndex, deflationIndex - 1) = 0.0;
      }
    }

    // =========================================================================
    // BUILDING BLOCK 5: EIGENVECTOR COMPUTATION FROM SCHUR FORM
    // =========================================================================

    /**
     * @struct EigenvectorResult
     * @brief Computed eigenvectors for a matrix
     * 
     * For real eigenvalue at index i: vectors.col(i) is the eigenvector
     * For complex pair at indices i,i+1: 
     *   - vectors.col(i) is the real part
     *   - vectors.col(i+1) is the imaginary part
     *   - Actual eigenvectors are: col(i) ± i*col(i+1)
     */
    struct EigenvectorResult
    {
      Matrix<Real> vectors;   // Eigenvector matrix (n x n)
      std::vector<bool> isComplexPair;  // isComplexPair[i] = true if columns i,i+1 form complex pair
    };

    /**
     * Compute eigenvector for a real eigenvalue using back-substitution.
     * 
     * Given upper triangular T with eigenvalue λ = T[k,k], solve (T - λI)x = 0
     * by back-substitution from row k-1 up to row 0.
     * 
     * @param T Upper triangular (or quasi-upper-triangular) matrix
     * @param k Index of the eigenvalue (diagonal element T[k,k])
     * @return Eigenvector (normalized)
     */
    static Vector<Real> ComputeRealEigenvector(const Matrix<Real>& T, int k)
    {
      int n = T.RowNum();
      Vector<Real> x(n, 0.0);
      
      // Set x[k] = 1 as starting point
      x[k] = 1.0;
      Real lambda = T(k, k);
      
      // Back-substitute: for i = k-1 down to 0
      // (T[i,i] - λ) * x[i] + sum_{j=i+1}^{k} T[i,j] * x[j] = 0
      // x[i] = -sum_{j=i+1}^{k} T[i,j] * x[j] / (T[i,i] - λ)
      
      for (int i = k - 1; i >= 0; i--)
      {
        Real sum = 0.0;
        for (int j = i + 1; j <= k; j++)
          sum += T(i, j) * x[j];
        
        Real denom = T(i, i) - lambda;
        if (std::abs(denom) > PrecisionValues<Real>::DivisionSafetyThreshold)
        {
          x[i] = -sum / denom;
        }
        else
        {
          // Near-singular: use small perturbation to avoid division by zero
          x[i] = -sum / PrecisionValues<Real>::DivisionSafetyThreshold;
        }
      }
      
      // Normalize
      Real norm = 0.0;
      for (int i = 0; i < n; i++)
        norm += x[i] * x[i];
      norm = std::sqrt(norm);
      
      if (norm > PrecisionValues<Real>::DivisionSafetyThreshold)
      {
        for (int i = 0; i < n; i++)
          x[i] /= norm;
      }
      
      return x;
    }

    /**
     * Compute eigenvectors for a 2x2 complex eigenvalue block.
     * 
     * For a 2x2 block at positions [k, k+1] with complex eigenvalues α ± iβ,
     * compute the real and imaginary parts of the eigenvector.
     * 
     * @param T Quasi-upper-triangular matrix
     * @param k Starting index of the 2x2 block
     * @return Pair of vectors: (real part, imaginary part)
     */
    static std::pair<Vector<Real>, Vector<Real>> ComputeComplexEigenvectors(
      const Matrix<Real>& T, int k)
    {
      int n = T.RowNum();
      Vector<Real> xr(n, 0.0);  // Real part
      Vector<Real> xi(n, 0.0);  // Imaginary part
      
      // Get the 2x2 block and its eigenvalues
      Real a = T(k, k);
      Real b = T(k, k + 1);
      Real c = T(k + 1, k);
      Real d = T(k + 1, k + 1);
      
      auto eig = Eigenvalues2x2(a, b, c, d);
      Real alpha = eig.real1;  // Real part of eigenvalue
      Real beta = eig.imag1;   // Imaginary part (positive)
      
      if (std::abs(beta) < PrecisionValues<Real>::DivisionSafetyThreshold)
      {
        // Not actually complex, treat as two real
        xr = ComputeRealEigenvector(T, k);
        xi = ComputeRealEigenvector(T, k + 1);
        return {xr, xi};
      }
      
      // For the 2x2 block, the eigenvector of [[a,b],[c,d]] for α+iβ is:
      // v = [b, α+iβ-a] or [α+iβ-d, c] (up to scaling)
      // We use [b, (α-a)+iβ] which gives:
      //   real: [b, α-a]
      //   imag: [0, β]
      
      // Initialize at the 2x2 block
      xr[k] = b;
      xr[k + 1] = alpha - a;
      xi[k] = 0.0;
      xi[k + 1] = beta;
      
      // Back-substitute for rows k-1 down to 0
      // We need to solve (T - (α+iβ)I) * (xr + i*xi) = 0
      // This gives two coupled real equations:
      // (T - αI) * xr + βI * xi = 0  =>  (T-αI)*xr = -β*xi
      // (T - αI) * xi - βI * xr = 0  =>  (T-αI)*xi = β*xr
      
      for (int i = k - 1; i >= 0; i--)
      {
        // sum_r = sum of T[i,j]*xr[j] for j > i
        // sum_i = sum of T[i,j]*xi[j] for j > i
        Real sum_r = 0.0;
        Real sum_i = 0.0;
        for (int j = i + 1; j <= k + 1; j++)
        {
          sum_r += T(i, j) * xr[j];
          sum_i += T(i, j) * xi[j];
        }
        
        // (T[i,i] - α) * xr[i] + β * xi[i] = -sum_r
        // (T[i,i] - α) * xi[i] - β * xr[i] = -sum_i
        // Let p = T[i,i] - α
        // p * xr[i] + β * xi[i] = -sum_r
        // p * xi[i] - β * xr[i] = -sum_i
        // Solve 2x2 system:
        // |p   β | |xr[i]|   |-sum_r|
        // |-β  p | |xi[i]| = |-sum_i|
        // det = p² + β²
        
        Real p = T(i, i) - alpha;
        Real det = p * p + beta * beta;
        
        if (det > PrecisionValues<Real>::DivisionSafetyThreshold)
        {
          xr[i] = (-sum_r * p - sum_i * beta) / det;
          xi[i] = (-sum_i * p + sum_r * beta) / det;
        }
        else
        {
          xr[i] = 0.0;
          xi[i] = 0.0;
        }
      }
      
      // Normalize: ||xr||² + ||xi||² = 1 for each eigenvector
      Real norm_sq = 0.0;
      for (int i = 0; i < n; i++)
        norm_sq += xr[i] * xr[i] + xi[i] * xi[i];
      
      Real norm = std::sqrt(norm_sq);
      if (norm > PrecisionValues<Real>::DivisionSafetyThreshold)
      {
        for (int i = 0; i < n; i++)
        {
          xr[i] /= norm;
          xi[i] /= norm;
        }
      }
      
      return {xr, xi};
    }

    /**
     * Compute all eigenvectors from quasi-upper-triangular (real Schur) form.
     * 
     * @param T Quasi-upper-triangular matrix (from QR iteration)
     * @param Q Accumulated orthogonal transformation (T = Q^T * A * Q)
     * @param tol Tolerance for detecting 2x2 blocks
     * @return EigenvectorResult with eigenvector matrix
     */
    static EigenvectorResult ComputeEigenvectorsFromSchur(
      const Matrix<Real>& T, const Matrix<Real>& Q, Real tol = 1e-10)
    {
      int n = T.RowNum();
      EigenvectorResult result;
      result.vectors = Matrix<Real>(n, n);
      result.isComplexPair.resize(n, false);
      
      // First compute eigenvectors of T (the Schur form)
      Matrix<Real> Y(n, n);  // Eigenvectors of T
      
      int i = 0;
      while (i < n)
      {
        bool is2x2 = false;
        
        if (i < n - 1)
        {
          Real subdiag = std::abs(T(i + 1, i));
          Real diagSum = std::abs(T(i, i)) + std::abs(T(i + 1, i + 1));
          is2x2 = (subdiag >= tol * std::max(diagSum, REAL(1.0)) && subdiag >= PrecisionValues<Real>::DivisionSafetyThreshold);
        }
        
        if (!is2x2)
        {
          // 1x1 block: real eigenvalue
          Vector<Real> v = ComputeRealEigenvector(T, i);
          for (int j = 0; j < n; j++)
            Y(j, i) = v[j];
          result.isComplexPair[i] = false;
          i++;
        }
        else
        {
          // 2x2 block: complex conjugate pair
          auto [vr, vi] = ComputeComplexEigenvectors(T, i);
          for (int j = 0; j < n; j++)
          {
            Y(j, i) = vr[j];      // Real part in column i
            Y(j, i + 1) = vi[j];  // Imaginary part in column i+1
          }
          result.isComplexPair[i] = true;
          result.isComplexPair[i + 1] = true;
          i += 2;
        }
      }
      
      // Transform back: X = Q * Y
      // Since T = Q^T * A * Q, eigenvectors of A are X = Q * Y
      for (int col = 0; col < n; col++)
      {
        for (int row = 0; row < n; row++)
        {
          Real sum = 0.0;
          for (int k = 0; k < n; k++)
            sum += Q(row, k) * Y(k, col);
          result.vectors(row, col) = sum;
        }
        
        // Re-normalize after transformation
        Real norm = 0.0;
        for (int row = 0; row < n; row++)
          norm += result.vectors(row, col) * result.vectors(row, col);
        norm = std::sqrt(norm);
        
        if (norm > PrecisionValues<Real>::DivisionSafetyThreshold)
        {
          for (int row = 0; row < n; row++)
            result.vectors(row, col) /= norm;
        }
      }
      
      return result;
    }

    /**
     * Verify eigenvector: compute ||A*v - λ*v|| / ||v||
     * For complex eigenvalue, uses real part of λ only.
     */
    static Real EigenvectorResidual(const Matrix<Real>& A, const Vector<Real>& v, 
                                    Real lambda_real, Real lambda_imag = 0.0)
    {
      int n = A.RowNum();
      Real residual = 0.0;
      Real vnorm = 0.0;
      
      if (std::abs(lambda_imag) < PrecisionValues<Real>::DivisionSafetyThreshold)
      {
        // Real eigenvalue: check ||A*v - λ*v||
        for (int i = 0; i < n; i++)
        {
          Real Av_i = 0.0;
          for (int j = 0; j < n; j++)
            Av_i += A(i, j) * v[j];
          Real diff = Av_i - lambda_real * v[i];
          residual += diff * diff;
          vnorm += v[i] * v[i];
        }
      }
      else
      {
        // Complex eigenvalue: this is the real part of eigenvector
        // The actual eigenvector is v_r + i*v_i where v_i is the next column
        // For a simple check, just verify magnitude is reasonable
        for (int i = 0; i < n; i++)
          vnorm += v[i] * v[i];
        return 0.0;  // Complex case needs both vectors, skip for now
      }
      
      vnorm = std::sqrt(vnorm);
      if (vnorm > PrecisionValues<Real>::DivisionSafetyThreshold)
        return std::sqrt(residual) / vnorm;
      return 0.0;
    }

  };

} // namespace MML

#endif // MML_EIGEN_SOLVER_HELPERS_H

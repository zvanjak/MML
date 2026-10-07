#ifndef MML_EIGEN_DETAIL_HESSENBERG_QR_ITERATION_H
#define MML_EIGEN_DETAIL_HESSENBERG_QR_ITERATION_H

#include <mml/MMLBase.h>
#include <mml/base/Vector/Vector.h>
#include <mml/base/Matrix/Matrix.h>
#include <mml/base/BaseUtils.h>
#include <mml/algorithms/Eigen/detail/RealSchurAnalysis.h>
#include <algorithm>
#include <cmath>
#include <complex>

namespace MML::detail {
	// =========================================================================
	// BUILDING BLOCK 2: SINGLE QR STEP FOR HESSENBERG MATRIX
	// =========================================================================

	/// @struct QRStepResult
	/// @brief Result of single QR step: H_new = Q^T * H * Q
	struct QRStepResult {
		Matrix<Real> H; // Updated Hessenberg matrix
		Matrix<Real> Q; // Orthogonal transformation (optional, for eigenvectors)
		Real shift;		// Shift used

		QRStepResult(int n)
			: H(n, n)
			, Q(n, n)
			, shift(0.0) {}
	};

	/// Compute Wilkinson shift from bottom 2x2 submatrix of Hessenberg H.
	///
	/// For 2x2 block:
	/// [ a   b ]
	/// [ c   d ]
	///
	/// The Wilkinson shift σ is the eigenvalue of this block closest to d.
	/// This choice provides cubic convergence for distinct eigenvalues.
	///
	/// Eigenvalues: λ = (a+d)/2 ± sqrt((a-d)²/4 + bc)
	/// We want the one closer to d.
	inline Real WilkinsonShift(Real a, Real b, Real c, Real d) {
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

		if (disc >= 0.0) {
			// Two real eigenvalues
			Real sqrtDisc = std::sqrt(disc);

			// λ₁ = (a+d)/2 + sqrt = d + δ + sqrt
			// λ₂ = (a+d)/2 - sqrt = d + δ - sqrt
			// Distance from d: |δ + sqrt| and |δ - sqrt|
			// The one closer is the one where δ and sqrt have opposite signs
			// i.e., |δ - sign(δ)*sqrt|

			if (delta >= 0.0) {
				// δ >= 0, so λ₂ = d + δ - sqrt is closer to d
				return d + delta - sqrtDisc;
			} else {
				// δ < 0, so λ₁ = d + δ + sqrt is closer to d
				return d + delta + sqrtDisc;
			}
		} else {
			// Complex eigenvalues - real part is (a+d)/2
			// For real QR, use the real part
			return (a + d) / 2.0;
		}
	}

	/// Apply single shifted QR step to upper Hessenberg matrix.
	///
	/// ALGORITHM:
	/// 1. Shift: H' = H - σI
	/// 2. QR factorize H' using Givens rotations: H' = QR (stored in place)
	/// 3. Reverse multiply: H_new = RQ + σI
	///
	/// For Hessenberg matrices, we only need n-1 Givens rotations,
	/// one for each subdiagonal element.
	///
	/// The key insight is:
	/// - Apply rotations from left: this zeros subdiagonal → gives R
	/// - Store the rotations in Q
	/// - After all left rotations, apply them from right: R → RQ
	/// - The result RQ is similar to H' and hence to H
	///
	/// @param H Upper Hessenberg matrix (modified in place)
	/// @param accumulateQ If true, returns Q; otherwise Q is identity
	/// @param shift If provided, uses this shift; otherwise computes Wilkinson shift
	/// @return QRStepResult with updated H and Q
	inline QRStepResult SingleQRStep(const Matrix<Real>& H, bool accumulateQ = true, Real* providedShift = nullptr) {
		int n = H.rows();
		QRStepResult result(n);
		result.H = H;
		result.Q = Matrix<Real>::Identity(n);

		if (n <= 1)
			return result;

		// Compute shift from bottom 2x2
		Real shift;
		if (providedShift != nullptr) {
			shift = *providedShift;
		} else {
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
		for (int i = 0; i < n - 1; i++) {
			Real a = result.H(i, i);
			Real b = result.H(i + 1, i);

			if (std::abs(b) < PrecisionValues<Real>::DivisionSafetyThreshold) {
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
			for (int j = i; j < n; j++) {
				Real t1 = result.H(i, j);
				Real t2 = result.H(i + 1, j);
				result.H(i, j) = c * t1 + s * t2;
				result.H(i + 1, j) = -s * t1 + c * t2;
			}
		}

		// At this point, result.H contains R (upper triangular)

		// Phase 2: Multiply R * Q by applying stored rotations from RIGHT
		// This gives RQ which is similar to H - σI
		for (int i = 0; i < n - 1; i++) {
			Real c = cosines[i];
			Real s = sines[i];

			if (c == 1.0 && s == 0.0)
				continue;

			// Apply rotation from right: H * G
			// For Hessenberg preservation: rotation G_{i,i+1} mixes columns i and i+1
			// The result is Hessenberg because each rotation only extends one row below diagonal
			// We affect rows 0 to i+1 (plus one row due to the implicit bulge chase)
			int rowMax = std::min(i + 2, n - 1);
			for (int j = 0; j <= rowMax; j++) {
				Real t1 = result.H(j, i);
				Real t2 = result.H(j, i + 1);
				result.H(j, i) = c * t1 + s * t2;
				result.H(j, i + 1) = -s * t1 + c * t2;
			}

			// Accumulate Q = Q * G for eigenvector computation
			if (accumulateQ) {
				for (int j = 0; j < n; j++) {
					Real t1 = result.Q(j, i);
					Real t2 = result.Q(j, i + 1);
					result.Q(j, i) = c * t1 + s * t2;
					result.Q(j, i + 1) = -s * t1 + c * t2;
				}
			}
		}

		// Unshift: H = H + σI
		for (int i = 0; i < n; i++)
			result.H(i, i) += shift;

		return result;
	}

	/// Apply multiple QR steps until subdiagonal element converges or max iterations.
	///
	/// @param H Upper Hessenberg matrix
	/// @param maxIter Maximum iterations
	/// @param tol Convergence tolerance for subdiagonal
	/// @return Number of iterations performed
	inline int MultipleQRSteps(Matrix<Real>& H, int maxIter = 30, Real tol = 1e-10) {
		int n = H.rows();
		if (n <= 1)
			return 0;

		for (int iter = 0; iter < maxIter; iter++) {
			// Check for convergence: is H(n-1, n-2) small?
			Real off = std::abs(H(n - 1, n - 2));
			Real diag = std::abs(H(n - 2, n - 2)) + std::abs(H(n - 1, n - 1));

			if (off < tol * diag || off < PrecisionValues<Real>::DivisionSafetyThreshold) {
				H(n - 1, n - 2) = 0.0; // Force exact zero
				return iter + 1;
			}

			// Apply single QR step
			auto result = SingleQRStep(H, false);
			H = result.H;
		}

		return maxIter; // Did not converge
	}

	// =========================================================================
	// BUILDING BLOCK 3: FRANCIS DOUBLE-SHIFT QR STEP
	// =========================================================================

	/// @struct DoubleShiftResult
	/// @brief Result of Francis double-shift QR step
	struct DoubleShiftResult {
		Matrix<Real> H;	  // Updated Hessenberg matrix
		Matrix<Real> Q;	  // Accumulated orthogonal transformation
		Real sigma1_real; // First shift (real part)
		Real sigma1_imag; // First shift (imaginary part)
		Real sigma2_real; // Second shift (real part)
		Real sigma2_imag; // Second shift (imaginary part)

		DoubleShiftResult(int n)
			: H(n, n)
			, Q(n, n)
			, sigma1_real(0)
			, sigma1_imag(0)
			, sigma2_real(0)
			, sigma2_imag(0) {}
	};

	/// Apply Francis implicit double-shift QR step.
	///
	/// For complex conjugate eigenvalue pairs σ ± iτ, the double shift
	/// implicitly applies two QR steps while keeping all arithmetic real.
	///
	/// ALGORITHM:
	/// 1. Compute shifts from bottom 2x2 block
	/// 2. Form first column of M = (H - σ₁I)(H - σ₂I) = H² - (σ₁+σ₂)H + σ₁σ₂I
	/// For complex conjugates: σ₁+σ₂ = 2*real, σ₁σ₂ = |σ|²
	/// 3. Apply Householder to zero elements 2,3 of first column → creates "bulge"
	/// 4. Chase bulge down the matrix with Householder reflections
	/// 5. Result: H undergoes implicit double QR step
	///
	/// This handles complex eigenvalues without complex arithmetic.
	///
	/// @param H Upper Hessenberg matrix (n >= 3)
	/// @param iLo Starting index of active submatrix (usually 0)
	/// @param iHi Ending index of active submatrix (usually n-1)
	/// @param accumulateQ If true, accumulates Q for eigenvector computation
	/// @return DoubleShiftResult with updated H and Q
	inline DoubleShiftResult FrancisDoubleShift(const Matrix<Real>& H, int iLo, int iHi, bool accumulateQ = true) {
		int n = H.rows();
		DoubleShiftResult result(n);
		result.H = H;
		result.Q = Matrix<Real>::Identity(n);

		int nn = iHi - iLo + 1; // Size of active block
		if (nn < 3) {
			// For 2x2, just use single shift
			if (nn == 2) {
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
		Real s = h11 + h22;				// trace of 2x2 = sum of shifts
		Real p = h11 * h22 - h12 * h21; // det of 2x2 = product of shifts

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
		for (int k = iLo; k <= iHi - 2; k++) {
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
			for (int j = col0; j < n; j++) {
				Real t0 = result.H(k, j);
				Real t1 = result.H(k + 1, j);
				Real t2 = result.H(k + 2, j);
				Real dot = v0 * t0 + v1 * t1 + v2 * t2;
				Real tau = 2.0 * dot;
				result.H(k, j) = t0 - tau * v0;
				result.H(k + 1, j) = t1 - tau * v1;
				result.H(k + 2, j) = t2 - tau * v2;
			}

			// Apply reflector from right: H * P
			// H[0:min(k+4,n), k:k+2] = H[...] * (I - 2vv^T)
			int row1 = std::min(k + 4, iHi + 1);
			for (int j = 0; j < row1; j++) {
				Real t0 = result.H(j, k);
				Real t1 = result.H(j, k + 1);
				Real t2 = result.H(j, k + 2);
				Real dot = v0 * t0 + v1 * t1 + v2 * t2;
				Real tau = 2.0 * dot;
				result.H(j, k) = t0 - tau * v0;
				result.H(j, k + 1) = t1 - tau * v1;
				result.H(j, k + 2) = t2 - tau * v2;
			}

			// Accumulate Q = Q * P
			if (accumulateQ) {
				for (int j = 0; j < n; j++) {
					Real t0 = result.Q(j, k);
					Real t1 = result.Q(j, k + 1);
					Real t2 = result.Q(j, k + 2);
					Real dot = v0 * t0 + v1 * t1 + v2 * t2;
					Real tau = 2.0 * dot;
					result.Q(j, k) = t0 - tau * v0;
					result.Q(j, k + 1) = t1 - tau * v1;
					result.Q(j, k + 2) = t2 - tau * v2;
				}
			}

			// Set up x, y, z for next iteration (bulge moves down)
			if (k < iHi - 2) {
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
		if (r > PrecisionValues<Real>::DivisionSafetyThreshold) {
			Real c = x / r;
			Real s_val = y / r;

			// Apply from left
			for (int j = k; j < n; j++) {
				Real t1 = result.H(k + 1, j);
				Real t2 = result.H(k + 2, j);
				result.H(k + 1, j) = c * t1 + s_val * t2;
				result.H(k + 2, j) = -s_val * t1 + c * t2;
			}

			// Apply from right
			for (int j = 0; j <= std::min(k + 3, iHi); j++) {
				Real t1 = result.H(j, k + 1);
				Real t2 = result.H(j, k + 2);
				result.H(j, k + 1) = c * t1 + s_val * t2;
				result.H(j, k + 2) = -s_val * t1 + c * t2;
			}

			// Accumulate
			if (accumulateQ) {
				for (int j = 0; j < n; j++) {
					Real t1 = result.Q(j, k + 1);
					Real t2 = result.Q(j, k + 2);
					result.Q(j, k + 1) = c * t1 + s_val * t2;
					result.Q(j, k + 2) = -s_val * t1 + c * t2;
				}
			}
		}

		// Clean up small subdiagonal elements
		for (int i = iLo + 1; i <= iHi; i++) {
			if (std::abs(result.H(i, i - 1)) < PrecisionValues<Real>::DivisionSafetyThreshold)
				result.H(i, i - 1) = 0.0;
		}

		return result;
	}

	/// Apply multiple double-shift QR iterations for complex eigenvalue convergence.
	/// Returns number of iterations used.
	inline int MultipleDoubleShiftSteps(Matrix<Real>& H, int iLo, int iHi, int maxIter = 30, Real tol = 1e-10) {
		for (int iter = 0; iter < maxIter; iter++) {
			// Check for deflation at bottom
			Real off = std::abs(H(iHi, iHi - 1));
			Real diag = std::abs(H(iHi - 1, iHi - 1)) + std::abs(H(iHi, iHi));

			if (off < tol * std::max(diag, REAL(1.0)) || off < PrecisionValues<Real>::DivisionSafetyThreshold) {
				H(iHi, iHi - 1) = 0.0;
				return iter + 1;
			}

			// Also check penultimate subdiagonal
			if (iHi >= iLo + 2) {
				Real off2 = std::abs(H(iHi - 1, iHi - 2));
				Real diag2 = std::abs(H(iHi - 2, iHi - 2)) + std::abs(H(iHi - 1, iHi - 1));
				if (off2 < tol * std::max(diag2, REAL(1.0)) || off2 < PrecisionValues<Real>::DivisionSafetyThreshold) {
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


} // namespace MML::detail

#endif // MML_EIGEN_DETAIL_HESSENBERG_QR_ITERATION_H

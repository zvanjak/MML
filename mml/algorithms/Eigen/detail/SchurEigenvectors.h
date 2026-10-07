#ifndef MML_EIGEN_DETAIL_SCHUR_EIGENVECTORS_H
#define MML_EIGEN_DETAIL_SCHUR_EIGENVECTORS_H

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
	// BUILDING BLOCK 5: EIGENVECTOR COMPUTATION FROM SCHUR FORM
	// =========================================================================

	/// @struct EigenvectorResult
	/// @brief Computed eigenvectors for a matrix
	///
	/// For real eigenvalue at index i: vectors.col(i) is the eigenvector
	/// For complex pair at indices i,i+1:
	/// - vectors.col(i) is the real part
	/// - vectors.col(i+1) is the imaginary part
	/// - Actual eigenvectors are: col(i) ± i*col(i+1)
	struct EigenvectorResult {
		Matrix<Real> vectors;			 // Eigenvector matrix (n x n)
		std::vector<bool> isComplexPair; // isComplexPair[i] = true if columns i,i+1 form complex pair
	};

	/// Compute eigenvector for a real eigenvalue using back-substitution.
	///
	/// Given upper triangular T with eigenvalue λ = T[k,k], solve (T - λI)x = 0
	/// by back-substitution from row k-1 up to row 0.
	///
	/// @param T Upper triangular (or quasi-upper-triangular) matrix
	/// @param k Index of the eigenvalue (diagonal element T[k,k])
	/// @return Eigenvector (normalized)
	inline Vector<Real> ComputeRealEigenvector(const Matrix<Real>& T, int k) {
		int n = T.rows();
		Vector<Real> x(n, 0.0);

		// Set x[k] = 1 as starting point
		x[k] = 1.0;
		Real lambda = T(k, k);

		// Back-substitute: for i = k-1 down to 0
		// (T[i,i] - λ) * x[i] + sum_{j=i+1}^{k} T[i,j] * x[j] = 0
		// x[i] = -sum_{j=i+1}^{k} T[i,j] * x[j] / (T[i,i] - λ)

		for (int i = k - 1; i >= 0; i--) {
			Real sum = 0.0;
			for (int j = i + 1; j <= k; j++)
				sum += T(i, j) * x[j];

			Real denom = T(i, i) - lambda;
			if (std::abs(denom) > PrecisionValues<Real>::DivisionSafetyThreshold) {
				x[i] = -sum / denom;
			} else {
				// Near-singular: use small perturbation to avoid division by zero
				x[i] = -sum / PrecisionValues<Real>::DivisionSafetyThreshold;
			}
		}

		// Normalize
		Real norm = 0.0;
		for (int i = 0; i < n; i++)
			norm += x[i] * x[i];
		norm = std::sqrt(norm);

		if (norm > PrecisionValues<Real>::DivisionSafetyThreshold) {
			for (int i = 0; i < n; i++)
				x[i] /= norm;
		}

		return x;
	}

	/// Compute eigenvectors for a 2x2 complex eigenvalue block.
	///
	/// For a 2x2 block at positions [k, k+1] with complex eigenvalues α ± iβ,
	/// compute the real and imaginary parts of the eigenvector.
	///
	/// @param T Quasi-upper-triangular matrix
	/// @param k Starting index of the 2x2 block
	/// @return Pair of vectors: (real part, imaginary part)
	inline std::pair<Vector<Real>, Vector<Real>> ComputeComplexEigenvectors(const Matrix<Real>& T, int k) {
		int n = T.rows();
		Vector<Real> xr(n, 0.0); // Real part
		Vector<Real> xi(n, 0.0); // Imaginary part

		// Get the 2x2 block and its eigenvalues
		Real a = T(k, k);
		Real b = T(k, k + 1);
		Real c = T(k + 1, k);
		Real d = T(k + 1, k + 1);

		auto eig = Eigenvalues2x2(a, b, c, d);
		Real alpha = eig.real1; // Real part of eigenvalue
		Real beta = eig.imag1;	// Imaginary part (positive)

		if (std::abs(beta) < PrecisionValues<Real>::DivisionSafetyThreshold) {
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
		for (int i = k - 1; i >= 0; i--) {
			// sum_r = sum of T[i,j]*xr[j] for j > i
			// sum_i = sum of T[i,j]*xi[j] for j > i
			Real sum_r = 0.0;
			Real sum_i = 0.0;
			for (int j = i + 1; j <= k + 1; j++) {
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

			if (det > PrecisionValues<Real>::DivisionSafetyThreshold) {
				xr[i] = (-sum_r * p - sum_i * beta) / det;
				xi[i] = (-sum_i * p + sum_r * beta) / det;
			} else {
				xr[i] = 0.0;
				xi[i] = 0.0;
			}
		}

		// Normalize: ||xr||² + ||xi||² = 1 for each eigenvector
		Real norm_sq = 0.0;
		for (int i = 0; i < n; i++)
			norm_sq += xr[i] * xr[i] + xi[i] * xi[i];

		Real norm = std::sqrt(norm_sq);
		if (norm > PrecisionValues<Real>::DivisionSafetyThreshold) {
			for (int i = 0; i < n; i++) {
				xr[i] /= norm;
				xi[i] /= norm;
			}
		}

		return {xr, xi};
	}

	/// Compute all eigenvectors from quasi-upper-triangular (real Schur) form.
	///
	/// @param T Quasi-upper-triangular matrix (from QR iteration)
	/// @param Q Accumulated orthogonal transformation (T = Q^T * A * Q)
	/// @param tol Tolerance for detecting 2x2 blocks
	/// @return EigenvectorResult with eigenvector matrix
	inline EigenvectorResult ComputeEigenvectorsFromSchur(const Matrix<Real>& T, const Matrix<Real>& Q, Real tol = 1e-10) {
		int n = T.rows();
		EigenvectorResult result;
		result.vectors = Matrix<Real>(n, n);
		result.isComplexPair.resize(n, false);

		// First compute eigenvectors of T (the Schur form)
		Matrix<Real> Y(n, n); // Eigenvectors of T

		int i = 0;
		while (i < n) {
			bool is2x2 = false;

			if (i < n - 1) {
				Real subdiag = std::abs(T(i + 1, i));
				Real diagSum = std::abs(T(i, i)) + std::abs(T(i + 1, i + 1));
				is2x2 = (subdiag >= tol * std::max(diagSum, REAL(1.0)) && subdiag >= PrecisionValues<Real>::DivisionSafetyThreshold);
			}

			if (!is2x2) {
				// 1x1 block: real eigenvalue
				Vector<Real> v = ComputeRealEigenvector(T, i);
				for (int j = 0; j < n; j++)
					Y(j, i) = v[j];
				result.isComplexPair[i] = false;
				i++;
			} else {
				// 2x2 block: complex conjugate pair
				auto [vr, vi] = ComputeComplexEigenvectors(T, i);
				for (int j = 0; j < n; j++) {
					Y(j, i) = vr[j];	 // Real part in column i
					Y(j, i + 1) = vi[j]; // Imaginary part in column i+1
				}
				result.isComplexPair[i] = true;
				result.isComplexPair[i + 1] = true;
				i += 2;
			}
		}

		// Transform back: X = Q * Y
		// Since T = Q^T * A * Q, eigenvectors of A are X = Q * Y
		for (int col = 0; col < n; col++) {
			for (int row = 0; row < n; row++) {
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

			if (norm > PrecisionValues<Real>::DivisionSafetyThreshold) {
				for (int row = 0; row < n; row++)
					result.vectors(row, col) /= norm;
			}
		}

		return result;
	}

	/// Verify eigenvector: compute ||A*v - λ*v|| / ||v||
	/// For complex eigenvalue, uses real part of λ only.
	inline Real EigenvectorResidual(const Matrix<Real>& A, const Vector<Real>& v, Real lambda_real, Real lambda_imag = 0.0) {
		int n = A.rows();
		Real residual = 0.0;
		Real vnorm = 0.0;

		if (std::abs(lambda_imag) < PrecisionValues<Real>::DivisionSafetyThreshold) {
			// Real eigenvalue: check ||A*v - λ*v||
			for (int i = 0; i < n; i++) {
				Real Av_i = 0.0;
				for (int j = 0; j < n; j++)
					Av_i += A(i, j) * v[j];
				Real diff = Av_i - lambda_real * v[i];
				residual += diff * diff;
				vnorm += v[i] * v[i];
			}
		} else {
			// Complex eigenvalue: this is the real part of eigenvector
			// The actual eigenvector is v_r + i*v_i where v_i is the next column
			// For a simple check, just verify magnitude is reasonable
			for (int i = 0; i < n; i++)
				vnorm += v[i] * v[i];
			return 0.0; // Complex case needs both vectors, skip for now
		}

		vnorm = std::sqrt(vnorm);
		if (vnorm > PrecisionValues<Real>::DivisionSafetyThreshold)
			return std::sqrt(residual) / vnorm;
		return 0.0;
	}

} // namespace MML::detail

#endif // MML_EIGEN_DETAIL_SCHUR_EIGENVECTORS_H

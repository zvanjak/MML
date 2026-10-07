#ifndef MML_EIGEN_HESSENBERG_REDUCTION_H
#define MML_EIGEN_HESSENBERG_REDUCTION_H

#include <mml/MMLBase.h>
#include <mml/base/Vector/Vector.h>
#include <mml/base/Matrix/Matrix.h>
#include <mml/base/BaseUtils.h>
#include <algorithm>
#include <cmath>

namespace MML {
	// =========================================================================
	// HESSENBERG REDUCTION - Canonical implementation
	// MatrixAlg::ReduceToHessenberg() delegates to this
	// =========================================================================

	/// @struct HessenbergResult
	/// @brief Result of Hessenberg reduction: H = Q^T * A * Q
	struct HessenbergResult {
		Matrix<Real> H; // Upper Hessenberg matrix
		Matrix<Real> Q; // Orthogonal transformation matrix

		HessenbergResult(int n)
			: H(n, n)
			, Q(n, n) {}
	};

	/// Reduce matrix A to upper Hessenberg form using Householder reflections.
	///
	/// This is the CANONICAL implementation used by all eigen solvers.
	/// MatrixAlg::ReduceToHessenberg() delegates to this function.
	inline HessenbergResult ReduceToHessenberg(const Matrix<Real>& A) {
		int n = A.rows();
		HessenbergResult result(n);

		if (n <= 2) {
			result.H = A;
			result.Q = Matrix<Real>::Identity(n);
			return result;
		}

		// Copy A to H (we'll transform H in place)
		result.H = A;
		// Initialize Q as identity
		result.Q = Matrix<Real>::Identity(n);

		// Temporary storage for Householder vector
		Vector<Real> v(n);

		// For each column k, zero out elements below the subdiagonal
		for (int k = 0; k < n - 2; k++) {
			// Compute the norm of the column below the diagonal
			Real sigma = 0.0;
			for (int i = k + 1; i < n; i++)
				sigma += result.H(i, k) * result.H(i, k);
			sigma = std::sqrt(sigma);

			if (sigma < PrecisionValues<Real>::DivisionSafetyThreshold)
				continue; // Column already zero below subdiagonal

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
			for (int j = k; j < n; j++) {
				Real dot = 0.0;
				for (int i = k + 1; i < n; i++)
					dot += v[i] * result.H(i, j);

				for (int i = k + 1; i < n; i++)
					result.H(i, j) -= beta * v[i] * dot;
			}

			// Apply Householder from right: H = H * (I - beta*v*v^T)
			// H(i,j) -= beta * H(i,m) * v(m) * v(j)
			for (int i = 0; i < n; i++) {
				Real dot = 0.0;
				for (int j = k + 1; j < n; j++)
					dot += result.H(i, j) * v[j];

				for (int j = k + 1; j < n; j++)
					result.H(i, j) -= beta * dot * v[j];
			}

			// Accumulate Q: Q = Q * (I - beta*v*v^T)
			// Q(i,j) -= beta * Q(i,m) * v(m) * v(j)
			for (int i = 0; i < n; i++) {
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


} // namespace MML

#endif // MML_EIGEN_HESSENBERG_REDUCTION_H

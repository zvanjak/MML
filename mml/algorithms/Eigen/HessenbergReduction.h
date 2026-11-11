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
		if (A.rows() != A.cols())
			throw MatrixDimensionError("ReduceToHessenberg: matrix must be square", A.rows(), A.cols(), -1, -1);
		int n = A.rows();
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < n; ++j)
				if (!std::isfinite(A(i, j)))
					throw MatrixNumericalError("ReduceToHessenberg: non-finite input matrix");
		HessenbergResult result(n);
		if (n == 0) return result;

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
			Real columnScale = 0;
			for (int i = k + 1; i < n; i++)
				columnScale = std::max(columnScale, std::abs(result.H(i, k)));

			if (columnScale == Real{0}) // tol-ok: exact structural zero; nonzero columns must be reduced at every scale.
				continue;

			// Dimensionless Householder construction (Golub & Van Loan, section 5.1).
			Real normSquared = 0;
			for (int i = k + 1; i < n; i++) {
				v[i] = result.H(i, k) / columnScale;
				normSquared += v[i] * v[i];
			}
			Real sigma = std::sqrt(normSquared);
			if (v[k + 1] > Real{0})
				sigma = -sigma;
			const Real subdiagonal = sigma * columnScale;
			if (!std::isfinite(subdiagonal))
				throw MatrixNumericalError("ReduceToHessenberg: transformed subdiagonal is not representable");

			v[k + 1] -= sigma;
			Real vNormSq = v[k + 1] * v[k + 1];
			for (int i = k + 2; i < n; i++)
				vNormSq += v[i] * v[i];
			const Real vNorm = std::sqrt(vNormSq);
			for (int i = k + 1; i < n; i++)
				v[i] /= vNorm;

			// Scale each projection and fuse the update so an intermediate need not fit in Real.
			for (int j = k + 1; j < n; j++) {
				Real scale = 0;
				for (int i = k + 1; i < n; i++)
					scale = std::max(scale, std::abs(result.H(i, j)));
				if (scale == Real{0}) // tol-ok: exact zero segment has no Householder projection.
					continue;
				Real dot = 0;
				for (int i = k + 1; i < n; i++)
					dot += v[i] * (result.H(i, j) / scale);

				for (int i = k + 1; i < n; i++) {
					result.H(i, j) = std::fma(-Real{2} * v[i] * dot, scale, result.H(i, j));
					if (!std::isfinite(result.H(i, j)))
						throw MatrixNumericalError("ReduceToHessenberg: left transformation is not representable");
				}
			}

			for (int i = 0; i < n; i++) {
				Real scale = 0;
				for (int j = k + 1; j < n; j++)
					scale = std::max(scale, std::abs(result.H(i, j)));
				if (scale == Real{0}) // tol-ok: exact zero segment has no Householder projection.
					continue;
				Real dot = 0;
				for (int j = k + 1; j < n; j++)
					dot += (result.H(i, j) / scale) * v[j];

				for (int j = k + 1; j < n; j++) {
					result.H(i, j) = std::fma(-Real{2} * dot * v[j], scale, result.H(i, j));
					if (!std::isfinite(result.H(i, j)))
						throw MatrixNumericalError("ReduceToHessenberg: right transformation is not representable");
				}
			}

			// Accumulate the same unit-vector reflection in Q.
			for (int i = 0; i < n; i++) {
				Real dot = 0.0;
				for (int j = k + 1; j < n; j++)
					dot += result.Q(i, j) * v[j];

				for (int j = k + 1; j < n; j++)
					result.Q(i, j) -= Real{2} * dot * v[j];
			}

			// After applying transformation, H(k+1, k) should be sigma
			// and H(k+2:n-1, k) should be 0
			// Clean up the zeros explicitly
			result.H(k + 1, k) = subdiagonal;
			for (int i = k + 2; i < n; i++)
				result.H(i, k) = 0.0;
		}

		return result;
	}


} // namespace MML

#endif // MML_EIGEN_HESSENBERG_REDUCTION_H

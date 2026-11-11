#ifndef MML_EIGEN_DETAIL_SELF_ADJOINT_SCALING_H
#define MML_EIGEN_DETAIL_SELF_ADJOINT_SCALING_H

#include <mml/MMLBase.h>
#include <mml/base/Vector/Vector.h>
#include <mml/base/Matrix/Matrix.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <iostream>

namespace MML::detail {

	inline void ValidateSelfAdjointControls(Real tolerance, int maxIterations) {
		if (!std::isfinite(tolerance) || tolerance < 0)
			throw DomainError("Self-adjoint eigensolver tolerance must be finite and nonnegative");
		if (maxIterations < 0)
			throw DomainError("Self-adjoint eigensolver iteration budget must be nonnegative");
	}

	template<class Result>
	inline void LogSelfAdjointResult(const Result& result) {
		std::clog << result.algorithm_name << ": " << to_string(result.status)
			<< ", iterations=" << result.iterations << ", residual=" << result.residual << '\n';
	}

	inline void ValidateSymmetricInput(const Matrix<Real>& matrix) {
		if (matrix.rows() != matrix.cols())
			throw MatrixDimensionError("Symmetric eigensolver requires a square matrix", matrix.rows(), matrix.cols(), -1, -1);
		for (int row = 0; row < matrix.rows(); ++row)
			for (int col = row + 1; col < matrix.cols(); ++col)
				if (std::isfinite(matrix(row, col)) && std::isfinite(matrix(col, row)) &&
					std::abs(matrix(row, col) - matrix(col, row)) > PrecisionValues<Real>::EigenSolverZeroThreshold)
					throw MatrixDimensionError("Symmetric eigensolver requires symmetric input; use SolveSymmetricPart for explicit averaging",
						matrix.rows(), matrix.cols(), -1, -1);
	}

	/// Overflow-safe average of two matrix entries used when symmetrizing general storage:
	/// 0.5*(a + b) overflows for |a|,|b| near the floating-point maximum, 0.5*a + 0.5*b does not.
	inline Real SymmetrizeEntries(Real a, Real b) {
		return MML::Real(0.5) * a + MML::Real(0.5) * b;
	}
	inline Complex HermitizeEntries(Complex a, Complex bConjugate) {
		return MML::Real(0.5) * a + MML::Real(0.5) * bConjugate;
	}

	/// Largest absolute entry of a real or complex square matrix (infinity if any entry is nonfinite).
	template<typename Scalar>
	inline Real MaxAbsEntry(const Matrix<Scalar>& A) {
		Real maxAbs = 0;
		for (int i = 0; i < A.rows(); i++)
			for (int j = 0; j < A.cols(); j++) {
				const Real magnitude = std::abs(A[i][j]);
				if (!std::isfinite(magnitude))
					return std::numeric_limits<Real>::infinity();
				maxAbs = std::max(maxAbs, magnitude);
			}
		return maxAbs;
	}

	/// Power-of-two scale that brings a matrix with maximal entry `maxAbs` into a range where
	/// squares and sums of entries cannot overflow or underflow. Returns 1 for ordinary data so
	/// the working arithmetic (and therefore all existing results) stays bit-identical; only
	/// extreme magnitudes (beyond roughly ±1e150 in double) are rescaled. Dividing by a power of
	/// two is exact, so the eigenvalues of A/scale are exactly the eigenvalues of A divided by scale.
	inline Real SelfAdjointScale(Real maxAbs) {
		if (!(maxAbs > 0) || !std::isfinite(maxAbs))
			return 1;
		const int safeExponent = std::numeric_limits<Real>::max_exponent / 2;
		int exponent = 0;
		std::frexp(maxAbs, &exponent); // maxAbs = m * 2^exponent, m in [0.5, 1)
		if (exponent > safeExponent || exponent < -safeExponent)
			return std::ldexp(Real(1), exponent - 1); // entries of A/scale lie in [-2, 2]
		return 1;
	}

	/// Off-diagonal Frobenius norm of the strict upper triangle without overflow/underflow:
	/// the sum of squares is accumulated relative to the largest off-diagonal magnitude.
	template<typename Scalar>
	inline Real OffDiagonalNormSafe(const Matrix<Scalar>& A) {
		const int n = A.rows();
		Real maxAbs = 0;
		for (int i = 0; i < n - 1; i++)
			for (int j = i + 1; j < n; j++)
				maxAbs = std::max(maxAbs, static_cast<Real>(std::abs(A[i][j])));
		if (!(maxAbs > 0) || !std::isfinite(maxAbs))
			return maxAbs; // 0, or propagate inf/NaN
		Real sum = 0;
		for (int i = 0; i < n - 1; i++)
			for (int j = i + 1; j < n; j++) {
				const Real ratio = static_cast<Real>(std::abs(A[i][j])) / maxAbs;
				sum += ratio * ratio;
			}
		return maxAbs * std::sqrt(sum);
	}

	/// Euclidean norm of a real vector without overflow/underflow.
	inline Real VectorNormSafe(const Vector<Real>& v) {
		Real maxAbs = 0;
		for (int i = 0; i < v.size(); i++)
			maxAbs = std::max(maxAbs, std::abs(v[i]));
		if (!(maxAbs > 0) || !std::isfinite(maxAbs))
			return maxAbs;
		Real sum = 0;
		for (int i = 0; i < v.size(); i++) {
			const Real ratio = v[i] / maxAbs;
			sum += ratio * ratio;
		}
		return maxAbs * std::sqrt(sum);
	}

	/// True when every eigenvalue and eigenvector entry is finite.
	template<typename Scalar>
	inline bool EigenResultIsFinite(const Vector<Real>& eigenvalues, const Matrix<Scalar>& eigenvectors) {
		for (int i = 0; i < eigenvalues.size(); i++)
			if (!std::isfinite(eigenvalues[i]))
				return false;
		for (int i = 0; i < eigenvectors.rows(); i++)
			for (int j = 0; j < eigenvectors.cols(); j++)
				if (!std::isfinite(std::abs(eigenvectors[i][j])))
					return false;
		return true;
	}

} // namespace MML::detail

#endif // MML_EIGEN_DETAIL_SELF_ADJOINT_SCALING_H

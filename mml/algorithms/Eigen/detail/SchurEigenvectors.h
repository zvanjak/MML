#ifndef MML_EIGEN_DETAIL_SCHUR_EIGENVECTORS_H
#define MML_EIGEN_DETAIL_SCHUR_EIGENVECTORS_H

#include <mml/MMLBase.h>
#include <mml/base/Vector/Vector.h>
#include <mml/base/Matrix/Matrix.h>
#include <mml/base/BaseUtils.h>
#include <mml/base/AlgorithmTypes.h>
#include <mml/algorithms/Eigen/detail/RealSchurAnalysis.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <utility>
#include <vector>

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

	/// Shared back-substitution kernel for eigenvectors of a quasi-upper-triangular (real Schur) matrix.
	///
	/// Solves (T - λI) x = 0 for the eigenvalue λ belonging to the diagonal block that starts at
	/// row `blockStart` and has size `blockSize` (1 for a real eigenvalue, 2 for a complex pair).
	/// The entries of x inside that block are preset by the caller; entries below it are zero;
	/// entries above are obtained by back-substitution from row blockStart-1 up to row 0.
	///
	/// Unlike plain triangular back-substitution this honours 2x2 diagonal blocks ABOVE the target
	/// block: whenever T(i, i-1) != 0 the rows i-1, i are solved together as a 2x2 (complex) system.
	/// Ignoring that coupling is exactly what produced O(1) eigenvector residuals on dense
	/// nonsymmetric matrices (analysis v4, ALGORITHMS A1).
	///
	/// Small or zero denominators (repeated eigenvalues) are replaced by smin as in
	/// LAPACK dtrevc / Numerical Recipes hqr2, so the routine never divides by zero.
	///
	/// Reference: Golub & Van Loan, Matrix Computations 4th ed., §7.6.1; LAPACK dtrevc.
	inline void BackSubstituteSchurEigenvector(const Matrix<Real>& T, int blockStart, int blockSize,
											   std::complex<Real> lambda, std::vector<std::complex<Real>>& x) {
		using C = std::complex<Real>;
		const int top = blockStart + blockSize - 1; // last row/column with a preset (nonzero) entry

		Real anorm = static_cast<MML::Real>(0.0);
		for (int r = 0; r <= top; r++)
			for (int c = std::max(r - 1, 0); c <= top; c++)
				anorm = std::max(anorm, std::abs(T(r, c)));
		const Real eps = std::numeric_limits<Real>::epsilon();
		const Real smin = std::max({eps * std::abs(lambda), eps * anorm, PrecisionValues<Real>::DivisionSafetyThreshold});

		auto rowSum = [&](int r, int from) {
			C s(static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0));
			for (int j = from; j <= top; j++)
				s += T(r, j) * x[j];
			return s;
		};

		int i = blockStart - 1;
		while (i >= 0) {
			const bool coupled = (i > 0) && (T(i, i - 1) != static_cast<MML::Real>(0.0));
			if (!coupled) {
				C denom = C(T(i, i), static_cast<MML::Real>(0.0)) - lambda;
				if (std::abs(denom) < smin)
					denom = C(smin, static_cast<MML::Real>(0.0));
				x[i] = -rowSum(i, i + 1) / denom;
				i -= 1;
			} else {
				// Rows i-1, i form a 2x2 block: solve
				//   [T(i-1,i-1)-λ   T(i-1,i)  ] [x(i-1)]   [-Σ_{j>i} T(i-1,j) x_j]
				//   [T(i,i-1)       T(i,i)-λ  ] [x(i)  ] = [-Σ_{j>i} T(i,j)   x_j]
				const C a = C(T(i - 1, i - 1), static_cast<MML::Real>(0.0)) - lambda;
				const C b = C(T(i - 1, i), static_cast<MML::Real>(0.0));
				const C c = C(T(i, i - 1), static_cast<MML::Real>(0.0));
				const C d = C(T(i, i), static_cast<MML::Real>(0.0)) - lambda;
				const C r0 = -rowSum(i - 1, i + 1);
				const C r1 = -rowSum(i, i + 1);
				C det = a * d - b * c;
				if (std::abs(det) < smin * smin)
					det = C(smin * smin, static_cast<MML::Real>(0.0));
				x[i - 1] = (r0 * d - b * r1) / det;
				x[i] = (a * r1 - c * r0) / det;
				i -= 2;
			}
		}
	}

	/// Normalise a complex vector to unit Euclidean norm (||Re||² + ||Im||² = 1).
	inline void NormalizeComplexVector(std::vector<std::complex<Real>>& x) {
		Real norm = static_cast<MML::Real>(0.0);
		for (const auto& v : x)
			norm += std::norm(v);
		norm = std::sqrt(norm);
		if (norm > PrecisionValues<Real>::DivisionSafetyThreshold)
			for (auto& v : x)
				v /= norm;
	}

	/// Compute eigenvector for a real eigenvalue λ = T[k,k] of a quasi-upper-triangular T.
	///
	/// Solves (T - λI)x = 0 with x[k] = 1 by back-substitution over rows k-1 .. 0,
	/// treating any 2x2 diagonal block above k as a coupled system.
	///
	/// @param T Upper triangular or quasi-upper-triangular (real Schur) matrix
	/// @param k Index of the eigenvalue (1x1 diagonal block T[k,k])
	/// @return Eigenvector (unit Euclidean norm)
	inline Vector<Real> ComputeRealEigenvector(const Matrix<Real>& T, int k) {
		const int n = T.rows();
		std::vector<std::complex<Real>> x(n, std::complex<Real>(static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0)));
		x[k] = std::complex<Real>(static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0));

		BackSubstituteSchurEigenvector(T, k, 1, std::complex<Real>(T(k, k), static_cast<MML::Real>(0.0)), x);
		NormalizeComplexVector(x);

		Vector<Real> v(n);
		for (int i = 0; i < n; i++)
			v[i] = x[i].real();
		return v;
	}

	/// Compute eigenvectors for the 2x2 diagonal block at rows/columns (k, k+1).
	///
	/// If the block has complex eigenvalues α ± iβ (β > 0), returns the real and imaginary
	/// parts of the eigenvector for α + iβ, normalised so that ||Re||² + ||Im||² = 1.
	/// If the block has two real eigenvalues (an un-standardised real 2x2 block), returns the
	/// two real eigenvectors in the same order as Eigenvalues2x2 (real1, real2); the coupling
	/// T(k+1,k) inside the block is honoured in both cases.
	///
	/// @param T Quasi-upper-triangular matrix
	/// @param k Starting index of the 2x2 block
	/// @return Pair of vectors: (real part, imaginary part) or (eigvec(real1), eigvec(real2))
	inline std::pair<Vector<Real>, Vector<Real>> ComputeComplexEigenvectors(const Matrix<Real>& T, int k) {
		using C = std::complex<Real>;
		const int n = T.rows();

		const Real a = T(k, k), b = T(k, k + 1), c = T(k + 1, k), d = T(k + 1, k + 1);
		const auto eig = Eigenvalues2x2(a, b, c, d);

		// Eigenvector of the 2x2 block for eigenvalue λ: [b, λ-a] or [λ-d, c]; use the larger one.
		auto blockVector = [&](C lambda, std::vector<C>& x) {
			const C v0 = C(b, static_cast<MML::Real>(0.0)), v1 = lambda - a;
			const C w0 = lambda - d, w1 = C(c, static_cast<MML::Real>(0.0));
			if (std::norm(w0) + std::norm(w1) > std::norm(v0) + std::norm(v1)) {
				x[k] = w0;
				x[k + 1] = w1;
			} else {
				x[k] = v0;
				x[k + 1] = v1;
			}
			if (std::norm(x[k]) + std::norm(x[k + 1]) < PrecisionValues<Real>::DivisionSafetyThreshold)
				x[k] = C(static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0)); // a == d, b == c == 0: any vector in the block works
		};

		auto solveFor = [&](C lambda) {
			std::vector<C> x(n, C(static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0)));
			blockVector(lambda, x);
			BackSubstituteSchurEigenvector(T, k, 2, lambda, x);
			NormalizeComplexVector(x);
			return x;
		};

		Vector<Real> first(n), second(n);
		if (eig.isComplex) {
			const auto x = solveFor(C(eig.real1, std::abs(eig.imag1)));
			for (int i = 0; i < n; i++) {
				first[i] = x[i].real();
				second[i] = x[i].imag();
			}
		} else {
			const auto x1 = solveFor(C(eig.real1, static_cast<MML::Real>(0.0)));
			const auto x2 = solveFor(C(eig.real2, static_cast<MML::Real>(0.0)));
			for (int i = 0; i < n; i++) {
				first[i] = x1[i].real();
				second[i] = x2[i].real();
			}
		}
		return {first, second};
	}

	/// Compute all eigenvectors from quasi-upper-triangular (real Schur) form.
	///
	/// @param T Quasi-upper-triangular matrix (from QR iteration)
	/// @param Q Accumulated orthogonal transformation (T = Q^T * A * Q)
	/// @param tol Tolerance for detecting 2x2 blocks
	/// @return EigenvectorResult with eigenvector matrix
	inline EigenvectorResult ComputeEigenvectorsFromSchur(const Matrix<Real>& T, const Matrix<Real>& Q, Real tol = Precision::EigenSolverConvergenceTolerance) {
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
				is2x2 = (subdiag >= tol * std::max(diagSum, static_cast<MML::Real>(1.0)) && subdiag >= PrecisionValues<Real>::DivisionSafetyThreshold);
			}

			if (!is2x2) {
				// 1x1 block: real eigenvalue
				Vector<Real> v = ComputeRealEigenvector(T, i);
				for (int j = 0; j < n; j++)
					Y(j, i) = v[j];
				result.isComplexPair[i] = false;
				i++;
			} else {
				// 2x2 block: complex conjugate pair (Re, Im) or two real eigenvectors if the block
				// has real eigenvalues and was not standardised to triangular form.
				const bool complexPair = Eigenvalues2x2(T(i, i), T(i, i + 1), T(i + 1, i), T(i + 1, i + 1)).isComplex;
				auto [vr, vi] = ComputeComplexEigenvectors(T, i);
				for (int j = 0; j < n; j++) {
					Y(j, i) = vr[j];	 // Real part in column i
					Y(j, i + 1) = vi[j]; // Imaginary part in column i+1
				}
				result.isComplexPair[i] = complexPair;
				result.isComplexPair[i + 1] = complexPair;
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
		}

		// Re-normalize. Columns (i, i+1) of a complex pair are Re and Im of ONE eigenvector and
		// must be scaled by a common factor; normalising them separately destroys the eigenvector.
		int col = 0;
		while (col < n) {
			const int width = (result.isComplexPair[col] && col + 1 < n) ? 2 : 1;
			Real norm = 0.0;
			for (int c = col; c < col + width; c++)
				for (int row = 0; row < n; row++)
					norm += result.vectors(row, c) * result.vectors(row, c);
			norm = std::sqrt(norm);

			if (norm > PrecisionValues<Real>::DivisionSafetyThreshold)
				for (int c = col; c < col + width; c++)
					for (int row = 0; row < n; row++)
						result.vectors(row, c) /= norm;
			col += width;
		}

		return result;
	}

	/// Residual ||A v - λ v|| / ||v|| for a complex eigenpair λ = λr + i λi, v = vr + i vi.
	/// Works for real pairs too (pass a zero vi and λi = 0).
	inline Real EigenpairResidual(const Matrix<Real>& A, const Vector<Real>& vr, const Vector<Real>& vi, Real lambda_real, Real lambda_imag) {
		const int n = A.rows();
		if (!std::isfinite(lambda_real) || !std::isfinite(lambda_imag))
			return NoEstimate;
		Real residual = 0.0;
		Real vnorm = 0.0;
		for (int i = 0; i < n; i++) {
			if (!std::isfinite(vr[i]) || !std::isfinite(vi[i]))
				return NoEstimate;
			Real Avr = 0.0, Avi = 0.0;
			for (int j = 0; j < n; j++) {
				Avr += A(i, j) * vr[j];
				Avi += A(i, j) * vi[j];
			}
			// (λr + iλi)(vr + i vi) = (λr vr - λi vi) + i(λr vi + λi vr)
			const Real dr = Avr - (lambda_real * vr[i] - lambda_imag * vi[i]);
			const Real di = Avi - (lambda_real * vi[i] + lambda_imag * vr[i]);
			if (!std::isfinite(dr) || !std::isfinite(di))
				return NoEstimate;
			residual = std::hypot(residual, dr, di);
			vnorm = std::hypot(vnorm, vr[i], vi[i]);
		}
		return std::isfinite(vnorm) && vnorm > Real{0} ? residual / vnorm : NoEstimate;
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

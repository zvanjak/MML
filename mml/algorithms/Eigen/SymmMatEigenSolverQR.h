#ifndef MML_SYMM_MAT_EIGEN_SOLVER_QR_H
#define MML_SYMM_MAT_EIGEN_SOLVER_QR_H

#include <mml/MMLBase.h>
#include <mml/base/Vector/Vector.h>
#include <mml/base/Matrix/Matrix.h>
#include <mml/base/Matrix/MatrixSym.h>
#include <mml/base/AlgorithmTypes.h>
#include <mml/algorithms/Eigen/EigenSolverConfig.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <utility>

namespace MML {

	/// @class SymmMatEigenSolverQR
	/// @brief QR algorithm for symmetric matrix eigenvalue decomposition
	///
	/// Algorithm: Tridiagonal reduction + Implicit QR iteration with Wilkinson shift
	/// 1. Reduce symmetric matrix to tridiagonal form using Householder reflections
	/// 2. Apply implicit QR algorithm with Wilkinson shift to tridiagonal matrix
	/// 3. Accumulate transformations to obtain eigenvectors
	///
	/// Complexity: O(n³) for reduction + O(n²) per QR iteration (typically O(n) iterations)
	/// Total: O(n³) - more efficient than Jacobi for large matrices
	///
	/// Pros: - Faster convergence than Jacobi (especially for larger matrices)
	/// - Cubic convergence rate with Wilkinson shift
	/// - Industry-standard algorithm
	/// Cons: - More complex implementation
	/// - Requires careful handling of deflation
	///
	/// REFERENCES:
	/// - Golub & Van Loan, "Matrix Computations", 4th ed., Sections 8.3-8.4
	/// - Numerical Recipes, 3rd ed., Section 11.3
	/// - Wilkinson, "The Algebraic Eigenvalue Problem"
	class SymmMatEigenSolverQR {
	public:
		/// Result structure for QR eigenvalue decomposition.
		///
		/// Contains eigenvalues, eigenvectors, and diagnostic information.
		/// Always check `converged` before using the results.
		struct Result {
			// === Primary Output ===
			Vector<Real> eigenvalues;  // Sorted eigenvalues (ascending)
			Matrix<Real> eigenvectors; // Corresponding eigenvectors (columns)

			// === Convergence Status ===
			bool converged = false;    // True if algorithm converged
			int iterations = 0;        // Number of QR iterations performed
			Real residual = 0.0;       // Final off-diagonal norm

			// === Diagnostics ===
			std::string algorithm_name = "QR";  // Algorithm identifier
			AlgorithmStatus status = AlgorithmStatus::Success;
			std::string error_message;          // Error description (empty on success)
			double elapsed_time_ms = 0.0;

			Result(int n)
				: eigenvalues(n)
				, eigenvectors(n, n) {}
			
			Result() = default;
		};

		// Main interface: Solve symmetric eigenvalue problem
		static Result Solve(const MatrixSym<Real>& A, Real tol = PrecisionValues<Real>::EigenSolverConvergenceTolerance, int maxIter = 1000) {
			int n = A.rows();
			Result result(n);

			if (n == 0) {
				result.converged = true;
				return result;
			}

			if (n == 1) {
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
			Matrix<Real> Q = Matrix<Real>::Identity(n);

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
			ImplicitQRAlgorithm(diag, subdiag, Q, n, tol, maxIter, result.iterations, result.converged);

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
		static Result Solve(const Matrix<Real>& A, Real tol = PrecisionValues<Real>::EigenSolverConvergenceTolerance, int maxIter = 1000) {
			int n = A.rows();
			MatrixSym<Real> symA(n);
			for (int i = 0; i < n; i++)
				for (int j = i; j < n; j++)
					symA(i, j) = 0.5 * (A[i][j] + A[j][i]);

			return Solve(symA, tol, maxIter);
		}

		/// Solve symmetric eigenvalue problem using configuration object.
		///
		/// @param A      Symmetric matrix
		/// @param config Solver configuration (tolerance, max iterations, etc.)
		/// @return Result with eigenvalues, eigenvectors, and full diagnostics
		static Result Solve(const MatrixSym<Real>& A, const EigenSolverConfig& config) {
			AlgorithmTimer timer;  // Starts automatically

			int maxIter = (config.max_iterations > 0) ? config.max_iterations : 1000;
			Result result = Solve(A, config.tolerance, maxIter);

			// Populate diagnostic fields
			result.elapsed_time_ms = timer.elapsed_ms();
			if (!result.converged) {
				result.status = AlgorithmStatus::MaxIterationsExceeded;
				result.error_message = "QR iteration did not converge within " + std::to_string(maxIter) + " iterations";
			}

			return result;
		}

		/// Solve for a regular Matrix using configuration object.
		static Result Solve(const Matrix<Real>& A, const EigenSolverConfig& config) {
			int n = A.rows();
			MatrixSym<Real> symA(n);
			for (int i = 0; i < n; i++)
				for (int j = i; j < n; j++)
					symA(i, j) = 0.5 * (A[i][j] + A[j][i]);

			return Solve(symA, config);
		}

	private:
		/// Reduce symmetric matrix to tridiagonal form using Householder reflections.
		///
		/// Algorithm: For k = 0, 1, ..., n-3:
		/// 1. Compute Householder vector v to zero out T(k+2:n, k)
		/// 2. Apply transformation: T = H * T * H where H = I - 2*v*v^T
		/// 3. Accumulate: Q = Q * H
		///
		/// Result: T becomes tridiagonal, Q accumulates the orthogonal transformations
		///
		/// Based on Golub & Van Loan, "Matrix Computations", Algorithm 8.3.1
		static void TridiagonalReduction(Matrix<Real>& T, Matrix<Real>& Q, int n) {
			for (int k = 0; k < n - 2; k++) {
				// Compute ||x||^2 where x = T(k+1:n-1, k)
				Real sigma = 0.0;
				for (int i = k + 1; i < n; i++)
					sigma += T(i, k) * T(i, k);

				if (sigma < PrecisionValues<Real>::DivisionSafetyThreshold) // Column already zeroed
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
				for (int i = k + 2; i < n; i++) {
					T(i, k) = 0.0;
					T(k, i) = 0.0;
				}

				// Accumulate transformation: Q = Q * H = Q * (I - 2*v*v^T)
				// Q_new[:, j] = Q[:, j] - 2*(Q*v)*v[j]
				for (int i = 0; i < n; i++) {
					Real dot = 0.0;
					for (int m = k + 1; m < n; m++)
						dot += Q(i, m) * v[m];

					for (int j = k + 1; j < n; j++)
						Q(i, j) -= 2.0 * dot * v[j];
				}
			}
		}

		/// Implicit QR algorithm with Wilkinson shift for tridiagonal matrices.
		///
		/// This is the workhorse of the eigenvalue computation. Uses:
		/// - Wilkinson shift for cubic convergence (default)
		/// - Francis double shift using both trailing 2x2 eigenvalues (when stagnating)
		/// - Exceptional shifts to break convergence cycles
		/// - Implicit QR step via bulge chasing with Givens rotations
		/// - Deflation when subdiagonal elements become negligible
		static void ImplicitQRAlgorithm(Vector<Real>& d, Vector<Real>& e, Matrix<Real>& Q, int n, Real tol, int maxIter, int& iterations,
										bool& converged) {
			iterations = 0;
			converged = false;

			// Work on unreduced portion [lo, hi]
			int lo = 0;
			int hi = n - 1;
			int iterSinceDeflation = 0;

			while (hi > 0 && iterations < maxIter) {
				// Deflate: Find the largest unreduced block
				// Check if e[hi-1] is negligible
				Real threshold = tol * (std::abs(d[hi - 1]) + std::abs(d[hi]));
				if (threshold < tol)
					threshold = tol;

				if (std::abs(e[hi - 1]) <= threshold) {
					e[hi - 1] = 0.0;
					hi--;
					iterSinceDeflation = 0;
					continue;
				}

				// Find lo: start of the unreduced block ending at hi
				lo = hi - 1;
				while (lo > 0) {
					threshold = tol * (std::abs(d[lo - 1]) + std::abs(d[lo]));
					if (threshold < tol)
						threshold = tol;

					if (std::abs(e[lo - 1]) <= threshold) {
						e[lo - 1] = 0.0;
						break;
					}
					lo--;
				}

				// Now work on the unreduced block [lo, hi]
				// Choose shift strategy based on convergence progress

				if (iterSinceDeflation >= 20) {
					// Exceptional shift: ad-hoc perturbation to break convergence cycles
					// Alternating sign prevents getting stuck in periodic orbits
					Real shift;
					if (iterSinceDeflation % 2 == 0)
						shift = d[hi] + 1.5 * std::abs(e[hi - 1]);
					else
						shift = d[hi] - 1.5 * std::abs(e[hi - 1]);

					ImplicitQRStep(d, e, Q, lo, hi, shift, n);
					iterations++;
					iterSinceDeflation++;
				}
				else if (iterSinceDeflation >= 10 && hi - lo >= 2) {
					// Francis double shift: use both eigenvalues of trailing 2x2
					// Helps with clustered or near-degenerate eigenvalue pairs
					Real s1, s2;
					TrailingEigenvalues(d[hi - 1], e[hi - 1], d[hi], s1, s2);

					ImplicitQRStep(d, e, Q, lo, hi, s1, n);
					iterations++;
					iterSinceDeflation++;

					// Quick deflation check between shifts
					threshold = tol * (std::abs(d[hi - 1]) + std::abs(d[hi]));
					if (threshold < tol) threshold = tol;
					if (std::abs(e[hi - 1]) <= threshold) {
						e[hi - 1] = 0.0;
						hi--;
						iterSinceDeflation = 0;
						continue;
					}

					ImplicitQRStep(d, e, Q, lo, hi, s2, n);
					iterations++;
					iterSinceDeflation++;
				}
				else {
					// Standard Wilkinson shift (cubic convergence for isolated eigenvalues)
					Real shift = WilkinsonShift(d[hi - 1], e[hi - 1], d[hi]);
					ImplicitQRStep(d, e, Q, lo, hi, shift, n);
					iterations++;
					iterSinceDeflation++;
				}
			}

			converged = (hi <= 0) || IsConverged(e, n, tol);
		}

		/// Compute Wilkinson shift: eigenvalue of 2x2 trailing submatrix closer to d[hi]
		///
		/// For matrix [[a, b], [b, c]], eigenvalues are:
		/// λ = (a+c)/2 ± sqrt(((a-c)/2)^2 + b^2)
		///
		/// We want the one closer to c (= d[hi])
		static Real WilkinsonShift(Real a, Real b, Real c) {
			Real d = (a - c) * 0.5;
			Real b_sq = b * b;

			if (std::abs(d) < PrecisionValues<Real>::DivisionSafetyThreshold)
				return c - std::abs(b);

			Real sign_d = (d >= 0.0) ? 1.0 : -1.0;
			Real r = std::sqrt(d * d + b_sq);

			// Shift = c - b^2 / (d + sign(d)*sqrt(d^2+b^2))
			return c - b_sq / (d + sign_d * r);
		}

		/// Compute both eigenvalues of the trailing 2x2 block [[a, b], [b, c]].
		///
		/// Returns both eigenvalues with s1 being the Wilkinson shift (closer to c = d[hi])
		/// and s2 being the other eigenvalue. Used for Francis double-shift strategy.
		static void TrailingEigenvalues(Real a, Real b, Real c, Real& s1, Real& s2) {
			Real avg = (a + c) * 0.5;
			Real delta = (a - c) * 0.5;
			Real r = std::sqrt(delta * delta + b * b);
			s1 = avg - r;
			s2 = avg + r;

			// Make s1 the Wilkinson shift (closer to c = d[hi])
			if (std::abs(s1 - c) > std::abs(s2 - c))
				std::swap(s1, s2);
		}

		/// Implicit QR step with shift using Givens rotations (bulge chasing).
		///
		/// Creates a bulge at position (1,0) by introducing shift, then chases
		/// it down the diagonal using Givens rotations.
		///
		/// For Givens G = [c s; -s c], the similarity transformation G*T*G^T
		/// on the 2x2 block [d_k e_k; e_k d_{k+1}] gives:
		/// d'_k     = c²*d_k + 2*c*s*e_k + s²*d_{k+1}
		/// d'_{k+1} = s²*d_k - 2*c*s*e_k + c²*d_{k+1}
		/// e'_k     = c*s*(d_{k+1} - d_k) + (c² - s²)*e_k
		static void ImplicitQRStep(Vector<Real>& d, Vector<Real>& e, Matrix<Real>& Q, int lo, int hi, Real shift, int n) {
			// Initial bulge: the first rotation is applied to [d[lo]-shift; e[lo]]
			Real x = d[lo] - shift;
			Real z = e[lo];

			for (int k = lo; k < hi; k++) {
				// Compute Givens rotation to eliminate z
				// G = [c s; -s c] such that G * [x; z] = [r; 0]
				Real c, s;
				ComputeGivens(x, z, c, s);

				// Update the previous off-diagonal (if not first iteration)
				if (k > lo) {
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
				if (k < hi - 1) {
					// The bulge element is s * e[k+1] (from applying G to column k+1)
					x = e[k];
					z = s * e[k + 1];
					e[k + 1] = c * e[k + 1];
				}

				// Accumulate rotation in eigenvector matrix
				// We're computing Q * G^T (to get eigenvectors of original matrix)
				// Q_new[:, k] = c*Q[:, k] + s*Q[:, k+1]
				// Q_new[:, k+1] = -s*Q[:, k] + c*Q[:, k+1]
				for (int i = 0; i < n; i++) {
					Real qik = Q(i, k);
					Real qik1 = Q(i, k + 1);
					Q(i, k) = c * qik + s * qik1;
					Q(i, k + 1) = -s * qik + c * qik1;
				}
			}
		}

		/// Compute Givens rotation coefficients c and s such that:
		/// [c  s] [a]   [r]
		/// [-s c] [b] = [0]
		///
		/// where r = sqrt(a^2 + b^2) (or -sqrt if needed for sign consistency)
		///
		/// Standard formulas: c = a/r, s = b/r
		/// This zeros out b: c*a + s*b = (a^2 + b^2)/r = r
		/// -s*a + c*b = (-ab + ab)/r = 0
		static void ComputeGivens(Real a, Real b, Real& c, Real& s) {
			if (std::abs(b) < PrecisionValues<Real>::DivisionSafetyThreshold) {
				c = (a >= 0.0) ? 1.0 : -1.0;
				s = 0.0;
			} else if (std::abs(a) < PrecisionValues<Real>::DivisionSafetyThreshold) {
				c = 0.0;
				s = (b >= 0.0) ? 1.0 : -1.0;
			} else if (std::abs(b) > std::abs(a)) {
				// Use t = a/b to avoid overflow
				Real t = a / b;
				Real u = std::sqrt(1.0 + t * t);
				if (b < 0.0)
					u = -u;
				s = 1.0 / u;
				c = s * t;
			} else {
				// Use t = b/a to avoid overflow
				Real t = b / a;
				Real u = std::sqrt(1.0 + t * t);
				if (a < 0.0)
					u = -u;
				c = 1.0 / u;
				s = c * t;
			}
		}

		/// Check if the tridiagonal matrix has converged (all subdiagonals negligible)
		static bool IsConverged(const Vector<Real>& e, int n, Real tol) {
			for (int i = 0; i < n - 1; i++) {
				if (std::abs(e[i]) > tol)
					return false;
			}
			return true;
		}

		/// Sort eigenvalues and corresponding eigenvectors in ascending order
		static void SortEigenpairs(Vector<Real>& eigenvalues, Matrix<Real>& eigenvectors, int n) {
			// Selection sort (n is typically small)
			for (int i = 0; i < n - 1; i++) {
				int minIdx = i;
				for (int j = i + 1; j < n; j++) {
					if (eigenvalues[j] < eigenvalues[minIdx])
						minIdx = j;
				}

				if (minIdx != i) {
					// Swap eigenvalues
					std::swap(eigenvalues[i], eigenvalues[minIdx]);

					// Swap corresponding eigenvectors (columns)
					for (int row = 0; row < n; row++)
						std::swap(eigenvectors(row, i), eigenvectors(row, minIdx));
				}
			}
		}
	};

} // namespace MML

#endif // MML_SYMM_MAT_EIGEN_SOLVER_QR_H

#ifndef MML_SYMM_MAT_EIGEN_SOLVER_JACOBI_H
#define MML_SYMM_MAT_EIGEN_SOLVER_JACOBI_H

#include <mml/MMLBase.h>
#include <mml/base/Vector/Vector.h>
#include <mml/base/Matrix/Matrix.h>
#include <mml/base/Matrix/MatrixSym.h>
#include <mml/base/AlgorithmTypes.h>
#include <mml/algorithms/Eigen/EigenSolverConfig.h>
#include <mml/algorithms/Eigen/detail/SelfAdjointScaling.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <string>
#include <utility>

namespace MML {

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
	class SymmMatEigenSolverJacobi {
	public:
		/// Result structure for Jacobi eigenvalue decomposition.
		///
		/// Contains eigenvalues, eigenvectors, and diagnostic information.
		/// Always check `converged` before using the results.
		struct Result {
			// === Primary Output ===
			Vector<Real> eigenvalues;  // Eigenvalues (sorted ascending by default)
			Matrix<Real> eigenvectors; // Column i = eigenvector for eigenvalue i

			// === Convergence Status ===
			bool converged = false;    // True if converged within tolerance
			int iterations = 0;        // Number of sweeps performed
			Real residual = 0.0;       // Final off-diagonal norm (achieved tolerance)

			// === Diagnostics ===
			std::string algorithm_name = "Jacobi";  // Algorithm identifier
			AlgorithmStatus status = AlgorithmStatus::Success;
			std::string error_message;              // Error description (empty on success)
			double elapsed_time_ms = 0.0;

			Result(int n, bool computeEigenvectors = true)
				: eigenvalues(n)
				, eigenvectors(computeEigenvectors ? n : 0, computeEigenvectors ? n : 0) {}
			
			Result() = default;
		};

		/// Solve symmetric eigenvalue problem: A*v = λ*v
		///
		/// @param A        Symmetric packed matrix
		/// @param tol      Convergence tolerance (default: 1e-10)
		/// @param maxIter  Maximum number of sweeps (default: 100)
		/// @return Result containing eigenvalues, eigenvectors, and convergence info
		///
		/// POSTCONDITIONS:
		/// - eigenvalues[i] <= eigenvalues[i+1] (ascending order)
		/// - eigenvectors.Column(i) corresponds to eigenvalues[i]
		/// - Eigenvectors are orthonormal: V^T * V = I
		/// - A * V = V * diag(λ)
		/// - `converged` is true only when every eigenvalue and eigenvector entry is finite
		///
		/// Matrices with extreme magnitudes (beyond roughly 1e150 or below 1e-150) are scaled by a
		/// power of two before iterating so that squares and sums of entries cannot overflow; the
		/// tolerance is interpreted in the units of the original matrix.
		static Result Solve(const MatrixSym<Real>& A, Real tol = PrecisionValues<Real>::EigenSolverConvergenceTolerance, int maxIter = 100) {
			return SolveImpl(A, tol, maxIter, true, true);
		}

	private:
		static Result SolveImpl(const MatrixSym<Real>& A, Real tol, int maxIter, bool sortEigenvalues, bool computeEigenvectors) {
			detail::ValidateSelfAdjointControls(tol, maxIter);
			int n = A.rows();
			Result result(n, computeEigenvectors);
			if (n == 0) {
				result.status = AlgorithmStatus::Completed;
				result.residual = NoEstimate;
				return result;
			}

			// Initialize: Copy A to working matrix, V to identity
			Matrix<Real> D(n, n);
			for (int i = 0; i < n; i++)
				for (int j = 0; j < n; j++)
					D[i][j] = (i <= j) ? A(i, j) : A(j, i);

			const Real maxAbs = detail::MaxAbsEntry(D);
			if (!std::isfinite(maxAbs)) {
				result.eigenvalues = Vector<Real>(n, std::numeric_limits<Real>::quiet_NaN());
				if (computeEigenvectors)
					result.eigenvectors = Matrix<Real>(n, n, std::numeric_limits<Real>::quiet_NaN());
				result.residual = std::numeric_limits<Real>::quiet_NaN();
				result.converged = false;
				result.status = AlgorithmStatus::NumericalInstability;
				result.error_message = "SymmMatEigenSolverJacobi: matrix contains nonfinite entries";
				return result;
			}

			// Overflow guard: work on D/scale (exact power-of-two scaling). Entries that are skipped
			// as "already zero" must stay within the convergence tolerance, so the skip threshold
			// is capped by the (scaled) tolerance; exact zeros are always skipped.
			const Real scale = detail::SelfAdjointScale(maxAbs);
			const Real tolScaled = tol / scale;
			const Real zeroThresholdScaled = std::min(PrecisionValues<Real>::EigenSolverZeroThreshold, tolScaled);
			if (scale != 1)
				for (int i = 0; i < n; i++)
					for (int j = 0; j < n; j++)
						D[i][j] /= scale;

			Matrix<Real> V;
			if (computeEigenvectors) V = Matrix<Real>::Identity(n);

			// Jacobi iterations
			int sweep = 0;
			Real offDiagNorm = detail::OffDiagonalNormSafe(D);

			while (offDiagNorm > tolScaled && sweep < maxIter) {
				// One sweep: process all off-diagonal elements
				for (int p = 0; p < n - 1; p++) {
					for (int q = p + 1; q < n; q++) {
						if (std::abs(D[p][q]) <= zeroThresholdScaled) // Skip if already zero
							continue;
						// Compute Givens rotation parameters
						Real theta = ComputeRotationAngle(D, p, q, zeroThresholdScaled);
						Real c = std::cos(theta);
						Real s = std::sin(theta);

						// Apply rotation
						ApplyRotation(D, V, p, q, c, s, n);
					}
				}

				sweep++;
				offDiagNorm = detail::OffDiagonalNormSafe(D);
			}

			// Extract eigenvalues from diagonal (undo scaling)
			for (int i = 0; i < n; i++)
				result.eigenvalues[i] = D[i][i] * scale;

			// Copy eigenvectors
			result.eigenvectors = V;

			// Sort eigenvalues and eigenvectors
			if (sortEigenvalues) SortEigenvalues(result.eigenvalues, result.eigenvectors);

			// Set convergence info
			result.iterations = sweep;
			result.residual = offDiagNorm * scale;
			result.converged = (offDiagNorm <= tolScaled);
			if (!detail::EigenResultIsFinite(result.eigenvalues, result.eigenvectors)) {
				result.converged = false;
				result.status = AlgorithmStatus::NumericalInstability;
				result.error_message = "SymmMatEigenSolverJacobi: iteration produced nonfinite eigenvalues or eigenvectors";
			} else if (!result.converged) {
				result.status = AlgorithmStatus::MaxIterationsExceeded;
				result.error_message = "Jacobi iteration did not converge within " + std::to_string(maxIter) + " sweeps";
			}

			return result;
		}

	public:
		/// Solve for a regular Matrix - validates square shape and symmetry.
		/// @throws MatrixDimensionError if A is not symmetric within EigenSolverZeroThreshold
		/// @see SolveSymmetricPart to solve for the symmetric part (A + A^T)/2 of a
		///      nearly-symmetric matrix without validation
		static Result Solve(const Matrix<Real>& A, Real tol = PrecisionValues<Real>::EigenSolverConvergenceTolerance, int maxIter = 100) {
			detail::ValidateSymmetricInput(A);
			return SolveSymmetricPart(A, tol, maxIter);
		}

		/// Solve for the symmetric part (A + A^T)/2 of a general matrix (no symmetry validation)
		static Result SolveSymmetricPart(const Matrix<Real>& A, Real tol = PrecisionValues<Real>::EigenSolverConvergenceTolerance, int maxIter = 100) {
			return Solve(SymmetricPart(A), tol, maxIter);
		}

		static Result SolveSymmetricPart(const Matrix<Real>& A, const EigenSolverConfig& config) {
			return Solve(SymmetricPart(A), config);
		}

		/// Solve symmetric eigenvalue problem using configuration object.
		///
		/// @param A      Symmetric matrix
		/// @param config Solver configuration (tolerance, max iterations, etc.)
		/// @return Result with eigenvalues, eigenvectors, and full diagnostics
		static Result Solve(const MatrixSym<Real>& A, const EigenSolverConfig& config) {
			AlgorithmTimer timer;  // Starts automatically

			detail::ValidateSelfAdjointControls(config.tolerance, config.max_iterations);
			int maxIter = config.max_iterations == 0 ? 100 : config.max_iterations;
			Result result = SolveImpl(A, config.tolerance, maxIter, config.sort_eigenvalues, config.compute_eigenvectors);

			// Populate diagnostic fields
			result.elapsed_time_ms = timer.elapsed_ms();
			if (!result.converged && result.status == AlgorithmStatus::Success) {
				result.status = AlgorithmStatus::MaxIterationsExceeded;
				result.error_message = "Jacobi iteration did not converge within " + std::to_string(maxIter) + " sweeps";
			}
			if (config.verbose) detail::LogSelfAdjointResult(result);

			return result;
		}

		/// Solve for a regular Matrix using configuration object.
		static Result Solve(const Matrix<Real>& A, const EigenSolverConfig& config) {
			detail::ValidateSymmetricInput(A);
			return Solve(SymmetricPart(A), config);
		}

	private:
		// Overflow-safe symmetric part (A + A^T)/2
		static MatrixSym<Real> SymmetricPart(const Matrix<Real>& A) {
			if (A.rows() != A.cols())
				throw MatrixDimensionError("SolveSymmetricPart requires a square matrix", A.rows(), A.cols(), -1, -1);
			int n = A.rows();
			MatrixSym<Real> symA(n);
			for (int i = 0; i < n; i++)
				for (int j = i; j < n; j++)
					symA(i, j) = detail::SymmetrizeEntries(A[i][j], A[j][i]);
			return symA;
		}

		// Compute rotation angle to zero out A[p][q]
		static Real ComputeRotationAngle(const Matrix<Real>& A, int p, int q, Real zeroThreshold) {
			Real theta;

			if (std::abs(A[p][q]) <= zeroThreshold) {
				theta = 0.0;
			} else {
				Real tau = (A[q][q] - A[p][p]) / (MML::Real(2.0) * A[p][q]);
				Real t;

				// Choose sign to avoid loss of significance
				if (tau >= 0)
					t = MML::Real(1.0) / (tau + std::sqrt(MML::Real(1.0) + tau * tau));
				else
					t = -MML::Real(1.0) / (-tau + std::sqrt(MML::Real(1.0) + tau * tau));

				theta = std::atan(t);
			}

			return theta;
		}

		// Apply Givens rotation J(p,q,θ) to A: A' = J^T * A * J
		static void ApplyRotation(Matrix<Real>& A, Matrix<Real>& V, int p, int q, Real c, Real s, int n) {
			// Update A
			for (int i = 0; i < n; i++) {
				if (i != p && i != q) {
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

			A[p][p] = c * c * App - MML::Real(2.0) * s * c * Apq + s * s * Aqq;
			A[q][q] = s * s * App + MML::Real(2.0) * s * c * Apq + c * c * Aqq;
			A[p][q] = 0.0;
			A[q][p] = 0.0;

			// Accumulate rotation in eigenvector matrix V
			for (int i = 0; i < V.rows(); i++) {
				Real Vip = V[i][p];
				Real Viq = V[i][q];
				V[i][p] = c * Vip - s * Viq;
				V[i][q] = s * Vip + c * Viq;
			}
		}

		// Sort eigenvalues and eigenvectors in ascending order
		static void SortEigenvalues(Vector<Real>& eigenvalues, Matrix<Real>& eigenvectors) {
			int n = eigenvalues.size();

			// Simple selection sort (eigenvalues should be small n)
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
					for (int row = 0; row < eigenvectors.rows(); row++)
						std::swap(eigenvectors[row][i], eigenvectors[row][minIdx]);
				}
			}
		}
	};

} // namespace MML

#endif // MML_SYMM_MAT_EIGEN_SOLVER_JACOBI_H

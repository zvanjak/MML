#ifndef MML_EIGEN_SOLVER_H
#define MML_EIGEN_SOLVER_H

#include <mml/MMLBase.h>
#include <mml/base/Vector/Vector.h>
#include <mml/base/Matrix/Matrix.h>
#include <mml/base/Matrix/MatrixSym.h>
#include <mml/base/AlgorithmTypes.h>
#include <mml/algorithms/Eigen/EigenSolverConfig.h>
#include <mml/algorithms/Eigen/HessenbergReduction.h>
#include <mml/algorithms/Eigen/detail/RealSchurAnalysis.h>
#include <mml/algorithms/Eigen/detail/HessenbergQRIteration.h>
#include <mml/algorithms/Eigen/detail/SchurEigenvectors.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <ostream>
#include <string>
#include <utility>

namespace MML {

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
	class EigenSolver {
	public:
		/// @struct ComplexEigenvalue
		/// @brief Represents a complex eigenvalue λ = real + imag*i
		struct ComplexEigenvalue {
			Real real;
			Real imag;

			ComplexEigenvalue(Real r = 0.0, Real i = 0.0)
				: real(r)
				, imag(i) {}

			bool isComplex(Real tol = 1e-12) const { return std::abs(imag) > tol; }
			Real magnitude() const { return std::sqrt(real * real + imag * imag); }

			friend std::ostream& operator<<(std::ostream& os, const ComplexEigenvalue& value) {
				return os << value.real << (value.imag < REAL(0.0) ? " - " : " + ") << std::abs(value.imag) << 'i';
			}
		};

		/// @struct Result
		/// @brief Complete eigensolution result with diagnostic information
		struct Result {
			// === Primary Output ===
			std::vector<ComplexEigenvalue> eigenvalues;
			Matrix<Real> eigenvectors;		 // Column i is eigenvector for eigenvalue i
			std::vector<bool> isComplexPair; // True if columns i,i+1 form complex pair

			// === Convergence Status ===
			bool converged = false;
			int iterations = 0;
			Real maxResidual = 0.0;

			// === Diagnostics ===
			std::string algorithm_name = "GeneralQR";
			AlgorithmStatus status = AlgorithmStatus::Success;
			std::string error_message;
			double elapsed_time_ms = 0.0;

			Result(int n) : eigenvectors(n, n) {}
			Result() = default;
		};

		/// Solve the complete eigenvalue problem for a general matrix.
		///
		/// @param A Input matrix (n x n)
		/// @param tol Convergence tolerance (default: 1e-10)
		/// @param maxIter Maximum total QR iterations (default: 1000)
		/// @return Result with all eigenvalues and eigenvectors
		///
		/// POSTCONDITIONS:
		/// - eigenvalues.size() == n
		/// - eigenvectors has n columns
		/// - For real eigenvalue λ: A*v ≈ λ*v
		/// - For complex pair λ±μi with eigenvector (vr ± i*vi):
		/// - eigenvalues[k] = {λ, μ}, eigenvalues[k+1] = {λ, -μ}
		/// - eigenvectors.Column(k) = vr, Column(k+1) = vi
		static Result Solve(const Matrix<Real>& A, Real tol = PrecisionValues<Real>::EigenSolverConvergenceTolerance, int maxIter = 1000) {
			int n = A.rows();
			Result result(n);

			if (n == 0) {
				result.converged = true;
				return result;
			}

			if (n == 1) {
				result.eigenvalues.push_back(ComplexEigenvalue(A(0, 0), 0.0));
				result.eigenvectors(0, 0) = 1.0;
				result.isComplexPair.push_back(false);
				result.converged = true;
				return result;
			}

			// Step 1: Reduce to Hessenberg form
			auto hess = ReduceToHessenberg(A);
			Matrix<Real> H = hess.H;
			Matrix<Real> Q = hess.Q; // Accumulated transformation

			// Step 2: QR iteration with deflation
			int iLo = 0;	 // Start of active region
			int iHi = n - 1; // End of active region
			int totalIter = 0;

			while (iHi > iLo && totalIter < maxIter) {
				// Check for deflation at bottom of active region
				auto defl = detail::CheckDeflation(H, iLo, iHi, tol);

				if (defl.canDeflate) {
					// Apply deflation
					detail::ApplyDeflation(H, defl.deflationIndex);

					// Shrink active region
					if (defl.deflationIndex == iHi) {
						iHi--;
					} else if (defl.deflationIndex == iHi - 1 && defl.blockSize == 2) {
						iHi -= 2; // 2x2 block deflated
					} else {
						// Middle deflation - just continue
						iHi = defl.deflationIndex - 1;
					}
					continue;
				}

				// Determine if we should use single or double shift
				// Check if bottom 2x2 has complex eigenvalues
				bool useDoubleShift = false;
				if (iHi >= iLo + 1) {
					Real a = H(iHi - 1, iHi - 1);
					Real b = H(iHi - 1, iHi);
					Real c = H(iHi, iHi - 1);
					Real d = H(iHi, iHi);
					auto eig = detail::Eigenvalues2x2(a, b, c, d);
					useDoubleShift = eig.isComplex;

					// Special case: if we have only a 2x2 block with complex eigenvalues,
					// it's already in final form - we're done with this block
					if (useDoubleShift && iHi == iLo + 1) {
						// 2x2 block with complex eigenvalues is irreducible - accept it
						break;
					}
				}

				if (useDoubleShift && iHi >= iLo + 2) {
					// Apply Francis double-shift
					auto step = detail::FrancisDoubleShift(H, iLo, iHi, true);
					H = step.H;

					// Update accumulated Q
					Matrix<Real> newQ(n, n);
					for (int i = 0; i < n; i++)
						for (int j = 0; j < n; j++) {
							newQ(i, j) = 0.0;
							for (int k = 0; k < n; k++)
								newQ(i, j) += Q(i, k) * step.Q(k, j);
						}
					Q = newQ;
					totalIter += 2; // Double shift counts as 2 iterations
				} else {
					// Apply single Wilkinson shift
					auto step = detail::SingleQRStep(H, true);
					H = step.H;

					// Update accumulated Q
					Matrix<Real> newQ(n, n);
					for (int i = 0; i < n; i++)
						for (int j = 0; j < n; j++) {
							newQ(i, j) = 0.0;
							for (int k = 0; k < n; k++)
								newQ(i, j) += Q(i, k) * step.Q(k, j);
						}
					Q = newQ;
					totalIter++;
				}

				// Apply exceptional shift if stuck (every 30 iterations without progress)
				if (totalIter % 30 == 29 && iHi >= 2) {
					// Exceptional shift to break cycles
					Real exceptionalShift = std::abs(H(iHi, iHi - 1)) + std::abs(H(iHi - 1, iHi - 2));
					auto step = detail::SingleQRStep(H, true, &exceptionalShift);
					H = step.H;

					Matrix<Real> newQ(n, n);
					for (int i = 0; i < n; i++)
						for (int j = 0; j < n; j++) {
							newQ(i, j) = 0.0;
							for (int k = 0; k < n; k++)
								newQ(i, j) += Q(i, k) * step.Q(k, j);
						}
					Q = newQ;
					totalIter++;
				}
			}

			// Step 3: Extract eigenvalues from quasi-triangular form
			auto eigResult = detail::ExtractEigenvalues(H, tol);

			// Convert to our ComplexEigenvalue format
			for (const auto& e : eigResult.eigenvalues)
				result.eigenvalues.push_back(ComplexEigenvalue(e.real, e.imag));

			// Step 4: Compute eigenvectors
			auto vecResult = detail::ComputeEigenvectorsFromSchur(H, Q, tol);
			result.eigenvectors = vecResult.vectors;
			result.isComplexPair = vecResult.isComplexPair;

			// Step 5: Compute residuals and check convergence
			result.iterations = totalIter;
			result.maxResidual = 0.0;
			result.converged = (totalIter < maxIter);

			// Verify eigenvector quality
			for (size_t i = 0; i < result.eigenvalues.size(); i++) {
				if (!result.isComplexPair[i] || (i > 0 && result.isComplexPair[i - 1])) {
					Vector<Real> v(n);
					for (int j = 0; j < n; j++)
						v[j] = result.eigenvectors(j, static_cast<int>(i));

					Real residual = detail::EigenvectorResidual(A, v, result.eigenvalues[i].real, result.eigenvalues[i].imag);
					result.maxResidual = std::max(result.maxResidual, residual);
				}
			}

			return result;
		}

		/// Solve the complete eigenvalue problem using configuration object.
		///
		/// @param A      General (possibly non-symmetric) matrix
		/// @param config Solver configuration (tolerance, max iterations, etc.)
		/// @return Result with eigenvalues, eigenvectors, and full diagnostics
		static Result Solve(const Matrix<Real>& A, const EigenSolverConfig& config) {
			AlgorithmTimer timer;  // Starts automatically

			int maxIter = (config.max_iterations > 0) ? config.max_iterations : 1000;
			Result result = Solve(A, config.tolerance, maxIter);

			// Populate diagnostic fields
			result.elapsed_time_ms = timer.elapsed_ms();
			if (!result.converged) {
				result.status = AlgorithmStatus::MaxIterationsExceeded;
				result.error_message = "General QR iteration did not converge within " + std::to_string(maxIter) + " iterations";
			}

			return result;
		}
	};

} // namespace MML

#endif // MML_EIGEN_SOLVER_H

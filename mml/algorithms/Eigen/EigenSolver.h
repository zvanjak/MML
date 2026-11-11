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
   * Complete eigensolver for general (non-symmetric) matrices using Francis implicit
   * double-shift QR iteration on the Hessenberg form (EISPACK hqr2 / LAPACK dlahqr structure).
   * 
   * ALGORITHM:
   * 1. Reduce matrix to upper Hessenberg form H = Q^T * A * Q
   * 2. Repeatedly locate the lowest unreduced diagonal block [l, iHi], deflate converged 1x1 / 2x2
   *    blocks from the bottom (real 2x2 blocks are triangularised), otherwise apply one implicit
   *    double-shift step to that block only, until H is quasi-upper-triangular (real Schur form)
   * 3. Extract eigenvalues from the 1x1 and 2x2 diagonal blocks
   * 4. Compute eigenvectors by quasi-triangular back-substitution on T, then X = Q * Y
   * 5. Report ||A v - λ v|| / ||v|| over all eigenpairs as maxResidual
   * 
   * FEATURES:
   * - Handles both real and complex eigenvalues
   * - Complex eigenvalues are stored as conjugate pairs (+β first)
   * - Shifts are always taken from the active block, so convergence is quadratic/cubic
   * - Exceptional shifts after 10 and 20 fruitless iterations on one eigenvalue
   * - converged is true only when every block was deflated within maxIter
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

			bool isComplex(Real tol = Precision::ZeroTolerance) const { return std::abs(imag) > tol; }
			Real magnitude() const { return std::sqrt(real * real + imag * imag); }

			friend std::ostream& operator<<(std::ostream& os, const ComplexEigenvalue& value) {
				return os << value.real << (value.imag < static_cast<MML::Real>(0.0) ? " - " : " + ") << std::abs(value.imag) << 'i';
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
			if (A.rows() != A.cols())
				throw MatrixDimensionError("EigenSolver::Solve - matrix must be square", A.rows(), A.cols(), -1, -1);
			if (!std::isfinite(tol) || tol < static_cast<MML::Real>(0.0))
				throw DomainError("EigenSolver::Solve - tolerance must be finite and nonnegative");

			int n = A.rows();
			Result result(n);

			if (n == 0) {
				// Empty problem: nothing to iterate or check (RESULT_POLICY: Completed)
				result.converged = false;
				result.status = AlgorithmStatus::Completed;
				result.maxResidual = NoEstimate;
				return result;
			}

			if (n == 1) {
				// Closed form: A v = A(0,0) v with v = 1, residual exactly 0 - the criterion is checked and met
				result.eigenvalues.push_back(ComplexEigenvalue(A(0, 0), 0.0));
				result.eigenvectors(0, 0) = 1.0;
				result.isComplexPair.push_back(false);
				result.maxResidual = 0.0;
				result.converged = result.maxResidual <= tol;
				return result;
			}

			// Step 1: Reduce to Hessenberg form  (H = Q^T A Q)
			auto hess = ReduceToHessenberg(A);
			Matrix<Real> H = hess.H;
			Matrix<Real> Q = hess.Q; // Accumulated transformation

			// Step 2: Francis double-shift QR iteration with deflation, in the structure of
			// EISPACK hqr / Numerical Recipes hqr2 / LAPACK dlahqr:
			//   - find the lowest unreduced block [l, iHi] (NOT the first small subdiagonal from the
			//     bottom: iterating on [0, iHi] with shifts from a converged block, or shrinking iHi
			//     past an interior split, both break convergence - analysis v4, ALGORITHMS A1);
			//   - deflate 1x1 and 2x2 blocks from the bottom, standardising real 2x2 blocks;
			//   - otherwise apply one implicit double-shift step to [l, iHi] only, with exceptional
			//     shifts after 10 and 20 fruitless iterations on the same eigenvalue.
			const Real anorm = detail::HessenbergAbsNorm(H);
			int iHi = n - 1;		// Bottom of the active region (everything below is converged)
			int its = 0;			// Iterations spent on the current eigenvalue
			int totalIter = 0;
			bool allDeflated = true;

			while (iHi >= 0) {
				const int l = detail::FindUnreducedBlockStart(H, 0, iHi, tol, anorm);
				if (l > 0)
					H(l, l - 1) = static_cast<MML::Real>(0.0); // make the split exact for extraction / back-substitution

				if (l == iHi) {
					// 1x1 block converged
					iHi--;
					its = 0;
					continue;
				}
				if (l == iHi - 1) {
					// 2x2 block converged: triangularise it if its eigenvalues are real
					detail::StandardizeSchurBlock(H, Q, l);
					iHi -= 2;
					its = 0;
					continue;
				}

				if (totalIter >= maxIter) {
					allDeflated = false;
					break;
				}

				// Exceptional shifts (LAPACK dlahqr form) to break cycles
				const Real* shiftTrace = nullptr;
				const Real* shiftProduct = nullptr;
				Real exTrace = static_cast<MML::Real>(0.0), exProduct = static_cast<MML::Real>(0.0);
				if (its == 10 || its == 20) {
					const Real s = std::abs(H(iHi, iHi - 1)) + std::abs(H(iHi - 1, iHi - 2));
					const Real h = static_cast<MML::Real>(0.75) * s + H(iHi, iHi);
					exTrace = static_cast<MML::Real>(2.0) * h;
					exProduct = h * h + static_cast<MML::Real>(0.4375) * s * s;
					shiftTrace = &exTrace;
					shiftProduct = &exProduct;
				}

				auto step = detail::FrancisDoubleShift(H, l, iHi, true, shiftTrace, shiftProduct);
				H = step.H;
				Q = Q * step.Q;
				its++;
				totalIter++;
			}

			// Step 3: Extract eigenvalues from the (now exactly) quasi-triangular form.
			// Deflated subdiagonals are exact zeros, so tol = 0 classifies blocks structurally.
			auto eigResult = detail::ExtractEigenvalues(H, static_cast<MML::Real>(0.0));
			for (const auto& e : eigResult.eigenvalues)
				result.eigenvalues.push_back(ComplexEigenvalue(e.real, e.imag));

			// Step 4: Compute eigenvectors (quasi-triangular back-substitution, then X = Q Y)
			auto vecResult = detail::ComputeEigenvectorsFromSchur(H, Q, static_cast<MML::Real>(0.0));
			result.eigenvectors = vecResult.vectors;
			result.isComplexPair = vecResult.isComplexPair;

			// Step 5: Residuals ||A v - λ v|| / ||v|| for every eigenpair (complex pairs included)
			result.iterations = totalIter;
			result.maxResidual = 0.0;
			{
				Vector<Real> vr(n), vi(n);
				size_t i = 0;
				while (i < result.eigenvalues.size()) {
					for (int j = 0; j < n; j++)
						vr[j] = result.eigenvectors(j, static_cast<int>(i));
					const bool pair = result.isComplexPair[i] && i + 1 < result.eigenvalues.size();
					for (int j = 0; j < n; j++)
						vi[j] = pair ? result.eigenvectors(j, static_cast<int>(i + 1)) : static_cast<MML::Real>(0.0);

					Real residual = detail::EigenpairResidual(A, vr, vi, result.eigenvalues[i].real,
															  pair ? result.eigenvalues[i].imag : static_cast<MML::Real>(0.0));
					result.maxResidual = std::max(result.maxResidual, residual);
					i += pair ? 2 : 1;
				}
			}

			result.converged = allDeflated && std::isfinite(result.maxResidual);
			if (!allDeflated) {
				result.status = AlgorithmStatus::MaxIterationsExceeded;
				result.error_message = "General QR iteration did not converge within " + std::to_string(maxIter) + " iterations";
			} else if (!std::isfinite(result.maxResidual)) {
				result.status = AlgorithmStatus::NumericalInstability;
				result.error_message = "Non-finite eigenvector residual";
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

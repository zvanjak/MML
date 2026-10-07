#ifndef MML_HERMITIAN_MAT_EIGEN_SOLVER_JACOBI_H
#define MML_HERMITIAN_MAT_EIGEN_SOLVER_JACOBI_H

#include <mml/MMLBase.h>
#include <mml/base/Vector/Vector.h>
#include <mml/base/Matrix/Matrix.h>
#include <mml/base/Matrix/MatrixSym.h>
#include <mml/base/AlgorithmTypes.h>
#include <mml/algorithms/Eigen/EigenSolverConfig.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <string>
#include <utility>

namespace MML {

	/// Jacobi eigensolver for complex Hermitian matrices.
	/// Produces real eigenvalues and unitary complex eigenvectors such that A = V diag(lambda) V*.
	class HermitianMatEigenSolverJacobi {
	public:
		struct Result {
			Vector<Real> eigenvalues;
			Matrix<Complex> eigenvectors;
			bool converged = false;
			int iterations = 0;
			Real residual = 0.0;
			std::string algorithm_name = "HermitianJacobi";
			AlgorithmStatus status = AlgorithmStatus::Success;
			std::string error_message;
			double elapsed_time_ms = 0.0;

			Result(int n)
				: eigenvalues(n)
				, eigenvectors(n, n) {}

			Result() = default;
		};

		static Result Solve(const Matrix<Complex>& A, Real tol = PrecisionValues<Real>::EigenSolverConvergenceTolerance, int maxIter = 100) {
			if (A.rows() != A.cols())
				throw MatrixDimensionError("HermitianMatEigenSolverJacobi::Solve - matrix must be square", A.rows(), A.cols(), A.rows(), A.rows());

			const int n = A.rows();
			const Real hermitianTolerance = PrecisionValues<Real>::EigenSolverZeroThreshold;
			for (int i = 0; i < n; i++) {
				if (std::abs(A[i][i].imag()) > hermitianTolerance)
					throw MatrixDimensionError("HermitianMatEigenSolverJacobi::Solve - matrix diagonal must be real", n, n, -1, -1);
				for (int j = i + 1; j < n; j++)
					if (std::abs(A[i][j] - std::conj(A[j][i])) > hermitianTolerance)
						throw MatrixDimensionError("HermitianMatEigenSolverJacobi::Solve - matrix is not Hermitian", n, n, -1, -1);
			}

			Result result(n);
			Matrix<Complex> D(A);
			Matrix<Complex> V = Matrix<Complex>::Identity(n);
			int sweep = 0;
			Real offDiagonalNorm = OffDiagonalNorm(D);

			while (offDiagonalNorm > tol && sweep < maxIter) {
				for (int p = 0; p < n - 1; p++) {
					for (int q = p + 1; q < n; q++) {
						const Real magnitude = std::abs(D[p][q]);
						if (magnitude < PrecisionValues<Real>::EigenSolverZeroThreshold)
							continue;

						const Real app = D[p][p].real();
						const Real aqq = D[q][q].real();
						const Real tau = (aqq - app) / (2.0 * magnitude);
						const Real t = tau >= 0.0
							? 1.0 / (tau + std::sqrt(1.0 + tau * tau))
							: -1.0 / (-tau + std::sqrt(1.0 + tau * tau));
						const Real c = 1.0 / std::sqrt(1.0 + t * t);
						const Real s = t * c;
						ApplyRotation(D, V, p, q, c, s, D[p][q] / magnitude);
					}
				}

				sweep++;
				offDiagonalNorm = OffDiagonalNorm(D);
			}

			for (int i = 0; i < n; i++)
				result.eigenvalues[i] = D[i][i].real();
			result.eigenvectors = V;
			SortEigenpairs(result.eigenvalues, result.eigenvectors);
			result.iterations = sweep;
			result.residual = offDiagonalNorm;
			result.converged = offDiagonalNorm <= tol;
			return result;
		}

		static Result Solve(const Matrix<Complex>& A, const EigenSolverConfig& config) {
			AlgorithmTimer timer;
			const int maxIter = config.max_iterations > 0 ? config.max_iterations : 100;
			Result result = Solve(A, config.tolerance, maxIter);
			result.elapsed_time_ms = timer.elapsed_ms();
			if (!result.converged) {
				result.status = AlgorithmStatus::MaxIterationsExceeded;
				result.error_message = "Hermitian Jacobi iteration did not converge within " + std::to_string(maxIter) + " sweeps";
			}
			return result;
		}

	private:
		static void ApplyRotation(Matrix<Complex>& A, Matrix<Complex>& V, int p, int q, Real c, Real s, Complex phase) {
			const int n = A.rows();
			for (int i = 0; i < n; i++) {
				if (i == p || i == q)
					continue;
				const Complex aip = A[i][p];
				const Complex aiq = A[i][q];
				A[i][p] = c * phase * aip - s * aiq;
				A[i][q] = s * phase * aip + c * aiq;
				A[p][i] = std::conj(A[i][p]);
				A[q][i] = std::conj(A[i][q]);
			}

			const Real app = A[p][p].real();
			const Real aqq = A[q][q].real();
			const Real apq = std::abs(A[p][q]);
			A[p][p] = Complex(c * c * app - 2.0 * s * c * apq + s * s * aqq, 0.0);
			A[q][q] = Complex(s * s * app + 2.0 * s * c * apq + c * c * aqq, 0.0);
			A[p][q] = A[q][p] = Complex{0.0, 0.0};

			for (int i = 0; i < n; i++) {
				const Complex vip = V[i][p];
				const Complex viq = V[i][q];
				V[i][p] = c * phase * vip - s * viq;
				V[i][q] = s * phase * vip + c * viq;
			}
		}

		static Real OffDiagonalNorm(const Matrix<Complex>& A) {
			Real sum = 0.0;
			for (int i = 0; i < A.rows() - 1; i++)
				for (int j = i + 1; j < A.cols(); j++)
					sum += std::norm(A[i][j]);
			return std::sqrt(sum);
		}

		static void SortEigenpairs(Vector<Real>& eigenvalues, Matrix<Complex>& eigenvectors) {
			const int n = eigenvalues.size();
			for (int i = 0; i < n - 1; i++) {
				int minIndex = i;
				for (int j = i + 1; j < n; j++)
					if (eigenvalues[j] < eigenvalues[minIndex])
						minIndex = j;

				if (minIndex != i) {
					std::swap(eigenvalues[i], eigenvalues[minIndex]);
					for (int row = 0; row < n; row++)
						std::swap(eigenvectors[row][i], eigenvectors[row][minIndex]);
				}
			}
		}
	};

} // namespace MML

#endif // MML_HERMITIAN_MAT_EIGEN_SOLVER_JACOBI_H

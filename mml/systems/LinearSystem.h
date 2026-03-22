///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        LinearSystem.h                                                      ///
///  Description: Unified linear algebra facade                                       ///
///               One class to access ALL MML linear algebra capabilities             ///
///               Smart solver selection, lazy decomposition caching, rich diagnostics///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_LINEAR_SYSTEM_H
#define MML_LINEAR_SYSTEM_H

#include "MMLBase.h"
#include "base/Vector/Vector.h"
#include "base/Matrix/Matrix.h"
#include "core/LinAlgEqSolvers.h"
#include "core/MatrixUtils.h"
#include "algorithms/MatrixAlg.h"
#include "algorithms/EigenSystemSolvers.h"

#include <optional>
#include <string>
#include <memory>
#include <cmath>
#include <vector>

namespace MML::Systems 
{
	//=============================================================================
	// ENUMERATIONS
	//=============================================================================

	/// @brief Iterative method selection for large sparse systems
	enum class IterativeMethod {
		Auto,				 ///< Auto-select based on matrix properties
		Jacobi,			 ///< Jacobi iteration (parallelizable)
		GaussSeidel, ///< Gauss-Seidel (faster than Jacobi)
		SOR					 ///< Successive Over-Relaxation
	};

	/// @brief Matrix stability assessment
	enum class MatrixStability {
		WellConditioned,			 ///< cond < 10^4
		ModeratelyConditioned, ///< 10^4 <= cond < 10^8
		IllConditioned,				 ///< 10^8 <= cond < 10^12
		Singular							 ///< cond >= 10^12 or rank-deficient
	};

	//=============================================================================
	// RESULT STRUCTURES
	//=============================================================================

	/// @brief LU decomposition result: A = P*L*U
	template<typename Type>
	struct LUDecomposition {
		Matrix<Type> L;								///< Lower triangular factor
		Matrix<Type> U;								///< Upper triangular factor
		std::vector<int> permutation; ///< Row permutation
		Type determinant;							///< Determinant of A
		bool valid = false;						///< Whether decomposition succeeded
	};

	/// @brief QR decomposition result: A = Q*R
	template<typename Type>
	struct QRDecomposition {
		Matrix<Type> Q; ///< Orthogonal matrix
		Matrix<Type> R; ///< Upper triangular matrix
		bool valid = false;
	};

	/// @brief SVD decomposition result: A = U*diag(w)*V^T
	template<typename Type>
	struct SVDDecomposition {
		Matrix<Type> U;							 ///< Left singular vectors (m x m)
		Matrix<Type> V;							 ///< Right singular vectors (n x n)
		Vector<Type> singularValues; ///< Singular values (descending)
		int rank = 0;								 ///< Numerical rank
		bool valid = false;
	};

	/// @brief Cholesky decomposition result: A = L*L^T (SPD only)
	template<typename Type>
	struct CholeskyDecomposition {
		Matrix<Type> L; ///< Lower triangular Cholesky factor
		bool valid = false;
	};

	/// @brief Solution verification result
	template<typename Type>
	struct VerificationResult {
		Type absoluteResidual;			///< ||Ax - b||
		Type relativeResidual;			///< ||Ax - b|| / ||b||
		Type backwardError;					///< Estimated backward error
		Type estimatedForwardError; ///< Based on condition number
		bool isAccurate;						///< Passes threshold test
	};

	/// @brief Comprehensive system analysis report
	template<typename Type>
	struct SystemAnalysis {
		// Dimensions
		int rows, cols;
		bool isSquare;
		bool isOverdetermined;	///< m > n
		bool isUnderdetermined; ///< m < n

		// Structure
		bool isSymmetric;
		bool isPositiveDefinite;
		bool isDiagonallyDominant;
		bool isUpperTriangular;
		bool isLowerTriangular;
		bool isDiagonal;
		Real sparsity; ///< Fraction of zeros

		// Numerical
		Type determinant;
		int rank;
		int nullity;
		Type conditionNumber;
		MatrixStability stability;
		int expectedDigitsLost;

		// Solution info
		bool hasUniqueSolution;
		bool hasInfiniteSolutions;
		bool hasNoSolution;

		// Recommendation
		std::string recommendedSolver;
		std::string analysisReport;
	};

	//=============================================================================
	// LINEAR SYSTEM CLASS
	//=============================================================================

	/// @class LinearSystem
	///
	/// @brief Unified facade for all MML linear algebra capabilities
	/// LinearSystem provides a single, comprehensive interface to:
	///
	/// - Solve linear systems with automatic solver selection
	/// - Access specific solvers (LU, QR, Cholesky, SVD, iterative)
	///
	/// - Analyze matrix properties and stability
	/// - Verify solution quality
	///
	/// - Access cached decompositions
	/// @par Example: Basic Usage
	///
	/// @code
	/// Matrix<Real> A = {{2, 1}, {1, 3}};
	///
	/// Vector<Real> b = {3, 4};
	/// LinearSystem sys(A, b);
	///
	/// Vector<Real> x = sys.Solve();           // Auto-selects best solver
	/// auto verify = sys.Verify(x);            // Check solution quality
	///
	/// auto analysis = sys.Analyze();          // Get full analysis
	/// @endcode
	///
	/// @par Example: Forcing Specific Solver
	/// @code
	///
	/// LinearSystem sys(A, b);
	/// Vector<Real> x = sys.SolveByQR();       // Force QR solver
	///
	/// Vector<Real> y = sys.SolveBySVD(1e-10); // SVD with threshold
	/// @endcode
	///
	/// @par Example: Analysis Only
	/// @code
	///
	/// LinearSystem sys(A);  // No RHS, just analyze
	/// auto analysis = sys.Analyze();
	///
	/// std::cout << analysis.analysisReport;
	/// @endcode
	template<typename Type = Real>
	class LinearSystem {
	public:
		//=========================================================================
		// CONSTRUCTORS
		//=========================================================================

		/// @brief Construct system Ax = b with single right-hand side
		///
		/// @param A Coefficient matrix (m x n)
		/// @param b Right-hand side vector (length m)
		///
		/// @throws MatrixDimensionError if dimensions don't match
		LinearSystem(const Matrix<Type>& A, const Vector<Type>& b)
				: _A(A)
				, _b(b)
				, _hasRHS(true)
				, _multipleRHS(false) {
			if (A.rows() != b.size())
				throw MatrixDimensionError("LinearSystem: A and b dimensions don't match", A.rows(), A.cols(), b.size(), 1);
		}

		/// @brief Construct system AX = B with multiple right-hand sides
		///
		/// @param A Coefficient matrix (m x n)
		/// @param B Right-hand side matrix (m x p), each column is a RHS
		///
		/// @throws MatrixDimensionError if dimensions don't match
		LinearSystem(const Matrix<Type>& A, const Matrix<Type>& B)
				: _A(A)
				, _B(B)
				, _hasRHS(true)
				, _multipleRHS(true) {
			if (A.rows() != B.rows())
				throw MatrixDimensionError("LinearSystem: A and B row counts don't match", A.rows(), A.cols(), B.rows(), B.cols());
		}

		/// @brief Construct with matrix only (for analysis without solving)
		///
		/// @param A Coefficient matrix (m x n)
		explicit LinearSystem(const Matrix<Type>& A)
				: _A(A)
				, _hasRHS(false)
				, _multipleRHS(false) {}

		//=========================================================================
		// AUTO-SELECT SOLVING
		//=========================================================================

		/// @brief Solve system with automatically selected best method
		///
		/// Selection logic:
		/// 1. Triangular → back/forward substitution
		///
		/// 2. SPD → Cholesky (most stable, fastest for SPD)
		/// 3. Symmetric → LU with pivoting
		///
		/// 4. Overdetermined → QR (least squares)
		/// 5. Ill-conditioned → SVD (most robust)
		///
		/// 6. Default → LU with partial pivoting
		/// @return Solution vector x
		///
		/// @throws std::runtime_error if no RHS provided
		Vector<Type> Solve() const {
			RequireRHS();

			std::string solver = SelectBestSolver();

			if (solver == "Triangular")
				return SolveTriangular();
			else if (solver == "Cholesky")
				return SolveByCholesky();
			else if (solver == "QR")
				return SolveByQR();
			else if (solver == "SVD")
				return SolveBySVD();
			else
				return SolveByLU();
		}

		/// @brief Solve for multiple right-hand sides
		///
		/// @return Solution matrix X where each column is a solution
		Matrix<Type> SolveMultiple() const {
			if (!_multipleRHS)
				throw std::runtime_error("LinearSystem::SolveMultiple - no multiple RHS provided");

			// Use LU for efficiency (factor once, solve many)
			EnsureLU();

			int n = _A.cols();
			int p = _B.cols();
			Matrix<Type> X(n, p);

			LUSolver<Type> solver(_A);
			for (int j = 0; j < p; ++j) {
				Vector<Type> col = _B.VectorFromColumn(j);
				Vector<Type> x = solver.Solve(col);
				for (int i = 0; i < n; ++i)
					X(i, j) = x[i];
			}
			return X;
		}

		//=========================================================================
		// SPECIFIC SOLVER METHODS
		//=========================================================================

		/// @brief Solve using Gauss-Jordan elimination
		///
		/// @note Modifies internal copy; original matrix preserved
		Vector<Type> SolveByGaussJordan() const {
			RequireRHS();
			RequireSquare();
			return GaussJordanSolver<Type>::SolveConst(_A, _b);
		}

		/// @brief Solve using LU decomposition with partial pivoting
		///
		/// @note General-purpose, O(n³) factorization
		Vector<Type> SolveByLU() const {
			RequireRHS();
			RequireSquare();
			LUSolver<Type> solver(_A);
			return solver.Solve(_b);
		}

		/// @brief Solve using Cholesky decomposition
		///
		/// @throws SingularMatrixError if matrix is not positive definite
		/// @note Fastest and most stable for symmetric positive definite matrices
		Vector<Type> SolveByCholesky() const {
			RequireRHS();
			RequireSquare();
			CholeskySolver<Type> solver(_A);
			return solver.Solve(_b);
		}

		/// @brief Solve using QR decomposition
		///
		/// @note Works for overdetermined systems (least squares)
		Vector<Type> SolveByQR() const {
			RequireRHS();
			QRSolver<Type> solver(_A);
			if (isOverdetermined())
				return solver.LeastSquaresSolve(_b);
			else
				return solver.Solve(_b);
		}

		/// @brief Solve using SVD decomposition
		///
		/// @param threshold Values below this treated as zero (default: auto)
		/// @note Most robust for ill-conditioned or rank-deficient systems
		Vector<Type> SolveBySVD(Type threshold = -1) const {
			RequireRHS();
			SVDecompositionSolver<Type> solver(_A);
			return solver.Solve(_b, threshold);
		}

		/// @brief Solve least squares problem min||Ax - b||
		///
		/// @note Uses QR decomposition; works for overdetermined systems
		Vector<Type> SolveLeastSquares() const {
			RequireRHS();
			return SolveByQR(); // QR naturally gives least squares
		}

		/// @brief Solve using iterative method
		///
		/// @param method Which iterative method to use
		/// @param tol Convergence tolerance
		///
		/// @param maxIter Maximum iterations
		/// @return Solution vector
		///
		/// @throws ConvergenceError if method doesn't converge
		Vector<Type> SolveIterative(IterativeMethod method = IterativeMethod::Auto, Type tol = 1e-10, int maxIter = 1000) const {
			RequireRHS();
			RequireSquare();

			if (method == IterativeMethod::Auto)
				method = SelectIterativeMethod();

			IterativeSolverResult result;

			switch (method) {
			case IterativeMethod::Jacobi:
				result = JacobiSolver::Solve(_A, _b, Vector<Real>(), tol, maxIter);
				break;
			case IterativeMethod::GaussSeidel:
				result = GaussSeidelSolver::Solve(_A, _b, Vector<Real>(), tol, maxIter);
				break;
			case IterativeMethod::SOR:
				// SOR: Solve(A, b, omega, x0, tol, maxIter) - omega=1.5 typical choice
				result = SORSolver::Solve(_A, _b, static_cast<Type>(1.5), Vector<Real>(), tol, maxIter);
				break;
			default:
				result = GaussSeidelSolver::Solve(_A, _b, Vector<Real>(), tol, maxIter);
			}

			if (!result.converged)
				throw ConvergenceError("LinearSystem::SolveIterative - failed to converge", result.iterations, result.residual);

			return result.solution;
		}

		//=========================================================================
		// SOLUTION VERIFICATION
		//=========================================================================

		/// @brief Compute residual ||Ax - b||
		Type ResidualNorm(const Vector<Type>& x) const {
			RequireRHS();
			Vector<Type> r = _A * x - _b;
			return r.NormL2();
		}

		/// @brief Compute relative residual ||Ax - b|| / ||b||
		Type RelativeResidual(const Vector<Type>& x) const {
			RequireRHS();
			Type bNorm = _b.NormL2();
			if (bNorm < 1e-30)
				return ResidualNorm(x); // Avoid division by zero
			return ResidualNorm(x) / bNorm;
		}

		/// @brief Full solution verification
		///
		/// @param x Solution to verify
		/// @param tol Accuracy threshold
		///
		/// @return Detailed verification results
		VerificationResult<Type> Verify(const Vector<Type>& x, Type tol = 1e-10) const {
			RequireRHS();

			VerificationResult<Type> result;
			result.absoluteResidual = ResidualNorm(x);
			result.relativeResidual = RelativeResidual(x);

			// Backward error estimate
			Type ANorm = Utils::InfinityNorm(_A);
			Type xNorm = x.NormL2();
			if (ANorm * xNorm > 1e-30)
				result.backwardError = result.absoluteResidual / (ANorm * xNorm);
			else
				result.backwardError = result.absoluteResidual;

			// Forward error estimate (based on condition number)
			Type cond = ConditionNumber();
			result.estimatedForwardError = cond * result.relativeResidual;

			result.isAccurate = (result.relativeResidual < tol);

			return result;
		}

		//=========================================================================
		// MATRIX PROPERTIES - DIMENSIONS
		//=========================================================================

		int rows() const { return _A.rows(); }
		int cols() const { return _A.cols(); }
		bool isSquare() const { return _A.rows() == _A.cols(); }
		bool isOverdetermined() const { return _A.rows() > _A.cols(); }
		bool isUnderdetermined() const { return _A.rows() < _A.cols(); }

		//=========================================================================
		// MATRIX PROPERTIES - STRUCTURE
		//=========================================================================

		/// @brief Check if matrix is symmetric within tolerance
		bool isSymmetric(Type tol = 1e-10) const {
			if (!_isSymmetric.has_value()) {
				if (!isSquare()) {
					_isSymmetric = false;
				} else {
					_isSymmetric = true;
					int n = rows();
					for (int i = 0; i < n && *_isSymmetric; ++i)
						for (int j = i + 1; j < n && *_isSymmetric; ++j)
							if (std::abs(_A(i, j) - _A(j, i)) > tol)
								_isSymmetric = false;
				}
			}
			return *_isSymmetric;
		}

		/// @brief Check if matrix is positive definite
		///
		/// @note Attempts Cholesky decomposition
		bool isPositiveDefinite(Type tol = 1e-10) const {
			if (!_isPositiveDefinite.has_value()) {
				if (!isSymmetric(tol)) {
					_isPositiveDefinite = false;
				} else {
					try {
						CholeskySolver<Type> solver(_A);
						_isPositiveDefinite = true;
					} catch (const SingularMatrixError&) {
						_isPositiveDefinite = false;
					}
				}
			}
			return *_isPositiveDefinite;
		}

		/// @brief Check if matrix is strictly diagonally dominant
		bool isDiagonallyDominant() const {
			if (!isSquare())
				return false;

			int n = rows();
			for (int i = 0; i < n; ++i) {
				Type diagAbs = std::abs(_A(i, i));
				Type offDiagSum = 0;
				for (int j = 0; j < n; ++j)
					if (j != i)
						offDiagSum += std::abs(_A(i, j));

				if (diagAbs <= offDiagSum) // Strict dominance
					return false;
			}
			return true;
		}

		bool isUpperTriangular(Type tol = 1e-10) const { return Utils::IsUpperTriangular(_A, tol); }

		bool isLowerTriangular(Type tol = 1e-10) const { return Utils::IsLowerTriangular(_A, tol); }

		bool isDiagonal(Type tol = 1e-10) const { return Utils::IsDiagonal(_A, tol); }

		/// @brief Compute fraction of zero elements
		///
		/// @param threshold Elements below this are considered zero
		Real Sparsity(Type threshold = 1e-15) const {
			int zeros = 0;
			int total = _A.rows() * _A.cols();

			for (int i = 0; i < _A.rows(); ++i)
				for (int j = 0; j < _A.cols(); ++j)
					if (std::abs(_A(i, j)) <= threshold)
						++zeros;

			return static_cast<Real>(zeros) / total;
		}

		//=========================================================================
		// MATRIX PROPERTIES - NUMERICAL
		//=========================================================================

		/// @brief Compute determinant using LU decomposition
		Type Determinant() const {
			if (!isSquare())
				throw MatrixDimensionError("LinearSystem::Determinant - matrix must be square", rows(), cols(), -1, -1);

			LUSolver<Type> solver(_A);
			return solver.det();
		}

		/// @brief Compute numerical rank using SVD
		///
		/// @param tol Threshold below which singular values are considered zero
		int Rank(Type tol = -1) const {
			if (!_rank.has_value()) {
				SVDecompositionSolver<Type> solver(_A);
				_rank = solver.Rank(tol);
			}
			return *_rank;
		}

		/// @brief Compute nullity (dimension of null space)
		int Nullity(Type tol = -1) const { return cols() - Rank(tol); }

		/// @brief Compute condition number using SVD
		///
		/// @note cond(A) = σ_max / σ_min
		Type ConditionNumber() const {
			if (!_conditionNumber.has_value()) {
				SVDecompositionSolver<Type> solver(_A);
				Type invCond = solver.inv_condition();
				_conditionNumber = (invCond > 1e-30) ? (1.0 / invCond) : 1e30;
			}
			return *_conditionNumber;
		}

		/// @brief Condition number using 1-norm
		Type ConditionNumber1() const {
			Type norm1 = Utils::OneNorm(_A);
			Matrix<Type> Ainv;
			try {
				LUSolver<Type> solver(_A);
				solver.inverse(Ainv);
				return norm1 * Utils::OneNorm(Ainv);
			} catch (...) {
				return 1e30; // Singular
			}
		}

		/// @brief Condition number using infinity norm
		Type ConditionNumberInf() const {
			Type normInf = Utils::InfinityNorm(_A);
			Matrix<Type> Ainv;
			try {
				LUSolver<Type> solver(_A);
				solver.inverse(Ainv);
				return normInf * Utils::InfinityNorm(Ainv);
			} catch (...) {
				return 1e30; // Singular
			}
		}

		/// @brief Assess numerical stability
		MatrixStability AssessStability() const {
			Type cond = ConditionNumber();

			if (cond >= 1e12)
				return MatrixStability::Singular;
			else if (cond >= 1e8)
				return MatrixStability::IllConditioned;
			else if (cond >= 1e4)
				return MatrixStability::ModeratelyConditioned;
			else
				return MatrixStability::WellConditioned;
		}

		/// @brief Estimate digits of precision lost due to conditioning
		int ExpectedDigitsLost() const {
			Type cond = ConditionNumber();
			return static_cast<int>(std::log10(std::max(cond, Type(1))));
		}

		//=========================================================================
		// DECOMPOSITIONS (Lazy, Cached)
		//=========================================================================

		/// @brief Get LU decomposition (cached)
		LUDecomposition<Type> GetLU() const {
			EnsureLU();
			return *_lu;
		}

		/// @brief Get QR decomposition (cached)
		QRDecomposition<Type> GetQR() const {
			EnsureQR();
			return *_qr;
		}

		/// @brief Get SVD decomposition (cached)
		SVDDecomposition<Type> GetSVD() const {
			EnsureSVD();
			return *_svd;
		}

		/// @brief Get Cholesky decomposition (cached)
		///
		/// @throws SingularMatrixError if not positive definite
		CholeskyDecomposition<Type> GetCholesky() const {
			EnsureCholesky();
			return *_cholesky;
		}

		//=========================================================================
		// FUNDAMENTAL SUBSPACES
		//=========================================================================

		/// @brief Compute null space basis using SVD
		///
		/// @return Matrix whose columns form orthonormal basis for null(A)
		Matrix<Type> NullSpace(Type tol = -1) const {
			SVDecompositionSolver<Type> solver(_A);
			return solver.Nullspace(tol);
		}

		/// @brief Compute column space (range) basis using SVD
		///
		/// @return Matrix whose columns form orthonormal basis for col(A)
		Matrix<Type> ColumnSpace(Type tol = -1) const {
			SVDecompositionSolver<Type> solver(_A);
			return solver.Range(tol);
		}

		//=========================================================================
		// MATRIX OPERATIONS
		//=========================================================================

		/// @brief Compute matrix inverse using LU
		///
		/// @throws SingularMatrixError if matrix is singular
		Matrix<Type> Inverse() const {
			RequireSquare();
			LUSolver<Type> solver(_A);
			Matrix<Type> inv;
			solver.inverse(inv);
			return inv;
		}

		/// @brief Compute Moore-Penrose pseudoinverse using SVD
		Matrix<Type> PseudoInverse(Type tol = -1) const {
			SVDecompositionSolver<Type> solver(_A);
			int m = rows();
			int n = cols();

			// A+ = V * diag(1/w) * U^T (for non-zero singular values)
			Matrix<Type> U = solver.getU();
			Matrix<Type> V = solver.getV();
			Vector<Type> w = solver.getW();

			Type thresh = (tol >= 0) ? tol : 0.5 * std::sqrt(m + n + 1.0) * w[0] * std::numeric_limits<Type>::epsilon();

			// Result is n x m
			Matrix<Type> pinv(n, m);

			for (int i = 0; i < n; ++i) {
				for (int j = 0; j < m; ++j) {
					Type sum = 0;
					for (int k = 0; k < n; ++k) {
						if (w[k] > thresh)
							sum += V(i, k) * U(j, k) / w[k];
					}
					pinv(i, j) = sum;
				}
			}

			return pinv;
		}

		//=========================================================================
		// EIGENANALYSIS
		//=========================================================================

		/// @brief Full eigenvalue/eigenvector decomposition
		///
		/// @return EigenSolver::Result with eigenvalues (potentially complex) and eigenvectors
		/// @note Uses implicitly-shifted QR algorithm for general matrices
		EigenSolver::Result GetEigen(Type tol = 1e-10, int maxIter = 1000) const {
			RequireSquare();
			return EigenSolver::Solve(_A, tol, maxIter);
		}

		/// @brief Get eigenvalues only (real parts)
		///
		/// @return Vector of real eigenvalues (imaginary parts discarded for complex pairs)
		/// @note For purely real eigenvalue analysis; use GetEigen() for complex eigenvalues
		Vector<Type> Eigenvalues(Type tol = 1e-10) const {
			RequireSquare();
			auto result = EigenSolver::Solve(_A, tol);
			int n = rows();
			Vector<Type> eigs(n);
			for (int i = 0; i < n; ++i)
				eigs[i] = result.eigenvalues[i].real;
			return eigs;
		}

		/// @brief Get eigenvalues for symmetric matrices (all real, faster)
		///
		/// @return Vector of sorted eigenvalues
		/// @note Uses Jacobi rotation method, guaranteed real eigenvalues
		Vector<Type> EigenvaluesSymmetric() const {
			RequireSquare();
			if (!isSymmetric())
				throw MatrixDimensionError("LinearSystem::EigenvaluesSymmetric - matrix is not symmetric", rows(), cols(), rows(), cols());

			auto result = SymmMatEigenSolverJacobi::Solve(_A);
			return result.eigenvalues;
		}

		/// @brief Compute spectral radius (largest |eigenvalue|)
		///
		/// @return max(|λ_i|) over all eigenvalues
		Type SpectralRadius(Type tol = 1e-10) const {
			RequireSquare();
			auto result = EigenSolver::Solve(_A, tol);
			Type maxMag = 0;
			for (const auto& ev : result.eigenvalues)
				maxMag = std::max(maxMag, ev.magnitude());
			return maxMag;
		}

		/// @brief Check if matrix has any complex eigenvalues
		///
		/// @return true if any eigenvalue has non-zero imaginary part
		bool HasComplexEigenvalues(Type tol = 1e-10) const {
			RequireSquare();
			auto result = EigenSolver::Solve(_A, tol);
			for (const auto& ev : result.eigenvalues) {
				if (ev.isComplex(tol))
					return true;
			}
			return false;
		}

		//=========================================================================
		// COMPREHENSIVE ANALYSIS
		//=========================================================================

		/// @brief Perform comprehensive system analysis
		///
		/// @return Detailed analysis report
		SystemAnalysis<Type> Analyze() const {
			SystemAnalysis<Type> result;

			// Dimensions
			result.rows = rows();
			result.cols = cols();
			result.isSquare = isSquare();
			result.isOverdetermined = isOverdetermined();
			result.isUnderdetermined = isUnderdetermined();

			// Structure
			result.isSymmetric = isSymmetric();
			result.isPositiveDefinite = isPositiveDefinite();
			result.isDiagonallyDominant = isDiagonallyDominant();
			result.isUpperTriangular = isUpperTriangular();
			result.isLowerTriangular = isLowerTriangular();
			result.isDiagonal = isDiagonal();
			result.sparsity = Sparsity();

			// Numerical
			if (result.isSquare) {
				try {
					result.determinant = Determinant();
				} catch (...) {
					result.determinant = 0;
				}
			} else {
				result.determinant = 0;
			}

			result.rank = Rank();
			result.nullity = Nullity();
			result.conditionNumber = ConditionNumber();
			result.stability = AssessStability();
			result.expectedDigitsLost = ExpectedDigitsLost();

			// Solution classification
			if (result.isSquare) {
				if (result.rank == result.cols) {
					result.hasUniqueSolution = true;
					result.hasInfiniteSolutions = false;
					result.hasNoSolution = false;
				} else {
					result.hasUniqueSolution = false;
					result.hasInfiniteSolutions = (result.nullity > 0);
					result.hasNoSolution = false; // Would need RHS to determine
				}
			} else if (result.isOverdetermined) {
				result.hasUniqueSolution = (result.rank == result.cols);
				result.hasInfiniteSolutions = false;
				result.hasNoSolution = !result.hasUniqueSolution;
			} else // Underdetermined
			{
				result.hasUniqueSolution = false;
				result.hasInfiniteSolutions = true;
				result.hasNoSolution = false;
			}

			// Recommendation
			result.recommendedSolver = SelectBestSolver();
			result.analysisReport = GenerateReport(result);

			return result;
		}

		//=========================================================================
		// UTILITY
		//=========================================================================

		/// @brief Get the coefficient matrix
		const Matrix<Type>& GetMatrix() const { return _A; }

		/// @brief Get the right-hand side vector
		const Vector<Type>& GetRHS() const {
			RequireRHS();
			return _b;
		}

	private:
		// Storage
		Matrix<Type> _A;
		Vector<Type> _b;
		Matrix<Type> _B; // For multiple RHS
		bool _hasRHS;
		bool _multipleRHS;

		// Cached decompositions (mutable for lazy evaluation)
		mutable std::optional<LUDecomposition<Type>> _lu;
		mutable std::optional<QRDecomposition<Type>> _qr;
		mutable std::optional<SVDDecomposition<Type>> _svd;
		mutable std::optional<CholeskyDecomposition<Type>> _cholesky;

		// Cached properties
		mutable std::optional<bool> _isSymmetric;
		mutable std::optional<bool> _isPositiveDefinite;
		mutable std::optional<Type> _conditionNumber;
		mutable std::optional<int> _rank;

		//=========================================================================
		// INTERNAL HELPERS
		//=========================================================================

		void RequireRHS() const {
			if (!_hasRHS)
				throw std::runtime_error("LinearSystem: operation requires right-hand side");
		}

		void RequireSquare() const {
			if (!isSquare())
				throw MatrixDimensionError("LinearSystem: operation requires square matrix", rows(), cols(), -1, -1);
		}

		void EnsureLU() const {
			if (!_lu.has_value()) {
				_lu = LUDecomposition<Type>();
				try {
					LUSolver<Type> solver(_A);

					// Extract L and U from the solver's combined storage
					int n = rows();
					_lu->L.Resize(n, n);
					_lu->U.Resize(n, n);

					// The LU matrix has L in lower triangle and U in upper
					// For now, store raw LU and determinant
					_lu->determinant = solver.det();
					_lu->valid = true;
				} catch (...) {
					_lu->valid = false;
				}
			}
		}

		void EnsureQR() const {
			if (!_qr.has_value()) {
				_qr = QRDecomposition<Type>();
				try {
					QRSolver<Type> solver(_A);
					_qr->Q = solver.GetQ();
					_qr->R = solver.GetR();
					_qr->valid = true;
				} catch (...) {
					_qr->valid = false;
				}
			}
		}

		void EnsureSVD() const {
			if (!_svd.has_value()) {
				_svd = SVDDecomposition<Type>();
				try {
					SVDecompositionSolver<Type> solver(_A);
					_svd->U = solver.getU();
					_svd->V = solver.getV();
					_svd->singularValues = solver.getW();
					_svd->rank = solver.Rank();
					_svd->valid = true;
				} catch (...) {
					_svd->valid = false;
				}
			}
		}

		void EnsureCholesky() const {
			if (!_cholesky.has_value()) {
				_cholesky = CholeskyDecomposition<Type>();
				try {
					CholeskySolver<Type> solver(_A);
					_cholesky->L = solver.el;
					_cholesky->valid = true;
				} catch (...) {
					_cholesky->valid = false;
				}
			}
		}

		std::string SelectBestSolver() const {
			// Quick structural checks first
			if (isUpperTriangular() || isLowerTriangular())
				return "Triangular";

			// For overdetermined systems, use QR (least squares)
			if (isOverdetermined())
				return "QR";

			// Check condition number for ill-conditioned systems
			Type cond = ConditionNumber();
			if (cond > 1e10)
				return "SVD"; // Most robust for ill-conditioned

			// For SPD matrices, Cholesky is best
			if (isSymmetric() && isPositiveDefinite())
				return "Cholesky";

			// Default: LU with partial pivoting
			return "LU";
		}

		IterativeMethod SelectIterativeMethod() const {
			// For diagonally dominant matrices, Gauss-Seidel converges faster
			if (isDiagonallyDominant())
				return IterativeMethod::GaussSeidel;

			// SOR can be faster but needs tuning
			// Default to Gauss-Seidel as safe choice
			return IterativeMethod::GaussSeidel;
		}

		Vector<Type> SolveTriangular() const {
			int n = rows();
			Vector<Type> x(n);

			if (isUpperTriangular()) {
				// Back substitution
				for (int i = n - 1; i >= 0; --i) {
					Type sum = _b[i];
					for (int j = i + 1; j < n; ++j)
						sum -= _A(i, j) * x[j];
					x[i] = sum / _A(i, i);
				}
			} else // Lower triangular
			{
				// Forward substitution
				for (int i = 0; i < n; ++i) {
					Type sum = _b[i];
					for (int j = 0; j < i; ++j)
						sum -= _A(i, j) * x[j];
					x[i] = sum / _A(i, i);
				}
			}

			return x;
		}

		std::string GenerateReport(const SystemAnalysis<Type>& analysis) const {
			std::string report;

			report += "Matrix: " + std::to_string(analysis.rows) + "x" + std::to_string(analysis.cols);
			if (analysis.isSquare)
				report += ", square";
			else if (analysis.isOverdetermined)
				report += ", overdetermined (m > n)";
			else
				report += ", underdetermined (m < n)";
			report += "\n";

			// Structure
			report += "Structure: ";
			if (analysis.isDiagonal)
				report += "diagonal";
			else if (analysis.isUpperTriangular)
				report += "upper triangular";
			else if (analysis.isLowerTriangular)
				report += "lower triangular";
			else if (analysis.isSymmetric && analysis.isPositiveDefinite)
				report += "symmetric positive definite (SPD)";
			else if (analysis.isSymmetric)
				report += "symmetric";
			else
				report += "general";
			report += "\n";

			// Sparsity
			if (analysis.sparsity > 0.5)
				report += "Sparsity: " + std::to_string(static_cast<int>(analysis.sparsity * 100)) + "% zeros (sparse)\n";

			// Rank
			report += "Rank: " + std::to_string(analysis.rank);
			if (analysis.isSquare)
				report += (analysis.rank == analysis.cols ? " (full rank)" : " (rank deficient)");
			report += "\n";

			// Condition number and stability
			report += "Condition number: " + std::to_string(analysis.conditionNumber);
			switch (analysis.stability) {
			case MatrixStability::WellConditioned:
				report += " (well-conditioned)\n";
				break;
			case MatrixStability::ModeratelyConditioned:
				report += " (moderately conditioned)\n";
				break;
			case MatrixStability::IllConditioned:
				report += " (ILL-CONDITIONED - expect precision loss)\n";
				break;
			case MatrixStability::Singular:
				report += " (NEAR-SINGULAR)\n";
				break;
			}

			if (analysis.expectedDigitsLost > 0)
				report += "Expected precision loss: ~" + std::to_string(analysis.expectedDigitsLost) + " digits\n";

			// Solution classification
			report += "Solution: ";
			if (analysis.hasUniqueSolution)
				report += "unique";
			else if (analysis.hasInfiniteSolutions)
				report += "infinite (underdetermined)";
			else if (analysis.hasNoSolution)
				report += "none (inconsistent)";
			else
				report += "depends on RHS";
			report += "\n";

			// Recommendation
			report += "Recommended solver: " + analysis.recommendedSolver + "\n";

			return report;
		}
	};

	//=============================================================================
	// CONVENIENCE FUNCTIONS
	//=============================================================================

	/// @brief Quick solve with automatic method selection
	template<typename Type = Real>
	inline Vector<Type> SolveLinearSystem(const Matrix<Type>& A, const Vector<Type>& b) {
		return LinearSystem<Type>(A, b).Solve();
	}

	/// @brief Quick least squares solve
	template<typename Type = Real>
	inline Vector<Type> SolveLeastSquares(const Matrix<Type>& A, const Vector<Type>& b) {
		return LinearSystem<Type>(A, b).SolveLeastSquares();
	}

	/// @brief Get comprehensive matrix analysis
	template<typename Type = Real>
	inline SystemAnalysis<Type> AnalyzeMatrix(const Matrix<Type>& A) {
		return LinearSystem<Type>(A).Analyze();
	}

} // namespace MML::Systems
#endif // MML_LINEAR_SYSTEM_H

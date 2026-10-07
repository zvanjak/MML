///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        LinearSystem.h                                                      ///
///  Description: RHS-aware linear-system solving facade                              ///
///               Smart solver selection, analyzer delegation, and rich diagnostics   ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_LINEAR_SYSTEM_H
#define MML_LINEAR_SYSTEM_H

#include <mml/MMLBase.h>
#include <mml/base/Vector/Vector.h>
#include <mml/base/Matrix/Matrix.h>
#include <mml/core/LinAlgEqSolvers.h>
#include <mml/algorithms/MatrixAlg.h>
#include <mml/algorithms/Analyzers/MatrixAnalyzer.h>

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

	//=============================================================================
	// RESULT STRUCTURES
	//=============================================================================

	enum class SolutionStatus { NotAnalyzed, Inconsistent, Unique, Infinite };
	enum class MultipleRHSAggregate { AllSame, Mixed };
	enum class LinearSolverRecommendation { Triangular, Cholesky, QR, SVD, LU };
	using MatrixStability = MatrixAlg::MatrixStability;

	/// @brief Solution verification result
	template<MMLReal Type>
	struct VerificationResult {
		MatrixAlg::MatrixMagnitude<Type> absoluteResidual{};
		MatrixAlg::MatrixMagnitude<Type> relativeResidual{};
		MatrixAlg::MatrixMagnitude<Type> backwardError{};
		MatrixAlg::MatrixMagnitude<Type> estimatedForwardError{};
		bool isAccurate = false;
	};

	/// @brief Comprehensive system analysis report
	template<MMLReal Type>
	struct SystemAnalysis {
		MatrixAlg::MatrixAnalysis<Type> matrix;
		std::vector<SolutionStatus> solutionStatuses;
		std::optional<MultipleRHSAggregate> multipleRHSAggregate;
		std::optional<LinearSolverRecommendation> recommendedSolver;
		std::string report;
	};

	//=============================================================================
	// LINEAR SYSTEM CLASS
	//=============================================================================

	/// @class LinearSystem
	///
	/// @brief RHS-aware facade for linear-system solving, analysis, and diagnostics
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
	/// std::cout << analysis.report;
	/// @endcode
	template<MMLReal Type = Real>
		requires std::same_as<std::remove_cvref_t<Type>, Real>
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
				, _matrixAnalyzer(A)
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
				, _matrixAnalyzer(A)
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
				, _matrixAnalyzer(A)
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

			switch (SelectBestSolver()) {
			case LinearSolverRecommendation::Triangular: return SolveTriangular();
			case LinearSolverRecommendation::Cholesky:   return SolveByCholesky();
			case LinearSolverRecommendation::QR:         return SolveByQR();
			case LinearSolverRecommendation::SVD:        return SolveBySVD();
			case LinearSolverRecommendation::LU:         return SolveByLU();
			}
			throw InvalidStateError("LinearSystem::Solve - invalid solver recommendation");
		}

		/// @brief Solve for multiple right-hand sides
		///
		/// @return Solution matrix X where each column is a solution
		Matrix<Type> SolveMultiple() const {
			if (!_multipleRHS)
				throw InvalidStateError("LinearSystem::SolveMultiple - no multiple RHS provided");

			// Use cached LU factorization (factor once, solve many)
			LUSolver<Type>& solver = EnsureLUSolver();

			int n = _A.cols();
			int p = _B.cols();
			Matrix<Type> X(n, p);

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
			if (IsTall())
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
		Vector<Type> SolveIterative(IterativeMethod method = IterativeMethod::Auto, Type tol = Precision::DefaultToleranceStrict, int maxIter = 1000) const {
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
			if (bNorm < Precision::DivisionSafetyThreshold)
				return ResidualNorm(x); // Avoid division by zero
			return ResidualNorm(x) / bNorm;
		}

		/// @brief Full solution verification
		///
		/// @param x Solution to verify
		/// @param tol Accuracy threshold
		///
		/// @return Detailed verification results
		VerificationResult<Type> Verify(const Vector<Type>& x, Type tol = Precision::DefaultToleranceStrict) const {
			RequireRHS();

			VerificationResult<Type> result;
			result.absoluteResidual = ResidualNorm(x);
			result.relativeResidual = RelativeResidual(x);

			// Backward error estimate
			Magnitude ANorm = _matrixAnalyzer.InfinityNorm();
			Type xNorm = x.NormL2();
			if (ANorm * xNorm > Precision::DivisionSafetyThreshold)
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
		using Magnitude = MatrixAlg::MatrixMagnitude<Type>;
		using ComparisonTolerance = MatrixAlg::MatrixComparisonTolerance<Magnitude>;
		using Threshold = std::optional<Magnitude>;

		int Rows() const noexcept { return _matrixAnalyzer.Rows(); }
		int Cols() const noexcept { return _matrixAnalyzer.Cols(); }
		bool IsSquare() const noexcept { return _matrixAnalyzer.IsSquare(); }
		bool IsTall() const noexcept { return _matrixAnalyzer.IsTall(); }
		bool IsWide() const noexcept { return _matrixAnalyzer.IsWide(); }

		bool IsSymmetric(ComparisonTolerance tolerance = MatrixAlg::Detail::DefaultSymmetryTolerance<Type>()) const {
			return _matrixAnalyzer.IsSymmetric(tolerance);
		}
		bool IsPositiveDefinite(Magnitude tolerance = PrecisionValues<Magnitude>::EigenSolverConvergenceTolerance) const {
			return _matrixAnalyzer.IsPositiveDefinite(tolerance);
		}
		bool IsDiagonallyDominant(ComparisonTolerance tolerance = MatrixAlg::Detail::DefaultDiagonalTolerance<Type>()) const {
			return _matrixAnalyzer.IsDiagonallyDominant(tolerance);
		}
		bool IsUpperTriangular(ComparisonTolerance tolerance = MatrixAlg::Detail::DefaultDiagonalTolerance<Type>()) const {
			return _matrixAnalyzer.IsUpperTriangular(tolerance);
		}
		bool IsLowerTriangular(ComparisonTolerance tolerance = MatrixAlg::Detail::DefaultDiagonalTolerance<Type>()) const {
			return _matrixAnalyzer.IsLowerTriangular(tolerance);
		}
		bool IsDiagonal(ComparisonTolerance tolerance = MatrixAlg::Detail::DefaultDiagonalTolerance<Type>()) const {
			return _matrixAnalyzer.IsDiagonal(tolerance);
		}

		/// @brief Compute fraction of zero elements
		///
		/// @param threshold Elements below this are considered zero
		Real Sparsity(Type threshold = Precision::EigenSolverZeroThreshold) const { return _matrixAnalyzer.Sparsity(threshold); }

		//=========================================================================
		// MATRIX PROPERTIES - NUMERICAL
		//=========================================================================

		/// @brief Compute determinant using LU decomposition
		Type Determinant() const {
			return _matrixAnalyzer.Determinant();
		}

		/// @brief Compute numerical rank using SVD
		///
		/// @param tol Threshold below which singular values are considered zero
		int Rank(Threshold threshold = std::nullopt) const { return _matrixAnalyzer.Rank(threshold); }

		/// @brief Compute nullity (dimension of null space)
		int Nullity(Threshold threshold = std::nullopt) const { return _matrixAnalyzer.Nullity(threshold); }

		/// @brief Compute condition number using SVD
		///
		/// @note cond(A) = σ_max / σ_min
		Magnitude ConditionNumber(Threshold threshold = std::nullopt) const { return _matrixAnalyzer.ConditionNumber(threshold); }

		/// @brief Condition number using 1-norm
		Magnitude ConditionNumber1(Threshold threshold = std::nullopt) const { return _matrixAnalyzer.ConditionNumber1(threshold); }

		/// @brief Condition number using infinity norm
		Magnitude ConditionNumberInfinity(Threshold threshold = std::nullopt) const {
			return _matrixAnalyzer.ConditionNumberInfinity(threshold);
		}

		/// @brief Assess numerical stability
		MatrixAlg::MatrixStability AssessStability(Threshold threshold = std::nullopt) const {
			return _matrixAnalyzer.AssessStability(threshold);
		}

		/// @brief Estimate digits of precision lost due to conditioning
		std::optional<int> ExpectedDigitsLost(Threshold threshold = std::nullopt) const {
			return _matrixAnalyzer.ExpectedDigitsLost(threshold);
		}

		//=========================================================================
		// DECOMPOSITIONS (Lazy, Cached)
		//=========================================================================

		/// @brief Get LU decomposition (cached)
		const MatrixAlg::LUDecomposition<Type>& LUDecompose() const { return _matrixAnalyzer.LUDecompose(); }

		/// @brief Get QR decomposition (cached)
		const MatrixAlg::QRDecomposition<Type>& QRDecompose() const { return _matrixAnalyzer.QRDecompose(); }

		/// @brief Get SVD decomposition (cached)
		const MatrixAlg::SVDDecomposition<Type>& SVDDecompose(Threshold threshold = std::nullopt) const {
			return _matrixAnalyzer.SVDDecompose(threshold);
		}

		/// @brief Get Cholesky decomposition (cached)
		///
		/// @throws SingularMatrixError if not positive definite
		const MatrixAlg::CholeskyDecomposition<Type>& CholeskyDecompose() const { return _matrixAnalyzer.CholeskyDecompose(); }

		//=========================================================================
		// FUNDAMENTAL SUBSPACES
		//=========================================================================

		/// @brief Compute null space basis using SVD
		///
		/// @return Matrix whose columns form orthonormal basis for null(A)
		Matrix<Type> NullSpace(Threshold threshold = std::nullopt) const { return _matrixAnalyzer.NullSpace(threshold); }

		/// @brief Compute column space (range) basis using SVD
		///
		/// @return Matrix whose columns form orthonormal basis for col(A)
		Matrix<Type> ColumnSpace(Threshold threshold = std::nullopt) const { return _matrixAnalyzer.ColumnSpace(threshold); }
		Matrix<Type> RowSpace(Threshold threshold = std::nullopt) const { return _matrixAnalyzer.RowSpace(threshold); }
		Matrix<Type> LeftNullSpace(Threshold threshold = std::nullopt) const { return _matrixAnalyzer.LeftNullSpace(threshold); }
		const MatrixAlg::FundamentalSubspaces<Type>& FundamentalSubspacesOf(Threshold threshold = std::nullopt) const {
			return _matrixAnalyzer.FundamentalSubspacesOf(threshold);
		}

		//=========================================================================
		// MATRIX OPERATIONS
		//=========================================================================

		/// @brief Compute matrix inverse using LU
		///
		/// @throws SingularMatrixError if matrix is singular
		Matrix<Type> Inverse() const { return _matrixAnalyzer.Inverse(); }

		/// @brief Compute Moore-Penrose pseudoinverse using SVD
		Matrix<Type> PseudoInverse(Threshold threshold = std::nullopt) const { return _matrixAnalyzer.PseudoInverse(threshold); }

		//=========================================================================
		// EIGENANALYSIS
		//=========================================================================

		/// @brief Canonical full eigenvalue/eigenvector decomposition
		///
		/// @return EigensystemResult with explicit complex eigenvalues and eigenvectors
		const MatrixAlg::EigensystemResult<Type>& Eigensystem(
				Type tol = Precision::DefaultToleranceStrict, int maxIter = 1000) const {
			return _matrixAnalyzer.Eigensystem(tol, maxIter);
		}

		/// @brief Get eigenvalues with explicit complex components
		///
		/// @return Vector of complex eigenvalues
		Vector<MatrixAlg::MatrixComplexScalar<Type>> Eigenvalues(Type tol = Precision::DefaultToleranceStrict) const {
			return _matrixAnalyzer.Eigenvalues(tol);
		}

		/// @brief Get eigenvalues for symmetric matrices (all real, faster)
		///
		/// @return Vector of sorted eigenvalues
		/// @note Uses Jacobi rotation method, guaranteed real eigenvalues
		Vector<Type> SymmetricEigenvalues() const {
			return _matrixAnalyzer.SymmetricEigenvalues();
		}

		/// @brief Compute spectral radius (largest |eigenvalue|)
		///
		/// @return max(|λ_i|) over all eigenvalues
		Type SpectralRadius(Type tol = Precision::DefaultToleranceStrict) const {
			return _matrixAnalyzer.SpectralRadius(tol);
		}

		/// @brief Check if matrix has any complex eigenvalues
		///
		/// @return true if any eigenvalue has non-zero imaginary part
		bool HasComplexEigenvalues(Type tol = Precision::DefaultToleranceStrict) const {
			return _matrixAnalyzer.HasComplexEigenvalues(tol);
		}

		//=========================================================================
		// COMPREHENSIVE ANALYSIS
		//=========================================================================

		/// @brief Perform comprehensive system analysis
		///
		/// @return Detailed analysis report
		SystemAnalysis<Type> Analyze(Threshold threshold = std::nullopt) const {
			SystemAnalysis<Type> result;
			const auto& coefficientSVD = _matrixAnalyzer.SVDDecompose(threshold);
			result.matrix = _matrixAnalyzer.Analyze(threshold);
			if (_hasRHS) {
				if (_multipleRHS) {
					for (int col = 0; col < _B.cols(); ++col)
						result.solutionStatuses.push_back(ClassifySolution(_B.VectorFromColumn(col), coefficientSVD.rank, coefficientSVD.threshold));
					if (!result.solutionStatuses.empty()) {
						const bool allSame = std::all_of(result.solutionStatuses.begin() + 1, result.solutionStatuses.end(),
							[&](SolutionStatus status) { return status == result.solutionStatuses.front(); });
						result.multipleRHSAggregate = allSame ? MultipleRHSAggregate::AllSame : MultipleRHSAggregate::Mixed;
					}
				}
				else {
					result.solutionStatuses.push_back(ClassifySolution(_b, coefficientSVD.rank, coefficientSVD.threshold));
				}
				result.recommendedSolver = SelectBestSolver(threshold);
			}
			result.report = GenerateReport(result);

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
		MatrixAnalyzer<Type> _matrixAnalyzer;
		Vector<Type> _b;
		Matrix<Type> _B; // For multiple RHS
		bool _hasRHS;
		bool _multipleRHS;

		// Cached solver used only for repeated RHS solves.
		mutable std::optional<LUSolver<Type>> _luSolver;

		//=========================================================================
		// INTERNAL HELPERS
		//=========================================================================

		void RequireRHS() const {
			if (!_hasRHS)
				throw InvalidStateError("LinearSystem: operation requires right-hand side");
		}

		void RequireSquare() const {
			if (!IsSquare())
				throw MatrixDimensionError("LinearSystem: operation requires square matrix", Rows(), Cols(), -1, -1);
		}

		/// @brief Get (or lazily create) the cached LU solver - factor once, solve many
		LUSolver<Type>& EnsureLUSolver() const {
			if (!_luSolver.has_value())
				_luSolver.emplace(_A);
			return *_luSolver;
		}

		LinearSolverRecommendation SelectBestSolver(Threshold threshold = std::nullopt) const {
			const auto stability = AssessStability(threshold);
			if (Rank(threshold) < Cols() || stability == MatrixAlg::MatrixStability::IllConditioned ||
				stability == MatrixAlg::MatrixStability::Singular)
				return LinearSolverRecommendation::SVD;

			if (IsSquare() && (IsUpperTriangular() || IsLowerTriangular()))
				return LinearSolverRecommendation::Triangular;

			if (IsTall())
				return LinearSolverRecommendation::QR;

			// For SPD matrices, Cholesky is best
			if (IsSymmetric() && IsPositiveDefinite())
				return LinearSolverRecommendation::Cholesky;

			// Default: LU with partial pivoting
			return LinearSolverRecommendation::LU;
		}

		SolutionStatus ClassifySolution(const Vector<Type>& rhs, int coefficientRank, Magnitude threshold) const {
			Matrix<Type> augmented(Rows(), Cols() + 1);
			for (int row = 0; row < Rows(); ++row) {
				for (int col = 0; col < Cols(); ++col)
					augmented(row, col) = _A(row, col);
				augmented(row, Cols()) = rhs[row];
			}
			if (MatrixAlg::Rank(augmented, threshold) > coefficientRank)
				return SolutionStatus::Inconsistent;
			return coefficientRank == Cols() ? SolutionStatus::Unique : SolutionStatus::Infinite;
		}

		IterativeMethod SelectIterativeMethod() const {
			// For diagonally dominant matrices, Gauss-Seidel converges faster
			if (IsDiagonallyDominant())
				return IterativeMethod::GaussSeidel;

			// SOR can be faster but needs tuning
			// Default to Gauss-Seidel as safe choice
			return IterativeMethod::GaussSeidel;
		}

		Vector<Type> SolveTriangular() const {
			int n = Rows();
			Vector<Type> x(n);

			if (IsUpperTriangular()) {
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
			report += analysis.matrix.report;
			for (size_t index = 0; index < analysis.solutionStatuses.size(); ++index) {
				report += "RHS " + std::to_string(index) + ": ";
				switch (analysis.solutionStatuses[index]) {
				case SolutionStatus::NotAnalyzed: report += "not analyzed"; break;
				case SolutionStatus::Inconsistent: report += "inconsistent"; break;
				case SolutionStatus::Unique: report += "unique"; break;
				case SolutionStatus::Infinite: report += "infinite"; break;
				}
				report += '\n';
			}
			if (analysis.recommendedSolver.has_value())
				report += "Recommended solver: " + SolverName(*analysis.recommendedSolver) + "\n";

			return report;
		}

		static std::string SolverName(LinearSolverRecommendation solver) {
			switch (solver) {
			case LinearSolverRecommendation::Triangular: return "Triangular";
			case LinearSolverRecommendation::Cholesky: return "Cholesky";
			case LinearSolverRecommendation::QR: return "QR";
			case LinearSolverRecommendation::SVD: return "SVD";
			case LinearSolverRecommendation::LU: return "LU";
			}
			return "Unknown";
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

} // namespace MML::Systems
#endif // MML_LINEAR_SYSTEM_H

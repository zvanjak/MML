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
#include <mml/MMLNumericValidation.h>
#include <algorithm>
#include <mml/base/Vector/Vector.h>
#include <mml/base/Matrix/Matrix.h>
#include <mml/core/LinAlgEqSolvers.h>
#include <mml/algorithms/MatrixAlg.h>
#include <mml/algorithms/Analyzers/MatrixAnalyzer.h>
#include <mml/tools/analyzers/MatrixAnalysisReports.h>

#include <optional>
#include <sstream>
#include <string>
#include <memory>
#include <cmath>
#include <limits>
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
				: _coefficientMatrix(A)
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
				: _coefficientMatrix(A)
				, _matrixAnalyzer(A)
				, _rhsMatrix(B)
				, _hasRHS(true)
				, _multipleRHS(true) {
			if (A.rows() != B.rows())
				throw MatrixDimensionError("LinearSystem: A and B row counts don't match", A.rows(), A.cols(), B.rows(), B.cols());
		}

		/// @brief Construct with matrix only (for analysis without solving)
		///
		/// @param A Coefficient matrix (m x n)
		explicit LinearSystem(const Matrix<Type>& A)
				: _coefficientMatrix(A)
				, _matrixAnalyzer(A)
				, _hasRHS(false)
				, _multipleRHS(false) {}

		//=========================================================================
		// AUTO-SELECT SOLVING
		//=========================================================================

		/// @brief Solve system with automatically selected best method
		///
		/// Selection logic:
		/// 1. Rank-deficient, singular or ill-conditioned → SVD
		/// 2. Exactly triangular square → back/forward substitution
		/// 3. Overdetermined → QR (least squares)
		/// 4. Exactly symmetric positive definite → Cholesky
		/// 5. Default → LU with partial pivoting
		/// @return Solution vector x
		///
		/// @throws InvalidStateError unless constructed with a vector RHS
		Vector<Type> Solve() const {
			RequireVectorRHS();

			switch (SelectBestSolver()) {
			case LinearSolverRecommendation::Triangular: return SolveTriangular(_b);
			case LinearSolverRecommendation::Cholesky:   return SolveByCholesky();
			case LinearSolverRecommendation::QR:         return SolveByQR();
			case LinearSolverRecommendation::SVD:        return SolveBySVD();
			case LinearSolverRecommendation::LU:         return SolveByLU();
			}
			throw InvalidStateError("LinearSystem::Solve - invalid solver recommendation");
		}

		/// @brief Solve for multiple right-hand sides using the same selection as Solve()
		///
		/// @note Reuses one solver across columns; general square LU is cached across calls.
		/// @throws InvalidStateError unless constructed with a matrix RHS
		/// @return Solution matrix X where each column is a solution
		Matrix<Type> SolveMultiple() const {
			RequireMatrixRHS();
			for (int row = 0; row < _rhsMatrix.rows(); ++row)
				for (int col = 0; col < _rhsMatrix.cols(); ++col)
					ValidateFinite(_rhsMatrix(row, col), "LinearSystem::SolveMultiple RHS");

			const auto solveColumns = [&](auto&& solve) {
				Matrix<Type> solutions(Cols(), _rhsMatrix.cols());
				for (int col = 0; col < _rhsMatrix.cols(); ++col) {
					const auto solution = solve(_rhsMatrix.VectorFromColumn(col));
					for (int row = 0; row < Cols(); ++row) {
						if (!std::isfinite(solution[row]))
							throw MatrixNumericalError("LinearSystem::SolveMultiple - non-finite solution");
						solutions(row, col) = solution[row];
					}
				}
				return solutions;
			};

			switch (SelectBestSolver()) {
			case LinearSolverRecommendation::Triangular:
				return solveColumns([&](const Vector<Type>& rhs) { return SolveTriangular(rhs); });
			case LinearSolverRecommendation::Cholesky: {
				CholeskySolver<Type> solver(_coefficientMatrix);
				return solveColumns([&](const Vector<Type>& rhs) { return solver.Solve(rhs); });
			}
			case LinearSolverRecommendation::QR: {
				QRSolver<Type> solver(_coefficientMatrix);
				return solveColumns([&](const Vector<Type>& rhs) { return solver.LeastSquaresSolve(rhs); });
			}
			case LinearSolverRecommendation::SVD: {
				SVDecompositionSolver<Type> solver(_coefficientMatrix);
				return solveColumns([&](const Vector<Type>& rhs) { return solver.Solve(rhs); });
			}
			case LinearSolverRecommendation::LU: {
				auto& solver = EnsureLUSolver();
				return solveColumns([&](const Vector<Type>& rhs) { return solver.Solve(rhs); });
			}
			}
			throw InvalidStateError("LinearSystem::SolveMultiple - invalid solver recommendation");
		}

		//=========================================================================
		// SPECIFIC SOLVER METHODS
		//=========================================================================

		/// @brief Solve using Gauss-Jordan elimination
		///
		/// @note Modifies internal copy; original matrix preserved
		Vector<Type> SolveByGaussJordan() const {
			RequireVectorRHS();
			RequireSquare();
			return GaussJordanSolver<Type>::SolveConst(_coefficientMatrix, _b);
		}

		/// @brief Solve using LU decomposition with partial pivoting
		///
		/// @note General-purpose, O(n³) factorization
		Vector<Type> SolveByLU() const {
			RequireVectorRHS();
			RequireSquare();
			LUSolver<Type> solver(_coefficientMatrix);
			return solver.Solve(_b);
		}

		/// @brief Solve using Cholesky decomposition
		///
		/// @throws SingularMatrixError if matrix is not positive definite
		/// @note Fastest and most stable for symmetric positive definite matrices
		Vector<Type> SolveByCholesky() const {
			RequireVectorRHS();
			RequireSquare();
			CholeskySolver<Type> solver(_coefficientMatrix);
			return solver.Solve(_b);
		}

		/// @brief Solve using QR decomposition
		///
		/// @note Works for overdetermined systems (least squares)
		Vector<Type> SolveByQR() const {
			RequireVectorRHS();
			QRSolver<Type> solver(_coefficientMatrix);
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
			RequireVectorRHS();
			SVDecompositionSolver<Type> solver(_coefficientMatrix);
			return solver.Solve(_b, threshold);
		}

		/// @brief Solve least squares problem min||Ax - b||
		///
		/// @note Uses QR decomposition; works for overdetermined systems
		Vector<Type> SolveLeastSquares() const {
			RequireVectorRHS();
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
			RequireVectorRHS();
			RequireSquare();

			if (method == IterativeMethod::Auto)
				method = SelectIterativeMethod();

			IterativeSolverResult result;

			switch (method) {
			case IterativeMethod::Jacobi:
				result = JacobiSolver::Solve(_coefficientMatrix, _b, Vector<Real>(), tol, maxIter);
				break;
			case IterativeMethod::GaussSeidel:
				result = GaussSeidelSolver::Solve(_coefficientMatrix, _b, Vector<Real>(), tol, maxIter);
				break;
			case IterativeMethod::SOR:
				// SOR: Solve(A, b, omega, x0, tol, maxIter) - omega=1.5 typical choice
				result = SORSolver::Solve(_coefficientMatrix, _b, static_cast<Type>(1.5), Vector<Real>(), tol, maxIter);
				break;
			default:
				result = GaussSeidelSolver::Solve(_coefficientMatrix, _b, Vector<Real>(), tol, maxIter);
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
			return MeasureVerification(x).residual.Value();
		}

		/// @brief Compute relative residual ||Ax - b|| / ||b||
		Type RelativeResidual(const Vector<Type>& x) const {
			const auto measures = MeasureVerification(x);
			return ScaledMagnitude::Ratio(measures.residual, measures.rhs);
		}

		/// @brief Full solution verification
		///
		/// @param x Solution to verify
		/// @param tol Accuracy threshold
		///
		/// @return Detailed verification results
		VerificationResult<Type> Verify(const Vector<Type>& x, Type tol = Precision::DefaultToleranceStrict) const {
			ValidateTolerance(tol, "LinearSystem::Verify");
			const auto measures = MeasureVerification(x);
			VerificationResult<Type> result;
			result.absoluteResidual = measures.residual.Value();
			result.relativeResidual = ScaledMagnitude::Ratio(measures.residual, measures.rhs);
			result.backwardError = ScaledMagnitude::Ratio(measures.residual,
				ScaledMagnitude::Product(measures.matrixInfinity, measures.solution));

			// Forward error estimate (based on condition number)
			Type cond = ConditionNumber();
			result.estimatedForwardError = std::isfinite(cond)
				? cond * result.relativeResidual : std::numeric_limits<Type>::quiet_NaN();

			result.isAccurate = std::isfinite(result.relativeResidual) && result.relativeResidual < tol;

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
					for (int col = 0; col < _rhsMatrix.cols(); ++col)
						result.solutionStatuses.push_back(ClassifySolution(_rhsMatrix.VectorFromColumn(col), coefficientSVD.rank, coefficientSVD.threshold));
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
		const Matrix<Type>& getMatrix() const { return _coefficientMatrix; }
		/// @deprecated Use getMatrix() (CONTRACT_POLICY section 11: accessors are getX()).
		[[deprecated("use getMatrix()")]]
		const Matrix<Type>& GetMatrix() const { return _coefficientMatrix; }

		/// @brief Get the right-hand side vector
		/// @throws InvalidStateError unless constructed with a vector RHS
		const Vector<Type>& getRHS() const {
			RequireVectorRHS();
			return _b;
		}
		/// @deprecated Use getRHS() (CONTRACT_POLICY section 11: accessors are getX()).
		[[deprecated("use getRHS()")]]
		const Vector<Type>& GetRHS() const { return getRHS(); }

		/// @brief Get the right-hand side matrix
		/// @throws InvalidStateError unless constructed with a matrix RHS
		const Matrix<Type>& getRHSMatrix() const {
			RequireMatrixRHS();
			return _rhsMatrix;
		}

	private:
		static constexpr int NoExponent = std::numeric_limits<int>::min(); // tol-ok: absent integer exponent, not a tolerance

		// Keep magnitudes normalized until the final ratio, even when a norm/product exceeds Type.
		struct ScaledMagnitude {
			Type fraction = 0;
			int exponent = 0;

			static ScaledMagnitude From(Type value, int scaleExponent = 0) {
				int exponent = 0;
				const Type fraction = std::frexp(std::abs(value), &exponent);
				return {fraction, fraction == Type(0) ? 0 : exponent + scaleExponent};
			}
			Type Value() const { return std::ldexp(fraction, exponent); }
			static ScaledMagnitude Product(const ScaledMagnitude& a, const ScaledMagnitude& b) {
				return From(a.fraction * b.fraction, a.exponent + b.exponent);
			}
			static Type Ratio(const ScaledMagnitude& numerator, const ScaledMagnitude& denominator) {
				if (denominator.fraction == Type(0))
					return numerator.fraction == Type(0) ? Type(0) : std::numeric_limits<Type>::infinity();
				return std::ldexp(numerator.fraction / denominator.fraction,
					numerator.exponent - denominator.exponent);
			}
		};

		struct VerificationMeasures {
			ScaledMagnitude residual, rhs, solution, matrixInfinity;
		};

		static ScaledMagnitude ScaledNorm(const Vector<Type>& values) {
			StableNormL2Accumulator norm;
			for (int i = 0; i < values.size(); ++i)
				norm.add(std::abs(values[i]));
			return ScaledMagnitude::Product(ScaledMagnitude::From(norm.scale),
				ScaledMagnitude::From(std::sqrt(norm.ssq)));
		}

		VerificationMeasures MeasureVerification(const Vector<Type>& x) const {
			RequireVectorRHS();
			if (x.size() != _coefficientMatrix.cols())
				throw VectorDimensionError("LinearSystem verification: solution size must match matrix columns",
					x.size(), _coefficientMatrix.cols());
			for (int j = 0; j < x.size(); ++j)
				ValidateFinite(x[j], "LinearSystem verification solution");
			for (int i = 0; i < _b.size(); ++i)
				ValidateFinite(_b[i], "LinearSystem verification RHS");

			std::vector<ScaledMagnitude> residuals(static_cast<size_t>(Rows()));
			ScaledMagnitude matrixInfinity;
			int residualExponent = NoExponent;
			for (int i = 0; i < Rows(); ++i) {
				int rowExponent = NoExponent;
				if (_b[i] != Type(0))
					rowExponent = ScaledMagnitude::From(_b[i]).exponent;
				Type matrixScale = 0;
				for (int j = 0; j < Cols(); ++j) {
					const Type coefficient = _coefficientMatrix(i, j);
					ValidateFinite(coefficient, "LinearSystem verification matrix");
					matrixScale = std::max(matrixScale, std::abs(coefficient));
					if (coefficient != Type(0) && x[j] != Type(0)) {
						const auto a = ScaledMagnitude::From(coefficient);
						const auto b = ScaledMagnitude::From(x[j]);
						rowExponent = std::max(rowExponent, a.exponent + b.exponent);
					}
				}
				Type rowSum = 0;
				if (matrixScale != Type(0)) {
					for (int j = 0; j < Cols(); ++j)
						rowSum += std::abs(_coefficientMatrix(i, j)) / matrixScale;
				}
				const auto rowNorm = ScaledMagnitude::Product(ScaledMagnitude::From(matrixScale),
					ScaledMagnitude::From(rowSum));
				if (rowNorm.fraction != Type(0) && (matrixInfinity.fraction == Type(0) ||
					rowNorm.exponent > matrixInfinity.exponent ||
					(rowNorm.exponent == matrixInfinity.exponent && rowNorm.fraction > matrixInfinity.fraction)))
					matrixInfinity = rowNorm;
				if (rowExponent == NoExponent)
					continue;

				Type sum = 0, correction = 0;
				auto accumulate = [&](Type term) {
					const Type next = sum + term;
					correction += std::abs(sum) >= std::abs(term) ? (sum - next) + term : (term - next) + sum;
					sum = next;
				};
				for (int j = 0; j < Cols(); ++j) {
					int aExponent = 0, xExponent = 0;
					const Type a = std::frexp(_coefficientMatrix(i, j), &aExponent);
					const Type value = std::frexp(x[j], &xExponent);
					const Type term = std::ldexp(a * value, aExponent + xExponent - rowExponent);
					if (term == Type(0) && a != Type(0) && value != Type(0))
						throw MatrixNumericalError("LinearSystem verification: residual term exceeds accumulation range");
					accumulate(term);
				}
				int bExponent = 0;
				const Type b = std::frexp(_b[i], &bExponent);
				const Type rhsTerm = std::ldexp(b, bExponent - rowExponent);
				if (rhsTerm == Type(0) && b != Type(0))
					throw MatrixNumericalError("LinearSystem verification: RHS exceeds residual accumulation range");
				accumulate(-rhsTerm);
				residuals[i] = ScaledMagnitude::From(sum + correction, rowExponent);
				if (residuals[i].fraction != Type(0))
					residualExponent = std::max(residualExponent, residuals[i].exponent);
			}
			StableNormL2Accumulator residual;
			for (const auto& component : residuals) {
				if (component.fraction != Type(0))
					residual.add(std::ldexp(component.fraction, component.exponent - residualExponent));
			}
			const auto residualNorm = ScaledMagnitude::From(residual.value(),
				residualExponent == NoExponent ? 0 : residualExponent);
			return {residualNorm, ScaledNorm(_b), ScaledNorm(x), matrixInfinity};
		}

		// Storage
		Matrix<Type> _coefficientMatrix;
		MatrixAnalyzer<Type> _matrixAnalyzer;
		Vector<Type> _b;
		Matrix<Type> _rhsMatrix; // For multiple RHS
		bool _hasRHS;
		bool _multipleRHS;

		// Cached solver used only for repeated RHS solves.
		mutable std::optional<LUSolver<Type>> _luSolver;

		//=========================================================================
		// INTERNAL HELPERS
		//=========================================================================

		void RequireVectorRHS() const {
			if (!_hasRHS || _multipleRHS)
				throw InvalidStateError("LinearSystem: operation requires a single vector right-hand side");
		}

		void RequireMatrixRHS() const {
			if (!_hasRHS || !_multipleRHS)
				throw InvalidStateError("LinearSystem: operation requires a matrix right-hand side");
		}

		void RequireSquare() const {
			if (!IsSquare())
				throw MatrixDimensionError("LinearSystem: operation requires square matrix", Rows(), Cols(), -1, -1);
		}

		/// @brief Get (or lazily create) the cached LU solver - factor once, solve many
		LUSolver<Type>& EnsureLUSolver() const {
			if (!_luSolver.has_value())
				_luSolver.emplace(_coefficientMatrix);
			return *_luSolver;
		}

		LinearSolverRecommendation SelectBestSolver(Threshold threshold = std::nullopt) const {
			const auto stability = AssessStability(threshold);
			if (Rank(threshold) < Cols() || stability == MatrixAlg::MatrixStability::IllConditioned ||
				stability == MatrixAlg::MatrixStability::Singular)
				return LinearSolverRecommendation::SVD;

			// Approximate property queries must not discard entries during a solve.
			constexpr ComparisonTolerance exactStructure{};
			if (IsSquare() && (IsUpperTriangular(exactStructure) || IsLowerTriangular(exactStructure)))
				return LinearSolverRecommendation::Triangular;

			if (IsTall())
				return LinearSolverRecommendation::QR;

			// For SPD matrices, Cholesky is best
			if (IsSymmetric(exactStructure) && IsPositiveDefinite())
				return LinearSolverRecommendation::Cholesky;

			// Default: LU with partial pivoting
			return LinearSolverRecommendation::LU;
		}

		SolutionStatus ClassifySolution(const Vector<Type>& rhs, int coefficientRank, Magnitude threshold) const {
			Matrix<Type> augmented(Rows(), Cols() + 1);
			for (int row = 0; row < Rows(); ++row) {
				for (int col = 0; col < Cols(); ++col)
					augmented(row, col) = _coefficientMatrix(row, col);
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

		Vector<Type> SolveTriangular(const Vector<Type>& rhs) const {
			int n = Rows();
			Vector<Type> x(n);

			if (IsUpperTriangular(ComparisonTolerance{})) {
				// Back substitution
				for (int i = n - 1; i >= 0; --i) {
					Type sum = rhs[i];
					for (int j = i + 1; j < n; ++j)
						sum -= _coefficientMatrix(i, j) * x[j];
					x[i] = sum / _coefficientMatrix(i, i);
				}
			} else // Lower triangular
			{
				// Forward substitution
				for (int i = 0; i < n; ++i) {
					Type sum = rhs[i];
					for (int j = 0; j < i; ++j)
						sum -= _coefficientMatrix(i, j) * x[j];
					x[i] = sum / _coefficientMatrix(i, i);
				}
			}

			return x;
		}

		std::string GenerateReport(const SystemAnalysis<Type>& analysis) const {
			std::ostringstream matrixReport;
			AnalysisReports::PrintMatrixAnalysis(matrixReport, analysis.matrix);
			std::string report = matrixReport.str();
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

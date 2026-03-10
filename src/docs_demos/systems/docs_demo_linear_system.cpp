#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "systems/LinearSystem.h"
#endif

#include "../test_beds/linear_alg_eq_systems_test_bed.h"

using namespace MML;
using namespace MML::Systems;

///////////////////////////////////////////////////////////////////////////////////////////
///                         BASIC USAGE                                                 ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_LinearSystem_BasicUsage()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: LinearSystem - Basic Usage\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nLinearSystem is MML's unified facade for ALL linear algebra operations.\n";
	std::cout << "One class gives you: solving, analysis, decompositions, eigenvalues.\n";
	
	// Simple 3x3 system
	Matrix<Real> A(3, 3, {4, 1, 0,
	                      1, 4, 1,
	                      0, 1, 4});
	Vector<Real> b({5, 6, 5});
	
	std::cout << "\n--- Simple 3x3 System ---\n";
	std::cout << "Matrix A (tridiagonal SPD):\n"; A.Print(std::cout, 8, 3);
	std::cout << "RHS b = "; b.Print(std::cout, 8, 3); std::cout << "\n";
	
	// Create system and solve
	LinearSystem<Real> sys(A, b);
	Vector<Real> x = sys.Solve();  // Auto-selects best solver
	
	std::cout << "\nSolution x = "; x.Print(std::cout, 10, 6); std::cout << "\n";
	
	// Quick verification
	std::cout << "Residual ||Ax - b|| = " << sys.ResidualNorm(x) << "\n";
	std::cout << "Relative residual = " << sys.RelativeResidual(x) << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         AUTO SOLVER SELECTION                                       ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_LinearSystem_AutoSolverSelection()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: LinearSystem - Automatic Solver Selection\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nLinearSystem.Solve() intelligently selects the best algorithm:\n";
	std::cout << "- SPD matrices → Cholesky (fastest, most stable)\n";
	std::cout << "- Triangular → Back/forward substitution (O(n²))\n";
	std::cout << "- Overdetermined → QR (least squares)\n";
	std::cout << "- Ill-conditioned → SVD (most robust)\n";
	std::cout << "- General → LU with partial pivoting\n";
	
	// Example 1: SPD matrix → will use Cholesky
	std::cout << "\n--- Symmetric Positive Definite ---\n";
	Matrix<Real> A_spd(3, 3, {4, 2, 0,
	                          2, 5, 2,
	                          0, 2, 4});
	Vector<Real> b1({8, 13, 8});
	LinearSystem<Real> sys1(A_spd, b1);
	
	std::cout << "Matrix type: SPD (will use Cholesky)\n";
	std::cout << "IsSymmetric: " << (sys1.IsSymmetric() ? "yes" : "no") << "\n";
	std::cout << "IsPositiveDefinite: " << (sys1.IsPositiveDefinite() ? "yes" : "no") << "\n";
	
	Vector<Real> x1 = sys1.Solve();
	std::cout << "Solution: "; x1.Print(std::cout, 10, 6); std::cout << "\n";
	std::cout << "Residual: " << sys1.ResidualNorm(x1) << "\n";
	
	// Example 2: Upper triangular → back substitution
	std::cout << "\n--- Upper Triangular ---\n";
	Matrix<Real> A_tri(3, 3, {2, 1, 3,
	                          0, 4, 2,
	                          0, 0, 5});
	Vector<Real> b2({13, 14, 10});
	LinearSystem<Real> sys2(A_tri, b2);
	
	std::cout << "Matrix type: Upper triangular (direct substitution)\n";
	std::cout << "IsUpperTriangular: " << (sys2.IsUpperTriangular() ? "yes" : "no") << "\n";
	
	Vector<Real> x2 = sys2.Solve();
	std::cout << "Solution: "; x2.Print(std::cout, 10, 6); std::cout << "\n";
	
	// Example 3: General matrix → LU
	std::cout << "\n--- General Matrix ---\n";
	Matrix<Real> A_gen(3, 3, {1, 2, 3,
	                          4, 5, 6,
	                          7, 8, 10});  // Not triangular, not symmetric
	Vector<Real> b3({14, 32, 53});
	LinearSystem<Real> sys3(A_gen, b3);
	
	std::cout << "Matrix type: General (will use LU)\n";
	Vector<Real> x3 = sys3.Solve();
	std::cout << "Solution: "; x3.Print(std::cout, 10, 6); std::cout << "\n";
	std::cout << "Residual: " << sys3.ResidualNorm(x3) << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         SPECIFIC SOLVERS                                            ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_LinearSystem_SpecificSolvers()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: LinearSystem - Specific Solver Methods\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nYou can force specific solvers for comparison or teaching purposes.\n";
	
	Matrix<Real> A = TestBeds::mat_5x5();
	Vector<Real> b = TestBeds::mat_5x5_rhs0();
	LinearSystem<Real> sys(A, b);
	
	std::cout << "\n--- Comparing Solvers on 5x5 System ---\n";
	std::cout << "Matrix A:\n"; A.Print(std::cout, 10, 4);
	std::cout << "RHS b = "; b.Print(std::cout, 10, 4); std::cout << "\n\n";
	
	// Solve with each method
	std::cout << std::setw(20) << "Method" << std::setw(18) << "Residual" << "\n";
	std::cout << std::string(38, '-') << "\n";
	
	Vector<Real> x_auto = sys.Solve();
	std::cout << std::setw(20) << "Auto (best)" << std::setw(18) << std::scientific 
	          << sys.ResidualNorm(x_auto) << "\n";
	
	Vector<Real> x_gj = sys.SolveByGaussJordan();
	std::cout << std::setw(20) << "Gauss-Jordan" << std::setw(18) << std::scientific 
	          << sys.ResidualNorm(x_gj) << "\n";
	
	Vector<Real> x_lu = sys.SolveByLU();
	std::cout << std::setw(20) << "LU" << std::setw(18) << std::scientific 
	          << sys.ResidualNorm(x_lu) << "\n";
	
	Vector<Real> x_qr = sys.SolveByQR();
	std::cout << std::setw(20) << "QR" << std::setw(18) << std::scientific 
	          << sys.ResidualNorm(x_qr) << "\n";
	
	Vector<Real> x_svd = sys.SolveBySVD();
	std::cout << std::setw(20) << "SVD" << std::setw(18) << std::scientific 
	          << sys.ResidualNorm(x_svd) << "\n";
	
	std::cout << std::fixed;  // Reset formatting
	
	std::cout << "\nAll methods give essentially the same solution for well-conditioned systems.\n";
	std::cout << "Differences become significant for ill-conditioned systems.\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         SOLUTION VERIFICATION                                       ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_LinearSystem_Verification()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: LinearSystem - Solution Verification\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nAlways verify your solutions! LinearSystem provides comprehensive checks.\n";
	
	Matrix<Real> A = TestBeds::mat_5x5();
	Vector<Real> b = TestBeds::mat_5x5_rhs0();
	LinearSystem<Real> sys(A, b);
	
	Vector<Real> x = sys.Solve();
	
	std::cout << "\n--- Quick Verification ---\n";
	std::cout << "Solution x = "; x.Print(std::cout, 10, 6); std::cout << "\n";
	std::cout << "Absolute residual ||Ax - b|| = " << sys.ResidualNorm(x) << "\n";
	std::cout << "Relative residual ||Ax - b|| / ||b|| = " << sys.RelativeResidual(x) << "\n";
	
	std::cout << "\n--- Full Verification Report ---\n";
	auto verify = sys.Verify(x);
	
	std::cout << "Absolute residual: " << verify.absoluteResidual << "\n";
	std::cout << "Relative residual: " << verify.relativeResidual << "\n";
	std::cout << "Backward error: " << verify.backwardError << "\n";
	std::cout << "Estimated forward error: " << verify.estimatedForwardError << "\n";
	std::cout << "Is accurate (rel < 1e-10): " << (verify.isAccurate ? "Yes" : "No") << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         MATRIX ANALYSIS                                             ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_LinearSystem_MatrixAnalysis()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: LinearSystem - Matrix Analysis\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nLinearSystem provides rich matrix analysis without solving.\n";
	
	Matrix<Real> A(4, 4, {10, 1, 2, 0,
	                       1, 8, 1, 2,
	                       2, 1, 7, 1,
	                       0, 2, 1, 6});
	
	LinearSystem<Real> sys(A);  // No RHS - just analysis
	
	std::cout << "\nMatrix A:\n"; A.Print(std::cout, 8, 3);
	
	std::cout << "\n--- Dimension Properties ---\n";
	std::cout << "Rows: " << sys.Rows() << "\n";
	std::cout << "Cols: " << sys.Cols() << "\n";
	std::cout << "IsSquare: " << (sys.IsSquare() ? "yes" : "no") << "\n";
	std::cout << "IsOverdetermined: " << (sys.IsOverdetermined() ? "yes" : "no") << "\n";
	std::cout << "IsUnderdetermined: " << (sys.IsUnderdetermined() ? "yes" : "no") << "\n";
	
	std::cout << "\n--- Structure Detection ---\n";
	std::cout << "IsSymmetric: " << (sys.IsSymmetric() ? "yes" : "no") << "\n";
	std::cout << "IsPositiveDefinite: " << (sys.IsPositiveDefinite() ? "yes" : "no") << "\n";
	std::cout << "IsDiagonallyDominant: " << (sys.IsDiagonallyDominant() ? "yes" : "no") << "\n";
	std::cout << "IsUpperTriangular: " << (sys.IsUpperTriangular() ? "yes" : "no") << "\n";
	std::cout << "IsLowerTriangular: " << (sys.IsLowerTriangular() ? "yes" : "no") << "\n";
	std::cout << "IsDiagonal: " << (sys.IsDiagonal() ? "yes" : "no") << "\n";
	std::cout << "Sparsity: " << std::fixed << std::setprecision(2) 
	          << (sys.Sparsity() * 100) << "%\n";
	
	std::cout << "\n--- Numerical Properties ---\n";
	std::cout << "Determinant: " << sys.Determinant() << "\n";
	std::cout << "Rank: " << sys.Rank() << "\n";
	std::cout << "Nullity: " << sys.Nullity() << "\n";
	std::cout << "Condition number: " << std::scientific << sys.ConditionNumber() << "\n";
	std::cout << "Expected digits lost: " << sys.ExpectedDigitsLost() << "\n";
	std::cout << std::fixed;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         STABILITY ASSESSMENT                                        ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_LinearSystem_StabilityAssessment()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: LinearSystem - Stability Assessment\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nCondition number determines how errors are amplified:\n";
	std::cout << "- Well-conditioned (κ < 10⁴): Accurate results expected\n";
	std::cout << "- Moderate (10⁴ ≤ κ < 10⁸): Some precision loss\n";
	std::cout << "- Ill-conditioned (10⁸ ≤ κ < 10¹²): Significant errors possible\n";
	std::cout << "- Singular (κ ≥ 10¹²): Results unreliable\n";
	
	// Well-conditioned matrix
	std::cout << "\n--- Well-Conditioned Matrix ---\n";
	Matrix<Real> A_good(3, 3, {4, 1, 1,
	                           1, 3, 1,
	                           1, 1, 4});
	LinearSystem<Real> sys_good(A_good);
	
	MatrixStability stability = sys_good.AssessStability();
	std::cout << "Condition number: " << std::scientific << sys_good.ConditionNumber() << "\n";
	std::cout << "Stability: ";
	switch (stability) {
		case MatrixStability::WellConditioned: std::cout << "Well-Conditioned"; break;
		case MatrixStability::ModeratelyConditioned: std::cout << "Moderately Conditioned"; break;
		case MatrixStability::IllConditioned: std::cout << "Ill-Conditioned"; break;
		case MatrixStability::Singular: std::cout << "Singular"; break;
	}
	std::cout << "\n";
	std::cout << "Expected digits lost: " << sys_good.ExpectedDigitsLost() << "\n";
	
	// Create an ill-conditioned Hilbert matrix
	std::cout << "\n--- Ill-Conditioned Matrix (Hilbert 5x5) ---\n";
	Matrix<Real> H(5, 5);
	for (int i = 0; i < 5; i++)
		for (int j = 0; j < 5; j++)
			H(i, j) = 1.0 / (i + j + 1);
	
	LinearSystem<Real> sys_hilbert(H);
	std::cout << "Hilbert matrix H(i,j) = 1/(i+j+1)\n";
	std::cout << "Condition number: " << std::scientific << sys_hilbert.ConditionNumber() << "\n";
	
	stability = sys_hilbert.AssessStability();
	std::cout << "Stability: ";
	switch (stability) {
		case MatrixStability::WellConditioned: std::cout << "Well-Conditioned"; break;
		case MatrixStability::ModeratelyConditioned: std::cout << "Moderately Conditioned"; break;
		case MatrixStability::IllConditioned: std::cout << "Ill-Conditioned"; break;
		case MatrixStability::Singular: std::cout << "Singular"; break;
	}
	std::cout << "\n";
	std::cout << "Expected digits lost: " << sys_hilbert.ExpectedDigitsLost() << "\n";
	std::cout << std::fixed;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         COMPREHENSIVE ANALYSIS                                      ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_LinearSystem_ComprehensiveAnalysis()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: LinearSystem - Comprehensive Analysis\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nThe Analyze() method provides a complete system report.\n";
	
	Matrix<Real> A(4, 4, {5, 1, 0, 0,
	                      1, 5, 1, 0,
	                      0, 1, 5, 1,
	                      0, 0, 1, 5});
	
	LinearSystem<Real> sys(A);
	auto analysis = sys.Analyze();
	
	std::cout << "\n" << analysis.analysisReport << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         DECOMPOSITIONS                                              ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_LinearSystem_Decompositions()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: LinearSystem - Cached Decompositions\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nDecompositions are computed lazily and cached for reuse.\n";
	
	Matrix<Real> A(3, 3, {4, 2, 0,
	                      2, 5, 2,
	                      0, 2, 4});
	LinearSystem<Real> sys(A);
	
	std::cout << "\nMatrix A:\n"; A.Print(std::cout, 8, 3);
	
	// QR decomposition
	std::cout << "\n--- QR Decomposition: A = Q*R ---\n";
	auto qr = sys.GetQR();
	if (qr.valid) {
		std::cout << "Q (orthogonal):\n"; qr.Q.Print(std::cout, 10, 6);
		std::cout << "R (upper triangular):\n"; qr.R.Print(std::cout, 10, 6);
		
		// Verify: Q*R should equal A
		Matrix<Real> QR = qr.Q * qr.R;
		std::cout << "Verification Q*R:\n"; QR.Print(std::cout, 10, 6);
	}
	
	// SVD decomposition
	std::cout << "\n--- SVD Decomposition: A = U*diag(w)*V^T ---\n";
	auto svd = sys.GetSVD();
	if (svd.valid) {
		std::cout << "Singular values: "; svd.singularValues.Print(std::cout, 10, 6); 
		std::cout << "\n";
		std::cout << "Rank: " << svd.rank << "\n";
	}
	
	// Cholesky (for SPD)
	std::cout << "\n--- Cholesky Decomposition: A = L*L^T (SPD only) ---\n";
	auto chol = sys.GetCholesky();
	if (chol.valid) {
		std::cout << "L (lower triangular):\n"; chol.L.Print(std::cout, 10, 6);
		
		// Verify: L*L^T should equal A
		Matrix<Real> LLT = chol.L * chol.L.transpose();
		std::cout << "Verification L*L^T:\n"; LLT.Print(std::cout, 10, 6);
	}
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         EIGENANALYSIS                                               ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_LinearSystem_Eigenanalysis()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: LinearSystem - Eigenanalysis\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nLinearSystem provides access to eigenvalues and eigenvectors.\n";
	
	// Symmetric matrix (guaranteed real eigenvalues)
	std::cout << "\n--- Symmetric Matrix (all real eigenvalues) ---\n";
	Matrix<Real> A_sym(3, 3, {4, 1, 0,
	                          1, 4, 1,
	                          0, 1, 4});
	LinearSystem<Real> sys_sym(A_sym);
	
	std::cout << "Matrix A:\n"; A_sym.Print(std::cout, 8, 3);
	
	Vector<Real> eigs = sys_sym.EigenvaluesSymmetric();
	std::cout << "Eigenvalues (symmetric solver): "; eigs.Print(std::cout, 10, 6);
	std::cout << "\n";
	
	std::cout << "Spectral radius (max |λ|): " << sys_sym.SpectralRadius() << "\n";
	std::cout << "Has complex eigenvalues: " << (sys_sym.HasComplexEigenvalues() ? "yes" : "no") << "\n";
	
	// General matrix (may have complex eigenvalues)
	std::cout << "\n--- General Matrix (may have complex eigenvalues) ---\n";
	Matrix<Real> A_gen(3, 3, {0, -1, 0,
	                          1,  0, 0,
	                          0,  0, 2});
	LinearSystem<Real> sys_gen(A_gen);
	
	std::cout << "Matrix A:\n"; A_gen.Print(std::cout, 8, 3);
	
	auto eigen = sys_gen.GetEigen();
	std::cout << "Eigenvalues:\n";
	for (size_t i = 0; i < eigen.eigenvalues.size(); i++) {
		std::cout << "  λ" << i << " = " << eigen.eigenvalues[i].real;
		if (std::abs(eigen.eigenvalues[i].imag) > 1e-10)
			std::cout << " + " << eigen.eigenvalues[i].imag << "i";
		std::cout << "\n";
	}
	
	std::cout << "Has complex eigenvalues: " << (sys_gen.HasComplexEigenvalues() ? "yes" : "no") << "\n";
	std::cout << "Spectral radius: " << sys_gen.SpectralRadius() << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         FUNDAMENTAL SUBSPACES                                       ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_LinearSystem_FundamentalSubspaces()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: LinearSystem - Fundamental Subspaces\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nCompute null space and column space bases using SVD.\n";
	
	// Rank-deficient matrix
	Matrix<Real> A(3, 3, {1, 2, 3,
	                      2, 4, 6,
	                      1, 1, 2});
	
	LinearSystem<Real> sys(A);
	
	std::cout << "\nMatrix A (rank-deficient):\n"; A.Print(std::cout, 8, 3);
	std::cout << "Rank: " << sys.Rank() << "\n";
	std::cout << "Nullity: " << sys.Nullity() << "\n";
	
	if (sys.Nullity() > 0) {
		std::cout << "\nNull space basis (orthonormal columns):\n";
		Matrix<Real> nullspace = sys.NullSpace();
		nullspace.Print(std::cout, 10, 6);
		
		std::cout << "\nVerification A*null_basis ≈ 0:\n";
		Matrix<Real> Az = A * nullspace;
		Az.Print(std::cout, 12, 8);
	}
	
	std::cout << "\nColumn space basis (orthonormal columns):\n";
	Matrix<Real> colspace = sys.ColumnSpace();
	colspace.Print(std::cout, 10, 6);
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         MATRIX INVERSE                                              ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_LinearSystem_MatrixOperations()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: LinearSystem - Matrix Operations\n";
	std::cout << "==========================================================================\n";
	
	Matrix<Real> A(3, 3, {4, 1, 0,
	                      1, 4, 1,
	                      0, 1, 4});
	LinearSystem<Real> sys(A);
	
	std::cout << "\nMatrix A:\n"; A.Print(std::cout, 8, 3);
	
	// Matrix inverse
	std::cout << "\n--- Matrix Inverse ---\n";
	Matrix<Real> A_inv = sys.Inverse();
	std::cout << "A^-1:\n"; A_inv.Print(std::cout, 10, 6);
	
	// Verify
	Matrix<Real> I = A * A_inv;
	std::cout << "Verification A*A^-1:\n"; I.Print(std::cout, 10, 6);
	
	// Pseudo-inverse for non-square
	std::cout << "\n--- Pseudo-Inverse (non-square) ---\n";
	Matrix<Real> B(4, 2, {1, 2,
	                      3, 4,
	                      5, 6,
	                      7, 8});
	LinearSystem<Real> sys_rect(B);
	
	std::cout << "Matrix B (4x2):\n"; B.Print(std::cout, 8, 3);
	
	Matrix<Real> B_pinv = sys_rect.PseudoInverse();
	std::cout << "B+ (pseudo-inverse, 2x4):\n"; B_pinv.Print(std::cout, 10, 6);
	
	// Verify: B+ * B should be close to identity (2x2)
	Matrix<Real> BpB = B_pinv * B;
	std::cout << "Verification B+*B:\n"; BpB.Print(std::cout, 10, 6);
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         OVERDETERMINED SYSTEMS (LEAST SQUARES)                      ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_LinearSystem_LeastSquares()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: LinearSystem - Least Squares (Overdetermined Systems)\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nOverdetermined systems (more equations than unknowns) are solved\n";
	std::cout << "using least squares: minimize ||Ax - b||²\n";
	
	// 5 data points, fit a line y = mx + c
	std::cout << "\n--- Linear Regression Example ---\n";
	std::cout << "Data points: (0,1), (1,2.5), (2,3.8), (3,5.2), (4,6.9)\n";
	std::cout << "Fitting y = c₀ + c₁*x\n";
	
	// Design matrix [1, x] for each point
	Matrix<Real> A(5, 2, {1, 0,
	                      1, 1,
	                      1, 2,
	                      1, 3,
	                      1, 4});
	Vector<Real> b({1.0, 2.5, 3.8, 5.2, 6.9});
	
	LinearSystem<Real> sys(A, b);
	
	std::cout << "Design matrix A:\n"; A.Print(std::cout, 8, 3);
	std::cout << "Observations b = "; b.Print(std::cout, 8, 3); std::cout << "\n";
	
	std::cout << "System is overdetermined: " << (sys.IsOverdetermined() ? "yes" : "no") << "\n";
	
	Vector<Real> x = sys.SolveLeastSquares();
	
	std::cout << "\nLeast squares solution:\n";
	std::cout << "  c₀ (intercept) = " << x[0] << "\n";
	std::cout << "  c₁ (slope) = " << x[1] << "\n";
	std::cout << "\nFitted line: y = " << x[0] << " + " << x[1] << "*x\n";
	
	// Show residuals
	std::cout << "\nResidual ||Ax - b|| = " << sys.ResidualNorm(x) << "\n";
	
	// Predictions
	std::cout << "\nPredicted vs Actual:\n";
	Vector<Real> pred = A * x;
	for (int i = 0; i < 5; i++) {
		std::cout << "  x=" << i << ": predicted=" << std::fixed << std::setprecision(3) 
		          << pred[i] << ", actual=" << b[i] 
		          << ", error=" << std::abs(pred[i] - b[i]) << "\n";
	}
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         ITERATIVE METHODS                                           ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_LinearSystem_IterativeMethods()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: LinearSystem - Iterative Methods\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nIterative methods are useful for large sparse systems:\n";
	std::cout << "- Jacobi: Simple, parallelizable\n";
	std::cout << "- Gauss-Seidel: Faster convergence\n";
	std::cout << "- SOR: Even faster with good ω parameter\n";
	
	// Diagonally dominant matrix (ensures convergence)
	Matrix<Real> A(4, 4, {10, 1, 1, 1,
	                       1, 10, 1, 1,
	                       1, 1, 10, 1,
	                       1, 1, 1, 10});
	Vector<Real> b({13, 13, 13, 13});
	
	LinearSystem<Real> sys(A, b);
	
	std::cout << "\nDiagonally dominant matrix (guarantees convergence):\n";
	A.Print(std::cout, 8, 3);
	std::cout << "IsDiagonallyDominant: " << (sys.IsDiagonallyDominant() ? "yes" : "no") << "\n";
	
	std::cout << "\n--- Comparing Iterative Methods ---\n";
	
	// Direct solution for reference
	Vector<Real> x_direct = sys.Solve();
	std::cout << "Direct solution: "; x_direct.Print(std::cout, 10, 6); std::cout << "\n";
	
	// Jacobi
	try {
		Vector<Real> x_jacobi = sys.SolveIterative(IterativeMethod::Jacobi);
		std::cout << "Jacobi solution: "; x_jacobi.Print(std::cout, 10, 6); std::cout << "\n";
	} catch (const ConvergenceError& e) {
		std::cout << "Jacobi: " << e.what() << "\n";
	}
	
	// Gauss-Seidel
	try {
		Vector<Real> x_gs = sys.SolveIterative(IterativeMethod::GaussSeidel);
		std::cout << "Gauss-Seidel: "; x_gs.Print(std::cout, 10, 6); std::cout << "\n";
	} catch (const ConvergenceError& e) {
		std::cout << "Gauss-Seidel: " << e.what() << "\n";
	}
	
	// SOR
	try {
		Vector<Real> x_sor = sys.SolveIterative(IterativeMethod::SOR);
		std::cout << "SOR (ω=1.5): "; x_sor.Print(std::cout, 10, 6); std::cout << "\n";
	} catch (const ConvergenceError& e) {
		std::cout << "SOR: " << e.what() << "\n";
	}
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         MULTIPLE RIGHT-HAND SIDES                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_LinearSystem_MultipleRHS()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: LinearSystem - Multiple Right-Hand Sides\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nEfficiently solve Ax₁=b₁, Ax₂=b₂, ... using one factorization.\n";
	
	Matrix<Real> A(3, 3, {4, 1, 0,
	                      1, 4, 1,
	                      0, 1, 4});
	
	// Multiple RHS as columns of B
	Matrix<Real> B(3, 3, {5, 1, 0,
	                      6, 0, 1,
	                      5, 0, 0});
	
	std::cout << "Matrix A:\n"; A.Print(std::cout, 8, 3);
	std::cout << "RHS matrix B (each column is a right-hand side):\n"; 
	B.Print(std::cout, 8, 3);
	
	LinearSystem<Real> sys(A, B);
	Matrix<Real> X = sys.SolveMultiple();
	
	std::cout << "\nSolution matrix X (each column is a solution):\n";
	X.Print(std::cout, 10, 6);
	
	// Verify
	std::cout << "\nVerification A*X:\n";
	Matrix<Real> AX = A * X;
	AX.Print(std::cout, 10, 6);
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         MAIN ENTRY POINT                                            ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_LinearSystem()
{
	std::cout << "\n";
	std::cout << "##########################################################################\n";
	std::cout << "###                  LinearSystem Documentation Demos                 ###\n";
	std::cout << "##########################################################################\n";
	std::cout << "\nLinearSystem is MML's unified facade for ALL linear algebra capabilities.\n";
	std::cout << "It provides: solving, analysis, decompositions, eigenvalues, and more.\n";
	
	Docs_Demo_LinearSystem_BasicUsage();
	Docs_Demo_LinearSystem_AutoSolverSelection();
	Docs_Demo_LinearSystem_SpecificSolvers();
	Docs_Demo_LinearSystem_Verification();
	Docs_Demo_LinearSystem_MatrixAnalysis();
	Docs_Demo_LinearSystem_StabilityAssessment();
	Docs_Demo_LinearSystem_ComprehensiveAnalysis();
	Docs_Demo_LinearSystem_Decompositions();
	Docs_Demo_LinearSystem_Eigenanalysis();
	Docs_Demo_LinearSystem_FundamentalSubspaces();
	Docs_Demo_LinearSystem_MatrixOperations();
	Docs_Demo_LinearSystem_LeastSquares();
	Docs_Demo_LinearSystem_IterativeMethods();
	Docs_Demo_LinearSystem_MultipleRHS();
	
	std::cout << "\n";
	std::cout << "##########################################################################\n";
	std::cout << "###              LinearSystem Demos Completed                         ###\n";
	std::cout << "##########################################################################\n";
}

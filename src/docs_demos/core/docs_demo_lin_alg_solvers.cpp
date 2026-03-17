#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "core/LinAlgEqSolvers.h"

#endif

#include "../test_beds/linear_alg_eq_systems_test_bed.h"

using namespace MML;

///////////////////////////////////////////////////////////////////////////////////////////
///                       GAUSS-JORDAN ELIMINATION                                     ///
///////////////////////////////////////////////////////////////////////////////////////////

// Demo: Gauss-Jordan solver for single and multiple RHS
void Docs_Demo_LinAlgSolvers_GaussJordan()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Gauss-Jordan Elimination\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nGaussJordanSolver - Full pivoting, computes A^-1 as byproduct\n";
	std::cout << "Best for: Small systems, when matrix inverse is also needed\n";
	
	// Example 1: Single RHS
	std::cout << "\n--- Single Right-Hand Side ---\n";
	Matrix<Real> A = TestBeds::mat_5x5();
	Matrix<Real> A_copy(A);
	Vector<Real> b = TestBeds::mat_5x5_rhs0();
	Vector<Real> b_copy(b);
	
	std::cout << "System A*x = b, where A is 5x5:\n";
	std::cout << "  A =\n"; A.Print(std::cout, 8, 3);
	std::cout << "  b = "; b.Print(std::cout, 8, 3); std::cout << std::endl;
	
	GaussJordanSolver<Real>::SolveInPlace(A_copy, b_copy);
	
	std::cout << "\nSolution x = "; b_copy.Print(std::cout, 10, 6); std::cout << std::endl;
	
	// Verify: A*x should equal b
	Vector<Real> Ax = A * b_copy;
	std::cout << "Verification A*x = "; Ax.Print(std::cout, 8, 3); std::cout << std::endl;
	std::cout << "Residual ||A*x - b|| = " << (Ax - b).NormL2() << std::endl;
	
	// Example 2: Matrix inversion (byproduct)
	std::cout << "\n--- Matrix Inversion (byproduct) ---\n";
	std::cout << "After Gauss-Jordan, A_copy contains A^-1:\n";
	std::cout << "  A^-1 =\n"; A_copy.Print(std::cout, 8, 3);
	
	Matrix<Real> AA_inv = A * A_copy;
	std::cout << "Verification A*A^-1 (approx I):\n"; AA_inv.Print(std::cout, 10, 2);
	
	// Example 3: Multiple RHS
	std::cout << "\n--- Multiple Right-Hand Sides ---\n";
	A_copy = A;
	Matrix<Real> B = TestBeds::mat_5x5_rhs_multi();
	std::cout << "Solving A*X = B where B has 2 columns:\n";
	
	GaussJordanSolver<Real>::SolveInPlace(A_copy, B);
	std::cout << "Solution X =\n"; B.Print(std::cout, 10, 4);
}

///////////////////////////////////////////////////////////////////////////////////////////
///                       LU DECOMPOSITION                                             ///
///////////////////////////////////////////////////////////////////////////////////////////

// Demo: LU decomposition solver
void Docs_Demo_LinAlgSolvers_LU()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: LU Decomposition Solver\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nLUSolver - Factor A = L*U once, solve for multiple RHS efficiently\n";
	std::cout << "Best for: Multiple RHS, general square systems\n";
	
	// Create solver (performs decomposition)
	Matrix<Real> A = TestBeds::mat_5x5();
	LUSolver<Real> lu(A);
	
	std::cout << "\n--- Solve single RHS ---\n";
	Vector<Real> b = TestBeds::mat_5x5_rhs0();
	Vector<Real> x = lu.Solve(b);
	
	std::cout << "Solution x = "; x.Print(std::cout, 10, 6); std::cout << std::endl;
	std::cout << "Residual ||A*x - b|| = " << (A*x - b).NormL2() << std::endl;
	
	// Multiple RHS (efficient - reuses decomposition)
	std::cout << "\n--- Solve multiple RHS (reusing decomposition) ---\n";
	Vector<Real> b2({1.0, 0.0, 0.0, 0.0, 0.0});
	Vector<Real> b3({0.0, 1.0, 0.0, 0.0, 0.0});
	
	Vector<Real> x2 = lu.Solve(b2);
	Vector<Real> x3 = lu.Solve(b3);
	
	std::cout << "Solve for e1: x = "; x2.Print(std::cout, 10, 4); std::cout << std::endl;
	std::cout << "Solve for e2: x = "; x3.Print(std::cout, 10, 4); std::cout << std::endl;
	
	// Determinant (available from LU)
	std::cout << "\n--- Determinant from LU ---\n";
	std::cout << "det(A) = " << lu.det() << std::endl;
	
	// Matrix inverse
	std::cout << "\n--- Matrix Inverse via LU ---\n";
	Matrix<Real> A_inv(5, 5);
	lu.inverse(A_inv);
	Matrix<Real> I = A * A_inv;
	Real err = 0;
	for (int i = 0; i < 5; i++)
		for (int j = 0; j < 5; j++)
			err += std::abs(I(i,j) - (i==j ? 1.0 : 0.0));
	std::cout << "A^-1 computed. Verification ||A*A^-1 - I||_1 = " << err << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                       QR DECOMPOSITION                                             ///
///////////////////////////////////////////////////////////////////////////////////////////

// Demo: QR decomposition solver
void Docs_Demo_LinAlgSolvers_QR()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: QR Decomposition Solver\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nQRSolver - Uses Householder reflections, excellent stability\n";
	std::cout << "Best for: Ill-conditioned systems, overdetermined systems (least squares)\n";
	
	// Square system
	std::cout << "\n--- Square System ---\n";
	Matrix<Real> A = TestBeds::mat_5x5();
	QRSolver<Real> qr(A);
	
	Vector<Real> b = TestBeds::mat_5x5_rhs0();
	Vector<Real> x = qr.Solve(b);
	
	std::cout << "Solution x = "; x.Print(std::cout, 10, 6); std::cout << std::endl;
	std::cout << "Residual ||A*x - b|| = " << (A*x - b).NormL2() << std::endl;
	
	// Overdetermined system (least squares)
	std::cout << "\n--- Overdetermined System (Least Squares) ---\n";
	// 6x4 system (more equations than unknowns)
	Matrix<Real> A_over(6, 4);
	A_over(0,0)=1; A_over(0,1)=2; A_over(0,2)=3; A_over(0,3)=4;
	A_over(1,0)=5; A_over(1,1)=6; A_over(1,2)=7; A_over(1,3)=8;
	A_over(2,0)=9; A_over(2,1)=10; A_over(2,2)=11; A_over(2,3)=12;
	A_over(3,0)=2; A_over(3,1)=1; A_over(3,2)=3; A_over(3,3)=2;
	A_over(4,0)=1; A_over(4,1)=3; A_over(4,2)=2; A_over(4,3)=1;
	A_over(5,0)=3; A_over(5,1)=2; A_over(5,2)=1; A_over(5,3)=3;
	
	Vector<Real> b_over({1.0, 2.0, 3.0, 1.5, 1.0, 2.0});
	
	std::cout << "System A (6x4)*x = b has no exact solution\n";
	std::cout << "QR finds least-squares solution minimizing ||A*x - b||_2\n";
	
	QRSolver<Real> qr_over(A_over);
	Vector<Real> x_ls = qr_over.LeastSquaresSolve(b_over);
	
	std::cout << "Least-squares solution x = "; x_ls.Print(std::cout, 10, 6); std::cout << std::endl;
	std::cout << "Residual norm ||A*x - b||_2 = " << (A_over*x_ls - b_over).NormL2() << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                       SVD - SINGULAR VALUE DECOMPOSITION                           ///
///////////////////////////////////////////////////////////////////////////////////////////

// Demo: SVD solver
void Docs_Demo_LinAlgSolvers_SVD()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: SVD (Singular Value Decomposition) Solver\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nSVDecompositionSolver - Most robust, handles singular matrices\n";
	std::cout << "Best for: Rank-deficient systems, singular matrices, pseudoinverse\n";
	
	// Regular system
	std::cout << "\n--- Regular System ---\n";
	Matrix<Real> A = TestBeds::mat_5x5();
	SVDecompositionSolver<Real> svd(A);
	
	Vector<Real> b = TestBeds::mat_5x5_rhs0();
	Vector<Real> x = svd.Solve(b);
	
	std::cout << "Solution x = "; x.Print(std::cout, 10, 6); std::cout << std::endl;
	std::cout << "Residual ||A*x - b|| = " << (A*x - b).NormL2() << std::endl;
	
	// Singular values
	std::cout << "\n--- Singular Values ---\n";
	Vector<Real> w = svd.getW();
	std::cout << "Singular values: "; w.Print(std::cout, 12, 6); std::cout << std::endl;
	
	Real cond = w[0] / w[w.size()-1];
	std::cout << "Condition number (sigma_max/sigma_min) = " << cond << std::endl;
	
	// Rank and nullspace
	std::cout << "\n--- Matrix Rank ---\n";
	int rank = svd.Rank();
	int nullity = svd.Nullity();
	std::cout << "Rank = " << rank << ", Nullity = " << nullity << std::endl;
	
	// Nearly singular matrix
	std::cout << "\n--- Nearly Singular Matrix ---\n";
	Matrix<Real> A_sing(3, 3);
	A_sing(0,0)=1; A_sing(0,1)=2; A_sing(0,2)=3;
	A_sing(1,0)=4; A_sing(1,1)=5; A_sing(1,2)=6;
	A_sing(2,0)=7; A_sing(2,1)=8; A_sing(2,2)=9.0001;  // Nearly rank-deficient
	
	SVDecompositionSolver<Real> svd_sing(A_sing);
	std::cout << "Singular values: "; svd_sing.getW().Print(std::cout, 15, 10); std::cout << std::endl;
	std::cout << "Note: Smallest singular value near 0, indicating near-singularity\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                       CHOLESKY DECOMPOSITION                                       ///
///////////////////////////////////////////////////////////////////////////////////////////

// Demo: Cholesky solver for symmetric positive definite matrices
void Docs_Demo_LinAlgSolvers_Cholesky()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Cholesky Decomposition Solver\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nCholeskySolver - For symmetric positive definite (SPD) matrices only\n";
	std::cout << "Best for: SPD matrices (2x faster than LU, guaranteed stability)\n";
	
	// Create a symmetric positive definite matrix
	// A = B^T * B guarantees SPD
	std::cout << "\n--- Creating SPD matrix: A = B^T*B ---\n";
	Matrix<Real> B(4, 4);
	B(0,0)=2; B(0,1)=1; B(0,2)=0; B(0,3)=0;
	B(1,0)=1; B(1,1)=3; B(1,2)=1; B(1,3)=0;
	B(2,0)=0; B(2,1)=1; B(2,2)=4; B(2,3)=1;
	B(3,0)=0; B(3,1)=0; B(3,2)=1; B(3,3)=5;
	Matrix<Real> A_spd = B.transpose() * B;
	
	std::cout << "A (SPD) =\n"; A_spd.Print(std::cout, 8, 3);
	
	// Solve
	CholeskySolver<Real> chol(A_spd);
	Vector<Real> b({1.0, 2.0, 3.0, 4.0});
	Vector<Real> x = chol.Solve(b);
	
	std::cout << "\nSolving A*x = b:\n";
	std::cout << "  b = "; b.Print(std::cout, 8, 3); std::cout << std::endl;
	std::cout << "  x = "; x.Print(std::cout, 10, 6); std::cout << std::endl;
	std::cout << "Residual ||A*x - b|| = " << (A_spd*x - b).NormL2() << std::endl;
	
	// Inverse and determinant
	std::cout << "\n--- Inverse and Determinant ---\n";
	Matrix<Real> A_inv(4, 4);
	chol.inverse(A_inv);
	Matrix<Real> I = A_spd * A_inv;
	Real err = 0;
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			err += std::abs(I(i,j) - (i==j ? 1.0 : 0.0));
	std::cout << "||A*A^-1 - I||_1 = " << err << std::endl;
	
	Real logdet = chol.logdet();
	std::cout << "log(det(A)) = " << logdet << " (det = " << std::exp(logdet) << ")\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                       ILL-CONDITIONED SYSTEM (HILBERT MATRIX)                      ///
///////////////////////////////////////////////////////////////////////////////////////////

// Demo: Compare solver accuracy on ill-conditioned Hilbert matrix
void Docs_Demo_LinAlgSolvers_IllConditioned()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Ill-Conditioned System (Hilbert Matrix)\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nHilbert matrix H[i,j] = 1/(i+j+1) is notoriously ill-conditioned\n";
	std::cout << "Demonstrates: QR and SVD provide much better accuracy than LU\n";
	
	// Create Hilbert matrix
	int n = 5;
	Matrix<Real> H(n, n);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			H[i][j] = 1.0 / (i + j + 1.0);
	
	std::cout << "\nHilbert matrix (5x5):\n"; H.Print(std::cout, 10, 6);
	
	// Known solution
	Vector<Real> x_exact({1.0, 2.0, 3.0, 4.0, 5.0});
	Vector<Real> b = H * x_exact;
	
	std::cout << "\nExact solution: "; x_exact.Print(std::cout, 8, 3); std::cout << std::endl;
	std::cout << "RHS b = H*x_exact: "; b.Print(std::cout, 10, 6); std::cout << std::endl;
	
	// Method 1: LU (may lose accuracy)
	std::cout << "\n--- LU Decomposition ---\n";
	LUSolver<Real> luSolver(H);
	Vector<Real> x_lu = luSolver.Solve(b);
	Real error_lu = (x_lu - x_exact).NormL2();
	std::cout << "LU solution: "; x_lu.Print(std::cout, 10, 6); std::cout << std::endl;
	std::cout << "Error ||x_LU - x_exact||_2 = " << error_lu << std::endl;
	
	// Method 2: QR (more stable)
	std::cout << "\n--- QR Decomposition ---\n";
	QRSolver<Real> qrSolver(H);
	Vector<Real> x_qr = qrSolver.Solve(b);
	Real error_qr = (x_qr - x_exact).NormL2();
	std::cout << "QR solution: "; x_qr.Print(std::cout, 10, 6); std::cout << std::endl;
	std::cout << "Error ||x_QR - x_exact||_2 = " << error_qr << std::endl;
	
	// Method 3: SVD (most stable)
	std::cout << "\n--- SVD Decomposition ---\n";
	SVDecompositionSolver<Real> svdSolver(H);
	Vector<Real> x_svd = svdSolver.Solve(b);
	Real error_svd = (x_svd - x_exact).NormL2();
	std::cout << "SVD solution: "; x_svd.Print(std::cout, 10, 6); std::cout << std::endl;
	std::cout << "Error ||x_SVD - x_exact||_2 = " << error_svd << std::endl;
	
	// Condition number
	Vector<Real> w = svdSolver.getW();
	Real cond = w[0] / w[w.size()-1];
	std::cout << "\nHilbert matrix condition number: " << cond << std::endl;
	
	// Summary
	std::cout << "\n--- Error Comparison ---\n";
	std::cout << "LU  error: " << error_lu << std::endl;
	std::cout << "QR  error: " << error_qr << "  (improvement: " << error_lu/error_qr << "x)\n";
	std::cout << "SVD error: " << error_svd << "  (improvement: " << error_lu/error_svd << "x)\n";
	std::cout << "\nConclusion: For ill-conditioned systems, QR and SVD provide\n";
	std::cout << "            significantly better accuracy (often 100-1000x better)!\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                       RANK-DEFICIENT SYSTEM WITH NULLSPACE                         ///
///////////////////////////////////////////////////////////////////////////////////////////

// Demo: Rank-deficient system with nullspace analysis
void Docs_Demo_LinAlgSolvers_RankDeficient()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Rank-Deficient System (SVD with Nullspace)\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nRank-deficient matrix: column 3 = 2*column 1 + column 2\n";
	std::cout << "Only SVD can handle this robustly!\n";
	
	// Rank-deficient matrix (column 3 = 2*col1 + col2)
	Matrix<Real> A(4, 3);
	A(0,0)=1.0; A(0,1)=2.0; A(0,2)=4.0;   // Row 1
	A(1,0)=2.0; A(1,1)=1.0; A(1,2)=5.0;   // Row 2
	A(2,0)=1.0; A(2,1)=3.0; A(2,2)=5.0;   // Row 3
	A(3,0)=3.0; A(3,1)=2.0; A(3,2)=8.0;   // Row 4
	
	std::cout << "\nMatrix A (4x3, rank-deficient):\n"; A.Print(std::cout, 8, 3);
	
	Vector<Real> b({1.0, 2.0, 3.0, 4.0});
	std::cout << "RHS b: "; b.Print(std::cout, 8, 3); std::cout << std::endl;
	
	// SVD can handle rank-deficient matrices
	SVDecompositionSolver<Real> svd(A);
	
	std::cout << "\n--- Singular Value Analysis ---\n";
	Vector<Real> w = svd.getW();
	std::cout << "Singular values: "; w.Print(std::cout, 15, 8); std::cout << std::endl;
	std::cout << "Rank: " << svd.Rank() << " (should be 2, not 3!)\n";
	std::cout << "Nullity: " << svd.Nullity() << " (dimension of null space)\n";
	
	// Solve with automatic thresholding
	std::cout << "\n--- Minimum-Norm Solution ---\n";
	Vector<Real> x = svd.Solve(b);
	std::cout << "SVD solution (minimum norm): "; x.Print(std::cout, 10, 6); std::cout << std::endl;
	
	// Verify residual
	Vector<Real> residual = A * x - b;
	std::cout << "Residual ||A*x - b||_2 = " << residual.NormL2() << std::endl;
	
	// Extract nullspace basis
	std::cout << "\n--- Null Space Basis ---\n";
	Matrix<Real> nullBasis = svd.Nullspace();
	std::cout << "Null space dimension: " << nullBasis.cols() << std::endl;
	
	if (nullBasis.cols() > 0) {
		std::cout << "Null space basis vectors:\n"; nullBasis.Print(std::cout, 10, 6);
		
		// Verify: A*v should be zero for null space vectors
		std::cout << "\nVerification A*v = 0 for null space vectors:\n";
		for (int i = 0; i < nullBasis.cols(); i++) {
			Vector<Real> v = nullBasis.VectorFromColumn(i);
			Vector<Real> Av = A * v;
			std::cout << "  ||A*v[" << i << "]||_2 = " << Av.NormL2() << " (should be ~0)\n";
		}
	}
	
	std::cout << "\nNote: LU/QR would fail with 'singular matrix' exception.\n";
	std::cout << "      Only SVD handles rank deficiency robustly!\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                       SOLVER COMPARISON                                            ///
///////////////////////////////////////////////////////////////////////////////////////////

// Demo: Compare solvers on the same system
void Docs_Demo_LinAlgSolvers_Comparison()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Solver Comparison\n";
	std::cout << "==========================================================================\n";
	
	Matrix<Real> A = TestBeds::mat_5x5();
	Vector<Real> b = TestBeds::mat_5x5_rhs0();
	
	std::cout << "Solving same 5x5 system with all methods:\n\n";
	
	// Gauss-Jordan
	Matrix<Real> A_gj(A);
	Vector<Real> b_gj(b);
	GaussJordanSolver<Real>::SolveInPlace(A_gj, b_gj);
	std::cout << "Gauss-Jordan residual: " << (A*b_gj - b).NormL2() << std::endl;
	
	// LU
	LUSolver<Real> lu(A);
	Vector<Real> x_lu = lu.Solve(b);
	std::cout << "LU residual:           " << (A*x_lu - b).NormL2() << std::endl;
	
	// QR
	QRSolver<Real> qr(A);
	Vector<Real> x_qr = qr.Solve(b);
	std::cout << "QR residual:           " << (A*x_qr - b).NormL2() << std::endl;
	
	// SVD
	SVDecompositionSolver<Real> svd(A);
	Vector<Real> x_svd = svd.Solve(b);
	std::cout << "SVD residual:          " << (A*x_svd - b).NormL2() << std::endl;
	
	std::cout << "\nAll methods produce essentially the same solution for well-conditioned systems.\n";
	std::cout << "Differences appear for ill-conditioned or singular matrices.\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                       MASTER DEMO FUNCTION                                         ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_LinAlgSolvers()
{
	std::cout << "\n\n";
	std::cout << "##########################################################################\n";
	std::cout << "#                LINEAR EQUATION SOLVER DEMOS                           #\n";
	std::cout << "##########################################################################\n";
	
	Docs_Demo_LinAlgSolvers_GaussJordan();
	Docs_Demo_LinAlgSolvers_LU();
	Docs_Demo_LinAlgSolvers_QR();
	Docs_Demo_LinAlgSolvers_SVD();
	Docs_Demo_LinAlgSolvers_Cholesky();
	Docs_Demo_LinAlgSolvers_IllConditioned();
	Docs_Demo_LinAlgSolvers_RankDeficient();
	Docs_Demo_LinAlgSolvers_Comparison();
}

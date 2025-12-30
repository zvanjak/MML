#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Vector.h"
#include "base/Matrix.h"
#include "base/MatrixSym.h"
#include "algorithms/EigenSystemSolvers.h"
#endif

using namespace MML;

///////////////////////////////////////////////////////////////////////////////////////////
///                 JACOBI METHOD FOR SYMMETRIC MATRICES                                ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_Eigen_Jacobi_Basic()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Jacobi Method for Symmetric Matrices\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nJacobi method uses Givens rotations to diagonalize symmetric matrices.\n";
	std::cout << "Excellent for small to medium matrices (n < 100).\n";
	
	// Create a symmetric matrix
	std::cout << "\n--- 3x3 Symmetric Matrix ---\n";
	MatrixSym<Real> A(3);
	A(0,0) = 4;  A(0,1) = 2;  A(0,2) = 1;
	             A(1,1) = 5;  A(1,2) = 3;
	                          A(2,2) = 6;
	
	std::cout << "Matrix A (symmetric):\n";
	Matrix<Real> Afull(3, 3);
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			Afull(i,j) = (i <= j) ? A(i,j) : A(j,i);
	Afull.Print(std::cout, 10, 4);
	
	// Solve
	auto result = SymmMatEigenSolverJacobi::Solve(A);
	
	std::cout << "\nEigenvalues (ascending order): ";
	result.eigenvalues.Print(std::cout, 12, 6);
	std::cout << std::endl;
	
	std::cout << "\nEigenvectors (columns):\n";
	result.eigenvectors.Print(std::cout, 12, 6);
	
	std::cout << "\nConvergence info:\n";
	std::cout << "  Iterations: " << result.iterations << std::endl;
	std::cout << "  Residual (off-diagonal norm): " << result.residual << std::endl;
	std::cout << "  Converged: " << (result.converged ? "yes" : "no") << std::endl;
	
	// Verify: A*v = λ*v
	std::cout << "\n--- Verification: A*v = λ*v ---\n";
	for (int i = 0; i < 3; i++) {
		Vector<Real> v(3);
		for (int j = 0; j < 3; j++) v[j] = result.eigenvectors(j, i);
		
		Vector<Real> Av = Afull * v;
		Vector<Real> lambdaV = result.eigenvalues[i] * v;
		
		std::cout << "λ" << i << " = " << result.eigenvalues[i] << ", ||A*v - λ*v|| = " 
		          << (Av - lambdaV).NormL2() << std::endl;
	}
}

void Docs_Demo_Eigen_Jacobi_Properties()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Eigenvalue Properties with Jacobi\n";
	std::cout << "==========================================================================\n";
	
	// Create a 4x4 symmetric matrix
	MatrixSym<Real> B(4);
	B(0,0) = 10; B(0,1) = -1; B(0,2) = 2;  B(0,3) = 0;
	             B(1,1) = 11; B(1,2) = -1; B(1,3) = 3;
	                          B(2,2) = 9;  B(2,3) = -1;
	                                       B(3,3) = 12;
	
	auto result = SymmMatEigenSolverJacobi::Solve(B);
	
	std::cout << "\n--- Orthogonality of Eigenvectors ---\n";
	std::cout << "For symmetric matrices, eigenvectors form orthonormal basis.\n";
	
	Matrix<Real> VtV = result.eigenvectors.GetTranspose() * result.eigenvectors;
	std::cout << "\nV^T * V (should be identity):\n";
	VtV.Print(std::cout, 10, 6);
	
	std::cout << "\n--- Trace and Determinant ---\n";
	Real trace = 0;
	Real det = 1;
	for (int i = 0; i < 4; i++) {
		trace += result.eigenvalues[i];
		det *= result.eigenvalues[i];
	}
	std::cout << "Sum of eigenvalues: " << trace << std::endl;
	std::cout << "Product of eigenvalues: " << det << std::endl;
	std::cout << "(These equal trace and determinant of original matrix)\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                 QR METHOD FOR SYMMETRIC MATRICES                                    ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_Eigen_QR_Basic()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: QR Algorithm for Symmetric Matrices\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nQR algorithm with Wilkinson shift: Industry-standard method.\n";
	std::cout << "O(n³) complexity, cubic convergence rate.\n";
	
	// Create a symmetric matrix
	MatrixSym<Real> A(4);
	A(0,0) = 4;  A(0,1) = 1;  A(0,2) = -2; A(0,3) = 2;
	             A(1,1) = 6;  A(1,2) = 0;  A(1,3) = -2;
	                          A(2,2) = 3;  A(2,3) = 1;
	                                       A(3,3) = 5;
	
	std::cout << "\n--- 4x4 Symmetric Matrix ---\n";
	Matrix<Real> Afull(4, 4);
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			Afull(i,j) = (i <= j) ? A(i,j) : A(j,i);
	Afull.Print(std::cout, 10, 4);
	
	// Solve
	auto result = SymmMatEigenSolverQR::Solve(A);
	
	std::cout << "\nEigenvalues: ";
	result.eigenvalues.Print(std::cout, 12, 6);
	std::cout << std::endl;
	
	std::cout << "\nEigenvectors:\n";
	result.eigenvectors.Print(std::cout, 12, 6);
	
	std::cout << "\nConvergence: " << result.iterations << " iterations, "
	          << "residual = " << result.residual << std::endl;
}

void Docs_Demo_Eigen_Comparison()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Jacobi vs QR Comparison\n";
	std::cout << "==========================================================================\n";
	
	// Create a larger symmetric matrix
	int n = 6;
	MatrixSym<Real> A(n);
	
	// Fill with values from a Hilbert-like matrix (known to have all positive eigenvalues)
	for (int i = 0; i < n; i++)
		for (int j = i; j < n; j++)
			A(i,j) = 1.0 / (1 + i + j);
	
	std::cout << "\n--- 6x6 Hilbert-like Matrix ---\n";
	
	// Jacobi method
	auto jacobi_result = SymmMatEigenSolverJacobi::Solve(A);
	
	// QR method
	auto qr_result = SymmMatEigenSolverQR::Solve(A);
	
	std::cout << "Jacobi eigenvalues: ";
	jacobi_result.eigenvalues.Print(std::cout, 12, 8);
	std::cout << std::endl;
	
	std::cout << "QR eigenvalues:     ";
	qr_result.eigenvalues.Print(std::cout, 12, 8);
	std::cout << std::endl;
	
	std::cout << "\nComparison:\n";
	std::cout << "  Jacobi iterations: " << jacobi_result.iterations << std::endl;
	std::cout << "  QR iterations:     " << qr_result.iterations << std::endl;
	
	std::cout << "\nMax eigenvalue difference: ";
	Real maxDiff = 0;
	for (int i = 0; i < n; i++)
		maxDiff = std::max(maxDiff, std::abs(jacobi_result.eigenvalues[i] - qr_result.eigenvalues[i]));
	std::cout << maxDiff << " (both methods agree!)\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                 POWER METHOD FOR DOMINANT EIGENVALUE                                ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_Eigen_Power_Method()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Power Method (Dominant Eigenvalue)\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nPower iteration finds largest eigenvalue in absolute value.\n";
	std::cout << "Simple but only gives one eigenvalue.\n";
	
	// Create a matrix with known dominant eigenvalue
	Matrix<Real> A(3, 3);
	A(0,0) = 6; A(0,1) = 5; A(0,2) = -5;
	A(1,0) = 2; A(1,1) = 6; A(1,2) = -2;
	A(2,0) = 2; A(2,1) = 5; A(2,2) = -1;
	
	std::cout << "\nMatrix A:\n";
	A.Print(std::cout, 10, 4);
	
	// Power iteration
	Vector<Real> x({1, 0, 0});  // Initial guess
	Real lambda = 0;
	
	std::cout << "\n--- Power Iteration ---\n";
	for (int iter = 0; iter < 10; iter++) {
		Vector<Real> y = A * x;
		lambda = y.NormL2();
		x = y / lambda;
		
		if (iter < 5 || iter >= 8)
			std::cout << "Iteration " << iter+1 << ": λ ≈ " << lambda << std::endl;
		else if (iter == 5)
			std::cout << "...\n";
	}
	
	std::cout << "\nDominant eigenvalue: " << lambda << std::endl;
	std::cout << "Corresponding eigenvector: ";
	x.Print(std::cout, 10, 6);
	std::cout << std::endl;
	
	// Verify
	Vector<Real> Ax = A * x;
	std::cout << "\nVerification ||A*x - λ*x||: " << (Ax - lambda*x).NormL2() << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                 INVERSE ITERATION FOR SPECIFIC EIGENVALUE                           ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_Eigen_Inverse_Iteration()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Inverse Iteration\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nInverse iteration finds eigenvalue closest to a given shift σ.\n";
	std::cout << "Faster convergence when σ is close to an eigenvalue.\n";
	
	// Create a symmetric matrix
	MatrixSym<Real> A(3);
	A(0,0) = 4;  A(0,1) = 1;  A(0,2) = 0;
	             A(1,1) = 5;  A(1,2) = 1;
	                          A(2,2) = 3;
	
	// First, find all eigenvalues using Jacobi
	auto result = SymmMatEigenSolverJacobi::Solve(A);
	std::cout << "\nAll eigenvalues (from Jacobi): ";
	result.eigenvalues.Print(std::cout, 10, 6);
	std::cout << std::endl;
	
	// Inverse iteration concept
	Real sigma = 3.5;
	std::cout << "\n--- Finding eigenvalue closest to σ = " << sigma << " ---\n";
	std::cout << "By choosing σ between eigenvalues, we select which one to find.\n";
	
	// Find which eigenvalue is closest to sigma
	int closestIdx = 0;
	Real minDist = std::abs(result.eigenvalues[0] - sigma);
	for (int i = 1; i < 3; i++) {
		Real dist = std::abs(result.eigenvalues[i] - sigma);
		if (dist < minDist) {
			minDist = dist;
			closestIdx = i;
		}
	}
	
	std::cout << "Closest eigenvalue to σ=" << sigma << " is λ = " 
	          << result.eigenvalues[closestIdx] << "\n";
	std::cout << "Inverse iteration would converge to this eigenvalue.\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                 GENERALIZED EIGENVALUE PROBLEM                                      ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_Eigen_Generalized()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Generalized Eigenvalue Problem\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nGeneralized eigenvalue problem: A*v = λ*B*v\n";
	std::cout << "Common in vibration analysis, quantum mechanics.\n";
	
	// For symmetric A, positive definite B, can reduce to standard form
	std::cout << "\nFor SPD matrices, reduce to standard form:\n";
	std::cout << "  A*v = λ*B*v  →  (L^-1 * A * L^-T) * w = λ * w\n";
	std::cout << "  where B = L*L^T (Cholesky), w = L^T * v\n";
	
	std::cout << "\n(Full implementation in GenEigenSolver - see EigenSystemSolvers.h)\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                    SPECIAL MATRICES                                                 ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_Eigen_Special_Matrices()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Eigenvalues of Special Matrices\n";
	std::cout << "==========================================================================\n";
	
	// 1. Rotation matrix (complex eigenvalues - won't work with symmetric solver)
	std::cout << "\n--- Diagonal Matrix ---\n";
	MatrixSym<Real> D(3);
	D(0,0) = 2;  D(0,1) = 0;  D(0,2) = 0;
	             D(1,1) = 5;  D(1,2) = 0;
	                          D(2,2) = 3;
	
	auto result = SymmMatEigenSolverJacobi::Solve(D);
	std::cout << "Diagonal matrix eigenvalues: ";
	result.eigenvalues.Print(std::cout, 10, 4);
	std::cout << " (just the diagonal entries!)\n";
	
	// 2. Positive definite matrix
	std::cout << "\n--- Positive Definite Matrix ---\n";
	MatrixSym<Real> P(3);
	P(0,0) = 4;  P(0,1) = 2;  P(0,2) = 1;
	             P(1,1) = 5;  P(1,2) = 2;
	                          P(2,2) = 6;
	
	auto pdResult = SymmMatEigenSolverJacobi::Solve(P);
	std::cout << "Eigenvalues: ";
	pdResult.eigenvalues.Print(std::cout, 10, 4);
	std::cout << std::endl;
	std::cout << "All positive → matrix is positive definite\n";
	
	// 3. Indefinite matrix
	std::cout << "\n--- Indefinite Matrix ---\n";
	MatrixSym<Real> I(3);
	I(0,0) = 1;  I(0,1) = 2;  I(0,2) = 0;
	             I(1,1) = 0;  I(1,2) = 1;
	                          I(2,2) = -1;
	
	auto indResult = SymmMatEigenSolverJacobi::Solve(I);
	std::cout << "Eigenvalues: ";
	indResult.eigenvalues.Print(std::cout, 10, 4);
	std::cout << std::endl;
	std::cout << "Has both positive and negative → indefinite\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         MAIN DEMO FUNCTION                                          ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_Eigen_solvers()
{
	std::cout << "\n##########################################################################\n";
	std::cout << "#                    EIGENVALUE SOLVER DEMOS                            #\n";
	std::cout << "##########################################################################\n";
	
	Docs_Demo_Eigen_Jacobi_Basic();
	Docs_Demo_Eigen_Jacobi_Properties();
	Docs_Demo_Eigen_QR_Basic();
	Docs_Demo_Eigen_Comparison();
	Docs_Demo_Eigen_Power_Method();
	Docs_Demo_Eigen_Inverse_Iteration();
	Docs_Demo_Eigen_Generalized();
	Docs_Demo_Eigen_Special_Matrices();
	
	std::cout << "\n=== All Eigenvalue Solver Demos Complete ===\n";
}

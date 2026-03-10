///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        linear_system_tests.cpp                                             ///
///  Description: Comprehensive tests for LinearSystem unified facade                 ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                    ///
///////////////////////////////////////////////////////////////////////////////////////////

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "MMLBase.h"
#include "systems/LinearSystem.h"

using namespace MML;
using namespace MML::Systems;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace MML::Tests::Systems::LinearSystemTests {

// Helper to create Hilbert matrix (notoriously ill-conditioned)
Matrix<Real> CreateHilbertMatrix(int n)
{
	Matrix<Real> H(n, n);
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			H(i, j) = 1.0 / (i + j + 1);
	return H;
}

// Helper to create SPD matrix
Matrix<Real> CreateSPDMatrix(int n)
{
	Matrix<Real> A(n, n);
	// Create A = diag dominant + small off-diagonal
	for (int i = 0; i < n; ++i)
	{
		A(i, i) = 10.0 + i;  // Large diagonal
		for (int j = 0; j < n; ++j)
			if (i != j)
				A(i, j) = 0.5;  // Small symmetric off-diagonal
	}
	return A;
}

//=============================================================================
// CONSTRUCTOR TESTS
//=============================================================================

TEST_CASE("LinearSystem - Constructors", "[LinearSystem][Constructor]")
{
	SECTION("Matrix + Vector constructor")
	{
		Matrix<Real> A(3, 3, {1, 2, 3, 4, 5, 6, 7, 8, 10});
		Vector<Real> b({1, 2, 3});
		
		LinearSystem<Real> sys(A, b);
		
		REQUIRE(sys.Rows() == 3);
		REQUIRE(sys.Cols() == 3);
		REQUIRE(sys.IsSquare() == true);
	}
	
	SECTION("Matrix + Matrix constructor (multiple RHS)")
	{
		Matrix<Real> A(3, 3, {2, 1, 0, 1, 3, 1, 0, 1, 2});
		Matrix<Real> B(3, 2, {1, 2, 3, 4, 5, 6});  // 2 right-hand sides
		
		LinearSystem<Real> sys(A, B);
		
		REQUIRE(sys.Rows() == 3);
		REQUIRE(sys.Cols() == 3);
	}
	
	SECTION("Matrix only constructor (analysis)")
	{
		Matrix<Real> A(4, 3, {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12});
		
		LinearSystem<Real> sys(A);
		
		REQUIRE(sys.Rows() == 4);
		REQUIRE(sys.Cols() == 3);
		REQUIRE(sys.IsOverdetermined() == true);
	}
	
	SECTION("Dimension mismatch throws")
	{
		Matrix<Real> A(3, 3, {1, 0, 0, 0, 1, 0, 0, 0, 1});
		Vector<Real> b({1, 2});  // Wrong size
		
		REQUIRE_THROWS_AS(LinearSystem<Real>(A, b), MatrixDimensionError);
	}
}

//=============================================================================
// BASIC SOLVING TESTS
//=============================================================================

TEST_CASE("LinearSystem - Auto Solve", "[LinearSystem][Solve]")
{
	SECTION("Simple 2x2 system")
	{
		Matrix<Real> A(2, 2, {2, 1, 1, 3});
		Vector<Real> b({5, 5});  // Solution: x = [2, 1]
		
		LinearSystem<Real> sys(A, b);
		Vector<Real> x = sys.Solve();
		
		REQUIRE_THAT(x[0], WithinAbs(2.0, 1e-10));
		REQUIRE_THAT(x[1], WithinAbs(1.0, 1e-10));
	}
	
	SECTION("3x3 system")
	{
		Matrix<Real> A(3, 3, {1, 2, 3, 4, 5, 6, 7, 8, 10});
		Vector<Real> b({14, 32, 53});  // Solution: x = [1, 2, 3]
		
		LinearSystem<Real> sys(A, b);
		Vector<Real> x = sys.Solve();
		
		REQUIRE_THAT(x[0], WithinAbs(1.0, 1e-9));
		REQUIRE_THAT(x[1], WithinAbs(2.0, 1e-9));
		REQUIRE_THAT(x[2], WithinAbs(3.0, 1e-9));
	}
}

TEST_CASE("LinearSystem - Specific Solvers", "[LinearSystem][Solvers]")
{
	Matrix<Real> A(3, 3, {4, 1, 0, 1, 4, 1, 0, 1, 4});  // Tridiagonal SPD
	Vector<Real> b({5, 6, 5});
	LinearSystem<Real> sys(A, b);
	
	SECTION("Gauss-Jordan")
	{
		Vector<Real> x = sys.SolveByGaussJordan();
		REQUIRE_THAT(sys.RelativeResidual(x), WithinAbs(0.0, 1e-12));
	}
	
	SECTION("LU Decomposition")
	{
		Vector<Real> x = sys.SolveByLU();
		REQUIRE_THAT(sys.RelativeResidual(x), WithinAbs(0.0, 1e-12));
	}
	
	SECTION("Cholesky (SPD)")
	{
		Vector<Real> x = sys.SolveByCholesky();
		REQUIRE_THAT(sys.RelativeResidual(x), WithinAbs(0.0, 1e-12));
	}
	
	SECTION("QR Decomposition")
	{
		Vector<Real> x = sys.SolveByQR();
		REQUIRE_THAT(sys.RelativeResidual(x), WithinAbs(0.0, 1e-10));
	}
	
	SECTION("SVD Decomposition")
	{
		Vector<Real> x = sys.SolveBySVD();
		REQUIRE_THAT(sys.RelativeResidual(x), WithinAbs(0.0, 1e-10));
	}
	
	SECTION("All solvers give same result")
	{
		Vector<Real> xGJ = sys.SolveByGaussJordan();
		Vector<Real> xLU = sys.SolveByLU();
		Vector<Real> xCh = sys.SolveByCholesky();
		Vector<Real> xQR = sys.SolveByQR();
		Vector<Real> xSVD = sys.SolveBySVD();
		
		for (int i = 0; i < 3; ++i)
		{
			REQUIRE_THAT(xLU[i], WithinAbs(xGJ[i], 1e-10));
			REQUIRE_THAT(xCh[i], WithinAbs(xGJ[i], 1e-10));
			REQUIRE_THAT(xQR[i], WithinAbs(xGJ[i], 1e-8));
			REQUIRE_THAT(xSVD[i], WithinAbs(xGJ[i], 1e-8));
		}
	}
}

TEST_CASE("LinearSystem - Iterative Solvers", "[LinearSystem][Iterative]")
{
	// Diagonally dominant matrix (convergence guaranteed)
	Matrix<Real> A(4, 4, {
		10, 1, 1, 1,
		1, 10, 1, 1,
		1, 1, 10, 1,
		1, 1, 1, 10
	});
	Vector<Real> b({13, 13, 13, 13});  // Solution: x = [1, 1, 1, 1]
	
	LinearSystem<Real> sys(A, b);
	
	SECTION("Jacobi iteration")
	{
		Vector<Real> x = sys.SolveIterative(IterativeMethod::Jacobi, 1e-10, 1000);
		REQUIRE_THAT(sys.RelativeResidual(x), WithinAbs(0.0, 1e-9));
	}
	
	SECTION("Gauss-Seidel iteration")
	{
		Vector<Real> x = sys.SolveIterative(IterativeMethod::GaussSeidel, 1e-10, 1000);
		REQUIRE_THAT(sys.RelativeResidual(x), WithinAbs(0.0, 1e-9));
	}
	
	SECTION("SOR iteration")
	{
		Vector<Real> x = sys.SolveIterative(IterativeMethod::SOR, 1e-10, 1000);
		REQUIRE_THAT(sys.RelativeResidual(x), WithinAbs(0.0, 1e-9));
	}
	
	SECTION("Auto-select iterative")
	{
		Vector<Real> x = sys.SolveIterative(IterativeMethod::Auto, 1e-10, 1000);
		REQUIRE_THAT(sys.RelativeResidual(x), WithinAbs(0.0, 1e-9));
	}
}

//=============================================================================
// OVERDETERMINED SYSTEMS (LEAST SQUARES)
//=============================================================================

TEST_CASE("LinearSystem - Overdetermined Systems", "[LinearSystem][LeastSquares]")
{
	SECTION("Simple least squares")
	{
		// Fit y = mx + c to points (0,1), (1,2), (2,4)
		// System: c = 1, m+c = 2, 2m+c = 4
		Matrix<Real> A(3, 2, {1, 0, 1, 1, 1, 2});
		Vector<Real> b({1, 2, 4});
		
		LinearSystem<Real> sys(A, b);
		REQUIRE(sys.IsOverdetermined() == true);
		
		Vector<Real> x = sys.SolveLeastSquares();
		
		// Least squares solution via A^T*A*x = A^T*b:
		// A^T*A = [[3,3],[3,5]], A^T*b = [7,10]
		// Solution: c = 5/6 ≈ 0.8333, m = 1.5
		REQUIRE_THAT(x[0], WithinAbs(0.8333, 0.01));  // c = 5/6
		REQUIRE_THAT(x[1], WithinAbs(1.5, 0.01));     // m
	}
	
	SECTION("Overdetermined with exact solution")
	{
		// Redundant equations with exact solution
		Matrix<Real> A(4, 2, {1, 0, 0, 1, 1, 0, 0, 1});
		Vector<Real> b({2, 3, 2, 3});
		
		LinearSystem<Real> sys(A, b);
		Vector<Real> x = sys.SolveLeastSquares();
		
		REQUIRE_THAT(x[0], WithinAbs(2.0, 1e-10));
		REQUIRE_THAT(x[1], WithinAbs(3.0, 1e-10));
	}
}

//=============================================================================
// SOLUTION VERIFICATION
//=============================================================================

TEST_CASE("LinearSystem - Verification", "[LinearSystem][Verify]")
{
	Matrix<Real> A(3, 3, {4, 1, 0, 1, 4, 1, 0, 1, 4});
	Vector<Real> b({5, 6, 5});
	LinearSystem<Real> sys(A, b);
	
	SECTION("Verify correct solution")
	{
		Vector<Real> x = sys.Solve();
		auto verify = sys.Verify(x);
		
		REQUIRE_THAT(verify.absoluteResidual, WithinAbs(0.0, 1e-12));
		REQUIRE_THAT(verify.relativeResidual, WithinAbs(0.0, 1e-12));
		REQUIRE(verify.isAccurate == true);
	}
	
	SECTION("Verify incorrect solution")
	{
		// x = {1,1,1} is actually the correct solution for this system!
		// Use a truly wrong solution
		Vector<Real> wrong({0, 0, 0});  // Zero vector is wrong
		auto verify = sys.Verify(wrong);
		
		// ||A*0 - b|| = ||b|| = sqrt(5^2 + 6^2 + 5^2) = sqrt(86) ≈ 9.27
		REQUIRE(verify.absoluteResidual > 5.0);
		REQUIRE(verify.isAccurate == false);
	}
	
	SECTION("Residual norm")
	{
		Vector<Real> x = sys.Solve();
		Real residual = sys.ResidualNorm(x);
		
		REQUIRE_THAT(residual, WithinAbs(0.0, 1e-12));
	}
	
	SECTION("Relative residual")
	{
		Vector<Real> x = sys.Solve();
		Real relResidual = sys.RelativeResidual(x);
		
		REQUIRE_THAT(relResidual, WithinAbs(0.0, 1e-12));
	}
}

//=============================================================================
// MATRIX PROPERTY TESTS
//=============================================================================

TEST_CASE("LinearSystem - Dimension Properties", "[LinearSystem][Properties]")
{
	SECTION("Square matrix")
	{
		Matrix<Real> A(3, 3, {1, 0, 0, 0, 1, 0, 0, 0, 1});
		LinearSystem<Real> sys(A);
		
		REQUIRE(sys.IsSquare() == true);
		REQUIRE(sys.IsOverdetermined() == false);
		REQUIRE(sys.IsUnderdetermined() == false);
	}
	
	SECTION("Overdetermined (tall)")
	{
		Matrix<Real> A(5, 3);
		LinearSystem<Real> sys(A);
		
		REQUIRE(sys.IsSquare() == false);
		REQUIRE(sys.IsOverdetermined() == true);
		REQUIRE(sys.IsUnderdetermined() == false);
	}
	
	SECTION("Underdetermined (wide)")
	{
		Matrix<Real> A(3, 5);
		LinearSystem<Real> sys(A);
		
		REQUIRE(sys.IsSquare() == false);
		REQUIRE(sys.IsOverdetermined() == false);
		REQUIRE(sys.IsUnderdetermined() == true);
	}
}

TEST_CASE("LinearSystem - Structure Detection", "[LinearSystem][Structure]")
{
	SECTION("Symmetric matrix")
	{
		Matrix<Real> A(3, 3, {2, 1, 0, 1, 3, 1, 0, 1, 2});
		LinearSystem<Real> sys(A);
		
		REQUIRE(sys.IsSymmetric() == true);
	}
	
	SECTION("Non-symmetric matrix")
	{
		Matrix<Real> A(3, 3, {2, 1, 0, 0, 3, 1, 0, 0, 2});
		LinearSystem<Real> sys(A);
		
		REQUIRE(sys.IsSymmetric() == false);
	}
	
	SECTION("Positive definite")
	{
		Matrix<Real> A = CreateSPDMatrix(4);
		LinearSystem<Real> sys(A);
		
		REQUIRE(sys.IsSymmetric() == true);
		REQUIRE(sys.IsPositiveDefinite() == true);
	}
	
	SECTION("Diagonally dominant")
	{
		Matrix<Real> A(3, 3, {10, 1, 2, 1, 10, 2, 1, 2, 10});
		LinearSystem<Real> sys(A);
		
		REQUIRE(sys.IsDiagonallyDominant() == true);
	}
	
	SECTION("Upper triangular")
	{
		Matrix<Real> A(3, 3, {1, 2, 3, 0, 4, 5, 0, 0, 6});
		LinearSystem<Real> sys(A);
		
		REQUIRE(sys.IsUpperTriangular() == true);
		REQUIRE(sys.IsLowerTriangular() == false);
	}
	
	SECTION("Lower triangular")
	{
		Matrix<Real> A(3, 3, {1, 0, 0, 2, 3, 0, 4, 5, 6});
		LinearSystem<Real> sys(A);
		
		REQUIRE(sys.IsUpperTriangular() == false);
		REQUIRE(sys.IsLowerTriangular() == true);
	}
	
	SECTION("Diagonal")
	{
		Matrix<Real> A(3, 3, {1, 0, 0, 0, 2, 0, 0, 0, 3});
		LinearSystem<Real> sys(A);
		
		REQUIRE(sys.IsDiagonal() == true);
	}
}

//=============================================================================
// NUMERICAL PROPERTIES
//=============================================================================

TEST_CASE("LinearSystem - Numerical Properties", "[LinearSystem][Numerical]")
{
	SECTION("Determinant")
	{
		Matrix<Real> A(3, 3, {1, 2, 3, 0, 4, 5, 0, 0, 6});
		LinearSystem<Real> sys(A);
		
		Real det = sys.Determinant();
		REQUIRE_THAT(det, WithinAbs(24.0, 1e-10));  // 1*4*6 = 24
	}
	
	SECTION("Rank of full rank matrix")
	{
		Matrix<Real> A(3, 3, {1, 0, 0, 0, 1, 0, 0, 0, 1});
		LinearSystem<Real> sys(A);
		
		REQUIRE(sys.Rank() == 3);
		REQUIRE(sys.Nullity() == 0);
	}
	
	SECTION("Rank of rank-deficient matrix")
	{
		Matrix<Real> A(3, 3, {1, 2, 3, 2, 4, 6, 3, 6, 9});  // Rows proportional
		LinearSystem<Real> sys(A);
		
		REQUIRE(sys.Rank() == 1);
		REQUIRE(sys.Nullity() == 2);
	}
	
	SECTION("Condition number of identity")
	{
		Matrix<Real> I(3, 3, {1, 0, 0, 0, 1, 0, 0, 0, 1});
		LinearSystem<Real> sys(I);
		
		Real cond = sys.ConditionNumber();
		REQUIRE_THAT(cond, WithinAbs(1.0, 1e-10));
	}
	
	SECTION("Ill-conditioned Hilbert matrix")
	{
		Matrix<Real> H = CreateHilbertMatrix(5);
		LinearSystem<Real> sys(H);
		
		Real cond = sys.ConditionNumber();
		REQUIRE(cond > 1e4);  // Hilbert(5) has cond ~4.8e5
	}
	
	SECTION("Stability assessment")
	{
		// Well-conditioned (identity)
		Matrix<Real> I(3, 3, {1, 0, 0, 0, 1, 0, 0, 0, 1});
		LinearSystem<Real> sysGood(I);
		REQUIRE(sysGood.AssessStability() == MatrixStability::WellConditioned);
		
		// Ill-conditioned (Hilbert)
		Matrix<Real> H = CreateHilbertMatrix(10);
		LinearSystem<Real> sysBad(H);
		REQUIRE(sysBad.AssessStability() != MatrixStability::WellConditioned);
	}
	
	SECTION("Expected digits lost")
	{
		Matrix<Real> A(3, 3, {1, 0, 0, 0, 1, 0, 0, 0, 1});
		LinearSystem<Real> sys(A);
		
		int digitsLost = sys.ExpectedDigitsLost();
		REQUIRE(digitsLost == 0);
	}
}

//=============================================================================
// MATRIX OPERATIONS
//=============================================================================

TEST_CASE("LinearSystem - Matrix Operations", "[LinearSystem][Operations]")
{
	SECTION("Matrix inverse")
	{
		Matrix<Real> A(3, 3, {1, 2, 0, 0, 1, 2, 2, 0, 1});
		LinearSystem<Real> sys(A);
		
		Matrix<Real> Ainv = sys.Inverse();
		
		// A * Ainv should be identity
		Matrix<Real> I = A * Ainv;
		for (int i = 0; i < 3; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				Real expected = (i == j) ? 1.0 : 0.0;
				REQUIRE_THAT(I(i, j), WithinAbs(expected, 1e-10));
			}
		}
	}
	
	SECTION("Pseudo-inverse")
	{
		Matrix<Real> A(3, 2, {1, 0, 0, 1, 1, 1});
		LinearSystem<Real> sys(A);
		
		Matrix<Real> Apinv = sys.PseudoInverse();
		
		// A * A+ * A should equal A
		Matrix<Real> result = A * Apinv * A;
		for (int i = 0; i < A.rows(); ++i)
			for (int j = 0; j < A.cols(); ++j)
				REQUIRE_THAT(result(i, j), WithinAbs(A(i, j), 1e-10));
	}
	
	SECTION("Null space")
	{
		// Rank-deficient matrix
		Matrix<Real> A(2, 3, {1, 2, 3, 2, 4, 6});  // Rows proportional
		LinearSystem<Real> sys(A);
		
		Matrix<Real> N = sys.NullSpace();
		
		// Verify A * N = 0
		if (N.cols() > 0)
		{
			Matrix<Real> AN = A * N;
			for (int i = 0; i < AN.rows(); ++i)
				for (int j = 0; j < AN.cols(); ++j)
					REQUIRE_THAT(AN(i, j), WithinAbs(0.0, 1e-10));
		}
	}
}

//=============================================================================
// COMPREHENSIVE ANALYSIS
//=============================================================================

TEST_CASE("LinearSystem - Comprehensive Analysis", "[LinearSystem][Analysis]")
{
	SECTION("Well-conditioned SPD system")
	{
		Matrix<Real> A = CreateSPDMatrix(4);
		LinearSystem<Real> sys(A);
		
		auto analysis = sys.Analyze();
		
		REQUIRE(analysis.rows == 4);
		REQUIRE(analysis.cols == 4);
		REQUIRE(analysis.isSquare == true);
		REQUIRE(analysis.isSymmetric == true);
		REQUIRE(analysis.isPositiveDefinite == true);
		REQUIRE(analysis.rank == 4);
		REQUIRE(analysis.hasUniqueSolution == true);
		REQUIRE(analysis.recommendedSolver == "Cholesky");
	}
	
	SECTION("Triangular system")
	{
		Matrix<Real> A(3, 3, {1, 2, 3, 0, 4, 5, 0, 0, 6});
		LinearSystem<Real> sys(A);
		
		auto analysis = sys.Analyze();
		
		REQUIRE(analysis.isUpperTriangular == true);
		REQUIRE(analysis.recommendedSolver == "Triangular");
	}
	
	SECTION("Overdetermined system")
	{
		Matrix<Real> A(5, 3);
		for (int i = 0; i < 5; ++i)
			for (int j = 0; j < 3; ++j)
				A(i, j) = i + j + 1;
		
		LinearSystem<Real> sys(A);
		
		auto analysis = sys.Analyze();
		
		REQUIRE(analysis.isOverdetermined == true);
		REQUIRE(analysis.recommendedSolver == "QR");
	}
	
	SECTION("Analysis report generated")
	{
		Matrix<Real> A = CreateSPDMatrix(3);
		LinearSystem<Real> sys(A);
		
		auto analysis = sys.Analyze();
		
		REQUIRE(analysis.analysisReport.size() > 0);
		REQUIRE(analysis.analysisReport.find("Matrix:") != std::string::npos);
	}
}

//=============================================================================
// SPECIAL MATRICES
//=============================================================================

TEST_CASE("LinearSystem - Triangular Solve", "[LinearSystem][Triangular]")
{
	SECTION("Upper triangular - back substitution")
	{
		Matrix<Real> A(3, 3, {2, 1, 1, 0, 3, 1, 0, 0, 4});
		Vector<Real> b({4, 4, 4});  // Solution: x = [1, 1, 1]
		
		LinearSystem<Real> sys(A, b);
		Vector<Real> x = sys.Solve();  // Should auto-select triangular
		
		REQUIRE_THAT(x[0], WithinAbs(1.0, 1e-10));
		REQUIRE_THAT(x[1], WithinAbs(1.0, 1e-10));
		REQUIRE_THAT(x[2], WithinAbs(1.0, 1e-10));
	}
	
	SECTION("Lower triangular - forward substitution")
	{
		Matrix<Real> A(3, 3, {2, 0, 0, 1, 3, 0, 1, 1, 4});
		Vector<Real> b({2, 4, 6});  // Solution: x = [1, 1, 1]
		
		LinearSystem<Real> sys(A, b);
		Vector<Real> x = sys.Solve();
		
		REQUIRE_THAT(x[0], WithinAbs(1.0, 1e-10));
		REQUIRE_THAT(x[1], WithinAbs(1.0, 1e-10));
		REQUIRE_THAT(x[2], WithinAbs(1.0, 1e-10));
	}
}

//=============================================================================
// CONVENIENCE FUNCTIONS
//=============================================================================

TEST_CASE("LinearSystem - Convenience Functions", "[LinearSystem][Convenience]")
{
	SECTION("SolveLinearSystem")
	{
		Matrix<Real> A(2, 2, {2, 1, 1, 3});
		Vector<Real> b({5, 5});
		
		Vector<Real> x = SolveLinearSystem(A, b);
		
		REQUIRE_THAT(x[0], WithinAbs(2.0, 1e-10));
		REQUIRE_THAT(x[1], WithinAbs(1.0, 1e-10));
	}
	
	SECTION("SolveLeastSquares convenience")
	{
		Matrix<Real> A(4, 2, {1, 0, 0, 1, 1, 0, 0, 1});
		Vector<Real> b({2, 3, 2, 3});
		
		Vector<Real> x = SolveLeastSquares(A, b);
		
		REQUIRE_THAT(x[0], WithinAbs(2.0, 1e-10));
		REQUIRE_THAT(x[1], WithinAbs(3.0, 1e-10));
	}
	
	SECTION("AnalyzeMatrix")
	{
		Matrix<Real> A(3, 3, {1, 0, 0, 0, 1, 0, 0, 0, 1});
		
		auto analysis = AnalyzeMatrix(A);
		
		REQUIRE(analysis.isSquare == true);
		REQUIRE(analysis.rank == 3);
		REQUIRE(analysis.isDiagonal == true);
	}
}

//=============================================================================
// ERROR HANDLING
//=============================================================================

TEST_CASE("LinearSystem - Error Handling", "[LinearSystem][Errors]")
{
	SECTION("Solve without RHS throws")
	{
		Matrix<Real> A(3, 3, {1, 0, 0, 0, 1, 0, 0, 0, 1});
		LinearSystem<Real> sys(A);  // No RHS
		
		REQUIRE_THROWS(sys.Solve());
	}
	
	SECTION("Square operation on non-square throws")
	{
		Matrix<Real> A(4, 3);  // Non-square
		Vector<Real> b(4);
		LinearSystem<Real> sys(A, b);
		
		REQUIRE_THROWS(sys.SolveByLU());
		REQUIRE_THROWS(sys.Determinant());
	}
	
	SECTION("Cholesky on non-SPD throws")
	{
		Matrix<Real> A(3, 3, {1, 2, 0, 2, 1, 0, 0, 0, -1});  // Not positive definite
		Vector<Real> b({1, 2, 3});
		LinearSystem<Real> sys(A, b);
		
		REQUIRE_THROWS(sys.SolveByCholesky());
	}
}

//=============================================================================
// EIGENANALYSIS
//=============================================================================

TEST_CASE("LinearSystem - Eigenanalysis", "[LinearSystem][Eigen]")
{
	SECTION("Eigenvalues of symmetric matrix")
	{
		// Symmetric matrix with known eigenvalues
		Matrix<Real> A(3, 3, {4, 1, 1, 1, 4, 1, 1, 1, 4});
		LinearSystem<Real> sys(A);
		
		Vector<Real> eigs = sys.EigenvaluesSymmetric();
		
		// For this matrix, eigenvalues are 6, 3, 3
		REQUIRE(eigs.size() == 3);
		
		// Check product of eigenvalues equals determinant
		Real eigProduct = eigs[0] * eigs[1] * eigs[2];
		Real det = sys.Determinant();
		REQUIRE_THAT(eigProduct, WithinAbs(det, 0.01));
	}
	
	SECTION("Full eigensystem")
	{
		// Simple 2x2 with eigenvalues 3, 1
		Matrix<Real> A(2, 2, {2, 1, 1, 2});
		LinearSystem<Real> sys(A);
		
		auto result = sys.GetEigen();
		
		REQUIRE(result.converged == true);
		REQUIRE(result.eigenvalues.size() == 2);
		
		// Eigenvalues should be 3 and 1
		Real e1 = result.eigenvalues[0].real;
		Real e2 = result.eigenvalues[1].real;
		
		// Sort for comparison
		if (e1 < e2) std::swap(e1, e2);
		
		REQUIRE_THAT(e1, WithinAbs(3.0, 0.01));
		REQUIRE_THAT(e2, WithinAbs(1.0, 0.01));
	}
	
	SECTION("Spectral radius")
	{
		// Matrix with eigenvalues 3, 1
		Matrix<Real> A(2, 2, {2, 1, 1, 2});
		LinearSystem<Real> sys(A);
		
		Real spectral = sys.SpectralRadius();
		
		REQUIRE_THAT(spectral, WithinAbs(3.0, 0.01));
	}
	
	SECTION("Real eigenvalues detection")
	{
		// Symmetric matrix - always real eigenvalues
		Matrix<Real> A(2, 2, {1, 2, 2, 4});
		LinearSystem<Real> sys(A);
		
		REQUIRE(sys.HasComplexEigenvalues() == false);
	}
	
	SECTION("Rotation matrix - complex eigenvalues")
	{
		// 90-degree rotation matrix has eigenvalues ±i
		Matrix<Real> A(2, 2, {0, -1, 1, 0});
		LinearSystem<Real> sys(A);
		
		REQUIRE(sys.HasComplexEigenvalues() == true);
	}
}

//=============================================================================
// MULTIPLE RHS
//=============================================================================

TEST_CASE("LinearSystem - Multiple RHS", "[LinearSystem][MultipleRHS]")
{
	SECTION("Solve multiple systems")
	{
		Matrix<Real> A(3, 3, {4, 1, 0, 1, 4, 1, 0, 1, 4});
		Matrix<Real> B(3, 2, {5, 10, 6, 12, 5, 10});  // Two RHS
		
		LinearSystem<Real> sys(A, B);
		Matrix<Real> X = sys.SolveMultiple();
		
		// Verify each column
		REQUIRE(X.cols() == 2);
		
		// First column should be same solution scaled
		Vector<Real> col1 = X.VectorFromColumn(0);
		Vector<Real> col2 = X.VectorFromColumn(1);
		
		// col2 should be 2 * col1
		for (int i = 0; i < 3; ++i)
			REQUIRE_THAT(col2[i], WithinAbs(2.0 * col1[i], 1e-10));
	}
}

} // namespace MML::Tests::Systems::LinearSystemTests

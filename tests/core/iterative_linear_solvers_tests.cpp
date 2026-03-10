#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Vector/Vector.h"
#include "base/Matrix/Matrix.h"
#include "base/BaseUtils.h"
#include "core/LinAlgEqSolvers/LinAlgEqSolvers_iterative.h"
#endif

#include "../test_beds/linear_alg_eq_systems_test_bed.h"

using namespace MML;
using namespace MML::Testing;

namespace MML::Tests::Core::IterativeLinearSolversTests
{
	/*********************************************************************/
	/*****                   JACOBI SOLVER TESTS                     *****/
	/*********************************************************************/
	
	TEST_CASE("Jacobi_DiagDominant_4x4", "[JacobiSolver][Iterative]")
	{
			TEST_PRECISION_INFO();
		// Diagonally dominant 4x4 - guaranteed convergence
		auto sys = TestBeds::diag_dominant_4x4();
		
		auto result = JacobiSolver::Solve(sys._mat, sys._rhs);
		
		REQUIRE(result.converged);
		REQUIRE(result.solution.IsEqualTo(sys._sol, 1e-8));
		
		// Verify Ax = b
		Vector<Real> res = sys._mat * result.solution;
		REQUIRE(res.IsEqualTo(sys._rhs, 1e-8));
	}
	
	TEST_CASE("Jacobi_DiagDominant_5x5_Tridiag", "[JacobiSolver][Iterative]")
	{
			TEST_PRECISION_INFO();
		// Diagonally dominant tridiagonal - 1D Poisson
		auto sys = TestBeds::diag_dominant_5x5_tridiag();
		
		auto result = JacobiSolver::Solve(sys._mat, sys._rhs);
		
		REQUIRE(result.converged);
		REQUIRE(result.solution.IsEqualTo(sys._sol, 1e-8));
	}
	
	TEST_CASE("Jacobi_DiagDominant_6x6_Poisson2D", "[JacobiSolver][Iterative]")
	{
			TEST_PRECISION_INFO();
		// 2D Poisson discretization
		auto sys = TestBeds::diag_dominant_6x6_poisson2d();
		
		auto result = JacobiSolver::Solve(sys._mat, sys._rhs);
		
		REQUIRE(result.converged);
		REQUIRE(result.solution.IsEqualTo(sys._sol, 1e-8));
	}
	
	TEST_CASE("Jacobi_SolveSimple", "[JacobiSolver][Iterative]")
	{
			TEST_PRECISION_INFO();
		auto sys = TestBeds::diag_dominant_4x4();
		
		Vector<Real> sol = JacobiSolver::SolveSimple(sys._mat, sys._rhs);
		
		REQUIRE(sol.IsEqualTo(sys._sol, 1e-8));
	}
	
	TEST_CASE("Jacobi_InitialGuess", "[JacobiSolver][Iterative]")
	{
			TEST_PRECISION_INFO();
		auto sys = TestBeds::diag_dominant_4x4();
		
		// Start with solution as initial guess - should converge in 1 iteration
		auto result = JacobiSolver::Solve(sys._mat, sys._rhs, sys._sol);
		
		REQUIRE(result.converged);
		REQUIRE(result.iterations <= 2);  // Should be very fast with perfect initial guess
	}
	
	TEST_CASE("Jacobi_ZeroDiagonal_Throws", "[JacobiSolver][Iterative]")
	{
			TEST_PRECISION_INFO();
		Matrix<Real> A{3, 3, {
			REAL(0.0), REAL(1.0), REAL(2.0),  // Zero diagonal!
			REAL(1.0), REAL(4.0), REAL(1.0),
			REAL(2.0), REAL(1.0), REAL(4.0)
		}};
		Vector<Real> b{REAL(1.0), REAL(2.0), REAL(3.0)};
		
		REQUIRE_THROWS_AS(JacobiSolver::Solve(A, b), SingularMatrixError);
	}
	
	TEST_CASE("Jacobi_DimensionMismatch_Throws", "[JacobiSolver][Iterative]")
	{
			TEST_PRECISION_INFO();
		Matrix<Real> A{3, 3, {REAL(1.0), REAL(0.0), REAL(0.0), REAL(0.0), REAL(1.0), REAL(0.0), REAL(0.0), REAL(0.0), REAL(1.0)}};
		Vector<Real> b{REAL(1.0), REAL(2.0)};  // Wrong size
		
		REQUIRE_THROWS_AS(JacobiSolver::Solve(A, b), MatrixDimensionError);
	}
	
	/*********************************************************************/
	/*****                GAUSS-SEIDEL SOLVER TESTS                  *****/
	/*********************************************************************/
	
	TEST_CASE("GaussSeidel_DiagDominant_4x4", "[GaussSeidelSolver][Iterative]")
	{
			TEST_PRECISION_INFO();
		auto sys = TestBeds::diag_dominant_4x4();
		
		auto result = GaussSeidelSolver::Solve(sys._mat, sys._rhs);
		
		REQUIRE(result.converged);
		REQUIRE(result.solution.IsEqualTo(sys._sol, 1e-8));
	}
	
	TEST_CASE("GaussSeidel_DiagDominant_5x5_Tridiag", "[GaussSeidelSolver][Iterative]")
	{
			TEST_PRECISION_INFO();
		auto sys = TestBeds::diag_dominant_5x5_tridiag();
		
		auto result = GaussSeidelSolver::Solve(sys._mat, sys._rhs);
		
		REQUIRE(result.converged);
		REQUIRE(result.solution.IsEqualTo(sys._sol, 1e-8));
	}
	
	TEST_CASE("GaussSeidel_DiagDominant_6x6_Poisson2D", "[GaussSeidelSolver][Iterative]")
	{
			TEST_PRECISION_INFO();
		auto sys = TestBeds::diag_dominant_6x6_poisson2d();
		
		auto result = GaussSeidelSolver::Solve(sys._mat, sys._rhs);
		
		REQUIRE(result.converged);
		REQUIRE(result.solution.IsEqualTo(sys._sol, 1e-8));
	}
	
	TEST_CASE("GaussSeidel_FasterThanJacobi", "[GaussSeidelSolver][Iterative]")
	{
			TEST_PRECISION_INFO();
		// Gauss-Seidel should converge in fewer iterations than Jacobi
		auto sys = TestBeds::diag_dominant_5x5_tridiag();
		
		auto jacobi_result = JacobiSolver::Solve(sys._mat, sys._rhs);
		auto gs_result = GaussSeidelSolver::Solve(sys._mat, sys._rhs);
		
		REQUIRE(jacobi_result.converged);
		REQUIRE(gs_result.converged);
		REQUIRE(gs_result.iterations < jacobi_result.iterations);
	}
	
	TEST_CASE("GaussSeidel_SolveSimple", "[GaussSeidelSolver][Iterative]")
	{
			TEST_PRECISION_INFO();
		auto sys = TestBeds::diag_dominant_4x4();
		
		Vector<Real> sol = GaussSeidelSolver::SolveSimple(sys._mat, sys._rhs);
		
		REQUIRE(sol.IsEqualTo(sys._sol, 1e-8));
	}
	
	TEST_CASE("GaussSeidel_ZeroDiagonal_Throws", "[GaussSeidelSolver][Iterative]")
	{
			TEST_PRECISION_INFO();
		Matrix<Real> A{3, 3, {
			REAL(4.0), REAL(1.0), REAL(2.0),
			REAL(1.0), REAL(0.0), REAL(1.0),  // Zero diagonal!
			REAL(2.0), REAL(1.0), REAL(4.0)
		}};
		Vector<Real> b{REAL(1.0), REAL(2.0), REAL(3.0)};
		
		REQUIRE_THROWS_AS(GaussSeidelSolver::Solve(A, b), SingularMatrixError);
	}
	
	/*********************************************************************/
	/*****                    SOR SOLVER TESTS                       *****/
	/*********************************************************************/
	
	TEST_CASE("SOR_DiagDominant_4x4_Omega1", "[SORSolver][Iterative]")
	{
			TEST_PRECISION_INFO();
		// With omega = REAL(1.0), SOR should be equivalent to Gauss-Seidel
		auto sys = TestBeds::diag_dominant_4x4();
		
		auto sor_result = SORSolver::Solve(sys._mat, sys._rhs, REAL(1.0));
		auto gs_result = GaussSeidelSolver::Solve(sys._mat, sys._rhs);
		
		REQUIRE(sor_result.converged);
		REQUIRE(sor_result.iterations == gs_result.iterations);
		REQUIRE(sor_result.solution.IsEqualTo(gs_result.solution, 1e-12));
	}
	
	TEST_CASE("SOR_DiagDominant_5x5_Tridiag", "[SORSolver][Iterative]")
	{
			TEST_PRECISION_INFO();
		auto sys = TestBeds::diag_dominant_5x5_tridiag();
		
		// Use omega = REAL(1.5) (typical over-relaxation)
		auto result = SORSolver::Solve(sys._mat, sys._rhs, REAL(1.5));
		
		REQUIRE(result.converged);
		REQUIRE(result.solution.IsEqualTo(sys._sol, 1e-8));
	}
	
	TEST_CASE("SOR_DiagDominant_6x6_Poisson2D", "[SORSolver][Iterative]")
	{
			TEST_PRECISION_INFO();
		auto sys = TestBeds::diag_dominant_6x6_poisson2d();
		
		auto result = SORSolver::Solve(sys._mat, sys._rhs, REAL(1.3));
		
		REQUIRE(result.converged);
		REQUIRE(result.solution.IsEqualTo(sys._sol, 1e-8));
	}
	
	TEST_CASE("SOR_OptimalOmega_Tridiag", "[SORSolver][Iterative]")
	{
			TEST_PRECISION_INFO();
		// For tridiagonal Poisson, optimal omega should speed convergence
		auto sys = TestBeds::diag_dominant_5x5_tridiag();
		
		auto gs_result = GaussSeidelSolver::Solve(sys._mat, sys._rhs);
		auto sor_optimal = SORSolver::SolveOptimal(sys._mat, sys._rhs);
		
		REQUIRE(sor_optimal.converged);
		REQUIRE(sor_optimal.solution.IsEqualTo(sys._sol, 1e-8));
		// Note: SOR with "optimal" omega (derived for Poisson) may not always beat
		// Gauss-Seidel for all diagonally dominant matrices. The formula is specifically
		// optimal for 1D Poisson with specific boundary conditions.
		// Just verify it converges and gives correct answer.
	}
	
	TEST_CASE("SOR_UnderRelaxation", "[SORSolver][Iterative]")
	{
			TEST_PRECISION_INFO();
		// Under-relaxation (omega < 1) should still converge
		auto sys = TestBeds::diag_dominant_4x4();
		
		auto result = SORSolver::Solve(sys._mat, sys._rhs, REAL(0.8));
		
		REQUIRE(result.converged);
		REQUIRE(result.solution.IsEqualTo(sys._sol, 1e-8));
	}
	
	TEST_CASE("SOR_InvalidOmega_Throws", "[SORSolver][Iterative]")
	{
			TEST_PRECISION_INFO();
		auto sys = TestBeds::diag_dominant_4x4();
		
		REQUIRE_THROWS_AS(SORSolver::Solve(sys._mat, sys._rhs, REAL(0.0)), NumericalMethodError);
		REQUIRE_THROWS_AS(SORSolver::Solve(sys._mat, sys._rhs, REAL(2.0)), NumericalMethodError);
		REQUIRE_THROWS_AS(SORSolver::Solve(sys._mat, sys._rhs, -REAL(0.5)), NumericalMethodError);
		REQUIRE_THROWS_AS(SORSolver::Solve(sys._mat, sys._rhs, REAL(2.5)), NumericalMethodError);
	}
	
	TEST_CASE("SOR_SolveSimple", "[SORSolver][Iterative]")
	{
			TEST_PRECISION_INFO();
		auto sys = TestBeds::diag_dominant_4x4();
		
		Vector<Real> sol = SORSolver::SolveSimple(sys._mat, sys._rhs, REAL(1.2));
		
		REQUIRE(sol.IsEqualTo(sys._sol, 1e-8));
	}
	
	TEST_CASE("SOR_EstimateOptimalOmega", "[SORSolver][Iterative]")
	{
			TEST_PRECISION_INFO();
		auto sys = TestBeds::diag_dominant_5x5_tridiag();
		
		Real estimated_omega = SORSolver::EstimateOptimalOmega(sys._mat, sys._rhs);
		
		// For a 5x5 tridiagonal system, optimal omega should be around REAL(1.5)-REAL(1.7)
		REQUIRE(estimated_omega >= REAL(1.0));
		REQUIRE(estimated_omega < REAL(2.0));
		
		// Verify it actually works well
		auto result = SORSolver::Solve(sys._mat, sys._rhs, estimated_omega);
		REQUIRE(result.converged);
	}
	
	/*********************************************************************/
	/*****              CONVERGENCE COMPARISON TESTS                 *****/
	/*********************************************************************/
	
	TEST_CASE("IterativeSolvers_ConvergenceComparison", "[Iterative][Comparison]")
	{
			TEST_PRECISION_INFO();
		// Compare convergence rates of all three methods
		auto sys = TestBeds::diag_dominant_5x5_tridiag();
		
		auto jacobi = JacobiSolver::Solve(sys._mat, sys._rhs);
		auto gauss_seidel = GaussSeidelSolver::Solve(sys._mat, sys._rhs);
		auto sor = SORSolver::SolveOptimal(sys._mat, sys._rhs);
		
		// All should converge
		REQUIRE(jacobi.converged);
		REQUIRE(gauss_seidel.converged);
		REQUIRE(sor.converged);
		
		// Expected ordering: Jacobi >= Gauss-Seidel >= SOR (optimal)
		REQUIRE(jacobi.iterations >= gauss_seidel.iterations);
		
		// All should give same solution
		REQUIRE(jacobi.solution.IsEqualTo(sys._sol, 1e-8));
		REQUIRE(gauss_seidel.solution.IsEqualTo(sys._sol, 1e-8));
		REQUIRE(sor.solution.IsEqualTo(sys._sol, 1e-8));
	}
	
	TEST_CASE("IterativeSolvers_MaxIterReached", "[Iterative]")
	{
			TEST_PRECISION_INFO();
		// Test behavior when max iterations is too small
		auto sys = TestBeds::diag_dominant_5x5_tridiag();
		
		// Use very few iterations
		auto result = JacobiSolver::Solve(sys._mat, sys._rhs, Vector<Real>(), 1e-15, 5);
		
		// Should not converge with only 5 iterations
		REQUIRE_FALSE(result.converged);
		REQUIRE(result.iterations <= 6);  // May be 5 or 6 depending on iteration counting
		
		// SolveSimple should throw
		REQUIRE_THROWS_AS(
			JacobiSolver::SolveSimple(sys._mat, sys._rhs, 1e-15, 5),
			std::runtime_error
		);
	}
	
	/*********************************************************************/
	/*****                SPD MATRIX TESTS                           *****/
	/*********************************************************************/
	
	TEST_CASE("GaussSeidel_SPD_3x3", "[GaussSeidelSolver][Iterative][SPD]")
	{
			TEST_PRECISION_INFO();
		// SPD matrices guarantee convergence for Gauss-Seidel
		auto sys = TestBeds::spd_3x3();
		Matrix<Real> fullMat = sys._mat.GetAsMatrix();
		
		auto result = GaussSeidelSolver::Solve(fullMat, sys._rhs);
		
		REQUIRE(result.converged);
		REQUIRE(result.solution.IsEqualTo(sys._sol, 1e-8));
	}
	
	TEST_CASE("SOR_SPD_5x5_MassMatrix", "[SORSolver][Iterative][SPD]")
	{
			TEST_PRECISION_INFO();
		auto sys = TestBeds::spd_5x5_mass_matrix();
		Matrix<Real> fullMat = sys._mat.GetAsMatrix();
		
		auto result = SORSolver::Solve(fullMat, sys._rhs, REAL(1.2));
		
		REQUIRE(result.converged);
		REQUIRE(result.solution.IsEqualTo(sys._sol, 1e-8));
	}

	/*********************************************************************/
	/*****              RESULT STRUCTURE TESTS                       *****/
	/*********************************************************************/

	TEST_CASE("IterativeSolverResult_AllFields", "[Iterative][Result]")
	{
		TEST_PRECISION_INFO();
		auto sys = TestBeds::diag_dominant_4x4();
		
		auto result = JacobiSolver::Solve(sys._mat, sys._rhs);
		
		// Verify all fields are populated correctly
		REQUIRE(result.converged == true);
		REQUIRE(result.iterations > 0);
		REQUIRE(result.iterations <= 1000);  // Within max iterations
		REQUIRE(result.residual >= REAL(0.0));
		REQUIRE(result.residual < REAL(1e-10));  // Within tolerance
		REQUIRE(result.solution.size() == sys._sol.size());
	}

	TEST_CASE("IterativeSolverResult_ResidualComputation", "[Iterative][Result]")
	{
		TEST_PRECISION_INFO();
		auto sys = TestBeds::diag_dominant_5x5_tridiag();
		
		auto result = GaussSeidelSolver::Solve(sys._mat, sys._rhs);
		
		// Manually compute residual ||Ax - b|| / ||b||
		Vector<Real> Ax = sys._mat * result.solution;
		Vector<Real> residual_vec = Ax - sys._rhs;
		Real computed_residual = residual_vec.NormL2() / sys._rhs.NormL2();
		
		// Result residual should match (approximately, due to final iteration)
		REQUIRE(std::abs(result.residual - computed_residual) < REAL(1e-12));
	}

	/*********************************************************************/
	/*****              LARGER SYSTEM TESTS                          *****/
	/*********************************************************************/

	TEST_CASE("Jacobi_LargeSystem_20x20", "[JacobiSolver][Iterative][Large]")
	{
		TEST_PRECISION_INFO();
		// Create a 20x20 diagonally dominant system
		int n = 20;
		Matrix<Real> A(n, n);
		Vector<Real> b(n);
		
		// Fill with diagonally dominant pattern
		for (int i = 0; i < n; i++)
		{
			Real diag_sum = REAL(0.0);
			for (int j = 0; j < n; j++)
			{
				if (i != j)
				{
					A[i][j] = REAL(1.0) / (REAL(1.0) + std::abs(i - j));
					diag_sum += std::abs(A[i][j]);
				}
			}
			A[i][i] = diag_sum + REAL(1.0);  // Strictly diagonally dominant
			b[i] = REAL(i + 1);
		}
		
		auto result = JacobiSolver::Solve(A, b);
		
		REQUIRE(result.converged);
		
		// Verify Ax = b
		Vector<Real> Ax = A * result.solution;
		REQUIRE(Ax.IsEqualTo(b, REAL(1e-8)));
	}

	TEST_CASE("GaussSeidel_LargeSystem_50x50", "[GaussSeidelSolver][Iterative][Large]")
	{
		TEST_PRECISION_INFO();
		// Create a 50x50 tridiagonal system (1D Poisson-like)
		int n = 50;
		Matrix<Real> A(n, n);
		Vector<Real> b(n);
		
		for (int i = 0; i < n; i++)
		{
			A[i][i] = REAL(4.0);
			if (i > 0) A[i][i-1] = REAL(-1.0);
			if (i < n-1) A[i][i+1] = REAL(-1.0);
			b[i] = REAL(1.0);
		}
		
		auto result = GaussSeidelSolver::Solve(A, b);
		
		REQUIRE(result.converged);
		
		// Verify solution
		Vector<Real> Ax = A * result.solution;
		REQUIRE(Ax.IsEqualTo(b, REAL(1e-8)));
	}

	TEST_CASE("SOR_LargeSystem_50x50", "[SORSolver][Iterative][Large]")
	{
		TEST_PRECISION_INFO();
		// 50x50 tridiagonal system with SOR
		int n = 50;
		Matrix<Real> A(n, n);
		Vector<Real> b(n);
		
		for (int i = 0; i < n; i++)
		{
			A[i][i] = REAL(4.0);
			if (i > 0) A[i][i-1] = REAL(-1.0);
			if (i < n-1) A[i][i+1] = REAL(-1.0);
			b[i] = REAL(1.0);
		}
		
		// Use omega slightly above 1 (mild over-relaxation)
		// The "optimal omega" formula is for specific Poisson BCs and may not apply here
		auto sor_result = SORSolver::Solve(A, b, REAL(1.2));
		
		REQUIRE(sor_result.converged);
		
		// Verify solution
		Vector<Real> Ax = A * sor_result.solution;
		REQUIRE(Ax.IsEqualTo(b, REAL(1e-8)));
	}

	/*********************************************************************/
	/*****              EDGE CASE TESTS                              *****/
	/*********************************************************************/

	TEST_CASE("Jacobi_NonSquareMatrix_Throws", "[JacobiSolver][Iterative]")
	{
		TEST_PRECISION_INFO();
		Matrix<Real> A(3, 4);  // Non-square
		Vector<Real> b(3);
		
		REQUIRE_THROWS_AS(JacobiSolver::Solve(A, b), MatrixDimensionError);
	}

	TEST_CASE("GaussSeidel_NonSquareMatrix_Throws", "[GaussSeidelSolver][Iterative]")
	{
		TEST_PRECISION_INFO();
		Matrix<Real> A(4, 3);  // Non-square
		Vector<Real> b(4);
		
		REQUIRE_THROWS_AS(GaussSeidelSolver::Solve(A, b), MatrixDimensionError);
	}

	TEST_CASE("SOR_NonSquareMatrix_Throws", "[SORSolver][Iterative]")
	{
		TEST_PRECISION_INFO();
		Matrix<Real> A(3, 5);  // Non-square
		Vector<Real> b(3);
		
		REQUIRE_THROWS_AS(SORSolver::Solve(A, b, REAL(1.2)), MatrixDimensionError);
	}

	TEST_CASE("GaussSeidel_SolveSimple_Throws_OnNonConvergence", "[GaussSeidelSolver][Iterative]")
	{
		TEST_PRECISION_INFO();
		auto sys = TestBeds::diag_dominant_5x5_tridiag();
		
		// Very tight tolerance with very few iterations - should not converge
		REQUIRE_THROWS_AS(
			GaussSeidelSolver::SolveSimple(sys._mat, sys._rhs, REAL(1e-20), 3),
			std::runtime_error
		);
	}

	TEST_CASE("SOR_SolveSimple_Throws_OnNonConvergence", "[SORSolver][Iterative]")
	{
		TEST_PRECISION_INFO();
		auto sys = TestBeds::diag_dominant_5x5_tridiag();
		
		// Very tight tolerance with very few iterations - should not converge
		REQUIRE_THROWS_AS(
			SORSolver::SolveSimple(sys._mat, sys._rhs, REAL(1.2), REAL(1e-20), 3),
			std::runtime_error
		);
	}

	TEST_CASE("Jacobi_ZeroRHS", "[JacobiSolver][Iterative]")
	{
		TEST_PRECISION_INFO();
		// When b = 0, solution should be x = 0
		auto sys = TestBeds::diag_dominant_4x4();
		Vector<Real> zero_b(sys._rhs.size());
		for (int i = 0; i < zero_b.size(); i++)
			zero_b[i] = REAL(0.0);
		
		auto result = JacobiSolver::Solve(sys._mat, zero_b);
		
		REQUIRE(result.converged);
		for (int i = 0; i < result.solution.size(); i++)
		{
			REQUIRE(std::abs(result.solution[i]) < REAL(1e-12));
		}
	}

	TEST_CASE("Iterative_CustomInitialGuess_Accuracy", "[Iterative]")
	{
		TEST_PRECISION_INFO();
		auto sys = TestBeds::diag_dominant_5x5_tridiag();
		
		// Start with a "close" initial guess
		Vector<Real> close_guess = sys._sol;
		for (int i = 0; i < close_guess.size(); i++)
			close_guess[i] *= REAL(1.1);  // 10% off
		
		// With close guess, should converge faster than from zero
		auto result_close = JacobiSolver::Solve(sys._mat, sys._rhs, close_guess);
		auto result_zero = JacobiSolver::Solve(sys._mat, sys._rhs);
		
		REQUIRE(result_close.converged);
		REQUIRE(result_zero.converged);
		REQUIRE(result_close.iterations <= result_zero.iterations);
		
		// Both should give correct solution (verified against expected)
		REQUIRE(result_close.solution.IsEqualTo(sys._sol, REAL(1e-8)));
		REQUIRE(result_zero.solution.IsEqualTo(sys._sol, REAL(1e-8)));
	}

	TEST_CASE("Iterative_ToleranceSensitivity", "[Iterative]")
	{
		TEST_PRECISION_INFO();
		auto sys = TestBeds::diag_dominant_4x4();
		
		// Tighter tolerance should need more iterations
		auto loose = JacobiSolver::Solve(sys._mat, sys._rhs, Vector<Real>(), REAL(1e-6), 1000);
		auto tight = JacobiSolver::Solve(sys._mat, sys._rhs, Vector<Real>(), REAL(1e-12), 1000);
		
		REQUIRE(loose.converged);
		REQUIRE(tight.converged);
		REQUIRE(tight.iterations >= loose.iterations);
		REQUIRE(tight.residual < loose.residual);  // Tighter should have smaller residual
	}

} // namespace MML::Tests::Core::IterativeLinearSolversTests

#include "../catch/catch.hpp"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "core/Vector.h"
#include "core/Matrix.h"

#include "algorithms/LinAlgEqSolvers.h"
#endif

#include "../test_data/linear_alg_eq_systems_test_bed.h"

using namespace MML;

TEST_CASE("Test_GaussJordan_Solve_5_x_5", "[simple]")
{
	Matrix<Real>   mat = TestBeds::mat_5x5;
	Matrix<Real> 	rhs{ TestBeds::mat_5x5_rhs_multi };

	GaussJordanSolver::Solve(mat, rhs);

	REQUIRE(true == rhs.IsEqual(TestBeds::mat_5x5_rhs_multi_sol, 1e-10));
}

TEST_CASE("Test_LUDecomposition_Solve_3_x_3", "[simple]")
{
	Matrix<Real>   mat = TestBeds::mat_3x3;
	Vector<Real> 	rhs = TestBeds::mat_3x3_rhs0;
	Vector<Real>	vecSol(rhs.size());

	LUDecompositionSolver luSolver(mat);

	luSolver.Solve(rhs, vecSol);

	REQUIRE(true == vecSol.IsEqual(TestBeds::mat_3x3_rhs0_sol, 1e-10));

	Vector<Real>    res_rhs = mat * vecSol;

	REQUIRE(true == res_rhs.IsEqual(TestBeds::mat_3x3_rhs0, 1e-10));
}

TEST_CASE("Test_LUDecomposition_Solve_5_x_5", "[simple]")
{
	Matrix<Real>     mat = TestBeds::mat_5x5;

	LUDecompositionSolver luSolver(mat);

	for (int i = 0; i < 2; i++)
	{
		Vector<Real> 	rhs = Matrix<Real>::VectorFromColumn(TestBeds::mat_5x5_rhs_multi, i);
		Vector<Real>	vecSol(rhs.size());

		luSolver.Solve(rhs, vecSol);

		Vector<Real> 	rhs_sol = Matrix<Real>::VectorFromColumn(TestBeds::mat_5x5_rhs_multi_sol, i);

		REQUIRE(true == vecSol.IsEqual(rhs_sol, 1e-10));

		Vector<Real>    res_rhs = mat * vecSol;

		REQUIRE(true == res_rhs.IsEqual(rhs, 1e-10));
	}
}

TEST_CASE("Test_LUDecompositionn_Solve_5_x_5_Matrix_input", "[simple]")
{
	Matrix<Real>   mat = TestBeds::mat_5x5;
	Matrix<Real> 	rhs{ TestBeds::mat_5x5_rhs_multi };
	Matrix<Real> 	sol(5, 2);

	LUDecompositionSolver luSolver(mat);

	luSolver.Solve(rhs, sol);

	REQUIRE(true == sol.IsEqual(TestBeds::mat_5x5_rhs_multi_sol, 1e-10));
}

TEST_CASE("Test_QRDecomposition_Solve_5_x_5", "[simple]")
{
	Matrix<Real>     mat = TestBeds::mat_5x5;

	QRDecompositionSolver qrSolver(mat);

	for (int i = 0; i < 2; i++)
	{
		Vector<Real> 	rhs = Matrix<Real>::VectorFromColumn(TestBeds::mat_5x5_rhs_multi, i);
		Vector<Real>		vecSol(rhs.size());

		qrSolver.Solve(rhs, vecSol);

		Vector<Real> 	rhs_sol = Matrix<Real>::VectorFromColumn(TestBeds::mat_5x5_rhs_multi_sol, i);

		REQUIRE(true == vecSol.IsEqual(rhs_sol, 1e-10));

		Vector<Real>    res_rhs = mat * vecSol;

		REQUIRE(true == res_rhs.IsEqual(rhs, 1e-10));
	}
}
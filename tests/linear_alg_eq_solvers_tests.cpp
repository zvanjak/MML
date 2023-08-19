#include "../catch/catch.hpp"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "core/Vector.h"
#include "core/Matrix.h"

#include "algorithms/LinAlgEqSolvers.h"
#endif

#include "../test_data/linear_alg_eq_systems_test_bed.h"


TEST_CASE("Test_GaussJordan_Solve_5_x_5", "[simple]")
{
	MML::Matrix<Real>   mat = MML::TestBeds::mat_5x5;
	MML::Matrix<Real> 	rhs{ MML::TestBeds::mat_5x5_rhs_multi };

	MML::GaussJordanSolver::Solve(mat, rhs);

	REQUIRE(true == rhs.IsEqual(MML::TestBeds::mat_5x5_rhs0_sol, 1e-10));
}

TEST_CASE("Test_LUDecomposition_Solve_3_x_3", "[simple]")
{
	MML::Matrix<Real>   mat = MML::TestBeds::mat_3x3;
	MML::Vector<Real> 	rhs = MML::TestBeds::mat_3x3_rhs0;
	MML::Vector<Real>	vecSol(rhs.size());

	MML::LUDecompositionSolver luSolver(mat);

	luSolver.Solve(rhs, vecSol);

	REQUIRE(true == vecSol.IsEqual(MML::TestBeds::mat_3x3_rhs0_sol, 1e-10));

	MML::Vector<Real>    res_rhs = mat * vecSol;

	REQUIRE(true == res_rhs.IsEqual(MML::TestBeds::mat_3x3_rhs0, 1e-10));
}

TEST_CASE("Test_LUDecomposition_Solve_5_x_5", "[simple]")
{
	MML::Matrix<Real>     mat = MML::TestBeds::mat_5x5;

	MML::LUDecompositionSolver luSolver(mat);

	for (int i = 0; i < 2; i++)
	{
		MML::Vector<Real> 	rhs = MML::Matrix<Real>::VectorFromColumn(MML::TestBeds::mat_5x5_rhs_multi, i);
		MML::Vector<Real>	vecSol(rhs.size());

		luSolver.Solve(rhs, vecSol);

		MML::Vector<Real> 	rhs_sol = MML::Matrix<Real>::VectorFromColumn(MML::TestBeds::mat_5x5_rhs0_sol, i);

		REQUIRE(true == vecSol.IsEqual(rhs_sol, 1e-10));

		MML::Vector<Real>    res_rhs = mat * vecSol;

		REQUIRE(true == res_rhs.IsEqual(rhs, 1e-10));
	}
}

TEST_CASE("Test_LUDecompositionn_Solve_5_x_5_Matrix_input", "[simple]")
{
	MML::Matrix<Real>   mat = MML::TestBeds::mat_5x5;
	MML::Matrix<Real> 	rhs{ MML::TestBeds::mat_5x5_rhs_multi };
	MML::Matrix<Real> 	sol(5, 2);

	MML::LUDecompositionSolver luSolver(mat);

	luSolver.Solve(rhs, sol);

	REQUIRE(true == sol.IsEqual(MML::TestBeds::mat_5x5_rhs0_sol, 1e-10));
}

TEST_CASE("Test_QRDecomposition_Solve_5_x_5", "[simple]")
{
	MML::Matrix<Real>     mat = MML::TestBeds::mat_5x5;

	MML::QRDecompositionSolver qrSolver(mat);

	for (int i = 0; i < 2; i++)
	{
		MML::Vector<Real> 	rhs = MML::Matrix<Real>::VectorFromColumn(MML::TestBeds::mat_5x5_rhs_multi, i);
		MML::Vector<Real>		vecSol(rhs.size());

		qrSolver.Solve(rhs, vecSol);

		MML::Vector<Real> 	rhs_sol = MML::Matrix<Real>::VectorFromColumn(MML::TestBeds::mat_5x5_rhs0_sol, i);

		REQUIRE(true == vecSol.IsEqual(rhs_sol, 1e-10));

		MML::Vector<Real>    res_rhs = mat * vecSol;

		REQUIRE(true == res_rhs.IsEqual(rhs, 1e-10));
	}
}
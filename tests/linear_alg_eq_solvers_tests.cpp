#include "../catch/catch.hpp"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "basic_types/Vector.h"
#include "basic_types/Matrix.h"
#include "algorithms/LinAlgEqSolvers.h"
#endif

#include "test_data/linear_alg_eq_solvers_test_bed.h"


TEST_CASE("Test_GaussJordan_Solve_5_x_5", "[simple]") 
{
    MML::Matrix     mat = MML::Tests::LinearAlgEqTestBed::mat4;
    MML::Matrix 	rhs{MML::Tests::LinearAlgEqTestBed::mat4_rhs4};

    MML::GaussJordanSolver::Solve(mat, rhs);

    REQUIRE( true == rhs.IsEqual(MML::Tests::LinearAlgEqTestBed::mat4_rhs0_sol, 1e-10));
}

TEST_CASE("Test_LUDecomposition_Solve_3_x_3", "[simple]") 
{
    MML::Matrix     mat = MML::Tests::LinearAlgEqTestBed::mat0;
	MML::Vector 	rhs = MML::Tests::LinearAlgEqTestBed::mat0_rhs0;
	MML::Vector		vecSol(rhs.size());

    MML::LUDecompositionSolver luSolver(mat);

    luSolver.Solve(rhs, vecSol);

	REQUIRE( true == vecSol.IsEqual(MML::Tests::LinearAlgEqTestBed::mat0_rhs0_sol, 1e-10));

    MML::Vector    res_rhs = mat * vecSol;

	REQUIRE( true == res_rhs.IsEqual(MML::Tests::LinearAlgEqTestBed::mat0_rhs0, 1e-10));
}

TEST_CASE("Test_LUDecomposition_Solve_5_x_5", "[simple]") 
{
    MML::Matrix     mat = MML::Tests::LinearAlgEqTestBed::mat4;

    MML::LUDecompositionSolver luSolver(mat);

    for(int i=0; i<2; i++ )
    {
        MML::Vector 	rhs = MML::Matrix::VectorFromColumn(MML::Tests::LinearAlgEqTestBed::mat4_rhs4, i);
        MML::Vector		vecSol(rhs.size());

        luSolver.Solve(rhs, vecSol);

        MML::Vector 	rhs_sol = MML::Matrix::VectorFromColumn(MML::Tests::LinearAlgEqTestBed::mat4_rhs0_sol, i);

        REQUIRE( true == vecSol.IsEqual(rhs_sol, 1e-10));

        MML::Vector    res_rhs = mat * vecSol;

        REQUIRE( true == res_rhs.IsEqual(rhs, 1e-10));
    }
}

TEST_CASE("Test_LUDecompositionn_Solve_5_x_5_Matrix_input", "[simple]") 
{
    MML::Matrix     mat = MML::Tests::LinearAlgEqTestBed::mat4;
    MML::Matrix 	rhs{MML::Tests::LinearAlgEqTestBed::mat4_rhs4};
    MML::Matrix 	sol(5,2);
    
    MML::LUDecompositionSolver luSolver(mat);

    luSolver.Solve(rhs, sol);

    REQUIRE( true == sol.IsEqual(MML::Tests::LinearAlgEqTestBed::mat4_rhs0_sol, 1e-10));
}

TEST_CASE("Test_QRDecomposition_Solve_5_x_5", "[simple]") 
{
    MML::Matrix     mat = MML::Tests::LinearAlgEqTestBed::mat4;

    MML::QRDecompositionSolver qrSolver(mat);

    for(int i=0; i<2; i++ )
    {
        MML::Vector 	rhs = MML::Matrix::VectorFromColumn(MML::Tests::LinearAlgEqTestBed::mat4_rhs4, i);
        MML::Vector		vecSol(rhs.size());

        qrSolver.Solve(rhs, vecSol);

        MML::Vector 	rhs_sol = MML::Matrix::VectorFromColumn(MML::Tests::LinearAlgEqTestBed::mat4_rhs0_sol, i);

        REQUIRE( true == vecSol.IsEqual(rhs_sol, 1e-10));

        MML::Vector    res_rhs = mat * vecSol;

        REQUIRE( true == res_rhs.IsEqual(rhs, 1e-10));
    }
}
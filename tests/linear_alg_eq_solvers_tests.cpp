#include "../catch/catch.hpp"

#ifdef MML_USE_SINGLE_HEADER
#include "MMLBasicTypes.h"
#else
#include "basic_types/Vector.h"
#include "basic_types/Matrix.h"
#include "algorithms/LinAlgEqSolvers.h"
#endif

#include "linear_alg_eq_solvers_test_bed.h"

TEST_CASE("Test_LUDecomposition_Solve", "[simple]") 
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
#include "../catch/catch.hpp"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Vector.h"
#include "base/Matrix.h"

#include "base/BaseUtils.h"

#include "core/LinAlgEqSolvers.h"
#endif

#include "../test_data/linear_alg_eq_systems_test_bed.h"

using namespace MML;

namespace MML::Tests::Core::LinearAlgSolversTests
{
	/*********************************************************************/
	/*****                 Gauss Jordan solver real                  *****/
	/*********************************************************************/

	TEST_CASE("Test_GaussJordan_Solve_5_x_5", "[simple]")
	{
		Matrix<Real>    mat = TestBeds::mat_5x5;
		Vector<Real>    rhs = TestBeds::mat_5x5_rhs0;

		GaussJordanSolver<Real>::Solve(mat, rhs);
		Vector<Real> vecSol = rhs;

		REQUIRE(true == rhs.IsEqualTo(TestBeds::mat_5x5_rhs0_sol, 1e-14));

		Vector<Real>    res_rhs = TestBeds::mat_5x5 * vecSol;
		REQUIRE(true == res_rhs.IsEqualTo(TestBeds::mat_5x5_rhs0, 1e-14));
	}
	TEST_CASE("Test_GaussJordan_Solve_5_x_5_multi", "[simple]")
	{
		Matrix<Real>    mat = TestBeds::mat_5x5;
		Matrix<Real>    rhs = TestBeds::mat_5x5_rhs_multi;

		GaussJordanSolver<Real>::Solve(mat, rhs);

		REQUIRE(true == rhs.IsEqualTo(TestBeds::mat_5x5_rhs_multi_sol, 1e-13));
	}

	/*********************************************************************/
	/*****                 Gauss Jordan solver complex               *****/
	/*********************************************************************/

	TEST_CASE("Test_GaussJordan_Solve_Complex_3_x_3", "[simple]")
	{
		auto    mat = TestBeds::mat_cmplx_1_3x3;
		auto    rhs = TestBeds::mat_cmplx_1_3x3_rhs0;

		GaussJordanSolver<Complex>::Solve(mat, rhs);
		Vector<Complex> vecSol = rhs;

		REQUIRE(true == rhs.IsEqualTo(TestBeds::mat_cmplx_1_3x3_rhs0_sol, 1e-14));

		Vector<Complex>   res_rhs = TestBeds::mat_cmplx_1_3x3 * vecSol;
		REQUIRE(true == res_rhs.IsEqualTo(TestBeds::mat_cmplx_1_3x3_rhs0, 1e-14));
	}

	/*********************************************************************/
	/*****                    Complete test beds                     *****/
	/*********************************************************************/

	TEST_CASE("Test_GaussJordan_COMPLETE_TEST_BED", "[simple]")
	{
		for (int i = 0; i < TestBeds::LinearAlgEqTestBed::numLinAlgEqSystems(); i++)
		{
			TestBeds::TestLinearSystem testBed = TestBeds::LinearAlgEqTestBed::getLinAlgEqSystem(i);

			Matrix<Real>     mat = testBed._mat;
			Vector<Real>     rhs = testBed._rhs;

			GaussJordanSolver<Real>::Solve(mat, rhs);
			REQUIRE(true == rhs.IsEqualTo(testBed._sol, 1e-13));

			Vector<Real>    res_rhs = testBed._mat * rhs;
			REQUIRE(true == res_rhs.IsEqualTo(testBed._rhs, 1e-13));
		}
	}

	// TODO - LOW, Gauss Jordan solver COMPLEX COMPLETE TEST BED

}
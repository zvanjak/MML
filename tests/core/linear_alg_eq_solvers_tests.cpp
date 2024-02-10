#include "../catch/catch.hpp"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Vector.h"
#include "base/Matrix.h"

#include "core/CoreUtils.h"

#include "core/LinAlgEqSolvers.h"
#endif

#include "../test_data/linear_alg_eq_systems_test_bed.h"

using namespace MML;

// TODO - HIGH, dodati još testova za kompleksne matrice

///////////////////////////////////////////////////////////////////////////////////////////
/*********                       Gauss Jordan solver REAL                         ********/
///////////////////////////////////////////////////////////////////////////////////////////
TEST_CASE("Test_GaussJordan_Solve_5_x_5", "[simple]")
{
	Matrix<Real>    mat = TestBeds::mat_5x5;
	Vector<Real>    rhs = TestBeds::mat_5x5_rhs0;

    GaussJordanSolver<Real>::Solve(mat, rhs);
	Vector<Real> vecSol = rhs;

	REQUIRE(true == rhs.IsEqual(TestBeds::mat_5x5_rhs0_sol, 1e-14));

	Vector<Real>    res_rhs = TestBeds::mat_5x5 * vecSol;
	REQUIRE(true == res_rhs.IsEqual(TestBeds::mat_5x5_rhs0, 1e-14));
}
TEST_CASE("Test_GaussJordan_Solve_5_x_5_multi", "[simple]")
{
	Matrix<Real>    mat = TestBeds::mat_5x5;
	Matrix<Real>    rhs = TestBeds::mat_5x5_rhs_multi;

	GaussJordanSolver<Real>::Solve(mat, rhs);

	REQUIRE(true == rhs.IsEqual(TestBeds::mat_5x5_rhs_multi_sol, 1e-13));   
}

///////////////////////////////////////////////////////////////////////////////////////////
/*********                     Gauss Jordan solver COMPLEX                        ********/
///////////////////////////////////////////////////////////////////////////////////////////
TEST_CASE("Test_GaussJordan_Solve_Complex_3_x_3", "[simple]")
{
	auto    mat = TestBeds::mat_cmplx_1_3x3;
	auto    rhs = TestBeds::mat_cmplx_1_3x3_rhs0;

	GaussJordanSolver<Complex>::Solve(mat, rhs);
	Vector<Complex> vecSol = rhs;

	REQUIRE(true == rhs.IsEqual(TestBeds::mat_cmplx_1_3x3_rhs0_sol, 1e-14));

	Vector<Complex>   res_rhs = TestBeds::mat_cmplx_1_3x3 * vecSol;
	REQUIRE(true == res_rhs.IsEqual(TestBeds::mat_cmplx_1_3x3_rhs0, 1e-14));    
}

///////////////////////////////////////////////////////////////////////////////////////////
/**********                       LU decomposition REAL                           ********/
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Test_LUDecomposition_Solve_3_x_3", "[simple]")
{
	Matrix<Real>    mat = TestBeds::mat_3x3;
	Vector<Real> 	rhs = TestBeds::mat_3x3_rhs0;
	Vector<Real>	vecSol(rhs.size());

	LUDecompositionSolver<Real> luSolver(mat);

	luSolver.Solve(rhs, vecSol);
	REQUIRE(true == vecSol.IsEqual(TestBeds::mat_3x3_rhs0_sol, 1e-15));

	Vector<Real>    res_rhs = mat * vecSol;
	REQUIRE(true == res_rhs.IsEqual(TestBeds::mat_3x3_rhs0, 1e-15));
}

TEST_CASE("Test_LUDecomposition_Solve_5_x_5", "[simple]")
{
	Matrix<Real>    mat = TestBeds::mat_5x5;
	Vector<Real> 	rhs = TestBeds::mat_5x5_rhs0;
	Vector<Real>	vecSol(rhs.size());

	LUDecompositionSolver<Real> luSolver(mat);

	luSolver.Solve(rhs, vecSol);
	REQUIRE(true == vecSol.IsEqual(TestBeds::mat_5x5_rhs0_sol, 1e-16));

	Vector<Real>    res_rhs = mat * vecSol;
	REQUIRE(true == res_rhs.IsEqual(TestBeds::mat_5x5_rhs0, 1e-14));
}

TEST_CASE("Test_LUDecompositionn_Solve_5_x_5_multi", "[simple]")
{
	Matrix<Real>    mat = TestBeds::mat_5x5;
	Matrix<Real> 	rhs{ TestBeds::mat_5x5_rhs_multi };
	Matrix<Real> 	sol(5, 2);

	LUDecompositionSolver<Real> luSolver(mat);

	luSolver.Solve(rhs, sol);
	REQUIRE(true == sol.IsEqual(TestBeds::mat_5x5_rhs_multi_sol, 1e-10));
}

TEST_CASE("Test_LUDecomposition_Solve_8_x_8", "[simple]")
{
	Matrix<Real>    mat = TestBeds::mat_8x8;
	Vector<Real> 	rhs = TestBeds::mat_8x8_rhs0;
	Vector<Real>	vecSol(rhs.size());

	LUDecompositionSolver<Real> luSolver(mat);

	luSolver.Solve(rhs, vecSol);
	REQUIRE(true == vecSol.IsEqual(TestBeds::mat_8x8_rhs0_sol, 1e-14));

	Vector<Real>    res_rhs = mat * vecSol;
	REQUIRE(true == res_rhs.IsEqual(TestBeds::mat_8x8_rhs0, 1e-14));
}

TEST_CASE("Test_LUDecomposition_Solve_10_x_10", "[simple]")
{
	Matrix<Real>    mat = TestBeds::mat_10x10;
	Vector<Real> 	rhs = TestBeds::mat_10x10_rhs0;
	Vector<Real>	vecSol(rhs.size());

	LUDecompositionSolver<Real> luSolver(mat);

	luSolver.Solve(rhs, vecSol);
	REQUIRE(true == vecSol.IsEqual(TestBeds::mat_10x10_rhs0_sol, 1e-14));

	Vector<Real>    res_rhs = mat * vecSol;
	REQUIRE(true == res_rhs.IsEqual(TestBeds::mat_10x10_rhs0, 1e-14));
}

TEST_CASE("Test_LUDecomposition_Solve_20_x_20", "[simple]")
{
	Matrix<Real>    mat = TestBeds::mat_20x20;
	Vector<Real> 	rhs = TestBeds::mat_20x20_rhs0;
	Vector<Real>	vecSol(rhs.size());

	LUDecompositionSolver<Real> luSolver(mat);

	luSolver.Solve(rhs, vecSol);
	REQUIRE(true == vecSol.IsEqual(TestBeds::mat_20x20_rhs0_sol, 1e-13));

	Vector<Real>    res_rhs = mat * vecSol;
	REQUIRE(true == res_rhs.IsEqual(TestBeds::mat_20x20_rhs0, 1e-13));
}

TEST_CASE("Test_LUDecomposition_Solve_50_x_50", "[simple]")
{
	Matrix<Real>    mat = TestBeds::mat_50x50;
	Vector<Real> 	rhs = TestBeds::mat_50x50_rhs0;
	Vector<Real>	vecSol(rhs.size());

	LUDecompositionSolver<Real> luSolver(mat);

	luSolver.Solve(rhs, vecSol);
	REQUIRE(true == vecSol.IsEqual(TestBeds::mat_50x50_rhs0_sol, 1e-13));

	Vector<Real>    res_rhs = mat * vecSol;
	REQUIRE(true == res_rhs.IsEqual(TestBeds::mat_50x50_rhs0, 1e-12));
}

///////////////////////////////////////////////////////////////////////////////////////////
/**********                      LU decomposition COMPLEX                         ********/
///////////////////////////////////////////////////////////////////////////////////////////
TEST_CASE("Test_LUDecomposition_Solve_3_x_3_complex", "[simple]")
{
	Matrix<Complex> mat = TestBeds::mat_cmplx_1_3x3;
	Vector<Complex> rhs = TestBeds::mat_cmplx_1_3x3_rhs0;
	Vector<Complex> vecSol(rhs.size());

	LUDecompositionSolver<Complex> luSolver(mat);

	luSolver.Solve(rhs, vecSol);
	REQUIRE(true == vecSol.IsEqual(TestBeds::mat_cmplx_1_3x3_rhs0_sol, 1e-14));

	Vector<Complex>    res_rhs = mat * vecSol;
	REQUIRE(true == res_rhs.IsEqual(TestBeds::mat_cmplx_1_3x3_rhs0, 1e-14));
}

// TODO - MED, dodati 5x5 i 8x8

///////////////////////////////////////////////////////////////////////////////////////////
/**********                         QR decomposition                              ********/
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Test_QRDecomposition_Solve_3_x_3", "[simple]")
{
	Matrix<Real>    mat = TestBeds::mat_3x3;
	Vector<Real> 	rhs = TestBeds::mat_3x3_rhs0;
	Vector<Real>	vecSol(rhs.size());

	QRDecompositionSolver qrSolver(mat);

	qrSolver.Solve(rhs, vecSol);
	REQUIRE(true == vecSol.IsEqual(TestBeds::mat_3x3_rhs0_sol, 1e-15));

	Vector<Real>    res_rhs = mat * vecSol;
	REQUIRE(true == res_rhs.IsEqual(TestBeds::mat_3x3_rhs0, 1e-15));
}

TEST_CASE("Test_QRDecomposition_Solve_5_x_5", "[simple]")
{
	Matrix<Real>    mat = TestBeds::mat_5x5;
	Vector<Real> 	rhs = TestBeds::mat_5x5_rhs0;
	Vector<Real>	vecSol(rhs.size());

	QRDecompositionSolver qrSolver(mat);

	qrSolver.Solve(rhs, vecSol);
	REQUIRE(true == vecSol.IsEqual(TestBeds::mat_5x5_rhs0_sol, 1e-13));

	Vector<Real>    res_rhs = mat * vecSol;
	REQUIRE(true == res_rhs.IsEqual(TestBeds::mat_5x5_rhs0, 1e-13));
}

TEST_CASE("Test_QRDecomposition_Solve_5_x_5_multi", "[simple]")
{
	Matrix<Real>     mat = TestBeds::mat_5x5;

	QRDecompositionSolver qrSolver(mat);

	for (int i = 0; i < 2; i++)
	{
		Vector<Real> 	rhs = TestBeds::mat_5x5_rhs_multi.VectorFromColumn(i);
		Vector<Real>	vecSol(rhs.size());

		qrSolver.Solve(rhs, vecSol);

		Vector<Real> 	rhs_sol = TestBeds::mat_5x5_rhs_multi_sol.VectorFromColumn(i);

		REQUIRE(true == vecSol.IsEqual(rhs_sol, 1e-10));

		Vector<Real>    res_rhs = mat * vecSol;

		REQUIRE(true == res_rhs.IsEqual(rhs, 1e-10));
	}
}

TEST_CASE("Test_QRDecomposition_Solve_8_x_8", "[simple]")
{
	Matrix<Real>    mat = TestBeds::mat_8x8;
	Vector<Real> 	rhs = TestBeds::mat_8x8_rhs0;
	Vector<Real>	vecSol(rhs.size());

	QRDecompositionSolver qrSolver(mat);

	qrSolver.Solve(rhs, vecSol);
	REQUIRE(true == vecSol.IsEqual(TestBeds::mat_8x8_rhs0_sol, 1e-14));

	Vector<Real>    res_rhs = mat * vecSol;
	REQUIRE(true == res_rhs.IsEqual(TestBeds::mat_8x8_rhs0, 1e-14));
}

TEST_CASE("Test_QRDecomposition_Solve_10_x_10", "[simple]")
{
	Matrix<Real>    mat = TestBeds::mat_10x10;
	Vector<Real> 	rhs = TestBeds::mat_10x10_rhs0;
	Vector<Real>	vecSol(rhs.size());

	QRDecompositionSolver qrSolver(mat);

	qrSolver.Solve(rhs, vecSol);
	REQUIRE(true == vecSol.IsEqual(TestBeds::mat_10x10_rhs0_sol, 1e-13));

	Vector<Real>    res_rhs = mat * vecSol;
	REQUIRE(true == res_rhs.IsEqual(TestBeds::mat_10x10_rhs0, 1e-14));
}

TEST_CASE("Test_QRDecomposition_Solve_20_x_20", "[simple]")
{
	Matrix<Real>    mat = TestBeds::mat_20x20;
	Vector<Real> 	rhs = TestBeds::mat_20x20_rhs0;
	Vector<Real>	vecSol(rhs.size());

	QRDecompositionSolver qrSolver(mat);

	qrSolver.Solve(rhs, vecSol);
	REQUIRE(true == vecSol.IsEqual(TestBeds::mat_20x20_rhs0_sol, 1e-11));

	Vector<Real>    res_rhs = mat * vecSol;
	REQUIRE(true == res_rhs.IsEqual(TestBeds::mat_20x20_rhs0, 1e-13));
}

TEST_CASE("Test_QRDecomposition_Solve_50_x_50", "[simple]")
{
	Matrix<Real>    mat = TestBeds::mat_50x50;
	Vector<Real> 	rhs = TestBeds::mat_50x50_rhs0;
	Vector<Real>	vecSol(rhs.size());

	QRDecompositionSolver qrSolver(mat);

	qrSolver.Solve(rhs, vecSol);
	REQUIRE(true == vecSol.IsEqual(TestBeds::mat_50x50_rhs0_sol, 1e-12));

	Vector<Real>    res_rhs = mat * vecSol;
	REQUIRE(true == res_rhs.IsEqual(TestBeds::mat_50x50_rhs0, 1e-13));
}

///////////////////////////////////////////////////////////////////////////////////////////
/**********                      Cholesky decomposition                            *******/
///////////////////////////////////////////////////////////////////////////////////////////
TEST_CASE("Test_CholeskyDecomposition_Solve_3_x_3", "[simple]")
{
    Matrix<Real>  mat{3, 3, { 2.0, -1.0,  0.0, 
                             -1.0,  3.0, -2.0, 
                              0.0, -2.0,  4.0 }};
    Vector<Real>     vecSol(3), rhs{1, 5, -2 }, solution{ 2.0, 3.0, 1.0 };
	
    CholeskyDecompositionSolver choleskySolver(mat);

	choleskySolver.Solve(rhs, vecSol);
	REQUIRE(true == vecSol.IsEqual(solution, 1e-14));

	Vector<Real>    res_rhs = mat * vecSol;
	REQUIRE(true == res_rhs.IsEqual(rhs, 1e-14));
}

///////////////////////////////////////////////////////////////////////////////////////////
/**********                         SVD decomposition                              *******/
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Test_SVDDecomposition_Solve_5_x_5", "[simple]")
{
	Matrix<Real>    mat = TestBeds::mat_5x5;
	Vector<Real> 	rhs = TestBeds::mat_5x5_rhs0;
	Vector<Real>	vecSol(rhs.size());

	SVDecompositionSolver svdSolver(mat);

	svdSolver.Solve(rhs, vecSol);
	REQUIRE(true == vecSol.IsEqual(TestBeds::mat_5x5_rhs0_sol, 1e-13));

	Vector<Real>    res_rhs = mat * vecSol;
	REQUIRE(true == res_rhs.IsEqual(TestBeds::mat_5x5_rhs0, 1e-13));
}

TEST_CASE("Test_SVDDecomposition_decomposition", "[simple]")
{
	Matrix<Real>    mat = TestBeds::mat_5x5;

	SVDecompositionSolver svdSolver(mat);

    Matrix<Real> u = svdSolver.getU();
    Matrix<Real> v = svdSolver.getV();
    Vector<Real> w = svdSolver.getW();

    Matrix<Real> wMat = MatrixUtils::DiagonalMatrixFromVector<Real>(w);
    
    Matrix<Real> b = u * wMat * v.GetTranspose();

    REQUIRE(true == b.IsEqual(mat, 1e-14));
}

///////////////////////////////////////////////////////////////////////////////////////////
/**********                         Complete test beds                             ********/
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Test_GaussJordan_COMPLETE_TEST_BED", "[simple]")
{
    for(int i=0; i<TestBeds::LinearAlgEqTestBed::numLinAlgEqSystems(); i++ )
    {
        TestBeds::TestLinearSystem testBed = TestBeds::LinearAlgEqTestBed::getLinAlgEqSystem(i);

        Matrix<Real>     mat = testBed._mat;
        Vector<Real>     rhs = testBed._rhs;

    	GaussJordanSolver<Real>::Solve(mat, rhs);
        REQUIRE(true == rhs.IsEqual(testBed._sol, 1e-13));

        Vector<Real>    res_rhs = testBed._mat * rhs;
        REQUIRE(true == res_rhs.IsEqual(testBed._rhs, 1e-13));        
    }
}

// TODO - LOW, Gauss Jordan solver COMPLEX COMPLETE TEST BED

TEST_CASE("Test_LUDecomposition_COMPLETE_TEST_BED", "[simple]")
{
    for(int i=0; i<TestBeds::LinearAlgEqTestBed::numLinAlgEqSystems(); i++ )
    {
        TestBeds::TestLinearSystem testBed = TestBeds::LinearAlgEqTestBed::getLinAlgEqSystem(i);

        Matrix<Real>     mat = testBed._mat;
        Vector<Real>     rhs = testBed._rhs;
        Vector<Real>     sol(rhs.size());

        LUDecompositionSolver<Real> luSolver(mat);

        luSolver.Solve(rhs, sol);
        REQUIRE(true == sol.IsEqual(testBed._sol, 1e-13));

        Vector<Real>    res_rhs = testBed._mat * sol;
        REQUIRE(true == res_rhs.IsEqual(testBed._rhs, 1e-13));
    }
}

// TODO - LOW, LU Decomposition COMPLEX COMPLETE TEST BED

TEST_CASE("Test_QRDecomposition_COMPLETE_TEST_BED", "[simple]")
{
    for(int i=0; i<TestBeds::LinearAlgEqTestBed::numLinAlgEqSystems(); i++ )
    {
        TestBeds::TestLinearSystem testBed = TestBeds::LinearAlgEqTestBed::getLinAlgEqSystem(i);

        Matrix<Real>     mat = testBed._mat;
        Vector<Real>     rhs = testBed._rhs;
        Vector<Real>     sol(rhs.size());

        QRDecompositionSolver luSolver(mat);

        luSolver.Solve(rhs, sol);
        REQUIRE(true == sol.IsEqual(testBed._sol, 1e-11));

        Vector<Real>    res_rhs = testBed._mat * sol;
        REQUIRE(true == res_rhs.IsEqual(testBed._rhs, 1e-12));
    }
}

TEST_CASE("Test_SVDDecomposition_COMPLETE_TEST_BED", "[simple]")
{
    for(int i=0; i<TestBeds::LinearAlgEqTestBed::numLinAlgEqSystems(); i++ )
    {
        TestBeds::TestLinearSystem testBed = TestBeds::LinearAlgEqTestBed::getLinAlgEqSystem(i);

        Matrix<Real>     mat = testBed._mat;
        Vector<Real>     rhs = testBed._rhs;
        Vector<Real>     sol(rhs.size());

        SVDecompositionSolver svdSolver(mat);

        svdSolver.Solve(rhs, sol);
        REQUIRE(true == sol.IsEqual(testBed._sol, 1e-12));

        Vector<Real>    res_rhs = testBed._mat * sol;
        REQUIRE(true == res_rhs.IsEqual(testBed._rhs, 1e-12));
    }
}
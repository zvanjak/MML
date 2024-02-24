#include "../catch/catch.hpp"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Vector.h"
#include "base/MatrixBandDiag.h"

#include "core/LinAlgEqSolvers.h"
#endif

using namespace MML;


// TODO 0.9 - Tridiag tests
namespace MML::Tests::TridiagMatrixTests
{
    TEST_CASE("MatrixTridiag_init_from_3_vectors", "[simple]") 
    {
        // initializing with 3 vectors
        TridiagonalMatrix<Real> a(4, {0.0, 4.5, 9.0, 10.0}, {4.0, 1.5, 6.0, 7.0}, {1.0, 2.0, 3.0, 0.0}  );

        // expected matrix
        Matrix<Real> b(4, 4, {  4.0, 1.0, 0.0, 0.0,
                                4.5, 1.5, 2.0, 0.0, 
                                0.0, 9.0, 6.0, 3.0,
                                0.0, 0.0, 10.0, 7.0 });    

        REQUIRE( 4 == a.RowNum());
        REQUIRE( 4 == a.ColNum());

        REQUIRE(a(0,0) == 4.0);
        REQUIRE(a(0,1) == 1.0);
        REQUIRE(a(1,0) == 4.5);
        REQUIRE(a(1,1) == 1.5);
        REQUIRE(a(1,2) == 2.0);
        REQUIRE(a(2,1) == 9.0);
        REQUIRE(a(2,2) == 6.0);
        REQUIRE(a(2,3) == 3.0);
        REQUIRE(a(3,2) == 10.0);
        REQUIRE(a(3,3) == 7.0);
        
        // these throw exception
        // REQUIRE(0.0 == a(0,2));
        // REQUIRE(a(0,3) == 0.0);
        // REQUIRE(a(1,3) == 0.0);
        // REQUIRE(a(2,0) == 0.0);
        // REQUIRE(a(3,0) == 0.0);
        // REQUIRE(a(3,1) == 0.0);
    }
    TEST_CASE("MatrixTridiag_init_from_initializer_list", "[simple]") 
    {
        // initializing with values in single initializer list
        TridiagonalMatrix<Real> a(4, {  4.0, 1.0,
                                        4.5, 1.5,  2.0, 
                                            9.0,  6.0, 3.0,
                                                10.0, 7.0 });
        
        // expected matrix
        Matrix<Real> b(4, 4, {  4.0, 1.0, 0.0, 0.0,
                                4.5, 1.5, 2.0, 0.0, 
                                0.0, 9.0, 6.0, 3.0,
                                0.0, 0.0, 10.0, 7.0 });    

        
        REQUIRE( 4 == a.RowNum());
        REQUIRE( 4 == a.ColNum());

        REQUIRE(a(0,0) == 4.0);
        REQUIRE(a(0,1) == 1.0);
        REQUIRE(a(1,0) == 4.5);
        REQUIRE(a(1,1) == 1.5);
        REQUIRE(a(1,2) == 2.0);
        REQUIRE(a(2,1) == 9.0);
        REQUIRE(a(2,2) == 6.0);
        REQUIRE(a(2,3) == 3.0);
        REQUIRE(a(3,2) == 10.0);
        REQUIRE(a(3,3) == 7.0);
    }

    TEST_CASE("MatrixTridiag_Solve", "[simple]") 
    {
        // initializing with values in single initializer list
        TridiagonalMatrix<Real> a(4, {  4.0, 1.0,
                                        4.5, 1.5,  2.0, 
                                            9.0,  6.0, 3.0,
                                                10.0, 7.0 });
        Matrix<Real> b(4, 4, {  4.0, 1.0, 0.0, 0.0,
                                4.5, 1.5, 2.0, 0.0, 
                                0.0, 9.0, 6.0, 3.0,
                                0.0, 0.0, 10.0, 7.0 });   
        Vector<Real> rhs{1.0, 2.0, 3.0, 4.0};

        Vector<Real> sol_a = a.Solve(rhs);

        REQUIRE( 4 == sol_a.size());
       
        LUDecompositionSolver<Real> solver(b); 
        Vector<Real> solLU = solver.Solve(rhs);     

        REQUIRE(sol_a[0] == Approx(solLU[0]));
        REQUIRE(sol_a[1] == Approx(solLU[1]));
        REQUIRE(sol_a[2] == Approx(solLU[2]));
        REQUIRE(sol_a[3] == Approx(solLU[3]));         
    }
}

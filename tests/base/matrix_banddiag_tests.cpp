#include "../catch/catch.hpp"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Vector.h"
#include "base/MatrixBandDiag.h"

#include "core/LinAlgEqSolvers.h"
#endif

using namespace MML;


namespace MML::Tests::BandDiagMatrixTests
{
    TEST_CASE("MatrixBandDiag_init_with_matrix", "[simple]") 
    {
        Matrix<Real> mat(7, 4,  { 0.0, 0.0, 3.0, 1.0,
                                0.0, 4.0, 1.0, 5.0,
                                9.0, 5.0, 6.0, 5.0,
                                4.0, 2.0, 2.0, 5.0,
                                5.0, 1.0, 2.0, 4.0,
                                3.0, 8.0, 4.0, 6.0,
                                2.0, 4.0, 4.0, 0.0 } );

        BandDiagonalMatrix band_diag(7, 2, 1, mat);

        // Expected layout
        // [      3,      1,      0,      0,      0,      0,      0,  ]
        // [      4,      1,      5,      0,      0,      0,      0,  ]
        // [      9,      5,      6,      5,      0,      0,      0,  ]
        // [      0,      4,      2,      2,      5,      0,      0,  ]
        // [      0,      0,      5,      1,      2,      4,      0,  ]
        // [      0,      0,      0,      3,      8,      4,      6,  ]
        // [      0,      0,      0,      0,      2,      4,      4,  ]   

        REQUIRE( 7 == band_diag.RowNum());
        REQUIRE( 7 == band_diag.ColNum());  

        REQUIRE(band_diag(0,0) == 3.0);
        REQUIRE(band_diag(0,1) == 1.0);
        
        REQUIRE(band_diag(1,0) == 4.0);
        REQUIRE(band_diag(1,1) == 1.0);
        REQUIRE(band_diag(1,2) == 5.0);

        REQUIRE(band_diag(2,0) == 9.0);
        REQUIRE(band_diag(2,1) == 5.0);
        REQUIRE(band_diag(2,2) == 6.0);

        REQUIRE(band_diag(3,1) == 4.0);
        REQUIRE(band_diag(3,2) == 2.0);
        REQUIRE(band_diag(3,3) == 2.0);
        REQUIRE(band_diag(3,4) == 5.0);

        REQUIRE(band_diag(4,2) == 5.0);
        REQUIRE(band_diag(4,3) == 1.0);
        REQUIRE(band_diag(4,4) == 2.0);
        REQUIRE(band_diag(4,5) == 4.0);

        REQUIRE(band_diag(5,3) == 3.0);
        REQUIRE(band_diag(5,4) == 8.0);
        REQUIRE(band_diag(5,5) == 4.0);
        REQUIRE(band_diag(5,6) == 6.0);

        REQUIRE(band_diag(6,4) == 2.0);
        REQUIRE(band_diag(6,5) == 4.0);
        REQUIRE(band_diag(6,6) == 4.0);
    }

    TEST_CASE("MatrixBandDiag_Solve", "[simple]") 
    {
        // matrix with 2 sub-diagonals and 1 super-diagonal
        Matrix<Real> mat(7, 4,  { 0.0, 0.0, 3.0, 1.0,
                                0.0, 4.0, 1.0, 5.0,
                                9.0, 5.0, 6.0, 5.0,
                                4.0, 2.0, 2.0, 5.0,
                                5.0, 1.0, 2.0, 4.0,
                                3.0, 8.0, 4.0, 6.0,
                                2.0, 4.0, 4.0, 0.0 } );
        BandDiagonalMatrix band_diag(7, 2, 1, mat);

        Matrix<Real> b(7, 7, {  3,      1,      0,      0,      0,      0,      0,  
                                4,      1,      5,      0,      0,      0,      0,  
                                9,      5,      6,      5,      0,      0,      0,  
                                0,      4,      2,      2,      5,      0,      0,  
                                0,      0,      5,      1,      2,      4,      0,  
                                0,      0,      0,      3,      8,      4,      6,  
                                0,      0,      0,      0,      2,      4,      4});  

        Vector<Real> rhs{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0};

        BandDiagLUSolver solverBand(mat, 2, 1);
        Vector<Real> solBand = solverBand.Solve(rhs);
       
        LUDecompositionSolver<Real> solverLU(b); 
        Vector<Real> solLU = solverLU.Solve(rhs);     

        REQUIRE( 7 == solBand.size());
        REQUIRE( 7 == solLU.size());
        
        REQUIRE( solBand.IsEqual(solLU,1e-15) );
    }

    TEST_CASE("MatrixBandDiag_operator*(Vector<Real> &)", "[simple]") 
    {
        // matrix with 2 sub-diagonals and 1 super-diagonal
        Matrix<Real> mat(7, 4,  { 0.0, 0.0, 3.0, 1.0,
                                0.0, 4.0, 1.0, 5.0,
                                9.0, 5.0, 6.0, 5.0,
                                4.0, 2.0, 2.0, 5.0,
                                5.0, 1.0, 2.0, 4.0,
                                3.0, 8.0, 4.0, 6.0,
                                2.0, 4.0, 4.0, 0.0 } );
        BandDiagonalMatrix band_diag(7, 2, 1, mat);

        Matrix<Real> mat_reg(7, 7, {  3,      1,      0,      0,      0,      0,      0,  
                                    4,      1,      5,      0,      0,      0,      0,  
                                    9,      5,      6,      5,      0,      0,      0,  
                                    0,      4,      2,      2,      5,      0,      0,  
                                    0,      0,      5,      1,      2,      4,      0,  
                                    0,      0,      0,      3,      8,      4,      6,  
                                    0,      0,      0,      0,      2,      4,      4});  

        Vector<Real> b{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0};

        Vector<Real> mulBand = band_diag * b;
        Vector<Real> mulReg  = mat_reg * b;

        REQUIRE( mulBand == mulReg );
    }
}

// Test that MML.h single header works standalone
// This file explicitly includes ONLY MML.h - no other MML headers

#include <catch2/catch_all.hpp>

// Force single header mode and include it directly
#define MML_USE_SINGLE_HEADER
#include "single_header/MML.h"

using namespace MML;

namespace MML::Tests::SingleHeader
{
    TEST_CASE("SingleHeader_Vector_Basic", "[single_header]")
    {
        Vector<Real> v1({ REAL(1.0), REAL(2.0), REAL(3.0) });
        Vector<Real> v2({ REAL(4.0), REAL(5.0), REAL(6.0) });
        
        // Vector addition
        auto v3 = v1 + v2;
        REQUIRE(v3.size() == 3);
        REQUIRE(v3[0] == REAL(5.0));
        REQUIRE(v3[1] == REAL(7.0));
        REQUIRE(v3[2] == REAL(9.0));
        
        // Scalar product (free function in VectorUtils)
        Real dot = Utils::ScalarProduct(v1, v2);
        REQUIRE(dot == REAL(32.0)); // 1*4 + 2*5 + 3*6 = 4 + 10 + 18 = 32
    }
    
    TEST_CASE("SingleHeader_Matrix_Basic", "[single_header]")
    {
        Matrix<Real> A(2, 3, {
            REAL(1.0), REAL(2.0), REAL(3.0),
            REAL(4.0), REAL(5.0), REAL(6.0)
        });
        
        Matrix<Real> B(3, 2, {
            REAL(7.0), REAL(8.0),
            REAL(9.0), REAL(10.0),
            REAL(11.0), REAL(12.0)
        });
        
        // Matrix multiplication: (2x3) * (3x2) = (2x2)
        auto C = A * B;
        
        REQUIRE(C.RowNum() == 2);
        REQUIRE(C.ColNum() == 2);
        
        // C[0][0] = 1*7 + 2*9 + 3*11 = 7 + 18 + 33 = 58
        REQUIRE(C[0][0] == REAL(58.0));
        // C[0][1] = 1*8 + 2*10 + 3*12 = 8 + 20 + 36 = 64
        REQUIRE(C[0][1] == REAL(64.0));
        // C[1][0] = 4*7 + 5*9 + 6*11 = 28 + 45 + 66 = 139
        REQUIRE(C[1][0] == REAL(139.0));
        // C[1][1] = 4*8 + 5*10 + 6*12 = 32 + 50 + 72 = 154
        REQUIRE(C[1][1] == REAL(154.0));
    }
    
    TEST_CASE("SingleHeader_MatrixVector_Multiply", "[single_header]")
    {
        Matrix<Real> A(2, 3, {
            REAL(1.0), REAL(2.0), REAL(3.0),
            REAL(4.0), REAL(5.0), REAL(6.0)
        });
        
        Vector<Real> v({ REAL(1.0), REAL(2.0), REAL(3.0) });
        
        // Matrix-vector multiplication: (2x3) * (3x1) = (2x1)
        auto result = A * v;
        
        REQUIRE(result.size() == 2);
        // result[0] = 1*1 + 2*2 + 3*3 = 1 + 4 + 9 = 14
        REQUIRE(result[0] == REAL(14.0));
        // result[1] = 4*1 + 5*2 + 6*3 = 4 + 10 + 18 = 32
        REQUIRE(result[1] == REAL(32.0));
    }
}

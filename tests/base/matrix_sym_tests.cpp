#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Vector.h"
#include "base/Matrix.h"
#include "base/MatrixSym.h"
#endif

using namespace MML;
using namespace MML::Testing;

namespace MML::Tests::Base::MatrixSymTests
{
    //=============================================================================
    // SECTION 1: CONSTRUCTORS (12 tests)
    //=============================================================================

    TEST_CASE("MatrixSym_default_ctor_creates_empty", "[MatrixSym][constructors]")
    {
        MatrixSym<Real> a;
        
        REQUIRE(0 == a.Dim());
        REQUIRE(0 == a.RowNum());
        REQUIRE(0 == a.ColNum());
        REQUIRE(0 == a.size());
        REQUIRE(a.empty());
    }

    TEST_CASE("MatrixSym_dimension_ctor_init_to_zero", "[MatrixSym][constructors]")
    {
        MatrixSym<Real> a(3);

        REQUIRE(3 == a.Dim());
        REQUIRE(3 == a.RowNum());
        REQUIRE(3 == a.ColNum());
        REQUIRE(6 == a.size());  // 3*(3+1)/2 = 6
        REQUIRE_FALSE(a.empty());

        // Verify all elements are zero
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                REQUIRE(REAL(0.0) == a(i, j));
            }
        }
    }

    TEST_CASE("MatrixSym_dimension_ctor_init_to_value", "[MatrixSym][constructors]")
    {
        MatrixSym<Real> a(2, REAL(5.0));

        REQUIRE(2 == a.Dim());
        REQUIRE(3 == a.size());  // 2*(2+1)/2 = 3

        REQUIRE(REAL(5.0) == a(0, 0));
        REQUIRE(REAL(5.0) == a(0, 1));
        REQUIRE(REAL(5.0) == a(1, 0));
        REQUIRE(REAL(5.0) == a(1, 1));
    }

    TEST_CASE("MatrixSym_initializer_list_2x2", "[MatrixSym][constructors]")
    {
        // Lower triangular: a00, a10, a11
        MatrixSym<Real> a(2, { REAL(1.0), REAL(2.0), REAL(3.0) });

        REQUIRE(2 == a.Dim());
        REQUIRE(REAL(1.0) == a(0, 0));
        REQUIRE(REAL(2.0) == a(0, 1));  // Symmetric: a(0,1) = a(1,0)
        REQUIRE(REAL(2.0) == a(1, 0));
        REQUIRE(REAL(3.0) == a(1, 1));
    }

    TEST_CASE("MatrixSym_initializer_list_3x3", "[MatrixSym][constructors]")
    {
        // Lower triangular: a00, a10, a11, a20, a21, a22
        MatrixSym<Real> a(3, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0), REAL(6.0) });

        REQUIRE(3 == a.Dim());
        REQUIRE(6 == a.size());

        // Verify symmetry
        REQUIRE(REAL(1.0) == a(0, 0));
        REQUIRE(REAL(2.0) == a(1, 0));
        REQUIRE(REAL(2.0) == a(0, 1));  // Symmetric
        REQUIRE(REAL(3.0) == a(1, 1));
        REQUIRE(REAL(4.0) == a(2, 0));
        REQUIRE(REAL(4.0) == a(0, 2));  // Symmetric
        REQUIRE(REAL(5.0) == a(2, 1));
        REQUIRE(REAL(5.0) == a(1, 2));  // Symmetric
        REQUIRE(REAL(6.0) == a(2, 2));
    }

    TEST_CASE("MatrixSym_initializer_list_wrong_size_throws", "[MatrixSym][constructors]")
    {
        // Too many elements
        REQUIRE_THROWS_AS(MatrixSym<Real>(2, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) }), MatrixDimensionError);
        
        // Too few elements
        REQUIRE_THROWS_AS(MatrixSym<Real>(3, { REAL(1.0), REAL(2.0), REAL(3.0) }), MatrixDimensionError);
    }

    TEST_CASE("MatrixSym_copy_constructor", "[MatrixSym][constructors]")
    {
        MatrixSym<Real> a(3, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0), REAL(6.0) });
        MatrixSym<Real> b(a);

        REQUIRE(a.Dim() == b.Dim());
        REQUIRE(a.IsEqual(b));

        // Verify deep copy - modifying b doesn't affect a
        b(0, 0) = REAL(999.0);
        REQUIRE(REAL(1.0) == a(0, 0));
        REQUIRE(REAL(999.0) == b(0, 0));
    }

    TEST_CASE("MatrixSym_move_constructor", "[MatrixSym][constructors]")
    {
        MatrixSym<Real> a(3, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0), REAL(6.0) });
        MatrixSym<Real> b(std::move(a));

        REQUIRE(3 == b.Dim());
        REQUIRE(REAL(6.0) == b(2, 2));
        
        // Source should be empty after move
        REQUIRE(0 == a.Dim());
        REQUIRE(a.empty());
    }

    TEST_CASE("MatrixSym_copy_empty_matrix", "[MatrixSym][constructors]")
    {
        MatrixSym<Real> a;
        MatrixSym<Real> b(a);

        REQUIRE(b.empty());
        REQUIRE(0 == b.Dim());
    }

    TEST_CASE("MatrixSym_negative_dimension_throws", "[MatrixSym][constructors]")
    {
        REQUIRE_THROWS_AS(MatrixSym<Real>(-1), MatrixDimensionError);
    }

    TEST_CASE("MatrixSym_1x1_matrix", "[MatrixSym][constructors][edge-cases]")
    {
        MatrixSym<Real> a(1, { REAL(42.0) });

        REQUIRE(1 == a.Dim());
        REQUIRE(1 == a.size());
        REQUIRE(REAL(42.0) == a(0, 0));
    }

    //=============================================================================
    // SECTION 2: ASSIGNMENT OPERATORS (6 tests)
    //=============================================================================

    TEST_CASE("MatrixSym_copy_assignment_same_size", "[MatrixSym][assignment]")
    {
        MatrixSym<Real> a(2, { REAL(1.0), REAL(2.0), REAL(3.0) });
        MatrixSym<Real> b(2, { REAL(4.0), REAL(5.0), REAL(6.0) });

        b = a;

        REQUIRE(a.IsEqual(b));
        
        // Verify deep copy
        b(0, 0) = REAL(999.0);
        REQUIRE(REAL(1.0) == a(0, 0));
    }

    TEST_CASE("MatrixSym_copy_assignment_different_size", "[MatrixSym][assignment]")
    {
        MatrixSym<Real> a(3, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0), REAL(6.0) });
        MatrixSym<Real> b(2, { REAL(7.0), REAL(8.0), REAL(9.0) });

        b = a;

        REQUIRE(3 == b.Dim());
        REQUIRE(a.IsEqual(b));
    }

    TEST_CASE("MatrixSym_move_assignment", "[MatrixSym][assignment]")
    {
        MatrixSym<Real> a(3, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0), REAL(6.0) });
        MatrixSym<Real> b(2);

        b = std::move(a);

        REQUIRE(3 == b.Dim());
        REQUIRE(REAL(6.0) == b(2, 2));
        REQUIRE(a.empty());
    }

    TEST_CASE("MatrixSym_self_assignment", "[MatrixSym][assignment]")
    {
        MatrixSym<Real> a(2, { REAL(1.0), REAL(2.0), REAL(3.0) });
        
        a = a;  // Self-assignment

        REQUIRE(2 == a.Dim());
        REQUIRE(REAL(1.0) == a(0, 0));
        REQUIRE(REAL(2.0) == a(1, 0));
        REQUIRE(REAL(3.0) == a(1, 1));
    }

    TEST_CASE("MatrixSym_move_self_assignment", "[MatrixSym][assignment]")
    {
        MatrixSym<Real> a(2, { REAL(1.0), REAL(2.0), REAL(3.0) });
        
        a = std::move(a);  // Self move-assignment

        REQUIRE(2 == a.Dim());
        REQUIRE(REAL(1.0) == a(0, 0));
    }

    TEST_CASE("MatrixSym_assign_empty_to_non_empty", "[MatrixSym][assignment]")
    {
        MatrixSym<Real> a;
        MatrixSym<Real> b(3, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0), REAL(6.0) });

        b = a;

        REQUIRE(b.empty());
        REQUIRE(0 == b.Dim());
    }

    //=============================================================================
    // SECTION 3: ACCESS OPERATORS (8 tests)
    //=============================================================================

    TEST_CASE("MatrixSym_access_operator_read", "[MatrixSym][access]")
    {
        const MatrixSym<Real> a(2, { REAL(1.0), REAL(2.0), REAL(3.0) });

        REQUIRE(REAL(1.0) == a(0, 0));
        REQUIRE(REAL(2.0) == a(0, 1));
        REQUIRE(REAL(2.0) == a(1, 0));
        REQUIRE(REAL(3.0) == a(1, 1));
    }

    TEST_CASE("MatrixSym_access_operator_write", "[MatrixSym][access]")
    {
        MatrixSym<Real> a(2);

        a(0, 0) = REAL(1.0);
        a(0, 1) = REAL(2.0);  // Writing to (0,1) should affect (1,0) due to symmetry
        a(1, 1) = REAL(3.0);

        REQUIRE(REAL(1.0) == a(0, 0));
        REQUIRE(REAL(2.0) == a(0, 1));
        REQUIRE(REAL(2.0) == a(1, 0));  // Same as (0,1)
        REQUIRE(REAL(3.0) == a(1, 1));
    }

    TEST_CASE("MatrixSym_symmetry_verification", "[MatrixSym][access]")
    {
        MatrixSym<Real> a(4, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0), REAL(6.0), REAL(7.0), REAL(8.0), REAL(9.0), REAL(10.0) });

        // Verify a(i,j) == a(j,i) for all i,j
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                REQUIRE(a(i, j) == a(j, i));
            }
        }
    }

    TEST_CASE("MatrixSym_at_bounds_checking", "[MatrixSym][access]")
    {
        MatrixSym<Real> a(3, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0), REAL(6.0) });

        // Valid access
        REQUIRE(REAL(1.0) == a.at(0, 0));
        REQUIRE(REAL(6.0) == a.at(2, 2));

        // Out of bounds
        REQUIRE_THROWS_AS(a.at(-1, 0), MatrixAccessBoundsError);
        REQUIRE_THROWS_AS(a.at(0, -1), MatrixAccessBoundsError);
        REQUIRE_THROWS_AS(a.at(3, 0), MatrixAccessBoundsError);
        REQUIRE_THROWS_AS(a.at(0, 3), MatrixAccessBoundsError);
    }

    TEST_CASE("MatrixSym_ElemAt_bounds_checking", "[MatrixSym][access]")
    {
        MatrixSym<Real> a(2, { REAL(1.0), REAL(2.0), REAL(3.0) });

        REQUIRE(REAL(1.0) == a.ElemAt(0, 0));
        REQUIRE(REAL(2.0) == a.ElemAt(0, 1));

        REQUIRE_THROWS_AS(a.ElemAt(2, 0), MatrixAccessBoundsError);
    }

    TEST_CASE("MatrixSym_data_accessor", "[MatrixSym][access]")
    {
        MatrixSym<Real> a(2, { REAL(1.0), REAL(2.0), REAL(3.0) });

        Real* ptr = a.data();
        REQUIRE(ptr != nullptr);

        // Data stored as lower triangular: [a00, a10, a11]
        REQUIRE(REAL(1.0) == ptr[0]);
        REQUIRE(REAL(2.0) == ptr[1]);
        REQUIRE(REAL(3.0) == ptr[2]);

        // Modify through raw pointer
        ptr[0] = REAL(10.0);
        REQUIRE(REAL(10.0) == a(0, 0));
    }

    TEST_CASE("MatrixSym_const_data_accessor", "[MatrixSym][access]")
    {
        const MatrixSym<Real> a(2, { REAL(1.0), REAL(2.0), REAL(3.0) });

        const Real* ptr = a.data();
        REQUIRE(ptr != nullptr);
        REQUIRE(REAL(1.0) == ptr[0]);
    }

    TEST_CASE("MatrixSym_empty_data_accessor", "[MatrixSym][access]")
    {
        MatrixSym<Real> a;
        
        // data() on empty should return nullptr or be safe
        Real* ptr = a.data();
        REQUIRE(a.size() == 0);
    }

    //=============================================================================
    // SECTION 4: ARITHMETIC OPERATORS (15 tests)
    //=============================================================================

    TEST_CASE("MatrixSym_addition", "[MatrixSym][arithmetic]")
    {
        MatrixSym<Real> a(2, { REAL(1.0), REAL(2.0), REAL(3.0) });
        MatrixSym<Real> b(2, { REAL(4.0), REAL(5.0), REAL(6.0) });

        auto c = a + b;

        REQUIRE(2 == c.Dim());
        REQUIRE(REAL(5.0) == c(0, 0));
        REQUIRE(REAL(7.0) == c(0, 1));
        REQUIRE(REAL(9.0) == c(1, 1));
    }

    TEST_CASE("MatrixSym_addition_compound", "[MatrixSym][arithmetic]")
    {
        MatrixSym<Real> a(2, { REAL(1.0), REAL(2.0), REAL(3.0) });
        MatrixSym<Real> b(2, { REAL(4.0), REAL(5.0), REAL(6.0) });

        a += b;

        REQUIRE(REAL(5.0) == a(0, 0));
        REQUIRE(REAL(7.0) == a(0, 1));
        REQUIRE(REAL(9.0) == a(1, 1));
    }

    TEST_CASE("MatrixSym_subtraction", "[MatrixSym][arithmetic]")
    {
        MatrixSym<Real> a(2, { REAL(5.0), REAL(7.0), REAL(9.0) });
        MatrixSym<Real> b(2, { REAL(1.0), REAL(2.0), REAL(3.0) });

        auto c = a - b;

        REQUIRE(REAL(4.0) == c(0, 0));
        REQUIRE(REAL(5.0) == c(0, 1));
        REQUIRE(REAL(6.0) == c(1, 1));
    }

    TEST_CASE("MatrixSym_subtraction_compound", "[MatrixSym][arithmetic]")
    {
        MatrixSym<Real> a(2, { REAL(5.0), REAL(7.0), REAL(9.0) });
        MatrixSym<Real> b(2, { REAL(1.0), REAL(2.0), REAL(3.0) });

        a -= b;

        REQUIRE(REAL(4.0) == a(0, 0));
        REQUIRE(REAL(5.0) == a(0, 1));
    }

    TEST_CASE("MatrixSym_unary_negation", "[MatrixSym][arithmetic]")
    {
        MatrixSym<Real> a(2, { REAL(1.0), -REAL(2.0), REAL(3.0) });

        auto b = -a;

        REQUIRE(-REAL(1.0) == b(0, 0));
        REQUIRE(REAL(2.0) == b(0, 1));
        REQUIRE(-REAL(3.0) == b(1, 1));
    }

    TEST_CASE("MatrixSym_matrix_multiplication", "[MatrixSym][arithmetic]")
    {
        // a = [1 2; 2 3], b = [1 2; 2 3]
        // a*b = [1*1+2*2  1*2+2*3; 2*1+3*2  2*2+3*3] = [5 8; 8 13]
        MatrixSym<Real> a(2, { REAL(1.0), REAL(2.0), REAL(3.0) });
        MatrixSym<Real> b(2, { REAL(1.0), REAL(2.0), REAL(3.0) });

        Matrix<Real> c = a * b;

        REQUIRE(REAL(5.0) == c(0, 0));
        REQUIRE(REAL(8.0) == c(0, 1));
        REQUIRE(REAL(8.0) == c(1, 0));
        REQUIRE(REAL(13.0) == c(1, 1));
    }

    TEST_CASE("MatrixSym_times_Matrix", "[MatrixSym][arithmetic]")
    {
        MatrixSym<Real> a(2, { REAL(1.0), REAL(2.0), REAL(3.0) });
        Matrix<Real> b(2, 3, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0), REAL(6.0) });

        Matrix<Real> c = a * b;

        REQUIRE(2 == c.RowNum());
        REQUIRE(3 == c.ColNum());

        // c[0][0] = 1*1 + 2*4 = 9
        REQUIRE(REAL(9.0) == c(0, 0));
        // c[0][1] = 1*2 + 2*5 = 12
        REQUIRE(REAL(12.0) == c(0, 1));
    }

    TEST_CASE("MatrixSym_scalar_multiplication_right", "[MatrixSym][arithmetic]")
    {
        MatrixSym<Real> a(2, { REAL(1.0), REAL(2.0), REAL(3.0) });

        auto b = a * REAL(2.0);

        REQUIRE(REAL(2.0) == b(0, 0));
        REQUIRE(REAL(4.0) == b(0, 1));
        REQUIRE(REAL(6.0) == b(1, 1));
    }

    TEST_CASE("MatrixSym_scalar_multiplication_left", "[MatrixSym][arithmetic]")
    {
        MatrixSym<Real> a(2, { REAL(1.0), REAL(2.0), REAL(3.0) });

        auto b = REAL(3.0) * a;

        REQUIRE(REAL(3.0) == b(0, 0));
        REQUIRE(REAL(6.0) == b(0, 1));
        REQUIRE(REAL(9.0) == b(1, 1));
    }

    TEST_CASE("MatrixSym_scalar_multiplication_compound", "[MatrixSym][arithmetic]")
    {
        MatrixSym<Real> a(2, { REAL(1.0), REAL(2.0), REAL(3.0) });

        a *= REAL(2.0);

        REQUIRE(REAL(2.0) == a(0, 0));
        REQUIRE(REAL(4.0) == a(0, 1));
        REQUIRE(REAL(6.0) == a(1, 1));
    }

    TEST_CASE("MatrixSym_scalar_division", "[MatrixSym][arithmetic]")
    {
        MatrixSym<Real> a(2, { REAL(2.0), REAL(4.0), REAL(6.0) });

        auto b = a / REAL(2.0);

        REQUIRE(REAL(1.0) == b(0, 0));
        REQUIRE(REAL(2.0) == b(0, 1));
        REQUIRE(REAL(3.0) == b(1, 1));
    }

    TEST_CASE("MatrixSym_scalar_division_compound", "[MatrixSym][arithmetic]")
    {
        MatrixSym<Real> a(2, { REAL(2.0), REAL(4.0), REAL(6.0) });

        a /= REAL(2.0);

        REQUIRE(REAL(1.0) == a(0, 0));
        REQUIRE(REAL(2.0) == a(0, 1));
    }

    TEST_CASE("MatrixSym_vector_multiplication_right", "[MatrixSym][arithmetic]")
    {
        // a = [1 10; 10 5], v = [1; 2]
        // a*v = [1*1+10*2; 10*1+5*2] = [21; 20]
        MatrixSym<Real> a(2, { REAL(1.0), REAL(10.0), REAL(5.0) });
        Vector<Real> v({ REAL(1.0), REAL(2.0) });

        Vector<Real> c = a * v;

        REQUIRE(REAL(21.0) == c[0]);
        REQUIRE(REAL(20.0) == c[1]);
    }

    TEST_CASE("MatrixSym_vector_multiplication_left", "[MatrixSym][arithmetic]")
    {
        MatrixSym<Real> a(2, { REAL(1.0), REAL(10.0), REAL(5.0) });
        Vector<Real> v({ REAL(1.0), REAL(2.0) });

        Vector<Real> c = v * a;

        REQUIRE(REAL(21.0) == c[0]);
        REQUIRE(REAL(20.0) == c[1]);
    }

    TEST_CASE("MatrixSym_dimension_mismatch_throws", "[MatrixSym][arithmetic]")
    {
        MatrixSym<Real> a(2);
        MatrixSym<Real> b(3);

        REQUIRE_THROWS_AS(a + b, MatrixDimensionError);
        REQUIRE_THROWS_AS(a - b, MatrixDimensionError);
        REQUIRE_THROWS_AS(a * b, MatrixDimensionError);
    }

    //=============================================================================
    // SECTION 5: PROPERTIES AND UTILITIES (12 tests)
    //=============================================================================

    TEST_CASE("MatrixSym_Trace", "[MatrixSym][properties]")
    {
        MatrixSym<Real> a(3, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0), REAL(6.0) });

        REQUIRE(REAL(10.0) == a.Trace());  // 1 + 3 + 6 = 10
    }

    TEST_CASE("MatrixSym_Trace_empty", "[MatrixSym][properties]")
    {
        MatrixSym<Real> a;
        REQUIRE(REAL(0.0) == a.Trace());
    }

    TEST_CASE("MatrixSym_NormFrobenius", "[MatrixSym][properties]")
    {
        // a = [1 2; 2 3]
        // ||a||_F = sqrt(1^2 + 2*2^2 + 3^2) = sqrt(1 + 8 + 9) = sqrt(18)
        MatrixSym<Real> a(2, { REAL(1.0), REAL(2.0), REAL(3.0) });

        Real norm = a.NormFrobenius();
        REQUIRE_THAT(norm, RealApprox(std::sqrt(REAL(18.0))));
    }

    TEST_CASE("MatrixSym_NormInf", "[MatrixSym][properties]")
    {
        // a = [1 2; 2 4]
        // row 0 sum: |1| + |2| = 3
        // row 1 sum: |2| + |4| = 6
        MatrixSym<Real> a(2, { REAL(1.0), REAL(2.0), REAL(4.0) });

        REQUIRE(REAL(6.0) == a.NormInf());
    }

    TEST_CASE("MatrixSym_Norm1_equals_NormInf_for_symmetric", "[MatrixSym][properties]")
    {
        MatrixSym<Real> a(3, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0), REAL(6.0) });

        REQUIRE(a.Norm1() == a.NormInf());
    }

    TEST_CASE("MatrixSym_GetAsMatrix", "[MatrixSym][properties]")
    {
        MatrixSym<Real> a(2, { REAL(1.0), REAL(2.0), REAL(3.0) });
        
        Matrix<Real> m = a.GetAsMatrix();

        REQUIRE(2 == m.RowNum());
        REQUIRE(2 == m.ColNum());
        REQUIRE(REAL(1.0) == m(0, 0));
        REQUIRE(REAL(2.0) == m(0, 1));
        REQUIRE(REAL(2.0) == m(1, 0));
        REQUIRE(REAL(3.0) == m(1, 1));
    }

    TEST_CASE("MatrixSym_VectorFromRow", "[MatrixSym][properties]")
    {
        MatrixSym<Real> a(3, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0), REAL(6.0) });

        Vector<Real> row1 = a.VectorFromRow(1);

        REQUIRE(3 == row1.size());
        REQUIRE(REAL(2.0) == row1[0]);  // a(1,0)
        REQUIRE(REAL(3.0) == row1[1]);  // a(1,1)
        REQUIRE(REAL(5.0) == row1[2]);  // a(1,2)
    }

    TEST_CASE("MatrixSym_VectorFromColumn", "[MatrixSym][properties]")
    {
        MatrixSym<Real> a(3, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0), REAL(6.0) });

        Vector<Real> col1 = a.VectorFromColumn(1);

        REQUIRE(3 == col1.size());
        REQUIRE(REAL(2.0) == col1[0]);  // a(0,1)
        REQUIRE(REAL(3.0) == col1[1]);  // a(1,1)
        REQUIRE(REAL(5.0) == col1[2]);  // a(2,1)
    }

    TEST_CASE("MatrixSym_VectorFromDiagonal", "[MatrixSym][properties]")
    {
        MatrixSym<Real> a(3, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0), REAL(6.0) });

        Vector<Real> diag = a.VectorFromDiagonal();

        REQUIRE(3 == diag.size());
        REQUIRE(REAL(1.0) == diag[0]);
        REQUIRE(REAL(3.0) == diag[1]);
        REQUIRE(REAL(6.0) == diag[2]);
    }

    TEST_CASE("MatrixSym_GetInverse", "[MatrixSym][properties]")
    {
        MatrixSym<Real> a(2, { REAL(1.0), REAL(2.0), REAL(3.0) });

        Matrix<Real> inv = a.GetInverse();
        Matrix<Real> product = a.GetAsMatrix() * inv;

        REQUIRE(product.IsEqualTo(Matrix<Real>::GetUnitMatrix(2)));
    }

    TEST_CASE("MatrixSym_IsEqual", "[MatrixSym][properties]")
    {
        MatrixSym<Real> a(2, { REAL(1.0), REAL(2.0), REAL(3.0) });
        MatrixSym<Real> b(2, { REAL(1.0), REAL(2.0), REAL(3.0) });
        MatrixSym<Real> c(2, { REAL(1.0), REAL(2.0), REAL(4.0) });

        REQUIRE(a.IsEqual(b));
        REQUIRE_FALSE(a.IsEqual(c));
    }

    TEST_CASE("MatrixSym_IsEqual_with_tolerance", "[MatrixSym][properties]")
    {
        MatrixSym<Real> a(2, { REAL(1.0), REAL(2.0), REAL(4.0) });
        MatrixSym<Real> b(2, { REAL(1.0), REAL(2.0), REAL(4.0001) });

        REQUIRE(a.IsEqual(b, 1e-3));
        REQUIRE_FALSE(a.IsEqual(b, 1e-5));
    }

    //=============================================================================
    // SECTION 6: STATIC FACTORY METHODS (6 tests)
    //=============================================================================

    TEST_CASE("MatrixSym_GetUnitMatrix", "[MatrixSym][static-factory]")
    {
        auto I = MatrixSym<Real>::GetUnitMatrix(3);

        REQUIRE(3 == I.Dim());
        
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                if (i == j) {
                    REQUIRE(REAL(1.0) == I(i, j));
                } else {
                    REQUIRE(REAL(0.0) == I(i, j));
                }
            }
        }
    }

    TEST_CASE("MatrixSym_GetDiagonalMatrix", "[MatrixSym][static-factory]")
    {
        Vector<Real> diag({ REAL(1.0), REAL(2.0), REAL(3.0) });
        auto D = MatrixSym<Real>::GetDiagonalMatrix(diag);

        REQUIRE(3 == D.Dim());
        REQUIRE(REAL(1.0) == D(0, 0));
        REQUIRE(REAL(0.0) == D(0, 1));
        REQUIRE(REAL(2.0) == D(1, 1));
        REQUIRE(REAL(3.0) == D(2, 2));
    }

    TEST_CASE("MatrixSym_FromLower", "[MatrixSym][static-factory]")
    {
        Matrix<Real> m(3, 3, {
            REAL(1.0), REAL(99.0), REAL(99.0),
            REAL(2.0), REAL(3.0), REAL(99.0),
            REAL(4.0), REAL(5.0), REAL(6.0)
        });

        auto s = MatrixSym<Real>::FromLower(m);

        REQUIRE(3 == s.Dim());
        REQUIRE(REAL(1.0) == s(0, 0));
        REQUIRE(REAL(2.0) == s(1, 0));
        REQUIRE(REAL(3.0) == s(1, 1));
        REQUIRE(REAL(4.0) == s(2, 0));
        REQUIRE(REAL(5.0) == s(2, 1));
        REQUIRE(REAL(6.0) == s(2, 2));
    }

    TEST_CASE("MatrixSym_FromUpper", "[MatrixSym][static-factory]")
    {
        Matrix<Real> m(3, 3, {
            REAL(1.0), REAL(2.0), REAL(4.0),
            REAL(99.0), REAL(3.0), REAL(5.0),
            REAL(99.0), REAL(99.0), REAL(6.0)
        });

        auto s = MatrixSym<Real>::FromUpper(m);

        REQUIRE(3 == s.Dim());
        REQUIRE(REAL(1.0) == s(0, 0));
        REQUIRE(REAL(2.0) == s(0, 1));
        REQUIRE(REAL(2.0) == s(1, 0));  // Symmetric
        REQUIRE(REAL(3.0) == s(1, 1));
    }

    TEST_CASE("MatrixSym_FromFullMatrix", "[MatrixSym][static-factory]")
    {
        // Non-symmetric matrix - should average
        Matrix<Real> m(2, 2, {
            REAL(1.0), REAL(3.0),
            REAL(5.0), REAL(4.0)
        });

        auto s = MatrixSym<Real>::FromFullMatrix(m);

        REQUIRE(2 == s.Dim());
        REQUIRE(REAL(1.0) == s(0, 0));
        REQUIRE(REAL(4.0) == s(0, 1));  // (3+5)/2 = 4
        REQUIRE(REAL(4.0) == s(1, 0));
        REQUIRE(REAL(4.0) == s(1, 1));
    }

    TEST_CASE("MatrixSym_AreEqual_static", "[MatrixSym][static-factory]")
    {
        MatrixSym<Real> a(2, { REAL(1.0), REAL(2.0), REAL(3.0) });
        MatrixSym<Real> b(2, { REAL(1.0), REAL(2.0), REAL(3.0) });

        REQUIRE(MatrixSym<Real>::AreEqual(a, b));
    }

    //=============================================================================
    // SECTION 7: RESIZE AND MODIFICATION (6 tests)
    //=============================================================================

    TEST_CASE("MatrixSym_Resize_larger", "[MatrixSym][modification]")
    {
        MatrixSym<Real> a(2, { REAL(1.0), REAL(2.0), REAL(3.0) });

        a.Resize(4);

        REQUIRE(4 == a.Dim());
        REQUIRE(10 == a.size());  // 4*5/2 = 10
        
        // All elements should be zero after resize
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                REQUIRE(REAL(0.0) == a(i, j));
            }
        }
    }

    TEST_CASE("MatrixSym_Resize_smaller", "[MatrixSym][modification]")
    {
        MatrixSym<Real> a(4);

        a.Resize(2);

        REQUIRE(2 == a.Dim());
        REQUIRE(3 == a.size());
    }

    TEST_CASE("MatrixSym_Resize_preserve_larger", "[MatrixSym][modification]")
    {
        MatrixSym<Real> a(2, { REAL(1.0), REAL(2.0), REAL(3.0) });

        a.Resize(4, true);

        REQUIRE(4 == a.Dim());
        
        // Original data preserved
        REQUIRE(REAL(1.0) == a(0, 0));
        REQUIRE(REAL(2.0) == a(0, 1));
        REQUIRE(REAL(3.0) == a(1, 1));
        
        // New elements are zero
        REQUIRE(REAL(0.0) == a(2, 2));
        REQUIRE(REAL(0.0) == a(3, 3));
    }

    TEST_CASE("MatrixSym_Resize_preserve_smaller", "[MatrixSym][modification]")
    {
        MatrixSym<Real> a(3, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0), REAL(6.0) });

        a.Resize(2, true);

        REQUIRE(2 == a.Dim());
        REQUIRE(REAL(1.0) == a(0, 0));
        REQUIRE(REAL(2.0) == a(0, 1));
        REQUIRE(REAL(3.0) == a(1, 1));
    }

    TEST_CASE("MatrixSym_SetToZero", "[MatrixSym][modification]")
    {
        MatrixSym<Real> a(2, { REAL(1.0), REAL(2.0), REAL(3.0) });

        a.SetToZero();

        REQUIRE(2 == a.Dim());
        REQUIRE(REAL(0.0) == a(0, 0));
        REQUIRE(REAL(0.0) == a(0, 1));
        REQUIRE(REAL(0.0) == a(1, 1));
    }

    TEST_CASE("MatrixSym_SetToValue", "[MatrixSym][modification]")
    {
        MatrixSym<Real> a(2);

        a.SetToValue(REAL(5.0));

        REQUIRE(REAL(5.0) == a(0, 0));
        REQUIRE(REAL(5.0) == a(0, 1));
        REQUIRE(REAL(5.0) == a(1, 1));
    }

    //=============================================================================
    // SECTION 8: EDGE CASES (6 tests)
    //=============================================================================

    TEST_CASE("MatrixSym_1x1_operations", "[MatrixSym][edge-cases]")
    {
        MatrixSym<Real> a(1, { REAL(5.0) });
        MatrixSym<Real> b(1, { REAL(3.0) });

        auto sum = a + b;
        REQUIRE(REAL(8.0) == sum(0, 0));

        auto prod = a * b;
        REQUIRE(REAL(15.0) == prod(0, 0));

        REQUIRE(REAL(5.0) == a.Trace());
        REQUIRE(REAL(5.0) == a.NormFrobenius());
    }

    TEST_CASE("MatrixSym_large_matrix", "[MatrixSym][edge-cases]")
    {
        const int dim = 50;
        MatrixSym<Real> a(dim, REAL(1.0));

        REQUIRE(dim == a.Dim());
        REQUIRE(dim * (dim + 1) / 2 == a.size());

        // Trace should be dim (all diagonal elements = 1)
        REQUIRE(static_cast<Real>(dim) == a.Trace());
    }

    TEST_CASE("MatrixSym_numerical_precision", "[MatrixSym][edge-cases]")
    {
        MatrixSym<Real> a(2, { 1e-15, 2e-15, 3e-15 });
        MatrixSym<Real> b(2, { 1e-15, 2e-15, 3e-15 });

        REQUIRE(a.IsEqual(b, 1e-16));
    }

    TEST_CASE("MatrixSym_storage_efficiency", "[MatrixSym][edge-cases]")
    {
        // Verify storage uses n*(n+1)/2 elements, not n*n
        const int dim = 10;
        MatrixSym<Real> a(dim);

        size_t symmetric_storage = dim * (dim + 1) / 2;  // 55
        size_t full_storage = dim * dim;                  // 100

        REQUIRE(symmetric_storage == a.size());
        REQUIRE(symmetric_storage < full_storage);
    }

    TEST_CASE("MatrixSym_row_column_equivalence_for_symmetric", "[MatrixSym][edge-cases]")
    {
        MatrixSym<Real> a(3, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0), REAL(6.0) });

        // For symmetric matrix, row i should equal column i
        for (int i = 0; i < 3; ++i) {
            Vector<Real> row = a.VectorFromRow(i);
            Vector<Real> col = a.VectorFromColumn(i);
            
            for (int j = 0; j < 3; ++j) {
                REQUIRE(row[j] == col[j]);
            }
        }
    }

    TEST_CASE("MatrixSym_to_string", "[MatrixSym][io]")
    {
        MatrixSym<Real> a(2, { REAL(1.0), REAL(2.0), REAL(3.0) });

        std::string str = a.to_string(5, 1);
        
        // Just verify it produces non-empty output
        REQUIRE_FALSE(str.empty());
    }

} // namespace MML::Tests::Base::MatrixSymTests

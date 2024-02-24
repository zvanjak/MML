#include "../catch/catch.hpp"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Vector.h"
#include "base/MatrixSym.h"
#endif

using namespace MML;

// TODO 0.9 - finish details
namespace MML::Tests::MatrixOtherTests
{
TEST_CASE("MatrixSym_default_ctor_init_to_zero", "[simple]") 
{
    MatrixSym<Real> a(2);

    REQUIRE( 2 == a.RowNum());
    REQUIRE( 2 == a.ColNum());

	REQUIRE(0.0 ==  a(0,0));
	REQUIRE(0.0 ==  a(0,1));
	REQUIRE(0.0 ==  a(1,0));
	REQUIRE(0.0 ==  a(1,1));
}

TEST_CASE("MatrixSym_default_ctor_init_to_value", "[simple]") {
    MatrixSym<Real> a(2, 5.0);

    REQUIRE( 2 == a.RowNum());
    REQUIRE( 2 == a.ColNum());

	REQUIRE(5.0 ==  a(0,0));
	REQUIRE(5.0 ==  a(0,1));
	REQUIRE(5.0 ==  a(1,0));
	REQUIRE(5.0 ==  a(1,1));
}
TEST_CASE("MatrixSym_initializer_list_ctor", "[simple]") {
    MatrixSym<Real> a(2, {1.0, 
                          2.0, 3.0});

	REQUIRE(2 == a.RowNum());
	REQUIRE(2 == a.ColNum());

	REQUIRE(1.0 ==  a(0,0));
	REQUIRE(2.0 ==  a(0,1));
	REQUIRE(2.0 ==  a(1,0));
	REQUIRE(3.0 ==  a(1,1));
}
TEST_CASE("MatrixSym_initializer_list_ctor2", "[simple]") {
    MatrixSym<Real> a(3, {1.0, 
                          2.0, 3.0, 
                          4.0, 5.0, 6.0});

	REQUIRE(3 == a.RowNum());
	REQUIRE(3 == a.ColNum());

	REQUIRE(1.0 ==  a(0,0));
	
    REQUIRE(2.0 ==  a(1,0));
	REQUIRE(3.0 ==  a(1,1));
	
    REQUIRE(4.0 ==  a(2,0));
    REQUIRE(5.0 ==  a(2,1));
    REQUIRE(6.0 ==  a(2,2));
}
TEST_CASE("MatrixSym_initializer_list_ctor_extra_term_throws", "[simple]") {
    
    REQUIRE_THROWS_AS(MatrixSym<Real>(2, {1.0, 2.0, 3.0, 4.0}), MatrixDimensionError); 
}

// IsEqual
TEST_CASE("MatrixSym_IsEqual_diff_matrices", "[simple]") {
    MatrixSym<Real> a(2), b(4);

    REQUIRE(false == a.IsEqual(b));
}
TEST_CASE("MatrixSym_IsEqual", "[simple]") {
    MatrixSym<Real> a(2, {1.0, 2.0, 3.0});
    MatrixSym<Real> b(2, {1.0, 2.0, 3.0});

    REQUIRE(true == a.IsEqual(b));
}
TEST_CASE("MatrixSym_IsEqual2", "[simple]") {
    MatrixSym<Real> a(2, {1.0, 2.0, 3.0});
    MatrixSym<Real> b(2, {1.0, 2.0, 4.0});

    REQUIRE(false == a.IsEqual(b));
}
TEST_CASE("MatrixSym_IsEqual3", "[simple]") {
    MatrixSym<Real> a(2, {1.0, 2.0, 4.0});
    MatrixSym<Real> b(2, {1.0, 2.0, 4.0001});

    REQUIRE(true == a.IsEqual(b, 1e-4));
    REQUIRE(false == a.IsEqual(b, 1e-5));
}

TEST_CASE("Test_MatrixSym_Op*", "[simple]") {
    MatrixSym<Real> a(2, { 1.0, 2.0, 
                                3.0});
    MatrixSym<Real> b(2, { 1.0, 2.0, 
                                3.0});
    auto c = a * b;

	REQUIRE( 5.0 ==  c[0][0]);
	REQUIRE( 8.0 ==  c[0][1]);
	REQUIRE( 8.0 ==  c[1][0]);
	REQUIRE(13.0 ==  c[1][1]);
}

// op. sa vektorom
TEST_CASE("Test_MatrixSym_mul_Vector_right", "[simple]") {
    MatrixSym<Real> a(2, {1.0, 10.0, 5.0});
    Vector<Real> b({1.0, 2.0});

	auto c = a * b;

	REQUIRE(21.0 ==  c[0]);
	REQUIRE(20.0 ==  c[1]);
}

TEST_CASE("Test_MatrixSym_mul_Vector_left", "[simple]") {
    MatrixSym<Real> a(2, {1.0, 10.0, 5.0});
    Vector<Real> b({1.0, 2.0});

	auto d = b * a;

	REQUIRE(21.0 ==  d[0]);
	REQUIRE(20.0 ==  d[1]);
}

TEST_CASE("MatrixSym_GetInverse", "[simple]") 
{
    MatrixSym<Real> mat(2, {1.0, 2.0, 3.0} );

	auto b = mat.GetInverse();
	auto c = mat.GetAsMatrix() * b;

	REQUIRE(c.IsEqual(Matrix<Real>::GetUnitMatrix(2)));
}

} // namespace MML::Tests::MatrixOtherTests
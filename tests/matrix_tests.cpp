#include "../catch/catch.hpp"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "core/Vector.h"
#include "core/Matrix.h"
#endif

using namespace MML;

TEST_CASE("Matrix_default_ctor_init_to_zero", "[simple]") {
    Matrix<Real> a(2,2);

    REQUIRE( 2 == a.RowNum());
    REQUIRE( 2 == a.ColNum());

	REQUIRE(0.0 ==  a(0,0));
	REQUIRE(0.0 ==  a(0,1));
	REQUIRE(0.0 ==  a(1,0));
	REQUIRE(0.0 ==  a(1,1));
}

TEST_CASE("Matrix_initializer_list_ctor", "[simple]") {
    Matrix<Real> a(2, 2, {1.0, 2.0, 3.0, 4.0});

	REQUIRE(2 == a.RowNum());
	REQUIRE(2 == a.ColNum());

	REQUIRE(1.0 ==  a[0][0]);
	REQUIRE(2.0 ==  a(0,1));
	REQUIRE(3.0 ==  a(1,0));
	REQUIRE(4.0 ==  a(1,1));
}

// TODO - matrix_tests - test access variants, const, pointers, references

TEST_CASE("Test_Matrix_Op+-", "[simple]") {
    Matrix<Real> a(2, 2, {1.0, 2.0, 3.0, 4.0});
    Matrix<Real> b(2, 2, {1.0, 2.0, 3.0, 4.0});

    auto c = a + b;
    auto d = a - b;

	REQUIRE(2.0 ==  c[0][0]);
	REQUIRE(4.0 ==  c[0][1]);

	REQUIRE(0.0 ==  d[0][0]);
	REQUIRE(0.0 ==  d[0][1]);
}

TEST_CASE("Test_Matrix_Op*", "[simple]") {
    Matrix<Real> a(2, 2, { 1.0, 2.0, 
                                3.0, 4.0});
    Matrix<Real> b(2, 2, { 1.0, 2.0, 
                                3.0, 4.0});
    auto c = a * b;

	REQUIRE( 7.0 ==  c[0][0]);
	REQUIRE(10.0 ==  c[0][1]);
	REQUIRE(15.0 ==  c[1][0]);
	REQUIRE(22.0 ==  c[1][1]);
}

// op. sa skalaraom
TEST_CASE("Test_Matrix_mul_double", "[simple]") {
    Matrix<Real> a(2, 2, {1.0, 100.0, 50.0, 100.0});

	auto b = a * 2.0;
	auto c = 2.0 * a;

	REQUIRE(2.0 ==  b[0][0]);
	REQUIRE(2.0 ==  c[0][0]);

	REQUIRE(200.0 ==  b[0][1]);
	REQUIRE(200.0 ==  c[0][1]);	
}

TEST_CASE("Test_Matrix_div_double", "[simple]") {
    Matrix<Real> a(2, 2, {4.0, 400.0, 1.0, 1});

	auto b = a / 2.0;

	REQUIRE(2.0 ==  b[0][0]);
	REQUIRE(200.0 ==  b[0][1]);
}

// op. sa vektorom
TEST_CASE("Test_Matrix_mul_Vector", "[simple]") {
    Matrix<Real> a(2, 2, {1.0, 10.0, 5.0, 2.0});
    Vector b({1.0, 2.0});

	auto c = a * b;
	auto d = b * a;

	REQUIRE(21.0 ==  c[0]);
	REQUIRE( 9.0 ==  c[1]);

	REQUIRE(11.0 ==  d[0]);
	REQUIRE(14.0 ==  d[1]);
}

TEST_CASE("Matrix_Transpose", "[simple]") 
{
    Matrix<Real> mat(2,2, {1.0, 2.0, 3.0, 4.0} );
    Matrix<Real> matTransp(2,2, {1.0, 3.0, 2.0, 4.0} );

	mat.Transpose();

	REQUIRE(mat.IsEqual(matTransp));
}

TEST_CASE("Matrix_GetTranspose", "[simple]") 
{
    Matrix<Real> mat(2,2, {1.0, 2.0, 3.0, 4.0} );
    Matrix<Real> matTransp(2,2, {1.0, 3.0, 2.0, 4.0} );

	auto trans = mat.GetTranspose();

	REQUIRE(trans.IsEqual(matTransp));
}

TEST_CASE("Matrix_GetInverse", "[simple]") 
{
    Matrix<Real> mat(2,2, {1.0, 2.0, 3.0, 4.0} );

	auto b = mat.GetInverse();

	auto c = mat * b;

	REQUIRE(c.IsEqual(Matrix<Real>::GetUnitMatrix(2)));
}

TEST_CASE("Test_Matrix_exceptions", "[simple]") 
{
    Matrix<Real> mat_1(2,3);
    REQUIRE_THROWS_AS(mat_1.MakeUnitMatrix(), MatrixDimensionError); 

    Matrix<Real> mat_21(3,3);
    Matrix<Real> mat_22(2,2);
    REQUIRE_THROWS_AS(mat_21 + mat_22, MatrixDimensionError); 
    REQUIRE_THROWS_AS(mat_21 - mat_22, MatrixDimensionError); 
    REQUIRE_THROWS_AS(mat_21 * mat_22, MatrixDimensionError); 
 
    Matrix<Real> mat(2,2, {1.0, 0.0, 0.0, 1.0} );
    Vector<Real> a{1.0, 0.0, 0.0};

    REQUIRE_THROWS_AS(a * mat, MatrixDimensionError); 
    REQUIRE_THROWS_AS(mat * a, MatrixDimensionError); 

    // invert i transpose
    Matrix<Real> mat_3(3,4);
    REQUIRE_THROWS_AS(mat_3.Invert(), MatrixDimensionError); 
    REQUIRE_THROWS_AS(mat_3.GetInverse(), MatrixDimensionError); 
    REQUIRE_THROWS_AS(mat_3.Transpose(), MatrixDimensionError); 
}
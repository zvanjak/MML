#include "../catch/catch.hpp"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/MatrixNM.h"
#endif

using namespace MML;

TEST_CASE("MatrixNM_default_ctor_init_to_zero", "[simple]") {
    MatrixNM<Real,2,2> a;

	REQUIRE(0.0 ==  a(0,0));
	REQUIRE(0.0 ==  a(0,1));
	REQUIRE(0.0 ==  a(1,0));
	REQUIRE(0.0 ==  a(1,1));
}

TEST_CASE("MatrixNM_initializer_list_ctor", "[simple]") {
    MatrixNM<Real,2,2> a({1.0, 2.0, 3.0, 4.0});

	REQUIRE(2 == a.RowNum());
	REQUIRE(2 == a.ColNum());

	REQUIRE(1.0 ==  a(0,0));
	REQUIRE(2.0 ==  a(0,1));
	REQUIRE(3.0 ==  a(1,0));
	REQUIRE(4.0 ==  a(1,1));
}

TEST_CASE("MatrixNM_MakeUnitMatrix", "[simple]") {
    MatrixNM<Real, 2, 2> a;

    a.MakeUnitMatrix();

    REQUIRE(1.0 ==  a(0,0));
    REQUIRE(0.0 ==  a(0,1));
    REQUIRE(0.0 ==  a(1,0));
    REQUIRE(1.0 ==  a(1,1));
}

TEST_CASE("MatrixNM_GetUnitMatrix", "[simple]") {
    auto a = MatrixNM<Real, 2, 2>::GetUnitMatrix();

    REQUIRE(1.0 ==  a[0][0]);
    REQUIRE(0.0 ==  a[0][1]);
    REQUIRE(0.0 ==  a[1][0]);
    REQUIRE(1.0 ==  a[1][1]);
}

TEST_CASE("MatrixNM_IsEqual", "[simple]") {
    MatrixNM<Real, 2, 2> a({1.0, 2.0, 3.0, 4.0});
    MatrixNM<Real, 2, 2> b({1.0, 2.0, 3.0, 4.0});

    REQUIRE(true == a.IsEqual(b));
}
TEST_CASE("MatrixNM_IsEqual2", "[simple]") {
    MatrixNM<Real, 2, 2> a({1.0, 2.0, 3.0, 4.0});
    MatrixNM<Real, 2, 2> b({1.0, 2.0, 3.0, 5.0});

    REQUIRE(false == a.IsEqual(b));
}
TEST_CASE("MatrixNM_IsEqual3", "[simple]") {
    MatrixNM<Real, 2, 2> a({1.0, 2.0, 3.0, 4.0});
    MatrixNM<Real, 2, 2> b({1.0, 2.0, 3.0, 4.0001});

    REQUIRE(true == a.IsEqual(b, 1e-4));
    REQUIRE(false == a.IsEqual(b, 1e-5));
}

TEST_CASE("Test_MatrixNM_Op+-", "[simple]") {
    MatrixNM<Real, 2, 2> a({1.0, 2.0, 3.0, 4.0});
    MatrixNM<Real, 2, 2> b({1.0, 2.0, 3.0, 4.0});

    auto c = a + b;
    auto d = a - b;

	REQUIRE(2.0 ==  c(0,0));
	REQUIRE(4.0 ==  c(0,1));

	REQUIRE(0.0 ==  d(0,0));
	REQUIRE(0.0 ==  d(0,1));
}

TEST_CASE("Test_MatrixNM_Op*", "[simple]") {
    MatrixNM<Real, 2, 2> a({1.0, 2.0, 3.0, 4.0});
    MatrixNM<Real, 2, 2> b({1.0, 2.0, 3.0, 4.0});

    auto c = a * b;

	REQUIRE( 7.0 ==  c(0,0));
	REQUIRE(10.0 ==  c(0,1));
	REQUIRE(15.0 ==  c(1,0));
	REQUIRE(22.0 ==  c(1,1));
}

TEST_CASE("Test_MatrixNM_mul_double", "[simple]") {
    MatrixNM<Real, 2, 2> a({1.0, 100.0, 50.0, 100.0});

	auto b = a * 2.0;
	auto c = 2.0 * a;

	REQUIRE(2.0 ==  b(0,0));
	REQUIRE(2.0 ==  c(0,0));

	REQUIRE(200.0 ==  b(0,1));
	REQUIRE(200.0 ==  c(0,1));	
}

TEST_CASE("Test_MatrixNM_div_double", "[simple]") {
    MatrixNM<Real, 2, 2> a({4.0, 400.0, 1.0, 1});

	auto b = a / 2.0;

	REQUIRE(2.0 ==  b(0,0));
	REQUIRE(200.0 ==  b(0,1));
}

TEST_CASE("Test_MatrixNM_mul_Vector", "[simple]") {
    MatrixNM<Real, 2, 2> a({1.0, 10.0, 
                                 5.0, 2.0});
    VectorN<Real, 2> b({1.0, 2.0});

	auto c = a * b;
	auto d = b * a;

	REQUIRE(21.0 ==  c[0]);
	REQUIRE( 9.0 ==  c[1]);

	REQUIRE(11.0 ==  d[0]);
	REQUIRE(14.0 ==  d[1]);
}

TEST_CASE("MatrixNM_Transpose", "[simple]") 
{
    MatrixNM<Real,2,2> mat({1.0, 2.0, 3.0, 4.0} );
    MatrixNM<Real,2,2> matTransp({1.0, 3.0, 2.0, 4.0} );

	mat.Transpose();

	REQUIRE(mat.IsEqual(matTransp));
}

TEST_CASE("MatrixNM_GetTranspose", "[simple]") 
{
    MatrixNM<Real,2,2> mat({1.0, 2.0, 3.0, 4.0} );
    MatrixNM<Real,2,2> matTransp({1.0, 3.0, 2.0, 4.0} );

	auto trans = mat.GetTranspose();

	REQUIRE(trans.IsEqual(matTransp));
}

TEST_CASE("MatrixNM_GetTranspose_nonrectangular", "[simple]") 
{
    MatrixNM<Real,2,3> mat({1.0, 2.0, 3.0, 4.0, 5.0, 6.0} );
    MatrixNM<Real,3,2> matTransp({1.0, 4.0, 2.0, 5.0, 3.0, 6.0} );

	auto trans = mat.GetTranspose();

	REQUIRE(trans.IsEqual(matTransp));
}

TEST_CASE("MatrixNM_GetInverse", "[simple]") 
{
    MatrixNM<Real,2,2> mat({1.0, 2.0, 3.0, 4.0} );

	auto b = mat.GetInverse();

	auto c = mat * b;

	MatrixNM<Real,2,2> d;
	d.MakeUnitMatrix();

	REQUIRE(c.IsEqual(d));
}

TEST_CASE("Test_MatrixNM_exceptions", "[simple]") 
{
    MatrixNM<Real, 3, 4> mat_3;

    REQUIRE_THROWS_AS(mat_3.Invert(), MatrixDimensionError); 
    REQUIRE_THROWS_AS(mat_3.GetInverse(), MatrixDimensionError); 
    REQUIRE_THROWS_AS(mat_3.Transpose(), MatrixDimensionError); 
}
#include "../catch/catch.hpp"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "basic_types/Vector.h"
#include "basic_types/Matrix.h"
#endif

TEST_CASE("Matrix_default_ctor_init_to_zero", "[simple]") {
    MML::Matrix<Real> a(2,2);

	REQUIRE(0.0 ==  a(0,0));
	REQUIRE(0.0 ==  a(0,1));
	REQUIRE(0.0 ==  a(1,0));
	REQUIRE(0.0 ==  a(1,1));
}

TEST_CASE("Matrix_initializer_list_ctor", "[simple]") {
    MML::Matrix<Real> a(2, 2, {1.0, 2.0, 3.0, 4.0});

	REQUIRE(2 == a.RowNum());
	REQUIRE(2 == a.ColNum());

	REQUIRE(1.0 ==  a(0,0));
	REQUIRE(2.0 ==  a(0,1));
	REQUIRE(3.0 ==  a(1,0));
	REQUIRE(4.0 ==  a(1,1));
}

TEST_CASE("Matrix_Transpose", "[simple]") 
{
    MML::Matrix<Real> mat(2,2, {1.0, 2.0, 3.0, 4.0} );
    MML::Matrix<Real> matTransp(2,2, {1.0, 3.0, 2.0, 4.0} );

	mat.Transpose();

	REQUIRE(mat.IsEqual(matTransp));
}

TEST_CASE("Matrix_GetTranspose", "[simple]") 
{
    MML::Matrix<Real> mat(2,2, {1.0, 2.0, 3.0, 4.0} );
    MML::Matrix<Real> matTransp(2,2, {1.0, 3.0, 2.0, 4.0} );

	auto trans = mat.GetTranspose();

	REQUIRE(trans.IsEqual(matTransp));
}

TEST_CASE("Matrix_GetInverse", "[simple]") 
{
    MML::Matrix<Real> mat(2,2, {1.0, 2.0, 3.0, 4.0} );

	auto b = mat.GetInverse();

	auto c = mat * b;

	MML::Matrix<Real> d(2,2);
	d.MakeUnitMatrix();

	REQUIRE(c.IsEqual(d));
}

TEST_CASE("Test_Matrix", "[simple]") 
{
    MML::Matrix<Real> mat(2,2, {1.0, 0.0, 0.0, 1.0} );
    MML::Vector<Real> a{1.0, 0.0, 0.0};

    REQUIRE_THROWS_AS(a * mat, MML::MatrixDimensionError); 
    
}
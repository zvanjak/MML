#include "../catch/catch.hpp"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "basic_types/MatrixNM.h"
#endif

TEST_CASE("MatrixNM_default_ctor_init_to_zero", "[simple]") {
    MML::MatrixNM<Real,2,2> a;

	REQUIRE(0.0 ==  a(0,0));
	REQUIRE(0.0 ==  a(0,1));
	REQUIRE(0.0 ==  a(1,0));
	REQUIRE(0.0 ==  a(1,1));
}

TEST_CASE("MatrixNM_initializer_list_ctor", "[simple]") {
    MML::MatrixNM<Real,2,2> a({1.0, 2.0, 3.0, 4.0});

	REQUIRE(2 == a.RowNum());
	REQUIRE(2 == a.ColNum());

	REQUIRE(1.0 ==  a(0,0));
	REQUIRE(2.0 ==  a(0,1));
	REQUIRE(3.0 ==  a(1,0));
	REQUIRE(4.0 ==  a(1,1));
}


TEST_CASE("MatrixNM_Transpose", "[simple]") 
{
    MML::MatrixNM<Real,2,2> mat({1.0, 2.0, 3.0, 4.0} );
    MML::MatrixNM<Real,2,2> matTransp({1.0, 3.0, 2.0, 4.0} );

	mat.Transpose();

	REQUIRE(mat.IsEqual(matTransp));
}

TEST_CASE("MatrixNM_GetTranspose", "[simple]") 
{
    MML::MatrixNM<Real,2,2> mat({1.0, 2.0, 3.0, 4.0} );
    MML::MatrixNM<Real,2,2> matTransp({1.0, 3.0, 2.0, 4.0} );

	auto trans = mat.GetTranspose();

	REQUIRE(trans.IsEqual(matTransp));
}

TEST_CASE("MatrixNM_GetInverse", "[simple]") 
{
    MML::MatrixNM<Real,2,2> mat({1.0, 2.0, 3.0, 4.0} );

	auto b = mat.GetInverse();

	auto c = mat * b;

	MML::MatrixNM<Real,2,2> d;
	d.MakeUnitMatrix();

	REQUIRE(c.IsEqual(d));
}

TEST_CASE("Test_MatrixNM", "[simple]") {
    MML::MatrixNM<Real,2,2> mat({1.0, 0.0, 0.0, 1.0} );

}
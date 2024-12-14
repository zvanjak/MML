#include "../catch/catch.hpp"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Vector.h"
#include "base/Matrix.h"
#include "base/BaseUtils.h"
#endif

using namespace MML;

namespace MML::Tests::Base::MatrixTests
{
	///////////////////////          Constructors and destructor       //////////////////////
	TEST_CASE("Matrix_default_ctor_init_to_zero", "[simple]") {
		Matrix<Real> a(2, 2);

		REQUIRE(2 == a.RowNum());
		REQUIRE(2 == a.ColNum());

		REQUIRE(0.0 == a(0, 0));
		REQUIRE(0.0 == a(0, 1));
		REQUIRE(0.0 == a(1, 0));
		REQUIRE(0.0 == a(1, 1));

		Matrix<Real> b(1000, 1000);

		REQUIRE(1000 == b.RowNum());
		REQUIRE(1000 == b.ColNum());

		REQUIRE(0.0 == b(0, 0));
		REQUIRE(0.0 == b(158, 738));
		REQUIRE(0.0 == b(35, 999));
		REQUIRE(0.0 == b(0, 100));
	}
	TEST_CASE("Matrix_default_ctor_init_to_value", "[simple]") {
		Matrix<Real> a(2, 2, 5.0);

		REQUIRE(2 == a.RowNum());
		REQUIRE(2 == a.ColNum());

		REQUIRE(5.0 == a(0, 0));
		REQUIRE(5.0 == a(0, 1));
		REQUIRE(5.0 == a(1, 0));
		REQUIRE(5.0 == a(1, 1));
	}
	TEST_CASE("Matrix_initializer_list_ctor", "[simple]") {
		Matrix<Real> a(2, 2, { 1.0, 2.0, 3.0, 4.0 });

		REQUIRE(2 == a.RowNum());
		REQUIRE(2 == a.ColNum());

		REQUIRE(1.0 == a[0][0]);
		REQUIRE(2.0 == a(0, 1));
		REQUIRE(3.0 == a(1, 0));
		REQUIRE(4.0 == a(1, 1));
	}
	TEST_CASE("Matrix_double_ptr_ctor", "[simple]") {
		Real vec[4] = { 1.0, 2.0, 3.0, 4.0 };
		Matrix<Real> a(2, 2, vec);

		REQUIRE(2 == a.RowNum());
		REQUIRE(2 == a.ColNum());

		REQUIRE(1.0 == a[0][0]);
		REQUIRE(2.0 == a(0, 1));
		REQUIRE(3.0 == a(1, 0));
		REQUIRE(4.0 == a(1, 1));
	}

	///////////////////////              Standard stuff                //////////////////////
	TEST_CASE("Matrix_Resize", "[simple]") {
		Matrix<Real> a(2, 2);

		a.Resize(3, 5);

		REQUIRE(3 == a.RowNum());
		REQUIRE(5 == a.ColNum());
	}
	TEST_CASE("Matrix_Resize_throws_for_negative_size", "[simple]") {
		Matrix<Real> a(2, 2);

		REQUIRE_THROWS_AS(a.Resize(3, -1), MatrixDimensionError);
	}
	TEST_CASE("Matrix_MakeUnitMatrix", "[simple]") {
		Matrix<Real> a(2, 2);

		a.MakeUnitMatrix();

		REQUIRE(1.0 == a[0][0]);
		REQUIRE(0.0 == a[0][1]);
		REQUIRE(0.0 == a[1][0]);
		REQUIRE(1.0 == a[1][1]);
	}
	TEST_CASE("Matrix_MakeUnitMatrix_throws_if_not_square_matrix", "[simple]") {
		Matrix<Real> mat_1(2, 3);

		REQUIRE_THROWS_AS(mat_1.MakeUnitMatrix(), MatrixDimensionError);
	}
	TEST_CASE("Matrix_GetUnitMatrix", "[simple]") {
		auto a = Matrix<Real>::GetUnitMatrix(2);

		REQUIRE(1.0 == a[0][0]);
		REQUIRE(0.0 == a[0][1]);
		REQUIRE(0.0 == a[1][0]);
		REQUIRE(1.0 == a[1][1]);
	}

	// TODO - test GetLower i GetUpper
	
	///////////////////////          Matrix to Vector conversions      //////////////////////
	TEST_CASE("Matrix_VectorFromRow", "[simple]") {
		Matrix<Real> a(1, 3, { 1.0, 2.0, 3.0 });

		auto b = a.VectorFromRow(0);

		REQUIRE(3 == b.size());

		REQUIRE(1.0 == b[0]);
		REQUIRE(2.0 == b[1]);
		REQUIRE(3.0 == b[2]);
	}
	TEST_CASE("Matrix_VectorFromRow_throws_for_wrong_index", "[simple]") {
		Matrix<Real> a(1, 3, { 1.0, 2.0, 3.0 });

		Vector<Real> b;
		REQUIRE_THROWS_AS(b = a.VectorFromRow(-1), MatrixAccessBoundsError);
		REQUIRE_THROWS_AS(b = a.VectorFromRow(3), MatrixAccessBoundsError);
	}
	TEST_CASE("Matrix_VectorFromColumn", "[simple]") {
		Matrix<Real> a(3, 1, { 1.0, 2.0, 3.0 });

		auto b = a.VectorFromColumn(0);

		REQUIRE(3 == b.size());

		REQUIRE(1.0 == b[0]);
		REQUIRE(2.0 == b[1]);
		REQUIRE(3.0 == b[2]);
	}
	TEST_CASE("Matrix_VectorFromColumn_throws_for_wrong_index", "[simple]") {
		Matrix<Real> a(3, 1, { 1.0, 2.0, 3.0 });

		Vector<Real> b;
		REQUIRE_THROWS_AS(b = a.VectorFromColumn(-1), MatrixAccessBoundsError);
		REQUIRE_THROWS_AS(b = a.VectorFromColumn(3), MatrixAccessBoundsError);
	}
	TEST_CASE("Matrix_VectorFromDiagonal", "[simple]") {
		Matrix<Real> a(3, 3, { 1.0, 2.0, 3.0,
													3.0, 2.0, 5.0,
													6.0, 1.0, 3.0 });

		auto b = a.VectorFromDiagonal();

		REQUIRE(3 == b.size());

		REQUIRE(1.0 == b[0]);
		REQUIRE(2.0 == b[1]);
		REQUIRE(3.0 == b[2]);
	}
	TEST_CASE("Matrix_VectorFromDiagonal_throws_for_non_square_matrix", "[simple]") {
		Matrix<Real> a(3, 2, { 1.0, 2.0, 3.0,
													3.0, 2.0, 5.0 });

		Vector<Real> b;
		REQUIRE_THROWS_AS(b = a.VectorFromDiagonal(), MatrixDimensionError);
	}

	///////////////////////               Matrix properties            //////////////////////
	TEST_CASE("Matrix_IsUnit", "[simple]") {
		Matrix<Real> a(2, 2, {1.0, 0.0, 0.0, 1.0});

		REQUIRE(true == a.IsUnit());

		Matrix<Real> b(2, 2, { 1.0, 0.0, 0.0, 2.0 });

		REQUIRE(false == b.IsUnit());
	}
	TEST_CASE("Matrix_IsDiagonal", "[simple]") {
		Matrix<Real> a(2, 2, { 1.0, 0.0, 0.0, 5.0 });

		REQUIRE(true == a.IsDiagonal());

		Matrix<Real> b(2, 2, { 1.0, 2.0, 0.0, 2.0 });

		REQUIRE(false == b.IsDiagonal());
	}

	///////////////////////             Assignment operators           //////////////////////

	///////////////////////               Access operators             //////////////////////
	
	///////////////////////              Equality operations           //////////////////////
	TEST_CASE("Matrix_IsEqual_diff_size_matrices", "[simple]") {
		Matrix<Real> a(2, 3), b(4, 5);

		REQUIRE(false == a.IsEqualTo(b));
	}
	TEST_CASE("Matrix_IsEqual", "[simple]") {
		Matrix<Real> a(2, 2, { 1.0, 2.0, 3.0, 4.0 });
		Matrix<Real> b(2, 2, { 1.0, 2.0, 3.0, 4.0 });

		REQUIRE(true == a.IsEqualTo(b));
	}
	TEST_CASE("Matrix_IsEqual2", "[simple]") {
		Matrix<Real> a(2, 2, { 1.0, 2.0, 3.0, 4.0 });
		Matrix<Real> b(2, 2, { 1.0, 2.0, 3.0, 5.0 });

		REQUIRE(false == a.IsEqualTo(b));
	}
	TEST_CASE("Matrix_IsEqual3", "[simple]") {
		Matrix<Real> a(2, 2, { 1.0, 2.0, 3.0, 4.0 });
		Matrix<Real> b(2, 2, { 1.0, 2.0, 3.0, 4.0001 });

		REQUIRE(true == a.IsEqualTo(b, 1e-4));
		REQUIRE(false == a.IsEqualTo(b, 1e-5));
	}

	///////////////////////              Arithmetic operators          //////////////////////
	TEST_CASE("Matrix_Op+-", "[simple]") {
		Matrix<Real> a(2, 2, { 1.0, 2.0, 3.0, 4.0 });
		Matrix<Real> b(2, 2, { 1.0, 2.0, 3.0, 4.0 });

		auto c = a + b;
		auto d = a - b;

		REQUIRE(2.0 == c[0][0]);
		REQUIRE(4.0 == c[0][1]);

		REQUIRE(0.0 == d[0][0]);
		REQUIRE(0.0 == d[0][1]);
	}
	TEST_CASE("Matrix_Op*", "[simple]") {
		Matrix<Real> a(2, 2, { 1.0, 2.0,
																3.0, 4.0 });
		Matrix<Real> b(2, 2, { 1.0, 2.0,
																3.0, 4.0 });
		auto c = a * b;

		REQUIRE(7.0 == c[0][0]);
		REQUIRE(10.0 == c[0][1]);
		REQUIRE(15.0 == c[1][0]);
		REQUIRE(22.0 == c[1][1]);
	}
	TEST_CASE("Matrix_mul_double", "[simple]") {
		Matrix<Real> a(2, 2, { 1.0, 100.0, 50.0, 100.0 });

		auto b = a * 2.0;
		auto c = 2.0 * a;

		REQUIRE(2.0 == b[0][0]);
		REQUIRE(2.0 == c[0][0]);

		REQUIRE(200.0 == b[0][1]);
		REQUIRE(200.0 == c[0][1]);
	}
	TEST_CASE("Matrix_div_double", "[simple]") {
		Matrix<Real> a(2, 2, { 4.0, 400.0, 1.0, 1 });

		auto b = a / 2.0;

		REQUIRE(2.0 == b[0][0]);
		REQUIRE(200.0 == b[0][1]);
	}
	TEST_CASE("Matrix_mul_Vector_right", "[simple]") {
		Matrix<Real> a(2, 2, { 1.0, 10.0, 5.0, 2.0 });
		Vector<Real> b({ 1.0, 2.0 });

		auto c = a * b;

		REQUIRE(21.0 == c[0]);
		REQUIRE(9.0 == c[1]);
	}
	TEST_CASE("Matrix_mul_Vector_left", "[simple]") {
		Matrix<Real> a(2, 2, { 1.0, 10.0, 5.0, 2.0 });
		Vector<Real> b({ 1.0, 2.0 });

		auto d = b * a;

		REQUIRE(11.0 == d[0]);
		REQUIRE(14.0 == d[1]);
	}

	///////////////////////            Trace, Inverse & Transpose      //////////////////////
	TEST_CASE("Matrix_Transpose", "[simple]")
	{
		Matrix<Real> mat(2, 2, { 1.0, 2.0, 3.0, 4.0 });
		Matrix<Real> matTransp(2, 2, { 1.0, 3.0, 2.0, 4.0 });

		mat.Transpose();

		REQUIRE(mat.IsEqualTo(matTransp));
	}
	TEST_CASE("Matrix_GetTranspose", "[simple]")
	{
		Matrix<Real> mat(2, 2, { 1.0, 2.0, 3.0, 4.0 });
		Matrix<Real> matTransp(2, 2, { 1.0, 3.0, 2.0, 4.0 });

		auto trans = mat.GetTranspose();

		REQUIRE(trans.IsEqualTo(matTransp));

		Matrix<Real> mat1(2, 3, { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 });
		Matrix<Real> matTransp1(3, 2, { 1.0, 4.0, 2.0, 5.0, 3.0, 6.0 });

		trans = mat1.GetTranspose();

		REQUIRE(trans.IsEqualTo(matTransp1));

		// test transposing row matrix
		Matrix<Real> mat2(1, 4, { 1.0, 2.0, 3.0, 4.0 });
		Matrix<Real> matTransp2(4, 1, { 1.0, 2.0, 3.0, 4.0 });

		auto mat22 = mat2.GetTranspose();

		REQUIRE(mat22.IsEqualTo(matTransp2));

		// test transposing column matrix
		Matrix<Real> mat3(4, 1, { 1.0, 2.0, 3.0, 4.0 });
		Matrix<Real> matTransp3(1, 4, { 1.0, 2.0, 3.0, 4.0 });

		auto mat33 = mat3.GetTranspose();

		REQUIRE(mat33.IsEqualTo(matTransp3));
	}
	TEST_CASE("Matrix_GetInverse", "[simple]")
	{
		Matrix<Real> mat(2, 2, { 1.0, 2.0, 3.0, 4.0 });

		auto b = mat.GetInverse();
		auto c = mat * b;

		REQUIRE(c.IsEqualTo(Matrix<Real>::GetUnitMatrix(2)));
	}
	TEST_CASE("Matrix_GetInverse_almost_signular", "[simple]")
	{
		Matrix<Real> mat(3, 3, { 1.0, 2.0, 1.0,
														 1.0, 2.0, 1.00001,
														-4.0, 6.0, 4.0});

		auto b = mat.GetInverse();
		auto c = mat * b;

		REQUIRE(c.IsEqualTo(Matrix<Real>::GetUnitMatrix(3), 1e-10));
		REQUIRE(!c.IsEqualTo(Matrix<Real>::GetUnitMatrix(3), 1e-11));
	}
	TEST_CASE("Matrix_GetInverse_throws_for_singular_matrix", "[simple]")
	{
		// test case for singular matrices with Invert
		Matrix<Real> mat2(2, 2, { 1.0, 2.0, 2.0, 4.0 });

		REQUIRE_THROWS_AS(mat2.GetInverse(), SingularMatrixError);
	}
	
	///////////////////////////               I/O                 ///////////////////////////
	TEST_CASE("Matrix_to_string", "[simple]")
	{
		Matrix<Real> mat(2, 2, { 1.0, 2.0, 3.0, 4.0 });

		REQUIRE("Rows: 2 Cols: 2\n[     1,     2 ]\n[     3,     4 ]" == mat.to_string(5, 3));
	}
	TEST_CASE("Matrix_Print_with_zerotreshold", "[simple]")
	{
		std::stringstream str;

		Matrix<Real> mat(2, 2, { 1.0, 2.0, 3.0, 1e-10 });
		mat.Print(str, 5, 3, 1e-9);

		REQUIRE("Rows: 2 Cols: 2\n[     1,     2 ]\n[     3,     0 ]" == str.str());
	}
	TEST_CASE("Test_Matrix_exceptions", "[simple]")
	{
		REQUIRE_THROWS_AS(Matrix<Real>(-1, 5), MatrixDimensionError);


		Matrix<Real> mat_21(3, 3);
		Matrix<Real> mat_22(2, 2);
		REQUIRE_THROWS_AS(mat_21 + mat_22, MatrixDimensionError);
		REQUIRE_THROWS_AS(mat_21 - mat_22, MatrixDimensionError);
		REQUIRE_THROWS_AS(mat_21 * mat_22, MatrixDimensionError);

		Matrix<Real> mat(2, 2, { 1.0, 0.0, 0.0, 1.0 });
		Vector<Real> a{ 1.0, 0.0, 0.0 };

		REQUIRE_THROWS_AS(a * mat, MatrixDimensionError);
		REQUIRE_THROWS_AS(mat * a, MatrixDimensionError);

		// invert i transpose
		Matrix<Real> mat_3(3, 4);
		REQUIRE_THROWS_AS(mat_3.Invert(), MatrixDimensionError);
		REQUIRE_THROWS_AS(mat_3.GetInverse(), MatrixDimensionError);
		REQUIRE_THROWS_AS(mat_3.Transpose(), MatrixDimensionError);
	}
} // namespace MML::Tests::Matrix
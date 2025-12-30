#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Vector.h"
#include "base/Matrix.h"
#include "base/BaseUtils.h"
#endif

using namespace MML;
using namespace MML::Testing;

namespace MML::Tests::Base::MatrixTests
{
	///////////////////////          Constructors and destructor       //////////////////////
	TEST_CASE("Matrix::default_ctor_init_to_zero", "[simple]") {
			TEST_PRECISION_INFO();
		Matrix<Real> a(2, 2);

		REQUIRE(2 == a.RowNum());
		REQUIRE(2 == a.ColNum());

		REQUIRE(REAL(0.0) == a(0, 0));
		REQUIRE(REAL(0.0) == a(0, 1));
		REQUIRE(REAL(0.0) == a(1, 0));
		REQUIRE(REAL(0.0) == a(1, 1));

		Matrix<Real> b(1000, 1000);

		REQUIRE(1000 == b.RowNum());
		REQUIRE(1000 == b.ColNum());

		REQUIRE(REAL(0.0) == b(0, 0));
		REQUIRE(REAL(0.0) == b(158, 738));
		REQUIRE(REAL(0.0) == b(35, 999));
		REQUIRE(REAL(0.0) == b(0, 100));
	}
	TEST_CASE("Matrix::default_ctor_init_to_value", "[simple]") {
			TEST_PRECISION_INFO();
		Matrix<Real> a(2, 2, REAL(5.0));

		REQUIRE(2 == a.RowNum());
		REQUIRE(2 == a.ColNum());

		REQUIRE(REAL(5.0) == a(0, 0));
		REQUIRE(REAL(5.0) == a(0, 1));
		REQUIRE(REAL(5.0) == a(1, 0));
		REQUIRE(REAL(5.0) == a(1, 1));
	}
	TEST_CASE("Matrix::initializer_list_ctor", "[simple]") {
			TEST_PRECISION_INFO();
		Matrix<Real> a(2, 2, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });

		REQUIRE(2 == a.RowNum());
		REQUIRE(2 == a.ColNum());

		REQUIRE(REAL(1.0) == a[0][0]);
		REQUIRE(REAL(2.0) == a(0, 1));
		REQUIRE(REAL(3.0) == a(1, 0));
		REQUIRE(REAL(4.0) == a(1, 1));
	}
	TEST_CASE("Matrix::double_ptr_ctor", "[simple]") {
			TEST_PRECISION_INFO();
		Real vec[4] = { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) };
		Matrix<Real> a(2, 2, vec);

		REQUIRE(2 == a.RowNum());
		REQUIRE(2 == a.ColNum());

		REQUIRE(REAL(1.0) == a[0][0]);
		REQUIRE(REAL(2.0) == a(0, 1));
		REQUIRE(REAL(3.0) == a(1, 0));
		REQUIRE(REAL(4.0) == a(1, 1));
	}

	///////////////////////              Standard stuff                //////////////////////
	TEST_CASE("Matrix::Resize", "[simple]") {
			TEST_PRECISION_INFO();
		Matrix<Real> a(2, 2);

		a.Resize(3, 5);

		REQUIRE(3 == a.RowNum());
		REQUIRE(5 == a.ColNum());
	}
	TEST_CASE("Matrix::Resize_throws_for_negative_size", "[simple]") {
			TEST_PRECISION_INFO();
		Matrix<Real> a(2, 2);

		REQUIRE_THROWS_AS(a.Resize(3, -1), MatrixDimensionError);
	}
	TEST_CASE("Matrix::MakeUnitMatrix", "[simple]") {
			TEST_PRECISION_INFO();
		Matrix<Real> a(2, 2);

		a.MakeUnitMatrix();

		REQUIRE(REAL(1.0) == a[0][0]);
		REQUIRE(REAL(0.0) == a[0][1]);
		REQUIRE(REAL(0.0) == a[1][0]);
		REQUIRE(REAL(1.0) == a[1][1]);
	}
	TEST_CASE("Matrix::MakeUnitMatrix throws_if_not_square_matrix", "[simple]") {
			TEST_PRECISION_INFO();
		Matrix<Real> mat_1(2, 3);

		REQUIRE_THROWS_AS(mat_1.MakeUnitMatrix(), MatrixDimensionError);
	}
	TEST_CASE("Matrix::GetUnitMatrix", "[simple]") {
			TEST_PRECISION_INFO();
		auto a = Matrix<Real>::GetUnitMatrix(2);

		REQUIRE(REAL(1.0) == a[0][0]);
		REQUIRE(REAL(0.0) == a[0][1]);
		REQUIRE(REAL(0.0) == a[1][0]);
		REQUIRE(REAL(1.0) == a[1][1]);
	}

	// TODO - test GetLower i GetUpper
	
	///////////////////////          Matrix to Vector conversions      //////////////////////
	TEST_CASE("Matrix::VectorFromRow", "[simple]") {
			TEST_PRECISION_INFO();
		Matrix<Real> a(1, 3, { REAL(1.0), REAL(2.0), REAL(3.0) });

		auto b = a.VectorFromRow(0);

		REQUIRE(3 == b.size());

		REQUIRE(REAL(1.0) == b[0]);
		REQUIRE(REAL(2.0) == b[1]);
		REQUIRE(REAL(3.0) == b[2]);
	}
	TEST_CASE("Matrix::VectorFromRow_throws_for_wrong_index", "[simple]") {
			TEST_PRECISION_INFO();
		Matrix<Real> a(1, 3, { REAL(1.0), REAL(2.0), REAL(3.0) });

		Vector<Real> b;
		REQUIRE_THROWS_AS(b = a.VectorFromRow(-1), MatrixAccessBoundsError);
		REQUIRE_THROWS_AS(b = a.VectorFromRow(3), MatrixAccessBoundsError);
	}
	TEST_CASE("Matrix::VectorFromColumn", "[simple]") {
			TEST_PRECISION_INFO();
		Matrix<Real> a(3, 1, { REAL(1.0), REAL(2.0), REAL(3.0) });

		auto b = a.VectorFromColumn(0);

		REQUIRE(3 == b.size());

		REQUIRE(REAL(1.0) == b[0]);
		REQUIRE(REAL(2.0) == b[1]);
		REQUIRE(REAL(3.0) == b[2]);
	}
	TEST_CASE("Matrix::VectorFromColumn_throws_for_wrong_index", "[simple]") {
			TEST_PRECISION_INFO();
		Matrix<Real> a(3, 1, { REAL(1.0), REAL(2.0), REAL(3.0) });

		Vector<Real> b;
		REQUIRE_THROWS_AS(b = a.VectorFromColumn(-1), MatrixAccessBoundsError);
		REQUIRE_THROWS_AS(b = a.VectorFromColumn(3), MatrixAccessBoundsError);
	}
	TEST_CASE("Matrix::VectorFromDiagonal", "[simple]") {
			TEST_PRECISION_INFO();
		Matrix<Real> a(3, 3, { REAL(1.0), REAL(2.0), REAL(3.0),
													REAL(3.0), REAL(2.0), REAL(5.0),
													REAL(6.0), REAL(1.0), REAL(3.0) });

		auto b = a.VectorFromDiagonal();

		REQUIRE(3 == b.size());

		REQUIRE(REAL(1.0) == b[0]);
		REQUIRE(REAL(2.0) == b[1]);
		REQUIRE(REAL(3.0) == b[2]);
	}
	TEST_CASE("Matrix::VectorFromDiagonal_throws_for_non_square_matrix", "[simple]") {
			TEST_PRECISION_INFO();
		Matrix<Real> a(3, 2, { REAL(1.0), REAL(2.0), REAL(3.0),
													REAL(3.0), REAL(2.0), REAL(5.0) });

		Vector<Real> b;
		REQUIRE_THROWS_AS(b = a.VectorFromDiagonal(), MatrixDimensionError);
	}

	///////////////////////               Matrix properties            //////////////////////
	TEST_CASE("Matrix::IsUnit", "[simple]") {
			TEST_PRECISION_INFO();
		Matrix<Real> a(2, 2, {REAL(1.0), REAL(0.0), REAL(0.0), REAL(1.0)});

		REQUIRE(true == a.IsUnit());

		Matrix<Real> b(2, 2, { REAL(1.0), REAL(0.0), REAL(0.0), REAL(2.0) });

		REQUIRE(false == b.IsUnit());
	}
	TEST_CASE("Matrix::IsDiagonal", "[simple]") {
			TEST_PRECISION_INFO();
		Matrix<Real> a(2, 2, { REAL(1.0), REAL(0.0), REAL(0.0), REAL(5.0) });

		REQUIRE(true == a.IsDiagonal());

		Matrix<Real> b(2, 2, { REAL(1.0), REAL(2.0), REAL(0.0), REAL(2.0) });

		REQUIRE(false == b.IsDiagonal());
	}

	///////////////////////             Assignment operators           //////////////////////

	///////////////////////               Access operators             //////////////////////
	
	///////////////////////              Equality operations           //////////////////////
	TEST_CASE("Matrix::IsEqual_diff_size_matrices", "[simple]") {
			TEST_PRECISION_INFO();
		Matrix<Real> a(2, 3), b(4, 5);

		REQUIRE(false == a.IsEqualTo(b));
	}
	TEST_CASE("Matrix::IsEqual", "[simple]") {
			TEST_PRECISION_INFO();
		Matrix<Real> a(2, 2, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });
		Matrix<Real> b(2, 2, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });

		REQUIRE(true == a.IsEqualTo(b));
	}
	TEST_CASE("Matrix::IsEqual2", "[simple]") {
			TEST_PRECISION_INFO();
		Matrix<Real> a(2, 2, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });
		Matrix<Real> b(2, 2, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(5.0) });

		REQUIRE(false == a.IsEqualTo(b));
	}
	TEST_CASE("Matrix::IsEqual3", "[simple]") {
			TEST_PRECISION_INFO();
		Matrix<Real> a(2, 2, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });
		Matrix<Real> b(2, 2, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0001) });

		// Difference is 1e-4, use 1e-3 (yes) and 1e-7 (no)
		// 1e-3 scales to 1e-1 (float), 1e-7 scales to 1e-5 (float)
		REQUIRE(true == a.IsEqualTo(b, ScaleTolerance(REAL(1e-3))));
		REQUIRE(false == a.IsEqualTo(b, ScaleTolerance(REAL(1e-7))));
	}

	///////////////////////              Arithmetic operators          //////////////////////
	TEST_CASE("Matrix::Op+-", "[simple]") {
			TEST_PRECISION_INFO();
		Matrix<Real> a(2, 2, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });
		Matrix<Real> b(2, 2, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });

		auto c = a + b;
		auto d = a - b;

		REQUIRE(REAL(2.0) == c[0][0]);
		REQUIRE(REAL(4.0) == c[0][1]);

		REQUIRE(REAL(0.0) == d[0][0]);
		REQUIRE(REAL(0.0) == d[0][1]);
	}
	TEST_CASE("Matrix::Op*", "[simple]") {
			TEST_PRECISION_INFO();
		Matrix<Real> a(2, 2, { REAL(1.0), REAL(2.0),
																REAL(3.0), REAL(4.0) });
		Matrix<Real> b(2, 2, { REAL(1.0), REAL(2.0),
																REAL(3.0), REAL(4.0) });
		auto c = a * b;

		REQUIRE(REAL(7.0) == c[0][0]);
		REQUIRE(REAL(10.0) == c[0][1]);
		REQUIRE(REAL(15.0) == c[1][0]);
		REQUIRE(REAL(22.0) == c[1][1]);
	}
	TEST_CASE("Matrix::mul_double", "[simple]") {
			TEST_PRECISION_INFO();
		Matrix<Real> a(2, 2, { REAL(1.0), REAL(100.0), REAL(50.0), REAL(100.0) });

		auto b = a * REAL(2.0);
		auto c = REAL(2.0) * a;

		REQUIRE(REAL(2.0) == b[0][0]);
		REQUIRE(REAL(2.0) == c[0][0]);

		REQUIRE(REAL(200.0) == b[0][1]);
		REQUIRE(REAL(200.0) == c[0][1]);
	}
	TEST_CASE("Matrix::div_double", "[simple]") {
			TEST_PRECISION_INFO();
		Matrix<Real> a(2, 2, { REAL(4.0), REAL(400.0), REAL(1.0), 1 });

		auto b = a / REAL(2.0);

		REQUIRE(REAL(2.0) == b[0][0]);
		REQUIRE(REAL(200.0) == b[0][1]);
	}
	TEST_CASE("Matrix::mul_Vector_right", "[simple]") {
			TEST_PRECISION_INFO();
		Matrix<Real> a(2, 2, { REAL(1.0), REAL(10.0), REAL(5.0), REAL(2.0) });
		Vector<Real> b({ REAL(1.0), REAL(2.0) });

		auto c = a * b;

		REQUIRE(REAL(21.0) == c[0]);
		REQUIRE(REAL(9.0) == c[1]);
	}
	TEST_CASE("Matrix::mul_Vector_left", "[simple]") {
			TEST_PRECISION_INFO();
		Matrix<Real> a(2, 2, { REAL(1.0), REAL(10.0), REAL(5.0), REAL(2.0) });
		Vector<Real> b({ REAL(1.0), REAL(2.0) });

		auto d = b * a;

		REQUIRE(REAL(11.0) == d[0]);
		REQUIRE(REAL(14.0) == d[1]);
	}

	///////////////////////            Trace, Inverse & Transpose      //////////////////////
	TEST_CASE("Matrix::Transpose", "[simple]")
	{
			TEST_PRECISION_INFO();
		Matrix<Real> mat(2, 2, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });
		Matrix<Real> matTransp(2, 2, { REAL(1.0), REAL(3.0), REAL(2.0), REAL(4.0) });

		mat.Transpose();

		REQUIRE(mat.IsEqualTo(matTransp));
	}
	TEST_CASE("Matrix::GetTranspose", "[simple]")
	{
			TEST_PRECISION_INFO();
		Matrix<Real> mat(2, 2, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });
		Matrix<Real> matTransp(2, 2, { REAL(1.0), REAL(3.0), REAL(2.0), REAL(4.0) });

		auto trans = mat.GetTranspose();

		REQUIRE(trans.IsEqualTo(matTransp));

		Matrix<Real> mat1(2, 3, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0), REAL(6.0) });
		Matrix<Real> matTransp1(3, 2, { REAL(1.0), REAL(4.0), REAL(2.0), REAL(5.0), REAL(3.0), REAL(6.0) });

		trans = mat1.GetTranspose();

		REQUIRE(trans.IsEqualTo(matTransp1));

		// test transposing row matrix
		Matrix<Real> mat2(1, 4, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });
		Matrix<Real> matTransp2(4, 1, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });

		auto mat22 = mat2.GetTranspose();

		REQUIRE(mat22.IsEqualTo(matTransp2));

		// test transposing column matrix
		Matrix<Real> mat3(4, 1, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });
		Matrix<Real> matTransp3(1, 4, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });

		auto mat33 = mat3.GetTranspose();

		REQUIRE(mat33.IsEqualTo(matTransp3));
	}
	TEST_CASE("Matrix::GetInverse", "[simple]")
	{
			TEST_PRECISION_INFO();
		Matrix<Real> mat(2, 2, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });

		auto b = mat.GetInverse();
		auto c = mat * b;

		REQUIRE(c.IsEqualTo(Matrix<Real>::GetUnitMatrix(2)));
	}
	TEST_CASE("Matrix::GetInverse_almost_signular", "[simple]")
	{
			TEST_PRECISION_INFO();
		Matrix<Real> mat(3, 3, { REAL(1.0), REAL(2.0), REAL(1.0),
														 REAL(1.0), REAL(2.0), REAL(1.00001),
														-REAL(4.0), REAL(6.0), REAL(4.0)});

		auto b = mat.GetInverse();
		auto c = mat * b;

		REQUIRE(c.IsEqualTo(Matrix<Real>::GetUnitMatrix(3), 1e-10));
		REQUIRE(!c.IsEqualTo(Matrix<Real>::GetUnitMatrix(3), 1e-11));
	}
	TEST_CASE("Matrix::GetInverse_throws_for_singular_matrix", "[simple]")
	{
			TEST_PRECISION_INFO();
		// test case for singular matrices with Invert
		Matrix<Real> mat2(2, 2, { REAL(1.0), REAL(2.0), REAL(2.0), REAL(4.0) });

		REQUIRE_THROWS_AS(mat2.GetInverse(), SingularMatrixError);
	}
	
	///////////////////////////               I/O                 ///////////////////////////
	TEST_CASE("Matrix::to_string", "[simple]")
	{
			TEST_PRECISION_INFO();
		Matrix<Real> mat(2, 2, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });

		REQUIRE("Rows: 2 Cols: 2\n[     1,     2 ]\n[     3,     4 ]" == mat.to_string(5, 3));
	}
	TEST_CASE("Matrix::Print_with_zerotreshold", "[simple]")
	{
			TEST_PRECISION_INFO();
		std::stringstream str;

		Matrix<Real> mat(2, 2, { REAL(1.0), REAL(2.0), REAL(3.0), 1e-10 });
		mat.Print(str, 5, 3, 1e-9);

		REQUIRE("Rows: 2 Cols: 2\n[     1,     2 ]\n[     3,     0 ]" == str.str());
	}
	TEST_CASE("Matrix::Test exceptions", "[simple]")
	{
			TEST_PRECISION_INFO();
		REQUIRE_THROWS_AS(Matrix<Real>(-1, 5), MatrixDimensionError);


		Matrix<Real> mat_21(3, 3);
		Matrix<Real> mat_22(2, 2);
		REQUIRE_THROWS_AS(mat_21 + mat_22, MatrixDimensionError);
		REQUIRE_THROWS_AS(mat_21 - mat_22, MatrixDimensionError);
		REQUIRE_THROWS_AS(mat_21 * mat_22, MatrixDimensionError);

		Matrix<Real> mat(2, 2, { REAL(1.0), REAL(0.0), REAL(0.0), REAL(1.0) });
		Vector<Real> a{ REAL(1.0), REAL(0.0), REAL(0.0) };

		REQUIRE_THROWS_AS(a * mat, MatrixDimensionError);
		REQUIRE_THROWS_AS(mat * a, MatrixDimensionError);

		// invert i transpose
		Matrix<Real> mat_3(3, 4);
		REQUIRE_THROWS_AS(mat_3.Invert(), MatrixDimensionError);
		REQUIRE_THROWS_AS(mat_3.GetInverse(), MatrixDimensionError);
		REQUIRE_THROWS_AS(mat_3.Transpose(), MatrixDimensionError);
	}

	///////////////////////          Formatted Output Tests       //////////////////////
	TEST_CASE("Matrix::Print_with_default_format", "[matrix_format]") {
			TEST_PRECISION_INFO();
		Matrix<Real> a(2, 3, { REAL(1.234), REAL(2.567), REAL(3.891),
		                       REAL(4.123), REAL(5.678), REAL(6.901) });

		std::stringstream ss;
		a.Print(ss, MatrixPrintFormat::Default());
		std::string output = ss.str();

		// Should contain header
		REQUIRE(output.find("Rows: 2") != std::string::npos);
		REQUIRE(output.find("Cols: 3") != std::string::npos);

		// Should contain brackets
		REQUIRE(output.find("[") != std::string::npos);
		REQUIRE(output.find("]") != std::string::npos);

		// Should contain values (check one)
		REQUIRE(output.find("1.234") != std::string::npos);
	}

	TEST_CASE("Matrix::Print_with_scientific_format", "[matrix_format]") {
			TEST_PRECISION_INFO();
		Matrix<Real> a(2, 2, { REAL(1234.567), REAL(0.000123),
		                       REAL(9876543.21), REAL(0.0000000987) });

		std::stringstream ss;
		a.Print(ss, MatrixPrintFormat::Scientific());
		std::string output = ss.str();

		// Should use scientific notation (contains 'e' for exponent)
		REQUIRE(output.find("e") != std::string::npos);
	}

	TEST_CASE("Matrix::Print_with_compact_format", "[matrix_format]") {
			TEST_PRECISION_INFO();
		Matrix<Real> a(2, 2, { REAL(1.1), REAL(2.2), REAL(3.3), REAL(4.4) });

		std::stringstream ss;
		a.Print(ss, MatrixPrintFormat::Compact());
		std::string output = ss.str();

		// Compact mode should be single line
		size_t newline_count = 0;
		for (char c : output) {
			if (c == '\n') newline_count++;
		}
		// Only one newline at the end
		REQUIRE(newline_count == 1);
	}

	TEST_CASE("Matrix::Print_with_high_precision", "[matrix_format]") {
			TEST_PRECISION_INFO();
		Matrix<Real> a(2, 2, { Constants::PI, Constants::E,
		                       std::sqrt(REAL(2.0)), std::sqrt(REAL(3.0)) });

		std::stringstream ss;
		a.Print(ss, MatrixPrintFormat::HighPrecision());
		std::string output = ss.str();

		// High precision should show more decimal places
		// Check for PI with at least 8 decimals
		REQUIRE(output.find("3.14159265") != std::string::npos);
	}

	TEST_CASE("Matrix::Print_with_no_delimiter", "[matrix_format]") {
			TEST_PRECISION_INFO();
		Matrix<Real> a(2, 2, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });

		std::stringstream ss;
		a.Print(ss, MatrixPrintFormat::NoDelimiter());
		std::string output = ss.str();

		// Should not contain commas
		REQUIRE(output.find(",") == std::string::npos);
	}

	TEST_CASE("Matrix::Print_custom_format", "[matrix_format]") {
			TEST_PRECISION_INFO();
		Matrix<Real> a(3, 3, { REAL(1.0), REAL(2.0), REAL(3.0),
		                       REAL(4.0), REAL(5.0), REAL(6.0),
		                       REAL(7.0), REAL(8.0), REAL(9.0) });

		MatrixPrintFormat fmt;
		fmt.width = 5;
		fmt.precision = 1;
		fmt.showHeader = false;
		fmt.showBrackets = false;
		fmt.delimiter = " | ";

		std::stringstream ss;
		a.Print(ss, fmt);
		std::string output = ss.str();

		// Should not have header
		REQUIRE(output.find("Rows:") == std::string::npos);

		// Should not have brackets
		REQUIRE(output.find("[") == std::string::npos);

		// Should have custom delimiter
		REQUIRE(output.find("|") != std::string::npos);
	}

	TEST_CASE("Matrix::Print_single_element", "[matrix_format]") {
			TEST_PRECISION_INFO();
		// Test 1x1 matrix (smallest valid matrix)
		Matrix<Real> a(1, 1, { REAL(42.0) });

		std::stringstream ss;
		a.Print(ss, MatrixPrintFormat::Default());
		std::string output = ss.str();

		// Should contain the value
		REQUIRE(output.find("42") != std::string::npos);
		
		// Should have header
		REQUIRE(output.find("Rows: 1") != std::string::npos);
		REQUIRE(output.find("Cols: 1") != std::string::npos);
	}

	TEST_CASE("Matrix::Print_backward_compatibility", "[matrix_format]") {
			TEST_PRECISION_INFO();
		Matrix<Real> a(2, 2, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });

		// Old-style Print should still work
		std::stringstream ss1;
		a.Print(ss1, 10, 3);

		std::stringstream ss2;
		a.Print(ss2, 10, 3, 1e-10);  // with zero threshold

		// Both should produce valid output
		REQUIRE(!ss1.str().empty());
		REQUIRE(!ss2.str().empty());

		// operator<< should still work
		std::stringstream ss3;
		ss3 << a;
		REQUIRE(!ss3.str().empty());
	}

	// ========== EDGE CASE TESTS ==========

	TEST_CASE("Matrix::EdgeCase_1x1_matrix", "[Matrix][edge_cases]")
	{
		TEST_PRECISION_INFO();
		Matrix<Real> m1(1, 1, {REAL(5.0)});
		Matrix<Real> m2(1, 1, {REAL(3.0)});
		
		REQUIRE(m1.RowNum() == 1);
		REQUIRE(m1.ColNum() == 1);
		REQUIRE(m1(0, 0) == REAL(5.0));
		
		// Arithmetic
		Matrix<Real> sum = m1 + m2;
		REQUIRE(sum(0, 0) == REAL(8.0));
		
		Matrix<Real> prod = m1 * m2;
		REQUIRE(prod(0, 0) == REAL(15.0));
		
		// Inverse of 1x1 [a] is [1/a]
		Matrix<Real> inv = m1.GetInverse();
		REQUIRE_THAT(inv(0, 0), Catch::Matchers::WithinAbs(REAL(0.2), 1e-10));
		
		// m1 * inv = I (1x1 identity)
		Matrix<Real> identity = m1 * inv;
		REQUIRE_THAT(identity(0, 0), Catch::Matchers::WithinAbs(REAL(1.0), 1e-10));
	}

	TEST_CASE("Matrix::EdgeCase_rectangular_multiply", "[Matrix][edge_cases]")
	{
		TEST_PRECISION_INFO();
		// (2x3) * (3x4) = (2x4)
		Matrix<Real> a(2, 3, {REAL(1.0), REAL(2.0), REAL(3.0),
		                      REAL(4.0), REAL(5.0), REAL(6.0)});
		Matrix<Real> b(3, 4, {REAL(1.0), REAL(0.0), REAL(0.0), REAL(0.0),
		                      REAL(0.0), REAL(1.0), REAL(0.0), REAL(0.0),
		                      REAL(0.0), REAL(0.0), REAL(1.0), REAL(0.0)});
		
		Matrix<Real> c = a * b;
		REQUIRE(c.RowNum() == 2);
		REQUIRE(c.ColNum() == 4);
		
		// First row of result
		REQUIRE(c(0, 0) == REAL(1.0));
		REQUIRE(c(0, 1) == REAL(2.0));
		REQUIRE(c(0, 2) == REAL(3.0));
		REQUIRE(c(0, 3) == REAL(0.0));
	}

	TEST_CASE("Matrix::EdgeCase_identity_properties", "[Matrix][edge_cases]")
	{
		TEST_PRECISION_INFO();
		Matrix<Real> I = Matrix<Real>::GetUnitMatrix(3);
		Matrix<Real> A(3, 3, {REAL(1.0), REAL(2.0), REAL(3.0),
		                      REAL(4.0), REAL(5.0), REAL(6.0),
		                      REAL(7.0), REAL(8.0), REAL(10.0)});
		
		// A * I = A
		Matrix<Real> AI = A * I;
		REQUIRE(A.IsEqualTo(AI, 1e-10));
		
		// I * A = A
		Matrix<Real> IA = I * A;
		REQUIRE(A.IsEqualTo(IA, 1e-10));
		
		// I^-1 = I
		Matrix<Real> Iinv = I.GetInverse();
		REQUIRE(I.IsEqualTo(Iinv, 1e-10));
	}

} // namespace MML::Tests::Matrix
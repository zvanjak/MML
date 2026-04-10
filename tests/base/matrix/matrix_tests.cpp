#include <catch2/catch_all.hpp>
#include "../../TestPrecision.h"
#include "../../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Vector/Vector.h"
#include "base/Matrix/Matrix.h"
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

		REQUIRE(2 == a.rows());
		REQUIRE(2 == a.cols());

		REQUIRE(REAL(0.0) == a(0, 0));
		REQUIRE(REAL(0.0) == a(0, 1));
		REQUIRE(REAL(0.0) == a(1, 0));
		REQUIRE(REAL(0.0) == a(1, 1));

		Matrix<Real> b(1000, 1000);

		REQUIRE(1000 == b.rows());
		REQUIRE(1000 == b.cols());

		REQUIRE(REAL(0.0) == b(0, 0));
		REQUIRE(REAL(0.0) == b(158, 738));
		REQUIRE(REAL(0.0) == b(35, 999));
		REQUIRE(REAL(0.0) == b(0, 100));
	}
	TEST_CASE("Matrix::default_ctor_init_to_value", "[simple]") {
			TEST_PRECISION_INFO();
		Matrix<Real> a(2, 2, REAL(5.0));

		REQUIRE(2 == a.rows());
		REQUIRE(2 == a.cols());

		REQUIRE(REAL(5.0) == a(0, 0));
		REQUIRE(REAL(5.0) == a(0, 1));
		REQUIRE(REAL(5.0) == a(1, 0));
		REQUIRE(REAL(5.0) == a(1, 1));
	}
	TEST_CASE("Matrix::initializer_list_ctor", "[simple]") {
			TEST_PRECISION_INFO();
		Matrix<Real> a(2, 2, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });

		REQUIRE(2 == a.rows());
		REQUIRE(2 == a.cols());

		REQUIRE(REAL(1.0) == a[0][0]);
		REQUIRE(REAL(2.0) == a(0, 1));
		REQUIRE(REAL(3.0) == a(1, 0));
		REQUIRE(REAL(4.0) == a(1, 1));
	}
	TEST_CASE("Matrix::double_ptr_ctor", "[simple]") {
			TEST_PRECISION_INFO();
		Real vec[4] = { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) };
		Matrix<Real> a(2, 2, vec);

		REQUIRE(2 == a.rows());
		REQUIRE(2 == a.cols());

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

		REQUIRE(3 == a.rows());
		REQUIRE(5 == a.cols());
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
	TEST_CASE("Matrix::Identity", "[simple]") {
			TEST_PRECISION_INFO();
		auto a = Matrix<Real>::Identity(2);

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

		auto b = a.diagonal();

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
		REQUIRE_THROWS_AS(b = a.diagonal(), MatrixDimensionError);
	}

	///////////////////////               Matrix properties            //////////////////////
	TEST_CASE("Matrix::IsUnit", "[simple]") {
			TEST_PRECISION_INFO();
		Matrix<Real> a(2, 2, {REAL(1.0), REAL(0.0), REAL(0.0), REAL(1.0)});

		REQUIRE(true == a.isIdentity());

		Matrix<Real> b(2, 2, { REAL(1.0), REAL(0.0), REAL(0.0), REAL(2.0) });

		REQUIRE(false == b.isIdentity());
	}
	TEST_CASE("Matrix::IsDiagonal", "[simple]") {
			TEST_PRECISION_INFO();
		Matrix<Real> a(2, 2, { REAL(1.0), REAL(0.0), REAL(0.0), REAL(5.0) });

		REQUIRE(true == a.isDiagonal());

		Matrix<Real> b(2, 2, { REAL(1.0), REAL(2.0), REAL(0.0), REAL(2.0) });

		REQUIRE(false == b.isDiagonal());
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
		REQUIRE(true == a.IsEqualTo(b, TOL3(1e-3, 1e-1, 1e-3)));
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

		auto trans = mat.transpose();

		REQUIRE(trans.IsEqualTo(matTransp));

		Matrix<Real> mat1(2, 3, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0), REAL(6.0) });
		Matrix<Real> matTransp1(3, 2, { REAL(1.0), REAL(4.0), REAL(2.0), REAL(5.0), REAL(3.0), REAL(6.0) });

		trans = mat1.transpose();

		REQUIRE(trans.IsEqualTo(matTransp1));

		// test transposing row matrix
		Matrix<Real> mat2(1, 4, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });
		Matrix<Real> matTransp2(4, 1, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });

		auto mat22 = mat2.transpose();

		REQUIRE(mat22.IsEqualTo(matTransp2));

		// test transposing column matrix
		Matrix<Real> mat3(4, 1, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });
		Matrix<Real> matTransp3(1, 4, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });

		auto mat33 = mat3.transpose();

		REQUIRE(mat33.IsEqualTo(matTransp3));
	}
	TEST_CASE("Matrix::GetInverse", "[simple]")
	{
			TEST_PRECISION_INFO();
		Matrix<Real> mat(2, 2, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });

		auto b = mat.inverse();
		auto c = mat * b;

		REQUIRE(c.IsEqualTo(Matrix<Real>::Identity(2)));
	}
	TEST_CASE("Matrix::GetInverse_almost_signular", "[simple]")
	{
			TEST_PRECISION_INFO();
		Matrix<Real> mat(3, 3, { REAL(1.0), REAL(2.0), REAL(1.0),
														 REAL(1.0), REAL(2.0), REAL(1.00001),
														-REAL(4.0), REAL(6.0), REAL(4.0)});

		auto b = mat.inverse();
		auto c = mat * b;

		REQUIRE(c.IsEqualTo(Matrix<Real>::Identity(3), TOL(1e-10, 5e-1)));
	}
	TEST_CASE("Matrix::GetInverse_throws_for_singular_matrix", "[simple]")
	{
			TEST_PRECISION_INFO();
		// test case for singular matrices with Invert
		Matrix<Real> mat2(2, 2, { REAL(1.0), REAL(2.0), REAL(2.0), REAL(4.0) });

		REQUIRE_THROWS_AS(mat2.inverse(), SingularMatrixError);
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

		Matrix<Real> mat(2, 2, { REAL(1.0), REAL(2.0), REAL(3.0), TOL(1e-10, 1e-5) });
		mat.Print(str, 5, 3, TOL(1e-9, 1e-4));

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
		REQUIRE_THROWS_AS(mat_3.inverse(), MatrixDimensionError);
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
		if constexpr (std::is_same_v<Real, float>) {
			REQUIRE(output.find("3.14159") != std::string::npos);
		} else {
			REQUIRE(output.find("3.14159265") != std::string::npos);
		}
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
		a.Print(ss2, 10, 3, TOL(1e-10, 1e-5));  // with zero threshold

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
		
		REQUIRE(m1.rows() == 1);
		REQUIRE(m1.cols() == 1);
		REQUIRE(m1(0, 0) == REAL(5.0));
		
		// Arithmetic
		Matrix<Real> sum = m1 + m2;
		REQUIRE(sum(0, 0) == REAL(8.0));
		
		Matrix<Real> prod = m1 * m2;
		REQUIRE(prod(0, 0) == REAL(15.0));
		
		// Inverse of 1x1 [a] is [1/a]
		Matrix<Real> inv = m1.inverse();
		REQUIRE_THAT(inv(0, 0), Catch::Matchers::WithinAbs(REAL(0.2), TOL(1e-10, 1e-5)));
		
		// m1 * inv = I (1x1 identity)
		Matrix<Real> identity = m1 * inv;
		REQUIRE_THAT(identity(0, 0), Catch::Matchers::WithinAbs(REAL(1.0), TOL(1e-10, 1e-5)));
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
		REQUIRE(c.rows() == 2);
		REQUIRE(c.cols() == 4);
		
		// First row of result
		REQUIRE(c(0, 0) == REAL(1.0));
		REQUIRE(c(0, 1) == REAL(2.0));
		REQUIRE(c(0, 2) == REAL(3.0));
		REQUIRE(c(0, 3) == REAL(0.0));
	}

	TEST_CASE("Matrix::EdgeCase_identity_properties", "[Matrix][edge_cases]")
	{
		TEST_PRECISION_INFO();
		Matrix<Real> I = Matrix<Real>::Identity(3);
		Matrix<Real> A(3, 3, {REAL(1.0), REAL(2.0), REAL(3.0),
		                      REAL(4.0), REAL(5.0), REAL(6.0),
		                      REAL(7.0), REAL(8.0), REAL(10.0)});
		
		// A * I = A
		Matrix<Real> AI = A * I;
		REQUIRE(A.IsEqualTo(AI, TOL(1e-10, 1e-5)));
		
		// I * A = A
		Matrix<Real> IA = I * A;
		REQUIRE(A.IsEqualTo(IA, TOL(1e-10, 1e-5)));
		
		// I^-1 = I
		Matrix<Real> Iinv = I.inverse();
		REQUIRE(I.IsEqualTo(Iinv, TOL(1e-10, 1e-5)));
	}

	///////////////////////          Copy and Move Semantics          //////////////////////
	TEST_CASE("Matrix::copy_constructor", "[Matrix][copy_move]")
	{
		TEST_PRECISION_INFO();
		Matrix<Real> original(3, 3, { REAL(1.0), REAL(2.0), REAL(3.0),
		                              REAL(4.0), REAL(5.0), REAL(6.0),
		                              REAL(7.0), REAL(8.0), REAL(9.0) });
		
		Matrix<Real> copy(original);
		
		REQUIRE(copy.rows() == 3);
		REQUIRE(copy.cols() == 3);
		REQUIRE(copy.IsEqualTo(original));
		
		// Modify copy, original should be unchanged
		copy(0, 0) = REAL(999.0);
		REQUIRE(original(0, 0) == REAL(1.0));
	}

	TEST_CASE("Matrix::move_constructor", "[Matrix][copy_move]")
	{
		TEST_PRECISION_INFO();
		Matrix<Real> original(3, 3, { REAL(1.0), REAL(2.0), REAL(3.0),
		                              REAL(4.0), REAL(5.0), REAL(6.0),
		                              REAL(7.0), REAL(8.0), REAL(9.0) });
		
		Matrix<Real> moved(std::move(original));
		
		REQUIRE(moved.rows() == 3);
		REQUIRE(moved.cols() == 3);
		REQUIRE(moved(0, 0) == REAL(1.0));
		REQUIRE(moved(2, 2) == REAL(9.0));
		
		// Original should be in valid empty state
		REQUIRE(original.rows() == 0);
		REQUIRE(original.cols() == 0);
	}

	TEST_CASE("Matrix::copy_assignment", "[Matrix][copy_move]")
	{
		TEST_PRECISION_INFO();
		Matrix<Real> original(2, 2, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });
		Matrix<Real> target(3, 3);  // Different size
		
		target = original;
		
		REQUIRE(target.rows() == 2);
		REQUIRE(target.cols() == 2);
		REQUIRE(target.IsEqualTo(original));
		
		// Modify target, original unchanged
		target(0, 0) = REAL(100.0);
		REQUIRE(original(0, 0) == REAL(1.0));
	}

	TEST_CASE("Matrix::move_assignment", "[Matrix][copy_move]")
	{
		TEST_PRECISION_INFO();
		Matrix<Real> original(2, 2, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });
		Matrix<Real> target(5, 5);  // Different size
		
		target = std::move(original);
		
		REQUIRE(target.rows() == 2);
		REQUIRE(target.cols() == 2);
		REQUIRE(target(1, 1) == REAL(4.0));
		
		// Original should be empty
		REQUIRE(original.rows() == 0);
		REQUIRE(original.cols() == 0);
	}

	///////////////////////          Additional Constructors           //////////////////////
	TEST_CASE("Matrix::nested_vector_constructor", "[Matrix][constructors]")
	{
		TEST_PRECISION_INFO();
		std::vector<std::vector<Real>> data = {
			{ REAL(1.0), REAL(2.0), REAL(3.0) },
			{ REAL(4.0), REAL(5.0), REAL(6.0) }
		};
		
		Matrix<Real> m(data);
		
		REQUIRE(m.rows() == 2);
		REQUIRE(m.cols() == 3);
		REQUIRE(m(0, 0) == REAL(1.0));
		REQUIRE(m(0, 2) == REAL(3.0));
		REQUIRE(m(1, 0) == REAL(4.0));
		REQUIRE(m(1, 2) == REAL(6.0));
	}

	TEST_CASE("Matrix::nested_vector_constructor_inconsistent_rows_throws", "[Matrix][constructors]")
	{
		TEST_PRECISION_INFO();
		std::vector<std::vector<Real>> badData = {
			{ REAL(1.0), REAL(2.0), REAL(3.0) },
			{ REAL(4.0), REAL(5.0) }  // Different size!
		};
		
		REQUIRE_THROWS_AS(Matrix<Real>(badData), MatrixDimensionError);
	}

	TEST_CASE("Matrix::std_array_constructor", "[Matrix][constructors]")
	{
		TEST_PRECISION_INFO();
		std::array<std::array<Real, 3>, 2> data = {{
			{{ REAL(1.0), REAL(2.0), REAL(3.0) }},
			{{ REAL(4.0), REAL(5.0), REAL(6.0) }}
		}};
		
		Matrix<Real> m(data);
		
		REQUIRE(m.rows() == 2);
		REQUIRE(m.cols() == 3);
		REQUIRE(m(0, 0) == REAL(1.0));
		REQUIRE(m(1, 2) == REAL(6.0));
	}

	TEST_CASE("Matrix::submatrix_constructor", "[Matrix][constructors]")
	{
		TEST_PRECISION_INFO();
		Matrix<Real> source(4, 4, { REAL(1.0),  REAL(2.0),  REAL(3.0),  REAL(4.0),
		                            REAL(5.0),  REAL(6.0),  REAL(7.0),  REAL(8.0),
		                            REAL(9.0),  REAL(10.0), REAL(11.0), REAL(12.0),
		                            REAL(13.0), REAL(14.0), REAL(15.0), REAL(16.0) });
		
		// Extract 2x2 submatrix starting at (1,1)
		Matrix<Real> sub(source, 1, 1, 2, 2);
		
		REQUIRE(sub.rows() == 2);
		REQUIRE(sub.cols() == 2);
		REQUIRE(sub(0, 0) == REAL(6.0));
		REQUIRE(sub(0, 1) == REAL(7.0));
		REQUIRE(sub(1, 0) == REAL(10.0));
		REQUIRE(sub(1, 1) == REAL(11.0));
	}

	TEST_CASE("Matrix::submatrix_constructor_out_of_bounds_throws", "[Matrix][constructors]")
	{
		TEST_PRECISION_INFO();
		Matrix<Real> source(3, 3, { REAL(1.0), REAL(2.0), REAL(3.0),
		                            REAL(4.0), REAL(5.0), REAL(6.0),
		                            REAL(7.0), REAL(8.0), REAL(9.0) });
		
		// Out of bounds extraction
		REQUIRE_THROWS_AS(Matrix<Real>(source, 2, 2, 3, 3), MatrixDimensionError);
		REQUIRE_THROWS_AS(Matrix<Real>(source, -1, 0, 2, 2), MatrixDimensionError);
	}

	///////////////////////          Static Factory Methods            //////////////////////
	TEST_CASE("Matrix::GetDiagonalMatrix", "[Matrix][factories]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> diag({ REAL(1.0), REAL(2.0), REAL(3.0) });
		
		Matrix<Real> m = Matrix<Real>::Diagonal(diag);
		
		REQUIRE(m.rows() == 3);
		REQUIRE(m.cols() == 3);
		REQUIRE(m(0, 0) == REAL(1.0));
		REQUIRE(m(1, 1) == REAL(2.0));
		REQUIRE(m(2, 2) == REAL(3.0));
		REQUIRE(m(0, 1) == REAL(0.0));
		REQUIRE(m(1, 2) == REAL(0.0));
		REQUIRE(m.isDiagonal());
	}

	///////////////////////          Triangular Extraction             //////////////////////
	TEST_CASE("Matrix::GetLower_with_diagonal", "[Matrix][extraction]")
	{
		TEST_PRECISION_INFO();
		Matrix<Real> m(3, 3, { REAL(1.0), REAL(2.0), REAL(3.0),
		                       REAL(4.0), REAL(5.0), REAL(6.0),
		                       REAL(7.0), REAL(8.0), REAL(9.0) });
		
		Matrix<Real> lower = m.GetLower(true);
		
		REQUIRE(lower(0, 0) == REAL(1.0));
		REQUIRE(lower(1, 0) == REAL(4.0));
		REQUIRE(lower(1, 1) == REAL(5.0));
		REQUIRE(lower(2, 0) == REAL(7.0));
		REQUIRE(lower(2, 1) == REAL(8.0));
		REQUIRE(lower(2, 2) == REAL(9.0));
		// Upper triangle should be zero
		REQUIRE(lower(0, 1) == REAL(0.0));
		REQUIRE(lower(0, 2) == REAL(0.0));
		REQUIRE(lower(1, 2) == REAL(0.0));
	}

	TEST_CASE("Matrix::GetLower_without_diagonal", "[Matrix][extraction]")
	{
		TEST_PRECISION_INFO();
		Matrix<Real> m(3, 3, { REAL(1.0), REAL(2.0), REAL(3.0),
		                       REAL(4.0), REAL(5.0), REAL(6.0),
		                       REAL(7.0), REAL(8.0), REAL(9.0) });
		
		Matrix<Real> lower = m.GetLower(false);  // Strictly lower
		
		// Diagonal should be zero
		REQUIRE(lower(0, 0) == REAL(0.0));
		REQUIRE(lower(1, 1) == REAL(0.0));
		REQUIRE(lower(2, 2) == REAL(0.0));
		// Lower elements preserved
		REQUIRE(lower(1, 0) == REAL(4.0));
		REQUIRE(lower(2, 0) == REAL(7.0));
		REQUIRE(lower(2, 1) == REAL(8.0));
	}

	TEST_CASE("Matrix::GetUpper_with_diagonal", "[Matrix][extraction]")
	{
		TEST_PRECISION_INFO();
		Matrix<Real> m(3, 3, { REAL(1.0), REAL(2.0), REAL(3.0),
		                       REAL(4.0), REAL(5.0), REAL(6.0),
		                       REAL(7.0), REAL(8.0), REAL(9.0) });
		
		Matrix<Real> upper = m.GetUpper(true);
		
		REQUIRE(upper(0, 0) == REAL(1.0));
		REQUIRE(upper(0, 1) == REAL(2.0));
		REQUIRE(upper(0, 2) == REAL(3.0));
		REQUIRE(upper(1, 1) == REAL(5.0));
		REQUIRE(upper(1, 2) == REAL(6.0));
		REQUIRE(upper(2, 2) == REAL(9.0));
		// Lower triangle should be zero
		REQUIRE(upper(1, 0) == REAL(0.0));
		REQUIRE(upper(2, 0) == REAL(0.0));
		REQUIRE(upper(2, 1) == REAL(0.0));
	}

	TEST_CASE("Matrix::GetUpper_without_diagonal", "[Matrix][extraction]")
	{
		TEST_PRECISION_INFO();
		Matrix<Real> m(3, 3, { REAL(1.0), REAL(2.0), REAL(3.0),
		                       REAL(4.0), REAL(5.0), REAL(6.0),
		                       REAL(7.0), REAL(8.0), REAL(9.0) });
		
		Matrix<Real> upper = m.GetUpper(false);  // Strictly upper
		
		// Diagonal should be zero
		REQUIRE(upper(0, 0) == REAL(0.0));
		REQUIRE(upper(1, 1) == REAL(0.0));
		REQUIRE(upper(2, 2) == REAL(0.0));
		// Upper elements preserved
		REQUIRE(upper(0, 1) == REAL(2.0));
		REQUIRE(upper(0, 2) == REAL(3.0));
		REQUIRE(upper(1, 2) == REAL(6.0));
	}

	TEST_CASE("Matrix::GetLower_GetUpper_non_square_throws", "[Matrix][extraction]")
	{
		TEST_PRECISION_INFO();
		Matrix<Real> m(2, 3);
		
		REQUIRE_THROWS_AS(m.GetLower(), MatrixDimensionError);
		REQUIRE_THROWS_AS(m.GetUpper(), MatrixDimensionError);
	}

	TEST_CASE("Matrix::GetSubmatrix", "[Matrix][extraction]")
	{
		TEST_PRECISION_INFO();
		Matrix<Real> m(4, 4, { REAL(1.0),  REAL(2.0),  REAL(3.0),  REAL(4.0),
		                       REAL(5.0),  REAL(6.0),  REAL(7.0),  REAL(8.0),
		                       REAL(9.0),  REAL(10.0), REAL(11.0), REAL(12.0),
		                       REAL(13.0), REAL(14.0), REAL(15.0), REAL(16.0) });
		
		Matrix<Real> sub = m.GetSubmatrix(1, 2, 2, 2);
		
		REQUIRE(sub.rows() == 2);
		REQUIRE(sub.cols() == 2);
		REQUIRE(sub(0, 0) == REAL(7.0));
		REQUIRE(sub(0, 1) == REAL(8.0));
		REQUIRE(sub(1, 0) == REAL(11.0));
		REQUIRE(sub(1, 1) == REAL(12.0));
	}

	TEST_CASE("Matrix::GetDiagonal_alias", "[Matrix][extraction]")
	{
		TEST_PRECISION_INFO();
		Matrix<Real> m(3, 3, { REAL(1.0), REAL(2.0), REAL(3.0),
		                       REAL(4.0), REAL(5.0), REAL(6.0),
		                       REAL(7.0), REAL(8.0), REAL(9.0) });
		
		Vector<Real> diag = m.diagonal();
		
		REQUIRE(diag.size() == 3);
		REQUIRE(diag[0] == REAL(1.0));
		REQUIRE(diag[1] == REAL(5.0));
		REQUIRE(diag[2] == REAL(9.0));
	}

	///////////////////////          Row/Column Operations             //////////////////////
	TEST_CASE("Matrix::InitRowWithVector", "[Matrix][row_col_ops]")
	{
		TEST_PRECISION_INFO();
		Matrix<Real> m(3, 3);
		Vector<Real> row({ REAL(1.0), REAL(2.0), REAL(3.0) });
		
		m.InitRowWithVector(1, row);
		
		REQUIRE(m(1, 0) == REAL(1.0));
		REQUIRE(m(1, 1) == REAL(2.0));
		REQUIRE(m(1, 2) == REAL(3.0));
		// Other rows should be zero
		REQUIRE(m(0, 0) == REAL(0.0));
		REQUIRE(m(2, 2) == REAL(0.0));
	}

	TEST_CASE("Matrix::InitRowWithVector_throws", "[Matrix][row_col_ops]")
	{
		TEST_PRECISION_INFO();
		Matrix<Real> m(3, 3);
		Vector<Real> wrongSize({ REAL(1.0), REAL(2.0) });
		
		REQUIRE_THROWS_AS(m.InitRowWithVector(5, wrongSize), MatrixAccessBoundsError);
		REQUIRE_THROWS_AS(m.InitRowWithVector(0, wrongSize), MatrixDimensionError);
	}

	TEST_CASE("Matrix::InitColWithVector", "[Matrix][row_col_ops]")
	{
		TEST_PRECISION_INFO();
		Matrix<Real> m(3, 3);
		Vector<Real> col({ REAL(7.0), REAL(8.0), REAL(9.0) });
		
		m.InitColWithVector(2, col);
		
		REQUIRE(m(0, 2) == REAL(7.0));
		REQUIRE(m(1, 2) == REAL(8.0));
		REQUIRE(m(2, 2) == REAL(9.0));
	}

	TEST_CASE("Matrix::SwapRows", "[Matrix][row_col_ops]")
	{
		TEST_PRECISION_INFO();
		Matrix<Real> m(3, 3, { REAL(1.0), REAL(2.0), REAL(3.0),
		                       REAL(4.0), REAL(5.0), REAL(6.0),
		                       REAL(7.0), REAL(8.0), REAL(9.0) });
		
		m.SwapRows(0, 2);
		
		// Row 0 now has [7, 8, 9]
		REQUIRE(m(0, 0) == REAL(7.0));
		REQUIRE(m(0, 1) == REAL(8.0));
		REQUIRE(m(0, 2) == REAL(9.0));
		// Row 2 now has [1, 2, 3]
		REQUIRE(m(2, 0) == REAL(1.0));
		REQUIRE(m(2, 1) == REAL(2.0));
		REQUIRE(m(2, 2) == REAL(3.0));
		// Row 1 unchanged
		REQUIRE(m(1, 0) == REAL(4.0));
	}

	TEST_CASE("Matrix::SwapCols", "[Matrix][row_col_ops]")
	{
		TEST_PRECISION_INFO();
		Matrix<Real> m(3, 3, { REAL(1.0), REAL(2.0), REAL(3.0),
		                       REAL(4.0), REAL(5.0), REAL(6.0),
		                       REAL(7.0), REAL(8.0), REAL(9.0) });
		
		m.SwapCols(0, 2);
		
		// Col 0 now has [3, 6, 9]
		REQUIRE(m(0, 0) == REAL(3.0));
		REQUIRE(m(1, 0) == REAL(6.0));
		REQUIRE(m(2, 0) == REAL(9.0));
		// Col 2 now has [1, 4, 7]
		REQUIRE(m(0, 2) == REAL(1.0));
		REQUIRE(m(1, 2) == REAL(4.0));
		REQUIRE(m(2, 2) == REAL(7.0));
	}

	TEST_CASE("Matrix::SwapRows_same_row", "[Matrix][row_col_ops]")
	{
		TEST_PRECISION_INFO();
		Matrix<Real> m(2, 2, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });
		Matrix<Real> original(m);
		
		m.SwapRows(0, 0);  // No-op
		
		REQUIRE(m.IsEqualTo(original));
	}

	///////////////////////          Matrix Properties                 //////////////////////
	TEST_CASE("Matrix::IsDiagDominant", "[Matrix][properties]")
	{
		TEST_PRECISION_INFO();
		// Diagonally dominant: |a_ii| >= sum of |a_ij| for j != i
		Matrix<Real> dominant(3, 3, { REAL(10.0), REAL(1.0), REAL(2.0),
		                              REAL(1.0),  REAL(8.0), REAL(1.0),
		                              REAL(2.0),  REAL(1.0), REAL(9.0) });
		
		REQUIRE(dominant.isDiagonallyDominant() == true);
		
		Matrix<Real> notDominant(3, 3, { REAL(1.0), REAL(5.0), REAL(5.0),
		                                 REAL(1.0), REAL(1.0), REAL(1.0),
		                                 REAL(1.0), REAL(1.0), REAL(1.0) });
		
		REQUIRE(notDominant.isDiagonallyDominant() == false);
	}

	TEST_CASE("Matrix::IsSymmetric", "[Matrix][properties]")
	{
		TEST_PRECISION_INFO();
		Matrix<Real> symmetric(3, 3, { REAL(1.0), REAL(2.0), REAL(3.0),
		                               REAL(2.0), REAL(5.0), REAL(6.0),
		                               REAL(3.0), REAL(6.0), REAL(9.0) });
		
		REQUIRE(symmetric.isSymmetric() == true);
		
		Matrix<Real> notSymmetric(3, 3, { REAL(1.0), REAL(2.0), REAL(3.0),
		                                  REAL(4.0), REAL(5.0), REAL(6.0),
		                                  REAL(7.0), REAL(8.0), REAL(9.0) });
		
		REQUIRE(notSymmetric.isSymmetric() == false);
		
		// Non-square matrix cannot be symmetric
		Matrix<Real> nonSquare(2, 3);
		REQUIRE(nonSquare.isSymmetric() == false);
	}

	TEST_CASE("Matrix::IsAntiSymmetric", "[Matrix][properties]")
	{
		TEST_PRECISION_INFO();
		// Anti-symmetric: M(i,j) = -M(j,i)
		Matrix<Real> antiSym(3, 3, { REAL(0.0),  REAL(2.0),  REAL(3.0),
		                            -REAL(2.0),  REAL(0.0),  REAL(5.0),
		                            -REAL(3.0), -REAL(5.0),  REAL(0.0) });
		
		REQUIRE(antiSym.isAntiSymmetric() == true);
		
		Matrix<Real> notAntiSym(3, 3, { REAL(1.0), REAL(2.0), REAL(3.0),
		                                REAL(4.0), REAL(5.0), REAL(6.0),
		                                REAL(7.0), REAL(8.0), REAL(9.0) });
		
		REQUIRE(notAntiSym.isAntiSymmetric() == false);
	}

	TEST_CASE("Matrix::IsEmpty", "[Matrix][properties]")
	{
		TEST_PRECISION_INFO();
		Matrix<Real> empty;
		REQUIRE(empty.isEmpty() == true);
		
		Matrix<Real> notEmpty(2, 2);
		REQUIRE(notEmpty.isEmpty() == false);
	}

	///////////////////////          Matrix Norms                      //////////////////////
	TEST_CASE("Matrix::NormL1", "[Matrix][norms]")
	{
		TEST_PRECISION_INFO();
		Matrix<Real> m(2, 2, { REAL(1.0), REAL(-2.0), REAL(3.0), REAL(-4.0) });
		
		// L1 norm = sum of absolute values = 1 + 2 + 3 + 4 = 10
		REQUIRE_THAT(m.NormL1(), RealWithinRel(REAL(10.0)));
	}

	TEST_CASE("Matrix::NormL2_Frobenius", "[Matrix][norms]")
	{
		TEST_PRECISION_INFO();
		Matrix<Real> m(2, 2, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });
		
		// Frobenius norm = sqrt(1^2 + 2^2 + 3^2 + 4^2) = sqrt(30)
		Real expected = std::sqrt(REAL(30.0));
		REQUIRE_THAT(m.NormL2(), RealWithinRel(expected));
	}

	TEST_CASE("Matrix::NormLInf", "[Matrix][norms]")
	{
		TEST_PRECISION_INFO();
		Matrix<Real> m(2, 2, { REAL(1.0), REAL(-7.0), REAL(3.0), REAL(4.0) });
		
		// L-infinity norm = max absolute value = 7
		REQUIRE_THAT(m.NormLInf(), RealWithinRel(REAL(7.0)));
	}

	///////////////////////          Access and Iterators              //////////////////////
	TEST_CASE("Matrix::at_bounds_checked", "[Matrix][access]")
	{
		TEST_PRECISION_INFO();
		Matrix<Real> m(2, 2, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });
		
		REQUIRE(m.at(0, 0) == REAL(1.0));
		REQUIRE(m.at(1, 1) == REAL(4.0));
		
		// Mutable access
		m.at(0, 1) = REAL(99.0);
		REQUIRE(m(0, 1) == REAL(99.0));
		
		// Out of bounds throws
		REQUIRE_THROWS_AS(m.at(-1, 0), MatrixAccessBoundsError);
		REQUIRE_THROWS_AS(m.at(0, 5), MatrixAccessBoundsError);
		REQUIRE_THROWS_AS(m.at(10, 0), MatrixAccessBoundsError);
	}

	TEST_CASE("Matrix::data_access", "[Matrix][access]")
	{
		TEST_PRECISION_INFO();
		Matrix<Real> m(2, 2, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });
		
		Real* ptr = m.data();
		const Real* cptr = static_cast<const Matrix<Real>&>(m).data();
		
		REQUIRE(ptr[0] == REAL(1.0));
		REQUIRE(ptr[3] == REAL(4.0));
		REQUIRE(cptr[1] == REAL(2.0));
		
		// Modify via pointer
		ptr[0] = REAL(100.0);
		REQUIRE(m(0, 0) == REAL(100.0));
	}

	TEST_CASE("Matrix::iterators_flat", "[Matrix][iterators]")
	{
		TEST_PRECISION_INFO();
		Matrix<Real> m(2, 2, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });
		
		// Sum all elements using iterators
		Real sum = REAL(0.0);
		for (auto it = m.begin(); it != m.end(); ++it) {
			sum += *it;
		}
		REQUIRE_THAT(sum, RealWithinRel(REAL(10.0)));
		
		// Range-based for
		Real sum2 = REAL(0.0);
		for (Real val : m) {
			sum2 += val;
		}
		REQUIRE_THAT(sum2, RealWithinRel(REAL(10.0)));
		
		// cbegin/cend
		Real sum3 = REAL(0.0);
		for (auto it = m.cbegin(); it != m.cend(); ++it) {
			sum3 += *it;
		}
		REQUIRE_THAT(sum3, RealWithinRel(REAL(10.0)));
	}

	TEST_CASE("Matrix::row_iterator", "[Matrix][iterators]")
	{
		TEST_PRECISION_INFO();
		Matrix<Real> m(3, 3, { REAL(1.0), REAL(2.0), REAL(3.0),
		                       REAL(4.0), REAL(5.0), REAL(6.0),
		                       REAL(7.0), REAL(8.0), REAL(9.0) });
		
		// Iterate over row 1
		Real sum = REAL(0.0);
		for (Real val : m.row(1)) {
			sum += val;
		}
		REQUIRE_THAT(sum, RealWithinRel(REAL(15.0)));  // 4 + 5 + 6
	}

	TEST_CASE("Matrix::col_iterator", "[Matrix][iterators]")
	{
		TEST_PRECISION_INFO();
		Matrix<Real> m(3, 3, { REAL(1.0), REAL(2.0), REAL(3.0),
		                       REAL(4.0), REAL(5.0), REAL(6.0),
		                       REAL(7.0), REAL(8.0), REAL(9.0) });
		
		// Iterate over column 1
		Real sum = REAL(0.0);
		for (Real val : m.col(1)) {
			sum += val;
		}
		REQUIRE_THAT(sum, RealWithinRel(REAL(15.0)));  // 2 + 5 + 8
	}

	TEST_CASE("Matrix::block_view", "[Matrix][views]")
	{
		TEST_PRECISION_INFO();
		Matrix<Real> m(4, 4, { REAL(1.0),  REAL(2.0),  REAL(3.0),  REAL(4.0),
		                       REAL(5.0),  REAL(6.0),  REAL(7.0),  REAL(8.0),
		                       REAL(9.0),  REAL(10.0), REAL(11.0), REAL(12.0),
		                       REAL(13.0), REAL(14.0), REAL(15.0), REAL(16.0) });
		
		auto view = m.block(1, 1, 2, 2);
		
		REQUIRE(view.rows() == 2);
		REQUIRE(view.cols() == 2);
		REQUIRE(view(0, 0) == REAL(6.0));
		REQUIRE(view(0, 1) == REAL(7.0));
		REQUIRE(view(1, 0) == REAL(10.0));
		REQUIRE(view(1, 1) == REAL(11.0));
		
		// Modify through view
		view(0, 0) = REAL(100.0);
		REQUIRE(m(1, 1) == REAL(100.0));  // Original modified!
	}

	///////////////////////          Equality Operators                //////////////////////
	TEST_CASE("Matrix::operator_equals", "[Matrix][equality]")
	{
		TEST_PRECISION_INFO();
		Matrix<Real> a(2, 2, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });
		Matrix<Real> b(2, 2, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });
		Matrix<Real> c(2, 2, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(5.0) });
		
		REQUIRE(a == b);
		REQUIRE(!(a == c));
	}

	TEST_CASE("Matrix::operator_not_equals", "[Matrix][equality]")
	{
		TEST_PRECISION_INFO();
		Matrix<Real> a(2, 2, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });
		Matrix<Real> b(2, 2, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(5.0) });
		
		REQUIRE(a != b);
	}

	TEST_CASE("Matrix::AreEqual_static", "[Matrix][equality]")
	{
		TEST_PRECISION_INFO();
		Matrix<Real> a(2, 2, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });
		Matrix<Real> b(2, 2, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0001) });
		
		REQUIRE(Matrix<Real>::AreEqual(a, b, REAL(0.001)) == true);
		REQUIRE(Matrix<Real>::AreEqual(a, b, REAL(0.00001)) == false);
	}

	///////////////////////          Compound Operators                //////////////////////
	TEST_CASE("Matrix::operator_plus_equals", "[Matrix][compound]")
	{
		TEST_PRECISION_INFO();
		Matrix<Real> a(2, 2, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });
		Matrix<Real> b(2, 2, { REAL(10.0), REAL(20.0), REAL(30.0), REAL(40.0) });
		
		a += b;
		
		REQUIRE(a(0, 0) == REAL(11.0));
		REQUIRE(a(0, 1) == REAL(22.0));
		REQUIRE(a(1, 0) == REAL(33.0));
		REQUIRE(a(1, 1) == REAL(44.0));
	}

	TEST_CASE("Matrix::operator_minus_equals", "[Matrix][compound]")
	{
		TEST_PRECISION_INFO();
		Matrix<Real> a(2, 2, { REAL(10.0), REAL(20.0), REAL(30.0), REAL(40.0) });
		Matrix<Real> b(2, 2, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });
		
		a -= b;
		
		REQUIRE(a(0, 0) == REAL(9.0));
		REQUIRE(a(1, 1) == REAL(36.0));
	}

	TEST_CASE("Matrix::operator_times_equals_scalar", "[Matrix][compound]")
	{
		TEST_PRECISION_INFO();
		Matrix<Real> a(2, 2, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });
		
		a *= REAL(3.0);
		
		REQUIRE(a(0, 0) == REAL(3.0));
		REQUIRE(a(0, 1) == REAL(6.0));
		REQUIRE(a(1, 0) == REAL(9.0));
		REQUIRE(a(1, 1) == REAL(12.0));
	}

	TEST_CASE("Matrix::operator_div_equals_scalar", "[Matrix][compound]")
	{
		TEST_PRECISION_INFO();
		Matrix<Real> a(2, 2, { REAL(10.0), REAL(20.0), REAL(30.0), REAL(40.0) });
		
		a /= REAL(10.0);
		
		REQUIRE(a(0, 0) == REAL(1.0));
		REQUIRE(a(1, 1) == REAL(4.0));
	}

	TEST_CASE("Matrix::unary_negation", "[Matrix][arithmetic]")
	{
		TEST_PRECISION_INFO();
		Matrix<Real> a(2, 2, { REAL(1.0), REAL(-2.0), REAL(3.0), REAL(-4.0) });
		
		Matrix<Real> neg = -a;
		
		REQUIRE(neg(0, 0) == -REAL(1.0));
		REQUIRE(neg(0, 1) == REAL(2.0));
		REQUIRE(neg(1, 0) == -REAL(3.0));
		REQUIRE(neg(1, 1) == REAL(4.0));
	}

	///////////////////////          Trace and Invert                  //////////////////////
	TEST_CASE("Matrix::Trace", "[Matrix][trace]")
	{
		TEST_PRECISION_INFO();
		Matrix<Real> m(3, 3, { REAL(1.0), REAL(2.0), REAL(3.0),
		                       REAL(4.0), REAL(5.0), REAL(6.0),
		                       REAL(7.0), REAL(8.0), REAL(9.0) });
		
		REQUIRE_THAT(m.trace(), RealWithinRel(REAL(15.0)));  // 1 + 5 + 9
	}

	TEST_CASE("Matrix::Trace_non_square_throws", "[Matrix][trace]")
	{
		TEST_PRECISION_INFO();
		Matrix<Real> m(2, 3);
		
		REQUIRE_THROWS_AS(m.trace(), MatrixDimensionError);
	}

	TEST_CASE("Matrix::Invert_in_place", "[Matrix][invert]")
	{
		TEST_PRECISION_INFO();
		Matrix<Real> original(2, 2, { REAL(4.0), REAL(7.0), REAL(2.0), REAL(6.0) });
		Matrix<Real> m(original);
		
		m.Invert();
		
		// m * original = I
		Matrix<Real> product = m * original;
		REQUIRE(product.IsEqualTo(Matrix<Real>::Identity(2), TOL(1e-10, 1e-5)));
	}

	///////////////////////          Type Aliases                      //////////////////////
	TEST_CASE("Matrix::type_aliases", "[Matrix][aliases]")
	{
		TEST_PRECISION_INFO();
		
		// Test that aliases compile and work
		MatrixInt mi(2, 2, { 1, 2, 3, 4 });
		REQUIRE(mi(0, 0) == 1);
		
		MatrixFlt mf(2, 2, { 1.0f, 2.0f, 3.0f, 4.0f });
		REQUIRE(mf(1, 1) == 4.0f);
		
		MatrixDbl md(2, 2, { 1.0, 2.0, 3.0, 4.0 });
		REQUIRE(md(0, 1) == 2.0);
		
		// Short aliases
		MatI mi2(2, 2, { 5, 6, 7, 8 });
		REQUIRE(mi2(0, 0) == 5);
		
		MatF mf2(2, 2, { 5.0f, 6.0f, 7.0f, 8.0f });
		REQUIRE(mf2(1, 0) == 7.0f);
		
		MatD md2(2, 2, { 5.0, 6.0, 7.0, 8.0 });
		REQUIRE(md2(1, 1) == 8.0);
	}

	TEST_CASE("Matrix::complex_type", "[Matrix][complex]")
	{
		TEST_PRECISION_INFO();
		
		MatrixComplex mc(2, 2, { Complex(1.0, 2.0), Complex(3.0, 4.0),
		                         Complex(5.0, 6.0), Complex(7.0, 8.0) });
		
		REQUIRE(mc(0, 0).real() == 1.0);
		REQUIRE(mc(0, 0).imag() == 2.0);
		REQUIRE(mc(1, 1).real() == 7.0);
		
		// MatC short alias
		MatC mc2(2, 2, { Complex(1.0, 0.0), Complex(0.0, 1.0),
		                 Complex(0.0, -1.0), Complex(1.0, 0.0) });
		REQUIRE(mc2(0, 1).imag() == 1.0);
	}

	TEST_CASE("Matrix::complex_operations", "[Matrix][complex]")
	{
		TEST_PRECISION_INFO();
		
		MatrixComplex a(2, 2, { Complex(1.0, 1.0), Complex(2.0, 2.0),
		                        Complex(3.0, 3.0), Complex(4.0, 4.0) });
		MatrixComplex b(2, 2, { Complex(1.0, -1.0), Complex(2.0, -2.0),
		                        Complex(3.0, -3.0), Complex(4.0, -4.0) });
		
		// Addition
		MatrixComplex sum = a + b;
		REQUIRE(sum(0, 0).real() == 2.0);
		REQUIRE(sum(0, 0).imag() == 0.0);
		
		// Frobenius norm for complex
		MatrixComplex c(2, 2, { Complex(3.0, 4.0), Complex(0.0, 0.0),
		                        Complex(0.0, 0.0), Complex(0.0, 0.0) });
		// |3+4i|^2 = 25, so norm = 5
		REQUIRE_THAT(c.NormL2(), Catch::Matchers::WithinAbs(5.0, TOL(1e-10, 1e-5)));
	}

	///////////////////////          Resize with Preserve              //////////////////////
	TEST_CASE("Matrix::Resize_preserve_elements", "[Matrix][resize]")
	{
		TEST_PRECISION_INFO();
		Matrix<Real> m(2, 2, { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });
		
		m.Resize(3, 3, true);  // Preserve elements
		
		REQUIRE(m.rows() == 3);
		REQUIRE(m.cols() == 3);
		// Original elements preserved
		REQUIRE(m(0, 0) == REAL(1.0));
		REQUIRE(m(0, 1) == REAL(2.0));
		REQUIRE(m(1, 0) == REAL(3.0));
		REQUIRE(m(1, 1) == REAL(4.0));
		// New elements are zero
		REQUIRE(m(2, 2) == REAL(0.0));
		REQUIRE(m(0, 2) == REAL(0.0));
	}

	TEST_CASE("Matrix::Resize_shrink_preserve", "[Matrix][resize]")
	{
		TEST_PRECISION_INFO();
		Matrix<Real> m(3, 3, { REAL(1.0), REAL(2.0), REAL(3.0),
		                       REAL(4.0), REAL(5.0), REAL(6.0),
		                       REAL(7.0), REAL(8.0), REAL(9.0) });
		
		m.Resize(2, 2, true);
		
		REQUIRE(m.rows() == 2);
		REQUIRE(m.cols() == 2);
		REQUIRE(m(0, 0) == REAL(1.0));
		REQUIRE(m(0, 1) == REAL(2.0));
		REQUIRE(m(1, 0) == REAL(4.0));
		REQUIRE(m(1, 1) == REAL(5.0));
	}

} // namespace MML::Tests::Matrix



#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/BaseUtils.h"
#include "core/MatrixUtils.h"
#endif

#include "../test_data/linear_alg_eq_systems_test_bed.h"

using namespace MML;
using namespace MML::Testing;

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace MML::Tests::Core::MatrixUtilsTests
{
	/***************************************************************************/
	/*******                   DET algorithm tests                      ********/
	/***************************************************************************/
	TEST_CASE("Determinant of identity matrix is 1", "[MatrixDet]") {
			TEST_PRECISION_INFO();
		Matrix<Real> mat = Matrix<Real>::GetUnitMatrix(3);
		REQUIRE_THAT(Utils::Det(mat) , WithinRel(REAL(1.0), REAL(1e-5)));
	}

	TEST_CASE("Determinant of zero matrix is 0", "[MatrixDet]") {
			TEST_PRECISION_INFO();
		Matrix<Real> mat(3, 3, REAL(0.0));
		REQUIRE_THAT(Utils::Det(mat) , WithinRel(REAL(0.0), REAL(1e-5)));
	}

	TEST_CASE("Determinant of diagonal matrix is product of diagonal", "[MatrixDet]") {
			TEST_PRECISION_INFO();
		Matrix<Real> mat(3, 3, REAL(0.0));
		mat(0, 0) = REAL(2.0); mat(1, 1) = REAL(3.0); mat(2, 2) = REAL(4.0);
		REQUIRE_THAT(Utils::Det(mat) , WithinRel(REAL(24.0), REAL(1e-5))); // 2*3*4
	}

	TEST_CASE("Determinant of singular matrix is 0", "[MatrixDet]") {
			TEST_PRECISION_INFO();
		Matrix<Real> mat(2, 2);
		
		mat(0, 0) = REAL(1.0); mat(0, 1) = REAL(2.0);
		mat(1, 0) = REAL(2.0); mat(1, 1) = REAL(4.0); // Rows are linearly dependent
		
		REQUIRE_THAT(Utils::Det(mat) , WithinRel(REAL(0.0), REAL(1e-5)));
	}

	TEST_CASE("Determinant of known 3x3 matrix", "[MatrixDet]") {
			TEST_PRECISION_INFO();
		Matrix<Real> mat(3, 3);
		mat(0, 0) = REAL(6.0); mat(0, 1) = REAL(1.0); mat(0, 2) = REAL(1.0);
		mat(1, 0) = REAL(4.0); mat(1, 1) = -REAL(2.0); mat(1, 2) = REAL(5.0);
		mat(2, 0) = REAL(2.0); mat(2, 1) = REAL(8.0); mat(2, 2) = REAL(7.0);
		REQUIRE_THAT(Utils::Det(mat) , WithinRel(-REAL(306.0), REAL(1e-5)));
	}

	TEST_CASE("Determinant of known 4x4 matrix", "[MatrixDet]") {
			TEST_PRECISION_INFO();
		Matrix<Real> mat(4, 4);
		mat(0, 0) = REAL(3.0); mat(0, 1) = REAL(2.0); mat(0, 2) = REAL(0.0); mat(0, 3) = REAL(1.0);
		mat(1, 0) = REAL(4.0); mat(1, 1) = REAL(0.0); mat(1, 2) = REAL(1.0); mat(1, 3) = REAL(2.0);
		mat(2, 0) = REAL(3.0); mat(2, 1) = REAL(0.0); mat(2, 2) = REAL(2.0); mat(2, 3) = REAL(1.0);
		mat(3, 0) = REAL(9.0); mat(3, 1) = REAL(2.0); mat(3, 2) = REAL(3.0); mat(3, 3) = REAL(1.0);
		REQUIRE_THAT(Utils::Det(mat) , WithinRel(REAL(24.0), REAL(1e-5)));
	}

	// Complex matrix tests
	TEST_CASE("Determinant of zero complex matrix is 0", "[MatrixDet][Complex]") {
			TEST_PRECISION_INFO();
		Matrix<Complex> mat(3, 3, Complex(REAL(0.0), REAL(0.0)));
		REQUIRE(Utils::Det(mat) == Complex(REAL(0.0), REAL(0.0)));
	}

	TEST_CASE("Determinant of identity complex matrix is 1", "[MatrixDet][Complex]") {
			TEST_PRECISION_INFO();
		Matrix<Complex> mat(3, 3);
		for (int i = 0; i < 3; ++i)
			mat(i, i) = Complex(REAL(1.0), REAL(0.0));
		REQUIRE(Utils::Det(mat) == Complex(REAL(1.0), REAL(0.0)));
	}

	TEST_CASE("Determinant of diagonal complex matrix", "[MatrixDet][Complex]") {
			TEST_PRECISION_INFO();
		Matrix<Complex> mat(2, 2);
		mat(0, 0) = Complex(REAL(2.0), REAL(1.0));
		mat(1, 1) = Complex(REAL(3.0), -REAL(2.0));
		mat(0, 1) = mat(1, 0) = Complex(REAL(0.0), REAL(0.0));
		REQUIRE(Utils::Det(mat) == Complex(REAL(2.0), REAL(1.0)) * Complex(REAL(3.0), -REAL(2.0)));
	}

	TEST_CASE("Determinant of singular complex matrix is 0", "[MatrixDet][Complex]") {
			TEST_PRECISION_INFO();
		Matrix<Complex> mat(2, 2);
		mat(0, 0) = Complex(REAL(1.0), REAL(1.0)); mat(0, 1) = Complex(REAL(2.0), REAL(2.0));
		mat(1, 0) = Complex(REAL(2.0), REAL(2.0)); mat(1, 1) = Complex(REAL(4.0), REAL(4.0)); // 2x row 0
		REQUIRE(Utils::Det(mat) == Complex(REAL(0.0), REAL(0.0)));
	}

	TEST_CASE("Determinant of matrix with purely imaginary entries", "[MatrixDet][Complex]") {
			TEST_PRECISION_INFO();
		Matrix<Complex> mat(2, 2);
		mat(0, 0) = Complex(REAL(0.0), REAL(1.0)); mat(0, 1) = Complex(REAL(0.0), REAL(2.0));
		mat(1, 0) = Complex(REAL(0.0), REAL(3.0)); mat(1, 1) = Complex(REAL(0.0), REAL(4.0));
		// Determinant: (0+i1)*(0+i4) - (0+i2)*(0+i3) = (-1*4) - (-2*3) = -4 + 6 = 2 (real)
		REQUIRE_THAT(Utils::Det(mat).real() , WithinRel(REAL(2.0), REAL(1e-5)));
	}

	/***************************************************************************/
	/*******                   RANK algorithm tests                     ********/
	/***************************************************************************/
	TEST_CASE("Rank of all-zero matrix", "[MatrixRank][Edge]") {
			TEST_PRECISION_INFO();
		Matrix<Real> mat(4, 4, REAL(0.0));
		REQUIRE(Utils::Rank(mat) == 0);
	}

	TEST_CASE("Rank of identity matrix is full", "[MatrixRank]") {
			TEST_PRECISION_INFO();
		Matrix<Real> mat = Matrix<Real>::GetUnitMatrix(4);

		REQUIRE(Utils::Rank(mat) == 4);
	}

	TEST_CASE("Rank of matrix with linearly dependent rows", "[MatrixRank]") {
			TEST_PRECISION_INFO();
		Matrix<Real> mat(3, 3);
		mat(0, 0) = REAL(1.0); mat(0, 1) = REAL(2.0); mat(0, 2) = REAL(3.0);
		mat(1, 0) = REAL(2.0); mat(1, 1) = REAL(4.0); mat(1, 2) = REAL(6.0); // 2x row 0
		mat(2, 0) = REAL(3.0); mat(2, 1) = REAL(6.0); mat(2, 2) = REAL(9.0); // 3x row 0

		REQUIRE(Utils::Rank(mat) == 1);
	}

	TEST_CASE("Rank of rectangular matrix", "[MatrixRank]") {
			TEST_PRECISION_INFO();
		Matrix<Real> mat(2, 3);
		mat(0, 0) = REAL(1.0); mat(0, 1) = REAL(2.0); mat(0, 2) = REAL(3.0);
		mat(1, 0) = REAL(4.0); mat(1, 1) = REAL(5.0); mat(1, 2) = REAL(6.0);

		REQUIRE(Utils::Rank(mat) == 2);
	}

	TEST_CASE("Rank of matrix with one nonzero row", "[MatrixRank][Edge]") {
			TEST_PRECISION_INFO();
		Matrix<Real> mat(3, 3, REAL(0.0));
		mat(1, 0) = REAL(1.0); mat(1, 1) = REAL(2.0); mat(1, 2) = REAL(3.0);
		REQUIRE(Utils::Rank(mat) == 1);
	}

	TEST_CASE("Rank of matrix with one nonzero column", "[MatrixRank][Edge]") {
			TEST_PRECISION_INFO();
		Matrix<Real> mat(3, 3, REAL(0.0));
		mat(0, 2) = REAL(5.0); mat(1, 2) = REAL(7.0); mat(2, 2) = REAL(9.0);
		REQUIRE(Utils::Rank(mat) == 1);
	}

	TEST_CASE("Rank of matrix with nearly linearly dependent rows", "[MatrixRank][Edge]") {
			TEST_PRECISION_INFO();
		Matrix<Real> mat(2, 2);
		mat(0, 0) = REAL(1.0); mat(0, 1) = REAL(1.0);
		mat(1, 0) = REAL(1.0); mat(1, 1) = REAL(1.0) + 1e-10;		// Very close to row 0
		
		REQUIRE(Utils::Rank(mat, 1e-10) == 2); 

		// with bigger EPS, we consider them dependent
		REQUIRE(Utils::Rank(mat, 1e-9) == 1);
	}

	TEST_CASE("Rank of rectangular matrix with more columns than rows", "[MatrixRank][Edge]") {
			TEST_PRECISION_INFO();
		Matrix<Real> mat(2, 5, REAL(0.0));
		mat(0, 0) = REAL(1.0); mat(0, 1) = REAL(2.0); mat(0, 2) = REAL(3.0); mat(0, 3) = REAL(4.0); mat(0, 4) = REAL(5.0);
		mat(1, 0) = REAL(2.0); mat(1, 1) = REAL(4.0); mat(1, 2) = REAL(6.0); mat(1, 3) = REAL(8.0); mat(1, 4) = REAL(10.0); // 2x row 0
		REQUIRE(Utils::Rank(mat) == 1);
	}

	TEST_CASE("Rank of rectangular matrix with more rows than columns", "[MatrixRank][Edge]") {
			TEST_PRECISION_INFO();
		Matrix<Real> mat(5, 2, REAL(0.0));
		mat(0, 0) = REAL(1.0); mat(0, 1) = REAL(2.0);
		mat(1, 0) = REAL(3.0); mat(1, 1) = REAL(4.0);
		mat(2, 0) = REAL(5.0); mat(2, 1) = REAL(6.0);
		mat(3, 0) = REAL(7.0); mat(3, 1) = REAL(8.0);
		mat(4, 0) = REAL(9.0); mat(4, 1) = REAL(10.0);
		REQUIRE(Utils::Rank(mat) == 2);
	}

	TEST_CASE("Test Rank on TestBeds matrices")
	{
			TEST_PRECISION_INFO();
		Matrix<Real> A(3, 3, { 1,2,3,4,5,6,7,8,9 });

		// interesting
		REQUIRE(2 == Utils::Rank(A));

		REQUIRE(3 == Utils::Rank(TestBeds::mat_3x3));
		REQUIRE(3 == Utils::Rank(TestBeds::mat_3x3_1));
		REQUIRE(3 == Utils::Rank(TestBeds::mat_3x3_2));
		REQUIRE(3 == Utils::Rank(TestBeds::mat_3x3_3));
		REQUIRE(3 == Utils::Rank(TestBeds::mat_3x3_4));

		REQUIRE(5 == Utils::Rank(TestBeds::mat_5x5));
		REQUIRE(5 == Utils::Rank(TestBeds::mat_5x5_2));
		REQUIRE(5 == Utils::Rank(TestBeds::mat_5x5_3));
		REQUIRE(5 == Utils::Rank(TestBeds::mat_5x5_4));

		REQUIRE(8 == Utils::Rank(TestBeds::mat_8x8));

		REQUIRE(10 == Utils::Rank(TestBeds::mat_10x10));
		REQUIRE(10 == Utils::Rank(TestBeds::mat_10x10_1));
		REQUIRE(10 == Utils::Rank(TestBeds::mat_10x10_2));

		REQUIRE(20 == Utils::Rank(TestBeds::get_mat_20x20()));
		REQUIRE(20 == Utils::Rank(TestBeds::get_mat_20x20_1()));
		REQUIRE(20 == Utils::Rank(TestBeds::get_mat_20x20_2()));

		REQUIRE(50 == Utils::Rank(TestBeds::get_mat_50x50()));
		REQUIRE(50 == Utils::Rank(TestBeds::get_mat_50x50_1()));
		REQUIRE(50 == Utils::Rank(TestBeds::get_mat_50x50_2()));
	}

	TEST_CASE("Rank of zero complex matrix is 0", "[MatrixRank][Complex]") {
			TEST_PRECISION_INFO();
		Matrix<Complex> mat(3, 3);
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j)
				mat(i, j) = Complex(REAL(0.0), REAL(0.0));
		REQUIRE(Utils::Rank(mat) == 0);
	}

	TEST_CASE("Rank of identity complex matrix is full", "[MatrixRank][Complex]") {
			TEST_PRECISION_INFO();
		Matrix<Complex> mat(3, 3);
		for (int i = 0; i < 3; ++i)
			mat(i, i) = Complex(REAL(1.0), REAL(0.0));
		REQUIRE(Utils::Rank(mat) == 3);
	}

	TEST_CASE("Rank of complex matrix with linearly dependent rows", "[MatrixRank][Complex]") {
			TEST_PRECISION_INFO();
		Matrix<Complex> mat(3, 3);
		mat(0, 0) = Complex(REAL(1.0), REAL(1.0)); mat(0, 1) = Complex(REAL(2.0), REAL(2.0)); mat(0, 2) = Complex(REAL(3.0), REAL(3.0));
		mat(1, 0) = Complex(REAL(2.0), REAL(2.0)); mat(1, 1) = Complex(REAL(4.0), REAL(4.0)); mat(1, 2) = Complex(REAL(6.0), REAL(6.0)); // 2x row 0
		mat(2, 0) = Complex(REAL(3.0), REAL(3.0)); mat(2, 1) = Complex(REAL(6.0), REAL(6.0)); mat(2, 2) = Complex(REAL(9.0), REAL(9.0)); // 3x row 0
		REQUIRE(Utils::Rank(mat) == 1);
	}

	TEST_CASE("Rank of rectangular complex matrix", "[MatrixRank][Complex]") {
			TEST_PRECISION_INFO();
		Matrix<Complex> mat(2, 3);
		mat(0, 0) = Complex(REAL(1.0), REAL(0.0)); mat(0, 1) = Complex(REAL(2.0), REAL(1.0)); mat(0, 2) = Complex(REAL(3.0), -REAL(1.0));
		mat(1, 0) = Complex(REAL(4.0), REAL(2.0)); mat(1, 1) = Complex(REAL(5.0), REAL(0.0)); mat(1, 2) = Complex(REAL(6.0), REAL(1.0));
		REQUIRE(Utils::Rank(mat) == 2);
	}

	TEST_CASE("Rank of complex matrix with real and imaginary parts", "[MatrixRank][Complex]") {
			TEST_PRECISION_INFO();
		Matrix<Complex> mat(2, 2);
		mat(0, 0) = Complex(REAL(1.0), REAL(2.0)); mat(0, 1) = Complex(REAL(3.0), REAL(4.0));
		mat(1, 0) = Complex(REAL(5.0), REAL(6.0)); mat(1, 1) = Complex(REAL(7.0), REAL(8.0));
		REQUIRE(Utils::Rank(mat) == 2);
	}

	TEST_CASE("Test Utils::IsOrthogonal")
	{
			TEST_PRECISION_INFO();
		REQUIRE(Utils::IsOrthogonal(Matrix<Real>(3, 3, { 1,0,0, 0,1,0, 0,0,1 })));

		Matrix<Real> B(3, 3, { 1,0,0,0,1,0,0,0,2 });

		REQUIRE_FALSE(Utils::IsOrthogonal(B));

		Matrix<Real> C(3, 3, { 1,0,0,0,1,0,0,0,0 });

		REQUIRE_FALSE(Utils::IsOrthogonal(C));

		Matrix<Real> D(3, 3, { 1,0,0,0,1,0,0,0,0 });

		REQUIRE_FALSE(Utils::IsOrthogonal(D));

		Matrix<Real> E(3, 3, { 1,0,0,0,1,0,0,0,0 });

		REQUIRE_FALSE(Utils::IsOrthogonal(E));

		Matrix<Real> F(3, 3, { 1,0,0,0,1,0,0,0,0 });

		REQUIRE_FALSE(Utils::IsOrthogonal(F));

		Matrix<Real> G(3, 3, { 1,0,0,0,1,0,0,0,0 });

		REQUIRE_FALSE(Utils::IsOrthogonal(G));

		Matrix<Real> H(3, 3, { 1,0,0,0,1,0,0,0,0 });

		REQUIRE_FALSE(Utils::IsOrthogonal(H));

		Matrix<Real> I(3, 3, { 1,0,0,0,1,0,0,0,0 });

		REQUIRE_FALSE(Utils::IsOrthogonal(I));

		Matrix<Real> J(3, 3, { 1,0,0,0,1,0,0,0,0 });

		REQUIRE_FALSE(Utils::IsOrthogonal(J));

		Matrix<Real> K(3, 3, { 1,0,0,0,1,0,0,0,0 });

		REQUIRE_FALSE(Utils::IsOrthogonal(K));

		Matrix<Real> L(3, 3, { 1,0,0,0,1,0,0,0,0 });

		REQUIRE_FALSE(Utils::IsOrthogonal(L));
	}

	TEST_CASE("Test Utils::Definitness")
	{

	}

	TEST_CASE("Test Utils::IsHermitian")
	{

	}

	TEST_CASE("Test Utils::IsUnitary")
	{

	}
}
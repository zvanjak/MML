#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/BaseUtils.h"
#include "core/MatrixUtils.h"
#endif

#include "../test_beds/linear_alg_eq_systems_test_bed.h"

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
		Matrix<Real> mat = Matrix<Real>::Identity(3);
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
		Matrix<Real> mat = Matrix<Real>::Identity(4);

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

		REQUIRE(3 == Utils::Rank(TestBeds::mat_3x3()));
		REQUIRE(3 == Utils::Rank(TestBeds::mat_3x3_1()));
		REQUIRE(3 == Utils::Rank(TestBeds::mat_3x3_2()));
		REQUIRE(3 == Utils::Rank(TestBeds::mat_3x3_3()));
		REQUIRE(3 == Utils::Rank(TestBeds::mat_3x3_4()));

		REQUIRE(5 == Utils::Rank(TestBeds::mat_5x5()));
		REQUIRE(5 == Utils::Rank(TestBeds::mat_5x5_2()));
		REQUIRE(5 == Utils::Rank(TestBeds::mat_5x5_3()));
		REQUIRE(5 == Utils::Rank(TestBeds::mat_5x5_4()));

		REQUIRE(8 == Utils::Rank(TestBeds::mat_8x8()));

		REQUIRE(10 == Utils::Rank(TestBeds::mat_10x10()));
		REQUIRE(10 == Utils::Rank(TestBeds::mat_10x10_1()));
		REQUIRE(10 == Utils::Rank(TestBeds::mat_10x10_2()));

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

	/***************************************************************************/
	/*******              MATRIX PROPERTY CHECKS                        ********/
	/***************************************************************************/
	
	TEST_CASE("IsUpperTriangular tests", "[MatrixUtils][Properties]")
	{
		TEST_PRECISION_INFO();
		
		SECTION("Identity is upper triangular")
		{
			auto I = Matrix<Real>::Identity(3);
			REQUIRE(Utils::IsUpperTriangular(I));
		}
		
		SECTION("True upper triangular")
		{
			Matrix<Real> U(3, 3, { 1,2,3, 0,4,5, 0,0,6 });
			REQUIRE(Utils::IsUpperTriangular(U));
		}
		
		SECTION("Not upper triangular")
		{
			Matrix<Real> A(3, 3, { 1,2,3, 1,4,5, 0,0,6 });
			REQUIRE_FALSE(Utils::IsUpperTriangular(A));
		}
		
		SECTION("Diagonal is upper triangular")
		{
			Matrix<Real> D(3, 3, { 1,0,0, 0,2,0, 0,0,3 });
			REQUIRE(Utils::IsUpperTriangular(D));
		}
	}
	
	TEST_CASE("IsLowerTriangular tests", "[MatrixUtils][Properties]")
	{
		TEST_PRECISION_INFO();
		
		SECTION("Identity is lower triangular")
		{
			auto I = Matrix<Real>::Identity(3);
			REQUIRE(Utils::IsLowerTriangular(I));
		}
		
		SECTION("True lower triangular")
		{
			Matrix<Real> L(3, 3, { 1,0,0, 2,4,0, 3,5,6 });
			REQUIRE(Utils::IsLowerTriangular(L));
		}
		
		SECTION("Not lower triangular")
		{
			Matrix<Real> A(3, 3, { 1,1,0, 2,4,0, 3,5,6 });
			REQUIRE_FALSE(Utils::IsLowerTriangular(A));
		}
	}
	
	TEST_CASE("IsDiagonal tests", "[MatrixUtils][Properties]")
	{
		TEST_PRECISION_INFO();
		
		SECTION("Identity is diagonal")
		{
			auto I = Matrix<Real>::Identity(3);
			REQUIRE(Utils::IsDiagonal(I));
		}
		
		SECTION("True diagonal")
		{
			Matrix<Real> D(3, 3, { 1,0,0, 0,2,0, 0,0,3 });
			REQUIRE(Utils::IsDiagonal(D));
		}
		
		SECTION("Zero matrix is diagonal")
		{
			Matrix<Real> Z(3, 3, REAL(0.0));
			REQUIRE(Utils::IsDiagonal(Z));
		}
		
		SECTION("Not diagonal")
		{
			Matrix<Real> A(3, 3, { 1,1,0, 0,2,0, 0,0,3 });
			REQUIRE_FALSE(Utils::IsDiagonal(A));
		}
	}
	
	TEST_CASE("IsSymmetric tests", "[MatrixUtils][Properties]")
	{
		TEST_PRECISION_INFO();
		
		SECTION("Identity is symmetric")
		{
			auto I = Matrix<Real>::Identity(3);
			REQUIRE(Utils::IsSymmetric(I));
		}
		
		SECTION("Symmetric matrix")
		{
			Matrix<Real> S(3, 3, { 1,2,3, 2,4,5, 3,5,6 });
			REQUIRE(Utils::IsSymmetric(S));
		}
		
		SECTION("Not symmetric")
		{
			Matrix<Real> A(3, 3, { 1,2,3, 4,5,6, 7,8,9 });
			REQUIRE_FALSE(Utils::IsSymmetric(A));
		}
		
		SECTION("Non-square is not symmetric")
		{
			Matrix<Real> R(2, 3, { 1,2,3, 4,5,6 });
			REQUIRE_FALSE(Utils::IsSymmetric(R));
		}
	}
	
	TEST_CASE("IsSkewSymmetric tests", "[MatrixUtils][Properties]")
	{
		TEST_PRECISION_INFO();
		
		SECTION("True skew-symmetric")
		{
			Matrix<Real> S(3, 3, { 0,2,-3, -2,0,5, 3,-5,0 });
			REQUIRE(Utils::IsSkewSymmetric(S));
		}
		
		SECTION("Not skew-symmetric - nonzero diagonal")
		{
			Matrix<Real> A(3, 3, { 1,2,-3, -2,0,5, 3,-5,0 });
			REQUIRE_FALSE(Utils::IsSkewSymmetric(A));
		}
		
		SECTION("Zero matrix is skew-symmetric")
		{
			Matrix<Real> Z(3, 3, REAL(0.0));
			REQUIRE(Utils::IsSkewSymmetric(Z));
		}
	}
	
	TEST_CASE("IsUpperHessenberg tests", "[MatrixUtils][Properties]")
	{
		TEST_PRECISION_INFO();
		
		SECTION("Upper triangular is upper Hessenberg")
		{
			Matrix<Real> U(4, 4, { 1,2,3,4, 0,5,6,7, 0,0,8,9, 0,0,0,10 });
			REQUIRE(Utils::IsUpperHessenberg(U));
		}
		
		SECTION("True upper Hessenberg")
		{
			Matrix<Real> H(4, 4, { 1,2,3,4, 5,6,7,8, 0,9,10,11, 0,0,12,13 });
			REQUIRE(Utils::IsUpperHessenberg(H));
		}
		
		SECTION("Not upper Hessenberg")
		{
			Matrix<Real> A(4, 4, { 1,2,3,4, 5,6,7,8, 1,9,10,11, 0,0,12,13 });  // A(2,0)=1 violates
			REQUIRE_FALSE(Utils::IsUpperHessenberg(A));
		}
	}
	
	/***************************************************************************/
	/*******              MATRIX NORMS AND SCALAR QUANTITIES            ********/
	/***************************************************************************/
	
	TEST_CASE("Trace tests", "[MatrixUtils][Norms]")
	{
		TEST_PRECISION_INFO();
		
		SECTION("Identity trace equals dimension")
		{
			auto I = Matrix<Real>::Identity(4);
			REQUIRE_THAT(Utils::Trace(I), WithinAbs(REAL(4.0), REAL(1e-14)));
		}
		
		SECTION("Diagonal matrix trace is sum of diagonal")
		{
			Matrix<Real> D(3, 3, { 1,0,0, 0,2,0, 0,0,3 });
			REQUIRE_THAT(Utils::Trace(D), WithinAbs(REAL(6.0), REAL(1e-14)));
		}
		
		SECTION("General matrix trace")
		{
			Matrix<Real> A(3, 3, { 1,2,3, 4,5,6, 7,8,9 });
			REQUIRE_THAT(Utils::Trace(A), WithinAbs(REAL(15.0), REAL(1e-14)));  // 1+5+9
		}
		
		SECTION("Zero matrix trace is zero")
		{
			Matrix<Real> Z(3, 3, REAL(0.0));
			REQUIRE_THAT(Utils::Trace(Z), WithinAbs(REAL(0.0), REAL(1e-14)));
		}
	}
	
	TEST_CASE("FrobeniusNorm tests", "[MatrixUtils][Norms]")
	{
		TEST_PRECISION_INFO();
		
		SECTION("Identity Frobenius norm")
		{
			auto I = Matrix<Real>::Identity(3);
			REQUIRE_THAT(Utils::FrobeniusNorm(I), WithinAbs(std::sqrt(REAL(3.0)), REAL(1e-14)));
		}
		
		SECTION("Simple matrix")
		{
			Matrix<Real> A(2, 2, { 1,2, 3,4 });
			// ||A||_F = sqrt(1+4+9+16) = sqrt(30)
			REQUIRE_THAT(Utils::FrobeniusNorm(A), WithinAbs(std::sqrt(REAL(30.0)), REAL(1e-14)));
		}
		
		SECTION("Zero matrix")
		{
			Matrix<Real> Z(3, 3, REAL(0.0));
			REQUIRE_THAT(Utils::FrobeniusNorm(Z), WithinAbs(REAL(0.0), REAL(1e-14)));
		}
	}
	
	TEST_CASE("InfinityNorm tests (max row sum)", "[MatrixUtils][Norms]")
	{
		TEST_PRECISION_INFO();
		
		SECTION("Identity infinity norm")
		{
			auto I = Matrix<Real>::Identity(3);
			REQUIRE_THAT(Utils::InfinityNorm(I), WithinAbs(REAL(1.0), REAL(1e-14)));
		}
		
		SECTION("Simple matrix")
		{
			Matrix<Real> A(2, 3, { 1,2,3, 4,5,6 });
			// Row sums: 6, 15 -> max = 15
			REQUIRE_THAT(Utils::InfinityNorm(A), WithinAbs(REAL(15.0), REAL(1e-14)));
		}
		
		SECTION("Matrix with negative entries")
		{
			Matrix<Real> A(2, 2, { -1,2, 3,-4 });
			// Row sums: |−1|+|2|=3, |3|+|−4|=7 -> max = 7
			REQUIRE_THAT(Utils::InfinityNorm(A), WithinAbs(REAL(7.0), REAL(1e-14)));
		}
	}
	
	TEST_CASE("OneNorm tests (max column sum)", "[MatrixUtils][Norms]")
	{
		TEST_PRECISION_INFO();
		
		SECTION("Identity one norm")
		{
			auto I = Matrix<Real>::Identity(3);
			REQUIRE_THAT(Utils::OneNorm(I), WithinAbs(REAL(1.0), REAL(1e-14)));
		}
		
		SECTION("Simple matrix")
		{
			Matrix<Real> A(3, 2, { 1,2, 3,4, 5,6 });
			// Column sums: 1+3+5=9, 2+4+6=12 -> max = 12
			REQUIRE_THAT(Utils::OneNorm(A), WithinAbs(REAL(12.0), REAL(1e-14)));
		}
	}
	
	TEST_CASE("MaxAbsDiff tests", "[MatrixUtils][Norms]")
	{
		TEST_PRECISION_INFO();
		
		SECTION("Same matrices")
		{
			Matrix<Real> A(2, 2, { 1,2, 3,4 });
			REQUIRE_THAT(Utils::MaxAbsDiff(A, A), WithinAbs(REAL(0.0), REAL(1e-14)));
		}
		
		SECTION("Different matrices")
		{
			Matrix<Real> A(2, 2, { 1,2, 3,4 });
			Matrix<Real> B(2, 2, { 1,2, 3,10 });
			REQUIRE_THAT(Utils::MaxAbsDiff(A, B), WithinAbs(REAL(6.0), REAL(1e-14)));
		}
	}
	
	/***************************************************************************/
	/*******              SIMILARITY TRANSFORM                          ********/
	/***************************************************************************/
	
	TEST_CASE("SimilarityTransform tests", "[MatrixUtils][Transform]")
	{
		TEST_PRECISION_INFO();
		
		SECTION("Identity transform preserves matrix")
		{
			Matrix<Real> A(3, 3, { 1,2,3, 4,5,6, 7,8,9 });
			auto I = Matrix<Real>::Identity(3);
			auto result = Utils::SimilarityTransform(I, A);
			REQUIRE_THAT(Utils::MaxAbsDiff(result, A), WithinAbs(REAL(0.0), REAL(1e-12)));
		}
		
		SECTION("Permutation matrix")
		{
			Matrix<Real> A(2, 2, { 1,2, 3,4 });
			Matrix<Real> P(2, 2, { 0,1, 1,0 });  // Swap rows/cols
			auto result = Utils::SimilarityTransform(P, A);
			// P^T * A * P swaps both rows and columns
			Matrix<Real> expected(2, 2, { 4,3, 2,1 });
			REQUIRE_THAT(Utils::MaxAbsDiff(result, expected), WithinAbs(REAL(0.0), REAL(1e-12)));
		}
	}
	
	/***************************************************************************/
	/*******              NILPOTENT AND UNIPOTENT                       ********/
	/***************************************************************************/
	
	TEST_CASE("IsNilpotent tests", "[MatrixUtils][Properties]")
	{
		TEST_PRECISION_INFO();
		
		SECTION("Zero matrix is nilpotent")
		{
			Matrix<Real> Z(3, 3, REAL(0.0));
			REQUIRE(Utils::IsNilpotent(Z));
		}
		
		SECTION("Strictly upper triangular is nilpotent")
		{
			Matrix<Real> N(3, 3, { 0,1,2, 0,0,3, 0,0,0 });
			REQUIRE(Utils::IsNilpotent(N));
		}
		
		SECTION("Identity is not nilpotent")
		{
			auto I = Matrix<Real>::Identity(3);
			REQUIRE_FALSE(Utils::IsNilpotent(I));
		}
		
		SECTION("General non-nilpotent matrix")
		{
			Matrix<Real> A(2, 2, { 1,1, 0,1 });
			REQUIRE_FALSE(Utils::IsNilpotent(A));
		}
	}
	
	TEST_CASE("IsUnipotent tests", "[MatrixUtils][Properties]")
	{
		TEST_PRECISION_INFO();
		
		SECTION("Identity is unipotent")
		{
			auto I = Matrix<Real>::Identity(3);
			REQUIRE(Utils::IsUnipotent(I));
		}
		
		SECTION("I + nilpotent is unipotent")
		{
			// U = I + N where N is strictly upper triangular
			Matrix<Real> U(3, 3, { 1,1,2, 0,1,3, 0,0,1 });
			REQUIRE(Utils::IsUnipotent(U));
		}
		
		SECTION("General matrix is not unipotent")
		{
			Matrix<Real> A(2, 2, { 2,1, 0,2 });
			REQUIRE_FALSE(Utils::IsUnipotent(A));
		}
	}
	
	/***************************************************************************/
	/*******              FADDEEV-LEVERRIER ALGORITHM                   ********/
	/***************************************************************************/
	
	TEST_CASE("FaddeevAlg tests", "[MatrixUtils][Faddeev]")
	{
		TEST_PRECISION_INFO();
		
		SECTION("2x2 matrix")
		{
			Matrix<Real> A(2, 2, { 1,2, 3,4 });
			PolynomRealFunc charPoly;
			Real det;
			Matrix<Real> inv(2, 2);
			
			Utils::FaddeevAlg(A, charPoly, det, inv);
			
			// Characteristic polynomial: λ² - tr(A)λ + det(A) = λ² - 5λ - 2
			REQUIRE_THAT(charPoly[2], WithinAbs(REAL(1.0), REAL(1e-10)));
			REQUIRE_THAT(charPoly[1], WithinAbs(-REAL(5.0), REAL(1e-10)));  // -trace
			REQUIRE_THAT(charPoly[0], WithinAbs(-REAL(2.0), REAL(1e-10)));  // det
			
			// Determinant: 1*4 - 2*3 = -2
			REQUIRE_THAT(det, WithinAbs(-REAL(2.0), REAL(1e-10)));
			
			// Verify inverse: A * inv = I
			auto product = A * inv;
			auto I = Matrix<Real>::Identity(2);
			REQUIRE_THAT(Utils::MaxAbsDiff(product, I), WithinAbs(REAL(0.0), REAL(1e-10)));
		}
		
		SECTION("3x3 identity matrix")
		{
			auto A = Matrix<Real>::Identity(3);
			PolynomRealFunc charPoly;
			Real det;
			Matrix<Real> inv(3, 3);
			
			Utils::FaddeevAlg(A, charPoly, det, inv);
			
			// Characteristic polynomial: (λ-1)³ = λ³ - 3λ² + 3λ - 1
			REQUIRE_THAT(charPoly[3], WithinAbs(REAL(1.0), REAL(1e-10)));
			REQUIRE_THAT(charPoly[2], WithinAbs(-REAL(3.0), REAL(1e-10)));
			REQUIRE_THAT(charPoly[1], WithinAbs(REAL(3.0), REAL(1e-10)));
			REQUIRE_THAT(charPoly[0], WithinAbs(-REAL(1.0), REAL(1e-10)));
			
			// det(I) = 1
			REQUIRE_THAT(det, WithinAbs(REAL(1.0), REAL(1e-10)));
			
			// inv(I) = I
			REQUIRE_THAT(Utils::MaxAbsDiff(inv, A), WithinAbs(REAL(0.0), REAL(1e-10)));
		}
	}
	
	/***************************************************************************/
	/*******              SVD-BASED UTILITIES                           ********/
	/***************************************************************************/
	
	TEST_CASE("ComputeSVD tests", "[MatrixUtils][SVD]")
	{
		TEST_PRECISION_INFO();
		
		SECTION("Identity matrix SVD")
		{
			auto I = Matrix<Real>::Identity(3);
			auto result = Utils::ComputeSVD(I);
			
			REQUIRE(result.rank == 3);
			REQUIRE_THAT(result.conditionNumber, WithinAbs(REAL(1.0), REAL(1e-10)));
			
			// All singular values should be 1
			for (int i = 0; i < 3; i++)
				REQUIRE_THAT(result.singularValues[i], WithinAbs(REAL(1.0), REAL(1e-10)));
		}
		
		SECTION("Rank-deficient matrix")
		{
			Matrix<Real> A(3, 3, { 1,2,3, 2,4,6, 3,6,9 });  // Rank 1
			auto result = Utils::ComputeSVD(A);
			
			REQUIRE(result.rank == 1);
		}
	}
	
	TEST_CASE("SingularValues tests", "[MatrixUtils][SVD]")
	{
		TEST_PRECISION_INFO();
		
		SECTION("Diagonal matrix")
		{
			Matrix<Real> D(3, 3, { 3,0,0, 0,2,0, 0,0,1 });
			auto sv = Utils::SingularValues(D);
			
			// Singular values should be 3, 2, 1 in descending order
			REQUIRE_THAT(sv[0], WithinAbs(REAL(3.0), REAL(1e-10)));
			REQUIRE_THAT(sv[1], WithinAbs(REAL(2.0), REAL(1e-10)));
			REQUIRE_THAT(sv[2], WithinAbs(REAL(1.0), REAL(1e-10)));
		}
	}
	
	TEST_CASE("RankSVD tests", "[MatrixUtils][SVD]")
	{
		TEST_PRECISION_INFO();
		
		SECTION("Full rank matrix")
		{
			Matrix<Real> A(3, 3, { 1,2,3, 4,5,6, 7,8,10 });  // Full rank
			REQUIRE(Utils::RankSVD(A) == 3);
		}
		
		SECTION("Rank 2 matrix")
		{
			Matrix<Real> A(3, 3, { 1,2,3, 4,5,6, 5,7,9 });  // row3 = row1 + row2
			REQUIRE(Utils::RankSVD(A) == 2);
		}
		
		SECTION("Consistency with Gaussian rank")
		{
			Matrix<Real> A(3, 3, { 1,2,3, 4,5,6, 7,8,9 });  // Rank 2
			REQUIRE(Utils::RankSVD(A) == Utils::RankGaussian(A));
		}
	}
	
	TEST_CASE("Nullity tests", "[MatrixUtils][SVD]")
	{
		TEST_PRECISION_INFO();
		
		SECTION("Full rank - zero nullity")
		{
			auto I = Matrix<Real>::Identity(3);
			REQUIRE(Utils::Nullity(I) == 0);
		}
		
		SECTION("Rank 2 in 3x3 - nullity 1")
		{
			Matrix<Real> A(3, 3, { 1,2,3, 4,5,6, 5,7,9 });
			REQUIRE(Utils::Nullity(A) == 1);
		}
		
		SECTION("Rank + Nullity = n")
		{
			Matrix<Real> A(3, 3, { 1,2,3, 2,4,6, 3,6,9 });  // Rank 1
			REQUIRE(Utils::RankSVD(A) + Utils::Nullity(A) == 3);
		}
	}
	
	TEST_CASE("ConditionNumber tests", "[MatrixUtils][SVD]")
	{
		TEST_PRECISION_INFO();
		
		SECTION("Identity has condition number 1")
		{
			auto I = Matrix<Real>::Identity(3);
			REQUIRE_THAT(Utils::ConditionNumber(I), WithinAbs(REAL(1.0), REAL(1e-10)));
		}
		
		SECTION("Well-conditioned matrix")
		{
			Matrix<Real> A(2, 2, { 2,0, 0,1 });
			REQUIRE_THAT(Utils::ConditionNumber(A), WithinAbs(REAL(2.0), REAL(1e-10)));
		}
		
		SECTION("Ill-conditioned matrix has large condition number")
		{
			Matrix<Real> A(2, 2, { 1,0, 0,1e-10 });
			REQUIRE(Utils::ConditionNumber(A) > 1e9);
		}
	}
	
	TEST_CASE("NullSpace tests", "[MatrixUtils][SVD]")
	{
		TEST_PRECISION_INFO();
		
		SECTION("Full rank has empty null space")
		{
			auto I = Matrix<Real>::Identity(3);
			auto ns = Utils::NullSpace(I);
			REQUIRE(ns.cols() == 0);  // Empty
		}
		
		SECTION("Rank-deficient matrix has non-trivial null space")
		{
			Matrix<Real> A(3, 3, { 1,2,3, 2,4,6, 3,6,9 });  // Rank 1, nullity 2
			auto ns = Utils::NullSpace(A);
			REQUIRE(ns.cols() == 2);  // Two basis vectors
			
			// Verify A * ns ≈ 0
			auto product = A * ns;
			REQUIRE_THAT(Utils::FrobeniusNorm(product), WithinAbs(REAL(0.0), REAL(1e-10)));
		}
	}
	
	TEST_CASE("ColumnSpace tests", "[MatrixUtils][SVD]")
	{
		TEST_PRECISION_INFO();
		
		SECTION("Full rank has full column space")
		{
			auto I = Matrix<Real>::Identity(3);
			auto cs = Utils::ColumnSpace(I);
			REQUIRE(cs.cols() == 3);
		}
		
		SECTION("Rank-deficient matrix")
		{
			Matrix<Real> A(3, 3, { 1,2,3, 2,4,6, 3,6,9 });  // Rank 1
			auto cs = Utils::ColumnSpace(A);
			REQUIRE(cs.cols() == 1);
		}
	}
	
	TEST_CASE("FundamentalSubspaces tests", "[MatrixUtils][SVD]")
	{
		TEST_PRECISION_INFO();
		
		SECTION("Dimension relationships")
		{
			Matrix<Real> A(3, 4, { 1,2,3,4, 2,4,6,8, 3,6,9,12 });  // 3x4, rank 1
			auto fs = Utils::ComputeFundamentalSubspaces(A);
			
			REQUIRE(fs.rank == 1);
			REQUIRE(fs.columnSpace.cols() == 1);      // dim = rank
			REQUIRE(fs.rowSpace.cols() == 1);         // dim = rank
			REQUIRE(fs.nullSpace.cols() == 3);        // dim = n - rank = 4 - 1
			REQUIRE(fs.leftNullSpace.cols() == 2);    // dim = m - rank = 3 - 1
		}
	}
	
	TEST_CASE("PseudoInverse tests", "[MatrixUtils][SVD]")
	{
		TEST_PRECISION_INFO();
		
		SECTION("Pseudoinverse of identity is identity")
		{
			auto I = Matrix<Real>::Identity(3);
			auto pinv = Utils::PseudoInverse(I);
			REQUIRE_THAT(Utils::MaxAbsDiff(pinv, I), WithinAbs(REAL(0.0), REAL(1e-10)));
		}
		
		SECTION("Pseudoinverse of full-rank square equals inverse")
		{
			Matrix<Real> A(2, 2, { 1,2, 3,4 });
			auto pinv = Utils::PseudoInverse(A);
			
			// A * pinv should be identity
			auto product = A * pinv;
			auto I = Matrix<Real>::Identity(2);
			REQUIRE_THAT(Utils::MaxAbsDiff(product, I), WithinAbs(REAL(0.0), REAL(1e-10)));
		}
		
		SECTION("Moore-Penrose property: A * A+ * A = A")
		{
			Matrix<Real> A(3, 2, { 1,2, 3,4, 5,6 });  // 3x2 rectangular
			auto pinv = Utils::PseudoInverse(A);
			
			auto AAA = A * pinv * A;
			REQUIRE_THAT(Utils::MaxAbsDiff(AAA, A), WithinAbs(REAL(0.0), REAL(1e-10)));
		}
		
		SECTION("Moore-Penrose property: A+ * A * A+ = A+")
		{
			Matrix<Real> A(2, 3, { 1,2,3, 4,5,6 });  // 2x3 rectangular
			auto pinv = Utils::PseudoInverse(A);
			
			auto pApAp = pinv * A * pinv;
			REQUIRE_THAT(Utils::MaxAbsDiff(pApAp, pinv), WithinAbs(REAL(0.0), REAL(1e-10)));
		}
	}
}


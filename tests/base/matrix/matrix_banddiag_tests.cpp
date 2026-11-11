#include <catch2/catch_all.hpp>
#include "../../TestPrecision.h"
#include "../../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/base/Vector/Vector.h>
#include <mml/base/Matrix/MatrixBandDiag.h>
#include <mml/base/Matrix/Matrix.h>
#endif

using namespace MML;
using namespace MML::Testing;

namespace MML::Tests::Base::BandDiagMatrixTests
{
	TEST_CASE("BandDiagonalMatrix_Basic_Construction", "[BandDiag]")
	{
			TEST_PRECISION_INFO();
		// Create a 5x5 band diagonal matrix with m1=1 (1 subdiagonal), m2=2 (2 superdiagonals)
		// Storage: [n][m1+m2+1] = [5][4]
		Matrix<Real> data(5, 4, {
			// Each row: [sub1, diag, super1, super2]
			static_cast<MML::Real>(0.0),  static_cast<MML::Real>(1.0),  static_cast<MML::Real>(2.0),  static_cast<MML::Real>(3.0),   // Row 0: diag=1, super1=2, super2=3
			static_cast<MML::Real>(4.0),  static_cast<MML::Real>(5.0),  static_cast<MML::Real>(6.0),  static_cast<MML::Real>(7.0),   // Row 1: sub1=4, diag=5, super1=6, super2=7
			static_cast<MML::Real>(8.0),  static_cast<MML::Real>(9.0), static_cast<MML::Real>(10.0), static_cast<MML::Real>(11.0),   // Row 2: sub1=8, diag=9, super1=10, super2=11
			static_cast<MML::Real>(12.0), static_cast<MML::Real>(13.0), static_cast<MML::Real>(14.0),  static_cast<MML::Real>(0.0),  // Row 3: sub1=12, diag=13, super1=14
			static_cast<MML::Real>(16.0), static_cast<MML::Real>(17.0),  static_cast<MML::Real>(0.0),  static_cast<MML::Real>(0.0)   // Row 4: sub1=16, diag=17
		});

		BandDiagonalMatrix mat(5, 1, 2, data);

		REQUIRE(mat.rows() == 5);
		REQUIRE(mat.cols() == 5);
		REQUIRE(mat.getLowerBandwidth() == 1);
		REQUIRE(mat.getUpperBandwidth() == 2);
		REQUIRE(mat.rows() == 5);
	}

	TEST_CASE("BandDiagonalMatrix_Element_Access", "[BandDiag]")
	{
			TEST_PRECISION_INFO();
		// Create a simple 5x5 band matrix with m1=1, m2=1 (tridiagonal)
		Matrix<Real> data(5, 3, {
			// [sub1, diag, super1]
			static_cast<MML::Real>(0.0),  static_cast<MML::Real>(2.0),  static_cast<MML::Real>(1.0),
			static_cast<MML::Real>(1.0),  static_cast<MML::Real>(2.0),  static_cast<MML::Real>(1.0),
			static_cast<MML::Real>(1.0),  static_cast<MML::Real>(2.0),  static_cast<MML::Real>(1.0),
			static_cast<MML::Real>(1.0),  static_cast<MML::Real>(2.0),  static_cast<MML::Real>(1.0),
			static_cast<MML::Real>(1.0),  static_cast<MML::Real>(2.0),  static_cast<MML::Real>(0.0)
		});

		BandDiagonalMatrix mat(5, 1, 1, data);

		// Test diagonal elements
		REQUIRE(mat(0, 0) == static_cast<MML::Real>(2.0));
		REQUIRE(mat(1, 1) == static_cast<MML::Real>(2.0));
		REQUIRE(mat(2, 2) == static_cast<MML::Real>(2.0));
		REQUIRE(mat(3, 3) == static_cast<MML::Real>(2.0));
		REQUIRE(mat(4, 4) == static_cast<MML::Real>(2.0));

		// Test superdiagonal
		REQUIRE(mat(0, 1) == static_cast<MML::Real>(1.0));
		REQUIRE(mat(1, 2) == static_cast<MML::Real>(1.0));
		REQUIRE(mat(2, 3) == static_cast<MML::Real>(1.0));
		REQUIRE(mat(3, 4) == static_cast<MML::Real>(1.0));

		// Test subdiagonal
		REQUIRE(mat(1, 0) == static_cast<MML::Real>(1.0));
		REQUIRE(mat(2, 1) == static_cast<MML::Real>(1.0));
		REQUIRE(mat(3, 2) == static_cast<MML::Real>(1.0));
		REQUIRE(mat(4, 3) == static_cast<MML::Real>(1.0));

		// Note: accessing elements outside band throws exception
		// (following TridiagonalMatrix pattern)
		// REQUIRE(mat(0, 2) == static_cast<MML::Real>(0.0));   // would throw
		// REQUIRE(mat(2, 0) == static_cast<MML::Real>(0.0));   // would throw
	}

	TEST_CASE("BandDiagonalMatrix_IsInBand", "[BandDiag]")
	{
			TEST_PRECISION_INFO();
		Matrix<Real> data(5, 3, {
			static_cast<MML::Real>(0.0),  static_cast<MML::Real>(1.0),  static_cast<MML::Real>(1.0),
			static_cast<MML::Real>(1.0),  static_cast<MML::Real>(1.0),  static_cast<MML::Real>(1.0),
			static_cast<MML::Real>(1.0),  static_cast<MML::Real>(1.0),  static_cast<MML::Real>(1.0),
			static_cast<MML::Real>(1.0),  static_cast<MML::Real>(1.0),  static_cast<MML::Real>(1.0),
			static_cast<MML::Real>(1.0),  static_cast<MML::Real>(1.0),  static_cast<MML::Real>(0.0)
		});

		BandDiagonalMatrix mat(5, 1, 1, data);

		// Elements on the band
		REQUIRE(mat.IsInBand(0, 0) == true);
		REQUIRE(mat.IsInBand(0, 1) == true);
		REQUIRE(mat.IsInBand(1, 0) == true);
		REQUIRE(mat.IsInBand(2, 2) == true);
		REQUIRE(mat.IsInBand(3, 4) == true);

		// Elements outside the band
		REQUIRE(mat.IsInBand(0, 2) == false);
		REQUIRE(mat.IsInBand(0, 3) == false);
		REQUIRE(mat.IsInBand(2, 0) == false);
		REQUIRE(mat.IsInBand(4, 1) == false);
	}

	TEST_CASE("BandDiagonalMatrix_MatrixVector_Multiplication", "[BandDiag]")
	{
			TEST_PRECISION_INFO();
		// Create a 4x4 tridiagonal matrix
		// [ 2  1  0  0 ]
		// [ 1  2  1  0 ]
		// [ 0  1  2  1 ]
		// [ 0  0  1  2 ]
		Matrix<Real> data(4, 3, {
			static_cast<MML::Real>(0.0),  static_cast<MML::Real>(2.0),  static_cast<MML::Real>(1.0),
			static_cast<MML::Real>(1.0),  static_cast<MML::Real>(2.0),  static_cast<MML::Real>(1.0),
			static_cast<MML::Real>(1.0),  static_cast<MML::Real>(2.0),  static_cast<MML::Real>(1.0),
			static_cast<MML::Real>(1.0),  static_cast<MML::Real>(2.0),  static_cast<MML::Real>(0.0)
		});

		BandDiagonalMatrix mat(4, 1, 1, data);
		Vector<Real> x{static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0), static_cast<MML::Real>(4.0)};

		Vector<Real> result = mat * x;

		// Expected result: [2*1 + 1*2, 1*1 + 2*2 + 1*3, 1*2 + 2*3 + 1*4, 1*3 + 2*4]
		//                = [4, 8, 12, 11]
		REQUIRE(result.size() == 4);
		REQUIRE(result[0] == static_cast<MML::Real>(4.0));
		REQUIRE(result[1] == static_cast<MML::Real>(8.0));
		REQUIRE(result[2] == static_cast<MML::Real>(12.0));
		REQUIRE(result[3] == static_cast<MML::Real>(11.0));
	}

	TEST_CASE("BandDiagonalMatrix_WiderBand_Multiplication", "[BandDiag]")
	{
			TEST_PRECISION_INFO();
		// Create a 5x5 matrix with m1=2, m2=1
		// [ 1  2  0  0  0 ]
		// [ 3  4  5  0  0 ]
		// [ 6  7  8  9  0 ]
		// [ 0 10 11 12 13 ]
		// [ 0  0 14 15 16 ]
		Matrix<Real> data(5, 4, {
			// [sub2, sub1, diag, super1]
			 static_cast<MML::Real>(0.0),  static_cast<MML::Real>(0.0),  static_cast<MML::Real>(1.0),  static_cast<MML::Real>(2.0),
			 static_cast<MML::Real>(0.0),  static_cast<MML::Real>(3.0),  static_cast<MML::Real>(4.0),  static_cast<MML::Real>(5.0),
			 static_cast<MML::Real>(6.0),  static_cast<MML::Real>(7.0),  static_cast<MML::Real>(8.0),  static_cast<MML::Real>(9.0),
			static_cast<MML::Real>(10.0), static_cast<MML::Real>(11.0), static_cast<MML::Real>(12.0), static_cast<MML::Real>(13.0),
			static_cast<MML::Real>(14.0), static_cast<MML::Real>(15.0), static_cast<MML::Real>(16.0),  static_cast<MML::Real>(0.0)
		});

		BandDiagonalMatrix mat(5, 2, 1, data);
		Vector<Real> x{static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0)};

		Vector<Real> result = mat * x;

		// Row 0: 1*1 + 2*1 = 3
		// Row 1: 3*1 + 4*1 + 5*1 = 12
		// Row 2: 6*1 + 7*1 + 8*1 + 9*1 = 30
		// Row 3: 10*1 + 11*1 + 12*1 + 13*1 = 46
		// Row 4: 14*1 + 15*1 + 16*1 = 45
		REQUIRE(result.size() == 5);
		REQUIRE(result[0] == static_cast<MML::Real>(3.0));
		REQUIRE(result[1] == static_cast<MML::Real>(12.0));
		REQUIRE(result[2] == static_cast<MML::Real>(30.0));
		REQUIRE(result[3] == static_cast<MML::Real>(46.0));
		REQUIRE(result[4] == static_cast<MML::Real>(45.0));
	}

	TEST_CASE("BandDiagonalMatrix_ToFullMatrix", "[BandDiag]")
	{
			TEST_PRECISION_INFO();
		// Create a 4x4 tridiagonal matrix
		Matrix<Real> data(4, 3, {
			static_cast<MML::Real>(0.0),  static_cast<MML::Real>(2.0),  static_cast<MML::Real>(1.0),
			static_cast<MML::Real>(1.0),  static_cast<MML::Real>(2.0),  static_cast<MML::Real>(1.0),
			static_cast<MML::Real>(1.0),  static_cast<MML::Real>(2.0),  static_cast<MML::Real>(1.0),
			static_cast<MML::Real>(1.0),  static_cast<MML::Real>(2.0),  static_cast<MML::Real>(0.0)
		});

		BandDiagonalMatrix mat(4, 1, 1, data);
		Matrix<Real> full = mat.ToFullMatrix();

		REQUIRE(full.rows() == 4);
		REQUIRE(full.cols() == 4);

		// Check all elements
		REQUIRE(full[0][0] == static_cast<MML::Real>(2.0));
		REQUIRE(full[0][1] == static_cast<MML::Real>(1.0));
		REQUIRE(full[0][2] == static_cast<MML::Real>(0.0));
		REQUIRE(full[0][3] == static_cast<MML::Real>(0.0));

		REQUIRE(full[1][0] == static_cast<MML::Real>(1.0));
		REQUIRE(full[1][1] == static_cast<MML::Real>(2.0));
		REQUIRE(full[1][2] == static_cast<MML::Real>(1.0));
		REQUIRE(full[1][3] == static_cast<MML::Real>(0.0));

		REQUIRE(full[2][0] == static_cast<MML::Real>(0.0));
		REQUIRE(full[2][1] == static_cast<MML::Real>(1.0));
		REQUIRE(full[2][2] == static_cast<MML::Real>(2.0));
		REQUIRE(full[2][3] == static_cast<MML::Real>(1.0));

		REQUIRE(full[3][0] == static_cast<MML::Real>(0.0));
		REQUIRE(full[3][1] == static_cast<MML::Real>(0.0));
		REQUIRE(full[3][2] == static_cast<MML::Real>(1.0));
		REQUIRE(full[3][3] == static_cast<MML::Real>(2.0));
	}

	TEST_CASE("BandDiagonalMatrix_Symmetric_Band", "[BandDiag]")
	{
			TEST_PRECISION_INFO();
		// Test with symmetric upper/lower bandwidths m1 = m2 = 2
		Matrix<Real> data(6, 5, {
			// [sub2, sub1, diag, super1, super2]
			 static_cast<MML::Real>(0.0),  static_cast<MML::Real>(0.0),  static_cast<MML::Real>(5.0),  static_cast<MML::Real>(1.0),  static_cast<MML::Real>(2.0),
			 static_cast<MML::Real>(0.0),  static_cast<MML::Real>(1.0),  static_cast<MML::Real>(5.0),  static_cast<MML::Real>(1.0),  static_cast<MML::Real>(2.0),
			 static_cast<MML::Real>(1.0),  static_cast<MML::Real>(1.0),  static_cast<MML::Real>(5.0),  static_cast<MML::Real>(1.0),  static_cast<MML::Real>(2.0),
			 static_cast<MML::Real>(1.0),  static_cast<MML::Real>(1.0),  static_cast<MML::Real>(5.0),  static_cast<MML::Real>(1.0),  static_cast<MML::Real>(2.0),
			 static_cast<MML::Real>(1.0),  static_cast<MML::Real>(1.0),  static_cast<MML::Real>(5.0),  static_cast<MML::Real>(1.0),  static_cast<MML::Real>(0.0),
			 static_cast<MML::Real>(1.0),  static_cast<MML::Real>(1.0),  static_cast<MML::Real>(5.0),  static_cast<MML::Real>(0.0),  static_cast<MML::Real>(0.0)
		});

		BandDiagonalMatrix mat(6, 2, 2, data);

		// Test diagonal
		for (int i = 0; i < 6; i++)
			REQUIRE(mat(i, i) == static_cast<MML::Real>(5.0));

		// Test first superdiagonal
		for (int i = 0; i < 5; i++)
			REQUIRE(mat(i, i + 1) == static_cast<MML::Real>(1.0));

		// Test second superdiagonal
		for (int i = 0; i < 4; i++)
			REQUIRE(mat(i, i + 2) == static_cast<MML::Real>(2.0));

		// Test first subdiagonal
		for (int i = 1; i < 6; i++)
			REQUIRE(mat(i, i - 1) == static_cast<MML::Real>(1.0));

		// Test second subdiagonal
		for (int i = 2; i < 6; i++)
			REQUIRE(mat(i, i - 2) == static_cast<MML::Real>(1.0));

		// Note: accessing elements outside band throws exception
		// REQUIRE(mat(0, 3) == static_cast<MML::Real>(0.0));   // would throw
		// REQUIRE(mat(5, 2) == static_cast<MML::Real>(0.0));   // would throw
	}

	TEST_CASE("BandDiagonalMatrix_Identity_Like", "[BandDiag]")
	{
			TEST_PRECISION_INFO();
		// Create a diagonal matrix (m1=0, m2=0)
		Matrix<Real> data(4, 1, {
			static_cast<MML::Real>(3.0),
			static_cast<MML::Real>(5.0),
			static_cast<MML::Real>(7.0),
			static_cast<MML::Real>(9.0)
		});

		BandDiagonalMatrix mat(4, 0, 0, data);

		// All diagonal elements should be non-zero
		REQUIRE(mat(0, 0) == static_cast<MML::Real>(3.0));
		REQUIRE(mat(1, 1) == static_cast<MML::Real>(5.0));
		REQUIRE(mat(2, 2) == static_cast<MML::Real>(7.0));
		REQUIRE(mat(3, 3) == static_cast<MML::Real>(9.0));

		// Note: For a diagonal-only matrix (m1=0, m2=0), all off-diagonal
		// accesses would throw exception, so we can't test them directly.
		// However, multiplication and ToFullMatrix() work correctly.

		// Test multiplication
		Vector<Real> x{static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0), static_cast<MML::Real>(4.0)};
		Vector<Real> result = mat * x;

		REQUIRE(result[0] == static_cast<MML::Real>(3.0));
		REQUIRE(result[1] == static_cast<MML::Real>(10.0));
		REQUIRE(result[2] == static_cast<MML::Real>(21.0));
		REQUIRE(result[3] == static_cast<MML::Real>(36.0));
	}

	TEST_CASE("BandDiagonalMatrix_Element_Modification", "[BandDiag]")
	{
			TEST_PRECISION_INFO();
		// Create a matrix and test element modification
		Matrix<Real> data(3, 3, {
			static_cast<MML::Real>(0.0),  static_cast<MML::Real>(1.0),  static_cast<MML::Real>(1.0),
			static_cast<MML::Real>(1.0),  static_cast<MML::Real>(1.0),  static_cast<MML::Real>(1.0),
			static_cast<MML::Real>(1.0),  static_cast<MML::Real>(1.0),  static_cast<MML::Real>(0.0)
		});

		BandDiagonalMatrix mat(3, 1, 1, data);

		// Modify elements within the band
		mat(0, 0) = static_cast<MML::Real>(10.0);
		mat(1, 2) = static_cast<MML::Real>(20.0);
		mat(2, 1) = static_cast<MML::Real>(30.0);

		REQUIRE(mat(0, 0) == static_cast<MML::Real>(10.0));
		REQUIRE(mat(1, 2) == static_cast<MML::Real>(20.0));
		REQUIRE(mat(2, 1) == static_cast<MML::Real>(30.0));

		// Original elements should remain
		REQUIRE(mat(0, 1) == static_cast<MML::Real>(1.0));
		REQUIRE(mat(1, 0) == static_cast<MML::Real>(1.0));
	}

	TEST_CASE("BandDiagonalMatrix_Exception_OutOfBand_Write", "[BandDiag]")
	{
			TEST_PRECISION_INFO();
		Matrix<Real> data(3, 3, {
			static_cast<MML::Real>(0.0),  static_cast<MML::Real>(1.0),  static_cast<MML::Real>(1.0),
			static_cast<MML::Real>(1.0),  static_cast<MML::Real>(1.0),  static_cast<MML::Real>(1.0),
			static_cast<MML::Real>(1.0),  static_cast<MML::Real>(1.0),  static_cast<MML::Real>(0.0)
		});

		BandDiagonalMatrix mat(3, 1, 1, data);

		// Attempting to write to elements outside the band should throw
		REQUIRE_THROWS_AS(mat(0, 2) = static_cast<MML::Real>(5.0), MatrixAccessBoundsError);
		REQUIRE_THROWS_AS(mat(2, 0) = static_cast<MML::Real>(5.0), MatrixAccessBoundsError);
	}

	TEST_CASE("BandDiagonalMatrix_Large_Sparse", "[BandDiag]")
	{
			TEST_PRECISION_INFO();
		// Create a larger 10x10 band matrix with m1=1, m2=1
		Matrix<Real> data(10, 3);
		for (int i = 0; i < 10; i++)
		{
			if (i > 0) data[i][0] = -static_cast<MML::Real>(1.0);  // subdiagonal
			data[i][1] = static_cast<MML::Real>(2.0);               // diagonal
			if (i < 9) data[i][2] = -static_cast<MML::Real>(1.0);  // superdiagonal
		}

		BandDiagonalMatrix mat(10, 1, 1, data);
		Vector<Real> x(10);
		for (int i = 0; i < 10; i++)
			x[i] = static_cast<MML::Real>(1.0);

		Vector<Real> result = mat * x;

		// First and last rows: 2*1 + (-1)*1 = 1
		// Middle rows: (-1)*1 + 2*1 + (-1)*1 = 0
		REQUIRE(result[0] == static_cast<MML::Real>(1.0));
		for (int i = 1; i < 9; i++)
			REQUIRE(result[i] == static_cast<MML::Real>(0.0));
		REQUIRE(result[9] == static_cast<MML::Real>(1.0));
	}

	///////////////////     ARITHMETIC OPERATIONS TESTS     ///////////////////

	TEST_CASE("BandDiagonalMatrix_Scalar_Multiplication", "[BandDiag][Arithmetic]")
	{
			TEST_PRECISION_INFO();
		Matrix<Real> data(4, 3, {
			static_cast<MML::Real>(0.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(1.0),
			static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(1.0),
			static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(1.0),
			static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(0.0)
		});
		BandDiagonalMatrix mat(4, 1, 1, data);

		// Test scalar multiplication
		BandDiagonalMatrix result1 = mat * static_cast<MML::Real>(3.0);
		REQUIRE(result1(0, 0) == static_cast<MML::Real>(6.0));
		REQUIRE(result1(1, 1) == static_cast<MML::Real>(6.0));
		REQUIRE(result1(0, 1) == static_cast<MML::Real>(3.0));
		REQUIRE(result1(1, 0) == static_cast<MML::Real>(3.0));

		// Test left scalar multiplication
		BandDiagonalMatrix result2 = static_cast<MML::Real>(2.0) * mat;
		REQUIRE(result2(0, 0) == static_cast<MML::Real>(4.0));
		REQUIRE(result2(1, 1) == static_cast<MML::Real>(4.0));
		REQUIRE(result2(0, 1) == static_cast<MML::Real>(2.0));

		// Test in-place scalar multiplication
		mat *= static_cast<MML::Real>(2.0);
		REQUIRE(mat(0, 0) == static_cast<MML::Real>(4.0));
		REQUIRE(mat(1, 1) == static_cast<MML::Real>(4.0));
	}

	TEST_CASE("BandDiagonalMatrix_Scalar_Division", "[BandDiag][Arithmetic]")
	{
			TEST_PRECISION_INFO();
		Matrix<Real> data(4, 3, {
			static_cast<MML::Real>(0.0), static_cast<MML::Real>(4.0), static_cast<MML::Real>(2.0),
			static_cast<MML::Real>(2.0), static_cast<MML::Real>(4.0), static_cast<MML::Real>(2.0),
			static_cast<MML::Real>(2.0), static_cast<MML::Real>(4.0), static_cast<MML::Real>(2.0),
			static_cast<MML::Real>(2.0), static_cast<MML::Real>(4.0), static_cast<MML::Real>(0.0)
		});
		BandDiagonalMatrix mat(4, 1, 1, data);

		BandDiagonalMatrix result = mat / static_cast<MML::Real>(2.0);
		REQUIRE(result(0, 0) == static_cast<MML::Real>(2.0));
		REQUIRE(result(1, 1) == static_cast<MML::Real>(2.0));
		REQUIRE(result(0, 1) == static_cast<MML::Real>(1.0));
		REQUIRE(result(1, 0) == static_cast<MML::Real>(1.0));

		// Test division by zero
		REQUIRE_THROWS_AS(mat / static_cast<MML::Real>(0.0), MML::DivisionByZeroError);

		// Test in-place division
		mat /= static_cast<MML::Real>(4.0);
		REQUIRE(mat(0, 0) == static_cast<MML::Real>(1.0));
		REQUIRE(mat(1, 1) == static_cast<MML::Real>(1.0));
	}

	TEST_CASE("BandDiagonalMatrix_Negation", "[BandDiag][Arithmetic]")
	{
			TEST_PRECISION_INFO();
		Matrix<Real> data(4, 3, {
			static_cast<MML::Real>(0.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(1.0),
			static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(1.0),
			static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(1.0),
			static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(0.0)
		});
		BandDiagonalMatrix mat(4, 1, 1, data);

		BandDiagonalMatrix neg = -mat;
		REQUIRE(neg(0, 0) == -static_cast<MML::Real>(2.0));
		REQUIRE(neg(1, 1) == -static_cast<MML::Real>(2.0));
		REQUIRE(neg(0, 1) == -static_cast<MML::Real>(1.0));
		REQUIRE(neg(1, 0) == -static_cast<MML::Real>(1.0));
	}

	TEST_CASE("BandDiagonalMatrix_Matrix_Addition", "[BandDiag][Arithmetic]")
	{
			TEST_PRECISION_INFO();
		Matrix<Real> data1(4, 3, {
			static_cast<MML::Real>(0.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(1.0),
			static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(1.0),
			static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(1.0),
			static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(0.0)
		});
		Matrix<Real> data2(4, 3, {
			static_cast<MML::Real>(0.0), static_cast<MML::Real>(3.0), static_cast<MML::Real>(2.0),
			static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0), static_cast<MML::Real>(2.0),
			static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0), static_cast<MML::Real>(2.0),
			static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0), static_cast<MML::Real>(0.0)
		});

		BandDiagonalMatrix mat1(4, 1, 1, data1);
		BandDiagonalMatrix mat2(4, 1, 1, data2);

		BandDiagonalMatrix sum = mat1 + mat2;
		REQUIRE(sum(0, 0) == static_cast<MML::Real>(5.0));
		REQUIRE(sum(1, 1) == static_cast<MML::Real>(5.0));
		REQUIRE(sum(0, 1) == static_cast<MML::Real>(3.0));
		REQUIRE(sum(1, 0) == static_cast<MML::Real>(3.0));

		// Test in-place addition
		mat1 += mat2;
		REQUIRE(mat1(0, 0) == static_cast<MML::Real>(5.0));
		REQUIRE(mat1(1, 1) == static_cast<MML::Real>(5.0));
	}

	TEST_CASE("BandDiagonalMatrix_Matrix_Subtraction", "[BandDiag][Arithmetic]")
	{
			TEST_PRECISION_INFO();
		Matrix<Real> data1(4, 3, {
			static_cast<MML::Real>(0.0), static_cast<MML::Real>(5.0), static_cast<MML::Real>(3.0),
			static_cast<MML::Real>(3.0), static_cast<MML::Real>(5.0), static_cast<MML::Real>(3.0),
			static_cast<MML::Real>(3.0), static_cast<MML::Real>(5.0), static_cast<MML::Real>(3.0),
			static_cast<MML::Real>(3.0), static_cast<MML::Real>(5.0), static_cast<MML::Real>(0.0)
		});
		Matrix<Real> data2(4, 3, {
			static_cast<MML::Real>(0.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(1.0),
			static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(1.0),
			static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(1.0),
			static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(0.0)
		});

		BandDiagonalMatrix mat1(4, 1, 1, data1);
		BandDiagonalMatrix mat2(4, 1, 1, data2);

		BandDiagonalMatrix diff = mat1 - mat2;
		REQUIRE(diff(0, 0) == static_cast<MML::Real>(3.0));
		REQUIRE(diff(1, 1) == static_cast<MML::Real>(3.0));
		REQUIRE(diff(0, 1) == static_cast<MML::Real>(2.0));
		REQUIRE(diff(1, 0) == static_cast<MML::Real>(2.0));

		// Test in-place subtraction
		mat1 -= mat2;
		REQUIRE(mat1(0, 0) == static_cast<MML::Real>(3.0));
		REQUIRE(mat1(1, 1) == static_cast<MML::Real>(3.0));
	}

	TEST_CASE("BandDiagonalMatrix_Incompatible_Operations", "[BandDiag][Arithmetic]")
	{
			TEST_PRECISION_INFO();
		Matrix<Real> data1(4, 3);  // m1=1, m2=1
		Matrix<Real> data2(4, 4);  // m1=1, m2=2 (different bandwidth)
		Matrix<Real> data3(5, 3);  // Different dimension

		BandDiagonalMatrix mat1(4, 1, 1, data1);
		BandDiagonalMatrix mat2(4, 1, 2, data2);
		BandDiagonalMatrix mat3(5, 1, 1, data3);

		// Different bandwidth
		REQUIRE_THROWS_AS(mat1 + mat2, std::invalid_argument);
		REQUIRE_THROWS_AS(mat1 - mat2, std::invalid_argument);

		// Different dimensions
		REQUIRE_THROWS_AS(mat1 + mat3, std::invalid_argument);
		REQUIRE_THROWS_AS(mat1 - mat3, std::invalid_argument);
	}

	TEST_CASE("BandDiagonalMatrix_Trace", "[BandDiag][Properties]")
	{
			TEST_PRECISION_INFO();
		Matrix<Real> data(4, 3, {
			static_cast<MML::Real>(0.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(1.0),
			static_cast<MML::Real>(1.0), static_cast<MML::Real>(3.0), static_cast<MML::Real>(1.0),
			static_cast<MML::Real>(1.0), static_cast<MML::Real>(5.0), static_cast<MML::Real>(1.0),
			static_cast<MML::Real>(1.0), static_cast<MML::Real>(7.0), static_cast<MML::Real>(0.0)
		});
		BandDiagonalMatrix mat(4, 1, 1, data);

		Real trace = mat.Trace();
		REQUIRE(trace == static_cast<MML::Real>(17.0));  // 2 + 3 + 5 + 7
	}

	TEST_CASE("BandDiagonalMatrix_Norms", "[BandDiag][Properties]")
	{
			TEST_PRECISION_INFO();
		Matrix<Real> data(3, 3, {
			static_cast<MML::Real>(0.0), static_cast<MML::Real>(3.0), static_cast<MML::Real>(4.0),
			static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(1.0),
			static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0)
		});
		BandDiagonalMatrix mat(3, 1, 1, data);

		// Frobenius norm: sqrt(3^2 + 4^2 + 1^2 + 2^2 + 1^2 + 1^2 + 1^2) = sqrt(9+16+1+4+1+1+1) = sqrt(33)
		Real frob = mat.NormFrobenius();
		REQUIRE(std::abs(frob - std::sqrt(static_cast<MML::Real>(33.0))) < TOL(1e-10, 1e-5));

		// Infinity norm: max row sum = max(7, 4, 2) = 7
		Real inf = mat.NormInf();
		REQUIRE(inf == static_cast<MML::Real>(7.0));
	}

	TEST_CASE("BandDiagonalMatrix_Transpose", "[BandDiag][Properties]")
	{
			TEST_PRECISION_INFO();
		// Create a non-symmetric band matrix
		Matrix<Real> data(4, 4, {
			// m1=1, m2=2: [sub1, diag, super1, super2]
			static_cast<MML::Real>(0.0),  static_cast<MML::Real>(1.0),  static_cast<MML::Real>(2.0),  static_cast<MML::Real>(3.0),
			static_cast<MML::Real>(4.0),  static_cast<MML::Real>(5.0),  static_cast<MML::Real>(6.0),  static_cast<MML::Real>(7.0),
			static_cast<MML::Real>(8.0),  static_cast<MML::Real>(9.0), static_cast<MML::Real>(10.0), static_cast<MML::Real>(11.0),
			static_cast<MML::Real>(12.0), static_cast<MML::Real>(13.0), static_cast<MML::Real>(14.0),  static_cast<MML::Real>(0.0)
		});
		BandDiagonalMatrix mat(4, 1, 2, data);

		BandDiagonalMatrix trans = mat.Transposed();

		// After transpose, m1 and m2 are swapped
		REQUIRE(trans.getLowerBandwidth() == 2);
		REQUIRE(trans.getUpperBandwidth() == 1);

		// Verify transpose property: A^T[i][j] = A[j][i] for elements in the band
		REQUIRE(trans(0, 0) == mat(0, 0));  // Diagonal
		REQUIRE(trans(0, 1) == mat(1, 0));  // trans superdiag = mat subdiag
		REQUIRE(trans(1, 0) == mat(0, 1));  // trans subdiag = mat superdiag  
		REQUIRE(trans(1, 2) == mat(2, 1));  // mat(2,1) is in band (m1=1)
		REQUIRE(trans(2, 1) == mat(1, 2));  // mat(1,2) is in band (m2=2)
		REQUIRE(trans(2, 0) == mat(0, 2));  // mat(0,2) is in band (m2=2), trans(2,0) in band (new m1=2)
	}

	TEST_CASE("BandDiagonalMatrix_IsSymmetric", "[BandDiag][Properties]")
	{
			TEST_PRECISION_INFO();
		// Symmetric tridiagonal matrix
		Matrix<Real> data_sym(4, 3, {
			static_cast<MML::Real>(0.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(1.0),
			static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(1.0),
			static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(1.0),
			static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(0.0)
		});
		BandDiagonalMatrix mat_sym(4, 1, 1, data_sym);
		REQUIRE(mat_sym.IsSymmetric());

		// Non-symmetric matrix
		Matrix<Real> data_nonsym(4, 3, {
			static_cast<MML::Real>(0.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(1.0),
			static_cast<MML::Real>(2.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(1.0),  // Different subdiagonal value
			static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(1.0),
			static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(0.0)
		});
		BandDiagonalMatrix mat_nonsym(4, 1, 1, data_nonsym);
		REQUIRE_FALSE(mat_nonsym.IsSymmetric());

		// Different bandwidths => not symmetric
		Matrix<Real> data_asym(4, 4);
		BandDiagonalMatrix mat_asym(4, 1, 2, data_asym);
		REQUIRE_FALSE(mat_asym.IsSymmetric());
	}

	TEST_CASE("BandDiagonalMatrix_Equality", "[BandDiag][Properties]")
	{
			TEST_PRECISION_INFO();
		Matrix<Real> data1(4, 3, {
			static_cast<MML::Real>(0.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(1.0),
			static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(1.0),
			static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(1.0),
			static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(0.0)
		});
		Matrix<Real> data2(4, 3, {
			static_cast<MML::Real>(0.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(1.0),
			static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(1.0),
			static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(1.0),
			static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(0.0)
		});

		BandDiagonalMatrix mat1(4, 1, 1, data1);
		BandDiagonalMatrix mat2(4, 1, 1, data2);

		REQUIRE(mat1.IsEqualTo(mat2));

		// Modify one element
		data2[1][1] = static_cast<MML::Real>(3.0);
		BandDiagonalMatrix mat3(4, 1, 1, data2);
		REQUIRE_FALSE(mat1.IsEqualTo(mat3));
	}

	TEST_CASE("BandDiagonalMatrix_GetDiagonal", "[BandDiag][Properties]")
	{
			TEST_PRECISION_INFO();
		Matrix<Real> data(5, 4, {
			// m1=1, m2=2
			static_cast<MML::Real>(0.0),  static_cast<MML::Real>(1.0),  static_cast<MML::Real>(2.0),  static_cast<MML::Real>(3.0),
			static_cast<MML::Real>(4.0),  static_cast<MML::Real>(5.0),  static_cast<MML::Real>(6.0),  static_cast<MML::Real>(7.0),
			static_cast<MML::Real>(8.0),  static_cast<MML::Real>(9.0), static_cast<MML::Real>(10.0), static_cast<MML::Real>(11.0),
			static_cast<MML::Real>(12.0), static_cast<MML::Real>(13.0), static_cast<MML::Real>(14.0),  static_cast<MML::Real>(0.0),
			static_cast<MML::Real>(16.0), static_cast<MML::Real>(17.0),  static_cast<MML::Real>(0.0),  static_cast<MML::Real>(0.0)
		});
		BandDiagonalMatrix mat(5, 1, 2, data);

		// Main diagonal
		Vector<Real> diag = mat.getDiagonal();
		REQUIRE(diag.size() == 5);
		REQUIRE(diag[0] == static_cast<MML::Real>(1.0));
		REQUIRE(diag[1] == static_cast<MML::Real>(5.0));
		REQUIRE(diag[2] == static_cast<MML::Real>(9.0));
		REQUIRE(diag[3] == static_cast<MML::Real>(13.0));
		REQUIRE(diag[4] == static_cast<MML::Real>(17.0));

		// First subdiagonal (k=1)
		Vector<Real> sub1 = mat.getDiagonal(1);
		REQUIRE(sub1.size() == 4);
		REQUIRE(sub1[0] == static_cast<MML::Real>(4.0));
		REQUIRE(sub1[1] == static_cast<MML::Real>(8.0));
		REQUIRE(sub1[2] == static_cast<MML::Real>(12.0));
		REQUIRE(sub1[3] == static_cast<MML::Real>(16.0));

		// First superdiagonal (k=-1)
		Vector<Real> super1 = mat.getDiagonal(-1);
		REQUIRE(super1.size() == 4);
		REQUIRE(super1[0] == static_cast<MML::Real>(2.0));
		REQUIRE(super1[1] == static_cast<MML::Real>(6.0));
		REQUIRE(super1[2] == static_cast<MML::Real>(10.0));
		REQUIRE(super1[3] == static_cast<MML::Real>(14.0));

		// Second superdiagonal (k=-2)
		Vector<Real> super2 = mat.getDiagonal(-2);
		REQUIRE(super2.size() == 3);
		REQUIRE(super2[0] == static_cast<MML::Real>(3.0));
		REQUIRE(super2[1] == static_cast<MML::Real>(7.0));
		REQUIRE(super2[2] == static_cast<MML::Real>(11.0));

		// Out of band diagonal should throw
		REQUIRE_THROWS_AS(mat.getDiagonal(2), MML::IndexError);
		REQUIRE_THROWS_AS(mat.getDiagonal(-3), MML::IndexError);
	}

	TEST_CASE("BandDiagonalMatrix_Combined_Operations", "[BandDiag][Arithmetic]")
	{
			TEST_PRECISION_INFO();
		// Test combining multiple operations
		Matrix<Real> data(4, 3, {
			static_cast<MML::Real>(0.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(1.0),
			static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(1.0),
			static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(1.0),
			static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(0.0)
		});
		BandDiagonalMatrix mat(4, 1, 1, data);

		// (2*A - A) + A = 2*A
		BandDiagonalMatrix result = (static_cast<MML::Real>(2.0) * mat - mat) + mat;
		REQUIRE(result(0, 0) == static_cast<MML::Real>(4.0));
		REQUIRE(result(1, 1) == static_cast<MML::Real>(4.0));
		REQUIRE(result(0, 1) == static_cast<MML::Real>(2.0));

		// Test trace after operations
		REQUIRE(result.Trace() == static_cast<MML::Real>(16.0));  // 4*4
	}

	TEST_CASE("BandDiagonalMatrix_operator()_rejects_indices_outside_the_matrix", "[BandDiag][wewi.1.5]")
	{
		// Only the band test existed: (i,j) outside the matrix but satisfying the band inequality
		// (e.g. (5,5) on a 5x5, or negative indices) read/wrote past the compact storage.
		Matrix<Real> data(5, 3, {
			static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0),
			static_cast<MML::Real>(3.0), static_cast<MML::Real>(4.0), static_cast<MML::Real>(5.0),
			static_cast<MML::Real>(6.0), static_cast<MML::Real>(7.0), static_cast<MML::Real>(8.0),
			static_cast<MML::Real>(9.0), static_cast<MML::Real>(10.0), static_cast<MML::Real>(11.0),
			static_cast<MML::Real>(12.0), static_cast<MML::Real>(13.0), static_cast<MML::Real>(0.0)
		});
		BandDiagonalMatrix mat(5, 1, 1, data);
		const BandDiagonalMatrix& cmat = mat;

		REQUIRE(cmat(0, 0) == static_cast<MML::Real>(1.0));
		REQUIRE(cmat(0, 3) == static_cast<MML::Real>(0.0));   // inside the matrix, outside the band
		REQUIRE_THROWS_AS(cmat(5, 5), MatrixAccessBoundsError);
		REQUIRE_THROWS_AS(cmat(-1, 0), MatrixAccessBoundsError);
		REQUIRE_THROWS_AS(cmat(0, 5), MatrixAccessBoundsError);
		REQUIRE_THROWS_AS(cmat(4, -1), MatrixAccessBoundsError);

		REQUIRE_THROWS_AS(mat(5, 5) = static_cast<MML::Real>(1.0), MatrixAccessBoundsError);
		REQUIRE_THROWS_AS(mat(-1, 0) = static_cast<MML::Real>(1.0), MatrixAccessBoundsError);
		REQUIRE_THROWS_AS(mat(0, 3) = static_cast<MML::Real>(1.0), MatrixAccessBoundsError);  // outside band, still rejected
		REQUIRE_NOTHROW(mat(4, 4) = static_cast<MML::Real>(2.0));
		REQUIRE(cmat(4, 4) == static_cast<MML::Real>(2.0));
	}
}


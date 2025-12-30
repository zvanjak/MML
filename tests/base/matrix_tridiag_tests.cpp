#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Vector.h"
#include "base/MatrixTriDiag.h"

#include "core/LinAlgEqSolvers.h"
#endif

using namespace MML;
using namespace MML::Testing;


// TODO REAL(0.9) - Tridiag tests
namespace MML::Tests::Base::TridiagMatrixTests
{
	TEST_CASE("MatrixTridiag_init_from_3_vectors", "[simple]")
	{
			TEST_PRECISION_INFO();
		// initializing with 3 vectors
		TridiagonalMatrix<Real> a(4, Vector<Real>{REAL(0.0), REAL(4.5), REAL(9.0), REAL(10.0)}, Vector<Real>{REAL(4.0), REAL(1.5), REAL(6.0), REAL(7.0)}, Vector<Real>{REAL(1.0), REAL(2.0), REAL(3.0), REAL(0.0)});

		// expected matrix
		Matrix<Real> b(4, 4, {  REAL(4.0), REAL(1.0), REAL(0.0), REAL(0.0),
														REAL(4.5), REAL(1.5), REAL(2.0), REAL(0.0),
														REAL(0.0), REAL(9.0), REAL(6.0), REAL(3.0),
														REAL(0.0), REAL(0.0), REAL(10.0), REAL(7.0) });

		REQUIRE(4 == a.RowNum());
		REQUIRE(4 == a.ColNum());

		REQUIRE(a(0, 0) == REAL(4.0));
		REQUIRE(a(0, 1) == REAL(1.0));
		REQUIRE(a(1, 0) == REAL(4.5));
		REQUIRE(a(1, 1) == REAL(1.5));
		REQUIRE(a(1, 2) == REAL(2.0));
		REQUIRE(a(2, 1) == REAL(9.0));
		REQUIRE(a(2, 2) == REAL(6.0));
		REQUIRE(a(2, 3) == REAL(3.0));
		REQUIRE(a(3, 2) == REAL(10.0));
		REQUIRE(a(3, 3) == REAL(7.0));

		// these throw exception
		// REQUIRE(REAL(0.0) == a(0,2));
		// REQUIRE(a(0,3) == REAL(0.0));
		// REQUIRE(a(1,3) == REAL(0.0));
		// REQUIRE(a(2,0) == REAL(0.0));
		// REQUIRE(a(3,0) == REAL(0.0));
		// REQUIRE(a(3,1) == REAL(0.0));
	}
	TEST_CASE("MatrixTridiag_init_from_initializer_list", "[simple]")
	{
			TEST_PRECISION_INFO();
		// initializing with values in single initializer list
		TridiagonalMatrix<Real> a(4, {  REAL(4.0), REAL(1.0),
																		REAL(4.5), REAL(1.5),  REAL(2.0),
																				 REAL(9.0),  REAL(6.0), REAL(3.0),
																						  REAL(10.0), REAL(7.0) });

		// expected matrix
		Matrix<Real> b(4, 4, {  REAL(4.0), REAL(1.0), REAL(0.0), REAL(0.0),
														REAL(4.5), REAL(1.5), REAL(2.0), REAL(0.0),
														REAL(0.0), REAL(9.0), REAL(6.0), REAL(3.0),
														REAL(0.0), REAL(0.0), REAL(10.0), REAL(7.0) });


		REQUIRE(4 == a.RowNum());
		REQUIRE(4 == a.ColNum());

		REQUIRE(a(0, 0) == REAL(4.0));
		REQUIRE(a(0, 1) == REAL(1.0));
		REQUIRE(a(1, 0) == REAL(4.5));
		REQUIRE(a(1, 1) == REAL(1.5));
		REQUIRE(a(1, 2) == REAL(2.0));
		REQUIRE(a(2, 1) == REAL(9.0));
		REQUIRE(a(2, 2) == REAL(6.0));
		REQUIRE(a(2, 3) == REAL(3.0));
		REQUIRE(a(3, 2) == REAL(10.0));
		REQUIRE(a(3, 3) == REAL(7.0));
	}

	TEST_CASE("MatrixTridiag_IsEqual", "[tridiag][equality]")
	{
			TEST_PRECISION_INFO();
		TridiagonalMatrix<Real> a(4, { REAL(4.0), REAL(1.0),
																	 REAL(4.5), REAL(1.5), REAL(2.0),
																			REAL(9.0), REAL(6.0), REAL(3.0),
																				 REAL(10.0), REAL(7.0) });

		TridiagonalMatrix<Real> b(4, { REAL(4.0), REAL(1.0),
																	 REAL(4.5), REAL(1.5), REAL(2.0),
																			REAL(9.0), REAL(6.0), REAL(3.0),
																				 REAL(10.0), REAL(7.0) });

		TridiagonalMatrix<Real> c(4, { REAL(4.0), REAL(1.0),
																	 REAL(4.5), REAL(1.5), REAL(2.1),  // different value
																			REAL(9.0), REAL(6.0), REAL(3.0),
																				 REAL(10.0), REAL(7.0) });

		REQUIRE(a.IsEqual(b));
		REQUIRE(b.IsEqual(a));
		REQUIRE_FALSE(a.IsEqual(c));
		REQUIRE_FALSE(c.IsEqual(a));
	}

	TEST_CASE("MatrixTridiag_IsEqual_with_tolerance", "[tridiag][equality]")
	{
			TEST_PRECISION_INFO();
		TridiagonalMatrix<Real> a(3, { REAL(2.0), REAL(1.0),
																	 REAL(1.0), REAL(2.0), REAL(1.0),
																			REAL(1.0), REAL(2.0) });

		TridiagonalMatrix<Real> b(3, { REAL(2.000001), REAL(1.000001),
																	 REAL(1.000001), REAL(2.000001), REAL(1.000001),
																			REAL(1.000001), REAL(2.000001) });

		REQUIRE(a.IsEqual(b, 1e-5));     // Within tolerance
		REQUIRE_FALSE(a.IsEqual(b, 1e-7));  // Outside tolerance
	}

	TEST_CASE("MatrixTridiag_GetTranspose", "[tridiag][transpose]")
	{
			TEST_PRECISION_INFO();
		// Create a tridiagonal matrix
		TridiagonalMatrix<Real> a(4, { REAL(4.0), REAL(1.0),
																	 REAL(4.5), REAL(1.5), REAL(2.0),
																			REAL(9.0), REAL(6.0), REAL(3.0),
																				 REAL(10.0), REAL(7.0) });

		// Transpose swaps above and below diagonals
		// Original:     below   diag   above
		//               [0, REAL(4.5), 9, 10] [4, REAL(1.5), 6, 7] [1, 2, 3, 0]
		// Transposed:   [1, 2, 3, 0]    [4, REAL(1.5), 6, 7] [0, REAL(4.5), 9, 10]

		TridiagonalMatrix<Real> at = a.GetTranspose();

		REQUIRE(at.RowNum() == 4);
		REQUIRE(at.ColNum() == 4);

		// Check diagonal (should be unchanged)
		REQUIRE(at(0, 0) == REAL(4.0));
		REQUIRE(at(1, 1) == REAL(1.5));
		REQUIRE(at(2, 2) == REAL(6.0));
		REQUIRE(at(3, 3) == REAL(7.0));

		// Check that transpose correctly swaps (i,j) with (j,i)
		REQUIRE(at(0, 1) == REAL(4.5));   // was a(1,0) = below[1] = REAL(4.5)
		REQUIRE(at(1, 0) == REAL(1.0));   // was a(0,1) = above[0] = 1
		REQUIRE(at(1, 2) == REAL(9.0));   // was a(2,1) = below[2] = 9
		REQUIRE(at(2, 1) == REAL(2.0));   // was a(1,2) = above[1] = 2
		REQUIRE(at(2, 3) == REAL(10.0));  // was a(3,2) = below[3] = 10
		REQUIRE(at(3, 2) == REAL(3.0));   // was a(2,3) = above[2] = 3
	}

	TEST_CASE("MatrixTridiag_GetTranspose_symmetric", "[tridiag][transpose]")
	{
			TEST_PRECISION_INFO();
		// Symmetric tridiagonal matrix (transpose equals itself)
		TridiagonalMatrix<Real> a(3, { REAL(2.0), REAL(1.0),
																	 REAL(1.0), REAL(2.0), REAL(1.0),
																			REAL(1.0), REAL(2.0) });

		TridiagonalMatrix<Real> at = a.GetTranspose();

		REQUIRE(a.IsEqual(at, 1e-10));
	}

	TEST_CASE("MatrixTridiag_GetInverse_throws", "[tridiag][inverse]")
	{
			TEST_PRECISION_INFO();
		TridiagonalMatrix<Real> a(3, { REAL(2.0), REAL(1.0),
																	 REAL(1.0), REAL(2.0), REAL(1.0),
																			REAL(1.0), REAL(2.0) });

		// Inverse of tridiagonal matrix is not tridiagonal, so should throw
		REQUIRE_THROWS(a.GetInverse());
	}

	//TEST_CASE("MatrixTridiag_Solve", "[simple]") 
	//{
	//    // initializing with values in single initializer list
	//    TridiagonalMatrix<Real> a(4, {  REAL(4.0), REAL(1.0),
	//                                    REAL(4.5), REAL(1.5),  REAL(2.0), 
	//                                        REAL(9.0),  REAL(6.0), REAL(3.0),
	//                                            REAL(10.0), REAL(7.0) });
	//    Matrix<Real> b(4, 4, {  REAL(4.0), REAL(1.0), REAL(0.0), REAL(0.0),
	//                            REAL(4.5), REAL(1.5), REAL(2.0), REAL(0.0), 
	//                            REAL(0.0), REAL(9.0), REAL(6.0), REAL(3.0),
	//                            REAL(0.0), REAL(0.0), REAL(10.0), REAL(7.0) });   
	//    Vector<Real> rhs{REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0)};

	//    Vector<Real> sol_a = a.Solve(rhs);

	//    REQUIRE( 4 == sol_a.size());
	//   
	//    LUDecompositionSolver<Real> solver(b); 
	//    Vector<Real> solLU = solver.Solve(rhs);     

	//    REQUIRE_THAT(sol_a[0] , WithinRel(solLU[0], REAL(1e-5)));
	//    REQUIRE_THAT(sol_a[1] , WithinRel(solLU[1], REAL(1e-5)));
	//    REQUIRE_THAT(sol_a[2] , WithinRel(solLU[2], REAL(1e-5)));
	//    REQUIRE_THAT(sol_a[3] , WithinRel(solLU[3], REAL(1e-5)));         
	//}
}

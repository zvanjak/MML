#include <catch2/catch_all.hpp>
#include "../../TestPrecision.h"
#include "../../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Vector/Vector.h"
#include "base/Matrix/MatrixTriDiag.h"

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

		REQUIRE(4 == a.rows());
		REQUIRE(4 == a.cols());

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


		REQUIRE(4 == a.rows());
		REQUIRE(4 == a.cols());

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

		REQUIRE(a.IsEqualTo(b));
		REQUIRE(b.IsEqualTo(a));
		REQUIRE_FALSE(a.IsEqualTo(c));
		REQUIRE_FALSE(c.IsEqualTo(a));
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

		REQUIRE(a.IsEqualTo(b, 1e-5));     // Within tolerance
		REQUIRE_FALSE(a.IsEqualTo(b, 1e-7));  // Outside tolerance
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

		REQUIRE(at.rows() == 4);
		REQUIRE(at.cols() == 4);

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

		REQUIRE(a.IsEqualTo(at, 1e-10));
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

	///////////////////////          Arithmetic Operators              //////////////////////
	TEST_CASE("MatrixTridiag_operator_plus", "[tridiag][arithmetic]")
	{
		TEST_PRECISION_INFO();
		TridiagonalMatrix<Real> a(3, { REAL(2.0), REAL(1.0),
		                               REAL(1.0), REAL(3.0), REAL(1.0),
		                                    REAL(1.0), REAL(4.0) });
		
		TridiagonalMatrix<Real> b(3, { REAL(1.0), REAL(2.0),
		                               REAL(2.0), REAL(1.0), REAL(3.0),
		                                    REAL(3.0), REAL(1.0) });
		
		auto c = a + b;
		
		REQUIRE(c(0, 0) == REAL(3.0));  // 2 + 1
		REQUIRE(c(0, 1) == REAL(3.0));  // 1 + 2
		REQUIRE(c(1, 0) == REAL(3.0));  // 1 + 2
		REQUIRE(c(1, 1) == REAL(4.0));  // 3 + 1
		REQUIRE(c(1, 2) == REAL(4.0));  // 1 + 3
		REQUIRE(c(2, 1) == REAL(4.0));  // 1 + 3
		REQUIRE(c(2, 2) == REAL(5.0));  // 4 + 1
	}

	TEST_CASE("MatrixTridiag_operator_minus", "[tridiag][arithmetic]")
	{
		TEST_PRECISION_INFO();
		TridiagonalMatrix<Real> a(3, { REAL(4.0), REAL(3.0),
		                               REAL(2.0), REAL(5.0), REAL(2.0),
		                                    REAL(1.0), REAL(6.0) });
		
		TridiagonalMatrix<Real> b(3, { REAL(1.0), REAL(1.0),
		                               REAL(1.0), REAL(1.0), REAL(1.0),
		                                    REAL(1.0), REAL(1.0) });
		
		auto c = a - b;
		
		REQUIRE(c(0, 0) == REAL(3.0));  // 4 - 1
		REQUIRE(c(0, 1) == REAL(2.0));  // 3 - 1
		REQUIRE(c(1, 0) == REAL(1.0));  // 2 - 1
		REQUIRE(c(1, 1) == REAL(4.0));  // 5 - 1
		REQUIRE(c(2, 2) == REAL(5.0));  // 6 - 1
	}

	TEST_CASE("MatrixTridiag_scalar_multiply", "[tridiag][arithmetic]")
	{
		TEST_PRECISION_INFO();
		TridiagonalMatrix<Real> a(3, { REAL(2.0), REAL(1.0),
		                               REAL(1.0), REAL(3.0), REAL(1.0),
		                                    REAL(1.0), REAL(4.0) });
		
		auto c1 = a * REAL(2.0);
		auto c2 = REAL(3.0) * a;
		
		REQUIRE(c1(0, 0) == REAL(4.0));  // 2 * 2
		REQUIRE(c1(1, 1) == REAL(6.0));  // 3 * 2
		REQUIRE(c1(2, 2) == REAL(8.0));  // 4 * 2
		
		REQUIRE(c2(0, 0) == REAL(6.0));  // 2 * 3
		REQUIRE(c2(1, 1) == REAL(9.0));  // 3 * 3
		REQUIRE(c2(2, 2) == REAL(12.0)); // 4 * 3
	}

	TEST_CASE("MatrixTridiag_scalar_divide", "[tridiag][arithmetic]")
	{
		TEST_PRECISION_INFO();
		TridiagonalMatrix<Real> a(3, { REAL(4.0), REAL(2.0),
		                               REAL(6.0), REAL(8.0), REAL(4.0),
		                                    REAL(10.0), REAL(12.0) });
		
		auto c = a / REAL(2.0);
		
		REQUIRE(c(0, 0) == REAL(2.0));   // 4 / 2
		REQUIRE(c(0, 1) == REAL(1.0));   // 2 / 2
		REQUIRE(c(1, 0) == REAL(3.0));   // 6 / 2
		REQUIRE(c(1, 1) == REAL(4.0));   // 8 / 2
		REQUIRE(c(2, 2) == REAL(6.0));   // 12 / 2
	}

	TEST_CASE("MatrixTridiag_dimension_mismatch_throws", "[tridiag][arithmetic]")
	{
		TEST_PRECISION_INFO();
		TridiagonalMatrix<Real> a(3, { REAL(2.0), REAL(1.0),
		                               REAL(1.0), REAL(3.0), REAL(1.0),
		                                    REAL(1.0), REAL(4.0) });
		
		TridiagonalMatrix<Real> b(4, { REAL(1.0), REAL(1.0),
		                               REAL(1.0), REAL(1.0), REAL(1.0),
		                                    REAL(1.0), REAL(1.0), REAL(1.0),
		                                         REAL(1.0), REAL(1.0) });
		
		REQUIRE_THROWS(a + b);
		REQUIRE_THROWS(a - b);
	}

	///////////////////////          Solve                             //////////////////////
	TEST_CASE("MatrixTridiag_Solve_simple", "[tridiag][solve]")
	{
		TEST_PRECISION_INFO();
		// Simple tridiagonal system
		// | 2  1  0 |   | x1 |   | 1 |
		// | 1  2  1 | * | x2 | = | 2 |
		// | 0  1  2 |   | x3 |   | 1 |
		TridiagonalMatrix<Real> a(3, { REAL(2.0), REAL(1.0),
		                               REAL(1.0), REAL(2.0), REAL(1.0),
		                                    REAL(1.0), REAL(2.0) });
		
		Vector<Real> rhs{ REAL(1.0), REAL(2.0), REAL(1.0) };
		
		Vector<Real> sol = a.Solve(rhs);
		
		// Verify A*x = b
		REQUIRE_THAT(a(0,0)*sol[0] + a(0,1)*sol[1], RealApprox(REAL(1.0)).margin(Tolerance::Strict));
		REQUIRE_THAT(a(1,0)*sol[0] + a(1,1)*sol[1] + a(1,2)*sol[2], RealApprox(REAL(2.0)).margin(Tolerance::Strict));
		REQUIRE_THAT(a(2,1)*sol[1] + a(2,2)*sol[2], RealApprox(REAL(1.0)).margin(Tolerance::Strict));
	}

	TEST_CASE("MatrixTridiag_Solve_4x4", "[tridiag][solve]")
	{
		TEST_PRECISION_INFO();
		// Larger tridiagonal system
		TridiagonalMatrix<Real> a(4, { REAL(4.0), REAL(1.0),
		                               REAL(1.0), REAL(4.0), REAL(1.0),
		                                    REAL(1.0), REAL(4.0), REAL(1.0),
		                                         REAL(1.0), REAL(4.0) });
		
		Vector<Real> rhs{ REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) };
		
		Vector<Real> sol = a.Solve(rhs);
		
		// Verify solution dimensions
		REQUIRE(sol.size() == 4);
		
		// Verify A*x = b for first and last rows
		REQUIRE_THAT(a(0,0)*sol[0] + a(0,1)*sol[1], RealApprox(REAL(1.0)).margin(Tolerance::Strict));
		REQUIRE_THAT(a(3,2)*sol[2] + a(3,3)*sol[3], RealApprox(REAL(4.0)).margin(Tolerance::Strict));
	}

	TEST_CASE("MatrixTridiag_Solve_in_place", "[tridiag][solve]")
	{
		TEST_PRECISION_INFO();
		TridiagonalMatrix<Real> a(3, { REAL(2.0), REAL(1.0),
		                               REAL(1.0), REAL(2.0), REAL(1.0),
		                                    REAL(1.0), REAL(2.0) });
		
		Vector<Real> rhs{ REAL(3.0), REAL(4.0), REAL(3.0) };
		Vector<Real> sol(3);
		
		a.Solve(rhs, sol);
		
		// Verify solution
		REQUIRE(sol.size() == 3);
		
		// The solution should be x = [1, 1, 1] for this system
		REQUIRE_THAT(a(0,0)*sol[0] + a(0,1)*sol[1], RealApprox(REAL(3.0)).margin(Tolerance::Strict));
	}

	///////////////////////          Print                             //////////////////////
	TEST_CASE("MatrixTridiag_Print_stream", "[tridiag][output]")
	{
		TEST_PRECISION_INFO();
		TridiagonalMatrix<Real> a(3, { REAL(2.0), REAL(1.0),
		                               REAL(1.0), REAL(3.0), REAL(1.0),
		                                    REAL(1.0), REAL(4.0) });
		
		std::stringstream ss;
		a.Print(ss, 5, 2);
		
		std::string output = ss.str();
		REQUIRE(!output.empty());
	}

	TEST_CASE("MatrixTridiag_Print_with_format", "[tridiag][output]")
	{
		TEST_PRECISION_INFO();
		TridiagonalMatrix<Real> a(3, { REAL(2.0), REAL(1.0),
		                               REAL(1.0), REAL(3.0), REAL(1.0),
		                                    REAL(1.0), REAL(4.0) });
		
		std::stringstream ss;
		MatrixPrintFormat fmt;
		fmt.showHeader = true;
		a.Print(ss, fmt);
		
		std::string output = ss.str();
		REQUIRE(output.find("Tridiagonal") != std::string::npos);
	}

	///////////////////////          Edge Cases                        //////////////////////
	TEST_CASE("MatrixTridiag_min_size_3", "[tridiag][edge_cases]")
	{
		TEST_PRECISION_INFO();
		// Minimum valid size is 3
		TridiagonalMatrix<Real> a(3, { REAL(2.0), REAL(1.0),
		                               REAL(1.0), REAL(2.0), REAL(1.0),
		                                    REAL(1.0), REAL(2.0) });
		
		REQUIRE(a.rows() == 3);
		REQUIRE(a.cols() == 3);
	}

	TEST_CASE("MatrixTridiag_zero_elements", "[tridiag][edge_cases]")
	{
		// Elements outside the band should be zero when accessing via const operator()
		const TridiagonalMatrix<Real> a(3, { REAL(2.0), REAL(1.0),
		                                     REAL(1.0), REAL(3.0), REAL(1.0),
		                                          REAL(1.0), REAL(4.0) });
		
		REQUIRE(a(0, 2) == REAL(0.0));  // Outside band: |0-2| > 1
		REQUIRE(a(2, 0) == REAL(0.0));  // Outside band: |2-0| > 1
	}

	TEST_CASE("MatrixTridiag_mutable_access", "[tridiag][access]")
	{
		TEST_PRECISION_INFO();
		TridiagonalMatrix<Real> a(3, { REAL(2.0), REAL(1.0),
		                               REAL(1.0), REAL(3.0), REAL(1.0),
		                                    REAL(1.0), REAL(4.0) });
		
		// Modify diagonal element
		a(1, 1) = REAL(10.0);
		REQUIRE(a(1, 1) == REAL(10.0));
		
		// Modify off-diagonal elements
		a(0, 1) = REAL(5.0);  // above diagonal
		REQUIRE(a(0, 1) == REAL(5.0));
		
		a(1, 0) = REAL(7.0);  // below diagonal
		REQUIRE(a(1, 0) == REAL(7.0));
	}

	TEST_CASE("MatrixTridiag_mutable_access_throws_outside_band", "[tridiag][access]")
	{
		TEST_PRECISION_INFO();
		TridiagonalMatrix<Real> a(3, { REAL(2.0), REAL(1.0),
		                               REAL(1.0), REAL(3.0), REAL(1.0),
		                                    REAL(1.0), REAL(4.0) });
		
		// Trying to modify element outside the band should throw
		REQUIRE_THROWS_AS(a(0, 2) = REAL(5.0), MatrixAccessBoundsError);
		REQUIRE_THROWS_AS(a(2, 0) = REAL(5.0), MatrixAccessBoundsError);
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

	TEST_CASE("MatrixTridiag_Solve_singular_detected", "[tridiag][solve]")
	{
		// Singular tridiagonal matrix (zero diagonal) should throw
		TridiagonalMatrix<Real> a(3,
			Vector<Real>{REAL(0.0), REAL(1.0), REAL(1.0)},  // below
			Vector<Real>{REAL(0.0), REAL(2.0), REAL(3.0)},  // diagonal (first element = 0)
			Vector<Real>{REAL(1.0), REAL(1.0), REAL(0.0)}); // above

		Vector<Real> rhs{REAL(1.0), REAL(2.0), REAL(3.0)};
		REQUIRE_THROWS_AS(a.Solve(rhs), SingularMatrixError);
	}

	TEST_CASE("MatrixTridiag_Solve_large_values_scaled_tolerance", "[tridiag][solve]")
	{
		// Large diagonal values - the scaled tolerance should still allow valid solves
		Real scale = REAL(1e12);
		TridiagonalMatrix<Real> a(3,
			Vector<Real>{REAL(0.0), -scale, -scale},           // below
			Vector<Real>{REAL(2.0)*scale, REAL(2.0)*scale, REAL(2.0)*scale},  // diagonal
			Vector<Real>{-scale, -scale, REAL(0.0)});          // above

		Vector<Real> rhs{scale, REAL(0.0), scale};
		Vector<Real> sol = a.Solve(rhs);
		REQUIRE(sol.size() == 3);
		// Solution should be [1, 1, 1] for this Laplacian-like system
		REQUIRE_THAT(sol[0], Catch::Matchers::WithinAbs(REAL(1.0), REAL(1e-3)));
		REQUIRE_THAT(sol[1], Catch::Matchers::WithinAbs(REAL(1.0), REAL(1e-3)));
		REQUIRE_THAT(sol[2], Catch::Matchers::WithinAbs(REAL(1.0), REAL(1e-3)));
	}
}


#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Vector.h"
#include "base/Matrix.h"
#include "base/MatrixSym.h"
#include "algorithms/MatrixAlg.h"
#endif

#include "../test_data/linear_alg_eq_systems_test_bed.h"

using namespace MML;
using namespace MML::Testing;
using namespace MML::MatrixAlg;

namespace MML::Tests::Algorithms::MatrixAlgTests
{
	/*********************************************************************/
	/*****               MATRIX SYMMETRY TESTS                       *****/
	/*********************************************************************/
	
	TEST_CASE("IsSymmetric_SymmetricMatrix", "[MatrixAlg][Symmetry]")
	{
			TEST_PRECISION_INFO();
		Matrix<Real> A{3, 3, {
			REAL(1.0), REAL(2.0), REAL(3.0),
			REAL(2.0), REAL(4.0), REAL(5.0),
			REAL(3.0), REAL(5.0), REAL(6.0)
		}};
		
		REQUIRE(IsSymmetric(A));
	}
	
	TEST_CASE("IsSymmetric_NonSymmetricMatrix", "[MatrixAlg][Symmetry]")
	{
			TEST_PRECISION_INFO();
		Matrix<Real> A{3, 3, {
			REAL(1.0), REAL(2.0), REAL(3.0),
			REAL(4.0), REAL(5.0), REAL(6.0),  // A[1][0] = 4 != 2 = A[0][1]
			REAL(7.0), REAL(8.0), REAL(9.0)
		}};
		
		REQUIRE_FALSE(IsSymmetric(A));
	}
	
	TEST_CASE("IsSkewSymmetric_True", "[MatrixAlg][Symmetry]")
	{
			TEST_PRECISION_INFO();
		Matrix<Real> A{3, 3, {
			 REAL(0.0),  REAL(2.0), -REAL(3.0),
			-REAL(2.0),  REAL(0.0),  REAL(5.0),
			 REAL(3.0), -REAL(5.0),  REAL(0.0)
		}};
		
		REQUIRE(IsSkewSymmetric(A));
	}
	
	TEST_CASE("IsSkewSymmetric_NonZeroDiagonal", "[MatrixAlg][Symmetry]")
	{
			TEST_PRECISION_INFO();
		Matrix<Real> A{3, 3, {
			 REAL(1.0),  REAL(2.0), -REAL(3.0),  // Non-zero diagonal
			-REAL(2.0),  REAL(0.0),  REAL(5.0),
			 REAL(3.0), -REAL(5.0),  REAL(0.0)
		}};
		
		REQUIRE_FALSE(IsSkewSymmetric(A));
	}
	
	/*********************************************************************/
	/*****               MATRIX DEFINITENESS TESTS                   *****/
	/*********************************************************************/
	
	TEST_CASE("IsPositiveDefinite_SPD_3x3", "[MatrixAlg][Definiteness]")
	{
			TEST_PRECISION_INFO();
		// SPD matrix from test bed
		auto sys = TestBeds::spd_3x3();
		
		REQUIRE(IsPositiveDefinite(sys._mat));
		REQUIRE(IsPositiveSemiDefinite(sys._mat));
		REQUIRE_FALSE(IsNegativeDefinite(sys._mat));
		REQUIRE_FALSE(IsIndefinite(sys._mat));
	}
	
	TEST_CASE("IsPositiveDefinite_SPD_CorrelationMatrix", "[MatrixAlg][Definiteness]")
	{
			TEST_PRECISION_INFO();
		auto sys = TestBeds::spd_4x4_correlation();
		
		REQUIRE(IsPositiveDefinite(sys._mat));
	}
	
	TEST_CASE("IsPositiveDefinite_SPD_MassMatrix", "[MatrixAlg][Definiteness]")
	{
			TEST_PRECISION_INFO();
		auto sys = TestBeds::spd_5x5_mass_matrix();
		
		REQUIRE(IsPositiveDefinite(sys._mat));
	}
	
	TEST_CASE("IsPositiveSemiDefinite_GraphLaplacian", "[MatrixAlg][Definiteness]")
	{
			TEST_PRECISION_INFO();
		// Graph Laplacian is PSD (has zero eigenvalue for constant vector)
		auto sys = TestBeds::spd_6x6_graph_laplacian();
		
		REQUIRE(IsPositiveSemiDefinite(sys._mat));
		// Graph Laplacians have a zero eigenvalue, so not strictly positive definite
		// (though with floating point, it might appear as very small positive)
	}
	
	TEST_CASE("IsNegativeDefinite_NegativeDiagonal", "[MatrixAlg][Definiteness]")
	{
			TEST_PRECISION_INFO();
		MatrixSym<Real> A{3, {
			-REAL(4.0),
			 REAL(0.0), -REAL(5.0),
			 REAL(0.0),  REAL(0.0), -REAL(6.0)
		}};
		
		REQUIRE(IsNegativeDefinite(A));
		REQUIRE(IsNegativeSemiDefinite(A));
		REQUIRE_FALSE(IsPositiveDefinite(A));
	}
	
	TEST_CASE("IsIndefinite_MixedEigenvalues", "[MatrixAlg][Definiteness]")
	{
			TEST_PRECISION_INFO();
		// Indefinite matrix: has both positive and negative eigenvalues
		MatrixSym<Real> A{3, {
			 REAL(2.0),
			 REAL(0.0), -REAL(1.0),
			 REAL(0.0),  REAL(0.0),  REAL(3.0)
		}};
		
		REQUIRE(IsIndefinite(A));
		REQUIRE_FALSE(IsPositiveDefinite(A));
		REQUIRE_FALSE(IsNegativeDefinite(A));
	}
	
	TEST_CASE("ClassifyDefiniteness_All", "[MatrixAlg][Definiteness]")
	{
			TEST_PRECISION_INFO();
		// Positive definite
		MatrixSym<Real> pd{2, {REAL(2.0), REAL(0.5), REAL(2.0)}};
		REQUIRE(ClassifyDefiniteness(pd) == Definiteness::PositiveDefinite);
		
		// Negative definite
		MatrixSym<Real> nd{2, {-REAL(2.0), REAL(0.5), -REAL(2.0)}};
		REQUIRE(ClassifyDefiniteness(nd) == Definiteness::NegativeDefinite);
		
		// Positive semi-definite (singular)
		MatrixSym<Real> psd{2, {REAL(1.0), REAL(1.0), REAL(1.0)}};  // rank 1 matrix
		REQUIRE(ClassifyDefiniteness(psd) == Definiteness::PositiveSemiDefinite);
		
		// Indefinite
		MatrixSym<Real> indef{2, {REAL(1.0), REAL(0.0), -REAL(1.0)}};
		REQUIRE(ClassifyDefiniteness(indef) == Definiteness::Indefinite);
	}
	
	TEST_CASE("Definiteness_GeneralMatrix", "[MatrixAlg][Definiteness]")
	{
			TEST_PRECISION_INFO();
		// Test with general (non-symmetric) matrix - uses symmetric part
		Matrix<Real> A{3, 3, {
			REAL(4.0), REAL(1.0), REAL(0.0),
			REAL(1.0), REAL(4.0), REAL(1.0),
			REAL(0.0), REAL(1.0), REAL(4.0)
		}};
		
		REQUIRE(IsPositiveDefinite(A));
	}
	
	/*********************************************************************/
	/*****                SVD-BASED UTILITIES                        *****/
	/*********************************************************************/
	
	TEST_CASE("SingularValues_DiagonalMatrix", "[MatrixAlg][SVD]")
	{
			TEST_PRECISION_INFO();
		Matrix<Real> A{3, 3, {
			REAL(5.0), REAL(0.0), REAL(0.0),
			REAL(0.0), REAL(3.0), REAL(0.0),
			REAL(0.0), REAL(0.0), REAL(1.0)
		}};
		
		Vector<Real> sv = SingularValues(A);
		
		// Singular values should be sorted in descending order
		REQUIRE(sv.size() == 3);
		REQUIRE(std::abs(sv[0] - REAL(5.0)) < 1e-10);
		REQUIRE(std::abs(sv[1] - REAL(3.0)) < 1e-10);
		REQUIRE(std::abs(sv[2] - REAL(1.0)) < 1e-10);
	}
	
	TEST_CASE("Rank_FullRank", "[MatrixAlg][SVD]")
	{
			TEST_PRECISION_INFO();
		auto sys = TestBeds::diag_dominant_4x4();
		
		int rank = Rank(sys._mat);
		REQUIRE(rank == 4);
	}
	
	TEST_CASE("Rank_RankDeficient", "[MatrixAlg][SVD]")
	{
			TEST_PRECISION_INFO();
		// Rank 2 matrix (3x3)
		Matrix<Real> A{3, 3, {
			REAL(1.0), REAL(2.0), REAL(3.0),
			REAL(2.0), REAL(4.0), REAL(6.0),   // Row 2 = 2 * Row 1
			REAL(3.0), REAL(6.0), REAL(9.0)    // Row 3 = 3 * Row 1
		}};
		
		int rank = Rank(A);
		REQUIRE(rank == 1);
	}
	
	TEST_CASE("Nullity_FullRank", "[MatrixAlg][SVD]")
	{
			TEST_PRECISION_INFO();
		auto sys = TestBeds::diag_dominant_4x4();
		
		int nullity = Nullity(sys._mat);
		REQUIRE(nullity == 0);
	}
	
	TEST_CASE("Nullity_RankDeficient", "[MatrixAlg][SVD]")
	{
			TEST_PRECISION_INFO();
		// Rank 1 matrix → nullity = 2
		Matrix<Real> A{3, 3, {
			REAL(1.0), REAL(2.0), REAL(3.0),
			REAL(2.0), REAL(4.0), REAL(6.0),
			REAL(3.0), REAL(6.0), REAL(9.0)
		}};
		
		int nullity = Nullity(A);
		REQUIRE(nullity == 2);
	}
	
	TEST_CASE("ConditionNumber_WellConditioned", "[MatrixAlg][SVD]")
	{
			TEST_PRECISION_INFO();
		// Diagonal matrix with condition number = 5/1 = 5
		Matrix<Real> A{3, 3, {
			REAL(5.0), REAL(0.0), REAL(0.0),
			REAL(0.0), REAL(3.0), REAL(0.0),
			REAL(0.0), REAL(0.0), REAL(1.0)
		}};
		
		Real cond = ConditionNumber(A);
		REQUIRE(std::abs(cond - REAL(5.0)) < 1e-10);
	}
	
	TEST_CASE("ConditionNumber_Identity", "[MatrixAlg][SVD]")
	{
			TEST_PRECISION_INFO();
		auto I = Matrix<Real>::GetUnitMatrix(5);
		
		Real cond = ConditionNumber(I);
		REQUIRE(std::abs(cond - REAL(1.0)) < 1e-10);
	}
	
	TEST_CASE("NullSpace_FullRank", "[MatrixAlg][SVD]")
	{
			TEST_PRECISION_INFO();
		auto sys = TestBeds::diag_dominant_4x4();
		
		Matrix<Real> ns = NullSpace(sys._mat);
		
		// Full rank matrix has trivial null space
		REQUIRE(ns.ColNum() == 0);
	}
	
	TEST_CASE("NullSpace_RankDeficient", "[MatrixAlg][SVD]")
	{
			TEST_PRECISION_INFO();
		// Rank 1 matrix
		Matrix<Real> A{3, 3, {
			REAL(1.0), REAL(2.0), REAL(3.0),
			REAL(2.0), REAL(4.0), REAL(6.0),
			REAL(3.0), REAL(6.0), REAL(9.0)
		}};
		
		Matrix<Real> ns = NullSpace(A);
		
		// Should have 2-dimensional null space
		REQUIRE(ns.ColNum() == 2);
		REQUIRE(ns.RowNum() == 3);
		
		// Verify: A * x should be zero for each null space vector
		for (int col = 0; col < ns.ColNum(); col++)
		{
			Vector<Real> x(3);
			for (int i = 0; i < 3; i++)
				x[i] = ns(i, col);
			
			Vector<Real> Ax = A * x;
			REQUIRE(Ax.NormL2() < 1e-10);
		}
	}
	
	TEST_CASE("Kernel_SameAsNullSpace", "[MatrixAlg][SVD]")
	{
			TEST_PRECISION_INFO();
		Matrix<Real> A{3, 3, {
			REAL(1.0), REAL(2.0), REAL(3.0),
			REAL(2.0), REAL(4.0), REAL(6.0),
			REAL(3.0), REAL(6.0), REAL(9.0)
		}};
		
		Matrix<Real> ns = NullSpace(A);
		Matrix<Real> ker = Kernel(A);
		
		REQUIRE(ns.RowNum() == ker.RowNum());
		REQUIRE(ns.ColNum() == ker.ColNum());
	}
	
	TEST_CASE("ColumnSpace_FullRank", "[MatrixAlg][SVD]")
	{
			TEST_PRECISION_INFO();
		auto sys = TestBeds::diag_dominant_4x4();
		
		Matrix<Real> cs = ColumnSpace(sys._mat);
		
		// Full rank 4x4 → column space is R^4
		REQUIRE(cs.ColNum() == 4);
		REQUIRE(cs.RowNum() == 4);
	}
	
	TEST_CASE("ColumnSpace_RankDeficient", "[MatrixAlg][SVD]")
	{
			TEST_PRECISION_INFO();
		// Rank 1 matrix
		Matrix<Real> A{3, 3, {
			REAL(1.0), REAL(2.0), REAL(3.0),
			REAL(2.0), REAL(4.0), REAL(6.0),
			REAL(3.0), REAL(6.0), REAL(9.0)
		}};
		
		Matrix<Real> cs = ColumnSpace(A);
		
		// Column space is 1-dimensional
		REQUIRE(cs.ColNum() == 1);
		REQUIRE(cs.RowNum() == 3);
	}
	
	TEST_CASE("FundamentalSubspaces_RankNullityTheorem", "[MatrixAlg][SVD]")
	{
			TEST_PRECISION_INFO();
		// For m×n matrix: rank + nullity = n
		Matrix<Real> A{4, 3, {
			REAL(1.0), REAL(2.0), REAL(3.0),
			REAL(2.0), REAL(4.0), REAL(6.0),
			REAL(3.0), REAL(6.0), REAL(9.0),
			REAL(4.0), REAL(8.0), REAL(12.0)
		}};
		
		auto fs = ComputeFundamentalSubspaces(A);
		
		int m = A.RowNum();  // 4
		int n = A.ColNum();  // 3
		
		// Rank-nullity theorem: rank + nullity = n
		REQUIRE(fs.rank + fs.nullSpace.ColNum() == n);
		
		// dim(col space) = rank
		REQUIRE(fs.columnSpace.ColNum() == fs.rank);
		
		// dim(row space) = rank
		REQUIRE(fs.rowSpace.ColNum() == fs.rank);
		
		// dim(left null space) = m - rank
		REQUIRE(fs.leftNullSpace.ColNum() == m - fs.rank);
	}
	
	TEST_CASE("PseudoInverse_FullRankSquare", "[MatrixAlg][SVD]")
	{
			TEST_PRECISION_INFO();
		// For invertible matrix, pseudoinverse = inverse
		Matrix<Real> A{3, 3, {
			REAL(2.0), REAL(0.0), REAL(0.0),
			REAL(0.0), REAL(3.0), REAL(0.0),
			REAL(0.0), REAL(0.0), REAL(4.0)
		}};
		
		Matrix<Real> Ainv = PseudoInverse(A);
		
		// A * A⁺ should be identity
		Matrix<Real> I = A * Ainv;
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
			{
				Real expected = (i == j) ? REAL(1.0) : REAL(0.0);
				REQUIRE(std::abs(I(i, j) - expected) < 1e-10);
			}
	}
	
	TEST_CASE("PseudoInverse_LeastSquares", "[MatrixAlg][SVD]")
	{
			TEST_PRECISION_INFO();
		// Overdetermined system: pseudoinverse gives least squares solution
		Matrix<Real> A{4, 2, {
			REAL(1.0), REAL(0.0),
			REAL(0.0), REAL(1.0),
			REAL(1.0), REAL(0.0),
			REAL(0.0), REAL(1.0)
		}};
		
		Vector<Real> b{REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0)};
		
		Matrix<Real> Aplus = PseudoInverse(A);
		
		// x = A⁺ * b is least squares solution
		Vector<Real> x(2);
		for (int i = 0; i < 2; i++)
		{
			x[i] = REAL(0.0);
			for (int j = 0; j < 4; j++)
				x[i] += Aplus(i, j) * b[j];
		}
		
		// Expected: x = [2, 3] (average of (1,3) and (2,4))
		REQUIRE(std::abs(x[0] - REAL(2.0)) < 1e-10);
		REQUIRE(std::abs(x[1] - REAL(3.0)) < 1e-10);
	}
	
	TEST_CASE("ComputeSVD_Basic", "[MatrixAlg][SVD]")
	{
			TEST_PRECISION_INFO();
		Matrix<Real> A{3, 3, {
			REAL(1.0), REAL(0.0), REAL(0.0),
			REAL(0.0), REAL(2.0), REAL(0.0),
			REAL(0.0), REAL(0.0), REAL(3.0)
		}};
		
		auto svd = ComputeSVD(A);
		
		REQUIRE(svd.rank == 3);
		REQUIRE(std::abs(svd.conditionNumber - REAL(3.0)) < 1e-10);
		REQUIRE(svd.singularValues.size() == 3);
	}

} // namespace MML::Tests::Algorithms::MatrixAlgTests

#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Vector.h"
#include "base/Matrix.h"
#include "core/LinAlgEqSolvers.h"
#endif

using namespace MML;
using namespace MML::Testing;

namespace MML::Tests::Core::SVDTests
{
	// Helper function to verify orthogonality
	bool IsOrthogonal(const Matrix<Real>& mat, Real tolerance = 1e-10)
	{
		int n = mat.ColNum();
		Matrix<Real> prod = mat.GetTranspose() * mat;
		
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				Real expected = (i == j) ? REAL(1.0) : REAL(0.0);
				if (std::abs(prod[i][j] - expected) > tolerance)
					return false;
			}
		}
		return true;
	}

	// Helper to verify A = U * W * V^T
	bool VerifyDecomposition(const Matrix<Real>& A, const Matrix<Real>& U, 
	                         const Vector<Real>& W, const Matrix<Real>& V, 
	                         Real tolerance = 1e-10)
	{
		int m = A.RowNum();
		int n = A.ColNum();
		
		// Reconstruct: A_reconstructed = U * W * V^T
		Matrix<Real> reconstructed(m, n);
		for (int i = 0; i < m; i++)
		{
			for (int j = 0; j < n; j++)
			{
				Real sum = REAL(0.0);
				for (int k = 0; k < n; k++)
					sum += U[i][k] * W[k] * V[j][k];
				reconstructed[i][j] = sum;
			}
		}
		
		// Compare with original
		for (int i = 0; i < m; i++)
		{
			for (int j = 0; j < n; j++)
			{
				if (std::abs(A[i][j] - reconstructed[i][j]) > tolerance)
					return false;
			}
		}
		return true;
	}

	///////////////// PHASE 1: BASIC DECOMPOSITION TESTS /////////////////

	TEST_CASE("SVD_Identity_3x3", "[SVD][Basic]")
	{
			TEST_PRECISION_INFO();
		Matrix<Real> A(3, 3, {
			REAL(1.0), REAL(0.0), REAL(0.0),
			REAL(0.0), REAL(1.0), REAL(0.0),
			REAL(0.0), REAL(0.0), REAL(1.0)
		});

		SVDecompositionSolver svd(A);
		
		auto U = svd.getU();
		auto W = svd.getW();
		auto V = svd.getV();

		// Singular values should all be REAL(1.0)
		REQUIRE(std::abs(W[0] - REAL(1.0)) < 1e-10);
		REQUIRE(std::abs(W[1] - REAL(1.0)) < 1e-10);
		REQUIRE(std::abs(W[2] - REAL(1.0)) < 1e-10);

		// Verify reconstruction
		REQUIRE(VerifyDecomposition(A, U, W, V, 1e-10));
		
		// U and V should be orthogonal
		REQUIRE(IsOrthogonal(U));
		REQUIRE(IsOrthogonal(V));
	}

	TEST_CASE("SVD_Diagonal_3x3", "[SVD][Basic]")
	{
			TEST_PRECISION_INFO();
		Matrix<Real> A(3, 3, {
			REAL(3.0), REAL(0.0), REAL(0.0),
			REAL(0.0), REAL(2.0), REAL(0.0),
			REAL(0.0), REAL(0.0), REAL(1.0)
		});

		SVDecompositionSolver svd(A);
		
		auto U = svd.getU();
		auto W = svd.getW();
		auto V = svd.getV();

		// Singular values should be 3, 2, 1 (in descending order)
		REQUIRE(std::abs(W[0] - REAL(3.0)) < 1e-10);
		REQUIRE(std::abs(W[1] - REAL(2.0)) < 1e-10);
		REQUIRE(std::abs(W[2] - REAL(1.0)) < 1e-10);

		// Verify reconstruction
		REQUIRE(VerifyDecomposition(A, U, W, V, 1e-10));
		
		// U and V should be orthogonal
		REQUIRE(IsOrthogonal(U));
		REQUIRE(IsOrthogonal(V));
	}

	TEST_CASE("SVD_Simple_2x2", "[SVD][Basic]")
	{
			TEST_PRECISION_INFO();
		Matrix<Real> A(2, 2, {
			REAL(4.0), REAL(0.0),
			REAL(3.0), -REAL(5.0)
		});

		SVDecompositionSolver svd(A);
		
		auto U = svd.getU();
		auto W = svd.getW();
		auto V = svd.getV();

		// Singular values should be positive and in descending order
		REQUIRE(W[0] > 0);
		REQUIRE(W[1] > 0);
		REQUIRE(W[0] >= W[1]);

		// Verify reconstruction
		REQUIRE(VerifyDecomposition(A, U, W, V, 1e-9));
		
		// U and V should be orthogonal
		REQUIRE(IsOrthogonal(U, 1e-9));
		REQUIRE(IsOrthogonal(V, 1e-9));
	}

	TEST_CASE("SVD_Rectangular_4x3", "[SVD][Basic]")
	{
			TEST_PRECISION_INFO();
		// Overdetermined system (m > n)
		Matrix<Real> A(4, 3, {
			REAL(1.0), REAL(2.0), REAL(3.0),
			REAL(4.0), REAL(5.0), REAL(6.0),
			REAL(7.0), REAL(8.0), REAL(9.0),
			REAL(10.0), REAL(11.0), REAL(12.0)
		});

		SVDecompositionSolver svd(A);
		
		auto U = svd.getU();
		auto W = svd.getW();
		auto V = svd.getV();

		// Check dimensions
		REQUIRE(U.RowNum() == 4);
		REQUIRE(U.ColNum() == 3);
		REQUIRE(V.RowNum() == 3);
		REQUIRE(V.ColNum() == 3);
		REQUIRE(W.size() == 3);

		// Singular values should be non-negative and descending
		for (int i = 0; i < 3; i++)
			REQUIRE(W[i] >= 0);
		REQUIRE(W[0] >= W[1]);
		REQUIRE(W[1] >= W[2]);

		// Verify reconstruction
		REQUIRE(VerifyDecomposition(A, U, W, V, 1e-8));
		
		// V should be orthogonal (U^T*U = I for overdetermined)
		REQUIRE(IsOrthogonal(V, 1e-9));
	}

	TEST_CASE("SVD_Rectangular_3x4", "[SVD][Basic]")
	{
			TEST_PRECISION_INFO();
		// Underdetermined system (m < n)
		Matrix<Real> A(3, 4, {
			REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0),
			REAL(5.0), REAL(6.0), REAL(7.0), REAL(8.0),
			REAL(9.0), REAL(10.0), REAL(11.0), REAL(12.0)
		});

		SVDecompositionSolver svd(A);
		
		auto U = svd.getU();
		auto W = svd.getW();
		auto V = svd.getV();

		// Check dimensions
		REQUIRE(U.RowNum() == 3);
		REQUIRE(U.ColNum() == 4);
		REQUIRE(V.RowNum() == 4);
		REQUIRE(V.ColNum() == 4);
		REQUIRE(W.size() == 4);

		// Singular values should be non-negative and descending
		for (int i = 0; i < 4; i++)
			REQUIRE(W[i] >= 0);
		REQUIRE(W[0] >= W[1]);
		REQUIRE(W[1] >= W[2]);
		REQUIRE(W[2] >= W[3]);

		// Verify reconstruction
		REQUIRE(VerifyDecomposition(A, U, W, V, 1e-8));
		
		// U and V should be orthogonal
		REQUIRE(IsOrthogonal(V, 1e-9));
	}

	TEST_CASE("SVD_Symmetric_3x3", "[SVD][Symmetric]")
	{
			TEST_PRECISION_INFO();
		// Symmetric matrix: eigenvalues = singular values
		Matrix<Real> A(3, 3, {
			REAL(4.0), REAL(1.0), REAL(2.0),
			REAL(1.0), REAL(3.0), REAL(1.0),
			REAL(2.0), REAL(1.0), REAL(5.0)
		});

		SVDecompositionSolver svd(A);
		
		auto U = svd.getU();
		auto W = svd.getW();
		auto V = svd.getV();

		// Verify reconstruction
		REQUIRE(VerifyDecomposition(A, U, W, V, 1e-9));
		
		// Both U and V should be orthogonal
		REQUIRE(IsOrthogonal(U, 1e-9));
		REQUIRE(IsOrthogonal(V, 1e-9));

		// For symmetric matrices, U and V should be similar (up to signs)
		// This is a weaker test - just verify they're both orthogonal
	}

	TEST_CASE("SVD_Condition_Number", "[SVD][Properties]")
	{
			TEST_PRECISION_INFO();
		// Well-conditioned matrix
		Matrix<Real> A(3, 3, {
			REAL(2.0), REAL(0.0), REAL(0.0),
			REAL(0.0), REAL(2.0), REAL(0.0),
			REAL(0.0), REAL(0.0), REAL(2.0)
		});

		SVDecompositionSolver svd(A);
		
		// Condition number should be REAL(1.0) (all singular values equal)
		Real cond = svd.inv_condition();
		REQUIRE(std::abs(cond - REAL(1.0)) < 1e-10);

		// Ill-conditioned matrix
		Matrix<Real> B(3, 3, {
			REAL(10.0), REAL(0.0), REAL(0.0),
			REAL(0.0), REAL(1.0), REAL(0.0),
			REAL(0.0), REAL(0.0), REAL(0.1)
		});

		SVDecompositionSolver svd2(B);
		Real cond2 = svd2.inv_condition();
		
		// Condition should be REAL(0.1)/10 = REAL(0.01)
		REQUIRE(std::abs(cond2 - REAL(0.01)) < 1e-10);
	}

	///////////////// PHASE 3: SOLVE TESTS /////////////////

	TEST_CASE("SVD_Solve_Square_Simple", "[SVD][Solve]")
	{
			TEST_PRECISION_INFO();
		// Simple 3x3 system
		Matrix<Real> A(3, 3, {
			REAL(1.0), REAL(2.0), REAL(3.0),
			REAL(4.0), REAL(5.0), REAL(6.0),
			REAL(7.0), REAL(8.0), REAL(10.0)  // Modified to make full rank
		});

		Vector<Real> b(3);
		b[0] = REAL(6.0); b[1] = REAL(15.0); b[2] = REAL(26.0);  // A * [1, 1, 1]^T = b (approximately)

		SVDecompositionSolver svd(A);
		Vector<Real> x = svd.Solve(b);

		// Check solution by computing A*x
		Vector<Real> result(3);
		for (int i = 0; i < 3; i++)
		{
			Real sum = REAL(0.0);
			for (int j = 0; j < 3; j++)
				sum += A[i][j] * x[j];
			result[i] = sum;
		}

		// Verify A*x ≈ b
		for (int i = 0; i < 3; i++)
			REQUIRE(std::abs(result[i] - b[i]) < 1e-8);
	}

	TEST_CASE("SVD_Solve_Overdetermined", "[SVD][Solve]")
	{
			TEST_PRECISION_INFO();
		// Overdetermined system (more equations than unknowns) - least squares solution
		Matrix<Real> A(4, 3, {
			REAL(1.0), REAL(0.0), REAL(0.0),
			REAL(0.0), REAL(1.0), REAL(0.0),
			REAL(0.0), REAL(0.0), REAL(1.0),
			REAL(1.0), REAL(1.0), REAL(1.0)
		});

		Vector<Real> b(4);
		b[0] = REAL(1.0); b[1] = REAL(2.0); b[2] = REAL(3.0); b[3] = REAL(5.0);  // Inconsistent system

		SVDecompositionSolver svd(A);
		Vector<Real> x = svd.Solve(b);

		// Check that x minimizes ||A*x - b||
		Vector<Real> residual(4);
		for (int i = 0; i < 4; i++)
		{
			Real sum = REAL(0.0);
			for (int j = 0; j < 3; j++)
				sum += A[i][j] * x[j];
			residual[i] = sum - b[i];
		}

		// Residual should be small (least squares solution)
		Real residualNorm = REAL(0.0);
		for (int i = 0; i < 4; i++)
			residualNorm += residual[i] * residual[i];
		residualNorm = std::sqrt(residualNorm);

		REQUIRE(residualNorm < REAL(2.0));  // Should find reasonable least squares solution
	}

	TEST_CASE("SVD_Solve_Identity", "[SVD][Solve]")
	{
			TEST_PRECISION_INFO();
		Matrix<Real> A(3, 3, {
			REAL(1.0), REAL(0.0), REAL(0.0),
			REAL(0.0), REAL(1.0), REAL(0.0),
			REAL(0.0), REAL(0.0), REAL(1.0)
		});

		Vector<Real> b(3);
		b[0] = REAL(5.0); b[1] = REAL(7.0); b[2] = REAL(9.0);

		SVDecompositionSolver svd(A);
		Vector<Real> x = svd.Solve(b);

		// For identity matrix, solution should equal b
		REQUIRE(std::abs(x[0] - REAL(5.0)) < 1e-10);
		REQUIRE(std::abs(x[1] - REAL(7.0)) < 1e-10);
		REQUIRE(std::abs(x[2] - REAL(9.0)) < 1e-10);
	}

	TEST_CASE("SVD_Solve_Matrix_RHS", "[SVD][Solve]")
	{
			TEST_PRECISION_INFO();
		// Solve multiple right-hand sides at once
		Matrix<Real> A(3, 3, {
			REAL(2.0), REAL(0.0), REAL(0.0),
			REAL(0.0), REAL(3.0), REAL(0.0),
			REAL(0.0), REAL(0.0), REAL(4.0)
		});

		Matrix<Real> B(3, 2, {
			REAL(4.0), REAL(6.0),
			REAL(9.0), REAL(12.0),
			REAL(16.0), REAL(20.0)
		});

		SVDecompositionSolver svd(A);
		Matrix<Real> X(3, 2);
		svd.Solve(B, X);

		// Expected solutions: X = [[2,3], [3,4], [4,5]]
		REQUIRE(std::abs(X[0][0] - REAL(2.0)) < 1e-10);
		REQUIRE(std::abs(X[0][1] - REAL(3.0)) < 1e-10);
		REQUIRE(std::abs(X[1][0] - REAL(3.0)) < 1e-10);
		REQUIRE(std::abs(X[1][1] - REAL(4.0)) < 1e-10);
		REQUIRE(std::abs(X[2][0] - REAL(4.0)) < 1e-10);
		REQUIRE(std::abs(X[2][1] - REAL(5.0)) < 1e-10);
	}

	///////////////// PHASE 4: PROPERTY TESTS /////////////////

	TEST_CASE("SVD_Rank_FullRank", "[SVD][Properties]")
	{
			TEST_PRECISION_INFO();
		Matrix<Real> A(3, 3, {
			REAL(1.0), REAL(0.0), REAL(0.0),
			REAL(0.0), REAL(1.0), REAL(0.0),
			REAL(0.0), REAL(0.0), REAL(1.0)
		});

		SVDecompositionSolver svd(A);
		REQUIRE(svd.Rank() == 3);
		REQUIRE(svd.Nullity() == 0);
	}

	TEST_CASE("SVD_Rank_Deficient", "[SVD][Properties]")
	{
			TEST_PRECISION_INFO();
		// Rank-deficient matrix (rank 2)
		Matrix<Real> A(3, 3, {
			REAL(1.0), REAL(2.0), REAL(3.0),
			REAL(2.0), REAL(4.0), REAL(6.0),
			REAL(3.0), REAL(6.0), REAL(9.0)
		});

		SVDecompositionSolver svd(A);
		
		int rank = svd.Rank();
		int nullity = svd.Nullity();

		REQUIRE(rank == 1);  // Only one independent row
		REQUIRE(nullity == 2);
		REQUIRE(rank + nullity == 3);
	}

	TEST_CASE("SVD_Range", "[SVD][Properties]")
	{
			TEST_PRECISION_INFO();
		Matrix<Real> A(3, 3, {
			REAL(1.0), REAL(0.0), REAL(0.0),
			REAL(0.0), REAL(2.0), REAL(0.0),
			REAL(0.0), REAL(0.0), REAL(0.0)  // Zero row - rank 2
		});

		SVDecompositionSolver svd(A);
		
		int rank = svd.Rank();
		REQUIRE(rank == 2);

		Matrix<Real> range = svd.Range();
		
		// Range should have correct dimensions
		REQUIRE(range.RowNum() == 3);
		REQUIRE(range.ColNum() == 2);

		// Columns of range should be orthonormal
		for (int i = 0; i < 2; i++)
		{
			Real norm = REAL(0.0);
			for (int j = 0; j < 3; j++)
				norm += range[j][i] * range[j][i];
			REQUIRE(std::abs(norm - REAL(1.0)) < 1e-9);
		}
	}

	TEST_CASE("SVD_Nullspace", "[SVD][Properties]")
	{
			TEST_PRECISION_INFO();
		// Matrix with known nullspace
		Matrix<Real> A(3, 3, {
			REAL(1.0), REAL(2.0), REAL(3.0),
			REAL(2.0), REAL(4.0), REAL(6.0),
			REAL(3.0), REAL(6.0), REAL(9.0)
		});

		SVDecompositionSolver svd(A);
		
		int nullity = svd.Nullity();
		REQUIRE(nullity == 2);

		Matrix<Real> nullspace = svd.Nullspace();
		
		// Nullspace should have correct dimensions
		REQUIRE(nullspace.RowNum() == 3);
		REQUIRE(nullspace.ColNum() == 2);

		// Verify A * v = 0 for each nullspace vector v
		for (int col = 0; col < 2; col++)
		{
			Vector<Real> v(3);
			for (int i = 0; i < 3; i++)
				v[i] = nullspace[i][col];

			// Compute A*v
			Vector<Real> result(3);
			for (int i = 0; i < 3; i++)
			{
				Real sum = REAL(0.0);
				for (int j = 0; j < 3; j++)
					sum += A[i][j] * v[j];
				result[i] = sum;
			}

			// Should be zero vector
			for (int i = 0; i < 3; i++)
				REQUIRE(std::abs(result[i]) < 1e-8);
		}
	}

	TEST_CASE("SVD_RankNullity_Theorem", "[SVD][Properties]")
	{
			TEST_PRECISION_INFO();
		// Verify rank-nullity theorem: rank + nullity = n for various matrices

		// Full rank
		Matrix<Real> A1(3, 3, {
			REAL(1.0), REAL(2.0), REAL(3.0),
			REAL(4.0), REAL(5.0), REAL(6.0),
			REAL(7.0), REAL(8.0), REAL(10.0)
		});
		SVDecompositionSolver svd1(A1);
		REQUIRE(svd1.Rank() + svd1.Nullity() == 3);

		// Rank deficient
		Matrix<Real> A2(4, 3, {
			REAL(1.0), REAL(0.0), REAL(0.0),
			REAL(0.0), REAL(1.0), REAL(0.0),
			REAL(0.0), REAL(0.0), REAL(0.0),
			REAL(0.0), REAL(0.0), REAL(0.0)
		});
		SVDecompositionSolver svd2(A2);
		REQUIRE(svd2.Rank() + svd2.Nullity() == 3);
	}
}

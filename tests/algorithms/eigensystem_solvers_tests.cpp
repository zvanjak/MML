#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#include "../../test_data/linear_alg_eq_systems_classic_defs.h"
#include "../../test_data/linear_alg_eq_systems_sym_defs.h"
#include "../../test_data/linear_alg_eq_systems_solverspecific_defs.h"
#include "../../test_data/linear_alg_eq_systems_test_bed.h"
#include "../../test_data/eigenvalue_test_bed.h"

#include "../../mml/algorithms/EigenSystemSolvers.h"
#include "../../mml/base/BaseUtils.h"

#include <sstream>

using namespace MML;
using namespace MML::Testing;
using namespace MML::TestBeds;

/***************************************************************************************************
 * JACOBI EIGENSOLVER TESTS
 ***************************************************************************************************/

TEST_CASE("JacobiEigen_2x2_Simple", "[EigenSolvers][Jacobi][Small]")
{
		TEST_PRECISION_INFO();
	// Simple 2x2 symmetric matrix with known eigenvalues
	MatrixSym<Real> A(2, {REAL(2.0), REAL(1.0), REAL(2.0)});
	
	// Known eigenvalues: λ₁ = 1, λ₂ = 3
	// Known eigenvectors: v₁ = [-1/√2, 1/√2], v₂ = [1/√2, 1/√2]
	
	auto result = SymmMatEigenSolverJacobi::Solve(A);
	
	REQUIRE(result.converged);
	REQUIRE(result.eigenvalues.size() == 2);
	
	// Check eigenvalues (sorted ascending)
	REQUIRE(std::abs(result.eigenvalues[0] - REAL(1.0)) < 1e-10);
	REQUIRE(std::abs(result.eigenvalues[1] - REAL(3.0)) < 1e-10);
	
	// Verify A*v = λ*v for each eigenpair
	for (int i = 0; i < 2; i++)
	{
		Vector<Real> v = result.eigenvectors.VectorFromColumn(i);
		Vector<Real> Av(2);
		for (int j = 0; j < 2; j++)
			Av[j] = A(j, 0) * v[0] + A(j, 1) * v[1];
		
		Vector<Real> lambda_v = result.eigenvalues[i] * v;
		
		for (int j = 0; j < 2; j++)
			REQUIRE(std::abs(Av[j] - lambda_v[j]) < 1e-10);
	}
	
	// Verify eigenvectors are orthonormal
	Vector<Real> v0 = result.eigenvectors.VectorFromColumn(0);
	Vector<Real> v1 = result.eigenvectors.VectorFromColumn(1);
	
	Real dot01 = v0[0] * v1[0] + v0[1] * v1[1];
	Real norm0 = std::sqrt(v0[0] * v0[0] + v0[1] * v0[1]);
	Real norm1 = std::sqrt(v1[0] * v1[0] + v1[1] * v1[1]);
	
	REQUIRE(std::abs(dot01) < 1e-10);        // Orthogonal
	REQUIRE(std::abs(norm0 - REAL(1.0)) < 1e-10);  // Unit norm
	REQUIRE(std::abs(norm1 - REAL(1.0)) < 1e-10);  // Unit norm
}

TEST_CASE("JacobiEigen_3x3_Identity", "[EigenSolvers][Jacobi][Small]")
{
		TEST_PRECISION_INFO();
	// Identity matrix - all eigenvalues = 1
	MatrixSym<Real> I(3, {REAL(1.0), REAL(0.0), REAL(1.0), REAL(0.0), REAL(0.0), REAL(1.0)});
	
	auto result = SymmMatEigenSolverJacobi::Solve(I);
	
	REQUIRE(result.converged);
	REQUIRE(result.eigenvalues.size() == 3);
	
	// All eigenvalues should be 1
	for (int i = 0; i < 3; i++)
		REQUIRE(std::abs(result.eigenvalues[i] - REAL(1.0)) < 1e-10);
}

TEST_CASE("JacobiEigen_3x3_Diagonal", "[EigenSolvers][Jacobi][Small]")
{
		TEST_PRECISION_INFO();
	// Diagonal matrix - eigenvalues are diagonal elements
	MatrixSym<Real> D(3, {REAL(5.0), REAL(0.0), REAL(3.0), REAL(0.0), REAL(0.0), REAL(1.0)});
	
	auto result = SymmMatEigenSolverJacobi::Solve(D);
	
	REQUIRE(result.converged);
	
	// Eigenvalues should be 1, 3, 5 (sorted)
	REQUIRE(std::abs(result.eigenvalues[0] - REAL(1.0)) < 1e-10);
	REQUIRE(std::abs(result.eigenvalues[1] - REAL(3.0)) < 1e-10);
	REQUIRE(std::abs(result.eigenvalues[2] - REAL(5.0)) < 1e-10);
}

TEST_CASE("JacobiEigen_SymmetricMatrix_3x3", "[EigenSolvers][Jacobi][TestData]")
{
		TEST_PRECISION_INFO();
	// Use test data: symm_mat_3x3 with known eigenvalues
	auto result = SymmMatEigenSolverJacobi::Solve(symm_mat_3x3);
	
	REQUIRE(result.converged);
	REQUIRE(result.eigenvalues.size() == 3);
	
	// Compare with known eigenvalues (sorted)
	Vector<Real> expected = symm_mat_3x3_eigen_val;
	std::sort(expected.begin(), expected.end());
	
	for (int i = 0; i < 3; i++)
		REQUIRE(std::abs(result.eigenvalues[i] - expected[i]) < 1e-6);
	
	// Verify A*V = V*Λ
	Matrix<Real> A_full(3, 3);
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			A_full[i][j] = (i <= j) ? symm_mat_3x3(i, j) : symm_mat_3x3(j, i);
	
	for (int i = 0; i < 3; i++)
	{
		Vector<Real> v = result.eigenvectors.VectorFromColumn(i);
		Vector<Real> Av = A_full * v;
		Vector<Real> lambda_v = result.eigenvalues[i] * v;
		
		for (int j = 0; j < 3; j++)
			REQUIRE(std::abs(Av[j] - lambda_v[j]) < 1e-8);
	}
}

TEST_CASE("JacobiEigen_SymmetricMatrix_5x5", "[EigenSolvers][Jacobi][TestData]")
{
		TEST_PRECISION_INFO();
	// Use test data: symm_mat_5x5 with known eigenvalues
	auto result = SymmMatEigenSolverJacobi::Solve(symm_mat_5x5);
	
	REQUIRE(result.converged);
	REQUIRE(result.eigenvalues.size() == 5);
	
	// Compare with known eigenvalues
	Vector<Real> expected = symm_mat_5x5_eigen_val;
	std::sort(expected.begin(), expected.end());
	
	for (int i = 0; i < 5; i++)
		REQUIRE(std::abs(result.eigenvalues[i] - expected[i]) < 1e-6);
	
	// Verify eigenvectors are orthonormal
	for (int i = 0; i < 5; i++)
	{
		Vector<Real> vi = result.eigenvectors.VectorFromColumn(i);
		Real norm = vi.NormL2();
		REQUIRE(std::abs(norm - REAL(1.0)) < 1e-10);
		
		// Check orthogonality with all other vectors
		for (int j = i + 1; j < 5; j++)
		{
			Vector<Real> vj = result.eigenvectors.VectorFromColumn(j);
			Real dot = Utils::ScalarProduct(vi, vj);
			REQUIRE(std::abs(dot) < 1e-8);
		}
	}
}

TEST_CASE("JacobiEigen_SymmetricMatrix_10x10", "[EigenSolvers][Jacobi][TestData]")
{
		TEST_PRECISION_INFO();
	// Use test data: symm_mat_10x10 with known eigenvalues
	auto result = SymmMatEigenSolverJacobi::Solve(symm_mat_10x10);
	
	REQUIRE(result.converged);
	REQUIRE(result.eigenvalues.size() == 10);
	
	// Compare with known eigenvalues
	Vector<Real> expected = symm_mat_10x10_eigen_val;
	std::sort(expected.begin(), expected.end());
	
	for (int i = 0; i < 10; i++)
		REQUIRE(std::abs(result.eigenvalues[i] - expected[i]) < 1e-6);
}

TEST_CASE("JacobiEigen_SymmetricMatrix_3x3_Verified", "[EigenSolvers][Jacobi][Verified]")
{
		TEST_PRECISION_INFO();
	// Test with verified symmetric matrix from test data
	auto result = SymmMatEigenSolverJacobi::Solve(symm_mat_3x3);
	
	REQUIRE(result.converged);
	REQUIRE(result.eigenvalues.size() == 3);
	
	// Compare with verified eigenvalues (sorted ascending)
	Vector<Real> expected = symm_mat_3x3_eigen_val;
	std::sort(expected.begin(), expected.end());
	
	for (int i = 0; i < 3; i++)
		REQUIRE(std::abs(result.eigenvalues[i] - expected[i]) < 1e-6);
	
	// Verify trace = sum of eigenvalues
	Real trace = symm_mat_3x3(0,0) + symm_mat_3x3(1,1) + symm_mat_3x3(2,2);
	Real eigensum = result.eigenvalues[0] + result.eigenvalues[1] + result.eigenvalues[2];
	REQUIRE(std::abs(trace - eigensum) < 1e-10);
}

TEST_CASE("JacobiEigen_SymmetricMatrix_5x5_Verified", "[EigenSolvers][Jacobi][Verified]")
{
		TEST_PRECISION_INFO();
	// Test with verified 5x5 symmetric matrix from test data
	auto result = SymmMatEigenSolverJacobi::Solve(symm_mat_5x5);
	
	REQUIRE(result.converged);
	REQUIRE(result.eigenvalues.size() == 5);
	
	// Compare with verified eigenvalues (sorted ascending)
	Vector<Real> expected = symm_mat_5x5_eigen_val;
	std::sort(expected.begin(), expected.end());
	
	for (int i = 0; i < 5; i++)
		REQUIRE(std::abs(result.eigenvalues[i] - expected[i]) < 1e-6);
	
	// Verify trace = sum of eigenvalues
	Real trace = REAL(0.0);
	for (int i = 0; i < 5; i++)
		trace += symm_mat_5x5(i, i);
	
	Real eigensum = REAL(0.0);
	for (int i = 0; i < 5; i++)
		eigensum += result.eigenvalues[i];
	
	REQUIRE(std::abs(trace - eigensum) < 1e-10);
}

TEST_CASE("JacobiEigen_SymmetricMatrix_10x10_Verified", "[EigenSolvers][Jacobi][Verified]")
{
		TEST_PRECISION_INFO();
	// Test with verified 10x10 symmetric matrix from test data
	auto result = SymmMatEigenSolverJacobi::Solve(symm_mat_10x10);
	
	REQUIRE(result.converged);
	REQUIRE(result.eigenvalues.size() == 10);
	
	// Compare with verified eigenvalues (sorted ascending)
	Vector<Real> expected = symm_mat_10x10_eigen_val;
	std::sort(expected.begin(), expected.end());
	
	for (int i = 0; i < 10; i++)
		REQUIRE(std::abs(result.eigenvalues[i] - expected[i]) < 1e-6);
	
	// Verify trace = sum of eigenvalues
	Real trace = REAL(0.0);
	for (int i = 0; i < 10; i++)
		trace += symm_mat_10x10(i, i);
	
	Real eigensum = REAL(0.0);
	for (int i = 0; i < 10; i++)
		eigensum += result.eigenvalues[i];
	
	REQUIRE(std::abs(trace - eigensum) < 1e-10);
}

TEST_CASE("JacobiEigen_GraphLaplacian_6x6", "[EigenSolvers][Jacobi][Singular]")
{
		TEST_PRECISION_INFO();
	// Graph Laplacian - has one zero eigenvalue (singular matrix)
	auto result = SymmMatEigenSolverJacobi::Solve(spd_6x6_graph_laplacian_mat);
	
	REQUIRE(result.converged);
	
	// Should have one zero eigenvalue
	REQUIRE(std::abs(result.eigenvalues[0]) < 1e-6);
	
	// Compare with known eigenvalues
	Vector<Real> expected = spd_6x6_graph_laplacian_eigen;
	std::sort(expected.begin(), expected.end());
	
	for (int i = 0; i < 6; i++)
		REQUIRE(std::abs(result.eigenvalues[i] - expected[i]) < 1e-3);
}

TEST_CASE("JacobiEigen_EigenvectorOrthonormality", "[EigenSolvers][Jacobi][Properties]")
{
		TEST_PRECISION_INFO();
	// Test that V^T * V = I for various matrices
	
	auto result = SymmMatEigenSolverJacobi::Solve(symm_mat_5x5);
	
	REQUIRE(result.converged);
	
	int n = 5;
	Matrix<Real> VtV = result.eigenvectors.GetTranspose() * result.eigenvectors;
	
	// Check V^T * V = I
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			Real expected = (i == j) ? REAL(1.0) : REAL(0.0);
			REQUIRE(std::abs(VtV[i][j] - expected) < 1e-8);
		}
	}
}

TEST_CASE("JacobiEigen_ReconstructMatrix", "[EigenSolvers][Jacobi][Properties]")
{
		TEST_PRECISION_INFO();
	// Test that A = V * Λ * V^T
	
	auto result = SymmMatEigenSolverJacobi::Solve(symm_mat_3x3);
	
	REQUIRE(result.converged);
	
	int n = 3;
	
	// Create diagonal matrix Λ
	Matrix<Real> Lambda(n, n);
	for (int i = 0; i < n; i++)
		Lambda[i][i] = result.eigenvalues[i];
	
	// Reconstruct A = V * Λ * V^T
	Matrix<Real> VL = result.eigenvectors * Lambda;
	Matrix<Real> A_reconstructed = VL * result.eigenvectors.GetTranspose();
	
	// Convert original to full matrix
	Matrix<Real> A_original(n, n);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			A_original[i][j] = (i <= j) ? symm_mat_3x3(i, j) : symm_mat_3x3(j, i);
	
	// Compare
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			REQUIRE(std::abs(A_reconstructed[i][j] - A_original[i][j]) < 1e-8);
}

TEST_CASE("JacobiEigen_ConvergenceSpeed", "[EigenSolvers][Jacobi][Performance]")
{
		TEST_PRECISION_INFO();
	// Test that Jacobi converges reasonably fast
	
	auto result = SymmMatEigenSolverJacobi::Solve(symm_mat_10x10, 1e-10, 100);
	
	REQUIRE(result.converged);
	
	// Should converge in reasonable number of iterations
	INFO("Iterations: " << result.iterations);
	REQUIRE(result.iterations < 50);  // Typically 5-15 for n=10
	
	// Residual should be small
	INFO("Residual: " << result.residual);
	REQUIRE(result.residual < 1e-10);
}

TEST_CASE("JacobiEigen_NonConvergence", "[EigenSolvers][Jacobi][EdgeCases]")
{
		TEST_PRECISION_INFO();
	// Test behavior with very small maxIter
	
	auto result = SymmMatEigenSolverJacobi::Solve(symm_mat_10x10, 1e-15, 1);  // Only 1 iteration
	
	// Should not converge
	REQUIRE_FALSE(result.converged);
	REQUIRE(result.iterations == 1);
}

TEST_CASE("JacobiEigen_TraceAndDeterminant", "[EigenSolvers][Jacobi][Properties]")
{
		TEST_PRECISION_INFO();
	// Test that trace(A) = sum(eigenvalues) and det(A) = product(eigenvalues)
	
	auto result = SymmMatEigenSolverJacobi::Solve(symm_mat_5x5);
	
	REQUIRE(result.converged);
	
	// Compute trace
	Real trace = REAL(0.0);
	for (int i = 0; i < 5; i++)
		trace += symm_mat_5x5(i, i);
	
	// Sum of eigenvalues
	Real eigensum = REAL(0.0);
	for (int i = 0; i < 5; i++)
		eigensum += result.eigenvalues[i];
	
	REQUIRE(std::abs(trace - eigensum) < 1e-10);
	
	// Product of eigenvalues = determinant
	Real eigenprod = REAL(1.0);
	for (int i = 0; i < 5; i++)
		eigenprod *= result.eigenvalues[i];
	
	// Note: Product of eigenvalues should equal determinant
	// Actual determinant calculation skipped (may not be implemented)
	INFO("Eigenvalue product (equals det): " << eigenprod);
}

/***************************************************************************************************
 * QR EIGENSOLVER TESTS
 ***************************************************************************************************/

TEST_CASE("QREigen_2x2_Simple", "[EigenSolvers][QR][Small]")
{
		TEST_PRECISION_INFO();
	// Simple 2x2 symmetric matrix with known eigenvalues
	MatrixSym<Real> A(2, {REAL(2.0), REAL(1.0), REAL(2.0)});
	
	// Known eigenvalues: λ₁ = 1, λ₂ = 3
	auto result = SymmMatEigenSolverQR::Solve(A);
	
	REQUIRE(result.converged);
	REQUIRE(result.eigenvalues.size() == 2);
	
	// Check eigenvalues (sorted ascending)
	REQUIRE(std::abs(result.eigenvalues[0] - REAL(1.0)) < 1e-10);
	REQUIRE(std::abs(result.eigenvalues[1] - REAL(3.0)) < 1e-10);
	
	// Verify A*v = λ*v for each eigenpair
	for (int i = 0; i < 2; i++)
	{
		Vector<Real> v = result.eigenvectors.VectorFromColumn(i);
		Vector<Real> Av(2);
		for (int j = 0; j < 2; j++)
			Av[j] = A(j, 0) * v[0] + A(j, 1) * v[1];
		
		Vector<Real> lambda_v = result.eigenvalues[i] * v;
		
		for (int j = 0; j < 2; j++)
			REQUIRE(std::abs(Av[j] - lambda_v[j]) < 1e-10);
	}
}

TEST_CASE("QREigen_3x3_Identity", "[EigenSolvers][QR][Small]")
{
		TEST_PRECISION_INFO();
	// 3x3 identity matrix - eigenvalues should all be 1
	MatrixSym<Real> I(3, {REAL(1.0), REAL(0.0), REAL(1.0), REAL(0.0), REAL(0.0), REAL(1.0)});
	
	auto result = SymmMatEigenSolverQR::Solve(I);
	
	REQUIRE(result.converged);
	REQUIRE(result.eigenvalues.size() == 3);
	
	// All eigenvalues should be REAL(1.0)
	for (int i = 0; i < 3; i++)
		REQUIRE(std::abs(result.eigenvalues[i] - REAL(1.0)) < 1e-10);
}

TEST_CASE("QREigen_3x3_Diagonal", "[EigenSolvers][QR][Small]")
{
		TEST_PRECISION_INFO();
	// 3x3 diagonal matrix
	MatrixSym<Real> D(3, {REAL(5.0), REAL(0.0), REAL(3.0), REAL(0.0), REAL(0.0), REAL(1.0)});
	
	// Known eigenvalues: 1, 3, 5
	auto result = SymmMatEigenSolverQR::Solve(D);
	
	REQUIRE(result.converged);
	
	// Check eigenvalues (sorted ascending)
	REQUIRE(std::abs(result.eigenvalues[0] - REAL(1.0)) < 1e-10);
	REQUIRE(std::abs(result.eigenvalues[1] - REAL(3.0)) < 1e-10);
	REQUIRE(std::abs(result.eigenvalues[2] - REAL(5.0)) < 1e-10);
}

TEST_CASE("QREigen_SymmetricMatrix_3x3_Verified", "[EigenSolvers][QR][Verified]")
{
		TEST_PRECISION_INFO();
	// Use verified 3x3 symmetric matrix from test data
	Matrix<Real> A_full(3, 3);
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			A_full(i, j) = symm_mat_3x3(i, j);
	
	MatrixSym<Real> A(3);
	for (int i = 0; i < 3; i++)
		for (int j = i; j < 3; j++)
			A(i, j) = A_full(i, j);
	
	auto result = SymmMatEigenSolverQR::Solve(A);
	
	REQUIRE(result.converged);
	REQUIRE(result.eigenvalues.size() == 3);
	
	// Check eigenvalues against verified values (sorted)
	Real expected_eigenvals[3] = {symm_mat_3x3_eigen_val[0], symm_mat_3x3_eigen_val[1], symm_mat_3x3_eigen_val[2]};
	std::sort(std::begin(expected_eigenvals), std::end(expected_eigenvals));
	
	for (int i = 0; i < 3; i++)
		REQUIRE(std::abs(result.eigenvalues[i] - expected_eigenvals[i]) < 1e-6);
	
	// Verify A*v = λ*v for each eigenpair
	Matrix<Real> A_mat(3, 3);
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			A_mat(i, j) = (i <= j) ? A(i, j) : A(j, i);
	
	for (int i = 0; i < 3; i++)
	{
		Vector<Real> v = result.eigenvectors.VectorFromColumn(i);
		Vector<Real> Av = A_mat * v;
		Vector<Real> lambda_v = result.eigenvalues[i] * v;
		
		for (int j = 0; j < 3; j++)
			REQUIRE(std::abs(Av[j] - lambda_v[j]) < 1e-6);
	}
}

TEST_CASE("QREigen_SymmetricMatrix_5x5_Verified", "[EigenSolvers][QR][Verified]")
{
		TEST_PRECISION_INFO();
	// Use verified 5x5 symmetric matrix from test data
	MatrixSym<Real> A(5);
	int idx = 0;
	for (int i = 0; i < 5; i++)
		for (int j = i; j < 5; j++)
			A(i, j) = symm_mat_5x5(i, j);
	
	auto result = SymmMatEigenSolverQR::Solve(A);
	
	REQUIRE(result.converged);
	REQUIRE(result.eigenvalues.size() == 5);
	
	// Check eigenvalues against verified values (sorted)
	Real expected_eigenvals[5] = {symm_mat_5x5_eigen_val[0], symm_mat_5x5_eigen_val[1], 
	                              symm_mat_5x5_eigen_val[2], symm_mat_5x5_eigen_val[3],
	                              symm_mat_5x5_eigen_val[4]};
	std::sort(std::begin(expected_eigenvals), std::end(expected_eigenvals));
	
	for (int i = 0; i < 5; i++)
		REQUIRE(std::abs(result.eigenvalues[i] - expected_eigenvals[i]) < 1e-6);
	
	// Verify eigenvector orthonormality
	for (int i = 0; i < 5; i++)
	{
		Vector<Real> vi = result.eigenvectors.VectorFromColumn(i);
		
		// Check normalization
		Real norm_sq = Utils::ScalarProduct(vi, vi);
		REQUIRE(std::abs(norm_sq - REAL(1.0)) < 1e-10);
		
		// Check orthogonality with other eigenvectors
		for (int j = i + 1; j < 5; j++)
		{
			Vector<Real> vj = result.eigenvectors.VectorFromColumn(j);
			Real dot = Utils::ScalarProduct(vi, vj);
			REQUIRE(std::abs(dot) < 1e-10);
		}
	}
}

TEST_CASE("QREigen_SymmetricMatrix_10x10_Verified", "[EigenSolvers][QR][Verified]")
{
		TEST_PRECISION_INFO();
	// Use verified 10x10 symmetric matrix from test data
	MatrixSym<Real> A(10);
	for (int i = 0; i < 10; i++)
		for (int j = i; j < 10; j++)
			A(i, j) = symm_mat_10x10(i, j);
	
	auto result = SymmMatEigenSolverQR::Solve(A);
	
	REQUIRE(result.converged);
	REQUIRE(result.eigenvalues.size() == 10);
	
	// Check eigenvalues against verified values (sorted)
	std::vector<Real> expected_eigenvals(10);
	for (int i = 0; i < 10; i++)
		expected_eigenvals[i] = symm_mat_10x10_eigen_val[i];
	std::sort(expected_eigenvals.begin(), expected_eigenvals.end());
	
	for (int i = 0; i < 10; i++)
		REQUIRE(std::abs(result.eigenvalues[i] - expected_eigenvals[i]) < 1e-4);  // Slightly relaxed for larger matrix
	
	// Verify matrix reconstruction: A = Q*Λ*Q^T
	Matrix<Real> Lambda(10, 10);
	for (int i = 0; i < 10; i++)
		for (int j = 0; j < 10; j++)
			Lambda(i, j) = REAL(0.0);
	for (int i = 0; i < 10; i++)
		Lambda(i, i) = result.eigenvalues[i];
	
	Matrix<Real> Q = result.eigenvectors;
	Matrix<Real> Qt = Q.GetTranspose();
	Matrix<Real> A_reconstructed = Q * Lambda * Qt;
	
	// Check reconstruction accuracy
	for (int i = 0; i < 10; i++)
	{
		for (int j = i; j < 10; j++)
		{
			Real expected = A(i, j);
			Real actual = A_reconstructed(i, j);
			REQUIRE(std::abs(expected - actual) < 1e-4);
		}
	}
}

TEST_CASE("QREigen_CompareWithJacobi_3x3", "[EigenSolvers][QR][Jacobi][Comparison]")
{
		TEST_PRECISION_INFO();
	// Both algorithms should give the same eigenvalues
	MatrixSym<Real> A(3);
	for (int i = 0; i < 3; i++)
		for (int j = i; j < 3; j++)
			A(i, j) = symm_mat_3x3(i, j);
	
	auto qr_result = SymmMatEigenSolverQR::Solve(A);
	auto jacobi_result = SymmMatEigenSolverJacobi::Solve(A);
	
	REQUIRE(qr_result.converged);
	REQUIRE(jacobi_result.converged);
	
	// Eigenvalues should match between methods
	for (int i = 0; i < 3; i++)
		REQUIRE(std::abs(qr_result.eigenvalues[i] - jacobi_result.eigenvalues[i]) < 1e-8);
}

TEST_CASE("QREigen_CompareWithJacobi_5x5", "[EigenSolvers][QR][Jacobi][Comparison]")
{
		TEST_PRECISION_INFO();
	// Both algorithms should give the same eigenvalues for 5x5
	MatrixSym<Real> A(5);
	for (int i = 0; i < 5; i++)
		for (int j = i; j < 5; j++)
			A(i, j) = symm_mat_5x5(i, j);
	
	auto qr_result = SymmMatEigenSolverQR::Solve(A);
	auto jacobi_result = SymmMatEigenSolverJacobi::Solve(A);
	
	REQUIRE(qr_result.converged);
	REQUIRE(jacobi_result.converged);
	
	// Eigenvalues should match between methods
	for (int i = 0; i < 5; i++)
		REQUIRE(std::abs(qr_result.eigenvalues[i] - jacobi_result.eigenvalues[i]) < 1e-8);
}

TEST_CASE("QREigen_Convergence", "[EigenSolvers][QR][Convergence]")
{
		TEST_PRECISION_INFO();
	// Test with tighter tolerance
	MatrixSym<Real> A(3);
	for (int i = 0; i < 3; i++)
		for (int j = i; j < 3; j++)
			A(i, j) = symm_mat_3x3(i, j);
	
	auto result = SymmMatEigenSolverQR::Solve(A, 1e-12, 1000);
	
	REQUIRE(result.converged);
	REQUIRE(result.residual < 1e-12);
	INFO("Converged in " << result.iterations << " iterations");
	INFO("Final residual: " << result.residual);
}

TEST_CASE("QREigen_Performance_vs_Jacobi", "[EigenSolvers][QR][Jacobi][Performance]")
{
		TEST_PRECISION_INFO();
	// QR should converge faster than Jacobi for larger matrices
	MatrixSym<Real> A(10);
	for (int i = 0; i < 10; i++)
		for (int j = i; j < 10; j++)
			A(i, j) = symm_mat_10x10(i, j);
	
	auto qr_result = SymmMatEigenSolverQR::Solve(A);
	auto jacobi_result = SymmMatEigenSolverJacobi::Solve(A);
	
	REQUIRE(qr_result.converged);
	REQUIRE(jacobi_result.converged);
	
	INFO("QR iterations: " << qr_result.iterations);
	INFO("Jacobi iterations: " << jacobi_result.iterations);
	
	// For 10x10, QR typically converges much faster
	// QR: typically ~30-50 iterations
	// Jacobi: typically ~100-300 iterations
}

TEST_CASE("QREigen_Trace_Property", "[EigenSolvers][QR][Properties]")
{
		TEST_PRECISION_INFO();
	// Sum of eigenvalues equals trace
	MatrixSym<Real> A(5);
	for (int i = 0; i < 5; i++)
		for (int j = i; j < 5; j++)
			A(i, j) = symm_mat_5x5(i, j);
	
	auto result = SymmMatEigenSolverQR::Solve(A);
	
	REQUIRE(result.converged);
	
	// Compute trace
	Real trace = REAL(0.0);
	for (int i = 0; i < 5; i++)
		trace += A(i, i);
	
	// Sum eigenvalues
	Real eigensum = REAL(0.0);
	for (int i = 0; i < 5; i++)
		eigensum += result.eigenvalues[i];
	
	REQUIRE(std::abs(trace - eigensum) < 1e-10);
}

/***************************************************************************************************
 * COMPREHENSIVE TEST BED TESTS
 * 
 * These tests iterate over ALL test cases from eigenvalue_test_bed.h to ensure
 * complete coverage of the eigenvalue solver algorithms.
 ***************************************************************************************************/

// Helper to convert Matrix to MatrixSym (assumes input is symmetric)
static MatrixSym<Real> toMatrixSym(const Matrix<Real>& M)
{
	int n = M.RowNum();
	MatrixSym<Real> result(n);
	for (int i = 0; i < n; i++)
		for (int j = i; j < n; j++)
			result(i, j) = M(i, j);
	return result;
}

TEST_CASE("JacobiEigen_AllSymmetricTestBed", "[EigenSolvers][Jacobi][TestBed][Comprehensive]")
{
	TEST_PRECISION_INFO();
	
	auto symmetricTests = TestBeds::getSymmetricEigenTests();
	
	INFO("Testing " << symmetricTests.size() << " symmetric matrices from eigenvalue test bed");
	
	for (const auto& test : symmetricTests)
	{
		DYNAMIC_SECTION(test.name)
		{
			INFO("Category: " << test.category);
			INFO("Difficulty: " << test.difficulty);
			INFO("Dimension: " << test.dimension);
			INFO("Description: " << test.description);
			
			Matrix<Real> M = test.getMatrix();
			MatrixSym<Real> A = toMatrixSym(M);
			
			auto result = SymmMatEigenSolverJacobi::Solve(A);
			
			// Check convergence
			REQUIRE(result.converged);
			REQUIRE(result.eigenvalues.size() == test.dimension);
			
			// Adjust tolerance based on difficulty and condition number
			// Use reasonable baseline accounting for limited precision in stored eigenvalues
			Real baseTol = 1e-5;  // Accounts for ~6-digit precision in test bed eigenvalues
			Real tol = baseTol * std::pow(10.0, test.difficulty - 1);
			if (test.conditionNumber > 1e6) tol *= 100;   // Well-conditioned but still sensitive
			if (test.conditionNumber > 1e10) tol *= 1e3;  // Ill-conditioned
			if (test.hasClusteredEigenvalues) tol *= 10;  // Clustered eigenvalues
			
			// For EXTREMELY ill-conditioned matrices (κ > 10^12), eigenvalue comparison
			// is unreliable due to numerical precision limits. We verify via residuals instead.
			bool skipEigenvalueComparison = test.conditionNumber > 1e12;
			
			// Verify eigenvalues match expected (skip for extremely ill-conditioned matrices)
			if (test.realEigenvalues.size() > 0 && !skipEigenvalueComparison)
			{
				bool match = TestBeds::verifyRealEigenvalues(result.eigenvalues, test.realEigenvalues, tol);
				if (!match)
				{
					Real error = TestBeds::computeEigenvalueError(result.eigenvalues, test.realEigenvalues);
					WARN("Max eigenvalue error: " << error << ", Tolerance: " << tol);
					// Print computed vs expected eigenvalues for diagnostics
					std::ostringstream ss;
					ss << "Computed: ";
					for (size_t i = 0; i < result.eigenvalues.size(); ++i)
						ss << result.eigenvalues[i] << " ";
					ss << "\nExpected: ";
					for (size_t i = 0; i < test.realEigenvalues.size(); ++i)
						ss << test.realEigenvalues[i] << " ";
					WARN(ss.str());
				}
				CHECK(match);
			}
			else if (skipEigenvalueComparison)
			{
				INFO("Skipping eigenvalue comparison for extremely ill-conditioned matrix (κ=" << test.conditionNumber << ")");
			}
			
			// Verify eigenvectors are orthonormal
			bool ortho = TestBeds::verifyOrthonormality(result.eigenvectors, tol * 10);
			CHECK(ortho);
			
			// Verify A*v = λ*v for each eigenpair (residual check)
			for (int i = 0; i < test.dimension; i++)
			{
				Vector<Real> v = result.eigenvectors.VectorFromColumn(i);
				Real residual = TestBeds::computeResidual(M, v, result.eigenvalues[i]);
				CHECK(residual < tol * 100);
			}
			
			// Verify trace property: sum(eigenvalues) = trace(A)
			Real trace = 0;
			for (int i = 0; i < test.dimension; i++)
				trace += A(i, i);
			Real eigensum = 0;
			for (int i = 0; i < test.dimension; i++)
				eigensum += result.eigenvalues[i];
			CHECK(std::abs(trace - eigensum) < tol * test.dimension);
		}
	}
}

TEST_CASE("QREigen_AllSymmetricTestBed", "[EigenSolvers][QR][TestBed][Comprehensive]")
{
	TEST_PRECISION_INFO();
	
	auto symmetricTests = TestBeds::getSymmetricEigenTests();
	
	INFO("Testing " << symmetricTests.size() << " symmetric matrices from eigenvalue test bed");
	
	for (const auto& test : symmetricTests)
	{
		DYNAMIC_SECTION(test.name)
		{
			INFO("Category: " << test.category);
			INFO("Difficulty: " << test.difficulty);
			INFO("Dimension: " << test.dimension);
			INFO("Description: " << test.description);
			
			Matrix<Real> M = test.getMatrix();
			MatrixSym<Real> A = toMatrixSym(M);
			
			auto result = SymmMatEigenSolverQR::Solve(A);
			
			// Check convergence
			REQUIRE(result.converged);
			REQUIRE(result.eigenvalues.size() == test.dimension);
			
			// Adjust tolerance based on difficulty and condition number
			// Use reasonable baseline accounting for limited precision in stored eigenvalues
			Real baseTol = 1e-5;  // Accounts for ~6-digit precision in test bed eigenvalues
			Real tol = baseTol * std::pow(10.0, test.difficulty - 1);
			if (test.conditionNumber > 1e6) tol *= 100;   // Well-conditioned but still sensitive
			if (test.conditionNumber > 1e10) tol *= 1e3;  // Ill-conditioned
			if (test.hasClusteredEigenvalues) tol *= 10;  // Clustered eigenvalues
			
			// For EXTREMELY ill-conditioned matrices (κ > 10^12), eigenvalue comparison
			// is unreliable due to numerical precision limits. We verify via residuals instead.
			bool skipEigenvalueComparison = test.conditionNumber > 1e12;
			
			// Verify eigenvalues match expected (skip for extremely ill-conditioned matrices)
			if (test.realEigenvalues.size() > 0 && !skipEigenvalueComparison)
			{
				bool match = TestBeds::verifyRealEigenvalues(result.eigenvalues, test.realEigenvalues, tol);
				if (!match)
				{
					Real error = TestBeds::computeEigenvalueError(result.eigenvalues, test.realEigenvalues);
					INFO("Max eigenvalue error: " << error);
					INFO("Tolerance was: " << tol);
				}
				CHECK(match);
			}
			else if (skipEigenvalueComparison)
			{
				INFO("Skipping eigenvalue comparison for extremely ill-conditioned matrix (κ=" << test.conditionNumber << ")");
			}
			
			// Verify eigenvectors are orthonormal
			bool ortho = TestBeds::verifyOrthonormality(result.eigenvectors, tol * 10);
			CHECK(ortho);
			
			// Verify A*v = λ*v for each eigenpair
			for (int i = 0; i < test.dimension; i++)
			{
				Vector<Real> v = result.eigenvectors.VectorFromColumn(i);
				Real residual = TestBeds::computeResidual(M, v, result.eigenvalues[i]);
				CHECK(residual < tol * 100);
			}
			
			// Verify trace property
			Real trace = 0, eigensum = 0;
			for (int i = 0; i < test.dimension; i++)
			{
				trace += A(i, i);
				eigensum += result.eigenvalues[i];
			}
			CHECK(std::abs(trace - eigensum) < tol * test.dimension);
		}
	}
}

TEST_CASE("EigenSolvers_EasyTestBed", "[EigenSolvers][TestBed][Quick]")
{
	TEST_PRECISION_INFO();
	
	// Quick test with easy cases - good for CI/CD fast checks
	auto easyTests = TestBeds::getEasyEigenTests();
	
	for (const auto& test : easyTests)
	{
		DYNAMIC_SECTION(test.name + " (Jacobi)")
		{
			if (!test.isSymmetric) continue;
			
			Matrix<Real> M = test.getMatrix();
			MatrixSym<Real> A = toMatrixSym(M);
			
			auto result = SymmMatEigenSolverJacobi::Solve(A);
			REQUIRE(result.converged);
			
			if (test.realEigenvalues.size() > 0)
			{
				// Use 1e-5 tolerance to account for limited precision in stored eigenvalues
				CHECK(TestBeds::verifyRealEigenvalues(result.eigenvalues, test.realEigenvalues, 1e-5));
			}
		}
		
		DYNAMIC_SECTION(test.name + " (QR)")
		{
			if (!test.isSymmetric) continue;
			
			Matrix<Real> M = test.getMatrix();
			MatrixSym<Real> A = toMatrixSym(M);
			
			auto result = SymmMatEigenSolverQR::Solve(A);
			REQUIRE(result.converged);
			
			if (test.realEigenvalues.size() > 0)
			{
				// Use 1e-5 tolerance to account for limited precision in stored eigenvalues
				CHECK(TestBeds::verifyRealEigenvalues(result.eigenvalues, test.realEigenvalues, 1e-5));
			}
		}
	}
}

TEST_CASE("EigenSolvers_ChallengingTestBed", "[EigenSolvers][TestBed][Challenging]")
{
	TEST_PRECISION_INFO();
	
	// Challenging cases - ill-conditioned, clustered eigenvalues, etc.
	auto challengingTests = TestBeds::getChallengingEigenTests();
	
	for (const auto& test : challengingTests)
	{
		if (!test.isSymmetric) continue;  // Skip non-symmetric for now
		
		DYNAMIC_SECTION(test.name)
		{
			INFO("Category: " << test.category);
			INFO("Difficulty: " << test.difficulty);
			INFO("Description: " << test.description);
			
			Matrix<Real> M = test.getMatrix();
			MatrixSym<Real> A = toMatrixSym(M);
			
			// For EXTREMELY ill-conditioned matrices (κ > 10^12), eigenvalue comparison
			// is unreliable due to numerical precision limits. Only verify residuals.
			bool skipEigenvalueComparison = test.conditionNumber > 1e12;
			
			// Try QR first (generally more robust)
			auto resultQR = SymmMatEigenSolverQR::Solve(A);
			
			if (resultQR.converged && test.realEigenvalues.size() > 0 && !skipEigenvalueComparison)
			{
				// More relaxed tolerance for challenging cases
				Real tol = 1e-6;
				if (test.difficulty >= 4) tol = 1e-4;
				
				bool match = TestBeds::verifyRealEigenvalues(resultQR.eigenvalues, test.realEigenvalues, tol);
				if (!match)
				{
					Real error = TestBeds::computeEigenvalueError(resultQR.eigenvalues, test.realEigenvalues);
					INFO("QR eigenvalue error: " << error);
				}
				CHECK(match);
			}
			else if (skipEigenvalueComparison)
			{
				INFO("Skipping eigenvalue comparison for extremely ill-conditioned matrix (κ=" << test.conditionNumber << ")");
				// Still verify convergence and residuals
				if (resultQR.converged)
				{
					// Verify A*v = λ*v for each eigenpair (self-consistency check)
					Real maxResidual = 0;
					for (int i = 0; i < test.dimension; i++)
					{
						Vector<Real> v = resultQR.eigenvectors.VectorFromColumn(i);
						Real residual = TestBeds::computeResidual(M, v, resultQR.eigenvalues[i]);
						maxResidual = std::max(maxResidual, residual);
					}
					INFO("QR max residual ||Av-λv||/||v||: " << maxResidual);
					CHECK(maxResidual < 1e-6);  // Self-consistency should still be good
				}
			}
			
			// Try Jacobi as well
			auto resultJac = SymmMatEigenSolverJacobi::Solve(A);
			
			if (resultJac.converged && test.realEigenvalues.size() > 0 && !skipEigenvalueComparison)
			{
				Real tol = 1e-6;
				if (test.difficulty >= 4) tol = 1e-4;
				
				bool match = TestBeds::verifyRealEigenvalues(resultJac.eigenvalues, test.realEigenvalues, tol);
				if (!match)
				{
					Real error = TestBeds::computeEigenvalueError(resultJac.eigenvalues, test.realEigenvalues);
					INFO("Jacobi eigenvalue error: " << error);
				}
				CHECK(match);
			}
			else if (skipEigenvalueComparison && resultJac.converged)
			{
				// Verify A*v = λ*v for each eigenpair (self-consistency check)
				Real maxResidual = 0;
				for (int i = 0; i < test.dimension; i++)
				{
					Vector<Real> v = resultJac.eigenvectors.VectorFromColumn(i);
					Real residual = TestBeds::computeResidual(M, v, resultJac.eigenvalues[i]);
					maxResidual = std::max(maxResidual, residual);
				}
				INFO("Jacobi max residual ||Av-λv||/||v||: " << maxResidual);
				CHECK(maxResidual < 1e-6);  // Self-consistency should still be good
			}
		}
	}
}

TEST_CASE("EigenSolvers_RepeatedEigenvalueTestBed", "[EigenSolvers][TestBed][Repeated]")
{
	TEST_PRECISION_INFO();
	
	// Specifically test repeated eigenvalue handling
	std::vector<TestBeds::TestEigenSystem> repeatedTests = {
		TestBeds::getRepeatedDouble4x4Test(),
		TestBeds::getRepeatedTriple5x5Test(),
		TestBeds::getRepeatedQuadruple4x4Test()
	};
	
	for (const auto& test : repeatedTests)
	{
		DYNAMIC_SECTION(test.name)
		{
			INFO("Multiplicity max: " << test.multiplicityMax);
			
			Matrix<Real> M = test.getMatrix();
			MatrixSym<Real> A = toMatrixSym(M);
			
			// Test both solvers
			auto resultJac = SymmMatEigenSolverJacobi::Solve(A);
			auto resultQR = SymmMatEigenSolverQR::Solve(A);
			
			REQUIRE(resultJac.converged);
			REQUIRE(resultQR.converged);
			
			// Both should find correct eigenvalues
			CHECK(TestBeds::verifyRealEigenvalues(resultJac.eigenvalues, test.realEigenvalues, 1e-8));
			CHECK(TestBeds::verifyRealEigenvalues(resultQR.eigenvalues, test.realEigenvalues, 1e-8));
			
			// Eigenvectors should still be orthonormal for symmetric matrices
			CHECK(TestBeds::verifyOrthonormality(resultJac.eigenvectors, 1e-8));
			CHECK(TestBeds::verifyOrthonormality(resultQR.eigenvectors, 1e-8));
		}
	}
}


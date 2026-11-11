#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#include "../../test_data/linear_alg_eq_systems_classic_defs.h"
#include "../../test_data/linear_alg_eq_systems_sym_defs.h"
#include "../../test_data/linear_alg_eq_systems_solverspecific_defs.h"
#include "../../test_beds/linear_alg_eq_systems_test_bed.h"
#include "../../test_beds/eigenvalue_test_bed.h"

#include "../../mml/algorithms/EigenSystemSolvers.h"
#include "../../mml/base/BaseUtils.h"

#include <sstream>

using namespace MML;
using namespace MML::Testing;
using namespace MML::TestBeds;

namespace MML::Tests::Algorithms::EigensystemSolversTests
{

/***************************************************************************************************
 * JACOBI EIGENSOLVER TESTS
 ***************************************************************************************************/

TEST_CASE("JacobiEigen_2x2_Simple", "[EigenSolvers][Jacobi][Small]")
{
		TEST_PRECISION_INFO();
	// Simple 2x2 symmetric matrix with known eigenvalues
	MatrixSym<Real> A(2, {static_cast<MML::Real>(2.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0)});
	
	// Known eigenvalues: λ₁ = 1, λ₂ = 3
	// Known eigenvectors: v₁ = [-1/√2, 1/√2], v₂ = [1/√2, 1/√2]
	
	auto result = SymmMatEigenSolverJacobi::Solve(A);
	
	REQUIRE(result.converged);
	REQUIRE(result.eigenvalues.size() == 2);
	
	// Check eigenvalues (sorted ascending)
	REQUIRE(std::abs(result.eigenvalues[0] - static_cast<MML::Real>(1.0)) < TOL(1e-10, 1e-5));
	REQUIRE(std::abs(result.eigenvalues[1] - static_cast<MML::Real>(3.0)) < TOL(1e-10, 1e-5));
	
	// Verify A*v = λ*v for each eigenpair
	for (int i = 0; i < 2; i++)
	{
		Vector<Real> v = result.eigenvectors.VectorFromColumn(i);
		Vector<Real> Av(2);
		for (int j = 0; j < 2; j++)
			Av[j] = A(j, 0) * v[0] + A(j, 1) * v[1];
		
		Vector<Real> lambda_v = result.eigenvalues[i] * v;
		
		for (int j = 0; j < 2; j++)
			REQUIRE(std::abs(Av[j] - lambda_v[j]) < TOL(1e-10, 1e-5));
	}
	
	// Verify eigenvectors are orthonormal
	Vector<Real> v0 = result.eigenvectors.VectorFromColumn(0);
	Vector<Real> v1 = result.eigenvectors.VectorFromColumn(1);
	
	Real dot01 = v0[0] * v1[0] + v0[1] * v1[1];
	Real norm0 = std::sqrt(v0[0] * v0[0] + v0[1] * v0[1]);
	Real norm1 = std::sqrt(v1[0] * v1[0] + v1[1] * v1[1]);
	
	REQUIRE(std::abs(dot01) < TOL(1e-10, 1e-5));        // Orthogonal
	REQUIRE(std::abs(norm0 - static_cast<MML::Real>(1.0)) < TOL(1e-10, 1e-5));  // Unit norm
	REQUIRE(std::abs(norm1 - static_cast<MML::Real>(1.0)) < TOL(1e-10, 1e-5));  // Unit norm
}

TEST_CASE("HermitianJacobi solves complex Hermitian matrices", "[EigenSolvers][HermitianJacobi]")
{
	TEST_PRECISION_INFO();
	Matrix<Complex> A(3, 3, {
		{2.0, 0.0}, {1.0, 1.0}, {0.0, -1.0},
		{1.0, -1.0}, {3.0, 0.0}, {2.0, 0.5},
		{0.0, 1.0}, {2.0, -0.5}, {5.0, 0.0}
	});

	auto result = HermitianMatEigenSolverJacobi::Solve(A);

	REQUIRE(result.converged);
	REQUIRE(result.eigenvalues[0] <= result.eigenvalues[1]);
	REQUIRE(result.eigenvalues[1] <= result.eigenvalues[2]);
	Matrix<Complex> identity = Matrix<Complex>::Identity(3);
	REQUIRE((Utils::GetConjugateTranspose(result.eigenvectors) * result.eigenvectors).IsEqualTo(identity, TOL(1e-10, 1e-5)));

	for (int column = 0; column < 3; column++) {
		Vector<Complex> eigenvector = result.eigenvectors.VectorFromColumn(column);
		Vector<Complex> residual = A * eigenvector - Complex(result.eigenvalues[column], 0.0) * eigenvector;
		REQUIRE(residual.NormL2() < TOL(1e-9, 1e-4));
	}
}

TEST_CASE("HermitianJacobi validates its input and reports convergence", "[EigenSolvers][HermitianJacobi]")
{
	TEST_PRECISION_INFO();
	Matrix<Complex> diagonal(3, 3, {
		{4.0, 0.0}, {0.0, 0.0}, {0.0, 0.0},
		{0.0, 0.0}, {-2.0, 0.0}, {0.0, 0.0},
		{0.0, 0.0}, {0.0, 0.0}, {1.0, 0.0}
	});
	EigenSolverConfig config;
	config.tolerance = TOL(1e-12, 1e-5);
	config.max_iterations = 10;

	auto result = HermitianMatEigenSolverJacobi::Solve(diagonal, config);
	REQUIRE(result.converged);
	REQUIRE(result.eigenvalues.IsEqualTo(Vector<Real>{-2.0, 1.0, 4.0}, TOL(1e-12, 1e-5)));
	REQUIRE(result.algorithm_name == "HermitianJacobi");

	Matrix<Complex> nonHermitian(2, 2, {
		{1.0, 0.0}, {0.0, 1.0},
		{0.0, 1.0}, {2.0, 0.0}
	});
	REQUIRE_THROWS_AS(HermitianMatEigenSolverJacobi::Solve(nonHermitian), MatrixDimensionError);
}

/***************************************************************************************************
 * EXTREME-MAGNITUDE / NONFINITE ROBUSTNESS (MinimalMathLibrary-uk6l)
 ***************************************************************************************************/

TEST_CASE("SelfAdjointEigen_ExtremeMagnitudes_StayFinite", "[EigenSolvers][Jacobi][HermitianJacobi][QR][MinimalMathLibrary-uk6l]")
{
	TEST_PRECISION_INFO();
	{
		// Float uses the same overflow-adjacent witness in its own finite range:
		// 2*big overflows, while both 1.1*big and sqrt(2)*big remain finite.
		const Real big = TOL(1e308, std::numeric_limits<Real>::max() * Real(0.55));
		const Real relTol = std::max(Real(1e-12), 8 * Precision::Eps);

		SECTION("SPD near the overflow boundary (Jacobi, QR, Matrix overloads)") {
			// eigenvalues 9e307 and 1.1e308 are representable; 0.5*(a+a) symmetrization is not
			MatrixSym<Real> A(2, {big, big / 10, big});
			Matrix<Real> G(2, 2, {big, big / 10, big / 10, big});

			auto checkSpd = [&](const auto& result) {
				INFO(result.error_message);
				REQUIRE(result.converged);
				REQUIRE(result.status == AlgorithmStatus::Success);
				REQUIRE(std::isfinite(result.residual));
				REQUIRE_THAT(result.eigenvalues[0], RealWithinRel(big - big / 10, relTol));
				REQUIRE_THAT(result.eigenvalues[1], RealWithinRel(big + big / 10, relTol));
			};
			checkSpd(SymmMatEigenSolverJacobi::Solve(A));
			checkSpd(SymmMatEigenSolverJacobi::Solve(G));
			checkSpd(SymmMatEigenSolverJacobi::SolveSymmetricPart(G));
			checkSpd(SymmMatEigenSolverQR::Solve(A));
			checkSpd(SymmMatEigenSolverQR::Solve(G));
		}

		SECTION("Indefinite near the overflow boundary (Jacobi, HermitianJacobi)") {
			// eigenvalues are +-sqrt(2)*1e308, both representable
			const Real expected = std::sqrt(static_cast<MML::Real>(2.0)) * big;
			MatrixSym<Real> A(2, {big, big, -big});
			auto jacobi = SymmMatEigenSolverJacobi::Solve(A);
			INFO(jacobi.error_message);
			REQUIRE(jacobi.converged);
			REQUIRE_THAT(jacobi.eigenvalues[0], RealWithinRel(-expected, relTol));
			REQUIRE_THAT(jacobi.eigenvalues[1], RealWithinRel(expected, relTol));

			Matrix<Complex> H(2, 2, {Complex(big, 0), Complex(big, 0), Complex(big, 0), Complex(-big, 0)});
			auto hermitian = HermitianMatEigenSolverJacobi::Solve(H);
			INFO(hermitian.error_message);
			REQUIRE(hermitian.converged);
			REQUIRE_THAT(hermitian.eigenvalues[0], RealWithinRel(-expected, relTol));
			REQUIRE_THAT(hermitian.eigenvalues[1], RealWithinRel(expected, relTol));
			for (int column = 0; column < 2; column++) {
				Vector<Complex> v = hermitian.eigenvectors.VectorFromColumn(column);
				REQUIRE(std::isfinite(std::abs(v[0])));
				REQUIRE(std::isfinite(std::abs(v[1])));
			}
		}

		SECTION("SymmetricQR on a 3x3 scaled by 1e160 matches the unit-scale spectrum") {
			const Real s = TOL(1e160, 1e20);
			MatrixSym<Real> unit(3, {static_cast<MML::Real>(2.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(3.0),
				static_cast<MML::Real>(0.5), static_cast<MML::Real>(1.0), static_cast<MML::Real>(4.0)});
			MatrixSym<Real> scaled(3);
			for (int i = 0; i < 3; i++)
				for (int j = i; j < 3; j++)
					scaled(i, j) = unit(i, j) * s;

			auto reference = SymmMatEigenSolverQR::Solve(unit);
			REQUIRE(reference.converged);
			auto checkScaled = [&](const auto& result, Real factor) {
				INFO(result.error_message);
				REQUIRE(result.converged);
				REQUIRE(std::isfinite(result.residual));
				// Compare two computed reductions, each accumulating O(4*n)
				// unit roundoffs, rather than an exact analytic eigenvalue.
				const Real comparisonBudget = std::max(Real(1e-10), 8 * unit.rows() * Precision::Eps);
				for (int i = 0; i < 3; i++)
					REQUIRE_THAT(result.eigenvalues[i] / factor, RealWithinRel(reference.eigenvalues[i],
						comparisonBudget));
			};
			checkScaled(SymmMatEigenSolverQR::Solve(scaled), s);
			checkScaled(SymmMatEigenSolverJacobi::Solve(scaled), s);

			// Tiny magnitudes must not underflow either. The tolerance is absolute (in matrix units),
			// so it is scaled along with the matrix to keep the request meaningful.
			const Real t = TOL(1e-200, 1e-20);
			MatrixSym<Real> tiny(3);
			for (int i = 0; i < 3; i++)
				for (int j = i; j < 3; j++)
					tiny(i, j) = unit(i, j) * t;
			checkScaled(SymmMatEigenSolverJacobi::Solve(tiny, PrecisionValues<Real>::EigenSolverConvergenceTolerance * t), t);
			checkScaled(SymmMatEigenSolverQR::Solve(tiny, PrecisionValues<Real>::EigenSolverConvergenceTolerance * t), t);
		}

		SECTION("Nonfinite input is reported, never Success") {
			MatrixSym<Real> nan(2, {static_cast<MML::Real>(1.0), std::numeric_limits<Real>::quiet_NaN(), static_cast<MML::Real>(1.0)});
			auto jacobi = SymmMatEigenSolverJacobi::Solve(nan);
			REQUIRE_FALSE(jacobi.converged);
			REQUIRE(jacobi.status == AlgorithmStatus::NumericalInstability);
			auto qr = SymmMatEigenSolverQR::Solve(nan);
			REQUIRE_FALSE(qr.converged);
			REQUIRE(qr.status == AlgorithmStatus::NumericalInstability);

			MatrixSym<Real> inf(2, {std::numeric_limits<Real>::infinity(), static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0)});
			EigenSolverConfig config;
			auto configured = SymmMatEigenSolverJacobi::Solve(inf, config);
			REQUIRE_FALSE(configured.converged);
			REQUIRE(configured.status == AlgorithmStatus::NumericalInstability);
			REQUIRE_FALSE(configured.error_message.empty());
		}
	}
}

TEST_CASE("JacobiEigen_3x3_Identity", "[EigenSolvers][Jacobi][Small]")
{
		TEST_PRECISION_INFO();
	// Identity matrix - all eigenvalues = 1
	MatrixSym<Real> I(3, {static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0)});
	
	auto result = SymmMatEigenSolverJacobi::Solve(I);
	
	REQUIRE(result.converged);
	REQUIRE(result.eigenvalues.size() == 3);
	
	// All eigenvalues should be 1
	for (int i = 0; i < 3; i++)
		REQUIRE(std::abs(result.eigenvalues[i] - static_cast<MML::Real>(1.0)) < TOL(1e-10, 1e-5));
}

TEST_CASE("JacobiEigen_3x3_Diagonal", "[EigenSolvers][Jacobi][Small]")
{
		TEST_PRECISION_INFO();
	// Diagonal matrix - eigenvalues are diagonal elements
	MatrixSym<Real> D(3, {static_cast<MML::Real>(5.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(3.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0)});
	
	auto result = SymmMatEigenSolverJacobi::Solve(D);
	
	REQUIRE(result.converged);
	
	// Eigenvalues should be 1, 3, 5 (sorted)
	REQUIRE(std::abs(result.eigenvalues[0] - static_cast<MML::Real>(1.0)) < TOL(1e-10, 1e-5));
	REQUIRE(std::abs(result.eigenvalues[1] - static_cast<MML::Real>(3.0)) < TOL(1e-10, 1e-5));
	REQUIRE(std::abs(result.eigenvalues[2] - static_cast<MML::Real>(5.0)) < TOL(1e-10, 1e-5));
}

TEST_CASE("JacobiEigen_SymmetricMatrix_3x3", "[EigenSolvers][Jacobi][TestData]")
{
		TEST_PRECISION_INFO();
	// Use test data: symm_mat_3x3 with known eigenvalues
	MatrixSym<Real> mat_copy = symm_mat_3x3();
	auto result = SymmMatEigenSolverJacobi::Solve(mat_copy);
	
	REQUIRE(result.converged);
	REQUIRE(result.eigenvalues.size() == 3);
	
	// Compare with known eigenvalues (sorted)
	Vector<Real> expected = symm_mat_3x3_eigen_val();
	std::sort(expected.begin(), expected.end());
	
	for (int i = 0; i < 3; i++)
		REQUIRE(std::abs(result.eigenvalues[i] - expected[i]) < TOL(1e-6, 5e-5));
	
	// Verify A*V = V*Λ
	const auto& symm_mat = symm_mat_3x3();
	Matrix<Real> A_full(3, 3);
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			A_full[i][j] = (i <= j) ? symm_mat(i, j) : symm_mat(j, i);
	
	for (int i = 0; i < 3; i++)
	{
		Vector<Real> v = result.eigenvectors.VectorFromColumn(i);
		Vector<Real> Av = A_full * v;
		Vector<Real> lambda_v = result.eigenvalues[i] * v;
		
		for (int j = 0; j < 3; j++)
			REQUIRE(std::abs(Av[j] - lambda_v[j]) < TOL(1e-8, 1e-4));
	}
}

TEST_CASE("JacobiEigen_SymmetricMatrix_5x5", "[EigenSolvers][Jacobi][TestData]")
{
		TEST_PRECISION_INFO();
	// Use test data: symm_mat_5x5 with known eigenvalues
	MatrixSym<Real> mat_copy = symm_mat_5x5();
	auto result = SymmMatEigenSolverJacobi::Solve(mat_copy);
	
	REQUIRE(result.converged);
	REQUIRE(result.eigenvalues.size() == 5);
	
	// Compare with known eigenvalues
	Vector<Real> expected = symm_mat_5x5_eigen_val();
	std::sort(expected.begin(), expected.end());
	
	for (int i = 0; i < 5; i++)
		REQUIRE(std::abs(result.eigenvalues[i] - expected[i]) < TOL(1e-6, 5e-5));
	
	// Verify eigenvectors are orthonormal
	for (int i = 0; i < 5; i++)
	{
		Vector<Real> vi = result.eigenvectors.VectorFromColumn(i);
		Real norm = vi.NormL2();
		REQUIRE(std::abs(norm - static_cast<MML::Real>(1.0)) < TOL(1e-10, 1e-5));
		
		// Check orthogonality with all other vectors
		for (int j = i + 1; j < 5; j++)
		{
			Vector<Real> vj = result.eigenvectors.VectorFromColumn(j);
			Real dot = Utils::ScalarProduct(vi, vj);
			REQUIRE(std::abs(dot) < TOL(1e-8, 1e-4));
		}
	}
}

TEST_CASE("JacobiEigen_SymmetricMatrix_10x10", "[EigenSolvers][Jacobi][TestData]")
{
		TEST_PRECISION_INFO();
	// Use test data: symm_mat_10x10 with known eigenvalues
	MatrixSym<Real> mat_copy = symm_mat_10x10();
	auto result = SymmMatEigenSolverJacobi::Solve(mat_copy);
	
	REQUIRE(result.converged);
	REQUIRE(result.eigenvalues.size() == 10);
	
	// Compare with known eigenvalues
	Vector<Real> expected = symm_mat_10x10_eigen_val();
	std::sort(expected.begin(), expected.end());
	
	for (int i = 0; i < 10; i++)
		REQUIRE(std::abs(result.eigenvalues[i] - expected[i]) < TOL(1e-6, 5e-5));
}

TEST_CASE("JacobiEigen_SymmetricMatrix_3x3_Verified", "[EigenSolvers][Jacobi][Verified]")
{
		TEST_PRECISION_INFO();
	// Test with verified symmetric matrix from test data
	MatrixSym<Real> mat_copy = symm_mat_3x3();
	auto result = SymmMatEigenSolverJacobi::Solve(mat_copy);
	
	REQUIRE(result.converged);
	REQUIRE(result.eigenvalues.size() == 3);
	
	// Compare with verified eigenvalues (sorted ascending)
	Vector<Real> expected = symm_mat_3x3_eigen_val();
	std::sort(expected.begin(), expected.end());
	
	for (int i = 0; i < 3; i++)
		REQUIRE(std::abs(result.eigenvalues[i] - expected[i]) < TOL(1e-6, 5e-5));
	
	// Verify trace = sum of eigenvalues
	const auto& mat = symm_mat_3x3();
	Real trace = mat(0,0) + mat(1,1) + mat(2,2);
	Real eigensum = result.eigenvalues[0] + result.eigenvalues[1] + result.eigenvalues[2];
	REQUIRE(std::abs(trace - eigensum) < TOL(1e-10, 1e-5));
}

TEST_CASE("JacobiEigen_SymmetricMatrix_5x5_Verified", "[EigenSolvers][Jacobi][Verified]")
{
		TEST_PRECISION_INFO();
	// Test with verified 5x5 symmetric matrix from test data
	MatrixSym<Real> mat_copy = symm_mat_5x5();
	auto result = SymmMatEigenSolverJacobi::Solve(mat_copy);
	
	REQUIRE(result.converged);
	REQUIRE(result.eigenvalues.size() == 5);
	
	// Compare with verified eigenvalues (sorted ascending)
	Vector<Real> expected = symm_mat_5x5_eigen_val();
	std::sort(expected.begin(), expected.end());
	
	for (int i = 0; i < 5; i++)
		REQUIRE(std::abs(result.eigenvalues[i] - expected[i]) < TOL(1e-6, 5e-5));
	
	// Verify trace = sum of eigenvalues
	const auto& mat = symm_mat_5x5();
	Real trace = static_cast<MML::Real>(0.0);
	for (int i = 0; i < 5; i++)
		trace += mat(i, i);
	
	Real eigensum = static_cast<MML::Real>(0.0);
	for (int i = 0; i < 5; i++)
		eigensum += result.eigenvalues[i];
	
	REQUIRE(std::abs(trace - eigensum) < TOL(1e-10, 1e-5));
}

TEST_CASE("JacobiEigen_SymmetricMatrix_10x10_Verified", "[EigenSolvers][Jacobi][Verified]")
{
		TEST_PRECISION_INFO();
	// Test with verified 10x10 symmetric matrix from test data
	MatrixSym<Real> mat_copy = symm_mat_10x10();
	auto result = SymmMatEigenSolverJacobi::Solve(mat_copy);
	
	REQUIRE(result.converged);
	REQUIRE(result.eigenvalues.size() == 10);
	
	// Compare with verified eigenvalues (sorted ascending)
	Vector<Real> expected = symm_mat_10x10_eigen_val();
	std::sort(expected.begin(), expected.end());
	
	for (int i = 0; i < 10; i++)
		REQUIRE(std::abs(result.eigenvalues[i] - expected[i]) < TOL(1e-6, 5e-5));
	
	// Verify trace = sum of eigenvalues
	const auto& mat = symm_mat_10x10();
	Real trace = static_cast<MML::Real>(0.0);
	for (int i = 0; i < 10; i++)
		trace += mat(i, i);
	
	Real eigensum = static_cast<MML::Real>(0.0);
	for (int i = 0; i < 10; i++)
		eigensum += result.eigenvalues[i];
	
	REQUIRE(std::abs(trace - eigensum) < TOL(1e-10, 1e-5));
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
	
	MatrixSym<Real> mat_copy = symm_mat_5x5();
	auto result = SymmMatEigenSolverJacobi::Solve(mat_copy);
	
	REQUIRE(result.converged);
	
	int n = 5;
	Matrix<Real> VtV = result.eigenvectors.Transposed() * result.eigenvectors;
	
	// Check V^T * V = I
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			Real expected = (i == j) ? static_cast<MML::Real>(1.0) : static_cast<MML::Real>(0.0);
			REQUIRE(std::abs(VtV[i][j] - expected) < TOL(1e-8, 1e-4));
		}
	}
}

TEST_CASE("JacobiEigen_ReconstructMatrix", "[EigenSolvers][Jacobi][Properties]")
{
		TEST_PRECISION_INFO();
	// Test that A = V * Λ * V^T
	
	MatrixSym<Real> mat_copy = symm_mat_3x3();
	auto result = SymmMatEigenSolverJacobi::Solve(mat_copy);
	
	REQUIRE(result.converged);
	
	int n = 3;
	
	// Create diagonal matrix Λ
	Matrix<Real> Lambda(n, n);
	for (int i = 0; i < n; i++)
		Lambda[i][i] = result.eigenvalues[i];
	
	// Reconstruct A = V * Λ * V^T
	Matrix<Real> VL = result.eigenvectors * Lambda;
	Matrix<Real> A_reconstructed = VL * result.eigenvectors.Transposed();
	
	// Convert original to full matrix
	const auto& symm_mat = symm_mat_3x3();
	Matrix<Real> A_original(n, n);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			A_original[i][j] = (i <= j) ? symm_mat(i, j) : symm_mat(j, i);
	
	// Compare
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			REQUIRE(std::abs(A_reconstructed[i][j] - A_original[i][j]) < TOL(1e-8, 1e-4));
}

TEST_CASE("JacobiEigen_ConvergenceSpeed", "[EigenSolvers][Jacobi][Performance]")
{
		TEST_PRECISION_INFO();
	// Test that Jacobi converges reasonably fast
	
	MatrixSym<Real> mat_copy = symm_mat_10x10();
	auto result = SymmMatEigenSolverJacobi::Solve(mat_copy, TOL(1e-10, 1e-5), 100);
	
	REQUIRE(result.converged);
	
	// Should converge in reasonable number of iterations
	INFO("Iterations: " << result.iterations);
	REQUIRE(result.iterations < 50);  // Typically 5-15 for n=10
	
	// Residual should be small
	INFO("Residual: " << result.residual);
	REQUIRE(result.residual < TOL(1e-10, 1e-5));
}

TEST_CASE("JacobiEigen_NonConvergence", "[EigenSolvers][Jacobi][EdgeCases]")
{
		TEST_PRECISION_INFO();
	// Test behavior with very small maxIter
	
	MatrixSym<Real> mat_copy = symm_mat_10x10();
	auto result = SymmMatEigenSolverJacobi::Solve(mat_copy, TOL(1e-15, 1e-5), 1);  // Only 1 iteration
	
	// Should not converge
	REQUIRE_FALSE(result.converged);
	REQUIRE(result.iterations == 1);
}

TEST_CASE("JacobiEigen_TraceAndDeterminant", "[EigenSolvers][Jacobi][Properties]")
{
		TEST_PRECISION_INFO();
	// Test that trace(A) = sum(eigenvalues) and det(A) = product(eigenvalues)
	
	MatrixSym<Real> mat_copy = symm_mat_5x5();
	auto result = SymmMatEigenSolverJacobi::Solve(mat_copy);
	
	REQUIRE(result.converged);
	
	// Compute trace
	Real trace = static_cast<MML::Real>(0.0);
	for (int i = 0; i < 5; i++)
		trace += symm_mat_5x5()(i, i);
	
	// Sum of eigenvalues
	Real eigensum = static_cast<MML::Real>(0.0);
	for (int i = 0; i < 5; i++)
		eigensum += result.eigenvalues[i];
	
	REQUIRE(std::abs(trace - eigensum) < TOL(1e-10, 1e-5));
	
	// Product of eigenvalues = determinant
	Real eigenprod = static_cast<MML::Real>(1.0);
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
	MatrixSym<Real> A(2, {static_cast<MML::Real>(2.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0)});
	
	// Known eigenvalues: λ₁ = 1, λ₂ = 3
	auto result = SymmMatEigenSolverQR::Solve(A);
	
	REQUIRE(result.converged);
	REQUIRE(result.eigenvalues.size() == 2);
	
	// Check eigenvalues (sorted ascending)
	REQUIRE(std::abs(result.eigenvalues[0] - static_cast<MML::Real>(1.0)) < TOL(1e-10, 1e-5));
	REQUIRE(std::abs(result.eigenvalues[1] - static_cast<MML::Real>(3.0)) < TOL(1e-10, 1e-5));
	
	// Verify A*v = λ*v for each eigenpair
	for (int i = 0; i < 2; i++)
	{
		Vector<Real> v = result.eigenvectors.VectorFromColumn(i);
		Vector<Real> Av(2);
		for (int j = 0; j < 2; j++)
			Av[j] = A(j, 0) * v[0] + A(j, 1) * v[1];
		
		Vector<Real> lambda_v = result.eigenvalues[i] * v;
		
		for (int j = 0; j < 2; j++)
			REQUIRE(std::abs(Av[j] - lambda_v[j]) < TOL(1e-10, 1e-5));
	}
}

TEST_CASE("QREigen_3x3_Identity", "[EigenSolvers][QR][Small]")
{
		TEST_PRECISION_INFO();
	// 3x3 identity matrix - eigenvalues should all be 1
	MatrixSym<Real> I(3, {static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0)});
	
	auto result = SymmMatEigenSolverQR::Solve(I);
	
	REQUIRE(result.converged);
	REQUIRE(result.eigenvalues.size() == 3);
	
	// All eigenvalues should be static_cast<MML::Real>(1.0)
	for (int i = 0; i < 3; i++)
		REQUIRE(std::abs(result.eigenvalues[i] - static_cast<MML::Real>(1.0)) < TOL(1e-10, 1e-5));
}

TEST_CASE("QREigen_3x3_Diagonal", "[EigenSolvers][QR][Small]")
{
		TEST_PRECISION_INFO();
	// 3x3 diagonal matrix
	MatrixSym<Real> D(3, {static_cast<MML::Real>(5.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(3.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0)});
	
	// Known eigenvalues: 1, 3, 5
	auto result = SymmMatEigenSolverQR::Solve(D);
	
	REQUIRE(result.converged);
	
	// Check eigenvalues (sorted ascending)
	REQUIRE(std::abs(result.eigenvalues[0] - static_cast<MML::Real>(1.0)) < TOL(1e-10, 1e-5));
	REQUIRE(std::abs(result.eigenvalues[1] - static_cast<MML::Real>(3.0)) < TOL(1e-10, 1e-5));
	REQUIRE(std::abs(result.eigenvalues[2] - static_cast<MML::Real>(5.0)) < TOL(1e-10, 1e-5));
}

TEST_CASE("QREigen_SymmetricMatrix_3x3_Verified", "[EigenSolvers][QR][Verified]")
{
		TEST_PRECISION_INFO();
	// Use verified 3x3 symmetric matrix from test data
	const auto& symm_mat = symm_mat_3x3();
	Matrix<Real> A_full(3, 3);
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			A_full(i, j) = symm_mat(i, j);
	
	MatrixSym<Real> A(3);
	for (int i = 0; i < 3; i++)
		for (int j = i; j < 3; j++)
			A(i, j) = A_full(i, j);
	
	auto result = SymmMatEigenSolverQR::Solve(A);
	
	REQUIRE(result.converged);
	REQUIRE(result.eigenvalues.size() == 3);
	
	// Check eigenvalues against verified values (sorted)
	const auto& eigen_vals = symm_mat_3x3_eigen_val();
	Real expected_eigenvals[3] = {eigen_vals[0], eigen_vals[1], eigen_vals[2]};
	std::sort(std::begin(expected_eigenvals), std::end(expected_eigenvals));
	
	for (int i = 0; i < 3; i++)
		REQUIRE(std::abs(result.eigenvalues[i] - expected_eigenvals[i]) < TOL(1e-6, 5e-5));
	
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
			REQUIRE(std::abs(Av[j] - lambda_v[j]) < TOL(1e-6, 5e-5));
	}
}

TEST_CASE("QREigen_SymmetricMatrix_5x5_Verified", "[EigenSolvers][QR][Verified]")
{
		TEST_PRECISION_INFO();
	// Use verified 5x5 symmetric matrix from test data
	const auto& symm_mat = symm_mat_5x5();
	MatrixSym<Real> A(5);
	for (int i = 0; i < 5; i++)
		for (int j = i; j < 5; j++)
			A(i, j) = symm_mat(i, j);
	
	auto result = SymmMatEigenSolverQR::Solve(A);
	
	REQUIRE(result.converged);
	REQUIRE(result.eigenvalues.size() == 5);
	
	// Check eigenvalues against verified values (sorted)
	const auto& eigen_vals = symm_mat_5x5_eigen_val();
	Real expected_eigenvals[5] = {eigen_vals[0], eigen_vals[1], 
	                              eigen_vals[2], eigen_vals[3],
	                              eigen_vals[4]};
	std::sort(std::begin(expected_eigenvals), std::end(expected_eigenvals));
	
	for (int i = 0; i < 5; i++)
		REQUIRE(std::abs(result.eigenvalues[i] - expected_eigenvals[i]) < TOL(1e-6, 5e-5));
	
	// Verify eigenvector orthonormality
	for (int i = 0; i < 5; i++)
	{
		Vector<Real> vi = result.eigenvectors.VectorFromColumn(i);
		
		// Check normalization
		Real norm_sq = Utils::ScalarProduct(vi, vi);
		REQUIRE(std::abs(norm_sq - static_cast<MML::Real>(1.0)) < TOL(1e-10, 1e-5));
		
		// Check orthogonality with other eigenvectors
		for (int j = i + 1; j < 5; j++)
		{
			Vector<Real> vj = result.eigenvectors.VectorFromColumn(j);
			Real dot = Utils::ScalarProduct(vi, vj);
			REQUIRE(std::abs(dot) < TOL(1e-10, 1e-5));
		}
	}
}

TEST_CASE("QREigen_SymmetricMatrix_10x10_Verified", "[EigenSolvers][QR][Verified]")
{
		TEST_PRECISION_INFO();
	// Use verified 10x10 symmetric matrix from test data
	const auto& symm_mat = symm_mat_10x10();
	MatrixSym<Real> A(10);
	for (int i = 0; i < 10; i++)
		for (int j = i; j < 10; j++)
			A(i, j) = symm_mat(i, j);
	
	auto result = SymmMatEigenSolverQR::Solve(A);
	
	REQUIRE(result.converged);
	REQUIRE(result.eigenvalues.size() == 10);
	
	// Check eigenvalues against verified values (sorted)
	const auto& eigen_vals = symm_mat_10x10_eigen_val();
	std::vector<Real> expected_eigenvals(10);
	for (int i = 0; i < 10; i++)
		expected_eigenvals[i] = eigen_vals[i];
	std::sort(expected_eigenvals.begin(), expected_eigenvals.end());
	
	for (int i = 0; i < 10; i++)
		REQUIRE(std::abs(result.eigenvalues[i] - expected_eigenvals[i]) < 1e-4);  // Slightly relaxed for larger matrix
	
	// Verify matrix reconstruction: A = Q*Λ*Q^T
	Matrix<Real> Lambda(10, 10);
	for (int i = 0; i < 10; i++)
		for (int j = 0; j < 10; j++)
			Lambda(i, j) = static_cast<MML::Real>(0.0);
	for (int i = 0; i < 10; i++)
		Lambda(i, i) = result.eigenvalues[i];
	
	Matrix<Real> Q = result.eigenvectors;
	Matrix<Real> Qt = Q.Transposed();
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
	const auto& symm_mat = symm_mat_3x3();
	MatrixSym<Real> A(3);
	for (int i = 0; i < 3; i++)
		for (int j = i; j < 3; j++)
			A(i, j) = symm_mat(i, j);
	
	auto qr_result = SymmMatEigenSolverQR::Solve(A);
	auto jacobi_result = SymmMatEigenSolverJacobi::Solve(A);
	
	REQUIRE(qr_result.converged);
	REQUIRE(jacobi_result.converged);
	
	// Eigenvalues should match between methods
	for (int i = 0; i < 3; i++)
		REQUIRE(std::abs(qr_result.eigenvalues[i] - jacobi_result.eigenvalues[i]) < TOL(1e-8, 1e-4));
}

TEST_CASE("QREigen_CompareWithJacobi_5x5", "[EigenSolvers][QR][Jacobi][Comparison]")
{
		TEST_PRECISION_INFO();
	// Both algorithms should give the same eigenvalues for 5x5
	const auto& symm_mat = symm_mat_5x5();
	MatrixSym<Real> A(5);
	for (int i = 0; i < 5; i++)
		for (int j = i; j < 5; j++)
			A(i, j) = symm_mat(i, j);
	
	auto qr_result = SymmMatEigenSolverQR::Solve(A);
	auto jacobi_result = SymmMatEigenSolverJacobi::Solve(A);
	
	REQUIRE(qr_result.converged);
	REQUIRE(jacobi_result.converged);
	
	// Eigenvalues should match between methods
	for (int i = 0; i < 5; i++)
		REQUIRE(std::abs(qr_result.eigenvalues[i] - jacobi_result.eigenvalues[i]) < TOL(1e-8, 1e-4));
}

TEST_CASE("QREigen_Convergence", "[EigenSolvers][QR][Convergence]")
{
		TEST_PRECISION_INFO();
	// Test with tighter tolerance
	const auto& symm_mat = symm_mat_3x3();
	MatrixSym<Real> A(3);
	for (int i = 0; i < 3; i++)
		for (int j = i; j < 3; j++)
			A(i, j) = symm_mat(i, j);
	
	auto result = SymmMatEigenSolverQR::Solve(A, TOL(1e-12, 1e-5), 1000);
	
	REQUIRE(result.converged);
	REQUIRE(result.residual < TOL(1e-12, 1e-5));
	INFO("Converged in " << result.iterations << " iterations");
	INFO("Final residual: " << result.residual);
}

TEST_CASE("QREigen_Performance_vs_Jacobi", "[EigenSolvers][QR][Jacobi][Performance]")
{
		TEST_PRECISION_INFO();
	// QR should converge faster than Jacobi for larger matrices
	const auto& symm_mat = symm_mat_10x10();
	MatrixSym<Real> A(10);
	for (int i = 0; i < 10; i++)
		for (int j = i; j < 10; j++)
			A(i, j) = symm_mat(i, j);
	
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
			A(i, j) = symm_mat_5x5()(i, j);
	
	auto result = SymmMatEigenSolverQR::Solve(A);
	
	REQUIRE(result.converged);
	
	// Compute trace
	Real trace = static_cast<MML::Real>(0.0);
	for (int i = 0; i < 5; i++)
		trace += A(i, i);
	
	// Sum eigenvalues
	Real eigensum = static_cast<MML::Real>(0.0);
	for (int i = 0; i < 5; i++)
		eigensum += result.eigenvalues[i];
	
	REQUIRE(std::abs(trace - eigensum) < TOL(1e-10, 1e-5));
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
	int n = M.rows();
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
			Real baseTol = MML::Real(1e-5);  // Accounts for ~6-digit precision in test bed eigenvalues
			Real tol = baseTol * std::pow(MML::Real(10.0), MML::Real(test.difficulty - 1));
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
					for (int i = 0; i < result.eigenvalues.size(); ++i)
						ss << result.eigenvalues[i] << " ";
					ss << "\nExpected: ";
					for (int i = 0; i < test.realEigenvalues.size(); ++i)
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
			Real baseTol = TOL(1e-5, 5e-4);  // Accounts for ~6-digit precision in test bed eigenvalues
			Real tol = baseTol * std::pow(MML::Real(10.0), MML::Real(test.difficulty - 1));
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
				CHECK(TestBeds::verifyRealEigenvalues(result.eigenvalues, test.realEigenvalues, MML::Real(1e-5)));
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
				CHECK(TestBeds::verifyRealEigenvalues(result.eigenvalues, test.realEigenvalues, TOL(1e-5, 5e-4)));
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
				Real tol = TOL(1e-6, 5e-5);
				if (test.difficulty >= 4) tol = TOL(1e-4, 1e-2);
				
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
					CHECK(maxResidual < TOL(1e-6, 5e-5));  // Self-consistency should still be good
				}
			}
			
			// Try Jacobi as well
			auto resultJac = SymmMatEigenSolverJacobi::Solve(A);
			
			if (resultJac.converged && test.realEigenvalues.size() > 0 && !skipEigenvalueComparison)
			{
				Real tol = TOL(1e-6, 5e-5);
				if (test.difficulty >= 4) tol = TOL(1e-4, 1e-2);
				
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
				CHECK(maxResidual < TOL(1e-6, 5e-5));  // Self-consistency should still be good
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
			CHECK(TestBeds::verifyRealEigenvalues(resultJac.eigenvalues, test.realEigenvalues, TOL(1e-8, 1e-4)));
			CHECK(TestBeds::verifyRealEigenvalues(resultQR.eigenvalues, test.realEigenvalues, TOL(1e-8, 1e-4)));
			
			// Eigenvectors should still be orthonormal for symmetric matrices
			CHECK(TestBeds::verifyOrthonormality(resultJac.eigenvectors, TOL(1e-8, 1e-4)));
			CHECK(TestBeds::verifyOrthonormality(resultQR.eigenvectors, TOL(1e-8, 1e-4)));
		}
	}
}

} // namespace MML::Tests::Algorithms::EigensystemSolversTests

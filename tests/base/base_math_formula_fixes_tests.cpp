///////////////////////////////////////////////////////////////////////////////////////////
///                         Base Math Formula Fixes Tests                             ///
///                                                                                   ///
///  Tests for P0: Fix incorrect base math formulas (MinimalMathLibrary-om8d.4)       ///
///  - VectorN::NormL2() complex handling                                             ///
///  - Matrix::NormL2() complex handling                                              ///
///  - DiracExp formula (correct /2 factor)                                           ///
///  - DiracSin singularity at x=0                                                    ///
///  - Random::UniformVecDirection3 uniform sphere sampling                           ///
///  - Vector3Cartesian::AngleToVector acos clamping                                  ///
///////////////////////////////////////////////////////////////////////////////////////////

#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"
#include <cmath>

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#include "base/Vector/VectorN.h"
#include "base/Matrix/Matrix.h"
#include "base/Matrix/MatrixSym.h"
#include "base/Matrix/MatrixTriDiag.h"
#include "base/DiracDeltaFunction.h"
#include "base/Random.h"
#include "base/Vector/VectorTypes.h"
#endif

using namespace MML;
using namespace MML::Testing;
using namespace Catch::Matchers;

namespace MML::Tests::Base::BaseMathFormulaFixesTests
{
	/*********************************************************************/
	/*****          VectorN::NormL2 Complex Tests                   *****/
	/*********************************************************************/
	TEST_CASE("VectorN_NormL2_Complex_correct_magnitude", "[VectorN][Complex][P0]")
	{
		TEST_PRECISION_INFO();
		// A complex vector with known norm
		// v = (3+4i, 0, 0) should have norm 5
		VectorN<Complex, 3> v{Complex(3.0, 4.0), Complex(0.0, 0.0), Complex(0.0, 0.0)};
		
		Real norm = v.NormL2();
		REQUIRE_THAT(norm, WithinAbs(REAL(5.0), TOL(1e-10, 1e-5)));
	}
	
	TEST_CASE("VectorN_NormL2_Complex_multiple_components", "[VectorN][Complex][P0]")
	{
		TEST_PRECISION_INFO();
		// v = (1+i, 1+i, 1+i) 
		// |1+i|² = 2, so ||v||² = 3*2 = 6, ||v|| = √6
		VectorN<Complex, 3> v{Complex(1.0, 1.0), Complex(1.0, 1.0), Complex(1.0, 1.0)};
		
		Real norm = v.NormL2();
		Real expected = std::sqrt(6.0);
		REQUIRE_THAT(norm, WithinAbs(expected, TOL(1e-10, 1e-5)));
	}
	
	TEST_CASE("VectorN_NormL2_Complex_pure_imaginary", "[VectorN][Complex][P0]")
	{
		TEST_PRECISION_INFO();
		// v = (0+3i, 0+4i) should have norm 5
		// |3i|² = 9, |4i|² = 16, total = 25, norm = 5
		VectorN<Complex, 2> v{Complex(0.0, 3.0), Complex(0.0, 4.0)};
		
		Real norm = v.NormL2();
		REQUIRE_THAT(norm, WithinAbs(REAL(5.0), TOL(1e-10, 1e-5)));
	}
	
	TEST_CASE("VectorN_NormL2_Real_unchanged", "[VectorN][P0]")
	{
		TEST_PRECISION_INFO();
		// Real vector should still work: v = (3, 4) -> norm = 5
		VectorN<Real, 2> v{REAL(3.0), REAL(4.0)};
		
		Real norm = v.NormL2();
		REQUIRE_THAT(norm, WithinAbs(REAL(5.0), TOL(1e-10, 1e-5)));
	}

	/*********************************************************************/
	/*****          Matrix::NormL2 Complex Tests                    *****/
	/*********************************************************************/
	TEST_CASE("Matrix_NormL2_Complex_correct_frobenius", "[Matrix][Complex][P0]")
	{
		TEST_PRECISION_INFO();
		// 2x2 matrix with one complex element
		// M = [[3+4i, 0], [0, 0]]
		// ||M||_F = |3+4i| = 5
		Matrix<Complex> M(2, 2);
		M(0, 0) = Complex(3.0, 4.0);
		M(0, 1) = Complex(0.0, 0.0);
		M(1, 0) = Complex(0.0, 0.0);
		M(1, 1) = Complex(0.0, 0.0);
		
		Real norm = M.NormL2();
		REQUIRE_THAT(norm, WithinAbs(REAL(5.0), TOL(1e-10, 1e-5)));
	}
	
	TEST_CASE("Matrix_NormL2_Complex_full_matrix", "[Matrix][Complex][P0]")
	{
		TEST_PRECISION_INFO();
		// M = [[1+i, 1+i], [1+i, 1+i]]
		// |1+i|² = 2, 4 elements, total = 8, norm = √8 = 2√2
		Matrix<Complex> M(2, 2);
		M(0, 0) = Complex(1.0, 1.0);
		M(0, 1) = Complex(1.0, 1.0);
		M(1, 0) = Complex(1.0, 1.0);
		M(1, 1) = Complex(1.0, 1.0);
		
		Real norm = M.NormL2();
		Real expected = std::sqrt(8.0);
		REQUIRE_THAT(norm, WithinAbs(expected, TOL(1e-10, 1e-5)));
	}
	
	TEST_CASE("Matrix_NormL2_Real_unchanged", "[Matrix][P0]")
	{
		TEST_PRECISION_INFO();
		// Real matrix: M = [[3, 0], [4, 0]] -> ||M||_F = 5
		Matrix<Real> M(2, 2);
		M(0, 0) = REAL(3.0);
		M(0, 1) = REAL(0.0);
		M(1, 0) = REAL(4.0);
		M(1, 1) = REAL(0.0);
		
		Real norm = M.NormL2();
		REQUIRE_THAT(norm, WithinAbs(REAL(5.0), TOL(1e-10, 1e-5)));
	}

	/*********************************************************************/
	/*****          DiracExp Integration Tests                      *****/
	/*********************************************************************/
	TEST_CASE("DiracExp_integrates_to_one", "[DiracDeltaFunction][DiracExp][P0]")
	{
		TEST_PRECISION_INFO();
		// The corrected DiracExp should integrate to 1
		// ∫ N/√(2π) exp(-N²x²/2) dx = 1
		DiracExp delta(10);
		
		// Numerical integration using trapezoidal rule
		Real sum = 0.0;
		Real dx = REAL(0.001);
		Real range = REAL(5.0);  // Well beyond the peak
		
		for (Real x = -range; x <= range; x += dx) {
			sum += delta(x) * dx;
		}
		
		// Should be very close to 1
		REQUIRE_THAT(sum, WithinAbs(REAL(1.0), REAL(0.01)));
	}
	
	TEST_CASE("DiracExp_formula_at_origin", "[DiracDeltaFunction][DiracExp][P0]")
	{
		TEST_PRECISION_INFO();
		// At x=0: δ(0) = N/√(2π) * exp(0) = N/√(2π)
		DiracExp delta(10);
		Real expected = REAL(10.0) / std::sqrt(2 * Constants::PI);
		
		REQUIRE_THAT(delta(REAL(0.0)), WithinAbs(expected, TOL(1e-10, 1e-5)));
	}

	/*********************************************************************/
	/*****          DiracSin Singularity Tests                      *****/
	/*********************************************************************/
	TEST_CASE("DiracSin_handles_x_equals_zero", "[DiracDeltaFunction][DiracSin][P0]")
	{
		TEST_PRECISION_INFO();
		// At x=0, the limit is N/π (previously would return NaN or inf)
		DiracSin delta(10);
		Real expected = REAL(10.0) / Constants::PI;
		
		// This should NOT throw or return NaN
		Real val = delta(REAL(0.0));
		REQUIRE_FALSE(std::isnan(val));
		REQUIRE_FALSE(std::isinf(val));
		REQUIRE_THAT(val, WithinAbs(expected, TOL(1e-10, 1e-5)));
	}
	
	TEST_CASE("DiracSin_handles_very_small_x", "[DiracDeltaFunction][DiracSin][P0]")
	{
		TEST_PRECISION_INFO();
		DiracSin delta(10);
		Real expected = REAL(10.0) / Constants::PI;
		
		// Very small x values should approach N/π
		REQUIRE_THAT(delta(REAL(1e-16)), WithinAbs(expected, TOL(1e-8, 1e-4)));
		REQUIRE_THAT(delta(-REAL(1e-16)), WithinAbs(expected, TOL(1e-8, 1e-4)));
	}
	
	TEST_CASE("DiracSin_continuous_near_origin", "[DiracDeltaFunction][DiracSin][P0]")
	{
		TEST_PRECISION_INFO();
		DiracSin delta(10);
		
		// Should be continuous: limit from both sides equals value at 0
		Real valAt0 = delta(REAL(0.0));
		Real valNearPos = delta(REAL(0.0001));
		Real valNearNeg = delta(-REAL(0.0001));
		
		// All should be close to N/π
		REQUIRE_THAT(valAt0, WithinAbs(valNearPos, REAL(0.01)));
		REQUIRE_THAT(valAt0, WithinAbs(valNearNeg, REAL(0.01)));
	}

	/*********************************************************************/
	/*****          Random::UniformVecDirection3 Tests              *****/
	/*********************************************************************/
	TEST_CASE("UniformVecDirection3_returns_correct_magnitude", "[Random][P0]")
	{
		TEST_PRECISION_INFO();
		Real vx, vy, vz;
		Real abs = REAL(5.0);
		
		for (int i = 0; i < 100; ++i) {
			Random::UniformVecDirection3(vx, vy, vz, abs);
			Real mag = std::sqrt(vx*vx + vy*vy + vz*vz);
			REQUIRE_THAT(mag, WithinAbs(abs, TOL(1e-10, 1e-5)));
		}
	}
	
	TEST_CASE("UniformVecDirection3_z_distribution_uniform", "[Random][P0]")
	{
		TEST_PRECISION_INFO();
		// For uniform sphere sampling, z = cos(θ) should be uniform in [-1, 1]
		// This is the key property that distinguishes correct from incorrect sampling
		Real vx, vy, vz;
		Real abs = REAL(1.0);
		
		// Count samples in different z regions
		int countPosZ = 0;
		int countNegZ = 0;
		int numSamples = 10000;
		
		for (int i = 0; i < numSamples; ++i) {
			Random::UniformVecDirection3(vx, vy, vz, abs);
			if (vz > 0) countPosZ++;
			else countNegZ++;
		}
		
		// Should be roughly equal (within statistical noise)
		Real ratio = static_cast<Real>(countPosZ) / numSamples;
		REQUIRE_THAT(ratio, WithinAbs(REAL(0.5), REAL(0.05)));  // 5% tolerance
	}
	
	TEST_CASE("UniformVecDirection3_polar_angle_distribution", "[Random][P0]")
	{
		TEST_PRECISION_INFO();
		// The old incorrect method would cluster points near the poles
		// The correct method should have uniform cos(θ) distribution
		Real vx, vy, vz;
		Real abs = REAL(1.0);
		
		// Count samples in equatorial band vs polar caps
		// Equatorial band: |z| < 0.5 (covers same solid angle as two polar caps |z| > 0.5)
		int countEquatorial = 0;
		int countPolar = 0;
		int numSamples = 10000;
		
		for (int i = 0; i < numSamples; ++i) {
			Random::UniformVecDirection3(vx, vy, vz, abs);
			if (std::abs(vz) < 0.5) countEquatorial++;
			else countPolar++;
		}
		
		// With uniform sampling, equatorial and polar should be equal
		Real ratio = static_cast<Real>(countEquatorial) / numSamples;
		REQUIRE_THAT(ratio, WithinAbs(REAL(0.5), REAL(0.05)));
	}

	/*********************************************************************/
	/*****          Vector3Cartesian::AngleToVector Tests           *****/
	/*********************************************************************/
	TEST_CASE("AngleToVector_parallel_vectors", "[Vector3Cartesian][P0]")
	{
		TEST_PRECISION_INFO();
		Vector3Cartesian a(REAL(1.0), REAL(0.0), REAL(0.0));
		Vector3Cartesian b(REAL(2.0), REAL(0.0), REAL(0.0));  // Parallel
		
		Real angle = a.AngleToVector(b);
		REQUIRE_THAT(angle, WithinAbs(REAL(0.0), TOL(1e-10, 1e-5)));
	}
	
	TEST_CASE("AngleToVector_antiparallel_vectors", "[Vector3Cartesian][P0]")
	{
		TEST_PRECISION_INFO();
		Vector3Cartesian a(REAL(1.0), REAL(0.0), REAL(0.0));
		Vector3Cartesian b(-REAL(1.0), REAL(0.0), REAL(0.0));  // Antiparallel
		
		Real angle = a.AngleToVector(b);
		REQUIRE_THAT(angle, WithinAbs(Constants::PI, TOL(1e-10, 1e-5)));
	}
	
	TEST_CASE("AngleToVector_perpendicular_vectors", "[Vector3Cartesian][P0]")
	{
		TEST_PRECISION_INFO();
		Vector3Cartesian a(REAL(1.0), REAL(0.0), REAL(0.0));
		Vector3Cartesian b(REAL(0.0), REAL(1.0), REAL(0.0));  // Perpendicular
		
		Real angle = a.AngleToVector(b);
		REQUIRE_THAT(angle, WithinAbs(Constants::PI / 2, TOL(1e-10, 1e-5)));
	}
	
	TEST_CASE("AngleToVector_near_boundary_no_NaN", "[Vector3Cartesian][P0]")
	{
		TEST_PRECISION_INFO();
		// Test case that could produce cos > 1 due to floating-point errors
		// if acos is not clamped
		Vector3Cartesian a(REAL(1.0), REAL(0.0), REAL(0.0));
		
		// Create a vector that's numerically almost identical
		// Floating-point arithmetic might produce dot product > 1
		Vector3Cartesian b(REAL(1.0) + TOL(1e-15, 1e-5), REAL(0.0), REAL(0.0));
		b = b.GetAsUnitVector();  // Normalize should make it unit
		
		Real angle = a.AngleToVector(b);
		
		// Should NOT be NaN
		REQUIRE_FALSE(std::isnan(angle));
		// Should be close to 0
		REQUIRE_THAT(angle, WithinAbs(REAL(0.0), REAL(1e-6)));
	}
	
	TEST_CASE("AngleToVector_large_vectors_no_overflow", "[Vector3Cartesian][P0]")
	{
		TEST_PRECISION_INFO();
		// Large magnitude vectors should still work correctly
		Real large = REAL(1e10);
		Vector3Cartesian a(large, REAL(0.0), REAL(0.0));
		Vector3Cartesian b(REAL(0.0), large, REAL(0.0));
		
		Real angle = a.AngleToVector(b);
		
		REQUIRE_FALSE(std::isnan(angle));
		REQUIRE_THAT(angle, WithinAbs(Constants::PI / 2, REAL(1e-6)));
	}

	/*********************************************************************/
	/*****     Matrix::IsDiagDominant Complex Tests (P0 om8d.2)     *****/
	/*********************************************************************/
	TEST_CASE("Matrix_IsDiagDominant_Complex_correct_magnitude", "[Matrix][Complex][P0]")
	{
		TEST_PRECISION_INFO();
		// Diagonally dominant complex matrix
		// Diagonal: (5+0i), off-diagonal sum: |1+1i| + |1-1i| = √2 + √2 ≈ 2.83
		// 5 > 2.83, so diagonally dominant
		Matrix<Complex> M(2, 2);
		M(0, 0) = Complex(5.0, 0.0);
		M(0, 1) = Complex(1.0, 1.0);   // |1+i| = √2 ≈ 1.41
		M(1, 0) = Complex(1.0, -1.0);  // |1-i| = √2 ≈ 1.41
		M(1, 1) = Complex(5.0, 0.0);
		
		REQUIRE(M.isDiagonallyDominant());
	}
	
	TEST_CASE("Matrix_IsDiagDominant_Complex_not_dominant", "[Matrix][Complex][P0]")
	{
		TEST_PRECISION_INFO();
		// Not diagonally dominant: diagonal smaller than off-diagonal sum
		// |1+0i| = 1, |2+2i| = √8 ≈ 2.83 > 1
		Matrix<Complex> M(2, 2);
		M(0, 0) = Complex(1.0, 0.0);
		M(0, 1) = Complex(2.0, 2.0);   // |2+2i| = √8 ≈ 2.83
		M(1, 0) = Complex(0.0, 0.0);
		M(1, 1) = Complex(1.0, 0.0);
		
		REQUIRE_FALSE(M.isDiagonallyDominant());
	}

	/*********************************************************************/
	/*****     TridiagonalMatrix::Solve Complex Tests (P0 om8d.2)   *****/
	/*********************************************************************/
	TEST_CASE("TridiagonalMatrix_Solve_Complex_basic", "[MatrixTriDiag][Complex][P0]")
	{
		TEST_PRECISION_INFO();
		// Simple complex tridiagonal system
		// Diagonal-dominant to ensure stability
		Vector<Complex> below{Complex(0.0, 0.0), Complex(0.5, 0.0), Complex(0.5, 0.0)};
		Vector<Complex> diag{Complex(2.0, 0.0), Complex(2.0, 0.0), Complex(2.0, 0.0)};
		Vector<Complex> above{Complex(0.5, 0.0), Complex(0.5, 0.0), Complex(0.0, 0.0)};
		
		TridiagonalMatrix<Complex> T(3, below, diag, above);
		
		// RHS with complex values
		Vector<Complex> rhs{Complex(1.0, 1.0), Complex(2.0, 0.0), Complex(1.0, -1.0)};
		
		Vector<Complex> sol = T.Solve(rhs);
		
		// Verify: T * sol ≈ rhs
		// Row 0: diag[0]*sol[0] + above[0]*sol[1]
		Complex row0 = diag[0] * sol[0] + above[0] * sol[1];
		Complex row1 = below[1] * sol[0] + diag[1] * sol[1] + above[1] * sol[2];
		Complex row2 = below[2] * sol[1] + diag[2] * sol[2];
		
		REQUIRE_THAT(std::abs(row0 - rhs[0]), WithinAbs(REAL(0.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(std::abs(row1 - rhs[1]), WithinAbs(REAL(0.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(std::abs(row2 - rhs[2]), WithinAbs(REAL(0.0), TOL(1e-10, 1e-5)));
	}
	
	TEST_CASE("TridiagonalMatrix_Solve_throws_on_singular", "[MatrixTriDiag][P0]")
	{
		TEST_PRECISION_INFO();
		// Singular system (zero on diagonal)
		Vector<Real> below{REAL(0.0), REAL(1.0), REAL(1.0)};
		Vector<Real> diag{REAL(0.0), REAL(2.0), REAL(2.0)};  // Zero at [0]!
		Vector<Real> above{REAL(1.0), REAL(1.0), REAL(0.0)};
		
		TridiagonalMatrix<Real> T(3, below, diag, above);
		Vector<Real> rhs{REAL(1.0), REAL(2.0), REAL(3.0)};
		
		REQUIRE_THROWS_AS(T.Solve(rhs), SingularMatrixError);
	}

	/*********************************************************************/
	/*****     MatrixSym::Norm Complex Tests (P0 om8d.2)            *****/
	/*********************************************************************/
	TEST_CASE("MatrixSym_NormFrobenius_Complex_returns_Real", "[MatrixSym][Complex][P0]")
	{
		TEST_PRECISION_INFO();
		// 2x2 symmetric complex matrix
		// Elements: (1+i), (2+0i), (1+i), (3+i)  (symmetric: (0,1) = (1,0))
		MatrixSym<Complex> M(2);
		M(0, 0) = Complex(1.0, 1.0);   // |1+i|² = 2
		M(0, 1) = Complex(2.0, 0.0);   // |2+0i|² = 4, appears twice (off-diagonal)
		M(1, 1) = Complex(3.0, 1.0);   // |3+i|² = 10
		
		Real norm = M.NormFrobenius();
		
		// Frobenius: sqrt(|a00|² + 2*|a01|² + |a11|²) = sqrt(2 + 2*4 + 10) = sqrt(20)
		Real expected = std::sqrt(REAL(2.0) + REAL(2.0) * REAL(4.0) + REAL(10.0));
		
		REQUIRE_THAT(norm, WithinAbs(expected, TOL(1e-10, 1e-5)));
	}
	
	TEST_CASE("MatrixSym_NormInf_Complex_returns_Real", "[MatrixSym][Complex][P0]")
	{
		TEST_PRECISION_INFO();
		// 2x2 symmetric complex matrix
		MatrixSym<Complex> M(2);
		M(0, 0) = Complex(1.0, 0.0);   // |1| = 1
		M(0, 1) = Complex(3.0, 4.0);   // |3+4i| = 5
		M(1, 1) = Complex(2.0, 0.0);   // |2| = 2
		
		Real norm = M.NormInf();
		
		// Row 0: |a00| + |a01| = 1 + 5 = 6
		// Row 1: |a10| + |a11| = 5 + 2 = 7 (max)
		REQUIRE_THAT(norm, WithinAbs(REAL(7.0), TOL(1e-10, 1e-5)));
	}
	
	TEST_CASE("MatrixSym_Norm1_Complex_returns_Real", "[MatrixSym][Complex][P0]")
	{
		TEST_PRECISION_INFO();
		// Norm1 = NormInf for symmetric matrices
		MatrixSym<Complex> M(2);
		M(0, 0) = Complex(1.0, 0.0);
		M(0, 1) = Complex(3.0, 4.0);
		M(1, 1) = Complex(2.0, 0.0);
		
		Real norm1 = M.Norm1();
		Real normInf = M.NormInf();
		
		// For symmetric matrices, Norm1 == NormInf
		REQUIRE_THAT(norm1, WithinAbs(normInf, TOL(1e-10, 1e-5)));
	}
}


#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/algorithms/MatrixAlg.h>
#endif

#include "../../test_data/linear_alg_eq_systems_complex_defs.h"

using namespace MML;
using namespace MML::Testing;

namespace MML::Tests::Core::ComplexSVDTests
{
	Matrix<Complex> Adjoint(const Matrix<Complex>& matrix)
	{
		Matrix<Complex> result(matrix.cols(), matrix.rows());
		for (int row = 0; row < matrix.rows(); ++row)
			for (int col = 0; col < matrix.cols(); ++col)
				result(col, row) = std::conj(matrix(row, col));
		return result;
	}

	Matrix<Complex> Sigma(int rows, int cols, const Vector<Real>& singularValues)
	{
		Matrix<Complex> result(rows, cols);
		for (int index = 0; index < singularValues.size(); ++index)
			result(index, index) = Complex{singularValues[index], REAL(0.0)};
		return result;
	}

	void RequireUnitary(const Matrix<Complex>& matrix, Real tolerance)
	{
		REQUIRE((Adjoint(matrix) * matrix).IsEqualTo(Matrix<Complex>::Identity(matrix.cols()), tolerance));
	}

	void RequireSVD(const Matrix<Complex>& matrix, int expectedRank, Real tolerance)
	{
		const auto svd = MatrixAlg::SVDDecompose(matrix);
		REQUIRE(svd.U.rows() == matrix.rows());
		REQUIRE(svd.U.cols() == matrix.rows());
		REQUIRE(svd.V.rows() == matrix.cols());
		REQUIRE(svd.V.cols() == matrix.cols());
		REQUIRE(svd.singularValues.size() == std::min(matrix.rows(), matrix.cols()));
		REQUIRE(svd.rank == expectedRank);
		for (int index = 0; index + 1 < svd.singularValues.size(); ++index)
			REQUIRE(svd.singularValues[index] >= svd.singularValues[index + 1]);
		RequireUnitary(svd.U, tolerance);
		RequireUnitary(svd.V, tolerance);
		REQUIRE((svd.U * Sigma(matrix.rows(), matrix.cols(), svd.singularValues) * Adjoint(svd.V))
			.IsEqualTo(matrix, tolerance));
	}

	TEST_CASE("Complex SVD matches structured fixture singular values", "[ComplexSVD][KnownValues][MML2]")
	{
		const auto& matrix = TestBeds::mat_cmplx_hermitian_2x2();
		const auto& expected = TestBeds::mat_cmplx_hermitian_2x2_singular_values();
		const auto svd = MatrixAlg::SVDDecompose(matrix);

		REQUIRE(svd.singularValues.size() == expected.size());
		for (int index = 0; index < expected.size(); ++index)
			REQUIRE(std::abs(svd.singularValues[index] - expected[index]) < TOL(1e-10, 1e-4));
		REQUIRE(std::abs(MatrixAlg::ConditionNumber(matrix) - TestBeds::mat_cmplx_hermitian_2x2_cond_2()) < TOL(1e-9, 1e-4));
	}

	TEST_CASE("Complex SVD reconstructs square tall and wide matrices", "[ComplexSVD][Reconstruction][MML2]")
	{
		RequireSVD(TestBeds::mat_cmplx_1_3x3(), 3, TOL(1e-9, 1e-4));
		RequireSVD(Matrix<Complex>{4, 2, {
			Complex{1.0, 1.0}, Complex{2.0, -1.0},
			Complex{0.0, 2.0}, Complex{1.0, 0.0},
			Complex{3.0, -1.0}, Complex{0.0, 1.0},
			Complex{2.0, 0.5}, Complex{-1.0, 2.0}
		}}, 2, TOL(1e-8, 1e-4));
		RequireSVD(Matrix<Complex>{2, 4, {
			Complex{1.0, 1.0}, Complex{0.0, 2.0}, Complex{3.0, -1.0}, Complex{2.0, 0.5},
			Complex{2.0, -1.0}, Complex{1.0, 0.0}, Complex{0.0, 1.0}, Complex{-1.0, 2.0}
		}}, 2, TOL(1e-8, 1e-4));
	}

	TEST_CASE("Complex SVD handles rank deficiency and repeated singular values", "[ComplexSVD][RankDeficient][MML2]")
	{
		const Matrix<Complex> rankOne{3, 2, {
			Complex{1.0, 1.0}, Complex{2.0, -1.0},
			Complex{2.0, 2.0}, Complex{4.0, -2.0},
			Complex{-1.0, -1.0}, Complex{-2.0, 1.0}
		}};
		RequireSVD(rankOne, 1, TOL(1e-8, 1e-4));
		const auto spaces = MatrixAlg::FundamentalSubspacesOf(rankOne);
		REQUIRE(spaces.nullSpace.rows() == 2);
		REQUIRE(spaces.nullSpace.cols() == 1);
		REQUIRE(spaces.leftNullSpace.rows() == 3);
		REQUIRE(spaces.leftNullSpace.cols() == 2);
		REQUIRE(MatrixAlg::FrobeniusNorm(rankOne * spaces.nullSpace) < TOL(1e-8, 1e-4));
		REQUIRE(MatrixAlg::FrobeniusNorm(Adjoint(rankOne) * spaces.leftNullSpace) < TOL(1e-8, 1e-4));

		const Real invSqrt2 = REAL(1.0) / std::sqrt(REAL(2.0));
		const Matrix<Complex> repeated{2, 2, {
			Complex{invSqrt2, 0.0}, Complex{0.0, invSqrt2},
			Complex{0.0, invSqrt2}, Complex{invSqrt2, 0.0}
		}};
		RequireSVD(repeated, 2, TOL(1e-9, 1e-4));
	}

	TEST_CASE("Complex pseudoinverse satisfies Moore Penrose identities", "[ComplexSVD][PseudoInverse][MML2]")
	{
		const Matrix<Complex> matrix{2, 3, {
			Complex{1.0, 1.0}, Complex{2.0, 0.0}, Complex{0.0, -1.0},
			Complex{2.0, 2.0}, Complex{4.0, 0.0}, Complex{0.0, -2.0}
		}};
		const auto pseudoinverse = MatrixAlg::PseudoInverse(matrix);
		REQUIRE((matrix * pseudoinverse * matrix).IsEqualTo(matrix, TOL(1e-8, 1e-4)));
		REQUIRE((pseudoinverse * matrix * pseudoinverse).IsEqualTo(pseudoinverse, TOL(1e-8, 1e-4)));
		REQUIRE((matrix * pseudoinverse).IsEqualTo(Adjoint(matrix * pseudoinverse), TOL(1e-8, 1e-4)));
		REQUIRE((pseudoinverse * matrix).IsEqualTo(Adjoint(pseudoinverse * matrix), TOL(1e-8, 1e-4)));
	}

	TEST_CASE("Complex SVD is scale aware", "[ComplexSVD][Scaling][MML2]")
	{
		const Matrix<Complex> matrix{2, 2, {
			Complex{2.0, 1.0}, Complex{1.0, -2.0},
			Complex{-1.0, 0.5}, Complex{3.0, 1.0}
		}};
		const auto baseline = MatrixAlg::SVDDecompose(matrix);
		const Real scale = std::sqrt(std::numeric_limits<Real>::max()) * REAL(0.01);
		const auto scaled = MatrixAlg::SVDDecompose(matrix * Complex{scale, 0.0});
		for (int index = 0; index < baseline.singularValues.size(); ++index)
			REQUIRE(std::abs(scaled.singularValues[index] / scale - baseline.singularValues[index]) < TOL(1e-8, 1e-4));
	}
}

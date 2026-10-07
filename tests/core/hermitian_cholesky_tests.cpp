#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/algorithms/MatrixAlg.h>
#endif

using namespace MML;
using namespace MML::Testing;

namespace MML::Tests::Core::HermitianCholeskyTests
{
	Matrix<Complex> Adjoint(const Matrix<Complex>& matrix)
	{
		Matrix<Complex> result(matrix.cols(), matrix.rows());
		for (int row = 0; row < matrix.rows(); ++row)
			for (int col = 0; col < matrix.cols(); ++col)
				result(col, row) = std::conj(matrix(row, col));
		return result;
	}

	TEST_CASE("Complex Cholesky reconstructs and solves Hermitian positive-definite systems",
				  "[HermitianCholesky][Cholesky][MML2]")
	{
		const Matrix<Complex> expectedLower{3, 3, {
			Complex{2.0, 0.0}, Complex{}, Complex{},
			Complex{1.0, 1.0}, Complex{3.0, 0.0}, Complex{},
			Complex{-0.5, 2.0}, Complex{0.25, -1.0}, Complex{1.5, 0.0}
		}};
		const Matrix<Complex> matrix = expectedLower * Adjoint(expectedLower);
		const Vector<Complex> expectedSolution{
			Complex{1.0, -1.0}, Complex{-2.0, 0.5}, Complex{0.25, 2.0}
		};

		CholeskySolver<Complex> solver(matrix);
		REQUIRE((solver.L() * Adjoint(solver.L())).IsEqualTo(matrix, TOL(1e-10, 1e-4)));
		REQUIRE(solver.Solve(matrix * expectedSolution).IsEqualTo(expectedSolution, TOL(1e-10, 1e-4)));

		Matrix<Complex> inverse;
		solver.inverse(inverse);
		REQUIRE((matrix * inverse).IsEqualTo(Matrix<Complex>::Identity(3), TOL(1e-9, 1e-4)));
	}

	TEST_CASE("Complex Cholesky enforces Hermitian positive definiteness",
				  "[HermitianCholesky][Cholesky][MML2]")
	{
		const Matrix<Complex> semidefinite{2, 2, {
			Complex{1.0, 0.0}, Complex{0.0, 1.0},
			Complex{0.0, -1.0}, Complex{1.0, 0.0}
		}};
		const Matrix<Complex> indefinite{2, 2, {
			Complex{1.0, 0.0}, Complex{}, Complex{}, Complex{-1.0, 0.0}
		}};
		const Matrix<Complex> nonHermitian{2, 2, {
			Complex{2.0, 0.0}, Complex{1.0, 1.0},
			Complex{1.0, 1.0}, Complex{2.0, 0.0}
		}};

		REQUIRE_THROWS_AS(CholeskySolver<Complex>(semidefinite), SingularMatrixError);
		REQUIRE_THROWS_AS(CholeskySolver<Complex>(indefinite), SingularMatrixError);
		REQUIRE_THROWS_AS(CholeskySolver<Complex>(nonHermitian), MatrixDimensionError);
	}

	TEST_CASE("MatrixAlg exposes complex Hermitian Cholesky and definiteness",
				  "[HermitianCholesky][MatrixAlg][Definiteness][MML2]")
	{
		const Matrix<Complex> positiveDefinite{2, 2, {
			Complex{4.0, 0.0}, Complex{1.0, 2.0},
			Complex{1.0, -2.0}, Complex{3.0, 0.0}
		}};
		const auto cholesky = MatrixAlg::CholeskyDecompose(positiveDefinite);
		REQUIRE((cholesky.L * Adjoint(cholesky.L)).IsEqualTo(positiveDefinite, TOL(1e-10, 1e-4)));
		REQUIRE(MatrixAlg::ClassifyDefiniteness(positiveDefinite) == MatrixAlg::Definiteness::PositiveDefinite);
		REQUIRE(MatrixAlg::IsPositiveDefinite(positiveDefinite));

		const Matrix<Complex> positiveSemidefinite{2, 2, {
			Complex{1.0, 0.0}, Complex{0.0, 1.0},
			Complex{0.0, -1.0}, Complex{1.0, 0.0}
		}};
		REQUIRE(MatrixAlg::ClassifyDefiniteness(positiveSemidefinite) == MatrixAlg::Definiteness::PositiveSemidefinite);
		REQUIRE(MatrixAlg::IsPositiveSemiDefinite(positiveSemidefinite));

		const Matrix<Complex> indefinite{2, 2, {
			Complex{2.0, 0.0}, Complex{}, Complex{}, Complex{-1.0, 0.0}
		}};
		REQUIRE(MatrixAlg::ClassifyDefiniteness(indefinite) == MatrixAlg::Definiteness::Indefinite);
		REQUIRE(MatrixAlg::IsIndefinite(indefinite));
	}

	TEST_CASE("Complex definiteness uses a Hermitian tolerance without analyzing the Hermitian part",
				  "[HermitianCholesky][MatrixAlg][Definiteness][MML2]")
	{
		const Real tolerance = TOL(1e-8, 1e-4);
		const Matrix<Complex> nearlyHermitian{2, 2, {
			Complex{3.0, tolerance * REAL(0.25)}, Complex{1.0, 2.0},
			Complex{1.0, REAL(-2.0) + tolerance * REAL(0.25)}, Complex{4.0, 0.0}
		}};
		REQUIRE(MatrixAlg::ClassifyDefiniteness(nearlyHermitian, tolerance) == MatrixAlg::Definiteness::PositiveDefinite);

		const Matrix<Complex> nonHermitian{2, 2, {
			Complex{3.0, 0.0}, Complex{1.0, 2.0},
			Complex{1.0, 2.0}, Complex{4.0, 0.0}
		}};
		REQUIRE_THROWS_AS(MatrixAlg::ClassifyDefiniteness(nonHermitian, tolerance), MatrixDimensionError);
		REQUIRE_THROWS_AS(MatrixAlg::IsPositiveDefinite(nonHermitian, tolerance), MatrixDimensionError);
	}
}
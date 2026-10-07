#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/algorithms/MatrixAlg.h>
#endif

#include <algorithm>
#include <concepts>

using namespace MML;
using namespace MML::Testing;

namespace MML::Tests::Algorithms::ComplexEigenanalysisTests
{
	Real EigenResidual(const Matrix<Complex>& matrix, const Vector<Complex>& vector, Complex eigenvalue)
	{
		return (matrix * vector - eigenvalue * vector).NormL2() /
			std::max(REAL(1.0), MatrixAlg::FrobeniusNorm(matrix) * vector.NormL2());
	}

	void RequireResiduals(const Matrix<Complex>& matrix,
			const MatrixAlg::EigensystemResult<Complex>& result, Real tolerance)
	{
		REQUIRE(result.eigenvalues.size() == matrix.rows());
		REQUIRE(result.eigenvectors.rows() == matrix.rows());
		REQUIRE(result.eigenvectors.cols() == matrix.cols());
		for (int index = 0; index < matrix.rows(); ++index)
			REQUIRE(EigenResidual(matrix, result.eigenvectors.VectorFromColumn(index),
				result.eigenvalues[index]) < tolerance);
	}

	TEST_CASE("Canonical general eigensystem preserves explicit complex types",
				  "[ComplexEigen][Contract][MML2]")
	{
		STATIC_REQUIRE(std::same_as<
			decltype(MatrixAlg::Eigensystem(std::declval<const Matrix<Real>&>())),
			MatrixAlg::EigensystemResult<Real>>);
		STATIC_REQUIRE(std::same_as<
			decltype(MatrixAlg::Eigensystem(std::declval<const Matrix<Complex>&>())),
			MatrixAlg::EigensystemResult<Complex>>);
		STATIC_REQUIRE(std::same_as<
			decltype(MatrixAlg::Eigenvalues(std::declval<const Matrix<Real>&>())),
			Vector<Complex>>);
		STATIC_REQUIRE(std::same_as<
			decltype(MatrixAlg::Eigenvalues(std::declval<const Matrix<Complex>&>())),
			Vector<Complex>>);
		STATIC_REQUIRE(std::same_as<
			decltype(MatrixAlg::SpectralRadius(std::declval<const Matrix<Complex>&>())), Real>);
		STATIC_REQUIRE(std::same_as<
			decltype(MatrixAlg::SymmetricEigenvalues(std::declval<const Matrix<Real>&>())), Vector<Real>>);
		STATIC_REQUIRE(std::same_as<
			decltype(MatrixAlg::HermitianEigenvalues(std::declval<const Matrix<Complex>&>())), Vector<Real>>);
	}

	TEST_CASE("Canonical real eigensystem expands rotation eigenpairs into complex columns",
				  "[ComplexEigen][RealInput][Rotation][MML2]")
	{
		const Real angle = REAL(0.7);
		const Matrix<Real> rotation{2, 2, {
			std::cos(angle), -std::sin(angle),
			std::sin(angle), std::cos(angle)
		}};
		const auto result = MatrixAlg::Eigensystem(rotation);

		REQUIRE(result.converged);
		REQUIRE(result.eigenvalues.size() == 2);
		REQUIRE(std::abs(std::abs(result.eigenvalues[0]) - REAL(1.0)) < TOL(1e-9, 1e-4));
		REQUIRE(std::abs(std::abs(result.eigenvalues[1]) - REAL(1.0)) < TOL(1e-9, 1e-4));
		REQUIRE(result.eigenvalues[0].imag() * result.eigenvalues[1].imag() < REAL(0.0));
		for (int index = 0; index < 2; ++index) {
			const Vector<Complex> vector = result.eigenvectors.VectorFromColumn(index);
			Vector<Complex> product(2);
			for (int row = 0; row < 2; ++row)
				for (int col = 0; col < 2; ++col)
					product[row] += rotation(row, col) * vector[col];
			REQUIRE((product - result.eigenvalues[index] * vector).NormL2() < TOL(1e-8, 1e-4));
		}
		REQUIRE(std::abs(MatrixAlg::SpectralRadius(rotation) - REAL(1.0)) < TOL(1e-9, 1e-4));
	}

	TEST_CASE("General complex eigensystem solves a nonnormal similarity transform",
				  "[ComplexEigen][General][Residual][MML2]")
	{
		const Matrix<Complex> similarity{2, 2, {
			Complex{1.0, 1.0}, Complex{2.0, -1.0},
			Complex{-0.5, 0.25}, Complex{1.0, 2.0}
		}};
		const Matrix<Complex> diagonal{2, 2, {
			Complex{1.0, 2.0}, Complex{}, Complex{}, Complex{-3.0, 0.5}
		}};
		const Matrix<Complex> matrix = similarity * diagonal * MatrixAlg::Inverse(similarity);
		const auto result = MatrixAlg::Eigensystem(matrix);

		REQUIRE(result.converged);
		REQUIRE(result.status == AlgorithmStatus::Success);
		REQUIRE(result.algorithmName == "ComplexQR");
		RequireResiduals(matrix, result, TOL(1e-9, 1e-4));

		std::vector<Complex> values{result.eigenvalues[0], result.eigenvalues[1]};
		std::sort(values.begin(), values.end(), [](Complex left, Complex right) { return left.real() < right.real(); });
		REQUIRE(std::abs(values[0] - Complex{-3.0, 0.5}) < TOL(1e-9, 1e-4));
		REQUIRE(std::abs(values[1] - Complex{1.0, 2.0}) < TOL(1e-9, 1e-4));
		REQUIRE(std::abs(MatrixAlg::SpectralRadius(matrix) - std::abs(Complex{-3.0, 0.5})) < TOL(1e-9, 1e-4));
	}

	TEST_CASE("General complex eigensystem handles triangular and defective matrices",
				  "[ComplexEigen][Triangular][Defective][MML2]")
	{
		const Matrix<Complex> triangular{3, 3, {
			Complex{1.0, 2.0}, Complex{4.0, -1.0}, Complex{-2.0, 0.5},
			Complex{}, Complex{-3.0, 0.5}, Complex{1.0, 3.0},
			Complex{}, Complex{}, Complex{2.0, -4.0}
		}};
		const auto triangularResult = MatrixAlg::Eigensystem(triangular);
		REQUIRE(triangularResult.converged);
		RequireResiduals(triangular, triangularResult, TOL(1e-9, 1e-4));

		const Matrix<Complex> jordan{3, 3, {
			Complex{2.0, 1.0}, Complex{1.0, 0.0}, Complex{},
			Complex{}, Complex{2.0, 1.0}, Complex{1.0, 0.0},
			Complex{}, Complex{}, Complex{2.0, 1.0}
		}};
		const auto jordanResult = MatrixAlg::Eigensystem(jordan);
		REQUIRE(jordanResult.converged);
		for (int index = 0; index < 3; ++index)
			REQUIRE(std::abs(jordanResult.eigenvalues[index] - Complex{2.0, 1.0}) < TOL(1e-9, 1e-4));
		REQUIRE(jordanResult.maxResidual < TOL(1e-8, 1e-3));
	}

	TEST_CASE("General complex QR performs multiple deflations on a dense matrix",
				  "[ComplexEigen][General][Deflation][MML2]")
	{
		const Matrix<Complex> similarity{3, 3, {
			Complex{1.0, 0.0}, Complex{1.0, -0.5}, Complex{0.0, 1.0},
			Complex{0.5, 1.0}, Complex{2.0, 0.0}, Complex{1.0, 0.25},
			Complex{1.0, -1.0}, Complex{0.0, 0.5}, Complex{1.5, 0.0}
		}};
		const Vector<Complex> expected{
			Complex{-2.0, 0.75}, Complex{0.5, -1.5}, Complex{3.0, 2.0}
		};
		Matrix<Complex> diagonal(3, 3);
		for (int index = 0; index < 3; ++index)
			diagonal(index, index) = expected[index];
		const Matrix<Complex> matrix = similarity * diagonal * MatrixAlg::Inverse(similarity);
		const auto result = MatrixAlg::Eigensystem(matrix);

		REQUIRE(result.converged);
		REQUIRE(result.iterations > 0);
		RequireResiduals(matrix, result, TOL(1e-8, 1e-3));
		for (int expectedIndex = 0; expectedIndex < expected.size(); ++expectedIndex) {
			Real closest = std::numeric_limits<Real>::infinity();
			for (int actualIndex = 0; actualIndex < result.eigenvalues.size(); ++actualIndex)
				closest = std::min(closest, static_cast<Real>(
					std::abs(expected[expectedIndex] - result.eigenvalues[actualIndex])));
			REQUIRE(closest < TOL(1e-8, 1e-3));
		}
	}

	TEST_CASE("Hermitian matrices use the Hermitian eigensolver path",
				  "[ComplexEigen][Hermitian][MML2]")
	{
		const Matrix<Complex> matrix{2, 2, {
			Complex{4.0, 0.0}, Complex{1.0, 2.0},
			Complex{1.0, -2.0}, Complex{3.0, 0.0}
		}};
		const auto result = MatrixAlg::Eigensystem(matrix);

		REQUIRE(result.converged);
		REQUIRE(result.algorithmName == "HermitianJacobi");
		REQUIRE(result.eigenvalues[0].imag() == REAL(0.0));
		REQUIRE(result.eigenvalues[1].imag() == REAL(0.0));
		REQUIRE(MatrixAlg::IsUnitary(result.eigenvectors));
		RequireResiduals(matrix, result, TOL(1e-9, 1e-4));
		REQUIRE(MatrixAlg::HermitianEigenvalues(matrix).IsEqualTo(
			Vector<Real>{REAL(1.2087121525220803), REAL(5.79128784747792)}, TOL(1e-9, 1e-4)));
	}

	TEST_CASE("General complex eigensystem reports invalid shape and nonconvergence",
				  "[ComplexEigen][Failure][MML2]")
	{
		REQUIRE_THROWS_AS(MatrixAlg::Eigensystem(Matrix<Complex>{2, 3}), MatrixDimensionError);
		REQUIRE_THROWS_AS(MatrixAlg::Eigensystem(Matrix<Complex>::Identity(2), REAL(-1.0)), DomainError);
		REQUIRE_THROWS_AS(MatrixAlg::Eigensystem(Matrix<Complex>::Identity(2), REAL(1e-8), -1), DomainError);

		const Matrix<Complex> matrix{3, 3, {
			Complex{1.0, 1.0}, Complex{2.0, -1.0}, Complex{3.0, 0.5},
			Complex{-2.0, 0.25}, Complex{0.5, -1.0}, Complex{1.0, 2.0},
			Complex{4.0, -0.5}, Complex{-1.0, 1.0}, Complex{2.0, 0.0}
		}};
		const auto result = MatrixAlg::Eigensystem(matrix, TOL(1e-12, 1e-5), 0);
		REQUIRE_FALSE(result.converged);
		REQUIRE(result.status == AlgorithmStatus::MaxIterationsExceeded);
		REQUIRE_FALSE(result.errorMessage.empty());
	}
}

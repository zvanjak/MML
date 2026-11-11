#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"

#if defined(MML_USE_SINGLE_HEADER) && !defined(MML_SINGLE_HEADER_INCLUDE)
#include <MML.h>
#elif !defined(MML_SINGLE_HEADER_INCLUDE)
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
		REQUIRE(std::isfinite(vector.NormL2()));
		REQUIRE(vector.NormL2() > Real{0});
		return (matrix * vector - eigenvalue * vector).NormL2() /
			std::max(static_cast<MML::Real>(1.0), MatrixAlg::FrobeniusNorm(matrix) * vector.NormL2());
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
		const Real angle = static_cast<MML::Real>(0.7);
		const Matrix<Real> rotation{2, 2, {
			std::cos(angle), -std::sin(angle),
			std::sin(angle), std::cos(angle)
		}};
		const auto result = MatrixAlg::Eigensystem(rotation);

		REQUIRE(result.converged);
		REQUIRE(result.eigenvalues.size() == 2);
		REQUIRE(std::abs(std::abs(result.eigenvalues[0]) - static_cast<MML::Real>(1.0)) < TOL(1e-9, 1e-4));
		REQUIRE(std::abs(std::abs(result.eigenvalues[1]) - static_cast<MML::Real>(1.0)) < TOL(1e-9, 1e-4));
		REQUIRE(result.eigenvalues[0].imag() * result.eigenvalues[1].imag() < static_cast<MML::Real>(0.0));
		for (int index = 0; index < 2; ++index) {
			const Vector<Complex> vector = result.eigenvectors.VectorFromColumn(index);
			Vector<Complex> product(2);
			for (int row = 0; row < 2; ++row)
				for (int col = 0; col < 2; ++col)
					product[row] += rotation(row, col) * vector[col];
			REQUIRE((product - result.eigenvalues[index] * vector).NormL2() < TOL(1e-8, 1e-4));
		}
		REQUIRE(std::abs(MatrixAlg::SpectralRadius(rotation) - static_cast<MML::Real>(1.0)) < TOL(1e-9, 1e-4));
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
		REQUIRE(result.eigenvalues[0].imag() == static_cast<MML::Real>(0.0));
		REQUIRE(result.eigenvalues[1].imag() == static_cast<MML::Real>(0.0));
		REQUIRE(MatrixAlg::IsUnitary(result.eigenvectors));
		RequireResiduals(matrix, result, TOL(1e-9, 1e-4));
		REQUIRE(MatrixAlg::HermitianEigenvalues(matrix).IsEqualTo(
			Vector<Real>{static_cast<MML::Real>(1.2087121525220803), static_cast<MML::Real>(5.79128784747792)}, TOL(1e-9, 1e-4)));
	}

	TEST_CASE("General complex eigensystem reports invalid shape and nonconvergence",
				  "[ComplexEigen][Failure][MML2]")
	{
		REQUIRE_THROWS_AS(MatrixAlg::Eigensystem(Matrix<Complex>{2, 3}), MatrixDimensionError);
		REQUIRE_THROWS_AS(MatrixAlg::Eigensystem(Matrix<Complex>::Identity(2), static_cast<MML::Real>(-1.0)), DomainError);
		REQUIRE_THROWS_AS(MatrixAlg::Eigensystem(Matrix<Complex>::Identity(2), static_cast<MML::Real>(1e-8), -1), DomainError);

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

	TEST_CASE("Canonical eigensystem reconstructs the tiny unbalanced conjugate pair",
				  "[ComplexEigen][StructuralPair][MML2]")
	{
		const Matrix<Real> matrix{2, 2, {0, Real(-1e-22), 1, 0}};
		const auto legacy = EigenSolver::Solve(matrix);
		const auto result = MatrixAlg::Eigensystem(matrix);
		REQUIRE(legacy.isComplexPair[0]);
		REQUIRE(legacy.isComplexPair[1]);
		REQUIRE(result.converged);
		REQUIRE(result.status == AlgorithmStatus::Success);
		REQUIRE(std::isfinite(result.maxResidual));
		for (int index = 0; index < 2; ++index) {
			REQUIRE(result.eigenvectors.VectorFromColumn(index).NormL2() > Real{0});
			for (int row = 0; row < 2; ++row) {
				const Complex expected{legacy.eigenvectors(row, 0), legacy.eigenvectors(row, 1)};
				CHECK(result.eigenvectors(row, index) == (index == 0 ? expected : std::conj(expected)));
			}
		}
		CHECK(IsWithinRel(result.eigenvalues[0].imag(), std::sqrt(-matrix(0, 1)),
			Precision::RoundoffZeroThreshold));
		CHECK(result.eigenvalues[0].imag() > Real{0});
		CHECK(result.eigenvalues[1] == std::conj(result.eigenvalues[0]));
	}

	TEST_CASE("Structural eigenpairs retain their columns across scales shifts and positions",
				  "[ComplexEigen][StructuralPair][MML2]")
	{
		const int exponent = GENERATE(0, 32, 36, 40);
		const Real scale = std::ldexp(Real{1}, GENERATE(-2, 2));
		const Real shift = Real(GENERATE(0, 2));
		const int pairStart = GENERATE(0, 1, 2);
		const Real omega = std::ldexp(Real{1}, -exponent);
		CAPTURE(exponent, scale, shift, pairStart);
		Matrix<Real> matrix(4, 4);
		for (int index = 0; index < 4; ++index)
			matrix(index, index) = scale * Real(5 + index);
		matrix(pairStart, pairStart) = matrix(pairStart + 1, pairStart + 1) = scale * shift;
		matrix(pairStart, pairStart + 1) = -scale * omega * omega;
		matrix(pairStart + 1, pairStart) = scale;
		const auto result = MatrixAlg::Eigensystem(matrix);
		REQUIRE(result.converged);
		REQUIRE(result.status == AlgorithmStatus::Success);
		CHECK(std::isfinite(result.maxResidual));
		CHECK(result.maxResidual <= Precision::RoundoffZeroThreshold);
		for (int index = 0; index < 4; ++index) {
			Real norm{};
			for (int row = 0; row < 4; ++row)
				norm = std::hypot(norm, std::abs(result.eigenvectors(row, index)));
			REQUIRE(std::isfinite(norm));
			REQUIRE(norm > Real{0});
		}
		CHECK(std::abs(result.eigenvalues[pairStart] - Complex{scale * shift, scale * omega}) <=
			Precision::RoundoffZeroThreshold * scale * omega);
		CHECK(result.eigenvalues[pairStart + 1] == std::conj(result.eigenvalues[pairStart]));
		const Complex lower = result.eigenvectors(pairStart + 1, pairStart);
		REQUIRE(std::abs(lower) > Real{0});
		CHECK(std::abs(result.eigenvectors(pairStart, pairStart) / lower - Complex{0, omega}) <=
			Precision::RoundoffZeroThreshold * omega);
		for (int row = 0; row < 4; ++row)
			CHECK(result.eigenvectors(row, pairStart + 1) == std::conj(result.eigenvectors(row, pairStart)));
	}

	TEST_CASE("Canonical residual validation rejects invalid eigenpairs and scaling",
				  "[ComplexEigen][ResidualValidation][MML2]")
	{
		Matrix<Real> matrix = Matrix<Real>::Identity(2);
		Vector<Complex> values{Complex{1}, Complex{1}};
		Matrix<Complex> vectors = Matrix<Complex>::Identity(2);
		SECTION("zero vector") { vectors(0, 0) = Complex{}; }
		SECTION("nonfinite vector") { vectors(0, 0) = Complex{std::numeric_limits<Real>::infinity()}; }
		SECTION("NaN vector") { vectors(0, 0) = Complex{0, std::numeric_limits<Real>::quiet_NaN()}; }
		SECTION("nonfinite eigenvalue") { values[0] = Complex{std::numeric_limits<Real>::infinity()}; }
		SECTION("NaN eigenvalue") { values[0] = Complex{std::numeric_limits<Real>::quiet_NaN()}; }
		SECTION("unrepresentable matrix norm") {
			for (int row = 0; row < 2; ++row)
				for (int col = 0; col < 2; ++col)
					matrix(row, col) = std::numeric_limits<Real>::max();
		}
		SECTION("unrepresentable denominator") {
			matrix(0, 0) = std::sqrt(std::numeric_limits<Real>::max());
			vectors(0, 0) = Complex{std::numeric_limits<Real>::max() / Real{2}};
		}
		SECTION("unrepresentable residual") {
			values[0] = Complex{std::numeric_limits<Real>::max()};
			vectors(1, 0) = Complex{1};
		}
		REQUIRE_THROWS_AS(MatrixAlg::Detail::MaximumEigenResidual(matrix, values, vectors), MatrixNumericalError);
	}

	TEST_CASE("Residual validation accepts nonzero subnormal vectors and the zero matrix",
				  "[ComplexEigen][ResidualValidation][MML2]")
	{
		const Real tiny = std::numeric_limits<Real>::denorm_min();
		const Vector<Complex> values{Complex{1}, Complex{1}};
		const Matrix<Complex> vectors{2, 2, {Complex{tiny}, Complex{}, Complex{}, Complex{tiny}}};
		CHECK(MatrixAlg::Detail::MaximumEigenResidual(Matrix<Real>::Identity(2), values, vectors) == Real{0});
		const auto zero = MatrixAlg::Eigensystem(Matrix<Real>(2, 2));
		REQUIRE(zero.converged);
		CHECK(zero.status == AlgorithmStatus::Success);
		CHECK(zero.maxResidual == Real{0});
		for (int index = 0; index < 2; ++index)
			CHECK(zero.eigenvectors.VectorFromColumn(index).NormL2() > Real{0});
	}

	TEST_CASE("Canonical eigensystems report unrepresentable residuals as failure data",
				  "[ComplexEigen][ResidualValidation][MML2]")
	{
		Matrix<Real> real = Matrix<Real>::Identity(9) * (std::numeric_limits<Real>::max() / Real{2});
		Matrix<Complex> complex(9, 9);
		for (int row = 0; row < 9; ++row)
			complex(row, row) = Complex{real(row, row)};
		auto requireFailure = [](const auto& result) {
			CHECK_FALSE(result.converged);
			CHECK(result.status == AlgorithmStatus::NumericalInstability);
			CHECK_FALSE(HasEstimate(result.maxResidual));
			CHECK_FALSE(result.errorMessage.empty());
		};
		requireFailure(MatrixAlg::SymmetricEigensystem(real));
		requireFailure(MatrixAlg::HermitianEigensystem(complex));
		requireFailure(MatrixAlg::Eigensystem(real));
		requireFailure(MatrixAlg::Eigensystem(complex));
		CHECK_THROWS_AS(MatrixAlg::SymmetricEigenvalues(real), MatrixNumericalError);
		CHECK_THROWS_AS(MatrixAlg::HermitianEigenvalues(complex), MatrixNumericalError);
		real(0, 1) = real(0, 0);
		complex(0, 1) = complex(0, 0);
		requireFailure(MatrixAlg::Eigensystem(real));
		requireFailure(MatrixAlg::Eigensystem(complex));
		CHECK_THROWS_AS(MatrixAlg::Eigenvalues(real), MatrixNumericalError);
		CHECK_THROWS_AS(MatrixAlg::SpectralRadius(complex), MatrixNumericalError);
	}

	TEST_CASE("Residual validation preserves iteration-limit results on all canonical paths",
				  "[ComplexEigen][ResidualValidation][MML2]")
	{
		const Real tolerance = Precision::EigenSolverConvergenceTolerance;
		const Matrix<Real> general{3, 3, {1, 2, 3, -2, 1, 1, 4, -1, 2}};
		const Matrix<Real> symmetric{3, 3, {1, 2, 3, 2, 1, 1, 3, 1, 2}};
		const Matrix<Complex> complex{3, 3, {
			Complex{1, 1}, Complex{2, -1}, Complex{3},
			Complex{-2}, Complex{1}, Complex{1, 2},
			Complex{4}, Complex{-1}, Complex{2}
		}};
		const Matrix<Complex> hermitian{3, 3, {
			Complex{1}, Complex{2, 1}, Complex{3},
			Complex{2, -1}, Complex{1}, Complex{1},
			Complex{3}, Complex{1}, Complex{2}
		}};
		auto requireLimit = [](const auto& result) {
			CHECK_FALSE(result.converged);
			CHECK(result.status == AlgorithmStatus::MaxIterationsExceeded);
			CHECK_FALSE(result.errorMessage.empty());
		};
		requireLimit(MatrixAlg::Eigensystem(general, tolerance, 0));
		requireLimit(MatrixAlg::Eigensystem(complex, tolerance, 0));
		requireLimit(MatrixAlg::SymmetricEigensystem(symmetric, tolerance, 0));
		requireLimit(MatrixAlg::HermitianEigensystem(hermitian, tolerance, 0));
	}
}

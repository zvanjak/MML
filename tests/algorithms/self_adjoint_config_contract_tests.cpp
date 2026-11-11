#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"

#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/algorithms/EigenSystemSolvers.h>
#include <mml/algorithms/MatrixAlg.h>
#endif

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>

using namespace MML;
using Catch::Matchers::WithinAbs;

namespace {
	template<class Solver>
	void CheckRealDomain() {
		const Matrix<Real> wide(2, 3, {1, 0, 7, 0, 2, 8});
		const Matrix<Real> nonsymmetric(2, 2, {1, 2, 0, 3});
		CHECK_THROWS_AS(Solver::Solve(wide), MatrixDimensionError);
		CHECK_THROWS_AS(Solver::Solve(wide, EigenSolverConfig{}), MatrixDimensionError);
		CHECK_THROWS_AS(Solver::Solve(nonsymmetric), MatrixDimensionError);
		CHECK_THROWS_AS(Solver::Solve(nonsymmetric, EigenSolverConfig{}), MatrixDimensionError);
		const Matrix<Real> narrow(3, 2, {1, 0, 0, 2, 7, 8});
		CHECK_THROWS_AS(Solver::Solve(narrow), MatrixDimensionError);
		CHECK_THROWS_AS(Solver::Solve(narrow, EigenSolverConfig{}), MatrixDimensionError);
		CHECK_THROWS_AS(Solver::SolveSymmetricPart(narrow), MatrixDimensionError);
		CHECK_THROWS_AS(Solver::SolveSymmetricPart(wide, EigenSolverConfig{}), MatrixDimensionError);
	}

	template<class Result, class Scalar>
	void CheckEigenpairs(const Result& result, const Matrix<Scalar>& matrix,
		Real residualBudget = Precision::DefaultToleranceStrict) {
		REQUIRE(result.converged);
		CHECK(result.status == AlgorithmStatus::Success);
		REQUIRE(result.eigenvectors.rows() == matrix.rows());
		REQUIRE(result.eigenvectors.cols() == matrix.rows());
		for (int col = 0; col < matrix.rows(); ++col) {
			for (int row = 0; row < matrix.rows(); ++row) {
				Scalar product{};
				for (int k = 0; k < matrix.rows(); ++k)
					product += matrix(row, k) * result.eigenvectors(k, col);
				CHECK(std::abs(product - result.eigenvalues[col] * result.eigenvectors(row, col)) <=
					residualBudget);
			}
			for (int other = 0; other < matrix.rows(); ++other) {
				Scalar dot{};
				for (int row = 0; row < matrix.rows(); ++row) {
					if constexpr (MMLComplex<Scalar>)
						dot += std::conj(result.eigenvectors(row, col)) * result.eigenvectors(row, other);
					else
						dot += result.eigenvectors(row, col) * result.eigenvectors(row, other);
				}
				CHECK(std::abs(dot - Scalar(col == other ? Real(1) : Real(0))) <= Precision::DefaultToleranceStrict);
			}
		}
	}

	template<class Solver>
	void CheckRealSwitches() {
		const Matrix<Real> diagonal(3, 3, {3, 0, 0, 0, 1, 0, 0, 0, 2});
		const Matrix<Real> dense(3, 3, {2, 1, 0, 1, 2, 1, 0, 1, 2});
		for (bool sorted : {false, true})
			for (bool vectors : {false, true}) {
				EigenSolverConfig config;
				config.sort_eigenvalues = sorted;
				config.compute_eigenvectors = vectors;
				const auto result = Solver::Solve(diagonal, config);
				REQUIRE(result.converged);
				CHECK(result.eigenvalues == (sorted ? Vector<Real>{1, 2, 3} : Vector<Real>{3, 1, 2}));
				if (vectors) CheckEigenpairs(result, diagonal);
				else {
					CHECK(result.eigenvectors.rows() == 0);
					CHECK(result.eigenvectors.cols() == 0);
				}
				const auto computed = Solver::Solve(dense, config);
				REQUIRE(computed.converged);
				auto values = computed.eigenvalues;
				std::sort(values.begin(), values.end());
				for (int i = 0; i < 3; ++i)
					CHECK_THAT(values[i], WithinAbs(2 + (i - 1) * std::sqrt(Real(2)), Precision::DefaultToleranceStrict));
				if (vectors) CheckEigenpairs(computed, dense);
				else CHECK(computed.eigenvectors.rows() == 0);
				const auto packed = Solver::Solve(MatrixSym<Real>(3, {2, 1, 2, 0, 1, 2}), config);
				REQUIRE(packed.converged);
				CHECK(packed.eigenvalues == computed.eigenvalues);
				CHECK(packed.eigenvectors.rows() == computed.eigenvectors.rows());
			}
	}

	template<class Solver>
	void CheckRealControls() {
		const Matrix<Real> diagonal(2, 2, {3, 0, 0, 1});
		const Matrix<Real> coupled(2, 2, {2, 1, 1, 2});
		const Real tolerance = PrecisionValues<Real>::EigenSolverConvergenceTolerance;
		const auto complete = Solver::Solve(diagonal, tolerance, 0);
		CHECK(complete.converged);
		CHECK(complete.iterations == 0);
		const auto stopped = Solver::Solve(coupled, tolerance, 0);
		CHECK_FALSE(stopped.converged);
		CHECK(stopped.iterations == 0);
		CHECK(stopped.status == AlgorithmStatus::MaxIterationsExceeded);
		CHECK_FALSE(stopped.error_message.empty());
		EigenSolverConfig config;
		config.max_iterations = 0;
		CHECK(Solver::Solve(coupled, config).converged);
		config.max_iterations = -1;
		CHECK_THROWS_AS(Solver::Solve(diagonal, config), DomainError);
		CHECK_THROWS_AS(Solver::Solve(diagonal, tolerance, -1), DomainError);
		for (Real invalid : {Real(-1), std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity()}) {
			CHECK_THROWS_AS(Solver::Solve(diagonal, invalid), DomainError);
			config = EigenSolverConfig{};
			config.tolerance = invalid;
			CHECK_THROWS_AS(Solver::Solve(diagonal, config), DomainError);
		}
	}
}

TEST_CASE("Symmetric eigensolvers validate shape and domain in every overload", "[EigenSolvers][self-adjoint-config]")
{
	CheckRealDomain<SymmMatEigenSolverJacobi>();
	CheckRealDomain<SymmMatEigenSolverQR>();
	const Matrix<Real> wide(2, 3, {1, 0, 7, 0, 2, 8});
	CHECK_THROWS_AS(SymmMatEigenSolverJacobi::SolveSymmetricPart(wide), MatrixDimensionError);
	CHECK_THROWS_AS(ReduceToHessenberg(wide), MatrixDimensionError);
	CHECK_THROWS_AS(MatrixAlg::ReduceToHessenberg(wide), MatrixDimensionError);
}

TEST_CASE("Jacobi honors eigenvalue ordering and eigenvector configuration", "[EigenSolvers][self-adjoint-config]")
{
	CheckRealSwitches<SymmMatEigenSolverJacobi>();
}

TEST_CASE("Symmetric QR honors eigenvalue ordering and eigenvector configuration", "[EigenSolvers][self-adjoint-config]")
{
	CheckRealSwitches<SymmMatEigenSolverQR>();
}

TEST_CASE("Symmetric iteration budgets and invalid controls are explicit", "[EigenSolvers][self-adjoint-config]")
{
	CheckRealControls<SymmMatEigenSolverJacobi>();
	CheckRealControls<SymmMatEigenSolverQR>();
}

TEST_CASE("Hermitian Jacobi honors sorting eigenvector and budget controls", "[EigenSolvers][self-adjoint-config]")
{
	const Matrix<Complex> diagonal(2, 2, {3, 0, 0, 1});
	const Matrix<Complex> coupled(2, 2, {2, Complex(0, 1), Complex(0, -1), 2});
	for (bool sorted : {false, true})
		for (bool vectors : {false, true}) {
			EigenSolverConfig config;
			config.sort_eigenvalues = sorted;
			config.compute_eigenvectors = vectors;
			const auto result = HermitianMatEigenSolverJacobi::Solve(diagonal, config);
			REQUIRE(result.converged);
			CHECK(result.eigenvalues == (sorted ? Vector<Real>{1, 3} : Vector<Real>{3, 1}));
			if (vectors) CheckEigenpairs(result, diagonal);
			else CHECK(result.eigenvectors.rows() == 0);
			const auto dense = HermitianMatEigenSolverJacobi::Solve(coupled, config);
			REQUIRE(dense.converged);
			auto values = dense.eigenvalues;
			std::sort(values.begin(), values.end());
			CHECK_THAT(values[0], WithinAbs(1, Precision::DefaultToleranceStrict));
			CHECK_THAT(values[1], WithinAbs(3, Precision::DefaultToleranceStrict));
			if (vectors) CheckEigenpairs(dense, coupled);
			else CHECK(dense.eigenvectors.rows() == 0);
		}
	const auto stopped = HermitianMatEigenSolverJacobi::Solve(coupled, Precision::DefaultToleranceStrict, 0);
	CHECK_FALSE(stopped.converged);
	CHECK(stopped.status == AlgorithmStatus::MaxIterationsExceeded);
	CHECK_THROWS_AS(HermitianMatEigenSolverJacobi::Solve(diagonal, Real(-1)), DomainError);
}

TEST_CASE("Empty and nonfinite self-adjoint results never claim successful convergence", "[EigenSolvers][self-adjoint-config]")
{
	const auto jacobi = SymmMatEigenSolverJacobi::Solve(MatrixSym<Real>{});
	CHECK_FALSE(jacobi.converged);
	CHECK(jacobi.status == AlgorithmStatus::Completed);
	CHECK(std::isnan(jacobi.residual));
	const auto hermitian = HermitianMatEigenSolverJacobi::Solve(Matrix<Complex>{});
	CHECK_FALSE(hermitian.converged);
	CHECK(hermitian.status == AlgorithmStatus::Completed);
	CHECK(SymmMatEigenSolverQR::Solve(Matrix<Real>{}).status == AlgorithmStatus::Completed);
	CHECK(ReduceToHessenberg(Matrix<Real>{}).H.rows() == 0);
	for (bool vectors : {false, true}) {
		EigenSolverConfig config;
		config.compute_eigenvectors = vectors;
		MatrixSym<Real> invalid(2);
		invalid(0, 0) = std::numeric_limits<Real>::quiet_NaN();
		const auto j = SymmMatEigenSolverJacobi::Solve(invalid, config);
		const auto q = SymmMatEigenSolverQR::Solve(invalid, config);
		CHECK(j.status == AlgorithmStatus::NumericalInstability);
		CHECK(q.status == AlgorithmStatus::NumericalInstability);
		CHECK_FALSE(j.converged);
		CHECK_FALSE(q.converged);
		if (!vectors) {
			CHECK(j.eigenvectors.rows() == 0);
			CHECK(q.eigenvectors.rows() == 0);
		}
	}
}

TEST_CASE("Explicit symmetric parts and focused matrix consumers retain their domains", "[EigenSolvers][self-adjoint-config]")
{
	const Matrix<Real> nonsymmetric(2, 2, {1, 2, 0, 3});
	EigenSolverConfig config;
	config.compute_eigenvectors = false;
	const auto j = SymmMatEigenSolverJacobi::SolveSymmetricPart(nonsymmetric, config);
	const auto q = SymmMatEigenSolverQR::SolveSymmetricPart(nonsymmetric, config);
	REQUIRE(j.converged);
	REQUIRE(q.converged);
	CHECK_THAT(j.eigenvalues[0], WithinAbs(2 - std::sqrt(Real(2)), Precision::DefaultToleranceStrict));
	CHECK_THAT(q.eigenvalues[1], WithinAbs(2 + std::sqrt(Real(2)), Precision::DefaultToleranceStrict));
	CHECK(j.eigenvectors.rows() == 0);
	CHECK(q.eigenvectors.rows() == 0);
	CHECK_THROWS_AS(MatrixAlg::SymmetricEigensystem(nonsymmetric), MatrixDimensionError);
	CHECK_THROWS_AS(MatrixAlg::ClassifyDefiniteness(nonsymmetric), MatrixDimensionError);
	const Matrix<Real> symmetric(2, 2, {2, 1, 1, 2});
	const auto result = MatrixAlg::SymmetricEigensystem(symmetric);
	CheckEigenpairs(result, symmetric);
	CHECK(MatrixAlg::IsPositiveDefinite(symmetric));
}

TEST_CASE("Self-adjoint switches preserve spectra and residuals across scales", "[EigenSolvers][self-adjoint-config]")
{
	const Matrix<Real> unit(4, 4, {Real(2.25), Real(0.25), Real(-0.75), Real(5.25),
		Real(0.25), Real(2.25), Real(5.25), Real(-0.75),
		Real(-0.75), Real(5.25), Real(2.25), Real(0.25),
		Real(5.25), Real(-0.75), Real(0.25), Real(2.25)});
	// Hadamard similarity with prescribed eigenvalues {-4, -2, 7, 8}.
	for (int exponent : {std::max(-600, std::numeric_limits<Real>::min_exponent / 2), 0,
		std::min(600, std::numeric_limits<Real>::max_exponent / 2)}) {
		const Real scale = std::ldexp(Real(1), exponent);
		for (bool sorted : {false, true})
			for (bool vectors : {false, true}) {
				EigenSolverConfig config;
				config.tolerance = PrecisionValues<Real>::EigenSolverConvergenceTolerance * scale;
				config.sort_eigenvalues = sorted;
				config.compute_eigenvectors = vectors;
				CAPTURE(exponent, sorted, vectors);
				const auto j = SymmMatEigenSolverJacobi::Solve(unit * scale, config);
				// QR also uses tolerance as a relative deflation factor; keep that factor small.
				auto qrConfig = config;
				qrConfig.tolerance = std::min(config.tolerance, PrecisionValues<Real>::EigenSolverConvergenceTolerance);
				const auto q = SymmMatEigenSolverQR::Solve(unit * scale, qrConfig);
				REQUIRE(j.converged);
				REQUIRE(q.converged);
				auto jvalues = j.eigenvalues, qvalues = q.eigenvalues;
				std::sort(jvalues.begin(), jvalues.end());
				std::sort(qvalues.begin(), qvalues.end());
				const Vector<Real> expected{-4, -2, 7, 8};
				// Similarity updates and length-four residual products accumulate
				// roundoff proportional to the spectral scale (8 here).
				const Real roundoff = 4 * unit.rows() * Precision::Eps;
				for (int i = 0; i < 4; ++i) {
					const Real accuracy = std::max(Precision::DefaultToleranceStrict, roundoff * std::abs(expected[i]));
					CHECK_THAT(jvalues[i] / scale, WithinAbs(expected[i], accuracy));
					CHECK_THAT(qvalues[i] / scale, WithinAbs(expected[i], accuracy));
				}
				if (vectors) {
					auto normalizedJ = j;
					auto normalizedQ = q;
					for (auto& value : normalizedJ.eigenvalues) value /= scale;
					for (auto& value : normalizedQ.eigenvalues) value /= scale;
					const Real residualBudget = std::max(Precision::DefaultToleranceStrict, roundoff * 8);
					CheckEigenpairs(normalizedJ, unit, residualBudget);
					CheckEigenpairs(normalizedQ, unit, residualBudget);
				} else {
					CHECK(j.eigenvectors.rows() == 0);
					CHECK(q.eigenvectors.rows() == 0);
				}
			}
	}
}

TEST_CASE("Self-adjoint verbose diagnostics are opt-in and iteration budgets are bounded", "[EigenSolvers][self-adjoint-config]")
{
	std::ostringstream output;
	struct RestoreBuffer {
		std::streambuf* original;
		~RestoreBuffer() { std::clog.rdbuf(original); }
	} restore{std::clog.rdbuf(output.rdbuf())};
	EigenSolverConfig config;
	const MatrixSym<Real> packed(2, {2, 1, 2});
	(void)SymmMatEigenSolverJacobi::Solve(packed, config);
	CHECK(output.str().empty());
	config.verbose = true;
	(void)SymmMatEigenSolverJacobi::Solve(packed, config);
	(void)SymmMatEigenSolverQR::Solve(packed, config);
	(void)HermitianMatEigenSolverJacobi::Solve(Matrix<Complex>(2, 2, {2, 1, 1, 2}), config);
	CHECK(output.str().find("Jacobi:") != std::string::npos);
	CHECK(output.str().find("QR:") != std::string::npos);
	CHECK(output.str().find("HermitianJacobi:") != std::string::npos);
	CHECK(output.str().find("iterations=") != std::string::npos);
	const MatrixSym<Real> dense(4, {4, 1, 3, 2, 1, 5, 1, 2, 1, 6});
	for (int budget : {0, 1, 2, 10, 11, 12}) {
		const auto result = SymmMatEigenSolverQR::Solve(dense, Precision::RoundoffZeroThreshold, budget);
		CHECK(result.iterations <= budget);
		CHECK(result.converged == (result.status == AlgorithmStatus::Success));
		if (!result.converged) CHECK_FALSE(result.error_message.empty());
	}
}

TEST_CASE("Self-adjoint shortcut and Hermitian invalid-input contracts honor configuration", "[EigenSolvers][self-adjoint-config]")
{
	for (bool vectors : {false, true})
		for (int n : {0, 1, 2}) {
			EigenSolverConfig config;
			config.compute_eigenvectors = vectors;
			MatrixSym<Real> real(n);
			Matrix<Complex> complex;
			if (n > 0) complex = Matrix<Complex>(n, n);
			for (int i = 0; i < n; ++i) {
				real(i, i) = 3;
				complex(i, i) = 3;
			}
			auto check = [&](const auto& result, bool finite) {
				CHECK(result.eigenvectors.rows() == (vectors ? n : 0));
				CHECK(result.eigenvectors.cols() == (vectors ? n : 0));
				CHECK(result.iterations == 0);
				if (n == 0) {
					CHECK_FALSE(result.converged);
					CHECK(result.status == AlgorithmStatus::Completed);
					CHECK(std::isnan(result.residual));
				} else if (finite) {
					REQUIRE(result.converged);
					CHECK(result.status == AlgorithmStatus::Success);
					CHECK(result.eigenvalues[0] == 3);
				} else {
					CHECK_FALSE(result.converged);
					CHECK(result.status == AlgorithmStatus::NumericalInstability);
					CHECK_FALSE(result.error_message.empty());
				}
			};
			check(SymmMatEigenSolverJacobi::Solve(real, config), true);
			check(SymmMatEigenSolverQR::Solve(real, config), true);
			check(HermitianMatEigenSolverJacobi::Solve(complex, config), true);
			if (n > 0) {
				real(0, 0) = std::numeric_limits<Real>::quiet_NaN();
				complex(0, 0) = Complex(std::numeric_limits<Real>::quiet_NaN(), 0);
				check(SymmMatEigenSolverJacobi::Solve(real, config), false);
				check(SymmMatEigenSolverQR::Solve(real, config), false);
				check(HermitianMatEigenSolverJacobi::Solve(complex, config), false);
			}
		}
	for (Real invalid : {Real(-1), std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity()}) {
		EigenSolverConfig config;
		config.tolerance = invalid;
		CHECK_THROWS_AS(HermitianMatEigenSolverJacobi::Solve(Matrix<Complex>{}, invalid), DomainError);
		CHECK_THROWS_AS(HermitianMatEigenSolverJacobi::Solve(Matrix<Complex>{}, config), DomainError);
		CHECK_THROWS_AS(SymmMatEigenSolverJacobi::Solve(MatrixSym<Real>{}, config), DomainError);
		CHECK_THROWS_AS(SymmMatEigenSolverQR::Solve(MatrixSym<Real>{}, config), DomainError);
	}
	EigenSolverConfig negative;
	negative.max_iterations = -1;
	CHECK_THROWS_AS(HermitianMatEigenSolverJacobi::Solve(Matrix<Complex>{}, negative), DomainError);
	CHECK_THROWS_AS(HermitianMatEigenSolverJacobi::Solve(Matrix<Complex>{}, Precision::DefaultToleranceStrict, -1), DomainError);
	for (const auto& rectangular : {Matrix<Complex>(2, 3), Matrix<Complex>(3, 2)}) {
		CHECK_THROWS_AS(HermitianMatEigenSolverJacobi::Solve(rectangular), MatrixDimensionError);
		CHECK_THROWS_AS(HermitianMatEigenSolverJacobi::Solve(rectangular, EigenSolverConfig{}), MatrixDimensionError);
	}
	const Matrix<Complex> coupled(2, 2, {2, Complex(0, 1), Complex(0, -1), 2});
	CHECK(HermitianMatEigenSolverJacobi::Solve(coupled, EigenSolverConfig{}).converged);
}

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "../TestPrecision.h"

#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/systems/LinearSystem.h>
#endif

#include <cmath>
#include <limits>
#include <type_traits>
#include <utility>

using namespace MML;
using namespace MML::Systems;
using Catch::Matchers::WithinAbs;

namespace {
#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable: 4996)
#elif defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif
	const Vector<Real>& LegacyRHS(const LinearSystem<Real>& system) {
		return system.GetRHS();
	}
#if defined(_MSC_VER)
#pragma warning(pop)
#elif defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

	void CheckVectorOperationsReject(const LinearSystem<Real>& system) {
		CHECK_THROWS_AS(system.Solve(), InvalidStateError);
		CHECK_THROWS_AS(system.SolveByGaussJordan(), InvalidStateError);
		CHECK_THROWS_AS(system.SolveByLU(), InvalidStateError);
		CHECK_THROWS_AS(system.SolveByCholesky(), InvalidStateError);
		CHECK_THROWS_AS(system.SolveByQR(), InvalidStateError);
		CHECK_THROWS_AS(system.SolveBySVD(), InvalidStateError);
		CHECK_THROWS_AS(system.SolveLeastSquares(), InvalidStateError);
		CHECK_THROWS_AS(system.SolveIterative(), InvalidStateError);
		CHECK_THROWS_AS(system.ResidualNorm(Vector<Real>(system.Cols())), InvalidStateError);
		CHECK_THROWS_AS(system.RelativeResidual(Vector<Real>(system.Cols())), InvalidStateError);
		CHECK_THROWS_AS(system.Verify(Vector<Real>(system.Cols())), InvalidStateError);
		CHECK_THROWS_AS(system.getRHS(), InvalidStateError);
		CHECK_THROWS_AS(LegacyRHS(system), InvalidStateError);
	}

	void CheckMultipleSolution(const Matrix<Real>& coefficient, const Matrix<Real>& rhs,
		const Matrix<Real>& expected, LinearSolverRecommendation recommendation, Real scale = Real(1)) {
		const LinearSystem<Real> system(coefficient, rhs);
		CHECK(system.Analyze().recommendedSolver == recommendation);
		const auto actual = system.SolveMultiple();
		const auto repeated = system.SolveMultiple();
		REQUIRE(actual.rows() == coefficient.cols());
		REQUIRE(actual.cols() == rhs.cols());
		REQUIRE(repeated.rows() == actual.rows());
		REQUIRE(repeated.cols() == actual.cols());
		for (int col = 0; col < rhs.cols(); ++col) {
			const auto single = LinearSystem<Real>(coefficient, rhs.VectorFromColumn(col)).Solve();
			for (int row = 0; row < expected.rows(); ++row) {
				CHECK_THAT(actual(row, col), WithinAbs(expected(row, col), Precision::EqualityTolerance));
				CHECK_THAT(repeated(row, col), WithinAbs(expected(row, col), Precision::EqualityTolerance));
				CHECK_THAT(actual(row, col), WithinAbs(single[row], Precision::EqualityTolerance));
			}
			for (int index = 0; index < coefficient.cols(); ++index) {
				Real normalResidual = 0;
				Real normalMagnitude = 0;
				for (int row = 0; row < coefficient.rows(); ++row) {
					Real residual = -rhs(row, col) / scale;
					Real residualMagnitude = std::abs(rhs(row, col) / scale);
					for (int inner = 0; inner < coefficient.cols(); ++inner) {
						const Real term = (coefficient(row, inner) / scale) * actual(inner, col);
						residual += term;
						residualMagnitude += std::abs(term);
					}
					normalResidual += (coefficient(row, index) / scale) * residual;
					normalMagnitude += std::abs(coefficient(row, index) / scale) * residualMagnitude;
				}
				Real tolerance = Precision::EqualityTolerance;
				if constexpr (std::is_same_v<Real, float>) {
					const Real operations = Real(coefficient.cols() + coefficient.rows() + 2);
					tolerance = std::max(tolerance, operations * Precision::Eps * normalMagnitude);
				}
				CHECK_THAT(normalResidual, WithinAbs(Real(0), tolerance));
			}
		}
	}
}

TEST_CASE("Matrix RHS cannot be read as an empty vector RHS", "[LinearSystem][rhs-contract]") {
	const LinearSystem<Real> system(Matrix<Real>::Identity(2), Matrix<Real>::Identity(2));
	CHECK_THROWS_AS(system.getRHS(), InvalidStateError);
	CHECK_THROWS_AS(LegacyRHS(system), InvalidStateError);
}

TEST_CASE("RHS mode guards cover every vector entry point before solving", "[LinearSystem][rhs-contract]") {
	for (const Matrix<Real>& coefficient : {
		Matrix<Real>(2, 2, {2, 1, 0, 3}),
		Matrix<Real>(2, 2, {2, 0, 1, 3}),
		Matrix<Real>(2, 2, {4, 1, 1, 3}),
		Matrix<Real>(3, 2, {1, 0, 0, 1, 1, 1})}) {
		INFO("rows " << coefficient.rows());
		CheckVectorOperationsReject(LinearSystem<Real>(coefficient, Matrix<Real>(coefficient.rows(), 1)));
		CheckVectorOperationsReject(LinearSystem<Real>(coefficient));
		CHECK_THROWS_AS(LinearSystem<Real>(coefficient).getRHSMatrix(), InvalidStateError);
		CHECK_THROWS_AS(LinearSystem<Real>(coefficient).SolveMultiple(), InvalidStateError);
	}
}

TEST_CASE("RHS accessors retain owned const data and reject the other mode", "[LinearSystem][rhs-contract]") {
	STATIC_REQUIRE(std::is_same_v<decltype(std::declval<const LinearSystem<Real>&>().getRHS()),
		const Vector<Real>&>);
	STATIC_REQUIRE(std::is_same_v<decltype(std::declval<const LinearSystem<Real>&>().getRHSMatrix()),
		const Matrix<Real>&>);
	const auto vectorSystem = [] {
		Vector<Real> rhs{1, 2};
		LinearSystem<Real> system(Matrix<Real>::Identity(2), rhs);
		rhs[0] = 99;
		return system;
	}();
	const auto matrixSystem = [] {
		Matrix<Real> rhs(2, 1, {3, 4});
		LinearSystem<Real> system(Matrix<Real>::Identity(2), rhs);
		rhs(0, 0) = 99;
		return system;
	}();
	CHECK(vectorSystem.getRHS()[0] == Real(1));
	CHECK(vectorSystem.getRHS()[1] == Real(2));
	CHECK(&LegacyRHS(vectorSystem) == &vectorSystem.getRHS());
	CHECK_THROWS_AS(vectorSystem.getRHSMatrix(), InvalidStateError);
	CHECK_THROWS_AS(vectorSystem.SolveMultiple(), InvalidStateError);
	REQUIRE(matrixSystem.getRHSMatrix().rows() == 2);
	REQUIRE(matrixSystem.getRHSMatrix().cols() == 1);
	CHECK(matrixSystem.getRHSMatrix()(0, 0) == Real(3));
	CHECK(matrixSystem.getRHSMatrix()(1, 0) == Real(4));
	CHECK(matrixSystem.SolveMultiple()(1, 0) == Real(4));
	CHECK_THROWS_AS(matrixSystem.getRHS(), InvalidStateError);
	CHECK_THROWS_AS(LegacyRHS(matrixSystem), InvalidStateError);
}

TEST_CASE("Multiple RHS preserves both triangular orientations across units", "[LinearSystem][rhs-contract]") {
	for (const int exponent : {-60, 0, 60}) {
		const Real scale = std::ldexp(Real(1), exponent);
		const Matrix<Real> expected(2, 3, {1, 2, 0, 2, -1, 0});
		CheckMultipleSolution(Matrix<Real>(2, 2, {2 * scale, scale, 0, 3 * scale}),
			Matrix<Real>(2, 3, {4 * scale, 3 * scale, 0, 6 * scale, -3 * scale, 0}),
			expected, LinearSolverRecommendation::Triangular, scale);
		CheckMultipleSolution(Matrix<Real>(2, 2, {2 * scale, 0, scale, 3 * scale}),
			Matrix<Real>(2, 3, {2 * scale, 4 * scale, 0, 7 * scale, -scale, 0}),
			expected, LinearSolverRecommendation::Triangular, scale);
	}
}

TEST_CASE("Multiple RHS preserves SPD and general square solves", "[LinearSystem][rhs-contract]") {
	const Matrix<Real> expected(2, 3, {1, 2, 0, 2, -1, 0});
	CheckMultipleSolution(Matrix<Real>(2, 2, {4, 1, 1, 3}),
		Matrix<Real>(2, 3, {6, 7, 0, 7, -1, 0}), expected, LinearSolverRecommendation::Cholesky);
	CheckMultipleSolution(Matrix<Real>(2, 2, {4, 1, 2, 3}),
		Matrix<Real>(2, 3, {6, 7, 0, 8, 1, 0}), expected, LinearSolverRecommendation::LU);
}

TEST_CASE("Multiple RHS solves tall exact and inconsistent least squares", "[LinearSystem][rhs-contract]") {
	for (const int exponent : {-60, 0, 60}) {
		const Real scale = std::ldexp(Real(1), exponent);
		CheckMultipleSolution(Matrix<Real>(3, 2, {scale, 0, 0, scale, scale, scale}),
			Matrix<Real>(3, 3, {scale, scale, 0, 2 * scale, 2 * scale, 0, 3 * scale, 4 * scale, 0}),
			Matrix<Real>(2, 3, {1, Real(4) / 3, 0, 2, Real(7) / 3, 0}),
			LinearSolverRecommendation::QR, scale);
	}
}

TEST_CASE("Multiple RHS returns wide minimum-norm solutions", "[LinearSystem][rhs-contract]") {
	const Matrix<Real> coefficient(2, 3, {1, 0, 1, 0, 1, 1});
	const Matrix<Real> rhs(2, 3, {3, 1, 0, 0, 4, 0});
	CheckMultipleSolution(coefficient, rhs,
		Matrix<Real>(3, 3, {2, -Real(2) / 3, 0, -1, Real(7) / 3, 0, 1, Real(5) / 3, 0}),
		LinearSolverRecommendation::SVD);
	const auto actual = LinearSystem<Real>(coefficient, rhs).SolveMultiple();
	for (int col = 0; col < rhs.cols(); ++col)
		CHECK_THAT(actual(0, col) + actual(1, col) - actual(2, col),
			WithinAbs(Real(0), Precision::EqualityTolerance));
}

TEST_CASE("Multiple RHS handles square and tall rank deficiency", "[LinearSystem][rhs-contract]") {
	SECTION("Square") {
		CheckMultipleSolution(Matrix<Real>(2, 2, {1, 2, 2, 4}),
			Matrix<Real>(2, 3, {3, 1, 0, 6, 0, 0}),
			Matrix<Real>(2, 3, {Real(3) / 5, Real(1) / 25, 0, Real(6) / 5, Real(2) / 25, 0}),
			LinearSolverRecommendation::SVD);
	}
	SECTION("Tall") {
		CheckMultipleSolution(Matrix<Real>(3, 2, {1, 2, 2, 4, 3, 6}),
			Matrix<Real>(3, 3, {3, 1, 0, 6, 0, 0, 9, 0, 0}),
			Matrix<Real>(2, 3, {Real(3) / 5, Real(1) / 70, 0, Real(6) / 5, Real(1) / 35, 0}),
			LinearSolverRecommendation::SVD);
	}
	SECTION("Zero operator") {
		CheckMultipleSolution(Matrix<Real>(2, 3), Matrix<Real>(2, 2, {1, 0, 2, 0}),
			Matrix<Real>(3, 2), LinearSolverRecommendation::SVD);
	}
}

TEST_CASE("Multiple RHS retains owned factors through copy move and assignment", "[LinearSystem][rhs-contract]") {
	const auto createSolvedSystem = [] {
		LinearSystem<Real> system(Matrix<Real>(2, 2, {4, 1, 2, 3}),
			Matrix<Real>(2, 2, {6, 7, 8, 1}));
		(void)system.SolveMultiple();
		return system;
	};
	auto original = createSolvedSystem();
	auto copied = original;
	auto moved = std::move(original);
	LinearSystem<Real> assigned(Matrix<Real>::Identity(1), Vector<Real>{1});
	assigned = copied;
	LinearSystem<Real> moveAssigned(Matrix<Real>::Identity(1));
	moveAssigned = createSolvedSystem();
	for (const auto* system : {&copied, &moved, &assigned, &moveAssigned}) {
		const auto actual = system->SolveMultiple();
		CHECK_THAT(actual(0, 0), WithinAbs(Real(1), Precision::EqualityTolerance));
		CHECK_THAT(actual(1, 0), WithinAbs(Real(2), Precision::EqualityTolerance));
		CHECK_THAT(actual(0, 1), WithinAbs(Real(2), Precision::EqualityTolerance));
		CHECK_THAT(actual(1, 1), WithinAbs(Real(-1), Precision::EqualityTolerance));
		CHECK_THROWS_AS(system->Solve(), InvalidStateError);
		REQUIRE(system->getRHSMatrix().cols() == 2);
	}
	assigned = LinearSystem<Real>(Matrix<Real>::Identity(1), Vector<Real>{7});
	CHECK(assigned.getRHS()[0] == Real(7));
	CHECK(assigned.Solve()[0] == Real(7));
	CHECK_THROWS_AS(assigned.SolveMultiple(), InvalidStateError);
	CHECK_THROWS_AS(assigned.getRHSMatrix(), InvalidStateError);
}

TEST_CASE("Multiple RHS preserves empty column shapes and constructor dimensions", "[LinearSystem][rhs-contract]") {
	for (const Matrix<Real>& coefficient : {
		Matrix<Real>(2, 2, {2, 1, 0, 3}),
		Matrix<Real>(2, 2, {4, 1, 1, 3}),
		Matrix<Real>(2, 2, {4, 1, 2, 3}),
		Matrix<Real>(3, 2, {1, 0, 0, 1, 1, 1}),
		Matrix<Real>(2, 3, {1, 0, 1, 0, 1, 1}),
		Matrix<Real>(2, 2, {1, 2, 2, 4})}) {
		const LinearSystem<Real> system(coefficient, Matrix<Real>(coefficient.rows(), 0));
		const auto actual = system.SolveMultiple();
		CHECK(actual.rows() == coefficient.cols());
		CHECK(actual.cols() == 0);
		CHECK(system.getRHSMatrix().rows() == coefficient.rows());
		CHECK(system.getRHSMatrix().cols() == 0);
		CheckVectorOperationsReject(system);
		const auto analysis = system.Analyze();
		CHECK(analysis.solutionStatuses.empty());
		CHECK_FALSE(analysis.multipleRHSAggregate.has_value());
		CHECK(analysis.recommendedSolver.has_value());
		CHECK_THROWS_AS(LinearSystem<Real>(coefficient, Matrix<Real>(coefficient.rows() + 1, 1)), MatrixDimensionError);
	}
	CHECK_THROWS_AS(LinearSystem<Real>(Matrix<Real>(), Matrix<Real>()).SolveMultiple(), MatrixDimensionError);
}

TEST_CASE("Multiple RHS rejects nonfinite data and unrepresentable solutions", "[LinearSystem][rhs-contract]") {
	for (const Matrix<Real>& coefficient : {
		Matrix<Real>(2, 2, {2, 1, 0, 3}),
		Matrix<Real>(2, 2, {4, 1, 1, 3}),
		Matrix<Real>(2, 2, {4, 1, 2, 3}),
		Matrix<Real>(3, 2, {1, 0, 0, 1, 1, 1}),
		Matrix<Real>(2, 3, {1, 0, 1, 0, 1, 1}),
		Matrix<Real>(2, 2, {1, 2, 2, 4})}) {
		for (const Real invalid : {std::numeric_limits<Real>::quiet_NaN(),
			std::numeric_limits<Real>::infinity(), -std::numeric_limits<Real>::infinity()}) {
			Matrix<Real> rhs(coefficient.rows(), 2);
			rhs(coefficient.rows() - 1, 1) = invalid;
			CHECK_THROWS_AS(LinearSystem<Real>(coefficient, rhs).SolveMultiple(), NumericInputError);
		}
	}
	const Real underflow = std::numeric_limits<Real>::min();
	const Real maximum = std::numeric_limits<Real>::max();
	CHECK_THROWS_AS(LinearSystem<Real>(Matrix<Real>(1, 1, underflow),
		Matrix<Real>(1, 1, maximum)).SolveMultiple(), MatrixNumericalError);
}

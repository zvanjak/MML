#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "../TestPrecision.h"

#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/systems/LinearSystem.h>
#endif

#include <array>
#include <cmath>
#include <limits>

using namespace MML;
using namespace MML::Systems;
using Catch::Matchers::WithinAbs;

namespace {
	std::array<Real, 4> DispatchScales() {
		constexpr int exponent = std::numeric_limits<Real>::max_exponent / 4;
		return {std::ldexp(Real(1), -exponent), Real(1e-20), Real(1), std::ldexp(Real(1), exponent)};
	}

	void CheckDispatchSolution(const Matrix<Real>& coefficient, const Vector<Real>& rhs,
		const Vector<Real>& expected, Real scale, LinearSolverRecommendation recommendation) {
		const LinearSystem<Real> system(coefficient, rhs);
		CHECK(system.Analyze().recommendedSolver == recommendation);
		const auto actual = system.Solve();
		const auto quick = SolveLinearSystem(coefficient, rhs);
		for (int col = 0; col < coefficient.cols(); ++col) {
			CHECK_THAT(actual[col], WithinAbs(expected[col], Precision::EqualityTolerance));
			CHECK_THAT(quick[col], WithinAbs(expected[col], Precision::EqualityTolerance));
		}

		Real residualNorm = 0;
		Real rhsNorm = 0;
		for (int row = 0; row < coefficient.rows(); ++row) {
			Real residual = -rhs[row] / scale;
			for (int col = 0; col < coefficient.cols(); ++col)
				residual += (coefficient(row, col) / scale) * actual[col];
			residualNorm = std::hypot(residualNorm, residual);
			rhsNorm = std::hypot(rhsNorm, rhs[row] / scale);
		}
		CHECK(residualNorm / rhsNorm <= Precision::EqualityTolerance);
		CHECK(system.RelativeResidual(actual) <= Precision::EqualityTolerance);
	}
}

TEST_CASE("Linear dispatch does not discard dense entries under uniform scaling", "[LinearSystem][dispatch-contract]") {
	for (const Real scale : DispatchScales()) {
		INFO("scale " << scale);
		const Matrix<Real> coefficient(2, 2, {2 * scale, scale, scale, 2 * scale});
		const Vector<Real> rhs{3 * scale, 3 * scale};
		const LinearSystem<Real> system(coefficient, rhs);
		CHECK(system.Analyze().recommendedSolver != LinearSolverRecommendation::Triangular);
		const auto actual = system.Solve();
		CHECK_THAT(actual[0], WithinAbs(Real(1), Precision::EqualityTolerance));
		CHECK_THAT(actual[1], WithinAbs(Real(1), Precision::EqualityTolerance));
		CHECK(system.RelativeResidual(actual) <= Precision::EqualityTolerance);
	}
}

TEST_CASE("Linear dispatch preserves nonsymmetric equations across units and signs", "[LinearSystem][dispatch-contract]") {
	for (const Real magnitude : DispatchScales()) {
		for (const Real sign : {Real(-1), Real(1)}) {
			const Real scale = sign * magnitude;
			INFO("scale " << scale);
			CheckDispatchSolution(Matrix<Real>(2, 2, {4 * scale, scale, 2 * scale, 3 * scale}),
				Vector<Real>{6 * scale, 8 * scale}, Vector<Real>{1, 2}, scale, LinearSolverRecommendation::LU);
			CheckDispatchSolution(Matrix<Real>(2, 2, {4 * scale, 2 * scale, scale, 3 * scale}),
				Vector<Real>{8 * scale, 7 * scale}, Vector<Real>{1, 2}, scale, LinearSolverRecommendation::LU);
		}
	}
}

TEST_CASE("Linear dispatch retains both exact triangular orientations across units", "[LinearSystem][dispatch-contract]") {
	for (const Real magnitude : DispatchScales()) {
		for (const Real sign : {Real(-1), Real(1)}) {
			const Real scale = sign * magnitude;
			INFO("scale " << scale);
			CheckDispatchSolution(Matrix<Real>(2, 2, {2 * scale, scale, 0, 3 * scale}),
				Vector<Real>{4 * scale, 6 * scale}, Vector<Real>{1, 2}, scale, LinearSolverRecommendation::Triangular);
			CheckDispatchSolution(Matrix<Real>(2, 2, {2 * scale, 0, scale, 3 * scale}),
				Vector<Real>{2 * scale, 7 * scale}, Vector<Real>{1, 2}, scale, LinearSolverRecommendation::Triangular);
			CheckDispatchSolution(Matrix<Real>(2, 2, {2 * scale, -Real(0), -Real(0), 3 * scale}),
				Vector<Real>{2 * scale, 6 * scale}, Vector<Real>{1, 2}, scale, LinearSolverRecommendation::Triangular);
		}
	}
}

TEST_CASE("Linear dispatch keeps approximate triangular queries separate from structure", "[LinearSystem][dispatch-contract]") {
	const Real nonzero = Precision::IsMatrixDiagonalTolerance / 2;
	for (const Matrix<Real>& coefficient : {
		Matrix<Real>(2, 2, {2, 1, nonzero, 3}),
		Matrix<Real>(2, 2, {2, nonzero, 1, 3})}) {
		const Vector<Real> expected{1, 2};
		const Vector<Real> rhs = coefficient * expected;
		const LinearSystem<Real> system(coefficient, rhs);
		CHECK((system.IsUpperTriangular() || system.IsLowerTriangular()));
		CHECK_FALSE(system.IsUpperTriangular(LinearSystem<Real>::ComparisonTolerance{}));
		CHECK_FALSE(system.IsLowerTriangular(LinearSystem<Real>::ComparisonTolerance{}));
		CheckDispatchSolution(coefficient, rhs, expected, Real(1), LinearSolverRecommendation::LU);
	}
}

TEST_CASE("Linear dispatch requires exact symmetry before selecting Cholesky", "[LinearSystem][dispatch-contract]") {
	const Real adjacent = std::nextafter(Real(1), Real(2));
	const Matrix<Real> coefficient(2, 2, {4, adjacent, 1, 3});
	const Vector<Real> expected{1, 2};
	const Vector<Real> rhs = coefficient * expected;
	const LinearSystem<Real> system(coefficient, rhs);
	CHECK(system.IsSymmetric());
	CHECK_FALSE(system.IsSymmetric(LinearSystem<Real>::ComparisonTolerance{}));
	CHECK(system.Analyze().matrix.isSymmetric);
	CheckDispatchSolution(coefficient, rhs, expected, Real(1), LinearSolverRecommendation::LU);
}

TEST_CASE("Linear dispatch preserves genuine SPD and indefinite recommendations", "[LinearSystem][dispatch-contract]") {
	CheckDispatchSolution(Matrix<Real>(2, 2, {4, 1, 1, 3}), Vector<Real>{6, 7},
		Vector<Real>{1, 2}, Real(1), LinearSolverRecommendation::Cholesky);
	for (const Real scale : DispatchScales()) {
		INFO("scale " << scale);
		CheckDispatchSolution(Matrix<Real>(2, 2, {2 * scale, scale, scale, -2 * scale}),
			Vector<Real>{4 * scale, -3 * scale}, Vector<Real>{1, 2}, scale, LinearSolverRecommendation::LU);
	}
}

TEST_CASE("Linear dispatch retains SVD priority over singular triangular structure", "[LinearSystem][dispatch-contract]") {
	for (const Real scale : DispatchScales()) {
		INFO("scale " << scale);
		CheckDispatchSolution(Matrix<Real>(2, 2, {scale, 2 * scale, 0, 0}),
			Vector<Real>{3 * scale, 0}, Vector<Real>{Real(3) / 5, Real(6) / 5},
			scale, LinearSolverRecommendation::SVD);
		CheckDispatchSolution(Matrix<Real>(1, 2, {scale, 2 * scale}),
			Vector<Real>{3 * scale}, Vector<Real>{Real(3) / 5, Real(6) / 5},
			scale, LinearSolverRecommendation::SVD);
	}
}

TEST_CASE("Linear dispatch retains tall full-rank QR across units", "[LinearSystem][dispatch-contract]") {
	for (const Real scale : DispatchScales()) {
		INFO("scale " << scale);
		CheckDispatchSolution(Matrix<Real>(3, 2, {scale, 0, 0, scale, scale, scale}),
			Vector<Real>{scale, 2 * scale, 3 * scale}, Vector<Real>{1, 2},
			scale, LinearSolverRecommendation::QR);
	}
}

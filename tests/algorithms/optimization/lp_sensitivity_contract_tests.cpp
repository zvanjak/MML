#include <catch2/catch_all.hpp>
#include "../../TestPrecision.h"

#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/algorithms/Optimization/LinearProgramming.h>
#endif

#include <array>
#include <cmath>
#include <limits>

using namespace MML;
using namespace MML::Optimization;
using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinAbs;

namespace {
	void CheckValue(Real actual, Real expected) {
		if (std::isinf(expected)) CHECK(actual == expected);
		else CHECK_THAT(actual, WithinAbs(expected, Precision::DefaultToleranceStrict *
			std::max(Real(1), std::abs(expected))));
	}

	void CheckRange(const std::pair<Real, Real>& actual, Real lower, Real upper) {
		CheckValue(actual.first, lower);
		CheckValue(actual.second, upper);
	}

	const Real infinity = std::numeric_limits<Real>::infinity();

	// For a 2x2 matrix, sigma_max^2 / |det| = sigma_max / sigma_min.
	Real Condition2(Real a, Real b, Real c, Real d) {
		const Real trace = a * a + b * b + c * c + d * d;
		const Real determinant = a * d - b * c;
		return (trace + std::sqrt(trace * trace - 4 * determinant * determinant)) /
			(2 * std::abs(determinant));
	}

	std::pair<Real, Real> PerturbationRange(const std::array<Real, 2>& values,
		const std::array<Real, 2>& direction) {
		Real lower = -infinity, upper = infinity;
		for (int i = 0; i < 2; ++i) {
			if (direction[i] > 0) lower = std::max(lower, -values[i] / direction[i]);
			if (direction[i] < 0) upper = std::min(upper, -values[i] / direction[i]);
		}
		return {lower, upper};
	}
}

TEST_CASE("LP sensitivity uses original rather than pivoted RHS values", "[LP][sensitivity][lp-sensitivity-contract]")
{
	LinearProgram lp(1);
	lp.SetObjective({2}, LPObjective::Maximize);
	lp.AddConstraint({1}, LPConstraintType::LessEqual, 3);
	lp.AddConstraint({1}, LPConstraintType::LessEqual, 5);
	const auto result = SolveLPWithSensitivity(lp);
	REQUIRE(result.IsOptimal());
	REQUIRE(result.rhsRanges.size() == 2);
	CheckRange(result.rhsRanges[0], 0, 5);
	CheckRange(result.rhsRanges[1], 3, infinity);
	CheckRange(result.objectiveRanges[0], 0, infinity);
	CheckValue(result.dualValues[0], 2);
	CheckValue(result.dualValues[1], 0);
	CheckValue(result.basisConditionNumber, Condition2(1, 0, 1, 1));

	for (Real bound : {Real(2), Real(3), Real(4), Real(9)}) {
		LinearProgram changed(1);
		changed.SetObjective({2}, LPObjective::Maximize);
		changed.AddConstraint({1}, LPConstraintType::LessEqual, 3);
		changed.AddConstraint({1}, LPConstraintType::LessEqual, bound);
		const auto solved = SolveLP(changed);
		REQUIRE(solved.IsOptimal());
		CheckValue(solved.objectiveValue, 2 * std::min(Real(3), bound));
		if (bound > 3) CHECK(solved.basisIndices == result.basisIndices);
		if (bound < 3) CHECK(solved.basisIndices != result.basisIndices);
	}
}

TEST_CASE("LP basis conditioning measures original diagonal and sheared bases", "[LP][sensitivity][lp-sensitivity-contract]")
{
	for (Real small : {Real(1e-6), Real(1) / 1024, Real(1)}) {
		for (Real scale : {Real(1), Real(1024)}) {
			for (bool shear : {false, true}) {
				LinearProgram lp(2);
				lp.SetObjective({1, 2}, LPObjective::Maximize);
				lp.AddConstraint({scale, shear ? scale : Real(0)}, LPConstraintType::LessEqual,
					scale * (shear ? 5 : 2));
				lp.AddConstraint({0, scale * small}, LPConstraintType::LessEqual, scale * small * 3);
				LPConfig config;
				// This witness asks to pivot a small, nonzero basis entry, not to
				// classify it using the ordinary unit-scale pivot threshold.
				config.pivotTolerance = std::min(config.pivotTolerance, scale * small / 2);
				const auto result = SolveLPWithSensitivity(lp, config);
				REQUIRE(result.IsOptimal());
				CheckValue(result.x[0], 2);
				CheckValue(result.x[1], 3);
				const Real expected = shear ? Condition2(1, 1, 0, small) : 1 / small;
				CHECK_THAT(result.basisConditionNumber, WithinRel(expected, Precision::DefaultToleranceStrict));
			}
		}
	}
}

TEST_CASE("LP diagnostics respect equality surplus and original row signs", "[LP][sensitivity][lp-sensitivity-contract]")
{
	struct Witness {
		LPObjective sense;
		LPConstraintType type;
		Real a, b, cost, shadow, rhsLower, rhsUpper, costLower, costUpper;
	};
	const std::array<Witness, 8> witnesses{{
		{LPObjective::Minimize, LPConstraintType::Equal, 1, 3, 2, 2, 0, infinity, -infinity, infinity},
		{LPObjective::Minimize, LPConstraintType::GreaterEqual, 1, 3, 2, 2, 0, infinity, 0, infinity},
		{LPObjective::Minimize, LPConstraintType::LessEqual, -1, -3, 2, -2, -infinity, 0, 0, infinity},
		{LPObjective::Maximize, LPConstraintType::GreaterEqual, -1, -3, 2, -2, -infinity, 0, 0, infinity},
		{LPObjective::Maximize, LPConstraintType::LessEqual, 1, 3, 2, 2, 0, infinity, 0, infinity},
		{LPObjective::Minimize, LPConstraintType::LessEqual, 1, 3, -2, -2, 0, infinity, -infinity, 0},
		{LPObjective::Minimize, LPConstraintType::LessEqual, 1, 3, 2, 0, 0, infinity, 0, infinity},
		{LPObjective::Maximize, LPConstraintType::GreaterEqual, 1, 3, -2, -2, 0, infinity, -infinity, 0}
	}};
	for (const auto& w : witnesses) {
		LinearProgram lp(1);
		lp.SetObjective({w.cost}, w.sense);
		lp.AddConstraint({w.a}, w.type, w.b);
		const auto result = SolveLPWithSensitivity(lp);
		REQUIRE(result.IsOptimal());
		CheckRange(result.rhsRanges[0], w.rhsLower, w.rhsUpper);
		CheckRange(result.objectiveRanges[0], w.costLower, w.costUpper);
		CheckValue(result.dualValues[0], w.shadow);
		CheckValue(result.basisConditionNumber, 1);
		const auto ordinary = SolveLP(lp);
		REQUIRE(ordinary.IsOptimal());
		CheckValue(ordinary.dualValues[0], w.shadow);
	}
}

TEST_CASE("LP objective ranges are absolute coefficients in the original objective sense", "[LP][sensitivity][lp-sensitivity-contract]")
{
	for (bool maximize : {false, true}) {
		const Real sign = maximize ? Real(1) : Real(-1);
		LinearProgram lp(3);
		lp.SetObjective({3 * sign, 2 * sign, -sign}, maximize ? LPObjective::Maximize : LPObjective::Minimize);
		lp.AddConstraint({1, 1, 1}, LPConstraintType::LessEqual, 4);
		lp.AddConstraint({2, 1, 1}, LPConstraintType::LessEqual, 5);
		const auto result = SolveLPWithSensitivity(lp);
		REQUIRE(result.IsOptimal());
		CheckValue(result.x[0], 1);
		CheckValue(result.x[1], 3);
		CheckValue(result.x[2], 0);
		CheckRange(result.objectiveRanges[0], maximize ? Real(2) : Real(-4), maximize ? Real(4) : Real(-2));
		CheckRange(result.objectiveRanges[1], maximize ? Real(1.5) : Real(-3), maximize ? Real(3) : Real(-1.5));
		CheckRange(result.objectiveRanges[2], maximize ? -infinity : Real(-2), maximize ? Real(2) : infinity);
		CheckRange(result.rhsRanges[0], Real(2.5), 5);
		CheckRange(result.rhsRanges[1], 4, 8);
		CheckValue(result.dualValues[0], sign);
		CheckValue(result.dualValues[1], sign);
		CheckValue(result.reducedCosts[2], -3 * sign);
	}
}

TEST_CASE("LP diagnostics match independent inverse and dual certificates across units", "[LP][sensitivity][lp-sensitivity-contract]")
{
	const std::array<std::array<Real, 4>, 3> matrices{{{1, 1, 2, 1}, {2, 1, 1, 3}, {1, 2, 3, 1}}};
	for (const auto& entries : matrices)
		for (bool reverseRows : {false, true})
			for (bool maximize : {false, true})
				for (int flips = 0; flips < 4; ++flips)
					for (Real scale : {Real(0.25), Real(1), Real(1024)}) {
						const int r0 = reverseRows ? 2 : 0, r1 = reverseRows ? 0 : 2;
						const Real a = entries[r0], b = entries[r0 + 1], c = entries[r1], d = entries[r1 + 1];
						const Real determinant = a * d - b * c;
						const std::array<std::array<Real, 2>, 2> inverse{{{d / determinant, -b / determinant},
							{-c / determinant, a / determinant}}};
						const std::array<Real, 2> x{2, 3}, dual{1, 2};
						const std::array<Real, 2> rhs{a * x[0] + b * x[1], c * x[0] + d * x[1]};
						const std::array<Real, 2> cost{a * dual[0] + c * dual[1], b * dual[0] + d * dual[1]};
						const Real objectiveSign = maximize ? Real(1) : Real(-1);
						const std::array<Real, 2> rowSign{flips & 1 ? Real(-1) : Real(1), flips & 2 ? Real(-1) : Real(1)};
						LinearProgram lp(2);
						lp.SetObjective({objectiveSign * cost[0], objectiveSign * cost[1]},
							maximize ? LPObjective::Maximize : LPObjective::Minimize);
						lp.AddConstraint({rowSign[0] * scale * a, rowSign[0] * scale * b},
							rowSign[0] > 0 ? LPConstraintType::LessEqual : LPConstraintType::GreaterEqual,
							rowSign[0] * scale * rhs[0]);
						lp.AddConstraint({rowSign[1] * scale * c, rowSign[1] * scale * d},
							rowSign[1] > 0 ? LPConstraintType::LessEqual : LPConstraintType::GreaterEqual,
							rowSign[1] * scale * rhs[1]);
						const auto result = SolveLPWithSensitivity(lp);
						REQUIRE(result.IsOptimal());
						for (int j = 0; j < 2; ++j) {
							CheckValue(result.x[j], x[j]);
							CheckValue(result.dualValues[j], objectiveSign * rowSign[j] * dual[j] / scale);
							const auto delta = PerturbationRange(x, {inverse[0][j], inverse[1][j]});
							const Real lower = scale * (rhs[j] + delta.first), upper = scale * (rhs[j] + delta.second);
							CheckRange(result.rhsRanges[j], rowSign[j] > 0 ? lower : -upper, rowSign[j] > 0 ? upper : -lower);
							const auto objectiveDelta = PerturbationRange(dual, inverse[j]);
							const Real costLower = cost[j] + objectiveDelta.first, costUpper = cost[j] + objectiveDelta.second;
							CheckRange(result.objectiveRanges[j], maximize ? costLower : -costUpper, maximize ? costUpper : -costLower);
						}
						CheckValue(result.basisConditionNumber, Condition2(a, b, c, d));
						CheckValue(result.objectiveValue, objectiveSign * (cost[0] * x[0] + cost[1] * x[1]));
					}
}

TEST_CASE("LP raw tableau diagnostics survive pivots objective changes and reinitialization", "[LP][sensitivity][lp-sensitivity-contract]")
{
	SimplexTableau tableau;
	tableau.Initialize(Matrix<Real>(2, 3, {1, 1, 0, 1, 0, 1}), Vector<Real>{3, 5}, Vector<Real>{-2, 0, 0}, {}, 1, 2, 0);
	tableau.Pivot(0, 0);
	CheckRange(tableau.GetRHSRange(0), 0, 5);
	CheckRange(tableau.GetRHSRange(1), 3, infinity);
	CheckValue(tableau.ComputeBasisCondition(), Condition2(1, 0, 1, 1));
	tableau.SetupPhase2(Vector<Real>{-4, 0, 0});
	CheckValue(tableau.GetShadowPrice(0), -4);
	CheckValue(tableau.GetShadowPrice(1), 0);
	CheckRange(tableau.GetObjectiveCoefficientRange(0), -infinity, 0);

	tableau.Initialize(Matrix<Real>(2, 4, {2, 0, 1, 0, 0, Real(0.5), 0, 1}),
		Vector<Real>{4, 1}, Vector<Real>{-1, -1, 0, 0}, {}, 2, 2, 0);
	tableau.Pivot(0, 0);
	tableau.Pivot(1, 1);
	CheckRange(tableau.GetRHSRange(0), 0, infinity);
	CheckRange(tableau.GetRHSRange(1), 0, infinity);
	CheckValue(tableau.GetShadowPrice(0), Real(-0.5));
	CheckValue(tableau.GetShadowPrice(1), -2);
	CheckValue(tableau.ComputeBasisCondition(), 4);
}

TEST_CASE("LP RHS sensitivity retains small inverse basis coefficients", "[LP][sensitivity][lp-sensitivity-contract]")
{
	for (Real scale : {Real(1), Real(1e12), Testing::Tol(Real(1e100), Real(1e30))}) {
		LinearProgram lp(1);
		lp.SetObjective({1}, LPObjective::Maximize);
		lp.AddConstraint({scale}, LPConstraintType::LessEqual, scale);
		const auto result = SolveLPWithSensitivity(lp);
		REQUIRE(result.IsOptimal());
		CheckRange(result.rhsRanges[0], 0, infinity);
		CHECK_THAT(result.dualValues[0], WithinRel(1 / scale, Precision::DefaultToleranceStrict));
		CheckValue(result.basisConditionNumber, 1);
	}
}

TEST_CASE("LP condition ratios do not require representable unscaled singular values", "[LP][sensitivity][lp-sensitivity-contract]")
{
	for (Real scale : {Testing::Tol(Real(1e-100), Real(1e-30)), Real(1), std::numeric_limits<Real>::max() * Real(0.75)}) {
		LPConfig config;
		config.pivotTolerance = std::numeric_limits<Real>::min();
		SimplexTableau tableau;
		tableau.Initialize(Matrix<Real>(2, 4, {scale, scale, 1, 0, scale, scale * (Real(2) / 3), 0, 1}),
			Vector<Real>{0, 0}, Vector<Real>{0, 0, 0, 0}, {}, 2, 2, 0, config);
		tableau.Pivot(0, 0);
		tableau.Pivot(1, 1);
		CheckValue(tableau.ComputeBasisCondition(), Condition2(1, 1, 1, Real(2) / 3));
	}
}

TEST_CASE("LP mixed equality surplus and slack rows retain independent sensitivity directions", "[LP][sensitivity][lp-sensitivity-contract]")
{
	for (int flips = 0; flips < 4; ++flips)
		for (bool maximize : {false, true}) {
			const Real equalitySign = flips & 1 ? Real(-1) : Real(1);
			const Real inequalitySign = flips & 2 ? Real(-1) : Real(1);
			const Real objectiveSign = maximize ? Real(-1) : Real(1);
			LinearProgram lp(2);
			lp.SetObjective({objectiveSign, 2 * objectiveSign}, maximize ? LPObjective::Maximize : LPObjective::Minimize);
			lp.AddConstraint({equalitySign, equalitySign}, LPConstraintType::Equal, 4 * equalitySign);
			lp.AddConstraint({inequalitySign, 0},
				inequalitySign > 0 ? LPConstraintType::GreaterEqual : LPConstraintType::LessEqual, inequalitySign);
			lp.AddConstraint({0, 1}, LPConstraintType::LessEqual, 2);
			const auto result = SolveLPWithSensitivity(lp);
			REQUIRE(result.IsOptimal());
			CheckValue(result.x[0], 4);
			CheckValue(result.x[1], 0);
			CheckRange(result.rhsRanges[0], equalitySign > 0 ? Real(1) : -infinity,
				equalitySign > 0 ? infinity : Real(-1));
			CheckRange(result.rhsRanges[1], inequalitySign > 0 ? -infinity : Real(-4),
				inequalitySign > 0 ? Real(4) : infinity);
			CheckRange(result.rhsRanges[2], 0, infinity);
			CheckRange(result.objectiveRanges[0], maximize ? Real(-2) : -infinity,
				maximize ? infinity : Real(2));
			CheckRange(result.objectiveRanges[1], maximize ? -infinity : Real(1),
				maximize ? Real(-1) : infinity);
			CheckValue(result.dualValues[0], equalitySign * objectiveSign);
			CheckValue(result.dualValues[1], 0);
			CheckValue(result.dualValues[2], 0);
			CheckValue(result.basisConditionNumber, Condition2(1, 0, 1, 1));
		}
}

TEST_CASE("LP redundant equalities cannot gain feasibility through zero basic artificials", "[LP][sensitivity][lp-sensitivity-contract]")
{
	for (Real sign : {Real(-1), Real(1)}) {
		LinearProgram lp(1);
		lp.SetObjective({1});
		lp.AddConstraint({1}, LPConstraintType::Equal, 1);
		lp.AddConstraint({2 * sign}, LPConstraintType::Equal, 2 * sign);
		const auto result = SolveLPWithSensitivity(lp);
		REQUIRE(result.IsOptimal());
		CheckRange(result.rhsRanges[0], 1, 1);
		CheckRange(result.rhsRanges[1], 2 * sign, 2 * sign);
		CheckRange(result.objectiveRanges[0], -infinity, infinity);
		CheckValue(result.dualValues[0] + 2 * sign * result.dualValues[1], 1);
	}
}

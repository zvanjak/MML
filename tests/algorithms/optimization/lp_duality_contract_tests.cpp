#include <catch2/catch_test_macros.hpp>
#include "../../TestPrecision.h"

#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/algorithms/Optimization/LinearProgramming.h>
#endif

using namespace MML;
using namespace MML::Optimization;

namespace {
	void CheckClose(Real actual, Real expected)
	{
		CHECK(std::abs(actual - expected) <= Precision::DefaultToleranceStrict * std::max(Real(1), std::abs(expected)));
	}

	Vector<Real> CheckOptimalPair(const LinearProgram& primal, Real expectedObjective)
	{
		const auto dual = primal.Dual();
		const auto p = SolveLP(primal);
		const auto d = SolveLP(dual);
		REQUIRE(p.IsOptimal());
		REQUIRE(d.IsOptimal());
		REQUIRE(p.status == AlgorithmStatus::Success);
		REQUIRE(d.status == AlgorithmStatus::Success);
		REQUIRE(p.x.size() == primal.NumVariables());
		REQUIRE(d.x.size() == dual.NumVariables());
		CheckClose(p.objectiveValue, expectedObjective);
		CheckClose(d.objectiveValue, expectedObjective);
		for (int j = 0; j < d.x.size(); ++j)
			CHECK(d.x[j] >= -Precision::DefaultToleranceStrict);

		Vector<Real> y(primal.NumConstraints());
		int negativePart = primal.NumConstraints();
		Real dualValue = 0;
		for (int i = 0; i < primal.NumConstraints(); ++i) {
			const auto& row = primal.GetConstraint(i);
			y[i] = d.x[i];
			if (row.type == LPConstraintType::Equal) {
				y[i] -= d.x[negativePart++];
			} else {
				const bool nonpositive = (primal.ObjectiveSense() == LPObjective::Minimize) ==
					(row.type == LPConstraintType::LessEqual);
				if (nonpositive)
					y[i] = -y[i];
				CHECK((nonpositive ? -y[i] : y[i]) >= -Precision::DefaultToleranceStrict);
			}
			Real lhs = 0;
			for (int j = 0; j < p.x.size(); ++j)
				lhs += row.coefficients[j] * p.x[j];
			const Real slack = row.rhs - lhs;
			if (row.type == LPConstraintType::Equal)
				CheckClose(slack, 0);
			else
				CHECK((row.type == LPConstraintType::LessEqual ? slack : -slack) >= -Precision::DefaultToleranceStrict);
			CheckClose(y[i] * slack, 0);
			dualValue += row.rhs * y[i];
		}
		Real primalValue = 0;
		for (int j = 0; j < p.x.size(); ++j) {
			CHECK(p.x[j] >= -Precision::DefaultToleranceStrict);
			Real transposeProduct = 0;
			for (int i = 0; i < y.size(); ++i)
				transposeProduct += primal.GetConstraint(i).coefficients[j] * y[i];
			const Real reduced = primal.ObjectiveCoeffs()[j] - transposeProduct;
			CHECK((primal.ObjectiveSense() == LPObjective::Minimize ? reduced : -reduced) >= -Precision::DefaultToleranceStrict);
			CheckClose(p.x[j] * reduced, 0);
			primalValue += primal.ObjectiveCoeffs()[j] * p.x[j];
		}
		CheckClose(primalValue, expectedObjective);
		CheckClose(dualValue, expectedObjective);
		CheckClose(primalValue, dualValue);
		return y;
	}
}

TEST_CASE("LP dual: positive RHS inequality witness and control", "[LP][dual][lp-duality]")
{
	for (const auto type : {LPConstraintType::LessEqual, LPConstraintType::GreaterEqual}) {
		CAPTURE(type);
		LinearProgram primal(1);
		primal.SetObjective({1});
		primal.AddConstraint({1}, type, 1);
		const auto primalResult = SolveLP(primal);
		const auto dualResult = SolveLP(primal.Dual());
		REQUIRE(primalResult.IsOptimal());
		REQUIRE(dualResult.IsOptimal());
		const Real expected = type == LPConstraintType::LessEqual ? Real(0) : Real(1);
		CHECK(std::abs(primalResult.objectiveValue - expected) <= Precision::EqualityTolerance);
		CHECK(std::abs(dualResult.objectiveValue - expected) <= Precision::EqualityTolerance);
	}
}

TEST_CASE("LP dual: inequality multipliers transform both objective senses", "[LP][dual][lp-duality]")
{
	for (const auto sense : {LPObjective::Minimize, LPObjective::Maximize}) {
		for (const auto type : {LPConstraintType::LessEqual, LPConstraintType::GreaterEqual}) {
			CAPTURE(sense, type);
			LinearProgram primal(1);
			primal.SetObjective({3}, sense);
			primal.AddConstraint({2}, type, 4);
			const auto dual = primal.Dual();
			REQUIRE(dual.NumVariables() == 1);
			REQUIRE(dual.NumConstraints() == 1);
			const bool negative = (sense == LPObjective::Minimize && type == LPConstraintType::LessEqual) ||
				(sense == LPObjective::Maximize && type == LPConstraintType::GreaterEqual);
			const Real sign = negative ? Real(-1) : Real(1);
			CHECK(dual.ObjectiveCoeffs()[0] == sign * 4);
			CHECK(dual.GetConstraint(0).coefficients[0] == sign * 2);
			CHECK(dual.GetConstraint(0).rhs == 3);
			CHECK(dual.GetConstraint(0).type == (sense == LPObjective::Minimize ?
				LPConstraintType::LessEqual : LPConstraintType::GreaterEqual));
			CHECK(dual.ObjectiveSense() != sense);
		}
	}
}

TEST_CASE("LP dual: equality multipliers have independent positive and negative parts", "[LP][dual][lp-duality]")
{
	for (const auto sense : {LPObjective::Minimize, LPObjective::Maximize}) {
		LinearProgram primal(1);
		primal.SetObjective({3}, sense);
		primal.AddConstraint({2}, LPConstraintType::Equal, 4);
		const auto dual = primal.Dual();
		REQUIRE(dual.NumVariables() == 2);
		REQUIRE(dual.NumConstraints() == 1);
		CHECK(dual.ObjectiveCoeffs()[0] == 4);
		CHECK(dual.ObjectiveCoeffs()[1] == -4);
		CHECK(dual.GetConstraint(0).coefficients[0] == 2);
		CHECK(dual.GetConstraint(0).coefficients[1] == -2);
	}
}

TEST_CASE("LP dual: independently certified optima cover active signs and free multipliers", "[LP][dual][lp-duality]")
{
	for (const auto sense : {LPObjective::Minimize, LPObjective::Maximize}) {
		for (const Real cost : {Real(-2), Real(2)}) {
			const bool increasing = (sense == LPObjective::Maximize) == (cost > 0);
			for (const Real rowSign : {Real(-1), Real(1)}) {
				const auto boundedType = increasing == (rowSign > 0) ?
					LPConstraintType::LessEqual : LPConstraintType::GreaterEqual;
				for (const auto type : {boundedType, LPConstraintType::Equal}) {
					CAPTURE(sense, cost, rowSign, type);
					LinearProgram primal(1);
					primal.SetObjective({cost}, sense);
					primal.AddConstraint({rowSign * 2}, type, rowSign * 6);
					const auto y = CheckOptimalPair(primal, cost * 3);
					CheckClose(y[0], cost / (rowSign * 2));
				}
			}
		}
	}
}

TEST_CASE("LP dual: mixed rows preserve primary order and append equality negative parts", "[LP][dual][lp-duality]")
{
	for (const auto sense : {LPObjective::Minimize, LPObjective::Maximize}) {
		for (const Real scale : {Real(1) / 4, Real(1), Real(4)}) {
			CAPTURE(sense, scale);
			LinearProgram primal(2, "mixed");
			primal.SetObjective({scale * 3, scale}, sense);
			primal.AddConstraint({1, 1}, LPConstraintType::Equal, 4);
			primal.AddConstraint({1, 0}, LPConstraintType::LessEqual, 3);
			primal.AddConstraint({1, 0}, LPConstraintType::GreaterEqual, 1);
			primal.AddConstraint({0, 1}, LPConstraintType::LessEqual, 5);
			const Real expected = scale * (sense == LPObjective::Minimize ? 6 : 10);
			const auto y = CheckOptimalPair(primal, expected);
			CheckClose(y[0], scale);
			CheckClose(y[1], sense == LPObjective::Maximize ? scale * 2 : Real(0));
			CheckClose(y[2], sense == LPObjective::Minimize ? scale * 2 : Real(0));
			CheckClose(y[3], 0);
			const auto dual = primal.Dual();
			REQUIRE(dual.NumVariables() == 5);
			CHECK(dual.Name() == "mixed_dual");
			CHECK(dual.GetVariable(0).name == "y1_plus");
			CHECK(dual.GetVariable(4).name == "y1_minus");
			CHECK(dual.ObjectiveCoeffs()[4] == -4);
			CHECK(dual.GetConstraint(0).coefficients[4] == -1);
			CHECK(dual.GetConstraint(1).coefficients[4] == -1);
		}
	}
}

TEST_CASE("LP dual: multiple equality multipliers retain independent negative parts", "[LP][dual][lp-duality]")
{
	LinearProgram primal(2);
	primal.SetObjective({-3, 2});
	primal.AddConstraint({1, 0}, LPConstraintType::Equal, 2);
	primal.AddConstraint({0, 1}, LPConstraintType::Equal, 3);
	const auto dual = primal.Dual();
	REQUIRE(dual.NumVariables() == 4);
	CHECK(dual.ObjectiveCoeffs()[2] == -2);
	CHECK(dual.ObjectiveCoeffs()[3] == -3);
	CHECK(dual.GetConstraint(0).coefficients[2] == -1);
	CHECK(dual.GetConstraint(0).coefficients[3] == 0);
	CHECK(dual.GetConstraint(1).coefficients[2] == 0);
	CHECK(dual.GetConstraint(1).coefficients[3] == -1);
	const auto y = CheckOptimalPair(primal, 0);
	CheckClose(y[0], -3);
	CheckClose(y[1], 2);
}

TEST_CASE("LP dual: mixed negative RHS rows and a zero RHS equality solve consistently", "[LP][dual][lp-duality]")
{
	for (const auto sense : {LPObjective::Minimize, LPObjective::Maximize}) {
		for (const Real cost : {Real(-1), Real(1)}) {
			CAPTURE(sense, cost);
			LinearProgram primal(2);
			primal.SetObjective({cost, cost}, sense);
			primal.AddConstraint({-1, 0}, LPConstraintType::LessEqual, -1);
			primal.AddConstraint({0, -1}, LPConstraintType::GreaterEqual, -3);
			primal.AddConstraint({1, 1}, LPConstraintType::LessEqual, 5);
			primal.AddConstraint({1, -1}, LPConstraintType::Equal, 0);
			const bool increasing = (sense == LPObjective::Maximize) == (cost > 0);
			CheckOptimalPair(primal, cost * (increasing ? 5 : 2));
		}
	}
}

TEST_CASE("LP standard form: size auxiliary columns after RHS normalization", "[LP][standard-form][lp-duality]")
{
	for (const auto type : {LPConstraintType::LessEqual, LPConstraintType::GreaterEqual, LPConstraintType::Equal}) {
		for (const Real sign : {Real(-1), Real(1)}) {
			CAPTURE(type, sign);
			LinearProgram primal(1);
			primal.SetObjective({2});
			primal.AddConstraint({sign * 2}, type, sign * 6);
			Matrix<Real> a;
			Vector<Real> b, c;
			int slacks, artificials;
			std::vector<int> indices;
			primal.ToStandardForm(a, b, c, slacks, artificials, indices);
			const bool equality = type == LPConstraintType::Equal;
			const bool greater = !equality && ((type == LPConstraintType::GreaterEqual) == (sign > 0));
			CHECK(slacks == (equality ? 0 : 1));
			CHECK(artificials == (equality || greater ? 1 : 0));
			REQUIRE(a.rows() == 1);
			REQUIRE(a.cols() == 1 + slacks + artificials);
			REQUIRE(c.size() == a.cols());
			CHECK(a(0, 0) == 2);
			CHECK(b[0] == 6);
			CHECK(c[0] == 2);
			if (!equality)
				CHECK(a(0, 1) == (greater ? -1 : 1));
			REQUIRE(indices.size() == static_cast<size_t>(artificials));
			if (artificials) {
				CHECK(indices[0] == 1 + slacks);
				CHECK(a(0, indices[0]) == 1);
			}
		}
	}
}

TEST_CASE("LP dual: unsupported primal bounds and invalid models fail explicitly", "[LP][dual][lp-duality]")
{
	const Real inf = std::numeric_limits<Real>::infinity();
	LinearProgram primal(1);
	primal.SetObjective({1});
	primal.AddConstraint({1}, LPConstraintType::LessEqual, 1);
	CHECK_THROWS_AS(primal.SetVariableBounds(0, -inf, inf), LinearProgrammingError);
	CHECK_THROWS_AS(primal.SetVariableBounds(0, -inf, 0), LinearProgrammingError);
	CHECK_THROWS_AS(primal.SetVariableBounds(0, 1, inf), LinearProgrammingError);
	CHECK_THROWS_AS(primal.SetVariableBounds(0, 0, 2), LinearProgrammingError);
	CHECK(primal.GetVariable(0).lowerBound == 0);
	CHECK(primal.GetVariable(0).upperBound == inf);
	CHECK_THROWS_AS(LinearProgram().Dual(), LinearProgrammingError);
	primal.SetObjective({1, 2});
	CHECK_THROWS_AS(primal.Dual(), LinearProgrammingError);
}

TEST_CASE("LP dual: non-finite data and invalid row or objective enums are rejected", "[LP][dual][lp-duality]")
{
	for (const Real invalid : {std::numeric_limits<Real>::infinity(), std::numeric_limits<Real>::quiet_NaN()}) {
		for (int field = 0; field < 3; ++field) {
			LinearProgram primal(1);
			primal.SetObjective({field == 0 ? invalid : Real(1)});
			primal.AddConstraint({field == 1 ? invalid : Real(1)}, LPConstraintType::LessEqual,
				field == 2 ? invalid : Real(1));
			CHECK_THROWS_AS(primal.Dual(), LinearProgrammingError);
		}
	}
	LinearProgram primal(1);
	primal.SetObjective({1});
	primal.AddConstraint({1}, static_cast<LPConstraintType>(-1), 1);
	CHECK_THROWS_AS(primal.Dual(), LinearProgrammingError);
	LinearProgram invalidSense(1);
	invalidSense.SetObjective({1}, static_cast<LPObjective>(-1));
	invalidSense.AddConstraint({1}, LPConstraintType::LessEqual, 1);
	CHECK_THROWS_AS(invalidSense.Dual(), LinearProgrammingError);
}

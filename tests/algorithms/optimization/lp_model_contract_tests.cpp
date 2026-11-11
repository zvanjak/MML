#include <catch2/catch_test_macros.hpp>
#include "../../TestPrecision.h"

#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/algorithms/Optimization/LinearProgramming.h>
#endif

using namespace MML;
using namespace MML::Optimization;

namespace {
	void CheckRejectedModel(const LinearProgram& lp)
	{
		std::string error;
		CHECK_FALSE(lp.Validate(error));
		CHECK_FALSE(error.empty());
		Matrix<Real> a(1, 1, {7});
		Vector<Real> b{8}, c{9};
		int slacks = 4, artificials = 5;
		std::vector<int> indices{6};
		CHECK_THROWS_AS(lp.ToStandardForm(a, b, c, slacks, artificials, indices), LinearProgrammingError);
		REQUIRE(a.rows() == 1);
		REQUIRE(a.cols() == 1);
		CHECK(a(0, 0) == 7);
		REQUIRE(b.size() == 1);
		REQUIRE(c.size() == 1);
		CHECK(b[0] == 8);
		CHECK(c[0] == 9);
		CHECK(slacks == 4);
		CHECK(artificials == 5);
		CHECK(indices == std::vector<int>{6});
		CHECK_THROWS_AS(SolveLP(lp), LinearProgrammingError);
		CHECK_THROWS_AS(SolveLPTableau(lp), LinearProgrammingError);
		CHECK_THROWS_AS(SolveLPDual(lp), LinearProgrammingError);
		CHECK_THROWS_AS(SolveLPWithSensitivity(lp), LinearProgrammingError);
		CHECK_THROWS_AS(lp.Dual(), LinearProgrammingError);
	}
}

TEST_CASE("LP model: validation rejects every nonfinite data field", "[LP][lp-model]")
{
	for (const Real invalid : {std::numeric_limits<Real>::quiet_NaN(),
		std::numeric_limits<Real>::infinity(), -std::numeric_limits<Real>::infinity()}) {
		for (int field = 0; field < 3; ++field) {
			CAPTURE(field, invalid);
			LinearProgram lp(1);
			lp.SetObjective({field == 0 ? invalid : Real(1)});
			lp.AddConstraint({field == 1 ? invalid : Real(1)}, LPConstraintType::LessEqual,
				field == 2 ? invalid : Real(1));
			CheckRejectedModel(lp);
		}
	}
}

TEST_CASE("LP model: NaN RHS cannot masquerade as an optimal feasible solution", "[LP][lp-model]")
{
	LinearProgram lp(1);
	lp.SetObjective({1});
	lp.AddConstraint({1}, LPConstraintType::LessEqual, std::numeric_limits<Real>::quiet_NaN());
	CHECK_THROWS_AS(SolveLP(lp), LinearProgrammingError);
	CHECK_THROWS_AS(SolveLPTableau(lp), LinearProgrammingError);
	CHECK_THROWS_AS(SolveLPDual(lp), LinearProgrammingError);
	CHECK_THROWS_AS(SolveLPWithSensitivity(lp), LinearProgrammingError);
}

TEST_CASE("LP model: direct conversion rejects nonfinite models", "[LP][lp-model]")
{
	LinearProgram lp(1);
	lp.SetObjective({1});
	lp.AddConstraint({1}, LPConstraintType::LessEqual, std::numeric_limits<Real>::quiet_NaN());
	Matrix<Real> a;
	Vector<Real> b, c;
	int slacks = 0, artificials = 0;
	std::vector<int> indices;
	CHECK_THROWS_AS(lp.ToStandardForm(a, b, c, slacks, artificials, indices), LinearProgrammingError);
}

TEST_CASE("LP model: conversion clears reused auxiliary objective coefficients", "[LP][lp-model]")
{
	LinearProgram lp(1);
	lp.SetObjective({1});
	lp.AddConstraint({1}, LPConstraintType::LessEqual, 2);
	Matrix<Real> a;
	Vector<Real> b, c{9, 9};
	int slacks = 0, artificials = 0;
	std::vector<int> indices;
	lp.ToStandardForm(a, b, c, slacks, artificials, indices);
	REQUIRE(c.size() == 2);
	CHECK(c[0] == 1);
	CHECK(c[1] == 0);
}

TEST_CASE("LP model: raw tableau rejects nonfinite standard form", "[LP][lp-model]")
{
	const Matrix<Real> a(1, 2, {std::numeric_limits<Real>::quiet_NaN(), 1});
	SimplexTableau tableau;
	CHECK_THROWS_AS(tableau.Initialize(a, Vector<Real>{1}, Vector<Real>{1, 0}, {}, 1, 1, 0),
		LinearProgrammingError);
}

TEST_CASE("LP model: matrix convenience paths reject nonfinite data", "[LP][lp-model]")
{
	for (const Real invalid : {std::numeric_limits<Real>::quiet_NaN(),
		std::numeric_limits<Real>::infinity(), -std::numeric_limits<Real>::infinity()}) {
		for (int field = 0; field < 3; ++field) {
			CAPTURE(field, invalid);
			const Matrix<Real> a(1, 1, {field == 1 ? invalid : Real(1)});
			const Vector<Real> b{field == 2 ? invalid : Real(1)};
			const Vector<Real> c{field == 0 ? invalid : Real(1)};
			CHECK_THROWS_AS(SolveLP(c, a, b), LinearProgrammingError);
			CHECK_THROWS_AS(SolveLPDual(c, a, b), LinearProgrammingError);
			CHECK_THROWS_AS(SolveLPWithSensitivity(c, a, b), LinearProgrammingError);
		}
	}
}

TEST_CASE("LP model: dimension and enum failures stop before conversion writes", "[LP][lp-model]")
{
	CheckRejectedModel(LinearProgram());
	CheckRejectedModel(LinearProgram(1));
	LinearProgram lp(1);
	lp.SetObjective({1});
	lp.AddConstraint({1}, LPConstraintType::LessEqual, 2);
	lp.SetObjective({1, 2});
	CheckRejectedModel(lp);
	lp.SetObjective({1}, static_cast<LPObjective>(-1));
	CheckRejectedModel(lp);
	lp.SetObjective({1});
	lp.AddConstraint({1}, static_cast<LPConstraintType>(-1), 3);
	CheckRejectedModel(lp);
}

TEST_CASE("LP model: successful revalidation clears stale diagnostics", "[LP][lp-model]")
{
	LinearProgram lp(1);
	lp.AddConstraint({1}, LPConstraintType::LessEqual, 2);
	lp.SetObjective({1, 2});
	std::string error;
	CHECK_FALSE(lp.Validate(error));
	CHECK_FALSE(error.empty());
	lp.SetObjective({1});
	CHECK(lp.Validate(error));
	CHECK(error.empty());
}

TEST_CASE("LP model: incompatible standard form output aliases are rejected", "[LP][lp-model]")
{
	LinearProgram lp(1);
	lp.SetObjective({1});
	lp.AddConstraint({1}, LPConstraintType::LessEqual, 2);
	Matrix<Real> a(1, 1, {7});
	Vector<Real> b{8}, c{9};
	int slacks = 4, artificials = 5;
	std::vector<int> indices{6};
	CHECK_THROWS_AS(lp.ToStandardForm(a, b, b, slacks, artificials, indices), LinearProgrammingError);
	CHECK_THROWS_AS(lp.ToStandardForm(a, b, c, slacks, slacks, indices), LinearProgrammingError);
	REQUIRE(a.rows() == 1);
	REQUIRE(a.cols() == 1);
	CHECK(a(0, 0) == 7);
	CHECK(b[0] == 8);
	CHECK(c[0] == 9);
	CHECK(slacks == 4);
	CHECK(artificials == 5);
	CHECK(indices == std::vector<int>{6});
}

TEST_CASE("LP model: normalized mixed rows match a checked canonical tableau", "[LP][lp-model]")
{
	for (const auto sense : {LPObjective::Minimize, LPObjective::Maximize}) {
		LinearProgram lp(2);
		lp.SetObjective({2, -1}, sense);
		lp.AddConstraint({-1, 0}, LPConstraintType::LessEqual, -2);
		lp.AddConstraint({0, -1}, LPConstraintType::GreaterEqual, -3);
		lp.AddConstraint({1, -1}, LPConstraintType::Equal, 0);
		Matrix<Real> a;
		Vector<Real> b, c;
		int slacks = 0, artificials = 0;
		std::vector<int> indices;
		lp.ToStandardForm(a, b, c, slacks, artificials, indices);
		REQUIRE(a.rows() == 3);
		REQUIRE(a.cols() == 6);
		REQUIRE(b.size() == 3);
		REQUIRE(c.size() == 6);
		CHECK(slacks == 2);
		CHECK(artificials == 2);
		REQUIRE(indices == std::vector<int>{4, 5});
		const Matrix<Real> expected(3, 6, {1, 0, -1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, -1, 0, 0, 0, 1});
		for (int i = 0; i < a.rows(); ++i) {
			for (int j = 0; j < a.cols(); ++j)
				CHECK(a(i, j) == expected(i, j));
		}
		CHECK(b[0] == 2);
		CHECK(b[1] == 3);
		CHECK(b[2] == 0);
		const Real sign = sense == LPObjective::Minimize ? Real(1) : Real(-1);
		CHECK(c[0] == sign * 2);
		CHECK(c[1] == -sign);
		for (int j = 2; j < c.size(); ++j)
			CHECK(c[j] == 0);
		SimplexTableau tableau;
		tableau.Initialize(a, b, c, indices, 2, slacks, artificials);
		CHECK(tableau.Basis() == std::vector<int>{4, 3, 5});

		for (const auto& result : {SolveLP(lp), SolveLPTableau(lp), SolveLPWithSensitivity(lp)}) {
			REQUIRE(result.IsOptimal());
			REQUIRE(result.x.size() == 2);
			const Real optimum = sense == LPObjective::Minimize ? Real(2) : Real(3);
			CHECK(std::abs(result.x[0] - optimum) <= Precision::DefaultToleranceStrict);
			CHECK(std::abs(result.x[1] - optimum) <= Precision::DefaultToleranceStrict);
			CHECK(std::abs(result.objectiveValue - optimum) <= Precision::DefaultToleranceStrict);
			CHECK(result.achieved_tolerance <= Precision::DefaultToleranceStrict);
		}
	}
}

TEST_CASE("LP model: standard form allocation obeys the matrix element budget", "[LP][lp-model]")
{
	LinearProgram lp(1);
	lp.SetObjective({1});
	const int rows = static_cast<int>(std::sqrt(static_cast<Real>(Matrix<Real>::MaxElementCount()))) + 1;
	for (int i = 0; i < rows; ++i)
		lp.AddConstraint({1}, LPConstraintType::LessEqual, 1);
	Matrix<Real> a;
	Vector<Real> b, c;
	int slacks = 0, artificials = 0;
	std::vector<int> indices;
	CHECK_THROWS_AS(lp.ToStandardForm(a, b, c, slacks, artificials, indices), MatrixAllocationError);
	CHECK(a.rows() == 0);
	CHECK(a.cols() == 0);
}

TEST_CASE("LP model: raw tableau shape and index rejection preserves installed state", "[LP][lp-model]")
{
	const Matrix<Real> a(1, 2, {2, 1});
	const Vector<Real> b{3}, c{1, 0};
	SimplexTableau tableau;
	tableau.Initialize(a, b, c, {}, 1, 1, 0);
	tableau.Pivot(0, 0);
	const auto basis = tableau.Basis();
	const Real objective = tableau.GetObjectiveValue();
	auto reject = [&](const Matrix<Real>& matrix, const Vector<Real>& rhs, const Vector<Real>& costs,
		const std::vector<int>& indices, int originals, int slacks, int artificials) {
		CHECK_THROWS_AS(tableau.Initialize(matrix, rhs, costs, indices, originals, slacks, artificials),
			LinearProgrammingError);
		CHECK(tableau.NumRows() == 1);
		CHECK(tableau.NumCols() == 3);
		CHECK(tableau.Basis() == basis);
		CHECK(tableau.Iterations() == 1);
		CHECK(tableau.GetObjectiveValue() == objective);
		CHECK(tableau.GetOriginalSolution()[0] == Real(3) / 2);
	};
	reject(a, {}, c, {}, 1, 1, 0);
	reject(a, b, Vector<Real>{1}, {}, 1, 1, 0);
	reject(a, b, Vector<Real>{1, 0, 0}, {}, 1, 1, 0);
	reject(a, b, c, {}, 0, 2, 0);
	reject(a, b, c, {}, -1, 3, 0);
	reject(a, b, c, {}, 3, 0, 0);
	reject(a, b, c, {}, 1, -1, 2);
	reject(a, b, c, {}, 1, 2, -1);
	reject(a, b, c, {}, 1, 0, 0);
	reject(a, b, c, {}, std::numeric_limits<int>::max(), std::numeric_limits<int>::max(), 0);
	reject(a, b, c, {}, 1, 1, std::numeric_limits<int>::max());
	reject(a, b, c, {1}, 1, 1, 0);
	reject(a, b, c, {}, 1, 0, 1);
	reject(a, b, c, {-1}, 1, 0, 1);
	reject(a, b, c, {0}, 1, 0, 1);
	reject(a, b, c, {2}, 1, 0, 1);
	reject(Matrix<Real>(1, 3, {2, 1, 1}), b, Vector<Real>{1, 0, 0}, {1, 1}, 1, 0, 2);
	reject(Matrix<Real>(1, 2, {1, 0}), b, c, {}, 1, 1, 0);
	reject(Matrix<Real>(), {}, {}, {}, 1, 0, 0);
	reject(Matrix<Real>(std::numeric_limits<int>::max(), 0), {}, {}, {}, 1, 0, 0);
	reject(Matrix<Real>(0, std::numeric_limits<int>::max()), {}, {}, {}, 1, 0, 0);

	for (const Real invalid : {std::numeric_limits<Real>::quiet_NaN(),
		std::numeric_limits<Real>::infinity(), -std::numeric_limits<Real>::infinity()}) {
		reject(Matrix<Real>(1, 2, {invalid, 1}), b, c, {}, 1, 1, 0);
		reject(a, Vector<Real>{invalid}, c, {}, 1, 1, 0);
		reject(a, b, Vector<Real>{invalid, 0}, {}, 1, 1, 0);
	}
	tableau.Initialize(a, Vector<Real>{-3}, c, {}, 1, 1, 0);
	CHECK(tableau.Iterations() == 0);
	CHECK(tableau.Basis() == std::vector<int>{1});
	CHECK(tableau.GetRHS(0) == -3);
}

TEST_CASE("LP model: derived slack overflow cannot report optimal success", "[LP][lp-model]")
{
	const Real largest = std::numeric_limits<Real>::max();
	LinearProgram lp(1);
	lp.SetObjective({-1});
	lp.AddConstraint({-1}, LPConstraintType::LessEqual, largest);
	lp.AddConstraint({1}, LPConstraintType::LessEqual, largest);
	for (const auto& result : {SolveLP(lp), SolveLPWithSensitivity(lp)}) {
		CHECK_FALSE(result.IsOptimal());
		CHECK_FALSE(result.IsSuccess());
		CHECK_FALSE(result.converged);
		CHECK(result.status == AlgorithmStatus::NumericalInstability);
		CHECK(std::isnan(result.achieved_tolerance));
		CHECK_FALSE(result.error_message.empty());
	}
}

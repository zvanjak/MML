#include <catch2/catch_test_macros.hpp>
#include <mml/base/ODESystemSolution.h>
#include <mml/base/DAESystem.h>
#include <limits>

using namespace MML;

TEST_CASE("Matrix permits long narrow shapes above the former axis cap", "[matrix-element-limits]")
{
	Matrix<Real> wide(1, 100001);
	wide(0, 100000) = 7;
	const auto tall = wide.Transposed();
	REQUIRE(tall.rows() == 100001);
	REQUIRE(tall.cols() == 1);
	REQUIRE(tall(100000, 0) == Real(7));
	wide.resize(1, 150000, true);
	REQUIRE(wide(0, 100000) == Real(7));
	REQUIRE(wide.cols() == 150000);
}

TEST_CASE("Matrix enforces its product budget before allocating", "[matrix-element-limits]")
{
	REQUIRE(Matrix<Real>::MaxElementCount() == size_t{100000000});
	const int over = static_cast<int>(Matrix<Real>::MaxElementCount()) + 1;
	REQUIRE_THROWS_AS(Matrix<Real>(1, over), MatrixAllocationError);
	REQUIRE_THROWS_AS(Matrix<Real>(over, 1), MatrixAllocationError);
	REQUIRE_THROWS_AS(Matrix<Real>(std::numeric_limits<int>::max(), std::numeric_limits<int>::max()),
		MatrixAllocationError);
	Matrix<Real> original(1, 2, Real(7));
	REQUIRE_THROWS_AS(original.resize(1, over, true), MatrixAllocationError);
	REQUIRE(original.rows() == 1);
	REQUIRE(original.cols() == 2);
	REQUIRE(original(0, 1) == Real(7));
	REQUIRE_THROWS_AS(Matrix<Real>(-1, 2), MatrixDimensionError);
}

TEST_CASE("Empty Matrix axes are bounded by index representation rather than memory", "[matrix-element-limits]")
{
	const int limit = std::numeric_limits<int>::max();
	Matrix<Real> noColumns(limit, 0);
	Matrix<Real> noRows(0, limit);
	REQUIRE(noColumns.rows() == limit);
	REQUIRE(noColumns.cols() == 0);
	REQUIRE(noRows.rows() == 0);
	REQUIRE(noRows.cols() == limit);
	REQUIRE(Matrix<Real>::MaxDimension() == limit);
}

TEST_CASE("ODE trajectory growth passes the former column cap", "[trajectory-storage][ODE]")
{
	ODESystemSolution solution(0, 1, 1, 0);
	solution.fillValues(0, 0, Vector<Real>{7});
	REQUIRE_NOTHROW(solution.fillValues(99999, 1, Vector<Real>{9}));
	REQUIRE(solution.size() == 100000);
	REQUIRE(solution.capacity() >= solution.size());
	REQUIRE(solution.getXValue(0, 0) == Real(7));
	REQUIRE(solution.getXValue(99999, 0) == Real(9));
	REQUIRE_NOTHROW(solution.fillValues(100000, 2, Vector<Real>{11}));
	REQUIRE(solution.size() == 100001);
	const int capacity = solution.capacity();
	REQUIRE_THROWS_AS(solution.reserveForPoint(solution.max_size()), ResourceLimitError);
	REQUIRE(solution.size() == 100001);
	REQUIRE(solution.capacity() == capacity);
	REQUIRE(solution.getTValue(99999) == Real(1));
}

TEST_CASE("DAE trajectory growth passes the former column cap", "[trajectory-storage][DAE]")
{
	DAESolution solution(0, 1, 1, 1, 0);
	solution.fillValues(0, 0, Vector<Real>{7}, Vector<Real>{8});
	REQUIRE_NOTHROW(solution.fillValues(99999, 1, Vector<Real>{9}, Vector<Real>{10}));
	REQUIRE(solution.size() == 100000);
	REQUIRE(solution.capacity() >= solution.size());
	REQUIRE(solution.getXValue(0, 0) == Real(7));
	REQUIRE(solution.getYValue(99999, 0) == Real(10));
	REQUIRE_NOTHROW(solution.fillValues(100000, 2, Vector<Real>{11}, Vector<Real>{12}));
	REQUIRE(solution.size() == 100001);
	const int capacity = solution.capacity();
	REQUIRE_THROWS_AS(solution.reserveForPoint(solution.max_size()), ResourceLimitError);
	REQUIRE(solution.size() == 100001);
	REQUIRE(solution.capacity() == capacity);
	REQUIRE(solution.getTValue(99999) == Real(1));
}

TEST_CASE("Failed ODE reserve preserves capacity and existing points", "[trajectory-storage][ODE]")
{
	ODESystemSolution solution(0, 1, 1, 1);
	solution.fillValues(0, 0, Vector<Real>{7});
	const int capacity = solution.capacity();
	REQUIRE_THROWS_AS(solution.reserve(solution.max_size() + 1), ResourceLimitError);
	REQUIRE(solution.capacity() == capacity);
	REQUIRE(solution.size() == 1);
	REQUIRE(solution.getXValue(0, 0) == Real(7));
}

TEST_CASE("Trajectory storage rejects oversized indices and final sizes transactionally", "[trajectory-storage]")
{
	ODESystemSolution ode(0, 1, 1, 0);
	DAESolution dae(0, 1, 1, 1, 0);
	ode.fillValues(0, 0, Vector<Real>{7});
	dae.fillValues(0, 0, Vector<Real>{7}, Vector<Real>{8});
	constexpr int huge = std::numeric_limits<int>::max();
	REQUIRE_THROWS_AS(ode.setTVal(huge, 1), ResourceLimitError);
	REQUIRE_THROWS_AS(ode.setXVal(huge, 0, 1), ResourceLimitError);
	REQUIRE_THROWS_AS(dae.setTVal(huge, 1), ResourceLimitError);
	REQUIRE_THROWS_AS(dae.setXVal(huge, 0, 1), ResourceLimitError);
	REQUIRE_THROWS_AS(dae.setYVal(huge, 0, 1), ResourceLimitError);
	REQUIRE_THROWS_AS(ode.setFinalSize(huge), ResourceLimitError);
	REQUIRE_THROWS_AS(dae.setFinalSize(huge), ResourceLimitError);
	REQUIRE(ode.size() == 1);
	REQUIRE(ode.capacity() == 1);
	REQUIRE(dae.size() == 1);
	REQUIRE(dae.capacity() == 1);
	REQUIRE(ode.getXValue(0, 0) == Real(7));
	REQUIRE(dae.getYValue(0, 0) == Real(8));
	REQUIRE_THROWS_AS(ODESystemSolution(0, 1, 1, huge), ResourceLimitError);
	REQUIRE_THROWS_AS(DAESolution(0, 1, 1, 1, huge), ResourceLimitError);
}

TEST_CASE("Trajectory point limits honor both DAE matrices and element budgets", "[trajectory-storage]")
{
	ODESystemSolution ode(0, 1, 2000, 0);
	DAESolution dae(0, 1, 1, 2000, 0);
	const int limit = static_cast<int>(Matrix<Real>::MaxElementCount() / 2000);
	REQUIRE(ode.max_size() == limit);
	REQUIRE(dae.max_size() == limit);
	REQUIRE_THROWS_AS(ode.reserveForPoint(limit), ResourceLimitError);
	REQUIRE_THROWS_AS(dae.reserveForPoint(limit), ResourceLimitError);
	REQUIRE(ode.capacity() == 1);
	REQUIRE(dae.capacity() == 1);
	REQUIRE(detail::TrajectoryCapacity(40000, limit, limit, true) == limit);
	REQUIRE(detail::TrajectoryMaxPoints({1}) == static_cast<int>(Matrix<Real>::MaxElementCount()));
	REQUIRE(detail::TrajectoryMaxPoints({0, 0}) == static_cast<int>(Matrix<Real>::MaxElementCount()));
	REQUIRE(detail::TrajectoryCapacity(80000000, 100000000, 100000000, true) == 100000000);
	const int indexLimit = std::numeric_limits<int>::max();
	REQUIRE(detail::TrajectoryCapacity(indexLimit - 1, indexLimit, indexLimit, true) == indexLimit);
	DAESolution noAlgebraicVariables(0, 1, 1, 0, 0);
	noAlgebraicVariables.fillValues(0, 0, Vector<Real>{3}, Vector<Real>{});
	noAlgebraicVariables.fillValues(1, 1, Vector<Real>{4}, Vector<Real>{});
	REQUIRE(noAlgebraicVariables.getXValue(1, 0) == Real(4));
}

TEST_CASE("Resource limits are explicit non-usable algorithm results", "[trajectory-storage]")
{
	REQUIRE(to_string(AlgorithmStatus::ResourceLimitExceeded) == "ResourceLimitExceeded");
	REQUIRE_FALSE(IsUsable(AlgorithmStatus::ResourceLimitExceeded));
}

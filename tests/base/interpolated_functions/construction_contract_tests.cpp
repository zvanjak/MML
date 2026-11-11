#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "../../TestPrecision.h"

#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/base/InterpolatedFunction.h>
#include <mml/base/RichardsonExtrapolation.h>
#include <mml/base/DAESystem.h>
#include <mml/algorithms/DAESolvers/DAEBackwardEuler.h>
#endif

#include <array>
#include <climits>
#include <cmath>
#include <limits>
#include <type_traits>

using namespace MML;

namespace {
	void CheckConstructionValue(Real actual, Real expected) {
		CHECK_THAT(actual, Catch::Matchers::WithinAbs(
			expected, Precision::EqualityTolerance * std::max(Real(1), std::abs(expected))));
	}

	std::array<int, 4> jacobianCalls{};
	void ContractDifferential(Real, const Vector<Real>& x, const Vector<Real>& y, Vector<Real>& output) {
		output[0] = -2 * x[0] + 3 * y[0];
	}
	void ContractAlgebraic(Real, const Vector<Real>& x, const Vector<Real>& y, Vector<Real>& output) {
		output[0] = 5 * x[0] + 7 * y[0] - 7;
	}
	void ContractFX(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& output) {
		++jacobianCalls[0];
		output(0, 0) = -2;
	}
	void ContractFY(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& output) {
		++jacobianCalls[1];
		output(0, 0) = 3;
	}
	void ContractGX(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& output) {
		++jacobianCalls[2];
		output(0, 0) = 5;
	}
	void ContractGY(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& output) {
		++jacobianCalls[3];
		output(0, 0) = 7;
	}

	template<typename Adapter>
	void CheckDAESolverConsumption(Adapter& adapter) {
		if constexpr (std::is_convertible_v<Adapter*, IODESystemDAEWithJacobian*>) {
			IODESystemDAEWithJacobian& analytic = adapter;
			IODESystemDAE& ordinary = adapter;
			CHECK(ordinary.getDiffDim() == 1);
			CHECK(ordinary.getAlgDim() == 1);
			CHECK(ordinary.getTotalDim() == 2);
			const Vector<Real> x{0}, y{1};
			Vector<Real> output(1);
			ordinary.diffEqs(0, x, y, output);
			CheckConstructionValue(output[0], 3);
			ordinary.algConstraints(0, x, y, output);
			CheckConstructionValue(output[0], 0);
			Matrix<Real> fx(1, 1), fy(1, 1), gx(1, 1), gy(1, 1);
			analytic.allJacobians(0, x, y, fx, fy, gx, gy);
			CheckConstructionValue(fx(0, 0), -2);
			CheckConstructionValue(fy(0, 0), 3);
			CheckConstructionValue(gx(0, 0), 5);
			CheckConstructionValue(gy(0, 0), 7);

			jacobianCalls.fill(0);
			DAESolverConfig config;
			config.step_size = Real(0.125);
			const auto result = SolveDAEBackwardEuler(adapter, 0, x, y, Real(0.5), config);
			REQUIRE(result.status == AlgorithmStatus::Success);
			CHECK(result.accepted_steps == 4);
			// Independent recurrence for x' = 3 - (29/7)x, using four fixed BE steps.
			const Real expectedX = Real(21) / 29 * (1 - std::pow(1 + Real(29) / 7 * config.step_size, MML::Real(-4)));
			CheckConstructionValue(result.solution.getXValuesAtEnd()[0], expectedX);
			CheckConstructionValue(result.solution.getYValuesAtEnd()[0], 1 - Real(5) / 7 * expectedX);
			for (int count : jacobianCalls)
				CHECK(count > 0);
		} else {
			FAIL("DAESystemWithJacobian cannot bind to the analytic DAE solver interface");
		}
	}
}

TEST_CASE("Rational barycentric construction rejects incomplete node data",
	"[interpolation][construction-contract][safe-constructor-witness]") {
	const Vector<Real> x{0, 1, 2};
	CHECK_THROWS_AS(BarycentricRationalInterp(x, Vector<Real>{2}, 1), RealFuncInterpInitError);
	CHECK_THROWS_AS(BarycentricRationalInterp(x, Vector<Real>{2, 3, 4, 5}, 1), RealFuncInterpInitError);
	CHECK_THROWS_AS(BarycentricRationalInterp(Vector<Real>{}, Vector<Real>{}, 0), RealFuncInterpInitError);
	for (const Vector<Real>& nodes : {Vector<Real>{0, 0, 1}, Vector<Real>{0, 2, 1}}) {
		CHECK_THROWS_AS(BarycentricRationalInterp(nodes, Vector<Real>{1, 2, 3}, 1), RealFuncInterpInitError);
	}
	for (Real bad : {std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity()})
		for (int index = 0; index < 3; ++index) {
			Vector<Real> nodes = x;
			nodes[index] = bad;
			CHECK_THROWS_AS(BarycentricRationalInterp(nodes, Vector<Real>{1, 2, 3}, 1), RealFuncInterpInitError);
		}
}

TEST_CASE("Rational minimum data preserves constant and affine interpolation", "[interpolation][construction-contract]") {
	BarycentricRationalInterp constant(Vector<Real>{2}, Vector<Real>{7}, 0);
	for (Real query : {Real(-3), Real(2), Real(9), std::numeric_limits<Real>::max()}) {
		CheckConstructionValue(constant(query), 7);
		const auto detailed = constant.EvaluateDetailed(query);
		CheckConstructionValue(detailed.value, 7);
		CHECK(detailed.interval_index == -1);
	}
	constant.setExtrapolationPolicy(ExtrapolationPolicy::Throw);
	CHECK_THROWS_AS(constant(3), RealFuncInterpRuntimeError);
	constant.setExtrapolationPolicy(ExtrapolationPolicy::Clamp);
	CheckConstructionValue(constant(3), 7);
	InterpolationConfig config;
	config.extrapolation_policy = ExtrapolationPolicy::Clamp;
	const auto detailed = constant.EvaluateDetailed(3, config);
	CheckConstructionValue(detailed.value, 7);
	CHECK(detailed.WasClamped());
	CHECK(detailed.interval_index == -1);
	for (bool descending : {false, true}) {
		const Vector<Real> x = descending ? Vector<Real>{3, 1} : Vector<Real>{1, 3};
		const Vector<Real> y = descending ? Vector<Real>{11, 5} : Vector<Real>{5, 11};
		BarycentricRationalInterp affine(x, y, 1);
		for (Real query : {Real(0), Real(1), Real(2), Real(3), Real(4)})
			CheckConstructionValue(affine(query), 3 * query + 2);
	}
}

TEST_CASE("Shared scalar interpolators reject nonfinite nodes", "[interpolation][construction-contract]") {
	for (Real bad : {std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity(),
			 -std::numeric_limits<Real>::infinity()})
		for (int index = 0; index < 3; ++index) {
			Vector<Real> x{0, 1, 2};
			x[index] = bad;
			const Vector<Real> y{1, 2, 3};
			CHECK_THROWS_AS(LinearInterpRealFunc(x, y), RealFuncInterpInitError);
			CHECK_THROWS_AS(SplineInterpRealFunc(x, y), RealFuncInterpInitError);
		}
}

TEST_CASE("Scalar splines validate small datasets before initializing derivatives",
	"[interpolation][spline][construction-contract]") {
	CHECK_THROWS_AS(SplineInterpRealFunc(Vector<Real>{}, Vector<Real>{}), RealFuncInterpInitError);
	CHECK_THROWS_AS(SplineInterpRealFunc(Vector<Real>{0}, Vector<Real>{1}), RealFuncInterpInitError);
	CHECK_THROWS_AS(SplineInterpRealFunc(Vector<Real>{0, 1}, Vector<Real>{1}), RealFuncInterpInitError);
	CHECK_THROWS_AS(SplineInterpRealFunc(Vector<Real>{0, 1}, Vector<Real>{1, 2, 3}), RealFuncInterpInitError);
	CHECK_THROWS_AS(SplineInterpRealFunc(Vector<Real>{0, 0}, Vector<Real>{1, 2}), RealFuncInterpInitError);
}

TEST_CASE("Two dimensional grid constructors validate complete axes",
	"[interpolation][construction-contract][safe-constructor-witness]") {
	const Vector<Real> axis{0, 1};
	for (const Vector<Real>& shortAxis : {Vector<Real>{}, Vector<Real>{0}}) {
		const Matrix<Real> rows(shortAxis.size(), 2), columns(2, shortAxis.size());
		CHECK_THROWS_AS(BilinearInterp2D(shortAxis, axis, rows), RealFuncInterpInitError);
		CHECK_THROWS_AS(BilinearInterp2D(axis, shortAxis, columns), RealFuncInterpInitError);
		CHECK_THROWS_AS(BicubicSplineInterp2D(shortAxis, axis, rows), RealFuncInterpInitError);
		CHECK_THROWS_AS(BicubicSplineInterp2D(axis, shortAxis, columns), RealFuncInterpInitError);
	}
	const Matrix<Real> values(3, 2);
	for (const Vector<Real>& nodes : {Vector<Real>{0, 0, 1}, Vector<Real>{0, 2, 1},
			 Vector<Real>{0, std::numeric_limits<Real>::quiet_NaN(), 2},
			 Vector<Real>{0, 1, std::numeric_limits<Real>::infinity()}}) {
		CHECK_THROWS_AS(BilinearInterp2D(nodes, axis, values), RealFuncInterpInitError);
		CHECK_THROWS_AS(BicubicSplineInterp2D(nodes, axis, values), RealFuncInterpInitError);
		const Matrix<Real> transposed(2, 3);
		CHECK_THROWS_AS(BilinearInterp2D(axis, nodes, transposed), RealFuncInterpInitError);
		CHECK_THROWS_AS(BicubicSplineInterp2D(axis, nodes, transposed), RealFuncInterpInitError);
	}
	CHECK_THROWS_AS(BilinearInterp2D(axis, axis, Matrix<Real>(2, 3)), RealFuncInterpInitError);
	CHECK_THROWS_AS(BicubicSplineInterp2D(axis, axis, Matrix<Real>(3, 2)), RealFuncInterpInitError);
}

TEST_CASE("Minimum bilinear and spline grids reproduce affine data in both directions", "[interpolation][construction-contract]") {
	for (bool reverseX : {false, true})
		for (bool reverseY : {false, true}) {
			const Vector<Real> x = reverseX ? Vector<Real>{3, 1} : Vector<Real>{1, 3};
			const Vector<Real> y = reverseY ? Vector<Real>{4, -2} : Vector<Real>{-2, 4};
			Matrix<Real> z(2, 2);
			for (int i = 0; i < 2; ++i)
				for (int j = 0; j < 2; ++j)
					z(i, j) = 2 * x[i] - 3 * y[j] + 4;
			const BilinearInterp2D bilinear(x, y, z);
			const BicubicSplineInterp2D spline(x, y, z);
			for (Real queryX : {Real(0), Real(1), Real(2), Real(3), Real(4)})
				for (Real queryY : {Real(-3), Real(-2), Real(1), Real(4), Real(5)}) {
					const Real expected = 2 * queryX - 3 * queryY + 4;
					CheckConstructionValue(bilinear(queryX, queryY), expected);
					CheckConstructionValue(spline(queryX, queryY), expected);
					Real value, dx, dy;
					bilinear.interpWithDerivatives(queryX, queryY, value, dx, dy);
					CheckConstructionValue(value, expected);
					CheckConstructionValue(dx, 2);
					CheckConstructionValue(dy, -3);
				}
		}
}

TEST_CASE("Parametric curves reject incomplete or nonfinite point matrices before indexing", "[interpolation][construction-contract]") {
	for (int rows : {0, 1, 2})
		for (int cols : {0, 1, 2, 3})
			if (rows < 2 || cols != 2) {
				const Matrix<Real> points(rows, cols);
				for (bool closed : {false, true}) {
					if (cols != 2) {
						CHECK_THROWS_AS(LinInterpParametricCurve<2>(points, closed), MatrixDimensionError);
						CHECK_THROWS_AS(SplineInterpParametricCurve<2>(points, closed), MatrixDimensionError);
					} else {
						CHECK_THROWS_AS(LinInterpParametricCurve<2>(points, closed), RealFuncInterpInitError);
						CHECK_THROWS_AS(SplineInterpParametricCurve<2>(points, closed), RealFuncInterpInitError);
					}
				}
			}
	for (Real bad : {std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity()})
		for (int row = 0; row < 2; ++row)
			for (int col = 0; col < 2; ++col) {
				Matrix<Real> points(2, 2, {0, 0, 3, 4});
				points(row, col) = bad;
				CHECK_THROWS_AS(LinInterpParametricCurve<2>(points), RealFuncInterpInitError);
				CHECK_THROWS_AS(SplineInterpParametricCurve<2>(points), RealFuncInterpInitError);
			}
	const Matrix<Real> points(2, 2, {0, 0, 3, 4});
	for (const auto bounds : {std::array<Real, 2>{0, 0}, std::array<Real, 2>{2, 1},
			 std::array<Real, 2>{std::numeric_limits<Real>::quiet_NaN(), 1},
			 std::array<Real, 2>{0, std::numeric_limits<Real>::infinity()},
			 std::array<Real, 2>{-std::numeric_limits<Real>::max(), std::numeric_limits<Real>::max()}})
		CHECK_THROWS_AS(SplineInterpParametricCurve<2>(bounds[0], bounds[1], points), RealFuncInterpInitError);
}

TEST_CASE("Minimum curves preserve endpoints and chord length interpolation", "[interpolation][construction-contract]") {
	const Matrix<Real> points(2, 2, {1, 2, 4, 6});
	const LinInterpParametricCurve<2> linear(points);
	const SplineInterpParametricCurve<2> spline(2, 6, points);
	for (Real t : {Real(0), Real(0.25), Real(0.5), Real(1)}) {
		const auto a = linear(t);
		const auto b = spline(2 + 4 * t);
		CheckConstructionValue(a[0], 1 + 3 * t);
		CheckConstructionValue(a[1], 2 + 4 * t);
		CheckConstructionValue(b[0], 1 + 3 * t);
		CheckConstructionValue(b[1], 2 + 4 * t);
	}
	const Matrix<Real> polyline(3, 2, {0, 0, 3, 4, 3, 14});
	const LinInterpParametricCurve<2> chord(polyline);
	CheckConstructionValue(chord(Real(1) / 3)[0], 3);
	CheckConstructionValue(chord(Real(1) / 3)[1], 4);
	const Matrix<Real> repeated(2, 2, {3, 4, 3, 4});
	const LinInterpParametricCurve<2> constant(repeated);
	CheckConstructionValue(constant(Real(0.5))[0], 3);
	CheckConstructionValue(constant(Real(0.5))[1], 4);
	const LinInterpParametricCurve<2> closedLinear(points, true);
	const SplineInterpParametricCurve<2> closedSpline(points, true);
	CheckConstructionValue(closedLinear(0)[0], 1);
	CheckConstructionValue(closedLinear(0)[1], 2);
	CheckConstructionValue(closedSpline(0)[0], 1);
	CheckConstructionValue(closedSpline(0)[1], 2);
}

TEST_CASE("Curve chord lengths avoid squared norm overflow and underflow", "[interpolation][construction-contract]") {
	const int exponent = 3 * std::numeric_limits<Real>::max_exponent / 4;
	for (Real scale : {std::ldexp(Real(1), -exponent), std::ldexp(Real(1), exponent)}) {
		const Matrix<Real> points(2, 2, {0, 0, 3 * scale, 4 * scale});
		const LinInterpParametricCurve<2> linear(points);
		const SplineInterpParametricCurve<2> spline(points);
		for (const auto result : {linear(Real(0.5)), spline(Real(0.5))}) {
			CheckConstructionValue(result[0] / scale, Real(1.5));
			CheckConstructionValue(result[1] / scale, 2);
		}
	}
	const Real largest = std::numeric_limits<Real>::max();
	const Matrix<Real> unrepresentable(2, 2, {largest, 0, -largest, 0});
	CHECK_THROWS_AS(LinInterpParametricCurve<2>(unrepresentable), RealFuncInterpInitError);
	CHECK_THROWS_AS(SplineInterpParametricCurve<2>(unrepresentable), RealFuncInterpInitError);
}

TEST_CASE("Richardson rejects invalid configuration before evaluating", "[richardson][construction-contract]") {
	int calls = 0;
	const auto evaluator = [&](Real) { ++calls; return Real(1); };
	for (int budget : {0, -1, INT_MIN})
		CHECK_THROWS_AS(Richardson::Extrapolate(evaluator, 1, 2, budget), ArgumentError);
	for (Real step : {Real(0), -Real(0), std::numeric_limits<Real>::quiet_NaN(),
			 std::numeric_limits<Real>::infinity(), -std::numeric_limits<Real>::infinity()})
		CHECK_THROWS_AS(Richardson::Extrapolate(evaluator, step), ArgumentError);
	for (Real factor : {Real(-1), Real(0), Real(1), std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity()}) {
		CHECK_THROWS_AS(Richardson::Extrapolate(evaluator, 1, factor), ArgumentError);
		CHECK_THROWS_AS(Richardson::Extrapolate(evaluator, 1, 2, 3, factor), ArgumentError);
	}
	CHECK(calls == 0);
}

TEST_CASE("Richardson minimum budgets and quadratic error elimination retain meaningful results", "[richardson][construction-contract]") {
	for (Real step : {Real(-0.5), Real(0.5), std::numeric_limits<Real>::denorm_min()}) {
		int calls = 0;
		const auto result = Richardson::Extrapolate([&](Real h) {
			CHECK(h == step);
			++calls;
			return Real(7);
		}, step, 2, 1);
		CheckConstructionValue(result.value, 7);
		CHECK(result.iterations == 1);
		CHECK_FALSE(result.converged);
		CHECK(calls == 1);
	}
	for (Real factor : {Real(1.4), Real(2), Real(3)}) {
		const auto result = Richardson::Extrapolate([](Real h) { return 3 + 5 * h * h + 7 * h * h * h * h; },
			Real(0.5), factor, 5);
		CHECK(result.converged);
		CheckConstructionValue(result.value, 3);
	}
	const Real large = std::numeric_limits<Real>::max() / 2;
	const auto result = Richardson::Extrapolate([=](Real) { return large; }, 1, 2, 3);
	CHECK(result.value == large);
	CHECK(result.converged);
	const auto largeRatio = Richardson::Extrapolate([](Real) { return Real(7); }, 1, large, 2);
	CheckConstructionValue(largeRatio.value, 7);
	CHECK(largeRatio.converged);
}

TEST_CASE("Richardson refuses nonfinite estimates and unrepresentable refinement steps", "[richardson][construction-contract]") {
	for (Real bad : {std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity()})
		CHECK_THROWS_AS(Richardson::Extrapolate([=](Real) { return bad; }, 1), NumericInputError);
	int calls = 0;
	CHECK_THROWS_AS(Richardson::Extrapolate([&](Real h) {
		CHECK(h != 0);
		++calls;
		return Real(1);
	}, std::numeric_limits<Real>::denorm_min(), 2, 2), NumericalMethodError);
	CHECK(calls == 1);
	calls = 0;
	CHECK_THROWS_AS(Richardson::Extrapolate([&](Real) {
		++calls;
		return Real(1);
	}, std::numeric_limits<Real>::denorm_min(), std::nextafter(Real(1), Real(2)), 2), NumericalMethodError);
	CHECK(calls == 1);
	calls = 0;
	CHECK_THROWS_AS(Richardson::Extrapolate([&](Real) {
		return ++calls == 1 ? Real(1) : std::numeric_limits<Real>::infinity();
	}, 1), NumericInputError);
	CHECK(calls == 2);
}

TEST_CASE("DAE Jacobian adapter binds both interfaces and drives an implicit solver", "[dae][construction-contract]") {
	DAESystemWithJacobian adapter(1, 1, ContractDifferential, ContractAlgebraic,
		ContractFX, ContractFY, ContractGX, ContractGY);
	CheckDAESolverConsumption(adapter);
}

TEST_CASE("DAE Jacobian adapter rejects missing callbacks and default evaluation", "[dae][construction-contract]") {
	using Jacobian = void (*)(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>&);
	for (int missing = 0; missing < 4; ++missing) {
		std::array<Jacobian, 4> callbacks{ContractFX, ContractFY, ContractGX, ContractGY};
		callbacks[missing] = nullptr;
		CHECK_THROWS_AS(DAESystemWithJacobian(1, 1, ContractDifferential, ContractAlgebraic,
			callbacks[0], callbacks[1], callbacks[2], callbacks[3]), ArgumentError);
	}
	CHECK_THROWS_AS(DAESystemWithJacobian(1, 1, nullptr, ContractAlgebraic, ContractFX, ContractFY, ContractGX, ContractGY),
		ArgumentError);
	CHECK_THROWS_AS(DAESystemWithJacobian(1, 1, ContractDifferential, nullptr, ContractFX, ContractFY, ContractGX, ContractGY),
		ArgumentError);
	const DAESystemWithJacobian empty;
	const Vector<Real> x{0}, y{1};
	Vector<Real> vectorOutput(1);
	Matrix<Real> output(1, 1, {42});
	CHECK_THROWS_AS(empty.diffEqs(0, x, y, vectorOutput), NotImplementedError);
	CHECK_THROWS_AS(empty.algConstraints(0, x, y, vectorOutput), NotImplementedError);
	CHECK_THROWS_AS(empty.jacobian_fx(0, x, y, output), NotImplementedError);
	CHECK_THROWS_AS(empty.jacobian_fy(0, x, y, output), NotImplementedError);
	CHECK_THROWS_AS(empty.jacobian_gx(0, x, y, output), NotImplementedError);
	CHECK_THROWS_AS(empty.jacobian_gy(0, x, y, output), NotImplementedError);
	CHECK(output(0, 0) == 42);
}

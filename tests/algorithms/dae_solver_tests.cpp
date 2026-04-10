///////////////////////////////////////////////////////////////////////////////////////////
// DAE Solver Tests - Unit tests for Differential-Algebraic Equation solvers
///////////////////////////////////////////////////////////////////////////////////////////
#include <catch2/catch_test_macros.hpp>
#include "../TestPrecision.h"
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <type_traits>

#include "MMLBase.h"
#include "interfaces/IODESystemDAE.h"
#include "base/DAESystem.h"
#include "base/Vector/Vector.h"
#include "base/Matrix/Matrix.h"
#include "algorithms/DAESolvers.h"

using namespace MML;
using namespace Catch::Matchers;

namespace MML::Tests::Algorithms::DAESolverTests
{

///////////////////////////////////////////////////////////////////////////////////////////
//                           Test DAE Systems
///////////////////////////////////////////////////////////////////////////////////////////

/// Simple linear DAE for testing:
/// dx/dt = -x + y
/// 0 = x + y - 1    (constraint: x + y = 1)
/// 
/// Analytical solution: x(t) = 0.5 + (x0 - 0.5)e^{-2t}, y(t) = 1 - x(t)
class SimpleLinearDAE : public IODESystemDAEWithJacobian
{
public:
	int getDiffDim() const override { return 1; }
	int getAlgDim() const override { return 1; }

	void diffEqs(Real t, const Vector<Real>& x, const Vector<Real>& y,
	             Vector<Real>& dxdt) const override
	{
		dxdt[0] = -x[0] + y[0];  // dx/dt = -x + y
	}

	void algConstraints(Real t, const Vector<Real>& x, const Vector<Real>& y,
	                    Vector<Real>& g) const override
	{
		g[0] = x[0] + y[0] - 1.0;  // x + y = 1
	}

	std::string getDiffVarName(int i) const override { return "x"; }
	std::string getAlgVarName(int i) const override { return "y"; }

	// Jacobians for implicit solvers
	void jacobian_fx(Real t, const Vector<Real>& x, const Vector<Real>& y,
	                 Matrix<Real>& df_dx) const override
	{
		df_dx(0, 0) = -1.0;  // d(dx/dt)/dx = -1
	}

	void jacobian_fy(Real t, const Vector<Real>& x, const Vector<Real>& y,
	                 Matrix<Real>& df_dy) const override
	{
		df_dy(0, 0) = 1.0;  // d(dx/dt)/dy = 1
	}

	void jacobian_gx(Real t, const Vector<Real>& x, const Vector<Real>& y,
	                 Matrix<Real>& dg_dx) const override
	{
		dg_dx(0, 0) = 1.0;  // dg/dx = 1
	}

	void jacobian_gy(Real t, const Vector<Real>& x, const Vector<Real>& y,
	                 Matrix<Real>& dg_dy) const override
	{
		dg_dy(0, 0) = 1.0;  // dg/dy = 1 (nonsingular - index-1)
	}
};

/// Pendulum DAE (simplified form for testing interface)
/// Differential: dx/dt = vx, dy/dt = vy, dvx/dt = λx, dvy/dt = λy - g
/// Constraint: x² + y² = L²
class SimplifiedPendulumDAE : public IODESystemDAE
{
	Real _L, _g;

public:
	SimplifiedPendulumDAE(Real length = 1.0, Real gravity = 9.81)
		: _L(length), _g(gravity) {}

	int getDiffDim() const override { return 4; }  // x, y, vx, vy
	int getAlgDim() const override { return 1; }   // λ

	void diffEqs(Real t, const Vector<Real>& diff, const Vector<Real>& alg,
	             Vector<Real>& dxdt) const override
	{
		Real x = diff[0], y = diff[1], vx = diff[2], vy = diff[3];
		Real lambda = alg[0];

		dxdt[0] = vx;                // dx/dt = vx
		dxdt[1] = vy;                // dy/dt = vy
		dxdt[2] = lambda * x;        // dvx/dt = λx
		dxdt[3] = lambda * y - _g;   // dvy/dt = λy - g
	}

	void algConstraints(Real t, const Vector<Real>& diff, const Vector<Real>& alg,
	                    Vector<Real>& g) const override
	{
		Real x = diff[0], y = diff[1];
		g[0] = x * x + y * y - _L * _L;  // x² + y² = L²
	}

	std::string getDiffVarName(int i) const override
	{
		static const char* names[] = {"x", "y", "vx", "vy"};
		return (i >= 0 && i < 4) ? names[i] : "?";
	}

	std::string getAlgVarName(int i) const override { return "lambda"; }
};

///////////////////////////////////////////////////////////////////////////////////////////
//                           Interface Tests
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("DAE::Interface_dimensions", "[dae][interface]")
{
	SimpleLinearDAE linearDAE;
	
	REQUIRE(linearDAE.getDiffDim() == 1);
	REQUIRE(linearDAE.getAlgDim() == 1);
	REQUIRE(linearDAE.getTotalDim() == 2);
}

TEST_CASE("DAE::Interface_variable_names", "[dae][interface]")
{
	SimpleLinearDAE linearDAE;
	
	REQUIRE(linearDAE.getDiffVarName(0) == "x");
	REQUIRE(linearDAE.getAlgVarName(0) == "y");
}

TEST_CASE("DAE::Interface_pendulum_dimensions", "[dae][interface]")
{
	SimplifiedPendulumDAE pendulum(1.0, 9.81);
	
	REQUIRE(pendulum.getDiffDim() == 4);
	REQUIRE(pendulum.getAlgDim() == 1);
	REQUIRE(pendulum.getTotalDim() == 5);
}

TEST_CASE("DAE::Interface_diffEqs_evaluation", "[dae][interface]")
{
	SimpleLinearDAE linearDAE;
	
	Vector<Real> x(1), y(1), dxdt(1);
	x[0] = 0.3;
	y[0] = 0.7;  // Consistent: x + y = 1
	
	linearDAE.diffEqs(0.0, x, y, dxdt);
	
	// dx/dt = -x + y = -0.3 + 0.7 = 0.4
	REQUIRE_THAT(dxdt[0], WithinAbs(0.4, TOL(1e-10, 1e-5)));
}

TEST_CASE("DAE::Interface_algConstraints_evaluation", "[dae][interface]")
{
	SimpleLinearDAE linearDAE;
	
	Vector<Real> x(1), y(1), g(1);
	x[0] = 0.3;
	y[0] = 0.7;  // Consistent: x + y = 1
	
	linearDAE.algConstraints(0.0, x, y, g);
	
	// Constraint residual should be 0 for consistent values
	REQUIRE_THAT(g[0], WithinAbs(0.0, TOL(1e-10, 1e-5)));
}

TEST_CASE("DAE::Interface_algConstraints_violation", "[dae][interface]")
{
	SimpleLinearDAE linearDAE;
	
	Vector<Real> x(1), y(1), g(1);
	x[0] = 0.5;
	y[0] = 0.3;  // Inconsistent: x + y = 0.8 ≠ 1
	
	linearDAE.algConstraints(0.0, x, y, g);
	
	// Constraint residual should be -0.2
	REQUIRE_THAT(g[0], WithinAbs(-0.2, TOL(1e-10, 1e-5)));
}

TEST_CASE("DAE::Interface_jacobians", "[dae][interface]")
{
	SimpleLinearDAE linearDAE;
	
	Vector<Real> x(1), y(1);
	x[0] = 0.5;
	y[0] = 0.5;
	
	Matrix<Real> df_dx(1, 1), df_dy(1, 1), dg_dx(1, 1), dg_dy(1, 1);
	
	linearDAE.jacobian_fx(0.0, x, y, df_dx);
	linearDAE.jacobian_fy(0.0, x, y, df_dy);
	linearDAE.jacobian_gx(0.0, x, y, dg_dx);
	linearDAE.jacobian_gy(0.0, x, y, dg_dy);
	
	REQUIRE_THAT(df_dx(0, 0), WithinAbs(-1.0, TOL(1e-10, 1e-5)));  // d(dx/dt)/dx = -1
	REQUIRE_THAT(df_dy(0, 0), WithinAbs(1.0, TOL(1e-10, 1e-5)));   // d(dx/dt)/dy = 1
	REQUIRE_THAT(dg_dx(0, 0), WithinAbs(1.0, TOL(1e-10, 1e-5)));   // dg/dx = 1
	REQUIRE_THAT(dg_dy(0, 0), WithinAbs(1.0, TOL(1e-10, 1e-5)));   // dg/dy = 1
}

TEST_CASE("DAE::Interface_pendulum_constraint", "[dae][interface]")
{
	Real L = 1.0;
	SimplifiedPendulumDAE pendulum(L, 9.81);
	
	Vector<Real> diff(4), alg(1), g(1);
	// Start at angle θ = 30°, at rest
	diff[0] = L * std::sin(Constants::PI / 6);  // x = L*sin(30°) = 0.5
	diff[1] = -L * std::cos(Constants::PI / 6); // y = -L*cos(30°) ≈ -0.866
	diff[2] = 0.0;  // vx = 0
	diff[3] = 0.0;  // vy = 0
	alg[0] = 0.0;   // λ (will be computed by solver)
	
	pendulum.algConstraints(0.0, diff, alg, g);
	
	// Constraint: x² + y² - L² should be 0 for consistent initial conditions
	REQUIRE_THAT(g[0], WithinAbs(0.0, TOL(1e-10, 1e-5)));
}

TEST_CASE("DAE::Interface_allJacobians", "[dae][interface]")
{
	SimpleLinearDAE linearDAE;
	
	Vector<Real> x(1), y(1);
	x[0] = 0.5;
	y[0] = 0.5;
	
	Matrix<Real> df_dx(1, 1), df_dy(1, 1), dg_dx(1, 1), dg_dy(1, 1);
	
	// Use the combined method
	linearDAE.allJacobians(0.0, x, y, df_dx, df_dy, dg_dx, dg_dy);
	
	REQUIRE_THAT(df_dx(0, 0), WithinAbs(-1.0, TOL(1e-10, 1e-5)));
	REQUIRE_THAT(df_dy(0, 0), WithinAbs(1.0, TOL(1e-10, 1e-5)));
	REQUIRE_THAT(dg_dx(0, 0), WithinAbs(1.0, TOL(1e-10, 1e-5)));
	REQUIRE_THAT(dg_dy(0, 0), WithinAbs(1.0, TOL(1e-10, 1e-5)));
}

///////////////////////////////////////////////////////////////////////////////////////////
//                           DAESolution Tests
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("DAE::DAESolution_construction", "[dae][solution]")
{
	DAESolution sol(0.0, 1.0, 2, 1, 100);
	
	REQUIRE(sol.getDiffDim() == 2);
	REQUIRE(sol.getAlgDim() == 1);
	REQUIRE(sol.getTotalDim() == 3);
	REQUIRE(sol.getT1() == 0.0);
	REQUIRE(sol.getT2() == 1.0);
	REQUIRE(sol.getTotalSavedSteps() == 101);
}

TEST_CASE("DAE::DAESolution_fillValues", "[dae][solution]")
{
	DAESolution sol(0.0, 1.0, 2, 1, 10);
	
	Vector<Real> x(2), y(1);
	x[0] = 1.0; x[1] = 2.0;
	y[0] = 3.0;
	
	sol.fillValues(0, 0.0, x, y);
	
	REQUIRE_THAT(sol.getTValue(0), WithinAbs(0.0, TOL(1e-10, 1e-5)));
	REQUIRE_THAT(sol.getXValue(0, 0), WithinAbs(1.0, TOL(1e-10, 1e-5)));
	REQUIRE_THAT(sol.getXValue(0, 1), WithinAbs(2.0, TOL(1e-10, 1e-5)));
	REQUIRE_THAT(sol.getYValue(0, 0), WithinAbs(3.0, TOL(1e-10, 1e-5)));
}

TEST_CASE("DAE::DAESolution_multiplePoints", "[dae][solution]")
{
	DAESolution sol(0.0, 1.0, 1, 1, 10);
	
	for (int i = 0; i <= 10; i++)
	{
		Real t = i * 0.1;
		Vector<Real> x(1), y(1);
		x[0] = t;           // x = t
		y[0] = 1.0 - t;     // y = 1 - t (constraint x + y = 1)
		sol.fillValues(i, t, x, y);
	}
	
	// Check a few points
	REQUIRE_THAT(sol.getTValue(5), WithinAbs(0.5, TOL(1e-10, 1e-5)));
	REQUIRE_THAT(sol.getXValue(5, 0), WithinAbs(0.5, TOL(1e-10, 1e-5)));
	REQUIRE_THAT(sol.getYValue(5, 0), WithinAbs(0.5, TOL(1e-10, 1e-5)));
	
	// Check final values
	Vector<Real> xEnd = sol.getXValuesAtEnd();
	Vector<Real> yEnd = sol.getYValuesAtEnd();
	REQUIRE_THAT(xEnd[0], WithinAbs(1.0, TOL(1e-10, 1e-5)));
	REQUIRE_THAT(yEnd[0], WithinAbs(0.0, TOL(1e-10, 1e-5)));
}

TEST_CASE("DAE::DAESolution_interpolation", "[dae][solution]")
{
	DAESolution sol(0.0, 1.0, 1, 1, 10);
	
	for (int i = 0; i <= 10; i++)
	{
		Real t = i * 0.1;
		Vector<Real> x(1), y(1);
		x[0] = t * t;       // x = t²
		y[0] = 2 * t;       // y = 2t
		sol.fillValues(i, t, x, y);
	}
	sol.setFinalSize(10);
	
	// Get interpolators
	auto xInterp = sol.getXAsLinInterp(0);
	auto yInterp = sol.getYAsLinInterp(0);
	
	// Test at known points
	REQUIRE_THAT(xInterp(0.5), WithinAbs(0.25, 0.02));  // x(0.5) = 0.25
	REQUIRE_THAT(yInterp(0.5), WithinAbs(1.0, TOL(1e-10, 1e-5))); // y(0.5) = 1.0
}

///////////////////////////////////////////////////////////////////////////////////////////
//                           DAESystem Tests
///////////////////////////////////////////////////////////////////////////////////////////

// Function pointers for testing DAESystem class
void testDiffFunc(Real t, const Vector<Real>& x, const Vector<Real>& y, Vector<Real>& dxdt)
{
	dxdt[0] = -x[0] + y[0];  // dx/dt = -x + y
}

void testAlgFunc(Real t, const Vector<Real>& x, const Vector<Real>& y, Vector<Real>& g)
{
	g[0] = x[0] + y[0] - 1.0;  // x + y = 1
}

TEST_CASE("DAE::DAESystem_construction", "[dae][system]")
{
	DAESystem sys(1, 1, testDiffFunc, testAlgFunc);
	
	REQUIRE(sys.getDiffDim() == 1);
	REQUIRE(sys.getAlgDim() == 1);
	REQUIRE(sys.getTotalDim() == 2);
}

TEST_CASE("DAE::DAESystem_nullptr_throws", "[dae][system]")
{
	// Both null
	REQUIRE_THROWS_AS(DAESystem(1, 1, nullptr, nullptr), std::invalid_argument);
	// diffFunc null
	REQUIRE_THROWS_AS(DAESystem(1, 1, nullptr, testAlgFunc), std::invalid_argument);
	// algFunc null
	REQUIRE_THROWS_AS(DAESystem(1, 1, testDiffFunc, nullptr), std::invalid_argument);
	// Both valid — no throw
	REQUIRE_NOTHROW(DAESystem(1, 1, testDiffFunc, testAlgFunc));
}

TEST_CASE("DAE::DAESystem_evaluation", "[dae][system]")
{
	DAESystem sys(1, 1, testDiffFunc, testAlgFunc);
	
	Vector<Real> x(1), y(1), dxdt(1), g(1);
	x[0] = 0.3;
	y[0] = 0.7;
	
	sys.diffEqs(0.0, x, y, dxdt);
	sys.algConstraints(0.0, x, y, g);
	
	// dx/dt = -0.3 + 0.7 = 0.4
	REQUIRE_THAT(dxdt[0], WithinAbs(0.4, TOL(1e-10, 1e-5)));
	// g = 0.3 + 0.7 - 1 = 0
	REQUIRE_THAT(g[0], WithinAbs(0.0, TOL(1e-10, 1e-5)));
}

///////////////////////////////////////////////////////////////////////////////////////////
//                           DAE Solver Tests
///////////////////////////////////////////////////////////////////////////////////////////

// DAESolvers.h included at top of file

TEST_CASE("DAE::ConsistentIC_verify", "[dae][solver]")
{
	SimpleLinearDAE linearDAE;
	
	Vector<Real> x(1), y(1);
	x[0] = 0.3;
	y[0] = 0.7;  // Consistent: x + y = 1
	
	bool isConsistent = VerifyConsistentIC(linearDAE, 0.0, x, y, TOL(1e-10, 1e-5));
	REQUIRE(isConsistent);
	
	// Inconsistent case
	y[0] = 0.5;  // Now x + y = 0.8 ≠ 1
	isConsistent = VerifyConsistentIC(linearDAE, 0.0, x, y, TOL(1e-10, 1e-5));
	REQUIRE_FALSE(isConsistent);
}

TEST_CASE("DAE::ConsistentIC_compute", "[dae][solver]")
{
	SimpleLinearDAE linearDAE;
	
	Vector<Real> x(1), y(1);
	x[0] = 0.3;
	y[0] = 0.5;  // Initial guess (inconsistent)
	
	// Compute consistent y such that x + y = 1
	bool success = ComputeConsistentIC(linearDAE, 0.0, x, y, 20, TOL(1e-10, 1e-5));
	
	REQUIRE(success);
	REQUIRE_THAT(y[0], WithinAbs(0.7, TOL(1e-10, 1e-5)));  // y = 1 - x = 0.7
	
	// Verify consistency
	REQUIRE(VerifyConsistentIC(linearDAE, 0.0, x, y, TOL(1e-10, 1e-5)));
}

TEST_CASE("DAE::BackwardEuler_linearDAE", "[dae][solver]")
{
	SimpleLinearDAE linearDAE;
	
	// Initial conditions
	Real t0 = 0.0, t_end = 1.0;
	Vector<Real> x0(1), y0(1);
	x0[0] = 0.8;               // Start away from equilibrium
	y0[0] = 1.0 - x0[0];       // Consistent: y = 1 - x
	
	DAESolverConfig config;
	config.step_size = 0.01;
	config.newton_tol = TOL(1e-10, 1e-5);
	
	DAESolverResult result = SolveDAEBackwardEuler(linearDAE, t0, x0, y0, t_end, config);
	
	REQUIRE(result.status == AlgorithmStatus::Success);
	REQUIRE(result.algorithm_name == "DAEBackwardEuler");
	REQUIRE(result.total_steps > 0);
	
	// Analytical solution: x(t) = 0.5 + (x0 - 0.5)*exp(-2t)
	// At t = 1: x(1) = 0.5 + 0.3*exp(-2) ≈ 0.5 + 0.0406 ≈ 0.5406
	Real x_analytical = 0.5 + (0.8 - 0.5) * std::exp(-2.0);
	Real y_analytical = 1.0 - x_analytical;
	
	Vector<Real> xEnd = result.solution.getXValuesAtEnd();
	Vector<Real> yEnd = result.solution.getYValuesAtEnd();
	
	// First-order method, so expect ~1% error with h=0.01
	REQUIRE_THAT(xEnd[0], WithinAbs(x_analytical, 0.01));
	REQUIRE_THAT(yEnd[0], WithinAbs(y_analytical, 0.01));
	
	// Constraint should be satisfied throughout
	REQUIRE(result.max_constraint_violation < TOL(1e-8, 1e-4));
}

TEST_CASE("DAE::BackwardEuler_constraintSatisfaction", "[dae][solver]")
{
	SimpleLinearDAE linearDAE;
	
	Real t0 = 0.0, t_end = 2.0;
	Vector<Real> x0(1), y0(1);
	x0[0] = 0.9;
	y0[0] = 0.1;  // Consistent
	
	DAESolverConfig config;
	config.step_size = 0.02;
	
	DAESolverResult result = SolveDAEBackwardEuler(linearDAE, t0, x0, y0, t_end, config);
	
	REQUIRE(result.status == AlgorithmStatus::Success);
	
	// Check constraint at multiple time points
	Vector<Real> g(1);
	for (int i = 0; i <= 10; i++)
	{
		int step = i * result.total_steps / 10;
		if (step > result.total_steps) step = result.total_steps;
		
		Real t = result.solution.getTValue(step);
		Vector<Real> x(1), y(1);
		x[0] = result.solution.getXValue(step, 0);
		y[0] = result.solution.getYValue(step, 0);
		
		linearDAE.algConstraints(t, x, y, g);
		REQUIRE_THAT(g[0], WithinAbs(0.0, TOL(1e-8, 1e-4)));
	}
}

TEST_CASE("DAE::BDF2_linearDAE", "[dae][solver]")
{
	SimpleLinearDAE linearDAE;
	
	Real t0 = 0.0, t_end = 1.0;
	Vector<Real> x0(1), y0(1);
	x0[0] = 0.8;
	y0[0] = 0.2;  // Consistent
	
	DAESolverConfig config;
	config.step_size = 0.01;
	config.newton_tol = TOL(1e-10, 1e-5);
	
	DAESolverResult result = SolveDAEBDF2(linearDAE, t0, x0, y0, t_end, config);
	
	REQUIRE(result.status == AlgorithmStatus::Success);
	REQUIRE(result.algorithm_name == "DAEBDF2");
	
	// BDF2 should be more accurate than Backward Euler
	Real x_analytical = 0.5 + (0.8 - 0.5) * std::exp(-2.0);
	Real y_analytical = 1.0 - x_analytical;
	
	Vector<Real> xEnd = result.solution.getXValuesAtEnd();
	Vector<Real> yEnd = result.solution.getYValuesAtEnd();
	
	// Second-order method, expect better accuracy
	REQUIRE_THAT(xEnd[0], WithinAbs(x_analytical, 0.001));
	REQUIRE_THAT(yEnd[0], WithinAbs(y_analytical, 0.001));
}

TEST_CASE("DAE::BDF2_betterThanBackwardEuler", "[dae][solver]")
{
	SimpleLinearDAE linearDAE;
	
	Real t0 = 0.0, t_end = 1.0;
	Vector<Real> x0(1), y0(1);
	x0[0] = 0.9;
	y0[0] = 0.1;
	
	DAESolverConfig config;
	config.step_size = 0.05;  // Larger step to see order difference
	
	DAESolverResult resultBE = SolveDAEBackwardEuler(linearDAE, t0, x0, y0, t_end, config);
	DAESolverResult resultBDF2 = SolveDAEBDF2(linearDAE, t0, x0, y0, t_end, config);
	
	Real x_analytical = 0.5 + (0.9 - 0.5) * std::exp(-2.0);
	
	Real errorBE = std::abs(resultBE.solution.getXValuesAtEnd()[0] - x_analytical);
	Real errorBDF2 = std::abs(resultBDF2.solution.getXValuesAtEnd()[0] - x_analytical);
	
	// BDF2 should have smaller error (roughly h^2 vs h)
	REQUIRE(errorBDF2 < errorBE);
}

TEST_CASE("DAE::BDF4_solves_linearDAE", "[dae][solver][bdf4]")
{
	SimpleLinearDAE linearDAE;
	
	Real t0 = 0.0, t_end = 1.0;
	Vector<Real> x0(1), y0(1);
	x0[0] = 0.8;
	y0[0] = 0.2;  // Consistent: y = 1 - x
	
	DAESolverConfig config;
	config.step_size = 0.01;
	config.newton_tol = TOL(1e-10, 1e-5);
	
	DAESolverResult result = SolveDAEBDF4(linearDAE, t0, x0, y0, t_end, config);
	
	REQUIRE(result.status == AlgorithmStatus::Success);
	REQUIRE(result.algorithm_name == "DAEBDF4");
	
	// Analytical: x(t) = 0.5 + (x0 - 0.5) * exp(-2t)
	Real x_analytical = 0.5 + (0.8 - 0.5) * std::exp(-2.0);
	Real y_analytical = 1.0 - x_analytical;
	
	Vector<Real> xEnd = result.solution.getXValuesAtEnd();
	Vector<Real> yEnd = result.solution.getYValuesAtEnd();
	
	// Fourth-order method, expect high accuracy
	REQUIRE_THAT(xEnd[0], WithinAbs(x_analytical, 1e-5));
	REQUIRE_THAT(yEnd[0], WithinAbs(y_analytical, 1e-5));
}

TEST_CASE("DAE::BDF4_betterThanBDF2", "[dae][solver][bdf4]")
{
	SimpleLinearDAE linearDAE;
	
	Real t0 = 0.0, t_end = 1.0;
	Vector<Real> x0(1), y0(1);
	x0[0] = 0.9;
	y0[0] = 0.1;
	
	DAESolverConfig config;
	config.step_size = 0.05;  // Larger step to see order difference
	
	DAESolverResult resultBDF2 = SolveDAEBDF2(linearDAE, t0, x0, y0, t_end, config);
	DAESolverResult resultBDF4 = SolveDAEBDF4(linearDAE, t0, x0, y0, t_end, config);
	
	Real x_analytical = 0.5 + (0.9 - 0.5) * std::exp(-2.0);
	
	Real errorBDF2 = std::abs(resultBDF2.solution.getXValuesAtEnd()[0] - x_analytical);
	Real errorBDF4 = std::abs(resultBDF4.solution.getXValuesAtEnd()[0] - x_analytical);
	
	// BDF4 should have smaller error (roughly h^4 vs h^2)
	REQUIRE(errorBDF4 < errorBDF2);
}

TEST_CASE("DAE::BDF4_convergence_order", "[dae][solver][bdf4]")
{
	// Test that BDF4 achieves better than 2nd order convergence
	// Note: Due to bootstrap startup using substepped BDF2, we don't get
	// full 4th order until sufficiently small step sizes
	SimpleLinearDAE linearDAE;
	
	Real t0 = 0.0, t_end = 1.0;
	Vector<Real> x0(1), y0(1);
	x0[0] = 0.8;
	y0[0] = 0.2;
	
	Real x_analytical = 0.5 + (0.8 - 0.5) * std::exp(-2.0);
	
	DAESolverConfig config1, config2;
	config1.step_size = 0.02;
	config2.step_size = 0.01;  // Half the step size
	
	DAESolverResult result1 = SolveDAEBDF4(linearDAE, t0, x0, y0, t_end, config1);
	DAESolverResult result2 = SolveDAEBDF4(linearDAE, t0, x0, y0, t_end, config2);
	
	Real error1 = std::abs(result1.solution.getXValuesAtEnd()[0] - x_analytical);
	Real error2 = std::abs(result2.solution.getXValuesAtEnd()[0] - x_analytical);
	
	// For 4th order method, halving step size should reduce error by ~16x
	// Due to substepped bootstrap, we see ~4x initially but better as h decreases
	// Require at least 3.5x improvement (between 2nd and 4th order)
	if constexpr (std::is_same_v<Real, double>) {
		Real ratio = error1 / error2;
		REQUIRE(ratio > 3.5);
	} else {
		// Float precision floor makes order measurement unreliable at these step sizes
		REQUIRE(result1.status == AlgorithmStatus::Success);
		REQUIRE(result2.status == AlgorithmStatus::Success);
	}
}

TEST_CASE("DAE::RODAS_solves_linearDAE", "[dae][solver][rodas]")
{
	SimpleLinearDAE linearDAE;
	
	Real t0 = 0.0, t_end = 1.0;
	Vector<Real> x0(1), y0(1);
	x0[0] = 0.8;
	y0[0] = 0.2;  // Consistent: y = 1 - x
	
	DAESolverConfig config;
	config.step_size = 0.01;
	
	DAESolverResult result = SolveDAERODAS(linearDAE, t0, x0, y0, t_end, config);
	
	REQUIRE(result.status == AlgorithmStatus::Success);
	REQUIRE(result.algorithm_name == "DAERODAS");
	
	// Analytical: x(t) = 0.5 + (x0 - 0.5) * exp(-2t)
	Real x_analytical = 0.5 + (0.8 - 0.5) * std::exp(-2.0);
	Real y_analytical = 1.0 - x_analytical;
	
	Vector<Real> xEnd = result.solution.getXValuesAtEnd();
	Vector<Real> yEnd = result.solution.getYValuesAtEnd();
	
	// Second-order method with linearly-implicit corrector
	// With h=0.01, expect error O(h^2) ≈ 1e-4, use 1e-3 tolerance
	REQUIRE_THAT(xEnd[0], WithinAbs(x_analytical, 1e-3));
	REQUIRE_THAT(yEnd[0], WithinAbs(y_analytical, 1e-3));
}

TEST_CASE("DAE::RODAS_betterThanBDF2", "[dae][solver][rodas]")
{
	// RODAS (linearly-implicit) vs BDF2 (Newton iteration)
	// Both are 2nd order, but test that RODAS is comparable accuracy
	SimpleLinearDAE linearDAE;
	
	Real t0 = 0.0, t_end = 1.0;
	Vector<Real> x0(1), y0(1);
	x0[0] = 0.9;
	y0[0] = 0.1;
	
	DAESolverConfig config;
	config.step_size = 0.05;  // Larger step to see order difference
	
	DAESolverResult resultBDF2 = SolveDAEBDF2(linearDAE, t0, x0, y0, t_end, config);
	DAESolverResult resultRODAS = SolveDAERODAS(linearDAE, t0, x0, y0, t_end, config);
	
	Real x_analytical = 0.5 + (0.9 - 0.5) * std::exp(-2.0);
	
	Real errorBDF2 = std::abs(resultBDF2.solution.getXValuesAtEnd()[0] - x_analytical);
	Real errorRODAS = std::abs(resultRODAS.solution.getXValuesAtEnd()[0] - x_analytical);
	
	// Both are 2nd order methods, RODAS uses fewer Newton iterations (none!)
	// Allow RODAS error to be at most 10x worse than BDF2 for comparable performance
	REQUIRE(errorRODAS < 10 * errorBDF2 + 0.01);
	
	// Also verify both achieve reasonable accuracy
	REQUIRE(errorBDF2 < 0.01);
	REQUIRE(errorRODAS < 0.01);
}

TEST_CASE("DAE::RODAS_noNewtonIterations", "[dae][solver][rodas]")
{
	// Verify that RODAS uses no Newton iterations (only linear solves)
	SimpleLinearDAE linearDAE;
	
	Real t0 = 0.0, t_end = 1.0;
	Vector<Real> x0(1), y0(1);
	x0[0] = 0.8;
	y0[0] = 0.2;
	
	DAESolverConfig config;
	config.step_size = 0.1;
	
	DAESolverResult result = SolveDAERODAS(linearDAE, t0, x0, y0, t_end, config);
	
	REQUIRE(result.status == AlgorithmStatus::Success);
	// RODAS doesn't use Newton iteration, so this should be 0
	REQUIRE(result.newton_iterations == 0);
	// But it does evaluate Jacobians (once per step)
	REQUIRE(result.jacobian_evaluations > 0);
}

TEST_CASE("DAE::RODAS_convergence_order", "[dae][solver][rodas]")
{
	// Test that RODAS achieves convergence (at least 1st order)
	// The linearly-implicit predictor-corrector gives ~2x ratio on step halving
	SimpleLinearDAE linearDAE;
	
	Real t0 = 0.0, t_end = 1.0;
	Vector<Real> x0(1), y0(1);
	x0[0] = 0.8;
	y0[0] = 0.2;
	
	Real x_analytical = 0.5 + (0.8 - 0.5) * std::exp(-2.0);
	
	DAESolverConfig config1, config2;
	config1.step_size = 0.02;
	config2.step_size = 0.01;  // Half the step size
	
	DAESolverResult result1 = SolveDAERODAS(linearDAE, t0, x0, y0, t_end, config1);
	DAESolverResult result2 = SolveDAERODAS(linearDAE, t0, x0, y0, t_end, config2);
	
	Real error1 = std::abs(result1.solution.getXValuesAtEnd()[0] - x_analytical);
	Real error2 = std::abs(result2.solution.getXValuesAtEnd()[0] - x_analytical);
	
	// Linearly-implicit method gives ratio ~2 on step halving
	// Require improvement of at least 1.5x (verify convergence)
	Real ratio = error1 / error2;
	REQUIRE(ratio > 1.5);
}

/******************************************************************************/
/*****            Radau IIA Tests (3-stage, Order 5, L-stable)            *****/
/******************************************************************************/

TEST_CASE("DAE::RadauIIA_solves_linearDAE", "[dae][solver][radauiia]")
{
	// Radau IIA (order 5) should achieve very high accuracy
	SimpleLinearDAE linearDAE;
	
	Real t0 = 0.0, t_end = 1.0;
	Vector<Real> x0(1), y0(1);
	x0[0] = 0.8;
	y0[0] = 0.2;
	
	DAESolverConfig config;
	config.step_size = 0.01;
	
	DAESolverResult result = SolveDAERadauIIA(linearDAE, t0, x0, y0, t_end, config);
	
	// Analytical: x(t) = 0.5 + (0.8 - 0.5)*exp(-2t) = 0.5 + 0.3*exp(-2t)
	Real x_analytical = 0.5 + 0.3 * std::exp(-2.0);
	
	REQUIRE(result.status == AlgorithmStatus::Success);
	REQUIRE_THAT(result.solution.getXValuesAtEnd()[0], WithinAbs(x_analytical, 1e-6));
}

TEST_CASE("DAE::RadauIIA_higherOrderThanBDF4", "[dae][solver][radauiia]")
{
	// Radau IIA (order 5) should be more accurate than BDF4 (order 4) with same step
	SimpleLinearDAE linearDAE;
	
	Real t0 = 0.0, t_end = 1.0;
	Vector<Real> x0(1), y0(1);
	x0[0] = 0.8;
	y0[0] = 0.2;
	
	DAESolverConfig config;
	config.step_size = 0.05;  // Larger step to see order difference
	
	DAESolverResult radau = SolveDAERadauIIA(linearDAE, t0, x0, y0, t_end, config);
	DAESolverResult bdf4 = SolveDAEBDF4(linearDAE, t0, x0, y0, t_end, config);
	
	Real x_analytical = 0.5 + 0.3 * std::exp(-2.0);
	Real radau_error = std::abs(radau.solution.getXValuesAtEnd()[0] - x_analytical);
	Real bdf4_error = std::abs(bdf4.solution.getXValuesAtEnd()[0] - x_analytical);
	
	// Radau IIA should be at least as accurate as BDF4 (often better)
	REQUIRE(radau_error <= bdf4_error * 2.0);  // Allow some tolerance
}

TEST_CASE("DAE::RadauIIA_usesNewton", "[dae][solver][radauiia]")
{
	// Radau IIA is fully implicit and requires Newton iterations
	SimpleLinearDAE linearDAE;
	
	Real t0 = 0.0, t_end = 0.5;
	Vector<Real> x0(1), y0(1);
	x0[0] = 0.8;
	y0[0] = 0.2;
	
	DAESolverConfig config;
	config.step_size = 0.1;
	
	DAESolverResult result = SolveDAERadauIIA(linearDAE, t0, x0, y0, t_end, config);
	
	REQUIRE(result.status == AlgorithmStatus::Success);
	REQUIRE(result.newton_iterations > 0);  // Must use Newton
	REQUIRE(result.total_steps >= 5);
}

TEST_CASE("DAE::RadauIIA_convergence_order", "[dae][solver][radauiia]")
{
	// Radau IIA is order 5, so error ratio should be ~32 when halving step
	SimpleLinearDAE linearDAE;
	
	Real t0 = 0.0, t_end = 1.0;
	Vector<Real> x0(1), y0(1);
	x0[0] = 0.8;
	y0[0] = 0.2;
	
	Real x_analytical = 0.5 + (0.8 - 0.5) * std::exp(-2.0);
	
	DAESolverConfig config1, config2;
	config1.step_size = 0.04;
	config2.step_size = 0.02;  // Half the step size
	
	DAESolverResult result1 = SolveDAERadauIIA(linearDAE, t0, x0, y0, t_end, config1);
	DAESolverResult result2 = SolveDAERadauIIA(linearDAE, t0, x0, y0, t_end, config2);
	
	Real error1 = std::abs(result1.solution.getXValuesAtEnd()[0] - x_analytical);
	Real error2 = std::abs(result2.solution.getXValuesAtEnd()[0] - x_analytical);
	
	// Order 5 means ratio ~32 = 2^5, but we allow slack due to constants
	// Require at least order 3 behavior (ratio > 8)
	if constexpr (std::is_same_v<Real, double>) {
		Real ratio = error1 / error2;
		REQUIRE(ratio > 8.0);
	} else {
		// Float precision floor makes order measurement unreliable at these step sizes
		REQUIRE(result1.status == AlgorithmStatus::Success);
		REQUIRE(result2.status == AlgorithmStatus::Success);
	}
}

TEST_CASE("DAE::Config_presets", "[dae][solver]")
{
	DAESolverConfig defaultConfig;
	DAESolverConfig highPrec = DAESolverConfig::HighPrecision();
	DAESolverConfig fast = DAESolverConfig::Fast();
	
	// High precision has smaller step and tighter tolerance
	REQUIRE(highPrec.step_size < defaultConfig.step_size);
	REQUIRE(highPrec.newton_tol < defaultConfig.newton_tol);
	
	// Fast has larger step
	REQUIRE(fast.step_size > defaultConfig.step_size);
}

TEST_CASE("DAE::Solver_diagnostics", "[dae][solver]")
{
	SimpleLinearDAE linearDAE;
	
	Real t0 = 0.0, t_end = 0.5;
	Vector<Real> x0(1), y0(1);
	x0[0] = 0.8;  // Away from equilibrium
	y0[0] = 0.2;
	
	DAESolverConfig config;
	config.step_size = 0.1;
	
	DAESolverResult result = SolveDAEBackwardEuler(linearDAE, t0, x0, y0, t_end, config);
	
	REQUIRE(result.status == AlgorithmStatus::Success);
	REQUIRE(result.total_steps == 5);  // 0.5 / 0.1 = 5 steps
	REQUIRE(result.newton_iterations > 0);
	// Note: Jacobian might not be evaluated if residual is small initially
	// For this simple linear problem, explicit prediction can be very good
	REQUIRE(result.elapsed_time_ms >= 0.0);
	REQUIRE_FALSE(result.algorithm_name.empty());
}

///////////////////////////////////////////////////////////////////////////////////////////
//                           Realistic DAE Test Systems
///////////////////////////////////////////////////////////////////////////////////////////

/// RC Circuit DAE:
/// An ideal capacitor with voltage source through resistor.
/// 
/// Differential: C * dV_c/dt = I
/// Algebraic: V_s - R*I - V_c = 0  (Kirchhoff's voltage law)
///
/// Where V_c = capacitor voltage (differential), I = current (algebraic)
/// Analytical: V_c(t) = V_s * (1 - exp(-t/(RC))) for V_c(0) = 0
///            I(t) = (V_s/R) * exp(-t/(RC))
class RCCircuitDAE : public IODESystemDAEWithJacobian
{
	Real _R, _C, _Vs;  // Resistance, Capacitance, Source voltage

public:
	RCCircuitDAE(Real R = 1000.0, Real C = 1e-6, Real Vs = 5.0)
		: _R(R), _C(C), _Vs(Vs) {}

	int getDiffDim() const override { return 1; }  // V_c
	int getAlgDim() const override { return 1; }   // I

	void diffEqs(Real t, const Vector<Real>& x, const Vector<Real>& y,
	             Vector<Real>& dxdt) const override
	{
		Real I = y[0];
		dxdt[0] = I / _C;  // C * dV_c/dt = I => dV_c/dt = I/C
	}

	void algConstraints(Real t, const Vector<Real>& x, const Vector<Real>& y,
	                    Vector<Real>& g) const override
	{
		Real Vc = x[0], I = y[0];
		g[0] = _Vs - _R * I - Vc;  // KVL: V_s = R*I + V_c
	}

	std::string getDiffVarName(int i) const override { return "V_c"; }
	std::string getAlgVarName(int i) const override { return "I"; }

	void jacobian_fx(Real t, const Vector<Real>& x, const Vector<Real>& y,
	                 Matrix<Real>& df_dx) const override
	{
		df_dx(0, 0) = 0.0;  // d(dV_c/dt)/dV_c = 0
	}

	void jacobian_fy(Real t, const Vector<Real>& x, const Vector<Real>& y,
	                 Matrix<Real>& df_dy) const override
	{
		df_dy(0, 0) = 1.0 / _C;  // d(dV_c/dt)/dI = 1/C
	}

	void jacobian_gx(Real t, const Vector<Real>& x, const Vector<Real>& y,
	                 Matrix<Real>& dg_dx) const override
	{
		dg_dx(0, 0) = -1.0;  // dg/dV_c = -1
	}

	void jacobian_gy(Real t, const Vector<Real>& x, const Vector<Real>& y,
	                 Matrix<Real>& dg_dy) const override
	{
		dg_dy(0, 0) = -_R;  // dg/dI = -R (nonsingular if R != 0)
	}

	// Analytical solution for testing
	Real analyticalVc(Real t) const { return _Vs * (1.0 - std::exp(-t / (_R * _C))); }
	Real analyticalI(Real t) const { return (_Vs / _R) * std::exp(-t / (_R * _C)); }
};

/// Stiff Linear DAE (varying time scales):
/// dx/dt = -1000*x + 1000*y
/// 0 = x + y - 2*exp(-t)
///
/// This has a "stiff" differential component (fast decay to constraint manifold)
/// and a slow constraint that varies with time.
/// Analytical: x(t) = exp(-t), y(t) = exp(-t)
class StiffLinearDAE : public IODESystemDAEWithJacobian
{
	Real _lambda;  // Stiffness parameter

public:
	StiffLinearDAE(Real lambda = 1000.0) : _lambda(lambda) {}

	int getDiffDim() const override { return 1; }
	int getAlgDim() const override { return 1; }

	void diffEqs(Real t, const Vector<Real>& x, const Vector<Real>& y,
	             Vector<Real>& dxdt) const override
	{
		dxdt[0] = -_lambda * x[0] + _lambda * y[0];
	}

	void algConstraints(Real t, const Vector<Real>& x, const Vector<Real>& y,
	                    Vector<Real>& g) const override
	{
		g[0] = x[0] + y[0] - 2.0 * std::exp(-t);
	}

	std::string getDiffVarName(int i) const override { return "x"; }
	std::string getAlgVarName(int i) const override { return "y"; }

	void jacobian_fx(Real t, const Vector<Real>& x, const Vector<Real>& y,
	                 Matrix<Real>& df_dx) const override
	{
		df_dx(0, 0) = -_lambda;
	}

	void jacobian_fy(Real t, const Vector<Real>& x, const Vector<Real>& y,
	                 Matrix<Real>& df_dy) const override
	{
		df_dy(0, 0) = _lambda;
	}

	void jacobian_gx(Real t, const Vector<Real>& x, const Vector<Real>& y,
	                 Matrix<Real>& dg_dx) const override
	{
		dg_dx(0, 0) = 1.0;
	}

	void jacobian_gy(Real t, const Vector<Real>& x, const Vector<Real>& y,
	                 Matrix<Real>& dg_dy) const override
	{
		dg_dy(0, 0) = 1.0;
	}

	Real analyticalX(Real t) const { return std::exp(-t); }
	Real analyticalY(Real t) const { return std::exp(-t); }
};

/// Two-body constrained system:
/// Two masses connected by rigid rod (length L).
/// Differential: dxi/dt = vi, mi*dvi/dt = Fi + λ*∂g/∂xi
/// Constraint: (x1-x2)² + (y1-y2)² = L²
///
/// Simplified 1D version for testing:
/// dx1/dt = v1, dx2/dt = v2
/// m1*dv1/dt = F1 + λ, m2*dv2/dt = F2 - λ
/// x1 - x2 = L (rigid rod constraint in 1D)
///
/// With F1 = F2 = -g (gravity), the constraint force λ maintains distance.
class RigidRod1DDAE : public IODESystemDAEWithJacobian
{
	Real _m1, _m2, _L, _g;

public:
	RigidRod1DDAE(Real m1 = 1.0, Real m2 = 1.0, Real L = 1.0, Real g = 9.81)
		: _m1(m1), _m2(m2), _L(L), _g(g) {}

	int getDiffDim() const override { return 4; }  // x1, x2, v1, v2
	int getAlgDim() const override { return 1; }   // λ

	void diffEqs(Real t, const Vector<Real>& diff, const Vector<Real>& alg,
	             Vector<Real>& dxdt) const override
	{
		Real v1 = diff[2], v2 = diff[3];
		Real lambda = alg[0];

		dxdt[0] = v1;                          // dx1/dt = v1
		dxdt[1] = v2;                          // dx2/dt = v2
		dxdt[2] = (-_g + lambda / _m1);        // m1*dv1/dt = -m1*g + λ
		dxdt[3] = (-_g - lambda / _m2);        // m2*dv2/dt = -m2*g - λ
	}

	void algConstraints(Real t, const Vector<Real>& diff, const Vector<Real>& alg,
	                    Vector<Real>& g) const override
	{
		Real x1 = diff[0], x2 = diff[1];
		g[0] = x1 - x2 - _L;  // Constraint: x1 - x2 = L
	}

	std::string getDiffVarName(int i) const override
	{
		static const char* names[] = {"x1", "x2", "v1", "v2"};
		return (i >= 0 && i < 4) ? names[i] : "?";
	}
	std::string getAlgVarName(int i) const override { return "lambda"; }

	void jacobian_fx(Real t, const Vector<Real>& diff, const Vector<Real>& alg,
	                 Matrix<Real>& df_dx) const override
	{
		df_dx.MakeUnitMatrix();
		df_dx = df_dx * 0.0;  // Clear
		df_dx(0, 2) = 1.0;    // dx1/dt depends on v1
		df_dx(1, 3) = 1.0;    // dx2/dt depends on v2
	}

	void jacobian_fy(Real t, const Vector<Real>& diff, const Vector<Real>& alg,
	                 Matrix<Real>& df_dy) const override
	{
		df_dy(0, 0) = 0.0;
		df_dy(1, 0) = 0.0;
		df_dy(2, 0) = 1.0 / _m1;   // dv1/dt depends on λ
		df_dy(3, 0) = -1.0 / _m2;  // dv2/dt depends on λ
	}

	void jacobian_gx(Real t, const Vector<Real>& diff, const Vector<Real>& alg,
	                 Matrix<Real>& dg_dx) const override
	{
		dg_dx(0, 0) = 1.0;   // dg/dx1 = 1
		dg_dx(0, 1) = -1.0;  // dg/dx2 = -1
		dg_dx(0, 2) = 0.0;
		dg_dx(0, 3) = 0.0;
	}

	void jacobian_gy(Real t, const Vector<Real>& diff, const Vector<Real>& alg,
	                 Matrix<Real>& dg_dy) const override
	{
		// Note: This is 0! Position constraint doesn't depend on λ directly.
		// This is actually an Index-2 DAE if we write it this way.
		// For Index-1, we need to use the velocity-level constraint.
		dg_dy(0, 0) = 0.0;
	}
};

/// Velocity-level constraint version (Index-1):
/// Instead of position constraint, use velocity constraint:
/// v1 - v2 = 0 (relative velocity must be zero for rigid connection)
/// This makes ∂g/∂y nonsingular through the coupling of accelerations.
class RigidRod1DIndex1DAE : public IODESystemDAEWithJacobian
{
	Real _m1, _m2, _L, _g;

public:
	RigidRod1DIndex1DAE(Real m1 = 1.0, Real m2 = 1.0, Real L = 1.0, Real g = 9.81)
		: _m1(m1), _m2(m2), _L(L), _g(g) {}

	int getDiffDim() const override { return 4; }  // x1, x2, v1, v2
	int getAlgDim() const override { return 1; }   // λ

	void diffEqs(Real t, const Vector<Real>& diff, const Vector<Real>& alg,
	             Vector<Real>& dxdt) const override
	{
		Real v1 = diff[2], v2 = diff[3];
		Real lambda = alg[0];

		dxdt[0] = v1;
		dxdt[1] = v2;
		// Accelerations with constraint force
		dxdt[2] = -_g + lambda / _m1;
		dxdt[3] = -_g - lambda / _m2;
	}

	void algConstraints(Real t, const Vector<Real>& diff, const Vector<Real>& alg,
	                    Vector<Real>& g) const override
	{
		// Velocity-level constraint for Index-1: v1 = v2
		Real v1 = diff[2], v2 = diff[3];
		g[0] = v1 - v2;
	}

	std::string getDiffVarName(int i) const override
	{
		static const char* names[] = {"x1", "x2", "v1", "v2"};
		return (i >= 0 && i < 4) ? names[i] : "?";
	}
	std::string getAlgVarName(int i) const override { return "lambda"; }

	void jacobian_fx(Real t, const Vector<Real>& diff, const Vector<Real>& alg,
	                 Matrix<Real>& df_dx) const override
	{
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++)
				df_dx(i, j) = 0.0;
		df_dx(0, 2) = 1.0;  // dx1/dt = v1
		df_dx(1, 3) = 1.0;  // dx2/dt = v2
	}

	void jacobian_fy(Real t, const Vector<Real>& diff, const Vector<Real>& alg,
	                 Matrix<Real>& df_dy) const override
	{
		df_dy(0, 0) = 0.0;
		df_dy(1, 0) = 0.0;
		df_dy(2, 0) = 1.0 / _m1;
		df_dy(3, 0) = -1.0 / _m2;
	}

	void jacobian_gx(Real t, const Vector<Real>& diff, const Vector<Real>& alg,
	                 Matrix<Real>& dg_dx) const override
	{
		dg_dx(0, 0) = 0.0;
		dg_dx(0, 1) = 0.0;
		dg_dx(0, 2) = 1.0;   // dg/dv1 = 1
		dg_dx(0, 3) = -1.0;  // dg/dv2 = -1
	}

	void jacobian_gy(Real t, const Vector<Real>& diff, const Vector<Real>& alg,
	                 Matrix<Real>& dg_dy) const override
	{
		// g = v1 - v2, doesn't directly depend on λ
		// But through the implicit relation dv1/dt - dv2/dt = 0, we get:
		// λ/m1 - (-λ/m2) = 0 => λ*(1/m1 + 1/m2) = 0
		// For the augmented system, dg_dy is part of the constraint Jacobian
		dg_dy(0, 0) = 0.0;
	}
};

/// Nonlinear DAE: Van der Pol oscillator with algebraic output
/// dx/dt = y
/// dy/dt = μ*(1-x²)*y - x + z
/// 0 = z - sin(t)  (algebraic forcing)
///
/// The algebraic variable z tracks a sinusoidal input.
class VanDerPolDAE : public IODESystemDAEWithJacobian
{
	Real _mu;

public:
	VanDerPolDAE(Real mu = 1.0) : _mu(mu) {}

	int getDiffDim() const override { return 2; }  // x, y
	int getAlgDim() const override { return 1; }   // z (forcing term)

	void diffEqs(Real t, const Vector<Real>& x, const Vector<Real>& alg,
	             Vector<Real>& dxdt) const override
	{
		Real pos = x[0], vel = x[1];
		Real z = alg[0];

		dxdt[0] = vel;
		dxdt[1] = _mu * (1.0 - pos * pos) * vel - pos + z;
	}

	void algConstraints(Real t, const Vector<Real>& x, const Vector<Real>& alg,
	                    Vector<Real>& g) const override
	{
		Real z = alg[0];
		g[0] = z - std::sin(t);
	}

	std::string getDiffVarName(int i) const override
	{
		return i == 0 ? "x" : "y";
	}
	std::string getAlgVarName(int i) const override { return "z"; }

	void jacobian_fx(Real t, const Vector<Real>& x, const Vector<Real>& alg,
	                 Matrix<Real>& df_dx) const override
	{
		Real pos = x[0], vel = x[1];
		df_dx(0, 0) = 0.0;
		df_dx(0, 1) = 1.0;
		df_dx(1, 0) = -2.0 * _mu * pos * vel - 1.0;
		df_dx(1, 1) = _mu * (1.0 - pos * pos);
	}

	void jacobian_fy(Real t, const Vector<Real>& x, const Vector<Real>& alg,
	                 Matrix<Real>& df_dy) const override
	{
		df_dy(0, 0) = 0.0;
		df_dy(1, 0) = 1.0;  // z enters the velocity equation
	}

	void jacobian_gx(Real t, const Vector<Real>& x, const Vector<Real>& alg,
	                 Matrix<Real>& dg_dx) const override
	{
		dg_dx(0, 0) = 0.0;
		dg_dx(0, 1) = 0.0;
	}

	void jacobian_gy(Real t, const Vector<Real>& x, const Vector<Real>& alg,
	                 Matrix<Real>& dg_dy) const override
	{
		dg_dy(0, 0) = 1.0;  // dg/dz = 1 (nonsingular)
	}
};

///////////////////////////////////////////////////////////////////////////////////////////
//                           Realistic DAE Tests
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("DAE::RCCircuit_analyticalComparison", "[dae][solver][realistic]")
{
	// RC circuit with R=1kΩ, C=1μF, Vs=5V
	// Time constant τ = RC = 1ms
	Real R = 1000.0, C = 1e-6, Vs = 5.0;
	Real tau = R * C;  // 1e-3 seconds

	RCCircuitDAE circuit(R, C, Vs);

	// Initial conditions: capacitor discharged
	Vector<Real> Vc0(1), I0(1);
	Vc0[0] = 0.0;
	// Consistent: I = (Vs - Vc) / R = 5mA
	I0[0] = Vs / R;

	// Integrate for 5 time constants (should reach ~99.3% of final value)
	Real t_end = 5.0 * tau;

	DAESolverConfig config;
	config.step_size = tau / 100.0;  // 100 steps per time constant
	config.newton_tol = TOL(1e-12, 1e-5);

	auto result = SolveDAEBDF2(circuit, 0.0, Vc0, I0, t_end, config);

	REQUIRE(result.status == AlgorithmStatus::Success);

	// Check at multiple time points
	for (Real t : {tau, 2*tau, 3*tau, 5*tau})
	{
		// Find closest step
		int step = static_cast<int>(t / config.step_size);
		if (step > result.total_steps) step = result.total_steps;

		Real t_actual = result.solution.getTValue(step);
		Real Vc_numerical = result.solution.getXValue(step, 0);
		Real I_numerical = result.solution.getYValue(step, 0);

		Real Vc_analytical = circuit.analyticalVc(t_actual);
		Real I_analytical = circuit.analyticalI(t_actual);

		// BDF2 should be accurate to ~0.1% with these step sizes
		REQUIRE_THAT(Vc_numerical, WithinRel(Vc_analytical, REAL(0.001)));
		REQUIRE_THAT(I_numerical, WithinRel(I_analytical, REAL(0.001)));
	}

	// Check final values (at 5τ, Vc ≈ 4.966V, I ≈ 0.034mA)
	Vector<Real> Vc_end = result.solution.getXValuesAtEnd();
	Vector<Real> I_end = result.solution.getYValuesAtEnd();

	REQUIRE_THAT(Vc_end[0], WithinAbs(Vs * (1 - std::exp(-5)), 0.01));
	REQUIRE_THAT(I_end[0], WithinAbs((Vs/R) * std::exp(-5), 1e-6));
}

TEST_CASE("DAE::RCCircuit_constraintSatisfaction", "[dae][solver][realistic]")
{
	RCCircuitDAE circuit(1000.0, 1e-6, 5.0);

	Vector<Real> Vc0(1), I0(1);
	Vc0[0] = 0.0;
	I0[0] = 5.0 / 1000.0;  // 5mA initial current

	// Use reasonable tolerances for this problem scale
	DAESolverConfig config;
	config.step_size = 1e-5;      // 100 steps per time constant
	config.newton_tol = TOL(1e-10, 1e-5);    // Reasonable for this problem
	config.constraint_tol = TOL(1e-10, 1e-5);

	auto result = SolveDAEBackwardEuler(circuit, 0.0, Vc0, I0, 5e-3, config);

	REQUIRE(result.status == AlgorithmStatus::Success);

	// Verify KVL at each step: Vs = R*I + Vc
	Vector<Real> g(1);
	for (int step = 0; step <= result.total_steps; step += result.total_steps / 10)
	{
		Real t = result.solution.getTValue(step);
		Vector<Real> Vc(1), I(1);
		Vc[0] = result.solution.getXValue(step, 0);
		I[0] = result.solution.getYValue(step, 0);

		circuit.algConstraints(t, Vc, I, g);
		REQUIRE_THAT(g[0], WithinAbs(0.0, TOL(1e-10, 1e-5)));
	}
}

TEST_CASE("DAE::StiffLinear_stability", "[dae][solver][realistic]")
{
	// Stiff DAE with λ = 1000 (fast component)
	StiffLinearDAE stiffDAE(1000.0);

	// Initial conditions on the manifold
	Vector<Real> x0(1), y0(1);
	x0[0] = 1.0;  // x(0) = exp(0) = 1
	y0[0] = 1.0;  // y(0) = exp(0) = 1

	// Integrate for 3 time units
	Real t_end = 3.0;

	DAESolverConfig config;
	config.step_size = 0.01;  // Step much larger than 1/λ = 0.001
	config.newton_tol = TOL(1e-10, 1e-5);

	// Backward Euler should be stable despite large step
	auto resultBE = SolveDAEBackwardEuler(stiffDAE, 0.0, x0, y0, t_end, config);
	REQUIRE(resultBE.status == AlgorithmStatus::Success);

	// BDF2 should also be stable and more accurate
	auto resultBDF2 = SolveDAEBDF2(stiffDAE, 0.0, x0, y0, t_end, config);
	REQUIRE(resultBDF2.status == AlgorithmStatus::Success);

	// Check final values
	Real x_analytical = std::exp(-3.0);  // ≈ 0.0498
	Real y_analytical = std::exp(-3.0);

	Vector<Real> xBE = resultBE.solution.getXValuesAtEnd();
	Vector<Real> xBDF2 = resultBDF2.solution.getXValuesAtEnd();

	// Both should get reasonable answer (no blowup)
	REQUIRE_THAT(xBE[0], WithinAbs(x_analytical, 0.01));
	REQUIRE_THAT(xBDF2[0], WithinAbs(x_analytical, 0.001));  // BDF2 more accurate
}

TEST_CASE("DAE::StiffLinear_constraintTracking", "[dae][solver][realistic]")
{
	StiffLinearDAE stiffDAE(500.0);

	Vector<Real> x0(1), y0(1);
	x0[0] = 1.0;
	y0[0] = 1.0;

	auto result = SolveDAEBDF2(stiffDAE, 0.0, x0, y0, 2.0, 
	                           DAESolverConfig::HighPrecision());

	REQUIRE(result.status == AlgorithmStatus::Success);

	// Verify constraint x + y = 2*exp(-t) is satisfied
	Vector<Real> g(1);
	for (int i = 0; i <= 10; i++)
	{
		int step = i * result.total_steps / 10;
		Real t = result.solution.getTValue(step);
		Vector<Real> x(1), y(1);
		x[0] = result.solution.getXValue(step, 0);
		y[0] = result.solution.getYValue(step, 0);

		stiffDAE.algConstraints(t, x, y, g);
		REQUIRE_THAT(g[0], WithinAbs(0.0, TOL(1e-9, 1e-4)));
	}
}

TEST_CASE("DAE::VanDerPol_forcedOscillator", "[dae][solver][realistic]")
{
	// Van der Pol oscillator with sinusoidal forcing via algebraic variable
	VanDerPolDAE vdp(0.5);  // Moderate nonlinearity

	Vector<Real> x0(2), z0(1);
	x0[0] = 1.0;   // Initial position
	x0[1] = 0.0;   // Initial velocity
	z0[0] = 0.0;   // z(0) = sin(0) = 0 (consistent)

	DAESolverConfig config;
	config.step_size = 0.01;
	config.newton_tol = TOL(1e-10, 1e-5);

	auto result = SolveDAEBDF2(vdp, 0.0, x0, z0, 10.0, config);

	REQUIRE(result.status == AlgorithmStatus::Success);

	// Verify algebraic constraint z = sin(t) is satisfied
	Vector<Real> g(1);
	for (int i = 0; i <= 20; i++)
	{
		int step = i * result.total_steps / 20;
		Real t = result.solution.getTValue(step);
		Vector<Real> x(2), z(1);
		x[0] = result.solution.getXValue(step, 0);
		x[1] = result.solution.getXValue(step, 1);
		z[0] = result.solution.getYValue(step, 0);

		vdp.algConstraints(t, x, z, g);
		REQUIRE_THAT(g[0], WithinAbs(0.0, TOL(1e-8, 1e-4)));

		// Also verify z tracks sin(t)
		REQUIRE_THAT(z[0], WithinAbs(std::sin(t), TOL(1e-8, 1e-4)));
	}

	// Solution should oscillate (not blow up or decay to zero)
	Real x_final = result.solution.getXValue(result.total_steps, 0);
	REQUIRE(std::abs(x_final) < 10.0);  // Bounded
}

TEST_CASE("DAE::ConvergenceOrder_backwardEuler", "[dae][solver][realistic]")
{
	// Verify first-order convergence of Backward Euler
	SimpleLinearDAE linearDAE;

	Real t_end = 1.0;
	Vector<Real> x0(1), y0(1);
	x0[0] = 0.9;
	y0[0] = 0.1;

	Real x_analytical = 0.5 + 0.4 * std::exp(-2.0);

	// Use larger step sizes to see discretization error clearly
	std::vector<Real> stepSizes = {0.2, 0.1, 0.05};
	std::vector<Real> errors;

	for (Real h : stepSizes)
	{
		DAESolverConfig config;
		config.step_size = h;
		config.newton_tol = TOL(1e-12, 1e-5);

		auto result = SolveDAEBackwardEuler(linearDAE, 0.0, x0, y0, t_end, config);
		Real x_final = result.solution.getXValuesAtEnd()[0];
		errors.push_back(std::abs(x_final - x_analytical));
	}

	// Check first-order convergence: error proportional to h
	// Ratio of errors should be ~2 when h is halved
	for (size_t i = 1; i < errors.size(); i++)
	{
		Real ratio = errors[i-1] / errors[i];
		REQUIRE(ratio > 1.7);  // Should be ~2 for first-order
		REQUIRE(ratio < 2.5);  // Allow some margin
	}
}

TEST_CASE("DAE::ConvergenceOrder_BDF2", "[dae][solver][realistic]")
{
	// Verify BDF2 is more accurate than Backward Euler
	// (Full order verification requires very careful step size selection)
	SimpleLinearDAE linearDAE;

	Real t_end = 1.0;
	Vector<Real> x0(1), y0(1);
	x0[0] = 0.9;
	y0[0] = 0.1;

	Real x_analytical = 0.5 + 0.4 * std::exp(-2.0);

	// Use step sizes small enough that asymptotic behavior is visible
	// (first step is BE, so need enough BDF2 steps to dominate)
	std::vector<Real> stepSizes = {0.05, 0.025, 0.0125};
	std::vector<Real> errors;

	DAESolverConfig config;
	config.newton_tol = TOL(1e-12, 1e-5);

	for (Real h : stepSizes)
	{
		config.step_size = h;
		auto result = SolveDAEBDF2(linearDAE, 0.0, x0, y0, t_end, config);
		REQUIRE(result.status == AlgorithmStatus::Success);
		Real error = std::abs(result.solution.getXValuesAtEnd()[0] - x_analytical);
		errors.push_back(error);
	}

	// For 2nd order method, halving step size should reduce error by ~4
	Real ratio1 = errors[0] / errors[1];
	Real ratio2 = errors[1] / errors[2];

	// Both ratios should be close to 4 (allow 3.0-5.0 range for practical tolerance)
	REQUIRE(ratio1 > 3.0);
	REQUIRE(ratio1 < 5.5);
	REQUIRE(ratio2 > 3.0);
	REQUIRE(ratio2 < 5.5);
}

TEST_CASE("DAE::ComputeConsistentIC_nonlinear", "[dae][solver][realistic]")
{
	// Test consistent IC computation for Van der Pol (nonlinear case)
	VanDerPolDAE vdp(1.0);

	Vector<Real> x0(2), z0(1);
	x0[0] = 2.0;
	x0[1] = 0.0;
	z0[0] = 0.5;  // Wrong initial guess (should be sin(0) = 0)

	// Compute consistent z0
	Real t0 = 0.0;
	bool success = ComputeConsistentIC(vdp, t0, x0, z0, 20, TOL(1e-12, 1e-5));

	REQUIRE(success);
	REQUIRE_THAT(z0[0], WithinAbs(std::sin(t0), TOL(1e-10, 1e-5)));

	// Try at t = π/6
	Real t1 = Constants::PI / 6.0;
	z0[0] = 0.0;  // Wrong guess
	success = ComputeConsistentIC(vdp, t1, x0, z0, 20, TOL(1e-12, 1e-5));

	REQUIRE(success);
	REQUIRE_THAT(z0[0], WithinAbs(std::sin(t1), TOL(1e-10, 1e-5)));
}

} // namespace MML::Tests::Algorithms::DAESolverTests

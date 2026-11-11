#include <catch2/catch_all.hpp>
#include "../../TestPrecision.h"

#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/core/LinAlgEqSolvers/LinAlgEqSolvers_iterative.h>
#include <mml/systems/LinearSystem.h>
#endif

#include <array>
#include <limits>

using namespace MML;

namespace MML::Tests::DenseIterativeContracts
{
	enum class Method { Jacobi, GaussSeidel, SOR };
	constexpr std::array methods{Method::Jacobi, Method::GaussSeidel, Method::SOR};

	IterativeSolverResult Run(Method method, bool configured, const Matrix<Real>& a, const Vector<Real>& b,
		const Vector<Real>& guess = {}, Real tolerance = Precision::IterativeSolverTolerance, int budget = 1)
	{
		IterativeSolverConfig config(tolerance, budget);
		config.omega = Real(1.25);
		switch (method) {
		case Method::Jacobi:
			return configured ? JacobiSolver::Solve(a, b, config, guess) : JacobiSolver::Solve(a, b, guess, tolerance, budget);
		case Method::GaussSeidel:
			return configured ? GaussSeidelSolver::Solve(a, b, config, guess) : GaussSeidelSolver::Solve(a, b, guess, tolerance, budget);
		default:
			return configured ? SORSolver::Solve(a, b, config, guess) : SORSolver::Solve(a, b, config.omega, guess, tolerance, budget);
		}
	}

	Real Residual2(const Matrix<Real>& a, const Vector<Real>& b, const Vector<Real>& x)
	{
		return std::hypot(a(0, 0) * x[0] + a(0, 1) * x[1] - b[0],
			a(1, 0) * x[0] + a(1, 1) * x[1] - b[1]) / std::hypot(b[0], b[1]);
	}

	TEST_CASE("Dense iteration uncomputed results never claim success or zero residual", "[dense-iteration-contract][Iterative]")
	{
		IterativeSolverResult empty;
		CHECK(empty.status == AlgorithmStatus::InvalidInput);
		CHECK_FALSE(empty.converged);
		CHECK_FALSE(HasEstimate(empty.residual));
		CHECK_FALSE(HasEstimate(empty.initial_residual));
		CHECK_FALSE(empty.error_message.empty());
		IterativeSolverResult exhausted(Vector<Real>{1}, 1, 2, false);
		CHECK(exhausted.status == AlgorithmStatus::MaxIterationsExceeded);
		CHECK_FALSE(exhausted.error_message.empty());
	}

	TEST_CASE("Dense iteration zero budget evaluates the actual initial guess", "[dense-iteration-contract][Iterative]")
	{
		const Matrix<Real> a(2, 2, {4, 1, 1, 3});
		const Vector<Real> b{1, 2};
		for (auto method : methods)
			for (bool configured : {false, true}) {
				auto result = Run(method, configured, a, b, {}, Precision::IterativeSolverTolerance, 0);
				CHECK_FALSE(result.converged);
				CHECK(result.status == AlgorithmStatus::MaxIterationsExceeded);
				CHECK(result.iterations == 0);
				CHECK(result.initial_residual == 1);
				CHECK(result.residual == 1);
				CHECK_FALSE(result.error_message.empty());
				CHECK_FALSE(result.algorithm_name.empty());
				CHECK(std::isfinite(result.elapsed_time_ms));
			}
	}

	TEST_CASE("Dense iteration checks convergence before spending any sweep budget", "[dense-iteration-contract][Iterative]")
	{
		const Matrix<Real> a(2, 2, {4, 1, 1, 3});
		const Vector<Real> exact{1, 2}, b{6, 7};
		for (auto method : methods)
			for (bool configured : {false, true})
				for (int budget : {0, 5}) {
					const auto result = Run(method, configured, a, b, exact, Precision::IterativeSolverTolerance, budget);
					CHECK(result.converged);
					CHECK(result.status == AlgorithmStatus::Success);
					CHECK(result.iterations == 0);
					CHECK(result.initial_residual == 0);
					CHECK(result.residual == 0);
					CHECK(result.error_message.empty());
					CHECK(result.solution == exact);
				}
	}

	TEST_CASE("Dense iteration one sweep reports exact work and matching normalized residuals", "[dense-iteration-contract][Iterative]")
	{
		const Matrix<Real> a(2, 2, {4, 1, 1, 3});
		const Vector<Real> b{1, 2};
		for (auto method : methods)
			for (bool configured : {false, true}) {
				const auto result = Run(method, configured, a, b);
				const Vector<Real> expected = method == Method::Jacobi ? Vector<Real>{Real(.25), Real(2) / 3}
					: method == Method::GaussSeidel ? Vector<Real>{Real(.25), Real(7) / 12}
					: Vector<Real>{Real(.3125), Real(.703125)};
				CHECK(result.iterations == 1);
				CHECK_FALSE(result.converged);
				CHECK(result.status == AlgorithmStatus::MaxIterationsExceeded);
				CHECK(result.initial_residual == 1);
				CHECK(result.solution.IsEqualTo(expected, Precision::EqualityTolerance));
				CHECK(result.residual == Catch::Approx(Residual2(a, b, expected)).epsilon(Precision::EqualityTolerance));
				CHECK_FALSE(result.algorithm_name.empty());
				CHECK_FALSE(result.error_message.empty());
				CHECK(result.omega_used == (method == Method::SOR ? Real(1.25) : Real(0)));
			}
	}

	TEST_CASE("Dense iteration nonzero guesses use the same initial and final residual units", "[dense-iteration-contract][Iterative]")
	{
		const Matrix<Real> a(2, 2, {4, 1, 1, 3});
		const Vector<Real> b{1, 2}, guess{Real(.125), Real(.25)};
		for (auto method : methods)
			for (bool configured : {false, true}) {
				const auto result = Run(method, configured, a, b, guess);
				CHECK(result.initial_residual == Catch::Approx(Residual2(a, b, guess)).epsilon(Precision::EqualityTolerance));
				CHECK(result.residual == Catch::Approx(Residual2(a, b, result.solution)).epsilon(Precision::EqualityTolerance));
			}
	}

	TEST_CASE("Dense iteration zero RHS uses absolute residuals and can succeed without sweeps", "[dense-iteration-contract][Iterative]")
	{
		const Matrix<Real> a(2, 2, {1, 0, 0, 1});
		const Vector<Real> b{0, 0}, guess{3, 4};
		for (auto method : methods)
			for (bool configured : {false, true}) {
				const auto result = Run(method, configured, a, b, guess, Precision::IterativeSolverTolerance, 0);
				CHECK(result.initial_residual == 5);
				CHECK(result.residual == 5);
				CHECK(result.iterations == 0);
				const auto zero = Run(method, configured, a, b);
				CHECK(zero.converged);
				CHECK(zero.iterations == 0);
			}
	}

	TEST_CASE("Dense Jacobi SPD is not a convergence guarantee and retains its best finite iterate", "[dense-iteration-contract][Iterative]")
	{
		// Eigenvalues are 0.1, 0.1, 2.8; the Jacobi error factor on (1,1,1) is -1.8.
		const Matrix<Real> a(3, 3, {1, Real(.9), Real(.9), Real(.9), 1, Real(.9), Real(.9), Real(.9), 1});
		const Vector<Real> b{1, 1, 1};
		for (bool configured : {false, true}) {
			const auto result = Run(Method::Jacobi, configured, a, b, {}, Precision::IterativeSolverTolerance, 5);
			CHECK_FALSE(result.converged);
			CHECK(result.status == AlgorithmStatus::MaxIterationsExceeded);
			CHECK(result.iterations == 5);
			CHECK(result.initial_residual == 1);
			CHECK(result.residual == 1);
			CHECK(result.solution == Vector<Real>{0, 0, 0});
		}
	}

	TEST_CASE("Dense iteration controls and nonempty guesses are validated uniformly", "[dense-iteration-contract][Iterative]")
	{
		const Matrix<Real> a(2, 2, {4, 1, 1, 3});
		const Vector<Real> b{1, 2};
		for (auto method : methods)
			for (bool configured : {false, true}) {
				for (Real tolerance : {Real(0), Real(-1), std::numeric_limits<Real>::infinity(), NoEstimate})
					CHECK_THROWS_AS(Run(method, configured, a, b, {}, tolerance), ArgumentError);
				CHECK_THROWS_AS(Run(method, configured, a, b, {}, Precision::IterativeSolverTolerance, -1), ArgumentError);
				CHECK_THROWS_AS(Run(method, configured, a, b, Vector<Real>{1}), MatrixDimensionError);
				CHECK_THROWS_AS(Run(method, configured, Matrix<Real>(), Vector<Real>()), MatrixDimensionError);
			}
	}

	TEST_CASE("Dense iteration rejects nonfinite matrix RHS and initial-guess entries", "[dense-iteration-contract][Iterative]")
	{
		for (auto method : methods)
			for (bool configured : {false, true})
				for (Real bad : {NoEstimate, std::numeric_limits<Real>::infinity()}) {
					Matrix<Real> a(2, 2, {4, 1, 1, 3});
					Vector<Real> b{1, 2}, guess{0, 0};
					a(0, 1) = bad;
					CHECK_THROWS_AS(Run(method, configured, a, b), ArgumentError);
					a(0, 1) = 1;
					b[0] = bad;
					CHECK_THROWS_AS(Run(method, configured, a, b), ArgumentError);
					b[0] = 1;
					guess[1] = bad;
					CHECK_THROWS_AS(Run(method, configured, a, b, guess), ArgumentError);
				}
	}

	TEST_CASE("Dense iteration numerical breakdown is not overwritten by config wrappers", "[dense-iteration-contract][Iterative]")
	{
		const Real large = std::numeric_limits<Real>::max() / 8;
		const Matrix<Real> a(2, 2, {1, large, large, 1});
		const Vector<Real> b{1, 1};
		for (auto method : methods)
			for (bool configured : {false, true}) {
				const auto result = Run(method, configured, a, b, {}, Precision::IterativeSolverTolerance, 10);
				CHECK_FALSE(result.converged);
				CHECK(result.status == AlgorithmStatus::NumericalInstability);
				CHECK(result.iterations >= 1);
				CHECK(result.iterations <= 2);
				CHECK(result.residual == 1);
				CHECK(result.solution == Vector<Real>{0, 0});
				CHECK_FALSE(result.error_message.empty());
			}
	}

	TEST_CASE("Dense iteration unrepresentable initial residuals do not fabricate estimates", "[dense-iteration-contract][Iterative]")
	{
		for (auto method : methods)
			for (bool configured : {false, true}) {
				const Real large = std::numeric_limits<Real>::max();
				const Matrix<Real> identity(2, 2, {1, 0, 0, 1});
				const auto rhsOverflow = Run(method, configured, identity, Vector<Real>{large, large});
				CHECK(rhsOverflow.status == AlgorithmStatus::NumericalInstability);
				CHECK(rhsOverflow.iterations == 0);
				CHECK_FALSE(HasEstimate(rhsOverflow.initial_residual));
				CHECK_FALSE(HasEstimate(rhsOverflow.residual));
				const Matrix<Real> a(2, 2, {1, large, large, 1});
				const auto residualOverflow = Run(method, configured, a, Vector<Real>{1, 1}, Vector<Real>{large, large});
				CHECK(residualOverflow.status == AlgorithmStatus::NumericalInstability);
				CHECK_FALSE(HasEstimate(residualOverflow.residual));
				CHECK_FALSE(residualOverflow.error_message.empty());
			}
	}

	TEST_CASE("Dense iteration residuals are invariant under representable common system scaling", "[dense-iteration-contract][Iterative]")
	{
		const Real small = std::sqrt(std::numeric_limits<Real>::min());
		const Real large = std::sqrt(std::numeric_limits<Real>::max()) / 16;
		for (auto method : methods)
			for (bool configured : {false, true})
				for (Real scale : {small, Real(1), large}) {
					const Matrix<Real> a(2, 2, {4 * scale, scale, scale, 3 * scale});
					const Vector<Real> b{scale, 2 * scale};
					const auto result = Run(method, configured, a, b);
					CHECK(result.initial_residual == 1);
					CHECK(result.residual == Catch::Approx(Residual2(a, b, result.solution)).epsilon(Precision::EqualityTolerance));
					CHECK(result.iterations == 1);
				}
	}

	TEST_CASE("Dense SOR convenience entry points do not hide invalid controls or failed trials", "[dense-iteration-contract][Iterative]")
	{
		const Matrix<Real> a(1, 1, {2});
		const Vector<Real> b{4};
		for (Real omega : {NoEstimate, std::numeric_limits<Real>::infinity(), Real(0), Real(2)})
			CHECK_THROWS_AS(SORSolver::Solve(a, b, omega), NumericalMethodError);
		CHECK(SORSolver::SolveOptimal(a, b).converged);
		for (int samples : {0, -1})
			CHECK_THROWS_AS(SORSolver::EstimateOptimalOmega(a, b, samples), ArgumentError);
		for (int budget : {0, -1})
			CHECK_THROWS_AS(SORSolver::EstimateOptimalOmega(a, b, 2, budget), ArgumentError);
		const Real large = std::numeric_limits<Real>::max() / 8;
		const Matrix<Real> unstable(2, 2, {1, large, large, 1});
		CHECK_THROWS_AS(SORSolver::EstimateOptimalOmega(unstable, Vector<Real>{1, 1}), NumericalMethodError);
	}

	TEST_CASE("Dense iteration inputs stay unchanged and returned solutions own their storage", "[dense-iteration-contract][Iterative]")
	{
		for (auto method : methods)
			for (bool configured : {false, true}) {
				const auto solve = [&] {
					const Matrix<Real> a(2, 2, {4, 1, 1, 3});
					const Vector<Real> both{1, 2};
					auto result = Run(method, configured, a, both, both);
					CHECK(both == Vector<Real>{1, 2});
					CHECK(a(0, 0) == 4);
					return result;
				};
				const auto result = solve();
				CHECK(result.solution.size() == 2);
				CHECK(std::isfinite(result.solution[0]));
			}
	}

	TEST_CASE("Dense simple solvers propagate truthful zero-budget exhaustion", "[dense-iteration-contract][Iterative]")
	{
		const Matrix<Real> a(1, 1, {1});
		const Vector<Real> b{2};
		CHECK_THROWS_WITH(JacobiSolver::SolveSimple(a, b, Precision::IterativeSolverTolerance, 0),
			Catch::Matchers::ContainsSubstring("0 iterations"));
		CHECK_THROWS_WITH(GaussSeidelSolver::SolveSimple(a, b, Precision::IterativeSolverTolerance, 0),
			Catch::Matchers::ContainsSubstring("0 iterations"));
		CHECK_THROWS_WITH(SORSolver::SolveSimple(a, b, 1, Precision::IterativeSolverTolerance, 0),
			Catch::Matchers::ContainsSubstring("0 iterations"));
	}

	TEST_CASE("LinearSystem consumers receive truthful dense iteration failures", "[dense-iteration-contract][LinearSystem]")
	{
		using namespace MML::Systems;
		const Matrix<Real> a(2, 2, {4, 1, 1, 3});
		const LinearSystem<Real> system(a, Vector<Real>{1, 2});
		for (auto method : {IterativeMethod::Jacobi, IterativeMethod::GaussSeidel, IterativeMethod::SOR}) {
			try {
				const auto solution = system.SolveIterative(method, Precision::IterativeSolverTolerance, 0);
				FAIL("A nonconverged zero-budget solve must throw");
			} catch (const ConvergenceError& error) {
				CHECK(error.iterations() == 0);
				CHECK(error.residual() == 1);
			}
			CHECK_THROWS_AS(system.SolveIterative(method, NoEstimate), ArgumentError);
		}
	}
}

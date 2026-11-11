///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DAEBDF4.h                                                           ///
///  Description: BDF4 method for semi-explicit index-1 DAE systems                   ///
///                                                                                   ///
///  BDF4 formula for differential equations:                                         ///
///    x_{n+1} = (48/25)x_n - (36/25)x_{n-1} + (16/25)x_{n-2} - (3/25)x_{n-3}         ///
///              + (12/25)h·f(t_{n+1}, x_{n+1}, y_{n+1})                               ///
///    0 = g(t_{n+1}, x_{n+1}, y_{n+1})                                               ///
///                                                                                   ///
///  Properties:                                                                      ///
///  - 4th order accurate                                                             ///
///  - A(α)-stable with α ≈ 73°                                                       ///
///  - Industry standard for stiff DAEs (DASSL, SUNDIALS IDA)                         ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DAE_BDF4_H
#define MML_DAE_BDF4_H

#include "DAESolverBase.h"
#include "DAERadauIIA.h"
#include <algorithm>
#include <vector>

namespace MML {

	/// @brief Fourth-order Backward Differentiation Formula solver for index-1 DAEs
	/// 
	/// Complexity: O(N³) per step where N = diffDim + algDim.
	///            Newton iteration with augmented Jacobian + LU solve each step.
	///            First 3 steps use fifth-order Radau IIA for accurate history.
	///
	/// BDF4 formula for differential equations:
	/// @f[
	///   x_{n+1} = \frac{48}{25}x_n - \frac{36}{25}x_{n-1} + \frac{16}{25}x_{n-2} - \frac{3}{25}x_{n-3} 
	///             + \frac{12}{25}h \cdot f(t_{n+1}, x_{n+1}, y_{n+1})
	/// @f]
	/// @f[
	///   0 = g(t_{n+1}, x_{n+1}, y_{n+1})
	/// @f]
	///
	/// Uses Radau IIA for the first 3 intervals, then BDF4 on the uniform grid.
	/// A clipped final interval also uses Radau IIA: a BE restart would introduce
	/// O(h^2) endpoint error and destroy fourth-order global accuracy.
	/// Radau's coupled stage solves are included in the Newton diagnostics, but
	/// each completed interval counts as one accepted step.
	/// Fourth-order accurate, A(α)-stable with α ≈ 73°.
	/// Industry standard for stiff DAEs (DASSL, SUNDIALS IDA).
	///
	/// @param system DAE system with Jacobian
	/// @param t0 Initial time
	/// @param x0 Initial differential state
	/// @param y0 Initial algebraic state (must satisfy constraints)
	/// @param t_end Final time
	/// @param config Solver configuration
	/// @return DAESolverResult with solution and diagnostics
	inline DAESolverResult SolveDAEBDF4(IODESystemDAEWithJacobian& system,
	                                     Real t0, const Vector<Real>& x0, const Vector<Real>& y0,
	                                     Real t_end, const DAESolverConfig& config = DAESolverConfig())
	{
		const int diffDim = system.getDiffDim(), algDim = system.getAlgDim();
		if (!std::isfinite(t0) || !std::isfinite(t_end) || !std::isfinite(t_end - t0) || t_end == t0 ||
			!std::isfinite(config.step_size) || config.step_size <= 0 ||
			!std::isfinite(config.newton_tol) || config.newton_tol <= 0 ||
			!std::isfinite(config.constraint_tol) || config.constraint_tol <= 0 ||
			config.max_steps <= 0 || config.max_newton_iter <= 0 ||
			diffDim <= 0 || algDim < 0 || x0.size() != diffDim || y0.size() != algDim)
			throw ArgumentError("SolveDAEBDF4: requires distinct finite endpoints, positive finite "
				"step and tolerances, positive budgets and matching state dimensions");

		AlgorithmTimer timer;

		DAESolverResult result(t0, t_end, diffDim, algDim, 0);
		result.algorithm_name = "DAEBDF4";

		// Working vectors
		Vector<Real> dxdt(diffDim);
		Vector<Real> g(algDim);
		Matrix<Real> df_dx(diffDim, diffDim);
		Matrix<Real> df_dy(diffDim, algDim);
		Matrix<Real> dg_dx(algDim, diffDim);
		Matrix<Real> dg_dy(algDim, algDim);

		// BDF4: x_new = (48/25)*x_n - (36/25)*x_{n-1} + (16/25)*x_{n-2} - (3/25)*x_{n-3} + (12/25)*h*f
		constexpr Real bdf4_c0 = MML::Real(48.0) / MML::Real(25.0);
		constexpr Real bdf4_c1 = -MML::Real(36.0) / MML::Real(25.0);
		constexpr Real bdf4_c2 = MML::Real(16.0) / MML::Real(25.0);
		constexpr Real bdf4_c3 = -MML::Real(3.0) / MML::Real(25.0);
		constexpr Real bdf4_h_coeff = MML::Real(12.0) / MML::Real(25.0);

		auto newtonSolve = [&](Real t_next, Real h, Vector<Real>& x_new, Vector<Real>& y_new,
		                       const Vector<Real>& xn, const Vector<Real>& xn1,
		                       const Vector<Real>& xn2, const Vector<Real>& xn3)
		{
			const Real derivativeCoefficient = bdf4_h_coeff * h;

			auto evaluateResidual = [&](const Vector<Real>& trialX, const Vector<Real>& trialY, Vector<Real>& residual) {
				system.diffEqs(t_next, trialX, trialY, dxdt);
				system.algConstraints(t_next, trialX, trialY, g);
				for (int i = 0; i < diffDim; ++i)
					residual[i] = trialX[i] - bdf4_c0 * xn[i] - bdf4_c1 * xn1[i]
						- bdf4_c2 * xn2[i] - bdf4_c3 * xn3[i] - derivativeCoefficient * dxdt[i];
				for (int i = 0; i < algDim; ++i)
					residual[diffDim + i] = g[i];
			};
			auto evaluateJacobian = [&](const Vector<Real>& trialX, const Vector<Real>& trialY, Matrix<Real>& jacobian) {
				system.allJacobians(t_next, trialX, trialY, df_dx, df_dy, dg_dx, dg_dy);
				AssembleDAEAugmentedJacobian(df_dx, df_dy, dg_dx, dg_dy, derivativeCoefficient, jacobian);
			};
			return SolveDAENewton(x_new, y_new, diffDim, algDim, config, evaluateResidual, evaluateJacobian);
		};

		const Real direction = t_end > t0 ? Real(1) : Real(-1);
		const Real h_main = direction * config.step_size;

		// We need values at t0, t0+h, t0+2h, t0+3h for BDF4
		// Store these in history array: hist[0]=newest, hist[3]=oldest
		std::vector<Vector<Real>> x_history(4, x0);
		std::vector<Vector<Real>> y_history(4, y0);
		Real t_current = t0;
		result.solution.fillValues(0, t0, x0, y0);
		if (!ValidateDAEInitialState(system, t0, x0, y0, config, result)) {
			result.solution.setFinalSize(0);
			result.elapsed_time_ms = timer.elapsed_ms();
			return result;
		}

		int step = 1;
		const auto rejectAttempt = [&]() {
			++result.rejected_steps;
			result.solution.incrementRejectedSteps();
		};
		while (direction * (t_end - t_current) > 0 && step <= config.max_steps &&
			result.status == AlgorithmStatus::Success)
		{
			if (!PrepareDAETrajectoryPoint(result, step)) break;
			const Real remaining = t_end - t_current;
			Real h = direction * std::min(std::abs(h_main), std::abs(remaining));
			// Snap grid-roundoff at the endpoint without discarding genuine short intervals.
			if (std::abs(remaining - h_main) <= std::abs(h_main) * Precision::MinimumStepRatio)
				h = remaining;
			const Real t_next = h == remaining ? t_end : t_current + h;
			if (!ValidateDAETimeStep(t_current, t_next, result, "BDF4")) break;

			Vector<Real> x_new = x_history[0];
			Vector<Real> y_new = y_history[0];
			if (step <= 3 || std::abs(h - h_main) > DAETimeTolerance(t_end, h_main)) {
				// High-order starting values preserve BDF4 convergence (Hairer-Wanner, Solving ODEs II).
				// Reuse the L-stable collocation solver, including its algebraic stage constraints.
				DAESolverConfig starterConfig = config;
				starterConfig.step_size = std::abs(h);
				starterConfig.max_steps = 1;
				auto starter = SolveDAERadauIIA(system, t_current, x_new, y_new, t_next, starterConfig);
				DAENewtonResult diagnostics;
				diagnostics.iterations = starter.newton_iterations;
				diagnostics.jacobian_evaluations = starter.jacobian_evaluations;
				diagnostics.lu_factorizations = starter.lu_factorizations;
				diagnostics.linear_solves = starter.linear_solves;
				diagnostics.residual_norm = starter.final_residual_norm;
				AccumulateDAENewtonDiagnostics(result, diagnostics);
				result.max_constraint_violation = std::max(result.max_constraint_violation, starter.max_constraint_violation);
				if (starter.status != AlgorithmStatus::Success) {
					result.status = starter.status;
					result.failure_reason = starter.failure_reason;
					result.error_message = "BDF4 Radau history/endpoint step: " + starter.error_message;
					rejectAttempt();
					break;
				}
				if (starter.accepted_steps != 1 || starter.solution.getTValue(starter.solution.size() - 1) != t_next) {
					result.status = AlgorithmStatus::Stalled;
					result.failure_reason = DAEFailureReason::StepSizeUnderflow;
					result.error_message = "BDF4 Radau history/endpoint step did not reach the requested time";
					rejectAttempt();
					break;
				}
				x_new = starter.solution.getXValuesAtEnd();
				y_new = starter.solution.getYValuesAtEnd();
			} else {
				system.diffEqs(t_current, x_history[0], y_history[0], dxdt);
				x_new += dxdt * h;
				auto newton = newtonSolve(t_next, h, x_new, y_new,
					x_history[0], x_history[1], x_history[2], x_history[3]);
				AccumulateDAENewtonDiagnostics(result, newton);
				if (!newton.converged) {
					SetDAENewtonFailure(result, newton, "BDF4 Newton failure at t=" + std::to_string(t_next));
					rejectAttempt();
					break;
				}
			}
			system.algConstraints(t_next, x_new, y_new, g);
			const Real previousConstraintNorm = result.final_constraint_norm;
			if (!ValidateAcceptedDAEState(x_new, y_new, g, config, result, "BDF4 step")) {
				result.final_constraint_norm = previousConstraintNorm;
				rejectAttempt();
				break;
			}

			// Shift history
			for (int i = 3; i > 0; --i)
			{
				x_history[i] = x_history[i-1];
				y_history[i] = y_history[i-1];
			}
			x_history[0] = x_new;
			y_history[0] = y_new;
			t_current = t_next;

			result.solution.fillValues(step, t_current, x_history[0], y_history[0]);
			result.solution.incrementSuccessfulSteps();
			result.accepted_steps++;
			++step;
		}

		// Finalize
		if (result.status != AlgorithmStatus::ResourceLimitExceeded)
			result.solution.setFinalSize(step - 1);
		result.total_steps = step - 1;
		if (direction * (t_end - t_current) > 0 && result.status == AlgorithmStatus::Success) {
			result.status = AlgorithmStatus::MaxIterationsExceeded;
			result.failure_reason = DAEFailureReason::MaxStepsExceeded;
			result.error_message = "Maximum BDF4 step count reached before t_end";
		}
		result.elapsed_time_ms = timer.elapsed_ms();

		return result;
	}

} // namespace MML

#endif // MML_DAE_BDF4_H

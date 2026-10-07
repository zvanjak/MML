///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DAEBDF2.h                                                           ///
///  Description: BDF2 method for semi-explicit index-1 DAE systems                   ///
///                                                                                   ///
///  BDF2 formula for differential equations:                                         ///
///    x_{n+1} = (4/3)x_n - (1/3)x_{n-1} + (2/3)h·f(t_{n+1}, x_{n+1}, y_{n+1})        ///
///    0 = g(t_{n+1}, x_{n+1}, y_{n+1})                                               ///
///                                                                                   ///
///  Properties:                                                                      ///
///  - 2nd order accurate                                                             ///
///  - A-stable                                                                       ///
///  - Industry standard for moderately stiff DAEs                                    ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DAE_BDF2_H
#define MML_DAE_BDF2_H

#include "DAESolverBase.h"

namespace MML {

	/// @brief Second-order Backward Differentiation Formula solver for index-1 DAEs
	/// 
	/// Complexity: O(N³) per step where N = diffDim + algDim.
	///            Newton iteration with augmented Jacobian + LU solve each step.
	///
	/// BDF2 formula for differential equations:
	/// @f[
	///   x_{n+1} = \frac{4}{3}x_n - \frac{1}{3}x_{n-1} + \frac{2}{3}h \cdot f(t_{n+1}, x_{n+1}, y_{n+1})
	/// @f]
	/// @f[
	///   0 = g(t_{n+1}, x_{n+1}, y_{n+1})
	/// @f]
	///
	/// Uses Backward Euler for the first step to bootstrap, then switches to BDF2.
	/// Second-order accurate, A-stable, industry standard for moderately stiff DAEs.
	///
	/// @param system DAE system with Jacobian
	/// @param t0 Initial time
	/// @param x0 Initial differential state
	/// @param y0 Initial algebraic state (must satisfy constraints)
	/// @param t_end Final time
	/// @param config Solver configuration
	/// @return DAESolverResult with solution and diagnostics
	inline DAESolverResult SolveDAEBDF2(IODESystemDAEWithJacobian& system,
	                                     Real t0, const Vector<Real>& x0, const Vector<Real>& y0,
	                                     Real t_end, const DAESolverConfig& config = DAESolverConfig())
	{
		if (t_end <= t0)
			throw ArgumentError("SolveDAEBDF2: t_end must be greater than t0 (reverse-time integration is not supported)");
		if (config.step_size <= 0)
			throw ArgumentError("SolveDAEBDF2: config.step_size must be positive");

		AlgorithmTimer timer;

		int diffDim = system.getDiffDim();
		int algDim = system.getAlgDim();
		int num_steps = static_cast<int>((t_end - t0) / config.step_size) + 1;

		DAESolverResult result(t0, t_end, diffDim, algDim, num_steps);
		result.algorithm_name = "DAEBDF2";

		// Current and previous state
		Real t = t0;
		Vector<Real> x = x0;
		Vector<Real> y = y0;
		Vector<Real> x_prev = x0;  // For BDF2, need previous step

		// Working vectors
		Vector<Real> dxdt(diffDim);
		Vector<Real> g(algDim);
		Matrix<Real> df_dx(diffDim, diffDim);
		Matrix<Real> df_dy(diffDim, algDim);
		Matrix<Real> dg_dx(algDim, diffDim);
		Matrix<Real> dg_dy(algDim, algDim);

		// Save initial condition
		result.solution.fillValues(0, t, x, y);
		if (!ValidateDAEInitialState(system, t, x, y, config, result)) {
			result.solution.setFinalSize(0);
			result.elapsed_time_ms = timer.elapsed_ms();
			return result;
		}
		int step = 1;

		auto newtonSolve = [&](Real t_next, Real h, Vector<Real>& x_new, Vector<Real>& y_new,
		                       bool useBDF2, const Vector<Real>& x_prev_bdf)
		{
			const Real derivativeCoefficient = useBDF2 ? REAL(2.0 / 3.0) * h : h;
			auto evaluateResidual = [&](const Vector<Real>& trialX, const Vector<Real>& trialY, Vector<Real>& residual) {
				system.diffEqs(t_next, trialX, trialY, dxdt);
				system.algConstraints(t_next, trialX, trialY, g);
				if (useBDF2)
				{
					for (int i = 0; i < diffDim; ++i)
						residual[i] = trialX[i] - REAL(4.0 / 3.0) * x[i]
							+ REAL(1.0 / 3.0) * x_prev_bdf[i] - derivativeCoefficient * dxdt[i];
				}
				else
				{
					for (int i = 0; i < diffDim; ++i)
						residual[i] = trialX[i] - x[i] - h * dxdt[i];
				}
				for (int i = 0; i < algDim; ++i)
					residual[diffDim + i] = g[i];
			};
			auto evaluateJacobian = [&](const Vector<Real>& trialX, const Vector<Real>& trialY, Matrix<Real>& jacobian) {
				system.allJacobians(t_next, trialX, trialY, df_dx, df_dy, dg_dx, dg_dy);
				AssembleDAEAugmentedJacobian(df_dx, df_dy, dg_dx, dg_dy, derivativeCoefficient, jacobian);
			};
			return SolveDAENewton(x_new, y_new, diffDim, algDim, config, evaluateResidual, evaluateJacobian);
		};

		// First step: Backward Euler
		if (step < config.max_steps && t < t_end)
		{
			Real h = std::min(config.step_size, t_end - t);
			Real t_next = t + h;

			// Initial guess
			system.diffEqs(t, x, y, dxdt);
			Vector<Real> x_new = x + dxdt * h;
			Vector<Real> y_new = y;

			DAENewtonResult newton = newtonSolve(t_next, h, x_new, y_new, false, x_prev);
			AccumulateDAENewtonDiagnostics(result, newton);
			if (!newton.converged)
			{
				SetDAENewtonFailure(result, newton, "BDF2 bootstrap Newton failure");
				result.elapsed_time_ms = timer.elapsed_ms();
				return result;
			}

			// Accept step
			x_prev = x;
			x = x_new;
			y = y_new;
			t = t_next;

			system.algConstraints(t, x, y, g);
			if (!ValidateAcceptedDAEState(x, y, g, config, result, "BDF2 bootstrap step"))
				return result;

			result.solution.fillValues(step, t, x, y);
			result.solution.incrementSuccessfulSteps();
			result.accepted_steps++;
			++step;
		}

		// Remaining steps: BDF2
		while (t < t_end && step < config.max_steps)
		{
			Real h = std::min(config.step_size, t_end - t);
			
			// Skip degenerate final steps (floating-point round-off)
			if (h < config.step_size * 1e-10)
				break;
			
			Real t_next = t + h;

			// For very small final steps (< 10% of nominal), use Backward Euler
			bool use_be_for_final = (h < config.step_size * 0.1);

			// BDF2 prediction: extrapolate
			system.diffEqs(t, x, y, dxdt);
			Vector<Real> x_new = x + dxdt * h;
			Vector<Real> y_new = y;

			DAENewtonResult newton = newtonSolve(t_next, h, x_new, y_new, !use_be_for_final, x_prev);
			AccumulateDAENewtonDiagnostics(result, newton);
			if (!newton.converged)
			{
				SetDAENewtonFailure(result, newton, "BDF2 Newton failure at t=" + std::to_string(t_next));
				break;
			}

			// Accept step
			x_prev = x;
			x = x_new;
			y = y_new;
			t = t_next;

			system.algConstraints(t, x, y, g);
			if (!ValidateAcceptedDAEState(x, y, g, config, result, "BDF2 step"))
				break;

			result.solution.fillValues(step, t, x, y);
			result.solution.incrementSuccessfulSteps();
			result.accepted_steps++;
			++step;
		}

		// Finalize
		result.solution.setFinalSize(step - 1);
		result.total_steps = step - 1;
		result.final_constraint_norm = g.NormL2();
		if (t_end - t > DAETimeTolerance(t_end, config.step_size) && result.status == AlgorithmStatus::Success) {
			result.status = AlgorithmStatus::MaxIterationsExceeded;
			result.failure_reason = DAEFailureReason::MaxStepsExceeded;
			result.error_message = "Maximum BDF2 step count reached before t_end";
		}
		result.elapsed_time_ms = timer.elapsed_ms();

		return result;
	}

} // namespace MML

#endif // MML_DAE_BDF2_H

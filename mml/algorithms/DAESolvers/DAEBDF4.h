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
#include <vector>

namespace MML {

	/// @brief Fourth-order Backward Differentiation Formula solver for index-1 DAEs
	/// 
	/// Complexity: O(N³) per step where N = diffDim + algDim.
	///            Newton iteration with augmented Jacobian + LU solve each step.
	///            First 3 steps use BDF2 sub-stepping for accurate history.
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
	/// Uses sub-stepping with BDF2 for first 3 steps to build accurate history,
	/// then switches to BDF4 for remaining steps.
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
		if (t_end <= t0)
			throw ArgumentError("SolveDAEBDF4: t_end must be greater than t0 (reverse-time integration is not supported)");
		if (config.step_size <= 0)
			throw ArgumentError("SolveDAEBDF4: config.step_size must be positive");

		AlgorithmTimer timer;

		int diffDim = system.getDiffDim();
		int algDim = system.getAlgDim();
		int num_steps = static_cast<int>((t_end - t0) / config.step_size) + 1;

		DAESolverResult result(t0, t_end, diffDim, algDim, num_steps);
		result.algorithm_name = "DAEBDF4";

		// Working vectors
		Vector<Real> dxdt(diffDim);
		Vector<Real> g(algDim);
		Matrix<Real> df_dx(diffDim, diffDim);
		Matrix<Real> df_dy(diffDim, algDim);
		Matrix<Real> dg_dx(algDim, diffDim);
		Matrix<Real> dg_dy(algDim, algDim);

		// BDF coefficients
		// BDF2: x_new = (4/3)*x_n - (1/3)*x_{n-1} + (2/3)*h*f
		constexpr Real bdf2_c0 = 4.0 / 3.0;
		constexpr Real bdf2_c1 = -1.0 / 3.0;
		constexpr Real bdf2_h_coeff = 2.0 / 3.0;

		// BDF4: x_new = (48/25)*x_n - (36/25)*x_{n-1} + (16/25)*x_{n-2} - (3/25)*x_{n-3} + (12/25)*h*f
		constexpr Real bdf4_c0 = 48.0 / 25.0;
		constexpr Real bdf4_c1 = -36.0 / 25.0;
		constexpr Real bdf4_c2 = 16.0 / 25.0;
		constexpr Real bdf4_c3 = -3.0 / 25.0;
		constexpr Real bdf4_h_coeff = 12.0 / 25.0;

		auto newtonSolve = [&](Real t_next, Real h, Vector<Real>& x_new, Vector<Real>& y_new,
		                       int bdf_order,
		                       const Vector<Real>& xn, const Vector<Real>& xn1,
		                       const Vector<Real>& xn2, const Vector<Real>& xn3)
		{
			Real derivativeCoefficient = h;
			if (bdf_order == 4) derivativeCoefficient = bdf4_h_coeff * h;
			else if (bdf_order == 2) derivativeCoefficient = bdf2_h_coeff * h;

			auto evaluateResidual = [&](const Vector<Real>& trialX, const Vector<Real>& trialY, Vector<Real>& residual) {
				system.diffEqs(t_next, trialX, trialY, dxdt);
				system.algConstraints(t_next, trialX, trialY, g);
				if (bdf_order == 4)
				{
					for (int i = 0; i < diffDim; ++i)
						residual[i] = trialX[i] - bdf4_c0 * xn[i] - bdf4_c1 * xn1[i]
							- bdf4_c2 * xn2[i] - bdf4_c3 * xn3[i] - derivativeCoefficient * dxdt[i];
				}
				else if (bdf_order == 2)
				{
					for (int i = 0; i < diffDim; ++i)
						residual[i] = trialX[i] - bdf2_c0 * xn[i] - bdf2_c1 * xn1[i]
							- derivativeCoefficient * dxdt[i];
				}
				else
				{
					for (int i = 0; i < diffDim; ++i)
						residual[i] = trialX[i] - xn[i] - h * dxdt[i];
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

		// ========================================================================
		// BOOTSTRAP PHASE: Build accurate history using sub-stepped BDF2
		// Use 4 substeps per main step to achieve O(h^4) local error during bootstrap
		// ========================================================================
		constexpr int bootstrap_substeps = 4;
		Real h_main = config.step_size;
		Real h_sub = h_main / bootstrap_substeps;

		// We need values at t0, t0+h, t0+2h, t0+3h for BDF4
		// Store these in history array: hist[0]=newest, hist[3]=oldest
		std::vector<Vector<Real>> x_history(4, x0);
		std::vector<Vector<Real>> y_history(4, y0);
		std::vector<Real> t_history(4, t0);

		Real t_current = t0;
		Vector<Real> x_current = x0;
		Vector<Real> y_current = y0;
		Vector<Real> x_sub_prev = x0;  // Previous substep for BDF2
		result.solution.fillValues(0, t0, x0, y0);
		if (!ValidateDAEInitialState(system, t0, x0, y0, config, result)) {
			result.solution.setFinalSize(0);
			result.elapsed_time_ms = timer.elapsed_ms();
			return result;
		}

		// Take 3 main steps using sub-stepped BDF2
		for (int main_step = 1; main_step <= 3 && t_current < t_end; ++main_step)
		{
			Real h_this = std::min(h_main, t_end - t_current);
			Real h_substep = h_this / bootstrap_substeps;

			for (int sub = 0; sub < bootstrap_substeps; ++sub)
			{
				Real t_next = t_current + h_substep;
				
				system.diffEqs(t_current, x_current, y_current, dxdt);
				Vector<Real> x_pred = x_current + dxdt * h_substep;
				Vector<Real> y_pred = y_current;

				int order = (main_step == 1 && sub == 0) ? 1 : 2;  // First substep uses BE
				DAENewtonResult newton = newtonSolve(t_next, h_substep, x_pred, y_pred, order,
					x_current, x_sub_prev, x_current, x_current);
				AccumulateDAENewtonDiagnostics(result, newton);
				if (!newton.converged)
				{
					SetDAENewtonFailure(result, newton, "BDF4 bootstrap Newton failure at t=" + std::to_string(t_next));
					result.elapsed_time_ms = timer.elapsed_ms();
					return result;
				}

				x_sub_prev = x_current;
				x_current = x_pred;
				y_current = y_pred;
				t_current = t_next;
				system.algConstraints(t_current, x_current, y_current, g);
				if (!ValidateAcceptedDAEState(x_current, y_current, g, config, result, "BDF4 bootstrap step")) {
					result.elapsed_time_ms = timer.elapsed_ms();
					return result;
				}
			}

			// Shift history and store main step result
			for (int i = 3; i > 0; --i)
			{
				x_history[i] = x_history[i-1];
				y_history[i] = y_history[i-1];
				t_history[i] = t_history[i-1];
			}
			x_history[0] = x_current;
			y_history[0] = y_current;
			t_history[0] = t_current;
		}

		// Save initial condition and bootstrap points to solution
		for (int i = 1; i <= 3; ++i)
		{
			result.solution.fillValues(i, t_history[3-i], x_history[3-i], y_history[3-i]);
			result.solution.incrementSuccessfulSteps();
			result.accepted_steps++;
		}
		int step = 4;

		// ========================================================================
		// MAIN PHASE: BDF4 with accurate history
		// ========================================================================
		while (t_current < t_end && step < config.max_steps)
		{
			Real h = std::min(h_main, t_end - t_current);
			
			// Skip degenerate final steps
			if (h < h_main * 1e-10)
				break;
			
			Real t_next = t_current + h;

			// For very small final steps, fall back to lower order
			int order = 4;
			if (h < h_main * 0.1)
				order = 1;  // Backward Euler for tiny final step

			// Prediction
			system.diffEqs(t_current, x_history[0], y_history[0], dxdt);
			Vector<Real> x_new = x_history[0] + dxdt * h;
			Vector<Real> y_new = y_history[0];

			DAENewtonResult newton = newtonSolve(t_next, h, x_new, y_new, order,
				x_history[0], x_history[1], x_history[2], x_history[3]);
			AccumulateDAENewtonDiagnostics(result, newton);
			if (!newton.converged)
			{
				SetDAENewtonFailure(result, newton, "BDF4 Newton failure at t=" + std::to_string(t_next));
				break;
			}

			// Shift history
			for (int i = 3; i > 0; --i)
			{
				x_history[i] = x_history[i-1];
				y_history[i] = y_history[i-1];
				t_history[i] = t_history[i-1];
			}
			x_history[0] = x_new;
			y_history[0] = y_new;
			t_history[0] = t_next;
			t_current = t_next;

			// Track constraint violation
			system.algConstraints(t_current, x_history[0], y_history[0], g);
			if (!ValidateAcceptedDAEState(x_history[0], y_history[0], g, config, result, "BDF4 step"))
				break;

			result.solution.fillValues(step, t_current, x_history[0], y_history[0]);
			result.solution.incrementSuccessfulSteps();
			result.accepted_steps++;
			++step;
		}

		// Finalize
		result.solution.setFinalSize(step - 1);
		result.total_steps = step - 1;
		result.final_constraint_norm = g.NormL2();
		if (t_end - t_current > h_main * REAL(1e-10) && result.status == AlgorithmStatus::Success) {
			result.status = AlgorithmStatus::MaxIterationsExceeded;
			result.failure_reason = DAEFailureReason::MaxStepsExceeded;
			result.error_message = "Maximum BDF4 step count reached before t_end";
		}
		result.elapsed_time_ms = timer.elapsed_ms();

		return result;
	}

} // namespace MML

#endif // MML_DAE_BDF4_H

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
		AlgorithmTimer timer;

		int diffDim = system.getDiffDim();
		int algDim = system.getAlgDim();
		int totalDim = diffDim + algDim;
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

		// Augmented system Jacobian and residual
		Matrix<Real> J_aug(totalDim, totalDim);
		Vector<Real> residual(totalDim);
		Vector<Real> delta(totalDim);

		// Save initial condition
		result.solution.fillValues(0, t, x, y);
		int step = 1;

		// Helper function for Newton solve (returns true if converged)
		auto newtonSolve = [&](Real t_next, Real h, const Vector<Real>& x_pred,
		                       Vector<Real>& x_new, Vector<Real>& y_new,
		                       bool useBDF2, const Vector<Real>& x_prev_bdf) -> bool
		{
			x_new = x_pred;

			for (int iter = 0; iter < config.max_newton_iter; ++iter)
			{
				result.newton_iterations++;

				// Evaluate differential equations and constraints
				system.diffEqs(t_next, x_new, y_new, dxdt);
				system.algConstraints(t_next, x_new, y_new, g);

				// Build residual based on method
				if (useBDF2)
				{
					// BDF2: x_new = (4/3)*x - (1/3)*x_prev + (2/3)*h*f
					Real c0 = 4.0 / 3.0;
					Real c1 = -1.0 / 3.0;
					Real c2 = 2.0 / 3.0;
					for (int i = 0; i < diffDim; ++i)
						residual[i] = x_new[i] - c0 * x[i] - c1 * x_prev_bdf[i] - c2 * h * dxdt[i];
				}
				else
				{
					// Backward Euler: x_new = x + h*f
					for (int i = 0; i < diffDim; ++i)
						residual[i] = x_new[i] - x[i] - h * dxdt[i];
				}
				for (int i = 0; i < algDim; ++i)
					residual[diffDim + i] = g[i];

				// Check convergence
				Real res_norm = residual.NormL2();
				if (res_norm < config.newton_tol)
					return true;

				// Evaluate Jacobians
				system.allJacobians(t_next, x_new, y_new, df_dx, df_dy, dg_dx, dg_dy);
				result.jacobian_evaluations++;

				// Build augmented Jacobian
				Real h_eff = useBDF2 ? (2.0 / 3.0) * h : h;
				for (int i = 0; i < diffDim; ++i)
				{
					for (int j = 0; j < diffDim; ++j)
						J_aug(i, j) = (i == j ? 1.0 : 0.0) - h_eff * df_dx(i, j);
					for (int j = 0; j < algDim; ++j)
						J_aug(i, diffDim + j) = -h_eff * df_dy(i, j);
				}
				for (int i = 0; i < algDim; ++i)
				{
					for (int j = 0; j < diffDim; ++j)
						J_aug(diffDim + i, j) = dg_dx(i, j);
					for (int j = 0; j < algDim; ++j)
						J_aug(diffDim + i, diffDim + j) = dg_dy(i, j);
				}

				// Solve J_aug * delta = -residual
				delta = residual * (-1.0);
				Matrix<Real> J_copy = J_aug;
				GaussJordanSolver<Real>::SolveInPlace(J_copy, delta);

				// Update solution
				for (int i = 0; i < diffDim; ++i)
					x_new[i] += delta[i];
				for (int i = 0; i < algDim; ++i)
					y_new[i] += delta[diffDim + i];
			}
			return false;
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

			if (!newtonSolve(t_next, h, x_new, x_new, y_new, false, x_prev))
			{
				result.status = AlgorithmStatus::NumericalInstability;
				result.error_message = "Newton iteration failed at first step (Backward Euler)";
				result.elapsed_time_ms = timer.elapsed_ms();
				return result;
			}

			// Accept step
			x_prev = x;
			x = x_new;
			y = y_new;
			t = t_next;

			system.algConstraints(t, x, y, g);
			Real g_norm = g.NormL2();
			if (g_norm > result.max_constraint_violation)
				result.max_constraint_violation = g_norm;

			result.solution.fillValues(step, t, x, y);
			result.solution.incrementSuccessfulSteps();
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

			if (!newtonSolve(t_next, h, x_new, x_new, y_new, !use_be_for_final, x_prev))
			{
				result.status = AlgorithmStatus::NumericalInstability;
				result.error_message = "Newton iteration failed at t=" + std::to_string(t_next);
				break;
			}

			// Accept step
			x_prev = x;
			x = x_new;
			y = y_new;
			t = t_next;

			system.algConstraints(t, x, y, g);
			Real g_norm = g.NormL2();
			if (g_norm > result.max_constraint_violation)
				result.max_constraint_violation = g_norm;

			result.solution.fillValues(step, t, x, y);
			result.solution.incrementSuccessfulSteps();
			++step;
		}

		// Finalize
		result.solution.setFinalSize(step - 1);
		result.total_steps = step - 1;
		result.elapsed_time_ms = timer.elapsed_ms();

		return result;
	}

} // namespace MML

#endif // MML_DAE_BDF2_H

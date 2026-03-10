///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DAEBackwardEuler.h                                                  ///
///  Description: Backward Euler method for semi-explicit index-1 DAE systems         ///
///                                                                                   ///
///  The method solves the coupled system at each step using Newton iteration:        ///
///    x_{n+1} = x_n + h·f(t_{n+1}, x_{n+1}, y_{n+1})                                 ///
///    0 = g(t_{n+1}, x_{n+1}, y_{n+1})                                               ///
///                                                                                   ///
///  Properties:                                                                      ///
///  - 1st order accurate                                                             ///
///  - A-stable (handles stiff problems)                                              ///
///  - Simple baseline for comparison                                                 ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DAE_BACKWARD_EULER_H
#define MML_DAE_BACKWARD_EULER_H

#include "DAESolverBase.h"

namespace MML {

	/// @brief Backward Euler implicit solver for index-1 DAEs
	/// 
	/// Complexity: O(N³) per step where N = diffDim + algDim.
	///            Each step: Jacobian eval O(N²) + LU solve O(N³) × Newton iterations.
	///
	/// Solves the coupled system at each step:
	/// @f[
	///   x_{n+1} = x_n + h \cdot f(t_{n+1}, x_{n+1}, y_{n+1})
	/// @f]
	/// @f[
	///   0 = g(t_{n+1}, x_{n+1}, y_{n+1})
	/// @f]
	///
	/// Using Newton iteration on the combined residual with the augmented Jacobian:
	/// @f[
	///   \begin{bmatrix}
	///     I - h\frac{\partial f}{\partial x} & -h\frac{\partial f}{\partial y} \\
	///     \frac{\partial g}{\partial x} & \frac{\partial g}{\partial y}
	///   \end{bmatrix}
	///   \begin{bmatrix} \Delta x \\ \Delta y \end{bmatrix}
	///   =
	///   \begin{bmatrix} -F \\ -G \end{bmatrix}
	/// @f]
	///
	/// where F = x_new - x_old - h*f, G = g.
	///
	/// @param system DAE system with Jacobian
	/// @param t0 Initial time
	/// @param x0 Initial differential state
	/// @param y0 Initial algebraic state (must satisfy constraints)
	/// @param t_end Final time
	/// @param config Solver configuration
	/// @return DAESolverResult with solution and diagnostics
	inline DAESolverResult SolveDAEBackwardEuler(IODESystemDAEWithJacobian& system,
	                                              Real t0, const Vector<Real>& x0, const Vector<Real>& y0,
	                                              Real t_end, const DAESolverConfig& config = DAESolverConfig())
	{
		AlgorithmTimer timer;

		int diffDim = system.getDiffDim();
		int algDim = system.getAlgDim();
		int totalDim = diffDim + algDim;
		int num_steps = static_cast<int>((t_end - t0) / config.step_size) + 1;

		DAESolverResult result(t0, t_end, diffDim, algDim, num_steps);
		result.algorithm_name = "DAEBackwardEuler";

		// Current state
		Real t = t0;
		Vector<Real> x = x0;
		Vector<Real> y = y0;

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

		while (t < t_end && step < config.max_steps)
		{
			Real h = std::min(config.step_size, t_end - t);
			Real t_next = t + h;

			// Initial guess for next step (explicit Euler prediction)
			system.diffEqs(t, x, y, dxdt);
			Vector<Real> x_new = x + dxdt * h;
			Vector<Real> y_new = y;

			// Newton iteration
			bool converged = false;
			for (int iter = 0; iter < config.max_newton_iter; ++iter)
			{
				result.newton_iterations++;

				// Evaluate differential equations and constraints at new point
				system.diffEqs(t_next, x_new, y_new, dxdt);
				system.algConstraints(t_next, x_new, y_new, g);

				// Build residual: F = x_new - x - h*f, G = g
				for (int i = 0; i < diffDim; ++i)
					residual[i] = x_new[i] - x[i] - h * dxdt[i];
				for (int i = 0; i < algDim; ++i)
					residual[diffDim + i] = g[i];

				// Check convergence
				Real res_norm = residual.NormL2();
				if (res_norm < config.newton_tol)
				{
					converged = true;
					break;
				}

				// Evaluate Jacobians
				system.allJacobians(t_next, x_new, y_new, df_dx, df_dy, dg_dx, dg_dy);
				result.jacobian_evaluations++;

				// Build augmented Jacobian:
				// [ I - h*df/dx,   -h*df/dy ]
				// [     dg/dx,       dg/dy  ]
				for (int i = 0; i < diffDim; ++i)
				{
					for (int j = 0; j < diffDim; ++j)
						J_aug(i, j) = (i == j ? 1.0 : 0.0) - h * df_dx(i, j);
					for (int j = 0; j < algDim; ++j)
						J_aug(i, diffDim + j) = -h * df_dy(i, j);
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
				Matrix<Real> J_copy = J_aug;  // SolveInPlace modifies matrix
				GaussJordanSolver<Real>::SolveInPlace(J_copy, delta);

				// Update solution
				for (int i = 0; i < diffDim; ++i)
					x_new[i] += delta[i];
				for (int i = 0; i < algDim; ++i)
					y_new[i] += delta[diffDim + i];
			}

			if (!converged)
			{
				result.status = AlgorithmStatus::NumericalInstability;
				result.error_message = "Newton iteration failed to converge at t=" + std::to_string(t_next);
				break;
			}

			// Accept step
			x = x_new;
			y = y_new;
			t = t_next;

			// Track constraint violation
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

#endif // MML_DAE_BACKWARD_EULER_H

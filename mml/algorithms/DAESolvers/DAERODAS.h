///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DAERODAS.h                                                          ///
///  Description: RODAS (Rosenbrock) method for semi-explicit index-1 DAE systems     ///
///                                                                                   ///
///  Uses linearly implicit approach - only linear system solves, no Newton iteration ///
///  This makes it efficient for problems where Jacobian computation is expensive.    ///
///                                                                                   ///
///  Properties:                                                                      ///
///  - L-stable (strong damping of stiff components)                                  ///
///  - No Newton iteration required                                                   ///
///  - ~1st-2nd order accuracy (linearly-implicit predictor-corrector)                ///
///                                                                                   ///
///  Reference: Hairer & Wanner, "Solving ODEs II"                                    ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DAE_RODAS_H
#define MML_DAE_RODAS_H

#include "DAESolverBase.h"

namespace MML {

	/// @brief RODAS - 2-stage Rosenbrock method for index-1 DAEs
	/// 
	/// Complexity: O(N³) per step where N = diffDim + algDim.
	///            Only 1-2 LU factorizations per step (no Newton iteration).
	///            Most efficient implicit method per step for small-to-medium N.
	///
	/// Uses the linearly implicit approach that requires only linear system solves,
	/// not iterative Newton methods. Each step solves a linear system once (or twice
	/// for 2-stage methods) using a frozen Jacobian.
	///
	/// This implementation uses a 2-stage method based on the trapezoidal rule:
	/// 1. Solve (I - h*J/2) * Δ = h*F(t,x,y) for increment Δ
	/// 2. Update x = x + Δ_x, y = y + Δ_y
	///
	/// Reference: Hairer & Wanner, "Solving ODEs II"
	///
	/// @param system DAE system with Jacobian
	/// @param t0 Initial time
	/// @param x0 Initial differential state
	/// @param y0 Initial algebraic state (must satisfy constraints)
	/// @param t_end Final time
	/// @param config Solver configuration
	/// @return DAESolverResult with solution and diagnostics
	inline DAESolverResult SolveDAERODAS(IODESystemDAEWithJacobian& system,
	                                      Real t0, const Vector<Real>& x0, const Vector<Real>& y0,
	                                      Real t_end, const DAESolverConfig& config = DAESolverConfig())
	{
		AlgorithmTimer timer;

		int diffDim = system.getDiffDim();
		int algDim = system.getAlgDim();
		int totalDim = diffDim + algDim;
		int num_steps = static_cast<int>((t_end - t0) / config.step_size) + 1;

		DAESolverResult result(t0, t_end, diffDim, algDim, num_steps);
		result.algorithm_name = "DAERODAS";

		// Working vectors
		Matrix<Real> df_dx(diffDim, diffDim);
		Matrix<Real> df_dy(diffDim, algDim);
		Matrix<Real> dg_dx(algDim, diffDim);
		Matrix<Real> dg_dy(algDim, algDim);

		// Augmented Jacobian matrix
		Matrix<Real> J_aug(totalDim, totalDim);

		// Vectors for predictor-corrector
		Vector<Real> delta(totalDim);
		Vector<Real> rhs(totalDim);

		// Current state
		Real t = t0;
		Vector<Real> x = x0;
		Vector<Real> y = y0;

		// Temporary states
		Vector<Real> x_pred(diffDim);
		Vector<Real> y_pred(algDim);
		Vector<Real> f_curr(diffDim);
		Vector<Real> f_pred(diffDim);
		Vector<Real> g_curr(algDim);
		Vector<Real> g_pred(algDim);

		// Save initial condition
		result.solution.fillValues(0, t, x, y);
		int step = 1;

		// Main time-stepping loop
		while (t < t_end && step < config.max_steps)
		{
			Real h = std::min(config.step_size, t_end - t);
			
			if (h < config.step_size * 1e-10)
				break;

			// Evaluate f, g and Jacobians at current point
			system.diffEqs(t, x, y, f_curr);
			system.algConstraints(t, x, y, g_curr);
			system.allJacobians(t, x, y, df_dx, df_dy, dg_dx, dg_dy);
			result.jacobian_evaluations++;

			// ================================================================
			// Rosenbrock step using linearly-implicit trapezoidal method
			// 
			// For x' = f(x,y), 0 = g(x,y), the trapezoidal rule gives:
			// x_{n+1} = x_n + (h/2)*(f_n + f_{n+1})
			// 0 = g_{n+1}
			//
			// Linearizing f_{n+1} ≈ f_n + J*(Δx, Δy):
			// x_{n+1} - x_n = (h/2)*(f_n + f_n + df/dx*(x_{n+1}-x_n) + df/dy*(y_{n+1}-y_n))
			// 
			// Let Δx = x_{n+1}-x_n, Δy = y_{n+1}-y_n:
			// (I - h/2*df/dx)*Δx - h/2*df/dy*Δy = h*f_n
			// dg/dx*Δx + dg/dy*Δy = -g_n  (project to constraint)
			// ================================================================
			
			Real hHalf = h * 0.5;
			
			// Build augmented system:
			// [I - h/2*df/dx,    -h/2*df/dy] [Δx]   [h*f_n   ]
			// [dg/dx,            dg/dy     ] [Δy] = [-g_n    ]
			for (int i = 0; i < diffDim; ++i)
			{
				for (int j = 0; j < diffDim; ++j)
					J_aug(i, j) = (i == j ? 1.0 : 0.0) - hHalf * df_dx(i, j);
				for (int j = 0; j < algDim; ++j)
					J_aug(i, diffDim + j) = -hHalf * df_dy(i, j);
			}
			for (int i = 0; i < algDim; ++i)
			{
				for (int j = 0; j < diffDim; ++j)
					J_aug(diffDim + i, j) = dg_dx(i, j);
				for (int j = 0; j < algDim; ++j)
					J_aug(diffDim + i, diffDim + j) = dg_dy(i, j);
			}

			// Build RHS
			for (int i = 0; i < diffDim; ++i)
				rhs[i] = h * f_curr[i];
			for (int i = 0; i < algDim; ++i)
				rhs[diffDim + i] = -g_curr[i];

			// Solve J_aug * delta = rhs using Gauss-Jordan
			delta = rhs;
			Matrix<Real> J_copy = J_aug;
			GaussJordanSolver<Real>::SolveInPlace(J_copy, delta);

			// First approximation (predictor)
			for (int i = 0; i < diffDim; ++i)
				x_pred[i] = x[i] + delta[i];
			for (int i = 0; i < algDim; ++i)
				y_pred[i] = y[i] + delta[diffDim + i];

			// ================================================================
			// Second stage: corrector using updated function values
			// Re-evaluate at the predictor point and do one more linear solve
			// ================================================================
			system.diffEqs(t + h, x_pred, y_pred, f_pred);
			system.algConstraints(t + h, x_pred, y_pred, g_pred);

			// Corrector: use average of f values
			// Build RHS: h/2*(f_n + f_{n+1}) for diff eqs
			for (int i = 0; i < diffDim; ++i)
				rhs[i] = hHalf * (f_curr[i] + f_pred[i]);
			for (int i = 0; i < algDim; ++i)
				rhs[diffDim + i] = -g_pred[i];

			// Solve again with same Jacobian (Rosenbrock approach)
			delta = rhs;
			J_copy = J_aug;
			GaussJordanSolver<Real>::SolveInPlace(J_copy, delta);

			// Update solution
			for (int i = 0; i < diffDim; ++i)
				x[i] += delta[i];
			for (int i = 0; i < algDim; ++i)
				y[i] += delta[diffDim + i];

			t += h;

			// Track constraint violation
			system.algConstraints(t, x, y, g_pred);
			Real g_norm = g_pred.NormL2();
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

#endif // MML_DAE_RODAS_H

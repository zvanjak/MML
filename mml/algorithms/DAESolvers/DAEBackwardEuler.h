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
		if (t_end <= t0)
			throw ArgumentError("SolveDAEBackwardEuler: t_end must be greater than t0 (reverse-time integration is not supported)");
		if (config.step_size <= 0)
			throw ArgumentError("SolveDAEBackwardEuler: config.step_size must be positive");

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

		// Augmented system Jacobian
		Matrix<Real> J_aug(totalDim, totalDim);

		// Save initial condition
		result.solution.fillValues(0, t, x, y);
		system.algConstraints(t, x, y, g);
		result.final_constraint_norm = g.NormL2();
		result.max_constraint_violation = result.final_constraint_norm;
		if (!IsFiniteDAEVector(x) || !IsFiniteDAEVector(y) || !IsFiniteDAEVector(g)) {
			result.status = AlgorithmStatus::NumericalInstability;
			result.failure_reason = DAEFailureReason::NonFiniteState;
			result.error_message = "Non-finite initial DAE state or constraint residual";
			result.solution.setFinalSize(0);
			result.elapsed_time_ms = timer.elapsed_ms();
			return result;
		}
		if (result.final_constraint_norm > config.constraint_tol) {
			result.status = AlgorithmStatus::InvalidInput;
			result.failure_reason = DAEFailureReason::InconsistentInitialConditions;
			result.error_message = "Initial algebraic constraints are inconsistent";
			result.solution.setFinalSize(0);
			result.elapsed_time_ms = timer.elapsed_ms();
			return result;
		}
		int step = 1;

		while (t < t_end && step < config.max_steps)
		{
			Real h = std::min(config.step_size, t_end - t);
			Real t_next = t + h;

			// Initial guess for next step (explicit Euler prediction)
			system.diffEqs(t, x, y, dxdt);
			Vector<Real> x_new = x + dxdt * h;
			Vector<Real> y_new = y;

			auto evaluateResidual = [&](const Vector<Real>& trialX, const Vector<Real>& trialY, Vector<Real>& residual) {
				system.diffEqs(t_next, trialX, trialY, dxdt);
				system.algConstraints(t_next, trialX, trialY, g);

				for (int i = 0; i < diffDim; ++i)
					residual[i] = trialX[i] - x[i] - h * dxdt[i];
				for (int i = 0; i < algDim; ++i)
					residual[diffDim + i] = g[i];
			};
			auto evaluateJacobian = [&](const Vector<Real>& trialX, const Vector<Real>& trialY, Matrix<Real>& jacobian) {
				system.allJacobians(t_next, trialX, trialY, df_dx, df_dy, dg_dx, dg_dy);
				AssembleDAEAugmentedJacobian(df_dx, df_dy, dg_dx, dg_dy, h, jacobian);
			};

			DAENewtonResult newton = SolveDAENewton(x_new, y_new, diffDim, algDim, config, evaluateResidual, evaluateJacobian);
			AccumulateDAENewtonDiagnostics(result, newton);

			if (!newton.converged)
			{
				result.status = newton.failure_reason == DAEFailureReason::SingularJacobian
					? AlgorithmStatus::SingularMatrix : AlgorithmStatus::NumericalInstability;
				result.failure_reason = newton.failure_reason;
				result.error_message = "DAE Newton failure at t=" + std::to_string(t_next) + ": " + newton.message;
				break;
			}

			// Accept step
			x = x_new;
			y = y_new;
			t = t_next;

			// Track constraint violation
			system.algConstraints(t, x, y, g);
			if (!ValidateAcceptedDAEState(x, y, g, config, result, "Backward Euler step"))
				break;

			result.solution.fillValues(step, t, x, y);
			result.solution.incrementSuccessfulSteps();
			result.accepted_steps++;
			++step;
		}

		// Finalize
		result.solution.setFinalSize(step - 1);
		result.total_steps = step - 1;
		if (t < t_end && result.status == AlgorithmStatus::Success) {
			result.status = AlgorithmStatus::MaxIterationsExceeded;
			result.failure_reason = DAEFailureReason::MaxStepsExceeded;
			result.error_message = "Maximum DAE step count reached before t_end";
		}
		result.elapsed_time_ms = timer.elapsed_ms();

		return result;
	}

} // namespace MML

#endif // MML_DAE_BACKWARD_EULER_H

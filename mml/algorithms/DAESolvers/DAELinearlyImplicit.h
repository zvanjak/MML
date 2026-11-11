///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DAELinearlyImplicit.h                                               ///
///  Description: Linearly-implicit trapezoid (Rosenbrock-type) for index-1 DAEs      ///
///                                                                                   ///
///  Uses linearly implicit approach - only linear system solves, no Newton iteration ///
///  This makes it efficient for problems where Jacobian computation is expensive.    ///
///                                                                                   ///
///  Properties:                                                                      ///
///  - Second order for smooth index-1 DAEs                                           ///
///  - A-stable on linear problems, NOT L-stable                                     ///
///  - No Newton iteration required                                                   ///
///  - One frozen-matrix predictor and residual correction                           ///
///                                                                                   ///
///  Reference: Hairer & Wanner, "Solving ODEs II"                                    ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DAE_LINEARLY_IMPLICIT_H
#define MML_DAE_LINEARLY_IMPLICIT_H

#include "DAESolverBase.h"
#include <algorithm>

namespace MML {

	/// @brief Linearly-implicit trapezoid (predictor-corrector, Rosenbrock-type) for index-1 DAEs
	///
	/// @note Until analysis v4 this function was called SolveDAERODAS. It is NOT the 4th-order
	///       6-stage RODAS of Hairer & Wanner. It is a second-order linearly-implicit
	///       trapezoidal scheme with a frozen Jacobian. The old name remains a deprecated
	///       alias. Its linear stability function is (1+z/2)/(1-z/2), tending to -1,
	///       not zero, for stiff decay. Use SolveDAERODAS4 for L-stable damping.
	/// 
	/// Complexity: O(N³) per step where N = diffDim + algDim.
	///            One LU factorization and two back-substitutions per step (no Newton iteration).
	///            Most efficient implicit method per step for small-to-medium N.
	///
	/// The predictor linearizes the trapezoidal equations at the current state. A
	/// second solve corrects their residual at the predicted endpoint, using the
	/// same matrix. Nonlinear constraints are checked, not iterated to convergence.
	///
	/// Reference: Hairer & Wanner, "Solving ODEs II", Section IV.7 (linearly implicit methods)
	///
	/// @param system DAE system with Jacobian
	/// @param t0 Initial time
	/// @param x0 Initial differential state
	/// @param y0 Initial algebraic state (must satisfy constraints)
	/// @param t_end Final time
	/// @param config Solver configuration
	/// @return DAESolverResult with solution and diagnostics
	inline DAESolverResult SolveDAELinearlyImplicit(IODESystemDAEWithJacobian& system,
	                                      Real t0, const Vector<Real>& x0, const Vector<Real>& y0,
	                                      Real t_end, const DAESolverConfig& config = DAESolverConfig())
	{
		const int diffDim = system.getDiffDim(), algDim = system.getAlgDim();
		if (!std::isfinite(t0) || !std::isfinite(t_end) || !std::isfinite(t_end - t0) || t_end == t0 ||
			!std::isfinite(config.step_size) || config.step_size <= 0 ||
			!std::isfinite(config.constraint_tol) || config.constraint_tol <= 0 || config.max_steps <= 0 ||
			diffDim <= 0 || algDim < 0 || x0.size() != diffDim || y0.size() != algDim)
			throw ArgumentError("SolveDAELinearlyImplicit: requires distinct finite endpoints, positive finite "
				"step and constraint tolerance, a positive step budget and matching state dimensions");

		AlgorithmTimer timer;

		int totalDim = diffDim + algDim;
		DAESolverResult result(t0, t_end, diffDim, algDim, 0);
		result.algorithm_name = "DAELinearlyImplicit";

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
		const Real direction = t_end > t0 ? Real(1) : Real(-1);
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
		if (!ValidateDAEInitialState(system, t, x, y, config, result)) {
			result.solution.setFinalSize(0);
			result.elapsed_time_ms = timer.elapsed_ms();
			return result;
		}
		int step = 1;
		const auto rejectAttempt = [&]() {
			++result.rejected_steps;
			result.solution.incrementRejectedSteps();
		};

		// Main time-stepping loop
		while (direction * (t_end - t) > 0 && step <= config.max_steps)
		{
			if (!PrepareDAETrajectoryPoint(result, step)) break;
			Real h = direction * std::min(config.step_size, std::abs(t_end - t));
			if (!ValidateDAETimeStep(t, t + h, result, "Linearly-implicit DAE")) break;
			
			if (std::abs(h) < config.step_size * Precision::MinimumStepRatio)
				break;

			// Evaluate f, g and Jacobians at current point
			system.diffEqs(t, x, y, f_curr);
			system.algConstraints(t, x, y, g_curr);
			system.allJacobians(t, x, y, df_dx, df_dy, dg_dx, dg_dy);
			result.jacobian_evaluations++;
			if (!IsFiniteDAEVector(f_curr) || !IsFiniteDAEVector(g_curr)) {
				result.status = AlgorithmStatus::NumericalInstability;
				result.failure_reason = DAEFailureReason::NonFiniteState;
				result.error_message = "Linearly-implicit step has a non-finite derivative or constraint";
				rejectAttempt();
				break;
			}

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
			
			Real hHalf = h * MML::Real(0.5);
			
			AssembleDAEAugmentedJacobian(df_dx, df_dy, dg_dx, dg_dy, hHalf, J_aug);

			// Build RHS
			for (int i = 0; i < diffDim; ++i)
				rhs[i] = h * f_curr[i];
			for (int i = 0; i < algDim; ++i)
				rhs[diffDim + i] = -g_curr[i];

			Detail::ImplicitNewtonMatrix iterationMatrix;
			try {
				iterationMatrix.Factor(J_aug);
				++result.lu_factorizations;
				delta = iterationMatrix.Solve(rhs);
				++result.linear_solves;
			}
			catch (const SingularMatrixError& error) {
				result.status = AlgorithmStatus::SingularMatrix;
				result.failure_reason = DAEFailureReason::SingularJacobian;
				result.error_message = std::string("Linearly-implicit predictor solve failed: ") + error.what();
				rejectAttempt();
				break;
			}
			catch (const MatrixNumericalError& error) {
				result.status = AlgorithmStatus::NumericalInstability;
				result.failure_reason = DAEFailureReason::NonFiniteState;
				result.error_message = error.what();
				rejectAttempt();
				break;
			}

			// First approximation (predictor)
			for (int i = 0; i < diffDim; ++i)
				x_pred[i] = x[i] + delta[i];
			for (int i = 0; i < algDim; ++i)
				y_pred[i] = y[i] + delta[diffDim + i];
			if (!IsFiniteDAEVector(x_pred) || !IsFiniteDAEVector(y_pred)) {
				result.status = AlgorithmStatus::NumericalInstability;
				result.failure_reason = DAEFailureReason::NonFiniteState;
				result.error_message = "Linearly-implicit predictor produced a non-finite state";
				rejectAttempt();
				break;
			}

			// ================================================================
			// Second stage: corrector using updated function values
			// Re-evaluate at the predictor point and do one more linear solve
			// ================================================================
			system.diffEqs(t + h, x_pred, y_pred, f_pred);
			system.algConstraints(t + h, x_pred, y_pred, g_pred);
			if (!IsFiniteDAEVector(f_pred) || !IsFiniteDAEVector(g_pred)) {
				result.status = AlgorithmStatus::NumericalInstability;
				result.failure_reason = DAEFailureReason::NonFiniteState;
				result.error_message = "Linearly-implicit corrector has a non-finite derivative or constraint";
				rejectAttempt();
				break;
			}

			// One frozen-Jacobian correction of the trapezoidal residual (Hairer-Wanner II, IV.7).
			// Subtract the predictor increment, then apply the correction to the predictor.
			for (int i = 0; i < diffDim; ++i)
				rhs[i] = hHalf * f_curr[i] + hHalf * f_pred[i] - (x_pred[i] - x[i]);
			for (int i = 0; i < algDim; ++i)
				rhs[diffDim + i] = -g_pred[i];

			// Solve again with same Jacobian (Rosenbrock approach)
			try {
				delta = iterationMatrix.Solve(rhs);
				++result.linear_solves;
			}
			catch (const SingularMatrixError& error) {
				result.status = AlgorithmStatus::SingularMatrix;
				result.failure_reason = DAEFailureReason::SingularJacobian;
				result.error_message = std::string("Linearly-implicit corrector solve failed: ") + error.what();
				rejectAttempt();
				break;
			}
			catch (const MatrixNumericalError& error) {
				result.status = AlgorithmStatus::NumericalInstability;
				result.failure_reason = DAEFailureReason::NonFiniteState;
				result.error_message = std::string("Linearly-implicit corrector solve failed: ") + error.what();
				rejectAttempt();
				break;
			}

			for (int i = 0; i < diffDim; ++i)
				x_pred[i] += delta[i];
			for (int i = 0; i < algDim; ++i)
				y_pred[i] += delta[diffDim + i];
			if (!IsFiniteDAEVector(x_pred) || !IsFiniteDAEVector(y_pred)) {
				result.status = AlgorithmStatus::NumericalInstability;
				result.failure_reason = DAEFailureReason::NonFiniteState;
				result.error_message = "Linearly-implicit step produced a non-finite state";
				rejectAttempt();
				break;
			}

			system.algConstraints(t + h, x_pred, y_pred, g_pred);
			const Real previousConstraintNorm = result.final_constraint_norm;
			if (!ValidateAcceptedDAEState(x_pred, y_pred, g_pred, config, result, "Linearly-implicit step")) {
				result.final_constraint_norm = previousConstraintNorm;
				rejectAttempt();
				break;
			}
			x = x_pred;
			y = y_pred;
			t += h;

			result.solution.fillValues(step, t, x, y);
			result.solution.incrementSuccessfulSteps();
			result.accepted_steps++;
			++step;
		}

		// Finalize
		if (result.status != AlgorithmStatus::ResourceLimitExceeded)
			result.solution.setFinalSize(step - 1);
		result.total_steps = step - 1;
		if (direction * (t_end - t) > DAETimeTolerance(t_end, config.step_size) && result.status == AlgorithmStatus::Success) {
			result.status = AlgorithmStatus::MaxIterationsExceeded;
			result.failure_reason = DAEFailureReason::MaxStepsExceeded;
			result.error_message = "Maximum step count reached before t_end";
		}
		result.elapsed_time_ms = timer.elapsed_ms();

		return result;
	}

	/// @deprecated Renamed in analysis v4: the method is a linearly-implicit trapezoid, not RODAS.
	[[deprecated("SolveDAERODAS is a linearly-implicit trapezoid, not RODAS; use SolveDAELinearlyImplicit")]]
	inline DAESolverResult SolveDAERODAS(IODESystemDAEWithJacobian& system,
	                                      Real t0, const Vector<Real>& x0, const Vector<Real>& y0,
	                                      Real t_end, const DAESolverConfig& config = DAESolverConfig())
	{
		return SolveDAELinearlyImplicit(system, t0, x0, y0, t_end, config);
	}

} // namespace MML

#endif // MML_DAE_LINEARLY_IMPLICIT_H

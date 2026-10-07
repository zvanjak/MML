///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DAESolverBase.h                                                     ///
///  Description: Base types and utilities for DAE solvers                            ///
///                                                                                   ///
///  Contents:                                                                        ///
///  - DAESolverConfig: Configuration parameters for all DAE solvers                  ///
///  - DAESolverResult: Result structure with solution and diagnostics                ///
///  - ComputeConsistentIC: Compute consistent initial algebraic variables            ///
///  - VerifyConsistentIC: Verify initial conditions satisfy constraints              ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DAE_SOLVER_BASE_H
#define MML_DAE_SOLVER_BASE_H

#include <mml/MMLBase.h>

#include <mml/base/Vector/Vector.h>
#include <mml/base/Matrix/Matrix.h>
#include <mml/base/DAESystem.h>
#include <mml/interfaces/IODESystemDAE.h>
#include <mml/base/AlgorithmTypes.h>
#include <mml/core/LinAlgEqSolvers.h>

#include <algorithm>
#include <cmath>
#include <string>

namespace MML {

	/******************************************************************************/
	/*****            DAE Solver Configuration                                *****/
	/******************************************************************************/

	/// Configuration for DAE solvers.
	///
	/// Provides user control over step size, Newton iteration parameters,
	/// and constraint satisfaction tolerance.
	///
	/// @example
	/// DAESolverConfig config;
	/// config.step_size = 0.001;
	/// config.constraint_tol = 1e-10;
	/// auto sol = SolveDAEBackwardEuler(system, t0, x0, y0, t_end, config);
	struct DAESolverConfig {
		/// Step size for fixed-step methods (default: 0.01)
		Real step_size = 0.01;

		/// Maximum Newton iterations per step (default: 20)
		int max_newton_iter = 20;

		/// Newton convergence tolerance (precision-aware default)
		Real newton_tol = PrecisionValues<Real>::NewtonTolerance;

		/// Tolerance for constraint satisfaction g(t,x,y) ≈ 0 (precision-aware default)
		Real constraint_tol = PrecisionValues<Real>::IntegrationTolerance;

		/// Maximum number of steps (default: 100000)
		int max_steps = 100000;

		/// Absolute and relative tolerances for adaptive methods
		Real abs_tolerance = PrecisionValues<Real>::IntegrationTolerance;
		Real rel_tolerance = PrecisionValues<Real>::IntegrationTolerance * REAL(100.0);

		/// Adaptive step-size bounds and controller safety factor
		Real min_step_size = REAL(1e-8);
		Real max_step_size = REAL(1.0);
		Real step_safety = REAL(0.9);

		/// Factory method for high-precision configuration
		static DAESolverConfig HighPrecision() {
			DAESolverConfig config;
			config.step_size = 0.001;
			config.newton_tol = PrecisionValues<Real>::IntegrationTolerance;
			config.constraint_tol = PrecisionValues<Real>::IntegrationTolerance;
			config.max_newton_iter = 30;
			return config;
		}

		/// Factory method for fast (lower precision) configuration
		static DAESolverConfig Fast() {
			DAESolverConfig config;
			config.step_size = 0.05;
			config.newton_tol = PrecisionValues<Real>::NewtonTolerance * 100;
			config.constraint_tol = PrecisionValues<Real>::IntegrationTolerance * 100;
			config.max_newton_iter = 10;
			return config;
		}
	};

	enum class DAEFailureReason {
		None,
		InconsistentInitialConditions,
		NewtonFailure,
		SingularJacobian,
		HigherIndexSuspected,
		MaxStepsExceeded,
		ConstraintViolation,
		NonFiniteState
	};

	/// @brief Result of DAE integration with diagnostics
	struct DAESolverResult {
		DAESolution solution;			///< Solution trajectory

		// === Diagnostics ===
		std::string algorithm_name;		///< Algorithm used
		AlgorithmStatus status = AlgorithmStatus::Success;
		std::string error_message;
		double elapsed_time_ms = 0.0;

		// === Statistics ===
		int total_steps = 0;
		int newton_iterations = 0;		///< Total Newton iterations
		int jacobian_evaluations = 0;	///< Number of Jacobian evaluations
		int accepted_steps = 0;
		int rejected_steps = 0;
		int last_step_newton_iterations = 0;
		int max_step_newton_iterations = 0;
		Real final_residual_norm = 0.0;
		Real max_residual_norm = 0.0;
		Real final_constraint_norm = 0.0;
		Real max_constraint_violation = 0.0;  ///< Maximum |g(t,x,y)| seen
		DAEFailureReason failure_reason = DAEFailureReason::None;

		DAESolverResult(Real t0, Real t_end, int diffDim, int algDim, int num_steps)
			: solution(t0, t_end, diffDim, algDim, num_steps) {}
	};

	inline bool IsFiniteDAEVector(const Vector<Real>& values) {
		for (int i = 0; i < values.size(); ++i)
			if (!std::isfinite(values[i])) return false;
		return true;
	}

	inline Real DAETimeTolerance(Real tEnd, Real stepSize) {
		return std::max(std::abs(stepSize) * REAL(1e-10),
			REAL(10.0) * std::numeric_limits<Real>::epsilon() * std::max(REAL(1.0), std::abs(tEnd)));
	}

	inline void AssembleDAEAugmentedJacobian(const Matrix<Real>& df_dx, const Matrix<Real>& df_dy,
		const Matrix<Real>& dg_dx, const Matrix<Real>& dg_dy, Real derivativeCoefficient,
		Matrix<Real>& jacobian) {
		const int diffDim = df_dx.rows();
		const int algDim = dg_dy.rows();
		for (int i = 0; i < diffDim; ++i) {
			for (int j = 0; j < diffDim; ++j)
				jacobian(i, j) = (i == j ? REAL(1.0) : REAL(0.0)) - derivativeCoefficient * df_dx(i, j);
			for (int j = 0; j < algDim; ++j)
				jacobian(i, diffDim + j) = -derivativeCoefficient * df_dy(i, j);
		}
		for (int i = 0; i < algDim; ++i) {
			for (int j = 0; j < diffDim; ++j)
				jacobian(diffDim + i, j) = dg_dx(i, j);
			for (int j = 0; j < algDim; ++j)
				jacobian(diffDim + i, diffDim + j) = dg_dy(i, j);
		}
	}

	struct DAENewtonResult {
		bool converged = false;
		int iterations = 0;
		int jacobian_evaluations = 0;
		Real residual_norm = 0.0;
		DAEFailureReason failure_reason = DAEFailureReason::None;
		std::string message;
	};

	template<class ResidualEvaluator, class JacobianEvaluator>
	DAENewtonResult SolveDAENewton(Vector<Real>& x, Vector<Real>& y, int diffDim, int algDim,
		const DAESolverConfig& config, ResidualEvaluator&& evaluateResidual,
		JacobianEvaluator&& evaluateJacobian) {
		DAENewtonResult result;
		Vector<Real> residual(diffDim + algDim);
		Matrix<Real> jacobian(diffDim + algDim, diffDim + algDim);

		for (int iteration = 0; iteration < config.max_newton_iter; ++iteration) {
			result.iterations = iteration + 1;
			evaluateResidual(x, y, residual);
			if (!IsFiniteDAEVector(residual) || !IsFiniteDAEVector(x) || !IsFiniteDAEVector(y)) {
				result.failure_reason = DAEFailureReason::NonFiniteState;
				result.message = "non-finite DAE state or residual";
				return result;
			}

			result.residual_norm = residual.NormL2();
			if (result.residual_norm < config.newton_tol) {
				result.converged = true;
				return result;
			}

			evaluateJacobian(x, y, jacobian);
			result.jacobian_evaluations++;
			Vector<Real> delta = residual * REAL(-1.0);
			try {
				Matrix<Real> jacobianCopy = jacobian;
				GaussJordanSolver<Real>::SolveInPlace(jacobianCopy, delta);
			}
			catch (const std::exception& error) {
				result.failure_reason = DAEFailureReason::SingularJacobian;
				result.message = error.what();
				return result;
			}

			for (int i = 0; i < diffDim; ++i) x[i] += delta[i];
			for (int i = 0; i < algDim; ++i) y[i] += delta[diffDim + i];
		}

		result.failure_reason = DAEFailureReason::NewtonFailure;
		result.message = "Newton iteration limit exceeded";
		return result;
	}

	inline void AccumulateDAENewtonDiagnostics(DAESolverResult& solverResult, const DAENewtonResult& newtonResult) {
		solverResult.newton_iterations += newtonResult.iterations;
		solverResult.jacobian_evaluations += newtonResult.jacobian_evaluations;
		solverResult.last_step_newton_iterations = newtonResult.iterations;
		solverResult.max_step_newton_iterations = std::max(solverResult.max_step_newton_iterations, newtonResult.iterations);
		solverResult.final_residual_norm = newtonResult.residual_norm;
		solverResult.max_residual_norm = std::max(solverResult.max_residual_norm, newtonResult.residual_norm);
	}

	inline bool ValidateDAEInitialState(const IODESystemDAE& system, Real t,
		const Vector<Real>& x, const Vector<Real>& y, const DAESolverConfig& config,
		DAESolverResult& result) {
		Vector<Real> constraints(system.getAlgDim());
		system.algConstraints(t, x, y, constraints);
		result.final_constraint_norm = constraints.NormL2();
		result.max_constraint_violation = result.final_constraint_norm;
		if (!IsFiniteDAEVector(x) || !IsFiniteDAEVector(y) || !IsFiniteDAEVector(constraints)) {
			result.status = AlgorithmStatus::NumericalInstability;
			result.failure_reason = DAEFailureReason::NonFiniteState;
			result.error_message = "Non-finite initial DAE state or constraint residual";
			return false;
		}
		if (result.final_constraint_norm > config.constraint_tol) {
			result.status = AlgorithmStatus::InvalidInput;
			result.failure_reason = DAEFailureReason::InconsistentInitialConditions;
			result.error_message = "Initial algebraic constraints are inconsistent";
			return false;
		}
		return true;
	}

	inline void SetDAENewtonFailure(DAESolverResult& result, const DAENewtonResult& newton,
		const std::string& context) {
		result.failure_reason = newton.failure_reason;
		result.status = newton.failure_reason == DAEFailureReason::SingularJacobian
			? AlgorithmStatus::SingularMatrix : AlgorithmStatus::NumericalInstability;
		result.error_message = context + ": " + newton.message;
	}

	inline bool ValidateAcceptedDAEState(const Vector<Real>& x, const Vector<Real>& y,
		const Vector<Real>& constraints, const DAESolverConfig& config, DAESolverResult& result,
		const std::string& context) {
		result.final_constraint_norm = constraints.NormL2();
		result.max_constraint_violation = std::max(result.max_constraint_violation, result.final_constraint_norm);
		if (!IsFiniteDAEVector(x) || !IsFiniteDAEVector(y) || !IsFiniteDAEVector(constraints)) {
			result.status = AlgorithmStatus::NumericalInstability;
			result.failure_reason = DAEFailureReason::NonFiniteState;
			result.error_message = context + ": non-finite state or constraint residual";
			return false;
		}
		if (result.final_constraint_norm > config.constraint_tol) {
			result.status = AlgorithmStatus::NumericalInstability;
			result.failure_reason = DAEFailureReason::ConstraintViolation;
			result.error_message = context + ": algebraic constraint tolerance exceeded";
			return false;
		}
		return true;
	}

	/******************************************************************************/
	/*****            Consistent Initial Condition Solver                     *****/
	/******************************************************************************/

	/// @brief Compute consistent initial algebraic variables from initial differential state.
	/// 
	/// Given x₀ and an initial guess y₀, solve for y such that g(t₀, x₀, y) = 0.
	/// Uses Newton iteration on the constraint equations.
	///
	/// @param system DAE system with Jacobian
	/// @param t0 Initial time
	/// @param x0 Initial differential state (fixed)
	/// @param[in,out] y0 Initial algebraic state guess; modified to satisfy constraints
	/// @param max_iter Maximum Newton iterations (default: 20)
	/// @param tol Tolerance for constraint residual (default: 1e-10)
	/// @return True if consistent initial conditions were found
	struct ConsistentICResult {
		bool converged = false;
		int iterations = 0;
		int line_search_reductions = 0;
		Real initial_residual_norm = 0.0;
		Real final_residual_norm = 0.0;
		DAEFailureReason failure_reason = DAEFailureReason::None;
		std::string message;
	};

	inline ConsistentICResult ComputeConsistentICDetailed(const IODESystemDAEWithJacobian& system,
	                                Real t0, const Vector<Real>& x0, Vector<Real>& y0,
	                                int max_iter = 20, Real tol = PrecisionValues<Real>::IntegrationTolerance)
	{
		ConsistentICResult result;
		int algDim = system.getAlgDim();
		Vector<Real> g(algDim);
		Matrix<Real> dg_dy(algDim, algDim);

		for (int iter = 0; iter < max_iter; ++iter)
		{
			result.iterations = iter + 1;
			// Evaluate constraint residual
			system.algConstraints(t0, x0, y0, g);
			if (!IsFiniteDAEVector(g) || !IsFiniteDAEVector(y0)) {
				result.failure_reason = DAEFailureReason::NonFiniteState;
				result.message = "non-finite consistent-IC state or residual";
				return result;
			}

			// Check convergence
			Real g_norm = g.NormL2();
			if (iter == 0) result.initial_residual_norm = g_norm;
			result.final_residual_norm = g_norm;
			if (g_norm < tol) {
				result.converged = true;
				return result;
			}

			// Compute Jacobian ∂g/∂y
			system.jacobian_gy(t0, x0, y0, dg_dy);

			// Newton step: solve dg_dy * Δy = -g
			Vector<Real> rhs = g * (-1.0);
			try {
				Matrix<Real> A = dg_dy;
				GaussJordanSolver<Real>::SolveInPlace(A, rhs);
			}
			catch (const std::exception& error) {
				result.failure_reason = DAEFailureReason::HigherIndexSuspected;
				result.message = std::string("singular dg/dy while computing consistent IC: ") + error.what();
				return result;
			}

			Real damping = REAL(1.0);
			bool accepted = false;
			while (damping >= REAL(1.0 / 1024.0)) {
				Vector<Real> trial = y0 + rhs * damping;
				Vector<Real> trialResidual(algDim);
				system.algConstraints(t0, x0, trial, trialResidual);
				if (IsFiniteDAEVector(trialResidual) && trialResidual.NormL2() < g_norm) {
					y0 = trial;
					accepted = true;
					break;
				}
				damping *= REAL(0.5);
				result.line_search_reductions++;
			}
			if (!accepted) {
				result.failure_reason = DAEFailureReason::NewtonFailure;
				result.message = "consistent-IC Newton line search stalled";
				return result;
			}
		}

		result.failure_reason = DAEFailureReason::NewtonFailure;
		result.message = "consistent-IC iteration limit exceeded";
		return result;
	}

	inline bool ComputeConsistentIC(const IODESystemDAEWithJacobian& system,
	                                Real t0, const Vector<Real>& x0, Vector<Real>& y0,
	                                int max_iter = 20, Real tol = PrecisionValues<Real>::IntegrationTolerance)
	{
		return ComputeConsistentICDetailed(system, t0, x0, y0, max_iter, tol).converged;
	}

	/// @brief Verify that initial conditions satisfy constraints
	/// @param system DAE system
	/// @param t0 Initial time
	/// @param x0 Initial differential state
	/// @param y0 Initial algebraic state
	/// @param tol Tolerance for constraint residual
	/// @return True if |g(t0, x0, y0)| < tol
	inline bool VerifyConsistentIC(const IODESystemDAE& system,
	                               Real t0, const Vector<Real>& x0, const Vector<Real>& y0,
	                               Real tol = PrecisionValues<Real>::IntegrationTolerance)
	{
		int algDim = system.getAlgDim();
		Vector<Real> g(algDim);
		system.algConstraints(t0, x0, y0, g);
		return g.NormL2() < tol;
	}

} // namespace MML

#endif // MML_DAE_SOLVER_BASE_H

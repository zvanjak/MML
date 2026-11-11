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
#include <mml/algorithms/ImplicitNewtonMatrix.h>

#include <algorithm>
#include <cmath>
#include <string>
#include <utility>

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
		/// Requested step size; initial magnitude for adaptive methods such as RODAS4.
		/// BE/BDF2 retain a smaller step after Newton retries.
		Real step_size = 0.01;

		/// Maximum Newton iterations per step (default: 20)
		int max_newton_iter = 20;

		/// Newton convergence tolerance (precision-aware default)
		Real newton_tol = PrecisionValues<Real>::NewtonTolerance;

		/// Tolerance for constraint satisfaction g(t,x,y) ≈ 0 (precision-aware default)
		Real constraint_tol = PrecisionValues<Real>::IntegrationTolerance;

		/// Maximum steps; BE/BDF2 and RODAS4 count both accepted and rejected attempts.
		int max_steps = 100000;

		/// Absolute and relative tolerances for adaptive methods
		Real abs_tolerance = PrecisionValues<Real>::IntegrationTolerance;
		Real rel_tolerance = PrecisionValues<Real>::IntegrationTolerance * static_cast<MML::Real>(100.0);

		/// Adaptive bounds; min_step_size also bounds BE/BDF2 Newton retries.
		Real min_step_size = Precision::DAEMinStepSize;
		Real max_step_size = static_cast<MML::Real>(1.0);
		Real step_safety = static_cast<MML::Real>(0.9);

		/// Step-local simplified Newton; false refreshes the Jacobian for every correction.
		bool reuse_jacobian = true;

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
		NonFiniteState,
		StepSizeUnderflow,
		ResourceLimitExceeded
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
		int lu_factorizations = 0;
		int linear_solves = 0;
		int accepted_steps = 0;
		int rejected_steps = 0;
		int last_step_newton_iterations = 0;
		int max_step_newton_iterations = 0;
		Real final_residual_norm = 0.0;
		Real max_residual_norm = 0.0;
		Real final_constraint_norm = 0.0;
		Real max_constraint_violation = 0.0;  ///< Maximum |g(t,x,y)| seen
		Real last_error_norm = 0.0;		///< RODAS4 scaled embedded error, last completed attempt
		Real max_error_norm = 0.0;		///< RODAS4 maximum, including rejected attempts
		DAEFailureReason failure_reason = DAEFailureReason::None;

		DAESolverResult(Real t0, Real t_end, int diffDim, int algDim, int num_steps)
			: solution(t0, t_end, diffDim, algDim, num_steps) {}
	};

	inline bool PrepareDAETrajectoryPoint(DAESolverResult& result, int index) {
		if (detail::PrepareTrajectoryPoint(result.solution, index, result.status, result.error_message))
			return true;
		result.failure_reason = DAEFailureReason::ResourceLimitExceeded;
		return false;
	}

	inline bool IsFiniteDAEVector(const Vector<Real>& values) {
		for (int i = 0; i < values.size(); ++i)
			if (!std::isfinite(values[i])) return false;
		return true;
	}

	inline Real DAETimeTolerance(Real tEnd, Real stepSize) {
		return std::max(std::abs(stepSize) * Precision::MinimumStepRatio,
			static_cast<MML::Real>(10.0) * std::numeric_limits<Real>::epsilon() * std::max(static_cast<MML::Real>(1.0), std::abs(tEnd)));
	}

	inline void AssembleDAEAugmentedJacobian(const Matrix<Real>& df_dx, const Matrix<Real>& df_dy,
		const Matrix<Real>& dg_dx, const Matrix<Real>& dg_dy, Real derivativeCoefficient,
		Matrix<Real>& jacobian) {
		const int diffDim = df_dx.rows();
		const int algDim = dg_dy.rows();
		for (int i = 0; i < diffDim; ++i) {
			for (int j = 0; j < diffDim; ++j)
				jacobian(i, j) = (i == j ? static_cast<MML::Real>(1.0) : static_cast<MML::Real>(0.0)) - derivativeCoefficient * df_dx(i, j);
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
		int lu_factorizations = 0;
		int linear_solves = 0;
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
		Detail::ImplicitNewtonMatrix iterationMatrix;

		for (int iteration = 0; iteration < config.max_newton_iter; ++iteration) {
			result.iterations = iteration + 1;
			evaluateResidual(x, y, residual);
			if (!IsFiniteDAEVector(residual) || !IsFiniteDAEVector(x) || !IsFiniteDAEVector(y)) {
				result.failure_reason = DAEFailureReason::NonFiniteState;
				result.message = "non-finite DAE state or residual";
				return result;
			}

			result.residual_norm = residual.NormL2();
			if (!std::isfinite(result.residual_norm)) {
				result.failure_reason = DAEFailureReason::NonFiniteState;
				result.message = "non-finite DAE residual norm";
				return result;
			}
			if (result.residual_norm < config.newton_tol) {
				result.converged = true;
				return result;
			}

			Vector<Real> delta = residual * static_cast<MML::Real>(-1.0);
			try {
				if (iterationMatrix.NeedsRefresh(result.residual_norm, config.newton_tol,
					config.max_newton_iter - iteration - 1, config.reuse_jacobian)) {
					evaluateJacobian(x, y, jacobian);
					++result.jacobian_evaluations;
					iterationMatrix.Factor(jacobian);
					++result.lu_factorizations;
				}
				delta = iterationMatrix.Solve(delta);
				++result.linear_solves;
			}
			catch (const SingularMatrixError& error) {
				result.failure_reason = DAEFailureReason::SingularJacobian;
				result.message = error.what();
				return result;
			}
			catch (const MatrixNumericalError& error) {
				result.failure_reason = DAEFailureReason::NonFiniteState;
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
		solverResult.lu_factorizations += newtonResult.lu_factorizations;
		solverResult.linear_solves += newtonResult.linear_solves;
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

	inline bool ValidateDAETimeStep(Real t, Real nextTime, DAESolverResult& result,
		const std::string& context) {
		if (nextTime != t) return true;
		result.status = AlgorithmStatus::Stalled;
		result.failure_reason = DAEFailureReason::StepSizeUnderflow;
		result.error_message = context + ": step cannot advance floating-point time";
		return false;
	}

	namespace DAEImplicitDetail {
		inline DAESolverResult Solve(IODESystemDAEWithJacobian& system, Real t0,
			const Vector<Real>& x0, const Vector<Real>& y0, Real t_end,
			const DAESolverConfig& config, bool bdf2) {
			const std::string name = bdf2 ? "DAEBDF2" : "DAEBackwardEuler";
			const int diffDim = system.getDiffDim(), algDim = system.getAlgDim();
			if (!std::isfinite(t0) || !std::isfinite(t_end) || !std::isfinite(t_end - t0) ||
				t_end == t0 || !std::isfinite(config.step_size) || config.step_size <= 0 ||
				!std::isfinite(config.min_step_size) || config.min_step_size <= 0 ||
				!std::isfinite(config.newton_tol) || config.newton_tol <= 0 ||
				!std::isfinite(config.constraint_tol) || config.constraint_tol <= 0 ||
				config.max_steps <= 0 || config.max_newton_iter <= 0 ||
				diffDim <= 0 || algDim < 0 || x0.size() != diffDim || y0.size() != algDim)
				throw ArgumentError(name +
					": requires distinct finite endpoints, positive step/tolerances/budgets and matching state dimensions");

			AlgorithmTimer timer;
			DAESolverResult result(t0, t_end, diffDim, algDim, 0);
			result.algorithm_name = name;
			result.solution.fillValues(0, t0, x0, y0);
			if (!ValidateDAEInitialState(system, t0, x0, y0, config, result)) {
				result.elapsed_time_ms = timer.elapsed_ms();
				return result;
			}
			const Real direction = t_end > t0 ? Real(1) : Real(-1);
			Real t = t0, h = direction * config.step_size, previousH = 0;
			Vector<Real> x = x0, y = y0, previous = x0;
			Vector<Real> dxdt(diffDim), constraints(algDim);
			Matrix<Real> df_dx(diffDim, diffDim), df_dy(diffDim, algDim);
			Matrix<Real> dg_dx(algDim, diffDim), dg_dy(algDim, algDim);
			while (direction * (t_end - t) > 0 && result.accepted_steps + result.rejected_steps < config.max_steps) {
				if (!PrepareDAETrajectoryPoint(result, result.accepted_steps + 1)) break;
				h = direction * std::min(std::abs(h), std::abs(t_end - t));
				Real nextTime = t + h;
				if (!ValidateDAETimeStep(t, nextTime, result, name)) break;
				Vector<Real> history = x;
				Real coefficient = h;
				if (bdf2 && result.accepted_steps > 0) {
					// Variable-step BDF2; rejected attempts must not change the accepted history.
					Real ratio = h / previousH;
					history += (x - previous) * (ratio * ratio / (1 + 2 * ratio));
					coefficient *= (1 + ratio) / (1 + 2 * ratio);
				}
				system.diffEqs(t, x, y, dxdt);
				Vector<Real> trialX = x + dxdt * h, trialY = y;
				auto evaluateResidual = [&](const Vector<Real>& xx, const Vector<Real>& yy,
					Vector<Real>& residual) {
					system.diffEqs(nextTime, xx, yy, dxdt);
					system.algConstraints(nextTime, xx, yy, constraints);
					for (int i = 0; i < diffDim; ++i)
						residual[i] = xx[i] - history[i] - coefficient * dxdt[i];
					for (int i = 0; i < algDim; ++i)
						residual[diffDim + i] = constraints[i];
				};
				auto evaluateJacobian = [&](const Vector<Real>& xx, const Vector<Real>& yy,
					Matrix<Real>& jacobian) {
					system.allJacobians(nextTime, xx, yy, df_dx, df_dy, dg_dx, dg_dy);
					AssembleDAEAugmentedJacobian(df_dx, df_dy, dg_dx, dg_dy, coefficient, jacobian);
				};
				auto newton = SolveDAENewton(trialX, trialY, diffDim, algDim, config,
					evaluateResidual, evaluateJacobian);
				AccumulateDAENewtonDiagnostics(result, newton);
				if (!newton.converged) {
					++result.rejected_steps;
					result.solution.incrementRejectedSteps();
					if (newton.failure_reason == DAEFailureReason::NewtonFailure &&
						std::abs(h) > config.min_step_size) {
						h = direction * std::max(config.min_step_size, std::abs(h) / 2);
						continue;
					}
					SetDAENewtonFailure(result, newton, name + " at t=" + std::to_string(nextTime) +
						(newton.failure_reason == DAEFailureReason::NewtonFailure
							? " (minimum retry step reached)" : ""));
					break;
				}
				system.algConstraints(nextTime, trialX, trialY, constraints);
				if (!ValidateAcceptedDAEState(trialX, trialY, constraints, config, result, name + " step")) {
					++result.rejected_steps;
					result.solution.incrementRejectedSteps();
					break;
				}
				previous = x;
				previousH = h;
				x = std::move(trialX);
				y = std::move(trialY);
				t = nextTime;
				result.solution.fillValues(++result.accepted_steps, t, x, y);
				result.solution.incrementSuccessfulSteps();
			}
			result.total_steps = result.accepted_steps;
			system.algConstraints(t, x, y, constraints);
			result.final_constraint_norm = constraints.NormL2();
			if (direction * (t_end - t) > 0 && result.status == AlgorithmStatus::Success) {
				result.status = AlgorithmStatus::MaxIterationsExceeded;
				result.failure_reason = DAEFailureReason::MaxStepsExceeded;
				result.error_message = name + ": attempted-step budget exhausted at t=" + std::to_string(t);
			}
			result.elapsed_time_ms = timer.elapsed_ms();
			return result;
		}
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

			Real damping = static_cast<MML::Real>(1.0);
			bool accepted = false;
			while (damping >= static_cast<MML::Real>(1.0 / 1024.0)) {
				Vector<Real> trial = y0 + rhs * damping;
				Vector<Real> trialResidual(algDim);
				system.algConstraints(t0, x0, trial, trialResidual);
				if (IsFiniteDAEVector(trialResidual) && trialResidual.NormL2() < g_norm) {
					y0 = trial;
					accepted = true;
					break;
				}
				damping *= static_cast<MML::Real>(0.5);
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

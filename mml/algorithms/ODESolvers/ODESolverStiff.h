///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        ODESolverStiff.h                                                   ///
///  Description: Implicit methods for stiff ODE systems                              ///
///                                                                                   ///
///  Methods included:                                                                ///
///  - Backward Euler (1st order, A-stable, simple baseline)                         ///
///  - BDF2 (2nd order Backward Differentiation Formula, industry standard)          ///
///  - ROS2 (2-stage 2nd-order L-stable Rosenbrock, adaptive stepping)                 ///
///                                                                                   ///
///  REFERENCES:                                                                      ///
///    [HW2]  Hairer & Wanner, Solving ODEs II, Ch. III-IV                           ///
///    [G71]  Gear, C.W. (1971). Numerical Initial Value Problems in ODEs            ///
///    [SW87] Shampine & Watts (1979). Math. Software III, pp. 257-275               ///
///                                                                                   ///
///  All methods require IODESystemWithJacobian interface for Jacobian matrix.       ///
///  See references/book_references.md and references/paperes_references.md          ///
///                                                                                   ///
///  COMPLEXITY (per step, N = system dimension)                                      ///
///  ============================================                                      ///
///    Method          Order   Per-step cost        Note                              ///
///    ------          -----   -------------        ----                              ///
///    Backward Euler    1     O(N³)                Newton + LU solve per step        ///
///    BDF2              2     O(N³)                Newton + LU solve per step        ///
///    ROS2             2      O(N³)                LU factorization per step         ///
///                                                                                   ///
///  WARNING: Unlike explicit methods (Euler, RK4) which cost O(N) per step,         ///
///  implicit methods cost O(N³) per step due to Jacobian evaluation + LU solve.     ///
///  For N=100, this is ~10⁶× more expensive per step than explicit Euler.           ///
///  The trade-off: stiff problems may need 100-1000× fewer steps with implicit      ///
///  methods, making them faster overall despite the per-step cost.                   ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_ODE_SOLVER_STIFF_H
#define MML_ODE_SOLVER_STIFF_H

#include <mml/MMLBase.h>

#include <mml/base/Vector/Vector.h>
#include <mml/base/Matrix/Matrix.h>
#include <mml/base/ODESystemSolution.h>
#include <mml/interfaces/IODESystem.h>
#include <mml/base/AlgorithmTypes.h>
#include <mml/core/LinAlgEqSolvers.h>
#include <mml/algorithms/ImplicitNewtonMatrix.h>

#include <optional>
#include <algorithm>
#include <cmath>
#include <utility>

namespace MML {

	/******************************************************************************/
	/*****            Stiff Solver Configuration                              *****/
	/******************************************************************************/

	/// Configuration for stiff ODE solvers (Backward Euler, BDF2, Rosenbrock).
	///
	/// Provides user control over step size, Newton iteration parameters,
	/// and adaptive stepping bounds.
	///
	/// @example
	/// StiffSolverConfig config;
	/// config.step_size = 0.01;
	/// config.newton_tol = 1e-10;
	/// auto sol = SolveBackwardEuler(system, t0, y0, t_end, config);
	struct StiffSolverConfig {
		/// Requested step size; BE/BDF2 retain a smaller step after Newton retries.
		Real step_size = 0.01;

		/// Maximum Newton iterations per step (default: 10)
		int max_newton_iter = 10;

		/// Newton convergence tolerance (default: type-adaptive)
		Real newton_tol = Precision::NewtonTolerance;

		/// Error tolerance for adaptive methods (default: type-adaptive)
		Real tolerance = Precision::ODEDefaultTolerance;

		/// Minimum retry step size (a shorter final endpoint step is allowed).
		Real min_step_size = Precision::StiffMinStepSize;

		/// Maximum step size for adaptive methods (default: 0.1)
		Real max_step_size = 0.1;

		/// Maximum steps; BE/BDF2 count both accepted and rejected attempts.
		int max_steps = 100000;

		/// Reuse a step-local LU, refreshing on slow Newton convergence; false uses full Newton.
		bool reuse_jacobian = true;

		/// Factory method for high-precision configuration
		static StiffSolverConfig HighPrecision() {
			StiffSolverConfig config;
			config.newton_tol = Precision::ZeroTolerance;
			config.tolerance = Precision::RoundoffZeroThreshold;
			config.max_newton_iter = 20;
			return config;
		}

		/// Factory method for fast (lower precision) configuration
		static StiffSolverConfig Fast() {
			StiffSolverConfig config;
			config.newton_tol = Precision::NewtonTolerance * Real(100);
			config.tolerance = Precision::ODEDefaultTolerance * Real(100);
			config.max_newton_iter = 5;
			return config;
		}
	};

	/// @brief Result of stiff ODE integration with diagnostics
	struct StiffSolverResult {
		ODESystemSolution solution;		///< Solution trajectory

		// === Diagnostics ===
		std::string algorithm_name;		///< Algorithm used (BackwardEuler, BDF2, ROS2)
		AlgorithmStatus status = AlgorithmStatus::Success;
		std::string error_message;
		double elapsed_time_ms = 0.0;

		// === Statistics ===
		int total_steps = 0;
		int newton_iterations = 0;		///< Total Newton iterations
		int jacobian_evaluations = 0;	///< Number of Jacobian evaluations
		int function_evaluations = 0;	///< RHS evaluations, including rejected attempts
		int lu_factorizations = 0;		///< Successfully factored iteration matrices
		int linear_solves = 0;			///< Completed LU back-substitutions

		StiffSolverResult(Real t0, Real t_end, int dim, int num_steps)
			: solution(t0, t_end, dim, num_steps) {}
	};

	namespace ODEStiffDetail {
		inline StiffSolverResult SolveImplicit(IODESystemWithJacobian& system, Real t0,
			const Vector<Real>& y0, Real t_end, const StiffSolverConfig& config, bool bdf2) {
			AlgorithmTimer timer;
			StiffSolverResult result(t0, t_end, y0.size(), 0);
			result.algorithm_name = bdf2 ? "BDF2" : "BackwardEuler";
			if (!std::isfinite(t0) || !std::isfinite(t_end) || !std::isfinite(t_end - t0) ||
				t_end == t0 || !std::isfinite(config.step_size) || config.step_size <= 0 ||
				!std::isfinite(config.min_step_size) || config.min_step_size <= 0 ||
				!std::isfinite(config.newton_tol) || config.newton_tol <= 0 ||
				config.max_steps <= 0 || config.max_newton_iter <= 0 ||
				y0.size() == 0 || y0.size() != system.getDim()) {
				result.status = AlgorithmStatus::InvalidInput;
				result.error_message = result.algorithm_name +
					": requires distinct finite endpoints, positive step/tolerance/budgets and matching state dimensions";
				result.elapsed_time_ms = timer.elapsed_ms();
				return result;
			}
			auto isFinite = [](const Vector<Real>& values) {
				for (int i = 0; i < values.size(); ++i)
					if (!std::isfinite(values[i])) return false;
				return true;
			};
			result.solution.fillValues(0, t0, y0);
			if (!isFinite(y0)) {
				result.status = AlgorithmStatus::NumericalInstability;
				result.error_message = result.algorithm_name + ": non-finite initial state";
				result.elapsed_time_ms = timer.elapsed_ms();
				return result;
			}

			const Real direction = t_end > t0 ? Real(1) : Real(-1);
			Real t = t0, h = direction * config.step_size, previousH = 0;
			Vector<Real> y = y0, previous = y0, f(y0.size());
			Matrix<Real> J(y0.size(), y0.size());
			while (direction * (t_end - t) > 0 && result.solution.getTotalNumSteps() < config.max_steps) {
				if (!detail::PrepareTrajectoryPoint(result.solution, result.total_steps + 1,
					result.status, result.error_message)) break;
				h = direction * std::min(std::abs(h), std::abs(t_end - t));
				Real nextTime = t + h;
				if (direction * (nextTime - t) <= 0) {
					result.status = AlgorithmStatus::Stalled;
					result.error_message = result.algorithm_name + ": step cannot advance floating-point time";
					break;
				}
				Vector<Real> history = y;
				Vector<Real> trial = y;
				Real coefficient = h;
				if (bdf2 && result.total_steps > 0) {
					// Variable-step BDF2 (Hairer-Wanner II): r=h_n/h_{n-1}.
					Real ratio = h / previousH;
					history += (y - previous) * (ratio * ratio / (1 + 2 * ratio));
					coefficient *= (1 + ratio) / (1 + 2 * ratio);
					trial += (y - previous) * ratio;
				}

				bool converged = false;
				AlgorithmStatus failure = AlgorithmStatus::NumericalInstability;
				std::string message = "Newton iteration limit exceeded";
				bool retry = true;
				Detail::ImplicitNewtonMatrix iterationMatrix;
				for (int iteration = 0; iteration < config.max_newton_iter; ++iteration) {
					++result.newton_iterations;
					system.derivs(nextTime, trial, f);
					++result.function_evaluations;
					Vector<Real> residual = trial - history - f * coefficient;
					Real norm = residual.NormL2();
					if (!isFinite(trial) || !isFinite(f) || !std::isfinite(norm)) {
						message = "non-finite Newton state or residual";
						retry = false;
						break;
					}
					if (norm < config.newton_tol) {
						converged = true;
						break;
					}
					Vector<Real> delta;
					try {
						if (iterationMatrix.NeedsRefresh(norm, config.newton_tol,
							config.max_newton_iter - iteration - 1, config.reuse_jacobian)) {
							system.jacobian(nextTime, trial, J);
							++result.jacobian_evaluations;
							iterationMatrix.Factor(Matrix<Real>::Identity(y0.size()) - J * coefficient);
							++result.lu_factorizations;
						}
						delta = iterationMatrix.Solve(residual * Real(-1));
						++result.linear_solves;
					}
					catch (const SingularMatrixError& error) {
						failure = AlgorithmStatus::SingularMatrix;
						message = error.what();
						retry = false;
						break;
					}
					catch (const MatrixNumericalError& error) {
						message = error.what();
						retry = false;
						break;
					}
					trial += delta;
				}
				if (!converged) {
					result.solution.incrementRejectedSteps();
					if (!retry || std::abs(h) <= config.min_step_size) {
						result.status = failure;
						result.error_message = result.algorithm_name + " at t=" + std::to_string(nextTime) +
							": " + message + (retry ? " (minimum retry step reached)" : "");
						break;
					}
					h = direction * std::max(config.min_step_size, std::abs(h) / 2);
					continue;
				}
				previous = y;
				previousH = h;
				y = std::move(trial);
				t = nextTime;
				result.solution.fillValues(++result.total_steps, t, y);
				result.solution.incrementSuccessfulSteps();
			}
			if (direction * (t_end - t) > 0 && result.status == AlgorithmStatus::Success) {
				result.status = AlgorithmStatus::MaxIterationsExceeded;
				result.error_message = result.algorithm_name + ": attempted-step budget exhausted at t=" +
					std::to_string(t);
			}
			result.elapsed_time_ms = timer.elapsed_ms();
			return result;
		}
	}

	///////////////////////////////////////////////////////////////////////////////////////////
	//                           BACKWARD EULER METHOD                                        //
	///////////////////////////////////////////////////////////////////////////////////////////
	/// @brief Backward Euler implicit solver for stiff ODEs
	/// @details First-order implicit method: y_{n+1} = y_n + h*f(t_{n+1}, y_{n+1})
	///          Unconditionally stable (A-stable), good for very stiff systems.
	///          Uses Newton iteration to solve the implicit equation.
	/// Complexity: O(N³) per refreshed iteration matrix, O(N²) per reused Newton solve.
	/// @param system ODE system with Jacobian (must implement IODESystemWithJacobian)
	/// @param t0 Initial time
	/// @param y0 Initial condition vector
	/// @param t_end Final time
	/// @param h Requested step size (halved on Newton non-convergence)
	/// @param max_newton_iter Maximum Newton iterations per step (default: 10)
	/// @param newton_tol Newton convergence tolerance (default: 1e-8)
	/// @return ODESystemSolution containing solution trajectory
	inline ODESystemSolution SolveBackwardEuler(IODESystemWithJacobian& system, Real t0, const Vector<Real>& y0, Real t_end, Real h,
												int max_newton_iter = 10, Real newton_tol = Precision::NewtonTolerance) {
		StiffSolverConfig config;
		config.step_size = h;
		config.max_newton_iter = max_newton_iter;
		config.newton_tol = newton_tol;
		auto result = ODEStiffDetail::SolveImplicit(system, t0, y0, t_end, config, false);
		if (result.status == AlgorithmStatus::ResourceLimitExceeded)
			throw ResourceLimitError(result.error_message);
		if (result.status != AlgorithmStatus::Success)
			throw ODESolverError(result.error_message);
		return std::move(result.solution);
	}

	/// @brief Backward Euler solver with configuration object
	/// @param system ODE system with Jacobian
	/// @param t0 Initial time
	/// @param y0 Initial condition vector
	/// @param t_end Final time
	/// @param config Solver configuration
	/// @return StiffSolverResult with solution and diagnostics
	inline StiffSolverResult SolveBackwardEuler(IODESystemWithJacobian& system, Real t0, const Vector<Real>& y0, 
												Real t_end, const StiffSolverConfig& config) {
		return ODEStiffDetail::SolveImplicit(system, t0, y0, t_end, config, false);
	}

	///////////////////////////////////////////////////////////////////////////////////////////
	//                            BDF2 METHOD                                                 //
	///////////////////////////////////////////////////////////////////////////////////////////
	/// @brief Second-order Backward Differentiation Formula solver
	/// @details BDF2 formula: y_{n+1} = (4/3)*y_n - (1/3)*y_{n-1} + (2/3)*h*f(t_{n+1}, y_{n+1})
	///          A-stable, industry standard for moderately stiff problems.
	///          Uses Backward Euler for first step, then switches to BDF2.
	/// Complexity: O(N³) per step — Newton iteration with Jacobian + LU solve.
	///            Total: O(totalSteps × newtonIters × N³).
	/// @param system ODE system with Jacobian
	/// @param t0 Initial time
	/// @param y0 Initial condition vector
	/// @param t_end Final time
	/// @param h Requested step size (halved on Newton non-convergence)
	/// @param max_newton_iter Maximum Newton iterations per step (default: 10)
	/// @param newton_tol Newton convergence tolerance (default: 1e-8)
	/// @return ODESystemSolution containing solution trajectory
	inline ODESystemSolution SolveBDF2(IODESystemWithJacobian& system, Real t0, const Vector<Real>& y0, Real t_end, Real h,
									   int max_newton_iter = 10, Real newton_tol = Precision::NewtonTolerance) {
		StiffSolverConfig config;
		config.step_size = h;
		config.max_newton_iter = max_newton_iter;
		config.newton_tol = newton_tol;
		auto result = ODEStiffDetail::SolveImplicit(system, t0, y0, t_end, config, true);
		if (result.status == AlgorithmStatus::ResourceLimitExceeded)
			throw ResourceLimitError(result.error_message);
		if (result.status != AlgorithmStatus::Success)
			throw ODESolverError(result.error_message);
		return std::move(result.solution);
	}

	/// @brief BDF2 solver with configuration object
	/// @param system ODE system with Jacobian
	/// @param t0 Initial time
	/// @param y0 Initial condition vector
	/// @param t_end Final time
	/// @param config Solver configuration
	/// @return StiffSolverResult with solution and diagnostics
	inline StiffSolverResult SolveBDF2(IODESystemWithJacobian& system, Real t0, const Vector<Real>& y0, 
									   Real t_end, const StiffSolverConfig& config) {
		return ODEStiffDetail::SolveImplicit(system, t0, y0, t_end, config, true);
	}

	///////////////////////////////////////////////////////////////////////////////////////////
	//                         ROSENBROCK METHOD                                                //
	///////////////////////////////////////////////////////////////////////////////////////////
	/// @brief ROS2: 2-stage, 2nd-order, L-stable Rosenbrock (linearly implicit) method with adaptive stepping
	/// @details Semi-implicit Runge-Kutta method that uses the exact Jacobian but no Newton
	///          iteration. This is the method of Verwer, Spee, Blom & Hundsdorfer (1999),
	///          γ = 1 + 1/√2, with the embedded 1st-order solution for error estimation.
	///          Well-suited for moderately stiff problems.
	///          (Until analysis v4 this class was called Rosenbrock23Solver; it never was the
	///          3-stage 3rd-order Rosenbrock23 of Shampine & Reichelt / MATLAB ode23s.)
	/// Complexity: O(N³) per step — one LU factorization O(N³) + two back-substitutions O(N²)
	///            + Jacobian evaluation O(N² × f-eval cost). More efficient than BDF
	///            methods for small-to-medium N because no Newton iteration needed.
	///          Formula: (I - h*γ*J)*k_i = h*f(...) + Σ(γ_ij * k_j) for each stage
	///
	/// The key advantage of Rosenbrock methods is that they only require ONE
	/// Jacobian evaluation and LU decomposition per step (not per stage).
	///
	/// References:
	/// - Verwer, Spee, Blom, Hundsdorfer, "A second-order Rosenbrock method applied to
	///   photochemical dispersion problems", SIAM J. Sci. Comput. 20 (1999)
	/// - Hairer & Wanner, "Solving ODEs II" (1996), Section IV.7
	class ROS2Solver {
	private:
		IODESystemWithJacobian& _system;
		Real _abs_tol;
		Real _rel_tol;
		int _functionEvaluations = 0, _jacobianEvaluations = 0;
		int _luFactorizations = 0, _linearSolves = 0;

		// ROS2 coefficients - 2-stage, 2nd order, L-stable method
		// γ = 1 + 1/sqrt(2) ≈ 1.707 for L-stability
		const Real gamma_ros = 1.7071067811865475; // 1.0 + 0.5 * sqrt(2)

	public:
		ROS2Solver(IODESystemWithJacobian& system, Real abs_tol = Precision::ODEDefaultTolerance, Real rel_tol = Precision::ODEDefaultTolerance, int /*max_newton_iter*/ = 10,
						   Real /*newton_tol*/ = Precision::NewtonTolerance)
			: _system(system)
			, _abs_tol(abs_tol)
			, _rel_tol(rel_tol) {}

		/// @brief Take one Rosenbrock step with error estimation
		/// @param t Current time
		/// @param y Current solution
		/// @param h Step size to attempt
		/// @param y_out Output: 2nd order solution
		/// @param y_err Output: Error estimate
		/// @return true if step successful, false if linear solve failed
		bool Step(Real t, const Vector<Real>& y, Real h, Vector<Real>& y_out, Vector<Real>& y_err) {
			const int dim = y.size();
			if (dim <= 0)
				throw ODESolverError("ROS2Solver::Step - empty state vector");

			Vector<Real> dydt(dim);
			Matrix<Real> J(dim, dim);

			// Evaluate f and Jacobian at (t, y)
			_system.derivs(t, y, dydt);
			++_functionEvaluations;
			_system.derivsAndJacobian(t, y, dydt, J);
			++_functionEvaluations;
			++_jacobianEvaluations;

			// Form W = I - h*γ*J once and factor it once; both stages reuse the factorization
			Matrix<Real> W = Matrix<Real>::Identity(dim) - J * (h * gamma_ros);
			Vector<Real> k1, k2;
			try {
				LUSolver<Real> lu(W);
				++_luFactorizations;

				// Stage 1: W*k1 = h*f(t, y)
				k1 = lu.Solve(dydt * h);
				++_linearSolves;

				// Stage 2: W*k2 = h*f(t + h, y + k1) - 2*k1
				Vector<Real> y2 = y + k1;
				Vector<Real> f2(dim);
				_system.derivs(t + h, y2, f2);
				++_functionEvaluations;
				k2 = lu.Solve(f2 * h - k1 * 2.0);
				++_linearSolves;
			} catch (const SingularMatrixError&) {
				return false; // W singular: let the driver shrink the step
			}

			// 2nd order solution: y_new = y + 1.5*k1 + 0.5*k2
			y_out = y + k1 * 1.5 + k2 * 0.5;

			// 1st order solution: y1 = y + k1
			// Error = y_new - y1 = 0.5*k1 + 0.5*k2
			y_err = (k1 + k2) * 0.5;

			return true;
		}

		/// @brief Solve ODE system with adaptive stepping
		/// @param t0 Initial time
		/// @param y0 Initial condition
		/// @param t_end Final time
		/// @param h_init Initial step size
		/// @return Solution trajectory
	private:
		bool SolveInto(ODESystemSolution& sol, Real t0, const Vector<Real>& y0,
			Real t_end, Real h_init, int max_steps) {
			_functionEvaluations = _jacobianEvaluations = _luFactorizations = _linearSolves = 0;
			if (!std::isfinite(t_end - t0) || t_end == t0)
				throw ODESolverError("ROS2Solver::Solve: requires distinct finite endpoints");
			if (h_init <= 0)
				throw ODESolverError("ROS2Solver::Solve: initial step size must be positive");

			int dim = y0.size();
			Real t = t0;
			Vector<Real> y = y0;
			const Real direction = t_end > t0 ? Real(1) : Real(-1);
			Real h = direction * h_init;

			// Safety factors for step size control
			const Real safety = 0.9;
			const Real h_min = Precision::StiffMinStepSize;		   // Minimum step size
			const Real h_max = h_init * 100.0; // Maximum step size
			const Real fac_max = 5.0;		   // Maximum step increase factor
			const Real fac_min = 0.2;		   // Minimum step decrease factor

			sol.fillValues(0, t, y);
			int step = 1;

			while (direction * (t_end - t) > 0 && sol.getTotalNumSteps() < max_steps) {
				sol.reserveForPoint(step);
				Real h_actual = direction * std::min(std::abs(h), std::abs(t_end - t));

				// Don't take tiny final steps
				if (direction * (t + h_actual * 1.01 - t_end) >= 0) {
					h_actual = t_end - t;
				}
				if (t + h_actual == t)
					throw ODESolverError("ROS2 step cannot advance floating-point time");

				Vector<Real> y_new(dim);
				Vector<Real> y_err(dim);

				bool success = Step(t, y, h_actual, y_new, y_err);

				if (!success) {
					// Linear solve failed, reduce step size
					h = h_actual * 0.5;
					if (std::abs(h) < h_min) {
						throw ODESolverError("Rosenbrock: Linear solve failed, system may be singular");
					}
					sol.incrementRejectedSteps();
					continue;
				}

				// Compute error norm (mixed absolute/relative)
				Real err_norm = 0.0;
				for (int i = 0; i < dim; ++i) {
					Real scale = _abs_tol + _rel_tol * std::max(std::abs(y[i]), std::abs(y_new[i]));
					Real err_scaled = std::abs(y_err[i]) / scale;
					err_norm = std::max(err_norm, err_scaled);
				}

				// Avoid division by zero
				err_norm = std::max<Real>(err_norm, Precision::ErrorNormFloor);

				// Accept or reject step
				if (err_norm <= 1.0) {
					// Accept step
					t += h_actual;
					y = y_new;
					sol.fillValues(step, t, y);
					sol.incrementSuccessfulSteps();
					++step;

					// Compute new step size (PI controller formula for 2nd order method)
					// h_new = h * safety * (1/err_norm)^(1/3)
					Real factor = safety * std::pow(1.0 / err_norm, 1.0 / 3.0);
					factor = std::min(fac_max, std::max(fac_min, factor));
					h = direction * std::min(h_max, std::abs(h_actual) * factor);
				} else {
					// Reject and retry with smaller step
					Real factor = safety * std::pow(1.0 / err_norm, 0.5); // More aggressive for rejection
					factor = std::max(fac_min, factor);
					h = h_actual * factor;

					if (std::abs(h) < h_min) {
						throw ODESolverError("Rosenbrock: Step size too small (" + std::to_string(h) +
												 "), system may be too stiff or tolerance too tight");
					}
					sol.incrementRejectedSteps();
				}
			}

			sol.setFinalSize(step - 1);
			return direction * (t_end - t) <= 0;
		}

	public:
		ODESystemSolution Solve(Real t0, const Vector<Real>& y0, Real t_end, Real h_init) {
			ODESystemSolution solution(t0, t_end, y0.size(), 0);
			if (!SolveInto(solution, t0, y0, t_end, h_init, 50000))
				throw ODESolverError("Rosenbrock: Maximum steps exceeded");
			return solution;
		}

		/// @brief Solve stiff ODE system using configuration object
		/// @param y0 Initial condition vector
		/// @param t0 Initial time
		/// @param t_end Final time
		/// @param config Solver configuration
		/// @return StiffSolverResult with solution and diagnostics
		StiffSolverResult Solve(const Vector<Real>& y0, Real t0, Real t_end, const StiffSolverConfig& config) {
			AlgorithmTimer timer;  // Starts automatically

			int dim = y0.size();
			StiffSolverResult result(t0, t_end, dim, 0);
			result.algorithm_name = "ROS2";
			if (config.max_steps <= 0) {
				result.status = AlgorithmStatus::InvalidInput;
				result.error_message = "ROS2 requires a positive attempted-step budget";
				return result;
			}

			// Update internal tolerances
			_abs_tol = config.tolerance;
			_rel_tol = config.tolerance;

			try {
				Real h_init = (config.step_size > 0) ? config.step_size : std::abs(t_end - t0) / 100.0;
				if (!SolveInto(result.solution, t0, y0, t_end, h_init, config.max_steps)) {
					result.status = AlgorithmStatus::MaxIterationsExceeded;
					result.error_message = "ROS2: attempted-step budget exhausted";
				}
			} catch (const ResourceLimitError& error) {
				result.status = AlgorithmStatus::ResourceLimitExceeded;
				result.error_message = error.what();
			} catch (const std::exception& e) {
				result.status = AlgorithmStatus::NumericalInstability;
				result.error_message = e.what();
			}

			result.elapsed_time_ms = timer.elapsed_ms();
			result.total_steps = result.solution.getNumStepsOK();
			result.function_evaluations = _functionEvaluations;
			result.jacobian_evaluations = _jacobianEvaluations;
			result.lu_factorizations = _luFactorizations;
			result.linear_solves = _linearSolves;
			return result;
		}
	};


	///////////////////////////////////////////////////////////////////////////////////////////
	// ODEStiffConfig - Configuration for stiff ODE detailed APIs
	///////////////////////////////////////////////////////////////////////////////////////////
	struct ODEStiffConfig : public EvaluationConfigBase {
		/// Embedded stiff solver configuration (step sizes, Newton params, etc.)
		StiffSolverConfig solver;
	};

	///////////////////////////////////////////////////////////////////////////////////////////
	// ODEStiffResult - Result type for stiff ODE detailed APIs
	///////////////////////////////////////////////////////////////////////////////////////////
	struct ODEStiffResult : public EvaluationResultBase {
		std::optional<ODESystemSolution> solution;

		int total_steps = 0;
		int newton_iterations = 0;
		int jacobian_evaluations = 0;
		int lu_factorizations = 0;
		int linear_solves = 0;
	};

	///////////////////////////////////////////////////////////////////////////////////////////
	// ODEStiffDetail - Internal helpers for Detailed API execution
	///////////////////////////////////////////////////////////////////////////////////////////
	namespace ODEStiffDetail
	{
		template<typename ResultType, typename ComputeFn>
		ResultType ExecuteODEStiffDetailed(const char* algorithm_name,
		                                   const EvaluationConfigBase& config,
		                                   ComputeFn&& compute)
		{
			auto execute = [&]() {
				AlgorithmTimer timer;
				ResultType result = MakeEvaluationSuccessResult<ResultType>(algorithm_name);
				compute(result);
				result.elapsed_time_ms = timer.elapsed_ms();
				return result;
			};

			if (config.exception_policy == EvaluationExceptionPolicy::Propagate)
				return execute();

			try {
				return execute();
			}
			catch (const ResourceLimitError& ex) {
				return MakeEvaluationFailureResult<ResultType>(
					AlgorithmStatus::ResourceLimitExceeded, ex.what(), algorithm_name);
			}
			catch (const ODESolverError& ex) {
				return MakeEvaluationFailureResult<ResultType>(
					AlgorithmStatus::NumericalInstability, ex.what(), algorithm_name);
			}
			catch (const std::exception& ex) {
				return MakeEvaluationFailureResult<ResultType>(
					AlgorithmStatus::AlgorithmSpecificFailure, ex.what(), algorithm_name);
			}
		}
	} // namespace ODEStiffDetail

	///////////////////////////////////////////////////////////////////////////////////////////
	// Detailed free function API — Stiff Solvers
	///////////////////////////////////////////////////////////////////////////////////////////

	/// Backward Euler solve with full instrumentation.
	inline ODEStiffResult SolveBackwardEulerDetailed(
		IODESystemWithJacobian& system, Real t0, const Vector<Real>& y0, Real t_end,
		const ODEStiffConfig& config = {})
	{
		return ODEStiffDetail::ExecuteODEStiffDetailed<ODEStiffResult>(
			"BackwardEuler", config,
			[&](ODEStiffResult& result) {
				auto inner = SolveBackwardEuler(system, t0, y0, t_end, config.solver);
				result.solution             = std::move(inner.solution);
				result.total_steps          = inner.total_steps;
				result.newton_iterations    = inner.newton_iterations;
				result.jacobian_evaluations = inner.jacobian_evaluations;
				result.function_evaluations = inner.function_evaluations;
				result.lu_factorizations    = inner.lu_factorizations;
				result.linear_solves        = inner.linear_solves;

				if (inner.status != AlgorithmStatus::Success) {
					result.status        = inner.status;
					result.error_message = inner.error_message;
				}
			});
	}

	/// BDF2 solve with full instrumentation.
	inline ODEStiffResult SolveBDF2Detailed(
		IODESystemWithJacobian& system, Real t0, const Vector<Real>& y0, Real t_end,
		const ODEStiffConfig& config = {})
	{
		return ODEStiffDetail::ExecuteODEStiffDetailed<ODEStiffResult>(
			"BDF2", config,
			[&](ODEStiffResult& result) {
				auto inner = SolveBDF2(system, t0, y0, t_end, config.solver);
				result.solution             = std::move(inner.solution);
				result.total_steps          = inner.total_steps;
				result.newton_iterations    = inner.newton_iterations;
				result.jacobian_evaluations = inner.jacobian_evaluations;
				result.function_evaluations = inner.function_evaluations;
				result.lu_factorizations    = inner.lu_factorizations;
				result.linear_solves        = inner.linear_solves;

				if (inner.status != AlgorithmStatus::Success) {
					result.status        = inner.status;
					result.error_message = inner.error_message;
				}
			});
	}

	/// ROS2 solve with full instrumentation.
	inline ODEStiffResult SolveROS2Detailed(
		IODESystemWithJacobian& system, const Vector<Real>& y0, Real t0, Real t_end,
		const ODEStiffConfig& config = {})
	{
		return ODEStiffDetail::ExecuteODEStiffDetailed<ODEStiffResult>(
			"ROS2", config,
			[&](ODEStiffResult& result) {
				ROS2Solver solver(system);
				auto inner = solver.Solve(y0, t0, t_end, config.solver);
				result.solution             = std::move(inner.solution);
				result.total_steps          = inner.total_steps;
				result.newton_iterations    = inner.newton_iterations;
				result.jacobian_evaluations = inner.jacobian_evaluations;
				result.function_evaluations = inner.function_evaluations;
				result.lu_factorizations    = inner.lu_factorizations;
				result.linear_solves        = inner.linear_solves;

				if (inner.status != AlgorithmStatus::Success) {
					result.status        = inner.status;
					result.error_message = inner.error_message;
				}
			});
	}

	/// @deprecated The class was renamed in analysis v4: it implements the 2-stage 2nd-order ROS2
	///             method, not the 3-stage Rosenbrock23 (ode23s). Use ROS2Solver.
	using Rosenbrock23Solver [[deprecated("Rosenbrock23Solver implements ROS2 (2-stage, 2nd order); use ROS2Solver")]] = ROS2Solver;

	/// @deprecated Use SolveROS2Detailed.
	[[deprecated("SolveRosenbrock23Detailed was renamed SolveROS2Detailed (the method is ROS2, 2nd order)")]]
	inline ODEStiffResult SolveRosenbrock23Detailed(
		IODESystemWithJacobian& system, const Vector<Real>& y0, Real t0, Real t_end,
		const ODEStiffConfig& config = {})
	{
		return SolveROS2Detailed(system, y0, t0, t_end, config);
	}

} // namespace MML

#endif // MML_ODE_SOLVER_STIFF_H

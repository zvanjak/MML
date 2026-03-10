///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        ODEStiffSolvers.h                                                   ///
///  Description: Implicit methods for stiff ODE systems                              ///
///                                                                                   ///
///  Methods included:                                                                ///
///  - Backward Euler (1st order, A-stable, simple baseline)                         ///
///  - BDF2 (2nd order Backward Differentiation Formula, industry standard)          ///
///  - Rosenbrock23 (2nd/3rd order semi-implicit with adaptive stepping)             ///
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
///    Rosenbrock23     2/3    O(N³)                LU factorization per step         ///
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
#if !defined MML_ODE_STIFF_SOLVERS_H
#define MML_ODE_STIFF_SOLVERS_H

#include "mml/MMLBase.h"

#include "mml/base/Vector/Vector.h"
#include "mml/base/Matrix/Matrix.h"
#include "mml/base/ODESystemSolution.h"
#include "mml/interfaces/IODESystem.h"
#include "mml/core/AlgorithmTypes.h"
#include "mml/core/LinAlgEqSolvers.h"

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
		/// Step size for fixed-step methods (default: 0.01)
		Real step_size = 0.01;

		/// Maximum Newton iterations per step (default: 10)
		int max_newton_iter = 10;

		/// Newton convergence tolerance (default: 1e-8)
		Real newton_tol = 1e-8;

		/// Error tolerance for adaptive methods (default: 1e-6)
		Real tolerance = 1e-6;

		/// Minimum step size for adaptive methods (default: 1e-12)
		Real min_step_size = 1e-12;

		/// Maximum step size for adaptive methods (default: 0.1)
		Real max_step_size = 0.1;

		/// Maximum number of steps (default: 100000)
		int max_steps = 100000;

		/// Factory method for high-precision configuration
		static StiffSolverConfig HighPrecision() {
			StiffSolverConfig config;
			config.newton_tol = 1e-12;
			config.tolerance = 1e-10;
			config.max_newton_iter = 20;
			return config;
		}

		/// Factory method for fast (lower precision) configuration
		static StiffSolverConfig Fast() {
			StiffSolverConfig config;
			config.newton_tol = 1e-6;
			config.tolerance = 1e-4;
			config.max_newton_iter = 5;
			return config;
		}
	};

	/// @brief Result of stiff ODE integration with diagnostics
	struct StiffSolverResult {
		ODESystemSolution solution;		///< Solution trajectory

		// === Diagnostics ===
		std::string algorithm_name;		///< Algorithm used (BackwardEuler, BDF2, Rosenbrock23)
		AlgorithmStatus status = AlgorithmStatus::Success;
		std::string error_message;
		double elapsed_time_ms = 0.0;

		// === Statistics ===
		int total_steps = 0;
		int newton_iterations = 0;		///< Total Newton iterations
		int jacobian_evaluations = 0;	///< Number of Jacobian evaluations

		StiffSolverResult(Real t0, Real t_end, int dim, int num_steps)
			: solution(t0, t_end, dim, num_steps) {}
	};

	///////////////////////////////////////////////////////////////////////////////////////////
	//                           BACKWARD EULER METHOD                                        //
	///////////////////////////////////////////////////////////////////////////////////////////
	/// @brief Backward Euler implicit solver for stiff ODEs
	/// @details First-order implicit method: y_{n+1} = y_n + h*f(t_{n+1}, y_{n+1})
	///          Unconditionally stable (A-stable), good for very stiff systems.
	///          Uses Newton iteration to solve the implicit equation.
	/// Complexity: O(N³) per step — each step requires Jacobian evaluation O(N²)
	///            + LU factorization O(N³) + back-substitution O(N²) per Newton iteration.
	/// @param system ODE system with Jacobian (must implement IODESystemWithJacobian)
	/// @param t0 Initial time
	/// @param y0 Initial condition vector
	/// @param t_end Final time
	/// @param h Step size (fixed)
	/// @param max_newton_iter Maximum Newton iterations per step (default: 10)
	/// @param newton_tol Newton convergence tolerance (default: 1e-8)
	/// @return ODESystemSolution containing solution trajectory
	inline ODESystemSolution SolveBackwardEuler(IODESystemWithJacobian& system, Real t0, const Vector<Real>& y0, Real t_end, Real h,
												int max_newton_iter = 10, Real newton_tol = 1e-8) {
		int dim = y0.size();
		int num_steps = static_cast<int>((t_end - t0) / h) + 1;

		ODESystemSolution sol(t0, t_end, dim, num_steps);

		Real t = t0;
		Vector<Real> y = y0;
		Vector<Real> dydt(dim);
		Matrix<Real> J(dim, dim);

		// Save initial condition
		sol.fillValues(0, t, y);
		int step = 1;

		while (t < t_end && step < num_steps) {
			Real h_actual = std::min(h, t_end - t);
			Real t_next = t + h_actual;

			// Newton iteration to solve: y_new - y - h*f(t_next, y_new) = 0
			Vector<Real> y_new = y; // Initial guess

			bool converged = false;
			for (int iter = 0; iter < max_newton_iter; ++iter) {
				// Evaluate f and Jacobian at current guess
				system.derivs(t_next, y_new, dydt);
				system.jacobian(t_next, y_new, dydt, J);

				// Residual: r = y_new - y - h*f(t_next, y_new)
				Vector<Real> residual = y_new - y - dydt * h_actual;

				// Check convergence
				Real res_norm = residual.NormL2();
				if (res_norm < newton_tol) {
					converged = true;
					break;
				}

				// Newton step: solve (I - h*J)*delta = -residual
				Matrix<Real> A = Matrix<Real>::Identity(dim) - J * h_actual;
				Vector<Real> rhs = residual * (-1.0);
				GaussJordanSolver<Real>::SolveInPlace(A, rhs);
				Vector<Real> delta = rhs; // Solution is in rhs after SolveInPlace

				// Update guess
				y_new = y_new + delta;
			}

			if (!converged) {
				throw std::runtime_error("Backward Euler: Newton iteration failed to converge at t=" + std::to_string(t_next));
			}

			// Accept step
			y = y_new;
			t = t_next;
			sol.fillValues(step, t, y);
			sol.incrementSuccessfulSteps();
			++step;
		}

		// step is now one past the last filled index
		sol.setFinalSize(step - 1);
		return sol;
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
		AlgorithmTimer timer;  // Starts automatically

		int dim = y0.size();
		int num_steps = static_cast<int>((t_end - t0) / config.step_size) + 1;
		StiffSolverResult result(t0, t_end, dim, num_steps);
		result.algorithm_name = "BackwardEuler";

		try {
			result.solution = SolveBackwardEuler(system, t0, y0, t_end, config.step_size, 
												 config.max_newton_iter, config.newton_tol);
			result.total_steps = result.solution.getNumStepsOK();
		} catch (const std::exception& e) {
			result.status = AlgorithmStatus::NumericalInstability;
			result.error_message = e.what();
		}

		result.elapsed_time_ms = timer.elapsed_ms();
		return result;
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
	/// @param h Step size (fixed)
	/// @param max_newton_iter Maximum Newton iterations per step (default: 10)
	/// @param newton_tol Newton convergence tolerance (default: 1e-8)
	/// @return ODESystemSolution containing solution trajectory
	inline ODESystemSolution SolveBDF2(IODESystemWithJacobian& system, Real t0, const Vector<Real>& y0, Real t_end, Real h,
									   int max_newton_iter = 10, Real newton_tol = 1e-8) {
		int dim = y0.size();
		int num_steps = static_cast<int>((t_end - t0) / h) + 1;

		ODESystemSolution sol(t0, t_end, dim, num_steps);

		Real t = t0;
		Vector<Real> y = y0;
		Vector<Real> y_prev(dim); // For BDF2, we need previous step
		Vector<Real> dydt(dim);
		Matrix<Real> J(dim, dim);

		// Save initial condition
		sol.fillValues(0, t, y);
		int step = 1;

		// First step: Use Backward Euler (BDF1) to bootstrap
		if (step < num_steps) {
			Real h_actual = std::min(h, t_end - t);
			Real t_next = t + h_actual;

			Vector<Real> y_new = y;
			bool converged = false;

			for (int iter = 0; iter < max_newton_iter; ++iter) {
				system.derivs(t_next, y_new, dydt);
				system.jacobian(t_next, y_new, dydt, J);

				Vector<Real> residual = y_new - y - dydt * h_actual;
				Real res_norm = residual.NormL2();

				if (res_norm < newton_tol) {
					converged = true;
					break;
				}

				Matrix<Real> A = Matrix<Real>::Identity(dim) - J * h_actual;
				Vector<Real> rhs = residual * (-1.0);
				GaussJordanSolver<Real>::SolveInPlace(A, rhs);
				Vector<Real> delta = rhs;
				y_new = y_new + delta;
			}

			if (!converged) {
				throw std::runtime_error("BDF2: Bootstrap step failed at t=" + std::to_string(t_next));
			}

			y_prev = y;
			y = y_new;
			t = t_next;
			sol.fillValues(step, t, y);
			sol.incrementSuccessfulSteps();
			++step;
		}

		// Subsequent steps: Use BDF2
		while (t < t_end && step < num_steps) {
			Real h_actual = std::min(h, t_end - t);
			Real t_next = t + h_actual;

			// BDF2: y_{n+1} = (4/3)*y_n - (1/3)*y_{n-1} + (2/3)*h*f(t_{n+1}, y_{n+1})
			// Rearrange: y_{n+1} - (2/3)*h*f(t_{n+1}, y_{n+1}) = (4/3)*y_n - (1/3)*y_{n-1}

			Vector<Real> y_new = y * 2.0 - y_prev; // Extrapolation as initial guess
			Vector<Real> rhs = y * (4.0 / 3.0) - y_prev * (1.0 / 3.0);

			bool converged = false;
			for (int iter = 0; iter < max_newton_iter; ++iter) {
				system.derivs(t_next, y_new, dydt);
				system.jacobian(t_next, y_new, dydt, J);

				// Residual: r = y_new - (2/3)*h*f - rhs
				Vector<Real> residual = y_new - dydt * (2.0 / 3.0 * h_actual) - rhs;
				Real res_norm = residual.NormL2();

				if (res_norm < newton_tol) {
					converged = true;
					break;
				}

				// Newton step: solve (I - (2/3)*h*J)*delta = -residual
				Matrix<Real> A = Matrix<Real>::Identity(dim) - J * (2.0 / 3.0 * h_actual);
				Vector<Real> rhs = residual * (-1.0);
				GaussJordanSolver<Real>::SolveInPlace(A, rhs);
				Vector<Real> delta = rhs;
				y_new = y_new + delta;
			}

			if (!converged) {
				throw std::runtime_error("BDF2: Newton iteration failed at t=" + std::to_string(t_next));
			}

			// Accept step
			y_prev = y;
			y = y_new;
			t = t_next;
			sol.fillValues(step, t, y);
			sol.incrementSuccessfulSteps();
			++step;
		}

		sol.setFinalSize(step - 1);
		return sol;
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
		AlgorithmTimer timer;  // Starts automatically

		int dim = y0.size();
		int num_steps = static_cast<int>((t_end - t0) / config.step_size) + 1;
		StiffSolverResult result(t0, t_end, dim, num_steps);
		result.algorithm_name = "BDF2";

		try {
			result.solution = SolveBDF2(system, t0, y0, t_end, config.step_size, 
										config.max_newton_iter, config.newton_tol);
			result.total_steps = result.solution.getNumStepsOK();
		} catch (const std::exception& e) {
			result.status = AlgorithmStatus::NumericalInstability;
			result.error_message = e.what();
		}

		result.elapsed_time_ms = timer.elapsed_ms();
		return result;
	}

	///////////////////////////////////////////////////////////////////////////////////////////
	//                         ROSENBROCK METHOD                                                //
	///////////////////////////////////////////////////////////////////////////////////////////
	/// @brief Rosenbrock semi-implicit method with adaptive stepping
	/// @details Semi-implicit Runge-Kutta method that uses approximate Jacobian.
	///          Uses a 2-stage L-stable Rosenbrock method (ROS2) with embedded
	///          error estimation. Well-suited for moderately stiff problems.
	/// Complexity: O(N³) per step — one LU factorization O(N³) + back-substitutions O(N²)
	///            + Jacobian evaluation O(N² × f-eval cost). More efficient than BDF
	///            methods for small-to-medium N because no Newton iteration needed.
	///          Formula: (I - h*γ*J)*k_i = h*f(...) + Σ(γ_ij * k_j) for each stage
	///
	/// The key advantage of Rosenbrock methods is that they only require ONE
	/// Jacobian evaluation and LU decomposition per step (not per stage).
	///
	/// References:
	/// - Shampine, "Implementation of Rosenbrock methods" (1982)
	/// - Hairer & Wanner, "Solving ODEs II" (1996), Section IV.7
	class Rosenbrock23Solver {
	private:
		IODESystemWithJacobian& _system;
		Real _abs_tol;
		Real _rel_tol;

		// ROS2 coefficients - 2-stage, 2nd order, L-stable method
		// γ = 1 + 1/sqrt(2) ≈ 1.707 for L-stability
		const Real gamma_ros = 1.7071067811865475; // 1.0 + 0.5 * sqrt(2)

	public:
		Rosenbrock23Solver(IODESystemWithJacobian& system, Real abs_tol = 1e-6, Real rel_tol = 1e-6, int /*max_newton_iter*/ = 10,
						   Real /*newton_tol*/ = 1e-8)
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

			// DEBUG: Check dimension
			if (dim <= 0 || dim > 1000) {
				throw std::runtime_error("Rosenbrock::Step - invalid dimension: " + std::to_string(dim));
			}

			Vector<Real> dydt(dim);
			Matrix<Real> J(dim, dim);

			// Evaluate f and Jacobian at (t, y)
			_system.derivs(t, y, dydt);
			_system.jacobian(t, y, dydt, J);

			// Form matrix W = I - h*γ*J (same for all stages)
			Matrix<Real> W = Matrix<Real>::Identity(dim) - J * (h * gamma_ros);

			// Stage 1: Solve W*k1 = h*f(t, y)
			Vector<Real> k1 = dydt * h;
			Matrix<Real> W_copy = W;
			try {
				GaussJordanSolver<Real>::SolveInPlace(W_copy, k1);
			} catch (...) {
				return false; // Linear solve failed
			}

			// Stage 2: Solve W*k2 = h*f(t + h, y + k1) - 2*k1
			Vector<Real> y2 = y + k1;
			Vector<Real> f2(dim);
			_system.derivs(t + h, y2, f2);
			Vector<Real> rhs2 = f2 * h - k1 * 2.0;
			Vector<Real> k2 = rhs2;
			W_copy = W;
			try {
				GaussJordanSolver<Real>::SolveInPlace(W_copy, k2);
			} catch (...) {
				return false;
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
		ODESystemSolution Solve(Real t0, const Vector<Real>& y0, Real t_end, Real h_init) {
			int dim = y0.size();
			int max_steps = 50000; // Safety limit (storage grows dynamically if needed)

			ODESystemSolution sol(t0, t_end, dim, max_steps);

			Real t = t0;
			Vector<Real> y = y0;
			Real h = h_init;

			// Safety factors for step size control
			const Real safety = 0.9;
			const Real h_min = 1e-12;		   // Minimum step size
			const Real h_max = h_init * 100.0; // Maximum step size
			const Real fac_max = 5.0;		   // Maximum step increase factor
			const Real fac_min = 0.2;		   // Minimum step decrease factor

			sol.fillValues(0, t, y);
			int step = 1;

			while (t < t_end && step < max_steps) {
				Real h_actual = std::min(h, t_end - t);

				// Don't take tiny final steps
				if (t + h_actual * 1.01 >= t_end) {
					h_actual = t_end - t;
				}

				Vector<Real> y_new(dim);
				Vector<Real> y_err(dim);

				bool success = Step(t, y, h_actual, y_new, y_err);

				if (!success) {
					// Linear solve failed, reduce step size
					h = h_actual * 0.5;
					if (h < h_min) {
						throw std::runtime_error("Rosenbrock: Linear solve failed, system may be singular");
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
				err_norm = std::max(err_norm, 1e-10);

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
					h = std::min(h_max, h_actual * factor);
				} else {
					// Reject and retry with smaller step
					Real factor = safety * std::pow(1.0 / err_norm, 0.5); // More aggressive for rejection
					factor = std::max(fac_min, factor);
					h = h_actual * factor;

					if (h < h_min) {
						throw std::runtime_error("Rosenbrock: Step size too small (" + std::to_string(h) +
												 "), system may be too stiff or tolerance too tight");
					}
					sol.incrementRejectedSteps();
				}
			}

			if (step >= max_steps) {
				throw std::runtime_error("Rosenbrock: Maximum steps exceeded");
			}

			sol.setFinalSize(step - 1);
			return sol;
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
			int max_steps = config.max_steps;
			StiffSolverResult result(t0, t_end, dim, max_steps);
			result.algorithm_name = "Rosenbrock23";

			// Update internal tolerances
			_abs_tol = config.tolerance;
			_rel_tol = config.tolerance;

			try {
				Real h_init = (config.step_size > 0) ? config.step_size : (t_end - t0) / 100.0;
				result.solution = Solve(t0, y0, t_end, h_init);
				result.total_steps = result.solution.getNumStepsOK();
			} catch (const std::exception& e) {
				result.status = AlgorithmStatus::NumericalInstability;
				result.error_message = e.what();
			}

			result.elapsed_time_ms = timer.elapsed_ms();
			return result;
		}
	};

} // namespace MML

#endif // MML_ODE_STIFF_SOLVERS_H

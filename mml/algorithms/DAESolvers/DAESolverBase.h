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

#include "mml/MMLBase.h"

#include "mml/base/Vector/Vector.h"
#include "mml/base/Matrix/Matrix.h"
#include "mml/base/DAESystem.h"
#include "mml/interfaces/IODESystemDAE.h"
#include "mml/core/AlgorithmTypes.h"
#include "mml/core/LinAlgEqSolvers.h"

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

		/// Newton convergence tolerance (default: 1e-8)
		Real newton_tol = 1e-8;

		/// Tolerance for constraint satisfaction g(t,x,y) ≈ 0 (default: 1e-10)
		Real constraint_tol = 1e-10;

		/// Maximum number of steps (default: 100000)
		int max_steps = 100000;

		/// Factory method for high-precision configuration
		static DAESolverConfig HighPrecision() {
			DAESolverConfig config;
			config.step_size = 0.001;
			config.newton_tol = 1e-12;
			config.constraint_tol = 1e-12;
			config.max_newton_iter = 30;
			return config;
		}

		/// Factory method for fast (lower precision) configuration
		static DAESolverConfig Fast() {
			DAESolverConfig config;
			config.step_size = 0.05;
			config.newton_tol = 1e-6;
			config.constraint_tol = 1e-8;
			config.max_newton_iter = 10;
			return config;
		}
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
		Real max_constraint_violation = 0.0;  ///< Maximum |g(t,x,y)| seen

		DAESolverResult(Real t0, Real t_end, int diffDim, int algDim, int num_steps)
			: solution(t0, t_end, diffDim, algDim, num_steps) {}
	};

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
	inline bool ComputeConsistentIC(const IODESystemDAEWithJacobian& system,
	                                Real t0, const Vector<Real>& x0, Vector<Real>& y0,
	                                int max_iter = 20, Real tol = 1e-10)
	{
		int algDim = system.getAlgDim();
		Vector<Real> g(algDim);
		Matrix<Real> dg_dy(algDim, algDim);

		for (int iter = 0; iter < max_iter; ++iter)
		{
			// Evaluate constraint residual
			system.algConstraints(t0, x0, y0, g);

			// Check convergence
			Real g_norm = g.NormL2();
			if (g_norm < tol)
				return true;

			// Compute Jacobian ∂g/∂y
			system.jacobian_gy(t0, x0, y0, dg_dy);

			// Newton step: solve dg_dy * Δy = -g
			Vector<Real> rhs = g * (-1.0);
			Matrix<Real> A = dg_dy;  // Copy for in-place solve
			GaussJordanSolver<Real>::SolveInPlace(A, rhs);

			// Update y
			y0 = y0 + rhs;
		}

		return false;  // Failed to converge
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
	                               Real tol = 1e-10)
	{
		int algDim = system.getAlgDim();
		Vector<Real> g(algDim);
		system.algConstraints(t0, x0, y0, g);
		return g.NormL2() < tol;
	}

} // namespace MML

#endif // MML_DAE_SOLVER_BASE_H

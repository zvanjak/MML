///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DAERadauIIA.h                                                       ///
///  Description: Radau IIA method for semi-explicit index-1 DAE systems              ///
///                                                                                   ///
///  3-stage implicit Runge-Kutta method with Radau nodes.                            ///
///  Gold standard for stiff DAEs - used in SUNDIALS IDA.                             ///
///                                                                                   ///
///  Properties:                                                                      ///
///  - Order 5 accuracy                                                               ///
///  - L-stable (strong damping of stiff components)                                  ///
///  - Excellent for index-1 and some index-2 problems                                ///
///                                                                                   ///
///  Reference: Hairer & Wanner, "Solving ODEs II", Section IV.8                      ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DAE_RADAU_IIA_H
#define MML_DAE_RADAU_IIA_H

#include "DAESolverBase.h"
#include <cmath>

namespace MML {

	/// @brief Radau IIA - 3-stage implicit Runge-Kutta method for index-1 DAEs
	/// 
	/// Complexity: O((3N)³) = O(27N³) per step where N = diffDim + algDim.
	///            All 3 stages are coupled, requiring Newton on a 3N-dimensional system.
	///            Most expensive DAE solver per step, but highest order (5th).
	///
	/// Radau IIA is the gold standard for stiff DAEs:
	/// - L-stable (strong damping of stiff components)
	/// - Order 5 accuracy
	/// - Excellent for index-1 and some index-2 problems
	/// - Used in SUNDIALS IDA and other production solvers
	///
	/// The 3-stage Radau IIA has coefficients:
	///   c1 = (4 - sqrt(6))/10 ≈ 0.1550510257
	///   c2 = (4 + sqrt(6))/10 ≈ 0.6449489743
	///   c3 = 1
	///
	/// All stages are implicit and coupled, requiring Newton iteration
	/// on a 3*(diffDim+algDim) dimensional system.
	///
	/// Reference: Hairer & Wanner, "Solving ODEs II", Section IV.8
	///
	/// @param system DAE system with Jacobian
	/// @param t0 Initial time
	/// @param x0 Initial differential state
	/// @param y0 Initial algebraic state (must satisfy constraints)
	/// @param t_end Final time
	/// @param config Solver configuration
	/// @return DAESolverResult with solution and diagnostics
	inline DAESolverResult SolveDAERadauIIA(IODESystemDAEWithJacobian& system,
	                                         Real t0, const Vector<Real>& x0, const Vector<Real>& y0,
	                                         Real t_end, const DAESolverConfig& config = DAESolverConfig())
	{
		AlgorithmTimer timer;

		int diffDim = system.getDiffDim();
		int algDim = system.getAlgDim();
		int totalDim = diffDim + algDim;
		int num_steps = static_cast<int>((t_end - t0) / config.step_size) + 1;

		DAESolverResult result(t0, t_end, diffDim, algDim, num_steps);
		result.algorithm_name = "DAERadauIIA";
		result.solution.fillValues(0, t0, x0, y0);
		result.solution.incrementSuccessfulSteps();

		// Radau IIA coefficients (3-stage, order 5)
		const Real sqrt6 = std::sqrt(Real(6));
		
		// Abscissae (nodes)
		const Real c1 = (4 - sqrt6) / 10;  // ≈ 0.1550510257
		const Real c2 = (4 + sqrt6) / 10;  // ≈ 0.6449489743
		const Real c3 = 1;

		// Butcher tableau A matrix (for Radau IIA order 5)
		const Real a11 = (88 - 7*sqrt6) / 360;
		const Real a12 = (296 - 169*sqrt6) / 1800;
		const Real a13 = (-2 + 3*sqrt6) / 225;
		
		const Real a21 = (296 + 169*sqrt6) / 1800;
		const Real a22 = (88 + 7*sqrt6) / 360;
		const Real a23 = (-2 - 3*sqrt6) / 225;
		
		const Real a31 = (16 - sqrt6) / 36;
		const Real a32 = (16 + sqrt6) / 36;
		const Real a33 = 1.0 / 9;

		// Weights (b = last row of A for Radau methods)
		// b1 = a31, b2 = a32, b3 = a33

		// Current state
		Vector<Real> x = x0;
		Vector<Real> y = y0;
		Real t = t0;
		Real h = config.step_size;

		// Stage values: K_i for differential, L_i for algebraic
		// K_i = f(t + c_i*h, X_i, Y_i), L_i = y-values at stage
		Vector<Real> X1(diffDim), X2(diffDim), X3(diffDim);  // Stage x-values
		Vector<Real> Y1(algDim), Y2(algDim), Y3(algDim);     // Stage y-values
		Vector<Real> K1(diffDim), K2(diffDim), K3(diffDim);  // f evaluations
		Vector<Real> G1(algDim), G2(algDim), G3(algDim);     // g evaluations

		// Newton iteration vectors
		int newtDim = 3 * totalDim;  // 3 stages × (diffDim + algDim)
		Vector<Real> Z(newtDim), delta(newtDim), residual(newtDim);
		Matrix<Real> Jac_newton(newtDim, newtDim);

		// Jacobian matrices from the DAE system
		Matrix<Real> dfDx(diffDim, diffDim), dfDy(diffDim, algDim);
		Matrix<Real> dgDx(algDim, diffDim), dgDy(algDim, algDim);

		int step = 1;
		while (t + h/2 <= t_end && step < config.max_steps) {

			// Clamp step size to reach t_end exactly
			Real h_actual = std::min(h, t_end - t);
			if (h_actual < 1e-15) break;

			// Initialize stage values to current solution (predictor)
			X1 = x; X2 = x; X3 = x;
			Y1 = y; Y2 = y; Y3 = y;

			// Get Jacobians at current point (frozen for all Newton iterations in this step)
			system.allJacobians(t, x, y, dfDx, dfDy, dgDx, dgDy);

			// Newton iteration to solve the coupled stage equations
			bool converged = false;
			for (int newt = 0; newt < config.max_newton_iter && !converged; ++newt) {
				result.newton_iterations++;

				// Evaluate f and g at each stage point
				system.diffEqs(t + c1*h_actual, X1, Y1, K1);
				system.diffEqs(t + c2*h_actual, X2, Y2, K2);
				system.diffEqs(t + c3*h_actual, X3, Y3, K3);
				system.algConstraints(t + c1*h_actual, X1, Y1, G1);
				system.algConstraints(t + c2*h_actual, X2, Y2, G2);
				system.algConstraints(t + c3*h_actual, X3, Y3, G3);

				// Build residual for the stage equations
				// Stage equations for differential: X_i = x + h_actual * sum_j(a_ij * K_j)
				// Residual: R_xi = X_i - x - h_actual*(a_i1*K1 + a_i2*K2 + a_i3*K3)
				// Stage equations for algebraic: 0 = g(t + c_i*h_actual, X_i, Y_i)
				// Residual: R_yi = G_i

				Real residNorm = 0;
				
				// Stage 1 residuals
				for (int i = 0; i < diffDim; ++i) {
					Real r = X1[i] - x[i] - h_actual*(a11*K1[i] + a12*K2[i] + a13*K3[i]);
					residual[i] = r;
					residNorm += r*r;
				}
				for (int i = 0; i < algDim; ++i) {
					Real r = G1[i];
					residual[diffDim + i] = r;
					residNorm += r*r;
				}

				// Stage 2 residuals
				int off2 = totalDim;
				for (int i = 0; i < diffDim; ++i) {
					Real r = X2[i] - x[i] - h_actual*(a21*K1[i] + a22*K2[i] + a23*K3[i]);
					residual[off2 + i] = r;
					residNorm += r*r;
				}
				for (int i = 0; i < algDim; ++i) {
					Real r = G2[i];
					residual[off2 + diffDim + i] = r;
					residNorm += r*r;
				}

				// Stage 3 residuals
				int off3 = 2 * totalDim;
				for (int i = 0; i < diffDim; ++i) {
					Real r = X3[i] - x[i] - h_actual*(a31*K1[i] + a32*K2[i] + a33*K3[i]);
					residual[off3 + i] = r;
					residNorm += r*r;
				}
				for (int i = 0; i < algDim; ++i) {
					Real r = G3[i];
					residual[off3 + diffDim + i] = r;
					residNorm += r*r;
				}

				residNorm = std::sqrt(residNorm);
				if (residNorm < config.newton_tol) {
					converged = true;
					break;
				}

				// Build Newton Jacobian
				// The unknowns are [X1, Y1, X2, Y2, X3, Y3] (3*totalDim vector)
				// Jacobian structure (block form):
				// 
				// For stage i, diff eq: dR_xi/dXj = delta_ij*I - h_actual*a_ij*df/dx
				//                       dR_xi/dYj = -h_actual*a_ij*df/dy
				// For stage i, alg eq:  dR_yi/dXj = dg/dx (only for j=i)
				//                       dR_yi/dYj = dg/dy (only for j=i)

				// Zero out the Jacobian matrix
				for (int i = 0; i < newtDim; ++i)
					for (int j = 0; j < newtDim; ++j)
						Jac_newton(i, j) = 0;

				// Fill the 3x3 block structure
				// Block (i,j) corresponds to derivatives of stage i residual w.r.t. stage j variables
				for (int si = 0; si < 3; ++si) {
					int row_off = si * totalDim;
					
					for (int sj = 0; sj < 3; ++sj) {
						int col_off = sj * totalDim;
						
						// Get a_ij coefficient
						Real aij = 0;
						if (si == 0) aij = (sj == 0) ? a11 : (sj == 1) ? a12 : a13;
						else if (si == 1) aij = (sj == 0) ? a21 : (sj == 1) ? a22 : a23;
						else aij = (sj == 0) ? a31 : (sj == 1) ? a32 : a33;

						// Differential equations block: dR_xi/d(Xj, Yj)
						for (int i = 0; i < diffDim; ++i) {
							for (int j = 0; j < diffDim; ++j) {
								Real val = (si == sj && i == j) ? 1.0 : 0.0;
								val -= h_actual * aij * dfDx(i, j);
								Jac_newton(row_off + i, col_off + j) = val;
							}
							for (int j = 0; j < algDim; ++j) {
								Jac_newton(row_off + i, col_off + diffDim + j) = -h_actual * aij * dfDy(i, j);
							}
						}

						// Algebraic equations block: dR_yi/d(Xj, Yj)
						// Only non-zero when si == sj (same stage)
						if (si == sj) {
							for (int i = 0; i < algDim; ++i) {
								for (int j = 0; j < diffDim; ++j) {
									Jac_newton(row_off + diffDim + i, col_off + j) = dgDx(i, j);
								}
								for (int j = 0; j < algDim; ++j) {
									Jac_newton(row_off + diffDim + i, col_off + diffDim + j) = dgDy(i, j);
								}
							}
						}
					}
				}

				// Solve Newton system: Jac_newton * delta = -residual
				delta = residual;
				delta *= -1;
				GaussJordanSolver<Real>::SolveInPlace(Jac_newton, delta);

				// Update stage values
				for (int i = 0; i < diffDim; ++i) {
					X1[i] += delta[i];
					X2[i] += delta[totalDim + i];
					X3[i] += delta[2*totalDim + i];
				}
				for (int i = 0; i < algDim; ++i) {
					Y1[i] += delta[diffDim + i];
					Y2[i] += delta[totalDim + diffDim + i];
					Y3[i] += delta[2*totalDim + diffDim + i];
				}
			}

			if (!converged) {
				result.status = AlgorithmStatus::NumericalInstability;
				result.error_message = "Newton iteration failed to converge in Radau IIA step";
				break;
			}

			// Update solution using stage 3 (for Radau IIA, c3 = 1, so X3, Y3 is the solution at t+h_actual)
			x = X3;
			y = Y3;
			t += h_actual;

			// Track constraint violation
			Vector<Real> g_check(algDim);
			system.algConstraints(t, x, y, g_check);
			Real g_norm = g_check.NormL2();
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

#endif // MML_DAE_RADAU_IIA_H

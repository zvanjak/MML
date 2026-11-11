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
	/** @f[
	 *   \begin{bmatrix}
	 *     I - h\frac{\partial f}{\partial x} & -h\frac{\partial f}{\partial y} \\
	 *     \frac{\partial g}{\partial x} & \frac{\partial g}{\partial y}
	 *   \end{bmatrix}
	 *   \begin{bmatrix} \Delta x \\ \Delta y \end{bmatrix}
	 *   =
	 *   \begin{bmatrix} -F \\ -G \end{bmatrix}
	 * @f]
	 */
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
		return DAEImplicitDetail::Solve(system, t0, x0, y0, t_end, config, false);
	}

} // namespace MML

#endif // MML_DAE_BACKWARD_EULER_H

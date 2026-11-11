///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DAEBDF2.h                                                           ///
///  Description: BDF2 method for semi-explicit index-1 DAE systems                   ///
///                                                                                   ///
///  BDF2 formula for differential equations:                                         ///
///    x_{n+1} = (4/3)x_n - (1/3)x_{n-1} + (2/3)h·f(t_{n+1}, x_{n+1}, y_{n+1})        ///
///    0 = g(t_{n+1}, x_{n+1}, y_{n+1})                                               ///
///                                                                                   ///
///  Properties:                                                                      ///
///  - 2nd order accurate                                                             ///
///  - A-stable                                                                       ///
///  - Industry standard for moderately stiff DAEs                                    ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DAE_BDF2_H
#define MML_DAE_BDF2_H

#include "DAESolverBase.h"

namespace MML {

	/// @brief Second-order Backward Differentiation Formula solver for index-1 DAEs
	/// 
	/// Complexity: O(N³) per step where N = diffDim + algDim.
	///            Newton iteration with augmented Jacobian + LU solve each step.
	///
	/// BDF2 formula for differential equations:
	/// @f[
	///   x_{n+1} = \frac{4}{3}x_n - \frac{1}{3}x_{n-1} + \frac{2}{3}h \cdot f(t_{n+1}, x_{n+1}, y_{n+1})
	/// @f]
	/// @f[
	///   0 = g(t_{n+1}, x_{n+1}, y_{n+1})
	/// @f]
	///
	/// Uses Backward Euler for the first step to bootstrap, then switches to BDF2.
	/// Second-order accurate, A-stable, industry standard for moderately stiff DAEs.
	///
	/// @param system DAE system with Jacobian
	/// @param t0 Initial time
	/// @param x0 Initial differential state
	/// @param y0 Initial algebraic state (must satisfy constraints)
	/// @param t_end Final time
	/// @param config Solver configuration
	/// @return DAESolverResult with solution and diagnostics
	inline DAESolverResult SolveDAEBDF2(IODESystemDAEWithJacobian& system,
	                                     Real t0, const Vector<Real>& x0, const Vector<Real>& y0,
	                                     Real t_end, const DAESolverConfig& config = DAESolverConfig())
	{
		return DAEImplicitDetail::Solve(system, t0, x0, y0, t_end, config, true);
	}

} // namespace MML

#endif // MML_DAE_BDF2_H

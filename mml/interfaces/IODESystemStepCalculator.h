///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        IODESystemStepCalculator.h                                          ///
///  Description: Interface for ODE step calculation strategies                       ///
///               Defines single-step computation protocol                            ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @file IODESystemStepCalculator.h
 * @brief Interface for ODE single-step calculation algorithms.
 * 
 * Defines the contract for step calculators that advance an ODE solution
 * by one step of size h. This is the core computational unit used by
 * ODE solvers to implement various numerical methods.
 * 
 * Implementations include Euler, Runge-Kutta (RK4, RK45), and other
 * single-step methods. The interface separates step computation from
 * step-size control, allowing flexible composition.
 * 
 * @see IODESystemStepper for adaptive step-size control
 * @see IODESystem for the ODE system interface
 */

#if !defined MML_IODE_SYSTEM_STEP_CALC_H
#define MML_IODE_SYSTEM_STEP_CALC_H

#include "MMLBase.h"

#include "interfaces/IODESystem.h"

#include "base/Vector/Vector.h"

namespace MML
{
	/**
	 * @brief Interface for single-step ODE integration methods.
	 * 
	 * Computes one integration step from state x at time t to an
	 * approximation at time t+h, along with an error estimate.
	 * 
	 * The error estimate xerr is used by adaptive step-size controllers
	 * to adjust h for optimal accuracy/efficiency trade-off.
	 * 
	 * **Implementation Requirements:**
	 * - Must be stateless (all state passed as parameters)
	 * - Should compute error estimate for adaptive stepping
	 * - Thread-safe for parallel integration of independent systems
	 */
	class IODESystemStepCalculator
	{
	public:
		/**
		 * @brief Compute a single integration step.
		 * 
		 * Advances the ODE solution from (t, x_start) to approximately (t+h, xout)
		 * using the provided derivative dxdt = f(t, x_start).
		 * 
		 * @param odeSystem The ODE system being integrated
		 * @param t Current time
		 * @param x_start Current state vector
		 * @param dxdt Derivative at (t, x_start), pre-computed by caller
		 * @param h Step size
		 * @param[out] xout State vector at t+h (approximate solution)
		 * @param[out] xerr Error estimate for adaptive step control
		 * 
		 * @note The caller provides dxdt to avoid redundant derivative evaluations
		 *       when the same derivative is needed for multiple purposes.
		 */
		virtual void calcStep(const IODESystem& odeSystem, 
													const Real t, const Vector<Real>& x_start, const Vector<Real>& dxdt,
													const Real h, Vector<Real>& xout, Vector<Real>& xerr) const = 0;
	};
}
#endif

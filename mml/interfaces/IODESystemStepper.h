///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        IODESystemStepper.h                                                 ///
///  Description: Interface for ODE stepping with error control                       ///
///               Adaptive step size management protocol                              ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @file IODESystemStepper.h
 * @brief Interface for ODE stepping with adaptive error control.
 * 
 * Defines the contract for ODE steppers that manage step-size adaptation
 * based on local error estimates. Unlike IODESystemStepCalculator which
 * performs a single computation, steppers maintain state for adaptive
 * control and may reject/retry steps.
 * 
 * The stepper combines:
 * - Step calculation (via an embedded step calculator or direct implementation)
 * - Error estimation and comparison to tolerances
 * - Step-size adjustment heuristics
 * 
 * @see IODESystemStepCalculator for stateless step computation
 * @see ODESolver for the high-level integration driver
 */

#if !defined MML_IODE_SYSTEM_STEPPER_H
#define MML_IODE_SYSTEM_STEPPER_H

#include "MMLBase.h"

#include "interfaces/IODESystem.h"

#include "base/Vector/Vector.h"
#include "base/Matrix/Matrix.h"

namespace MML
{
	/**
	 * @brief Interface for adaptive ODE steppers.
	 * 
	 * Performs integration steps with error control, potentially maintaining
	 * internal state for step-size adaptation, FSAL (First Same As Last)
	 * optimizations, or dense output.
	 * 
	 * Unlike the stateless IODESystemStepCalculator, a stepper may:
	 * - Cache intermediate results between steps
	 * - Track step-size history for optimal adaptation
	 * - Provide dense output between step endpoints
	 * 
	 * @note Non-const doStep() allows implementations to update internal state.
	 */
	class IODESystemStepper	{
	public:
		/**
		 * @brief Perform an adaptive integration step.
		 * 
		 * Advances from (t, x_start) to (t+h, xout) with error estimation.
		 * Implementations may adjust internal parameters based on the
		 * computed error for subsequent step-size recommendations.
		 * 
		 * @param odeSystem The ODE system being integrated
		 * @param t Current time
		 * @param x_start Current state vector
		 * @param dxdt Derivative at (t, x_start)
		 * @param h Step size (may be suggested size; actual may differ)
		 * @param[out] xout State vector at t+h
		 * @param[out] xerr Error estimate for this step
		 */
		virtual void doStep(const IODESystem& odeSystem, 
												const Real t, const Vector<Real>& x_start, const Vector<Real>& dxdt,
												const Real h, Vector<Real>& xout, Vector<Real>& xerr) = 0;
	};
}
#endif

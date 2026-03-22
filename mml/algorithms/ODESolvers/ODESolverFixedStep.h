///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        ODESolverFixedStep.h                                           ///
///  Description: Fixed-step ODE system integrators                                   ///
///                                                                                   ///
///  Contents:    ODESystemFixedStepSolver - Fixed-step integration with pluggable    ///
///                                          step calculators (Euler, RK4, etc.)      ///
///                                                                                   ///
///  For adaptive step-size integration, see ODESolverAdaptive.h                  ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_ODE_SOLVER_FIXED_STEP_H
#define MML_ODE_SOLVER_FIXED_STEP_H

#include "mml/MMLBase.h"

#include "mml/interfaces/IODESystem.h"
#include "mml/interfaces/IODESystemStepCalculator.h"

#include "mml/base/ODESystem.h"
#include "mml/base/ODESystemSolution.h"

namespace MML {
	/// @brief Fixed-step ODE system solver
	///
	/// This class provides a simple fixed-step integration method using
	/// pluggable step calculators (Euler, RK4, etc.)
	///
	/// @note For adaptive step-size control, use ODEAdaptiveIntegrator instead.
	///
	/// @warning LIFETIME REQUIREMENT: The IODESystem and IODESystemStepCalculator references
	///          must remain valid for the entire lifetime of this solver object.
	class ODESystemFixedStepSolver {
		const IODESystem& _odeSys;
		const IODESystemStepCalculator& _stepCalc;

	public:
		/// @brief Construct solver with system and step calculator
		/// @warning The references must outlive this solver instance
		ODESystemFixedStepSolver(const IODESystem& inOdeSys, const IODESystemStepCalculator& inStepCalc)
			: _odeSys(inOdeSys)
			, _stepCalc(inStepCalc) {}

		/// @brief Integrate the ODE system with fixed step size
		/// @param initCond Initial conditions vector
		/// @param t1 Start time
		/// @param t2 End time
		/// @param numSteps Number of steps to take
		/// @return ODESystemSolution with numSteps+1 values
		ODESystemSolution integrate(const Vector<Real>& initCond, Real t1, Real t2, int numSteps) {
			int dim = _odeSys.getDim();

			ODESystemSolution sol(t1, t2, dim, numSteps);
			Vector<Real> x(initCond), x_out(dim), dxdt(dim), x_err(dim);

			sol.fillValues(0, t1, x); // store initial values

			Real t = t1;
			Real h = (t2 - t1) / numSteps;

			for (int k = 1; k <= numSteps; k++) {
				_odeSys.derivs(t, x, dxdt);
				_stepCalc.calcStep(_odeSys, t, x, dxdt, h, x_out, x_err);
				t += h;
				x = x_out;
				sol.fillValues(k, t, x);
			}

			return sol;
		}
	};

} // namespace MML

#endif // MML_ODE_SOLVER_FIXED_STEP_H

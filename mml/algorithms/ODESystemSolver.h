///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        ODESystemSolver.h                                                   ///
///  Description: ODE system integrator with adaptive step control                    ///
///                                                                                   ///
///  NOTE:        This file provides backward compatibility.                          ///
///               The implementation has been moved to ODEAdaptiveIntegrator.h        ///
///               Prefer using the new ODEAdaptiveIntegrator directly.                ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_ODE_SYSTEM_SOLVERS_H
#define MML_ODE_SYSTEM_SOLVERS_H

#include "MMLBase.h"

#include "interfaces/IODESystem.h"
#include "interfaces/IODESystemStepCalculator.h"

#include "base/ODESystem.h"
#include "base/ODESystemSolution.h"

#include "algorithms/ODEAdaptiveIntegrator.h"
#include "algorithms/ODESystemSteppers.h"

namespace MML
{
	///////////////////////////////////////////////////////////////////////////////////////////
	/// @brief      LEGACY: Adaptive-step ODE system solver with pluggable stepper
	///
	/// @deprecated Prefer using CashKarpIntegrator, DormandPrince5Integrator, or
	///             DormandPrince8Integrator directly from ODEAdaptiveIntegrator.h
	///
	/// This class provides backward compatibility for code using the old interface.
	/// It delegates to the new ODEAdaptiveIntegrator infrastructure.
	///
	/// @tparam Stepper The stepper type (e.g., RK5_CashKarp_Stepper)
	///
	/// @example Legacy usage (still works):
	/// @code
	///   ODESystemSolver<RK5_CashKarp_Stepper> solver(system);
	///   auto solution = solver.integrate(initCond, t1, t2, saveInterval, eps, h1, hmin);
	/// @endcode
	///
	/// @example New preferred usage:
	/// @code
	///   CashKarpIntegrator integrator(system);
	///   auto solution = integrator.integrate(initCond, t1, t2, saveInterval, eps);
	/// @endcode
	///////////////////////////////////////////////////////////////////////////////////////////
	template<class Stepper>
	class ODESystemSolver
	{
		const IODESystem& _sys;

	public:
		/// @brief Construct solver for given ODE system
		/// @param sys The ODE system to solve
		ODESystemSolver(const IODESystem& sys) : _sys(sys) {}

		int getDim() { return _sys.getDim(); }

		/// @brief Integrate ODE system with adaptive step control (legacy interface)
		/// @param initCond Initial conditions vector
		/// @param t1 Start time
		/// @param t2 End time
		/// @param minSaveInterval Minimum interval between saved points
		/// @param eps Error tolerance
		/// @param h1 Initial step size guess (ignored in new implementation - auto-estimated)
		/// @param hmin Minimum step size (ignored in new implementation)
		/// @return ODESystemSolution with integration results
		ODESystemSolution integrate(const Vector<Real>& initCond,
		                            Real t1, Real t2, Real minSaveInterval,
		                            Real eps, Real h1, Real hmin = 0)
		{
			// Use the new adaptive integrator with Cash-Karp stepper
			// h1 and hmin are ignored - new implementation auto-estimates initial step
			// and has better step-size control
			ODEAdaptiveIntegrator<CashKarp_Stepper> integrator(_sys);
			return integrator.integrate(initCond, t1, t2, minSaveInterval, eps, h1);
		}
	};

	// Convenience type alias for most common usage
	using AdaptiveODESolver = ODESystemSolver<RK5_CashKarp_Stepper>;
	/// @brief Fixed-step ODE system solver
	/// 
	/// This class provides a simple fixed-step integration method using
	/// pluggable step calculators (Euler, RK4, etc.)
	/// 
	/// @note For adaptive step-size control, use ODEAdaptiveIntegrator instead.
	/// 
	/// @warning LIFETIME REQUIREMENT: The IODESystem and IODESystemStepCalculator references
	///          must remain valid for the entire lifetime of this solver object.
	class ODESystemFixedStepSolver
	{
		const IODESystem& _odeSys;
		const IODESystemStepCalculator& _stepCalc;

	public:
		/// @brief Construct solver with system and step calculator
		/// @warning The references must outlive this solver instance
		ODESystemFixedStepSolver(const IODESystem& inOdeSys, const IODESystemStepCalculator& inStepCalc)
			: _odeSys(inOdeSys), _stepCalc(inStepCalc) {
		}

		/// @brief Integrate the ODE system with fixed step size
		/// @param initCond Initial conditions vector
		/// @param t1 Start time
		/// @param t2 End time
		/// @param numSteps Number of steps to take
		/// @return ODESystemSolution with numSteps+1 values
		ODESystemSolution integrate(const Vector<Real>& initCond, Real t1, Real t2, int numSteps)
		{
			int dim = _odeSys.getDim();

			ODESystemSolution sol(t1, t2, dim, numSteps);
			Vector<Real> x(initCond), x_out(dim), dxdt(dim), x_err(dim);

			sol.fillValues(0, t1, x);		// store initial values

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

	/// @brief Leapfrog (Velocity Verlet) integrator for Hamiltonian systems
	/// 
	/// Symplectic method that conserves energy for long-time integration.
	/// Expects state vector split as [positions..., velocities...].
	class ODESystemLeapfrogSolver
	{
		IODESystem& _odeSys;

	public:
		ODESystemLeapfrogSolver(IODESystem& inOdeSys) : _odeSys(inOdeSys) {}

		ODESystemSolution integrate(const Vector<Real>& initCond, Real t1, Real t2, int numSteps)
		{
			int dim = _odeSys.getDim();
			int half_n = dim / 2;
			ODESystemSolution sol(t1, t2, dim, numSteps);

			Vector<Real> x(initCond);
			Vector<Real> dxdt(dim);
			Vector<Real> x_out(dim);

			sol.fillValues(0, t1, x);

			Real t = t1;
			Real h = (t2 - t1) / numSteps;

			for (int k = 1; k <= numSteps; k++) {
				_odeSys.derivs(t, x, dxdt);

				// Half-step velocity
				Vector<Real> v_half(half_n);
				for (int i = 0; i < half_n; ++i)
					v_half[i] = x[half_n + i] + 0.5 * h * dxdt[half_n + i];

				// Full-step position
				for (int i = 0; i < half_n; ++i)
					x_out[i] = x[i] + h * v_half[i];

				// Compute new acceleration
				Vector<Real> x_next = x_out;
				for (int i = 0; i < half_n; ++i)
					x_next[half_n + i] = v_half[i];

				Vector<Real> dxdt_next(dim);
				_odeSys.derivs(t + h, x_next, dxdt_next);

				// Full-step velocity
				for (int i = 0; i < half_n; ++i)
					x_out[half_n + i] = v_half[i] + 0.5 * h * dxdt_next[half_n + i];

				t += h;
				x = x_out;
				sol.fillValues(k, t, x);
			}

			return sol;
		}
	};

} // namespace MML

#endif // MML_ODE_SYSTEM_SOLVERS_H

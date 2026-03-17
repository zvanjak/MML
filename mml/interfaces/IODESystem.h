///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        IODESystem.h                                                        ///
///  Description: Interface for ordinary differential equation systems                ///
///               Defines evaluation protocol for ODE solvers                         ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @file IODESystem.h
 * @brief Interfaces for ordinary differential equation (ODE) systems.
 * 
 * Defines the contract for ODE systems that can be solved by MML's ODE solvers.
 * The hierarchy supports:
 * 
 * - **IODESystem** - Basic first-order ODE system: dx/dt = f(t, x)
 * - **IODESystemWithJacobian** - System with analytic Jacobian for implicit methods
 * - **IODESystemParametrized** - System with adjustable parameters
 * 
 * All ODE systems are expressed in first-order form. Higher-order ODEs must be
 * converted to equivalent first-order systems (e.g., y'' = f becomes y1' = y2, y2' = f).
 * 
 * @see ODESolver, ODESystemSolution, RungeKutta methods
 */

#if !defined MML_IODESYSTEM_H
#define MML_IODESYSTEM_H

#include "MMLBase.h"

#include "base/Vector/Vector.h"
#include "base/Matrix/Matrix.h"
#include "interfaces/IParametrized.h"

namespace MML
{
	/**
	 * @brief Interface for first-order ODE systems.
	 * 
	 * Represents a system of first-order ordinary differential equations:
	 * @f[
	 *   \frac{d\mathbf{x}}{dt} = \mathbf{f}(t, \mathbf{x})
	 * @f]
	 * 
	 * where x is an N-dimensional state vector and f is the derivative function.
	 * 
	 * **Implementation Requirements:**
	 * - derivs() must be thread-safe for parallel integration
	 * - getDim() returns the system dimension (number of equations)
	 * - getVarName() provides human-readable variable names for output
	 * 
	 * @note This is the base interface used by all MML ODE solvers including
	 *       Euler, Runge-Kutta, and adaptive step-size methods.
	 */
	class IODESystem
	{
	public:
		virtual ~IODESystem() = default;
		
		/**
		 * @brief Get the dimension (number of equations) of the system.
		 * @return Number of first-order equations in the system
		 */
		virtual int   getDim() const = 0;
		
		/**
		 * @brief Compute the derivatives at a given state.
		 * 
		 * Evaluates the right-hand side f(t, x) of the ODE system.
		 * 
		 * @param t Current time/independent variable
		 * @param x Current state vector (size = getDim())
		 * @param[out] dxdt Output derivative vector dx/dt (size = getDim())
		 */
		virtual void  derivs(const Real t, const Vector<Real> &x, Vector<Real> &dxdt) const = 0;

		/**
		 * @brief Get a human-readable name for a state variable.
		 * 
		 * Override this method to provide meaningful names for output
		 * and debugging (e.g., "position", "velocity", "theta").
		 * 
		 * @param ind Index of the variable (0 to getDim()-1)
		 * @return Name string for the variable
		 */
		virtual std::string getVarName(int ind) const
		{
			if (ind < 0 || ind >= getDim())
				return "var" + std::to_string(ind);

			return "var" + std::to_string(ind + 1); // 1-based index 
		}
	};

	/**
	 * @brief ODE system with analytic Jacobian matrix.
	 * 
	 * Extends IODESystem to provide the Jacobian matrix ∂f/∂x, which is
	 * required by implicit integration methods (e.g., backward Euler,
	 * implicit Runge-Kutta) and stiff ODE solvers.
	 * 
	 * The Jacobian matrix J has elements:
	 * @f[
	 *   J_{ij} = \frac{\partial f_i}{\partial x_j}
	 * @f]
	 * 
	 * @note Providing an analytic Jacobian improves both accuracy and
	 *       performance compared to numerical differentiation.
	 */
	class IODESystemWithJacobian : public IODESystem
	{
	public:
		/**
		 * @brief Compute derivatives and Jacobian matrix simultaneously.
		 * 
		 * @param t Current time
		 * @param x Current state vector
		 * @param[out] dxdt Output derivative vector
		 * @param[out] dydx Output Jacobian matrix (getDim() × getDim())
		 */
		virtual void jacobian(const Real t, const Vector<Real>& x, Vector<Real>& dxdt, Matrix<Real>& dydx) const = 0;
	};

	/**
	 * @brief ODE system with adjustable parameters.
	 * 
	 * Extends IODESystem to support systems where coefficients or parameters
	 * can be modified at runtime. Useful for parameter sweeps, sensitivity
	 * analysis, and fitting ODE models to data.
	 * 
	 * Example: A damped oscillator with adjustable damping and stiffness:
	 * @f[
	 *   \ddot{x} + c\dot{x} + kx = 0
	 * @f]
	 * where c and k are parameters.
	 */
	class IODESystemParametrized : public IODESystem, public IParametrized
	{
	public:
		virtual ~IODESystemParametrized() = default;
	};

}
#endif
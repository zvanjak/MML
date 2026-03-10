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

	//===================================================================================
	//                           Event Detection Support
	//===================================================================================

	/**
	 * @brief Direction of zero-crossing to detect.
	 * 
	 * Controls which direction of event function zero-crossings trigger events:
	 * - Increasing: g(t,x) crosses zero from negative to positive
	 * - Decreasing: g(t,x) crosses zero from positive to negative  
	 * - Both: Any zero-crossing direction triggers event
	 */
	enum class EventDirection {
		Increasing,  ///< Trigger when g crosses from - to +
		Decreasing,  ///< Trigger when g crosses from + to -
		Both         ///< Trigger on any crossing
	};

	/**
	 * @brief Action to take when an event is detected.
	 * 
	 * - Continue: Record event and continue integration
	 * - Stop: Record event and terminate integration
	 * - Restart: Record event, call state update handler, restart with new state
	 */
	enum class EventAction {
		Continue,  ///< Continue integration after event
		Stop,      ///< Stop integration at event
		Restart    ///< Update state and restart integration
	};

	/**
	 * @brief Detailed information about a detected event.
	 * 
	 * Returned by event detection routines to describe when and where
	 * an event occurred during ODE integration.
	 */
	struct EventInfo {
		int eventIndex;         ///< Which event function triggered (0-based)
		Real time;              ///< Exact time of the event
		Vector<Real> state;     ///< State vector at event time
		EventDirection direction; ///< Actual crossing direction detected
		Real eventValue;        ///< Value of event function (should be ~0)
		
		EventInfo() : eventIndex(-1), time(0), direction(EventDirection::Both), eventValue(0) {}
		
		EventInfo(int idx, Real t, const Vector<Real>& x, EventDirection dir, Real val)
			: eventIndex(idx), time(t), state(x), direction(dir), eventValue(val) {}
	};

	/**
	 * @brief ODE system with event detection capability.
	 * 
	 * Extends IODESystem to support event detection during integration.
	 * Events are defined by event functions g_i(t, x) that trigger when they
	 * cross zero. This enables simulation of hybrid systems with discontinuities:
	 * 
	 * - Bouncing ball (floor contact)
	 * - Switching control systems
	 * - State-dependent termination conditions
	 * - Periodic event logging
	 * 
	 * **Implementation Requirements:**
	 * 1. Override getNumEvents() to return number of event functions
	 * 2. Override eventFunction() to evaluate each event function g_i(t, x)
	 * 3. Override getEventDirection() for each event's trigger direction
	 * 4. Override getEventAction() to specify behavior when event occurs
	 * 5. Optionally override handleEvent() for state updates (e.g., velocity reversal)
	 * 
	 * **Example: Bouncing Ball**
	 * @code
	 * class BouncingBall : public IODESystemWithEvents {
	 *     int getNumEvents() const override { return 1; }
	 *     
	 *     Real eventFunction(int i, Real t, const Vector<Real>& x) const override {
	 *         return x[0];  // Position = 0 is ground contact
	 *     }
	 *     
	 *     EventDirection getEventDirection(int i) const override {
	 *         return EventDirection::Decreasing;  // Falling toward ground
	 *     }
	 *     
	 *     EventAction getEventAction(int i) const override {
	 *         return EventAction::Restart;  // Need to reverse velocity
	 *     }
	 *     
	 *     void handleEvent(int i, Real t, Vector<Real>& x) const override {
	 *         x[1] = -0.9 * x[1];  // Reverse velocity with 10% loss
	 *     }
	 * };
	 * @endcode
	 */
	class IODESystemWithEvents : public IODESystem
	{
	public:
		/**
		 * @brief Get the number of event functions.
		 * @return Number of event functions g_i to monitor
		 */
		virtual int getNumEvents() const = 0;

		/**
		 * @brief Evaluate an event function.
		 * 
		 * The event triggers when g_i(t, x) crosses zero in the specified direction.
		 * Event functions should be smooth and continuous for accurate detection.
		 * 
		 * @param eventIndex Which event function to evaluate (0 to getNumEvents()-1)
		 * @param t Current time
		 * @param x Current state vector
		 * @return Value of event function g_i(t, x)
		 */
		virtual Real eventFunction(int eventIndex, Real t, const Vector<Real>& x) const = 0;

		/**
		 * @brief Get the direction for event detection.
		 * @param eventIndex Event index
		 * @return Direction that triggers this event
		 */
		virtual EventDirection getEventDirection(int eventIndex) const {
			return EventDirection::Both;  // Default: detect both directions
		}

		/**
		 * @brief Get the action to take when event is detected.
		 * @param eventIndex Event index  
		 * @return Action to perform (Continue, Stop, or Restart)
		 */
		virtual EventAction getEventAction(int eventIndex) const {
			return EventAction::Continue;  // Default: continue integration
		}

		/**
		 * @brief Handle an event by modifying the state.
		 * 
		 * Called when EventAction::Restart is specified. Override to implement
		 * state discontinuities such as velocity reversals or mode switches.
		 * 
		 * @param eventIndex Which event triggered
		 * @param t Time of event
		 * @param[in,out] x State vector to modify
		 */
		virtual void handleEvent(int eventIndex, Real t, Vector<Real>& x) const {
			// Default: no state modification
		}

		/**
		 * @brief Evaluate all event functions at once.
		 * 
		 * Default implementation calls eventFunction() for each index.
		 * Override for efficiency if evaluating multiple events shares computation.
		 * 
		 * @param t Current time
		 * @param x Current state vector
		 * @param[out] g Vector to fill with event function values (size = getNumEvents())
		 */
		virtual void eventFunctions(Real t, const Vector<Real>& x, Vector<Real>& g) const {
			for (int i = 0; i < getNumEvents(); ++i) {
				g[i] = eventFunction(i, t, x);
			}
		}
	};

}
#endif
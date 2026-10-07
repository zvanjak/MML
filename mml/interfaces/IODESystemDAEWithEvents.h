///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        IODESystemDAEWithEvents.h                                           ///
///  Description: Interface for DAE systems with zero-crossing event detection       ///
///               Supports direction filtering, terminal events, and state updates   ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @file IODESystemDAEWithEvents.h
 * @brief Interface for differential-algebraic systems with event detection.
 *
 * Extends the analytic-Jacobian DAE interface with scalar event functions
 * `g_i(t, x, y)`. An event is detected when one of these functions crosses zero
 * in its configured direction. Each event may continue integration, stop at the
 * root, or modify the differential and algebraic state before restarting.
 *
 * @see IODESystemDAEWithJacobian
 * @see IODESystemWithEvents
 * @see SolveDAEBackwardEulerWithEvents
 */

#if !defined MML_IODESYSTEM_DAE_WITH_EVENTS_H
#define MML_IODESYSTEM_DAE_WITH_EVENTS_H

#include <mml/interfaces/IODESystemDAE.h>
#include <mml/interfaces/IODESystemWithEvents.h>

namespace MML {

	/**
	 * @brief Analytic-Jacobian DAE system with zero-crossing events.
	 *
	 * Event functions may depend on time, the differential state `x`, and the
	 * algebraic state `y`. Implementations provide the number of event functions
	 * and evaluate each function independently. Direction, action, and state-update
	 * hooks have defaults and may be overridden per event.
	 *
	 * The DAE equations and all four Jacobian blocks remain defined by
	 * IODESystemDAEWithJacobian. Initial states and states produced by handleEvent()
	 * must be suitable for projection back onto the algebraic constraint manifold.
	 *
	 * @note Event functions should be continuous near their roots so the solver can
	 *       locate zero crossings reliably.
	 */
	class IODESystemDAEWithEvents : public IODESystemDAEWithJacobian {
	public:
		/**
		 * @brief Get the number of event functions monitored by the solver.
		 * @return Number of scalar event functions `g_i(t, x, y)`
		 */
		virtual int getNumEvents() const = 0;

		/**
		 * @brief Evaluate one event function.
		 *
		 * An event is eligible to trigger when this value crosses zero in the
		 * direction returned by getEventDirection().
		 *
		 * @param eventIndex Event function index in `[0, getNumEvents())`
		 * @param t Current time
		 * @param x Current differential state (size = getDiffDim())
		 * @param y Current algebraic state (size = getAlgDim())
		 * @return Value of the selected event function
		 */
		virtual Real eventFunction(int eventIndex, Real t, const Vector<Real>& x,
			const Vector<Real>& y) const = 0;

		/**
		 * @brief Get the zero-crossing direction that triggers an event.
		 * @param eventIndex Event function index
		 * @return Crossing direction to detect; defaults to EventDirection::Both
		 */
		virtual EventDirection getEventDirection(int eventIndex) const {
			return EventDirection::Both;
		}

		/**
		 * @brief Get the action performed after an event is located.
		 * @param eventIndex Event function index
		 * @return Event action; defaults to EventAction::Continue
		 */
		virtual EventAction getEventAction(int eventIndex) const {
			return EventAction::Continue;
		}

		/**
		 * @brief Modify the DAE state for a restart event.
		 *
		 * Called for EventAction::Restart after the event root is located. The
		 * solver subsequently recomputes a consistent algebraic state, so `y` may
		 * be used as an initial guess for that projection.
		 *
		 * @param eventIndex Event function index
		 * @param t Located event time
		 * @param[in,out] x Differential state to update
		 * @param[in,out] y Algebraic state or consistency-projection initial guess
		 */
		virtual void handleEvent(int eventIndex, Real t, Vector<Real>& x,
			Vector<Real>& y) const {}

		/**
		 * @brief Evaluate all event functions at once.
		 *
		 * The default implementation calls eventFunction() for every event index.
		 * Override when several event functions can share intermediate work.
		 *
		 * @param t Current time
		 * @param x Current differential state (size = getDiffDim())
		 * @param y Current algebraic state (size = getAlgDim())
		 * @param[out] values Pre-sized output vector (size = getNumEvents())
		 */
		virtual void eventFunctions(Real t, const Vector<Real>& x, const Vector<Real>& y,
			Vector<Real>& values) const {
			for (int i = 0; i < getNumEvents(); ++i)
				values[i] = eventFunction(i, t, x, y);
		}
	};

} // namespace MML

#endif // MML_IODESYSTEM_DAE_WITH_EVENTS_H

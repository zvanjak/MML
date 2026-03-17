///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        IODESystemWithEvents.h                                              ///
///  Description: Interface for ODE systems with event detection (zero-crossing)      ///
///               Supports direction filtering, terminal events, state modification   ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_IODE_SYSTEM_WITH_EVENTS_H
#define MML_IODE_SYSTEM_WITH_EVENTS_H

#include "mml/interfaces/IODESystem.h"

namespace MML {

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

} // namespace MML

#endif // MML_IODE_SYSTEM_WITH_EVENTS_H

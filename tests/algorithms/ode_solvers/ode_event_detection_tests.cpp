///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        ode_event_detection_tests.cpp                                       ///
///  Description: Tests for ODE event detection (IODESystemWithEvents)                ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

#include <catch2/catch_all.hpp>
#include "../../TestPrecision.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/ODESystem.h"
#include "mml/algorithms/ODESolvers/ODESolverEventDetection.h"
#endif

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

using namespace MML;

namespace MML::Tests::Algorithms::ODEEventDetectionTests
{
	//===================================================================================
	//                           Test ODE Systems with Events
	//===================================================================================

	/**
	 * @brief Bouncing ball - classic event detection example
	 * 
	 * Equations:
	 *   y' = v      (position derivative)
	 *   v' = -g     (velocity derivative, g = 9.81)
	 * 
	 * Event: y = 0 (ground contact) while falling (v < 0)
	 * Action: Reverse velocity with coefficient of restitution
	 */
	class BouncingBall : public IODESystemWithEvents {
	private:
		Real _g;           // Gravity
		Real _restitution; // Coefficient of restitution (0-1)
		
	public:
		BouncingBall(Real gravity = 9.81, Real restitution = 0.8)
			: _g(gravity), _restitution(restitution) {}
		
		int getDim() const override { return 2; }
		
		void derivs(Real t, const Vector<Real>& x, Vector<Real>& dxdt) const override {
			dxdt[0] = x[1];   // dy/dt = v
			dxdt[1] = -_g;    // dv/dt = -g
		}
		
		// Event functions
		int getNumEvents() const override { return 1; }
		
		Real eventFunction(int eventIndex, Real t, const Vector<Real>& x) const override {
			return x[0];  // y = 0 is ground
		}
		
		EventDirection getEventDirection(int eventIndex) const override {
			return EventDirection::Decreasing;  // Only when falling (y going from + to -)
		}
		
		EventAction getEventAction(int eventIndex) const override {
			return EventAction::Restart;  // Update velocity and continue
		}
		
		void handleEvent(int eventIndex, Real t, Vector<Real>& x) const override {
			x[0] = 0;  // Ensure exactly at ground
			x[1] = -_restitution * x[1];  // Reverse velocity with energy loss
		}
	};

	/**
	 * @brief Simple harmonic oscillator with zero-crossing detection
	 * 
	 * x' = v
	 * v' = -omega^2 * x
	 * 
	 * Events: Detect when x = 0 (crossing equilibrium)
	 */
	class HarmonicOscillatorWithEvents : public IODESystemWithEvents {
	private:
		Real _omega;  // Angular frequency
		
	public:
		explicit HarmonicOscillatorWithEvents(Real omega = 1.0) : _omega(omega) {}
		
		int getDim() const override { return 2; }
		
		void derivs(Real t, const Vector<Real>& x, Vector<Real>& dxdt) const override {
			dxdt[0] = x[1];              // dx/dt = v
			dxdt[1] = -_omega * _omega * x[0]; // dv/dt = -omega^2 * x
		}
		
		int getNumEvents() const override { return 1; }
		
		Real eventFunction(int eventIndex, Real t, const Vector<Real>& x) const override {
			return x[0];  // x = 0
		}
		
		EventDirection getEventDirection(int eventIndex) const override {
			return EventDirection::Both;  // Detect crossings in both directions
		}
		
		EventAction getEventAction(int eventIndex) const override {
			return EventAction::Continue;  // Just record, don't stop
		}
	};

	/**
	 * @brief System with terminal event
	 * 
	 * Exponential decay: x' = -x
	 * Event: Stop when x < threshold
	 */
	class ExponentialDecayWithTermination : public IODESystemWithEvents {
	private:
		Real _threshold;
		
	public:
		explicit ExponentialDecayWithTermination(Real threshold = 0.1) 
			: _threshold(threshold) {}
		
		int getDim() const override { return 1; }
		
		void derivs(Real t, const Vector<Real>& x, Vector<Real>& dxdt) const override {
			dxdt[0] = -x[0];
		}
		
		int getNumEvents() const override { return 1; }
		
		Real eventFunction(int eventIndex, Real t, const Vector<Real>& x) const override {
			return x[0] - _threshold;  // g = 0 when x = threshold
		}
		
		EventDirection getEventDirection(int eventIndex) const override {
			return EventDirection::Decreasing;  // x decreasing through threshold
		}
		
		EventAction getEventAction(int eventIndex) const override {
			return EventAction::Stop;  // Terminate integration
		}
	};

	/**
	 * @brief System with multiple events
	 * 
	 * 2D projectile: x' = vx, y' = vy, vx' = 0, vy' = -g
	 * Events:
	 *   0: y = 0 (ground impact)
	 *   1: y = y_max (apex detection when going up)
	 */
	class ProjectileWithMultipleEvents : public IODESystemWithEvents {
	private:
		Real _g;
		Real _yApex;
		
	public:
		ProjectileWithMultipleEvents(Real gravity = 9.81, Real apexHeight = 5.0)
			: _g(gravity), _yApex(apexHeight) {}
		
		int getDim() const override { return 4; }
		
		void derivs(Real t, const Vector<Real>& x, Vector<Real>& dxdt) const override {
			dxdt[0] = x[2];   // dx/dt = vx
			dxdt[1] = x[3];   // dy/dt = vy
			dxdt[2] = 0;      // dvx/dt = 0
			dxdt[3] = -_g;    // dvy/dt = -g
		}
		
		int getNumEvents() const override { return 2; }
		
		Real eventFunction(int eventIndex, Real t, const Vector<Real>& x) const override {
			if (eventIndex == 0) {
				return x[1];  // y = 0 (ground)
			} else {
				return x[1] - _yApex;  // y = y_apex
			}
		}
		
		EventDirection getEventDirection(int eventIndex) const override {
			if (eventIndex == 0) {
				return EventDirection::Decreasing;  // Falling toward ground
			} else {
				return EventDirection::Increasing;  // Rising through apex height
			}
		}
		
		EventAction getEventAction(int eventIndex) const override {
			if (eventIndex == 0) {
				return EventAction::Stop;  // Stop at ground
			} else {
				return EventAction::Continue;  // Just record apex crossing
			}
		}
	};

	//===================================================================================
	//                                    Tests
	//===================================================================================

	TEST_CASE("Event Detection - Bouncing Ball", "[ODEEventDetection]") {
		TEST_PRECISION_INFO();
		
		// Ball dropped from height 10 with zero initial velocity
		BouncingBall ball(9.81, 0.9);  // 10% energy loss per bounce
		Vector<Real> x0(std::vector<Real>{ 10.0, 0.0 });  // y=10, v=0
		
		DormandPrince5EventIntegrator integrator(ball);
		auto result = integrator.integrateWithEvents(ball, x0, 0.0, 5.0, 0.01, 1e-10, 1e-12);
		
		// Should have at least 2 bounces in 5 seconds (first bounce at ~1.43s)
		REQUIRE(result.events.size() >= 2);
		
		// Check first bounce
		SECTION("First bounce is accurate") {
			const auto& firstEvent = result.events[0];
			
			// For free fall from h=10: t = sqrt(2h/g) ≈ 1.428 s
			Real expectedTime = std::sqrt(2.0 * 10.0 / 9.81);
			REQUIRE_THAT(firstEvent.time, WithinRel((double)expectedTime, 1e-4));
			
			// Position should be at ground (≈ 0)
			REQUIRE_THAT(firstEvent.state[0], WithinAbs(0.0, 1e-10));
			
			// Velocity should be ≈ sqrt(2gh) = 14.0 m/s downward
			Real expectedVelocity = -std::sqrt(2.0 * 9.81 * 10.0);
			REQUIRE_THAT(firstEvent.state[1], WithinRel((double)expectedVelocity, 1e-4));
			
			// Direction should be decreasing
			REQUIRE(firstEvent.direction == EventDirection::Decreasing);
		}
		
		SECTION("Bounces decrease in height") {
			// After each bounce, max height should decrease by restitution^2
			// Since we're recording ground contact, check times between bounces
			if (result.events.size() >= 3) {
				Real t1 = result.events[0].time;
				Real t2 = result.events[1].time;
				Real t3 = result.events[2].time;
				
				Real period1 = t2 - t1;  // First bounce to second
				Real period2 = t3 - t2;  // Second to third
				
				// Period ratio should be approximately restitution coefficient
				Real ratio = period2 / period1;
				REQUIRE_THAT(ratio, WithinRel(0.9, 0.05));
			}
		}
	}

	TEST_CASE("Event Detection - Harmonic Oscillator Zero Crossings", "[ODEEventDetection]") {
		TEST_PRECISION_INFO();
		
		// omega = 2*pi => period = 1 second
		Real omega = 2.0 * Constants::PI;
		HarmonicOscillatorWithEvents osc(omega);
		Vector<Real> x0(std::vector<Real>{ 1.0, 0.0 });  // Start at max displacement
		
		DormandPrince5EventIntegrator integrator(osc);
		auto result = integrator.integrateWithEvents(osc, x0, 0.0, 2.0, 0.01, 1e-10, 1e-12);
		
		// In 2 periods, should have 4 zero crossings
		REQUIRE(result.events.size() == 4);
		
		SECTION("Zero crossings at expected times") {
			// Starting at x=1, first zero is at t = T/4 = 0.25
			// Then at 3T/4 = 0.75, 5T/4 = 1.25, 7T/4 = 1.75
			Real T = 1.0;  // Period
			
			REQUIRE_THAT(result.events[0].time, WithinRel(T / 4.0, 1e-4));
			REQUIRE_THAT(result.events[1].time, WithinRel(3.0 * T / 4.0, 1e-4));
			REQUIRE_THAT(result.events[2].time, WithinRel(5.0 * T / 4.0, 1e-4));
			REQUIRE_THAT(result.events[3].time, WithinRel(7.0 * T / 4.0, 1e-4));
		}
		
		SECTION("Alternating crossing directions") {
			// First crossing is decreasing (going from + to -)
			REQUIRE(result.events[0].direction == EventDirection::Decreasing);
			REQUIRE(result.events[1].direction == EventDirection::Increasing);
			REQUIRE(result.events[2].direction == EventDirection::Decreasing);
			REQUIRE(result.events[3].direction == EventDirection::Increasing);
		}
		
		SECTION("Position near zero at events") {
			for (const auto& evt : result.events) {
				REQUIRE_THAT(evt.state[0], WithinAbs(0.0, 1e-10));
			}
		}
	}

	TEST_CASE("Event Detection - Terminal Event (Stop)", "[ODEEventDetection]") {
		TEST_PRECISION_INFO();
		
		// x' = -x, x(0) = 1
		// Solution: x(t) = e^(-t)
		// Stop when x = 0.1 => t = -ln(0.1) ≈ 2.303
		Real threshold = 0.1;
		ExponentialDecayWithTermination decay(threshold);
		Vector<Real> x0(std::vector<Real>{ 1.0 });
		
		DormandPrince5EventIntegrator integrator(decay);
		auto result = integrator.integrateWithEvents(decay, x0, 0.0, 10.0, 0.1, 1e-10, 1e-12);
		
		// Should terminate before tEnd
		REQUIRE(result.terminatedByEvent == true);
		REQUIRE(result.events.size() == 1);
		
		SECTION("Termination time is accurate") {
			Real expectedTime = -std::log(threshold);  // ln(10) ≈ 2.303
			REQUIRE_THAT(result.events[0].time, WithinRel((double)expectedTime, 1e-4));
			REQUIRE_THAT(result.finalTime, WithinRel((double)expectedTime, 1e-4));
		}
		
		SECTION("State at termination") {
			REQUIRE_THAT(result.finalState[0], WithinRel((double)threshold, 1e-4));
		}
	}

	TEST_CASE("Event Detection - Multiple Events (Projectile)", "[ODEEventDetection]") {
		TEST_PRECISION_INFO();
		
		// Projectile launched at 45 degrees with v0 = 20 m/s
		Real v0 = 20.0;
		Real angle = Constants::PI / 4.0;  // 45 degrees
		Real g = 9.81;
		Real apexHeight = 5.0;  // Detect when y crosses 5m
		
		ProjectileWithMultipleEvents proj(g, apexHeight);
		Vector<Real> x0(std::vector<Real>{ 
			0.0,                      // x = 0
			0.0,                      // y = 0
			v0 * std::cos(angle),     // vx
			v0 * std::sin(angle)      // vy
		});
		
		DormandPrince5EventIntegrator integrator(proj);
		auto result = integrator.integrateWithEvents(proj, x0, 0.0, 10.0, 0.01, 1e-10, 1e-12);
		
		// Should detect apex crossing and ground impact
		REQUIRE(result.events.size() >= 2);
		REQUIRE(result.terminatedByEvent == true);  // Stops at ground
		
		SECTION("Apex height crossing detected") {
			// Find event 1 (apex crossing)
			auto apexEvent = std::find_if(result.events.begin(), result.events.end(),
				[](const EventInfo& e) { return e.eventIndex == 1; });
			
			REQUIRE(apexEvent != result.events.end());
			REQUIRE_THAT(apexEvent->state[1], WithinAbs(apexHeight, 1e-6));
			REQUIRE(apexEvent->direction == EventDirection::Increasing);
		}
		
		SECTION("Ground impact is final event") {
			const auto& lastEvent = result.events.back();
			REQUIRE(lastEvent.eventIndex == 0);  // Ground event
			REQUIRE_THAT(lastEvent.state[1], WithinAbs(0.0, 1e-10));
			
			// Flight time: t = 2 * v0 * sin(angle) / g
			Real expectedFlightTime = 2.0 * v0 * std::sin(angle) / g;
			REQUIRE_THAT(lastEvent.time, WithinRel((double)expectedFlightTime, 1e-3));
		}
	}

	TEST_CASE("Event Detection - Direction Filtering", "[ODEEventDetection]") {
		TEST_PRECISION_INFO();
		
		// Test that direction filtering works correctly
		// Use harmonic oscillator but only detect increasing crossings
		class HarmonicOscillatorIncreasingOnly : public HarmonicOscillatorWithEvents {
		public:
			HarmonicOscillatorIncreasingOnly() : HarmonicOscillatorWithEvents(2.0 * Constants::PI) {}
			
			EventDirection getEventDirection(int) const override {
				return EventDirection::Increasing;  // Only detect when x goes from - to +
			}
		};
		
		HarmonicOscillatorIncreasingOnly osc;
		Vector<Real> x0(std::vector<Real>{ 1.0, 0.0 });
		
		DormandPrince5EventIntegrator integrator(osc);
		auto result = integrator.integrateWithEvents(osc, x0, 0.0, 2.0, 0.01, 1e-10, 1e-12);
		
		// Should only have 2 events (one per period, increasing only)
		REQUIRE(result.events.size() == 2);
		
		// All events should be increasing direction
		for (const auto& evt : result.events) {
			REQUIRE(evt.direction == EventDirection::Increasing);
		}
		
		// Times should be at 3T/4 and 7T/4 (when x goes from - to +)
		Real T = 1.0;
		REQUIRE_THAT(result.events[0].time, WithinRel(3.0 * T / 4.0, 1e-4));
		REQUIRE_THAT(result.events[1].time, WithinRel(7.0 * T / 4.0, 1e-4));
	}

	TEST_CASE("Event Detection - No Events", "[ODEEventDetection]") {
		TEST_PRECISION_INFO();
		
		// Exponential decay that never reaches threshold
		ExponentialDecayWithTermination decay(0.001);  // Very low threshold
		Vector<Real> x0(std::vector<Real>{ 1.0 });
		
		DormandPrince5EventIntegrator integrator(decay);
		auto result = integrator.integrateWithEvents(decay, x0, 0.0, 2.0, 0.1, 1e-10, 1e-12);
		
		// x(2) = e^(-2) ≈ 0.135, still above threshold 0.001
		REQUIRE(result.events.empty());
		REQUIRE(result.terminatedByEvent == false);
		REQUIRE_THAT(result.finalTime, WithinAbs(2.0, 1e-10));
	}

	TEST_CASE("Event Detection - Different Steppers", "[ODEEventDetection]") {
		TEST_PRECISION_INFO();
		
		// Test that event detection works with different steppers
		BouncingBall ball(9.81, 0.9);
		Vector<Real> x0(std::vector<Real>{ 5.0, 0.0 });
		
		SECTION("DormandPrince5") {
			DormandPrince5EventIntegrator integrator(ball);
			auto result = integrator.integrateWithEvents(ball, x0, 0.0, 3.0, 0.01, 1e-10, 1e-12);
			REQUIRE(result.events.size() >= 2);
		}
		
		SECTION("DormandPrince8") {
			DormandPrince8EventIntegrator integrator(ball);
			auto result = integrator.integrateWithEvents(ball, x0, 0.0, 3.0, 0.01, 1e-10, 1e-12);
			REQUIRE(result.events.size() >= 2);
		}
		
		SECTION("CashKarp") {
			CashKarpEventIntegrator integrator(ball);
			auto result = integrator.integrateWithEvents(ball, x0, 0.0, 3.0, 0.01, 1e-10, 1e-12);
			REQUIRE(result.events.size() >= 2);
		}
	}

} // namespace MML::Tests::Algorithms::ODEEventDetectionTests

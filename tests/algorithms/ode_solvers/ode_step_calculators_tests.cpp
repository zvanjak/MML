#include <catch2/catch_all.hpp>
#include "../../TestPrecision.h"

#include <cmath>
#include <limits>

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "mml/algorithms/ODESolvers/ODESystemStepCalculators.h"
#include "interfaces/IODESystem.h"
#endif

using namespace MML;
using namespace MML::Testing;

namespace MML::Tests::Algorithms::ODEStepCalculatorTests
{
	/////////////////////////////////////////////////////////////////////////////////////
	///                         TEST ODE SYSTEMS                                      ///
	/////////////////////////////////////////////////////////////////////////////////////

	// Free particle: dx/dt = v, dv/dt = 0
	// Solution: x(t) = x0 + v0*t, v(t) = v0
	class FreeParticle1D final : public IODESystem
	{
	public:
		int getDim() const override { return 2; }

		void derivs(const Real /*t*/, const Vector<Real>& x, Vector<Real>& dxdt) const override
		{
			dxdt[0] = x[1];       // dx/dt = v
			dxdt[1] = REAL(0.0);  // dv/dt = 0
		}
	};

	// Simple harmonic oscillator: dx/dt = v, dv/dt = -x
	// Solution: x(t) = A*cos(t + φ), v(t) = -A*sin(t + φ)
	class SimpleHarmonicOscillator final : public IODESystem
	{
	public:
		int getDim() const override { return 2; }

		void derivs(const Real /*t*/, const Vector<Real>& x, Vector<Real>& dxdt) const override
		{
			dxdt[0] = x[1];        // dx/dt = v
			dxdt[1] = -x[0];       // dv/dt = -x (omega = 1)
		}
	};

	// Exponential decay: dy/dt = -y
	// Solution: y(t) = y0 * exp(-t)
	class ExponentialDecay final : public IODESystem
	{
	public:
		int getDim() const override { return 1; }

		void derivs(const Real /*t*/, const Vector<Real>& y, Vector<Real>& dydt) const override
		{
			dydt[0] = -y[0];
		}
	};

	// Exponential growth: dy/dt = y
	// Solution: y(t) = y0 * exp(t)
	class ExponentialGrowth final : public IODESystem
	{
	public:
		int getDim() const override { return 1; }

		void derivs(const Real /*t*/, const Vector<Real>& y, Vector<Real>& dydt) const override
		{
			dydt[0] = y[0];
		}
	};

	// Polynomial ODE: dy/dt = t^2
	// Solution: y(t) = y0 + t^3/3
	class PolynomialODE final : public IODESystem
	{
	public:
		int getDim() const override { return 1; }

		void derivs(const Real t, const Vector<Real>& /*y*/, Vector<Real>& dydt) const override
		{
			dydt[0] = t * t;
		}
	};

	// Coupled linear system: dx/dt = y, dy/dt = -x
	// (Same as SHO but framed differently for testing)
	class CoupledLinearSystem final : public IODESystem
	{
	public:
		int getDim() const override { return 2; }

		void derivs(const Real /*t*/, const Vector<Real>& x, Vector<Real>& dxdt) const override
		{
			dxdt[0] = x[1];
			dxdt[1] = -x[0];
		}
	};

	// Falling body with gravity: dx/dt = v, dv/dt = -g
	// Solution: x(t) = x0 + v0*t - 0.5*g*t^2, v(t) = v0 - g*t
	class FallingBody final : public IODESystem
	{
		Real _g;
	public:
		FallingBody(Real g = REAL(9.81)) : _g(g) {}
		int getDim() const override { return 2; }

		void derivs(const Real /*t*/, const Vector<Real>& x, Vector<Real>& dxdt) const override
		{
			dxdt[0] = x[1];    // dx/dt = v
			dxdt[1] = -_g;     // dv/dt = -g
		}
	};

	/////////////////////////////////////////////////////////////////////////////////////
	///                    HELPER: GENERIC STEP TESTER                                ///
	/////////////////////////////////////////////////////////////////////////////////////

	template <typename StepCalc>
	void TestStepCalculatorContract(const StepCalc& stepCalc, const char* name)
	{
		DYNAMIC_SECTION(name << " - basic contract")
		{
			TEST_PRECISION_INFO();

			FreeParticle1D sys;
			const Real t = REAL(1.23);
			const Real h = REAL(0.1);

			Vector<Real> x_start(2);
			x_start[0] = REAL(1.25);
			x_start[1] = REAL(-0.5);

			Vector<Real> dxdt(2);
			sys.derivs(t, x_start, dxdt);

			Vector<Real> x_out(2, std::numeric_limits<Real>::quiet_NaN());
			Vector<Real> x_err(2, std::numeric_limits<Real>::quiet_NaN());

			stepCalc.calcStep(sys, t, x_start, dxdt, h, x_out, x_err);

			const Vector<Real> expected{ x_start[0] + x_start[1] * h, x_start[1] };

			for (int i = 0; i < expected.size(); ++i)
			{
				REQUIRE(std::isfinite(x_out[i]));
				REQUIRE(std::isfinite(x_err[i]));
				REQUIRE_THAT(x_out[i], Catch::Matchers::WithinAbs(expected[i], REAL(1e-12)));
				REQUIRE_THAT(x_err[i], Catch::Matchers::WithinAbs(REAL(0.0), REAL(1e-12)));
			}
		}
	}

	template <typename StepCalc>
	void TestExponentialDecay(const StepCalc& stepCalc, const char* name, Real expectedTol)
	{
		DYNAMIC_SECTION(name << " - exponential decay")
		{
			TEST_PRECISION_INFO();

			ExponentialDecay sys;
			Real t = REAL(0.0);
			Real h = REAL(0.01);
			Real y0 = REAL(1.0);

			Vector<Real> y(1);
			y[0] = y0;
			Vector<Real> dydt(1);
			Vector<Real> y_out(1);
			Vector<Real> y_err(1);

			// Integrate for 100 steps to t=1.0
			for (int step = 0; step < 100; ++step)
			{
				sys.derivs(t, y, dydt);
				stepCalc.calcStep(sys, t, y, dydt, h, y_out, y_err);
				y = y_out;
				t += h;
			}

			Real expected = y0 * std::exp(-REAL(1.0));
			REQUIRE_THAT(y[0], Catch::Matchers::WithinRel(expected, expectedTol));
		}
	}

	template <typename StepCalc>
	void TestHarmonicOscillator(const StepCalc& stepCalc, const char* name, Real expectedTol)
	{
		DYNAMIC_SECTION(name << " - harmonic oscillator")
		{
			TEST_PRECISION_INFO();

			SimpleHarmonicOscillator sys;
			Real t = REAL(0.0);
			Real h = REAL(0.01);

			Vector<Real> x(2);
			x[0] = REAL(1.0);  // x(0) = 1
			x[1] = REAL(0.0);  // v(0) = 0, so A=1, φ=0

			Vector<Real> dxdt(2);
			Vector<Real> x_out(2);
			Vector<Real> x_err(2);

			// Integrate for one period (2*pi)
			int steps = static_cast<int>(2.0 * Constants::PI / h);
			for (int step = 0; step < steps; ++step)
			{
				sys.derivs(t, x, dxdt);
				stepCalc.calcStep(sys, t, x, dxdt, h, x_out, x_err);
				x = x_out;
				t += h;
			}

			// After one period, should return to initial conditions
			// Velocity tolerance is larger because it accumulates phase error differently
			REQUIRE_THAT(x[0], Catch::Matchers::WithinRel(REAL(1.0), expectedTol));
			REQUIRE_THAT(x[1], Catch::Matchers::WithinAbs(REAL(0.0), REAL(10.0) * expectedTol));
		}
	}

	template <typename StepCalc>
	void TestPolynomialODE(const StepCalc& stepCalc, const char* name, Real expectedTol)
	{
		DYNAMIC_SECTION(name << " - polynomial dy/dt = t^2")
		{
			TEST_PRECISION_INFO();

			PolynomialODE sys;
			Real t = REAL(0.0);
			Real h = REAL(0.01);

			Vector<Real> y(1);
			y[0] = REAL(0.0);  // y(0) = 0

			Vector<Real> dydt(1);
			Vector<Real> y_out(1);
			Vector<Real> y_err(1);

			// Integrate to t=1.0
			int steps = static_cast<int>(REAL(1.0) / h);
			for (int step = 0; step < steps; ++step)
			{
				sys.derivs(t, y, dydt);
				stepCalc.calcStep(sys, t, y, dydt, h, y_out, y_err);
				y = y_out;
				t += h;
			}

			// y(1) = 1^3/3 = 1/3
			Real expected = REAL(1.0) / REAL(3.0);
			REQUIRE_THAT(y[0], Catch::Matchers::WithinRel(expected, expectedTol));
		}
	}

	template <typename StepCalc>
	void TestFallingBody(const StepCalc& stepCalc, const char* name, Real expectedTol)
	{
		DYNAMIC_SECTION(name << " - falling body")
		{
			TEST_PRECISION_INFO();

			Real g = REAL(10.0);  // Simple gravity for easy calculation
			FallingBody sys(g);
			Real t = REAL(0.0);
			Real h = REAL(0.001);  // Small step for accuracy

			Vector<Real> x(2);
			x[0] = REAL(100.0);  // x0 = 100
			x[1] = REAL(0.0);    // v0 = 0

			Vector<Real> dxdt(2);
			Vector<Real> x_out(2);
			Vector<Real> x_err(2);

			// Integrate for 1 second
			int steps = static_cast<int>(REAL(1.0) / h);
			for (int step = 0; step < steps; ++step)
			{
				sys.derivs(t, x, dxdt);
				stepCalc.calcStep(sys, t, x, dxdt, h, x_out, x_err);
				x = x_out;
				t += h;
			}

			// x(1) = 100 + 0*1 - 0.5*10*1^2 = 100 - 5 = 95
			// v(1) = 0 - 10*1 = -10
			REQUIRE_THAT(x[0], Catch::Matchers::WithinRel(REAL(95.0), expectedTol));
			REQUIRE_THAT(x[1], Catch::Matchers::WithinRel(REAL(-10.0), expectedTol));
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////
	///                         EULER STEP CALCULATOR                                 ///
	/////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("Euler step calculator - basic contract", "[ODE][stepcalc][Euler]")
	{
		TestStepCalculatorContract(StepCalculators::EulerStepCalc, "Euler");
	}

	TEST_CASE("Euler step calculator - exponential decay", "[ODE][stepcalc][Euler]")
	{
		TestExponentialDecay(StepCalculators::EulerStepCalc, "Euler", REAL(0.01)); // 1% tolerance (1st order)
	}

	TEST_CASE("Euler step calculator - polynomial", "[ODE][stepcalc][Euler]")
	{
		TestPolynomialODE(StepCalculators::EulerStepCalc, "Euler", REAL(0.05)); // 5% tolerance (1st order, accumulated)
	}

	TEST_CASE("Euler step calculator - falling body", "[ODE][stepcalc][Euler]")
	{
		TestFallingBody(StepCalculators::EulerStepCalc, "Euler", REAL(0.001));
	}

	/////////////////////////////////////////////////////////////////////////////////////
	///                      EULER-CROMER STEP CALCULATOR                             ///
	/////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("Euler-Cromer step calculator - basic contract", "[ODE][stepcalc][EulerCromer]")
	{
		TestStepCalculatorContract(StepCalculators::EulerCromerStepCalc, "Euler-Cromer");
	}

	TEST_CASE("Euler-Cromer step calculator - exponential decay", "[ODE][stepcalc][EulerCromer]")
	{
		TestExponentialDecay(StepCalculators::EulerCromerStepCalc, "Euler-Cromer", REAL(0.001)); // Better than Euler
	}

	TEST_CASE("Euler-Cromer step calculator - provides error estimate", "[ODE][stepcalc][EulerCromer]")
	{
		TEST_PRECISION_INFO();

		ExponentialDecay sys;
		Vector<Real> y(1);
		y[0] = REAL(1.0);
		Vector<Real> dydt(1);
		Vector<Real> y_out(1);
		Vector<Real> y_err(1);

		sys.derivs(REAL(0.0), y, dydt);
		StepCalculators::EulerCromerStepCalc.calcStep(sys, REAL(0.0), y, dydt, REAL(0.1), y_out, y_err);

		// Euler-Cromer should provide non-zero error estimate
		REQUIRE(std::abs(y_err[0]) > REAL(0.0));
	}

	/////////////////////////////////////////////////////////////////////////////////////
	///                     VELOCITY VERLET STEP CALCULATOR                           ///
	/////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("Velocity Verlet step calculator - basic contract", "[ODE][stepcalc][VelocityVerlet]")
	{
		TestStepCalculatorContract(StepCalculators::VelocityVerletStepCalc, "VelocityVerlet");
	}

	TEST_CASE("Velocity Verlet step calculator - harmonic oscillator energy conservation", "[ODE][stepcalc][VelocityVerlet][symplectic]")
	{
		TEST_PRECISION_INFO();

		SimpleHarmonicOscillator sys;
		Real t = REAL(0.0);
		Real h = REAL(0.01);

		Vector<Real> x(2);
		x[0] = REAL(1.0);
		x[1] = REAL(0.0);

		Vector<Real> dxdt(2);
		Vector<Real> x_out(2);
		Vector<Real> x_err(2);

		// Initial energy: E = 0.5*(v^2 + x^2) = 0.5*(0 + 1) = 0.5
		Real E0 = REAL(0.5) * (x[0] * x[0] + x[1] * x[1]);

		// Integrate for many periods
		int steps = static_cast<int>(10.0 * 2.0 * Constants::PI / h);
		Real maxEnergyDeviation = REAL(0.0);

		for (int step = 0; step < steps; ++step)
		{
			sys.derivs(t, x, dxdt);
			StepCalculators::VelocityVerletStepCalc.calcStep(sys, t, x, dxdt, h, x_out, x_err);
			x = x_out;
			t += h;

			Real E = REAL(0.5) * (x[0] * x[0] + x[1] * x[1]);
			maxEnergyDeviation = std::max(maxEnergyDeviation, std::abs(E - E0));
		}

		// Symplectic integrators conserve energy well
		REQUIRE(maxEnergyDeviation < REAL(0.01));
	}

	/////////////////////////////////////////////////////////////////////////////////////
	///                      LEAPFROG STEP CALCULATOR                                 ///
	/////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("Leapfrog step calculator - basic contract", "[ODE][stepcalc][Leapfrog]")
	{
		TestStepCalculatorContract(StepCalculators::LeapfrogStepCalc, "Leapfrog");
	}

	TEST_CASE("Leapfrog step calculator - harmonic oscillator energy conservation", "[ODE][stepcalc][Leapfrog][symplectic]")
	{
		TEST_PRECISION_INFO();

		SimpleHarmonicOscillator sys;
		Real t = REAL(0.0);
		Real h = REAL(0.01);

		Vector<Real> x(2);
		x[0] = REAL(1.0);
		x[1] = REAL(0.0);

		Vector<Real> dxdt(2);
		Vector<Real> x_out(2);
		Vector<Real> x_err(2);

		Real E0 = REAL(0.5) * (x[0] * x[0] + x[1] * x[1]);

		int steps = static_cast<int>(10.0 * 2.0 * Constants::PI / h);
		Real maxEnergyDeviation = REAL(0.0);

		for (int step = 0; step < steps; ++step)
		{
			sys.derivs(t, x, dxdt);
			StepCalculators::LeapfrogStepCalc.calcStep(sys, t, x, dxdt, h, x_out, x_err);
			x = x_out;
			t += h;

			Real E = REAL(0.5) * (x[0] * x[0] + x[1] * x[1]);
			maxEnergyDeviation = std::max(maxEnergyDeviation, std::abs(E - E0));
		}

		REQUIRE(maxEnergyDeviation < REAL(0.01));
	}

	/////////////////////////////////////////////////////////////////////////////////////
	///                      MIDPOINT STEP CALCULATOR                                 ///
	/////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("Midpoint step calculator - basic contract", "[ODE][stepcalc][Midpoint]")
	{
		TestStepCalculatorContract(StepCalculators::MidpointStepCalc, "Midpoint");
	}

	TEST_CASE("Midpoint step calculator - exponential decay", "[ODE][stepcalc][Midpoint]")
	{
		TestExponentialDecay(StepCalculators::MidpointStepCalc, "Midpoint", REAL(0.0001)); // 2nd order
	}

	TEST_CASE("Midpoint step calculator - polynomial", "[ODE][stepcalc][Midpoint]")
	{
		TestPolynomialODE(StepCalculators::MidpointStepCalc, "Midpoint", REAL(0.0001));
	}

	/////////////////////////////////////////////////////////////////////////////////////
	///                         RK4 STEP CALCULATOR                                   ///
	/////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("RK4 step calculator - basic contract", "[ODE][stepcalc][RK4]")
	{
		TestStepCalculatorContract(StepCalculators::RK4_Basic, "RK4");
	}

	TEST_CASE("RK4 step calculator - exponential decay", "[ODE][stepcalc][RK4]")
	{
		TestExponentialDecay(StepCalculators::RK4_Basic, "RK4", REAL(1e-6)); // 4th order - very accurate
	}

	TEST_CASE("RK4 step calculator - harmonic oscillator", "[ODE][stepcalc][RK4]")
	{
		TestHarmonicOscillator(StepCalculators::RK4_Basic, "RK4", REAL(1e-2)); // Accumulated over 628 steps
	}

	TEST_CASE("RK4 step calculator - polynomial", "[ODE][stepcalc][RK4]")
	{
		TestPolynomialODE(StepCalculators::RK4_Basic, "RK4", REAL(1e-8)); // Exact for polynomials degree <= 4
	}

	TEST_CASE("RK4 step calculator - falling body", "[ODE][stepcalc][RK4]")
	{
		TestFallingBody(StepCalculators::RK4_Basic, "RK4", REAL(1e-8));
	}

	/////////////////////////////////////////////////////////////////////////////////////
	///                     RK5 CASH-KARP STEP CALCULATOR                             ///
	/////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("RK5 Cash-Karp step calculator - basic contract", "[ODE][stepcalc][RK5CashKarp]")
	{
		TestStepCalculatorContract(StepCalculators::RK5_CashKarp, "RK5-CashKarp");
	}

	TEST_CASE("RK5 Cash-Karp step calculator - exponential decay", "[ODE][stepcalc][RK5CashKarp]")
	{
		TestExponentialDecay(StepCalculators::RK5_CashKarp, "RK5-CashKarp", REAL(1e-8)); // 5th order
	}

	TEST_CASE("RK5 Cash-Karp step calculator - harmonic oscillator", "[ODE][stepcalc][RK5CashKarp]")
	{
		TestHarmonicOscillator(StepCalculators::RK5_CashKarp, "RK5-CashKarp", REAL(1e-2)); // Accumulated over 628 steps
	}

	TEST_CASE("RK5 Cash-Karp step calculator - provides error estimate", "[ODE][stepcalc][RK5CashKarp]")
	{
		TEST_PRECISION_INFO();

		ExponentialDecay sys;
		Vector<Real> y(1);
		y[0] = REAL(1.0);
		Vector<Real> dydt(1);
		Vector<Real> y_out(1);
		Vector<Real> y_err(1);

		sys.derivs(REAL(0.0), y, dydt);
		StepCalculators::RK5_CashKarp.calcStep(sys, REAL(0.0), y, dydt, REAL(0.1), y_out, y_err);

		// Should provide non-zero error estimate
		REQUIRE(std::abs(y_err[0]) > REAL(0.0));
		// Error estimate should be O(h^5) so relatively small
		REQUIRE(std::abs(y_err[0]) < REAL(1e-4));
	}

	/////////////////////////////////////////////////////////////////////////////////////
	///                   DORMAND-PRINCE 5 STEP CALCULATOR                            ///
	/////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("Dormand-Prince 5 step calculator - basic contract", "[ODE][stepcalc][DormandPrince5]")
	{
		TestStepCalculatorContract(StepCalculators::DormandPrince5StepCalc, "DormandPrince5");
	}

	TEST_CASE("Dormand-Prince 5 step calculator - exponential decay", "[ODE][stepcalc][DormandPrince5]")
	{
		TestExponentialDecay(StepCalculators::DormandPrince5StepCalc, "DormandPrince5", REAL(1e-8));
	}

	TEST_CASE("Dormand-Prince 5 step calculator - harmonic oscillator", "[ODE][stepcalc][DormandPrince5]")
	{
		TestHarmonicOscillator(StepCalculators::DormandPrince5StepCalc, "DormandPrince5", REAL(1e-2)); // Accumulated over 628 steps
	}

	TEST_CASE("Dormand-Prince 5 step calculator - provides error estimate", "[ODE][stepcalc][DormandPrince5]")
	{
		TEST_PRECISION_INFO();

		SimpleHarmonicOscillator sys;
		Vector<Real> x(2);
		x[0] = REAL(1.0);
		x[1] = REAL(0.0);
		Vector<Real> dxdt(2);
		Vector<Real> x_out(2);
		Vector<Real> x_err(2);

		sys.derivs(REAL(0.0), x, dxdt);
		StepCalculators::DormandPrince5StepCalc.calcStep(sys, REAL(0.0), x, dxdt, REAL(0.1), x_out, x_err);

		// Should provide non-zero error estimate
		REQUIRE((std::abs(x_err[0]) > REAL(0.0) || std::abs(x_err[1]) > REAL(0.0)));
	}

	/////////////////////////////////////////////////////////////////////////////////////
	///                   DORMAND-PRINCE 8 STEP CALCULATOR                            ///
	/////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("Dormand-Prince 8 step calculator - basic contract", "[ODE][stepcalc][DormandPrince8]")
	{
		TestStepCalculatorContract(StepCalculators::DormandPrince8StepCalc, "DormandPrince8");
	}

	TEST_CASE("Dormand-Prince 8 step calculator - exponential decay", "[ODE][stepcalc][DormandPrince8]")
	{
		TestExponentialDecay(StepCalculators::DormandPrince8StepCalc, "DormandPrince8", REAL(1e-10)); // 8th order - very accurate
	}

	TEST_CASE("Dormand-Prince 8 step calculator - harmonic oscillator", "[ODE][stepcalc][DormandPrince8]")
	{
		TestHarmonicOscillator(StepCalculators::DormandPrince8StepCalc, "DormandPrince8", REAL(1e-2)); // Accumulated over 628 steps
	}

	TEST_CASE("Dormand-Prince 8 step calculator - polynomial", "[ODE][stepcalc][DormandPrince8]")
	{
		TestPolynomialODE(StepCalculators::DormandPrince8StepCalc, "DormandPrince8", REAL(1e-12));
	}

	TEST_CASE("Dormand-Prince 8 step calculator - provides error estimate", "[ODE][stepcalc][DormandPrince8]")
	{
		TEST_PRECISION_INFO();

		ExponentialDecay sys;
		Vector<Real> y(1);
		y[0] = REAL(1.0);
		Vector<Real> dydt(1);
		Vector<Real> y_out(1);
		Vector<Real> y_err(1);

		sys.derivs(REAL(0.0), y, dydt);
		StepCalculators::DormandPrince8StepCalc.calcStep(sys, REAL(0.0), y, dydt, REAL(0.1), y_out, y_err);

		// Should provide non-zero error estimate
		REQUIRE(std::abs(y_err[0]) > REAL(0.0));
		// Error should be very small for 8th order
		REQUIRE(std::abs(y_err[0]) < REAL(1e-8));
	}

	/////////////////////////////////////////////////////////////////////////////////////
	///                        ORDER VERIFICATION TESTS                               ///
	/////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("Step calculators - order verification via step halving", "[ODE][stepcalc][order]")
	{
		TEST_PRECISION_INFO();

		ExponentialDecay sys;
		Real t0 = REAL(0.0);
		Real y0 = REAL(1.0);
		Real h = REAL(0.1);

		// Expected solution at t=h
		Real exact = y0 * std::exp(-h);

		// Test that higher-order methods have smaller errors
		Vector<Real> y(1), dydt(1), y_out(1), y_err(1);

		// Euler (1st order)
		y[0] = y0;
		sys.derivs(t0, y, dydt);
		StepCalculators::EulerStepCalc.calcStep(sys, t0, y, dydt, h, y_out, y_err);
		Real eulerErr = std::abs(y_out[0] - exact);

		// Midpoint (2nd order)
		y[0] = y0;
		sys.derivs(t0, y, dydt);
		StepCalculators::MidpointStepCalc.calcStep(sys, t0, y, dydt, h, y_out, y_err);
		Real midpointErr = std::abs(y_out[0] - exact);

		// RK4 (4th order)
		y[0] = y0;
		sys.derivs(t0, y, dydt);
		StepCalculators::RK4_Basic.calcStep(sys, t0, y, dydt, h, y_out, y_err);
		Real rk4Err = std::abs(y_out[0] - exact);

		// Dormand-Prince 5 (5th order)
		y[0] = y0;
		sys.derivs(t0, y, dydt);
		StepCalculators::DormandPrince5StepCalc.calcStep(sys, t0, y, dydt, h, y_out, y_err);
		Real dp5Err = std::abs(y_out[0] - exact);

		// Dormand-Prince 8 (8th order)
		y[0] = y0;
		sys.derivs(t0, y, dydt);
		StepCalculators::DormandPrince8StepCalc.calcStep(sys, t0, y, dydt, h, y_out, y_err);
		Real dp8Err = std::abs(y_out[0] - exact);

		// Higher order should mean smaller error
		REQUIRE(midpointErr < eulerErr);
		REQUIRE(rk4Err < midpointErr);
		REQUIRE(dp5Err < rk4Err);
		REQUIRE(dp8Err < dp5Err);
	}

	/////////////////////////////////////////////////////////////////////////////////////
	///                     ALL STEPPERS CONTRACT TEST                                ///
	/////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("All step calculators - unified contract test", "[ODE][stepcalc][AllSteppers]")
	{
		TestStepCalculatorContract(StepCalculators::EulerStepCalc, "Euler");
		TestStepCalculatorContract(StepCalculators::EulerCromerStepCalc, "Euler-Cromer");
		TestStepCalculatorContract(StepCalculators::VelocityVerletStepCalc, "VelocityVerlet");
		TestStepCalculatorContract(StepCalculators::LeapfrogStepCalc, "Leapfrog");
		TestStepCalculatorContract(StepCalculators::MidpointStepCalc, "Midpoint");
		TestStepCalculatorContract(StepCalculators::RK4_Basic, "RK4");
		TestStepCalculatorContract(StepCalculators::RK5_CashKarp, "RK5 Cash-Karp");
		TestStepCalculatorContract(StepCalculators::DormandPrince5StepCalc, "Dormand-Prince 5");
		TestStepCalculatorContract(StepCalculators::DormandPrince8StepCalc, "Dormand-Prince 8");
	}

} // namespace MML::Tests::Algorithms::ODEStepCalculatorTests

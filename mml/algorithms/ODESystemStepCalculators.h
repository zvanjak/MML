///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        ODESystemStepCalculators.h                                          ///
///  Description: ODE step calculators (Euler, RK4, RK5, Cash-Karp, Dormand-Prince)   ///
///               Single-step methods for ODE integration                             ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_ODE_SYSTEM_STEP_CALCULATORS_H
#define MML_ODE_SYSTEM_STEP_CALCULATORS_H

#include "MMLBase.h"

#include "interfaces/IODESystem.h"
#include "interfaces/IODESystemStepCalculator.h"

#include "base/ODESystem.h"

// NOTE: These calculators are stateless single-step implementations used by
// `ODESystemFixedStepSolver` for fixed-step integration. They intentionally
// implement the simple one-step interface (`IODESystemStepCalculator`) and
// therefore duplicate some coefficients that are also present in the
// adaptive steppers (implemented in `ODEAdaptiveIntegrator.h`).
//
// If you change coefficients here, ensure the adaptive steppers remain
// consistent where intended. The duplication is deliberate to keep the
// fixed-step and adaptive implementations independent and easy to reason
// about.

namespace MML
{
	// For a given IODESystem of dimension n, and given initial values for the variables x_start[0..n-1] 
	// and their derivatives dxdt[0..n-1] known at t, uses the Euler method to advance the solution 
	// over an interval h and return the incremented variables as x_out[0..n-1].
	class EulerStep_Calculator : public IODESystemStepCalculator
	{
	public:
		void calcStep(const IODESystem& odeSystem,
									const Real t, const Vector<Real>& x_start, const Vector<Real>& dxdt,
									const Real h, Vector<Real>& x_out, Vector<Real>& x_err_out) const override
		{
			int i, n = odeSystem.getDim();
			
			for (i = 0; i < n; i++)
				x_out[i] = x_start[i] + h * dxdt[i];
		}
	};

	// For a given IODESystem of dimension n, and given initial values for the variables x_start[0..n-1] 
	// and their derivatives dxdt[0..n-1] known at t, uses the Euler-Cromer method to advance the solution 
	// over an interval h and return the incremented variables as x_out[0..n-1].
	class EulerCromer_StepCalculator : public IODESystemStepCalculator
	{
	public:
		void calcStep(const IODESystem& odeSystem,
									const Real t, const Vector<Real>& x_start, const Vector<Real>& dxdt,
									const Real h, Vector<Real>& x_out, Vector<Real>& x_err_out) const override
		{
			int i, n = odeSystem.getDim();
			Vector<Real> x_new(n);

			for (i = 0; i < n; i++)
				x_new[i] = x_start[i] + h * dxdt[i];

			odeSystem.derivs(t + h, x_new, x_out);
		}
	};

	// For a given IODESystem of dimension n, and given initial values for the variables x_start[0..n-1] 
	// and their derivatives dxdt[0..n-1] known at t, uses the Velocity Verlet method to advance the solution 
	// over an interval h and return the incremented variables as x_out[0..n-1].
	class VelocityVerlet_StepCalculator : public IODESystemStepCalculator
	{
	public:
		void calcStep(const IODESystem& odeSystem,
									const Real t, const Vector<Real>& x_start, const Vector<Real>& dxdt,
									const Real h, Vector<Real>& x_out, Vector<Real>& x_err_out) const override
		{
			int n = odeSystem.getDim();
			int half_n = n / 2;

			// Split state
			// x_start[0..half_n-1] = positions
			// x_start[half_n..n-1] = velocities
			// dxdt[0..half_n-1] = velocities
			// dxdt[half_n..n-1] = accelerations

			// 1. Update positions
			for (int i = 0; i < half_n; ++i)
				x_out[i] = x_start[i] + x_start[half_n + i] * h + 0.5 * dxdt[half_n + i] * h * h;

			// 2. Compute new acceleration at new position
			Vector<Real> x_temp = x_out;
			for (int i = 0; i < half_n; ++i)
				x_temp[half_n + i] = x_start[half_n + i]; // velocities (will be updated)
			Vector<Real> dxdt_temp(n);
			odeSystem.derivs(t + h, x_temp, dxdt_temp);

			// 3. Update velocities
			for (int i = 0; i < half_n; ++i)
				x_out[half_n + i] = x_start[half_n + i] + 0.5 * (dxdt[half_n + i] + dxdt_temp[half_n + i]) * h;

			// No error estimate
			for (int i = 0; i < n; ++i)
				x_err_out[i] = 0.0;
		}
	};

	// Leapfrog (Velocity Verlet) step calculator for symplectic integration.
	// Best suited for Hamiltonian systems where state is split into positions and velocities.
	// Conserves energy for long-time integration of oscillatory systems.
	class Leapfrog_StepCalculator : public IODESystemStepCalculator
	{
	public:
		// For a given IODESystem of dimension n, and given initial values for the variables x_start[0..n-1] 
		// and their derivatives dxdt[0..n-1] known at t, uses the Leapfrog method to advance the solution 
		// over an interval h and return the incremented variables as x_out[0..n-1].
		void calcStep(const IODESystem& odeSystem,
									const Real t, const Vector<Real>& x_start, const Vector<Real>& dxdt,
									const Real h, Vector<Real>& x_out, Vector<Real>& x_err_out) const override
		{
			int n = odeSystem.getDim();
			int half_n = n / 2;

			// Split x_start into positions and velocities
			// x_start[0..half_n-1] = positions
			// x_start[half_n..n-1] = velocities

			// 1. Compute velocity at half step: v_{n+1/2} = v_n + (h/2) * a(x_n)
			Vector<Real> x_temp = x_start;
			Vector<Real> dxdt_temp(n);
			odeSystem.derivs(t, x_start, dxdt_temp);

			// positions: x_start[0..half_n-1]
			// velocities: x_start[half_n..n-1]
			// accelerations: dxdt_temp[half_n..n-1]

			Vector<Real> v_half(half_n);
			for (int i = 0; i < half_n; ++i)
				v_half[i] = x_start[half_n + i] + 0.5 * h * dxdt_temp[half_n + i];

			// 2. Update positions: x_{n+1} = x_n + h * v_{n+1/2}
			for (int i = 0; i < half_n; ++i)
				x_out[i] = x_start[i] + h * v_half[i];

			// 3. Compute new acceleration at x_{n+1}
			Vector<Real> x_next = x_out;
			for (int i = 0; i < half_n; ++i)
				x_next[half_n + i] = v_half[i]; // temporary velocities for acceleration calculation

			Vector<Real> dxdt_next(n);
			odeSystem.derivs(t + h, x_next, dxdt_next);

			// 4. Update velocities: v_{n+1} = v_{n+1/2} + (h/2) * a(x_{n+1})
			for (int i = 0; i < half_n; ++i)
				x_out[half_n + i] = v_half[i] + 0.5 * h * dxdt_next[half_n + i];

			// No error estimate for Leapfrog
			for (int i = 0; i < n; ++i)
				x_err_out[i] = 0.0;
		}
	};

	// For a given ODESystem, of dimension n, and given initial values for the variables x_start[0..n-1] 
	// and their derivatives dxdt[0..n-1] known at t, uses the Midpoint method to advance the solution 
	// over an interval h and return the incremented variables as xout[0..n-1].
	class Midpoint_StepCalculator : public IODESystemStepCalculator
	{
	public:
		void calcStep(const IODESystem& odeSystem,
									const Real t, const Vector<Real>& x_start, const Vector<Real>& dxdt,
									const Real h, Vector<Real>& x_out, Vector<Real>& x_err_out) const override
		{
			int i, n = odeSystem.getDim();
			Vector<Real> x_mid(n);

			for (i = 0; i < n; i++)
				x_mid[i] = x_start[i] + 0.5 * h * dxdt[i];

			odeSystem.derivs(t + 0.5 * h, x_mid, x_out);
		}
	};

	// For a given ODESystem, of dimension n, and given initial values for the variables x_start[0..n-1] 
	// and their derivatives dxdt[0..n-1] known at t, uses the fourth-order Runge-Kutta method 
	// to advance the solution over an interval h and return the incremented variables as xout[0..n-1].
	class RungeKutta4_StepCalculator : public IODESystemStepCalculator
	{
	public:
		void calcStep(const IODESystem& odeSystem,
			const Real t, const Vector<Real>& x_start, const Vector<Real>& dxdt,
			const Real h, Vector<Real>& x_out, Vector<Real>& x_err_out) const override
		{
			int i, n = odeSystem.getDim();
			Vector<Real> dx_mid(n), dx_temp(n), x_temp(n);

			Real xh, hh, h6;
			hh = h * 0.5;
			h6 = h / 6.0;
			xh = t + hh;

			for (i = 0; i < n; i++)												// First step
				x_temp[i] = x_start[i] + hh * dxdt[i];

			odeSystem.derivs(xh, x_temp, dx_temp);				// Second step

			for (i = 0; i < n; i++)
				x_temp[i] = x_start[i] + hh * dx_temp[i];

			odeSystem.derivs(xh, x_temp, dx_mid);					// Third step

			for (i = 0; i < n; i++) {
				x_temp[i] = x_start[i] + h * dx_mid[i];
				dx_mid[i] += dx_temp[i];
			}

			odeSystem.derivs(t + h, x_temp, dx_temp);			// Fourth step	

			for (i = 0; i < n; i++)
				x_out[i] = x_start[i] + h6 * (dxdt[i] + dx_temp[i] + 2.0 * dx_mid[i]);
		}
	};
	
	// Given initial values for n variables x[0..n-1] and their derivatives dxdt[0..n-1] known at t,
	// uses	the fifth-order Cash-Karp Runge-Kutta method to advance the solution over an interval h
	// and return the incremented variables as x_out[0..n-1].
	// Also returns an estimate of the local truncation error in x_out using the embedded fourth-order method.
	// VERIFIED: Coefficients match Numerical Recipes 2nd ed. rkck() implementation exactly.
	class RK5_CashKarp_Calculator : public IODESystemStepCalculator
	{
	public:

		void calcStep(const IODESystem& sys,
									Real t, const Vector<Real>& x, const Vector<Real>& dxdt,
									Real h, Vector<Real>& x_out, Vector<Real>& x_err) const override
		{
			static const Real a2 = 0.2, a3 = 0.3, a4 = 0.6, a5 = 1.0, a6 = 0.875,
				b21 = 0.2, b31 = 3.0 / 40.0, b32 = 9.0 / 40.0, b41 = 0.3, b42 = -0.9,
				b43 = 1.2, b51 = -11.0 / 54.0, b52 = 2.5, b53 = -70.0 / 27.0,
				b54 = 35.0 / 27.0, b61 = 1631.0 / 55296.0, b62 = 175.0 / 512.0,
				b63 = 575.0 / 13824.0, b64 = 44275.0 / 110592.0, b65 = 253.0 / 4096.0,
				c1 = 37.0 / 378.0, c3 = 250.0 / 621.0, c4 = 125.0 / 594.0, c6 = 512.0 / 1771.0,
				dc1 = c1 - 2825.0 / 27648.0, dc3 = c3 - 18575.0 / 48384.0,
				dc4 = c4 - 13525.0 / 55296.0, dc5 = -277.00 / 14336.0, dc6 = c6 - 0.25;

			int i, n = x.size();
			Vector<Real> ak2(n), ak3(n), ak4(n), ak5(n), ak6(n), xtemp(n);

			for (i = 0; i < n; i++)
				xtemp[i] = x[i] + b21 * h * dxdt[i];

			sys.derivs(t + a2 * h, xtemp, ak2);
			for (i = 0; i < n; i++)
				xtemp[i] = x[i] + h * (b31 * dxdt[i] + b32 * ak2[i]);

			sys.derivs(t + a3 * h, xtemp, ak3);
			for (i = 0; i < n; i++)
				xtemp[i] = x[i] + h * (b41 * dxdt[i] + b42 * ak2[i] + b43 * ak3[i]);

			sys.derivs(t + a4 * h, xtemp, ak4);
			for (i = 0; i < n; i++)
				xtemp[i] = x[i] + h * (b51 * dxdt[i] + b52 * ak2[i] + b53 * ak3[i] + b54 * ak4[i]);

			sys.derivs(t + a5 * h, xtemp, ak5);
			for (i = 0; i < n; i++)
				xtemp[i] = x[i] + h * (b61 * dxdt[i] + b62 * ak2[i] + b63 * ak3[i] + b64 * ak4[i] + b65 * ak5[i]);

			sys.derivs(t + a6 * h, xtemp, ak6);

			for (i = 0; i < n; i++)
				x_out[i] = x[i] + h * (c1 * dxdt[i] + c3 * ak3[i] + c4 * ak4[i] + c6 * ak6[i]);
			for (i = 0; i < n; i++)
				x_err[i] = h * (dc1 * dxdt[i] + dc3 * ak3[i] + dc4 * ak4[i] + dc5 * ak5[i] + dc6 * ak6[i]);
		}
	};

	// VERIFIED: Implementation matches standard Dormand-Prince (4)5 Butcher tableau.
	// Coefficients verified against multiple sources (Wikipedia, Hairer et al.).
	class DormandPrince5_StepCalculator : public IODESystemStepCalculator
	{
	public:
		// For a given ODESystem, of dimension n, and given initial values for the variables x_start[0..n-1] 
		// and their derivatives dxdt[0..n-1] known at t, uses the Dormand-Prince 5th-order Runge-Kutta method 
		// to advance the solution over an interval h and return the incremented variables as xout[0..n-1].
		void calcStep(const IODESystem& odeSystem,
									const Real t, const Vector<Real>& x_start, const Vector<Real>& dxdt,
									const Real h, Vector<Real>& x_out, Vector<Real>& x_err_out) const override
		{
			// Dormand-Prince coefficients (Butcher tableau)
			static const Real a2 = 1.0 / 5.0;
			static const Real a3 = 3.0 / 10.0;
			static const Real a4 = 4.0 / 5.0;
			static const Real a5 = 8.0 / 9.0;
			static const Real a6 = 1.0;
			static const Real a7 = 1.0;

			static const Real b21 = 1.0 / 5.0;

			static const Real b31 = 3.0 / 40.0;
			static const Real b32 = 9.0 / 40.0;

			static const Real b41 = 44.0 / 45.0;
			static const Real b42 = -56.0 / 15.0;
			static const Real b43 = 32.0 / 9.0;

			static const Real b51 = 19372.0 / 6561.0;
			static const Real b52 = -25360.0 / 2187.0;
			static const Real b53 = 64448.0 / 6561.0;
			static const Real b54 = -212.0 / 729.0;

			static const Real b61 = 9017.0 / 3168.0;
			static const Real b62 = -355.0 / 33.0;
			static const Real b63 = 46732.0 / 5247.0;
			static const Real b64 = 49.0 / 176.0;
			static const Real b65 = -5103.0 / 18656.0;

			static const Real b71 = 35.0 / 384.0;
			static const Real b72 = 0.0;
			static const Real b73 = 500.0 / 1113.0;
			static const Real b74 = 125.0 / 192.0;
			static const Real b75 = -2187.0 / 6784.0;
			static const Real b76 = 11.0 / 84.0;

			// 5th order solution weights
			static const Real c1 = 35.0 / 384.0;
			static const Real c2 = 0.0;
			static const Real c3 = 500.0 / 1113.0;
			static const Real c4 = 125.0 / 192.0;
			static const Real c5 = -2187.0 / 6784.0;
			static const Real c6 = 11.0 / 84.0;
			static const Real c7 = 0.0;

			// 4th order solution weights (for error estimate)
			static const Real d1 = 5179.0 / 57600.0;
			static const Real d2 = 0.0;
			static const Real d3 = 7571.0 / 16695.0;
			static const Real d4 = 393.0 / 640.0;
			static const Real d5 = -92097.0 / 339200.0;
			static const Real d6 = 187.0 / 2100.0;
			static const Real d7 = 1.0 / 40.0;

			int n = x_start.size();
			Vector<Real> k1(n), k2(n), k3(n), k4(n), k5(n), k6(n), k7(n), xtemp(n);

			// k1 = f(t, x)
			for (int i = 0; i < n; ++i)
				k1[i] = dxdt[i];

			// k2 = f(t + a2*h, x + h*b21*k1)
			for (int i = 0; i < n; ++i)
				xtemp[i] = x_start[i] + h * b21 * k1[i];
			odeSystem.derivs(t + a2 * h, xtemp, k2);

			// k3 = f(t + a3*h, x + h*(b31*k1 + b32*k2))
			for (int i = 0; i < n; ++i)
				xtemp[i] = x_start[i] + h * (b31 * k1[i] + b32 * k2[i]);
			odeSystem.derivs(t + a3 * h, xtemp, k3);

			// k4 = f(t + a4*h, x + h*(b41*k1 + b42*k2 + b43*k3))
			for (int i = 0; i < n; ++i)
				xtemp[i] = x_start[i] + h * (b41 * k1[i] + b42 * k2[i] + b43 * k3[i]);
			odeSystem.derivs(t + a4 * h, xtemp, k4);

			// k5 = f(t + a5*h, x + h*(b51*k1 + b52*k2 + b53*k3 + b54*k4))
			for (int i = 0; i < n; ++i)
				xtemp[i] = x_start[i] + h * (b51 * k1[i] + b52 * k2[i] + b53 * k3[i] + b54 * k4[i]);
			odeSystem.derivs(t + a5 * h, xtemp, k5);

			// k6 = f(t + a6*h, x + h*(b61*k1 + b62*k2 + b63*k3 + b64*k4 + b65*k5))
			for (int i = 0; i < n; ++i)
				xtemp[i] = x_start[i] + h * (b61 * k1[i] + b62 * k2[i] + b63 * k3[i] + b64 * k4[i] + b65 * k5[i]);
			odeSystem.derivs(t + a6 * h, xtemp, k6);

			// k7 = f(t + a7*h, x + h*(b71*k1 + b72*k2 + b73*k3 + b74*k4 + b75*k5 + b76*k6))
			for (int i = 0; i < n; ++i)
				xtemp[i] = x_start[i] + h * (b71 * k1[i] + b72 * k2[i] + b73 * k3[i] + b74 * k4[i] + b75 * k5[i] + b76 * k6[i]);
			odeSystem.derivs(t + a7 * h, xtemp, k7);

			// 5th order solution
			for (int i = 0; i < n; ++i)
				x_out[i] = x_start[i] + h * (c1 * k1[i] + c2 * k2[i] + c3 * k3[i] + c4 * k4[i] + c5 * k5[i] + c6 * k6[i] + c7 * k7[i]);

			// 4th order solution (for error estimate)
			for (int i = 0; i < n; ++i)
				x_err_out[i] = h * ((c1 - d1) * k1[i] + (c2 - d2) * k2[i] + (c3 - d3) * k3[i] +
					(c4 - d4) * k4[i] + (c5 - d5) * k5[i] + (c6 - d6) * k6[i] + (c7 - d7) * k7[i]);
		}
	};

	// VERIFIED: Implementation matches Dormand-Prince (7)8 Butcher tableau.
	// High-order method with 14 stages for problems requiring tight error control.
	class DormandPrince8_StepCalculator : public IODESystemStepCalculator
	{
	public:
		// For a given ODESystem, of dimension n, and given initial values for the variables x_start[0..n-1] 
		// and their derivatives dxdt[0..n-1] known at t, uses the Dormand-Prince 8th-order Runge-Kutta method 
		// to advance the solution over an interval h and return the incremented variables as xout[0..n-1].
		void calcStep(const IODESystem& odeSystem,
			const Real t, const Vector<Real>& x_start, const Vector<Real>& dxdt,
			const Real h, Vector<Real>& x_out, Vector<Real>& x_err_out) const override
		{
			constexpr int s = 14; // number of stages

			// c: nodes
			static const Real c[s] = {
					0.0,
					1.0 / 18.0,
					1.0 / 12.0,
					1.0 / 8.0,
					5.0 / 16.0,
					3.0 / 8.0,
					59.0 / 400.0,
					93.0 / 200.0,
					5490023248.0 / 9719169821.0,
					13.0 / 20.0,
					1201146811.0 / 1299019798.0,
					1.0,
					1.0
			};

			// a: stage coefficients
			static const Real a[s][s] = {
					{0},
					{1.0 / 18.0},
					{1.0 / 48.0, 1.0 / 16.0},
					{1.0 / 32.0, 0.0, 3.0 / 32.0},
					{5.0 / 16.0, 0.0, -75.0 / 64.0, 75.0 / 64.0},
					{3.0 / 80.0, 0.0, 0.0, 3.0 / 16.0, 3.0 / 20.0},
					{29443841.0 / 614563906.0, 0.0, 0.0, 77736538.0 / 692538347.0, -28693883.0 / 1125000000.0, 23124283.0 / 1800000000.0},
					{16016141.0 / 946692911.0, 0.0, 0.0, 61564180.0 / 158732637.0, 22789713.0 / 633445777.0, 545815736.0 / 2771057229.0, -180193667.0 / 1043307555.0},
					{39632708.0 / 573591083.0, 0.0, 0.0, -433636366.0 / 683701615.0, -421739975.0 / 2616292301.0, 100302831.0 / 723423059.0, 790204164.0 / 839813087.0, 800635310.0 / 3783071287.0},
					{246121993.0 / 1340847787.0, 0.0, 0.0, -37695042795.0 / 15268766246.0, -309121744.0 / 1061227803.0, -12992083.0 / 490766935.0, 6005943493.0 / 2108947869.0, 393006217.0 / 1396673457.0, 123872331.0 / 1001029789.0},
					{-1028468189.0 / 846180014.0, 0.0, 0.0, 8478235783.0 / 508512852.0, 1311729495.0 / 1432422823.0, -10304129995.0 / 1701304382.0, -48777925059.0 / 3047939560.0, 15336726248.0 / 1032824649.0, -45442868181.0 / 3398467696.0, 3065993473.0 / 597172653.0},
					{185892177.0 / 718116043.0, 0.0, 0.0, -3185094517.0 / 667107341.0, -477755414.0 / 1098053517.0, -703635378.0 / 230739211.0, 5731566787.0 / 1027545527.0, 5232866602.0 / 850066563.0, -4093664535.0 / 808688257.0, 3962137247.0 / 1805957418.0, 65686358.0 / 487910083.0},
					{403863854.0 / 491063109.0, 0.0, 0.0, -5068492393.0 / 434740067.0, -411421997.0 / 543043805.0, 652783627.0 / 914296604.0, 11173962825.0 / 925320556.0, -13158990841.0 / 6184727034.0, 3936647629.0 / 1978049680.0, -160528059.0 / 685178525.0, 248638103.0 / 1413531060.0, 0.0},
					{14005451.0 / 335480064.0, 0.0, 0.0, 0.0, 0.0, -59238493.0 / 1068277825.0, 181606767.0 / 758867731.0, 561292985.0 / 797845732.0, -1041891430.0 / 1371343529.0, 760417239.0 / 1151165299.0, 118820643.0 / 751138087.0, -528747749.0 / 2220607170.0, 1.0 / 4.0}
			};

			// b8: 8th order weights, b7: 7th order weights (for error estimate)
			static const Real b8[s] = {
					14005451.0 / 335480064.0,
					0.0,
					0.0,
					0.0,
					0.0,
					-59238493.0 / 1068277825.0,
					181606767.0 / 758867731.0,
					561292985.0 / 797845732.0,
					-1041891430.0 / 1371343529.0,
					760417239.0 / 1151165299.0,
					118820643.0 / 751138087.0,
					-528747749.0 / 2220607170.0,
					1.0 / 4.0
			};
			static const Real b7[s] = {
					13451932.0 / 455176623.0,
					0.0,
					0.0,
					0.0,
					0.0,
					-808719846.0 / 976000145.0,
					1757004468.0 / 5645159321.0,
					656045339.0 / 265891186.0,
					-3867574721.0 / 1518517206.0,
					465885868.0 / 322736535.0,
					53011238.0 / 667516719.0,
					2.0 / 45.0,
					0.0
			};

			int n = x_start.size();
			Vector<Vector<Real>> k(s, Vector<Real>(n));
			Vector<Real> xtemp(n);

			// k1 = dxdt
			for (int i = 0; i < n; ++i)
				k[0][i] = dxdt[i];

			// Compute all stages
			for (int stage = 1; stage < s; ++stage) 
			{
				for (int i = 0; i < n; ++i) {
					Real sum = 0.0;
					for (int j = 0; j < stage; ++j)
						sum += a[stage][j] * k[j][i];
					xtemp[i] = x_start[i] + h * sum;
				}
				odeSystem.derivs(t + c[stage] * h, xtemp, k[stage]);
			}

			// 8th order solution
			for (int i = 0; i < n; ++i) 
			{
				Real sum = 0.0;
				for (int stage = 0; stage < s; ++stage)
					sum += b8[stage] * k[stage][i];
				x_out[i] = x_start[i] + h * sum;
			}

			// Error estimate (8th - 7th order)
			for (int i = 0; i < n; ++i) 
			{
				Real sum = 0.0;
				for (int stage = 0; stage < s; ++stage)
					sum += (b8[stage] - b7[stage]) * k[stage][i];
				x_err_out[i] = h * sum;
			}
		}
	};

	class StepCalculators
	{
	public:
		static inline EulerStep_Calculator					EulerStepCalc;
		static inline EulerCromer_StepCalculator		EulerCromerStepCalc;
		static inline VelocityVerlet_StepCalculator	VelocityVerletStepCalc;
		static inline Leapfrog_StepCalculator				LeapfrogStepCalc;
		static inline Midpoint_StepCalculator				MidpointStepCalc;
		static inline RungeKutta4_StepCalculator		RK4_Basic;
		static inline RK5_CashKarp_Calculator				RK5_CashKarp;
		static inline DormandPrince5_StepCalculator	DormandPrince5StepCalc;
		static inline DormandPrince8_StepCalculator	DormandPrince8StepCalc;
	};
}

#endif // MML_ODE_SYSTEM_STEP_CALCULATORS_H

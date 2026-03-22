///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        ODEStepCalculators.h                                          ///
///  Description: ODE step calculators (Euler, RK4, RK5, Cash-Karp, Dormand-Prince)   ///
///               Single-step methods for ODE integration                             ///
///                                                                                   ///
///  REFERENCES:                                                                      ///
///    [NR3]  Press et al., Numerical Recipes 3rd ed., Ch. 17                         ///
///    [HNW1] Hairer et al., Solving ODEs I, Ch. II                                   ///
///    [DP80] Dormand & Prince (1980), J. Comp. Appl. Math. 6(1), pp. 19-26          ///
///    [CK90] Cash & Karp (1990), ACM TOMS 16(3), pp. 201-222                        ///
///                                                                                   ///
///  See references/book_references.md and references/paperes_references.md          ///
///                                                                                   ///
///  COMPLEXITY SUMMARY (per step, N = system dimension)                              ///
///  ====================================================                              ///
///    Stepper              Order  Stages  f-evals/step  Memory                       ///
///    -------              -----  ------  ------------  ------                       ///
///    Euler                  1      1          1         O(N)                         ///
///    Euler-Cromer           2      1          2         O(N)                         ///
///    Velocity Verlet        2      1          2         O(N)   (symplectic)          ///
///    Leapfrog               2      1          2         O(N)   (symplectic)          ///
///    Midpoint               2      1          2         O(N)                         ///
///    RK4 (classic)          4      4          4         O(N)                         ///
///    RK5 Cash-Karp        5(4)     6          6         O(N)   (with error est.)     ///
///    Dormand-Prince 5     5(4)     7          7         O(N)   (FSAL, with error)    ///
///    Dormand-Prince 8     8(7)    13         13         O(N)   (with error est.)     ///
///                                                                                   ///
///  Each f-eval costs O(N) for the function evaluation itself, so total per-step    ///
///  cost is O(stages × N). Choose stepper based on accuracy needs vs cost.          ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_ODE_STEP_CALCULATORS_H
#define MML_ODE_STEP_CALCULATORS_H

#include "mml/MMLBase.h"

#include "mml/core/AlgorithmTypes.h"
#include "mml/interfaces/IODESystem.h"
#include "mml/interfaces/IODESystemStepCalculator.h"
#include "ODERKCoefficients.h"

#include "mml/base/ODESystem.h"

// NOTE: These calculators are stateless single-step implementations used by
// `ODESystemFixedStepSolver` for fixed-step integration.
// RK coefficients are centralized in ODERKCoefficients.h for consistency
// with the adaptive steppers in ODESteppers.h.

namespace MML {
	// For a given IODESystem of dimension n, and given initial values for the variables x_start[0..n-1]
	// and their derivatives dxdt[0..n-1] known at t, uses the Euler method to advance the solution
	// over an interval h and return the incremented variables as x_out[0..n-1].
	// Complexity: O(N) per step, 1 f-eval (provided via dxdt). No error estimate.
	class EulerStep_Calculator : public IODESystemStepCalculator {
	public:
		void calcStep(const IODESystem& odeSystem, const Real t, const Vector<Real>& x_start, const Vector<Real>& dxdt, const Real h,
					  Vector<Real>& x_out, Vector<Real>& x_err_out) const override {
			int i, n = odeSystem.getDim();

			for (i = 0; i < n; i++)
				x_out[i] = x_start[i] + h * dxdt[i];

			// No error estimate
			for (i = 0; i < n; i++)
				x_err_out[i] = 0.0;
		}
	};

	// For a given IODESystem of dimension n, and given initial values for the variables x_start[0..n-1]
	// and their derivatives dxdt[0..n-1] known at t, uses the Euler-Cromer method to advance the solution
	// over an interval h and return the incremented variables as x_out[0..n-1].
	// Complexity: O(N) per step, 2 f-evals. Provides local error estimate.
	class EulerCromer_StepCalculator : public IODESystemStepCalculator {
	public:
		void calcStep(const IODESystem& odeSystem, const Real t, const Vector<Real>& x_start, const Vector<Real>& dxdt, const Real h,
					        Vector<Real>& x_out, Vector<Real>& x_err_out) const override 
{
			int i, n = odeSystem.getDim();
			Vector<Real> x_new(n);
			Vector<Real> dxdt_end(n);

			for (i = 0; i < n; i++)
				x_new[i] = x_start[i] + h * dxdt[i];

			// Use derivative at the end of the step to form a next-state update.
			// This keeps x_out as the state at (t+h) and provides a simple local error estimate.
			odeSystem.derivs(t + h, x_new, dxdt_end);
			for (i = 0; i < n; i++)
				x_out[i] = x_start[i] + 0.5 * h * (dxdt[i] + dxdt_end[i]);

			// Error estimate: difference between improved-Euler (2nd order) and Euler (1st order)
			for (i = 0; i < n; i++)
				x_err_out[i] = x_out[i] - x_new[i];
		}
	};

	// For a given IODESystem of dimension n, and given initial values for the variables x_start[0..n-1]
	// and their derivatives dxdt[0..n-1] known at t, uses the Velocity Verlet method to advance the solution
	// over an interval h and return the incremented variables as x_out[0..n-1].
	// Complexity: O(N) per step, 2 f-evals. Symplectic (conserves energy for Hamiltonian systems).
	class VelocityVerlet_StepCalculator : public IODESystemStepCalculator {
	public:
		void calcStep(const IODESystem& odeSystem, const Real t, const Vector<Real>& x_start, const Vector<Real>& dxdt, const Real h,
					  Vector<Real>& x_out, Vector<Real>& x_err_out) const override {
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
	// Complexity: O(N) per step, 2 f-evals. Symplectic (time-reversible, energy-conserving).
	class Leapfrog_StepCalculator : public IODESystemStepCalculator {
	public:
		// For a given IODESystem of dimension n, and given initial values for the variables x_start[0..n-1]
		// and their derivatives dxdt[0..n-1] known at t, uses the Leapfrog method to advance the solution
		// over an interval h and return the incremented variables as x_out[0..n-1].
		void calcStep(const IODESystem& odeSystem, const Real t, const Vector<Real>& x_start, const Vector<Real>& dxdt, const Real h,
					  Vector<Real>& x_out, Vector<Real>& x_err_out) const override {
			int n = odeSystem.getDim();
			int half_n = n / 2;

			// Split x_start into positions and velocities
			// x_start[0..half_n-1] = positions
			// x_start[half_n..n-1] = velocities
			// dxdt[half_n..n-1] = accelerations (passed in, computed by caller)

			// 1. Compute velocity at half step: v_{n+1/2} = v_n + (h/2) * a(x_n)
			// Uses the passed-in dxdt (already computed by the solver)
			Vector<Real> v_half(half_n);
			for (int i = 0; i < half_n; ++i)
				v_half[i] = x_start[half_n + i] + 0.5 * h * dxdt[half_n + i];

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
	// Complexity: O(N) per step, 2 f-evals. No error estimate.
	class Midpoint_StepCalculator : public IODESystemStepCalculator {
	public:
		void calcStep(const IODESystem& odeSystem, const Real t, const Vector<Real>& x_start, const Vector<Real>& dxdt, const Real h,
					  Vector<Real>& x_out, Vector<Real>& x_err_out) const override {
			int i, n = odeSystem.getDim();
			Vector<Real> x_mid(n);
			Vector<Real> dx_mid(n);

			for (i = 0; i < n; i++)
				x_mid[i] = x_start[i] + 0.5 * h * dxdt[i];

			odeSystem.derivs(t + 0.5 * h, x_mid, dx_mid);
			
      for (i = 0; i < n; i++)
				x_out[i] = x_start[i] + h * dx_mid[i];

			// No error estimate
			for (i = 0; i < n; i++)
				x_err_out[i] = 0.0;
		}
	};

	// For a given ODESystem, of dimension n, and given initial values for the variables x_start[0..n-1]
	// and their derivatives dxdt[0..n-1] known at t, uses the fourth-order Runge-Kutta method
	// to advance the solution over an interval h and return the incremented variables as xout[0..n-1].
	// Complexity: O(N) per step, 4 f-evals. No error estimate.
	class RungeKutta4_StepCalculator : public IODESystemStepCalculator {
	public:
		void calcStep(const IODESystem& odeSystem, const Real t, const Vector<Real>& x_start, const Vector<Real>& dxdt, const Real h,
					  Vector<Real>& x_out, Vector<Real>& x_err_out) const override {
			int i, n = odeSystem.getDim();
			Vector<Real> dx_mid(n), dx_temp(n), x_temp(n);

			Real xh, hh, h6;
			hh = h * 0.5;
			h6 = h / 6.0;
			xh = t + hh;

			for (i = 0; i < n; i++) // First step
				x_temp[i] = x_start[i] + hh * dxdt[i];

			odeSystem.derivs(xh, x_temp, dx_temp); // Second step

			for (i = 0; i < n; i++)
				x_temp[i] = x_start[i] + hh * dx_temp[i];

			odeSystem.derivs(xh, x_temp, dx_mid); // Third step

			for (i = 0; i < n; i++) {
				x_temp[i] = x_start[i] + h * dx_mid[i];
				dx_mid[i] += dx_temp[i];
			}

			odeSystem.derivs(t + h, x_temp, dx_temp); // Fourth step

			for (i = 0; i < n; i++)
				x_out[i] = x_start[i] + h6 * (dxdt[i] + dx_temp[i] + 2.0 * dx_mid[i]);

			// No error estimate
			for (i = 0; i < n; i++)
				x_err_out[i] = 0.0;
		}
	};

	/******************************************************************************
	 * CASH-KARP RK5(4) EMBEDDED METHOD
	 *
	 * Fifth-order Runge-Kutta with embedded fourth-order error estimate.
	 * 6 stages, FSAL property not used.
	 * Complexity: O(N) per step, 6 f-evals. Provides embedded error estimate.
	 *
	 * REFERENCES:
	 * - [CK90] Cash, J.R., & Karp, A.H. (1990). A variable order Runge-Kutta
	 *          method for initial value problems with rapidly varying right-hand
	 *          sides. ACM Trans. Math. Softw. 16(3), pp. 201-222.
	 * - [NR3]  Press et al., Numerical Recipes 3rd ed., Section 17.2
	 *
	 * VERIFIED: Coefficients match Numerical Recipes 2nd ed. rkck() exactly.
	 ******************************************************************************/
	class RK5_CashKarp_Calculator : public IODESystemStepCalculator {
	public:
		void calcStep(const IODESystem& sys, Real t, const Vector<Real>& x, const Vector<Real>& dxdt, Real h, Vector<Real>& x_out,
					  Vector<Real>& x_err) const override {
			using CK = RKCoeff::CashKarp5;

			int i, n = x.size();
			Vector<Real> ak2(n), ak3(n), ak4(n), ak5(n), ak6(n), xtemp(n);

			for (i = 0; i < n; i++)
				xtemp[i] = x[i] + CK::a21 * h * dxdt[i];

			sys.derivs(t + CK::c2 * h, xtemp, ak2);
			for (i = 0; i < n; i++)
				xtemp[i] = x[i] + h * (CK::a31 * dxdt[i] + CK::a32 * ak2[i]);

			sys.derivs(t + CK::c3 * h, xtemp, ak3);
			for (i = 0; i < n; i++)
				xtemp[i] = x[i] + h * (CK::a41 * dxdt[i] + CK::a42 * ak2[i] + CK::a43 * ak3[i]);

			sys.derivs(t + CK::c4 * h, xtemp, ak4);
			for (i = 0; i < n; i++)
				xtemp[i] = x[i] + h * (CK::a51 * dxdt[i] + CK::a52 * ak2[i] + CK::a53 * ak3[i] + CK::a54 * ak4[i]);

			sys.derivs(t + CK::c5 * h, xtemp, ak5);
			for (i = 0; i < n; i++)
				xtemp[i] = x[i] + h * (CK::a61 * dxdt[i] + CK::a62 * ak2[i] + CK::a63 * ak3[i] + CK::a64 * ak4[i] + CK::a65 * ak5[i]);

			sys.derivs(t + CK::c6 * h, xtemp, ak6);

			for (i = 0; i < n; i++)
				x_out[i] = x[i] + h * (CK::b1 * dxdt[i] + CK::b3 * ak3[i] + CK::b4 * ak4[i] + CK::b6 * ak6[i]);
			for (i = 0; i < n; i++)
				x_err[i] = h * (CK::e1 * dxdt[i] + CK::e3 * ak3[i] + CK::e4 * ak4[i] + CK::e5 * ak5[i] + CK::e6 * ak6[i]);
		}
	};

	/******************************************************************************
	 * DORMAND-PRINCE 5(4) METHOD
	 *
	 * Fifth-order Runge-Kutta with embedded fourth-order error estimate.
	 * 7 stages with FSAL (First Same As Last) property.
	 * The standard method in MATLAB's ode45.
	 * Complexity: O(N) per step, 7 f-evals (6 effective with FSAL).
	 *             Provides embedded error estimate.
	 *
	 * REFERENCES:
	 * - [DP80] Dormand, J.R., & Prince, P.J. (1980). A family of embedded
	 *          Runge-Kutta formulae. J. Comp. Appl. Math. 6(1), pp. 19-26.
	 * - [HNW1] Hairer et al., Solving ODEs I, Chapter II.5
	 * - [NR3]  Press et al., Numerical Recipes 3rd ed., Section 17.2
	 *
	 * VERIFIED: Coefficients match standard Butcher tableau.
	 ******************************************************************************/
	class DormandPrince5_StepCalculator : public IODESystemStepCalculator {
	public:
		// For a given ODESystem, of dimension n, and given initial values for the variables x_start[0..n-1]
		// and their derivatives dxdt[0..n-1] known at t, uses the Dormand-Prince 5th-order Runge-Kutta method
		// to advance the solution over an interval h and return the incremented variables as xout[0..n-1].
		void calcStep(const IODESystem& odeSystem, const Real t, const Vector<Real>& x_start, const Vector<Real>& dxdt, const Real h,
					  Vector<Real>& x_out, Vector<Real>& x_err_out) const override {
			using DP = RKCoeff::DormandPrince5;

			int n = x_start.size();
			Vector<Real> k1(n), k2(n), k3(n), k4(n), k5(n), k6(n), k7(n), xtemp(n);

			// k1 = f(t, x)
			for (int i = 0; i < n; ++i)
				k1[i] = dxdt[i];

			// k2 = f(t + c2*h, x + h*a21*k1)
			for (int i = 0; i < n; ++i)
				xtemp[i] = x_start[i] + h * DP::a21 * k1[i];
			odeSystem.derivs(t + DP::c2 * h, xtemp, k2);

			// k3 = f(t + c3*h, x + h*(a31*k1 + a32*k2))
			for (int i = 0; i < n; ++i)
				xtemp[i] = x_start[i] + h * (DP::a31 * k1[i] + DP::a32 * k2[i]);
			odeSystem.derivs(t + DP::c3 * h, xtemp, k3);

			// k4 = f(t + c4*h, x + h*(a41*k1 + a42*k2 + a43*k3))
			for (int i = 0; i < n; ++i)
				xtemp[i] = x_start[i] + h * (DP::a41 * k1[i] + DP::a42 * k2[i] + DP::a43 * k3[i]);
			odeSystem.derivs(t + DP::c4 * h, xtemp, k4);

			// k5 = f(t + c5*h, x + h*(a51*k1 + a52*k2 + a53*k3 + a54*k4))
			for (int i = 0; i < n; ++i)
				xtemp[i] = x_start[i] + h * (DP::a51 * k1[i] + DP::a52 * k2[i] + DP::a53 * k3[i] + DP::a54 * k4[i]);
			odeSystem.derivs(t + DP::c5 * h, xtemp, k5);

			// k6 = f(t + c6*h, x + h*(a61*k1 + a62*k2 + a63*k3 + a64*k4 + a65*k5))
			for (int i = 0; i < n; ++i)
				xtemp[i] = x_start[i] + h * (DP::a61 * k1[i] + DP::a62 * k2[i] + DP::a63 * k3[i] + DP::a64 * k4[i] + DP::a65 * k5[i]);
			odeSystem.derivs(t + DP::c6 * h, xtemp, k6);

			// k7 = f(t + c7*h, x + h*(a71*k1 + a72*k2 + a73*k3 + a74*k4 + a75*k5 + a76*k6))
			for (int i = 0; i < n; ++i)
				xtemp[i] = x_start[i] + h * (DP::a71 * k1[i] + DP::a72 * k2[i] + DP::a73 * k3[i] + DP::a74 * k4[i] + DP::a75 * k5[i] + DP::a76 * k6[i]);
			odeSystem.derivs(t + DP::c7 * h, xtemp, k7);

			// 5th order solution
			for (int i = 0; i < n; ++i)
				x_out[i] = x_start[i] + h * (DP::b1 * k1[i] + DP::b2 * k2[i] + DP::b3 * k3[i] + DP::b4 * k4[i] + DP::b5 * k5[i] + DP::b6 * k6[i] + DP::b7 * k7[i]);

			// Error estimate (5th - 4th order)
			for (int i = 0; i < n; ++i)
				x_err_out[i] = h * (DP::e1 * k1[i] + DP::e2 * k2[i] + DP::e3 * k3[i] + DP::e4 * k4[i] + DP::e5 * k5[i] + DP::e6 * k6[i] + DP::e7 * k7[i]);
		}
	};

	/******************************************************************************
	 * DORMAND-PRINCE 8(7) METHOD
	 *
	 * Eighth-order Runge-Kutta with embedded seventh-order error estimate.
	 * 13 stages. For high-precision problems requiring tight error control.
	 * Complexity: O(N) per step, 13 f-evals. Provides embedded error estimate.
	 *             13× more expensive per step than Euler — use only when high
	 *             order allows much larger step sizes to compensate.
	 *
	 * REFERENCES:
	 * - [PD81] Prince, P.J., & Dormand, J.R. (1981). High order embedded
	 *          Runge-Kutta formulae. J. Comp. Appl. Math. 7(1), pp. 67-75.
	 * - [HNW1] Hairer et al., Solving ODEs I, Chapter II.5
	 *
	 * VERIFIED: Coefficients match standard DOP853 Butcher tableau.
	 ******************************************************************************/
	class DormandPrince8_StepCalculator : public IODESystemStepCalculator {
	public:
		// For a given ODESystem, of dimension n, and given initial values for the variables x_start[0..n-1]
		// and their derivatives dxdt[0..n-1] known at t, uses the Dormand-Prince 8th-order Runge-Kutta method
		// to advance the solution over an interval h and return the incremented variables as xout[0..n-1].
		void calcStep(const IODESystem& odeSystem, const Real t, const Vector<Real>& x_start, const Vector<Real>& dxdt, const Real h,
					  Vector<Real>& x_out, Vector<Real>& x_err_out) const override {
			using DP8 = RKCoeff::DormandPrince8;
			constexpr int s = DP8::stages;

			int n = x_start.size();
			Vector<Vector<Real>> k(s, Vector<Real>(n));
			Vector<Real> xtemp(n);

			// k[0] = dxdt
			for (int i = 0; i < n; ++i)
				k[0][i] = dxdt[i];

			// Compute stages 2-13 using centralized coefficients
			// Stage 2
			for (int i = 0; i < n; ++i)
				xtemp[i] = x_start[i] + h * DP8::a2[0] * k[0][i];
			odeSystem.derivs(t + DP8::c[1] * h, xtemp, k[1]);

			// Stage 3
			for (int i = 0; i < n; ++i)
				xtemp[i] = x_start[i] + h * (DP8::a3[0] * k[0][i] + DP8::a3[1] * k[1][i]);
			odeSystem.derivs(t + DP8::c[2] * h, xtemp, k[2]);

			// Stage 4
			for (int i = 0; i < n; ++i)
				xtemp[i] = x_start[i] + h * (DP8::a4[0] * k[0][i] + DP8::a4[1] * k[1][i] + DP8::a4[2] * k[2][i]);
			odeSystem.derivs(t + DP8::c[3] * h, xtemp, k[3]);

			// Stage 5
			for (int i = 0; i < n; ++i)
				xtemp[i] = x_start[i] + h * (DP8::a5[0] * k[0][i] + DP8::a5[1] * k[1][i] + DP8::a5[2] * k[2][i] + DP8::a5[3] * k[3][i]);
			odeSystem.derivs(t + DP8::c[4] * h, xtemp, k[4]);

			// Stage 6
			for (int i = 0; i < n; ++i)
				xtemp[i] = x_start[i] + h * (DP8::a6[0] * k[0][i] + DP8::a6[1] * k[1][i] + DP8::a6[2] * k[2][i] + DP8::a6[3] * k[3][i] + DP8::a6[4] * k[4][i]);
			odeSystem.derivs(t + DP8::c[5] * h, xtemp, k[5]);

			// Stage 7
			for (int i = 0; i < n; ++i)
				xtemp[i] = x_start[i] + h * (DP8::a7[0] * k[0][i] + DP8::a7[1] * k[1][i] + DP8::a7[2] * k[2][i] + DP8::a7[3] * k[3][i] + DP8::a7[4] * k[4][i] + DP8::a7[5] * k[5][i]);
			odeSystem.derivs(t + DP8::c[6] * h, xtemp, k[6]);

			// Stage 8
			for (int i = 0; i < n; ++i)
				xtemp[i] = x_start[i] + h * (DP8::a8[0] * k[0][i] + DP8::a8[1] * k[1][i] + DP8::a8[2] * k[2][i] + DP8::a8[3] * k[3][i] + DP8::a8[4] * k[4][i] + DP8::a8[5] * k[5][i] + DP8::a8[6] * k[6][i]);
			odeSystem.derivs(t + DP8::c[7] * h, xtemp, k[7]);

			// Stage 9
			for (int i = 0; i < n; ++i)
				xtemp[i] = x_start[i] + h * (DP8::a9[0] * k[0][i] + DP8::a9[1] * k[1][i] + DP8::a9[2] * k[2][i] + DP8::a9[3] * k[3][i] + DP8::a9[4] * k[4][i] + DP8::a9[5] * k[5][i] + DP8::a9[6] * k[6][i] + DP8::a9[7] * k[7][i]);
			odeSystem.derivs(t + DP8::c[8] * h, xtemp, k[8]);

			// Stage 10
			for (int i = 0; i < n; ++i)
				xtemp[i] = x_start[i] + h * (DP8::a10[0] * k[0][i] + DP8::a10[1] * k[1][i] + DP8::a10[2] * k[2][i] + DP8::a10[3] * k[3][i] + DP8::a10[4] * k[4][i] + DP8::a10[5] * k[5][i] + DP8::a10[6] * k[6][i] + DP8::a10[7] * k[7][i] + DP8::a10[8] * k[8][i]);
			odeSystem.derivs(t + DP8::c[9] * h, xtemp, k[9]);

			// Stage 11
			for (int i = 0; i < n; ++i)
				xtemp[i] = x_start[i] + h * (DP8::a11[0] * k[0][i] + DP8::a11[1] * k[1][i] + DP8::a11[2] * k[2][i] + DP8::a11[3] * k[3][i] + DP8::a11[4] * k[4][i] + DP8::a11[5] * k[5][i] + DP8::a11[6] * k[6][i] + DP8::a11[7] * k[7][i] + DP8::a11[8] * k[8][i] + DP8::a11[9] * k[9][i]);
			odeSystem.derivs(t + DP8::c[10] * h, xtemp, k[10]);

			// Stage 12
			for (int i = 0; i < n; ++i)
				xtemp[i] = x_start[i] + h * (DP8::a12[0] * k[0][i] + DP8::a12[1] * k[1][i] + DP8::a12[2] * k[2][i] + DP8::a12[3] * k[3][i] + DP8::a12[4] * k[4][i] + DP8::a12[5] * k[5][i] + DP8::a12[6] * k[6][i] + DP8::a12[7] * k[7][i] + DP8::a12[8] * k[8][i] + DP8::a12[9] * k[9][i] + DP8::a12[10] * k[10][i]);
			odeSystem.derivs(t + DP8::c[11] * h, xtemp, k[11]);

			// Stage 13
			for (int i = 0; i < n; ++i)
				xtemp[i] = x_start[i] + h * (DP8::a13[0] * k[0][i] + DP8::a13[1] * k[1][i] + DP8::a13[2] * k[2][i] + DP8::a13[3] * k[3][i] + DP8::a13[4] * k[4][i] + DP8::a13[5] * k[5][i] + DP8::a13[6] * k[6][i] + DP8::a13[7] * k[7][i] + DP8::a13[8] * k[8][i] + DP8::a13[9] * k[9][i] + DP8::a13[10] * k[10][i] + DP8::a13[11] * k[11][i]);
			odeSystem.derivs(t + DP8::c[12] * h, xtemp, k[12]);

			// 8th order solution
			for (int i = 0; i < n; ++i) {
				Real sum = 0.0;
				for (int stage = 0; stage < s; ++stage)
					sum += DP8::b8[stage] * k[stage][i];
				x_out[i] = x_start[i] + h * sum;
			}

			// Error estimate (8th - 7th order)
			for (int i = 0; i < n; ++i) {
				Real sum = 0.0;
				for (int stage = 0; stage < s; ++stage)
					sum += (DP8::b8[stage] - DP8::b7[stage]) * k[stage][i];
				x_err_out[i] = h * sum;
			}
		}
	};

	class StepCalculators {
	public:
		static inline EulerStep_Calculator EulerStepCalc;
		static inline EulerCromer_StepCalculator EulerCromerStepCalc;
		static inline VelocityVerlet_StepCalculator VelocityVerletStepCalc;
		static inline Leapfrog_StepCalculator LeapfrogStepCalc;
		static inline Midpoint_StepCalculator MidpointStepCalc;
		static inline RungeKutta4_StepCalculator RK4_Basic;
		static inline RK5_CashKarp_Calculator RK5_CashKarp;
		static inline DormandPrince5_StepCalculator DormandPrince5StepCalc;
		static inline DormandPrince8_StepCalculator DormandPrince8StepCalc;
	};
} // namespace MML

#endif // MML_ODE_STEP_CALCULATORS_H

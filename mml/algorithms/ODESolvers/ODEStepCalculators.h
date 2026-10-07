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
///    Leapfrog               2      1          2         O(N)   (compatibility name)  ///
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

#include <mml/MMLBase.h>

#include <mml/base/AlgorithmTypes.h>
#include <mml/interfaces/IODESystem.h>
#include <mml/interfaces/IODESystemStepCalculator.h>
#include "ODERKCoefficients.h"
#include "ODEStepperInfrastructure.h"

#include <mml/base/ODESystem.h>

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

	// Velocity Verlet, equivalently synchronized kick-drift-kick Leapfrog, for systems
	// q' = v, v' = a(t, q). State layout must be [positions..., velocities...] and
	// the dimension must be even. Velocity-dependent acceleration is not supported.
	// Complexity: O(N) per step, 2 f-evals including the derivative supplied by the caller.
	// Second-order and symplectic for separable Hamiltonian systems with constant step size.
	class VelocityVerlet_StepCalculator : public IODESystemStepCalculator {
	public:
		void calcStep(const IODESystem& odeSystem, const Real t, const Vector<Real>& x_start, const Vector<Real>& dxdt, const Real h,
					  Vector<Real>& x_out, Vector<Real>& x_err_out) const override {
			int n = odeSystem.getDim();
			if (n % 2 != 0)
				throw ODESolverError("Velocity Verlet requires an even-dimensional [positions, velocities] state");

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

	// Compatibility name for synchronized kick-drift-kick Velocity Verlet.
	// This is not a staggered-state Leapfrog API; it returns q and v at the same time.
	class Leapfrog_StepCalculator : public VelocityVerlet_StepCalculator {};

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
			// Fixed-step calculators are stateless, so their stage storage is local to the call.
			std::vector<Vector<Real>> stages(CK::stages, Vector<Real>(x.size()));
			Vector<Real> workspace(x.size());
			ExplicitRKStageEvaluator<CK>::evaluate(sys, t, x, h, dxdt, stages, workspace);
			ExplicitRKStageEvaluator<CK>::combineSolution(x, h, stages, x_out);
			ExplicitRKStageEvaluator<CK>::combineError(h, stages, x_err);
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
			std::vector<Vector<Real>> stages(DP::stages, Vector<Real>(x_start.size()));
			Vector<Real> workspace(x_start.size());
			ExplicitRKStageEvaluator<DP>::evaluate(odeSystem, t, x_start, h, dxdt, stages, workspace);
			ExplicitRKStageEvaluator<DP>::combineSolution(x_start, h, stages, x_out);
			ExplicitRKStageEvaluator<DP>::combineError(h, stages, x_err_out);
		}
	};

	/******************************************************************************
	 * DOP853 8(5,3) METHOD
	 *
	 * Eighth-order Runge-Kutta with embedded fifth- and third-order estimators.
	 * 12 propagation stages. For high-precision problems requiring tight error
	 * control. Complexity: O(N) per step, 12 f-evals when the initial derivative
	 * is supplied. Use only when high
	 *             order allows much larger step sizes to compensate.
	 *
	 * REFERENCES:
	 * - [HNW1] Hairer, Norsett & Wanner, Solving ODEs I, Section II.5
	 *
	 * VERIFIED: Coefficients match Hairer's DOP853 implementation.
	 ******************************************************************************/
	class DormandPrince8_StepCalculator : public IODESystemStepCalculator {
	public:
		// For a given ODESystem, of dimension n, and given initial values for the variables x_start[0..n-1]
		// and their derivatives dxdt[0..n-1] known at t, uses the Dormand-Prince 8th-order Runge-Kutta method
		// to advance the solution over an interval h and return the incremented variables as xout[0..n-1].
		void calcStep(const IODESystem& odeSystem, const Real t, const Vector<Real>& x_start, const Vector<Real>& dxdt, const Real h,
					  Vector<Real>& x_out, Vector<Real>& x_err_out) const override {
			using DP8 = RKCoeff::DormandPrince8;
			std::vector<Vector<Real>> stages(DP8::stages, Vector<Real>(x_start.size()));
			Vector<Real> workspace(x_start.size());
			ExplicitRKStageEvaluator<DP8>::evaluate(odeSystem, t, x_start, h, dxdt, stages, workspace);
			ExplicitRKStageEvaluator<DP8>::combineSolution(x_start, h, stages, x_out);
			ExplicitRKStageEvaluator<DP8>::combineError(h, stages, x_err_out);
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

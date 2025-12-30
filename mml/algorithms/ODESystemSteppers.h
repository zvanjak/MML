///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        ODESystemSteppers.h                                                 ///
///  Description: ODE steppers with error estimation and adaptive control             ///
///                                                                                   ///
///  NOTE:        This file provides backward compatibility.                          ///
///               The implementation has been moved to ODEAdaptiveIntegrator.h        ///
///               Prefer using the new adaptive steppers directly:                    ///
///               - DormandPrince5_Stepper (FSAL, dense output)                       ///
///               - CashKarp_Stepper (classic 5(4) embedded pair)                     ///
///               - DormandPrince8_Stepper (high-order 8(7) method)                   ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_ODE_SYSTEM_STEPPERS_H
#define MML_ODE_SYSTEM_STEPPERS_H

#include "MMLBase.h"

#include "interfaces/IODESystem.h"
#include "interfaces/IODESystemStepCalculator.h"

#include "base/ODESystem.h"

#include "algorithms/ODESystemStepCalculators.h"
#include "algorithms/ODEAdaptiveIntegrator.h"

namespace MML
{
	/// @brief Legacy base stepper class
	/// @deprecated Use IAdaptiveStepper from ODEAdaptiveIntegrator.h instead
	class StepperBase 
	{
	protected:
		const IODESystem& _sys;

		Real& _t;
		Vector<Real>& _x;
		Vector<Real>& _dxdt;

		Real _absTol, _relTol;
		Real _eps;

		Real _tOld;
		Real _hDone, _hNext;

		Vector<Real> _xout, _xerr;

	public:
		StepperBase(const IODESystem& sys, Real& t, Vector<Real>& x, Vector<Real>& dxdt)
			: _sys(sys), _t(t), _x(x), _dxdt(dxdt) {
		}

		virtual void doStep(Real htry, Real eps) = 0;

		Real hDone() const { return _hDone; }
		Real& hDone() { return _hDone; }

		Real hNext() const { return _hNext; }
		Real& hNext() { return _hNext; }
	};

	/// @brief Legacy Cash-Karp stepper (reference-based interface)
	/// @deprecated Use CashKarp_Stepper with ODEAdaptiveIntegrator instead
	class RK5_CashKarp_Stepper : public StepperBase
	{
	private:
		RK5_CashKarp_Calculator _stepCalc;

	public:
		RK5_CashKarp_Stepper(const IODESystem& sys, Real& t, Vector<Real>& x, Vector<Real>& dxdt)
			: StepperBase(sys, t, x, dxdt) {}

		void doStep(Real htry, Real eps) override
		{
			const Real SAFETY = 0.9, PGROW = -0.2, PSHRNK = -0.25, ERRCON = 1.89e-4;
			Real errmax, h, htemp;
			int n = _sys.getDim();
			Vector<Real> xerr(n), xscale(n), xtemp(n);

			h = htry;
			
			for (int i = 0; i < n; i++)
				xscale[i] = std::abs(_x[i]) + std::abs(_dxdt[i] * h) + 1e-25;

			for (;;) {
				_stepCalc.calcStep(_sys, _t, _x, _dxdt, h, xtemp, xerr);

				for (int i = 0; i < n; i++) {
					if (std::isnan(xtemp[i]) || std::isinf(xtemp[i]))
						throw ODESolverError("Non-finite value in RK5 step calculation");
				}

				errmax = 0.0;
				for (int i = 0; i < n; i++)
					errmax = std::max(errmax, std::abs(xerr[i] / xscale[i]));
				errmax /= eps;

				if (errmax <= 1.0) break;

				htemp = SAFETY * h * std::pow(errmax, PSHRNK);
				h = (h >= 0.0) ? std::max(htemp, 0.1 * h) : std::min(htemp, 0.1 * h);

				if (std::abs(h) < Constants::Eps)
					throw ODESolverError("Stepsize underflow in RK5_CashKarp_Stepper");
			}

			if (errmax > ERRCON)
				_hNext = SAFETY * h * std::pow(errmax, PGROW);
			else
				_hNext = 5.0 * h;

			_t += h;
			_hDone = h;
			_x = xtemp;
		}
	};

} // namespace MML

#endif // MML_ODE_SYSTEM_STEPPERS_H

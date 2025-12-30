#if !defined MPL_PENDULUMS_H
#define MPL_PENDULUMS_H

#include "MMLBase.h"

#include "interfaces/IODESystem.h"

using namespace MML;

namespace MPL
{
	class PendulumODE : public IODESystem
	{
		Real _length;
	public:
		PendulumODE(Real length) : _length(length) {}

		int  getDim() const override { return 2; }
		void derivs(const Real t, const MML::Vector<Real>& x, MML::Vector<Real>& dxdt) const override
		{
			dxdt[0] = x[1];
			dxdt[1] = -9.81 / _length * sin(x[0]);
		}

		// providing explicit variable names for the system
		std::string getVarName(int ind) const override
		{
			switch (ind)
			{
			case 0: return "angle";			// angle of pendulum
			case 1: return "ang_vel";		// angular velocity
			default: return "unknown_var_" + std::to_string(ind);
			}
		}

		double calcPeriodLinearized()
		{
			const double g = 9.81;
			double period = 2 * Constants::PI * std::sqrt(_length / g);
			return period;
		}

		double calcExactPeriod(double initialAngle)
		{
			const double g = 9.81;
			double k = std::sin(initialAngle / 2);
			double ellipticIntegral = std::comp_ellint_1(k);

			double period = 4 * std::sqrt(_length / g) * ellipticIntegral;

			return period;
		}

		// check if the pendulum is in a closed orbit, 
		// i.e. if it will return to the initial angle and angular velocity after some time
		bool isClosedOrbit(Real initAngle, Real initAngVel) const
		{
			const Real g = 9.81;
			
			// Total energy of the pendulum
			Real E = 0.5 * _length * _length * initAngVel * initAngVel + g * _length * (1 - std::cos(initAngle));
			
			// Energy required to reach the upright position (theta = pi)
			Real E_critical = 2 * g * _length;
			
			// Closed orbit if energy is less than critical (no full rotation)
			return E < E_critical;
		}
	};

	class DampedPendulumODE : public IODESystem
	{
		Real _length, _damp_coeff;
	public:
		DampedPendulumODE(Real length, Real damp_coeff)
			: _length(length), _damp_coeff(damp_coeff) {
		}

		int  getDim() const override { return 2; }
		void derivs(const Real t, const MML::Vector<Real>& x,
								MML::Vector<Real>& dxdt) const override
		{
			dxdt[0] = x[1];
			dxdt[1] = -9.81 / _length * sin(x[0]) - _damp_coeff * x[1];
		}

		// calculate critical damping coefficient
		Real calcCriticalDampingCoeff() const
		{
			return 2 * std::sqrt(9.81 / _length);
		}

		// providing explicit variable names for the system
		std::string getVarName(int ind) const override
		{
			switch (ind)
			{
			case 0: return "angle";			// angle of pendulum
			case 1: return "ang_vel";		// angular velocity
			default: return "unknown_var_" + std::to_string(ind);
			}
		}
	};

	class ForcedPendulumODE : public IODESystem
	{
		Real _l, _amp, _period;
	public:
		ForcedPendulumODE(Real length, Real forceAmp, Real forcePeriod)
			: _l(length), _amp(forceAmp), _period(forcePeriod) {
		}
		int  getDim() const override { return 2; }
		void derivs(const Real t, const MML::Vector<Real>& x,
								MML::Vector<Real>& dxdt) const override
		{
			dxdt[0] = x[1];
			dxdt[1] = -9.81 / _l * sin(x[0]) + _amp * cos(Constants::PI * t / _period);
		}

		// providing explicit variable names for the system
		std::string getVarName(int ind) const override
		{
			switch (ind)
			{
			case 0: return "angle";			// angle of pendulum
			case 1: return "ang_vel";		// angular velocity
			default: return "unknown_var_" + std::to_string(ind);
			}
		}
	};

	class DampedForcedPendulumODE : public IODESystem
	{
		Real _l, _damp_coeff, _amp, _period;
	public:
		DampedForcedPendulumODE(Real length, Real damp_coeff, Real forceAmp, Real forcePeriod)
			: _l(length), _damp_coeff(damp_coeff), _amp(forceAmp), _period(forcePeriod) {
		}
		int  getDim() const override { return 2; }
		void derivs(const Real t, const MML::Vector<Real>& x,
								MML::Vector<Real>& dxdt) const override
		{
			dxdt[0] = x[1];
			dxdt[1] = -9.81 / _l * sin(x[0]) - _damp_coeff * x[1] + _amp * cos(Constants::PI * t / _period);
		}

		// calculate critical damping coefficient
		Real calcCriticalDampingCoeff() const
		{
			return 2 * std::sqrt(9.81 / _l);
		}

		// providing explicit variable names for the system
		std::string getVarName(int ind) const override
		{
			switch (ind)
			{
			case 0: return "angle";			// angle of pendulum
			case 1: return "ang_vel";		// angular velocity
			default: return "unknown_var_" + std::to_string(ind);
			}
		}
	};

	class DoublePendulumODE : public IODESystem
	{
		Real _m1, _m2;
		Real _l1, _l2;
	public:
		DoublePendulumODE(Real m1, Real m2, Real l1, Real l2) : _m1(m1), _m2(m2), _l1(l1), _l2(l2) {}

		int  getDim() const override { return 4; }
		void derivs(const Real t, const MML::Vector<Real>& x, MML::Vector<Real>& dxdt) const override
		{
			Real th1 = x[0];
			Real w1 = x[1];
			Real th2 = x[2];
			Real w2 = x[3];

			Real g = 9.81;
			Real divisor = (2 * _m1 + _m2 - _m2 * cos(2 * th1 - 2 * th2));

			dxdt[0] = w1;
			dxdt[1] = (-g * (2 * _m1 + _m2) * sin(th1) - _m2 * g * sin(th1 - 2 * th2) -
								  2 * sin(th1 - th2) * _m2 * (POW2(w2) * _l2 + POW2(w1) * _l1 * cos(th1 - th2))
								) / (_l1 * divisor);
			dxdt[2] = w2;
			dxdt[3] = (2 * sin(th1 - th2) * (POW2(w1) * _l1 * (_m1 + _m2) +
								 g * (_m1 + _m2) * cos(th1) +
								 POW2(w2) * _l2 * _m2 * cos(th1 - th2))
								) / (_l2 * divisor);
		}
	};

	class SphericalPendulumODE : public IODESystem
	{
		Real _l;
	public:
		SphericalPendulumODE() : _l(1.0) {}
		SphericalPendulumODE(Real l) : _l(l) {}

		int  getDim() const override { return 4; }
		void derivs(const Real t, const MML::Vector<Real>& x, MML::Vector<Real>& dxdt) const override
		{
			Real tht = x[0];
			Real w1 = x[1];
			Real phi = x[2];
			Real w2 = x[3];

			Real g = 9.81;

			dxdt[0] = w1;
			dxdt[1] = sin(tht) * cos(tht) * POW2(w2) - g / _l * sin(tht);
			dxdt[2] = w2;
			dxdt[3] = -2.0 * w1 * w2 / tan(tht);
		}
	};

	class SphericalPendulumHamiltonODE : public IODESystem
	{
		Real _m, _l;
	public:
		SphericalPendulumHamiltonODE() : _l(1.0), _m(1.0) {}
		SphericalPendulumHamiltonODE(Real l) : _l(l), _m(1.0) {}
		SphericalPendulumHamiltonODE(Real m, Real l) : _m(m), _l(l) {}

		int  getDim() const override { return 4; }
		void derivs(const Real t, const MML::Vector<Real>& x, MML::Vector<Real>& dxdt) const override
		{
			Real tht = x[0];
			Real phi = x[1];
			Real P_tht = x[2];
			Real P_phi = x[3];

			Real g = 9.81;

			dxdt[0] = P_tht / (_m * POW2(_l));
			dxdt[1] = P_phi / (_m * POW2(_l * sin(tht)));
			dxdt[2] = POW2(P_phi) * cos(tht) / (_m * POW2(_l) * POW3(sin(tht))) - _m * g * _l * sin(tht);
			dxdt[3] = 0.0;
		}
	};
} // namespace MPL

#endif // MPL_PENDULUMS_H
#if !defined MPL_HARMONIC_OSCILATOR_H
#define MPL_HARMONIC_OSCILATOR_H

#include "MMLBase.h"

#include "interfaces/IODESystem.h"

using namespace MML;

namespace MPL
{
	// class representing a general harmonic oscillator
	class HarmonicOscillator : public IODESystem
	{
	private:
		Real _m;		
		Real _k;		

	public:
		HarmonicOscillator(Real mass, Real stiffness)
			: _m(mass), _k(stiffness) {}

		int getDim() const override { return 2; }
		void derivs(const Real t, const Vector<Real>& x, Vector<Real>& dxdt) const override
		{
			dxdt[0] = x[1];									// velocity
			dxdt[1] = -(_k / _m) * x[0];		// acceleration
		}
		std::string getVarName(int ind) const override
		{
			switch (ind)
			{
			case 0: return "Position";
			case 1: return "Velocity";
			default: return IODESystem::getVarName(ind);
			}
		}
		Vector<Real> getInitCond(Real position, Real velocity) const
		{
			Vector<Real> initCond(2);
			initCond[0] = position; // initial position
			initCond[1] = velocity; // initial velocity
			return initCond;
		}

		Real getMass() const { return _m; }
		Real getStiffness() const { return _k; }

		Real getNaturalFrequency() const
		{
			return sqrt(_k / _m); // natural frequency of the harmonic oscillator
		}
		Real getPeriod() const
		{
			return 2 * Constants::PI / getNaturalFrequency(); // period of the harmonic oscillator
		}

		// get RealFunction representing the position as a function of time
		RealFunctionFromStdFunc getPositionFunction(Real initialPosition, Real initialVelocity) const
		{
			Real omega = getNaturalFrequency();
			return RealFunctionFromStdFunc([=](Real t) {
				return initialPosition * cos(omega * t) + (initialVelocity / omega) * sin(omega * t);
			});
		}
		// get RealFunction representing the velocity as a function of time
		RealFunctionFromStdFunc getVelocityFunction(Real initialPosition, Real initialVelocity) const
		{
			Real omega = getNaturalFrequency();
			return RealFunctionFromStdFunc([=](Real t) {
				return -initialPosition * omega * sin(omega * t) + initialVelocity * cos(omega * t);
			});
		}
	};

	class DampedHarmonicOscillator : public IODESystem
	{
		private:
		Real _m;			// mass
		Real _k;			// stiffness
		Real _c;			// damping coefficient
	public:
		DampedHarmonicOscillator(Real mass, Real stiffness, Real damping)
			: _m(mass), _k(stiffness), _c(damping) {}
		
		int getDim() const override { return 2; }
		void derivs(const Real t, const Vector<Real>& x, Vector<Real>& dxdt) const override
		{
			dxdt[0] = x[1];									// velocity
			dxdt[1] = -(_k / _m) * x[0] - (_c / _m) * x[1]; // acceleration with damping
		}
		std::string getVarName(int ind) const override
		{
			switch (ind)
			{
			case 0: return "Position";
			case 1: return "Velocity";
			default: return IODESystem::getVarName(ind);
			}
		}
		Vector<Real> getInitCond(Real position, Real velocity) const
		{
			Vector<Real> initCond(2);
			initCond[0] = position; // initial position
			initCond[1] = velocity; // initial velocity
			return initCond;
		}

		Real getMass() const { return _m; }
		Real getStiffness() const { return _k; }
		Real getDampingCoefficient() const { return _c; }

		Real getNaturalFrequency() const
		{
			return sqrt(_k / _m); // natural frequency of the damped harmonic oscillator
		}
		Real getDampedNaturalFrequency() const
		{
			return sqrt(_k / _m - (_c * _c) / (4 * _m * _m)); // damped natural frequency
		}
		Real getPeriod() const
		{
			return 2 * Constants::PI / getDampedNaturalFrequency(); // period of the damped harmonic oscillator
		}
		Real getDampingRatio() const
		{
			return _c / (2 * sqrt(_m * _k)); // damping ratio
		}
		Real getCriticalDamping() const
		{
			return 2 * sqrt(_m * _k); // critical damping coefficient
		}
		
		bool isOverDamped() const
		{
			return _c > getCriticalDamping(); 
		}
		bool isUnderDamped() const
		{
			return _c < getCriticalDamping(); 
		}
		bool isCriticallyDamped() const
		{
			return _c == getCriticalDamping(); 
		}

		// get RealFunction representing the position as a function of time
		RealFunctionFromStdFunc getPositionFunction(Real initialPosition, Real initialVelocity) const
		{
			Real omega_n = sqrt(_k / _m);
			Real omega = getDampedNaturalFrequency();
			Real zeta = _c / (2 * sqrt(_m * _k));

			return RealFunctionFromStdFunc([=](Real t) {
				return exp(-zeta * omega_n * t) *
					(initialPosition * cos(omega * t) +
						(initialVelocity + zeta * omega_n * initialPosition) / omega * sin(omega * t));
				});
		}
		// get RealFunction representing the velocity as a function of time
		RealFunctionFromStdFunc getVelocityFunction(Real initialPosition, Real initialVelocity) const
		{
			Real omega_n = sqrt(_k / _m);
			Real zeta = _c / (2 * sqrt(_m * _k));
			Real omega_d = sqrt(omega_n * omega_n - zeta * zeta * omega_n * omega_n);

			return RealFunctionFromStdFunc([=](Real t) {
				return exp(-zeta * omega_n * t) *
					(initialVelocity * cos(omega_d * t)
						+ (-initialPosition * omega_d - zeta * omega_n * initialVelocity) * sin(omega_d * t));
				});
		}
	};

	// representing a forced harmonic oscillator, but without damping
	class ForcedHarmonicOscillator : public IODESystem
	{
		private:
		Real _m;			// mass
		Real _k;			// stiffness
		Real _F0;			// amplitude of the forcing function
		Real _omegaF; // frequency of the forcing function
	public:
		ForcedHarmonicOscillator(Real mass, Real stiffness, Real F0, Real omegaF)
			: _m(mass), _k(stiffness), _F0(F0), _omegaF(omegaF) {}
		
		int getDim() const override { return 2; }
		void derivs(const Real t, const Vector<Real>& x, Vector<Real>& dxdt) const override
		{
			dxdt[0] = x[1];									// velocity
			dxdt[1] = -(_k / _m) * x[0] + (_F0 / _m) * cos(_omegaF * t); // acceleration with forcing
		}
		std::string getVarName(int ind) const override
		{
			switch (ind)
			{
			case 0: return "Position";
			case 1: return "Velocity";
			default: return IODESystem::getVarName(ind);
			}
		}
		Vector<Real> getInitCond(Real position, Real velocity) const
		{
			Vector<Real> initCond(2);
			initCond[0] = position; // initial position
			initCond[1] = velocity; // initial velocity
			return initCond;
		}
		
		Real getMass() const { return _m; }
		Real getStiffness() const { return _k; }
		Real getForceAmplitude() const { return _F0; }
		Real getForceFrequency() const { return _omegaF; }
		
		Real getNaturalFrequency() const
		{
			return sqrt(_k / _m); // natural frequency of the forced harmonic oscillator
		}
		Real getPeriod() const
		{
			return 2 * Constants::PI / getNaturalFrequency(); // period of the forced harmonic oscillator
		}
		Real getResonanceFrequency() const
		{
			return _omegaF; // resonance frequency is the frequency of the forcing function
		}

		// get RealFunction representing the position as a function of time
		RealFunctionFromStdFunc getPositionFunction(Real initialPosition, Real initialVelocity) const
		{
			Real omega = getNaturalFrequency();
			Real delta = omega * omega - _omegaF * _omegaF;
			return RealFunctionFromStdFunc([=](Real t) {
				return initialPosition * cos(omega * t)
					+ (initialVelocity / omega) * sin(omega * t)
					+ (_F0 / (_m * delta)) * (cos(_omegaF * t) - cos(omega * t));
				});
		}
		// get RealFunction representing the velocity as a function of time
		RealFunctionFromStdFunc getVelocityFunction(Real initialPosition, Real initialVelocity) const
		{
			Real omega = getNaturalFrequency();
			Real delta = omega * omega - _omegaF * _omegaF;
			return RealFunctionFromStdFunc([=](Real t) {
				return -initialPosition * omega * sin(omega * t)
					+ initialVelocity * cos(omega * t)
					+ (_F0 / _m) * (omega * sin(omega * t) - _omegaF * sin(_omegaF * t)) / delta;
				});
		}
	};

	// representing a forced harmonic oscillator with damping
	class DampedForcedHarmonicOscillator : public IODESystem
	{
		private:
		Real _m;			// mass
		Real _k;			// stiffness
		Real _c;			// damping coefficient
		Real _F0;			// amplitude of the forcing function
		Real _omegaF; // frequency of the forcing function

	public:
		DampedForcedHarmonicOscillator(Real mass, Real stiffness, Real damping, Real F0, Real omegaF)
			: _m(mass), _k(stiffness), _c(damping), _F0(F0), _omegaF(omegaF) {}
		
		int getDim() const override { return 2; }
		void derivs(const Real t, const Vector<Real>& x, Vector<Real>& dxdt) const override
		{
			dxdt[0] = x[1];									// velocity
			dxdt[1] = -(_k / _m) * x[0] - (_c / _m) * x[1] + (_F0 / _m) * cos(_omegaF * t); // acceleration with forcing
		}
		std::string getVarName(int ind) const override
		{
			switch (ind)
			{
			case 0: return "Position";
			case 1: return "Velocity";
			default: return IODESystem::getVarName(ind);
			}
		}
		Vector<Real> getInitCond(Real position, Real velocity) const
		{
			Vector<Real> initCond(2);
			initCond[0] = position; // initial position
			initCond[1] = velocity; // initial velocity
			return initCond;
		}

		Real getMass() const { return _m; }
		Real getStiffness() const { return _k; }
		Real getDampingCoefficient() const { return _c; }
		Real getForceAmplitude() const { return _F0; }
		Real getForceFrequency() const { return _omegaF; }
		
		Real getNaturalFrequency() const
		{
			return sqrt(_k / _m); // natural frequency of the forced harmonic oscillator
		}
		Real getDampedNaturalFrequency() const
		{
			return sqrt(_k / _m - (_c * _c) / (4 * _m * _m)); // damped natural frequency
		}
		Real getPeriod() const
		{
			return 2 * Constants::PI / getDampedNaturalFrequency(); // period of the forced harmonic oscillator
		}
		Real getDampingRatio() const
		{
			return _c / (2 * sqrt(_m * _k)); // damping ratio
		}
		Real getCriticalDamping() const
		{
			return 2 * sqrt(_m * _k); // critical damping coefficient
		}
		
		bool isOverDamped() const
		{
			return _c > getCriticalDamping(); 
		}
		bool isUnderDamped() const
		{
			return _c < getCriticalDamping(); 
		}
		bool isCriticallyDamped() const
		{
			return _c == getCriticalDamping(); 
		}

		// get RealFunction representing the position as a function of time
		RealFunctionFromStdFunc getPositionFunction(Real initialPosition, Real initialVelocity) const
		{
			Real omega_n = sqrt(_k / _m); // natural frequency
			Real zeta = _c / (2 * sqrt(_m * _k)); // damping ratio
			Real omega_d = omega_n * sqrt(1 - zeta * zeta); // damped natural frequency

			// Transient (homogeneous) solution coefficients
			Real A = initialPosition;
			Real B = (initialVelocity + zeta * omega_n * initialPosition) / omega_d;

			// Steady-state (particular) solution amplitude and phase
			Real denom = sqrt(pow(omega_n * omega_n - _omegaF * _omegaF, 2) + pow(2 * zeta * omega_n * _omegaF, 2));
			Real Xp = (_F0 / _m) / denom;
			Real phi = atan2(2 * zeta * omega_n * _omegaF, omega_n * omega_n - _omegaF * _omegaF);

			return RealFunctionFromStdFunc([=](Real t) {
				Real transient = exp(-zeta * omega_n * t) * (A * cos(omega_d * t) + B * sin(omega_d * t));
				Real steady = Xp * cos(_omegaF * t - phi);
				return transient + steady;
				});
		}
		// get RealFunction representing the velocity as a function of time
		RealFunctionFromStdFunc getVelocityFunction(Real initialPosition, Real initialVelocity) const
		{
			Real omega_n = sqrt(_k / _m);
			Real zeta = _c / (2 * sqrt(_m * _k));
			Real omega_d = omega_n * sqrt(1 - zeta * zeta);
			Real A = initialPosition;
			Real B = (initialVelocity + zeta * omega_n * initialPosition) / omega_d;
			Real denom = sqrt(pow(omega_n * omega_n - _omegaF * _omegaF, 2) + pow(2 * zeta * omega_n * _omegaF, 2));
			Real phi = atan2(2 * zeta * omega_n * _omegaF, omega_n * omega_n - _omegaF * _omegaF);

			return RealFunctionFromStdFunc([=](Real t) {
				Real exp_term = exp(-zeta * omega_n * t);
				Real d_transient =
					exp_term * (
						-A * zeta * omega_n * cos(omega_d * t)
						- A * omega_d * sin(omega_d * t)
						- B * zeta * omega_n * sin(omega_d * t)
						+ B * omega_d * cos(omega_d * t)
						);
				Real d_steady = -(_F0 / _m) / denom * _omegaF * sin(_omegaF * t - phi);
				return d_transient + d_steady;
				});
		}
	};
}

#endif
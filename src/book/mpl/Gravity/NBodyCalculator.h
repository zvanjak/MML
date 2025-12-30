#if !defined MPL_NBODY_CALCULATOR_H
#define MPL_NBODY_CALCULATOR_H

#include "MMLBase.h"

#include "interfaces/IODESystem.h"

#include "base/BaseUtils.h"
#include "base/Function.h"
#include "base/VectorTypes.h"

#include "algorithms/ODESystemSolver.h"
#include "algorithms/ODESystemStepCalculators.h"
#include "algorithms/ODESystemSteppers.h"

#include "GravityBase.h"
#include "NBodyBase.h"

using namespace MML;

namespace MPL
{
	class NBodyGravityCalculator
	{
	private:
		Real _G;
		NBodyState _state;

	public:
		NBodyGravityCalculator(NBodyState state) : _G(6.67430e-11), _state(state) {}
		NBodyGravityCalculator(Real G, NBodyState state) : _G(G), _state(state) {}
		NBodyGravityCalculator(NBodyGravitySimConfig config) : _G(config.G()), _state(config.InitState()) {}
		NBodyGravityCalculator(Real G, NBodyGravitySimConfig config) : _G(G), _state(config.InitState()) {}

		Vec3Cart ForceOnBody(int i, const Vec3Cart& x)
		{
			Vec3Cart force;
			for (int j = 0; j < _state.NumBodies(); j++)
			{
				if (i != j)
				{
					Vec3Cart radialVec = x - _state.Pos(j);
					force = force + _G * _state.Mass(i) * _state.Mass(j) / POW3(radialVec.NormL2()) * radialVec;
				}
			}
			return force;
		}

		Vec3Cart ForceAtPoint(const Vec3Cart& x)
		{
			Vec3Cart force;
			for (int i = 0; i < _state.NumBodies(); i++)
			{
				Vec3Cart radialVec = x - _state.Pos(i);
				force = force + _G * _state.Mass(i) / POW3(radialVec.NormL2()) * radialVec;
			}
			return force;
		}

		Real PotentialAtPoint(const Vec3Cart& x)
		{
			Real pot = 0.0;
			for (int i = 0; i < _state.NumBodies(); i++)
			{
				Vec3Cart radialVec = x - _state.Pos(i);
				pot += _G * _state.Mass(i) / radialVec.NormL2();
			}
			return pot;
		}
	};

	class NBodyGravityPotentialField : public IScalarFunction<3>
	{
	private:
		Real _G;
		NBodyGravitySimConfig _config;
	public:
		NBodyGravityPotentialField(NBodyGravitySimConfig config) : _G(6.67430e-11), _config(config) {}
		NBodyGravityPotentialField(Real G, NBodyGravitySimConfig config) : _G(G), _config(config) {}

		Real operator()(const VectorN<Real, 3>& x) const
		{
			Real pot = 0.0;
			for (int i = 0; i < _config.NumBodies(); i++)
			{
				// FIX!
				// pot += Fields::InverseRadialPotentialFieldCart(_G * _config.Mass(i), _config.Position(i) - x);
			}
			return pot;
		}
	};
	class NBodyGravityForceField : public IVectorFunction<3>
	{
	private:
		Real _G;
		NBodyGravitySimConfig _config;
	public:
		NBodyGravityForceField(NBodyGravitySimConfig config) : _G(6.67430e-11), _config(config) {}
		NBodyGravityForceField(Real G, NBodyGravitySimConfig config) : _G(G), _config(config) {}

		VectorN<Real, 3> operator()(const VectorN<Real, 3>& x) const
		{
			VectorN<Real, 3> force(0.0);
			for (int i = 0; i < _config.NumBodies(); i++)
			{
				// FIX!
				// force = force + Fields::InverseRadialPotentialForceFieldCart(_G * _config.Mass(i), _config.Position(i) - x);
			}
			return force;
		}
	};

}

#endif // MPL_NBODY_CALCULATOR_H
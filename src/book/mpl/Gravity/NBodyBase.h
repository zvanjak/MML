#if !defined MPL_NBODY_BASE_H
#define MPL_NBODY_BASE_H

#include "MMLBase.h"

#include "base/BaseUtils.h"
#include "base/Function.h"
#include "base/VectorTypes.h"

#include "GravityBase.h"

using namespace MML;

namespace MPL
{
	// representing complete state of N-body system
	class NBodyState
	{
	public:
		std::vector<GravityBodyState> _bodies;

	public:
		void addBody(const GravityBodyState& body)
		{
			_bodies.push_back(body);
		}
		void addBody(Real mass, Vec3Cart position)
		{
			_bodies.push_back(GravityBodyState(mass, position, Vec3Cart(0, 0, 0)));
		}
		void addBody(Real mass, Vec3Cart position, Vec3Cart velocity)
		{
			_bodies.push_back(GravityBodyState(mass, position, velocity));
		}

		int NumBodies() const { return static_cast<int>(_bodies.size()); }

		Real		 Mass(int i)	const { return _bodies[i].Mass(); }
		Vec3Cart Pos(int i)		const { return _bodies[i].R(); }
		Vec3Cart Vel(int i)		const { return _bodies[i].V(); }

		void SetPosition(int i, Vec3Cart pos) { _bodies[i].R() = pos; }
		void SetVelocity(int i, Vec3Cart vel) { _bodies[i].V() = vel; }

		Vec3Cart CentreOfMassPos() const
		{
			Vec3Cart com(0, 0, 0);
			Real totalMass = 0.0;
			for (const auto& body : _bodies)
			{
				com += body.Mass() * body.R();
				totalMass += body.Mass();
			}
			return com / totalMass;
		}
		Vec3Cart CentreOfMassVel() const
		{
			Vec3Cart com(0, 0, 0);
			Real totalMass = 0.0;
			for (const auto& body : _bodies)
			{
				com += body.Mass() * body.V();
				totalMass += body.Mass();
			}
			return com / totalMass;
		}

		Vec3Cart LinearMomentum() const
		{
			Vec3Cart lm(0, 0, 0);
			for (const auto& body : _bodies)
			{
				lm += body.Mass() * body.V();
			}
			return lm;
		}
		Vec3Cart AngularMomentum(Vec3Cart origin) const
		{
			Vec3Cart am(0, 0, 0);
			for (const auto& body : _bodies)
			{
				Vec3Cart r = body.R() - origin;
				am += VectorProduct(r, body.Mass() * body.V());
			}
			return am;
		}
		Vec3Cart AngularMomentumCM() const
		{
			Vec3Cart am(0, 0, 0);
			Vec3Cart com = CentreOfMassPos();
			for (const auto& body : _bodies)
			{
				Vec3Cart r = body.R() - com;
				am += VectorProduct(r, body.Mass() * body.V());
			}
			return am;
		}
		Vec3Cart AngularMomentumOrigin() const
		{
			Vec3Cart am(0, 0, 0);
			for (const auto& body : _bodies)
			{
				Vec3Cart r = body.R();
				am += VectorProduct(r, body.Mass() * body.V());
			}
			return am;
		}

		Real TotalKineticEnergy() const
		{
			Real ke = 0.0;
			for (const auto& body : _bodies)
			{
				ke += 0.5 * body.Mass() * POW2(body.V().NormL2());
			}
			return ke;
		}
		Real TotalPotentialEnergy() const
		{
			Real pe = 0.0;
			for (size_t i = 0; i < _bodies.size(); ++i)
			{
				for (size_t j = i + 1; j < _bodies.size(); ++j)
				{
					Vec3Cart r_ij = _bodies[j].R() - _bodies[i].R();
					pe -= (_bodies[i].Mass() * _bodies[j].Mass()) / r_ij.NormL2();
				}
			}
			return pe;
		}
	};

	// initial configuration for N-body gravity simulation
	class NBodyGravitySimConfig
	{
		Real _G;
		NBodyState _initState;

	public:
		NBodyGravitySimConfig() : _G(6.67430e-11) {}
		NBodyGravitySimConfig(Real G) : _G(G) {}

		Real G() const { return _G; }

		int NumBodies() const { return static_cast<int>(_initState._bodies.size()); }

		const NBodyState& InitState() const { return _initState; }
		NBodyState& InitStateAcc() { return _initState; }

		void AddBody(Real mass, Vec3Cart position, Vec3Cart velocity)
		{
			_initState.addBody(mass, position, velocity);
		}
		void AddBody(Real mass, Vec3Cart position, Vec3Cart velocity, std::string color, Real radius)
		{
			_initState._bodies.push_back(GravityBodyState(mass, position, velocity, color, radius));
		}

		Real Mass(int i) const { return _initState.Mass(i); }
		Real Radius(int i) const { return _initState._bodies[i].Radius(); }
		std::string Color(int i) const { return _initState._bodies[i].Color(); }

		Vec3Cart Position(int i) const { return _initState.Pos(i); }
		Vec3Cart Velocity(int i) const { return _initState.Vel(i); }

		void SetPosition(int i, Vec3Cart pos) { _initState.SetPosition(i, pos); }
		void SetVelocity(int i, Vec3Cart vel) { _initState.SetVelocity(i, vel); }

		Vector<Real> getInitCond()
		{
			Vector<Real> initCond(6 * NumBodies());
			for (int i = 0; i < NumBodies(); i++)
			{
				initCond[3 * i]			= Position(i).X();
				initCond[3 * i + 1] = Position(i).Y();
				initCond[3 * i + 2] = Position(i).Z();
				initCond[3 * NumBodies() + 3 * i]			= Velocity(i).X();
				initCond[3 * NumBodies() + 3 * i + 1] = Velocity(i).Y();
				initCond[3 * NumBodies() + 3 * i + 2] = Velocity(i).Z();
			}
			return initCond;
		}
	};

}

#endif 
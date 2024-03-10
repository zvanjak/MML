#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "interfaces/IODESystem.h"

#include "base/Geometry3D.h"

#include "core/Derivation.h"
#include "core/Fields.h"
#include "core/Function.h"
#include "core/InterpolatedFunction.h"
#include "core/Serializer.h"
#include "core/Visualizer.h"
#endif


using namespace MML;

class GravityMass
{
public:
	Real _mass;
	Vector3Cartesian _position;
	Vector3Cartesian _velocity;

	GravityMass(const Real& mass, const Vector3Cartesian& position) : _mass(mass), _position(position) {}
	GravityMass(const Real& mass, const Vector3Cartesian& position, const Vector3Cartesian& velocity) : _mass(mass), _position(position), _velocity(velocity) {}
};

class GravityMultibodyConfigCart
{
	std::vector<GravityMass> _masses;
public:
	int NumBodies() const { return (int)_masses.size(); }

	void AddBody(Real mass, Vector3Cartesian position, Vector3Cartesian velocity)
	{
		_masses.push_back(GravityMass(mass, position, velocity));
	}

	Real Mass(int i) const { return _masses[i]._mass; }
	Vector3Cartesian Position(int i) const { return _masses[i]._position; }
	Vector3Cartesian Velocity(int i) const { return _masses[i]._velocity; }

	void SetPosition(int i, Vector3Cartesian pos) { _masses[i]._position = pos; }
	void SetVelocity(int i, Vector3Cartesian vel) { _masses[i]._velocity = vel; }
};

class MultibodyGravityPotentialFieldCart : public IScalarFunction<3>
{
private:
	Real _G;
	GravityMultibodyConfigCart _config;
public:
	MultibodyGravityPotentialFieldCart(GravityMultibodyConfigCart config) : _G(6.67430e-11), _config(config) {}
	MultibodyGravityPotentialFieldCart(Real G, GravityMultibodyConfigCart config) : _G(G), _config(config) {}

	Real operator()(const VectorN<Real, 3>& x) const
	{
		Real pot = 0.0;
		for (int i = 0; i < _config.NumBodies(); i++)
		{
			pot += Fields::InverseRadialPotentialFieldCart(_G * _config.Mass(i), _config.Position(i) - x);
		}
		return pot;
	}
};
class MultibodyGravityForceFieldCart : public IVectorFunction<3>
{
private:
	Real _G;
	GravityMultibodyConfigCart _config;
public:
	MultibodyGravityForceFieldCart(GravityMultibodyConfigCart config) : _G(6.67430e-11), _config(config) {}
	MultibodyGravityForceFieldCart(Real G, GravityMultibodyConfigCart config) : _G(G), _config(config) {}

	VectorN<Real, 3> operator()(const VectorN<Real, 3>& x) const
	{
		VectorN<Real, 3> force(0.0);
		for (int i = 0; i < _config.NumBodies(); i++)
		{
			force = force + Fields::InverseRadialPotentialForceFieldCart(_G * _config.Mass(i), _config.Position(i) - x);
		}
		return force;
	}
};

// set up ODESystem za simuliranje sustava tijela u gravitacijskom polju
// a ima i druga opcija - u DANOM sustavu tijela, simulirati neko trece tijelo!
class GravitySystemMotionSolverSimple
{
protected:
	GravityMultibodyConfigCart _config;

public:
	GravitySystemMotionSolverSimple(GravityMultibodyConfigCart inConfig) : _config(inConfig) { }

	int getDim() const { return 3 * _config.NumBodies(); }

	void simulateOneStep(const Real dt)
	{
		std::vector<Vector3Cartesian> force(_config.NumBodies(), Vector3Cartesian(0.0, 0.0, 0.0));

		for (int i = 0; i < _config.NumBodies(); i++) {
			// calculate force on body i
			for (int j = 0; j < _config.NumBodies(); j++) {
				if (i != j)
					force[i] = force[i] + Fields::InverseRadialPotentialForceFieldCart(10000 * _config.Mass(j), _config.Position(i) - _config.Position(j));
			}
		}
		// after we calculated forces on all bodies, we can update velocities and positions
		for (int i = 0; i < _config.NumBodies(); i++)
		{
			_config.SetVelocity(i, _config.Velocity(i) + force[i] * dt / _config.Mass(i));
			_config.SetPosition(i, _config.Position(i) + _config.Velocity(i) * dt);
		}
	}
	std::vector<std::vector<VectorN< Real, 3 >>> simulate(const Real dt, const int steps)
	{
		std::vector<std::vector<VectorN<Real,3>>> trajectories(_config.NumBodies());

		// saving initial positions
		for (int i = 0; i < _config.NumBodies(); i++)
			trajectories[i].push_back(_config.Position(i));

		for (int i = 0; i < steps; i++)
		{
			simulateOneStep(dt);

			for (int i = 0; i < _config.NumBodies(); i++)
				trajectories[i].push_back(_config.Position(i));
		}
		return trajectories;
	}
};

class GravitySystemMotionSolver : public IODESystem
{
protected:
	GravityMultibodyConfigCart _initialConfig;

public:
	GravitySystemMotionSolver(GravityMultibodyConfigCart inConfig) : _initialConfig(inConfig) { }

	int getDim() const { return 3 * _initialConfig.NumBodies(); }
	void derivs(const Real t, const Vector<Real>& x, Vector<Real>& dxdt) const
	{
		// moramo popuniti dxdt vektor s N * 6 derivacija nasih varijabli (x, y, z, vx, vy, vz)
		for (int i = 0; i < _initialConfig.NumBodies(); i++)
		{
			Vector3Cartesian force(0, 0, 0);
			// calculating force on body i
			for (int j = 0; j < _initialConfig.NumBodies(); j++)
			{
				if (i != j)      // self-force would be infinite :)
				{
					Vector3Cartesian vec_dist(dxdt[6 * j] - dxdt[6 * i], dxdt[6 * j + 2] - dxdt[6 * i + 2], dxdt[6 * j + 4] - dxdt[6 * i + 4]);
					Vector3Cartesian f = Fields::InverseRadialPotentialForceFieldCart(6.67430e-11 * _initialConfig.Mass(i) * _initialConfig.Mass(j), vec_dist);
					force = force + f;
				}
			}

			// x coord
			dxdt[6 * i] = x[6 * i + 1];
			dxdt[6 * i + 1] = 1 / _initialConfig.Mass(i) * force.X();
			// y coord
			dxdt[6 * i + 2] = x[6 * i + 3];
			dxdt[6 * i + 3] = 1 / _initialConfig.Mass(i) * force.Y();
			// z coord
			dxdt[6 * i + 4] = x[6 * i + 5];
			dxdt[6 * i + 5] = 1 / _initialConfig.Mass(i) * force.Z();
		}
	}
};

void Example4_Gravity_field_visualization()
{
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "****            EXAMPLE 4 - gravity field investigations           ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	// prvo, simulacija 5 tijela

	GravityMultibodyConfigCart config;
	config.AddBody(1000, Vector3Cartesian{ 0.0,    0.0,    0.0 }, Vector3Cartesian{ 0.0,   0.0,   0.0 });
	config.AddBody(20, Vector3Cartesian{ -110.0,  -50.0,   10.0 }, Vector3Cartesian{ 0.0,   50,   0.0 });
	config.AddBody(10, Vector3Cartesian{ 130.0,   50.0,   70.0 }, Vector3Cartesian{ 0.0,   -50,   0 });
	config.AddBody(20, Vector3Cartesian{ -20.0,  100.0, -110.0 }, Vector3Cartesian{ 50,   0.0,   0.0 });
	config.AddBody(10, Vector3Cartesian{ 70.0, -110.0,   70.0 }, Vector3Cartesian{ -50,   50,   50.0 });

	GravitySystemMotionSolverSimple solver(config);

	const Real dt = 0.1;
	const int steps = 2000;

	auto res = solver.simulate(dt, steps);

	for (int i = 0; i < config.NumBodies(); i++)
	{
		Serializer::SaveAsParamCurve<3>(res[i], "PARAMETRIC_CURVE_CARTESIAN_3D", "Body" + std::to_string(i+1), 0.0, dt * steps, steps + 1,
																		GLOB_PATH_ResultFiles + "body" + std::to_string(i) + ".txt");
	}

	Visualizer::VisualizeMultiParamCurve3D({"body0.txt", "body1.txt", "body2.txt", "body3.txt", "body4.txt" });

	// drugo, Voyager kroz solarni sustav
}
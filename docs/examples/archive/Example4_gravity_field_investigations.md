# Gravity field investigations and visualizations

[C++ file with implementation](/docs_examples/examples/example4_gravity_field_investigations.cpp)

Classes for defining gravity field configuration.
~~~C++
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
~~~
 
Gravity solver class that calculates forces and updates positions and velocities of bodies.

~~~C++
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
 ~~~

Example of usage:
~~~C++
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
Serializer::SaveAsParamCurve<3>(res[i], "PARAMETRIC_CURVE_CARTESIAN_3D", 0.0, dt*steps, steps+1, 
																GLOB_PATH_ResultFiles + "body" + std::to_string(i) + ".txt");
}

Visualizer::VisualizeMultiParamCurve3D({"body0.txt", "body1.txt", "body2.txt", "body3.txt", "body4.txt" });
~~~

Result of the simulation is shown in the figure below showing the trajectories of the bodies in the gravity field.

![My Image](/docs/images/example4_gravity_trajectories.png)

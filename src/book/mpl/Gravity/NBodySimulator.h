#if !defined MPL_NBODY_SIMULATOR_H
#define MPL_NBODY_SIMULATOR_H

#include "MMLBase.h"

#include "interfaces/IODESystem.h"

#include "base/BaseUtils.h"
#include "base/Random.h"
#include "base/Function.h"
#include "base/VectorTypes.h"

#include "algorithms/ODESystemSolver.h"
#include "algorithms/ODESystemStepCalculators.h"
#include "algorithms/ODESystemSteppers.h"

#include "tools/Serializer.h"
#include "tools/Visualizer.h"

#include "../Base/PhysicalConstants.h"
#include "../Base/SolarSystem.h"

#include "GravityBase.h"
#include "NBodyBase.h"
#include "NBodyCalculator.h"

using namespace MML;

namespace MPL
{
	/*********************************************************************************************/
	class NBodyGravityConfigGenerator
	{
	public:
		// configuration with five bodies
		static NBodyGravitySimConfig Config1_five_bodies()
		{
			NBodyGravitySimConfig config(1.0);

			config.AddBody(10000, Vec3Cart{   0.0,    0.0,    0.0 }, Vec3Cart{  0.55, -0.7,  0.9 });
			config.AddBody(   20, Vec3Cart{-110.0,  -50.0,   10.0 }, Vec3Cart{  2.0,   5,   -3.0 });
			config.AddBody(   10, Vec3Cart{ 130.0,   50.0,   70.0 }, Vec3Cart{  2.0,  -5,   -2 });
			config.AddBody(   20, Vec3Cart{ -20.0, -100.0, -110.0 }, Vec3Cart{ -5,    -2.0,  3.0 });
			config.AddBody(   10, Vec3Cart{  70.0, -110.0,   70.0 }, Vec3Cart{ -3,     3,    3.0 });

			return config;
		}

		static NBodyGravitySimConfig Config2_Solar_system()
		{
			// Physical constants
			constexpr double G_SI = PhyConst::GravityConstant; // m^3 kg^-1 s^-2
			constexpr double million_km = 1.0e9;		// m
			constexpr double year = 3.15576e7;			// s

			double M_jup = SolarSystem::GetBodyInfo("Jupiter").Mass_kg; // kg
			double M_sun = SolarSystem::GetBodyInfo("Sun").Mass_kg; // kg

			// G in units: [ (million km)^3 / (Jupiter mass * year^2) ]
			double G = G_SI * (M_jup) * (year * year) / (million_km * million_km * million_km);

			NBodyGravitySimConfig config(G);

			// Sun
			config.AddBody(M_sun / M_jup, Vec3Cart{ 0, 0, 0 }, Vec3Cart{ 0, 0, 0 }, "Yellow", 30);

			// Planets: mass (in Jupiter masses), semi-major axis (in million km)
			struct Planet {
				const char* name;
				double mass_Mjup;
				double dist_million_km;
				std::string color;
				double radius;
			} planets[] = {
					{"Mercury", 1.6601e-7, 57.91,  "Black", 2.5},
					{"Venus",   2.447e-6,  108.21, "Orange", 6},
					{"Earth",   3.003e-6,  149.60, "Blue", 6},
					{"Mars",    3.227e-7,  227.92, "Red", 4},
					{"Jupiter", 1.0,       778.57, "Brown", 50},
					{"Saturn",  0.299,     1433.53, "Brown", 40},
					{"Uranus",  0.0457,    2872.46, "LightBlue", 20},
					{"Neptune", 0.0539,    4495.06, "Green", 20}
			};

			double Mtot = M_sun / M_jup;
			for (const auto& p : planets) Mtot += p.mass_Mjup;

			for (const auto& p : planets)
			{
				// Circular orbit velocity (in million km/year)
				// v = sqrt(G * M_sun/M_jup / r)
				double v = std::sqrt(G * (M_sun / M_jup) / p.dist_million_km);

				// Place planet at (a, 0, 0), velocity (0, v, 0)
				config.AddBody(
					p.mass_Mjup,
					Vec3Cart{ Real(p.dist_million_km), Real(0), Real(0) },
					Vec3Cart{ Real(0), Real(v), Real(10) },
					p.color,
					p.radius
				);
			}

			return config;
		}

		// configuration of N bodies around massive object
		static NBodyGravitySimConfig Config3_N_bodies_around_massive_object()
		{
			NBodyGravitySimConfig config(1.0);

			// adding central body
			config.AddBody(10000, Vec3Cart{ 0.0, 0.0, 0.0 }, Vec3Cart{ 0.0, 0.0, 0.0 }, "Black", 10);

			// create randomly 10 bodies within cube of side 100
			for (int i = 0; i < 100; i++)
			{
				Real mass = 10 + Random::UniformReal(0.0, 20.0);
				Real radius = 10 / std::pow(10000 / mass, 1. / 3);

				Real rad = Random::UniformReal(100.0, 300.0);
				Real theta = Random::UniformReal(0.0, 2 * Constants::PI);
				Real phi = Random::UniformReal(0.0, Constants::PI);
				Vec3Cart pos = Vec3Cart{ rad * sin(phi) * cos(theta), rad * sin(phi) * sin(theta), rad * cos(phi) };

				//Vec3Cart pos = Vec3Cart{ Random::UniformReal(-150.0, 150.0), Random::UniformReal(-150.0, 150.0), Random::UniformReal(-150.0, 150.0) };
				// we need velocity that is tangential to the radius vPector
				Vec3Cart normal = pos.GetAsUnitVector();

				// random velocity vector
				Real maxVel = 30;
				Vec3Cart vel = Vec3Cart{ Random::UniformReal(-maxVel, maxVel), Random::UniformReal(-maxVel, maxVel), Random::UniformReal(-maxVel, maxVel) };

				// let's get rid of velocity in the direction of radius vector
				Vec3Cart vel2 = vel - normal * (vel * normal) / POW2(normal.NormL2());

				Vec3Cart vel_final = vel2 * vel.NormL2() / vel2.NormL2();

				config.AddBody(mass, pos, vel_final, "Red", radius);
			}
			return config;
		}
	};

	/*********************************************************************************************/
	class NBodyGravitySimulationResults
	{
	public:
		const NBodyGravitySimConfig& _config;

		Real _duration;
		Vector<Real>				_vecTimes;
		Vector<NBodyState> _vecStates;

		NBodyGravitySimulationResults(const NBodyGravitySimConfig& config)
			: _config(config), _duration(0)
		{
		}
		NBodyGravitySimulationResults(const NBodyGravitySimConfig& config, const Real& duration)
			: _config(config), _duration(duration)
		{
		}

		int NumBodies() const { return _config.NumBodies(); }
		int NumSteps() const	{ return _vecStates.size(); }

		Real				Time(int i)  const { return _vecTimes[i]; }
		NBodyState	State(int i) const { return _vecStates[i]; }

		Vector<Real>			getTimes() const { return _vecTimes; }
		std::vector<Vec3Cart> getBodyPos(int i) const
		{
			std::vector<Vec3Cart> res(_vecStates.size());
			for (int j = 0; j < _vecStates.size(); j++)
				res[j] = _vecStates[j].Pos(i);
			return res;
		}
		std::vector<Vec3Cart> getBodyVel(int i) const
		{
			std::vector<Vec3Cart> res(_vecStates.size());
			for (int j = 0; j < _vecStates.size(); j++)
				res[j] = _vecStates[j].Vel(i);
			return res;
		}

		Vector<Real> getPosX(int i) const
		{
			Vector<Real> res(_vecStates.size());
			for (int j = 0; j < _vecStates.size(); j++)
				res[j] = _vecStates[j].Pos(i).X();
			return res;
		}
		Vector<Real> getPosY(int i) const
		{
			Vector<Real> res(_vecStates.size());
			for (int j = 0; j < _vecStates.size(); j++)
				res[j] = _vecStates[j].Pos(i).Y();
			return res;
		}
		Vector<Real> getPosZ(int i) const
		{
			Vector<Real> res(_vecStates.size());
			for (int j = 0; j < _vecStates.size(); j++)
				res[j] = _vecStates[j].Pos(i).Z();
			return res;
		}

		Vector<Real> getVelX(int i) const
		{
			Vector<Real> res(_vecStates.size());
			for (int j = 0; j < _vecStates.size(); j++)
				res[j] = _vecStates[j].Vel(i).X();
			return res;
		}
		Vector<Real> getVelY(int i) const
		{
			Vector<Real> res(_vecStates.size());
			for (int j = 0; j < _vecStates.size(); j++)
				res[j] = _vecStates[j].Vel(i).Y();
			return res;
		}
		Vector<Real> getVelZ(int i) const
		{
			Vector<Real> res(_vecStates.size());
			for (int j = 0; j < _vecStates.size(); j++)
				res[j] = _vecStates[j].Vel(i).Z();
			return res;
		}

		Vector<Real> getCMPosX() const
		{
			Vector<Real> res(_vecStates.size());
			for (int j = 0; j < _vecStates.size(); j++)
				res[j] = _vecStates[j].CentreOfMassPos().X();
			return res;
		}
		Vector<Real> getCMPosY() const
		{
			Vector<Real> res(_vecStates.size());
			for (int j = 0; j < _vecStates.size(); j++)
				res[j] = _vecStates[j].CentreOfMassPos().Y();
			return res;
		}
		Vector<Real> getCMPosZ() const
		{
			Vector<Real> res(_vecStates.size());
			for (int j = 0; j < _vecStates.size(); j++)
				res[j] = _vecStates[j].CentreOfMassPos().Z();
			return res;
		}

		Vector<Real> getCMVX() const
		{
			Vector<Real> res(_vecStates.size());
			for (int j = 0; j < _vecStates.size(); j++)
				res[j] = _vecStates[j].CentreOfMassVel().X();
			return res;
		}
		Vector<Real> getCMVY() const
		{
			Vector<Real> res(_vecStates.size());
			for (int j = 0; j < _vecStates.size(); j++)
				res[j] = _vecStates[j].CentreOfMassVel().Y();
			return res;
		}
		Vector<Real> getCMVZ() const
		{
			Vector<Real> res(_vecStates.size());
			for (int j = 0; j < _vecStates.size(); j++)
				res[j] = _vecStates[j].CentreOfMassVel().Z();
			return res;
		}

		Vector<Real> getTotalKineticEnergy() const
		{
			Vector<Real> res(_vecStates.size());
			for (int j = 0; j < _vecStates.size(); j++)
				res[j] = _vecStates[j].TotalKineticEnergy();
			return res;
		}
		Vector<Real> getTotalPotentialEnergy() const
		{
			Vector<Real> res(_vecStates.size());
			for (int j = 0; j < _vecStates.size(); j++)
				res[j] = _vecStates[j].TotalPotentialEnergy();
			return res;
		}
		Vector<Real> getTotalEnergy() const
		{
			Vector<Real> res(_vecStates.size());
			for (int j = 0; j < _vecStates.size(); j++)
				res[j] = _vecStates[j].TotalKineticEnergy() + _vecStates[j].TotalPotentialEnergy();
			return res;
		}

		Vector<Real> getLinearMomentumX() const
		{
			Vector<Real> res(_vecStates.size());
			for (int j = 0; j < _vecStates.size(); j++)
				res[j] = _vecStates[j].LinearMomentum().X();
			return res;
		}
		Vector<Real> getLinearMomentumY() const
		{
			Vector<Real> res(_vecStates.size());
			for (int j = 0; j < _vecStates.size(); j++)
				res[j] = _vecStates[j].LinearMomentum().Y();
			return res;
		}
		Vector<Real> getLinearMomentumZ() const
		{
			Vector<Real> res(_vecStates.size());
			for (int j = 0; j < _vecStates.size(); j++)
				res[j] = _vecStates[j].LinearMomentum().Z();
			return res;
		}

		Vector<Real> getAngularMomentumCMX() const
		{
			Vector<Real> res(_vecStates.size());
			for (int j = 0; j < _vecStates.size(); j++)
				res[j] = _vecStates[j].AngularMomentumCM().X();
			return res;
		}
		Vector<Real> getAngularMomentumCMY() const
		{
			Vector<Real> res(_vecStates.size());
			for (int j = 0; j < _vecStates.size(); j++)
				res[j] = _vecStates[j].AngularMomentumCM().Y();
			return res;
		}
		Vector<Real> getAngularMomentumCMZ() const
		{
			Vector<Real> res(_vecStates.size());
			for (int j = 0; j < _vecStates.size(); j++)
				res[j] = _vecStates[j].AngularMomentumCM().Z();
			return res;
		}

		Vector<Real> getAngularMomentumX(Vec3Cart origin) const
		{
			Vector<Real> res(_vecStates.size());
			for (int i = 0; i < _vecStates.size(); i++)
				res[i] = _vecStates[i].AngularMomentum(origin).X();
			return res;
		}
		Vector<Real> getAngularMomentumY(Vec3Cart origin) const
		{
			Vector<Real> res(_vecStates.size());
			for (int i = 0; i < _vecStates.size(); i++)
				res[i] = _vecStates[i].AngularMomentum(origin).Y();
			return res;
		}
		Vector<Real> getAngularMomentumZ(Vec3Cart origin) const
		{
			Vector<Real> res(_vecStates.size());
			for (int i = 0; i < _vecStates.size(); i++)
				res[i] = _vecStates[i].AngularMomentum(origin).Z();
			return res;
		}

		void VisualizeAsParamCurve(std::string baseFileName, Vector<int> vecBodiesIndexToVisualize)
		{
			Vector<Real> t_vals = getTimes();

			std::vector<std::string> fileNames;

			for (int i : vecBodiesIndexToVisualize)
			{
				Vector<Real> body_x_vals = getPosX(i);
				Vector<Real> body_y_vals = getPosY(i);
				Vector<Real> body_z_vals = getPosZ(i);

				// form ParametricCurve from these 3 vectors
				std::vector<VectorN<Real, 3>> res;
				for (int i = 0; i < t_vals.size(); i++)
					res.push_back(VectorN<Real, 3>{body_x_vals[i], body_y_vals[i], body_z_vals[i]});

				std::string fileName = baseFileName + std::to_string(i) + ".txt";
				fileNames.push_back(fileName);

				Real t1 = t_vals[0];
				Real t2 = t_vals[t_vals.size() - 1];

				Serializer::SaveAsParamCurve<3>(res, "PARAMETRIC_CURVE_CARTESIAN_3D", baseFileName + std::to_string(i),
																				t1, t2, t_vals.size(), GetResultFilesPath() + fileName);
			}

			Visualizer::VisualizeMultiParamCurve3D(fileNames);
		}

		void VisualizeAsParticleSimulation(std::string baseFileName, Vector<int> vecBodiesIndexToVisualize, double dT)
		{
			std::vector<std::string> vecColors(NumBodies());
			std::vector<Real> vecRad(NumBodies());
			for (int i = 0; i < NumBodies(); i++)
			{
				vecColors[i] = _config.Color(i);
				vecRad[i] = _config.Radius(i);
			}

			std::vector<std::vector<Pnt3Cart>> res2;
			// copy elements
			res2.resize(NumBodies());
			for (int i = 0; i < NumBodies(); i++)
			{
				res2[i].resize(NumSteps());
				for (int j = 0; j < NumSteps(); j++)
					res2[i][j] = Pnt3Cart(_vecStates[j].Pos(i).X(), _vecStates[j].Pos(i).Y(), _vecStates[j].Pos(i).Z());
			}

			Serializer::SaveParticleSimulation3D(GetResultFilesPath() + baseFileName, NumBodies(), -1, -1, -1, res2,
																					 vecColors, vecRad, dT);

			Visualizer::VisualizeParticleSimulation3D(baseFileName);
		}
	};

	/*********************************************************************************************/
	class NBodyGravitySystemODE : public IODESystem
	{
		NBodyGravitySimConfig _config;

	public:
		NBodyGravitySystemODE(NBodyGravitySimConfig inConfig) : _config(inConfig) {}

		int getDim() const { return 6 * _config.NumBodies(); }
		void derivs(const Real t, const Vector<Real>& x, Vector<Real>& dxdt) const
		{
			// filling in dxdt vector with N * 6 derivations of our variables (x, y, z, vx, vy, vz)
			for (int i = 0; i < _config.NumBodies(); i++)
			{
				Vec3Cart force(0, 0, 0);				// calculating force on body i
				for (int j = 0; j < _config.NumBodies(); j++)
				{
					if (i != j)      // checking for self-force 
					{
						Vec3Cart vec_dist(x[3 * j] - x[3 * i],
															x[3 * j + 1] - x[3 * i + 1],
															x[3 * j + 2] - x[3 * i + 2]);

						force = force + _config.G() * _config.Mass(i) * _config.Mass(j) / POW3(vec_dist.NormL2()) * vec_dist;
					}
				}

				// x, y, z coord
				dxdt[3 * i]			= x[3 * 2 * i + 1];
				dxdt[3 * i + 1] = x[6 * i + 2];
				dxdt[3 * i + 2] = x[6 * i + 3];

				// vx, vy, vz velocity components
				int half = _config.NumBodies() * 3;
				dxdt[half + 3 * i]		 = 1 / _config.Mass(i) * force.X();
				dxdt[half + 3 * i + 1] = 1 / _config.Mass(i) * force.Y();
				dxdt[half + 3 * i + 2] = 1 / _config.Mass(i) * force.Z();
			}
		}
	};

	/*********************************************************************************************/
	class NBodyGravitySimulator
	{
	protected:
		NBodyGravitySimConfig _config;

	public:
		NBodyGravitySimulator(const NBodyGravitySimConfig& inConfig) : _config(inConfig) {}

		NBodyGravitySimulationResults SolveEuler(Real dT, int numSteps)
		{
			NBodyGravitySimulationResults results(_config, dT * numSteps);

			results._vecTimes.Resize(numSteps);
			results._vecStates.Resize(numSteps);

			// getting initil system state from given configuraton
			NBodyState sysState = _config.InitState();

			// save initial positions
			results._vecTimes[0] = 0.0;
			for (int i = 0; i < sysState.NumBodies(); i++)
				results._vecStates[0].addBody(sysState.Mass(i), sysState.Pos(i), sysState.Vel(i));

			// performing N Euler steps
			for (int i = 1; i < numSteps; i++)
			{
				// save time
				results._vecTimes[i] = results._vecTimes[i - 1] + dT;

				// perform one Euler step
				Vector<Vec3Cart> force(sysState.NumBodies(), Vec3Cart(0.0, 0.0, 0.0));
				for (int j = 0; j < sysState.NumBodies()-1; j++) {
					for (int k = j + 1; k < sysState.NumBodies(); k++) {
						Vec3Cart radialVec = sysState.Pos(j) - sysState.Pos(k);
						force[j] = force[j] - _config.G() * sysState.Mass(j) * sysState.Mass(k) / POW3(radialVec.NormL2()) * radialVec;
						force[k] = force[k] + _config.G() * sysState.Mass(j) * sysState.Mass(k) / POW3(radialVec.NormL2()) * radialVec; // opposite force on body j
					}
				}

				// advancing velocities and positions,  using the most simple Euler method of first order
				for (int j = 0; j < _config.NumBodies(); j++) {
					sysState.SetVelocity(j, sysState.Vel(j) + force[j] * dT / sysState.Mass(j));
					sysState.SetPosition(j, sysState.Pos(j) + sysState.Vel(j) * dT);
				}

				// save state
				for (int j = 0; j < _config.NumBodies(); j++)
				{
					results._vecStates[i].addBody(sysState.Mass(j), sysState.Pos(j), sysState.Vel(j));
				}
			}

			return results;
		}

		NBodyGravitySimulationResults SolveRK5(Real duration, Real eps, Real minSaveInterval, Real hStart)
		{
			NBodyGravitySystemODE ode(_config);

			ODESystemSolver<RK5_CashKarp_Stepper> rk5solver(ode);
			ODESystemSolution sol = rk5solver.integrate(_config.getInitCond(), 0, duration, minSaveInterval, eps, hStart);

			NBodyGravitySimulationResults results(_config, duration);

			results._duration = duration;
			results._vecTimes = sol.getTValues();

			int numBodies = _config.NumBodies();
			int numSteps = sol.getTValues().size();
			results._vecStates.Resize(numSteps);
			for (int i = 0; i < numSteps; i++)
			{
				NBodyState state;
				for (int j = 0; j < numBodies; j++)
				{
					Vec3Cart pos(sol.getXValues()[3 * j][i], sol.getXValues()[3 * j + 1][i], sol.getXValues()[3 * j + 2][i]);
					Vec3Cart vel(sol.getXValues()[numBodies * 3 + 3 * j][i], sol.getXValues()[numBodies * 3 + 3 * j + 1][i], sol.getXValues()[numBodies * 3 + 3 * j + 2][i]);
					state.addBody(_config.Mass(j), pos, vel);
				}
				results._vecStates[i] = state;
			}

			return results;
		}
	};

	// SolarSystem gravity simulator
	class SolarSystemGravitySimulator
	{
		// ima definiranu putanju za sve planete
		ParametricCurve<3> _planetPaths[9];
	};
} // namespace MPL

#endif
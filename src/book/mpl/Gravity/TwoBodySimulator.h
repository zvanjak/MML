#if !defined MPL_TWO_BODY_SIMULATOR_H
#define MPL_TWO_BODY_SIMULATOR_H

#include "MMLBase.h"

#include "interfaces/IODESystem.h"

#include "base/VectorTypes.h"
#include "base/Geometry3D.h"

#include "algorithms/ODESystemSolver.h"
#include "algorithms/ODESystemStepCalculators.h"
#include "algorithms/ODESystemSteppers.h"

#include "GravityBase.h"
#include "TwoBodyCalculator.h"

using namespace MML;

namespace MPL
{
	/*********************************************************************************************/
	class TwoBodyState
	{
	private:
		GravityBodyState _body1;
		GravityBodyState _body2;
	
	public:
		TwoBodyState() { }
		TwoBodyState(const GravityBodyState& m1, const GravityBodyState& m2)
			: _body1(m1), _body2(m2) { }

		const GravityBodyState& Body1() const { return _body1; }
		const GravityBodyState& Body2() const { return _body2; }

		GravityBodyState& BodyAcc1() { return _body1; }
		GravityBodyState& BodyAcc2() { return _body2; }

		Vec3Cart Pos1() const { return _body1.R(); }
		Vec3Cart Pos2() const { return _body2.R(); }
		Vec3Cart V1()		const { return _body1.V(); }
		Vec3Cart V2()		const { return _body2.V(); }

		Real Dist() const
		{
			return (_body2.R() - _body1.R()).NormL2();
		}

		Vec3Cart CenterOfMassPos() const
		{
			return (_body1.Mass() * _body1.R() + _body2.Mass() * _body2.R()) / (_body1.Mass() + _body2.Mass());
		}
		Vec3Cart CenterOfMassVelocity() const
		{
			return (_body1.Mass() * _body1.V() + _body2.Mass() * _body2.V()) / (_body1.Mass() + _body2.Mass());
		}

		Real ReducedMass() const { return _body1.Mass() * _body2.Mass() / (_body1.Mass() + _body2.Mass()); }

		Real KineticEnergy1()	const { return 0.5 * _body1.Mass() * POW2(_body1.V().NormL2()); }
		Real KineticEnergy2()	const { return 0.5 * _body2.Mass() * POW2(_body2.V().NormL2()); }
		
		Real TotalKineticEnergy() const { return KineticEnergy1() + KineticEnergy2(); }
		Real TotalPotentialEnergy(Real G) const
		{
			return -G * _body1.Mass() * _body2.Mass() / Dist();
		}

		Vec3Cart TotalMomentum() const
		{
			return _body1.Mass() * _body1.V() + _body2.Mass() * _body2.V();
		}
		Vec3Cart AngularMomentumCM() const
		{
			Vec3Cart r1 = _body1.R() - CenterOfMassPos();
			Vec3Cart r2 = _body2.R() - CenterOfMassPos();
			Vec3Cart L1 = VectorProduct(r1, _body1.V()) * _body1.Mass();
			Vec3Cart L2 = VectorProduct(r2, _body2.V()) * _body2.Mass();
			return L1 + L2;
		}
		Vec3Cart AngularMomentum(Vec3Cart origin) const
		{
			Vec3Cart r1 = _body1.R() - origin;
			Vec3Cart r2 = _body2.R() - origin;
			Vec3Cart L1 = VectorProduct(r1, _body1.V()) * _body1.Mass();
			Vec3Cart L2 = VectorProduct(r2, _body2.V()) * _body2.Mass();
			return L1 + L2;
		}
	};

	/*********************************************************************************************/
	class TwoBodyGravitySimConfig
	{
	private:
		Real _G = 100;
		TwoBodyState _initState;

	public:
		TwoBodyGravitySimConfig(const GravityBodyState& m1, const GravityBodyState& m2)
			: _initState(m1, m2) { }
		TwoBodyGravitySimConfig(const Real& G, const GravityBodyState& m1, const GravityBodyState& m2)
			: _G(G), _initState(m1, m2)	{	}

		Real G() const { return _G; }

		const TwoBodyState& InitState() const { return _initState; }
		TwoBodyState& InitStateAcc() { return _initState; }

		Real Mass1() const { return _initState.Body1().Mass(); }
		Real Mass2() const { return _initState.Body2().Mass(); }

		Vector<Real> getInitCond()
		{
			return Vector<Real>{_initState.Pos1()[0], _initState.Pos1()[1], _initState.Pos1()[2],
													_initState.Pos2()[0], _initState.Pos2()[1], _initState.Pos2()[2],
													_initState.V1()[0], _initState.V1()[1], _initState.V1()[2],
													_initState.V2()[0], _initState.V2()[1], _initState.V2()[2] };
		}
	};

	/*********************************************************************************************/
	class TwoBodyGravityConfigGenerator
	{
	public:
		static TwoBodyGravitySimConfig Config1_same_bodies_elliptic_CM_static()
		{
			Real G = 1;
			GravityBodyState m1(10000, Vec3Cart{ -100, 100, 50 }, Vec3Cart{ 4, 0, -2 });
			GravityBodyState m2(10000, Vec3Cart{ 50, -100, 50 }, Vec3Cart{ -4, 0, 2 });

			return TwoBodyGravitySimConfig(G, m1, m2);
		}
		static TwoBodyGravitySimConfig Config2_same_bodies_elliptic_CM_moving()
		{
			Real G = 1;
			GravityBodyState m1(10000, Vec3Cart{ -100, 100, 50 }, Vec3Cart{ 4, 3, -2 });
			GravityBodyState m2(10000, Vec3Cart{ 50, -100, 50 }, Vec3Cart{ -2, -1, 4 });

			return TwoBodyGravitySimConfig(G, m1, m2);
		}
	};

	/*********************************************************************************************/
	struct TwoBodyGravitySimulationResults
	{
		const TwoBodyGravitySimConfig& _config;

		Real _duration;
		Vector<TwoBodyState> _vecStates;
		Vector<Real> _vecTimes;

	public:
		TwoBodyGravitySimulationResults(const TwoBodyGravitySimConfig& config)
			: _config(config), _duration(0)
		{
		}
		TwoBodyGravitySimulationResults(const TwoBodyGravitySimConfig& config, const Real& duration)
			: _config(config), _duration(duration)
		{
		}

		const TwoBodyGravitySimConfig& Config() const { return _config; }
		Real Duration() const { return _duration; }

		Real				 Time(int i)  const { return _vecTimes[i]; }
		TwoBodyState State(int i) const { return _vecStates[i]; }

		Vector<Real>		 getTimes() const { return _vecTimes; }
		Vector<Vec3Cart> getBody1Pos() const
		{
			Vector<Vec3Cart> res(_vecStates.size());
			for (int i = 0; i < _vecStates.size(); i++)
				res[i] = _vecStates[i].Body1().R();
			return res;
		}
		Vector<Vec3Cart> getBody2Pos() const
		{
			Vector<Vec3Cart> res(_vecStates.size());
			for (int i = 0; i < _vecStates.size(); i++)
				res[i] = _vecStates[i].Body2().R();
			return res;
		}
		Vector<Vec3Cart> getBody1Velocity() const
		{
			Vector<Vec3Cart> res(_vecStates.size());
			for (int i = 0; i < _vecStates.size(); i++)
				res[i] = _vecStates[i].Body1().V();
			return res;
		}
		Vector<Vec3Cart> getBody2Velocity() const
		{
			Vector<Vec3Cart> res(_vecStates.size());
			for (int i = 0; i < _vecStates.size(); i++)
				res[i] = _vecStates[i].Body2().V();
			return res;
		}

		Vector<Real> getPos1X() const
		{
			Vector<Real> res(_vecStates.size());
			for (int i = 0; i < _vecStates.size(); i++)
				res[i] = _vecStates[i].Body1().R().X();
			return res;
		}
		Vector<Real> getPos1Y() const
		{
			Vector<Real> res(_vecStates.size());
			for (int i = 0; i < _vecStates.size(); i++)
				res[i] = _vecStates[i].Body1().R().Y();
			return res;
		}
		Vector<Real> getPos1Z() const
		{
			Vector<Real> res(_vecStates.size());
			for (int i = 0; i < _vecStates.size(); i++)
				res[i] = _vecStates[i].Body1().R().Z();
			return res;
		}
		Vector<Real> getPos2X() const
		{
			Vector<Real> res(_vecStates.size());
			for (int i = 0; i < _vecStates.size(); i++)
				res[i] = _vecStates[i].Body2().R().X();
			return res;
		}
		Vector<Real> getPos2Y() const
		{
			Vector<Real> res(_vecStates.size());
			for (int i = 0; i < _vecStates.size(); i++)
				res[i] = _vecStates[i].Body2().R().Y();
			return res;
		}
		Vector<Real> getPos2Z() const
		{
			Vector<Real> res(_vecStates.size());
			for (int i = 0; i < _vecStates.size(); i++)
				res[i] = _vecStates[i].Body2().R().Z();
			return res;
		}

		Vector<Real> getV1X() const
		{
			Vector<Real> res(_vecStates.size());
			for (int i = 0; i < _vecStates.size(); i++)
				res[i] = _vecStates[i].Body1().V().X();
			return res;
		}
		Vector<Real> getV1Y() const
		{
			Vector<Real> res(_vecStates.size());
			for (int i = 0; i < _vecStates.size(); i++)
				res[i] = _vecStates[i].Body1().V().Y();
			return res;
		}
		Vector<Real> getV1Z() const
		{
			Vector<Real> res(_vecStates.size());
			for (int i = 0; i < _vecStates.size(); i++)
				res[i] = _vecStates[i].Body1().V().Z();
			return res;
		}
		Vector<Real> getV2X() const
		{
			Vector<Real> res(_vecStates.size());
			for (int i = 0; i < _vecStates.size(); i++)
				res[i] = _vecStates[i].Body2().V().X();
			return res;
		}
		Vector<Real> getV2Y() const
		{
			Vector<Real> res(_vecStates.size());
			for (int i = 0; i < _vecStates.size(); i++)
				res[i] = _vecStates[i].Body2().V().Y();
			return res;
		}
		Vector<Real> getV2Z() const
		{
			Vector<Real> res(_vecStates.size());
			for (int i = 0; i < _vecStates.size(); i++)
				res[i] = _vecStates[i].Body2().V().Z();
			return res;
		}

		Vector<Real> getCMPosX() const
		{
			Vector<Real> res(_vecStates.size());
			for (int i = 0; i < _vecStates.size(); i++)
				res[i] = _vecStates[i].CenterOfMassPos().X();
			return res;
		}
		Vector<Real> getCMPosY() const
		{
			Vector<Real> res(_vecStates.size());
			for (int i = 0; i < _vecStates.size(); i++)
				res[i] = _vecStates[i].CenterOfMassPos().Y();
			return res;
		}
		Vector<Real> getCMPosZ() const
		{
			Vector<Real> res(_vecStates.size());
			for (int i = 0; i < _vecStates.size(); i++)
				res[i] = _vecStates[i].CenterOfMassPos().Z();
			return res;
		}

		Vector<Real> getCMVX() const
		{
			Vector<Real> res(_vecStates.size());
			for (int i = 0; i < _vecStates.size(); i++)
				res[i] = _vecStates[i].CenterOfMassVelocity().X();
			return res;
		}
		Vector<Real> getCMVY() const
		{
			Vector<Real> res(_vecStates.size());
			for (int i = 0; i < _vecStates.size(); i++)
				res[i] = _vecStates[i].CenterOfMassVelocity().Y();
			return res;
		}
		Vector<Real> getCMVZ() const
		{
			Vector<Real> res(_vecStates.size());
			for (int i = 0; i < _vecStates.size(); i++)
				res[i] = _vecStates[i].CenterOfMassVelocity().Z();
			return res;
		}

		Vector<Real> getTotalKineticEnergy() const
		{
			Vector<Real> res(_vecStates.size());
			for (int i = 0; i < _vecStates.size(); i++)
				res[i] = _vecStates[i].TotalKineticEnergy();
			return res;
		}
		Vector<Real> getTotalPotentialEnergy() const
		{
			Vector<Real> res(_vecStates.size());
			for (int i = 0; i < _vecStates.size(); i++)
				res[i] = _vecStates[i].TotalPotentialEnergy(_config.G());
			return res;
		}
		Vector<Real> getTotalEnergy() const
		{
			Vector<Real> res(_vecStates.size());
			for (int i = 0; i < _vecStates.size(); i++)
				res[i] = _vecStates[i].TotalKineticEnergy() + _vecStates[i].TotalPotentialEnergy(_config.G());
			return res;
		}

		Vector<Real> getTotalMomentumX() const
		{
			Vector<Real> res(_vecStates.size());
			for (int i = 0; i < _vecStates.size(); i++)
				res[i] = _vecStates[i].TotalMomentum().X();
			return res;
		}
		Vector<Real> getTotalMomentumY() const
		{
			Vector<Real> res(_vecStates.size());
			for (int i = 0; i < _vecStates.size(); i++)
				res[i] = _vecStates[i].TotalMomentum().Y();
			return res;
		}
		Vector<Real> getTotalMomentumZ() const
		{
			Vector<Real> res(_vecStates.size());
			for (int i = 0; i < _vecStates.size(); i++)
				res[i] = _vecStates[i].TotalMomentum().Z();
			return res;
		}

		Vector<Real> getAngularMomentumCMX() const
		{
			Vector<Real> res(_vecStates.size());
			for (int i = 0; i < _vecStates.size(); i++)
				res[i] = _vecStates[i].AngularMomentumCM().X();
			return res;
		}
		Vector<Real> getAngularMomentumCMY() const
		{
			Vector<Real> res(_vecStates.size());
			for (int i = 0; i < _vecStates.size(); i++)
				res[i] = _vecStates[i].AngularMomentumCM().Y();
			return res;
		}
		Vector<Real> getAngularMomentumCMZ() const
		{
			Vector<Real> res(_vecStates.size());
			for (int i = 0; i < _vecStates.size(); i++)
				res[i] = _vecStates[i].AngularMomentumCM().Z();
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

		// get parametric curve

		// visualize as parametric curve
		void VisualizeSolution(std::string baseFileName)
		{
			// first body
			std::vector<VectorN<Real, 3>> res1;
			for (int i = 0; i < _vecStates.size(); i++)
				res1.push_back(VectorN<Real, 3>{_vecStates[i].Body1().R().X(), _vecStates[i].Body1().R().Y(), _vecStates[i].Body1().R().Z()});

			Serializer::SaveAsParamCurve<3>(res1, "PARAMETRIC_CURVE_CARTESIAN_3D", baseFileName + "_body1",
				0, _duration, _vecTimes.size(),
				GetResultFilesPath() + baseFileName + "_body1.txt");

			// second body
			std::vector<VectorN<Real, 3>> res2;
			for (int i = 0; i < _vecStates.size(); i++)
				res2.push_back(VectorN<Real, 3>{_vecStates[i].Body2().R().X(), _vecStates[i].Body2().R().Y(), _vecStates[i].Body2().R().Z()});

			Serializer::SaveAsParamCurve<3>(res2, "PARAMETRIC_CURVE_CARTESIAN_3D", baseFileName + "_body2",
				0, _duration, _vecTimes.size(),
				GetResultFilesPath() + baseFileName + "_body2.txt");

			Visualizer::VisualizeMultiParamCurve3D({ baseFileName + "_body1.txt", baseFileName + "_body2.txt" });
		}
	};

	/*********************************************************************************************/
	class TwoBodyGravitySystemODE : public IODESystem
	{
		TwoBodyGravitySimConfig _config;

	public:
		TwoBodyGravitySystemODE(const TwoBodyGravitySimConfig& config)
			: _config(config)	{	}

		int getDim() const override { return 12; }
		void derivs(const Real t, const Vector<Real>& x, Vector<Real>& dxdt) const override
		{
			Vec3Cart r1{ x[0], x[1], x[2] };
			Vec3Cart r2{ x[3], x[4], x[5] };
			Vec3Cart v1{ x[6], x[7], x[8] };
			Vec3Cart v2{ x[9], x[10], x[11] };

			Vec3Cart r12 = r2 - r1;
			Real dist = r12.NormL2();

			Vec3Cart f12 = _config.G() * _config.Mass1() * _config.Mass2() / POW3(dist) * r12;
			Vec3Cart f21 = -f12;

			dxdt[0] = v1[0];
			dxdt[1] = v1[1];
			dxdt[2] = v1[2];
			dxdt[3] = v2[0];
			dxdt[4] = v2[1];
			dxdt[5] = v2[2];
			dxdt[6] = f12[0] / _config.Mass1();
			dxdt[7] = f12[1] / _config.Mass1();
			dxdt[8] = f12[2] / _config.Mass1();
			dxdt[9] = f21[0] / _config.Mass2();
			dxdt[10] = f21[1] / _config.Mass2();
			dxdt[11] = f21[2] / _config.Mass2();
		}
	};

	/*********************************************************************************************/
	class TwoBodiesGravitySimulator
	{
		TwoBodyGravitySimConfig _config;

	public:
		TwoBodiesGravitySimulator(const TwoBodyGravitySimConfig& config)
			: _config(config)	{	}

		TwoBodyGravitySimulationResults SolveRK5(Real duration, Real eps, Real minSaveInterval, Real hStart)
		{
			TwoBodyGravitySystemODE ode(_config);

			ODESystemSolver<RK5_CashKarp_Stepper> rk5solver(ode);
			ODESystemSolution sol = rk5solver.integrate(_config.getInitCond(), 0, duration, minSaveInterval, eps, hStart);

			TwoBodyGravitySimulationResults results(_config);

			results._vecTimes = sol.getTValues();
			results._duration = duration;

			int numSteps = sol.getTValues().size();
			results._vecStates.Resize(numSteps);
			for (int i = 0; i < numSteps; i++)
			{
				results._vecStates[i].BodyAcc1() = GravityBodyState(_config.Mass1(),
					Vec3Cart{ sol.getXValues(0)[i], sol.getXValues(1)[i], sol.getXValues(2)[i] },
					Vec3Cart{ sol.getXValues(6)[i], sol.getXValues(7)[i], sol.getXValues(8)[i] });
				results._vecStates[i].BodyAcc2() = GravityBodyState(_config.Mass2(),
					Vec3Cart{ sol.getXValues(3)[i], sol.getXValues(4)[i], sol.getXValues(5)[i] },
					Vec3Cart{ sol.getXValues(9)[i], sol.getXValues(10)[i], sol.getXValues(11)[i] });
			}

			return results;
		}
		TwoBodyGravitySimulationResults SolveRK5(Real duration, Real eps, Real minSaveInterval)
		{
			return SolveRK5(duration, eps, minSaveInterval, 0.01);
		}
		TwoBodyGravitySimulationResults SolveRK5(Real duration, Real eps)
		{
			return SolveRK5(duration, eps, 0.01);
		}
		TwoBodyGravitySimulationResults SolveRK5(Real duration)
		{
			return SolveRK5(duration, 1e-06);
		}

		// get trajectory type - Ellipse, Parabolic, Hyperbolic
		std::string GetTrajectoryType(const TwoBodyGravitySimConfig &sysConfig)
		{
			return TwoBodyGravityCalculator::GetTrajectoryType(sysConfig.G(), sysConfig.InitState().Body1(), sysConfig.InitState().Body2());
		}
		
		// get plane of motion
		Plane3D GetPlaneOfMotion(const TwoBodyGravitySimConfig &sysConfig)
		{
			return TwoBodyGravityCalculator::GetPlaneOfMotion(sysConfig.G(), sysConfig.InitState().Body1(), sysConfig.InitState().Body2());
		}

		// solve and visualize trajectory
		void SolveAndShowTrajectories(Real t)
		{
			TwoBodyGravitySystemODE ode(_config);

			ODESystemSolver<RK5_CashKarp_Stepper> adaptSolver(ode);
			ODESystemSolution solAdapt = adaptSolver.integrate(_config.getInitCond(), 0, t, 0.01, 1e-06, 0.01);

			Vector<Real> t_vals = solAdapt.getTValues();

			Vector<Real> r1_x_vals = solAdapt.getXValues(0);
			Vector<Real> r1_y_vals = solAdapt.getXValues(1);
			Vector<Real> r1_z_vals = solAdapt.getXValues(2);
			Vector<Real> r2_x_vals = solAdapt.getXValues(3);
			Vector<Real> r2_y_vals = solAdapt.getXValues(4);
			Vector<Real> r2_z_vals = solAdapt.getXValues(5);

			// form ParametricCurve from these 3 vectors
			std::vector<VectorN<Real, 3>> res1;
			for (int i = 0; i < t_vals.size(); i++)
				res1.push_back(VectorN<Real, 3>{r1_x_vals[i], r1_y_vals[i], r1_z_vals[i]});
			Serializer::SaveAsParamCurve<3>(res1, "PARAMETRIC_CURVE_CARTESIAN_3D", "gravity_example_body1",
				0, t, t_vals.size(),
				GetResultFilesPath() + "gravity_example_body1.txt");

			std::vector<VectorN<Real, 3>> res2;
			for (int i = 0; i < t_vals.size(); i++)
				res2.push_back(VectorN<Real, 3>{r2_x_vals[i], r2_y_vals[i], r2_z_vals[i]});
			Serializer::SaveAsParamCurve<3>(res2, "PARAMETRIC_CURVE_CARTESIAN_3D", "gravity_example_body2",
				0, t, t_vals.size(),
				GetResultFilesPath() + "gravity_example_body2.txt");

			Visualizer::VisualizeMultiParamCurve3D({ "gravity_example_body1.txt", "gravity_example_body2.txt" });
		}
	};
} // namespace MPL

#endif
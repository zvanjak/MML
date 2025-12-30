#if !defined MPL_TWO_BODY_SIMULATOR_2D_H
#define MPL_TWO_BODY_SIMULATOR_2D_H

#include "MMLBase.h"

#include "interfaces/IODESystem.h"

#include "base/VectorTypes.h"
#include "base/Geometry2D.h"

#include "algorithms/ODESystemSolver.h"
#include "algorithms/ODESystemStepCalculators.h"
#include "algorithms/ODESystemSteppers.h"

#include "GravityBase.h"
#include "TwoBodyCalculator.h"

using namespace MML;

namespace MPL
{
	class TwoBodyState2D
	{
	private:
		GravityBodyState2D _body1;
		GravityBodyState2D _body2;
	
	public:
		TwoBodyState2D() { }
		TwoBodyState2D(const GravityBodyState2D& m1, const GravityBodyState2D& m2)
			: _body1(m1), _body2(m2) { }

		const GravityBodyState2D& Body1() const { return _body1; }
		const GravityBodyState2D& Body2() const { return _body2; }

		GravityBodyState2D& BodyAcc1() { return _body1; }
		GravityBodyState2D& BodyAcc2() { return _body2; }

		Vec2Cart Pos1() const { return _body1.R(); }
		Vec2Cart Pos2() const { return _body2.R(); }
		Vec2Cart V1()		const { return _body1.V(); }
		Vec2Cart V2()		const { return _body2.V(); }

		Real Dist() const
		{
			return (_body2.R() - _body1.R()).NormL2();
		}

		Vec2Cart CenterOfMassPos() const
		{
			return (_body1.Mass() * _body1.R() + _body2.Mass() * _body2.R()) / (_body1.Mass() + _body2.Mass());
		}
		Vec2Cart CenterOfMassVelocity() const
		{
			return (_body1.Mass() * _body1.V() + _body2.Mass() * _body2.V()) / (_body1.Mass() + _body2.Mass());
		}

		Real ReducedMass() const { return _body1.Mass() * _body2.Mass() / (_body1.Mass() + _body2.Mass()); }

		Real KineticEnergy1()			const { return 0.5 * _body1.Mass() * POW2(_body1.V().NormL2()); }
		Real KineticEnergy2()			const { return 0.5 * _body2.Mass() * POW2(_body2.V().NormL2()); }
		Real TotalKineticEnergy() const { return KineticEnergy1() + KineticEnergy2(); }

		Real TotalPotentialEnergy(Real G) const
		{
			return -G * _body1.Mass() * _body2.Mass() / Dist();
		}
	};

	class TwoBodyGravitySimConfig2D
	{
	private:
		Real _G = 100;
		TwoBodyState2D _initState;

	public:
		TwoBodyGravitySimConfig2D(const GravityBodyState2D& m1, const GravityBodyState2D& m2)
			: _initState(m1, m2) { }
		TwoBodyGravitySimConfig2D(const Real& G, const GravityBodyState2D& m1, const GravityBodyState2D& m2)
			: _G(G), _initState(m1, m2)	{	}

		Real G() const { return _G; }

		const TwoBodyState2D& InitState() const { return _initState; }
		TwoBodyState2D& InitStateAcc() { return _initState; }

		Real Mass1() const { return _initState.Body1().Mass(); }
		Real Mass2() const { return _initState.Body2().Mass(); }

		Vector<Real> getInitCond()
		{
			return Vector<Real>{_initState.Pos1()[0], _initState.Pos1()[1],
													_initState.Pos2()[0], _initState.Pos2()[1],
													_initState.V1()[0], _initState.V1()[1],
													_initState.V2()[0], _initState.V2()[1] };
		}
	};

	class TwoBodyGravityConfigGenerator2D
	{
	public:
		// config for two same bodies, in elliptic orbit, center of mass static
		static TwoBodyGravitySimConfig2D Config1_same_bodies_elliptic_CM_static(Real G = 1)
		{
			GravityBodyState2D body1(10000, Vec2Cart{ -100, 100 }, Vec2Cart{ 4, 0 }, 10, "Black");
			GravityBodyState2D body2(10000, Vec2Cart{ 100, -100 }, Vec2Cart{ -4, 0 }, 10, "Blue");

			TwoBodyGravitySimConfig2D demo1(G, body1, body2);

			return demo1;
		}

		// config for two same bodies, in elliptic orbit, center of mass moving
		static TwoBodyGravitySimConfig2D Config2_same_bodies_elliptic_CM_moving(Real G = 1)
		{
			GravityBodyState2D body1(10000, Vec2Cart{ -100, 100 }, Vec2Cart{ 6, 0 }, 10, "Black");
			GravityBodyState2D body2(10000, Vec2Cart{ 100, -100 }, Vec2Cart{ -4, 0 }, 10, "Blue");

			TwoBodyGravitySimConfig2D demo1(G, body1, body2);

			return demo1;
		}

		// config for two different bodies, in elliptic orbit, center of mass static
		static TwoBodyGravitySimConfig2D Config3_diff_bodies_elliptic_CM_static(Real G = 1)
		{
			GravityBodyState2D body1(10000, Vec2Cart{ -100, 100 }, Vec2Cart{ 2, 0 }, 10, "Black");
			GravityBodyState2D body2(5000, Vec2Cart{ 100, -100 }, Vec2Cart{ -4, 0.0 }, 5, "Blue");

			TwoBodyGravitySimConfig2D demo1(G, body1, body2);

			return demo1;
		}

		// config for two different bodies, in elliptic orbit, center of mass moving
		static TwoBodyGravitySimConfig2D Config4_diff_bodies_elliptic_CM_moving(Real G = 1)
		{
			GravityBodyState2D body1(10000, Vec2Cart{ -100, 100 }, Vec2Cart{ 4, 0 }, 10, "Black");
			GravityBodyState2D body2(5000, Vec2Cart{ 100, -100 }, Vec2Cart{ -4, 0.0 }, 5, "Blue");
			
			TwoBodyGravitySimConfig2D demo1(G, body1, body2);
			
			return demo1;
		}

		// config for two same bodies, in hyperbolic orbit, center of mass static+
		static TwoBodyGravitySimConfig2D Config5_same_bodies_hyperbolic_CM_static(Real G = 1)
		{
			GravityBodyState2D m1(10000, Vec2Cart{ -100, 100 }, Vec2Cart{ 5.9, 0 });
			GravityBodyState2D m2(10000, Vec2Cart{ 100, -100 }, Vec2Cart{ -5.9, 0 });
			TwoBodyGravitySimConfig2D demo1(G, m1, m2);	
			return demo1;
		}

		// config for two same bodies, in hyperbolic orbit, center of mass static+
		static TwoBodyGravitySimConfig2D Config6_same_bodies_hyperbolic_CM_moving(Real G = 1)
		{
			GravityBodyState2D m1(10000, Vec2Cart{ -100, 100 }, Vec2Cart{ 10, 0 });
			GravityBodyState2D m2(10000, Vec2Cart{ 100, -100 }, Vec2Cart{ -6, 0 });
			TwoBodyGravitySimConfig2D demo1(G, m1, m2);
			return demo1;
		}
	};

	class TwoBodyGravitySimulationResults2D
	{
	public:
		const TwoBodyGravitySimConfig2D& _config;

		Real _duration;
		Vector<TwoBodyState2D> _vecStates;
		Vector<Real> _vecTimes;

		TwoBodyGravitySimulationResults2D(const TwoBodyGravitySimConfig2D& config)
			: _config(config), _duration(0)
		{
		}
		TwoBodyGravitySimulationResults2D(const TwoBodyGravitySimConfig2D& config, const Real& duration)
			: _config(config), _duration(duration)
		{
		}

		int   NumSteps() const { return (int)_vecStates.size(); }
		Real	Duration() const { return _duration; }

		Real						Time(int i)  const { return _vecTimes[i]; }
		TwoBodyState2D	State(int i) const { return _vecStates[i]; }

		Vector<Real>		 getTimes() const { return _vecTimes; }
		Vector<Vec2Cart> getBody1Pos() const
		{
			Vector<Vec2Cart> res((int) _vecStates.size());
			for (int i = 0; i < _vecStates.size(); i++)
				res[i] = _vecStates[i].Body1().R();
			return res;
		}
		Vector<Vec2Cart> getBody2Pos() const
		{
			Vector<Vec2Cart> res(_vecStates.size());
			for (int i = 0; i < _vecStates.size(); i++)
				res[i] = _vecStates[i].Body2().R();
			return res;
		}
		Vector<Vec2Cart> getBody1Velocity() const
		{
			Vector<Vec2Cart> res(_vecStates.size());
			for (int i = 0; i < _vecStates.size(); i++)
				res[i] = _vecStates[i].Body1().V();
			return res;
		}
		Vector<Vec2Cart> getBody2Velocity() const
		{
			Vector<Vec2Cart> res(_vecStates.size());
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

		Vector<Real> getTotalMomentumX() const
		{
			Vector<Real> res(_vecStates.size());
			for (int i = 0; i < _vecStates.size(); i++)
				res[i] = _vecStates[i].Body1().Mass() * _vecStates[i].Body1().V().X() +
								 _vecStates[i].Body2().Mass() * _vecStates[i].Body2().V().X();
			return res;
		}
		Vector<Real> getTotalMomentumY() const
		{
			Vector<Real> res(_vecStates.size());
			for (int i = 0; i < _vecStates.size(); i++)
				res[i] = _vecStates[i].Body1().Mass() * _vecStates[i].Body1().V().Y() +
								 _vecStates[i].Body2().Mass() * _vecStates[i].Body2().V().Y();
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

		// visualize as parametric curve
		void VisualizeSolution(std::string baseFileName)
		{
			// first body
			std::vector<VectorN<Real, 2>> res1;
			for (int i = 0; i < _vecStates.size(); i++)
				res1.push_back(VectorN<Real, 2>{_vecStates[i].Body1().R().X(), _vecStates[i].Body1().R().Y()});

			Serializer::SaveAsParamCurve<2>(res1, "PARAMETRIC_CURVE_CARTESIAN_2D", baseFileName + "_body1",
				0, _duration, _vecTimes.size(),
				GetResultFilesPath() + baseFileName + "_body1.txt");

			// second body
			std::vector<VectorN<Real, 2>> res2;
			for (int i = 0; i < _vecStates.size(); i++)
				res2.push_back(VectorN<Real, 2>{_vecStates[i].Body2().R().X(), _vecStates[i].Body2().R().Y()});

			Serializer::SaveAsParamCurve<2>(res2, "PARAMETRIC_CURVE_CARTESIAN_2D", baseFileName + "_body2",
				0, _duration, _vecTimes.size(),
				GetResultFilesPath() + baseFileName + "_body2.txt");

			Visualizer::VisualizeMultiParamCurve2D({ baseFileName + "_body1.txt", baseFileName + "_body2.txt" });
		}
	};

	class TwoBodyGravitySystemODE2D : public IODESystem
	{
		TwoBodyGravitySimConfig2D _config;

	public:
		TwoBodyGravitySystemODE2D(const TwoBodyGravitySimConfig2D& config)
			: _config(config)	{	}

		int getDim() const override { return 8; }
		void derivs(const Real t, const Vector<Real>& x, Vector<Real>& dxdt) const override
		{
			Vec2Cart r1{ x[0], x[1] }, r2{ x[2], x[3] };
			Vec2Cart v1{ x[4], x[5] }, v2{ x[6], x[7] };

			Vec2Cart r12 = r2 - r1;
			Real dist = r12.NormL2();

			Vec2Cart f12 = _config.G() * _config.Mass1() * _config.Mass2() / POW3(dist) * r12;
			Vec2Cart f21 = -f12;

			dxdt[0] = v1[0];
			dxdt[1] = v1[1];
			dxdt[2] = v2[0];
			dxdt[3] = v2[1];
			dxdt[4] = f12[0] / _config.Mass1();
			dxdt[5] = f12[1] / _config.Mass1();
			dxdt[6] = f21[0] / _config.Mass2();
			dxdt[7] = f21[1] / _config.Mass2();
		}
	};

	class TwoBodiesGravitySimulator2D
	{
		TwoBodyGravitySimConfig2D _config;

	public:
		TwoBodiesGravitySimulator2D(const TwoBodyGravitySimConfig2D& config)
			: _config(config)	{	}

		TwoBodyGravitySimulationResults2D SolveRK5(Real duration, Real eps, Real minSaveInterval, Real hStart)
		{
			TwoBodyGravitySystemODE2D ode(_config);

			ODESystemSolver<RK5_CashKarp_Stepper> rk5solver(ode);
			ODESystemSolution sol = rk5solver.integrate(_config.getInitCond(), 0, duration, minSaveInterval, eps, hStart);

			TwoBodyGravitySimulationResults2D results(_config);

			results._vecTimes = sol.getTValues();
			results._duration = duration;

			int numSteps = sol.getTValues().size();
			results._vecStates.Resize(numSteps);
			for (int i = 0; i < numSteps; i++)
			{
				results._vecStates[i].BodyAcc1() = GravityBodyState2D(_config.Mass1(),
					Vec2Cart{ sol.getXValues(0)[i], sol.getXValues(1)[i] },
					Vec2Cart{ sol.getXValues(4)[i], sol.getXValues(5)[i] });
				results._vecStates[i].BodyAcc2() = GravityBodyState2D(_config.Mass2(),
					Vec2Cart{ sol.getXValues(2)[i], sol.getXValues(3)[i] },
					Vec2Cart{ sol.getXValues(6)[i], sol.getXValues(7)[i] });
			}

			return results;
		}
		TwoBodyGravitySimulationResults2D SolveRK5(Real duration, Real eps, Real minSaveInterval)
		{
			return SolveRK5(duration, eps, minSaveInterval, 0.01);
		}
		TwoBodyGravitySimulationResults2D SolveRK5(Real duration, Real eps)
		{
			return SolveRK5(duration, eps, 0.01);
		}
		TwoBodyGravitySimulationResults2D SolveRK5(Real duration)
		{
			return SolveRK5(duration, 1e-06);
		}

		// get trajectory type - Ellipse, Parabolic, Hyperbolic
		std::string GetTrajectoryType(const TwoBodyGravitySimConfig2D &sysConfig)
		{
			return TwoBodyGravityCalculator2D::GetTrajectoryType(sysConfig.G(), sysConfig.InitState().Body1(), sysConfig.InitState().Body2());
		}

		// solve and visualize trajectory
		void SolveAndShowTrajectories(Real t)
		{
			TwoBodyGravitySystemODE2D ode(_config);

			ODESystemSolver<RK5_CashKarp_Stepper> adaptSolver(ode);
			ODESystemSolution solAdapt = adaptSolver.integrate(_config.getInitCond(), 0, t, 0.01, 1e-06, 0.01);

			Vector<Real> t_vals = solAdapt.getTValues();

			Vector<Real> r1_x_vals = solAdapt.getXValues(0);
			Vector<Real> r1_y_vals = solAdapt.getXValues(1);
			Vector<Real> r2_x_vals = solAdapt.getXValues(2);
			Vector<Real> r2_y_vals = solAdapt.getXValues(3);

			// form ParametricCurve from these 2 vectors
			std::vector<VectorN<Real, 2>> res1;
			for (int i = 0; i < t_vals.size(); i++)
				res1.push_back(VectorN<Real, 2>{r1_x_vals[i], r1_y_vals[i]});
			Serializer::SaveAsParamCurve<2>(res1, "PARAMETRIC_CURVE_CARTESIAN_2D", "gravity_example_body1",
				0, t, t_vals.size(),
				GetResultFilesPath() + "gravity_example_body1.txt");

			std::vector<VectorN<Real, 2>> res2;
			for (int i = 0; i < t_vals.size(); i++)
				res2.push_back(VectorN<Real, 2>{r2_x_vals[i], r2_y_vals[i]});
			Serializer::SaveAsParamCurve<2>(res2, "PARAMETRIC_CURVE_CARTESIAN_2D", "gravity_example_body2",
				0, t, t_vals.size(),
				GetResultFilesPath() + "gravity_example_body2.txt");

			Visualizer::VisualizeMultiParamCurve2D({ "gravity_example_body1.txt", "gravity_example_body2.txt" });
		}
	};
} // namespace MPL

#endif
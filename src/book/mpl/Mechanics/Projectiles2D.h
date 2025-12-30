#if !defined MPL_PROJECTILES_2D_H
#define MPL_PROJECTILES_2D_H

#include "MMLBase.h"

#include "interfaces/IODESystem.h"

#include "base/InterpolatedFunction.h"

#include "algorithms/ODESystemSolver.h"
#include "algorithms/ODESystemStepCalculators.h"
#include "algorithms/ODESystemSteppers.h"
using namespace MML;

namespace MPL
{
	class IProjectile2DODE : public IODESystem
	{
	public:
		virtual std::string getVarName(int ind) const override
		{
			switch (ind)
			{
			case 0: return "x";		// x position
			case 1: return "y";		// y position
			case 2: return "vx";	// x velocity
			case 3: return "vy";	// y velocity
			default: return IODESystem::getVarName(ind);
			}
		}
		Vector<Real> getInitCond(Real angle, Real initHeight, Real velocity) const
		{
			Vector<Real> initCond(4);
			initCond[0] = 0;											// initial x position
			initCond[1] = initHeight;							// initial y position
			initCond[2] = velocity * cos(angle);	// initial x velocity (m/s)
			initCond[3] = velocity * sin(angle);	// initial y velocity (m/s)
			return initCond;
		}
	};

	// representing ODE system for projectile motion in vacuum
	class Projectile2DInVacuumODE : public IProjectile2DODE
	{
	public:
		int  getDim() const override { return 4; }
		void derivs(const Real t, const Vector<Real>& x, Vector<Real>& dxdt) const override
		{
			dxdt[0] = x[2];
			dxdt[1] = x[3];
			dxdt[2] = 0;
			dxdt[3] = -9.81;

			if (x[1] < 0) // if the projectile hits the ground, stop the motion
			{
				dxdt[0] = 0; // stop x motion
				dxdt[1] = 0; // stop y motion
				dxdt[2] = 0; // stop x velocity
				dxdt[3] = 0; // stop y velocity

				//std::cout << "Projectile hit the ground at t = " << t << std::endl;
			}
		}

		// calculating range of the projectile motion in vacuum
		Real CalcRange(Real angle, Real initHeight, Real velocity) const
		{
			Real g = 9.81; 
			return (velocity * cos(angle) / g) * (velocity * sin(angle) + sqrt(velocity * velocity * sin(angle) * sin(angle) + 2 * g * initHeight));
		}
		// calculating time of flight for the projectile motion in vacuum
		Real TimeOfFlight(Real angle, Real initHeight, Real velocity) const
		{
			Real g = 9.81;
			Real v0y = velocity * sin(angle);
			return (v0y + sqrt(v0y * v0y + 2 * g * initHeight)) / g;
		}
		// calculating maximum height of the projectile motion in vacuum
		Real MaxHeight(Real angle, Real initHeight, Real velocity) const
		{
			Real g = 9.81;
			Real v0y = velocity * sin(angle);
			return initHeight + (v0y * v0y) / (2 * g);
		}
	};

	class Projectile2DWithAirResistanceODE : public IProjectile2DODE
	{
		Real _dragCoefficient; // air resistance coefficient
	public:
		Projectile2DWithAirResistanceODE(Real dragCoefficient = 0.1) : _dragCoefficient(dragCoefficient) {}

		int  getDim() const override { return 4; }
		void derivs(const Real t, const Vector<Real>& x, Vector<Real>& dxdt) const override
		{
			Real v = sqrt(x[2] * x[2] + x[3] * x[3]);			// speed of the projectile

			dxdt[0] = x[2];
			dxdt[1] = x[3];
			dxdt[2] = -_dragCoefficient * v * x[2];
			dxdt[3] = -9.81 - _dragCoefficient * v * x[3];

			if (x[1] < 0) // if the projectile hits the ground, stop the motion
			{
				dxdt[0] = 0; // stop x motion
				dxdt[1] = 0; // stop y motion
				dxdt[2] = 0; // stop x velocity
				dxdt[3] = 0; // stop y velocity
				//std::cout << "Projectile hit the ground at t = " << t << std::endl;
			}
		}
	};

	// representing ODE system for projectile motion with air resistance, 
	// but with included changes of air density with altitude
	enum class AirDensityModel
	{
		Isothermal,		// isothermal model of air density
		Adiabatic			// adiabatic model of air density
	};

	class Projectile2DChangingAirDensityODE : public IProjectile2DODE
	{
		Real _dragCoefficient0;					// air resistance coefficient at sea level
		Real _airDensity0 = 1.225;			// air density at sea level in kg/m^3

		AirDensityModel _airDensityModel = AirDensityModel::Isothermal; // model of air density

	public:
		Projectile2DChangingAirDensityODE(Real dragCoefficient0 = 0.1) : _dragCoefficient0(dragCoefficient0) {}
		Projectile2DChangingAirDensityODE(Real dragCoefficient0, AirDensityModel airDensityModel) 
			: _dragCoefficient0(dragCoefficient0), _airDensityModel(airDensityModel) {
		}

		int  getDim() const override { return 4; }
		void derivs(const Real t, const Vector<Real>& x, Vector<Real>& dxdt) const override
		{
			Real v = sqrt(x[2] * x[2] + x[3] * x[3]);				// speed of the projectile

			dxdt[0] = x[2];
			dxdt[1] = x[3];
			dxdt[2] = -dragCoefficient(x[3], v) * v * x[2];
			dxdt[3] = -9.81 - dragCoefficient(x[3], v) * v * x[3];

			if (x[1] < 0) // if the projectile hits the ground, stop the motion
			{
				dxdt[0] = 0; // stop x motion
				dxdt[1] = 0; // stop y motion
				dxdt[2] = 0; // stop x velocity
				dxdt[3] = 0; // stop y velocity
			}
		}

		// function to calculate air density based on altitude (isothermal)
		Real airDensity(Real altitude) const
		{
			// Simplified model: air density decreases exponentially with altitude
			// at sea level, air density is approximately 1.225 kg/m^3
			return 1.225 * exp(-altitude / 8000); // 8000 m is a rough scale height for the atmosphere
		}

		// function to calculate air density based on altitude (adiabatic)
		Real airDensityAdiabatic(Real altitude) const
		{
			// Adiabatic model: air density decreases with altitude, but with a different rate
			// at sea level, air density is approximately 1.225 kg/m^3
			return 1.225 * pow((1 - (0.0065 * altitude / 288.15)), 5.2561); // using the adiabatic lapse rate
		}

		// function to calculate drag coefficient based on air density and speed
		Real dragCoefficient(Real altitude, Real speed) const
		{
			if (_airDensityModel == AirDensityModel::Adiabatic)
			{
				Real density = airDensityAdiabatic(altitude);
				return _dragCoefficient0 * density / _airDensity0; // simplified drag coefficient formula
			}

			// Default to isothermal model
			Real density = airDensity(altitude);
			return _dragCoefficient0 * density / _airDensity0; // simplified drag coefficient formula
		}
	};

	enum class BaseballType
	{
		Smooth,	
		Normal,
		Rough,
		NormalApprox
	};
	// representing ODE system for projectile motion with air resistance,
	// but with drag coefficient dependent on speed and type of baseball
	class BaseballWithDragCoeffDependentOnSpeedODE : public IProjectile2DODE
	{
	private:
		BaseballType _baseballType; // type of baseball

		PolynomInterpRealFunc *_dragCoefficientSmooth;			// drag coefficient for soft baseball
		PolynomInterpRealFunc *_dragCoefficientNormal;			// drag coefficient for normal baseball
		PolynomInterpRealFunc *_dragCoefficientRough;				// drag coefficient for rough baseball

	public:
		BaseballWithDragCoeffDependentOnSpeedODE(BaseballType baseballType) : _baseballType(baseballType) 
		{
			// Initialize drag coefficients for different types of baseballs
			Vector<Real> smoothXVal{ 0.0, 50.0, 100.0, 150.0 };
			Vector<Real> smoothYVal{ 0.005, 0.0045, 0.004, 0.0035 }; 
			_dragCoefficientSmooth = new PolynomInterpRealFunc(smoothXVal, smoothYVal, 3);

			Vector<Real> normalXVal{ 0.0, 50.0, 100.0, 150.0 };
			Vector<Real> normalYVal{ 0.0039, 0.0035, 0.0032, 0.003 }; // typical values for a normal baseball
			_dragCoefficientNormal = new PolynomInterpRealFunc(normalXVal, normalYVal, 3);

			Vector<Real> roughXVal{ 0.0, 50.0, 100.0, 150.0 };
			Vector<Real> roughYVal{ 0.0045, 0.004, 0.0035, 0.003 }; // typical values for a rough baseball
			_dragCoefficientRough = new PolynomInterpRealFunc(roughXVal, roughYVal, 3);
		}
		BaseballWithDragCoeffDependentOnSpeedODE() : BaseballWithDragCoeffDependentOnSpeedODE(BaseballType::NormalApprox) {}
		~BaseballWithDragCoeffDependentOnSpeedODE()
		{
			delete _dragCoefficientSmooth;
			delete _dragCoefficientNormal;
			delete _dragCoefficientRough;
		}

		int  getDim() const override { return 4; }
		void derivs(const Real t, const Vector<Real>& x, Vector<Real>& dxdt) const override
		{
			Real v = sqrt(x[2] * x[2] + x[3] * x[3]);				// speed of the projectile

			dxdt[0] = x[2];
			dxdt[1] = x[3];
			dxdt[2] = -dragCoefficient(v) * v * x[2];
			dxdt[3] = -9.81 - dragCoefficient(v) * v * x[3];

			if (x[1] < 0) // if the projectile hits the ground, stop the motion
			{
				dxdt[0] = 0; // stop x motion
				dxdt[1] = 0; // stop y motion
				dxdt[2] = 0; // stop x velocity
				dxdt[3] = 0; // stop y velocity
			}
		}

		Real dragCoefficient(Real v) const
		{
			switch (_baseballType)
			{
			case BaseballType::Smooth:
				return dragCoefficientSmooth(v);
			case BaseballType::Normal:
				return dragCoefficientNormal(v);
			case BaseballType::Rough:
				return dragCoefficientRough(v);
			case BaseballType::NormalApprox:
				return dragCoefficientNormalApprox(v);
			default:
				return dragCoefficientNormal(v); // default to normal baseball
			}
		}
		Real dragCoefficientNormal(Real v) const
		{
			return (*_dragCoefficientNormal)(v); // use interpolated function for normal baseball
		}
		Real dragCoefficientSmooth(Real v) const
		{
			return (*_dragCoefficientSmooth)(v); // use interpolated function for smooth baseball
		}
		Real dragCoefficientRough(Real v) const
		{
			return (*_dragCoefficientRough)(v); // use interpolated function for rough baseball
		}
		Real dragCoefficientNormalApprox(Real v) const
		{
			Real vd = 35;
			Real delta = 5;

			return 0.0039 + 0.0058 / (1 + exp((v - vd) / delta));
		}
	};

	// Struct representing the trajectory of a projectile
	struct ProjectileTrajectory2D
	{
		Real _angle;						// launch angle in radians
		Real _initHeight;				// initial height in meters
		Real _velocity;					// initial velocity in m/s

		Real _timeOfFlight;			// maximum time of flight in seconds
		Real _range;						// range of the projectile in meters

		Vector<Real> _tValues;		// time values for the trajectory points
		Vector<Real> _xValues;		// x coordinates of the trajectory points
		Vector<Real> _yValues;		// y coordinates of the trajectory points

		Vector<Real> _vxValues;		// x velocities of the trajectory points
		Vector<Real> _vyValues;		// y velocities of the trajectory points

		ProjectileTrajectory2D() 
			: _angle(0.0), _initHeight(0.0), _velocity(0.0), _timeOfFlight(0.0), _range(0.0)
		{	}
		ProjectileTrajectory2D(Real angle, Real initHeight, Real velocity)
			: _angle(angle), _initHeight(initHeight), _velocity(velocity), _timeOfFlight(0.0), _range(0.0)
		{	}
		
		LinearInterpRealFunc getXOfT() const
		{
			if (_tValues.size() != _xValues.size())
				throw std::runtime_error("t and x values must have the same size");
			
			return LinearInterpRealFunc(_tValues, _xValues); // linear interpolation for x(t)
		}
		LinearInterpRealFunc getYOfT() const
		{
			if (_tValues.size() != _yValues.size())
				throw std::runtime_error("t and y values must have the same size");
			
			return LinearInterpRealFunc(_tValues, _yValues); // linear interpolation for y(t)
		}
		LinearInterpRealFunc getVxOfT() const
		{
			if (_tValues.size() != _vxValues.size())
				throw std::runtime_error("t and vx values must have the same size");
			
			return LinearInterpRealFunc(_tValues, _vxValues); // linear interpolation for vx(t)
		}
		LinearInterpRealFunc getVyOfT() const
		{
			if (_tValues.size() != _vyValues.size())
				throw std::runtime_error("t and vy values must have the same size");
			
			return LinearInterpRealFunc(_tValues, _vyValues); // linear interpolation for vy(t)
		}
		LinearInterpRealFunc getYOfX() const
		{
			if (_xValues.size() != _yValues.size())
				throw std::runtime_error("x and y values must have the same size");
			return LinearInterpRealFunc(_xValues, _yValues);
		}
		PolynomInterpRealFunc getYOfXPoly() const
		{
			if (_xValues.size() != _yValues.size())
				throw std::runtime_error("x and y values must have the same size");
			
			return PolynomInterpRealFunc(_xValues, _yValues, 3); // cubic interpolation for y(x)
		}
	};

	// Solver for projectile motion
	// Main responsibility is computing trajectory, and correctly handling the case when the projectile hits the ground.
	class ProjectileMotionSolver2D
	{
	private:
		IProjectile2DODE& _projectileODESystem;				// ODE system defining projectile equations of motion
		ODESystemSolution _lastSol;

	public:
		ProjectileMotionSolver2D(IProjectile2DODE& odeSystem) : _projectileODESystem(odeSystem), _lastSol(0.0, 1.0, 4) {}

		Real getRange(Real angle, Real initHeight, Real velocity) const
		{
			// if IProjectileODE is ProjectileInVacuumODE, use its CalcRange method
			if (dynamic_cast<const Projectile2DInVacuumODE*>(&_projectileODESystem) != nullptr)
			{
				return static_cast<const Projectile2DInVacuumODE&>(_projectileODESystem).CalcRange(angle, initHeight, velocity);
			}

			// if not, we have to use numerical integration to find the range

			return 0.0; 
		}
		Real getTimeOfFlight(Real angle, Real initHeight, Real velocity) const
		{
			// if IProjectileODE is ProjectileInVacuumODE, use its TimeOfFlight method
			if (dynamic_cast<const Projectile2DInVacuumODE*>(&_projectileODESystem) != nullptr)
			{
				return static_cast<const Projectile2DInVacuumODE&>(_projectileODESystem).TimeOfFlight(angle, initHeight, velocity);
			}

			// if not, we have to use numerical integration to find the range

			return 0.0; 
		}
		Real getTimeOfFlightVacuum(Real angle, Real initHeight, Real velocity) const
		{
			Real g = 9.81;
			Real v0y = velocity * sin(angle);

			return (v0y + sqrt(v0y * v0y + 2 * g * initHeight)) / g;
		}

		// Function to get the last solution of the ODE system
		ODESystemSolution& getLastODESolution()
		{
			return _lastSol;
		}
		// Helper function to prepare the results after integration and ground hit interpolation
		static ProjectileTrajectory2D prepareResults(Real angle, Real initHeight,	Real velocity,
																							const Vector<Real>& sol_t, const Vector<Real>& sol_x,	const Vector<Real>& sol_y,
																							const Vector<Real>& sol_vx,	const Vector<Real>& sol_vy)
		{
			ProjectileTrajectory2D result(angle, initHeight, velocity);

			int numPoints = sol_t.size();	

			// if the projectile hits the ground, we need to truncate the solution
			if (sol_y.back() < 0)
			{
				// find the last time when y >= 0
				int lastIndex = sol_y.size() - 1;
				while (lastIndex > 0 && sol_y[lastIndex] < 0)
					lastIndex--;

				// Make sure we have valid indices
				if (lastIndex < 0 || lastIndex + 1 >= sol_y.size())
				{
					// Something went wrong - just use the raw data
					result._tValues = sol_t;
					result._xValues = sol_x;
					result._yValues = sol_y;
					result._vxValues = sol_vx;
					result._vyValues = sol_vy;
					result._timeOfFlight = sol_t.back();
					result._range = sol_x.back();
					return result;
				}

				// interpolate the last point to find the exact time of flight
				Real t1 = sol_t[lastIndex];
				Real x1 = sol_x[lastIndex];
				Real y1 = sol_y[lastIndex];
				Real t2 = sol_t[lastIndex + 1];
				Real x2 = sol_x[lastIndex + 1];
				Real y2 = sol_y[lastIndex + 1];

				Real alpha = (0 - y1) / (y2 - y1);
			// Manual linear interpolation (std::lerp requires C++20, we use C++17)
			Real tHitGround = t1 + alpha * (t2 - t1);
			Real xHitGround = x1 + alpha * (x2 - x1);
			Real vxHitGround = sol_vx[lastIndex] + alpha * (sol_vx[lastIndex + 1] - sol_vx[lastIndex]);
			Real vyHitGround = sol_vy[lastIndex] + alpha * (sol_vy[lastIndex + 1] - sol_vy[lastIndex]);

				// Set the computed range and time of flight
				result._timeOfFlight = tHitGround;
				result._range = xHitGround;

				result._tValues.Resize(numPoints);
				result._xValues.Resize(numPoints);
				result._yValues.Resize(numPoints);
				result._vxValues.Resize(numPoints);
				result._vyValues.Resize(numPoints);

				// Fill up to lastIndex
				for (int i = 0; i <= lastIndex; i++)
				{
					result._tValues[i] = sol_t[i];
					result._xValues[i] = sol_x[i];
					result._yValues[i] = sol_y[i];
					result._vxValues[i] = sol_vx[i];
					result._vyValues[i] = sol_vy[i];
				}

				// Add the exact hit ground point
				result._tValues[lastIndex + 1] = tHitGround;
				result._xValues[lastIndex + 1] = xHitGround;
				result._yValues[lastIndex + 1] = 0;
				result._vxValues[lastIndex + 1] = vxHitGround;
				result._vyValues[lastIndex + 1] = vyHitGround;

				// Fill the rest with zeros
				for (int i = lastIndex + 2; i < numPoints; i++)
				{
					result._tValues[i] = sol_t[i];
					result._xValues[i] = xHitGround;
					result._yValues[i] = 0;			// y is zero at the ground
					result._vxValues[i] = 0;		// vx is zero at the ground
					result._vyValues[i] = 0;		// vy is zero at the ground
				}
			}
			else
			{
				result._tValues = sol_t;
				result._xValues = sol_x;
				result._yValues = sol_y;
				result._vxValues = sol_vx;
				result._vyValues = sol_vy;
				
				// Projectile didn't hit the ground - use last values
				result._timeOfFlight = sol_t.back();
				result._range = sol_x.back();
			}

			return result;
		}

		Vector<ProjectileTrajectory2D> solveForAnglesEuler(const Vector<Real>& angles, Real initHeight, Real velocity, Real tMax, Real dt)
		{
			Vector<ProjectileTrajectory2D> results;

			for (Real angle : angles)
			{
				// solve the system of ODEs using the specified ODE system
				ProjectileTrajectory2D result = solveEuler(angle, initHeight, velocity, tMax, dt);
				results.push_back(result);
			}
			return results;
		}
			
		Vector<ProjectileTrajectory2D> solveForAnglesEuler(const Vector<Real>& angles, Real initHeight, Real velocity, Real dt)
		{
			Real tMax = 0.0;
			// we need to find the maximum time of flight for all angles
			for (Real angle : angles)
			{
				Real timeOfFlight = getTimeOfFlightVacuum(angle, initHeight, velocity);
				if (timeOfFlight > tMax)
					tMax = timeOfFlight;
			}
			tMax *= 1.1; // add some margin to the maximum time of flight

			return solveForAnglesEuler(angles, initHeight, velocity, tMax, dt);
		}

		ProjectileTrajectory2D solveEuler(Real angle, Real initHeight, Real velocity, Real tMax, Real dt)
		{
			Vector<Real> initCond = _projectileODESystem.getInitCond(angle, initHeight, velocity);

			int numSteps = tMax / dt;			// expected number of steps in the solution

			ODESystemFixedStepSolver	odeSolver(_projectileODESystem, StepCalculators::EulerStepCalc);
			ODESystemSolution sol = odeSolver.integrate(initCond, 0.0, tMax, numSteps);

			_lastSol = sol; // store the last solution for later use

			Vector<Real> sol_t{ sol.getTValues() }, sol_x{ sol.getXValues(0) }, sol_y{ sol.getXValues(1) }, sol_vx{ sol.getXValues(2) }, sol_vy{ sol.getXValues(3) };

			ProjectileTrajectory2D result = prepareResults(angle, initHeight, velocity, sol_t, sol_x, sol_y, sol_vx, sol_vy);

			return result;
		}

		ProjectileTrajectory2D solveRK4(Real angle, Real initHeight, Real velocity, Real tMax, Real dt)
		{
			ProjectileTrajectory2D result;

			result._angle = angle;
			result._initHeight = initHeight;
			result._velocity = velocity;


			// solve the system of ODEs using the specified ODE system

			// analyse soluton in y direction, and find the time of flight

			// create representation of the trajectory

			return result;
		}

		ProjectileTrajectory2D solveRK5(Real angle, Real initHeight, Real velocity, Real tMax, Real eps=1e-7)
		{
			Vector<Real> initCond = _projectileODESystem.getInitCond(angle, initHeight, velocity);

			ODESystemSolver<RK5_CashKarp_Stepper> odeSolver(_projectileODESystem);
			ODESystemSolution sol = odeSolver.integrate(initCond, 0.0, tMax, 0.01, eps, 0.01);

			_lastSol = sol; // store the last solution for later use

			Vector<Real> sol_t{ sol.getTValues() }, sol_x{ sol.getXValues(0) }, sol_y{ sol.getXValues(1) }, sol_vx{ sol.getXValues(2) }, sol_vy{ sol.getXValues(3) };

			ProjectileTrajectory2D result = prepareResults(angle, initHeight, velocity, sol_t, sol_x, sol_y, sol_vx, sol_vy);

			return result;
		}
	};
}

#endif
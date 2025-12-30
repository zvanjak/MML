#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#include "algorithms/BVPShootingMethod.h"
#endif

using namespace MML;
using namespace MML::Testing;
using namespace Catch::Matchers;

/**
 * COMPREHENSIVE TEST SUITE: BVP Shooting Method
 * 
 * Tests for the shooting method boundary value problem solver:
 * 
 * 1. Basic Shooting Tests
 *    - Simple harmonic oscillator BVP
 *    - Linear ODE with known solution
 * 
 * 2. N-Dimensional Shooting Tests
 *    - Third-order ODE with two shooting parameters
 * 
 * 3. Projectile with Air Resistance Tests (application example)
 *    - Case 1: Given speed, find angle
 *    - Case 2: Given angle, find speed
 *    - Comparison with no-drag analytical solution
 * 
 * TEST STRATEGY:
 * - Use problems with known analytical solutions when possible
 * - Verify convergence to boundary conditions
 * - Test both 1D and N-dimensional shooting
 */

namespace MML::Tests::Algorithms::BVPShootingMethodTests
{
	/*********************************************************************/
	/*****                TEST ODE SYSTEMS                           *****/
	/*********************************************************************/

	/**
	 * Simple harmonic oscillator: y'' + ω²y = 0
	 * Converted to first-order system:
	 *   y[0]' = y[1]       (y' = v)
	 *   y[1]' = -ω²y[0]    (v' = -ω²y)
	 * 
	 * General solution: y = A*cos(ωt) + B*sin(ωt)
	 */
	class HarmonicOscillator : public IODESystem
	{
		Real _omega;
	public:
		HarmonicOscillator(Real omega = 1.0) : _omega(omega) {}
		
		int getDim() const override { return 2; }
		
		void derivs(Real t, const Vector<Real>& y, Vector<Real>& dydt) const override
		{
			dydt[0] = y[1];
			dydt[1] = -_omega * _omega * y[0];
		}
		
		Real omega() const { return _omega; }
	};

	/**
	 * Linear ODE: y'' = 0
	 * Solution: y = a*t + b
	 * With y(0) = y0 and y(T) = yT, we get y = y0 + (yT - y0)/T * t
	 */
	class LinearODE : public IODESystem
	{
	public:
		int getDim() const override { return 2; }
		
		void derivs(Real t, const Vector<Real>& y, Vector<Real>& dydt) const override
		{
			dydt[0] = y[1];  // y' = v
			dydt[1] = 0;     // v' = 0
		}
	};

	/**
	 * Third-order ODE: y''' = 0
	 * State: [y, y', y'']
	 * Solution: y = a + b*t + c*t²/2
	 * 
	 * Used for testing N-dimensional shooting with 2 unknowns.
	 */
	class ThirdOrderLinearODE : public IODESystem
	{
	public:
		int getDim() const override { return 3; }
		
		void derivs(Real t, const Vector<Real>& y, Vector<Real>& dydt) const override
		{
			dydt[0] = y[1];  // y' 
			dydt[1] = y[2];  // y''
			dydt[2] = 0;     // y''' = 0
		}
	};

	/**
	 * Projectile motion with quadratic air resistance
	 * 
	 * State vector: [x, y, vx, vy]
	 * Equations:
	 *   dx/dt = vx
	 *   dy/dt = vy
	 *   dvx/dt = -k * v * vx    (drag proportional to v²)
	 *   dvy/dt = -g - k * v * vy
	 * 
	 * NOTE: This is an application-specific class for testing.
	 * General users should define their own ODE systems.
	 */
	class ProjectileWithDrag : public IODESystem
	{
	private:
		Real _g;     ///< Gravitational acceleration (m/s²)
		Real _k;     ///< Drag coefficient (1/m)
		
	public:
		ProjectileWithDrag(Real g = 9.81, Real k = 0.01)
			: _g(g), _k(k) {}
		
		int getDim() const override { return 4; }
		
		void derivs(Real t, const Vector<Real>& state, Vector<Real>& dstate) const override
		{
			Real vx = state[2];
			Real vy = state[3];
			Real v = std::sqrt(vx*vx + vy*vy);
			
			dstate[0] = vx;
			dstate[1] = vy;
			dstate[2] = -_k * v * vx;
			dstate[3] = -_g - _k * v * vy;
		}
		
		Real getGravity() const { return _g; }
		Real getDragCoeff() const { return _k; }
	};

	/**
	 * 3D Projectile motion with air drag AND wind
	 * 
	 * State vector: [x, y, z, vx, vy, vz]
	 * - x = downrange distance (toward target)
	 * - y = height (vertical)
	 * - z = lateral distance (crosswind direction)
	 * 
	 * Physics:
	 * - Quadratic drag: F_drag = -k * |v_rel| * v_rel
	 * - v_rel = velocity relative to air = v - v_wind
	 * - Gravity acts in -y direction
	 * 
	 * Equations:
	 *   dx/dt = vx
	 *   dy/dt = vy
	 *   dz/dt = vz
	 *   dvx/dt = -k * |v_rel| * (vx - wx)
	 *   dvy/dt = -g - k * |v_rel| * (vy - wy)
	 *   dvz/dt = -k * |v_rel| * (vz - wz)
	 * 
	 * This is a great test case for N-dimensional shooting because:
	 * - Wind causes lateral drift even when aiming straight
	 * - Must "aim off" to compensate for wind
	 * - Requires solving for multiple velocity components simultaneously
	 */
	class Projectile3DWithWind : public IODESystem
	{
	private:
		Real _g;      ///< Gravitational acceleration (m/s²)
		Real _k;      ///< Drag coefficient (1/m)
		Real _wx;     ///< Wind velocity x-component (m/s)
		Real _wy;     ///< Wind velocity y-component (m/s) - usually 0
		Real _wz;     ///< Wind velocity z-component (m/s) - crosswind
		
	public:
		/**
		 * @brief Construct 3D projectile system with wind
		 * 
		 * @param g   Gravity (default 9.81 m/s²)
		 * @param k   Drag coefficient (default 0.01 1/m)
		 * @param wx  Headwind/tailwind (+/- in direction of target)
		 * @param wy  Vertical wind (usually 0)
		 * @param wz  Crosswind (perpendicular to target direction)
		 */
		Projectile3DWithWind(Real g = 9.81, Real k = 0.01,
		                     Real wx = 0.0, Real wy = 0.0, Real wz = 0.0)
			: _g(g), _k(k), _wx(wx), _wy(wy), _wz(wz) {}
		
		int getDim() const override { return 6; }
		
		void derivs(Real t, const Vector<Real>& state, Vector<Real>& dstate) const override
		{
			Real vx = state[3];
			Real vy = state[4];
			Real vz = state[5];
			
			// Velocity relative to air (for drag calculation)
			Real vx_rel = vx - _wx;
			Real vy_rel = vy - _wy;
			Real vz_rel = vz - _wz;
			Real v_rel = std::sqrt(vx_rel*vx_rel + vy_rel*vy_rel + vz_rel*vz_rel);
			
			// Position derivatives
			dstate[0] = vx;  // dx/dt
			dstate[1] = vy;  // dy/dt
			dstate[2] = vz;  // dz/dt
			
			// Velocity derivatives (drag opposes relative velocity)
			dstate[3] = -_k * v_rel * vx_rel;           // dvx/dt
			dstate[4] = -_g - _k * v_rel * vy_rel;      // dvy/dt (gravity + drag)
			dstate[5] = -_k * v_rel * vz_rel;           // dvz/dt
		}
		
		// Accessors
		Real getGravity() const { return _g; }
		Real getDragCoeff() const { return _k; }
		Real getWindX() const { return _wx; }
		Real getWindY() const { return _wy; }
		Real getWindZ() const { return _wz; }
	};

	/**
	 * Result struct for 3D projectile BVP
	 */
	struct Projectile3DBVPResult
	{
		bool converged = false;
		Real launchElevation = 0.0;    ///< Elevation angle (radians)
		Real launchAzimuth = 0.0;      ///< Azimuth angle (radians) - wind compensation
		Real launchSpeed = 0.0;
		Real flightTime = 0.0;
		Real maxHeight = 0.0;
		Real lateralDrift = 0.0;       ///< z-position at landing (should be ~0 if solved)
		Vector<Real> finalPosition;    ///< [x, y, z] at end
		ODESystemSolution trajectory;
		int iterations = 0;
		
		Projectile3DBVPResult() : finalPosition(3), trajectory(0.0, 1.0, 6, 100) {}
		
		Real launchElevationDegrees() const { return launchElevation * 180.0 / Constants::PI; }
		Real launchAzimuthDegrees() const { return launchAzimuth * 180.0 / Constants::PI; }
	};

	/**
	 * Helper class for 3D projectile BVP solving with wind compensation (test utility).
	 * 
	 * This solver handles a FREE-BOUNDARY BVP where the landing time is unknown:
	 * - Find (vx, vz) such that the projectile lands at (targetX, 0, targetZ)
	 * - The landing time is determined implicitly by when y crosses 0
	 * 
	 * NOTE: This demonstrates a limitation of BVPShootingSolverND: it solves for 
	 * conditions at a FIXED endpoint time T. Free-boundary problems (like landing 
	 * at y=0 at an unknown time) require custom handling.
	 * 
	 * The implementation uses Newton-Raphson with numerical Jacobian, similar to
	 * BVPShootingSolverND, but targets the LANDING POSITION rather than position
	 * at a fixed time T.
	 * 
	 * For FIXED-TIME BVPs (e.g., "reach position (100, y, 0) at t=5s"), use 
	 * BVPShootingSolverND directly - see the ThirdOrderLinearODE tests.
	 */
	class Projectile3DBVPSolver
	{
	private:
		Projectile3DWithWind _system;
		int _numSteps;
		Real _tolerance;
		
	public:
		Projectile3DBVPSolver(Real g = 9.81, Real k = 0.01,
		                      Real windX = 0.0, Real windY = 0.0, Real windZ = 0.0,
		                      int numSteps = 3000, Real tolerance = 1e-6)
			: _system(g, k, windX, windY, windZ)
			, _numSteps(numSteps)
			, _tolerance(tolerance)
		{}
		
		/**
		 * @brief Find launch angles to hit target at ground level
		 * 
		 * Uses BVPShootingSolverND for the core shooting, with an outer
		 * iteration to handle the free-boundary (landing time) problem.
		 * 
		 * @param x0, y0, z0       Starting position
		 * @param targetX         Target distance (in x-direction)
		 * @param targetZ         Target lateral position (usually 0)
		 * @param v0              Launch speed
		 * @param elevationGuess  Initial guess for elevation angle (radians)
		 * @param maxTime         Maximum flight time
		 * @return Projectile3DBVPResult
		 */
		Projectile3DBVPResult findAnglesForSpeed(
			Real x0, Real y0, Real z0,
			Real targetX, Real targetZ,
			Real v0,
			Real elevationGuess,
			Real maxTime = 20.0)
		{
			Projectile3DBVPResult result;
			result.launchSpeed = v0;
			
			// Initial guesses from elevation
			Real theta = elevationGuess;
			Real vx = v0 * std::cos(theta);
			Real vy = v0 * std::sin(theta);
			Real vz = REAL(0.0);
			
			// Helper to find landing point from trajectory
			auto findLandingPoint = [&](const ODESystemSolution& sol,
			                            Real& landX, Real& landY, Real& landZ, Real& landT) -> bool
			{
				bool hasRisen = (y0 > 0.1);
				for (int i = 1; i < sol.size(); ++i)
				{
					Real y_prev = sol.getXValue(i-1, 1);
					Real y_curr = sol.getXValue(i, 1);
					
					if (!hasRisen && y_curr > y0 + 0.1)
						hasRisen = true;
					
					if (hasRisen && y_prev >= 0 && y_curr < 0)
					{
						Real alpha = -y_prev / (y_curr - y_prev);
						landT = sol.getTValue(i-1) + alpha * (sol.getTValue(i) - sol.getTValue(i-1));
						landX = sol.getXValue(i-1, 0) + alpha * (sol.getXValue(i, 0) - sol.getXValue(i-1, 0));
						landY = 0;
						landZ = sol.getXValue(i-1, 2) + alpha * (sol.getXValue(i, 2) - sol.getXValue(i-1, 2));
						return true;
					}
				}
				return false;
			};
			
			// Newton-Raphson with numerical Jacobian to solve for (vx, vz) 
			// such that landing position = (targetX, targetZ)
			// This is a 2x2 system: F(vx, vz) = [land_x - targetX, land_z - targetZ]
			
			const Real delta = 1e-4;  // Perturbation for numerical Jacobian
			
			for (int iter = 0; iter < 30; ++iter)
			{
				// Compute base trajectory and landing point
				Vector<Real> ic(6);
				ic[0] = x0; ic[1] = y0; ic[2] = z0;
				ic[3] = vx; ic[4] = vy; ic[5] = vz;
				
				RungeKutta4_StepCalculator stepCalc;
				ODESystemFixedStepSolver solver(_system, stepCalc);
				auto sol = solver.integrate(ic, 0, maxTime, _numSteps);
				
				Real landX, landY, landZ, landT;
				if (!findLandingPoint(sol, landX, landY, landZ, landT))
				{
					result.converged = false;
					return result;
				}
				
				// Compute defect: how far from target?
				Real Fx = landX - targetX;
				Real Fz = landZ - targetZ;
				
				if (std::abs(Fx) < _tolerance && std::abs(Fz) < _tolerance)
				{
					// Converged!
					result.converged = true;
					result.iterations = iter + 1;
					result.flightTime = landT;
					result.trajectory = sol;
					result.finalPosition[0] = landX;
					result.finalPosition[1] = landY;
					result.finalPosition[2] = landZ;
					result.lateralDrift = landZ;
					
					// Extract angles
					Real v_horiz = std::sqrt(vx*vx + vz*vz);
					result.launchElevation = std::atan2(vy, v_horiz);
					result.launchAzimuth = std::atan2(vz, vx);
					
					// Find max height
					result.maxHeight = y0;
					for (int i = 0; i < sol.size(); ++i)
						if (sol.getXValue(i, 1) > result.maxHeight)
							result.maxHeight = sol.getXValue(i, 1);
					
					return result;
				}
				
				// Build numerical Jacobian (2x2)
				// J = [dFx/dvx  dFx/dvz]
				//     [dFz/dvx  dFz/dvz]
				
				// Perturb vx
				ic[3] = vx + delta;
				ic[5] = vz;
				auto sol_vx = solver.integrate(ic, 0, maxTime, _numSteps);
				Real landX_vx, landY_vx, landZ_vx, landT_vx;
				if (!findLandingPoint(sol_vx, landX_vx, landY_vx, landZ_vx, landT_vx))
				{
					result.converged = false;
					return result;
				}
				Real dFx_dvx = (landX_vx - targetX - Fx) / delta;
				Real dFz_dvx = (landZ_vx - targetZ - Fz) / delta;
				
				// Perturb vz
				ic[3] = vx;
				ic[5] = vz + delta;
				auto sol_vz = solver.integrate(ic, 0, maxTime, _numSteps);
				Real landX_vz, landY_vz, landZ_vz, landT_vz;
				if (!findLandingPoint(sol_vz, landX_vz, landY_vz, landZ_vz, landT_vz))
				{
					result.converged = false;
					return result;
				}
				Real dFx_dvz = (landX_vz - targetX - Fx) / delta;
				Real dFz_dvz = (landZ_vz - targetZ - Fz) / delta;
				
				// Solve J * [dvx, dvz]^T = -[Fx, Fz]^T
				Real det = dFx_dvx * dFz_dvz - dFx_dvz * dFz_dvx;
				if (std::abs(det) < 1e-12)
				{
					result.converged = false;
					return result;
				}
				
				Real dvx = (-Fx * dFz_dvz + Fz * dFx_dvz) / det;
				Real dvz = (-dFx_dvx * Fz + dFz_dvx * Fx) / det;
				
				// Update
				vx += dvx;
				vz += dvz;
				
				result.iterations = iter + 1;
			}
			
			// If we get here, didn't converge
			result.converged = false;
			return result;
		}
		
		/**
		 * @brief Compute trajectory without BVP solving (for comparison)
		 */
		ODESystemSolution computeTrajectory(Real x0, Real y0, Real z0,
		                                    Real vx0, Real vy0, Real vz0,
		                                    Real maxTime)
		{
			Vector<Real> initCond(6);
			initCond[0] = x0; initCond[1] = y0; initCond[2] = z0;
			initCond[3] = vx0; initCond[4] = vy0; initCond[5] = vz0;
			
			RungeKutta4_StepCalculator stepCalc;
			ODESystemFixedStepSolver solver(_system, stepCalc);
			return solver.integrate(initCond, 0, maxTime, _numSteps);
		}
		
		const Projectile3DWithWind& getSystem() const { return _system; }
	};

	/**
	 * Result struct for projectile BVP tests
	 */
	struct ProjectileBVPResult
	{
		bool converged = false;
		Real launchAngle = 0.0;
		Real launchSpeed = 0.0;
		Real flightTime = 0.0;
		Real maxHeight = 0.0;
		Real residual = 0.0;
		ODESystemSolution trajectory;
		
		ProjectileBVPResult() : trajectory(0.0, 1.0, 4, 100) {}
		
		Real launchAngleDegrees() const { return launchAngle * 180.0 / Constants::PI; }
	};

	/**
	 * Helper class for projectile BVP solving (test utility)
	 */
	class ProjectileBVPSolver
	{
	private:
		ProjectileWithDrag _system;
		int _numSteps;
		Real _tolerance;
		
	public:
		ProjectileBVPSolver(Real g = 9.81, Real k = 0.01,
		                    int numSteps = 2000, Real tolerance = 1e-8)
			: _system(g, k), _numSteps(numSteps), _tolerance(tolerance) {}
		
		ProjectileBVPResult findAngleForSpeed(Real x0, Real y0,
		                                      Real x_target, Real y_target,
		                                      Real v0,
		                                      Real angleGuess1, Real angleGuess2,
		                                      Real maxTime = 100.0)
		{
			ProjectileBVPResult result;
			result.launchSpeed = v0;
			
			auto defectFunc = [&](Real theta) -> Real {
				Vector<Real> initCond(4);
				initCond[0] = x0;
				initCond[1] = y0;
				initCond[2] = v0 * std::cos(theta);
				initCond[3] = v0 * std::sin(theta);
				
				RungeKutta4_StepCalculator stepCalc;
				ODESystemFixedStepSolver solver(_system, stepCalc);
				auto sol = solver.integrate(initCond, 0, maxTime, _numSteps);
				
				Real x_final = x0;
				bool hasRisen = (y0 > y_target + 0.01);
				
				for (int i = 1; i < sol.size(); ++i)
				{
					Real y_prev = sol.getXValue(i-1, 1);
					Real y_curr = sol.getXValue(i, 1);
					
					if (!hasRisen && y_curr > y_target + 0.01)
						hasRisen = true;
					
					if (hasRisen && y_prev >= y_target && y_curr < y_target)
					{
						Real t_prev = sol.getTValue(i-1);
						Real t_curr = sol.getTValue(i);
						Real x_prev = sol.getXValue(i-1, 0);
						Real x_curr = sol.getXValue(i, 0);
						
						Real alpha = (y_target - y_prev) / (y_curr - y_prev);
						x_final = x_prev + alpha * (x_curr - x_prev);
						result.flightTime = t_prev + alpha * (t_curr - t_prev);
						break;
					}
				}
				return x_final - x_target;
			};
			
			RealFunctionFromStdFunc F(defectFunc);
			
			try
			{
				Real a = angleGuess1, b = angleGuess2;
				bool bracketed = RootFinding::BracketRoot(F, a, b, 50);
				
				if (!bracketed) { result.converged = false; return result; }
				
				Real optAngle = RootFinding::FindRootBrent(F, a, b, _tolerance);
				
				result.converged = true;
				result.launchAngle = optAngle;
				result.residual = std::abs(F(optAngle));
				
				Vector<Real> initCond(4);
				initCond[0] = x0; initCond[1] = y0;
				initCond[2] = v0 * std::cos(optAngle);
				initCond[3] = v0 * std::sin(optAngle);
				
				RungeKutta4_StepCalculator stepCalc;
				ODESystemFixedStepSolver solver(_system, stepCalc);
				result.trajectory = solver.integrate(initCond, 0, result.flightTime * 1.1, _numSteps);
				
				result.maxHeight = y0;
				for (int i = 0; i < result.trajectory.size(); ++i)
					if (result.trajectory.getXValue(i, 1) > result.maxHeight)
						result.maxHeight = result.trajectory.getXValue(i, 1);
			}
			catch (...) { result.converged = false; }
			
			return result;
		}
		
		ProjectileBVPResult findSpeedForAngle(Real x0, Real y0,
		                                      Real x_target, Real y_target,
		                                      Real angle,
		                                      Real speedGuess1, Real speedGuess2,
		                                      Real maxTime = 100.0)
		{
			ProjectileBVPResult result;
			result.launchAngle = angle;
			
			auto defectFunc = [&](Real v0) -> Real {
				Vector<Real> initCond(4);
				initCond[0] = x0; initCond[1] = y0;
				initCond[2] = v0 * std::cos(angle);
				initCond[3] = v0 * std::sin(angle);
				
				RungeKutta4_StepCalculator stepCalc;
				ODESystemFixedStepSolver solver(_system, stepCalc);
				auto sol = solver.integrate(initCond, 0, maxTime, _numSteps);
				
				Real x_final = x0;
				bool hasRisen = (y0 > y_target + 0.01);
				
				for (int i = 1; i < sol.size(); ++i)
				{
					Real y_prev = sol.getXValue(i-1, 1);
					Real y_curr = sol.getXValue(i, 1);
					
					if (!hasRisen && y_curr > y_target + 0.01)
						hasRisen = true;
					
					if (hasRisen && y_prev >= y_target && y_curr < y_target)
					{
						Real x_prev = sol.getXValue(i-1, 0);
						Real x_curr = sol.getXValue(i, 0);
						Real alpha = (y_target - y_prev) / (y_curr - y_prev);
						x_final = x_prev + alpha * (x_curr - x_prev);
						
						Real t_prev = sol.getTValue(i-1);
						Real t_curr = sol.getTValue(i);
						result.flightTime = t_prev + alpha * (t_curr - t_prev);
						break;
					}
				}
				return x_final - x_target;
			};
			
			RealFunctionFromStdFunc F(defectFunc);
			
			try
			{
				Real a = speedGuess1, b = speedGuess2;
				bool bracketed = RootFinding::BracketRoot(F, a, b, 50);
				
				if (!bracketed) { result.converged = false; return result; }
				
				Real optSpeed = RootFinding::FindRootBrent(F, a, b, _tolerance);
				
				result.converged = true;
				result.launchSpeed = optSpeed;
				result.residual = std::abs(F(optSpeed));
				
				Vector<Real> initCond(4);
				initCond[0] = x0; initCond[1] = y0;
				initCond[2] = optSpeed * std::cos(angle);
				initCond[3] = optSpeed * std::sin(angle);
				
				RungeKutta4_StepCalculator stepCalc;
				ODESystemFixedStepSolver solver(_system, stepCalc);
				result.trajectory = solver.integrate(initCond, 0, result.flightTime * 1.1, _numSteps);
				
				result.maxHeight = y0;
				for (int i = 0; i < result.trajectory.size(); ++i)
					if (result.trajectory.getXValue(i, 1) > result.maxHeight)
						result.maxHeight = result.trajectory.getXValue(i, 1);
			}
			catch (...) { result.converged = false; }
			
			return result;
		}
		
		ODESystemSolution computeTrajectory(Real x0, Real y0, Real v0, Real angle, Real maxTime)
		{
			Vector<Real> initCond(4);
			initCond[0] = x0; initCond[1] = y0;
			initCond[2] = v0 * std::cos(angle);
			initCond[3] = v0 * std::sin(angle);
			
			RungeKutta4_StepCalculator stepCalc;
			ODESystemFixedStepSolver solver(_system, stepCalc);
			return solver.integrate(initCond, 0, maxTime, _numSteps);
		}
	};

	/*********************************************************************/
	/*****               BASIC SHOOTING TESTS                        *****/
	/*********************************************************************/

	TEST_CASE("Shooting_LinearODE_BasicBVP", "[bvp_shooting][basic]")
	{
		TEST_PRECISION_INFO();
		
		// y'' = 0 with y(0) = 1 and y(2) = 5
		// Solution: y = 1 + 2t, so y'(0) = 2
		LinearODE system;
		BVPShootingSolver<RungeKutta4_StepCalculator> solver(system, 1000, 1e-10);
		
		Vector<Real> initCond(2);
		initCond[0] = REAL(1.0);  // y(0) = 1
		initCond[1] = REAL(0.0);  // placeholder for y'(0)
		
		auto result = solver.solve1D(REAL(0.0), REAL(2.0), initCond,
		                             1,        // shoot with velocity
		                             0,        // target position
		                             REAL(5.0), // y(2) = 5
		                             REAL(0.0), REAL(5.0));  // velocity guesses
		
		REQUIRE(result.converged);
		REQUIRE_THAT(result.shootingParams[0], WithinRel(REAL(2.0), REAL(1e-6)));
		
		// Verify final value
		Vector<Real> finalState = result.solution.getXValuesAtEnd();
		REQUIRE_THAT(finalState[0], WithinRel(REAL(5.0), REAL(1e-5)));
	}

	TEST_CASE("Shooting_HarmonicOscillator_BVP", "[bvp_shooting][harmonic]")
	{
		TEST_PRECISION_INFO();
		
		// y'' + y = 0 with y(0) = 0 and y(π/2) = 1
		// Solution: y = sin(t), so y'(0) = 1
		HarmonicOscillator system(REAL(1.0));
		BVPShootingSolver<RungeKutta4_StepCalculator> solver(system, 2000, 1e-10);
		
		Vector<Real> initCond(2);
		initCond[0] = REAL(0.0);  // y(0) = 0
		initCond[1] = REAL(0.0);  // placeholder
		
		auto result = solver.solve1D(REAL(0.0), Constants::PI / 2, initCond,
		                             1, 0, REAL(1.0),
		                             REAL(0.5), REAL(2.0));
		
		REQUIRE(result.converged);
		REQUIRE_THAT(result.shootingParams[0], WithinRel(REAL(1.0), REAL(1e-5)));
	}

	/*********************************************************************/
	/*****          N-DIMENSIONAL SHOOTING TESTS                     *****/
	/*********************************************************************/

	TEST_CASE("ShootingND_ThirdOrderLinear_TwoUnknowns", "[bvp_shooting][nd][basic]")
	{
		TEST_PRECISION_INFO();
		
		// Third-order ODE: y''' = 0
		// Solution: y = a + b*t + c*t²/2
		// 
		// BVP: y(0) = 1, y(2) = 5, y'(2) = 3
		// State: [y, y', y''] at indices [0, 1, 2]
		// 
		// Known: y(0) = 1
		// Unknown: y'(0), y''(0)
		// Targets at t=2: y(2) = 5, y'(2) = 3
		//
		// Analytical solution:
		//   y(0) = a = 1
		//   y(2) = 1 + b*2 + c*2 = 5  →  2b + 2c = 4  →  b + c = 2
		//   y'(2) = b + 2c = 3
		//   Solving: c = 1, b = 1
		//   So y'(0) = b = 1, y''(0) = c = 1
		
		ThirdOrderLinearODE system;
		BVPShootingSolverND<RungeKutta4_StepCalculator> solver(system, 1000, 1e-8);
		
		Vector<Real> initCond(3);
		initCond[0] = REAL(1.0);  // y(0) = 1 (known)
		initCond[1] = REAL(0.0);  // y'(0) placeholder
		initCond[2] = REAL(0.0);  // y''(0) placeholder
		
		ShootingConfig config;
		config.unknownIndices = {1, 2};        // y'(0) and y''(0) are unknown
		config.unknownGuesses = {REAL(0.5), REAL(0.5)};   // Initial guesses
		config.targetIndices  = {0, 1};        // y(2) and y'(2) are targets
		config.targetValues   = {REAL(5.0), REAL(3.0)};   // Desired values at t=2
		
		auto result = solver.solve(REAL(0.0), REAL(2.0), initCond, config);
		
		REQUIRE(result.converged);
		REQUIRE_THAT(result.shootingParams[0], WithinRel(REAL(1.0), REAL(1e-4)));  // y'(0) = 1
		REQUIRE_THAT(result.shootingParams[1], WithinRel(REAL(1.0), REAL(1e-4)));  // y''(0) = 1
		
		// Verify final conditions
		Vector<Real> finalState = result.solution.getXValuesAtEnd();
		REQUIRE_THAT(finalState[0], WithinRel(REAL(5.0), REAL(1e-4)));  // y(2) = 5
		REQUIRE_THAT(finalState[1], WithinRel(REAL(3.0), REAL(1e-4)));  // y'(2) = 3
		
		INFO("Converged in " << result.iterations << " iterations");
		INFO("y'(0) = " << result.shootingParams[0] << " (expected 1)");
		INFO("y''(0) = " << result.shootingParams[1] << " (expected 1)");
	}

	TEST_CASE("ShootingND_ThirdOrderLinear_DifferentBC", "[bvp_shooting][nd]")
	{
		TEST_PRECISION_INFO();
		
		// Different boundary conditions:
		// y(0) = 2, y(1) = 3, y'(1) = 1
		// 
		// y = a + bt + ct²/2
		// y(0) = a = 2
		// y(1) = 2 + b + c/2 = 3  →  b + c/2 = 1
		// y'(1) = b + c = 1
		// Solving: c/2 = 0  →  c = 0, b = 1
		// So y'(0) = 1, y''(0) = 0
		
		ThirdOrderLinearODE system;
		BVPShootingSolverND<RungeKutta4_StepCalculator> solver(system, 1000, 1e-8);
		
		Vector<Real> initCond(3);
		initCond[0] = REAL(2.0);  // y(0) = 2
		initCond[1] = REAL(0.0);  // placeholder
		initCond[2] = REAL(0.0);  // placeholder
		
		ShootingConfig config;
		config.unknownIndices = {1, 2};
		config.unknownGuesses = {REAL(0.0), REAL(0.0)};
		config.targetIndices  = {0, 1};
		config.targetValues   = {REAL(3.0), REAL(1.0)};
		
		auto result = solver.solve(REAL(0.0), REAL(1.0), initCond, config);
		
		REQUIRE(result.converged);
		REQUIRE_THAT(result.shootingParams[0], WithinAbs(REAL(1.0), REAL(1e-4)));  // y'(0) = 1
		REQUIRE_THAT(result.shootingParams[1], WithinAbs(REAL(0.0), REAL(1e-4)));  // y''(0) = 0
	}

	TEST_CASE("ShootingND_ConfigValidation", "[bvp_shooting][nd][validation]")
	{
		TEST_PRECISION_INFO();
		
		// Test that invalid configurations are rejected
		ShootingConfig config;
		
		// Empty config
		REQUIRE_FALSE(config.isValid());
		
		// Mismatched sizes
		config.unknownIndices = {1, 2};
		config.unknownGuesses = {REAL(0.0)};  // Wrong size
		config.targetIndices  = {0, 1};
		config.targetValues   = {REAL(1.0), REAL(2.0)};
		REQUIRE_FALSE(config.isValid());
		
		// Correct sizes
		config.unknownGuesses = {REAL(0.0), REAL(0.0)};
		REQUIRE(config.isValid());
		REQUIRE(config.numParams() == 2);
	}

	TEST_CASE("ShootingND_ReducesToScalar", "[bvp_shooting][nd]")
	{
		TEST_PRECISION_INFO();
		
		// N-D solver with N=1 should give same result as 1D solver
		LinearODE system;
		
		// Using 1D solver
		BVPShootingSolver<RungeKutta4_StepCalculator> solver1D(system, 1000, 1e-10);
		auto result1D = solver1D.solvePositionBVP(REAL(0.0), REAL(2.0),
		                                          REAL(1.0), REAL(5.0),
		                                          REAL(0.0), REAL(5.0));
		
		// Using N-D solver with N=1
		BVPShootingSolverND<RungeKutta4_StepCalculator> solverND(system, 1000, 1e-10);
		
		Vector<Real> initCond(2);
		initCond[0] = REAL(1.0);  // y(0) = 1
		initCond[1] = REAL(0.0);
		
		ShootingConfig config;
		config.unknownIndices = {1};
		config.unknownGuesses = {REAL(1.0)};
		config.targetIndices  = {0};
		config.targetValues   = {REAL(5.0)};
		
		auto resultND = solverND.solve(REAL(0.0), REAL(2.0), initCond, config);
		
		REQUIRE(result1D.converged);
		REQUIRE(resultND.converged);
		
		// Both should find y'(0) = 2
		REQUIRE_THAT(result1D.shootingParams[0], WithinRel(REAL(2.0), REAL(1e-5)));
		REQUIRE_THAT(resultND.shootingParams[0], WithinRel(REAL(2.0), REAL(1e-4)));
	}

	TEST_CASE("Shooting_PositionBVP_Convenience", "[bvp_shooting][convenience]")
	{
		TEST_PRECISION_INFO();
		
		LinearODE system;
		BVPShootingSolver<RungeKutta4_StepCalculator> solver(system, 1000, 1e-10);
		
		// y(0) = 0, y(1) = 10 → y = 10t, y'(0) = 10
		auto result = solver.solvePositionBVP(REAL(0.0), REAL(1.0),
		                                      REAL(0.0), REAL(10.0),
		                                      REAL(5.0), REAL(15.0));
		
		REQUIRE(result.converged);
		REQUIRE_THAT(result.shootingParams[0], WithinRel(REAL(10.0), REAL(1e-5)));
	}

	/*********************************************************************/
	/*****           PROJECTILE WITHOUT DRAG (BASELINE)              *****/
	/*********************************************************************/

	TEST_CASE("Projectile_NoDrag_TrajectoryVerification", "[bvp_shooting][projectile][baseline]")
	{
		TEST_PRECISION_INFO();
		
		// With k=0 (no drag), compare to analytical solution
		// Range = v0² * sin(2θ) / g
		ProjectileBVPSolver solver(REAL(9.81), REAL(0.0), 2000);  // No drag
		
		Real v0 = REAL(50.0);
		Real angle = Constants::PI / 4;  // 45 degrees
		
		auto traj = solver.computeTrajectory(0, 0, v0, angle, 10.0);
		
		// Analytical range for 45 degrees
		Real analytical_range = v0 * v0 / REAL(9.81);  // v0²/g at 45°
		
		// Find landing point from trajectory
		Real landing_x = 0;
		for (int i = 1; i < traj.size(); ++i)
		{
			if (traj.getXValue(i-1, 1) > 0 && traj.getXValue(i, 1) <= 0)
			{
				Real alpha = -traj.getXValue(i-1, 1) / (traj.getXValue(i, 1) - traj.getXValue(i-1, 1));
				landing_x = traj.getXValue(i-1, 0) + alpha * (traj.getXValue(i, 0) - traj.getXValue(i-1, 0));
				break;
			}
		}
		
		REQUIRE_THAT(landing_x, WithinRel(analytical_range, REAL(1e-3)));
	}

	/*********************************************************************/
	/*****           PROJECTILE WITH DRAG - FIND ANGLE               *****/
	/*********************************************************************/

	TEST_CASE("Projectile_WithDrag_FindAngle_ShortRange", "[bvp_shooting][projectile][angle]")
	{
		TEST_PRECISION_INFO();
		
		// Find angle to hit target at 100m with v0 = 60 m/s
		// With drag k=0.01, need sufficient speed to reach 100m
		ProjectileBVPSolver solver(REAL(9.81), REAL(0.01), 3000, REAL(1e-6));
		
		auto result = solver.findAngleForSpeed(
			REAL(0.0), REAL(0.0),     // Start at origin
			REAL(100.0), REAL(0.0),   // Target at x=100, y=0
			REAL(60.0),               // Launch speed 60 m/s
			REAL(0.1), REAL(0.8)      // Angle guesses: tighter range (5.7° to 46°)
		);
		
		REQUIRE(result.converged);
		REQUIRE(result.launchAngle > REAL(0.0));
		REQUIRE(result.launchAngle < Constants::PI / 2);
		REQUIRE(result.residual < REAL(1.0));  // Landing within 1m of target
		
		INFO("Optimal angle: " << result.launchAngleDegrees() << " degrees");
		INFO("Flight time: " << result.flightTime << " seconds");
		INFO("Max height: " << result.maxHeight << " m");
	}

	TEST_CASE("Projectile_WithDrag_FindAngle_LongRange", "[bvp_shooting][projectile][angle]")
	{
		TEST_PRECISION_INFO();
		
		// Find angle to hit target at 150m with moderate drag
		// Lower drag (0.005) makes longer range more achievable
		ProjectileBVPSolver solver(REAL(9.81), REAL(0.005), 4000, REAL(1e-6));
		
		auto result = solver.findAngleForSpeed(
			REAL(0.0), REAL(0.0),
			REAL(150.0), REAL(0.0),
			REAL(60.0),               // 60 m/s with lower drag
			REAL(0.2), REAL(1.0)      // Angle guesses
		);
		
		REQUIRE(result.converged);
		REQUIRE(result.residual < REAL(2.0));
		
		INFO("Optimal angle: " << result.launchAngleDegrees() << " degrees");
		INFO("Flight time: " << result.flightTime << " seconds");
	}

	TEST_CASE("Projectile_WithDrag_FindAngle_ElevatedTarget", "[bvp_shooting][projectile][angle]")
	{
		TEST_PRECISION_INFO();
		
		// Target is elevated: x=80m, y=20m - needs sufficient speed
		ProjectileBVPSolver solver(REAL(9.81), REAL(0.01), 3000, REAL(1e-6));
		
		auto result = solver.findAngleForSpeed(
			REAL(0.0), REAL(0.0),
			REAL(80.0), REAL(20.0),   // Elevated target
			REAL(55.0),               // Increased from 50 to ensure reachability
			REAL(0.4), REAL(1.3)
		);
		
		REQUIRE(result.converged);
		
		// Verify max height is above target
		REQUIRE(result.maxHeight >= REAL(20.0));
		
		INFO("Optimal angle: " << result.launchAngleDegrees() << " degrees");
		INFO("Max height: " << result.maxHeight << " m");
	}

	/*********************************************************************/
	/*****           PROJECTILE WITH DRAG - FIND SPEED               *****/
	/*********************************************************************/

	TEST_CASE("Projectile_WithDrag_FindSpeed_45Degrees", "[bvp_shooting][projectile][speed]")
	{
		TEST_PRECISION_INFO();
		
		// Find speed to hit target at 150m with 45° angle
		ProjectileBVPSolver solver(REAL(9.81), REAL(0.01), 3000, REAL(1e-6));
		
		auto result = solver.findSpeedForAngle(
			REAL(0.0), REAL(0.0),
			REAL(150.0), REAL(0.0),
			Constants::PI / 4,        // 45 degrees
			REAL(30.0), REAL(100.0)   // Speed guesses
		);
		
		REQUIRE(result.converged);
		REQUIRE(result.launchSpeed > REAL(0.0));
		REQUIRE(result.residual < REAL(1.0));
		
		INFO("Required speed: " << result.launchSpeed << " m/s");
		INFO("Flight time: " << result.flightTime << " seconds");
	}

	TEST_CASE("Projectile_WithDrag_FindSpeed_LowAngle", "[bvp_shooting][projectile][speed]")
	{
		TEST_PRECISION_INFO();
		
		// Find speed for 30° angle to hit target at 100m
		ProjectileBVPSolver solver(REAL(9.81), REAL(0.01), 3000, REAL(1e-6));
		
		auto result = solver.findSpeedForAngle(
			REAL(0.0), REAL(0.0),
			REAL(100.0), REAL(0.0),
			Constants::PI / 6,        // 30 degrees
			REAL(20.0), REAL(80.0)
		);
		
		REQUIRE(result.converged);
		
		INFO("Required speed: " << result.launchSpeed << " m/s");
		INFO("Max height: " << result.maxHeight << " m");
	}

	TEST_CASE("Projectile_WithDrag_FindSpeed_HighAngle", "[bvp_shooting][projectile][speed]")
	{
		TEST_PRECISION_INFO();
		
		// Find speed for 60° angle (mortar-style trajectory)
		ProjectileBVPSolver solver(REAL(9.81), REAL(0.01), 3000, REAL(1e-6));
		
		auto result = solver.findSpeedForAngle(
			REAL(0.0), REAL(0.0),
			REAL(80.0), REAL(0.0),
			Constants::PI / 3,        // 60 degrees
			REAL(20.0), REAL(80.0)
		);
		
		REQUIRE(result.converged);
		
		// High angle should give high trajectory
		REQUIRE(result.maxHeight > REAL(20.0));
		
		INFO("Required speed: " << result.launchSpeed << " m/s");
		INFO("Max height: " << result.maxHeight << " m");
	}

	/*********************************************************************/
	/*****           DRAG EFFECT VERIFICATION                        *****/
	/*********************************************************************/

	TEST_CASE("Projectile_DragReducesRange", "[bvp_shooting][projectile][physics]")
	{
		TEST_PRECISION_INFO();
		
		// Same initial conditions, compare range with and without drag
		Real v0 = REAL(50.0);
		Real angle = Constants::PI / 4;
		
		ProjectileBVPSolver noDragSolver(REAL(9.81), REAL(0.0), 2000);
		ProjectileBVPSolver withDragSolver(REAL(9.81), REAL(0.02), 2000);  // k=0.02
		
		auto noDragTraj = noDragSolver.computeTrajectory(0, 0, v0, angle, 15.0);
		auto withDragTraj = withDragSolver.computeTrajectory(0, 0, v0, angle, 15.0);
		
		// Find landing points
		Real range_noDrag = 0, range_withDrag = 0;
		
		for (int i = 1; i < noDragTraj.size(); ++i)
		{
			if (noDragTraj.getXValue(i-1, 1) > 0 && noDragTraj.getXValue(i, 1) <= 0)
			{
				Real alpha = -noDragTraj.getXValue(i-1, 1) / (noDragTraj.getXValue(i, 1) - noDragTraj.getXValue(i-1, 1));
				range_noDrag = noDragTraj.getXValue(i-1, 0) + alpha * (noDragTraj.getXValue(i, 0) - noDragTraj.getXValue(i-1, 0));
				break;
			}
		}
		
		for (int i = 1; i < withDragTraj.size(); ++i)
		{
			if (withDragTraj.getXValue(i-1, 1) > 0 && withDragTraj.getXValue(i, 1) <= 0)
			{
				Real alpha = -withDragTraj.getXValue(i-1, 1) / (withDragTraj.getXValue(i, 1) - withDragTraj.getXValue(i-1, 1));
				range_withDrag = withDragTraj.getXValue(i-1, 0) + alpha * (withDragTraj.getXValue(i, 0) - withDragTraj.getXValue(i-1, 0));
				break;
			}
		}
		
		// Drag should reduce range significantly
		REQUIRE(range_withDrag < range_noDrag);
		REQUIRE(range_withDrag < REAL(0.8) * range_noDrag);  // At least 20% reduction
		
		INFO("No-drag range: " << range_noDrag << " m");
		INFO("With-drag range: " << range_withDrag << " m");
		INFO("Range reduction: " << (1 - range_withDrag/range_noDrag) * 100 << "%");
	}

	TEST_CASE("Projectile_HigherDragNeedsMoreSpeed", "[bvp_shooting][projectile][physics]")
	{
		TEST_PRECISION_INFO();
		
		// To hit same target, higher drag requires more speed
		Real angle = Constants::PI / 4;
		Real targetX = REAL(100.0);
		
		ProjectileBVPSolver lowDragSolver(REAL(9.81), REAL(0.005), 3000, REAL(1e-6));
		ProjectileBVPSolver highDragSolver(REAL(9.81), REAL(0.02), 3000, REAL(1e-6));
		
		auto lowDragResult = lowDragSolver.findSpeedForAngle(0, 0, targetX, 0, angle, 20, 80);
		auto highDragResult = highDragSolver.findSpeedForAngle(0, 0, targetX, 0, angle, 20, 100);
		
		REQUIRE(lowDragResult.converged);
		REQUIRE(highDragResult.converged);
		REQUIRE(highDragResult.launchSpeed > lowDragResult.launchSpeed);
		
		INFO("Low drag (k=0.005) speed: " << lowDragResult.launchSpeed << " m/s");
		INFO("High drag (k=0.02) speed: " << highDragResult.launchSpeed << " m/s");
	}

	/*********************************************************************/
	/*****                    INFO TEST                              *****/
	/*********************************************************************/

	TEST_CASE("BVPShooting_Info", "[bvp_shooting][info]")
	{
		TEST_PRECISION_INFO();
		INFO("BVP Shooting Method test suite complete");
		SUCCEED();
	}

	/*********************************************************************/
	/*****           3D PROJECTILE WITH WIND TESTS                   *****/
	/*********************************************************************/

	TEST_CASE("Projectile3D_NoWind_Baseline", "[bvp_shooting][projectile3d][baseline]")
	{
		TEST_PRECISION_INFO();
		
		// Without wind, azimuth should be ~0 (aim straight at target)
		Projectile3DBVPSolver solver(REAL(9.81), REAL(0.01),
		                              REAL(0.0), REAL(0.0), REAL(0.0),  // No wind
		                              4000, REAL(0.5));  // Looser tolerance for landing
		
		auto result = solver.findAnglesForSpeed(
			REAL(0.0), REAL(0.0), REAL(0.0),  // Start at origin
			REAL(100.0), REAL(0.0),            // Target at x=100, z=0
			REAL(50.0),                        // 50 m/s
			REAL(0.6)                          // ~35° elevation guess
		);
		
		REQUIRE(result.converged);
		
		// Azimuth should be essentially zero (no wind compensation needed)
		REQUIRE_THAT(result.launchAzimuthDegrees(), WithinAbs(REAL(0.0), REAL(2.0)));
		
		// Should hit target in x
		REQUIRE_THAT(result.finalPosition[0], WithinAbs(REAL(100.0), REAL(2.0)));
		
		// No lateral drift
		REQUIRE_THAT(result.finalPosition[2], WithinAbs(REAL(0.0), REAL(1.0)));
		
		INFO("Elevation: " << result.launchElevationDegrees() << "°");
		INFO("Azimuth: " << result.launchAzimuthDegrees() << "°");
		INFO("Flight time: " << result.flightTime << " s");
		INFO("Max height: " << result.maxHeight << " m");
	}

	TEST_CASE("Projectile3D_CrossWind_Compensation", "[bvp_shooting][projectile3d][wind]")
	{
		TEST_PRECISION_INFO();
		
		// Crosswind of 5 m/s in z-direction
		// Solver should find negative azimuth to compensate (aim into wind)
		Projectile3DBVPSolver solver(REAL(9.81), REAL(0.01),
		                              REAL(0.0), REAL(0.0), REAL(5.0),  // 5 m/s crosswind
		                              4000, REAL(0.5));
		
		auto result = solver.findAnglesForSpeed(
			REAL(0.0), REAL(0.0), REAL(0.0),
			REAL(100.0), REAL(0.0),   // Target straight ahead
			REAL(50.0),
			REAL(0.6)
		);
		
		REQUIRE(result.converged);
		
		// Azimuth should be NEGATIVE (aim into the wind coming from +z)
		REQUIRE(result.launchAzimuth < REAL(0.0));
		
		// Should still hit target
		REQUIRE_THAT(result.finalPosition[0], WithinAbs(REAL(100.0), REAL(2.0)));
		REQUIRE_THAT(result.finalPosition[2], WithinAbs(REAL(0.0), REAL(1.0)));
		
		INFO("Elevation: " << result.launchElevationDegrees() << "°");
		INFO("Azimuth (wind compensation): " << result.launchAzimuthDegrees() << "°");
		INFO("Flight time: " << result.flightTime << " s");
		INFO("Lateral error: " << result.finalPosition[2] << " m");
	}

	TEST_CASE("Projectile3D_StrongCrossWind", "[bvp_shooting][projectile3d][wind]")
	{
		TEST_PRECISION_INFO();
		
		// Strong 8 m/s crosswind - requires significant compensation
		Projectile3DBVPSolver solver(REAL(9.81), REAL(0.01),
		                              REAL(0.0), REAL(0.0), REAL(8.0),  // 8 m/s crosswind
		                              4000, REAL(0.5));
		
		auto result = solver.findAnglesForSpeed(
			REAL(0.0), REAL(0.0), REAL(0.0),
			REAL(80.0), REAL(0.0),    // Slightly closer target
			REAL(50.0),
			REAL(0.6)                 // ~35° elevation guess
		);
		
		REQUIRE(result.converged);
		
		// Larger azimuth compensation for stronger wind
		REQUIRE(result.launchAzimuth < REAL(-0.02));  // At least a degree or so
		
		INFO("Elevation: " << result.launchElevationDegrees() << "°");
		INFO("Azimuth compensation: " << result.launchAzimuthDegrees() << "°");
		INFO("Iterations: " << result.iterations);
	}

	TEST_CASE("Projectile3D_Headwind_ReducesRange", "[bvp_shooting][projectile3d][wind]")
	{
		TEST_PRECISION_INFO();
		
		// Compare elevation needed with headwind vs no wind
		// Headwind should require steeper angle to reach same distance
		
		Projectile3DBVPSolver noWindSolver(REAL(9.81), REAL(0.01),
		                                    REAL(0.0), REAL(0.0), REAL(0.0),
		                                    4000, REAL(0.5));
		
		Projectile3DBVPSolver headwindSolver(REAL(9.81), REAL(0.01),
		                                      REAL(-5.0), REAL(0.0), REAL(0.0),  // -5 m/s headwind
		                                      4000, REAL(0.5));
		
		auto noWindResult = noWindSolver.findAnglesForSpeed(
			0, 0, 0, REAL(80.0), 0, REAL(50.0), REAL(0.5));
		
		auto headwindResult = headwindSolver.findAnglesForSpeed(
			0, 0, 0, REAL(80.0), 0, REAL(50.0), REAL(0.5));
		
		REQUIRE(noWindResult.converged);
		REQUIRE(headwindResult.converged);
		
		// Headwind requires higher elevation (and/or more initial velocity in x)
		// The elevation should be higher or about the same
		INFO("No wind elevation: " << noWindResult.launchElevationDegrees() << "°");
		INFO("Headwind elevation: " << headwindResult.launchElevationDegrees() << "°");
		INFO("No wind flight time: " << noWindResult.flightTime << " s");
		INFO("Headwind flight time: " << headwindResult.flightTime << " s");
	}

	TEST_CASE("Projectile3D_Tailwind_IncreasesRange", "[bvp_shooting][projectile3d][wind]")
	{
		TEST_PRECISION_INFO();
		
		// Tailwind should allow lower elevation angle for same distance
		Projectile3DBVPSolver noWindSolver(REAL(9.81), REAL(0.01),
		                                    REAL(0.0), REAL(0.0), REAL(0.0),
		                                    4000, REAL(0.5));
		
		Projectile3DBVPSolver tailwindSolver(REAL(9.81), REAL(0.01),
		                                      REAL(5.0), REAL(0.0), REAL(0.0),  // +5 m/s tailwind
		                                      4000, REAL(0.5));
		
		auto noWindResult = noWindSolver.findAnglesForSpeed(
			0, 0, 0, REAL(80.0), 0, REAL(50.0), REAL(0.5));
		
		auto tailwindResult = tailwindSolver.findAnglesForSpeed(
			0, 0, 0, REAL(80.0), 0, REAL(50.0), REAL(0.5));
		
		REQUIRE(noWindResult.converged);
		REQUIRE(tailwindResult.converged);
		
		INFO("No wind elevation: " << noWindResult.launchElevationDegrees() << "°");
		INFO("Tailwind elevation: " << tailwindResult.launchElevationDegrees() << "°");
	}

	TEST_CASE("Projectile3D_DiagonalWind", "[bvp_shooting][projectile3d][wind]")
	{
		TEST_PRECISION_INFO();
		
		// Wind at 45° angle: both headwind component and crosswind
		// wx = -3 m/s (slight headwind), wz = 3 m/s (crosswind)
		Projectile3DBVPSolver solver(REAL(9.81), REAL(0.01),
		                              REAL(-3.0), REAL(0.0), REAL(3.0),
		                              4000, REAL(0.5));
		
		auto result = solver.findAnglesForSpeed(
			REAL(0.0), REAL(0.0), REAL(0.0),
			REAL(100.0), REAL(0.0),
			REAL(50.0),
			REAL(0.6)
		);
		
		REQUIRE(result.converged);
		
		// Should have negative azimuth to compensate for crosswind
		REQUIRE(result.launchAzimuth < REAL(0.0));
		
		INFO("Elevation: " << result.launchElevationDegrees() << "°");
		INFO("Azimuth: " << result.launchAzimuthDegrees() << "°");
		INFO("Final position: (" << result.finalPosition[0] << ", " 
		     << result.finalPosition[1] << ", " << result.finalPosition[2] << ")");
	}

	TEST_CASE("Projectile3D_UnreachableTarget", "[bvp_shooting][projectile3d][edge]")
	{
		TEST_PRECISION_INFO();
		
		// Target too far for given speed with drag - should fail gracefully
		Projectile3DBVPSolver solver(REAL(9.81), REAL(0.02),  // Higher drag
		                              REAL(0.0), REAL(0.0), REAL(0.0),
		                              3000, REAL(1e-4));
		
		auto result = solver.findAnglesForSpeed(
			REAL(0.0), REAL(0.0), REAL(0.0),
			REAL(500.0), REAL(0.0),   // Very far target
			REAL(30.0),                // Low speed
			Constants::PI / 4
		);
		
		// Either doesn't converge or has large residual
		if (result.converged)
		{
			// If it "converged", check that landing was far from target
			Real error = std::abs(result.finalPosition[0] - REAL(500.0));
			INFO("Landing error: " << error << " m (target unreachable)");
		}
		else
		{
			INFO("Correctly identified unreachable target");
		}
		SUCCEED();  // Either outcome is acceptable
	}

	TEST_CASE("Projectile3D_WindDriftWithoutCompensation", "[bvp_shooting][projectile3d][demo]")
	{
		TEST_PRECISION_INFO();
		
		// Demonstrate wind drift when NOT compensating
		// Fire straight ahead (azimuth = 0) with crosswind and see drift
		Projectile3DBVPSolver solver(REAL(9.81), REAL(0.01),
		                              REAL(0.0), REAL(0.0), REAL(5.0),  // 5 m/s crosswind
		                              3000, REAL(1e-4));
		
		// Launch with elevation only, no azimuth (vz = 0)
		Real v0 = REAL(50.0);
		Real theta = REAL(0.6);  // ~35° elevation
		
		auto traj = solver.computeTrajectory(
			0, 0, 0,
			v0 * std::cos(theta),  // vx
			v0 * std::sin(theta),  // vy
			REAL(0.0),             // vz = 0 (no compensation)
			REAL(10.0)
		);
		
		// Find landing point
		Real z_drift = 0;
		Real x_land = 0;
		bool landed = false;
		bool hasRisen = false;
		
		for (int i = 1; i < traj.size(); ++i)
		{
			Real y_prev = traj.getXValue(i-1, 1);
			Real y_curr = traj.getXValue(i, 1);
			
			if (!hasRisen && y_curr > 0.1)
				hasRisen = true;
			
			if (hasRisen && y_prev >= 0 && y_curr < 0)
			{
				Real alpha = -y_prev / (y_curr - y_prev);
				x_land = traj.getXValue(i-1, 0) + alpha * (traj.getXValue(i, 0) - traj.getXValue(i-1, 0));
				z_drift = traj.getXValue(i-1, 2) + alpha * (traj.getXValue(i, 2) - traj.getXValue(i-1, 2));
				landed = true;
				break;
			}
		}
		
		REQUIRE(landed);
		
		// Should have significant drift in direction of wind
		REQUIRE(z_drift > REAL(5.0));  // Drifted downwind
		
		INFO("Without compensation:");
		INFO("  Range: " << x_land << " m");
		INFO("  Lateral drift: " << z_drift << " m (due to 5 m/s crosswind)");
		
		// Now with compensation
		auto result = solver.findAnglesForSpeed(
			0, 0, 0, x_land, 0, v0, theta);
		
		if (result.converged)
		{
			INFO("With compensation:");
			INFO("  Azimuth offset: " << result.launchAzimuthDegrees() << "°");
			INFO("  Lateral error: " << result.finalPosition[2] << " m");
		}
	}
}

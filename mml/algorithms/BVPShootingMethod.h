///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        BVPShootingMethod.h                                                 ///
///  Description: Boundary Value Problem solver using the shooting method             ///
///                                                                                   ///
///  The shooting method converts a two-point boundary value problem into an         ///
///  initial value problem by "shooting" from one boundary and adjusting initial     ///
///  conditions until the other boundary condition is satisfied.                      ///
///                                                                                   ///
///  For a BVP: y'' = f(t, y, y') with y(a) = α and y(b) = β                         ///
///  We solve the IVP: y'' = f(t, y, y') with y(a) = α and y'(a) = s (unknown)       ///
///  Then find s such that y(b) = β using root-finding.                              ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_BVP_SHOOTING_METHOD_H
#define MML_BVP_SHOOTING_METHOD_H

#include "MMLBase.h"

#include "interfaces/IODESystem.h"

#include "base/Vector.h"
#include "base/Matrix.h"
#include "base/ODESystemSolution.h"

#include <vector>

#include "algorithms/ODESystemSolver.h"
#include "algorithms/ODESystemStepCalculators.h"
#include "algorithms/RootFinding.h"

namespace MML
{
	///////////////////////////////////////////////////////////////////////////////////////////
	///                         BVP SHOOTING METHOD RESULT                                  ///
	///////////////////////////////////////////////////////////////////////////////////////////
	
	/**
	 * @brief Result of a shooting method BVP solve
	 * 
	 * Contains the solution trajectory, optimal shooting parameter(s),
	 * and convergence information.
	 */
	struct BVPShootingResult
	{
		bool converged = false;          ///< Did the solver converge?
		int iterations = 0;              ///< Number of shooting iterations
		Real residual = 0.0;             ///< Final boundary condition residual
		Vector<Real> shootingParams;     ///< Optimal shooting parameter(s)
		ODESystemSolution solution;      ///< Full trajectory solution
		
		BVPShootingResult() : solution(0.0, 1.0, 2, 100) {}  // Default: 2D system, 100 steps
	};

	///////////////////////////////////////////////////////////////////////////////////////////
	///                         SINGLE SHOOTING SOLVER                                      ///
	///////////////////////////////////////////////////////////////////////////////////////////
	
	/**
	 * @brief Solve two-point BVP using single shooting method
	 * 
	 * ALGORITHM:
	 * Given an n-th order ODE as a first-order system of n equations:
	 *   dy/dt = f(t, y), with y ∈ ℝⁿ
	 * 
	 * Boundary conditions:
	 *   - Some components of y(t0) are known (initial boundary)
	 *   - Some components of y(t1) are specified (final boundary)
	 * 
	 * The shooting method:
	 * 1. Guess values for unknown initial conditions
	 * 2. Integrate the IVP from t0 to t1
	 * 3. Compare final state with desired boundary condition
	 * 4. Adjust guess using root-finding (Brent's method for 1D)
	 * 5. Repeat until boundary condition is satisfied
	 * 
	 * COMMON USE CASE: Second-order ODE
	 * For y'' = f(t, y, y'), we have state [y, y'] = [y0, y1]
	 * - Given y(t0) = α and y(t1) = β
	 * - Shoot with unknown y'(t0) = s
	 * - Find s such that y(t1) = β
	 * 
	 * @tparam StepCalc ODE step calculator type (default: RK4)
	 */
	template<typename StepCalc = RungeKutta4_StepCalculator>
	class BVPShootingSolver
	{
	private:
		const IODESystem& _system;
		StepCalc _stepCalc;
		int _numSteps;
		Real _tolerance;
		int _maxIterations;

	public:
		/**
		 * @brief Construct shooting solver for given ODE system
		 * 
		 * @param system      ODE system to solve (first-order system form)
		 * @param numSteps    Number of integration steps (more = more accurate)
		 * @param tolerance   Root-finding tolerance for boundary condition
		 * @param maxIter     Maximum shooting iterations
		 */
		BVPShootingSolver(const IODESystem& system,
		                  int numSteps = 1000,
		                  Real tolerance = 1e-10,
		                  int maxIter = 100)
			: _system(system)
			, _numSteps(numSteps)
			, _tolerance(tolerance)
			, _maxIterations(maxIter)
		{}

		/**
		 * @brief Solve 1D shooting problem (one unknown initial condition)
		 * 
		 * This is the most common case: second-order ODE with position BC at both ends.
		 * 
		 * @param t0                Start time
		 * @param t1                End time
		 * @param knownInitCond     Known initial conditions (partial)
		 * @param shootingIndex     Index of the unknown component to shoot with
		 * @param targetIndex       Index of the component with end boundary condition
		 * @param targetValue       Desired value at t1 for targetIndex component
		 * @param shootGuess1       First guess for shooting parameter
		 * @param shootGuess2       Second guess for shooting parameter (for bracketing)
		 * @return BVPShootingResult with solution and convergence info
		 */
		BVPShootingResult solve1D(Real t0, Real t1,
		                          const Vector<Real>& knownInitCond,
		                          int shootingIndex,
		                          int targetIndex,
		                          Real targetValue,
		                          Real shootGuess1,
		                          Real shootGuess2)
		{
			BVPShootingResult result;
			
			// Create the "defect function" that measures boundary condition error
			// F(s) = y[targetIndex](t1) - targetValue, where s is the shooting parameter
			auto defectFunc = [&](Real s) -> Real {
				Vector<Real> initCond = knownInitCond;
				initCond[shootingIndex] = s;
				
				ODESystemFixedStepSolver solver(_system, _stepCalc);
				auto sol = solver.integrate(initCond, t0, t1, _numSteps);
				
				Vector<Real> finalState = sol.getXValuesAtEnd();
				return finalState[targetIndex] - targetValue;
			};
			
			RealFunctionFromStdFunc F(defectFunc);
			
			try
			{
				// Try to bracket the root
				Real a = shootGuess1, b = shootGuess2;
				bool bracketed = RootFinding::BracketRoot(F, a, b, 50);
				
				if (!bracketed)
				{
					// If bracketing fails, try wider search
					a = shootGuess1 - 10 * std::abs(shootGuess2 - shootGuess1);
					b = shootGuess2 + 10 * std::abs(shootGuess2 - shootGuess1);
					bracketed = RootFinding::BracketRoot(F, a, b, 100);
				}
				
				if (!bracketed)
				{
					result.converged = false;
					result.residual = std::numeric_limits<Real>::max();
					return result;
				}
				
				// Use Brent's method to find the root
				Real shootingParam = RootFinding::FindRootBrent(F, a, b, _tolerance);
				
				// Compute final solution with optimal shooting parameter
				Vector<Real> initCond = knownInitCond;
				initCond[shootingIndex] = shootingParam;
				
				ODESystemFixedStepSolver solver(_system, _stepCalc);
				result.solution = solver.integrate(initCond, t0, t1, _numSteps);
				
				result.converged = true;
				result.shootingParams = Vector<Real>(1);
				result.shootingParams[0] = shootingParam;
				
				Vector<Real> finalState = result.solution.getXValuesAtEnd();
				result.residual = std::abs(finalState[targetIndex] - targetValue);
			}
			catch (const std::exception& e)
			{
				result.converged = false;
				result.residual = std::numeric_limits<Real>::max();
			}
			
			return result;
		}
		
		/**
		 * @brief Convenience method for common case: y(t0)=y0, y(t1)=y1, find y'(t0)
		 * 
		 * For second-order ODE converted to first-order system:
		 *   state = [y, y'] = [position, velocity]
		 * 
		 * @param t0           Start time
		 * @param t1           End time
		 * @param y0           Value of y at t0 (position at start)
		 * @param y1           Value of y at t1 (position at end)
		 * @param velocityGuess1  First guess for initial velocity
		 * @param velocityGuess2  Second guess for initial velocity
		 * @return BVPShootingResult
		 */
		BVPShootingResult solvePositionBVP(Real t0, Real t1,
		                                   Real y0, Real y1,
		                                   Real velocityGuess1, Real velocityGuess2)
		{
			Vector<Real> initCond(2);
			initCond[0] = y0;  // Known position at t0
			initCond[1] = 0;   // Placeholder for velocity (shooting parameter)
			
			return solve1D(t0, t1, initCond, 
			               1,    // Shooting index: velocity (y')
			               0,    // Target index: position (y)
			               y1,   // Target value at t1
			               velocityGuess1, velocityGuess2);
		}
	};

	///////////////////////////////////////////////////////////////////////////////////////////
	///                    N-DIMENSIONAL SHOOTING SOLVER                                   ///
	///////////////////////////////////////////////////////////////////////////////////////////
	
	/**
	 * @brief Configuration for N-dimensional shooting method
	 * 
	 * Specifies which initial conditions are unknown (shooting parameters)
	 * and which final conditions are targets (boundary conditions at t1).
	 * 
	 * EXAMPLE: Third-order ODE y''' = f(t, y, y', y'')
	 * State: [y, y', y''] with indices [0, 1, 2]
	 * If we know y(t0), y(t1), y'(t1), we shoot with y'(t0) and y''(t0):
	 *   unknownIndices = {1, 2}      // y'(t0) and y''(t0) are unknown
	 *   targetIndices  = {0, 1}      // y(t1) and y'(t1) are targets
	 */
	struct ShootingConfig
	{
		std::vector<int> unknownIndices;     ///< Indices of unknown initial conditions
		std::vector<Real> unknownGuesses;    ///< Initial guesses for unknowns
		std::vector<int> targetIndices;      ///< Indices of target final conditions  
		std::vector<Real> targetValues;      ///< Desired values at t1
		
		/// Number of shooting parameters (must equal number of targets)
		int numParams() const { return static_cast<int>(unknownIndices.size()); }
		
		/// Validate configuration
		bool isValid() const {
			return unknownIndices.size() == targetIndices.size() &&
			       unknownIndices.size() == unknownGuesses.size() &&
			       unknownIndices.size() == targetValues.size() &&
			       !unknownIndices.empty();
		}
	};
	
	/**
	 * @brief Result of N-dimensional shooting method
	 */
	struct ShootingResultND
	{
		bool converged = false;          ///< Did the solver converge?
		int iterations = 0;              ///< Number of Newton iterations
		Real residualNorm = 0.0;         ///< Final ||F(s)||₂
		Vector<Real> shootingParams;     ///< Optimal shooting parameters
		Vector<Real> residuals;          ///< Final residual vector
		ODESystemSolution solution;      ///< Full trajectory
		
		ShootingResultND() : solution(0.0, 1.0, 2, 100) {}
	};
	
	/**
	 * @brief N-dimensional shooting method solver
	 * 
	 * Solves BVPs with multiple unknown initial conditions using
	 * Newton-Raphson iteration with numerical Jacobian.
	 * 
	 * ALGORITHM:
	 * Given: dy/dt = f(t,y), some y_i(t0) known, some y_j(t1) specified
	 * 
	 * 1. Form shooting parameter vector s = [y_k(t0) for k in unknownIndices]
	 * 2. Define defect F(s) = [y_j(t1;s) - target_j for j in targetIndices]
	 * 3. Solve F(s) = 0 using Newton-Raphson:
	 *    - Compute Jacobian J_ij = ∂F_i/∂s_j numerically
	 *    - Update: s_{n+1} = s_n - J^{-1} F(s_n)
	 * 4. Iterate until ||F(s)|| < tolerance
	 * 
	 * @tparam StepCalc ODE step calculator type
	 */
	template<typename StepCalc = RungeKutta4_StepCalculator>
	class BVPShootingSolverND
	{
	private:
		const IODESystem& _system;
		StepCalc _stepCalc;
		int _numSteps;
		Real _tolerance;
		int _maxIterations;
		Real _jacobianDelta;    ///< Step size for numerical Jacobian
		
	public:
		/**
		 * @brief Construct N-dimensional shooting solver
		 * 
		 * @param system        ODE system
		 * @param numSteps      Integration steps
		 * @param tolerance     Convergence tolerance for ||F(s)||
		 * @param maxIter       Maximum Newton iterations
		 * @param jacobianDelta Step size for numerical derivatives
		 */
		BVPShootingSolverND(const IODESystem& system,
		                    int numSteps = 1000,
		                    Real tolerance = 1e-10,
		                    int maxIter = 50,
		                    Real jacobianDelta = 1e-8)
			: _system(system)
			, _numSteps(numSteps)
			, _tolerance(tolerance)
			, _maxIterations(maxIter)
			, _jacobianDelta(jacobianDelta)
		{}
		
		/**
		 * @brief Solve BVP with given configuration
		 * 
		 * @param t0           Start time
		 * @param t1           End time  
		 * @param knownInitCond Vector with known initial conditions
		 *                      (values at unknownIndices will be overwritten)
		 * @param config       Shooting configuration
		 * @return ShootingResultND
		 */
		ShootingResultND solve(Real t0, Real t1,
		                       const Vector<Real>& knownInitCond,
		                       const ShootingConfig& config)
		{
			ShootingResultND result;
			
			if (!config.isValid())
			{
				result.converged = false;
				return result;
			}
			
			int n = config.numParams();
			
			// Initialize shooting parameters from guesses
			Vector<Real> s(n);
			for (int i = 0; i < n; ++i)
				s[i] = config.unknownGuesses[i];
			
			// Defect function: computes F(s) = final_state[targets] - target_values
			auto computeDefect = [&](const Vector<Real>& params) -> Vector<Real> {
				Vector<Real> initCond = knownInitCond;
				for (int i = 0; i < n; ++i)
					initCond[config.unknownIndices[i]] = params[i];
				
				ODESystemFixedStepSolver solver(_system, _stepCalc);
				auto sol = solver.integrate(initCond, t0, t1, _numSteps);
				Vector<Real> finalState = sol.getXValuesAtEnd();
				
				Vector<Real> F(n);
				for (int i = 0; i < n; ++i)
					F[i] = finalState[config.targetIndices[i]] - config.targetValues[i];
				
				return F;
			};
			
			// Newton-Raphson iteration
			Vector<Real> F = computeDefect(s);
			Real residualNorm = F.NormL2();
			
			for (int iter = 0; iter < _maxIterations; ++iter)
			{
				result.iterations = iter + 1;
				
				if (residualNorm < _tolerance)
				{
					result.converged = true;
					break;
				}
				
				// Compute Jacobian numerically: J_ij = ∂F_i/∂s_j
				Matrix<Real> J(n, n);
				for (int j = 0; j < n; ++j)
				{
					Vector<Real> s_plus = s;
					s_plus[j] += _jacobianDelta;
					
					Vector<Real> F_plus = computeDefect(s_plus);
					
					for (int i = 0; i < n; ++i)
						J(i, j) = (F_plus[i] - F[i]) / _jacobianDelta;
				}
				
				// Solve J * delta_s = -F for delta_s
				// Using simple Gaussian elimination for small systems
				Vector<Real> delta_s = solveLinearSystem(J, F * (-1.0));
				
				// Update shooting parameters
				s = s + delta_s;
				
				// Recompute defect
				F = computeDefect(s);
				residualNorm = F.NormL2();
			}
			
			// Store results
			result.shootingParams = s;
			result.residuals = F;
			result.residualNorm = residualNorm;
			
			// Compute final solution
			Vector<Real> initCond = knownInitCond;
			for (int i = 0; i < n; ++i)
				initCond[config.unknownIndices[i]] = s[i];
			
			ODESystemFixedStepSolver solver(_system, _stepCalc);
			result.solution = solver.integrate(initCond, t0, t1, _numSteps);
			
			if (!result.converged && residualNorm < _tolerance * 100)
				result.converged = true;  // Close enough
			
			return result;
		}
		
	private:
		/**
		 * @brief Solve small linear system Ax = b using Gaussian elimination
		 */
		Vector<Real> solveLinearSystem(const Matrix<Real>& A, const Vector<Real>& b)
		{
			int n = b.size();
			Matrix<Real> Aug(n, n + 1);
			
			// Build augmented matrix [A|b]
			for (int i = 0; i < n; ++i)
			{
				for (int j = 0; j < n; ++j)
					Aug(i, j) = A(i, j);
				Aug(i, n) = b[i];
			}
			
			// Forward elimination with partial pivoting
			for (int k = 0; k < n; ++k)
			{
				// Find pivot
				int maxRow = k;
				Real maxVal = std::abs(Aug(k, k));
				for (int i = k + 1; i < n; ++i)
				{
					if (std::abs(Aug(i, k)) > maxVal)
					{
						maxVal = std::abs(Aug(i, k));
						maxRow = i;
					}
				}
				
				// Swap rows
				if (maxRow != k)
				{
					for (int j = k; j <= n; ++j)
						std::swap(Aug(k, j), Aug(maxRow, j));
				}
				
				// Check for singularity
				if (std::abs(Aug(k, k)) < 1e-14)
					return Vector<Real>(n);  // Return zeros if singular
				
				// Eliminate below
				for (int i = k + 1; i < n; ++i)
				{
					Real factor = Aug(i, k) / Aug(k, k);
					for (int j = k; j <= n; ++j)
						Aug(i, j) -= factor * Aug(k, j);
				}
			}
			
			// Back substitution
			Vector<Real> x(n);
			for (int i = n - 1; i >= 0; --i)
			{
				x[i] = Aug(i, n);
				for (int j = i + 1; j < n; ++j)
					x[i] -= Aug(i, j) * x[j];
				x[i] /= Aug(i, i);
			}
			
			return x;
		}
	};

}  // namespace MML

#endif  // MML_BVP_SHOOTING_METHOD_H

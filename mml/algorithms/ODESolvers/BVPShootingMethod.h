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
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_BVP_SHOOTING_METHOD_H
#define MML_BVP_SHOOTING_METHOD_H

#include "mml/MMLBase.h"

#include "mml/core/AlgorithmTypes.h"
#include "mml/interfaces/IODESystem.h"

#include "mml/base/Vector/Vector.h"
#include "mml/base/Matrix/Matrix.h"
#include "mml/base/ODESystemSolution.h"

#include <vector>

#include "mml/algorithms/ODESolvers/ODEFixedStepIntegrators.h"
#include "mml/algorithms/ODESolvers/ODESystemStepCalculators.h"
#include "mml/algorithms/RootFinding.h"

namespace MML {
	///////////////////////////////////////////////////////////////////////////////////////////
	///                         BVP SHOOTING METHOD FAILURE REASONS                          ///
	///////////////////////////////////////////////////////////////////////////////////////////
	/// @brief Enumeration of possible failure reasons for BVP shooting method
	enum class BVPFailureReason {
		None,                      ///< No failure (converged successfully)
		BracketingFailed,          ///< Could not bracket the root
		RootFindingFailed,         ///< Root finding algorithm failed
		MaxIterationsExceeded,     ///< Exceeded maximum Newton iterations
		SingularJacobian,          ///< Jacobian matrix became singular
		IntegrationFailed,         ///< ODE integration failed
		InvalidConfiguration,      ///< Invalid shooting configuration
		NumericalInstability       ///< Numerical instability detected (NaN/Inf)
	};

	/// @brief Convert failure reason to human-readable string
	inline const char* BVPFailureReasonToString(BVPFailureReason reason) {
		switch (reason) {
			case BVPFailureReason::None:                  return "None";
			case BVPFailureReason::BracketingFailed:      return "BracketingFailed";
			case BVPFailureReason::RootFindingFailed:     return "RootFindingFailed";
			case BVPFailureReason::MaxIterationsExceeded: return "MaxIterationsExceeded";
			case BVPFailureReason::SingularJacobian:      return "SingularJacobian";
			case BVPFailureReason::IntegrationFailed:     return "IntegrationFailed";
			case BVPFailureReason::InvalidConfiguration:  return "InvalidConfiguration";
			case BVPFailureReason::NumericalInstability:  return "NumericalInstability";
			default:                                      return "Unknown";
		}
	}

	///////////////////////////////////////////////////////////////////////////////////////////
	///                         BVP SHOOTING METHOD RESULT                                  ///
	///////////////////////////////////////////////////////////////////////////////////////////
	/// @brief Result of a shooting method BVP solve
	///
	/// Contains the solution trajectory, optimal shooting parameter(s),
	/// convergence information, and detailed diagnostics on failure.
	struct BVPShootingResult {
		bool converged = false;		              ///< Did the solver converge?
		BVPFailureReason failureReason = BVPFailureReason::None; ///< Reason for failure (if !converged)
		std::string failureDetails;               ///< Additional failure details/message
		int iterations = 0;			              ///< Number of shooting iterations
		Real residual = 0.0;		              ///< Final boundary condition residual
		Vector<Real> shootingParams;              ///< Optimal shooting parameter(s)
		ODESystemSolution solution;	              ///< Full trajectory solution

		// Enhanced diagnostic fields (API Standardization Phase 6)
		std::string algorithm_name;               ///< Name of the algorithm used
		AlgorithmStatus status = AlgorithmStatus::Success;  ///< Algorithm termination status
		double elapsed_time_ms = 0;               ///< Execution time in milliseconds

		BVPShootingResult()
			: solution(0.0, 1.0, 2, 100) {} // Default: 2D system, 100 steps
		
		/// @brief Check if result represents a successful solve
		bool success() const { return converged && failureReason == BVPFailureReason::None; }
		
		/// @brief Get human-readable failure description
		std::string getFailureDescription() const {
			if (converged) return "Success";
			std::string desc = BVPFailureReasonToString(failureReason);
			if (!failureDetails.empty()) {
				desc += ": ";
				desc += failureDetails;
			}
			return desc;
		}
	};

	///////////////////////////////////////////////////////////////////////////////////////////
	///                         SINGLE SHOOTING SOLVER                                      ///
	///////////////////////////////////////////////////////////////////////////////////////////
	/// @brief Solve two-point BVP using single shooting method
	///
	/// ALGORITHM:
	/// Given an n-th order ODE as a first-order system of n equations:
	/// dy/dt = f(t, y), with y ∈ ℝⁿ
	///
	/// Boundary conditions:
	/// - Some components of y(t0) are known (initial boundary)
	/// - Some components of y(t1) are specified (final boundary)
	///
	/// The shooting method:
	/// 1. Guess values for unknown initial conditions
	/// 2. Integrate the IVP from t0 to t1
	/// 3. Compare final state with desired boundary condition
	/// 4. Adjust guess using root-finding (Brent's method for 1D)
	/// 5. Repeat until boundary condition is satisfied
	///
	/// COMMON USE CASE: Second-order ODE
	/// For y'' = f(t, y, y'), we have state [y, y'] = [y0, y1]
	/// - Given y(t0) = α and y(t1) = β
	/// - Shoot with unknown y'(t0) = s
	/// - Find s such that y(t1) = β
	///
	/// @tparam StepCalc ODE step calculator type (default: RK4)
	template<typename StepCalc = RungeKutta4_StepCalculator>
	class BVPShootingSolver {
	private:
		const IODESystem& _system;
		StepCalc _stepCalc;
		int _numSteps;
		Real _tolerance;
		int _maxIterations;

	public:
		/// @brief Construct shooting solver for given ODE system
		///
		/// @param system      ODE system to solve (first-order system form)
		/// @param numSteps    Number of integration steps (more = more accurate)
		/// @param tolerance   Root-finding tolerance for boundary condition
		/// @param maxIter     Maximum shooting iterations
		BVPShootingSolver(const IODESystem& system, int numSteps = 1000, Real tolerance = 1e-10, int maxIter = 100)
			: _system(system)
			, _numSteps(numSteps)
			, _tolerance(tolerance)
			, _maxIterations(maxIter) {}

		/// @brief Solve 1D shooting problem (one unknown initial condition)
		///
		/// This is the most common case: second-order ODE with position BC at both ends.
		///
		/// @param t0                Start time
		/// @param t1                End time
		/// @param knownInitCond     Known initial conditions (partial)
		/// @param shootingIndex     Index of the unknown component to shoot with
		/// @param targetIndex       Index of the component with end boundary condition
		/// @param targetValue       Desired value at t1 for targetIndex component
		/// @param shootGuess1       First guess for shooting parameter
		/// @param shootGuess2       Second guess for shooting parameter (for bracketing)
		/// @return BVPShootingResult with solution and convergence info
		BVPShootingResult solve1D(Real t0, Real t1, const Vector<Real>& knownInitCond, int shootingIndex, int targetIndex, Real targetValue,
								  Real shootGuess1, Real shootGuess2) {
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

			try {
				// Try to bracket the root
				Real a = shootGuess1, b = shootGuess2;
				bool bracketed = RootFinding::BracketRoot(F, a, b, 50);

				if (!bracketed) {
					// If bracketing fails, try wider search
					a = shootGuess1 - 10 * std::abs(shootGuess2 - shootGuess1);
					b = shootGuess2 + 10 * std::abs(shootGuess2 - shootGuess1);
					bracketed = RootFinding::BracketRoot(F, a, b, 100);
				}

				if (!bracketed) {
					result.converged = false;
					result.failureReason = BVPFailureReason::BracketingFailed;
					result.failureDetails = "Could not bracket root with initial guesses [" + 
						std::to_string(shootGuess1) + ", " + std::to_string(shootGuess2) + 
						"] even after expanding range";
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
				result.failureReason = BVPFailureReason::None;
			} catch (const std::exception& e) {
				result.converged = false;
				result.failureReason = BVPFailureReason::RootFindingFailed;
				result.failureDetails = e.what();
				result.residual = std::numeric_limits<Real>::max();
			}

			return result;
		}

		/// @brief Convenience method for common case: y(t0)=y0, y(t1)=y1, find y'(t0)
		///
		/// For second-order ODE converted to first-order system:
		/// state = [y, y'] = [position, velocity]
		///
		/// @param t0           Start time
		/// @param t1           End time
		/// @param y0           Value of y at t0 (position at start)
		/// @param y1           Value of y at t1 (position at end)
		/// @param velocityGuess1  First guess for initial velocity
		/// @param velocityGuess2  Second guess for initial velocity
		/// @return BVPShootingResult
		BVPShootingResult solvePositionBVP(Real t0, Real t1, Real y0, Real y1, Real velocityGuess1, Real velocityGuess2) {
			Vector<Real> initCond(2);
			initCond[0] = y0; // Known position at t0
			initCond[1] = 0;  // Placeholder for velocity (shooting parameter)

			return solve1D(t0, t1, initCond,
						   1,  // Shooting index: velocity (y')
						   0,  // Target index: position (y)
						   y1, // Target value at t1
						   velocityGuess1, velocityGuess2);
		}
	};

	///////////////////////////////////////////////////////////////////////////////////////////
	///                    N-DIMENSIONAL SHOOTING SOLVER                                   ///
	///////////////////////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////////
	///                       BVPSolverConfig                               ///
	///////////////////////////////////////////////////////////////////////////
	/// @brief Configuration for BVP solver algorithm parameters
	///
	/// Standardized config object following API Standardization Phase 6 pattern.
	struct BVPSolverConfig {
		Real tolerance = 1e-10;          ///< Convergence tolerance
		int max_iterations = 50;         ///< Maximum Newton iterations
		int num_integration_steps = 1000; ///< Number of ODE integration steps
		Real jacobian_delta = 1e-8;      ///< Step size for numerical Jacobian
		bool verbose = false;            ///< Enable verbose output

		/// Default constructor
		BVPSolverConfig() = default;

		/// Constructor with tolerance
		BVPSolverConfig(Real tol, int max_iter = 50)
			: tolerance(tol)
			, max_iterations(max_iter) {}

		/// Factory: High precision configuration
		static BVPSolverConfig HighPrecision() {
			BVPSolverConfig cfg;
			cfg.tolerance = 1e-14;
			cfg.max_iterations = 100;
			cfg.num_integration_steps = 5000;
			cfg.jacobian_delta = 1e-10;
			return cfg;
		}

		/// Factory: Fast configuration
		static BVPSolverConfig Fast() {
			BVPSolverConfig cfg;
			cfg.tolerance = 1e-6;
			cfg.max_iterations = 20;
			cfg.num_integration_steps = 200;
			return cfg;
		}
	};

	/// @brief Configuration for N-dimensional shooting method
	///
	/// Specifies which initial conditions are unknown (shooting parameters)
	/// and which final conditions are targets (boundary conditions at t1).
	///
	/// EXAMPLE: Third-order ODE y''' = f(t, y, y', y'')
	/// State: [y, y', y''] with indices [0, 1, 2]
	/// If we know y(t0), y(t1), y'(t1), we shoot with y'(t0) and y''(t0):
	/// unknownIndices = {1, 2}      // y'(t0) and y''(t0) are unknown
	/// targetIndices  = {0, 1}      // y(t1) and y'(t1) are targets
	struct ShootingConfig {
		std::vector<int> unknownIndices;  ///< Indices of unknown initial conditions
		std::vector<Real> unknownGuesses; ///< Initial guesses for unknowns
		std::vector<int> targetIndices;	  ///< Indices of target final conditions
		std::vector<Real> targetValues;	  ///< Desired values at t1
		/// Number of shooting parameters (must equal number of targets)
		int numParams() const { return static_cast<int>(unknownIndices.size()); }

		/// Validate configuration
		bool isValid() const {
			return unknownIndices.size() == targetIndices.size() && unknownIndices.size() == unknownGuesses.size() &&
				   unknownIndices.size() == targetValues.size() && !unknownIndices.empty();
		}
	};

	/// @brief Result of N-dimensional shooting method
	struct ShootingResultND {
		bool converged = false;		              ///< Did the solver converge?
		BVPFailureReason failureReason = BVPFailureReason::None; ///< Reason for failure (if !converged)
		std::string failureDetails;               ///< Additional failure details/message
		int iterations = 0;			              ///< Number of Newton iterations
		Real residualNorm = 0.0;	              ///< Final ||F(s)||₂
		Vector<Real> shootingParams;              ///< Optimal shooting parameters
		Vector<Real> residuals;		              ///< Final residual vector
		ODESystemSolution solution;	              ///< Full trajectory

		// Enhanced diagnostic fields (API Standardization Phase 6)
		std::string algorithm_name;               ///< Name of the algorithm used
		AlgorithmStatus status = AlgorithmStatus::Success;  ///< Algorithm termination status
		double elapsed_time_ms = 0;               ///< Execution time in milliseconds

		ShootingResultND()
			: solution(0.0, 1.0, 2, 100) {}
		
		/// @brief Check if result represents a successful solve
		bool success() const { return converged && failureReason == BVPFailureReason::None; }
		
		/// @brief Get human-readable failure description
		std::string getFailureDescription() const {
			if (converged) return "Success";
			std::string desc = BVPFailureReasonToString(failureReason);
			if (!failureDetails.empty()) {
				desc += ": ";
				desc += failureDetails;
			}
			return desc;
		}
	};

	/// @brief N-dimensional shooting method solver
	///
	/// Solves BVPs with multiple unknown initial conditions using
	/// Newton-Raphson iteration with numerical Jacobian.
	///
	/// ALGORITHM:
	/// Given: dy/dt = f(t,y), some y_i(t0) known, some y_j(t1) specified
	///
	/// 1. Form shooting parameter vector s = [y_k(t0) for k in unknownIndices]
	/// 2. Define defect F(s) = [y_j(t1;s) - target_j for j in targetIndices]
	/// 3. Solve F(s) = 0 using Newton-Raphson:
	/// - Compute Jacobian J_ij = ∂F_i/∂s_j numerically
	/// - Update: s_{n+1} = s_n - J^{-1} F(s_n)
	/// 4. Iterate until ||F(s)|| < tolerance
	///
	/// @tparam StepCalc ODE step calculator type
	template<typename StepCalc = RungeKutta4_StepCalculator>
	class BVPShootingSolverND {
	private:
		const IODESystem& _system;
		StepCalc _stepCalc;
		int _numSteps;
		Real _tolerance;
		int _maxIterations;
		Real _jacobianDelta; ///< Step size for numerical Jacobian

	public:
		/// @brief Construct N-dimensional shooting solver
		///
		/// @param system        ODE system
		/// @param numSteps      Integration steps
		/// @param tolerance     Convergence tolerance for ||F(s)||
		/// @param maxIter       Maximum Newton iterations
		/// @param jacobianDelta Step size for numerical derivatives
		BVPShootingSolverND(const IODESystem& system, int numSteps = 1000, Real tolerance = 1e-10, int maxIter = 50,
							Real jacobianDelta = 1e-8)
			: _system(system)
			, _numSteps(numSteps)
			, _tolerance(tolerance)
			, _maxIterations(maxIter)
			, _jacobianDelta(jacobianDelta) {}

		/// @brief Construct with standardized config (Phase 6 API)
		/// @param system ODE system to solve
		/// @param config Solver configuration
		BVPShootingSolverND(const IODESystem& system, const BVPSolverConfig& config)
			: _system(system)
			, _numSteps(config.num_integration_steps)
			, _tolerance(config.tolerance)
			, _maxIterations(config.max_iterations)
			, _jacobianDelta(config.jacobian_delta) {}

		/// @brief Solve BVP with given configuration
		/// @param t0           Start time
		/// @param t1           End time
		/// @param knownInitCond Vector with known initial conditions
		/// (values at unknownIndices will be overwritten)
		/// @param config       Shooting configuration
		/// @return ShootingResultND
		ShootingResultND solve(Real t0, Real t1, const Vector<Real>& knownInitCond, const ShootingConfig& config) {
			ShootingResultND result;

			if (!config.isValid()) {
				result.converged = false;
				result.failureReason = BVPFailureReason::InvalidConfiguration;
				result.failureDetails = "Shooting configuration invalid: unknownIndices, targetIndices, "
					"unknownGuesses, and targetValues must all have same non-zero size";
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

			for (int iter = 0; iter < _maxIterations; ++iter) {
				result.iterations = iter + 1;

				if (residualNorm < _tolerance) {
					result.converged = true;
					result.failureReason = BVPFailureReason::None;
					break;
				}
				
				// Check for numerical instability
				if (!std::isfinite(residualNorm)) {
					result.converged = false;
					result.failureReason = BVPFailureReason::NumericalInstability;
					result.failureDetails = "Residual norm became NaN/Inf at iteration " + 
						std::to_string(iter);
					break;
				}

				// Compute Jacobian numerically: J_ij = ∂F_i/∂s_j
				Matrix<Real> J(n, n);
				for (int j = 0; j < n; ++j) {
					Vector<Real> s_plus = s;
					s_plus[j] += _jacobianDelta;

					Vector<Real> F_plus = computeDefect(s_plus);

					for (int i = 0; i < n; ++i)
						J(i, j) = (F_plus[i] - F[i]) / _jacobianDelta;
				}

				// Solve J * delta_s = -F for delta_s
				// Using simple Gaussian elimination for small systems
				auto linearResult = solveLinearSystem(J, F * (-1.0));
				
				if (!linearResult.success) {
					result.converged = false;
					result.failureReason = BVPFailureReason::SingularJacobian;
					result.failureDetails = linearResult.errorMessage;
					break;
				}

				// Update shooting parameters
				s = s + linearResult.solution;

				// Recompute defect
				F = computeDefect(s);
				residualNorm = F.NormL2();
			}
			
			// Check if we exceeded max iterations without converging
			if (!result.converged && result.failureReason == BVPFailureReason::None) {
				result.failureReason = BVPFailureReason::MaxIterationsExceeded;
				result.failureDetails = "Did not converge within " + std::to_string(_maxIterations) + 
					" iterations. Final residual norm: " + std::to_string(residualNorm) + 
					" (tolerance: " + std::to_string(_tolerance) + ")";
			}

			// Store results
			result.shootingParams = s;
			result.residuals = F;
			result.residualNorm = residualNorm;

			// Compute final solution (even for failed solves, provide best-effort trajectory)
			Vector<Real> initCond = knownInitCond;
			for (int i = 0; i < n; ++i)
				initCond[config.unknownIndices[i]] = s[i];

			ODESystemFixedStepSolver solver(_system, _stepCalc);
			result.solution = solver.integrate(initCond, t0, t1, _numSteps);

			return result;
		}

	private:
		/// @brief Result of linear system solve
		struct LinearSolveResult {
			bool success = false;      ///< Was the solve successful?
			Vector<Real> solution;     ///< Solution vector (valid only if success)
			std::string errorMessage;  ///< Error message if !success
		};

		/// @brief Solve small linear system Ax = b using Gaussian elimination
		/// @return LinearSolveResult with success flag and solution or error message
		LinearSolveResult solveLinearSystem(const Matrix<Real>& A, const Vector<Real>& b) {
			LinearSolveResult result;
			int n = b.size();
			Matrix<Real> Aug(n, n + 1);

			// Build augmented matrix [A|b]
			for (int i = 0; i < n; ++i) {
				for (int j = 0; j < n; ++j)
					Aug(i, j) = A(i, j);
				Aug(i, n) = b[i];
			}

			// Forward elimination with partial pivoting
			for (int k = 0; k < n; ++k) {
				// Find pivot
				int maxRow = k;
				Real maxVal = std::abs(Aug(k, k));
				for (int i = k + 1; i < n; ++i) {
					if (std::abs(Aug(i, k)) > maxVal) {
						maxVal = std::abs(Aug(i, k));
						maxRow = i;
					}
				}

				// Swap rows
				if (maxRow != k) {
					for (int j = k; j <= n; ++j)
						std::swap(Aug(k, j), Aug(maxRow, j));
				}

				// Check for singularity
				if (std::abs(Aug(k, k)) < 1e-14) {
					result.success = false;
					result.errorMessage = "Jacobian matrix is singular at row " + std::to_string(k) + 
						" (pivot = " + std::to_string(Aug(k, k)) + ")";
					return result;
				}

				// Eliminate below
				for (int i = k + 1; i < n; ++i) {
					Real factor = Aug(i, k) / Aug(k, k);
					for (int j = k; j <= n; ++j)
						Aug(i, j) -= factor * Aug(k, j);
				}
			}

			// Back substitution
			Vector<Real> x(n);
			for (int i = n - 1; i >= 0; --i) {
				x[i] = Aug(i, n);
				for (int j = i + 1; j < n; ++j)
					x[i] -= Aug(i, j) * x[j];
				x[i] /= Aug(i, i);
				
				// Check for NaN/Inf
				if (!std::isfinite(x[i])) {
					result.success = false;
					result.errorMessage = "Numerical instability during back substitution at index " + 
						std::to_string(i);
					return result;
				}
			}

			result.success = true;
			result.solution = x;
			return result;
		}
	};

} // namespace MML

#endif // MML_BVP_SHOOTING_METHOD_H

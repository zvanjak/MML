#if !defined MML_OPTIMIZATION_COMMON_H
#define MML_OPTIMIZATION_COMMON_H

#include "../../MMLBase.h"
#include "../../base/Vector.h"

#include <vector>
#include <string>
#include <limits>
#include <memory>
#include <functional>

///////////////////////////////////////////////////////////////////////////////////////////
// OptimizationCommon.h - Base infrastructure for flexible optimization framework
// 
// Provides:
//   - OptimizationState: Encapsulates current state during optimization
//   - OptimizationResult: Enhanced result structure with trajectory and diagnostics
//   - ITerminationCriterion: Interface for pluggable termination criteria
//   - IOptimizationObserver: Observer pattern for progress monitoring
//
// This is the foundation for the new optimization framework, designed to be:
//   - Extensible (add new criteria/observers without modifying optimizers)
//   - Flexible (compose criteria, chain observers)
//   - Backward compatible (old APIs still work)
//   - Zero overhead (no cost if features not used)
///////////////////////////////////////////////////////////////////////////////////////////

namespace MML
{

///////////////////////////////////////////////////////////////////////////////////////////
// OPTIMIZATION STATE - Encapsulates current state during optimization
///////////////////////////////////////////////////////////////////////////////////////////

/// @brief Encapsulates the current state of an optimization algorithm
/// @details Contains all relevant information about the current iteration,
///          allowing termination criteria and observers to make decisions
template <typename Real>
struct OptimizationState
{
    // Iteration counters
    int iteration;              ///< Current iteration number (0-based)
    int funcEvals;              ///< Total function evaluations so far
    int gradEvals;              ///< Total gradient evaluations (for gradient-based methods)
    
    // Current solution
    Vector<Real> xCurrent;      ///< Current solution vector
    Real fCurrent;              ///< Current function value
    
    // Best solution so far
    Vector<Real> xBest;         ///< Best solution found so far
    Real fBest;                 ///< Best function value found so far
    
    // Convergence metrics (method-specific)
    Real gradNorm;              ///< Gradient norm (for gradient-based methods)
    Real stepSize;              ///< Last step size taken
    Real simplexSize;           ///< Simplex size (for Nelder-Mead)
    Real directionChange;       ///< Direction change angle (for Powell, CG)
    
    // Stochastic method metrics
    Real temperature;           ///< Temperature (for Simulated Annealing)
    int acceptedMoves;          ///< Number of accepted moves (for SA/MCMC)
    int rejectedMoves;          ///< Number of rejected moves (for SA/MCMC)
    
    // Stagnation tracking
    int iterSinceImprovement;   ///< Iterations without improvement in fBest
    Real fBestPrevious;         ///< Previous best function value (for tracking improvement)
    
    // Timing
    double elapsedTime;         ///< Wall-clock time since start (seconds)
    
    /// @brief Default constructor initializes all fields to safe defaults
    OptimizationState()
        : iteration(0), funcEvals(0), gradEvals(0),
          fCurrent(std::numeric_limits<Real>::max()),
          fBest(std::numeric_limits<Real>::max()),
          gradNorm(0), stepSize(0), simplexSize(0), directionChange(0),
          temperature(0), acceptedMoves(0), rejectedMoves(0),
          iterSinceImprovement(0), 
          fBestPrevious(std::numeric_limits<Real>::max()),
          elapsedTime(0.0)
    {}
    
    /// @brief Check if there was improvement in this iteration
    bool HasImproved() const {
        return fBest < fBestPrevious - std::numeric_limits<Real>::epsilon();
    }
    
    /// @brief Get acceptance rate (for stochastic methods)
    Real GetAcceptanceRate() const {
        int total = acceptedMoves + rejectedMoves;
        return total > 0 ? static_cast<Real>(acceptedMoves) / total : Real(0);
    }
};

///////////////////////////////////////////////////////////////////////////////////////////
// OPTIMIZATION RESULT - Enhanced result structure with trajectory and diagnostics
///////////////////////////////////////////////////////////////////////////////////////////

/// @brief Enhanced optimization result with trajectory and diagnostics
/// @details Replaces the old simple result structs with richer information
template <typename Real>
struct OptimizationResult
{
    // Solution
    Vector<Real> xBest;         ///< Best solution found
    Real fBest;                 ///< Best function value
    
    // Convergence information
    bool converged;             ///< Whether optimization converged successfully
    std::string terminationReason;  ///< Human-readable reason for termination
    
    // Statistics
    int iterations;             ///< Total iterations performed
    int funcEvals;              ///< Total function evaluations
    int gradEvals;              ///< Total gradient evaluations (gradient-based methods)
    int acceptedMoves;          ///< Accepted moves (stochastic methods)
    int rejectedMoves;          ///< Rejected moves (stochastic methods)
    double elapsedTime;         ///< Total wall-clock time (seconds)
    
    // Final state diagnostics
    Real finalGradNorm;         ///< Final gradient norm (gradient methods)
    Real finalStepSize;         ///< Final step size
    Real finalTemperature;      ///< Final temperature (Simulated Annealing)
    Real finalSimplexSize;      ///< Final simplex size (Nelder-Mead)
    
    // Optional trajectory (populated if TrajectoryObserver used)
    std::vector<OptimizationState<Real>> trajectory;
    
    /// @brief Default constructor
    OptimizationResult()
        : fBest(std::numeric_limits<Real>::max()),
          converged(false),
          terminationReason("Not started"),
          iterations(0), funcEvals(0), gradEvals(0),
          acceptedMoves(0), rejectedMoves(0),
          elapsedTime(0.0),
          finalGradNorm(0), finalStepSize(0),
          finalTemperature(0), finalSimplexSize(0)
    {}
    
    /// @brief Get acceptance rate (for stochastic methods)
    Real GetAcceptanceRate() const {
        int total = acceptedMoves + rejectedMoves;
        return total > 0 ? static_cast<Real>(acceptedMoves) / total : Real(0);
    }
    
    /// @brief Check if trajectory was recorded
    bool HasTrajectory() const {
        return !trajectory.empty();
    }
};

///////////////////////////////////////////////////////////////////////////////////////////
// TERMINATION CRITERION INTERFACE
///////////////////////////////////////////////////////////////////////////////////////////

/// @brief Interface for termination criteria
/// @details Allows pluggable, composable stopping conditions for optimization
template <typename Real>
class ITerminationCriterion
{
public:
    virtual ~ITerminationCriterion() = default;
    
    /// @brief Check if optimization should terminate
    /// @param state Current optimization state
    /// @return true if optimization should stop, false to continue
    virtual bool ShouldTerminate(const OptimizationState<Real>& state) const = 0;
    
    /// @brief Get human-readable reason for termination
    /// @return String describing why termination occurred (if applicable)
    virtual std::string GetReason() const = 0;
    
    /// @brief Reset criterion for new optimization run
    /// @details Called before starting a new optimization to clear any state
    virtual void Reset() = 0;
};

///////////////////////////////////////////////////////////////////////////////////////////
// OPTIMIZATION OBSERVER INTERFACE
///////////////////////////////////////////////////////////////////////////////////////////

/// @brief Observer interface for monitoring optimization progress
/// @details Allows callbacks at key points during optimization for:
///          - Progress monitoring and logging
///          - Intermediate result saving (checkpointing)
///          - Custom early stopping logic
///          - GUI updates, plotting, etc.
template <typename Real>
class IOptimizationObserver
{
public:
    virtual ~IOptimizationObserver() = default;
    
    /// @brief Called at the start of optimization
    /// @param state Initial state (iteration 0)
    virtual void OnStart(const OptimizationState<Real>& state) = 0;
    
    /// @brief Called after each iteration
    /// @param state Current state after iteration
    /// @return true to continue optimization, false to abort
    /// @note Returning false allows observers to implement custom early stopping
    virtual bool OnIteration(const OptimizationState<Real>& state) = 0;
    
    /// @brief Called when optimization completes (success or abort)
    /// @param state Final state
    /// @param reason Termination reason string
    virtual void OnComplete(const OptimizationState<Real>& state, 
                           const std::string& reason) = 0;
};

///////////////////////////////////////////////////////////////////////////////////////////
// HELPER FUNCTIONS
///////////////////////////////////////////////////////////////////////////////////////////

/// @brief Create an OptimizationResult from final OptimizationState
/// @param state Final state
/// @param reason Termination reason
/// @param converged Whether optimization converged successfully
/// @return Populated OptimizationResult
template <typename Real>
OptimizationResult<Real> CreateResult(const OptimizationState<Real>& state,
                                       const std::string& reason,
                                       bool converged)
{
    OptimizationResult<Real> result;
    
    result.xBest = state.xBest;
    result.fBest = state.fBest;
    result.converged = converged;
    result.terminationReason = reason;
    
    result.iterations = state.iteration;
    result.funcEvals = state.funcEvals;
    result.gradEvals = state.gradEvals;
    result.acceptedMoves = state.acceptedMoves;
    result.rejectedMoves = state.rejectedMoves;
    result.elapsedTime = state.elapsedTime;
    
    result.finalGradNorm = state.gradNorm;
    result.finalStepSize = state.stepSize;
    result.finalTemperature = state.temperature;
    result.finalSimplexSize = state.simplexSize;
    
    return result;
}

} // namespace MML

#endif // MML_OPTIMIZATION_COMMON_H

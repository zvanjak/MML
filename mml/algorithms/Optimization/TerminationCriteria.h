#if !defined MML_TERMINATION_CRITERIA_H
#define MML_TERMINATION_CRITERIA_H

#include "OptimizationCommon.h"
#include <cmath>
#include <algorithm>

///////////////////////////////////////////////////////////////////////////////////////////
// TerminationCriteria.h - Standard termination criteria for optimization
//
// Provides concrete implementations of ITerminationCriterion:
//   - MaxIterationsCriterion: Stop after N iterations
//   - FunctionToleranceCriterion: Converged when f change is small
//   - GradientNormCriterion: Stop when gradient norm is small
//   - StagnationCriterion: Stop after N iterations without improvement
//   - TimeLimitCriterion: Wall-clock time limit
//   - TargetValueCriterion: Stop when target function value reached
//   - CompositeCriterion: Combine multiple criteria with AND/OR logic
//
// All criteria are:
//   - Stateless or properly reset-able
//   - Thread-safe for const methods
//   - Provide clear termination reasons
///////////////////////////////////////////////////////////////////////////////////////////

namespace MML
{

///////////////////////////////////////////////////////////////////////////////////////////
// MAX ITERATIONS CRITERION
///////////////////////////////////////////////////////////////////////////////////////////

/// @brief Terminate after maximum number of iterations
template <typename Real>
class MaxIterationsCriterion : public ITerminationCriterion<Real>
{
private:
    int _maxIter;
    
public:
    explicit MaxIterationsCriterion(int maxIter) 
        : _maxIter(maxIter) 
    {
        if (maxIter <= 0) {
            throw std::invalid_argument("MaxIterationsCriterion: maxIter must be positive");
        }
    }
    
    bool ShouldTerminate(const OptimizationState<Real>& state) const override {
        return state.iteration >= _maxIter;
    }
    
    std::string GetReason() const override {
        return "Maximum iterations reached (" + std::to_string(_maxIter) + ")";
    }
    
    void Reset() override {
        // Stateless - nothing to reset
    }
    
    int GetMaxIterations() const { return _maxIter; }
};

///////////////////////////////////////////////////////////////////////////////////////////
// FUNCTION TOLERANCE CRITERION
///////////////////////////////////////////////////////////////////////////////////////////

/// @brief Terminate when function value change is below tolerance
/// @details Checks relative change: |f_new - f_old| / (|f_old| + |f_new| + eps) < ftol
template <typename Real>
class FunctionToleranceCriterion : public ITerminationCriterion<Real>
{
private:
    Real _ftol;
    mutable Real _fBestPrev;
    mutable bool _initialized;
    
public:
    explicit FunctionToleranceCriterion(Real ftol) 
        : _ftol(ftol), _fBestPrev(0), _initialized(false)
    {
        if (ftol <= 0) {
            throw std::invalid_argument("FunctionToleranceCriterion: ftol must be positive");
        }
    }
    
    bool ShouldTerminate(const OptimizationState<Real>& state) const override {
        if (!_initialized) {
            _fBestPrev = state.fBest;
            _initialized = true;
            return false;  // Need at least one iteration
        }
        
        // Relative change in function value
        Real fChange = std::abs(_fBestPrev - state.fBest);
        Real fScale = std::abs(_fBestPrev) + std::abs(state.fBest) + Real(1e-10);
        Real relativeChange = fChange / fScale;
        
        bool converged = relativeChange < _ftol;
        
        // Update for next iteration
        _fBestPrev = state.fBest;
        
        return converged;
    }
    
    std::string GetReason() const override {
        return "Function value converged (ftol = " + std::to_string(_ftol) + ")";
    }
    
    void Reset() override {
        _initialized = false;
        _fBestPrev = 0;
    }
    
    Real GetTolerance() const { return _ftol; }
};

///////////////////////////////////////////////////////////////////////////////////////////
// GRADIENT NORM CRITERION
///////////////////////////////////////////////////////////////////////////////////////////

/// @brief Terminate when gradient norm is below threshold
/// @details For gradient-based optimization methods (BFGS, CG, etc.)
template <typename Real>
class GradientNormCriterion : public ITerminationCriterion<Real>
{
private:
    Real _gtol;
    
public:
    explicit GradientNormCriterion(Real gtol) 
        : _gtol(gtol)
    {
        if (gtol <= 0) {
            throw std::invalid_argument("GradientNormCriterion: gtol must be positive");
        }
    }
    
    bool ShouldTerminate(const OptimizationState<Real>& state) const override {
        return state.gradNorm < _gtol;
    }
    
    std::string GetReason() const override {
        return "Gradient norm below threshold (gtol = " + std::to_string(_gtol) + ")";
    }
    
    void Reset() override {
        // Stateless
    }
    
    Real GetTolerance() const { return _gtol; }
};

///////////////////////////////////////////////////////////////////////////////////////////
// STAGNATION CRITERION
///////////////////////////////////////////////////////////////////////////////////////////

/// @brief Terminate when no improvement for N iterations
/// @details Detects when optimization is stuck in a local minimum or plateau
template <typename Real>
class StagnationCriterion : public ITerminationCriterion<Real>
{
private:
    int _maxStagnation;
    
public:
    explicit StagnationCriterion(int maxStagnation) 
        : _maxStagnation(maxStagnation)
    {
        if (maxStagnation <= 0) {
            throw std::invalid_argument("StagnationCriterion: maxStagnation must be positive");
        }
    }
    
    bool ShouldTerminate(const OptimizationState<Real>& state) const override {
        return state.iterSinceImprovement >= _maxStagnation;
    }
    
    std::string GetReason() const override {
        return "Stagnation detected (no improvement for " + 
               std::to_string(_maxStagnation) + " iterations)";
    }
    
    void Reset() override {
        // State tracked in OptimizationState
    }
    
    int GetMaxStagnation() const { return _maxStagnation; }
};

///////////////////////////////////////////////////////////////////////////////////////////
// TIME LIMIT CRITERION
///////////////////////////////////////////////////////////////////////////////////////////

/// @brief Terminate after wall-clock time limit
/// @details Useful for real-time applications or benchmarking
template <typename Real>
class TimeLimitCriterion : public ITerminationCriterion<Real>
{
private:
    double _maxSeconds;
    
public:
    explicit TimeLimitCriterion(double maxSeconds) 
        : _maxSeconds(maxSeconds)
    {
        if (maxSeconds <= 0) {
            throw std::invalid_argument("TimeLimitCriterion: maxSeconds must be positive");
        }
    }
    
    bool ShouldTerminate(const OptimizationState<Real>& state) const override {
        return state.elapsedTime >= _maxSeconds;
    }
    
    std::string GetReason() const override {
        return "Time limit exceeded (" + std::to_string(_maxSeconds) + " seconds)";
    }
    
    void Reset() override {
        // Time tracked in OptimizationState
    }
    
    double GetMaxSeconds() const { return _maxSeconds; }
    double GetMaxMinutes() const { return _maxSeconds / 60.0; }
};

///////////////////////////////////////////////////////////////////////////////////////////
// TARGET VALUE CRITERION
///////////////////////////////////////////////////////////////////////////////////////////

/// @brief Terminate when target function value reached
/// @details Useful when you know the global optimum or a satisfactory value
template <typename Real>
class TargetValueCriterion : public ITerminationCriterion<Real>
{
private:
    Real _targetValue;
    Real _tolerance;
    
public:
    TargetValueCriterion(Real targetValue, Real tolerance = Real(1e-6)) 
        : _targetValue(targetValue), _tolerance(tolerance)
    {
        if (tolerance < 0) {
            throw std::invalid_argument("TargetValueCriterion: tolerance must be non-negative");
        }
    }
    
    bool ShouldTerminate(const OptimizationState<Real>& state) const override {
        return state.fBest <= _targetValue + _tolerance;
    }
    
    std::string GetReason() const override {
        return "Target function value reached (target = " + 
               std::to_string(_targetValue) + " ± " + std::to_string(_tolerance) + ")";
    }
    
    void Reset() override {
        // Stateless
    }
    
    Real GetTargetValue() const { return _targetValue; }
    Real GetTolerance() const { return _tolerance; }
};

///////////////////////////////////////////////////////////////////////////////////////////
// SIMPLEX SIZE CRITERION
///////////////////////////////////////////////////////////////////////////////////////////

/// @brief Terminate when simplex size is below threshold (Nelder-Mead specific)
template <typename Real>
class SimplexSizeCriterion : public ITerminationCriterion<Real>
{
private:
    Real _threshold;
    
public:
    explicit SimplexSizeCriterion(Real threshold) 
        : _threshold(threshold)
    {
        if (threshold <= 0) {
            throw std::invalid_argument("SimplexSizeCriterion: threshold must be positive");
        }
    }
    
    bool ShouldTerminate(const OptimizationState<Real>& state) const override {
        return state.simplexSize < _threshold;
    }
    
    std::string GetReason() const override {
        return "Simplex converged (size = " + std::to_string(_threshold) + ")";
    }
    
    void Reset() override {
        // Stateless
    }
    
    Real GetThreshold() const { return _threshold; }
};

///////////////////////////////////////////////////////////////////////////////////////////
// COMPOSITE CRITERION
///////////////////////////////////////////////////////////////////////////////////////////

/// @brief Combine multiple criteria with AND or OR logic
/// @details Allows complex stopping conditions like:
///          "Stop if (maxIter OR stagnation) AND (ftol OR gtol)"
template <typename Real>
class CompositeCriterion : public ITerminationCriterion<Real>
{
public:
    enum class Mode { 
        ANY,  ///< OR: Stop if ANY criterion met
        ALL   ///< AND: Stop only if ALL criteria met
    };
    
private:
    std::vector<std::unique_ptr<ITerminationCriterion<Real>>> _criteria;
    Mode _mode;
    mutable std::string _reason;
    
public:
    explicit CompositeCriterion(Mode mode = Mode::ANY) 
        : _mode(mode) 
    {}
    
    /// @brief Add a criterion to the composite
    void AddCriterion(std::unique_ptr<ITerminationCriterion<Real>> criterion) {
        if (!criterion) {
            throw std::invalid_argument("CompositeCriterion: Cannot add null criterion");
        }
        _criteria.push_back(std::move(criterion));
    }
    
    /// @brief Convenience: Add criterion by constructing in-place
    template<typename CriterionType, typename... Args>
    void AddCriterion(Args&&... args) {
        _criteria.push_back(
            std::make_unique<CriterionType>(std::forward<Args>(args)...)
        );
    }
    
    bool ShouldTerminate(const OptimizationState<Real>& state) const override {
        if (_criteria.empty()) {
            return false;  // No criteria = never stop
        }
        
        if (_mode == Mode::ANY) {
            // OR logic: stop if ANY criterion met
            for (const auto& crit : _criteria) {
                if (crit->ShouldTerminate(state)) {
                    _reason = crit->GetReason();
                    return true;
                }
            }
            return false;
        }
        else {
            // AND logic: stop only if ALL criteria met
            std::vector<std::string> reasons;
            bool allMet = true;
            
            for (const auto& crit : _criteria) {
                if (crit->ShouldTerminate(state)) {
                    reasons.push_back(crit->GetReason());
                }
                else {
                    allMet = false;
                }
            }
            
            if (allMet && !reasons.empty()) {
                // Combine all reasons
                _reason = "Multiple criteria met: ";
                for (size_t i = 0; i < reasons.size(); ++i) {
                    if (i > 0) _reason += "; ";
                    _reason += reasons[i];
                }
                return true;
            }
            return false;
        }
    }
    
    std::string GetReason() const override {
        return _reason.empty() ? "Composite criterion not yet evaluated" : _reason;
    }
    
    void Reset() override {
        for (auto& crit : _criteria) {
            crit->Reset();
        }
        _reason.clear();
    }
    
    Mode GetMode() const { return _mode; }
    size_t GetCriterionCount() const { return _criteria.size(); }
};

///////////////////////////////////////////////////////////////////////////////////////////
// FACTORY HELPERS
///////////////////////////////////////////////////////////////////////////////////////////

/// @brief Create default termination criterion (common pattern)
/// @details Returns: maxIter OR stagnation (whichever comes first)
template <typename Real>
std::unique_ptr<ITerminationCriterion<Real>> CreateDefaultCriterion(
    int maxIter = 1000, 
    int maxStagnation = 100)
{
    auto composite = std::make_unique<CompositeCriterion<Real>>(
        CompositeCriterion<Real>::Mode::ANY
    );
    composite->AddCriterion(std::make_unique<MaxIterationsCriterion<Real>>(maxIter));
    composite->AddCriterion(std::make_unique<StagnationCriterion<Real>>(maxStagnation));
    return composite;
}

/// @brief Create tolerance-based criterion (for precise optimization)
/// @details Returns: (ftol AND gtol) OR maxIter (safety net)
template <typename Real>
std::unique_ptr<ITerminationCriterion<Real>> CreateToleranceCriterion(
    Real ftol = Real(1e-8),
    Real gtol = Real(1e-6),
    int maxIter = 10000)
{
    auto composite = std::make_unique<CompositeCriterion<Real>>(
        CompositeCriterion<Real>::Mode::ANY
    );
    
    // ftol AND gtol
    auto innerComposite = std::make_unique<CompositeCriterion<Real>>(
        CompositeCriterion<Real>::Mode::ALL
    );
    innerComposite->AddCriterion(std::make_unique<FunctionToleranceCriterion<Real>>(ftol));
    innerComposite->AddCriterion(std::make_unique<GradientNormCriterion<Real>>(gtol));
    
    composite->AddCriterion(std::move(innerComposite));
    composite->AddCriterion(std::make_unique<MaxIterationsCriterion<Real>>(maxIter));
    
    return composite;
}

} // namespace MML

#endif // MML_TERMINATION_CRITERIA_H

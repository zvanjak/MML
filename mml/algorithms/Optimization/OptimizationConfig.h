#if !defined MML_OPTIMIZATION_CONFIG_H
#define MML_OPTIMIZATION_CONFIG_H

#include "OptimizationCommon.h"
#include "TerminationCriteria.h"
#include "OptimizationObservers.h"

///////////////////////////////////////////////////////////////////////////////////////////
// OptimizationConfig.h - Fluent builder API for optimization configuration
//
// Provides OptimizationConfig class with fluent API for configuring:
//   - Termination criteria (pluggable, composable)
//   - Observers (console, trajectory, custom callbacks)
//   - Convenient defaults and presets
//
// Design goals:
//   - User-friendly: Readable, self-documenting API
//   - Flexible: Mix and match criteria and observers
//   - Safe: Sensible defaults, validated inputs
//   - Zero overhead: Only allocate what's used
//
// Usage example:
//   OptimizationConfig<double> config;
//   config.WithMaxIterations(5000)
//         .WithTolerance(1e-8)
//         .WithConsoleOutput(50)
//         .WithTrajectory(10);
//
//   SimulatedAnnealing<double> sa(config);
//   auto result = sa.Minimize(func, x0);
///////////////////////////////////////////////////////////////////////////////////////////

namespace MML
{

///////////////////////////////////////////////////////////////////////////////////////////
// OPTIMIZATION CONFIG - Fluent Builder
///////////////////////////////////////////////////////////////////////////////////////////

/// @brief Fluent builder for optimization configuration
/// @details Provides convenient API for configuring termination and observation
template <typename Real>
class OptimizationConfig
{
private:
    std::unique_ptr<ITerminationCriterion<Real>> _criterion;
    std::vector<std::shared_ptr<IOptimizationObserver<Real>>> _observers;
    
public:
    /// @brief Default constructor with sensible defaults
    /// @details Default: 1000 iterations OR 100 stagnation (whichever first)
    OptimizationConfig() {
        // Default criterion: maxIter OR stagnation
        _criterion = CreateDefaultCriterion<Real>(1000, 100);
    }
    
    /// @brief Copy constructor (deep copy of criterion, shared observers)
    OptimizationConfig(const OptimizationConfig&) = delete;  // Disable for now
    
    /// @brief Move constructor
    OptimizationConfig(OptimizationConfig&&) = default;
    
    /// @brief Assignment operators
    OptimizationConfig& operator=(const OptimizationConfig&) = delete;
    OptimizationConfig& operator=(OptimizationConfig&&) = default;
    
    //=================================================================================
    // TERMINATION CRITERIA CONFIGURATION
    //=================================================================================
    
    /// @brief Set maximum iterations (replaces current criterion)
    OptimizationConfig& WithMaxIterations(int maxIter) {
        _criterion = std::make_unique<MaxIterationsCriterion<Real>>(maxIter);
        return *this;
    }
    
    /// @brief Set function tolerance (replaces current criterion)
    OptimizationConfig& WithTolerance(Real ftol) {
        _criterion = std::make_unique<FunctionToleranceCriterion<Real>>(ftol);
        return *this;
    }
    
    /// @brief Set gradient norm tolerance (for gradient-based methods)
    OptimizationConfig& WithGradientTolerance(Real gtol) {
        _criterion = std::make_unique<GradientNormCriterion<Real>>(gtol);
        return *this;
    }
    
    /// @brief Set stagnation detection (replaces current criterion)
    OptimizationConfig& WithStagnation(int maxStagnation) {
        _criterion = std::make_unique<StagnationCriterion<Real>>(maxStagnation);
        return *this;
    }
    
    /// @brief Set time limit in seconds (replaces current criterion)
    OptimizationConfig& WithTimeLimit(double maxSeconds) {
        _criterion = std::make_unique<TimeLimitCriterion<Real>>(maxSeconds);
        return *this;
    }
    
    /// @brief Set time limit in minutes (convenience)
    OptimizationConfig& WithTimeLimitMinutes(double maxMinutes) {
        return WithTimeLimit(maxMinutes * 60.0);
    }
    
    /// @brief Set target value (replaces current criterion)
    OptimizationConfig& WithTargetValue(Real target, Real tolerance = Real(1e-6)) {
        _criterion = std::make_unique<TargetValueCriterion<Real>>(target, tolerance);
        return *this;
    }
    
    /// @brief Set custom termination criterion (replaces current)
    OptimizationConfig& WithCustomCriterion(
        std::unique_ptr<ITerminationCriterion<Real>> criterion) {
        if (!criterion) {
            throw std::invalid_argument("OptimizationConfig: Cannot set null criterion");
        }
        _criterion = std::move(criterion);
        return *this;
    }
    
    //=================================================================================
    // COMPOSITE CRITERIA (Advanced)
    //=================================================================================
    
    /// @brief Start building composite criterion with OR logic
    /// @details After calling this, use AndCriterion() or OrCriterion() to add more
    OptimizationConfig& WithAnyCriterion() {
        _criterion = std::make_unique<CompositeCriterion<Real>>(
            CompositeCriterion<Real>::Mode::ANY
        );
        return *this;
    }
    
    /// @brief Start building composite criterion with AND logic
    OptimizationConfig& WithAllCriteria() {
        _criterion = std::make_unique<CompositeCriterion<Real>>(
            CompositeCriterion<Real>::Mode::ALL
        );
        return *this;
    }
    
    /// @brief Add criterion to composite with OR logic (must call WithAnyCriterion first)
    template<typename CriterionType, typename... Args>
    OptimizationConfig& OrCriterion(Args&&... args) {
        auto* composite = dynamic_cast<CompositeCriterion<Real>*>(_criterion.get());
        if (!composite) {
            throw std::logic_error("OrCriterion requires composite criterion");
        }
        composite->template AddCriterion<CriterionType>(std::forward<Args>(args)...);
        return *this;
    }
    
    /// @brief Add criterion to composite with AND logic (must call WithAllCriteria first)
    template<typename CriterionType, typename... Args>
    OptimizationConfig& AndCriterion(Args&&... args) {
        auto* composite = dynamic_cast<CompositeCriterion<Real>*>(_criterion.get());
        if (!composite) {
            throw std::logic_error("AndCriterion requires composite criterion");
        }
        composite->template AddCriterion<CriterionType>(std::forward<Args>(args)...);
        return *this;
    }
    
    //=================================================================================
    // COMMON PATTERNS (Convenience)
    //=================================================================================
    
    /// @brief Standard pattern: maxIter OR stagnation
    OptimizationConfig& WithStandardTermination(int maxIter = 1000, 
                                                 int maxStagnation = 100) {
        _criterion = CreateDefaultCriterion<Real>(maxIter, maxStagnation);
        return *this;
    }
    
    /// @brief Precise pattern: (ftol AND gtol) OR maxIter
    OptimizationConfig& WithPreciseTermination(Real ftol = Real(1e-8),
                                                Real gtol = Real(1e-6),
                                                int maxIter = 10000) {
        _criterion = CreateToleranceCriterion<Real>(ftol, gtol, maxIter);
        return *this;
    }
    
    /// @brief Timed pattern: maxTime OR maxIter (safety net)
    OptimizationConfig& WithTimedTermination(double maxSeconds, 
                                             int maxIter = 100000) {
        auto composite = std::make_unique<CompositeCriterion<Real>>(
            CompositeCriterion<Real>::Mode::ANY
        );
        composite->template AddCriterion<TimeLimitCriterion<Real>>(maxSeconds);
        composite->template AddCriterion<MaxIterationsCriterion<Real>>(maxIter);
        _criterion = std::move(composite);
        return *this;
    }
    
    //=================================================================================
    // OBSERVER CONFIGURATION
    //=================================================================================
    
    /// @brief Add console output observer
    OptimizationConfig& WithConsoleOutput(int printEvery = 1, bool verbose = false) {
        _observers.push_back(
            std::make_shared<ConsoleObserver<Real>>(printEvery, verbose)
        );
        return *this;
    }
    
    /// @brief Add trajectory recording observer
    OptimizationConfig& WithTrajectory(int saveEvery = 1, size_t maxSize = 0) {
        _observers.push_back(
            std::make_shared<TrajectoryObserver<Real>>(saveEvery, maxSize)
        );
        return *this;
    }
    
    /// @brief Add custom observer
    OptimizationConfig& AddObserver(std::shared_ptr<IOptimizationObserver<Real>> observer) {
        if (!observer) {
            throw std::invalid_argument("OptimizationConfig: Cannot add null observer");
        }
        _observers.push_back(observer);
        return *this;
    }
    
    /// @brief Add observer by constructing in-place
    template<typename ObserverType, typename... Args>
    OptimizationConfig& AddObserver(Args&&... args) {
        _observers.push_back(
            std::make_shared<ObserverType>(std::forward<Args>(args)...)
        );
        return *this;
    }
    
    /// @brief Add callback observer with lambda
    OptimizationConfig& WithCallback(
        std::function<bool(const OptimizationState<Real>&)> onIteration) {
        auto callback = std::make_shared<CallbackObserver<Real>>();
        callback->SetOnIteration(onIteration);
        _observers.push_back(callback);
        return *this;
    }
    
    /// @brief Add full callback observer (start, iteration, complete)
    OptimizationConfig& WithCallbacks(
        std::function<void(const OptimizationState<Real>&)> onStart,
        std::function<bool(const OptimizationState<Real>&)> onIteration,
        std::function<void(const OptimizationState<Real>&, const std::string&)> onComplete) {
        auto callback = std::make_shared<CallbackObserver<Real>>();
        callback->SetOnStart(onStart);
        callback->SetOnIteration(onIteration);
        callback->SetOnComplete(onComplete);
        _observers.push_back(callback);
        return *this;
    }
    
    //=================================================================================
    // ACCESSORS (for optimizers)
    //=================================================================================
    
    /// @brief Get termination criterion (const)
    const ITerminationCriterion<Real>* GetCriterion() const {
        return _criterion.get();
    }
    
    /// @brief Get termination criterion (mutable, for reset)
    ITerminationCriterion<Real>* GetCriterion() {
        return _criterion.get();
    }
    
    /// @brief Get observers (const)
    const std::vector<std::shared_ptr<IOptimizationObserver<Real>>>& 
    GetObservers() const {
        return _observers;
    }
    
    /// @brief Get observers (mutable)
    std::vector<std::shared_ptr<IOptimizationObserver<Real>>>& 
    GetObservers() {
        return _observers;
    }
    
    /// @brief Check if any observers configured
    bool HasObservers() const {
        return !_observers.empty();
    }
    
    /// @brief Get specific observer by type (returns first match)
    template<typename ObserverType>
    std::shared_ptr<ObserverType> GetObserver() const {
        for (const auto& obs : _observers) {
            auto typed = std::dynamic_pointer_cast<ObserverType>(obs);
            if (typed) return typed;
        }
        return nullptr;
    }
    
    /// @brief Get trajectory observer (convenience for accessing trajectory)
    std::shared_ptr<TrajectoryObserver<Real>> GetTrajectoryObserver() const {
        return GetObserver<TrajectoryObserver<Real>>();
    }
    
    //=================================================================================
    // UTILITY
    //=================================================================================
    
    /// @brief Clear all observers
    OptimizationConfig& ClearObservers() {
        _observers.clear();
        return *this;
    }
    
    /// @brief Reset criterion and observers for new optimization run
    void Reset() {
        if (_criterion) {
            _criterion->Reset();
        }
        // Observers are reset by optimizer calling OnStart
    }
    
    /// @brief Create a copy with new instances (deep copy)
    OptimizationConfig Clone() const {
        // Note: This is tricky with unique_ptr criterion
        // For now, throw if cloning needed - users should build new config
        throw std::logic_error("OptimizationConfig cloning not yet implemented");
    }
};

///////////////////////////////////////////////////////////////////////////////////////////
// PRESET CONFIGS (Factory Functions)
///////////////////////////////////////////////////////////////////////////////////////////

/// @brief Create quick config for fast prototyping
template <typename Real>
OptimizationConfig<Real> QuickConfig(int maxIter = 100) {
    OptimizationConfig<Real> config;
    config.WithMaxIterations(maxIter)
          .WithConsoleOutput(10);
    return config;
}

/// @brief Create production config (precise, monitored)
template <typename Real>
OptimizationConfig<Real> ProductionConfig(Real ftol = Real(1e-8),
                                          Real gtol = Real(1e-6)) {
    OptimizationConfig<Real> config;
    config.WithPreciseTermination(ftol, gtol, 10000)
          .WithConsoleOutput(100, true)
          .WithTrajectory(50);
    return config;
}

/// @brief Create debug config (verbose output, full trajectory)
template <typename Real>
OptimizationConfig<Real> DebugConfig(int maxIter = 1000) {
    OptimizationConfig<Real> config;
    config.WithMaxIterations(maxIter)
          .WithConsoleOutput(1, true)   // Every iteration, verbose
          .WithTrajectory(1);            // Full trajectory
    return config;
}

/// @brief Create benchmark config (timed, minimal output)
template <typename Real>
OptimizationConfig<Real> BenchmarkConfig(double maxSeconds) {
    OptimizationConfig<Real> config;
    config.WithTimeLimit(maxSeconds);
    // No observers for minimal overhead
    return config;
}

} // namespace MML

#endif // MML_OPTIMIZATION_CONFIG_H

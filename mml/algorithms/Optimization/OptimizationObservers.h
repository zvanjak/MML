#if !defined MML_OPTIMIZATION_OBSERVERS_H
#define MML_OPTIMIZATION_OBSERVERS_H

#include "OptimizationCommon.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

///////////////////////////////////////////////////////////////////////////////////////////
// OptimizationObservers.h - Standard observers for optimization progress monitoring
//
// Provides concrete implementations of IOptimizationObserver:
//   - ConsoleObserver: Prints progress to console
//   - TrajectoryObserver: Saves OptimizationState snapshots with CSV export
//   - CallbackObserver: Lambda-friendly custom callbacks
//   - NullObserver: No-op observer (for disabling observation)
//
// Observers enable:
//   - Progress monitoring and logging
//   - Intermediate result saving (checkpointing)
//   - Custom early stopping logic
//   - GUI updates, plotting, analysis
//
// Usage patterns:
//   - Single observer: direct use
//   - Multiple observers: wrap in MultiObserver
//   - Custom behavior: use CallbackObserver with lambdas
///////////////////////////////////////////////////////////////////////////////////////////

namespace MML
{

///////////////////////////////////////////////////////////////////////////////////////////
// NULL OBSERVER (for disabling observation)
///////////////////////////////////////////////////////////////////////////////////////////

/// @brief No-op observer that does nothing
/// @details Useful as a default when no observation needed
template <typename Real>
class NullObserver : public IOptimizationObserver<Real>
{
public:
    void OnStart(const OptimizationState<Real>&) override {}
    bool OnIteration(const OptimizationState<Real>&) override { return true; }
    void OnComplete(const OptimizationState<Real>&, const std::string&) override {}
};

///////////////////////////////////////////////////////////////////////////////////////////
// CONSOLE OBSERVER
///////////////////////////////////////////////////////////////////////////////////////////

/// @brief Prints optimization progress to console
/// @details Configurable output frequency and verbosity
template <typename Real>
class ConsoleObserver : public IOptimizationObserver<Real>
{
private:
    int _printEvery;        ///< Print every N iterations (1 = every iteration)
    bool _verbose;          ///< Verbose output (more details)
    std::ostream& _out;     ///< Output stream (default: std::cout)
    
public:
    /// @brief Constructor
    /// @param printEvery Print every N iterations (default: 1)
    /// @param verbose Verbose output with more details (default: false)
    /// @param out Output stream (default: std::cout)
    explicit ConsoleObserver(int printEvery = 1, bool verbose = false, 
                            std::ostream& out = std::cout)
        : _printEvery(printEvery), _verbose(verbose), _out(out)
    {
        if (printEvery <= 0) {
            throw std::invalid_argument("ConsoleObserver: printEvery must be positive");
        }
    }
    
    void OnStart(const OptimizationState<Real>& state) override {
        _out << "\n=== Optimization Started ===" << std::endl;
        _out << "Initial f = " << state.fCurrent << std::endl;
        if (_verbose) {
            _out << "Dimension = " << state.xCurrent.size() << std::endl;
        }
        _out << std::endl;
        
        // Print header
        _out << std::setw(8) << "Iter" 
             << std::setw(15) << "fBest"
             << std::setw(12) << "fEvals";
        
        if (_verbose) {
            _out << std::setw(12) << "Accepted"
                 << std::setw(12) << "AcceptRate"
                 << std::setw(12) << "Time(s)";
        }
        _out << std::endl;
        _out << std::string(80, '-') << std::endl;
    }
    
    bool OnIteration(const OptimizationState<Real>& state) override {
        if (state.iteration % _printEvery == 0) {
            _out << std::setw(8) << state.iteration
                 << std::setw(15) << std::scientific << std::setprecision(6) 
                 << state.fBest
                 << std::setw(12) << state.funcEvals;
            
            if (_verbose) {
                _out << std::setw(12) << state.acceptedMoves
                     << std::setw(12) << std::fixed << std::setprecision(3) 
                     << state.GetAcceptanceRate()
                     << std::setw(12) << std::fixed << std::setprecision(2) 
                     << state.elapsedTime;
            }
            _out << std::endl;
        }
        return true;  // Continue optimization
    }
    
    void OnComplete(const OptimizationState<Real>& state, 
                   const std::string& reason) override {
        _out << std::string(80, '-') << std::endl;
        _out << "\n=== Optimization Complete ===" << std::endl;
        _out << "Reason: " << reason << std::endl;
        _out << "Final fBest = " << state.fBest << std::endl;
        _out << "Iterations: " << state.iteration << std::endl;
        _out << "Function evaluations: " << state.funcEvals << std::endl;
        
        if (_verbose) {
            _out << "Elapsed time: " << state.elapsedTime << " seconds" << std::endl;
            if (state.acceptedMoves > 0 || state.rejectedMoves > 0) {
                _out << "Acceptance rate: " << state.GetAcceptanceRate() << std::endl;
            }
            if (state.gradNorm > 0) {
                _out << "Final gradient norm: " << state.gradNorm << std::endl;
            }
        }
        _out << std::endl;
    }
    
    void SetPrintEvery(int printEvery) { 
        if (printEvery <= 0) {
            throw std::invalid_argument("printEvery must be positive");
        }
        _printEvery = printEvery; 
    }
    void SetVerbose(bool verbose) { _verbose = verbose; }
};

///////////////////////////////////////////////////////////////////////////////////////////
// TRAJECTORY OBSERVER
///////////////////////////////////////////////////////////////////////////////////////////

/// @brief Saves OptimizationState snapshots for post-analysis
/// @details Allows trajectory visualization, convergence analysis, etc.
template <typename Real>
class TrajectoryObserver : public IOptimizationObserver<Real>
{
private:
    std::vector<OptimizationState<Real>> _trajectory;
    int _saveEvery;     ///< Save every N iterations (1 = every iteration)
    size_t _maxSize;    ///< Maximum trajectory size (0 = unlimited)
    
public:
    /// @brief Constructor
    /// @param saveEvery Save state every N iterations (default: 1)
    /// @param maxSize Maximum trajectory size, 0 for unlimited (default: 0)
    explicit TrajectoryObserver(int saveEvery = 1, size_t maxSize = 0)
        : _saveEvery(saveEvery), _maxSize(maxSize)
    {
        if (saveEvery <= 0) {
            throw std::invalid_argument("TrajectoryObserver: saveEvery must be positive");
        }
    }
    
    void OnStart(const OptimizationState<Real>& state) override {
        _trajectory.clear();
        _trajectory.push_back(state);  // Save initial state
    }
    
    bool OnIteration(const OptimizationState<Real>& state) override {
        if (state.iteration % _saveEvery == 0) {
            // Check size limit
            if (_maxSize > 0 && _trajectory.size() >= _maxSize) {
                // Keep first and last, downsample middle
                DownsampleTrajectory();
            }
            _trajectory.push_back(state);
        }
        return true;  // Continue
    }
    
    void OnComplete(const OptimizationState<Real>& state, 
                   const std::string&) override {
        // Always save final state if not already saved
        if (_trajectory.empty() || _trajectory.back().iteration != state.iteration) {
            _trajectory.push_back(state);
        }
    }
    
    /// @brief Get the recorded trajectory
    const std::vector<OptimizationState<Real>>& GetTrajectory() const {
        return _trajectory;
    }
    
    /// @brief Get trajectory size
    size_t GetSize() const {
        return _trajectory.size();
    }
    
    /// @brief Clear trajectory (free memory)
    void Clear() {
        _trajectory.clear();
    }
    
    /// @brief Export trajectory to CSV file
    /// @param filename Output file path
    /// @param includeX Include solution vector components (default: false, can be large)
    void ExportCSV(const std::string& filename, bool includeX = false) const {
        std::ofstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("TrajectoryObserver: Cannot open file " + filename);
        }
        
        // Header
        file << "iteration,funcEvals,fCurrent,fBest,gradNorm,stepSize,temperature,"
             << "acceptedMoves,rejectedMoves,acceptRate,iterSinceImprovement,elapsedTime";
        
        if (includeX && !_trajectory.empty()) {
            size_t dim = _trajectory[0].xBest.size();
            for (size_t i = 0; i < dim; ++i) {
                file << ",x" << i;
            }
        }
        file << std::endl;
        
        // Data
        for (const auto& state : _trajectory) {
            file << state.iteration << ","
                 << state.funcEvals << ","
                 << state.fCurrent << ","
                 << state.fBest << ","
                 << state.gradNorm << ","
                 << state.stepSize << ","
                 << state.temperature << ","
                 << state.acceptedMoves << ","
                 << state.rejectedMoves << ","
                 << state.GetAcceptanceRate() << ","
                 << state.iterSinceImprovement << ","
                 << state.elapsedTime;
            
            if (includeX) {
                for (size_t i = 0; i < state.xBest.size(); ++i) {
                    file << "," << state.xBest[i];
                }
            }
            file << std::endl;
        }
        
        file.close();
    }
    
    /// @brief Export convergence plot data (iteration vs fBest)
    void ExportConvergencePlot(const std::string& filename) const {
        std::ofstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("TrajectoryObserver: Cannot open file " + filename);
        }
        
        file << "iteration,fBest" << std::endl;
        for (const auto& state : _trajectory) {
            file << state.iteration << "," << state.fBest << std::endl;
        }
        
        file.close();
    }
    
private:
    /// @brief Downsample trajectory when size limit reached
    void DownsampleTrajectory() {
        if (_trajectory.size() <= 2) return;
        
        // Keep every other element (simple decimation)
        std::vector<OptimizationState<Real>> downsampled;
        downsampled.reserve(_trajectory.size() / 2 + 1);
        
        for (size_t i = 0; i < _trajectory.size(); i += 2) {
            downsampled.push_back(_trajectory[i]);
        }
        
        _trajectory = std::move(downsampled);
        _saveEvery *= 2;  // Adjust save frequency
    }
};

///////////////////////////////////////////////////////////////////////////////////////////
// CALLBACK OBSERVER
///////////////////////////////////////////////////////////////////////////////////////////

/// @brief Lambda-friendly observer for custom callbacks
/// @details Allows inline custom logic without creating new observer classes
template <typename Real>
class CallbackObserver : public IOptimizationObserver<Real>
{
private:
    std::function<void(const OptimizationState<Real>&)> _onStartCallback;
    std::function<bool(const OptimizationState<Real>&)> _onIterCallback;
    std::function<void(const OptimizationState<Real>&, const std::string&)> _onCompleteCallback;
    
public:
    CallbackObserver() = default;
    
    /// @brief Set callback for optimization start
    void SetOnStart(std::function<void(const OptimizationState<Real>&)> callback) {
        _onStartCallback = callback;
    }
    
    /// @brief Set callback for each iteration
    /// @note Callback should return true to continue, false to abort
    void SetOnIteration(std::function<bool(const OptimizationState<Real>&)> callback) {
        _onIterCallback = callback;
    }
    
    /// @brief Set callback for optimization completion
    void SetOnComplete(
        std::function<void(const OptimizationState<Real>&, const std::string&)> callback) {
        _onCompleteCallback = callback;
    }
    
    void OnStart(const OptimizationState<Real>& state) override {
        if (_onStartCallback) {
            _onStartCallback(state);
        }
    }
    
    bool OnIteration(const OptimizationState<Real>& state) override {
        if (_onIterCallback) {
            return _onIterCallback(state);
        }
        return true;  // Continue by default
    }
    
    void OnComplete(const OptimizationState<Real>& state, 
                   const std::string& reason) override {
        if (_onCompleteCallback) {
            _onCompleteCallback(state, reason);
        }
    }
};

///////////////////////////////////////////////////////////////////////////////////////////
// MULTI OBSERVER (Composite)
///////////////////////////////////////////////////////////////////////////////////////////

/// @brief Combines multiple observers into one
/// @details Calls all observers in sequence
template <typename Real>
class MultiObserver : public IOptimizationObserver<Real>
{
private:
    std::vector<std::shared_ptr<IOptimizationObserver<Real>>> _observers;
    
public:
    MultiObserver() = default;
    
    /// @brief Add an observer to the collection
    void AddObserver(std::shared_ptr<IOptimizationObserver<Real>> observer) {
        if (!observer) {
            throw std::invalid_argument("MultiObserver: Cannot add null observer");
        }
        _observers.push_back(observer);
    }
    
    /// @brief Add observer by constructing in-place
    template<typename ObserverType, typename... Args>
    void AddObserver(Args&&... args) {
        _observers.push_back(
            std::make_shared<ObserverType>(std::forward<Args>(args)...)
        );
    }
    
    void OnStart(const OptimizationState<Real>& state) override {
        for (auto& obs : _observers) {
            obs->OnStart(state);
        }
    }
    
    bool OnIteration(const OptimizationState<Real>& state) override {
        // If any observer returns false, abort
        for (auto& obs : _observers) {
            if (!obs->OnIteration(state)) {
                return false;
            }
        }
        return true;
    }
    
    void OnComplete(const OptimizationState<Real>& state, 
                   const std::string& reason) override {
        for (auto& obs : _observers) {
            obs->OnComplete(state, reason);
        }
    }
    
    size_t GetObserverCount() const { return _observers.size(); }
    void Clear() { _observers.clear(); }
};

///////////////////////////////////////////////////////////////////////////////////////////
// EXAMPLE CALLBACK PATTERNS
///////////////////////////////////////////////////////////////////////////////////////////

/// @brief Example: Checkpointing callback (saves state every N iterations)
template <typename Real>
std::function<bool(const OptimizationState<Real>&)> 
MakeCheckpointCallback(int saveEvery, const std::string& baseFilename) {
    return [saveEvery, baseFilename](const OptimizationState<Real>& state) {
        if (state.iteration % saveEvery == 0) {
            std::string filename = baseFilename + "_iter" + 
                                 std::to_string(state.iteration) + ".txt";
            std::ofstream file(filename);
            if (file.is_open()) {
                file << "iteration: " << state.iteration << std::endl;
                file << "fBest: " << state.fBest << std::endl;
                file << "xBest:";
                for (size_t i = 0; i < state.xBest.size(); ++i) {
                    file << " " << state.xBest[i];
                }
                file << std::endl;
                file.close();
            }
        }
        return true;  // Continue
    };
}

/// @brief Example: Early stopping if function value too good (sanity check)
template <typename Real>
std::function<bool(const OptimizationState<Real>&)> 
MakeEarlyStoppingCallback(Real stopIfBelow) {
    return [stopIfBelow](const OptimizationState<Real>& state) {
        return state.fBest >= stopIfBelow;  // Stop if too good
    };
}

/// @brief Example: Alert if stagnation detected
template <typename Real>
std::function<bool(const OptimizationState<Real>&)> 
MakeStagnationAlertCallback(int alertAfter) {
    return [alertAfter](const OptimizationState<Real>& state) {
        if (state.iterSinceImprovement == alertAfter) {
            std::cout << "⚠ ALERT: No improvement for " << alertAfter 
                     << " iterations!" << std::endl;
        }
        return true;  // Continue
    };
}

} // namespace MML

#endif // MML_OPTIMIZATION_OBSERVERS_H

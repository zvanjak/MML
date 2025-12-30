///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        SimulatedAnnealing.h                                                ///
///  Description: Simulated Annealing optimization algorithm                          ///
///               Global optimization for multimodal and non-convex problems          ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_SIMULATED_ANNEALING_H
#define MML_SIMULATED_ANNEALING_H

#include "MMLBase.h"
#include "MMLExceptions.h"

#include "interfaces/IFunction.h"
#include "base/Vector.h"

// New optimization framework (optional)
#include "OptimizationCommon.h"
#include "OptimizationConfig.h"

#include <random>
#include <functional>
#include <cmath>
#include <limits>
#include <chrono>

namespace MML
{
    /////////////////////////////////////////////////////////////////////
    ///                   HEURISTIC OPTIMIZATION                      ///
    /////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////
    ///                   HeuristicOptimizationError                        ///
    ///////////////////////////////////////////////////////////////////////////
    class HeuristicOptimizationError : public std::runtime_error
    {
    public:
        explicit HeuristicOptimizationError(const std::string& message)
            : std::runtime_error("HeuristicOptimizationError: " + message) {}
    };

    ///////////////////////////////////////////////////////////////////////////
    ///                  HeuristicOptimizationResult                        ///
    ///////////////////////////////////////////////////////////////////////////
    /**
     * @brief Result structure for heuristic optimization methods
     */
    struct HeuristicOptimizationResult
    {
        Vector<Real> xbest;        ///< Best solution found
        Real         fbest;         ///< Function value at best solution
        int          iterations;    ///< Number of iterations performed
        int          funcEvals;     ///< Total function evaluations
        int          acceptedMoves; ///< Number of accepted moves (for SA)
        bool         converged;     ///< True if convergence criterion met

        HeuristicOptimizationResult()
            : fbest(std::numeric_limits<Real>::max()), iterations(0), 
              funcEvals(0), acceptedMoves(0), converged(false) {}

        HeuristicOptimizationResult(const Vector<Real>& x, Real f, int iter, 
                                    int fEvals, int accepted, bool conv)
            : xbest(x), fbest(f), iterations(iter), funcEvals(fEvals),
              acceptedMoves(accepted), converged(conv) {}
    };

    ///////////////////////////////////////////////////////////////////////////
    ///                     Cooling Schedule Interface                      ///
    ///////////////////////////////////////////////////////////////////////////
    /**
     * @brief Interface for temperature cooling schedules
     */
    class ICoolingSchedule
    {
    public:
        virtual ~ICoolingSchedule() = default;
        
        /**
         * @brief Get the temperature at a given iteration
         * @param iteration Current iteration number (0-based)
         * @param T0 Initial temperature
         * @return Temperature at this iteration
         */
        virtual Real Temperature(int iteration, Real T0) const = 0;
    };

    ///////////////////////////////////////////////////////////////////////////
    ///                     Exponential Cooling Schedule                    ///
    ///////////////////////////////////////////////////////////////////////////
    /**
     * @brief Exponential (geometric) cooling: T(k) = T0 * alpha^k
     * 
     * This is the most common cooling schedule. Temperature decreases
     * geometrically with each iteration.
     */
    class ExponentialCooling : public ICoolingSchedule
    {
    private:
        Real _alpha;  ///< Cooling rate (0 < alpha < 1, typically 0.9-0.99)

    public:
        /**
         * @brief Construct exponential cooling schedule
         * @param alpha Cooling rate (default 0.95)
         */
        explicit ExponentialCooling(Real alpha = 0.95) : _alpha(alpha)
        {
            if (alpha <= 0.0 || alpha >= 1.0)
                throw HeuristicOptimizationError("Cooling rate alpha must be in (0, 1)");
        }

        Real Temperature(int iteration, Real T0) const override
        {
            return T0 * std::pow(_alpha, iteration);
        }

        Real getAlpha() const { return _alpha; }
        void setAlpha(Real alpha)
        {
            if (alpha <= 0.0 || alpha >= 1.0)
                throw HeuristicOptimizationError("Cooling rate alpha must be in (0, 1)");
            _alpha = alpha;
        }
    };

    ///////////////////////////////////////////////////////////////////////////
    ///                     Linear Cooling Schedule                         ///
    ///////////////////////////////////////////////////////////////////////////
    /**
     * @brief Linear cooling: T(k) = T0 * (1 - k/maxIter)
     * 
     * Temperature decreases linearly to zero at maxIter.
     */
    class LinearCooling : public ICoolingSchedule
    {
    private:
        int _maxIter;  ///< Maximum iterations (temperature reaches 0)

    public:
        /**
         * @brief Construct linear cooling schedule
         * @param maxIter Maximum iterations
         */
        explicit LinearCooling(int maxIter) : _maxIter(maxIter)
        {
            if (maxIter <= 0)
                throw HeuristicOptimizationError("maxIter must be positive");
        }

        Real Temperature(int iteration, Real T0) const override
        {
            Real ratio = 1.0 - static_cast<Real>(iteration) / _maxIter;
            return T0 * std::max(REAL(0.0), ratio);
        }

        int getMaxIter() const { return _maxIter; }
        void setMaxIter(int maxIter)
        {
            if (maxIter <= 0)
                throw HeuristicOptimizationError("maxIter must be positive");
            _maxIter = maxIter;
        }
    };

    ///////////////////////////////////////////////////////////////////////////
    ///                   Logarithmic Cooling Schedule                      ///
    ///////////////////////////////////////////////////////////////////////////
    /**
     * @brief Logarithmic (slow) cooling: T(k) = T0 / (1 + c * ln(1 + k))
     * 
     * Very slow cooling that theoretically guarantees convergence to
     * the global minimum, but may be impractically slow.
     */
    class LogarithmicCooling : public ICoolingSchedule
    {
    private:
        Real _c;  ///< Cooling constant

    public:
        /**
         * @brief Construct logarithmic cooling schedule
         * @param c Cooling constant (default 1.0)
         */
        explicit LogarithmicCooling(Real c = 1.0) : _c(c)
        {
            if (c <= 0.0)
                throw HeuristicOptimizationError("Cooling constant c must be positive");
        }

        Real Temperature(int iteration, Real T0) const override
        {
            return T0 / (1.0 + _c * std::log(1.0 + iteration));
        }

        Real getC() const { return _c; }
        void setC(Real c)
        {
            if (c <= 0.0)
                throw HeuristicOptimizationError("Cooling constant c must be positive");
            _c = c;
        }
    };

    ///////////////////////////////////////////////////////////////////////////
    ///                      Adaptive Cooling Schedule                      ///
    ///////////////////////////////////////////////////////////////////////////
    /**
     * @brief Adaptive cooling that adjusts based on acceptance rate
     * 
     * Slows cooling when acceptance rate is low (stuck),
     * speeds up when acceptance rate is high (too hot).
     */
    class AdaptiveCooling : public ICoolingSchedule
    {
    private:
        Real _alphaFast;   ///< Fast cooling rate
        Real _alphaSlow;   ///< Slow cooling rate
        Real _targetAcceptRate; ///< Target acceptance rate

        mutable Real _currentAlpha;
        mutable Real _currentTemp;
        mutable int _lastIter;

    public:
        /**
         * @brief Construct adaptive cooling schedule
         * @param alphaFast Fast cooling rate (default 0.99)
         * @param alphaSlow Slow cooling rate (default 0.8)
         * @param targetAcceptRate Target acceptance rate (default 0.3)
         */
        AdaptiveCooling(Real alphaFast = 0.99, Real alphaSlow = 0.8, 
                        Real targetAcceptRate = 0.3)
            : _alphaFast(alphaFast), _alphaSlow(alphaSlow), 
              _targetAcceptRate(targetAcceptRate),
              _currentAlpha(alphaFast), _currentTemp(0), _lastIter(-1) {}

        Real Temperature(int iteration, Real T0) const override
        {
            if (iteration == 0 || _lastIter < 0)
            {
                _currentTemp = T0;
                _lastIter = 0;
            }
            
            // Simple exponential decay with current alpha
            while (_lastIter < iteration)
            {
                _currentTemp *= _currentAlpha;
                _lastIter++;
            }
            
            return _currentTemp;
        }

        /**
         * @brief Update cooling rate based on actual acceptance rate
         * @param acceptRate Current acceptance rate
         */
        void UpdateRate(Real acceptRate)
        {
            if (acceptRate > _targetAcceptRate)
                _currentAlpha = _alphaFast;  // Cool faster
            else
                _currentAlpha = _alphaSlow;  // Cool slower
        }

        void reset()
        {
            _currentAlpha = _alphaFast;
            _currentTemp = 0;
            _lastIter = -1;
        }
    };

    ///////////////////////////////////////////////////////////////////////////
    ///                     Neighbor Generator Interface                    ///
    ///////////////////////////////////////////////////////////////////////////
    /**
     * @brief Interface for generating neighbor solutions
     */
    class INeighborGenerator
    {
    public:
        virtual ~INeighborGenerator() = default;
        
        /**
         * @brief Generate a neighbor of the current solution
         * @param current Current solution
         * @param temperature Current temperature (may affect step size)
         * @return Neighbor solution
         */
        virtual Vector<Real> Generate(const Vector<Real>& current, Real temperature) = 0;
    };

    ///////////////////////////////////////////////////////////////////////////
    ///                    Uniform Random Neighbor Generator                ///
    ///////////////////////////////////////////////////////////////////////////
    /**
     * @brief Generate neighbors by uniform random perturbation
     * 
     * Each component is perturbed by a uniform random value in [-delta, delta]
     * where delta can optionally scale with temperature.
     */
    class UniformNeighborGenerator : public INeighborGenerator
    {
    private:
        Real _delta;          ///< Maximum perturbation per component
        bool _scaleWithTemp;  ///< Whether to scale delta with temperature
        Real _T0;             ///< Reference temperature for scaling
        
        std::mt19937 _rng;
        std::uniform_real_distribution<Real> _dist;

        Vector<Real> _lowerBounds;  ///< Lower bounds (optional)
        Vector<Real> _upperBounds;  ///< Upper bounds (optional)
        bool _hasBounds;

    public:
        /**
         * @brief Construct uniform neighbor generator
         * @param delta Maximum perturbation magnitude
         * @param scaleWithTemp Scale delta proportionally to temperature
         * @param T0 Reference temperature (initial temperature)
         * @param seed Random seed (0 for random)
         */
        UniformNeighborGenerator(Real delta = 1.0, bool scaleWithTemp = false, 
                                  Real T0 = 1.0, unsigned int seed = 0)
            : _delta(delta), _scaleWithTemp(scaleWithTemp), _T0(T0),
              _dist(-1.0, 1.0), _hasBounds(false)
        {
            if (seed == 0)
                _rng.seed(std::random_device{}());
            else
                _rng.seed(seed);
        }

        /**
         * @brief Set bounds for the search space
         */
        void SetBounds(const Vector<Real>& lower, const Vector<Real>& upper)
        {
            _lowerBounds = lower;
            _upperBounds = upper;
            _hasBounds = true;
        }

        void ClearBounds() { _hasBounds = false; }

        Vector<Real> Generate(const Vector<Real>& current, Real temperature) override
        {
            Vector<Real> neighbor(current.size());
            Real scale = _scaleWithTemp ? _delta * (temperature / _T0) : _delta;
            
            for (int i = 0; i < current.size(); ++i)
            {
                neighbor[i] = current[i] + scale * _dist(_rng);
                
                // Apply bounds if set
                if (_hasBounds)
                {
                    if (neighbor[i] < _lowerBounds[i])
                        neighbor[i] = _lowerBounds[i];
                    if (neighbor[i] > _upperBounds[i])
                        neighbor[i] = _upperBounds[i];
                }
            }
            
            return neighbor;
        }

        Real getDelta() const { return _delta; }
        void setDelta(Real delta) { _delta = delta; }
        void setT0(Real T0) { _T0 = T0; }
    };

    ///////////////////////////////////////////////////////////////////////////
    ///                  Gaussian Neighbor Generator                        ///
    ///////////////////////////////////////////////////////////////////////////
    /**
     * @brief Generate neighbors by Gaussian perturbation
     * 
     * Each component is perturbed by a Gaussian random value with 
     * standard deviation sigma.
     */
    class GaussianNeighborGenerator : public INeighborGenerator
    {
    private:
        Real _sigma;          ///< Standard deviation of perturbation
        bool _scaleWithTemp;  ///< Whether to scale sigma with temperature
        Real _T0;             ///< Reference temperature for scaling
        
        std::mt19937 _rng;
        std::normal_distribution<Real> _dist;

        Vector<Real> _lowerBounds;
        Vector<Real> _upperBounds;
        bool _hasBounds;

    public:
        /**
         * @brief Construct Gaussian neighbor generator
         * @param sigma Standard deviation of perturbation
         * @param scaleWithTemp Scale sigma proportionally to sqrt(temperature)
         * @param T0 Reference temperature
         * @param seed Random seed (0 for random)
         */
        GaussianNeighborGenerator(Real sigma = 1.0, bool scaleWithTemp = false,
                                   Real T0 = 1.0, unsigned int seed = 0)
            : _sigma(sigma), _scaleWithTemp(scaleWithTemp), _T0(T0),
              _dist(0.0, 1.0), _hasBounds(false)
        {
            if (seed == 0)
                _rng.seed(std::random_device{}());
            else
                _rng.seed(seed);
        }

        void SetBounds(const Vector<Real>& lower, const Vector<Real>& upper)
        {
            _lowerBounds = lower;
            _upperBounds = upper;
            _hasBounds = true;
        }

        void ClearBounds() { _hasBounds = false; }

        Vector<Real> Generate(const Vector<Real>& current, Real temperature) override
        {
            Vector<Real> neighbor(current.size());
            Real scale = _scaleWithTemp ? _sigma * std::sqrt(temperature / _T0) : _sigma;
            
            for (int i = 0; i < current.size(); ++i)
            {
                neighbor[i] = current[i] + scale * _dist(_rng);
                
                if (_hasBounds)
                {
                    if (neighbor[i] < _lowerBounds[i])
                        neighbor[i] = _lowerBounds[i];
                    if (neighbor[i] > _upperBounds[i])
                        neighbor[i] = _upperBounds[i];
                }
            }
            
            return neighbor;
        }

        Real getSigma() const { return _sigma; }
        void setSigma(Real sigma) { _sigma = sigma; }
        void setT0(Real T0) { _T0 = T0; }
    };

    ///////////////////////////////////////////////////////////////////////////
    ///                      Simulated Annealing                            ///
    ///////////////////////////////////////////////////////////////////////////
    /**
     * @brief Simulated Annealing optimization algorithm
     * 
     * Simulated Annealing is a probabilistic metaheuristic for global optimization.
     * It mimics the physical process of heating and slowly cooling a material
     * to decrease defects (minimizing energy).
     * 
     * Algorithm:
     * 1. Start with initial solution x and temperature T
     * 2. Generate neighbor x' of current solution
     * 3. If f(x') < f(x), accept x' (downhill move)
     * 4. If f(x') >= f(x), accept with probability exp(-(f(x')-f(x))/T)
     * 5. Reduce temperature according to cooling schedule
     * 6. Repeat until stopping criterion met
     * 
     * The key insight is that at high temperatures, uphill moves are likely
     * accepted (allowing escape from local minima), while at low temperatures
     * the algorithm behaves like hill-climbing.
     * 
     * Reference: Kirkpatrick, Gelatt, Vecchi (1983)
     */
    class SimulatedAnnealing
    {
    public:
        /// Stopping criteria
        enum class StopCriteria
        {
            MaxIterations,   ///< Stop after maxIter iterations
            MinTemperature,  ///< Stop when temperature drops below threshold
            NoImprovement,   ///< Stop after N iterations without improvement
            Combined         ///< Use all criteria
        };

    private:
        Real _T0;                ///< Initial temperature
        Real _Tmin;              ///< Minimum temperature (stopping criterion)
        int  _maxIter;           ///< Maximum iterations
        int  _stagnationLimit;   ///< Iterations without improvement to stop
        StopCriteria _stopCrit;  ///< Stopping criterion to use

        // Components
        std::unique_ptr<ICoolingSchedule> _coolingSchedule;
        std::unique_ptr<INeighborGenerator> _neighborGen;

        // Random number generation for acceptance
        mutable std::mt19937 _rng;
        mutable std::uniform_real_distribution<Real> _acceptDist;

    public:
        /**
         * @brief Construct Simulated Annealing optimizer with default components
         * @param T0 Initial temperature (default 100.0)
         * @param Tmin Minimum temperature (default 1e-8)
         * @param maxIter Maximum iterations (default 10000)
         * @param stagnationLimit Iterations without improvement (default 1000)
         * @param stopCrit Stopping criterion (default Combined)
         * @param seed Random seed (0 for random)
         */
        SimulatedAnnealing(Real T0 = 100.0, Real Tmin = 1e-8, int maxIter = 10000,
                           int stagnationLimit = 1000, 
                           StopCriteria stopCrit = StopCriteria::Combined,
                           unsigned int seed = 0)
            : _T0(T0), _Tmin(Tmin), _maxIter(maxIter), 
              _stagnationLimit(stagnationLimit), _stopCrit(stopCrit),
              _acceptDist(0.0, 1.0)
        {
            // Default cooling schedule: exponential with alpha=0.995 (slow cooling)
            _coolingSchedule = std::make_unique<ExponentialCooling>(0.995);
            
            // Default neighbor generator: Gaussian with moderate sigma
            // Pass through the provided seed so neighbor generation is reproducible
            auto gen = std::make_unique<GaussianNeighborGenerator>(0.5, true, T0, seed);
            _neighborGen = std::move(gen);

            if (seed == 0)
                _rng.seed(std::random_device{}());
            else
                _rng.seed(seed);
        }

        /**
         * @brief Construct with custom cooling schedule and neighbor generator
         */
        SimulatedAnnealing(std::unique_ptr<ICoolingSchedule> cooling,
                           std::unique_ptr<INeighborGenerator> neighbor,
                           Real T0 = 100.0, Real Tmin = 1e-8, int maxIter = 10000,
                           int stagnationLimit = 1000,
                           StopCriteria stopCrit = StopCriteria::Combined,
                           unsigned int seed = 0)
            : _T0(T0), _Tmin(Tmin), _maxIter(maxIter),
              _stagnationLimit(stagnationLimit), _stopCrit(stopCrit),
              _coolingSchedule(std::move(cooling)),
              _neighborGen(std::move(neighbor)),
              _acceptDist(0.0, 1.0)
        {
            if (seed == 0)
                _rng.seed(std::random_device{}());
            else
                _rng.seed(seed);
        }

        /**
         * @brief Minimize a function using simulated annealing
         * @param func Function to minimize (callable: Vector<Real> -> Real)
         * @param x0 Initial solution
         * @return Optimization result
         */
        template<typename Func>
        HeuristicOptimizationResult Minimize(Func& func, const Vector<Real>& x0) const
        {
            int n = x0.size();
            
            // Initialize
            Vector<Real> x = x0;
            Real fx = func(x);
            
            Vector<Real> xbest = x;
            Real fbest = fx;
            
            int funcEvals = 1;
            int acceptedMoves = 0;
            int iterSinceImprovement = 0;
            
            // Main loop
            for (int iter = 0; iter < _maxIter; ++iter)
            {
                // Get current temperature
                Real T = _coolingSchedule->Temperature(iter, _T0);
                
                // Check temperature stopping criterion
                if ((_stopCrit == StopCriteria::MinTemperature || 
                     _stopCrit == StopCriteria::Combined) && T < _Tmin)
                {
                    return HeuristicOptimizationResult(xbest, fbest, iter, 
                                                        funcEvals, acceptedMoves, true);
                }
                
                // Generate neighbor
                Vector<Real> xnew = _neighborGen->Generate(x, T);
                Real fxnew = func(xnew);
                funcEvals++;
                
                // Compute acceptance probability
                Real deltaE = fxnew - fx;
                bool accept = false;
                
                if (deltaE < 0)
                {
                    // Always accept improvement
                    accept = true;
                }
                else if (T > 0)
                {
                    // Accept with Boltzmann probability
                    Real prob = std::exp(-deltaE / T);
                    accept = (_acceptDist(_rng) < prob);
                }
                
                if (accept)
                {
                    x = xnew;
                    fx = fxnew;
                    acceptedMoves++;
                    
                    // Update best solution
                    if (fx < fbest)
                    {
                        xbest = x;
                        fbest = fx;
                        iterSinceImprovement = 0;
                    }
                    else
                    {
                        iterSinceImprovement++;
                    }
                }
                else
                {
                    iterSinceImprovement++;
                }
                
                // Check stagnation stopping criterion
                if ((_stopCrit == StopCriteria::NoImprovement || 
                     _stopCrit == StopCriteria::Combined) && 
                    iterSinceImprovement >= _stagnationLimit)
                {
                    return HeuristicOptimizationResult(xbest, fbest, iter + 1,
                                                        funcEvals, acceptedMoves, true);
                }
            }
            
            // Reached max iterations
            bool converged = (_stopCrit == StopCriteria::MaxIterations);
            return HeuristicOptimizationResult(xbest, fbest, _maxIter, 
                                                funcEvals, acceptedMoves, converged);
        }

        /**
         * @brief Minimize using NEW optimization framework (with config)
         * @param func Function to minimize
         * @param x0 Initial solution
         * @param config Optimization configuration (criteria, observers, etc.)
         * @return Enhanced OptimizationResult with trajectory if configured
         */
        template<typename Func>
        OptimizationResult<Real> Minimize(Func& func, const Vector<Real>& x0,
                                           OptimizationConfig<Real>& config) const
        {
            using namespace std::chrono;
            auto startTime = steady_clock::now();
            
            int n = x0.size();
            
            // Initialize optimization state
            OptimizationState<Real> state;
            state.iteration = 0;
            state.funcEvals = 0;
            state.xCurrent = x0;
            state.fCurrent = func(x0);
            state.xBest = x0;
            state.fBest = state.fCurrent;
            state.fBestPrevious = state.fBest;
            state.iterSinceImprovement = 0;
            state.temperature = _T0;
            state.acceptedMoves = 0;
            state.rejectedMoves = 0;
            state.elapsedTime = 0.0;
            state.funcEvals = 1;
            
            // Get config components
            auto* criterion = config.GetCriterion();
            auto& observers = config.GetObservers();
            
            if (!criterion) {
                throw HeuristicOptimizationError("No termination criterion provided");
            }
            
            // Reset criterion for new run
            criterion->Reset();
            
            // Notify observers: optimization starting
            for (auto& obs : observers) {
                obs->OnStart(state);
            }
            
            // Main optimization loop
            while (true)
            {
                // Get current temperature
                Real T = _coolingSchedule->Temperature(state.iteration, _T0);
                state.temperature = T;
                
                // Check termination criterion
                if (criterion->ShouldTerminate(state)) {
                    // Notify observers: optimization complete
                    std::string reason = criterion->GetReason();
                    for (auto& obs : observers) {
                        obs->OnComplete(state, reason);
                    }
                    
                    // Create result
                    OptimizationResult<Real> result = CreateResult(state, reason, true);
                    
                    // Copy trajectory if TrajectoryObserver used
                    auto trajObs = config.GetTrajectoryObserver();
                    if (trajObs) {
                        result.trajectory = trajObs->GetTrajectory();
                    }
                    
                    return result;
                }
                
                // Generate neighbor
                Vector<Real> xnew = _neighborGen->Generate(state.xCurrent, T);
                Real fxnew = func(xnew);
                state.funcEvals++;
                
                // Compute acceptance probability
                Real deltaE = fxnew - state.fCurrent;
                bool accept = false;
                
                if (deltaE < 0)
                {
                    // Always accept improvement
                    accept = true;
                }
                else if (T > 0)
                {
                    // Accept with Boltzmann probability
                    Real prob = std::exp(-deltaE / T);
                    accept = (_acceptDist(_rng) < prob);
                }
                
                if (accept)
                {
                    state.xCurrent = xnew;
                    state.fCurrent = fxnew;
                    state.acceptedMoves++;
                    
                    // Update best solution
                    if (state.fCurrent < state.fBest)
                    {
                        state.xBest = state.xCurrent;
                        state.fBestPrevious = state.fBest;
                        state.fBest = state.fCurrent;
                        state.iterSinceImprovement = 0;
                    }
                    else
                    {
                        state.iterSinceImprovement++;
                    }
                }
                else
                {
                    state.rejectedMoves++;
                    state.iterSinceImprovement++;
                }
                
                // Update timing
                auto now = steady_clock::now();
                state.elapsedTime = duration_cast<duration<double>>(now - startTime).count();
                
                // Increment iteration
                state.iteration++;
                
                // Notify observers: iteration complete (can abort if returns false)
                bool continueOpt = true;
                for (auto& obs : observers) {
                    if (!obs->OnIteration(state)) {
                        continueOpt = false;
                        break;
                    }
                }
                
                if (!continueOpt) {
                    // Observer requested abort
                    std::string reason = "Aborted by observer";
                    for (auto& obs : observers) {
                        obs->OnComplete(state, reason);
                    }
                    
                    OptimizationResult<Real> result = CreateResult(state, reason, false);
                    
                    auto trajObs = config.GetTrajectoryObserver();
                    if (trajObs) {
                        result.trajectory = trajObs->GetTrajectory();
                    }
                    
                    return result;
                }
            }
        }

        /**
         * @brief Minimize a function with bounds
         */
        template<typename Func>
        HeuristicOptimizationResult Minimize(Func& func, const Vector<Real>& x0,
                                              const Vector<Real>& lowerBounds,
                                              const Vector<Real>& upperBounds) const
        {
            // Set bounds on neighbor generator
            if (auto* uniform = dynamic_cast<UniformNeighborGenerator*>(_neighborGen.get()))
            {
                uniform->SetBounds(lowerBounds, upperBounds);
            }
            else if (auto* gaussian = dynamic_cast<GaussianNeighborGenerator*>(_neighborGen.get()))
            {
                gaussian->SetBounds(lowerBounds, upperBounds);
            }
            
            return Minimize(func, x0);
        }

        // Accessors
        Real getT0() const { return _T0; }
        void setT0(Real T0) { _T0 = T0; }
        
        Real getTmin() const { return _Tmin; }
        void setTmin(Real Tmin) { _Tmin = Tmin; }
        
        int getMaxIter() const { return _maxIter; }
        void setMaxIter(int maxIter) { _maxIter = maxIter; }
        
        int getStagnationLimit() const { return _stagnationLimit; }
        void setStagnationLimit(int limit) { _stagnationLimit = limit; }
        
        StopCriteria getStopCriteria() const { return _stopCrit; }
        void setStopCriteria(StopCriteria crit) { _stopCrit = crit; }

        /**
         * @brief Set custom cooling schedule
         */
        void setCoolingSchedule(std::unique_ptr<ICoolingSchedule> schedule)
        {
            _coolingSchedule = std::move(schedule);
        }

        /**
         * @brief Set custom neighbor generator
         */
        void setNeighborGenerator(std::unique_ptr<INeighborGenerator> gen)
        {
            _neighborGen = std::move(gen);
        }
    };

    ///////////////////////////////////////////////////////////////////////////
    ///                     Convenience Functions                           ///
    ///////////////////////////////////////////////////////////////////////////

    /**
     * @brief Simple simulated annealing minimization
     * @param func Function to minimize
     * @param x0 Initial solution
     * @param T0 Initial temperature (default 100)
     * @param maxIter Maximum iterations (default 10000)
     * @return Optimization result
     */
    template<typename Func>
    HeuristicOptimizationResult SimulatedAnnealingMinimize(
        Func& func, 
        const Vector<Real>& x0,
        Real T0 = 100.0,
        int maxIter = 10000)
    {
        SimulatedAnnealing sa(T0, 1e-8, maxIter);
        return sa.Minimize(func, x0);
    }

    /**
     * @brief Simulated annealing with bounds
     */
    template<typename Func>
    HeuristicOptimizationResult SimulatedAnnealingMinimize(
        Func& func,
        const Vector<Real>& x0,
        const Vector<Real>& lowerBounds,
        const Vector<Real>& upperBounds,
        Real T0 = 100.0,
        int maxIter = 10000)
    {
        SimulatedAnnealing sa(T0, 1e-8, maxIter);
        return sa.Minimize(func, x0, lowerBounds, upperBounds);
    }

    /////////////////////////////////////////////////////////////////////
    ///                     FUTURE: Genetic Algorithms                ///
    /////////////////////////////////////////////////////////////////////
    
    // TODO: Implement genetic algorithms
    // - Binary GA
    // - Real-coded GA
    // - Selection strategies (roulette, tournament, rank)
    // - Crossover operators (single-point, multi-point, uniform)
    // - Mutation operators
    // - Elitism

} // namespace MML

#endif // MML_OPTIMIZATION_HEURISTIC_H
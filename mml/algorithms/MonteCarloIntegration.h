///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        MonteCarloIntegration.h                                             ///
///  Description: Monte Carlo integration for high-dimensional integrals              ///
///               Plain, stratified sampling, and importance sampling methods         ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_MONTE_CARLO_INTEGRATION_H
#define MML_MONTE_CARLO_INTEGRATION_H

#include "MMLBase.h"
#include "interfaces/IFunction.h"
#include "base/VectorN.h"

#include <random>
#include <numeric>
#include <cmath>

namespace MML
{
    //////////////////////////////////////////////////////////////////////////
    /// @brief Result of a Monte Carlo integration
    /// @note For production code, check error_estimate, variance, and converged!
    //////////////////////////////////////////////////////////////////////////
    struct MonteCarloResult
    {
        Real value;              ///< Estimated integral value
        Real error_estimate;     ///< Statistical error estimate (standard error)
        Real variance;           ///< Sample variance
        size_t samples_used;     ///< Number of samples actually used
        bool converged;          ///< Whether error is below requested tolerance
        
        MonteCarloResult(Real val = 0.0, Real err = 0.0, Real var = 0.0, 
                         size_t n = 0, bool conv = false)
            : value(val), error_estimate(err), variance(var), 
              samples_used(n), converged(conv) {}
        
        /// Implicit conversion to Real for backward compatibility
        /// @warning Silently discards error_estimate, variance, samples_used, and converged!
        /// @deprecated Prefer explicit .value access in new code
        operator Real() const { return value; }
    };
    
    //////////////////////////////////////////////////////////////////////////
    /// @brief Configuration for Monte Carlo integration
    //////////////////////////////////////////////////////////////////////////
    struct MonteCarloConfig
    {
        size_t num_samples = 100000;     ///< Number of random samples
        Real target_error = 1e-3;         ///< Target relative error (for adaptive)
        unsigned int seed = 0;            ///< Random seed (0 = use random_device)
        bool use_antithetic = false;      ///< Use antithetic variates for variance reduction
        
        MonteCarloConfig& samples(size_t n) { num_samples = n; return *this; }
        MonteCarloConfig& error(Real e) { target_error = e; return *this; }
        MonteCarloConfig& randomSeed(unsigned int s) { seed = s; return *this; }
        MonteCarloConfig& antithetic(bool a) { use_antithetic = a; return *this; }
    };

    //////////////////////////////////////////////////////////////////////////
    /// @brief Plain Monte Carlo integrator for N-dimensional integrals
    /// 
    /// Estimates ∫∫...∫ f(x) dx over a hyperrectangular domain [a,b]^N
    /// using simple random sampling.
    /// 
    /// Error scales as O(1/√N) regardless of dimension - this is the key
    /// advantage over deterministic methods for high dimensions.
    //////////////////////////////////////////////////////////////////////////
    template<int N>
    class MonteCarloIntegrator
    {
    private:
        std::mt19937_64 _rng;
        std::uniform_real_distribution<Real> _uniform{0.0, 1.0};
        
    public:
        /// Initialize with optional seed (0 = random)
        explicit MonteCarloIntegrator(unsigned int seed = 0)
        {
            if (seed == 0) {
                std::random_device rd;
                _rng.seed(rd());
            } else {
                _rng.seed(seed);
            }
        }
        
        /// Set the random seed
        void seed(unsigned int s) { _rng.seed(s); }
        
        /// @brief Integrate f over hyperrectangle [lower, upper]
        /// @param func The N-dimensional scalar function to integrate
        /// @param lower Lower bounds for each dimension
        /// @param upper Upper bounds for each dimension
        /// @param config Integration configuration
        /// @return MonteCarloResult with value, error estimate, and diagnostics
        MonteCarloResult integrate(
            const IScalarFunction<N>& func,
            const VectorN<Real, N>& lower,
            const VectorN<Real, N>& upper,
            const MonteCarloConfig& config = MonteCarloConfig())
        {
            // Compute volume of integration domain
            Real volume = 1.0;
            for (int i = 0; i < N; ++i) {
                volume *= (upper[i] - lower[i]);
            }
            
            if (volume == 0.0) {
                return MonteCarloResult(0.0, 0.0, 0.0, 0, true);
            }
            
            // Reseed if specified
            if (config.seed != 0) {
                _rng.seed(config.seed);
            }
            
            // Accumulate sum and sum of squares for variance estimation
            Real sum = 0.0;
            Real sum_sq = 0.0;
            size_t n = config.num_samples;
            
            VectorN<Real, N> point;
            
            if (config.use_antithetic) {
                // Antithetic variates: use both x and (1-x) to reduce variance
                n = (n / 2) * 2;  // Ensure even number
                for (size_t i = 0; i < n; i += 2) {
                    // Generate random point
                    for (int d = 0; d < N; ++d) {
                        Real u = _uniform(_rng);
                        point[d] = lower[d] + u * (upper[d] - lower[d]);
                    }
                    Real f1 = func(point);
                    
                    // Antithetic point
                    for (int d = 0; d < N; ++d) {
                        point[d] = lower[d] + upper[d] - point[d];  // Reflect
                    }
                    Real f2 = func(point);
                    
                    // Average the pair
                    Real f_avg = 0.5 * (f1 + f2);
                    sum += f_avg;
                    sum_sq += f_avg * f_avg;
                }
                n /= 2;  // We have n/2 averaged samples
            } else {
                // Standard Monte Carlo
                for (size_t i = 0; i < n; ++i) {
                    // Generate uniform random point in [lower, upper]
                    for (int d = 0; d < N; ++d) {
                        Real u = _uniform(_rng);
                        point[d] = lower[d] + u * (upper[d] - lower[d]);
                    }
                    
                    Real f_val = func(point);
                    sum += f_val;
                    sum_sq += f_val * f_val;
                }
            }
            
            // Compute statistics
            Real mean = sum / n;
            Real mean_sq = sum_sq / n;
            Real variance = mean_sq - mean * mean;
            
            // Standard error of the mean
            Real std_error = std::sqrt(variance / n);
            
            // Integral estimate = volume * mean
            Real integral = volume * mean;
            Real error = volume * std_error;
            
            // Check convergence against target
            bool converged = (std::abs(integral) < PrecisionValues<Real>::DivisionSafetyThreshold) || 
                             (error / std::abs(integral) < config.target_error);
            
            return MonteCarloResult(integral, error, variance, 
                                    config.use_antithetic ? n * 2 : n, converged);
        }
        
        /// @brief Convenience overload for unit hypercube [0,1]^N
        MonteCarloResult integrate(
            const IScalarFunction<N>& func,
            const MonteCarloConfig& config = MonteCarloConfig())
        {
            VectorN<Real, N> lower, upper;
            for (int i = 0; i < N; ++i) {
                lower[i] = 0.0;
                upper[i] = 1.0;
            }
            return integrate(func, lower, upper, config);
        }
    };
    
    //////////////////////////////////////////////////////////////////////////
    /// @brief Convenience function for 1D Monte Carlo integration
    //////////////////////////////////////////////////////////////////////////
    inline MonteCarloResult IntegrateMonteCarlo1D(
        const IRealFunction& func,
        Real a, Real b,
        const MonteCarloConfig& config = MonteCarloConfig())
    {
        // Wrap IRealFunction in IScalarFunction<1>
        class Wrapper : public IScalarFunction<1> {
            const IRealFunction& _f;
        public:
            Wrapper(const IRealFunction& f) : _f(f) {}
            Real operator()(const VectorN<Real, 1>& x) const override {
                return _f(x[0]);
            }
        };
        
        Wrapper wrapper(func);
        MonteCarloIntegrator<1> integrator(config.seed);
        
        VectorN<Real, 1> lower, upper;
        lower[0] = a;
        upper[0] = b;
        
        return integrator.integrate(wrapper, lower, upper, config);
    }
    
    //////////////////////////////////////////////////////////////////////////
    /// @brief Stratified sampling Monte Carlo for variance reduction
    /// 
    /// Divides the domain into strata and samples proportionally from each.
    /// Reduces variance when the function varies differently across regions.
    //////////////////////////////////////////////////////////////////////////
    template<int N>
    class StratifiedMonteCarloIntegrator
    {
    private:
        std::mt19937_64 _rng;
        std::uniform_real_distribution<Real> _uniform{0.0, 1.0};
        
    public:
        explicit StratifiedMonteCarloIntegrator(unsigned int seed = 0)
        {
            if (seed == 0) {
                std::random_device rd;
                _rng.seed(rd());
            } else {
                _rng.seed(seed);
            }
        }
        
        void seed(unsigned int s) { _rng.seed(s); }
        
        /// @brief Integrate using stratified sampling
        /// @param func Function to integrate
        /// @param lower Lower bounds
        /// @param upper Upper bounds
        /// @param strata_per_dim Number of strata per dimension (total = strata^N)
        /// @param samples_per_stratum Samples per stratum
        MonteCarloResult integrate(
            const IScalarFunction<N>& func,
            const VectorN<Real, N>& lower,
            const VectorN<Real, N>& upper,
            int strata_per_dim = 10,
            int samples_per_stratum = 10,
            unsigned int seed = 0)
        {
            if (seed != 0) {
                _rng.seed(seed);
            }
            
            // Compute total volume
            Real volume = 1.0;
            VectorN<Real, N> delta;
            for (int d = 0; d < N; ++d) {
                delta[d] = (upper[d] - lower[d]) / strata_per_dim;
                volume *= (upper[d] - lower[d]);
            }
            
            if (volume == 0.0) {
                return MonteCarloResult(0.0, 0.0, 0.0, 0, true);
            }
            
            // Volume per stratum
            Real stratum_volume = volume;
            for (int d = 0; d < N; ++d) {
                stratum_volume /= strata_per_dim;
            }
            
            // Total number of strata
            size_t total_strata = 1;
            for (int d = 0; d < N; ++d) {
                total_strata *= strata_per_dim;
            }
            
            Real sum = 0.0;
            Real sum_sq = 0.0;
            size_t total_samples = 0;
            
            // Iterate over all strata using N-digit base-strata_per_dim counter
            std::vector<int> stratum_idx(N, 0);
            
            for (size_t s = 0; s < total_strata; ++s) {
                // Compute stratum bounds
                VectorN<Real, N> stratum_lower, stratum_upper;
                for (int d = 0; d < N; ++d) {
                    stratum_lower[d] = lower[d] + stratum_idx[d] * delta[d];
                    stratum_upper[d] = stratum_lower[d] + delta[d];
                }
                
                // Sample within stratum
                Real stratum_sum = 0.0;
                VectorN<Real, N> point;
                
                for (int i = 0; i < samples_per_stratum; ++i) {
                    for (int d = 0; d < N; ++d) {
                        Real u = _uniform(_rng);
                        point[d] = stratum_lower[d] + u * delta[d];
                    }
                    stratum_sum += func(point);
                }
                
                Real stratum_mean = stratum_sum / samples_per_stratum;
                sum += stratum_mean;
                sum_sq += stratum_mean * stratum_mean;
                total_samples += samples_per_stratum;
                
                // Increment stratum index (N-digit counter)
                for (int d = 0; d < N; ++d) {
                    stratum_idx[d]++;
                    if (stratum_idx[d] < strata_per_dim) break;
                    stratum_idx[d] = 0;
                }
            }
            
            // Compute result
            Real mean = sum / total_strata;
            Real mean_sq = sum_sq / total_strata;
            Real variance = (mean_sq - mean * mean) / total_strata;
            Real std_error = std::sqrt(variance);
            
            Real integral = volume * mean;
            Real error = volume * std_error;
            
            bool converged = (std::abs(integral) < PrecisionValues<Real>::DivisionSafetyThreshold) || 
                             (error / std::abs(integral) < PrecisionValues<Real>::DefaultToleranceRelaxed);
            
            return MonteCarloResult(integral, error, variance, total_samples, converged);
        }
    };

    //////////////////////////////////////////////////////////////////////////
    /// @brief Hit-or-miss Monte Carlo for computing area/volume
    /// 
    /// Classic method: count fraction of random points that fall inside
    /// a region defined by an indicator function.
    //////////////////////////////////////////////////////////////////////////
    template<int N>
    class HitOrMissIntegrator
    {
    private:
        std::mt19937_64 _rng;
        std::uniform_real_distribution<Real> _uniform{0.0, 1.0};
        
    public:
        explicit HitOrMissIntegrator(unsigned int seed = 0)
        {
            if (seed == 0) {
                std::random_device rd;
                _rng.seed(rd());
            } else {
                _rng.seed(seed);
            }
        }
        
        void seed(unsigned int s) { _rng.seed(s); }
        
        /// @brief Estimate volume of region where indicator(x) returns true
        /// @param indicator Function returning true if point is inside region
        /// @param lower Lower bounds of bounding box
        /// @param upper Upper bounds of bounding box
        /// @param num_samples Number of random samples
        /// @return MonteCarloResult with estimated volume
        MonteCarloResult estimateVolume(
            std::function<bool(const VectorN<Real, N>&)> indicator,
            const VectorN<Real, N>& lower,
            const VectorN<Real, N>& upper,
            size_t num_samples = 100000,
            unsigned int seed = 0)
        {
            if (seed != 0) {
                _rng.seed(seed);
            }
            
            Real bounding_volume = 1.0;
            for (int d = 0; d < N; ++d) {
                bounding_volume *= (upper[d] - lower[d]);
            }
            
            size_t hits = 0;
            VectorN<Real, N> point;
            
            for (size_t i = 0; i < num_samples; ++i) {
                for (int d = 0; d < N; ++d) {
                    Real u = _uniform(_rng);
                    point[d] = lower[d] + u * (upper[d] - lower[d]);
                }
                if (indicator(point)) {
                    ++hits;
                }
            }
            
            Real p = static_cast<Real>(hits) / num_samples;
            Real volume = bounding_volume * p;
            
            // Binomial standard error
            Real variance = p * (1.0 - p);
            Real std_error = std::sqrt(variance / num_samples);
            Real error = bounding_volume * std_error;
            
            return MonteCarloResult(volume, error, variance, num_samples, true);
        }
    };

    //////////////////////////////////////////////////////////////////////////
    /// @brief Classic example: Estimate π using Monte Carlo
    /// 
    /// Textbook application: estimate π by computing area of unit circle
    /// inscribed in [-1,1]² square.
    //////////////////////////////////////////////////////////////////////////
    inline MonteCarloResult EstimatePi(size_t num_samples = 100000, unsigned int seed = 0)
    {
        HitOrMissIntegrator<2> integrator(seed);
        
        VectorN<Real, 2> lower, upper;
        lower[0] = lower[1] = -1.0;
        upper[0] = upper[1] = 1.0;
        
        auto inside_circle = [](const VectorN<Real, 2>& p) {
            return p[0]*p[0] + p[1]*p[1] <= 1.0;
        };
        
        auto result = integrator.estimateVolume(inside_circle, lower, upper, num_samples, seed);
        
        // Area of circle = π, so result directly gives π estimate
        return result;
    }
    
    //////////////////////////////////////////////////////////////////////////
    /// @brief Estimate volume of N-dimensional unit ball
    //////////////////////////////////////////////////////////////////////////
    template<int N>
    MonteCarloResult EstimateUnitBallVolume(size_t num_samples = 100000, unsigned int seed = 0)
    {
        HitOrMissIntegrator<N> integrator(seed);
        
        VectorN<Real, N> lower, upper;
        for (int d = 0; d < N; ++d) {
            lower[d] = -1.0;
            upper[d] = 1.0;
        }
        
        auto inside_ball = [](const VectorN<Real, N>& p) {
            Real r_sq = 0.0;
            for (int d = 0; d < N; ++d) {
                r_sq += p[d] * p[d];
            }
            return r_sq <= 1.0;
        };
        
        return integrator.estimateVolume(inside_ball, lower, upper, num_samples, seed);
    }

}  // namespace MML

#endif  // MML_MONTE_CARLO_INTEGRATION_H
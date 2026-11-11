///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        MonteCarloIntegration.h                                             ///
///  Description: Monte Carlo integration for high-dimensional integrals              ///
///               Plain, stratified sampling, and importance sampling methods         ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_MONTE_CARLO_INTEGRATION_H
#define MML_MONTE_CARLO_INTEGRATION_H

#include <mml/MMLBase.h>
#include <algorithm>
#include <mml/base/AlgorithmTypes.h>
#include <mml/core/Integration/IntegrationBase.h>
#include <mml/base/Random.h>
#include <mml/base/QuasiRandom.h>
#include <mml/interfaces/IFunction.h>
#include <mml/base/Vector/VectorN.h>

#include <numeric>
#include <cmath>
#include <cstddef>
#include <limits>
#include <type_traits>

namespace MML 
{
	namespace Detail
	{
		template<int N>
		inline Real ValidateMonteCarloBounds(const VectorN<Real, N>& lower, const VectorN<Real, N>& upper, const char* algorithmName)
		{
			Real volume = 1.0;
			bool zeroWidth = false;
			for (int d = 0; d < N; ++d) {
				if (!std::isfinite(lower[d]) || !std::isfinite(upper[d])) {
					throw ArgumentError(std::string(algorithmName) + ": bounds must be finite");
				}
				if (upper[d] < lower[d]) {
					throw ArgumentError(std::string(algorithmName) + ": upper bounds must be greater than or equal to lower bounds");
				}
				const Real width = upper[d] - lower[d];
				if (!std::isfinite(width))
					throw NumericalMethodError(std::string(algorithmName) + ": domain width is not representable");
				zeroWidth = zeroWidth || width == 0;
			}
			if (zeroWidth) return 0;
			for (int d = 0; d < N; ++d) {
				volume *= upper[d] - lower[d];
				if (!std::isfinite(volume) || volume == 0)
					throw NumericalMethodError(std::string(algorithmName) + ": domain volume is not representable");
			}
			return volume;
		}

		inline Real ClampMonteCarloVariance(Real variance)
		{
			return std::max(Real(0.0), variance);
		}

		/// Welford's online mean/variance (Knuth TAOCP vol. 2, 4.2.2). Numerically stable where
		/// E[f^2] - E[f]^2 suffers catastrophic cancellation for large |mean| relative to spread
		/// (analysis v4, CORE C7: the clamp to zero was hiding exactly that failure mode).
		struct WelfordAccumulator
		{
			// Double can represent every squared float deviation, including subnormals.
			using Moment = std::conditional_t<std::is_same_v<Real, float>, double, Real>;
			std::size_t count = 0;
			Real mean = 0.0;
			Moment m2 = 0.0;   ///< Sum of squared deviations from the running mean

			void add(Real x)
			{
				if (!std::isfinite(x))
					throw NumericalMethodError("Monte Carlo: non-finite integrand value");
				++count;
				const Real delta = x - mean;
				mean += delta / static_cast<Real>(count);
				m2 += Moment(delta) * Moment(x - mean);
				if (!std::isfinite(mean) || !std::isfinite(m2))
					throw NumericalMethodError("Monte Carlo: sample statistics overflow");
			}

			/// Chan-Golub-LeVeque parallel variance merge; includes between-block variation.
			void merge(const WelfordAccumulator& other) {
				if (other.count == 0) return;
				if (count == 0) { *this = other; return; }
				const size_t total = count + other.count;
				const Real delta = other.mean - mean;
				m2 += other.m2 + Moment(delta) * Moment(delta) * (Moment(count) * (Moment(other.count) / Moment(total)));
				mean += delta * (Real(other.count) / Real(total));
				count = total;
				if (!std::isfinite(mean) || !std::isfinite(m2))
					throw NumericalMethodError("Monte Carlo: merged sample statistics overflow");
			}

			/// Unbiased sample variance (n - 1); 0 for fewer than two samples
			Real sampleVariance() const {
				const Real variance = count > 1 ? Real(m2 / Moment(count - 1)) : Real(0);
				if (!std::isfinite(variance))
					throw NumericalMethodError("Monte Carlo: sample variance overflow");
				return variance;
			}
			/// Standard error of the mean, sqrt(s^2 / n)
			Real standardError() const { return count > 1 ? Real(std::sqrt(m2 / Moment(count - 1) / Moment(count))) : Real(0.0); }
		};
	}

	//////////////////////////////////////////////////////////////////////////
	/// @brief Result of a Monte Carlo integration
	/// @note For production code, check error_estimate, variance, and converged!
	//////////////////////////////////////////////////////////////////////////
	struct MonteCarloResult : public IntegrationResult
	{
		// value, error_estimate (standard error), converged, status, function_evaluations: IntegrationResult
		Real variance = 0.0;			 ///< Sample variance; stratified: variance of the unscaled mean estimator
		size_t samples_used = 0; ///< Number of samples actually used

		MonteCarloResult() { error_estimated = true; iterations = 1; }
	};

	//////////////////////////////////////////////////////////////////////////
	/// @brief Configuration for Monte Carlo integration
	//////////////////////////////////////////////////////////////////////////
	struct MonteCarloConfig {
		size_t num_samples = 100000; ///< Number of random samples
		Real target_error = Precision::MonteCarloTargetRelativeError;		 ///< Target relative error (for adaptive)
		unsigned int seed = 0;			 ///< Random seed (0 = use random_device)
		bool use_antithetic = false; ///< Use antithetic variates for variance reduction
		Real absolute_error = 0; ///< Absolute standard-error target; zero requests relative-only testing.

		MonteCarloConfig& samples(size_t n) {
			num_samples = n;
			return *this;
		}
		MonteCarloConfig& error(Real e) {
			target_error = e;
			return *this;
		}
		MonteCarloConfig& randomSeed(unsigned int s) {
			seed = s;
			return *this;
		}
		MonteCarloConfig& antithetic(bool a) {
			use_antithetic = a;
			return *this;
		}
	};

	namespace Detail {
		inline void ValidateMonteCarloConfig(const MonteCarloConfig& config) {
			if (config.num_samples > static_cast<size_t>(std::numeric_limits<int>::max()))
				throw ArgumentError("Monte Carlo sample count exceeds diagnostic counter range");
			if (config.num_samples < (config.use_antithetic ? size_t(4) : size_t(2)) ||
				(config.use_antithetic && config.num_samples % 2 != 0))
				throw ArgumentError("Monte Carlo requires two independent samples, or an even budget of at least four antithetic evaluations");
			if (!std::isfinite(config.target_error) || config.target_error <= 0 ||
				!std::isfinite(config.absolute_error) || config.absolute_error < 0)
				throw ArgumentError("Monte Carlo requires a positive finite relative target and nonnegative finite absolute target");
		}

		inline MonteCarloResult FinishMonteCarlo(const WelfordAccumulator& stats, Real volume,
			size_t evaluations, const MonteCarloConfig& config) {
			MonteCarloResult result;
			result.value = volume * stats.mean;
			result.error_estimate = volume * stats.standardError();
			result.variance = stats.sampleVariance();
			result.samples_used = evaluations;
			result.function_evaluations = static_cast<int>(evaluations);
			if (!std::isfinite(result.value) || !std::isfinite(result.error_estimate))
				throw NumericalMethodError("Monte Carlo: integral or standard error overflow");
			result.converged = result.error_estimate <= std::max(config.absolute_error,
				config.target_error * std::abs(result.value));
			result.status = result.converged ? AlgorithmStatus::Success : AlgorithmStatus::MaxIterationsExceeded;
			if (!result.converged)
				result.error_message = "Monte Carlo sample budget exhausted before the standard-error target was met";
			return result;
		}
	}

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
	class MonteCarloIntegrator {
	public:
		/// Initialize with optional seed (0 = random)
		explicit MonteCarloIntegrator(unsigned int seed = 0) {
			if (seed != 0) {
				Random::SetSeed(seed);
			}
		}

		/// Set the random seed
		void seed(unsigned int s) { Random::SetSeed(s); }

		/// @brief Integrate f over hyperrectangle [lower, upper]
		/// @param func The N-dimensional scalar function to integrate
		/// @param lower Lower bounds for each dimension
		/// @param upper Upper bounds for each dimension
		/// @param config Integration configuration
		/// @return MonteCarloResult with value, error estimate, and diagnostics
		MonteCarloResult integrate(const IScalarFunction<N>& func, const VectorN<Real, N>& lower, const VectorN<Real, N>& upper,
															 const MonteCarloConfig& config = MonteCarloConfig()) {
			AlgorithmTimer timer;

			Detail::ValidateMonteCarloConfig(config);

			// Compute volume of integration domain
			Real volume = Detail::ValidateMonteCarloBounds(lower, upper, "MonteCarloIntegrator");

			if (volume == 0.0) {
				MonteCarloResult r;
				// Measure-zero domain: the integral is exactly 0 by construction, nothing was sampled or checked
				r.converged = false;
				r.error_estimate = 0.0;   // exact, not a placeholder
				r.status = AlgorithmStatus::Completed;
				r.algorithm_name = "MonteCarloIntegrator";
				r.elapsed_time_ms = timer.elapsed_ms();
				return r;
			}

			// Reseed if specified
			if (config.seed != 0) {
				Random::SetSeed(config.seed);
			}

			// Welford running mean / variance (stable; see Detail::WelfordAccumulator)
			Detail::WelfordAccumulator stats;
			size_t n = config.num_samples;

			VectorN<Real, N> point;

			if (config.use_antithetic) {
				VectorN<Real, N> reflected;
				// Antithetic variates: use both x and (1-x) to reduce variance
				for (size_t i = 0; i < n; i += 2) {
					// Generate random point
					for (int d = 0; d < N; ++d) {
						Real u = Random::UniformReal(0.0, 1.0);
						point[d] = std::lerp(lower[d], upper[d], u);
						reflected[d] = std::lerp(lower[d], upper[d], Real(1) - u);
					}
					Real f1 = func(point);
					IntegrationDetail::RequireFinite(f1, "Monte Carlo integrand");
					Real f2 = func(reflected);
					IntegrationDetail::RequireFinite(f2, "Monte Carlo antithetic integrand");

					// Average the pair
					stats.add(f1 / Real(2) + f2 / Real(2));
				}
			} else {
				// Standard Monte Carlo
				for (size_t i = 0; i < n; ++i) {
					// Generate uniform random point in [lower, upper]
					for (int d = 0; d < N; ++d) {
						Real u = Random::UniformReal(0.0, 1.0);
						point[d] = lower[d] + u * (upper[d] - lower[d]);
					}

					stats.add(func(point));
				}
			}

			auto result = Detail::FinishMonteCarlo(stats, volume, config.num_samples, config);
			result.algorithm_name = "MonteCarloIntegrator";
			result.elapsed_time_ms = timer.elapsed_ms();
			return result;
		}

		/// @brief Convenience overload for unit hypercube [0,1]^N
		MonteCarloResult integrate(const IScalarFunction<N>& func, const MonteCarloConfig& config = MonteCarloConfig()) {
			VectorN<Real, N> lower, upper;
			for (int i = 0; i < N; ++i) {
				lower[i] = 0.0;
				upper[i] = 1.0;
			}
			return integrate(func, lower, upper, config);
		}

		/// @brief Quasi-Monte Carlo integration using a low-discrepancy sequence.
		/// @details Uses a Sobol sequence for N <= 6, otherwise a Halton sequence. Deterministic
		///          (no RNG). For smooth integrands this converges close to O(1/n) rather than the
		///          O(1/sqrt(n)) of plain Monte Carlo. The reported error_estimate is a heuristic
		///          sample-based value, not a rigorous QMC bound.
		MonteCarloResult integrateQuasiRandom(const IScalarFunction<N>& func, const VectorN<Real, N>& lower,
		                                       const VectorN<Real, N>& upper, size_t numSamples, bool useSobol = true) {
			AlgorithmTimer timer;
			if (numSamples == 0)
				throw ArgumentError("MonteCarloIntegrator: numSamples must be positive");
			if (numSamples > static_cast<size_t>(std::numeric_limits<int>::max()))
				throw ResourceLimitError("QuasiMonteCarloIntegrator: sample count exceeds diagnostic counter range");

			Real volume = Detail::ValidateMonteCarloBounds(lower, upper, "QuasiMonteCarloIntegrator");
			if (volume == 0.0) {
				MonteCarloResult r;
				// Measure-zero domain: the integral is exactly 0 by construction, nothing was sampled or checked
				r.converged = false;
				r.error_estimate = 0.0;   // exact, not a placeholder
				r.status = AlgorithmStatus::Completed;
				r.algorithm_name = "QuasiMonteCarloIntegrator";
				r.elapsed_time_ms = timer.elapsed_ms();
				return r;
			}

			Detail::WelfordAccumulator stats;
			VectorN<Real, N> u, point;

			auto accumulate = [&](auto& seq) {
				for (size_t i = 0; i < numSamples; ++i) {
					seq.Next(&u[0]);
					for (int d = 0; d < N; ++d)
						point[d] = std::lerp(lower[d], upper[d], u[d]);
					stats.add(func(point));
				}
			};

			if (useSobol && N <= SobolSequence::MaxDim) {
				SobolSequence seq(N);
				accumulate(seq);
			} else {
				HaltonSequence seq(N);
				accumulate(seq);
			}

			Real mean = stats.mean;
			Real variance = stats.sampleVariance();
			Real std_error = stats.standardError();

			MonteCarloResult result;
			result.value = volume * mean;
			result.error_estimate = volume * std_error;   // heuristic, not a rigorous QMC bound
			IntegrationDetail::RequireFinite(result.value, "QuasiMonteCarloIntegrator");
			if (numSamples > 1) {
				IntegrationDetail::RequireFinite(result.error_estimate, "QuasiMonteCarloIntegrator error");
				IntegrationDetail::RequireFinite(variance, "QuasiMonteCarloIntegrator variance");
			} else {
				result.error_estimated = false;
				result.error_estimate = NoEstimate;
			}
			result.variance = variance;
			result.samples_used = numSamples;
			result.function_evaluations = static_cast<int>(numSamples);
			// Fixed sample count: no tolerance was checked (RESULT_POLICY: Completed with a statistical estimate)
			result.converged = false;
			result.status = AlgorithmStatus::Completed;
			result.algorithm_name = "QuasiMonteCarloIntegrator";
			result.elapsed_time_ms = timer.elapsed_ms();
			return result;
		}
	};

	//////////////////////////////////////////////////////////////////////////
	/// @brief Convenience function for 1D Monte Carlo integration
	//////////////////////////////////////////////////////////////////////////
	inline MonteCarloResult IntegrateMonteCarlo1D(const IRealFunction& func, Real a, Real b,
																								const MonteCarloConfig& config = MonteCarloConfig()) {
		// Wrap IRealFunction in IScalarFunction<1>
		class Wrapper : public IScalarFunction<1> {
			const IRealFunction& _f;

		public:
			Wrapper(const IRealFunction& f)
					: _f(f) {}
			Real operator()(const VectorN<Real, 1>& x) const override { return _f(x[0]); }
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
	class StratifiedMonteCarloIntegrator {
	public:
		explicit StratifiedMonteCarloIntegrator(unsigned int seed = 0) {
			if (seed != 0) {
				Random::SetSeed(seed);
			}
		}

		void seed(unsigned int s) { Random::SetSeed(s); }

		/// @brief Integrate using stratified sampling
		/// @param func Function to integrate
		/// @param lower Lower bounds
		/// @param upper Upper bounds
		/// @param strata_per_dim Number of strata per dimension (total = strata^N)
		/// @param samples_per_stratum Samples per stratum
		MonteCarloResult integrate(const IScalarFunction<N>& func, const VectorN<Real, N>& lower, const VectorN<Real, N>& upper,
															 int strata_per_dim = 10, int samples_per_stratum = 10, unsigned int seed = 0) {
			AlgorithmTimer timer;

			if (strata_per_dim <= 0) {
				throw ArgumentError("StratifiedMonteCarloIntegrator: strata_per_dim must be positive");
			}
			if (samples_per_stratum <= 0) {
				throw ArgumentError("StratifiedMonteCarloIntegrator: samples_per_stratum must be positive");
			}

			// Compute total volume
			Real volume = Detail::ValidateMonteCarloBounds(lower, upper, "StratifiedMonteCarloIntegrator");

			if (volume == 0.0) {
				MonteCarloResult r;
				// Measure-zero domain: the integral is exactly 0 by construction, nothing was sampled or checked
				r.converged = false;
				r.error_estimate = 0.0;   // exact, not a placeholder
				r.status = AlgorithmStatus::Completed;
				r.algorithm_name = "StratifiedMonteCarloIntegrator";
				r.elapsed_time_ms = timer.elapsed_ms();
				return r;
			}

			// Total number of strata
			const size_t counterLimit = static_cast<size_t>(std::numeric_limits<int>::max());
			size_t total_strata = 1;
			for (int d = 0; d < N; ++d) {
				if (total_strata > counterLimit / static_cast<size_t>(strata_per_dim))
					throw ResourceLimitError("StratifiedMonteCarloIntegrator: stratum count exceeds diagnostic counter range");
				total_strata *= static_cast<size_t>(strata_per_dim);
			}
			if (total_strata > counterLimit / static_cast<size_t>(samples_per_stratum))
				throw ResourceLimitError("StratifiedMonteCarloIntegrator: sample count exceeds diagnostic counter range");
			const size_t total_samples = total_strata * static_cast<size_t>(samples_per_stratum);

			if (seed != 0) {
				Random::SetSeed(seed);
			}

			Real mean = 0;
			Real standardError = 0;

			// Iterate over all strata using N-digit base-strata_per_dim counter
			std::vector<int> stratum_idx(N, 0);

			for (size_t s = 0; s < total_strata; ++s) {
				// Compute stratum bounds
				VectorN<Real, N> stratum_lower, stratum_upper;
				for (int d = 0; d < N; ++d) {
					stratum_lower[d] = std::lerp(lower[d], upper[d], Real(stratum_idx[d]) / Real(strata_per_dim));
					stratum_upper[d] = std::lerp(lower[d], upper[d], Real(stratum_idx[d] + 1) / Real(strata_per_dim));
				}

				// Sample within stratum
				Detail::WelfordAccumulator stats;
				VectorN<Real, N> point;

				for (int i = 0; i < samples_per_stratum; ++i) {
					for (int d = 0; d < N; ++d) {
						Real u = Random::UniformReal(0.0, 1.0);
						point[d] = std::lerp(stratum_lower[d], stratum_upper[d], u);
					}
					stats.add(func(point));
				}

				mean = std::lerp(mean, stats.mean, Real(1) / Real(s + 1));
				// Var(mean) = sum_h s_h^2 / (n_h * H^2), for H equal-volume strata.
				// Combine standard errors before squaring to avoid premature range loss.
				if (samples_per_stratum > 1) {
					const Real stratumError = stats.standardError() / Real(total_strata);
					standardError = std::hypot(standardError, stratumError);
				}

				// Increment stratum index (N-digit counter)
				for (int d = 0; d < N; ++d) {
					stratum_idx[d]++;
					if (stratum_idx[d] < strata_per_dim)
						break;
					stratum_idx[d] = 0;
				}
			}

			MonteCarloResult result;
			result.value = volume * mean;
			IntegrationDetail::RequireFinite(result.value, "StratifiedMonteCarloIntegrator");
			result.samples_used = total_samples;
			result.function_evaluations = static_cast<int>(total_samples);
			result.algorithm_name = "StratifiedMonteCarloIntegrator";
			result.elapsed_time_ms = timer.elapsed_ms();
			if (samples_per_stratum == 1) {
				result.MarkFixedRule(static_cast<int>(total_samples));
				result.variance = NoEstimate;
				return result;
			}
			result.error_estimate = volume * standardError;
			result.variance = standardError * standardError;
			IntegrationDetail::RequireFinite(result.error_estimate, "StratifiedMonteCarloIntegrator error");
			IntegrationDetail::RequireFinite(result.variance, "StratifiedMonteCarloIntegrator variance");
			result.converged = result.error_estimate <= Precision::DefaultToleranceStrict * std::abs(result.value);
			result.status = result.converged ? AlgorithmStatus::Success : AlgorithmStatus::MaxIterationsExceeded;
			if (!result.converged)
				result.error_message = "Stratified Monte Carlo sample budget exhausted before the standard-error target was met";
			return result;
		}
	};

	//////////////////////////////////////////////////////////////////////////
	/// @brief Hit-or-miss Monte Carlo for computing area/volume
	///
	/// Classic method: count fraction of random points that fall inside
	/// a region defined by an indicator function.
	//////////////////////////////////////////////////////////////////////////
	template<int N>
	class HitOrMissIntegrator {
	public:
		explicit HitOrMissIntegrator(unsigned int seed = 0) {
			if (seed != 0) {
				Random::SetSeed(seed);
			}
		}

		void seed(unsigned int s) { Random::SetSeed(s); }

		/// @brief Estimate volume of region where indicator(x) returns true
		/// @param indicator Function returning true if point is inside region
		/// @param lower Lower bounds of bounding box
		/// @param upper Upper bounds of bounding box
		/// @param num_samples Number of random samples
		/// @return MonteCarloResult with estimated volume
		MonteCarloResult estimateVolume(std::function<bool(const VectorN<Real, N>&)> indicator, const VectorN<Real, N>& lower,
																		const VectorN<Real, N>& upper, size_t num_samples = 100000, unsigned int seed = 0) {
			AlgorithmTimer timer;

			if (num_samples == 0) {
				throw ArgumentError("HitOrMissIntegrator: num_samples must be positive");
			}

			if (seed != 0) {
				Random::SetSeed(seed);
			}

			Real bounding_volume = Detail::ValidateMonteCarloBounds(lower, upper, "HitOrMissIntegrator");

			size_t hits = 0;
			VectorN<Real, N> point;

			for (size_t i = 0; i < num_samples; ++i) {
				for (int d = 0; d < N; ++d) {
					Real u = Random::UniformReal(0.0, 1.0);
					point[d] = lower[d] + u * (upper[d] - lower[d]);
				}
				if (indicator(point)) {
					++hits;
				}
			}

			Real p = static_cast<Real>(hits) / num_samples;
			Real volume = bounding_volume * p;

			// Binomial standard error
			Real variance = Detail::ClampMonteCarloVariance(p * (MML::Real(1.0) - p));
			Real std_error = std::sqrt(variance / num_samples);
			Real error = bounding_volume * std_error;

			MonteCarloResult result;
			result.value = volume;
			result.error_estimate = error;
			result.variance = variance;
			result.samples_used = num_samples;
			// Fixed sample count: no tolerance was checked (RESULT_POLICY: Completed with a statistical estimate)
			result.converged = false;
			result.function_evaluations = static_cast<int>(num_samples);
			result.status = AlgorithmStatus::Completed;
			result.algorithm_name = "HitOrMissIntegrator";
			result.elapsed_time_ms = timer.elapsed_ms();
			return result;
		}
	};

	//////////////////////////////////////////////////////////////////////////
	/// @brief Classic example: Estimate π using Monte Carlo
	///
	/// Textbook application: estimate π by computing area of unit circle
	/// inscribed in [-1,1]² square.
	//////////////////////////////////////////////////////////////////////////
	inline MonteCarloResult EstimatePi(size_t num_samples = 100000, unsigned int seed = 0) {
		HitOrMissIntegrator<2> integrator(seed);

		VectorN<Real, 2> lower, upper;
		lower[0] = lower[1] = -1.0;
		upper[0] = upper[1] = 1.0;

		auto inside_circle = [](const VectorN<Real, 2>& p) { return p[0] * p[0] + p[1] * p[1] <= 1.0; };

		auto result = integrator.estimateVolume(inside_circle, lower, upper, num_samples, seed);

		// Area of circle = π, so result directly gives π estimate
		return result;
	}

	//////////////////////////////////////////////////////////////////////////
	/// @brief Estimate volume of N-dimensional unit ball
	//////////////////////////////////////////////////////////////////////////
	template<int N>
	MonteCarloResult EstimateUnitBallVolume(size_t num_samples = 100000, unsigned int seed = 0) {
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

} // namespace MML

#endif // MML_MONTE_CARLO_INTEGRATION_H
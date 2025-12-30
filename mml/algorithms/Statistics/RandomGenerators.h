///////////////////////////////////////////////////////////////////////////////////////////
// RandomGenerators.h
// 
// Random number generators (deviates) for MML (MinimalMathLibrary)
// Extracted from Statistics.h as part of refactoring to improve organization
//
// All generators use Mersenne Twister (std::mt19937_64) for high-quality randomness
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_RANDOM_GENERATORS_H
#define MML_RANDOM_GENERATORS_H

#include "../Statistics.h"
#include <random>

namespace MML
{
	namespace Statistics
	{
		/**
		 * @brief Exponential random deviate generator
		 * 
		 * Generates random numbers from an exponential distribution using the
		 * inverse transform method: X = -ln(U)/λ where U ~ Uniform(0,1)
		 * 
		 * Uses C++20 std::mt19937_64 for random number generation.
		 */
		class ExponentialDeviate
		{
		private:
			std::mt19937_64 _gen;                   // Mersenne Twister 64-bit generator
			std::uniform_real_distribution<Real> _uniform;  // Uniform(0,1) distribution
			Real _lambda;                            // Rate parameter

		public:
			/**
			 * @brief Construct an exponential deviate generator
			 * @param rate Rate parameter λ (must be > 0)
			 * @param seed Random seed (default: random_device)
			 */
			ExponentialDeviate(Real rate, uint64_t seed = std::random_device{}())
				: _gen(seed), _uniform(0.0, 1.0), _lambda(rate)
			{
				if (rate <= 0.0)
					throw StatisticsError("Rate parameter must be positive in ExponentialDeviate");
			}

			/**
			 * @brief Generate a random exponential deviate
			 * @return Random value from exponential distribution
			 */
			Real generate()
			{
				Real u;
				do {
					u = _uniform(_gen);
				} while (u == 0.0);  // Avoid log(0)
				return -std::log(u) / _lambda;
			}

			/// Alias for generate() to match common usage
			Real operator()() { return generate(); }
		};

		/**
		 * @brief Logistic random deviate generator
		 * 
		 * Generates random numbers from a logistic distribution using the
		 * inverse CDF method: X = μ + σ·ln(U/(1-U))
		 * 
		 * Uses C++20 std::mt19937_64 for random number generation.
		 */
		class LogisticDeviate
		{
		private:
			std::mt19937_64 _gen;
			std::uniform_real_distribution<Real> _uniform;
			Real _mu, _sigma;

		public:
			/**
			 * @brief Construct a logistic deviate generator
			 * @param location Location parameter μ
			 * @param scale Scale parameter σ (must be > 0)
			 * @param seed Random seed
			 */
			LogisticDeviate(Real location, Real scale, uint64_t seed = std::random_device{}())
				: _gen(seed), _uniform(0.0, 1.0), _mu(location), _sigma(scale)
			{
				if (scale <= 0.0)
					throw StatisticsError("Scale parameter must be positive in LogisticDeviate");
			}

			/**
			 * @brief Generate a random logistic deviate
			 * @return Random value from logistic distribution
			 */
			Real generate()
			{
				Real u;
				do {
					u = _uniform(_gen);
				} while (u == 0.0 || u == 1.0);  // Avoid log(0) and division by 0
				return _mu + _sigma * std::log(u / (1.0 - u));
			}

			Real operator()() { return generate(); }
		};

		/**
		 * @brief Normal (Gaussian) random deviate generator using Box-Muller transform
		 * 
		 * Box-Muller method generates pairs of independent normal variates from
		 * uniform random numbers. This implementation caches one value for efficiency.
		 * 
		 * Algorithm: If U1, U2 ~ Uniform(0,1), then
		 *   X = √(-2ln(R²)) · V1,  Y = √(-2ln(R²)) · V2
		 * where V1, V2 are uniform on unit circle, R² = V1² + V2²
		 */
		class NormalDeviateBoxMuller
		{
		private:
			std::mt19937_64 _gen;
			std::uniform_real_distribution<Real> _uniform;
			Real _mu, _sigma;
			Real _storedValue;  // Cached value from Box-Muller pair
			bool _hasStored;    // Whether cached value is available

		public:
			/**
			 * @brief Construct a normal deviate generator (Box-Muller method)
			 * @param mean Mean μ
			 * @param stddev Standard deviation σ (must be > 0)
			 * @param seed Random seed
			 */
			NormalDeviateBoxMuller(Real mean, Real stddev, uint64_t seed = std::random_device{}())
				: _gen(seed), _uniform(-1.0, 1.0), _mu(mean), _sigma(stddev),
				  _storedValue(0.0), _hasStored(false)
			{
				if (stddev <= 0.0)
					throw StatisticsError("Standard deviation must be positive in NormalDeviateBoxMuller");
			}

			/**
			 * @brief Generate a random normal deviate
			 * @return Random value from normal distribution
			 */
			Real generate()
			{
				if (_hasStored) {
					_hasStored = false;
					return _mu + _sigma * _storedValue;
				}

				Real v1, v2, rsq, fac;
				do {
					v1 = _uniform(_gen);
					v2 = _uniform(_gen);
					rsq = v1 * v1 + v2 * v2;
				} while (rsq >= 1.0 || rsq == 0.0);

				fac = std::sqrt(-2.0 * std::log(rsq) / rsq);
				_storedValue = v1 * fac;
				_hasStored = true;
				return _mu + _sigma * v2 * fac;
			}

			Real operator()() { return generate(); }
		};

		/**
		 * @brief Cauchy random deviate generator
		 * 
		 * Generates Cauchy variates using the ratio of two uniform random numbers.
		 * The Cauchy distribution has no defined mean or variance.
		 * 
		 * Algorithm: If V1, V2 uniform on unit circle, then X = V1/V2 ~ Cauchy(0,1)
		 */
		class CauchyDeviate
		{
		private:
			std::mt19937_64 _gen;
			std::uniform_real_distribution<Real> _uniform;
			Real _mu, _sigma;

		public:
			/**
			 * @brief Construct a Cauchy deviate generator
			 * @param location Location parameter μ
			 * @param scale Scale parameter σ (must be > 0)
			 * @param seed Random seed
			 */
			CauchyDeviate(Real location, Real scale, uint64_t seed = std::random_device{}())
				: _gen(seed), _uniform(-1.0, 1.0), _mu(location), _sigma(scale)
			{
				if (scale <= 0.0)
					throw StatisticsError("Scale parameter must be positive in CauchyDeviate");
			}

			/**
			 * @brief Generate a random Cauchy deviate
			 * @return Random value from Cauchy distribution
			 */
			Real generate()
			{
				Real v1, v2;
				do {
					v1 = _uniform(_gen);
					v2 = _uniform(_gen);
				} while (v1 * v1 + v2 * v2 >= 1.0 || v2 == 0.0);
				return _mu + _sigma * v1 / v2;
			}

			Real operator()() { return generate(); }
		};

		/**
		 * @brief Normal random deviate generator using Leva's ratio-of-uniforms method
		 * 
		 * Fast and accurate method for generating normal deviates. More efficient
		 * than Box-Muller for single values (no caching needed).
		 * 
		 * Uses a ratio-of-uniforms method with quick acceptance tests.
		 */
		class NormalDeviate
		{
		private:
			std::mt19937_64 _gen;
			std::uniform_real_distribution<Real> _uniform;
			Real _mu, _sigma;

		public:
			/**
			 * @brief Construct a normal deviate generator (Leva's method)
			 * @param mean Mean μ
			 * @param stddev Standard deviation σ (must be > 0)
			 * @param seed Random seed
			 */
			NormalDeviate(Real mean, Real stddev, uint64_t seed = std::random_device{}())
				: _gen(seed), _uniform(0.0, 1.0), _mu(mean), _sigma(stddev)
			{
				if (stddev <= 0.0)
					throw StatisticsError("Standard deviation must be positive in NormalDeviate");
			}

			/**
			 * @brief Generate a random normal deviate
			 * @return Random value from normal distribution
			 */
			Real generate()
			{
				Real u, v, x, y, q;
				do {
					u = _uniform(_gen);
					v = 1.7156 * (_uniform(_gen) - 0.5);
					x = u - 0.449871;
					y = std::abs(v) + 0.386595;
					q = x * x + y * (0.19600 * y - 0.25472 * x);
				} while (q > 0.27597 && 
				        (q > 0.27846 || v * v > -4.0 * std::log(u) * u * u));
				return _mu + _sigma * v / u;
			}

			Real operator()() { return generate(); }
		};

		/**
		 * @brief Gamma random deviate generator
		 * 
		 * Generates gamma distributed random numbers using Marsaglia and Tsang's method.
		 * For shape parameter α < 1, uses a transformation.
		 * 
		 * The gamma distribution is used in Bayesian statistics, queuing theory,
		 * and reliability analysis.
		 */
		class GammaDeviate
		{
		private:
			NormalDeviate _normalGen;
			std::mt19937_64 _gen;
			std::uniform_real_distribution<Real> _uniform;
			Real _alpha;      // Shape parameter
			Real _originalAlpha;  // Original shape (before transformation if α < 1)
			Real _beta;       // Scale parameter (1/rate)
			Real _a1, _a2;    // Precomputed constants

		public:
			/**
			 * @brief Construct a gamma deviate generator
			 * @param shape Shape parameter α (must be > 0)
			 * @param scale Scale parameter β (must be > 0, β = 1/rate)
			 * @param seed Random seed
			 */
			GammaDeviate(Real shape, Real scale, uint64_t seed = std::random_device{}())
				: _normalGen(0.0, 1.0, seed), _gen(seed), _uniform(0.0, 1.0),
				  _alpha(shape), _originalAlpha(shape), _beta(scale)
			{
				if (shape <= 0.0)
					throw StatisticsError("Shape parameter must be positive in GammaDeviate");
				if (scale <= 0.0)
					throw StatisticsError("Scale parameter must be positive in GammaDeviate");

				// For α < 1, use α' = α + 1 and later transform
				if (_alpha < 1.0)
					_alpha += 1.0;

				_a1 = _alpha - 1.0 / 3.0;
				_a2 = 1.0 / std::sqrt(9.0 * _a1);
			}

			/**
			 * @brief Generate a random gamma deviate
			 * @return Random value from gamma distribution
			 */
			Real generate()
			{
				Real u, v, x;
				do {
					do {
						x = _normalGen.generate();
						v = 1.0 + _a2 * x;
					} while (v <= 0.0);
					
					v = v * v * v;
					u = _uniform(_gen);
				} while (u > 1.0 - 0.331 * x * x * x * x &&
				         std::log(u) > 0.5 * x * x + _a1 * (1.0 - v + std::log(v)));

				// If original α < 1, apply transformation
				if (_alpha == _originalAlpha) {
					return _a1 * v / _beta;
				}
				else {
					do {
						u = _uniform(_gen);
					} while (u == 0.0);
					return std::pow(u, 1.0 / _originalAlpha) * _a1 * v / _beta;
				}
			}

			Real operator()() { return generate(); }
		};

		/**
		 * @brief Poisson random deviate generator
		 * 
		 * Generates Poisson distributed random integers. Uses different algorithms
		 * depending on mean λ:
		 * - λ < 5: Direct method (multiplicative algorithm)
		 * - λ ≥ 5: Ratio-of-uniforms method (faster for large λ)
		 * 
		 * The Poisson distribution models the number of events in a fixed interval.
		 */
		class PoissonDeviate
		{
		private:
			std::mt19937_64 _gen;
			std::uniform_real_distribution<Real> _uniform;
			Real _lambda;     // Mean parameter
			Real _sqrtLambda, _logLambda, _lambdaExp;
			Real _previousLambda;  // For caching
			std::vector<Real> _logFactorial;  // Cache of log(k!)

		public:
			/**
			 * @brief Construct a Poisson deviate generator
			 * @param mean Mean parameter λ (must be > 0)
			 * @param seed Random seed
			 */
			PoissonDeviate(Real mean, uint64_t seed = std::random_device{}())
				: _gen(seed), _uniform(0.0, 1.0), _lambda(mean),
				  _previousLambda(-1.0), _logFactorial(1024, -1.0)
			{
				if (mean <= 0.0)
					throw StatisticsError("Mean parameter must be positive in PoissonDeviate");
			}

			/**
			 * @brief Generate a random Poisson deviate
			 * @return Random integer from Poisson distribution
			 */
			int generate()
			{
				Real u, u2, v, v2, p, t, lfac;
				int k;

				if (_lambda < 5.0) {
					// Direct method for small λ
					if (_lambda != _previousLambda)
						_lambdaExp = std::exp(-_lambda);
					k = -1;
					t = 1.0;
					do {
						++k;
						t *= _uniform(_gen);
					} while (t > _lambdaExp);
				}
				else {
					// Ratio-of-uniforms method for large λ
					if (_lambda != _previousLambda) {
						_sqrtLambda = std::sqrt(_lambda);
						_logLambda = std::log(_lambda);
					}

					for (;;) {
						u = 0.64 * _uniform(_gen);
						v = -0.68 + 1.28 * _uniform(_gen);

						if (_lambda > 13.5) {
							v2 = v * v;
							if (v >= 0.0) {
								if (v2 > 6.5 * u * (0.64 - u) * (u + 0.2)) continue;
							}
							else {
								if (v2 > 9.6 * u * (0.66 - u) * (u + 0.07)) continue;
							}
						}

						k = static_cast<int>(std::floor(_sqrtLambda * (v / u) + _lambda + 0.5));
						if (k < 0) continue;

						u2 = u * u;
						if (_lambda > 13.5) {
							if (v >= 0.0) {
								if (v2 < 15.2 * u2 * (0.61 - u) * (0.8 - u)) break;
							}
							else {
								if (v2 < 6.76 * u2 * (0.62 - u) * (1.4 - u)) break;
							}
						}

						// Compute log(k!) using cache or lgamma
						if (k < 1024) {
							if (_logFactorial[k] < 0.0)
								_logFactorial[k] = std::lgamma(k + 1.0);
							lfac = _logFactorial[k];
						}
						else {
							lfac = std::lgamma(k + 1.0);
						}

						p = _sqrtLambda * std::exp(-_lambda + k * _logLambda - lfac);
						if (u2 < p) break;
					}
				}

				_previousLambda = _lambda;
				return k;
			}

			/**
			 * @brief Generate a Poisson deviate with different mean
			 * @param mean New mean parameter λ
			 * @return Random integer from Poisson(λ)
			 */
			int generate(Real mean)
			{
				_lambda = mean;
				return generate();
			}

			int operator()() { return generate(); }
		};

		/**
		 * @brief Binomial random deviate generator
		 * 
		 * Generates binomial distributed random integers B(n,p).
		 * Uses different algorithms based on n and p:
		 * - Small n: Bit comparison method
		 * - Medium n·p: Inverse CDF method  
		 * - Large n·p: Ratio-of-uniforms rejection method
		 * 
		 * Models the number of successes in n independent Bernoulli trials.
		 */
		class BinomialDeviate
		{
		private:
			std::mt19937_64 _gen;
			std::uniform_real_distribution<Real> _uniform;
			Real _pp, _p, _pb;          // Probability parameters
			Real _expNP, _np, _glnp, _plog, _pclog, _sq;
			int _n;                      // Number of trials
			int _method;                 // Which algorithm to use (0, 1, or 2)
			
			// For method 0 (small n)
			uint64_t _uz, _uo, _unfin, _diff, _rltp;
			int _pbits[5];
			
			// For method 1 (medium np)
			Real _cdf[64];
			
			// For method 2 (large np)
			Real _logFactorial[1024];

		public:
			/**
			 * @brief Construct a binomial deviate generator
			 * @param trials Number of trials n (must be > 0)
			 * @param probability Success probability p (must be in [0,1])
			 * @param seed Random seed
			 */
			BinomialDeviate(int trials, Real probability, uint64_t seed = std::random_device{}())
				: _gen(seed), _uniform(0.0, 1.0), _pp(probability), _n(trials)
			{
				if (trials <= 0)
					throw StatisticsError("Number of trials must be positive in BinomialDeviate");
				if (probability < 0.0 || probability > 1.0)
					throw StatisticsError("Probability must be in [0,1] in BinomialDeviate");

				// Work with p ≤ 0.5 for efficiency, flip at end if needed
				_pb = _p = (probability <= 0.5 ? probability : 1.0 - probability);

				if (_n <= 64) {
					// Method 0: Bit comparison for small n
					_uz = 0;
					_uo = 0xFFFFFFFFFFFFFFFFULL;
					_rltp = 0;
					for (int j = 0; j < 5; j++) {
						_pbits[j] = 1 & static_cast<int>(_pb *= 2.0);
					}
					_pb -= std::floor(_pb);
					_method = 0;
				}
				else if (_n * _p < 30.0) {
					// Method 1: Inverse CDF for medium np
					_cdf[0] = std::exp(_n * std::log(1.0 - _p));
					for (int j = 1; j < 64; j++) {
						_cdf[j] = _cdf[j - 1] + std::exp(
							std::lgamma(_n + 1.0) - std::lgamma(j + 1.0) - std::lgamma(_n - j + 1.0) +
							j * std::log(_p) + (_n - j) * std::log(1.0 - _p));
					}
					_method = 1;
				}
				else {
					// Method 2: Ratio-of-uniforms for large np
					_np = _n * _p;
					_glnp = std::lgamma(_n + 1.0);
					_plog = std::log(_p);
					_pclog = std::log(1.0 - _p);
					_sq = std::sqrt(_np * (1.0 - _p));
					if (_n < 1024) {
						for (int j = 0; j <= _n; j++)
							_logFactorial[j] = std::lgamma(j + 1.0);
					}
					_method = 2;
				}
			}

			/**
			 * @brief Generate a random binomial deviate
			 * @return Random integer from binomial distribution B(n,p)
			 */
			int generate()
			{
				int j, k, kl, km;
				Real y, u, v, u2, v2, b;

				if (_method == 0) {
					// Bit comparison method
					_unfin = _uo;
					for (j = 0; j < 5; j++) {
						_diff = _unfin & (_gen() ^ (_pbits[j] ? _uo : _uz));
						if (_pbits[j])
							_rltp |= _diff;
						else
							_rltp = _rltp & ~_diff;
						_unfin = _unfin & ~_diff;
					}
					k = 0;
					for (j = 0; j < _n; j++) {
						if (_unfin & 1) {
							if (_uniform(_gen) < _pb) ++k;
						}
						else {
							if (_rltp & 1) ++k;
						}
						_unfin >>= 1;
						_rltp >>= 1;
					}
				}
				else if (_method == 1) {
					// Inverse CDF method
					y = _uniform(_gen);
					kl = -1;
					k = 64;
					while (k - kl > 1) {
						km = (kl + k) / 2;
						if (y < _cdf[km])
							k = km;
						else
							kl = km;
					}
				}
				else {
					// Ratio-of-uniforms method
					for (;;) {
						u = 0.645 * _uniform(_gen);
						v = -0.63 + 1.25 * _uniform(_gen);
						v2 = v * v;

						if (v >= 0.0) {
							if (v2 > 6.5 * u * (0.645 - u) * (u + 0.2)) continue;
						}
						else {
							if (v2 > 8.4 * u * (0.645 - u) * (u + 0.1)) continue;
						}

						k = static_cast<int>(std::floor(_sq * (v / u) + _np + 0.5));
						if (k < 0) continue;

						u2 = u * u;
						if (v >= 0.0) {
							if (v2 < 12.25 * u2 * (0.615 - u) * (0.92 - u)) break;
						}
						else {
							if (v2 < 7.84 * u2 * (0.615 - u) * (1.2 - u)) break;
						}

						b = _sq * std::exp(_glnp + k * _plog + (_n - k) * _pclog -
							(_n < 1024 ? _logFactorial[k] + _logFactorial[_n - k]
							           : std::lgamma(k + 1.0) + std::lgamma(_n - k + 1.0)));
						if (u2 < b) break;
					}
				}

				// If we flipped p, flip the result
				if (_p != _pp)
					k = _n - k;

				return k;
			}

			int operator()() { return generate(); }
		};

	} // namespace Statistics
}  // namespace MML

#endif // MML_RANDOM_GENERATORS_H

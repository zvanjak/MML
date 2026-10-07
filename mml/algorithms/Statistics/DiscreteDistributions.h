///////////////////////////////////////////////////////////////////////////////////////////
// DiscreteDistributions.h
//
// Core discrete probability distributions for MML (MinimalMathLibrary).
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_STATISTICS_DISCRETE_DISTRIBUTIONS_H
#define MML_STATISTICS_DISCRETE_DISTRIBUTIONS_H

#include <mml/algorithms/Statistics.h>
#include <mml/base/SpecialFunctions.h>
#include <mml/base/Random.h>

#include <algorithm>
#include <cmath>
#include <limits>

namespace MML
{
	namespace Statistics
	{
		/// Bernoulli distribution for a single success/failure trial.
		struct BernoulliDistribution
		{
			Real p;

			/// Construct with success probability in [0, 1].
			BernoulliDistribution(Real probability) : p(probability)
			{
				if (p < 0.0 || p > 1.0)
					throw StatisticsError("Probability must be in [0, 1] in BernoulliDistribution");
			}

			/// @brief Draw a random 0/1 variate using the given URBG.
			template<class URBG> int sample(URBG& gen) const { return (std::uniform_real_distribution<Real>(0.0, 1.0)(gen) < p) ? 1 : 0; }
			/// @brief Draw a random variate using MML's thread-local engine.
			int sample() const { return sample(Random::Engine()); }

			/// Probability mass function P(X = k), where support is {0, 1}.
			Real pmf(int k) const
			{
				if (k == 0) return 1.0 - p;
				if (k == 1) return p;
				return 0.0;
			}

			/// Cumulative distribution function P(X <= k).
			Real cdf(int k) const
			{
				if (k < 0) return 0.0;
				if (k < 1) return 1.0 - p;
				return 1.0;
			}

			/// Return the smallest k such that CDF(k) >= probability.
			int inverseCdf(Real probability) const
			{
				if (probability <= 1.0 - p) return 0;
				return 1;
			}

			Real mean() const { return p; }
			Real variance() const { return p * (1.0 - p); }
		};

		/// Binomial distribution for the number of successes in independent Bernoulli trials.
		struct BinomialDistribution
		{
			int n;
			Real p;

			/// Construct with trial count n >= 0 and success probability in [0, 1].
			BinomialDistribution(int trials, Real probability)
				: n(trials), p(probability)
			{
				if (n < 0)
					throw StatisticsError("Number of trials must be non-negative in BinomialDistribution");
				if (p < 0.0 || p > 1.0)
					throw StatisticsError("Probability must be in [0, 1] in BinomialDistribution");
			}

			/// @brief Draw a random variate (sum of n Bernoulli trials) using the given URBG.
			template<class URBG> int sample(URBG& gen) const
			{
				std::uniform_real_distribution<Real> U(0.0, 1.0);
				int count = 0;
				for (int i = 0; i < n; ++i) if (U(gen) < p) ++count;
				return count;
			}
			/// @brief Draw a random variate using MML's thread-local engine.
			int sample() const { return sample(Random::Engine()); }

			/// Probability mass function P(X = k), where support is {0, ..., n}.
			Real pmf(int k) const
			{
				if (k < 0 || k > n) return 0.0;
				if (p == 0.0) return (k == 0) ? 1.0 : 0.0;
				if (p == 1.0) return (k == n) ? 1.0 : 0.0;

				Real logPmf = logBinomialCoeff(n, k)
					+ k * std::log(p)
					+ (n - k) * std::log(1.0 - p);
				return std::exp(logPmf);
			}

			/// Cumulative distribution function P(X <= k).
			Real cdf(int k) const
			{
				if (k < 0) return 0.0;
				if (k >= n) return 1.0;

				Real cumulative = 0.0;
				for (int i = 0; i <= k; ++i)
					cumulative += pmf(i);
				return cumulative;
			}

			/// Return the smallest k such that CDF(k) >= probability.
			int inverseCdf(Real probability) const
			{
				if (probability <= 0.0) return 0;
				if (probability >= 1.0) return n;

				Real cumulative = 0.0;
				for (int k = 0; k <= n; ++k) {
					cumulative += pmf(k);
					if (cumulative >= probability) return k;
				}
				return n;
			}

			Real mean() const { return n * p; }
			Real variance() const { return n * p * (1.0 - p); }
			Real stddev() const { return std::sqrt(variance()); }

		private:
			static Real logBinomialCoeff(int n, int k)
			{
				return std::lgamma(n + 1) - std::lgamma(k + 1) - std::lgamma(n - k + 1);
			}
		};

		/// Poisson distribution for counts with a positive event rate.
		struct PoissonDistribution
		{
			Real lambda;

			/// Construct with rate lambda > 0.
			PoissonDistribution(Real rate) : lambda(rate)
			{
				if (lambda <= 0.0)
					throw StatisticsError("Rate must be positive in PoissonDistribution");
			}

			/// @brief Draw a random variate (Knuth's algorithm) using the given URBG.
			template<class URBG> int sample(URBG& gen) const
			{
				Real L = std::exp(-lambda);
				int k = 0;
				Real prod = 1.0;
				std::uniform_real_distribution<Real> U(0.0, 1.0);
				do { ++k; prod *= U(gen); } while (prod > L);
				return k - 1;
			}
			/// @brief Draw a random variate using MML's thread-local engine.
			int sample() const { return sample(Random::Engine()); }

			/// Probability mass function P(X = k), where support is {0, 1, ...}.
			Real pmf(int k) const
			{
				if (k < 0) return 0.0;
				Real logPmf = k * std::log(lambda) - lambda - std::lgamma(k + 1);
				return std::exp(logPmf);
			}

			/// Cumulative distribution function P(X <= k).
			Real cdf(int k) const
			{
				if (k < 0) return 0.0;
				return 1.0 - SpecialFunctions::RegularizedGammaP(k + 1, lambda);
			}

			/// Return the smallest k such that CDF(k) >= probability.
			int inverseCdf(Real probability) const
			{
				if (probability <= 0.0) return 0;
				if (probability >= 1.0) return std::numeric_limits<int>::max();

				Real cumulative = 0.0;
				int maxIter = static_cast<int>(lambda + 10.0 * std::sqrt(lambda)) + 100;
				for (int k = 0; k <= maxIter; ++k) {
					cumulative += pmf(k);
					if (cumulative >= probability) return k;
				}
				return maxIter;
			}

			Real mean() const { return lambda; }
			Real variance() const { return lambda; }
			Real stddev() const { return std::sqrt(lambda); }
		};

		/// Geometric distribution for failures before the first success.
		struct GeometricDistribution
		{
			Real p;

			/// Construct with success probability in (0, 1].
			GeometricDistribution(Real probability) : p(probability)
			{
				if (p <= 0.0 || p > 1.0)
					throw StatisticsError("Probability must be in (0, 1] in GeometricDistribution");
			}

			/// @brief Draw a random variate (number of failures before first success) using the given URBG.
			template<class URBG> int sample(URBG& gen) const
			{
				if (p >= 1.0) return 0;
				Real u = std::uniform_real_distribution<Real>(0.0, 1.0)(gen);
				if (u < REAL(1e-300)) u = REAL(1e-300);
				return static_cast<int>(std::floor(std::log(u) / std::log(1.0 - p)));
			}
			/// @brief Draw a random variate using MML's thread-local engine.
			int sample() const { return sample(Random::Engine()); }

			/// Probability mass function P(X = k), where support is {0, 1, ...}.
			Real pmf(int k) const
			{
				if (k < 0) return 0.0;
				if (p == 1.0) return (k == 0) ? 1.0 : 0.0;
				return std::pow(1.0 - p, k) * p;
			}

			/// Cumulative distribution function P(X <= k).
			Real cdf(int k) const
			{
				if (k < 0) return 0.0;
				if (p == 1.0) return 1.0;
				return 1.0 - std::pow(1.0 - p, k + 1);
			}

			/// Return the smallest k such that CDF(k) >= probability.
			int inverseCdf(Real probability) const
			{
				if (probability <= 0.0) return 0;
				if (probability >= 1.0) return std::numeric_limits<int>::max();
				if (p == 1.0) return 0;

				return static_cast<int>(std::ceil(std::log(1.0 - probability) / std::log(1.0 - p) - 1.0));
			}

			Real mean() const { return (1.0 - p) / p; }
			Real variance() const { return (1.0 - p) / (p * p); }
			Real stddev() const { return std::sqrt(variance()); }
		};

		/// Negative binomial distribution for failures before the r-th success.
		struct NegativeBinomialDistribution
		{
			int r;
			Real p;

			/// Construct with required successes r > 0 and success probability in (0, 1].
			NegativeBinomialDistribution(int successes, Real probability)
				: r(successes), p(probability)
			{
				if (r <= 0)
					throw StatisticsError("Number of successes must be positive in NegativeBinomialDistribution");
				if (p <= 0.0 || p > 1.0)
					throw StatisticsError("Probability must be in (0, 1] in NegativeBinomialDistribution");
			}

			/// @brief Draw a random variate (sum of r geometric draws) using the given URBG.
			template<class URBG> int sample(URBG& gen) const
			{
				if (p >= 1.0) return 0;
				std::uniform_real_distribution<Real> U(0.0, 1.0);
				int total = 0;
				for (int i = 0; i < r; ++i) {
					Real u = U(gen);
					if (u < REAL(1e-300)) u = REAL(1e-300);
					total += static_cast<int>(std::floor(std::log(u) / std::log(1.0 - p)));
				}
				return total;
			}
			/// @brief Draw a random variate using MML's thread-local engine.
			int sample() const { return sample(Random::Engine()); }

			/// Probability mass function P(X = k), where support is {0, 1, ...}.
			Real pmf(int k) const
			{
				if (k < 0) return 0.0;
				if (p == 1.0) return (k == 0) ? 1.0 : 0.0;

				Real logPmf = std::lgamma(k + r) - std::lgamma(k + 1) - std::lgamma(r)
					+ r * std::log(p) + k * std::log(1.0 - p);
				return std::exp(logPmf);
			}

			/// Cumulative distribution function P(X <= k).
			Real cdf(int k) const
			{
				if (k < 0) return 0.0;

				Real cumulative = 0.0;
				for (int i = 0; i <= k; ++i)
					cumulative += pmf(i);
				return cumulative;
			}

			/// Return the smallest k such that CDF(k) >= probability.
			int inverseCdf(Real probability) const
			{
				if (probability <= 0.0) return 0;
				if (probability >= 1.0) return std::numeric_limits<int>::max();

				Real cumulative = 0.0;
				int maxIter = static_cast<int>(mean() + 10.0 * stddev()) + 100;
				for (int k = 0; k <= maxIter; ++k) {
					cumulative += pmf(k);
					if (cumulative >= probability) return k;
				}
				return maxIter;
			}

			Real mean() const { return r * (1.0 - p) / p; }
			Real variance() const { return r * (1.0 - p) / (p * p); }
			Real stddev() const { return std::sqrt(variance()); }
		};

		/// Hypergeometric distribution for successes in draws without replacement.
		struct HypergeometricDistribution
		{
			int N;
			int K;
			int n;

			/// Construct with population N, success states K, and draws n.
			HypergeometricDistribution(int population, int successStates, int draws)
				: N(population), K(successStates), n(draws)
			{
				if (N < 0)
					throw StatisticsError("Population must be non-negative in HypergeometricDistribution");
				if (K < 0 || K > N)
					throw StatisticsError("Success states must be in [0, N] in HypergeometricDistribution");
				if (n < 0 || n > N)
					throw StatisticsError("Draws must be in [0, N] in HypergeometricDistribution");
			}

			/// @brief Draw a random variate (sampling without replacement) using the given URBG.
			template<class URBG> int sample(URBG& gen) const
			{
				int successes = K, remaining = N, count = 0;
				std::uniform_real_distribution<Real> U(0.0, 1.0);
				for (int i = 0; i < n; ++i) {
					if (remaining <= 0) break;
					Real prob = static_cast<Real>(successes) / static_cast<Real>(remaining);
					if (U(gen) < prob) { ++count; --successes; }
					--remaining;
				}
				return count;
			}
			/// @brief Draw a random variate using MML's thread-local engine.
			int sample() const { return sample(Random::Engine()); }

			/// Probability mass function P(X = k), where support is max(0, n-N+K) to min(n, K).
			Real pmf(int k) const
			{
				int kMin = supportMin();
				int kMax = supportMax();
				if (k < kMin || k > kMax) return 0.0;

				Real logPmf = logBinomialCoeff(K, k)
					+ logBinomialCoeff(N - K, n - k)
					- logBinomialCoeff(N, n);
				return std::exp(logPmf);
			}

			/// Cumulative distribution function P(X <= k).
			Real cdf(int k) const
			{
				int kMin = supportMin();
				if (k < kMin) return 0.0;

				int kMax = supportMax();
				if (k >= kMax) return 1.0;

				Real cumulative = 0.0;
				for (int i = kMin; i <= k; ++i)
					cumulative += pmf(i);
				return cumulative;
			}

			/// Return the smallest k such that CDF(k) >= probability.
			int inverseCdf(Real probability) const
			{
				int kMin = supportMin();
				int kMax = supportMax();
				if (probability <= 0.0) return kMin;
				if (probability >= 1.0) return kMax;

				Real cumulative = 0.0;
				for (int k = kMin; k <= kMax; ++k) {
					cumulative += pmf(k);
					if (cumulative >= probability) return k;
				}
				return kMax;
			}

			Real mean() const
			{
				if (N == 0) return 0.0;
				return static_cast<Real>(n) * static_cast<Real>(K) / static_cast<Real>(N);
			}

			Real variance() const
			{
				if (N <= 1) return 0.0;
				Real population = static_cast<Real>(N);
				Real successes = static_cast<Real>(K);
				Real draws = static_cast<Real>(n);
				return draws * (successes / population) * (1.0 - successes / population)
					* ((population - draws) / (population - 1.0));
			}

			Real stddev() const { return std::sqrt(variance()); }

		private:
			int supportMin() const { return std::max(0, n - N + K); }
			int supportMax() const { return std::min(n, K); }

			static Real logBinomialCoeff(int n, int k)
			{
				if (k < 0 || k > n) return -std::numeric_limits<Real>::infinity();
				return std::lgamma(n + 1) - std::lgamma(k + 1) - std::lgamma(n - k + 1);
			}
		};
	} // namespace Statistics
} // namespace MML

#endif // MML_STATISTICS_DISCRETE_DISTRIBUTIONS_H
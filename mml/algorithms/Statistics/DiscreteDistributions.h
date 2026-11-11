///////////////////////////////////////////////////////////////////////////////////////////
// DiscreteDistributions.h
//
// Core discrete probability distributions for MML (MinimalMathLibrary).
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_STATISTICS_DISCRETE_DISTRIBUTIONS_H
#define MML_STATISTICS_DISCRETE_DISTRIBUTIONS_H

#include <mml/algorithms/Statistics.h>
#include <mml/MMLExceptions.h>
#include <mml/base/SpecialFunctions.h>
#include <mml/base/Random.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <type_traits>

namespace MML
{
	namespace Statistics
	{
		namespace StatisticsDetail
		{
			using DiscreteSamplingReal = std::common_type_t<Real, double>;

			inline int CheckedDiscreteSample(DiscreteSamplingReal value) {
				const auto count = std::floor(value);
				// The exclusive power-of-two bound stays exact even when Real is float.
				if (!std::isfinite(count) || count < 0
					|| count >= std::ldexp(DiscreteSamplingReal(1), std::numeric_limits<int>::digits))
					throw NumericalMethodError("Discrete sample exceeds the representable integer support");
				return static_cast<int>(count);
			}

			inline int CheckedGeometricSample(Real uniform, Real probability) {
				const Real fpmin = std::numeric_limits<Real>::min();
				const Real u = std::max(uniform, fpmin);
				return CheckedDiscreteSample(std::log(DiscreteSamplingReal(u))
					/ std::log1p(-DiscreteSamplingReal(probability)));
			}

			inline DiscreteSamplingReal PoissonSamplingLogPmf(DiscreteSamplingReal count, DiscreteSamplingReal rate) {
				using S = DiscreteSamplingReal;
				if (count < Precision::PoissonSamplingStirlingThreshold)
					return count * std::log(rate) - rate - std::lgamma(count + 1);
				const S relative = (count - rate) / rate;
				S deviance;
				if (std::abs(relative) <= Precision::PoissonSamplingDevianceLimit) {
					// (1+r)log(1+r)-r = sum_{n>=2} (-r)^n/[n(n-1)] near the mode.
					S power = relative * relative;
					S sum = 0;
					for (int n = 2; n <= Precision::PoissonSamplingDevianceTerms; ++n) {
						const S next = sum + power / (S(n) * (n - 1));
						if (next == sum) break;
						sum = next;
						power *= -relative;
					}
					deviance = rate * sum;
				} else {
					deviance = count * std::log(count / rate) + rate - count;
				}
				// Stirling's log-factorial correction through k^-11 avoids O(rate*log(rate)) cancellation.
				const S inv = 1 / count;
				const S inv2 = inv * inv;
				const S correction = inv * (S(1) / 12 + inv2 * (-S(1) / 360
					+ inv2 * (S(1) / 1260 + inv2 * (-S(1) / 1680
					+ inv2 * (S(1) / 1188 - inv2 * (S(691) / 360360))))));
				return -deviance - (std::log(S(2) * Constants::PI) + std::log(count)) / 2 - correction;
			}

			inline void RequireDiscreteQuantileInput(Real probability) {
				if (std::isnan(probability))
					throw StatisticsError("Probability must not be NaN in discrete inverse CDF");
			}

			inline Real CheckedDiscreteProbability(Real value) {
				if (!std::isfinite(value) || value < -Precision::EqualityTolerance || value > 1 + Precision::EqualityTolerance)
					throw NumericalMethodError("Invalid probability in discrete distribution");
				return std::clamp(value, Real(0), Real(1));
			}

			inline DiscreteSamplingReal LogBinomialCoefficient(int n, int k) {
				using Wide = DiscreteSamplingReal;
				if (k < 0 || k > n) return -std::numeric_limits<Wide>::infinity();
				k = std::min(k, n - k);
				if (k == 0) return 0;
				if (k == 1) return std::log(Wide(n));
				return std::lgamma(Wide(n) + 1) - std::lgamma(Wide(k) + 1) - std::lgamma(Wide(n - k) + 1);
			}

			template<class Cdf>
			int DiscreteQuantile(Real probability, Cdf cdf, int minimum, int maximum,
				int maxIterations = Precision::DiscreteQuantileMaxIterations) {
				auto evaluate = [&](int k) { return CheckedDiscreteProbability(cdf(k)); };
				if (evaluate(minimum) >= probability) return minimum;
				Real maximumProbability = evaluate(maximum);
				if (maximumProbability < probability)
					throw ConvergenceError("Discrete quantile exceeds the representable integer support", 0,
						probability - maximumProbability);
				int low = minimum;
				int high = maximum;
				int iteration = 0;
				for (; iteration < maxIterations && high - low > 1; ++iteration) {
					int mid = low + (high - low) / 2;
					if (evaluate(mid) >= probability)
						high = mid;
					else
						low = mid;
				}
				if (high - low > 1)
					throw ConvergenceError("Discrete quantile exhausted its integer search budget", iteration, Real(high - low));
				return high;
			}

			template<class Pmf>
			Real SumDiscreteProbabilities(int first, int last, Pmf pmf) {
				if (last < first) return 0;
				const int distance = last - first;
				if (distance >= Precision::DiscreteCdfMaxTerms)
					throw ConvergenceError("Discrete CDF exceeds its accumulation budget", 0, Real(distance) + 1);
				CompensatedSum sum;
				for (int k = first;; ++k) {
					sum.Add(CheckedDiscreteProbability(pmf(k)));
					if (k == last) break;
				}
				return CheckedDiscreteProbability(sum.Value());
			}
		}

		/// Bernoulli distribution for a single success/failure trial.
		struct BernoulliDistribution
		{
			Real p;

			/// Construct with success probability in [0, 1].
			BernoulliDistribution(Real probability) : p(probability)
			{
				ValidateParameters();
			}

			/// @brief Draw a random 0/1 variate using the given URBG.
			template<class URBG> int sample(URBG& gen) const { ValidateParameters(); return (std::uniform_real_distribution<Real>(0.0, 1.0)(gen) < p) ? 1 : 0; }
			/// @brief Draw a random variate using MML's thread-local engine.
			int sample() const { ValidateParameters(); return sample(Random::Engine()); }

			/// Probability mass function P(X = k), where support is {0, 1}.
			Real pmf(int k) const
			{
				ValidateParameters();
				if (k == 0) return MML::Real(1.0) - p;
				if (k == 1) return p;
				return 0.0;
			}

			/// Cumulative distribution function P(X <= k).
			Real cdf(int k) const
			{
				ValidateParameters();
				if (k < 0) return 0.0;
				if (k < 1) return MML::Real(1.0) - p;
				return 1.0;
			}

			/// Return the smallest k such that CDF(k) >= probability.
			int inverseCdf(Real probability) const
			{
				ValidateParameters();
				StatisticsDetail::RequireDiscreteQuantileInput(probability);
				if (probability <= 0) return 0;
				if (probability >= 1) return p == 0 ? 0 : 1;
				if (probability <= 1.0 - p) return 0;
				return 1;
			}

			Real mean() const { ValidateParameters(); return p; }
			Real variance() const { ValidateParameters(); return p * (MML::Real(1.0) - p); }

		private:
			void ValidateParameters() const {
				if (!(p >= 0 && p <= 1))
					throw StatisticsError("Invalid parameters in BernoulliDistribution");
			}
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
				ValidateParameters();
			}

			/// @brief Draw a random variate (sum of n Bernoulli trials) using the given URBG.
			template<class URBG> int sample(URBG& gen) const
			{
				ValidateParameters();
				std::uniform_real_distribution<Real> U(0.0, 1.0);
				int count = 0;
				for (int i = 0; i < n; ++i) if (U(gen) < p) ++count;
				return count;
			}
			/// @brief Draw a random variate using MML's thread-local engine.
			int sample() const { ValidateParameters(); return sample(Random::Engine()); }

			/// Probability mass function P(X = k), where support is {0, ..., n}.
			Real pmf(int k) const
			{
				ValidateParameters();
				if (k < 0 || k > n) return 0.0;
				if (p == MML::Real(0.0)) return (k == 0) ? MML::Real(1.0) : MML::Real(0.0);
				if (p == MML::Real(1.0)) return (k == n) ? MML::Real(1.0) : MML::Real(0.0);

				using Wide = StatisticsDetail::DiscreteSamplingReal;
				const Wide probability = p;
				const Wide logPmf = StatisticsDetail::LogBinomialCoefficient(n, k)
					+ Wide(k) * std::log(probability)
					+ Wide(n - k) * std::log1p(-probability);
				return static_cast<Real>(std::exp(logPmf));
			}

			/// Cumulative distribution function P(X <= k).
			Real cdf(int k) const
			{
				ValidateParameters();
				if (k < 0) return 0.0;
				if (k >= n) return 1.0;
				if (p == 0) return 1;
				if (p == 1) return 0;
				if (k == 0) return std::exp(Real(n) * std::log1p(-p));
				// Evaluate the upper beta tail when the CDF is near one, retaining tiny p.
				if (Real(k) + 1 > Real(n) * p)
					return StatisticsDetail::CheckedDiscreteProbability(1 - SpecialFunctions::RegularizedBetaI(p, Real(k) + 1, Real(n - k)));
				return SpecialFunctions::RegularizedBetaI(1 - p, Real(n - k), Real(k) + 1);
			}

			/// Return the smallest k such that CDF(k) >= probability.
			int inverseCdf(Real probability) const
			{
				ValidateParameters();
				StatisticsDetail::RequireDiscreteQuantileInput(probability);
				if (probability <= 0 || p == 0) return 0;
				if (probability >= 1 || p == 1) return n;
				return StatisticsDetail::DiscreteQuantile(probability, [this](int k) { return cdf(k); }, 0, n);
			}

			Real mean() const { ValidateParameters(); return n * p; }
			Real variance() const { ValidateParameters(); return n * p * (MML::Real(1.0) - p); }
			Real stddev() const { ValidateParameters(); return std::sqrt(variance()); }

		private:
			void ValidateParameters() const {
				if (!(n >= 0 && p >= 0 && p <= 1))
					throw StatisticsError("Invalid parameters in BinomialDistribution");
			}
		};

		/// Poisson distribution for counts with a positive event rate.
		struct PoissonDistribution
		{
			Real lambda;

			/// Construct with rate lambda > 0.
			PoissonDistribution(Real rate) : lambda(rate)
			{
				ValidateParameters();
			}

			/// @brief Draw a checked integer variate using bounded Knuth/PTRS sampling.
			template<class URBG> int sample(URBG& gen) const
			{
				ValidateSamplingParameters();
				using S = StatisticsDetail::DiscreteSamplingReal;
				std::uniform_real_distribution<Real> U(0.0, 1.0);
				if (lambda < Precision::PoissonSamplingSwitchRate) {
					const Real limit = std::exp(-lambda);
					Real product = 1;
					for (int k = 0; k < Precision::PoissonSamplingMaxIterations; ++k) {
						product *= U(gen);
						if (product <= limit) return k;
					}
				} else {
					// Hormann (1993), "The transformed rejection method for generating Poisson random variables" (PTRS).
					const S rate = lambda;
					const S b = S(0.931) + S(2.53) * std::sqrt(rate);
					const S a = -S(0.059) + S(0.02483) * b;
					const S inverseAlpha = S(1.1239) + S(1.1328) / (b - S(3.4));
					const S squeeze = S(0.9277) - S(3.6224) / (b - 2);
					for (int iteration = 0; iteration < Precision::PoissonSamplingMaxIterations; ++iteration) {
						const S u = S(U(gen)) - S(0.5);
						const S v = U(gen);
						const S width = S(0.5) - std::abs(u);
						if (width <= 0) continue;
						const S candidate = std::floor((2 * a / width + b) * u + rate + S(0.43));
						if (candidate < 0) continue;
						if (width >= Precision::PoissonSamplingSqueezeWidth && v <= squeeze)
							return StatisticsDetail::CheckedDiscreteSample(candidate);
						if (width < Precision::PoissonSamplingTailWidth && v > width) continue;
						const S logAcceptance = std::log(v) + std::log(inverseAlpha)
							- std::log(a / (width * width) + b);
						if (logAcceptance <= StatisticsDetail::PoissonSamplingLogPmf(candidate, rate))
							return StatisticsDetail::CheckedDiscreteSample(candidate);
					}
				}
				throw ConvergenceError("Poisson sampling exhausted its iteration budget",
					Precision::PoissonSamplingMaxIterations, 0);
			}
			/// @brief Draw a random variate using MML's thread-local engine.
			int sample() const { ValidateSamplingParameters(); return sample(Random::Engine()); }

			/// Probability mass function P(X = k), where support is {0, 1, ...}.
			Real pmf(int k) const
			{
				ValidateParameters();
				if (k < 0) return 0.0;
				Real logPmf = k * std::log(lambda) - lambda - std::lgamma(Real(k) + 1);
				return std::exp(logPmf);
			}

			/// Cumulative distribution function P(X <= k).
			Real cdf(int k) const
			{
				ValidateParameters();
				if (k < 0) return 0.0;
				return SpecialFunctions::RegularizedGammaQ(Real(k) + 1, lambda);
			}

			/// Return the smallest k such that CDF(k) >= probability.
			int inverseCdf(Real probability) const
			{
				ValidateParameters();
				StatisticsDetail::RequireDiscreteQuantileInput(probability);
				if (probability <= 0.0) return 0;
				if (probability >= 1.0) return std::numeric_limits<int>::max();
				return StatisticsDetail::DiscreteQuantile(probability, [this](int k) { return cdf(k); },
					0, std::numeric_limits<int>::max());
			}

			Real mean() const { ValidateParameters(); return lambda; }
			Real variance() const { ValidateParameters(); return lambda; }
			Real stddev() const { ValidateParameters(); return std::sqrt(lambda); }

		private:
			void ValidateSamplingParameters() const {
				ValidateParameters();
				if (static_cast<long double>(lambda) > static_cast<long double>(std::numeric_limits<int>::max()))
					throw NumericalMethodError("Poisson sampling rate exceeds INT_MAX");
			}

			void ValidateParameters() const {
				if (!(std::isfinite(lambda) && lambda > 0))
					throw StatisticsError("Invalid parameters in PoissonDistribution");
			}
		};

		/// Geometric distribution for failures before the first success.
		struct GeometricDistribution
		{
			Real p;

			/// Construct with success probability in (0, 1].
			GeometricDistribution(Real probability) : p(probability)
			{
				ValidateParameters();
			}

			/// @brief Draw a random variate (number of failures before first success) using the given URBG.
			template<class URBG> int sample(URBG& gen) const
			{
				ValidateParameters();
				if (p >= 1.0) return 0;
				Real u = std::uniform_real_distribution<Real>(0.0, 1.0)(gen);
				return StatisticsDetail::CheckedGeometricSample(u, p);
			}
			/// @brief Draw a random variate using MML's thread-local engine.
			int sample() const { ValidateParameters(); return sample(Random::Engine()); }

			/// Probability mass function P(X = k), where support is {0, 1, ...}.
			Real pmf(int k) const
			{
				ValidateParameters();
				if (k < 0) return 0.0;
				if (p == MML::Real(1.0)) return (k == 0) ? MML::Real(1.0) : MML::Real(0.0);
				return std::exp(Real(k) * std::log1p(-p)) * p;
			}

			/// Cumulative distribution function P(X <= k).
			Real cdf(int k) const
			{
				ValidateParameters();
				if (k < 0) return 0.0;
				if (p == 1.0) return 1.0;
				return -std::expm1((Real(k) + 1) * std::log1p(-p));
			}

			/// Return the smallest k such that CDF(k) >= probability.
			int inverseCdf(Real probability) const
			{
				ValidateParameters();
				StatisticsDetail::RequireDiscreteQuantileInput(probability);
				if (probability <= 0 || p == 1) return 0;
				if (probability >= 1.0) return std::numeric_limits<int>::max();
				return StatisticsDetail::DiscreteQuantile(probability, [this](int k) { return cdf(k); },
					0, std::numeric_limits<int>::max());
			}

			Real mean() const { ValidateParameters(); return (MML::Real(1.0) - p) / p; }
			Real variance() const { ValidateParameters(); return (MML::Real(1.0) - p) / (p * p); }
			Real stddev() const { ValidateParameters(); return std::sqrt(variance()); }

		private:
			void ValidateParameters() const {
				if (!(p > 0 && p <= 1))
					throw StatisticsError("Invalid parameters in GeometricDistribution");
			}
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
				ValidateParameters();
			}

			/// @brief Draw a random variate (sum of r geometric draws) using the given URBG.
			template<class URBG> int sample(URBG& gen) const
			{
				ValidateParameters();
				if (p >= 1.0) return 0;
				std::uniform_real_distribution<Real> U(0.0, 1.0);
				int total = 0;
				for (int i = 0; i < r; ++i) {
					Real u = U(gen);
					const int count = StatisticsDetail::CheckedGeometricSample(u, p);
					if (count > std::numeric_limits<int>::max() - total)
						throw NumericalMethodError("Negative-binomial sample sum exceeds INT_MAX");
					total += count;
				}
				return total;
			}
			/// @brief Draw a random variate using MML's thread-local engine.
			int sample() const { ValidateParameters(); return sample(Random::Engine()); }

			/// Probability mass function P(X = k), where support is {0, 1, ...}.
			Real pmf(int k) const
			{
				ValidateParameters();
				if (k < 0) return 0.0;
				if (p == MML::Real(1.0)) return (k == 0) ? MML::Real(1.0) : MML::Real(0.0);

				using Wide = StatisticsDetail::DiscreteSamplingReal;
				const Wide successes = r, failures = k, total = successes + failures;
				// Conditional Poisson identity avoids cancellation between O((k+r) log(k+r)) terms.
				const Wide logPmf = StatisticsDetail::PoissonSamplingLogPmf(successes, total * Wide(p))
					+ StatisticsDetail::PoissonSamplingLogPmf(failures, total * (Wide(1) - Wide(p)))
					- StatisticsDetail::PoissonSamplingLogPmf(total, total) + std::log(successes / total);
				return Real(std::exp(logPmf));
			}

			/// Cumulative distribution function P(X <= k).
			Real cdf(int k) const
			{
				ValidateParameters();
				if (k < 0) return 0.0;
				if (p == 1) return 1;
				if (r == 1) return -std::expm1((Real(k) + 1) * std::log1p(-p));
				return SpecialFunctions::RegularizedBetaI(p, Real(r), Real(k) + 1);
			}

			/// Return the smallest k such that CDF(k) >= probability.
			int inverseCdf(Real probability) const
			{
				ValidateParameters();
				StatisticsDetail::RequireDiscreteQuantileInput(probability);
				if (probability <= 0 || p == 1) return 0;
				if (probability >= 1.0) return std::numeric_limits<int>::max();
				return StatisticsDetail::DiscreteQuantile(probability, [this](int k) { return cdf(k); },
					0, std::numeric_limits<int>::max());
			}

			Real mean() const { ValidateParameters(); return r * (MML::Real(1.0) - p) / p; }
			Real variance() const { ValidateParameters(); return r * (MML::Real(1.0) - p) / (p * p); }
			Real stddev() const { ValidateParameters(); return std::sqrt(variance()); }

		private:
			void ValidateParameters() const {
				if (!(r > 0 && p > 0 && p <= 1))
					throw StatisticsError("Invalid parameters in NegativeBinomialDistribution");
			}
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
				ValidateParameters();
			}

			/// @brief Draw a random variate (sampling without replacement) using the given URBG.
			template<class URBG> int sample(URBG& gen) const
			{
				ValidateParameters();
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
			int sample() const { ValidateParameters(); return sample(Random::Engine()); }

			/// Probability mass function P(X = k), where support is max(0, n-N+K) to min(n, K).
			Real pmf(int k) const
			{
				ValidateParameters();
				int kMin = supportMin();
				int kMax = supportMax();
				if (k < kMin || k > kMax) return 0.0;

				const auto logPmf = StatisticsDetail::LogBinomialCoefficient(K, k)
					+ StatisticsDetail::LogBinomialCoefficient(N - K, n - k)
					- StatisticsDetail::LogBinomialCoefficient(N, n);
				return static_cast<Real>(std::exp(logPmf));
			}

			/// Cumulative distribution function P(X <= k).
			Real cdf(int k) const
			{
				ValidateParameters();
				int kMin = supportMin();
				if (k < kMin) return 0.0;

				int kMax = supportMax();
				if (k >= kMax) return 1.0;

				if (Real(k) < mean())
					return StatisticsDetail::SumDiscreteProbabilities(kMin, k, [this](int i) { return pmf(i); });
				return StatisticsDetail::CheckedDiscreteProbability(1
					- StatisticsDetail::SumDiscreteProbabilities(k + 1, kMax, [this](int i) { return pmf(i); }));
			}

			/// Return the smallest k such that CDF(k) >= probability.
			int inverseCdf(Real probability) const
			{
				ValidateParameters();
				StatisticsDetail::RequireDiscreteQuantileInput(probability);
				int kMin = supportMin();
				int kMax = supportMax();
				if (probability <= 0.0) return kMin;
				if (probability >= 1.0) return kMax;

				return StatisticsDetail::DiscreteQuantile(probability, [this](int k) { return cdf(k); }, kMin, kMax);
			}

			Real mean() const
			{
				ValidateParameters();
				if (N == 0) return 0.0;
				return static_cast<Real>(n) * static_cast<Real>(K) / static_cast<Real>(N);
			}

			Real variance() const
			{
				ValidateParameters();
				if (N <= 1) return 0.0;
				Real population = static_cast<Real>(N);
				Real successes = static_cast<Real>(K);
				Real draws = static_cast<Real>(n);
				return draws * (successes / population) * (Real(1) - successes / population)
					* ((population - draws) / (population - MML::Real(1.0)));
			}

			Real stddev() const { ValidateParameters(); return std::sqrt(variance()); }

		private:
			int supportMin() const { return std::max(0, n - N + K); }
			int supportMax() const { return std::min(n, K); }

			void ValidateParameters() const {
				if (!(N >= 0 && K >= 0 && K <= N && n >= 0 && n <= N))
					throw StatisticsError("Invalid parameters in HypergeometricDistribution");
			}
		};
	} // namespace Statistics
} // namespace MML

#endif // MML_STATISTICS_DISCRETE_DISTRIBUTIONS_H
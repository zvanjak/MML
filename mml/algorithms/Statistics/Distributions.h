///////////////////////////////////////////////////////////////////////////////////////////
// Distributions.h
// 
// Specialized probability distributions for MML (MinimalMathLibrary)
// Extracted from Statistics.h as part of refactoring to improve organization
//
// Contains: Cauchy, Exponential, and Logistic distributions
// Core distributions (Normal, T, ChiSquare, F) remain in Statistics.h
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_DISTRIBUTIONS_H
#define MML_DISTRIBUTIONS_H

#include "../Statistics.h"

namespace MML
{
	namespace Statistics
	{
		/**
		 * @brief Cauchy distribution (also known as Lorentz distribution)
		 * 
		 * The Cauchy distribution is a continuous probability distribution with
		 * no defined mean or variance. It has heavy tails and is used in physics
		 * and statistical mechanics.
		 * 
		 * PDF: f(x) = 1 / (π·σ·(1 + ((x-μ)/σ)²))
		 * where μ is the location parameter and σ is the scale parameter
		 */
		struct CauchyDistribution
		{
			Real mu;   // Location parameter (peak location)
			Real sigma; // Scale parameter (half-width at half-maximum)

			/**
			 * @brief Construct a Cauchy distribution
			 * @param location Location parameter (default: 0)
			 * @param scale Scale parameter (default: 1, must be > 0)
			 */
			CauchyDistribution(Real location = 0.0, Real scale = 1.0) 
				: mu(location), sigma(scale)
			{
				if (sigma <= 0.0)
					throw StatisticsError("Scale parameter must be positive in CauchyDistribution");
			}

			/**
			 * @brief Probability density function (PDF)
			 * @param x Value at which to evaluate the PDF
			 * @return Probability density at x
			 */
			Real pdf(Real x) const
			{
				Real z = (x - mu) / sigma;
				return 1.0 / (Constants::PI * sigma * (1.0 + z * z));
			}

			/**
			 * @brief Cumulative distribution function (CDF)
			 * @param x Value at which to evaluate the CDF
			 * @return Probability that X <= x
			 */
			Real cdf(Real x) const
			{
				return 0.5 + std::atan2(x - mu, sigma) / Constants::PI;
			}

			/**
			 * @brief Inverse cumulative distribution function (quantile function)
			 * @param p Probability (must be in (0, 1))
			 * @return Value x such that P(X <= x) = p
			 */
			Real inverseCdf(Real p) const
			{
				if (p <= 0.0 || p >= 1.0)
					throw StatisticsError("Probability must be in (0, 1) in CauchyDistribution::inverseCdf");
				return mu + sigma * std::tan(Constants::PI * (p - 0.5));
			}
		};

		/**
		 * @brief Exponential distribution
		 * 
		 * The exponential distribution models the time between events in a Poisson
		 * point process. It has the memoryless property.
		 * 
		 * PDF: f(x) = λ·exp(-λ·x) for x >= 0
		 * where λ is the rate parameter (λ = 1/mean)
		 */
		struct ExponentialDistribution
		{
			Real lambda; // Rate parameter (inverse of mean)

			/**
			 * @brief Construct an exponential distribution
			 * @param rate Rate parameter λ (must be > 0)
			 */
			ExponentialDistribution(Real rate) : lambda(rate)
			{
				if (lambda <= 0.0)
					throw StatisticsError("Rate parameter must be positive in ExponentialDistribution");
			}

			/**
			 * @brief Probability density function (PDF)
			 * @param x Value at which to evaluate the PDF (must be >= 0)
			 * @return Probability density at x
			 */
			Real pdf(Real x) const
			{
				if (x < 0.0)
					throw StatisticsError("x must be non-negative in ExponentialDistribution::pdf");
				return lambda * std::exp(-lambda * x);
			}

			/**
			 * @brief Cumulative distribution function (CDF)
			 * @param x Value at which to evaluate the CDF (must be >= 0)
			 * @return Probability that X <= x
			 */
			Real cdf(Real x) const
			{
				if (x < 0.0)
					throw StatisticsError("x must be non-negative in ExponentialDistribution::cdf");
				return 1.0 - std::exp(-lambda * x);
			}

			/**
			 * @brief Inverse cumulative distribution function (quantile function)
			 * @param p Probability (must be in [0, 1))
			 * @return Value x such that P(X <= x) = p
			 */
			Real inverseCdf(Real p) const
			{
				if (p < 0.0 || p >= 1.0)
					throw StatisticsError("Probability must be in [0, 1) in ExponentialDistribution::inverseCdf");
				return -std::log(1.0 - p) / lambda;
			}

			/**
			 * @brief Get the mean of the distribution
			 * @return Mean = 1/λ
			 */
			Real mean() const { return 1.0 / lambda; }

			/**
			 * @brief Get the variance of the distribution
			 * @return Variance = 1/λ²
			 */
			Real variance() const { return 1.0 / (lambda * lambda); }
		};

		/**
		 * @brief Logistic distribution
		 * 
		 * The logistic distribution is used in logistic regression and neural networks.
		 * It resembles the normal distribution but has heavier tails.
		 * 
		 * PDF: f(x) = exp(-z) / (σ·(1 + exp(-z))²)
		 * where z = (x - μ)/σ, μ is location, and σ is scale
		 */
		struct LogisticDistribution
		{
			Real mu;    // Location parameter (mean and median)
			Real sigma; // Scale parameter (related to variance by σ² = 3·s²/π²)

			/**
			 * @brief Construct a logistic distribution
			 * @param location Location parameter (default: 0)
			 * @param scale Scale parameter (default: 1, must be > 0)
			 */
			LogisticDistribution(Real location = 0.0, Real scale = 1.0)
				: mu(location), sigma(scale)
			{
				if (sigma <= 0.0)
					throw StatisticsError("Scale parameter must be positive in LogisticDistribution");
			}

			/**
			 * @brief Probability density function (PDF)
			 * @param x Value at which to evaluate the PDF
			 * @return Probability density at x
			 */
			Real pdf(Real x) const
			{
				// Standard logistic PDF: f(x) = e^(-z) / (σ(1 + e^(-z))²)
				// where z = (x - μ) / σ
				Real z = (x - mu) / sigma;
				Real exp_neg_z = std::exp(-z);
				Real denom = 1.0 + exp_neg_z;
				return exp_neg_z / (sigma * denom * denom);
			}

			/**
			 * @brief Cumulative distribution function (CDF)
			 * @param x Value at which to evaluate the CDF
			 * @return Probability that X <= x
			 */
			Real cdf(Real x) const
			{
				// Numerically stable computation using the logistic function
				Real z = (x - mu) / sigma;
				Real exp_z = std::exp(-std::abs(z));
				
				if (z >= 0.0)
					return 1.0 / (1.0 + exp_z);
				else
					return exp_z / (1.0 + exp_z);
			}

			/**
			 * @brief Inverse cumulative distribution function (quantile function)
			 * @param p Probability (must be in (0, 1))
			 * @return Value x such that P(X <= x) = p
			 */
			Real inverseCdf(Real p) const
			{
				if (p <= 0.0 || p >= 1.0)
					throw StatisticsError("Probability must be in (0, 1) in LogisticDistribution::inverseCdf");
				return mu + sigma * std::log(p / (1.0 - p));
			}

			/**
			 * @brief Get the mean of the distribution
			 * @return Mean = μ
			 */
			Real mean() const { return mu; }

			/**
			 * @brief Get the variance of the distribution
			 * @return Variance = (π·σ)²/3
			 */
			Real variance() const 
			{ 
				return (Constants::PI * sigma) * (Constants::PI * sigma) / 3.0;
			}
		};

	} // namespace Statistics
}  // namespace MML

#endif // MML_DISTRIBUTIONS_H

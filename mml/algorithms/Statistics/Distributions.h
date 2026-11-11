///////////////////////////////////////////////////////////////////////////////////////////
// Distributions.h
// 
// Probability distributions for MML (MinimalMathLibrary)
//
// Core:     Normal, Student's T, Chi-Square, F distributions
// Extended: Cauchy, Exponential, Logistic, Uniform, Gamma, Beta,
//           Weibull, Pareto, Log-Normal distributions
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_DISTRIBUTIONS_H
#define MML_DISTRIBUTIONS_H

#include <mml/algorithms/Statistics.h>
#include <mml/MMLExceptions.h>
#include <mml/base/SpecialFunctions.h>
#include <mml/base/Random.h>

namespace MML
{
	namespace Statistics
	{
		namespace StatisticsDetail
		{
			inline void RequireGammaParameters(Real shape, Real scale) {
				if (!(std::isfinite(shape) && shape > 0 && std::isfinite(scale) && scale > 0))
					throw StatisticsError("Gamma shape and scale must be finite and positive");
			}

			inline Real NormalCdf(Real z) {
				return Real(0.5) * std::erfc(-z / std::sqrt(Real(2)));
			}

			inline Real StandardizeNormal(Real x, Real mean, Real sigma) {
				Real difference = x - mean;
				if (std::isfinite(x) && std::isfinite(mean) && !std::isfinite(difference))
					return x / sigma - mean / sigma;
				return difference / sigma;
			}

			// Unlike general root solvers, require both relative tail accuracy and a relative x bracket.
			// Logarithmic bisection covers tiny positive roots without an arbitrary support clamp.
			template<class Probability>
			Real PositiveQuantile(Real target, Probability probability, bool increasing, Real high = 1,
				bool expand = true, int maxIterations = Precision::DistributionQuantileMaxIterations,
				int bracketIterations = Precision::DistributionQuantileBracketIterations) {
				auto evaluate = [&](Real x) {
					Real value = probability(x);
					if (!(std::isfinite(value) && value >= 0 && value <= 1))
						throw NumericalMethodError("Invalid probability during distribution quantile search");
					return value;
				};
				auto belowRoot = [&](Real value) { return increasing ? value < target : value > target; };
				const Real minimum = std::numeric_limits<Real>::denorm_min() > 0
					? std::numeric_limits<Real>::denorm_min() : std::numeric_limits<Real>::min();
				Real low = minimum;
				Real lowValue = evaluate(low);
				Real highValue = evaluate(high);
				if (belowRoot(highValue) && expand) {
					int iteration = 0;
					for (; iteration < bracketIterations && belowRoot(highValue); ++iteration) {
						if (high == std::numeric_limits<Real>::max())
							break;
						low = high;
						lowValue = highValue;
						high = high > std::numeric_limits<Real>::max() / 2
							? std::numeric_limits<Real>::max() : 2 * high;
						highValue = evaluate(high);
					}
					if (belowRoot(highValue))
						throw ConvergenceError("Cannot bracket a finite distribution quantile", iteration,
							std::abs(highValue / target - 1));
				}
				if (belowRoot(highValue) || (lowValue != target && !belowRoot(lowValue)))
					throw ConvergenceError("Distribution quantile lies outside the resolvable positive range", 0,
						std::abs(lowValue / target - 1));

				Real logLow = std::log(low);
				Real logHigh = std::log(high);
				Real residual = std::numeric_limits<Real>::infinity();
				for (int iteration = 0; iteration < maxIterations; ++iteration) {
					Real logX = std::midpoint(logLow, logHigh);
					Real x = std::clamp(std::exp(logX), low, high);
					Real value = evaluate(x);
					residual = std::abs(value / target - 1);
					if (residual <= Precision::DistributionQuantileProbabilityTolerance
						&& logHigh - logLow <= std::log1p(Precision::DistributionQuantileRelativeTolerance)) {
						Real probabilityResolution = (std::nextafter(value, std::numeric_limits<Real>::infinity()) - value) / target;
						if (probabilityResolution > Precision::DistributionQuantileProbabilityTolerance)
							throw ConvergenceError("Distribution quantile tail has insufficient probability resolution",
								iteration + 1, probabilityResolution);
						return x;
					}
					if (logX == logLow || logX == logHigh)
						throw ConvergenceError("Distribution quantile stagnated before reaching tail accuracy",
							iteration + 1, residual);
					if (belowRoot(value))
						logLow = logX;
					else
						logHigh = logX;
				}
				throw ConvergenceError("Distribution quantile exhausted its iteration budget", maxIterations, residual);
			}

			inline Real StudentCentralProbability(Real x, Real nu) {
				Real ratio = x / std::sqrt(nu);
				Real square = ratio * ratio;
				return SpecialFunctions::RegularizedBetaI(square / (1 + square), Real(0.5), nu / 2);
			}

			inline Real StudentTwoSidedTail(Real x, Real nu) {
				if (nu == 1)
					return 2 * std::atan(1 / x) / Constants::PI;
				if (x <= 1)
					return 1 - StudentCentralProbability(x, nu);
				Real ratio = std::sqrt(nu) / x;
				Real square = ratio * ratio;
				return SpecialFunctions::RegularizedBetaI(square / (1 + square), nu / 2, Real(0.5));
			}

			inline Real StudentCriticalValue(Real alpha, Real nu) {
				if (alpha > Real(0.5))
					return PositiveQuantile(1 - alpha,
						[=](Real x) { return StudentCentralProbability(x, nu); }, true);
				return PositiveQuantile(alpha,
					[=](Real x) { return StudentTwoSidedTail(x, nu); }, false);
			}
		}

		/// @brief Inverse of the standard normal CDF (probit): z such that Phi(z)=p.
		/// @details Uses the smaller probability directly via the shared complementary error-function inverse.
		inline Real InverseStandardNormalCdf(Real p)
		{
			if (!(p > 0.0 && p < 1.0))
				throw StatisticsError("Probability must be in (0, 1) in InverseStandardNormalCdf");
			Real tail = p < Real(0.5) ? p : 1 - p;
			Real magnitude = std::sqrt(Real(2)) * SpecialFunctions::ErfcInv(2 * tail);
			return p < Real(0.5) ? -magnitude : magnitude;
		}

		/// @brief Draw U ~ Uniform(0,1) from a uniform random bit generator.
		template<class URBG>
		inline Real SampleUniform01(URBG& gen)
		{
			return std::uniform_real_distribution<Real>(0.0, 1.0)(gen);
		}

		/// @brief Draw Z ~ N(0,1) via the Box-Muller transform.
		template<class URBG>
		inline Real SampleStandardNormal(URBG& gen)
		{
			Real u1;
			do { u1 = SampleUniform01(gen); } while (u1 == static_cast<MML::Real>(0.0));
			Real u2 = SampleUniform01(gen);
			return std::sqrt(-MML::Real(2.0) * std::log(u1)) * std::cos(MML::Real(2.0) * Constants::PI * u2);
		}

		/// @brief Draw X ~ Gamma(shape, scale) via the Marsaglia-Tsang method.
		template<class URBG>
		inline Real SampleGamma(URBG& gen, Real shape, Real scale)
		{
			StatisticsDetail::RequireGammaParameters(shape, scale);
			if (shape < 1.0)
			{
				Real u;
				do { u = SampleUniform01(gen); } while (u == static_cast<MML::Real>(0.0));
				return SampleGamma(gen, shape + MML::Real(1.0), scale) * std::pow(u, MML::Real(1.0) / shape);
			}
			Real d = shape - MML::Real(1.0) / MML::Real(3.0);
			Real c = MML::Real(1.0) / std::sqrt(MML::Real(9.0) * d);
			while (true)
			{
				Real x, v;
				do { x = SampleStandardNormal(gen); v = MML::Real(1.0) + c * x; } while (v <= MML::Real(0.0));
				v = v * v * v;
				Real u = SampleUniform01(gen);
				if (u < 1.0 - static_cast<MML::Real>(0.0331) * x * x * x * x) return d * v * scale;
				if (std::log(u) < 0.5 * x * x + d * (1.0 - v + std::log(v))) return d * v * scale;
			}
		}

		/// @brief Normal (Gaussian) distribution
		/// The normal distribution is the most important continuous distribution,
		/// central to the Central Limit Theorem and statistical inference.
		/// PDF: f(x) = (1 / (σ√(2π))) · exp(-½((x-μ)/σ)²)
		/// where μ is the mean and σ is the standard deviation
		/// The standard normal has μ=0 and σ=1.
		struct NormalDistribution
		{
			Real mu;    // Mean (location parameter)
			Real sigma; // Standard deviation (scale parameter)

			/// @brief Construct a normal distribution
			///
			/// @param mean Mean parameter μ (default: 0)
			/// @param stddev Standard deviation σ (default: 1, must be > 0)
			NormalDistribution(Real mean = 0.0, Real stddev = 1.0)
				: mu(mean), sigma(stddev)
			{
				ValidateParameters();
			}

			/// @brief Draw a random variate using the given URBG (e.g. std::mt19937).
			template<class URBG> Real sample(URBG& gen) const { ValidateParameters(); return mu + sigma * SampleStandardNormal(gen); }
			/// @brief Draw a random variate using MML's thread-local engine.
			Real sample() const { ValidateParameters(); return sample(Random::Engine()); }

			/// @brief Probability density function (PDF)
			///
			/// @param x Value at which to evaluate the PDF
			/// @return Probability density at x
			Real pdf(Real x) const
			{
				ValidateParameters();
				Real z = (x - mu) / sigma;
				return std::exp(-MML::Real(0.5) * z * z) / (sigma * std::sqrt(MML::Real(2.0) * Constants::PI));
			}

			/// @brief Cumulative distribution function (CDF)
			///
			/// Uses the complementary error function to preserve the small lower tail.
			/// @param x Value at which to evaluate the CDF
			///
			/// @return Probability that X <= x
			Real cdf(Real x) const
			{
				ValidateParameters();
				return StatisticsDetail::NormalCdf(StatisticsDetail::StandardizeNormal(x, mu, sigma));
			}

			/// @brief Inverse cumulative distribution function (quantile/probit function)
			///
			/// Uses the shared tail-aware inverse standard normal CDF.
			/// @param p Probability (must be in (0, 1))
			///
			/// @return Value x such that P(X <= x) = p
			Real inverseCdf(Real p) const
			{
				ValidateParameters();
				if (!(p > 0.0 && p < 1.0))
					throw StatisticsError("Probability must be in (0, 1) in NormalDistribution::inverseCdf");

				Real z = InverseStandardNormalCdf(p);
				return mu + sigma * z;
			}

			/// @brief Z-score: standardize a value
			///
			/// @param x Raw value
			/// @return Standardized value (z-score)
			Real zScore(Real x) const
			{
				ValidateParameters();
				return StatisticsDetail::StandardizeNormal(x, mu, sigma);
			}

			/// @brief Convert z-score back to raw value
			///
			/// @param z Standardized value
			/// @return Raw value
			Real fromZScore(Real z) const
			{
				ValidateParameters();
				return mu + sigma * z;
			}

			/// @brief Mean of the distribution
			Real mean() const { ValidateParameters(); return mu; }

			/// @brief Variance of the distribution
			Real variance() const {
				ValidateParameters();
				int exponent = 0;
				const Real mantissa = std::frexp(sigma, &exponent);
				return std::scalbn(mantissa * mantissa, 2 * exponent);
			}

			/// @brief Standard deviation of the distribution
			Real stddev() const { ValidateParameters(); return sigma; }

		private:
			void ValidateParameters() const {
				if (!(std::isfinite(mu) && std::isfinite(sigma) && sigma > 0))
					throw StatisticsError("Invalid parameters in NormalDistribution");
			}
		};

		/// @brief Standard normal distribution (mean=0, stddev=1)
		///
		/// Convenience functions for the standard normal distribution.
		struct StandardNormal
		{
			/// @brief PDF of standard normal
			static Real pdf(Real z)
			{
				return std::exp(-MML::Real(0.5) * z * z) / std::sqrt(MML::Real(2.0) * Constants::PI);
			}

			/// @brief CDF of standard normal (Φ function)
			static Real cdf(Real z)
			{
				return StatisticsDetail::NormalCdf(z);
			}

			/// @brief Inverse CDF (probit function)
			static Real inverseCdf(Real p)
			{
				NormalDistribution standard(0.0, 1.0);
				return standard.inverseCdf(p);
			}

			/// @brief Survival function: P(Z > z) = 1 - Φ(z)
			static Real sf(Real z)
			{
				return StatisticsDetail::NormalCdf(-z);
			}

			/// @brief Two-tailed p-value: P(|Z| > |z|)
			static Real twoTailedPValue(Real z)
			{
				return MML::Real(2.0) * sf(std::abs(z));
			}
		};

		/// @brief Student's t-distribution
		///
		/// The t-distribution is used for inference about means when the sample size
		/// is small and/or the population standard deviation is unknown. It has heavier
		///
		/// tails than the normal distribution, with the tail weight controlled by
		/// degrees of freedom (df).
		///
		/// As df → ∞, the t-distribution approaches the standard normal distribution.
		/// PDF: f(x) = Γ((ν+1)/2) / (√(νπ) · Γ(ν/2)) · (1 + x²/ν)^(-(ν+1)/2)
		///
		/// where ν is the degrees of freedom
		struct TDistribution
		{
			int df;  // Degrees of freedom

			/// @brief Construct a Student's t-distribution
			///
			/// @param degreesOfFreedom Degrees of freedom (must be >= 1)
			TDistribution(int degreesOfFreedom) : df(degreesOfFreedom)
			{
				ValidateParameters();
			}

			/// @brief Draw a random variate (Z / sqrt(chi2_df / df)) using the given URBG.
			template<class URBG> Real sample(URBG& gen) const
			{
				ValidateParameters();
				Real z = SampleStandardNormal(gen);
				Real v = SampleGamma(gen, df / MML::Real(2.0), MML::Real(2.0));
				return z / std::sqrt(v / df);
			}
			/// @brief Draw a random variate using MML's thread-local engine.
			Real sample() const { ValidateParameters(); return sample(Random::Engine()); }

			/// @brief Probability density function (PDF)
			///
			/// @param x Value at which to evaluate the PDF
			/// @return Probability density at x
			Real pdf(Real x) const
			{
				ValidateParameters();
				Real nu = static_cast<Real>(df);
				
				// Log-space computation for numerical stability
				Real logNumerator = std::lgamma((nu + MML::Real(1.0)) / MML::Real(2.0));
				Real logDenominator = MML::Real(0.5) * std::log(nu * Constants::PI) + std::lgamma(nu / MML::Real(2.0));
				Real logBase = -(nu + MML::Real(1.0)) / MML::Real(2.0) * std::log(MML::Real(1.0) + x * x / nu);
				
				return std::exp(logNumerator - logDenominator + logBase);
			}

			/// @brief Cumulative distribution function (CDF)
			///
			/// Uses the regularized incomplete beta function.
			/// @param x Value at which to evaluate the CDF
			///
			/// @return Probability that T <= x
			Real cdf(Real x) const
			{
				ValidateParameters();
				Real tail = Real(0.5) * StatisticsDetail::StudentTwoSidedTail(std::abs(x), Real(df));
				return x < 0 ? tail : 1 - tail;
			}

			/// @brief Inverse cumulative distribution function (quantile function)
			///
			/// Uses iterative root finding for the inverse.
			/// @param p Probability (must be in (0, 1))
			///
			/// @return Value t such that P(T <= t) = p
			Real inverseCdf(Real p) const
			{
				ValidateParameters();
				if (!(p > 0 && p < 1))
					throw StatisticsError("Probability must be in (0, 1) in TDistribution::inverseCdf");

				if (p == 0.5)
					return 0.0;

				Real tail = p < Real(0.5) ? p : 1 - p;
				Real magnitude = StatisticsDetail::StudentCriticalValue(2 * tail, Real(df));
				return p < Real(0.5) ? -magnitude : magnitude;
			}

			/// @brief Survival function: P(T > t)
			Real sf(Real x) const
			{
				ValidateParameters();
				Real tail = Real(0.5) * StatisticsDetail::StudentTwoSidedTail(std::abs(x), Real(df));
				return x < 0 ? 1 - tail : tail;
			}

			/// @brief Two-tailed p-value: P(|T| > |t|)
			///
			/// This is the standard p-value used in t-tests.
			Real twoTailedPValue(Real t) const
			{
				ValidateParameters();
				return StatisticsDetail::StudentTwoSidedTail(std::abs(t), Real(df));
			}

			/// @brief Critical value for given significance level (two-tailed)
			///
			/// Returns t* such that P(|T| > t*) = alpha
			/// @param alpha Significance level (e.g., 0.05 for 95% confidence)
			///
			/// @return Critical value t*
			Real criticalValue(Real alpha) const
			{
				ValidateParameters();
				if (!(alpha > 0 && alpha < 1))
					throw StatisticsError("Alpha must be in (0, 1) in TDistribution::criticalValue");
				
				return StatisticsDetail::StudentCriticalValue(alpha, Real(df));
			}

			/// @brief Mean of the distribution (if it exists)
			///
			/// Mean exists only for df >= 2, and equals 0.
			Real mean() const
			{
				ValidateParameters();
				if (df < 2)
					throw StatisticsError("Mean does not exist for df < 2 in TDistribution");
				return 0.0;
			}

			/// @brief Variance of the distribution (if it exists)
			///
			/// Variance exists only for df >= 3.
			/// Var(T) = ν/(ν-2) for ν > 2
			Real variance() const
			{
				ValidateParameters();
				if (df <= 2)
					throw StatisticsError("Variance does not exist for df <= 2 in TDistribution");
				
				Real nu = static_cast<Real>(df);
				return nu / (nu - MML::Real(2.0));
			}

		private:
			void ValidateParameters() const {
				if (!(df > 0))
					throw StatisticsError("Invalid parameters in TDistribution");
			}
		};

		/// @brief Chi-square (χ²) distribution
		///
		/// The chi-square distribution is the distribution of a sum
		/// of a sum of the squares of k independent standard normal random variables.
		///
		/// PDF: f(x; k) = (1 / (2^(k/2) * Γ(k/2))) * x^(k/2-1) * e^(-x/2) for x ≥ 0
		/// Mean: k
		///
		/// Variance: 2k
		/// Used extensively in hypothesis testing (goodness-of-fit, independence tests)
		///
		/// @param df Degrees of freedom (must be positive)
		class ChiSquareDistribution
		{
		private:
			int _df;  // degrees of freedom

		public:
			ChiSquareDistribution(int df) : _df(df)
			{
				if (df <= 0)
					throw StatisticsError("Chi-square distribution requires positive degrees of freedom");
			}

			/// @brief Draw a random variate (chi2 = Gamma(df/2, 2)) using the given URBG.
			template<class URBG> Real sample(URBG& gen) const { return SampleGamma(gen, _df / MML::Real(2.0), MML::Real(2.0)); }
			/// @brief Draw a random variate using MML's thread-local engine.
			Real sample() const { return sample(Random::Engine()); }

			int df() const { return _df; }
			Real mean() const { return static_cast<Real>(_df); }
			Real variance() const { return MML::Real(2.0) * _df; }
			Real stddev() const { return std::sqrt(variance()); }

			/// @brief Probability Density Function
			Real pdf(Real x) const
			{
				if (x < 0.0) return 0.0;
				if (x == MML::Real(0.0)) return (_df == 2) ? MML::Real(0.5) : MML::Real(0.0);

				Real k2 = _df / MML::Real(2.0);
				return std::exp((k2 - MML::Real(1.0)) * std::log(x) - x / MML::Real(2.0) - k2 * std::log(MML::Real(2.0)) - std::lgamma(k2));
			}

			/// @brief Cumulative Distribution Function (CDF)
			///
			/// P(X ≤ x) using the incomplete gamma function
			Real cdf(Real x) const
			{
				if (x <= 0.0) return 0.0;
				
				// CDF = P(k/2, x/2) where P is regularized lower incomplete gamma
				return SpecialFunctions::RegularizedGammaP(_df / MML::Real(2.0), x / MML::Real(2.0));
			}

			/// @brief Right-tail probability P(X ≥ x)
			///
			/// This is the p-value for chi-square tests
			Real rightTailPValue(Real chiSquare) const
			{
				if (chiSquare <= 0) return 1;
				return SpecialFunctions::RegularizedGammaQ(_df / MML::Real(2.0), chiSquare / 2);
			}

			/// @brief Critical value for given significance level (right-tail test)
			///
			/// Returns x such that P(X ≥ x) = alpha
			Real criticalValue(Real alpha) const
			{
				if (!(alpha > 0 && alpha < 1))
					throw StatisticsError("Alpha must be in (0, 1)");
				if (alpha > Real(0.5))
					return StatisticsDetail::PositiveQuantile(1 - alpha, [this](Real x) { return cdf(x); }, true, Real(_df));
				return StatisticsDetail::PositiveQuantile(alpha, [this](Real x) { return rightTailPValue(x); }, false, Real(_df));
			}
		};

    /// @brief F-Distribution (Fisher-Snedecor Distribution)
    ///
    /// The F-distribution is the ratio of two scaled chi-square distributions.
    /// F = (χ²₁/df₁) / (χ²₂/df₂)
    ///
    /// Used extensively in ANOVA, regression analysis, and comparing variances.
    /// PDF: f(x; d₁, d₂) = [complex formula involving beta function]
    ///
    /// Mean: d₂/(d₂-2) for d₂ > 2
    /// Variance: 2d₂²(d₁+d₂-2) / [d₁(d₂-2)²(d₂-4)] for d₂ > 4
    ///
    /// @param df1 Numerator degrees of freedom
    /// @param df2 Denominator degrees of freedom
    class FDistribution
    {
    private:
      int _df1;  // numerator df
      int _df2;  // denominator df

    public:
      FDistribution(int df1, int df2) : _df1(df1), _df2(df2)
      {
        if (df1 <= 0 || df2 <= 0)
          throw StatisticsError("F-distribution requires positive degrees of freedom");
      }

      /// @brief Draw a random variate (ratio of scaled chi-squares) using the given URBG.
      template<class URBG> Real sample(URBG& gen) const
      {
        Real u1 = SampleGamma(gen, _df1 / MML::Real(2.0), MML::Real(2.0));
        Real u2 = SampleGamma(gen, _df2 / MML::Real(2.0), MML::Real(2.0));
        return (u1 / _df1) / (u2 / _df2);
      }
      /// @brief Draw a random variate using MML's thread-local engine.
      Real sample() const { return sample(Random::Engine()); }

      int df1() const { return _df1; }
      int df2() const { return _df2; }

      Real mean() const
      {
        if (_df2 <= 2)
          throw StatisticsError("F-distribution mean undefined for df2 <= 2");
        return static_cast<Real>(_df2) / (_df2 - 2);
      }

      Real variance() const
      {
        if (_df2 <= 4)
          throw StatisticsError("F-distribution variance undefined for df2 <= 4");
        Real d1 = static_cast<Real>(_df1);
        Real d2 = static_cast<Real>(_df2);
        return (MML::Real(2.0) * d2 * d2 * (d1 + d2 - MML::Real(2.0))) /
              (d1 * (d2 - MML::Real(2.0)) * (d2 - MML::Real(2.0)) * (d2 - MML::Real(4.0)));
      }

      /// @brief Probability Density Function
      Real pdf(Real x) const
      {
        if (x < 0.0) return 0.0;
        if (x == MML::Real(0.0)) return (_df1 == 2) ? MML::Real(1.0) : MML::Real(0.0);

        Real d1 = _df1 / MML::Real(2.0);
        Real d2 = _df2 / MML::Real(2.0);

        // log PDF to avoid overflow
        Real logPdf = d1 * std::log(Real(_df1)) + d2 * std::log(Real(_df2)) +
                      (d1 - MML::Real(1.0)) * std::log(x) -
                      (d1 + d2) * std::log(_df2 + _df1 * x) +
                      std::lgamma(d1 + d2) - std::lgamma(d1) - std::lgamma(d2);

        return std::exp(logPdf);
      }

      /// @brief Cumulative Distribution Function
      ///
      /// Uses relationship with regularized incomplete beta function
      Real cdf(Real x) const
      {
        if (x <= 0.0) return 0.0;

        Real ratio = Real(_df2) / Real(_df1);
        if (x > ratio)
          return 1 - rightTailPValue(x);
        Real y = x / ratio;
        return SpecialFunctions::RegularizedBetaI(y / (1 + y), _df1 / MML::Real(2.0), _df2 / MML::Real(2.0));
      }

      /// @brief Right-tail probability P(F ≥ f)
      ///
      /// This is the p-value for F-tests (ANOVA, regression)
      Real rightTailPValue(Real f) const
      {
        if (f <= 0) return 1;
        Real ratio = Real(_df2) / Real(_df1);
        if (f < ratio)
          return 1 - cdf(f);
        Real y = ratio / f;
        return SpecialFunctions::RegularizedBetaI(y / (1 + y), _df2 / MML::Real(2.0), _df1 / MML::Real(2.0));
      }

      /// @brief Critical value for given significance level (right-tail test)
      ///
      /// Returns f such that P(F ≥ f) = alpha
      Real criticalValue(Real alpha) const
      {
        if (!(alpha > 0 && alpha < 1))
          throw StatisticsError("Alpha must be in (0, 1)");

        if (alpha > Real(0.5))
          return StatisticsDetail::PositiveQuantile(1 - alpha, [this](Real x) { return cdf(x); }, true);
        return StatisticsDetail::PositiveQuantile(alpha, [this](Real x) { return rightTailPValue(x); }, false);
      }
    };

/// @brief Cauchy distribution (Lorentz distribution)
		/// @note Heavy tails, undefined mean/variance, PDF: f(x) = 1/(π·σ·(1+((x-μ)/σ)²))
		struct CauchyDistribution
		{
			Real mu;   // Location parameter (peak location)
			Real sigma; // Scale parameter (half-width at half-maximum)

			/// @brief Construct a Cauchy distribution
			/// @param location Location parameter (default: 0)
			/// @param scale Scale parameter (default: 1, must be > 0)
			CauchyDistribution(Real location = 0.0, Real scale = 1.0) 
				: mu(location), sigma(scale)
			{
				ValidateParameters();
			}

			/// @brief Draw a random variate using the given URBG.
			template<class URBG> Real sample(URBG& gen) const { ValidateParameters(); return mu + sigma * std::tan(Constants::PI * (SampleUniform01(gen) - MML::Real(0.5))); }
			/// @brief Draw a random variate using MML's thread-local engine.
			Real sample() const { ValidateParameters(); return sample(Random::Engine()); }

			/// @brief Probability density function (PDF)
			///
			/// @param x Value at which to evaluate the PDF
			/// @return Probability density at x
			Real pdf(Real x) const
			{
				ValidateParameters();
				Real z = (x - mu) / sigma;
				return MML::Real(1.0) / (Constants::PI * sigma * (MML::Real(1.0) + z * z));
			}

			/// @brief Cumulative distribution function (CDF)
			///
			/// @param x Value at which to evaluate the CDF
			/// @return Probability that X <= x
			Real cdf(Real x) const
			{
				ValidateParameters();
				return MML::Real(0.5) + std::atan2(x - mu, sigma) / Constants::PI;
			}

			/// @brief Inverse cumulative distribution function (quantile function)
			///
			/// @param p Probability (must be in (0, 1))
			/// @return Value x such that P(X <= x) = p
			Real inverseCdf(Real p) const
			{
				ValidateParameters();
				if (!(p > 0.0 && p < 1.0))
					throw StatisticsError("Probability must be in (0, 1) in CauchyDistribution::inverseCdf");
				return mu + sigma * std::tan(Constants::PI * (p - MML::Real(0.5)));
			}

		private:
			void ValidateParameters() const {
				if (!(std::isfinite(mu) && std::isfinite(sigma) && sigma > 0))
					throw StatisticsError("Invalid parameters in CauchyDistribution");
			}
		};

		/// @brief Exponential distribution (memoryless property)
		/// @note Models time between Poisson events, PDF: f(x) = λ·exp(-λ·x), λ = rate = 1/mean
		struct ExponentialDistribution
		{
			Real lambda; // Rate parameter (inverse of mean)

			/// @brief Construct an exponential distribution
			/// @param rate Rate parameter λ (must be > 0)
			ExponentialDistribution(Real rate) : lambda(rate)
			{
				ValidateParameters();
			}

			/// @brief Draw a random variate (inverse-CDF) using the given URBG.
			template<class URBG> Real sample(URBG& gen) const { ValidateParameters(); return -std::log(MML::Real(1.0) - SampleUniform01(gen)) / lambda; }
			/// @brief Draw a random variate using MML's thread-local engine.
			Real sample() const { ValidateParameters(); return sample(Random::Engine()); }

			/// @brief Probability density function (PDF)
			/// @param x Value at which to evaluate the PDF (must be >= 0)
			/// @return Probability density at x
			Real pdf(Real x) const
			{
				ValidateParameters();
				if (x < 0.0)
					throw StatisticsError("x must be non-negative in ExponentialDistribution::pdf");
				return lambda * std::exp(-lambda * x);
			}

			/// @brief Cumulative distribution function (CDF)
			/// @param x Value at which to evaluate the CDF (must be >= 0)
			/// @return Probability that X <= x
			Real cdf(Real x) const
			{
				ValidateParameters();
				if (x < 0.0)
					throw StatisticsError("x must be non-negative in ExponentialDistribution::cdf");
				return MML::Real(1.0) - std::exp(-lambda * x);
			}

			/// @brief Inverse CDF (quantile function)
			/// @param p Probability (must be in [0, 1))
			/// @return Value x such that P(X <= x) = p
			Real inverseCdf(Real p) const
			{
				ValidateParameters();
				if (!(p >= 0.0 && p < 1.0))
					throw StatisticsError("Probability must be in [0, 1) in ExponentialDistribution::inverseCdf");
				return -std::log(MML::Real(1.0) - p) / lambda;
			}

			/// @brief Get the mean of the distribution
			/// @return Mean = 1/λ
			Real mean() const { ValidateParameters(); return MML::Real(1.0) / lambda; }

			/// @brief Get the variance of the distribution
			/// @return Variance = 1/λ²
			Real variance() const { ValidateParameters(); return MML::Real(1.0) / (lambda * lambda); }

		private:
			void ValidateParameters() const {
				if (!(std::isfinite(lambda) && lambda > 0))
					throw StatisticsError("Invalid parameters in ExponentialDistribution");
			}
		};

		/// @brief Logistic distribution (used in logistic regression, neural networks)
		/// @note Resembles normal with heavier tails, PDF: f(x) = exp(-z)/(σ·(1+exp(-z))²), z=(x-μ)/σ
		struct LogisticDistribution
		{
			Real mu;    // Location parameter (mean and median)
			Real sigma; // Scale parameter (related to variance by σ² = 3·s²/π²)

			/// @brief Construct a logistic distribution
			/// @param location Location parameter (default: 0)
			/// @param scale Scale parameter (default: 1, must be > 0)
			LogisticDistribution(Real location = 0.0, Real scale = 1.0)
				: mu(location), sigma(scale)
			{
				ValidateParameters();
			}

			/// @brief Draw a random variate (inverse-CDF) using the given URBG.
			template<class URBG> Real sample(URBG& gen) const
			{
				ValidateParameters();
				Real u = SampleUniform01(gen);
				if (u < static_cast<MML::Real>(1e-300)) u = static_cast<MML::Real>(1e-300);
				if (u > MML::Real(1.0) - static_cast<MML::Real>(1e-16)) u = MML::Real(1.0) - static_cast<MML::Real>(1e-16);
				return mu + sigma * std::log(u / (MML::Real(1.0) - u));
			}
			/// @brief Draw a random variate using MML's thread-local engine.
			Real sample() const { ValidateParameters(); return sample(Random::Engine()); }

			/// @brief Probability density function (PDF)
			///
			/// @param x Value at which to evaluate the PDF
			/// @return Probability density at x
			Real pdf(Real x) const
			{
				ValidateParameters();
				// Standard logistic PDF: f(x) = e^(-z) / (σ(1 + e^(-z))²)
				// where z = (x - μ) / σ
				Real z = (x - mu) / sigma;
				Real exp_neg_z = std::exp(-z);
				Real denom = MML::Real(1.0) + exp_neg_z;
				return exp_neg_z / (sigma * denom * denom);
			}

			/// @brief Cumulative distribution function (CDF)
			///
			/// @param x Value at which to evaluate the CDF
			/// @return Probability that X <= x
			Real cdf(Real x) const
			{
				ValidateParameters();
				// Numerically stable computation using the logistic function
				Real z = (x - mu) / sigma;
				Real exp_z = std::exp(-std::abs(z));
				
				if (z >= 0.0)
					return MML::Real(1.0) / (MML::Real(1.0) + exp_z);
				else
					return exp_z / (MML::Real(1.0) + exp_z);
			}

			/// @brief Inverse cumulative distribution function (quantile function)
			///
			/// @param p Probability (must be in (0, 1))
			/// @return Value x such that P(X <= x) = p
			Real inverseCdf(Real p) const
			{
				ValidateParameters();
				if (!(p > 0.0 && p < 1.0))
					throw StatisticsError("Probability must be in (0, 1) in LogisticDistribution::inverseCdf");
				return mu + sigma * std::log(p / (MML::Real(1.0) - p));
			}

			/// @brief Get the mean of the distribution
			/// @return Mean = μ
			Real mean() const { ValidateParameters(); return mu; }

			/// @brief Get the variance of the distribution
			/// @return Variance = (π·σ)²/3
			Real variance() const 
			{
				ValidateParameters();
				return (Constants::PI * sigma) * (Constants::PI * sigma) / MML::Real(3.0);
			}

		private:
			void ValidateParameters() const {
				if (!(std::isfinite(mu) && std::isfinite(sigma) && sigma > 0))
					throw StatisticsError("Invalid parameters in LogisticDistribution");
			}
		};

		//=========================================================================
		// UNIFORM DISTRIBUTION
		//=========================================================================

		/// @brief Continuous Uniform distribution
		/// @note Constant probability over [a, b], PDF: f(x) = 1/(b-a) for x ∈ [a,b]
		struct UniformDistribution
		{
			Real a;  // Lower bound
			Real b;  // Upper bound

			/// @brief Construct a uniform distribution
			/// @param lower Lower bound (default: 0)
			/// @param upper Upper bound (default: 1, must be > lower)
			UniformDistribution(Real lower = 0.0, Real upper = 1.0)
				: a(lower), b(upper)
			{
				ValidateParameters();
			}

			/// @brief Draw a random variate using the given URBG.
			template<class URBG> Real sample(URBG& gen) const { ValidateParameters(); return a + (b - a) * SampleUniform01(gen); }
			/// @brief Draw a random variate using MML's thread-local engine.
			Real sample() const { ValidateParameters(); return sample(Random::Engine()); }

			/// @brief Probability density function (PDF)
			Real pdf(Real x) const
			{
				ValidateParameters();
				if (x < a || x > b)
					return 0.0;
				return MML::Real(1.0) / (b - a);
			}

			/// @brief Cumulative distribution function (CDF)
			Real cdf(Real x) const
			{
				ValidateParameters();
				if (x < a) return 0.0;
				if (x > b) return 1.0;
				return (x - a) / (b - a);
			}

			/// @brief Inverse CDF (quantile function)
			Real inverseCdf(Real p) const
			{
				ValidateParameters();
				if (!(p >= 0.0 && p <= 1.0))
					throw StatisticsError("Probability must be in [0, 1] in UniformDistribution::inverseCdf");
				return a + p * (b - a);
			}

			Real mean() const { ValidateParameters(); return (a + b) / MML::Real(2.0); }
			Real variance() const {
				ValidateParameters();
				const Real halfWidth = b / 2 - a / 2;
				int exponent = 0;
				const Real mantissa = std::frexp(halfWidth, &exponent);
				return std::scalbn(mantissa * mantissa / Real(3), 2 * exponent);
			}

		private:
			void ValidateParameters() const {
				if (!(std::isfinite(a) && std::isfinite(b) && a < b))
					throw StatisticsError("Invalid parameters in UniformDistribution");
			}
		};

		//=========================================================================
		// GAMMA DISTRIBUTION
		//=========================================================================

		/// @brief Gamma distribution
		/// @note Generalizes exponential and chi-square distributions
		/// PDF: f(x; k, θ) = x^(k-1) * exp(-x/θ) / (θ^k * Γ(k))
		/// where k = shape, θ = scale
		struct GammaDistribution
		{
			Real shape;  // Shape parameter k (also called α)
			Real scale;  // Scale parameter θ (also called β)

			/// @brief Construct a gamma distribution
			/// @param k Shape parameter (must be > 0)
			/// @param theta Scale parameter (must be > 0)
			GammaDistribution(Real k, Real theta = 1.0)
				: shape(k), scale(theta)
			{
				ValidateParameters();
			}

			/// @brief Draw a random variate (Marsaglia-Tsang) using the given URBG.
			template<class URBG> Real sample(URBG& gen) const { ValidateParameters(); return SampleGamma(gen, shape, scale); }
			/// @brief Draw a random variate using MML's thread-local engine.
			Real sample() const { ValidateParameters(); return sample(Random::Engine()); }

			/// @brief Probability density function (PDF)
			Real pdf(Real x) const
			{
				ValidateParameters();
				if (x < 0.0) return 0.0;
				if (x == MML::Real(0.0)) return (shape < MML::Real(1.0)) ? std::numeric_limits<Real>::infinity()
				                                  : (shape == MML::Real(1.0) ? MML::Real(1.0)/scale : MML::Real(0.0));
				
				// Use log-space for numerical stability
				Real logPdf = (shape - MML::Real(1.0)) * std::log(x) - x / scale
				            - shape * std::log(scale) - std::lgamma(shape);
				return std::exp(logPdf);
			}

			/// @brief Cumulative distribution function (CDF)
			/// Uses regularized incomplete gamma function
			Real cdf(Real x) const
			{
				ValidateParameters();
				if (x <= 0.0) return 0.0;
				return SpecialFunctions::RegularizedGammaP(shape, x / scale);
			}

			/// @brief Inverse CDF with verified bracketing and relative tail/position accuracy.
			/// @throws ConvergenceError if no resolvable finite interior quantile is found.
			Real inverseCdf(Real p) const
			{
				ValidateParameters();
				if (std::isnan(p))
					throw StatisticsError("Probability must not be NaN in GammaDistribution::inverseCdf");
				if (p <= 0.0) return 0.0;
				if (p >= 1.0) return std::numeric_limits<Real>::infinity();
				if (p <= Real(0.5))
					return StatisticsDetail::PositiveQuantile(p, [this](Real x) { return cdf(x); }, true, scale);
				return StatisticsDetail::PositiveQuantile(1 - p,
					[this](Real x) { return SpecialFunctions::RegularizedGammaQ(shape, x / scale); }, false, scale);
			}

			Real mean() const { ValidateParameters(); return shape * scale; }
			Real variance() const { ValidateParameters(); return shape * scale * scale; }

		private:
			void ValidateParameters() const {
				StatisticsDetail::RequireGammaParameters(shape, scale);
			}
		};

		//=========================================================================
		// BETA DISTRIBUTION
		//=========================================================================

		/// @brief Beta distribution
		/// @note Defined on [0, 1], models proportions and probabilities
		/// PDF: f(x; α, β) = x^(α-1) * (1-x)^(β-1) / B(α, β)
		struct BetaDistribution
		{
			Real alpha;  // Shape parameter α
			Real beta;   // Shape parameter β

			/// @brief Construct a beta distribution
			/// @param a Shape parameter α (must be > 0)
			/// @param b Shape parameter β (must be > 0)
			BetaDistribution(Real a, Real b)
				: alpha(a), beta(b)
			{
				ValidateParameters();
			}

			/// @brief Draw a random variate (two-gamma method) using the given URBG.
			template<class URBG> Real sample(URBG& gen) const
			{
				ValidateParameters();
				Real x = SampleGamma(gen, alpha, 1.0);
				Real y = SampleGamma(gen, beta, 1.0);
				return x / (x + y);
			}
			/// @brief Draw a random variate using MML's thread-local engine.
			Real sample() const { ValidateParameters(); return sample(Random::Engine()); }

			/// @brief Probability density function (PDF)
			Real pdf(Real x) const
			{
				ValidateParameters();
				if (x < 0.0 || x > 1.0) return 0.0;
				if (x == MML::Real(0.0)) return (alpha < MML::Real(1.0)) ? std::numeric_limits<Real>::infinity()
				                                  : (alpha == MML::Real(1.0) ? beta : MML::Real(0.0));
				if (x == MML::Real(1.0)) return (beta < MML::Real(1.0)) ? std::numeric_limits<Real>::infinity()
				                                 : (beta == MML::Real(1.0) ? alpha : MML::Real(0.0));

				// Log-space for numerical stability
				Real logPdf = (alpha - MML::Real(1.0)) * std::log(x) + (beta - MML::Real(1.0)) * std::log(MML::Real(1.0) - x)
				            - std::lgamma(alpha) - std::lgamma(beta) + std::lgamma(alpha + beta);
				return std::exp(logPdf);
			}

			/// @brief Cumulative distribution function (CDF)
			/// Uses regularized incomplete beta function
			Real cdf(Real x) const
			{
				ValidateParameters();
				if (x <= 0.0) return 0.0;
				if (x >= 1.0) return 1.0;
				return SpecialFunctions::RegularizedBetaI(x, alpha, beta);
			}

			/// @brief Inverse CDF with verified bracketing and relative tail/position accuracy.
			/// @throws ConvergenceError if no resolvable interior quantile is found.
			Real inverseCdf(Real p) const
			{
				ValidateParameters();
				if (std::isnan(p))
					throw StatisticsError("Probability must not be NaN in BetaDistribution::inverseCdf");
				if (p <= 0.0) return 0.0;
				if (p >= 1.0) return 1.0;
				if (p <= Real(0.5))
					return StatisticsDetail::PositiveQuantile(p, [this](Real x) { return cdf(x); }, true, 1, false);
				return StatisticsDetail::PositiveQuantile(1 - p,
					[this](Real x) { return SpecialFunctions::RegularizedBetaI(1 - x, beta, alpha); }, false, 1, false);
			}

			Real mean() const { ValidateParameters(); return alpha / (alpha + beta); }
			Real variance() const 
			{
				ValidateParameters();
				Real ab = alpha + beta;
				return (alpha * beta) / (ab * ab * (ab + MML::Real(1.0)));
			}

		private:
			void ValidateParameters() const {
				if (!(std::isfinite(alpha) && alpha > 0 && std::isfinite(beta) && beta > 0))
					throw StatisticsError("Invalid parameters in BetaDistribution");
			}
		};

		//=========================================================================
		// WEIBULL DISTRIBUTION
		//=========================================================================

		/// @brief Weibull distribution
		/// @note Used in reliability analysis, failure modeling
		/// PDF: f(x; k, λ) = (k/λ) * (x/λ)^(k-1) * exp(-(x/λ)^k)
		struct WeibullDistribution
		{
			Real shape;  // Shape parameter k
			Real scale;  // Scale parameter λ

			/// @brief Construct a Weibull distribution
			/// @param k Shape parameter (must be > 0)
			/// @param lambda Scale parameter (must be > 0)
			WeibullDistribution(Real k, Real lambda = 1.0)
				: shape(k), scale(lambda)
			{
				ValidateParameters();
			}

			/// @brief Draw a random variate (inverse-CDF) using the given URBG.
			template<class URBG> Real sample(URBG& gen) const { ValidateParameters(); return scale * std::pow(-std::log(MML::Real(1.0) - SampleUniform01(gen)), MML::Real(1.0) / shape); }
			/// @brief Draw a random variate using MML's thread-local engine.
			Real sample() const { ValidateParameters(); return sample(Random::Engine()); }

			/// @brief Probability density function (PDF)
			Real pdf(Real x) const
			{
				ValidateParameters();
				if (x < 0.0) return 0.0;
				if (x == MML::Real(0.0)) return (shape < MML::Real(1.0)) ? std::numeric_limits<Real>::infinity()
				                                  : (shape == MML::Real(1.0) ? MML::Real(1.0)/scale : MML::Real(0.0));
				Real z = x / scale;
				return (shape / scale) * std::pow(z, shape - MML::Real(1.0)) * std::exp(-std::pow(z, shape));
			}

			/// @brief Cumulative distribution function (CDF)
			Real cdf(Real x) const
			{
				ValidateParameters();
				if (x <= 0.0) return 0.0;
				return MML::Real(1.0) - std::exp(-std::pow(x / scale, shape));
			}

			/// @brief Inverse CDF (quantile function) - closed form
			Real inverseCdf(Real p) const
			{
				ValidateParameters();
				if (std::isnan(p))
					throw StatisticsError("Probability must not be NaN in WeibullDistribution::inverseCdf");
				if (p <= 0.0) return 0.0;
				if (p >= 1.0) return std::numeric_limits<Real>::infinity();
				return scale * std::pow(-std::log(MML::Real(1.0) - p), MML::Real(1.0) / shape);
			}

			Real mean() const { ValidateParameters(); return scale * std::tgamma(MML::Real(1.0) + MML::Real(1.0) / shape); }
			Real variance() const 
			{
				ValidateParameters();
				Real g1 = std::tgamma(MML::Real(1.0) + MML::Real(1.0) / shape);
				Real g2 = std::tgamma(MML::Real(1.0) + MML::Real(2.0) / shape);
				return scale * scale * (g2 - g1 * g1);
			}

		private:
			void ValidateParameters() const {
				if (!(std::isfinite(shape) && shape > 0 && std::isfinite(scale) && scale > 0))
					throw StatisticsError("Invalid parameters in WeibullDistribution");
			}
		};

		//=========================================================================
		// PARETO DISTRIBUTION
		//=========================================================================

		/// @brief Pareto distribution (Type I)
		/// @note Power-law distribution, models wealth, file sizes, etc.
		/// PDF: f(x; α, x_m) = α * x_m^α / x^(α+1) for x >= x_m
		struct ParetoDistribution
		{
			Real alpha;  // Shape parameter (tail index)
			Real xm;     // Scale parameter (minimum value)

			/// @brief Construct a Pareto distribution
			/// @param shape Shape parameter α (must be > 0)
			/// @param scale Minimum value x_m (must be > 0)
			ParetoDistribution(Real shape, Real scale = 1.0)
				: alpha(shape), xm(scale)
			{
				ValidateParameters();
			}

			/// @brief Draw a random variate (inverse-CDF) using the given URBG.
			template<class URBG> Real sample(URBG& gen) const { ValidateParameters(); return xm / std::pow(MML::Real(1.0) - SampleUniform01(gen), MML::Real(1.0) / alpha); }
			/// @brief Draw a random variate using MML's thread-local engine.
			Real sample() const { ValidateParameters(); return sample(Random::Engine()); }

			/// @brief Probability density function (PDF)
			Real pdf(Real x) const
			{
				ValidateParameters();
				if (x < xm) return 0.0;
				return alpha * std::pow(xm, alpha) / std::pow(x, alpha + MML::Real(1.0));
			}

			/// @brief Cumulative distribution function (CDF)
			Real cdf(Real x) const
			{
				ValidateParameters();
				if (x < xm) return 0.0;
				return MML::Real(1.0) - std::pow(xm / x, alpha);
			}

			/// @brief Inverse CDF (quantile function) - closed form
			Real inverseCdf(Real p) const
			{
				ValidateParameters();
				if (std::isnan(p))
					throw StatisticsError("Probability must not be NaN in ParetoDistribution::inverseCdf");
				if (p <= 0.0) return xm;
				if (p >= 1.0) return std::numeric_limits<Real>::infinity();
				return xm / std::pow(MML::Real(1.0) - p, MML::Real(1.0) / alpha);
			}

			/// @brief Mean (only defined for α > 1)
			Real mean() const
			{
				ValidateParameters();
				if (alpha <= 1.0)
					throw StatisticsError("Mean undefined for alpha <= 1 in ParetoDistribution");
				return alpha * xm / (alpha - MML::Real(1.0));
			}

			/// @brief Variance (only defined for α > 2)
			Real variance() const
			{
				ValidateParameters();
				if (alpha <= 2.0)
					throw StatisticsError("Variance undefined for alpha <= 2 in ParetoDistribution");
				return (xm * xm * alpha) / ((alpha - MML::Real(1.0)) * (alpha - MML::Real(1.0)) * (alpha - MML::Real(2.0)));
			}

		private:
			void ValidateParameters() const {
				if (!(std::isfinite(alpha) && alpha > 0 && std::isfinite(xm) && xm > 0))
					throw StatisticsError("Invalid parameters in ParetoDistribution");
			}
		};

		//=========================================================================
		// LOG-NORMAL DISTRIBUTION
		//=========================================================================

		/// @brief Log-Normal distribution
		/// @note X is log-normal if ln(X) is normal; models multiplicative processes
		/// PDF: f(x; μ, σ) = (1/(x·σ·√(2π))) * exp(-((ln(x)-μ)²)/(2σ²))
		struct LogNormalDistribution
		{
			Real mu;     // Mean of ln(X)
			Real sigma;  // Std dev of ln(X)

			/// @brief Construct a log-normal distribution
			/// @param logMean Mean of the log (μ)
			/// @param logStdDev Standard deviation of the log (σ, must be > 0)
			LogNormalDistribution(Real logMean = 0.0, Real logStdDev = 1.0)
				: mu(logMean), sigma(logStdDev)
			{
				ValidateParameters();
			}

			/// @brief Draw a random variate (exp of a normal draw) using the given URBG.
			template<class URBG> Real sample(URBG& gen) const { ValidateParameters(); return std::exp(mu + sigma * SampleStandardNormal(gen)); }
			/// @brief Draw a random variate using MML's thread-local engine.
			Real sample() const { ValidateParameters(); return sample(Random::Engine()); }

			/// @brief Probability density function (PDF)
			Real pdf(Real x) const
			{
				ValidateParameters();
				if (x <= 0.0) return 0.0;
				Real z = (std::log(x) - mu) / sigma;
				return std::exp(-MML::Real(0.5) * z * z) / (x * sigma * std::sqrt(MML::Real(2.0) * Constants::PI));
			}

			/// @brief Cumulative distribution function (CDF)
			Real cdf(Real x) const
			{
				ValidateParameters();
				if (x <= 0.0) return 0.0;
				return StatisticsDetail::NormalCdf(StatisticsDetail::StandardizeNormal(std::log(x), mu, sigma));
			}

			/// @brief Inverse CDF (quantile function)
			Real inverseCdf(Real p) const
			{
				ValidateParameters();
				if (std::isnan(p))
					throw StatisticsError("Probability must not be NaN in LogNormalDistribution::inverseCdf");
				if (p <= 0.0) return 0.0;
				if (p >= 1.0) return std::numeric_limits<Real>::infinity();
				// Use inverse normal for ln(x)
				Real z = InverseStandardNormalCdf(p);
				return std::exp(mu + sigma * z);
			}

			Real mean() const { ValidateParameters(); return std::exp(mu + MML::Real(0.5) * sigma * sigma); }
			Real variance() const 
			{
				ValidateParameters();
				Real expSigma2 = std::exp(sigma * sigma);
				return std::exp(MML::Real(2.0) * mu + sigma * sigma) * (expSigma2 - MML::Real(1.0));
			}
			Real median() const { ValidateParameters(); return std::exp(mu); }

		private:
			void ValidateParameters() const {
				if (!(std::isfinite(mu) && std::isfinite(sigma) && sigma > 0))
					throw StatisticsError("Invalid parameters in LogNormalDistribution");
			}
		};

  }
}

#endif // MML_DISTRIBUTIONS_H

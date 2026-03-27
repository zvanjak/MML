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

#include "algorithms/Statistics.h"

namespace MML
{
	namespace Statistics
	{
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
				if (sigma <= 0.0)
					throw StatisticsError("Standard deviation must be positive in NormalDistribution");
			}

			/// @brief Probability density function (PDF)
			///
			/// @param x Value at which to evaluate the PDF
			/// @return Probability density at x
			Real pdf(Real x) const
			{
				Real z = (x - mu) / sigma;
				return std::exp(-0.5 * z * z) / (sigma * std::sqrt(2.0 * Constants::PI));
			}

			/// @brief Cumulative distribution function (CDF)
			///
			/// Uses the error function: Φ(x) = 0.5 * (1 + erf(x/√2))
			/// @param x Value at which to evaluate the CDF
			///
			/// @return Probability that X <= x
			Real cdf(Real x) const
			{
				Real z = (x - mu) / (sigma * std::sqrt(2.0));
				return 0.5 * (1.0 + std::erf(z));
			}

			/// @brief Inverse cumulative distribution function (quantile/probit function)
			///
			/// Uses Abramowitz and Stegun approximation for the inverse error function.
			/// @param p Probability (must be in (0, 1))
			///
			/// @return Value x such that P(X <= x) = p
			Real inverseCdf(Real p) const
			{
				if (p <= 0.0 || p >= 1.0)
					throw StatisticsError("Probability must be in (0, 1) in NormalDistribution::inverseCdf");

				// Use rational approximation for inverse normal CDF
				// Based on Abramowitz and Stegun approximation 26.2.23
				Real z = inverseStandardNormalCdf(p);
				return mu + sigma * z;
			}

			/// @brief Z-score: standardize a value
			///
			/// @param x Raw value
			/// @return Standardized value (z-score)
			Real zScore(Real x) const
			{
				return (x - mu) / sigma;
			}

			/// @brief Convert z-score back to raw value
			///
			/// @param z Standardized value
			/// @return Raw value
			Real fromZScore(Real z) const
			{
				return mu + sigma * z;
			}

			/// @brief Mean of the distribution
			Real mean() const { return mu; }

			/// @brief Variance of the distribution
			Real variance() const { return sigma * sigma; }

			/// @brief Standard deviation of the distribution
			Real stddev() const { return sigma; }

		private:
			/// @brief Inverse of the standard normal CDF (probit function)
			/// Uses rational approximation from Abramowitz and Stegun (26.2.23)
			/// with refinement for tails. Accurate to ~1e-9.
			static Real inverseStandardNormalCdf(Real p)
			{
				if (p <= 0.0 || p >= 1.0)
					throw StatisticsError("Probability must be in (0, 1) in inverseStandardNormalCdf");

				// Coefficients for rational approximation
				const Real a[] = {
					-3.969683028665376e+01,
					 2.209460984245205e+02,
					-2.759285104469687e+02,
					 1.383577518672690e+02,
					-3.066479806614716e+01,
					 2.506628277459239e+00
				};
				const Real b[] = {
					-5.447609879822406e+01,
					 1.615858368580409e+02,
					-1.556989798598866e+02,
					 6.680131188771972e+01,
					-1.328068155288572e+01
				};
				const Real c[] = {
					-7.784894002430293e-03,
					-3.223964580411365e-01,
					-2.400758277161838e+00,
					-2.549732539343734e+00,
					 4.374664141464968e+00,
					 2.938163982698783e+00
				};
				const Real d[] = {
					 7.784695709041462e-03,
					 3.224671290700398e-01,
					 2.445134137142996e+00,
					 3.754408661907416e+00
				};

				const Real pLow = 0.02425;
				const Real pHigh = 1.0 - pLow;

				Real q, r;

				if (p < pLow) {
					// Lower tail
					q = std::sqrt(-2.0 * std::log(p));
					return (((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) /
					        ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1.0);
				} else if (p <= pHigh) {
					// Central region
					q = p - 0.5;
					r = q * q;
					return (((((a[0]*r + a[1])*r + a[2])*r + a[3])*r + a[4])*r + a[5]) * q /
					       (((((b[0]*r + b[1])*r + b[2])*r + b[3])*r + b[4])*r + 1.0);
				} else {
					// Upper tail
					q = std::sqrt(-2.0 * std::log(1.0 - p));
					return -(((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) /
					         ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1.0);
				}
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
				return std::exp(-0.5 * z * z) / std::sqrt(2.0 * Constants::PI);
			}

			/// @brief CDF of standard normal (Φ function)
			static Real cdf(Real z)
			{
				return 0.5 * (1.0 + std::erf(z / std::sqrt(2.0)));
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
				return 1.0 - cdf(z);
			}

			/// @brief Two-tailed p-value: P(|Z| > |z|)
			static Real twoTailedPValue(Real z)
			{
				return 2.0 * sf(std::abs(z));
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
				if (df < 1)
					throw StatisticsError("Degrees of freedom must be at least 1 in TDistribution");
			}

			/// @brief Probability density function (PDF)
			///
			/// @param x Value at which to evaluate the PDF
			/// @return Probability density at x
			Real pdf(Real x) const
			{
				Real nu = static_cast<Real>(df);
				
				// Log-space computation for numerical stability
				Real logNumerator = std::lgamma((nu + 1.0) / 2.0);
				Real logDenominator = 0.5 * std::log(nu * Constants::PI) + std::lgamma(nu / 2.0);
				Real logBase = -(nu + 1.0) / 2.0 * std::log(1.0 + x * x / nu);
				
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
				if (df == 1) {
					// Special case: df=1 is Cauchy distribution centered at 0
					return 0.5 + std::atan(x) / Constants::PI;
				}

				Real nu = static_cast<Real>(df);
				
				if (x == 0.0)
					return 0.5;
				
				// Standard formula for t-distribution CDF:
				// For t > 0: P(T <= t) = 1 - 0.5 * I_x(df/2, 0.5)
				// where x = df/(df + t²)
				// For t < 0: Use symmetry P(T <= -|t|) = 1 - P(T <= |t|)
				
				Real absX = std::abs(x);
				Real x2 = absX * absX;
				Real xBeta = nu / (nu + x2);
				Real betaReg = incompleteBetaRegularized(xBeta, nu / 2.0, 0.5);
				
				Real result = 1.0 - 0.5 * betaReg;
				
				// If x < 0, use symmetry
				if (x < 0.0)
					result = 1.0 - result;
				
				return result;
			}

			/// @brief Inverse cumulative distribution function (quantile function)
			///
			/// Uses iterative root finding for the inverse.
			/// @param p Probability (must be in (0, 1))
			///
			/// @return Value t such that P(T <= t) = p
			Real inverseCdf(Real p) const
			{
				if (p <= 0.0 || p >= 1.0)
					throw StatisticsError("Probability must be in (0, 1) in TDistribution::inverseCdf");

				if (p == 0.5)
					return 0.0;

				if (df == 1) {
					// Special case: df=1 is Cauchy
					return std::tan(Constants::PI * (p - 0.5));
				}

				// Use Newton-Raphson iteration
				// Start with normal approximation
				Real t = StandardNormal::inverseCdf(p);
				
				// Refine with Newton-Raphson
				for (int iter = 0; iter < 10; iter++) {
					Real f = cdf(t) - p;
					Real fprime = pdf(t);
					
					if (std::abs(fprime) < 1e-14)
						break;
					
					Real delta = f / fprime;
					t -= delta;
					
					if (std::abs(delta) < Precision::DefaultToleranceStrict)
						break;
				}
				
				return t;
			}

			/// @brief Survival function: P(T > t)
			Real sf(Real x) const
			{
				return 1.0 - cdf(x);
			}

			/// @brief Two-tailed p-value: P(|T| > |t|)
			///
			/// This is the standard p-value used in t-tests.
			Real twoTailedPValue(Real t) const
			{
				return 2.0 * sf(std::abs(t));
			}

			/// @brief Critical value for given significance level (two-tailed)
			///
			/// Returns t* such that P(|T| > t*) = alpha
			/// @param alpha Significance level (e.g., 0.05 for 95% confidence)
			///
			/// @return Critical value t*
			Real criticalValue(Real alpha) const
			{
				if (alpha <= 0.0 || alpha >= 1.0)
					throw StatisticsError("Alpha must be in (0, 1) in TDistribution::criticalValue");
				
				// Two-tailed: find t* where P(T > t*) = alpha/2
				return inverseCdf(1.0 - alpha / 2.0);
			}

			/// @brief Mean of the distribution (if it exists)
			///
			/// Mean exists only for df >= 2, and equals 0.
			Real mean() const
			{
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
				if (df <= 2)
					throw StatisticsError("Variance does not exist for df <= 2 in TDistribution");
				
				Real nu = static_cast<Real>(df);
				return nu / (nu - 2.0);
			}

		private:
			/// @brief Regularized incomplete beta function I_x(a, b)
			/// Uses continued fraction expansion (Lentz's method).
			/// Based on Numerical Recipes formula 6.4.5
			static Real incompleteBetaRegularized(Real x, Real a, Real b)
			{
				if (x < 0.0 || x > 1.0)
					return 0.0;
				if (x == 0.0)
					return 0.0;
				if (x == 1.0)
					return 1.0;

				// Use symmetry relation if x > (a+1)/(a+b+2)
				bool useSymmetry = (x > (a + 1.0) / (a + b + 2.0));
				if (useSymmetry) {
					return 1.0 - incompleteBetaRegularized(1.0 - x, b, a);
				}

				// Compute log of beta function B(a, b) = Γ(a)Γ(b)/Γ(a+b)
				Real logBeta = std::lgamma(a) + std::lgamma(b) - std::lgamma(a + b);
				
				// Front factor: exp(a*ln(x) + b*ln(1-x) - ln(B(a,b))) / a
				Real front = std::exp(a * std::log(x) + b * std::log(1.0 - x) - logBeta) / a;

				// Modified Lentz's method for continued fraction
				// Formula from Numerical Recipes (Press et al.)
				const Real fpmin = 1.0e-30;
				const Real eps = 1.0e-15;
				const int maxIter = 200;
				
				Real qab = a + b;
				Real qap = a + 1.0;
				Real qam = a - 1.0;
				Real c = 1.0;
				Real d = 1.0 - qab * x / qap;
				if (std::abs(d) < fpmin) d = fpmin;
				d = 1.0 / d;
				Real h = d;
				
				for (int m = 1; m <= maxIter; m++) {
					int m2 = 2 * m;
					Real aa = m * (b - m) * x / ((qam + m2) * (a + m2));
					d = 1.0 + aa * d;
					if (std::abs(d) < fpmin) d = fpmin;
					c = 1.0 + aa / c;
					if (std::abs(c) < fpmin) c = fpmin;
					d = 1.0 / d;
					h *= d * c;
					
					aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
					d = 1.0 + aa * d;
					if (std::abs(d) < fpmin) d = fpmin;
					c = 1.0 + aa / c;
					if (std::abs(c) < fpmin) c = fpmin;
					d = 1.0 / d;
					Real del = d * c;
					h *= del;
					
					if (std::abs(del - 1.0) < eps)
						break;
				}

				return front * h;
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

			int df() const { return _df; }
			Real mean() const { return static_cast<Real>(_df); }
			Real variance() const { return 2.0 * _df; }
			Real stddev() const { return std::sqrt(variance()); }

			/// @brief Probability Density Function
			Real pdf(Real x) const
			{
				if (x < 0.0) return 0.0;
				if (x == 0.0) return (_df == 2) ? 0.5 : 0.0;

				Real k2 = _df / 2.0;
				return std::exp((k2 - 1.0) * std::log(x) - x / 2.0 - k2 * std::log(2.0) - std::lgamma(k2));
			}

			/// @brief Cumulative Distribution Function (CDF)
			///
			/// P(X ≤ x) using the incomplete gamma function
			Real cdf(Real x) const
			{
				if (x <= 0.0) return 0.0;
				
				// CDF = P(k/2, x/2) where P is regularized lower incomplete gamma
				return incompleteGammaP(_df / 2.0, x / 2.0);
			}

			/// @brief Right-tail probability P(X ≥ x)
			///
			/// This is the p-value for chi-square tests
			Real rightTailPValue(Real chiSquare) const
			{
				return 1.0 - cdf(chiSquare);
			}

			/// @brief Critical value for given significance level (right-tail test)
			///
			/// Returns x such that P(X ≥ x) = alpha
			Real criticalValue(Real alpha) const
			{
				if (alpha <= 0.0 || alpha >= 1.0)
					throw StatisticsError("Alpha must be in (0, 1)");

				// Find x where CDF(x) = 1 - alpha using binary search
				Real low = 0.0;
				Real high = mean() + 10.0 * stddev();  // Start with reasonable upper bound
				Real tolerance = 1e-9;
				int maxIterations = 100;

				for (int iter = 0; iter < maxIterations; iter++) {
					Real mid = (low + high) / 2.0;
					Real cdf_mid = cdf(mid);
					Real target = 1.0 - alpha;

					if (std::abs(cdf_mid - target) < tolerance)
						return mid;

					if (cdf_mid < target)
						low = mid;
					else
						high = mid;
				}

				return (low + high) / 2.0;
			}

		private:
			/// @brief Regularized lower incomplete gamma function P(a,x)
			/// P(a,x) = γ(a,x) / Γ(a)
			Real incompleteGammaP(Real a, Real x) const
			{
				if (x < 0.0 || a <= 0.0)
					return 0.0;

				if (x == 0.0)
					return 0.0;

				// Use series expansion for x < a+1
				if (x < a + 1.0) {
					return gammaSeriesExpansion(a, x);
				}
				// Use continued fraction for x >= a+1
				else {
					return 1.0 - gammaContinuedFraction(a, x);
				}
			}

			/// @brief Series expansion for incomplete gamma
			Real gammaSeriesExpansion(Real a, Real x) const
			{
				Real sum = 1.0 / a;
				Real term = 1.0 / a;
				Real tolerance = 1e-12;

				for (int n = 1; n <= 1000; n++) {
					term *= x / (a + n);
					sum += term;

					if (std::abs(term) < std::abs(sum) * tolerance)
						break;
				}

				return sum * std::exp(-x + a * std::log(x) - std::lgamma(a));
			}

			/// @brief Continued fraction for incomplete gamma
			Real gammaContinuedFraction(Real a, Real x) const
			{
				Real tolerance = 1e-12;
				Real fpmin = 1e-30;

				Real b = x + 1.0 - a;
				Real c = 1.0 / fpmin;
				Real d = 1.0 / b;
				Real h = d;

				for (int i = 1; i <= 1000; i++) {
					Real an = -i * (i - a);
					b += 2.0;
					d = an * d + b;
					if (std::abs(d) < fpmin) d = fpmin;
					c = b + an / c;
					if (std::abs(c) < fpmin) c = fpmin;
					d = 1.0 / d;
					Real delta = d * c;
					h *= delta;

					if (std::abs(delta - 1.0) < tolerance)
						break;
				}

				return h * std::exp(-x + a * std::log(x) - std::lgamma(a));
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
        return (2.0 * d2 * d2 * (d1 + d2 - 2.0)) /
              (d1 * (d2 - 2.0) * (d2 - 2.0) * (d2 - 4.0));
      }

      /// @brief Probability Density Function
      Real pdf(Real x) const
      {
        if (x < 0.0) return 0.0;
        if (x == 0.0) return (_df1 == 2) ? 1.0 : 0.0;

        Real d1 = _df1 / 2.0;
        Real d2 = _df2 / 2.0;

        // log PDF to avoid overflow
        Real logPdf = d1 * std::log(_df1) + d2 * std::log(_df2) +
                      (d1 - 1.0) * std::log(x) -
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

        // F-CDF = I_y(df1/2, df2/2) where y = (df1*x)/(df1*x + df2)
        // I is regularized incomplete beta function
        Real y = (_df1 * x) / (_df1 * x + _df2);
        return incompleteBeta(_df1 / 2.0, _df2 / 2.0, y);
      }

      /// @brief Right-tail probability P(F ≥ f)
      ///
      /// This is the p-value for F-tests (ANOVA, regression)
      Real rightTailPValue(Real f) const
      {
        return 1.0 - cdf(f);
      }

      /// @brief Critical value for given significance level (right-tail test)
      ///
      /// Returns f such that P(F ≥ f) = alpha
      Real criticalValue(Real alpha) const
      {
        if (alpha <= 0.0 || alpha >= 1.0)
          throw StatisticsError("Alpha must be in (0, 1)");

        // Binary search for critical value
        Real low = 0.0;
        Real high = 100.0;  // Start with reasonable upper bound
        Real tolerance = 1e-9;
        int maxIterations = 100;

        for (int iter = 0; iter < maxIterations; iter++) {
          Real mid = (low + high) / 2.0;
          Real pValue = rightTailPValue(mid);

          if (std::abs(pValue - alpha) < tolerance)
            return mid;

          if (pValue > alpha)  // Need larger F
            low = mid;
          else  // Need smaller F
            high = mid;
        }

        return (low + high) / 2.0;
      }

    private:
      /// @brief Regularized incomplete beta function I_x(a,b)
      Real incompleteBeta(Real a, Real b, Real x) const
      {
        if (x < 0.0 || x > 1.0)
          return 0.0;
        if (x == 0.0) return 0.0;
        if (x == 1.0) return 1.0;

        // Use symmetry relation if needed
        if (x > (a + 1.0) / (a + b + 2.0)) {
          return 1.0 - incompleteBeta(b, a, 1.0 - x);
        }

        // Continued fraction evaluation
        Real bt = std::exp(std::lgamma(a + b) - std::lgamma(a) - std::lgamma(b) +
                          a * std::log(x) + b * std::log(1.0 - x));

        if (x < (a + 1.0) / (a + b + 2.0))
          return bt * betaContinuedFraction(a, b, x) / a;
        else
          return 1.0 - bt * betaContinuedFraction(b, a, 1.0 - x) / b;
      }

      /// @brief Continued fraction for incomplete beta
      Real betaContinuedFraction(Real a, Real b, Real x) const
      {
        Real tolerance = 1e-12;
        Real fpmin = 1e-30;
        Real qab = a + b;
        Real qap = a + 1.0;
        Real qam = a - 1.0;
        Real c = 1.0;
        Real d = 1.0 - qab * x / qap;

        if (std::abs(d) < fpmin) d = fpmin;
        d = 1.0 / d;
        Real h = d;

        for (int m = 1; m <= 1000; m++) {
          int m2 = 2 * m;
          Real aa = m * (b - m) * x / ((qam + m2) * (a + m2));
          d = 1.0 + aa * d;
          if (std::abs(d) < fpmin) d = fpmin;
          c = 1.0 + aa / c;
          if (std::abs(c) < fpmin) c = fpmin;
          d = 1.0 / d;
          h *= d * c;

          aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
          d = 1.0 + aa * d;
          if (std::abs(d) < fpmin) d = fpmin;
          c = 1.0 + aa / c;
          if (std::abs(c) < fpmin) c = fpmin;
          d = 1.0 / d;
          Real delta = d * c;
          h *= delta;

          if (std::abs(delta - 1.0) < tolerance)
            break;
        }

        return h;
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
				if (sigma <= 0.0)
					throw StatisticsError("Scale parameter must be positive in CauchyDistribution");
			}

			/// @brief Probability density function (PDF)
			///
			/// @param x Value at which to evaluate the PDF
			/// @return Probability density at x
			Real pdf(Real x) const
			{
				Real z = (x - mu) / sigma;
				return 1.0 / (Constants::PI * sigma * (1.0 + z * z));
			}

			/// @brief Cumulative distribution function (CDF)
			///
			/// @param x Value at which to evaluate the CDF
			/// @return Probability that X <= x
			Real cdf(Real x) const
			{
				return 0.5 + std::atan2(x - mu, sigma) / Constants::PI;
			}

			/// @brief Inverse cumulative distribution function (quantile function)
			///
			/// @param p Probability (must be in (0, 1))
			/// @return Value x such that P(X <= x) = p
			Real inverseCdf(Real p) const
			{
				if (p <= 0.0 || p >= 1.0)
					throw StatisticsError("Probability must be in (0, 1) in CauchyDistribution::inverseCdf");
				return mu + sigma * std::tan(Constants::PI * (p - 0.5));
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
				if (lambda <= 0.0)
					throw StatisticsError("Rate parameter must be positive in ExponentialDistribution");
			}

			/// @brief Probability density function (PDF)
			/// @param x Value at which to evaluate the PDF (must be >= 0)
			/// @return Probability density at x
			Real pdf(Real x) const
			{
				if (x < 0.0)
					throw StatisticsError("x must be non-negative in ExponentialDistribution::pdf");
				return lambda * std::exp(-lambda * x);
			}

			/// @brief Cumulative distribution function (CDF)
			/// @param x Value at which to evaluate the CDF (must be >= 0)
			/// @return Probability that X <= x
			Real cdf(Real x) const
			{
				if (x < 0.0)
					throw StatisticsError("x must be non-negative in ExponentialDistribution::cdf");
				return 1.0 - std::exp(-lambda * x);
			}

			/// @brief Inverse CDF (quantile function)
			/// @param p Probability (must be in [0, 1))
			/// @return Value x such that P(X <= x) = p
			Real inverseCdf(Real p) const
			{
				if (p < 0.0 || p >= 1.0)
					throw StatisticsError("Probability must be in [0, 1) in ExponentialDistribution::inverseCdf");
				return -std::log(1.0 - p) / lambda;
			}

			/// @brief Get the mean of the distribution
			/// @return Mean = 1/λ
			Real mean() const { return 1.0 / lambda; }

			/// @brief Get the variance of the distribution
			/// @return Variance = 1/λ²
			Real variance() const { return 1.0 / (lambda * lambda); }
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
				if (sigma <= 0.0)
					throw StatisticsError("Scale parameter must be positive in LogisticDistribution");
			}

			/// @brief Probability density function (PDF)
			///
			/// @param x Value at which to evaluate the PDF
			/// @return Probability density at x
			Real pdf(Real x) const
			{
				// Standard logistic PDF: f(x) = e^(-z) / (σ(1 + e^(-z))²)
				// where z = (x - μ) / σ
				Real z = (x - mu) / sigma;
				Real exp_neg_z = std::exp(-z);
				Real denom = 1.0 + exp_neg_z;
				return exp_neg_z / (sigma * denom * denom);
			}

			/// @brief Cumulative distribution function (CDF)
			///
			/// @param x Value at which to evaluate the CDF
			/// @return Probability that X <= x
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

			/// @brief Inverse cumulative distribution function (quantile function)
			///
			/// @param p Probability (must be in (0, 1))
			/// @return Value x such that P(X <= x) = p
			Real inverseCdf(Real p) const
			{
				if (p <= 0.0 || p >= 1.0)
					throw StatisticsError("Probability must be in (0, 1) in LogisticDistribution::inverseCdf");
				return mu + sigma * std::log(p / (1.0 - p));
			}

			/// @brief Get the mean of the distribution
			/// @return Mean = μ
			Real mean() const { return mu; }

			/// @brief Get the variance of the distribution
			/// @return Variance = (π·σ)²/3
			Real variance() const 
			{ 
				return (Constants::PI * sigma) * (Constants::PI * sigma) / 3.0;
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
				if (b <= a)
					throw StatisticsError("Upper bound must be greater than lower bound in UniformDistribution");
			}

			/// @brief Probability density function (PDF)
			Real pdf(Real x) const
			{
				if (x < a || x > b)
					return 0.0;
				return 1.0 / (b - a);
			}

			/// @brief Cumulative distribution function (CDF)
			Real cdf(Real x) const
			{
				if (x < a) return 0.0;
				if (x > b) return 1.0;
				return (x - a) / (b - a);
			}

			/// @brief Inverse CDF (quantile function)
			Real inverseCdf(Real p) const
			{
				if (p < 0.0 || p > 1.0)
					throw StatisticsError("Probability must be in [0, 1] in UniformDistribution::inverseCdf");
				return a + p * (b - a);
			}

			Real mean() const { return (a + b) / 2.0; }
			Real variance() const { return (b - a) * (b - a) / 12.0; }
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
				if (shape <= 0.0)
					throw StatisticsError("Shape parameter must be positive in GammaDistribution");
				if (scale <= 0.0)
					throw StatisticsError("Scale parameter must be positive in GammaDistribution");
			}

			/// @brief Probability density function (PDF)
			Real pdf(Real x) const
			{
				if (x < 0.0) return 0.0;
				if (x == 0.0) return (shape < 1.0) ? std::numeric_limits<Real>::infinity() 
				                                  : (shape == 1.0 ? 1.0/scale : 0.0);
				
				// Use log-space for numerical stability
				Real logPdf = (shape - 1.0) * std::log(x) - x / scale 
				            - shape * std::log(scale) - std::lgamma(shape);
				return std::exp(logPdf);
			}

			/// @brief Cumulative distribution function (CDF)
			/// Uses regularized incomplete gamma function
			Real cdf(Real x) const
			{
				if (x <= 0.0) return 0.0;
				return incompleteGammaP(shape, x / scale);
			}

			/// @brief Inverse CDF (quantile function) - uses Newton-Raphson
			Real inverseCdf(Real p) const
			{
				if (p <= 0.0) return 0.0;
				if (p >= 1.0) return std::numeric_limits<Real>::infinity();

				// Initial guess using Wilson-Hilferty approximation for large shape
				Real x;
				if (shape >= 1.0) {
					Real z = inverseStandardNormalCdf(p);
					Real h = 1.0 / (9.0 * shape);
					x = shape * scale * std::pow(1.0 - h + z * std::sqrt(h), 3);
					if (x <= 0.0) x = 0.5 * shape * scale;
				} else {
					x = 0.5 * shape * scale;
				}

				// Newton-Raphson refinement
				for (int iter = 0; iter < 50; iter++) {
					Real f = cdf(x) - p;
					Real fp = pdf(x);
					if (std::abs(fp) < 1e-30) break;
					Real dx = f / fp;
					x -= dx;
					if (x <= 0.0) x = 1e-10;
					if (std::abs(dx) < 1e-12 * x) break;
				}
				return x;
			}

			Real mean() const { return shape * scale; }
			Real variance() const { return shape * scale * scale; }

		private:
			/// @brief Regularized lower incomplete gamma function P(a,x)
			static Real incompleteGammaP(Real a, Real x)
			{
				if (x < 0.0) return 0.0;
				if (x == 0.0) return 0.0;

				if (x < a + 1.0) {
					// Series representation
					Real ap = a;
					Real sum = 1.0 / a;
					Real del = sum;
					for (int n = 1; n <= 200; n++) {
						ap += 1.0;
						del *= x / ap;
						sum += del;
						if (std::abs(del) < std::abs(sum) * 1e-15) break;
					}
					return sum * std::exp(-x + a * std::log(x) - std::lgamma(a));
				} else {
					// Continued fraction representation
					Real b = x + 1.0 - a;
					Real c = 1.0 / 1e-30;
					Real d = 1.0 / b;
					Real h = d;
					for (int n = 1; n <= 200; n++) {
						Real an = -n * (n - a);
						b += 2.0;
						d = an * d + b;
						if (std::abs(d) < 1e-30) d = 1e-30;
						c = b + an / c;
						if (std::abs(c) < 1e-30) c = 1e-30;
						d = 1.0 / d;
						Real del = d * c;
						h *= del;
						if (std::abs(del - 1.0) < 1e-15) break;
					}
					return 1.0 - h * std::exp(-x + a * std::log(x) - std::lgamma(a));
				}
			}

			/// @brief Inverse standard normal CDF (probit)
			static Real inverseStandardNormalCdf(Real p)
			{
				// Rational approximation from Abramowitz and Stegun
				const Real a[] = { -3.969683028665376e+01, 2.209460984245205e+02,
				                   -2.759285104469687e+02, 1.383577518672690e+02,
				                   -3.066479806614716e+01, 2.506628277459239e+00 };
				const Real b[] = { -5.447609879822406e+01, 1.615858368580409e+02,
				                   -1.556989798598866e+02, 6.680131188771972e+01,
				                   -1.328068155288572e+01 };
				const Real c[] = { -7.784894002430293e-03, -3.223964580411365e-01,
				                   -2.400758277161838e+00, -2.549732539343734e+00,
				                    4.374664141464968e+00,  2.938163982698783e+00 };
				const Real d[] = { 7.784695709041462e-03, 3.224671290700398e-01,
				                   2.445134137142996e+00, 3.754408661907416e+00 };

				const Real pLow = 0.02425, pHigh = 1.0 - pLow;
				Real q, r;

				if (p < pLow) {
					q = std::sqrt(-2.0 * std::log(p));
					return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
					        ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1.0);
				} else if (p <= pHigh) {
					q = p - 0.5;
					r = q * q;
					return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
					       (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1.0);
				} else {
					q = std::sqrt(-2.0 * std::log(1.0 - p));
					return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
					         ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1.0);
				}
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
				if (alpha <= 0.0)
					throw StatisticsError("Alpha parameter must be positive in BetaDistribution");
				if (beta <= 0.0)
					throw StatisticsError("Beta parameter must be positive in BetaDistribution");
			}

			/// @brief Probability density function (PDF)
			Real pdf(Real x) const
			{
				if (x < 0.0 || x > 1.0) return 0.0;
				if (x == 0.0) return (alpha < 1.0) ? std::numeric_limits<Real>::infinity()
				                                  : (alpha == 1.0 ? beta : 0.0);
				if (x == 1.0) return (beta < 1.0) ? std::numeric_limits<Real>::infinity()
				                                 : (beta == 1.0 ? alpha : 0.0);

				// Log-space for numerical stability
				Real logPdf = (alpha - 1.0) * std::log(x) + (beta - 1.0) * std::log(1.0 - x)
				            - std::lgamma(alpha) - std::lgamma(beta) + std::lgamma(alpha + beta);
				return std::exp(logPdf);
			}

			/// @brief Cumulative distribution function (CDF)
			/// Uses regularized incomplete beta function
			Real cdf(Real x) const
			{
				if (x <= 0.0) return 0.0;
				if (x >= 1.0) return 1.0;
				return incompleteBetaRegularized(x, alpha, beta);
			}

			/// @brief Inverse CDF (quantile function) - uses Newton-Raphson
			Real inverseCdf(Real p) const
			{
				if (p <= 0.0) return 0.0;
				if (p >= 1.0) return 1.0;

				// Initial guess
				Real x = 0.5;
				if (alpha > 1.0 && beta > 1.0) {
					// Mode-based initial guess
					x = (alpha - 1.0) / (alpha + beta - 2.0);
				}

				// Newton-Raphson
				for (int iter = 0; iter < 50; iter++) {
					Real f = cdf(x) - p;
					Real fp = pdf(x);
					if (std::abs(fp) < 1e-30) break;
					Real dx = f / fp;
					x -= dx;
					if (x <= 0.0) x = 1e-10;
					if (x >= 1.0) x = 1.0 - 1e-10;
					if (std::abs(dx) < 1e-12) break;
				}
				return x;
			}

			Real mean() const { return alpha / (alpha + beta); }
			Real variance() const 
			{ 
				Real ab = alpha + beta;
				return (alpha * beta) / (ab * ab * (ab + 1.0));
			}

		private:
			/// @brief Regularized incomplete beta function I_x(a,b)
			static Real incompleteBetaRegularized(Real x, Real a, Real b)
			{
				if (x <= 0.0) return 0.0;
				if (x >= 1.0) return 1.0;

				// Use symmetry if needed for convergence
				if (x > (a + 1.0) / (a + b + 2.0)) {
					return 1.0 - incompleteBetaRegularized(1.0 - x, b, a);
				}

				Real bt = std::exp(std::lgamma(a + b) - std::lgamma(a) - std::lgamma(b)
				                 + a * std::log(x) + b * std::log(1.0 - x));

				// Continued fraction (Lentz's method)
				Real qab = a + b, qap = a + 1.0, qam = a - 1.0;
				Real c = 1.0, d = 1.0 - qab * x / qap;
				if (std::abs(d) < 1e-30) d = 1e-30;
				d = 1.0 / d;
				Real h = d;

				for (int m = 1; m <= 200; m++) {
					int m2 = 2 * m;
					Real aa = m * (b - m) * x / ((qam + m2) * (a + m2));
					d = 1.0 + aa * d; if (std::abs(d) < 1e-30) d = 1e-30;
					c = 1.0 + aa / c; if (std::abs(c) < 1e-30) c = 1e-30;
					d = 1.0 / d;
					h *= d * c;

					aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
					d = 1.0 + aa * d; if (std::abs(d) < 1e-30) d = 1e-30;
					c = 1.0 + aa / c; if (std::abs(c) < 1e-30) c = 1e-30;
					d = 1.0 / d;
					Real del = d * c;
					h *= del;
					if (std::abs(del - 1.0) < 1e-15) break;
				}
				return bt * h / a;
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
				if (shape <= 0.0)
					throw StatisticsError("Shape parameter must be positive in WeibullDistribution");
				if (scale <= 0.0)
					throw StatisticsError("Scale parameter must be positive in WeibullDistribution");
			}

			/// @brief Probability density function (PDF)
			Real pdf(Real x) const
			{
				if (x < 0.0) return 0.0;
				if (x == 0.0) return (shape < 1.0) ? std::numeric_limits<Real>::infinity()
				                                  : (shape == 1.0 ? 1.0/scale : 0.0);
				Real z = x / scale;
				return (shape / scale) * std::pow(z, shape - 1.0) * std::exp(-std::pow(z, shape));
			}

			/// @brief Cumulative distribution function (CDF)
			Real cdf(Real x) const
			{
				if (x <= 0.0) return 0.0;
				return 1.0 - std::exp(-std::pow(x / scale, shape));
			}

			/// @brief Inverse CDF (quantile function) - closed form
			Real inverseCdf(Real p) const
			{
				if (p <= 0.0) return 0.0;
				if (p >= 1.0) return std::numeric_limits<Real>::infinity();
				return scale * std::pow(-std::log(1.0 - p), 1.0 / shape);
			}

			Real mean() const { return scale * std::tgamma(1.0 + 1.0 / shape); }
			Real variance() const 
			{
				Real g1 = std::tgamma(1.0 + 1.0 / shape);
				Real g2 = std::tgamma(1.0 + 2.0 / shape);
				return scale * scale * (g2 - g1 * g1);
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
				if (alpha <= 0.0)
					throw StatisticsError("Shape parameter must be positive in ParetoDistribution");
				if (xm <= 0.0)
					throw StatisticsError("Scale parameter must be positive in ParetoDistribution");
			}

			/// @brief Probability density function (PDF)
			Real pdf(Real x) const
			{
				if (x < xm) return 0.0;
				return alpha * std::pow(xm, alpha) / std::pow(x, alpha + 1.0);
			}

			/// @brief Cumulative distribution function (CDF)
			Real cdf(Real x) const
			{
				if (x < xm) return 0.0;
				return 1.0 - std::pow(xm / x, alpha);
			}

			/// @brief Inverse CDF (quantile function) - closed form
			Real inverseCdf(Real p) const
			{
				if (p <= 0.0) return xm;
				if (p >= 1.0) return std::numeric_limits<Real>::infinity();
				return xm / std::pow(1.0 - p, 1.0 / alpha);
			}

			/// @brief Mean (only defined for α > 1)
			Real mean() const
			{
				if (alpha <= 1.0)
					throw StatisticsError("Mean undefined for alpha <= 1 in ParetoDistribution");
				return alpha * xm / (alpha - 1.0);
			}

			/// @brief Variance (only defined for α > 2)
			Real variance() const
			{
				if (alpha <= 2.0)
					throw StatisticsError("Variance undefined for alpha <= 2 in ParetoDistribution");
				return (xm * xm * alpha) / ((alpha - 1.0) * (alpha - 1.0) * (alpha - 2.0));
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
				if (sigma <= 0.0)
					throw StatisticsError("Sigma must be positive in LogNormalDistribution");
			}

			/// @brief Probability density function (PDF)
			Real pdf(Real x) const
			{
				if (x <= 0.0) return 0.0;
				Real z = (std::log(x) - mu) / sigma;
				return std::exp(-0.5 * z * z) / (x * sigma * std::sqrt(2.0 * Constants::PI));
			}

			/// @brief Cumulative distribution function (CDF)
			Real cdf(Real x) const
			{
				if (x <= 0.0) return 0.0;
				Real z = (std::log(x) - mu) / (sigma * std::sqrt(2.0));
				return 0.5 * (1.0 + std::erf(z));
			}

			/// @brief Inverse CDF (quantile function)
			Real inverseCdf(Real p) const
			{
				if (p <= 0.0) return 0.0;
				if (p >= 1.0) return std::numeric_limits<Real>::infinity();
				// Use inverse normal for ln(x)
				Real z = inverseStandardNormalCdf(p);
				return std::exp(mu + sigma * z);
			}

			Real mean() const { return std::exp(mu + 0.5 * sigma * sigma); }
			Real variance() const 
			{
				Real expSigma2 = std::exp(sigma * sigma);
				return std::exp(2.0 * mu + sigma * sigma) * (expSigma2 - 1.0);
			}
			Real median() const { return std::exp(mu); }

		private:
			static Real inverseStandardNormalCdf(Real p)
			{
				const Real a[] = { -3.969683028665376e+01, 2.209460984245205e+02,
				                   -2.759285104469687e+02, 1.383577518672690e+02,
				                   -3.066479806614716e+01, 2.506628277459239e+00 };
				const Real b[] = { -5.447609879822406e+01, 1.615858368580409e+02,
				                   -1.556989798598866e+02, 6.680131188771972e+01,
				                   -1.328068155288572e+01 };
				const Real c[] = { -7.784894002430293e-03, -3.223964580411365e-01,
				                   -2.400758277161838e+00, -2.549732539343734e+00,
				                    4.374664141464968e+00,  2.938163982698783e+00 };
				const Real d[] = { 7.784695709041462e-03, 3.224671290700398e-01,
				                   2.445134137142996e+00, 3.754408661907416e+00 };

				const Real pLow = 0.02425, pHigh = 1.0 - pLow;
				Real q, r;

				if (p < pLow) {
					q = std::sqrt(-2.0 * std::log(p));
					return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
					        ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1.0);
				} else if (p <= pHigh) {
					q = p - 0.5;
					r = q * q;
					return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
					       (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1.0);
				} else {
					q = std::sqrt(-2.0 * std::log(1.0 - p));
					return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
					         ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1.0);
				}
			}
		};

  }
}

#endif // MML_DISTRIBUTIONS_H

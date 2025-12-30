///////////////////////////////////////////////////////////////////////////////////////////
// CoreDistributions.h
// 
// Core probability distributions for MML (MinimalMathLibrary)
// Extracted from Statistics.h as part of refactoring to improve organization
//
// Contains: Normal, Student's T, Chi-Square, and F distributions
// These are fundamental distributions used throughout statistical inference
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_CORE_DISTRIBUTIONS_H
#define MML_CORE_DISTRIBUTIONS_H

#include "../Statistics.h"

namespace MML
{
	namespace Statistics
	{
		/**
		 * @brief Normal (Gaussian) distribution
		 * 
		 * The normal distribution is the most important continuous distribution,
		 * central to the Central Limit Theorem and statistical inference.
		 * 
		 * PDF: f(x) = (1 / (σ√(2π))) · exp(-½((x-μ)/σ)²)
		 * where μ is the mean and σ is the standard deviation
		 * 
		 * The standard normal has μ=0 and σ=1.
		 */
		struct NormalDistribution
		{
			Real mu;    // Mean (location parameter)
			Real sigma; // Standard deviation (scale parameter)

			/**
			 * @brief Construct a normal distribution
			 * @param mean Mean parameter μ (default: 0)
			 * @param stddev Standard deviation σ (default: 1, must be > 0)
			 */
			NormalDistribution(Real mean = 0.0, Real stddev = 1.0)
				: mu(mean), sigma(stddev)
			{
				if (sigma <= 0.0)
					throw StatisticsError("Standard deviation must be positive in NormalDistribution");
			}

			/**
			 * @brief Probability density function (PDF)
			 * @param x Value at which to evaluate the PDF
			 * @return Probability density at x
			 */
			Real pdf(Real x) const
			{
				Real z = (x - mu) / sigma;
				return std::exp(-0.5 * z * z) / (sigma * std::sqrt(2.0 * Constants::PI));
			}

			/**
			 * @brief Cumulative distribution function (CDF)
			 * 
			 * Uses the error function: Φ(x) = 0.5 * (1 + erf(x/√2))
			 * 
			 * @param x Value at which to evaluate the CDF
			 * @return Probability that X <= x
			 */
			Real cdf(Real x) const
			{
				Real z = (x - mu) / (sigma * std::sqrt(2.0));
				return 0.5 * (1.0 + std::erf(z));
			}

			/**
			 * @brief Inverse cumulative distribution function (quantile/probit function)
			 * 
			 * Uses Abramowitz and Stegun approximation for the inverse error function.
			 * 
			 * @param p Probability (must be in (0, 1))
			 * @return Value x such that P(X <= x) = p
			 */
			Real inverseCdf(Real p) const
			{
				if (p <= 0.0 || p >= 1.0)
					throw StatisticsError("Probability must be in (0, 1) in NormalDistribution::inverseCdf");

				// Use rational approximation for inverse normal CDF
				// Based on Abramowitz and Stegun approximation 26.2.23
				Real z = inverseStandardNormalCdf(p);
				return mu + sigma * z;
			}

			/**
			 * @brief Z-score: standardize a value
			 * @param x Raw value
			 * @return Standardized value (z-score)
			 */
			Real zScore(Real x) const
			{
				return (x - mu) / sigma;
			}

			/**
			 * @brief Convert z-score back to raw value
			 * @param z Standardized value
			 * @return Raw value
			 */
			Real fromZScore(Real z) const
			{
				return mu + sigma * z;
			}

			/**
			 * @brief Mean of the distribution
			 */
			Real mean() const { return mu; }

			/**
			 * @brief Variance of the distribution
			 */
			Real variance() const { return sigma * sigma; }

			/**
			 * @brief Standard deviation of the distribution
			 */
			Real stddev() const { return sigma; }

		private:
			/**
			 * @brief Inverse of the standard normal CDF (probit function)
			 * 
			 * Uses rational approximation from Abramowitz and Stegun (26.2.23)
			 * with refinement for tails. Accurate to ~1e-9.
			 */
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

		/**
		 * @brief Standard normal distribution (mean=0, stddev=1)
		 * 
		 * Convenience functions for the standard normal distribution.
		 */
		struct StandardNormal
		{
			/**
			 * @brief PDF of standard normal
			 */
			static Real pdf(Real z)
			{
				return std::exp(-0.5 * z * z) / std::sqrt(2.0 * Constants::PI);
			}

			/**
			 * @brief CDF of standard normal (Φ function)
			 */
			static Real cdf(Real z)
			{
				return 0.5 * (1.0 + std::erf(z / std::sqrt(2.0)));
			}

			/**
			 * @brief Inverse CDF (probit function)
			 */
			static Real inverseCdf(Real p)
			{
				NormalDistribution standard(0.0, 1.0);
				return standard.inverseCdf(p);
			}

			/**
			 * @brief Survival function: P(Z > z) = 1 - Φ(z)
			 */
			static Real sf(Real z)
			{
				return 1.0 - cdf(z);
			}

			/**
			 * @brief Two-tailed p-value: P(|Z| > |z|)
			 */
			static Real twoTailedPValue(Real z)
			{
				return 2.0 * sf(std::abs(z));
			}
		};

		/**
		 * @brief Student's t-distribution
		 * 
		 * The t-distribution is used for inference about means when the sample size
		 * is small and/or the population standard deviation is unknown. It has heavier
		 * tails than the normal distribution, with the tail weight controlled by
		 * degrees of freedom (df).
		 * 
		 * As df → ∞, the t-distribution approaches the standard normal distribution.
		 * 
		 * PDF: f(x) = Γ((ν+1)/2) / (√(νπ) · Γ(ν/2)) · (1 + x²/ν)^(-(ν+1)/2)
		 * where ν is the degrees of freedom
		 */
		struct TDistribution
		{
			int df;  // Degrees of freedom

			/**
			 * @brief Construct a Student's t-distribution
			 * @param degreesOfFreedom Degrees of freedom (must be >= 1)
			 */
			TDistribution(int degreesOfFreedom) : df(degreesOfFreedom)
			{
				if (df < 1)
					throw StatisticsError("Degrees of freedom must be at least 1 in TDistribution");
			}

			/**
			 * @brief Probability density function (PDF)
			 * @param x Value at which to evaluate the PDF
			 * @return Probability density at x
			 */
			Real pdf(Real x) const
			{
				Real nu = static_cast<Real>(df);
				
				// Log-space computation for numerical stability
				Real logNumerator = std::lgamma((nu + 1.0) / 2.0);
				Real logDenominator = 0.5 * std::log(nu * Constants::PI) + std::lgamma(nu / 2.0);
				Real logBase = -(nu + 1.0) / 2.0 * std::log(1.0 + x * x / nu);
				
				return std::exp(logNumerator - logDenominator + logBase);
			}

			/**
			 * @brief Cumulative distribution function (CDF)
			 * 
			 * Uses the regularized incomplete beta function.
			 * 
			 * @param x Value at which to evaluate the CDF
			 * @return Probability that T <= x
			 */
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

			/**
			 * @brief Inverse cumulative distribution function (quantile function)
			 * 
			 * Uses iterative root finding for the inverse.
			 * 
			 * @param p Probability (must be in (0, 1))
			 * @return Value t such that P(T <= t) = p
			 */
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
					
					if (std::abs(delta) < 1e-10)
						break;
				}
				
				return t;
			}

			/**
			 * @brief Survival function: P(T > t)
			 */
			Real sf(Real x) const
			{
				return 1.0 - cdf(x);
			}

			/**
			 * @brief Two-tailed p-value: P(|T| > |t|)
			 * 
			 * This is the standard p-value used in t-tests.
			 */
			Real twoTailedPValue(Real t) const
			{
				return 2.0 * sf(std::abs(t));
			}

			/**
			 * @brief Critical value for given significance level (two-tailed)
			 * 
			 * Returns t* such that P(|T| > t*) = alpha
			 * 
			 * @param alpha Significance level (e.g., 0.05 for 95% confidence)
			 * @return Critical value t*
			 */
			Real criticalValue(Real alpha) const
			{
				if (alpha <= 0.0 || alpha >= 1.0)
					throw StatisticsError("Alpha must be in (0, 1) in TDistribution::criticalValue");
				
				// Two-tailed: find t* where P(T > t*) = alpha/2
				return inverseCdf(1.0 - alpha / 2.0);
			}

			/**
			 * @brief Mean of the distribution (if it exists)
			 * 
			 * Mean exists only for df >= 2, and equals 0.
			 */
			Real mean() const
			{
				if (df < 2)
					throw StatisticsError("Mean does not exist for df < 2 in TDistribution");
				return 0.0;
			}

			/**
			 * @brief Variance of the distribution (if it exists)
			 * 
			 * Variance exists only for df >= 3.
			 * Var(T) = ν/(ν-2) for ν > 2
			 */
			Real variance() const
			{
				if (df <= 2)
					throw StatisticsError("Variance does not exist for df <= 2 in TDistribution");
				
				Real nu = static_cast<Real>(df);
				return nu / (nu - 2.0);
			}

		private:
			/**
			 * @brief Regularized incomplete beta function I_x(a, b)
			 * 
			 * Uses continued fraction expansion (Lentz's method).
			 * Based on Numerical Recipes formula 6.4.5
			 */
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

		/**
		 * @brief Chi-square (χ²) distribution
		 * 
		 * The chi-square distribution is the distribution of a sum
		 * of a sum of the squares of k independent standard normal random variables.
		 * 
		 * PDF: f(x; k) = (1 / (2^(k/2) * Γ(k/2))) * x^(k/2-1) * e^(-x/2) for x ≥ 0
		 * Mean: k
		 * Variance: 2k
		 * 
		 * Used extensively in hypothesis testing (goodness-of-fit, independence tests)
		 * 
		 * @param df Degrees of freedom (must be positive)
		 */
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

			/**
			 * @brief Probability Density Function
			 */
			Real pdf(Real x) const
			{
				if (x < 0.0) return 0.0;
				if (x == 0.0) return (_df == 2) ? 0.5 : 0.0;

				Real k2 = _df / 2.0;
				return std::exp((k2 - 1.0) * std::log(x) - x / 2.0 - k2 * std::log(2.0) - std::lgamma(k2));
			}

			/**
			 * @brief Cumulative Distribution Function (CDF)
			 * P(X ≤ x) using the incomplete gamma function
			 */
			Real cdf(Real x) const
			{
				if (x <= 0.0) return 0.0;
				
				// CDF = P(k/2, x/2) where P is regularized lower incomplete gamma
				return incompleteGammaP(_df / 2.0, x / 2.0);
			}

			/**
			 * @brief Right-tail probability P(X ≥ x)
			 * This is the p-value for chi-square tests
			 */
			Real rightTailPValue(Real chiSquare) const
			{
				return 1.0 - cdf(chiSquare);
			}

			/**
			 * @brief Critical value for given significance level (right-tail test)
			 * Returns x such that P(X ≥ x) = alpha
			 */
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
			/**
			 * @brief Regularized lower incomplete gamma function P(a,x)
			 * P(a,x) = γ(a,x) / Γ(a)
			 */
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

			/**
			 * @brief Series expansion for incomplete gamma
			 */
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

			/**
			 * @brief Continued fraction for incomplete gamma
			 */
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

    /**
     * @brief F-Distribution (Fisher-Snedecor Distribution)
     * 
     * The F-distribution is the ratio of two scaled chi-square distributions.
     * F = (χ²₁/df₁) / (χ²₂/df₂)
     * 
     * Used extensively in ANOVA, regression analysis, and comparing variances.
     * 
     * PDF: f(x; d₁, d₂) = [complex formula involving beta function]
     * Mean: d₂/(d₂-2) for d₂ > 2
     * Variance: 2d₂²(d₁+d₂-2) / [d₁(d₂-2)²(d₂-4)] for d₂ > 4
     * 
     * @param df1 Numerator degrees of freedom
     * @param df2 Denominator degrees of freedom
     */
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

      /**
       * @brief Probability Density Function
       */
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

      /**
       * @brief Cumulative Distribution Function
       * Uses relationship with regularized incomplete beta function
       */
      Real cdf(Real x) const
      {
        if (x <= 0.0) return 0.0;

        // F-CDF = I_y(df1/2, df2/2) where y = (df1*x)/(df1*x + df2)
        // I is regularized incomplete beta function
        Real y = (_df1 * x) / (_df1 * x + _df2);
        return incompleteBeta(_df1 / 2.0, _df2 / 2.0, y);
      }

      /**
       * @brief Right-tail probability P(F ≥ f)
       * This is the p-value for F-tests (ANOVA, regression)
       */
      Real rightTailPValue(Real f) const
      {
        return 1.0 - cdf(f);
      }

      /**
       * @brief Critical value for given significance level (right-tail test)
       * Returns f such that P(F ≥ f) = alpha
       */
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
      /**
       * @brief Regularized incomplete beta function I_x(a,b)
       */
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

      /**
       * @brief Continued fraction for incomplete beta
       */
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
    }
}

#endif // MML_CORE_DISTRIBUTIONS_H

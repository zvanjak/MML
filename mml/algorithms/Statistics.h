///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Statistics.h                                                        ///
///  Description: Statistical functions (mean, variance, correlation, distributions)  ///
///               Hypothesis testing and probability distributions                    ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_STATISTICS_H
#define MML_STATISTICS_H

#include "MMLBase.h"
#include "MMLExceptions.h"

#include "base/Vector.h"
#include "base/Matrix.h"

#include <random>

namespace MML
{
	namespace Statistics
	{
		static Real Avg(const Vector<Real>& data)
		{
			Real outAvg = 0.0;
			int n = data.size();
			
			if (n <= 0)
				throw StatisticsError("Vector size must be greater than 0 in Avg");

			for (int j = 0; j < n; j++)
				outAvg += data[j];
			outAvg /= n;

			return outAvg;
		}

		static void AvgVar(const Vector<Real>& data, Real& outAvg, Real& outVar)
		{
			Real s, ep;
			int j, n = data.size();

			if (n <= 0)
				throw StatisticsError("Vector size must be greater than 0 in AvgVar");

			outAvg = Avg(data);

			outVar = ep = 0.0;
			for (j = 0; j < n; j++) {
				s = data[j] - outAvg;
				ep += s;
				outVar += s * s;
			}
			outVar = (outVar - ep * ep / n) / (n - 1);
		}
		static void AvgStdDev(const Vector<Real>& data, Real& outAvg, Real& outStdDev)
		{
			Real var;
			int n = data.size();
			if (n <= 0)
				throw StatisticsError("Vector size must be greater than 0 in AvgStdDev");
			AvgVar(data, outAvg, var);
			outStdDev = sqrt(var);
		}

		// Helper functions for t-tests (aliases for clarity)
		static Real Mean(const Vector<Real>& data) { return Avg(data); }
		
		static Real Variance(const Vector<Real>& data)
		{
			Real avg, var;
			AvgVar(data, avg, var);
			return var;
		}
		
		static Real StdDev(const Vector<Real>& data)
		{
			Real avg, stddev;
			AvgStdDev(data, avg, stddev);
			return stddev;
		}

		static void Moments(const Vector<Real>& data, Real& ave, Real& adev, Real& sdev, Real& var, Real& skew, Real& curt)
		{
			int j, n = data.size();
			Real ep = 0.0, s, p;

			if (n <= 1)
				throw StatisticsError("Vector size must be greater than 1 in Moments");

			s = 0.0;
			for (j = 0; j < n; j++)
				s += data[j];
			ave = s / n;

			adev = var = skew = curt = 0.0;
			for (j = 0; j < n; j++) {
				adev += std::abs(s = data[j] - ave);
				ep += s;
				var += (p = s * s);
				skew += (p *= s);
				curt += (p *= s);
			}
			adev /= n;
			var = (var - ep * ep / n) / (n - 1);
			sdev = sqrt(var);

			if (var != 0.0) {
				skew /= (n * var * sdev);
				curt = curt / (n * var * var) - 3.0;
			}
			else
				throw StatisticsError("No skew/kurtosis when variance = 0 (in Moments)");
		}

		/*********************************************************************/
		/*****                  Order Statistics                         *****/
		/*********************************************************************/

		/**
		 * @brief Compute the median of a dataset
		 * 
		 * The median is the middle value when data is sorted. For even-sized
		 * datasets, returns the average of the two middle values.
		 * 
		 * @param data Input vector (will be copied and sorted internally)
		 * @return The median value
		 * @throws StatisticsError if data is empty
		 * 
		 * Complexity: O(n log n) due to sorting
		 */
		static Real Median(const Vector<Real>& data)
		{
			int n = data.size();
			if (n <= 0)
				throw StatisticsError("Vector size must be greater than 0 in Median");

			// Create a sorted copy
			std::vector<Real> sorted(n);
			for (int i = 0; i < n; i++)
				sorted[i] = data[i];
			std::sort(sorted.begin(), sorted.end());

			if (n % 2 == 0)
				return (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0;
			else
				return sorted[n / 2];
		}

		/**
		 * @brief Compute a percentile of a dataset
		 * 
		 * Uses linear interpolation between data points (similar to Excel PERCENTILE.INC).
		 * The p-th percentile is the value below which p percent of the data falls.
		 * 
		 * @param data Input vector (will be copied and sorted internally)
		 * @param p Percentile to compute (must be in [0, 100])
		 * @return The p-th percentile value
		 * @throws StatisticsError if data is empty or p is out of range
		 * 
		 * Complexity: O(n log n) due to sorting
		 */
		static Real Percentile(const Vector<Real>& data, Real p)
		{
			int n = data.size();
			if (n <= 0)
				throw StatisticsError("Vector size must be greater than 0 in Percentile");
			if (p < 0.0 || p > 100.0)
				throw StatisticsError("Percentile must be in range [0, 100]");

			// Create a sorted copy
			std::vector<Real> sorted(n);
			for (int i = 0; i < n; i++)
				sorted[i] = data[i];
			std::sort(sorted.begin(), sorted.end());

			// Handle edge cases
			if (p == 0.0) return sorted[0];
			if (p == 100.0) return sorted[n - 1];

			// Linear interpolation (Excel PERCENTILE.INC method)
			Real rank = (p / 100.0) * (n - 1);
			int lower = static_cast<int>(std::floor(rank));
			int upper = static_cast<int>(std::ceil(rank));
			
			if (lower == upper)
				return sorted[lower];

			Real fraction = rank - lower;
			return sorted[lower] + fraction * (sorted[upper] - sorted[lower]);
		}

		/**
		 * @brief Compute quartiles (Q1, Q2/Median, Q3) of a dataset
		 * 
		 * Q1 = 25th percentile (first quartile)
		 * Q2 = 50th percentile (median)
		 * Q3 = 75th percentile (third quartile)
		 * 
		 * @param data Input vector
		 * @param[out] q1 First quartile (25th percentile)
		 * @param[out] median Second quartile (50th percentile)
		 * @param[out] q3 Third quartile (75th percentile)
		 * @throws StatisticsError if data is empty
		 */
		static void Quartiles(const Vector<Real>& data, Real& q1, Real& median, Real& q3)
		{
			int n = data.size();
			if (n <= 0)
				throw StatisticsError("Vector size must be greater than 0 in Quartiles");

			q1 = Percentile(data, 25.0);
			median = Percentile(data, 50.0);
			q3 = Percentile(data, 75.0);
		}

		/**
		 * @brief Compute the range of a dataset
		 * 
		 * Range = max(data) - min(data)
		 * 
		 * @param data Input vector
		 * @return The range (maximum minus minimum)
		 * @throws StatisticsError if data is empty
		 * 
		 * Complexity: O(n)
		 */
		static Real Range(const Vector<Real>& data)
		{
			int n = data.size();
			if (n <= 0)
				throw StatisticsError("Vector size must be greater than 0 in Range");

			Real minVal = data[0];
			Real maxVal = data[0];
			for (int i = 1; i < n; i++) {
				if (data[i] < minVal) minVal = data[i];
				if (data[i] > maxVal) maxVal = data[i];
			}
			return maxVal - minVal;
		}

		/**
		 * @brief Compute the Interquartile Range (IQR)
		 * 
		 * IQR = Q3 - Q1 = 75th percentile - 25th percentile
		 * 
		 * The IQR is a robust measure of spread, less sensitive to outliers
		 * than the range or standard deviation.
		 * 
		 * @param data Input vector
		 * @return The interquartile range
		 * @throws StatisticsError if data is empty
		 */
		static Real IQR(const Vector<Real>& data)
		{
			int n = data.size();
			if (n <= 0)
				throw StatisticsError("Vector size must be greater than 0 in IQR");

			return Percentile(data, 75.0) - Percentile(data, 25.0);
		}

		/**
		 * @brief Compute minimum and maximum values
		 * 
		 * @param data Input vector
		 * @param[out] minVal Minimum value in the data
		 * @param[out] maxVal Maximum value in the data
		 * @throws StatisticsError if data is empty
		 * 
		 * Complexity: O(n)
		 */
		static void MinMax(const Vector<Real>& data, Real& minVal, Real& maxVal)
		{
			int n = data.size();
			if (n <= 0)
				throw StatisticsError("Vector size must be greater than 0 in MinMax");

			minVal = data[0];
			maxVal = data[0];
			for (int i = 1; i < n; i++) {
				if (data[i] < minVal) minVal = data[i];
				if (data[i] > maxVal) maxVal = data[i];
			}
		}

		/*********************************************************************/
		/*****               Robust Statistics                           *****/
		/*********************************************************************/

		/**
		 * @brief Compute the mode (most frequent value) of a dataset
		 * 
		 * For continuous data, finds the value that appears most frequently.
		 * If multiple values have the same highest frequency, returns the smallest.
		 * For truly continuous data where no values repeat, returns the first value.
		 * 
		 * @param data Input vector
		 * @return The mode (most frequent value)
		 * @throws StatisticsError if data is empty
		 * 
		 * Complexity: O(n log n) due to sorting
		 * 
		 * @note For continuous distributions, consider using kernel density estimation
		 *       or histogram-based methods instead.
		 */
		static Real Mode(const Vector<Real>& data)
		{
			int n = data.size();
			if (n <= 0)
				throw StatisticsError("Vector size must be greater than 0 in Mode");
			if (n == 1)
				return data[0];

			// Create a sorted copy
			std::vector<Real> sorted(n);
			for (int i = 0; i < n; i++)
				sorted[i] = data[i];
			std::sort(sorted.begin(), sorted.end());

			Real mode = sorted[0];
			int maxCount = 1;
			int currentCount = 1;

			for (int i = 1; i < n; i++) {
				if (sorted[i] == sorted[i - 1]) {
					currentCount++;
				} else {
					if (currentCount > maxCount) {
						maxCount = currentCount;
						mode = sorted[i - 1];
					}
					currentCount = 1;
				}
			}
			// Check the last run
			if (currentCount > maxCount) {
				mode = sorted[n - 1];
			}

			return mode;
		}

		/**
		 * @brief Compute the trimmed mean (truncated mean)
		 * 
		 * The trimmed mean excludes a percentage of the smallest and largest values
		 * before computing the mean. This provides robustness against outliers.
		 * 
		 * @param data Input vector
		 * @param trimPercent Percentage to trim from each end (0-50)
		 *                    e.g., 10 means remove bottom 10% and top 10%
		 * @return The trimmed mean
		 * @throws StatisticsError if data is empty or trimPercent is invalid
		 * 
		 * Complexity: O(n log n) due to sorting
		 */
		static Real TrimmedMean(const Vector<Real>& data, Real trimPercent)
		{
			int n = data.size();
			if (n <= 0)
				throw StatisticsError("Vector size must be greater than 0 in TrimmedMean");
			if (trimPercent < 0.0 || trimPercent >= 50.0)
				throw StatisticsError("trimPercent must be in range [0, 50) in TrimmedMean");

			// Create a sorted copy
			std::vector<Real> sorted(n);
			for (int i = 0; i < n; i++)
				sorted[i] = data[i];
			std::sort(sorted.begin(), sorted.end());

			// Calculate how many to trim from each end
			int trimCount = static_cast<int>(std::floor(n * trimPercent / 100.0));
			
			// Calculate mean of remaining values
			Real sum = 0.0;
			int count = 0;
			for (int i = trimCount; i < n - trimCount; i++) {
				sum += sorted[i];
				count++;
			}

			if (count == 0)
				throw StatisticsError("No values remain after trimming in TrimmedMean");

			return sum / count;
		}

		/**
		 * @brief Compute the Median Absolute Deviation (MAD)
		 * 
		 * MAD = median(|x_i - median(x)|)
		 * 
		 * MAD is a robust measure of variability. It is less sensitive to
		 * outliers than standard deviation. To estimate standard deviation
		 * for normal data, multiply MAD by 1.4826.
		 * 
		 * @param data Input vector
		 * @return The median absolute deviation
		 * @throws StatisticsError if data is empty
		 * 
		 * Complexity: O(n log n) due to sorting
		 */
		static Real MAD(const Vector<Real>& data)
		{
			int n = data.size();
			if (n <= 0)
				throw StatisticsError("Vector size must be greater than 0 in MAD");

			Real med = Median(data);

			// Compute absolute deviations
			Vector<Real> absDeviations(n);
			for (int i = 0; i < n; i++)
				absDeviations[i] = std::abs(data[i] - med);

			return Median(absDeviations);
		}

		/**
		 * @brief Compute the geometric mean
		 * 
		 * GeometricMean = (x_1 * x_2 * ... * x_n)^(1/n) = exp(mean(log(x_i)))
		 * 
		 * The geometric mean is appropriate for data that is multiplicative
		 * in nature (e.g., growth rates, ratios). All values must be positive.
		 * 
		 * @param data Input vector (all values must be positive)
		 * @return The geometric mean
		 * @throws StatisticsError if data is empty or contains non-positive values
		 * 
		 * Complexity: O(n)
		 */
		static Real GeometricMean(const Vector<Real>& data)
		{
			int n = data.size();
			if (n <= 0)
				throw StatisticsError("Vector size must be greater than 0 in GeometricMean");

			Real logSum = 0.0;
			for (int i = 0; i < n; i++) {
				if (data[i] <= 0.0)
					throw StatisticsError("All values must be positive in GeometricMean");
				logSum += std::log(data[i]);
			}

			return std::exp(logSum / n);
		}

		/**
		 * @brief Compute the harmonic mean
		 * 
		 * HarmonicMean = n / (1/x_1 + 1/x_2 + ... + 1/x_n)
		 * 
		 * The harmonic mean is appropriate for averaging rates or ratios
		 * (e.g., speeds, P/E ratios). All values must be positive.
		 * 
		 * @param data Input vector (all values must be positive)
		 * @return The harmonic mean
		 * @throws StatisticsError if data is empty or contains non-positive values
		 * 
		 * Complexity: O(n)
		 */
		static Real HarmonicMean(const Vector<Real>& data)
		{
			int n = data.size();
			if (n <= 0)
				throw StatisticsError("Vector size must be greater than 0 in HarmonicMean");

			Real reciprocalSum = 0.0;
			for (int i = 0; i < n; i++) {
				if (data[i] <= 0.0)
					throw StatisticsError("All values must be positive in HarmonicMean");
				reciprocalSum += 1.0 / data[i];
			}

			return n / reciprocalSum;
		}

		/**
		 * @brief Compute the weighted mean
		 * 
		 * WeightedMean = sum(w_i * x_i) / sum(w_i)
		 * 
		 * @param data Input vector of values
		 * @param weights Vector of weights (must have same size as data)
		 * @return The weighted mean
		 * @throws StatisticsError if data is empty, sizes don't match,
		 *         or total weight is zero/negative
		 * 
		 * Complexity: O(n)
		 */
		static Real WeightedMean(const Vector<Real>& data, const Vector<Real>& weights)
		{
			int n = data.size();
			if (n <= 0)
				throw StatisticsError("Vector size must be greater than 0 in WeightedMean");
			if (weights.size() != n)
				throw StatisticsError("Data and weights must have the same size in WeightedMean");

			Real weightedSum = 0.0;
			Real totalWeight = 0.0;
			for (int i = 0; i < n; i++) {
				if (weights[i] < 0.0)
					throw StatisticsError("Weights must be non-negative in WeightedMean");
				weightedSum += weights[i] * data[i];
				totalWeight += weights[i];
			}

			if (totalWeight <= 0.0)
				throw StatisticsError("Total weight must be positive in WeightedMean");

			return weightedSum / totalWeight;
		}

		/**
		 * @brief Compute the weighted variance
		 * 
		 * Uses reliability weights (frequency weights) formula:
		 * WeightedVar = sum(w_i * (x_i - weighted_mean)^2) / (sum(w_i) - 1)
		 * 
		 * @param data Input vector of values
		 * @param weights Vector of weights (must have same size as data)
		 * @return The weighted variance
		 * @throws StatisticsError if data has fewer than 2 elements, sizes don't match,
		 *         or weights are invalid
		 * 
		 * Complexity: O(n)
		 */
		static Real WeightedVariance(const Vector<Real>& data, const Vector<Real>& weights)
		{
			int n = data.size();
			if (n < 2)
				throw StatisticsError("Vector size must be at least 2 in WeightedVariance");
			if (weights.size() != n)
				throw StatisticsError("Data and weights must have the same size in WeightedVariance");

			Real wMean = WeightedMean(data, weights);

			Real weightedSqSum = 0.0;
			Real totalWeight = 0.0;
			for (int i = 0; i < n; i++) {
				if (weights[i] < 0.0)
					throw StatisticsError("Weights must be non-negative in WeightedVariance");
				Real diff = data[i] - wMean;
				weightedSqSum += weights[i] * diff * diff;
				totalWeight += weights[i];
			}

			if (totalWeight <= 1.0)
				throw StatisticsError("Total weight must be greater than 1 in WeightedVariance");

			return weightedSqSum / (totalWeight - 1.0);
		}

		/*********************************************************************/
		/*****               Correlation & Covariance                    *****/
		/*********************************************************************/

		/**
		 * @brief Compute the sample covariance between two variables
		 * 
		 * Cov(X,Y) = sum((x_i - mean_x)(y_i - mean_y)) / (n - 1)
		 * 
		 * Uses a two-pass algorithm for numerical stability.
		 * 
		 * @param x First variable (vector of values)
		 * @param y Second variable (must have same size as x)
		 * @return Sample covariance
		 * @throws StatisticsError if vectors are empty or have different sizes
		 * 
		 * Complexity: O(n)
		 */
		static Real Covariance(const Vector<Real>& x, const Vector<Real>& y)
		{
			int n = x.size();
			if (n <= 1)
				throw StatisticsError("Vector size must be greater than 1 in Covariance");
			if (y.size() != n)
				throw StatisticsError("Vectors must have the same size in Covariance");

			// Two-pass algorithm for numerical stability
			Real meanX = Avg(x);
			Real meanY = Avg(y);

			Real cov = 0.0;
			for (int i = 0; i < n; i++) {
				cov += (x[i] - meanX) * (y[i] - meanY);
			}

			return cov / (n - 1);
		}

		/**
		 * @brief Compute the Pearson correlation coefficient
		 * 
		 * r = Cov(X,Y) / (StdDev(X) * StdDev(Y))
		 * 
		 * The Pearson correlation measures linear relationship between variables.
		 * Values range from -1 (perfect negative) to +1 (perfect positive).
		 * A value of 0 indicates no linear correlation.
		 * 
		 * @param x First variable (vector of values)
		 * @param y Second variable (must have same size as x)
		 * @return Pearson correlation coefficient r in [-1, 1]
		 * @throws StatisticsError if vectors are empty, have different sizes,
		 *         or either variable has zero variance
		 * 
		 * Complexity: O(n)
		 */
		static Real PearsonCorrelation(const Vector<Real>& x, const Vector<Real>& y)
		{
			int n = x.size();
			if (n <= 1)
				throw StatisticsError("Vector size must be greater than 1 in PearsonCorrelation");
			if (y.size() != n)
				throw StatisticsError("Vectors must have the same size in PearsonCorrelation");

			// Compute means
			Real meanX = Avg(x);
			Real meanY = Avg(y);

			// Compute sums for correlation
			Real sumXY = 0.0;
			Real sumX2 = 0.0;
			Real sumY2 = 0.0;

			for (int i = 0; i < n; i++) {
				Real dx = x[i] - meanX;
				Real dy = y[i] - meanY;
				sumXY += dx * dy;
				sumX2 += dx * dx;
				sumY2 += dy * dy;
			}

			if (sumX2 == 0.0)
				throw StatisticsError("First variable has zero variance in PearsonCorrelation");
			if (sumY2 == 0.0)
				throw StatisticsError("Second variable has zero variance in PearsonCorrelation");

			return sumXY / std::sqrt(sumX2 * sumY2);
		}

		/**
		 * @brief Result structure for correlation with significance test
		 */
		struct CorrelationResult
		{
			Real r;                ///< Correlation coefficient
			Real tStatistic;       ///< t-statistic for significance test
			int degreesOfFreedom;  ///< Degrees of freedom (n-2)
			
			CorrelationResult(Real correlation, Real tStat, int df)
				: r(correlation), tStatistic(tStat), degreesOfFreedom(df) {}
		};

		/**
		 * @brief Compute Pearson correlation with t-test for significance
		 * 
		 * Tests the null hypothesis H0: ρ = 0 (no population correlation)
		 * The t-statistic follows a t-distribution with n-2 degrees of freedom.
		 * 
		 * t = r * sqrt((n-2) / (1-r²))
		 * 
		 * @param x First variable
		 * @param y Second variable
		 * @return CorrelationResult containing r, t-statistic, and degrees of freedom
		 * @throws StatisticsError if inputs are invalid
		 * 
		 * @note Use t-distribution critical values or p-value calculation to assess significance
		 */
		static CorrelationResult PearsonCorrelationWithTest(const Vector<Real>& x, const Vector<Real>& y)
		{
			int n = x.size();
			if (n <= 2)
				throw StatisticsError("Vector size must be greater than 2 for significance test");
			
			Real r = PearsonCorrelation(x, y);
			int df = n - 2;
			
			// Compute t-statistic: t = r * sqrt((n-2)/(1-r²))
			Real r2 = r * r;
			Real tStat;
			if (r2 >= 1.0) {
				// Perfect correlation - t approaches infinity
				tStat = (r > 0) ? std::numeric_limits<Real>::infinity() 
				                : -std::numeric_limits<Real>::infinity();
			} else {
				tStat = r * std::sqrt(static_cast<Real>(df) / (1.0 - r2));
			}
			
			return CorrelationResult(r, tStat, df);
		}

		/**
		 * @brief Compute the coefficient of determination (R²)
		 * 
		 * R² = r² where r is the Pearson correlation coefficient.
		 * Represents the proportion of variance in y explained by x.
		 * 
		 * @param x Independent variable
		 * @param y Dependent variable
		 * @return R² value in [0, 1]
		 */
		static Real RSquared(const Vector<Real>& x, const Vector<Real>& y)
		{
			Real r = PearsonCorrelation(x, y);
			return r * r;
		}

		/**
		 * @brief Compute the covariance matrix for multivariate data
		 * 
		 * For a data matrix where each column is a variable and each row is
		 * an observation, computes the pairwise covariance between all variables.
		 * 
		 * @param data Matrix of size (n_observations x n_variables)
		 * @return Square covariance matrix of size (n_variables x n_variables)
		 * @throws StatisticsError if matrix has insufficient observations
		 * 
		 * @note The diagonal contains variances, off-diagonal contains covariances
		 */
		static Matrix<Real> CovarianceMatrix(const Matrix<Real>& data)
		{
			int n = data.RowNum();  // observations
			int p = data.ColNum();  // variables

			if (n <= 1)
				throw StatisticsError("Need at least 2 observations for CovarianceMatrix");
			if (p <= 0)
				throw StatisticsError("Need at least 1 variable for CovarianceMatrix");

			// Compute means for each variable
			Vector<Real> means(p);
			for (int j = 0; j < p; j++) {
				Real sum = 0.0;
				for (int i = 0; i < n; i++)
					sum += data(i, j);
				means[j] = sum / n;
			}

			// Compute covariance matrix
			Matrix<Real> covMatrix(p, p);
			for (int j1 = 0; j1 < p; j1++) {
				for (int j2 = j1; j2 < p; j2++) {
					Real cov = 0.0;
					for (int i = 0; i < n; i++) {
						cov += (data(i, j1) - means[j1]) * (data(i, j2) - means[j2]);
					}
					cov /= (n - 1);
					covMatrix(j1, j2) = cov;
					covMatrix(j2, j1) = cov;  // Symmetric
				}
			}

			return covMatrix;
		}

		/**
		 * @brief Compute the correlation matrix for multivariate data
		 * 
		 * Similar to covariance matrix, but normalized to [-1, 1] range.
		 * The diagonal is always 1 (variable correlated with itself).
		 * 
		 * @param data Matrix of size (n_observations x n_variables)
		 * @return Square correlation matrix of size (n_variables x n_variables)
		 * @throws StatisticsError if matrix has insufficient observations or zero variance
		 */
		static Matrix<Real> CorrelationMatrix(const Matrix<Real>& data)
		{
			int n = data.RowNum();
			int p = data.ColNum();

			if (n <= 1)
				throw StatisticsError("Need at least 2 observations for CorrelationMatrix");
			if (p <= 0)
				throw StatisticsError("Need at least 1 variable for CorrelationMatrix");

			// Compute covariance matrix first
			Matrix<Real> covMatrix = CovarianceMatrix(data);

			// Extract standard deviations from diagonal
			Vector<Real> stdDevs(p);
			for (int j = 0; j < p; j++) {
				if (covMatrix(j, j) <= 0.0)
					throw StatisticsError("Variable has zero or negative variance in CorrelationMatrix");
				stdDevs[j] = std::sqrt(covMatrix(j, j));
			}

			// Normalize to correlation matrix
			Matrix<Real> corrMatrix(p, p);
			for (int j1 = 0; j1 < p; j1++) {
				for (int j2 = 0; j2 < p; j2++) {
					corrMatrix(j1, j2) = covMatrix(j1, j2) / (stdDevs[j1] * stdDevs[j2]);
				}
			}

			return corrMatrix;
		}
  } // namespace Statistics
}  // namespace MML
#endif

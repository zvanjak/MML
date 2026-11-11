///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Statistics.h                                                        ///
///  Description: Statistical functions (mean, variance, correlation, distributions)  ///
///               Hypothesis testing and probability distributions                    ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_STATISTICS_H
#define MML_STATISTICS_H

#include <mml/MMLBase.h>
#include <algorithm>
#include <numeric>
#include <mml/MMLExceptions.h>

#include <mml/base/Vector/Vector.h>
#include <mml/base/Matrix/Matrix.h>

namespace MML {
	namespace Statistics {
		namespace StatisticsDetail {
			inline void CompensatedAdd(Real value, Real& sum, Real& correction) {
				Real updated = sum + value;
				if (std::abs(sum) >= std::abs(value))
					correction += (sum - updated) + value;
				else
					correction += (value - updated) + sum;
				sum = updated;
			}

			struct CompensatedSum {
				Real sum = 0;
				Real correction = 0;
				void Add(Real value) { CompensatedAdd(value, sum, correction); }
				Real Value() const { return sum + correction; }
			};

			inline void RequireFiniteData(Real value) {
				if (!std::isfinite(value))
					throw StatisticsError("Data must be finite in unweighted statistics");
			}

			inline Real RequireFiniteResult(Real value) {
				if (!std::isfinite(value))
					throw StatisticsError("Result overflow in unweighted statistics");
				return value;
			}

			template<typename ValueAt>
			Real ArithmeticMean(int n, ValueAt valueAt) {
				CompensatedSum raw;
				Real maximum = 0;
				bool finiteSum = true;
				for (int i = 0; i < n; ++i) {
					Real value = valueAt(i);
					RequireFiniteData(value);
					maximum = std::max(maximum, std::abs(value));
					if (finiteSum) {
						finiteSum = std::isfinite(raw.sum + value);
						if (finiteSum)
							raw.Add(value);
					}
				}
				if (finiteSum) {
					Real total = raw.Value();
					return RequireFiniteResult(std::isfinite(total) ? total / n : raw.sum / n + raw.correction / n);
				}

				int exponent;
				std::frexp(maximum, &exponent);
				CompensatedSum scaled, small;
				for (int i = 0; i < n; ++i) {
					Real value = valueAt(i);
					Real normalized = std::scalbn(value, -exponent);
					// Keep values that would lose precision through subnormal scaling in their original units.
					if (value != 0 && !std::isnormal(normalized))
						small.Add(value);
					else
						scaled.Add(normalized);
				}
				return RequireFiniteResult(std::scalbn(scaled.Value() / n, exponent) + small.Value() / n);
			}

			struct CenteredData {
				int exponent = 0;
				Real origin = 0;
				Real meanOffset = 0;
				Real residual = 0;
				Real squareSum = 0;

				Real Center(Real value) const {
					return (std::scalbn(value, -exponent) - origin) - meanOffset;
				}
			};

			template<typename ValueAt>
			CenteredData CenterData(int n, ValueAt valueAt) {
				CenteredData result;
				Real maximum = 0;
				for (int i = 0; i < n; ++i) {
					RequireFiniteData(valueAt(i));
					maximum = std::max(maximum, std::abs(valueAt(i)));
				}
				std::frexp(maximum, &result.exponent);
				result.origin = std::scalbn(valueAt(0), -result.exponent);
				CompensatedSum offsets;
				for (int i = 0; i < n; ++i)
					offsets.Add(std::scalbn(valueAt(i), -result.exponent) - result.origin);
				result.meanOffset = offsets.Value() / n;

				// Corrected two-pass moments (Numerical Recipes), in binary-scaled, origin-relative coordinates.
				CompensatedSum residual, squares;
				for (int i = 0; i < n; ++i) {
					Real centered = result.Center(valueAt(i));
					residual.Add(centered);
					squares.Add(centered * centered);
				}
				result.residual = residual.Value();
				result.squareSum = squares.Value() - result.residual * (result.residual / n);
				if (result.squareSum < 0)
					throw StatisticsError("Negative centered square sum in unweighted statistics");
				return result;
			}

			template<typename XAt, typename YAt>
			Real CrossSum(int n, XAt xAt, YAt yAt, const CenteredData& x, const CenteredData& y) {
				CompensatedSum products;
				for (int i = 0; i < n; ++i)
					products.Add(x.Center(xAt(i)) * y.Center(yAt(i)));
				return products.Value() - x.residual * (y.residual / n);
			}

			inline Real RestoreUnits(Real value, int exponent) {
				return RequireFiniteResult(std::scalbn(value, exponent));
			}

			inline Real Correlation(Real crossSum, const CenteredData& x, const CenteredData& y) {
				if (x.squareSum == 0 || y.squareSum == 0)
					throw StatisticsError("Variable has zero variance in correlation");
				Real correlation = RequireFiniteResult(crossSum / std::sqrt(x.squareSum) / std::sqrt(y.squareSum));
				return std::clamp(correlation, Real(-1), Real(1));
			}

			inline std::vector<CenteredData> CenterColumns(const Matrix<Real>& data) {
				std::vector<CenteredData> columns(data.cols());
				for (int j = 0; j < data.cols(); ++j)
					columns[j] = CenterData(data.rows(), [&](int i) { return data(i, j); });
				return columns;
			}

			inline Real SelectMedian(std::vector<Real>& values) {
				const auto middle = values.begin() + values.size() / 2;
				std::nth_element(values.begin(), middle, values.end());
				if (values.size() % 2 != 0)
					return *middle;
				return std::midpoint(*std::max_element(values.begin(), middle), *middle);
			}
		} // namespace StatisticsDetail

		/// Compute arithmetic mean. Complexity: O(n)
		static Real Avg(const Vector<Real>& data) {
			int n = data.size();

			if (n <= 0)
				throw StatisticsError("Vector size must be greater than 0 in Avg");

			return StatisticsDetail::ArithmeticMean(n, [&](int i) { return data[i]; });
		}

		/// Compute mean and sample variance using denominator n-1. Complexity: O(n)
		static void AvgVar(const Vector<Real>& data, Real& outAvg, Real& outVar) {
			int n = data.size();

			if (n <= 1)
				throw StatisticsError("Vector size must be greater than 1 in AvgVar");

			Real average = Avg(data);
			auto centered = StatisticsDetail::CenterData(n, [&](int i) { return data[i]; });
			Real variance = StatisticsDetail::RestoreUnits(centered.squareSum / (n - 1), 2 * centered.exponent);
			outAvg = average;
			outVar = variance;
		}
		static inline void AvgStdDev(const Vector<Real>& data, Real& outAvg, Real& outStdDev) {
			int n = data.size();
			if (n <= 1)
				throw StatisticsError("Vector size must be greater than 1 in AvgStdDev");
			Real average = Avg(data);
			auto centered = StatisticsDetail::CenterData(n, [&](int i) { return data[i]; });
			Real stddev = StatisticsDetail::RestoreUnits(std::sqrt(centered.squareSum / (n - 1)), centered.exponent);
			outAvg = average;
			outStdDev = stddev;
		}

		// Helper functions for t-tests (aliases for clarity)
		static inline Real Mean(const Vector<Real>& data) { return Avg(data); }

		/// Compute sample variance using Bessel's correction (denominator n-1).
		static Real SampleVariance(const Vector<Real>& data) {
			Real avg, var;
			AvgVar(data, avg, var);
			return var;
		}

		/// Compute population variance using denominator n.
		static inline Real PopulationVariance(const Vector<Real>& data) {
			int n = data.size();
			if (n <= 0)
				throw StatisticsError("Vector size must be greater than 0 in PopulationVariance");

			auto centered = StatisticsDetail::CenterData(n, [&](int i) { return data[i]; });
			return StatisticsDetail::RestoreUnits(centered.squareSum / n, 2 * centered.exponent);
		}

		/// Compute sample standard deviation using denominator n-1.
		static Real SampleStdDev(const Vector<Real>& data) {
			int n = data.size();
			if (n <= 1)
				throw StatisticsError("Vector size must be greater than 1 in SampleStdDev");
			auto centered = StatisticsDetail::CenterData(n, [&](int i) { return data[i]; });
			return StatisticsDetail::RestoreUnits(std::sqrt(centered.squareSum / (n - 1)), centered.exponent);
		}

		/// Compute population standard deviation using denominator n.
		static inline Real PopulationStdDev(const Vector<Real>& data) {
			int n = data.size();
			if (n <= 0)
				throw StatisticsError("Vector size must be greater than 0 in PopulationStdDev");
			auto centered = StatisticsDetail::CenterData(n, [&](int i) { return data[i]; });
			return StatisticsDetail::RestoreUnits(std::sqrt(centered.squareSum / n), centered.exponent);
		}

		/// Backward-compatible alias for SampleVariance.
		static inline Real Variance(const Vector<Real>& data) { return SampleVariance(data); }

		/// Backward-compatible alias for SampleStdDev.
		static inline Real StdDev(const Vector<Real>& data) { return SampleStdDev(data); }

		/// Compute mean, variance, skewness, kurtosis. Complexity: O(n)
		static inline void Moments(const Vector<Real>& data, Real& ave, Real& adev, Real& sdev, Real& var, Real& skew, Real& curt) {
			int n = data.size();
			if (n <= 1)
				throw StatisticsError("Vector size must be greater than 1 in Moments");

			Real average = Avg(data);
			auto centered = StatisticsDetail::CenterData(n, [&](int i) { return data[i]; });
			if (centered.squareSum == 0)
				throw StatisticsError("No skew/kurtosis when variance = 0 (in Moments)");

			StatisticsDetail::CompensatedSum absolute, third, fourth;
			for (int i = 0; i < n; ++i) {
				Real d = centered.Center(data[i]);
				Real square = d * d;
				absolute.Add(std::abs(d));
				third.Add(square * d);
				fourth.Add(square * square);
			}
			Real normalizedVariance = centered.squareSum / (n - 1);
			Real normalizedStdDev = std::sqrt(normalizedVariance);
			Real variance = StatisticsDetail::RestoreUnits(normalizedVariance, 2 * centered.exponent);
			Real stddev = StatisticsDetail::RestoreUnits(normalizedStdDev, centered.exponent);
			Real deviation = StatisticsDetail::RestoreUnits(absolute.Value() / n, centered.exponent);
			Real skewness = StatisticsDetail::RequireFiniteResult(third.Value() / n / normalizedVariance / normalizedStdDev);
			Real kurtosis = StatisticsDetail::RequireFiniteResult(fourth.Value() / n / normalizedVariance / normalizedVariance - 3);
			ave = average;
			adev = deviation;
			sdev = stddev;
			var = variance;
			skew = skewness;
			curt = kurtosis;
		}

		/*********************************************************************/
		/*****                  Order Statistics                         *****/
		/*********************************************************************/

		/// @brief Compute the median of a dataset
		///
		/// The median is the middle value when data is sorted. For even-sized
		/// datasets, returns the average of the two middle values.
		///
		/// @param data Input vector (will be copied and sorted internally)
		/// @return The median value
		/// @throws StatisticsError if data is empty
		///
		/// Complexity: O(n) average using nth_element
		static Real Median(const Vector<Real>& data) {
			int n = data.size();
			if (n <= 0)
				throw StatisticsError("Vector size must be greater than 0 in Median");

			std::vector<Real> v(n);
			for (int i = 0; i < n; i++) {
				StatisticsDetail::RequireFiniteData(data[i]);
				v[i] = data[i];
			}
			return StatisticsDetail::SelectMedian(v);
		}

		/// @brief Compute a percentile of a dataset
		///
		/// Uses linear interpolation between data points (similar to Excel PERCENTILE.INC).
		/// The p-th percentile is the value below which p percent of the data falls.
		///
		/// @param data Input vector (will be copied and sorted internally)
		/// @param p Percentile to compute (must be in [0, 100])
		/// @return The p-th percentile value
		/// @throws StatisticsError if data is empty or p is out of range
		///
		/// Complexity: O(n log n) due to sorting
		static Real Percentile(const Vector<Real>& data, Real p) {
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
			if (p == 0.0)
				return sorted[0];
			if (p == 100.0)
				return sorted[n - 1];

			// Linear interpolation (Excel PERCENTILE.INC method)
			Real rank = (p / MML::Real(100.0)) * (n - 1);
			int lower = static_cast<int>(std::floor(rank));
			int upper = static_cast<int>(std::ceil(rank));

			if (lower == upper)
				return sorted[lower];

			Real fraction = rank - lower;
			return std::lerp(sorted[lower], sorted[upper], fraction);
		}

		/// @brief Compute quartiles (Q1, Q2/Median, Q3) of a dataset
		///
		/// Q1 = 25th percentile (first quartile)
		/// Q2 = 50th percentile (median)
		/// Q3 = 75th percentile (third quartile)
		///
		/// @param data Input vector
		/// @param[out] q1 First quartile (25th percentile)
		/// @param[out] median Second quartile (50th percentile)
		/// @param[out] q3 Third quartile (75th percentile)
		/// @throws StatisticsError if data is empty
		static inline void Quartiles(const Vector<Real>& data, Real& q1, Real& median, Real& q3) {
			int n = data.size();
			if (n <= 0)
				throw StatisticsError("Vector size must be greater than 0 in Quartiles");

			// Sort once instead of calling Percentile three times (each of which sorts)
			std::vector<Real> sorted(n);
			for (int i = 0; i < n; i++)
				sorted[i] = data[i];
			std::sort(sorted.begin(), sorted.end());

			auto percentile_sorted = [&](Real p) -> Real {
				if (p == 0.0) return sorted[0];
				if (p == 100.0) return sorted[n - 1];
				Real rank = (p / MML::Real(100.0)) * (n - 1);
				int lower = static_cast<int>(std::floor(rank));
				int upper = static_cast<int>(std::ceil(rank));
				if (lower == upper) return sorted[lower];
				Real fraction = rank - lower;
				return std::lerp(sorted[lower], sorted[upper], fraction);
			};

			q1 = percentile_sorted(25.0);
			median = percentile_sorted(50.0);
			q3 = percentile_sorted(75.0);
		}

		/// @brief Compute the range of a dataset
		///
		/// Range = max(data) - min(data)
		///
		/// @param data Input vector
		/// @return The range (maximum minus minimum)
		/// @throws StatisticsError if data is empty
		///
		/// Complexity: O(n)
		static inline Real Range(const Vector<Real>& data) {
			int n = data.size();
			if (n <= 0)
				throw StatisticsError("Vector size must be greater than 0 in Range");

			Real minVal = data[0];
			Real maxVal = data[0];
			for (int i = 1; i < n; i++) {
				if (data[i] < minVal)
					minVal = data[i];
				if (data[i] > maxVal)
					maxVal = data[i];
			}
			return maxVal - minVal;
		}

		/// @brief Compute the Interquartile Range (IQR)
		///
		/// IQR = Q3 - Q1 = 75th percentile - 25th percentile
		///
		/// The IQR is a robust measure of spread, less sensitive to outliers
		/// than the range or standard deviation.
		///
		/// @param data Input vector
		/// @return The interquartile range
		/// @throws StatisticsError if data is empty
		static inline Real IQR(const Vector<Real>& data) {
			int n = data.size();
			if (n <= 0)
				throw StatisticsError("Vector size must be greater than 0 in IQR");

			return Percentile(data, 75.0) - Percentile(data, 25.0);
		}

		/// @brief Compute minimum and maximum values
		///
		/// @param data Input vector
		/// @param[out] minVal Minimum value in the data
		/// @param[out] maxVal Maximum value in the data
		/// @throws StatisticsError if data is empty
		///
		/// Complexity: O(n)
		static inline void MinMax(const Vector<Real>& data, Real& minVal, Real& maxVal) {
			int n = data.size();
			if (n <= 0)
				throw StatisticsError("Vector size must be greater than 0 in MinMax");

			minVal = data[0];
			maxVal = data[0];
			for (int i = 1; i < n; i++) {
				if (data[i] < minVal)
					minVal = data[i];
				if (data[i] > maxVal)
					maxVal = data[i];
			}
		}

		/*********************************************************************/
		/*****               Robust Statistics                           *****/
		/*********************************************************************/

		/// @brief Compute the mode (most frequent value) of a dataset
		///
		/// For continuous data, finds the value that appears most frequently.
		/// If multiple values have the same highest frequency, returns the smallest.
		/// For truly continuous data where no values repeat, returns the first value.
		///
		/// @param data Input vector
		/// @return The mode (most frequent value)
		/// @throws StatisticsError if data is empty
		///
		/// Complexity: O(n log n) due to sorting
		///
		/// @note For continuous distributions, consider using kernel density estimation
		/// or histogram-based methods instead.
		static inline Real Mode(const Vector<Real>& data) {
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

		/// @brief Compute the trimmed mean (truncated mean)
		///
		/// The trimmed mean excludes a percentage of the smallest and largest values
		/// before computing the mean. This provides robustness against outliers.
		///
		/// @param data Input vector
		/// @param trimPercent Percentage to trim from each end (0-50)
		/// e.g., 10 means remove bottom 10% and top 10%
		/// @return The trimmed mean
		/// @throws StatisticsError if data is empty or trimPercent is invalid
		///
		/// Complexity: O(n log n) due to sorting
		static inline Real TrimmedMean(const Vector<Real>& data, Real trimPercent) {
			int n = data.size();
			if (n <= 0)
				throw StatisticsError("Vector size must be greater than 0 in TrimmedMean");
			if (!(trimPercent >= 0.0 && trimPercent < 50.0))
				throw StatisticsError("trimPercent must be in range [0, 50) in TrimmedMean");

			// Create a sorted copy
			std::vector<Real> sorted(n);
			for (int i = 0; i < n; i++) {
				StatisticsDetail::RequireFiniteData(data[i]);
				sorted[i] = data[i];
			}
			std::sort(sorted.begin(), sorted.end());

			// Calculate how many to trim from each end
			int trimCount = static_cast<int>(std::floor(n * trimPercent / 100.0));

			// Calculate mean of remaining values
			int count = n - 2 * trimCount;
			if (count == 0)
				throw StatisticsError("No values remain after trimming in TrimmedMean");

			return StatisticsDetail::ArithmeticMean(count, [&](int i) { return sorted[trimCount + i]; });
		}

		/// @brief Compute the Median Absolute Deviation (MAD)
		///
		/// MAD = median(|x_i - median(x)|)
		///
		/// MAD is a robust measure of variability. It is less sensitive to
		/// outliers than standard deviation. To estimate standard deviation
		/// for normal data, multiply MAD by 1.4826.
		///
		/// @param data Input vector
		/// @return The median absolute deviation
		/// @throws StatisticsError if data is empty
		///
		/// Complexity: O(n log n) due to sorting
		static inline Real MAD(const Vector<Real>& data) {
			int n = data.size();
			if (n <= 0)
				throw StatisticsError("Vector size must be greater than 0 in MAD");

			Real med = Median(data);

			// Compute absolute deviations
			std::vector<Real> absDeviations(n);
			for (int i = 0; i < n; i++)
				absDeviations[i] = std::abs(data[i] - med);

			// Overflowing outlier deviations can be ordered and discarded without invalidating a finite MAD.
			return StatisticsDetail::RequireFiniteResult(StatisticsDetail::SelectMedian(absDeviations));
		}

		/// @brief Compute the geometric mean
		///
		/// GeometricMean = (x_1 * x_2 * ... * x_n)^(1/n) = exp(mean(log(x_i)))
		///
		/// The geometric mean is appropriate for data that is multiplicative
		/// in nature (e.g., growth rates, ratios). All values must be positive.
		///
		/// @param data Input vector (all values must be positive)
		/// @return The geometric mean
		/// @throws StatisticsError if data is empty or contains non-positive values
		///
		/// Complexity: O(n)
		static inline Real GeometricMean(const Vector<Real>& data) {
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

		/// @brief Compute the harmonic mean
		///
		/// HarmonicMean = n / (1/x_1 + 1/x_2 + ... + 1/x_n)
		///
		/// The harmonic mean is appropriate for averaging rates or ratios
		/// (e.g., speeds, P/E ratios). All values must be positive.
		///
		/// @param data Input vector (all values must be positive)
		/// @return The harmonic mean
		/// @throws StatisticsError if data is empty or contains non-positive values
		///
		/// Complexity: O(n)
		static inline Real HarmonicMean(const Vector<Real>& data) {
			int n = data.size();
			if (n <= 0)
				throw StatisticsError("Vector size must be greater than 0 in HarmonicMean");

			Real reciprocalSum = 0.0;
			for (int i = 0; i < n; i++) {
				if (data[i] <= 0.0)
					throw StatisticsError("All values must be positive in HarmonicMean");
				reciprocalSum += MML::Real(1.0) / data[i];
			}

			return n / reciprocalSum;
		}

		/// @brief Compute the weighted mean
		///
		/// WeightedMean = sum(w_i * x_i) / sum(w_i)
		///
		/// @param data Input vector of values
		/// @param weights Vector of weights (must have same size as data)
		/// @return The weighted mean
		/// @throws StatisticsError if data is empty, sizes don't match,
		/// or total weight is zero/negative
		///
		/// Complexity: O(n)
		static Real WeightedMean(const Vector<Real>& data, const Vector<Real>& weights) {
			int n = data.size();
			if (n <= 0)
				throw StatisticsError("Vector size must be greater than 0 in WeightedMean");
			if (weights.size() != n)
				throw StatisticsError("Data and weights must have the same size in WeightedMean");

			Real weightedSum = 0.0;
			Real weightedCorrection = 0.0;
			Real totalWeight = 0.0;
			Real weightCorrection = 0.0;
			for (int i = 0; i < n; i++) {
				if (!std::isfinite(data[i]))
					throw StatisticsError("Data must be finite in WeightedMean");
				if (!std::isfinite(weights[i]) || weights[i] < 0.0)
					throw StatisticsError("Weights must be finite and non-negative in WeightedMean");

				Real weightedValue = weights[i] * data[i];
				if (!std::isfinite(weightedValue))
					throw StatisticsError("Weighted value overflow in WeightedMean");
				StatisticsDetail::CompensatedAdd(weightedValue, weightedSum, weightedCorrection);
				StatisticsDetail::CompensatedAdd(weights[i], totalWeight, weightCorrection);
			}

			totalWeight += weightCorrection;
			weightedSum += weightedCorrection;
			if (totalWeight <= 0.0)
				throw StatisticsError("Total weight must be positive in WeightedMean");
			if (!std::isfinite(totalWeight) || !std::isfinite(weightedSum))
				throw StatisticsError("Weighted accumulation overflow in WeightedMean");

			return weightedSum / totalWeight;
		}

		/// @brief Compute the weighted variance
		///
		/// Uses frequency weights, equivalent to repeating each value w_i times:
		/// WeightedVar = sum(w_i * (x_i - weighted_mean)^2) / (sum(w_i) - 1)
		///
		/// @param data Input vector of values
		/// @param weights Vector of weights (must have same size as data)
		/// @return The weighted variance
		/// @throws StatisticsError if data has fewer than 2 elements, sizes don't match,
		/// or weights are invalid
		///
		/// Complexity: O(n)
		static Real WeightedVariance(const Vector<Real>& data, const Vector<Real>& weights) {
			int n = data.size();
			if (n <= 0)
				throw StatisticsError("Vector size must be greater than 0 in WeightedVariance");
			if (weights.size() != n)
				throw StatisticsError("Data and weights must have the same size in WeightedVariance");

			Real wMean = WeightedMean(data, weights);

			Real weightedSqSum = 0.0;
			Real squaredCorrection = 0.0;
			Real totalWeight = 0.0;
			Real weightCorrection = 0.0;
			for (int i = 0; i < n; i++) {
				Real diff = data[i] - wMean;
				Real weightedSquare = weights[i] * diff * diff;
				if (!std::isfinite(weightedSquare))
					throw StatisticsError("Weighted squared deviation overflow in WeightedVariance");
				StatisticsDetail::CompensatedAdd(weightedSquare, weightedSqSum, squaredCorrection);
				StatisticsDetail::CompensatedAdd(weights[i], totalWeight, weightCorrection);
			}

			totalWeight += weightCorrection;
			weightedSqSum += squaredCorrection;
			if (totalWeight <= 1.0)
				throw StatisticsError("Total weight must be greater than 1 in WeightedVariance");
			if (!std::isfinite(totalWeight) || !std::isfinite(weightedSqSum))
				throw StatisticsError("Weighted accumulation overflow in WeightedVariance");

			return weightedSqSum / (totalWeight - MML::Real(1.0));
		}

		/// @brief Compute frequency-weighted sample standard deviation.
		static inline Real WeightedStdDev(const Vector<Real>& data, const Vector<Real>& weights) {
			return std::sqrt(WeightedVariance(data, weights));
		}

		/*********************************************************************/
		/*****               Correlation & Covariance                    *****/
		/*********************************************************************/

		/// @brief Compute the sample covariance between two variables
		///
		/// Cov(X,Y) = sum((x_i - mean_x)(y_i - mean_y)) / (n - 1)
		///
		/// Uses a two-pass algorithm for numerical stability.
		///
		/// @param x First variable (vector of values)
		/// @param y Second variable (must have same size as x)
		/// @return Sample covariance
		/// @throws StatisticsError if vectors are empty or have different sizes
		///
		/// Complexity: O(n)
		static inline Real Covariance(const Vector<Real>& x, const Vector<Real>& y) {
			int n = x.size();
			if (n <= 1)
				throw StatisticsError("Vector size must be greater than 1 in Covariance");
			if (y.size() != n)
				throw StatisticsError("Vectors must have the same size in Covariance");

			auto xAt = [&](int i) { return x[i]; };
			auto yAt = [&](int i) { return y[i]; };
			auto centeredX = StatisticsDetail::CenterData(n, xAt);
			auto centeredY = StatisticsDetail::CenterData(n, yAt);
			Real crossSum = StatisticsDetail::CrossSum(n, xAt, yAt, centeredX, centeredY);
			return StatisticsDetail::RestoreUnits(crossSum / (n - 1), centeredX.exponent + centeredY.exponent);
		}

		/// @brief Compute frequency-weighted sample covariance.
		///
		/// Equivalent to repeating each (x_i, y_i) pair w_i times and computing
		/// sample covariance with denominator sum(weights) - 1.
		static Real WeightedCovariance(const Vector<Real>& x, const Vector<Real>& y,
		                               const Vector<Real>& weights) {
			int n = x.size();
			if (n <= 0)
				throw StatisticsError("Vector size must be greater than 0 in WeightedCovariance");
			if (y.size() != n || weights.size() != n)
				throw StatisticsError("Vectors and weights must have the same size in WeightedCovariance");

			Real meanX = WeightedMean(x, weights);
			Real meanY = WeightedMean(y, weights);
			Real crossSum = 0.0;
			Real crossCorrection = 0.0;
			Real totalWeight = 0.0;
			Real weightCorrection = 0.0;
			for (int i = 0; i < n; i++) {
				Real weightedProduct = weights[i] * (x[i] - meanX) * (y[i] - meanY);
				if (!std::isfinite(weightedProduct))
					throw StatisticsError("Weighted cross-product overflow in WeightedCovariance");
				StatisticsDetail::CompensatedAdd(weightedProduct, crossSum, crossCorrection);
				StatisticsDetail::CompensatedAdd(weights[i], totalWeight, weightCorrection);
			}

			totalWeight += weightCorrection;
			crossSum += crossCorrection;
			if (totalWeight <= 1.0)
				throw StatisticsError("Total weight must be greater than 1 in WeightedCovariance");
			if (!std::isfinite(totalWeight) || !std::isfinite(crossSum))
				throw StatisticsError("Weighted accumulation overflow in WeightedCovariance");

			return crossSum / (totalWeight - MML::Real(1.0));
		}

		/// @brief Compute the Pearson correlation coefficient
		///
		/// r = Cov(X,Y) / (StdDev(X) * StdDev(Y))
		///
		/// The Pearson correlation measures linear relationship between variables.
		/// Values range from -1 (perfect negative) to +1 (perfect positive).
		/// A value of 0 indicates no linear correlation.
		///
		/// @param x First variable (vector of values)
		/// @param y Second variable (must have same size as x)
		/// @return Pearson correlation coefficient r in [-1, 1]
		/// @throws StatisticsError if vectors are empty, have different sizes,
		/// or either variable has zero variance
		///
		/// Complexity: O(n)
		static Real PearsonCorrelation(const Vector<Real>& x, const Vector<Real>& y) {
			int n = x.size();
			if (n <= 1)
				throw StatisticsError("Vector size must be greater than 1 in PearsonCorrelation");
			if (y.size() != n)
				throw StatisticsError("Vectors must have the same size in PearsonCorrelation");

			auto xAt = [&](int i) { return x[i]; };
			auto yAt = [&](int i) { return y[i]; };
			auto centeredX = StatisticsDetail::CenterData(n, xAt);
			auto centeredY = StatisticsDetail::CenterData(n, yAt);
			return StatisticsDetail::Correlation(
				StatisticsDetail::CrossSum(n, xAt, yAt, centeredX, centeredY), centeredX, centeredY);
		}

		/// @brief Compute frequency-weighted Pearson correlation.
		static inline Real WeightedPearsonCorrelation(const Vector<Real>& x, const Vector<Real>& y,
		                                       const Vector<Real>& weights) {
			Real covariance = WeightedCovariance(x, y, weights);
			Real varianceX = WeightedVariance(x, weights);
			Real varianceY = WeightedVariance(y, weights);
			if (varianceX <= 0.0)
				throw StatisticsError("First variable has zero variance in WeightedPearsonCorrelation");
			if (varianceY <= 0.0)
				throw StatisticsError("Second variable has zero variance in WeightedPearsonCorrelation");

			Real correlation = covariance / std::sqrt(varianceX) / std::sqrt(varianceY);
			if (!std::isfinite(correlation))
				throw StatisticsError("Non-finite result in WeightedPearsonCorrelation");

			return std::clamp(correlation, static_cast<MML::Real>(-1.0), static_cast<MML::Real>(1.0));
		}

		/// @brief Result structure for correlation with significance test
		struct CorrelationResult {
			Real r;								///< Correlation coefficient
			Real tStatistic;			///< t-statistic for significance test
			int degreesOfFreedom; ///< Degrees of freedom (n-2)

			CorrelationResult(Real correlation, Real tStat, int df)
					: r(correlation)
					, tStatistic(tStat)
					, degreesOfFreedom(df) {}
		};

		/// @brief Compute Pearson correlation with t-test for significance
		///
		/// Tests the null hypothesis H0: ρ = 0 (no population correlation)
		/// The t-statistic follows a t-distribution with n-2 degrees of freedom.
		///
		/// t = r * sqrt((n-2) / (1-r²))
		///
		/// @param x First variable
		/// @param y Second variable
		/// @return CorrelationResult containing r, t-statistic, and degrees of freedom
		/// @throws StatisticsError if inputs are invalid
		///
		/// @note Use t-distribution critical values or p-value calculation to assess significance
		static inline CorrelationResult PearsonCorrelationWithTest(const Vector<Real>& x, const Vector<Real>& y) {
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
				tStat = (r > 0) ? std::numeric_limits<Real>::infinity() : -std::numeric_limits<Real>::infinity();
			} else {
				tStat = r * std::sqrt(static_cast<Real>(df) / (MML::Real(1.0) - r2));
			}

			return CorrelationResult(r, tStat, df);
		}

		/// @brief Compute the coefficient of determination (R²)
		///
		/// R² = r² where r is the Pearson correlation coefficient.
		/// Represents the proportion of variance in y explained by x.
		///
		/// @param x Independent variable
		/// @param y Dependent variable
		/// @return R² value in [0, 1]
		static inline Real RSquared(const Vector<Real>& x, const Vector<Real>& y) {
			Real r = PearsonCorrelation(x, y);
			return r * r;
		}

		/// @brief Compute the covariance matrix for multivariate data
		///
		/// For a data matrix where each column is a variable and each row is
		/// an observation, computes the pairwise covariance between all variables.
		///
		/// @param data Matrix of size (n_observations x n_variables)
		/// @return Square covariance matrix of size (n_variables x n_variables)
		/// @throws StatisticsError if matrix has insufficient observations
		///
		/// @note The diagonal contains variances, off-diagonal contains covariances
		static inline Matrix<Real> CovarianceMatrix(const Matrix<Real>& data) {
			int n = data.rows(); // observations
			int p = data.cols(); // variables

			if (n <= 1)
				throw StatisticsError("Need at least 2 observations for CovarianceMatrix");
			if (p <= 0)
				throw StatisticsError("Need at least 1 variable for CovarianceMatrix");

			auto columns = StatisticsDetail::CenterColumns(data);
			Matrix<Real> covMatrix(p, p);
			for (int j1 = 0; j1 < p; j1++) {
				for (int j2 = j1; j2 < p; j2++) {
					Real crossSum = j1 == j2 ? columns[j1].squareSum : StatisticsDetail::CrossSum(n,
						[&](int i) { return data(i, j1); }, [&](int i) { return data(i, j2); }, columns[j1], columns[j2]);
					Real cov = StatisticsDetail::RestoreUnits(crossSum / (n - 1),
						columns[j1].exponent + columns[j2].exponent);
					covMatrix(j1, j2) = cov;
					covMatrix(j2, j1) = cov; // Symmetric
				}
			}

			return covMatrix;
		}

		/// @brief Compute the correlation matrix for multivariate data
		///
		/// Similar to covariance matrix, but normalized to [-1, 1] range.
		/// The diagonal is always 1 (variable correlated with itself).
		///
		/// @param data Matrix of size (n_observations x n_variables)
		/// @return Square correlation matrix of size (n_variables x n_variables)
		/// @throws StatisticsError if matrix has insufficient observations or zero variance
		static inline Matrix<Real> CorrelationMatrix(const Matrix<Real>& data) {
			int n = data.rows();
			int p = data.cols();

			if (n <= 1)
				throw StatisticsError("Need at least 2 observations for CorrelationMatrix");
			if (p <= 0)
				throw StatisticsError("Need at least 1 variable for CorrelationMatrix");

			auto columns = StatisticsDetail::CenterColumns(data);
			for (int j = 0; j < p; j++) {
				if (columns[j].squareSum == 0)
					throw StatisticsError("Variable has zero variance in CorrelationMatrix");
			}

			// Normalize upper triangle, then mirror to preserve exact symmetry.
			Matrix<Real> corrMatrix(p, p);
			for (int j = 0; j < p; j++)
				corrMatrix(j, j) = 1.0;
			for (int j1 = 0; j1 < p; j1++) {
				for (int j2 = j1 + 1; j2 < p; j2++) {
					Real crossSum = StatisticsDetail::CrossSum(n,
						[&](int i) { return data(i, j1); }, [&](int i) { return data(i, j2); }, columns[j1], columns[j2]);
					Real correlation = StatisticsDetail::Correlation(crossSum, columns[j1], columns[j2]);
					corrMatrix(j1, j2) = correlation;
					corrMatrix(j2, j1) = correlation;
				}
			}

			return corrMatrix;
		}
	} // namespace Statistics
} // namespace MML
#endif

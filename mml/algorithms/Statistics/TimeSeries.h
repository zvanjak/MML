///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        TimeSeries.h                                                        ///
///  Description: Time series analysis functions including moving averages,           ///
///               rolling statistics, autocorrelation, and differencing operations    ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_TIME_SERIES_H
#define MML_TIME_SERIES_H

#include "MMLBase.h"
#include "MMLExceptions.h"
#include "base/Vector.h"
#include "base/Matrix.h"

#include <deque>
#include <algorithm>
#include <cmath>

namespace MML {
	namespace Statistics {
		namespace TimeSeries {
			/////////////////////////////////////////////////////////////////////////////////////
			///                           MOVING AVERAGES                                     ///
			/////////////////////////////////////////////////////////////////////////////////////

			/// @brief Simple Moving Average (SMA)
			/// @details Computes the unweighted mean of the previous windowSize data points.
			///          Formula: SMA_t = (1/windowSize) * sum(x_{t-windowSize+1} to x_t)
			/// @param data Input time series data
			/// @param windowSize Number of points in the moving window
			/// @return Vector of size (n - windowSize + 1) containing smoothed values
			/// @throws StatisticsError if windowSize > data.size() or windowSize < 1 or data is empty
			static Vector<Real> SimpleMovingAverage(const Vector<Real>& data, int windowSize) {
				int n = data.size();

				if (n == 0)
					throw StatisticsError("SimpleMovingAverage: Input data cannot be empty");
				if (windowSize < 1)
					throw StatisticsError("SimpleMovingAverage: Window size must be at least 1");
				if (windowSize > n)
					throw StatisticsError("SimpleMovingAverage: Window size cannot exceed data length");

				int resultSize = n - windowSize + 1;
				Vector<Real> result(resultSize);

				// Compute first window sum
				Real windowSum = 0.0;
				for (int i = 0; i < windowSize; ++i)
					windowSum += data[i];

				result[0] = windowSum / windowSize;

				// Slide the window using running sum (O(n) complexity)
				for (int i = 1; i < resultSize; ++i) {
					windowSum = windowSum - data[i - 1] + data[i + windowSize - 1];
					result[i] = windowSum / windowSize;
				}

				return result;
			}

			/// @brief Exponential Moving Average (EMA) with smoothing factor alpha
			/// @details Applies exponential weighting where recent values have more influence.
			///          Formula: EMA_t = alpha * x_t + (1 - alpha) * EMA_{t-1}
			///          Initial EMA_0 = x_0
			/// @param data Input time series data
			/// @param alpha Smoothing factor in range (0, 1). Higher alpha = more weight on recent values.
			/// @return Vector of same size as input containing EMA values
			/// @throws StatisticsError if alpha not in (0,1) or data is empty
			static Vector<Real> ExponentialMovingAverage(const Vector<Real>& data, Real alpha) {
				int n = data.size();

				if (n == 0)
					throw StatisticsError("ExponentialMovingAverage: Input data cannot be empty");
				if (alpha <= 0.0 || alpha >= 1.0)
					throw StatisticsError("ExponentialMovingAverage: Alpha must be in range (0, 1)");

				Vector<Real> result(n);
				result[0] = data[0]; // Initialize with first value

				Real oneMinusAlpha = 1.0 - alpha;
				for (int i = 1; i < n; ++i)
					result[i] = alpha * data[i] + oneMinusAlpha * result[i - 1];

				return result;
			}

			/// @brief Exponential Moving Average (EMA) with span parameter
			/// @details Computes alpha from span: alpha = 2 / (span + 1)
			///          This is the common convention used in pandas.
			/// @param data Input time series data
			/// @param span The span parameter (must be >= 1)
			/// @return Vector of same size as input containing EMA values
			/// @throws StatisticsError if span < 1 or data is empty
			static Vector<Real> ExponentialMovingAverageSpan(const Vector<Real>& data, int span) {
				if (span < 1)
					throw StatisticsError("ExponentialMovingAverageSpan: Span must be at least 1");

				Real alpha = 2.0 / (span + 1.0);
				return ExponentialMovingAverage(data, alpha);
			}

			/// @brief Weighted Moving Average (WMA) with custom weights
			/// @details Computes weighted average over sliding window using provided weights.
			///          Formula: WMA_t = sum(w_i * x_{t-n+1+i}) / sum(w_i)
			/// @param data Input time series data
			/// @param weights Weight vector (window size = weights.size())
			/// @return Vector of size (n - windowSize + 1) containing weighted averages
			/// @throws StatisticsError if weights empty, window > data length, or all weights zero
			static Vector<Real> WeightedMovingAverage(const Vector<Real>& data, const Vector<Real>& weights) {
				int n = data.size();
				int windowSize = weights.size();

				if (n == 0)
					throw StatisticsError("WeightedMovingAverage: Input data cannot be empty");
				if (windowSize == 0)
					throw StatisticsError("WeightedMovingAverage: Weights cannot be empty");
				if (windowSize > n)
					throw StatisticsError("WeightedMovingAverage: Window size cannot exceed data length");

				// Compute weight sum once
				Real weightSum = 0.0;
				for (int i = 0; i < windowSize; ++i)
					weightSum += weights[i];

				if (std::abs(weightSum) < Constants::Eps)
					throw StatisticsError("WeightedMovingAverage: Sum of weights cannot be zero");

				int resultSize = n - windowSize + 1;
				Vector<Real> result(resultSize);

				for (int i = 0; i < resultSize; ++i) {
					Real weighted = 0.0;
					for (int j = 0; j < windowSize; ++j)
						weighted += weights[j] * data[i + j];
					result[i] = weighted / weightSum;
				}

				return result;
			}

			/// @brief Weighted Moving Average (WMA) with linear weights
			/// @details Uses linearly increasing weights: [1, 2, 3, ..., windowSize]
			///          Most recent value gets highest weight.
			/// @param data Input time series data
			/// @param windowSize Number of points in the moving window
			/// @return Vector of size (n - windowSize + 1) containing weighted averages
			static Vector<Real> WeightedMovingAverageLinear(const Vector<Real>& data, int windowSize) {
				if (windowSize < 1)
					throw StatisticsError("WeightedMovingAverageLinear: Window size must be at least 1");

				// Create linear weights [1, 2, 3, ..., windowSize]
				Vector<Real> weights(windowSize);
				for (int i = 0; i < windowSize; ++i)
					weights[i] = static_cast<Real>(i + 1);

				return WeightedMovingAverage(data, weights);
			}


			/////////////////////////////////////////////////////////////////////////////////////
			///                          ROLLING STATISTICS                                   ///
			/////////////////////////////////////////////////////////////////////////////////////

			/// @brief Rolling Mean (equivalent to SMA)
			/// @param data Input time series data
			/// @param windowSize Number of points in the rolling window
			/// @return Vector of size (n - windowSize + 1) containing rolling means
			static Vector<Real> RollingMean(const Vector<Real>& data, int windowSize) { return SimpleMovingAverage(data, windowSize); }

			/// @brief Rolling Sum over a sliding window
			/// @param data Input time series data
			/// @param windowSize Number of points in the rolling window
			/// @return Vector of size (n - windowSize + 1) containing rolling sums
			static Vector<Real> RollingSum(const Vector<Real>& data, int windowSize) {
				int n = data.size();

				if (n == 0)
					throw StatisticsError("RollingSum: Input data cannot be empty");
				if (windowSize < 1)
					throw StatisticsError("RollingSum: Window size must be at least 1");
				if (windowSize > n)
					throw StatisticsError("RollingSum: Window size cannot exceed data length");

				int resultSize = n - windowSize + 1;
				Vector<Real> result(resultSize);

				// Compute first window sum
				Real windowSum = 0.0;
				for (int i = 0; i < windowSize; ++i)
					windowSum += data[i];

				result[0] = windowSum;

				// Slide using running sum
				for (int i = 1; i < resultSize; ++i) {
					windowSum = windowSum - data[i - 1] + data[i + windowSize - 1];
					result[i] = windowSum;
				}

				return result;
			}

			/// @brief Rolling Variance using Welford's online algorithm for numerical stability
			/// @details Computes sample variance (n-1 denominator) over sliding window.
			/// @param data Input time series data
			/// @param windowSize Number of points in the rolling window (must be >= 2)
			/// @return Vector of size (n - windowSize + 1) containing rolling variances
			static Vector<Real> RollingVariance(const Vector<Real>& data, int windowSize) {
				int n = data.size();

				if (n == 0)
					throw StatisticsError("RollingVariance: Input data cannot be empty");
				if (windowSize < 2)
					throw StatisticsError("RollingVariance: Window size must be at least 2 for variance");
				if (windowSize > n)
					throw StatisticsError("RollingVariance: Window size cannot exceed data length");

				int resultSize = n - windowSize + 1;
				Vector<Real> result(resultSize);

				// For each window, compute variance directly
				// Using compensated summation for numerical stability
				for (int i = 0; i < resultSize; ++i) {
					// Compute mean
					Real mean = 0.0;
					for (int j = 0; j < windowSize; ++j)
						mean += data[i + j];
					mean /= windowSize;

					// Compute variance with compensation (Welford-like correction)
					Real sumSq = 0.0;
					Real compensation = 0.0;
					for (int j = 0; j < windowSize; ++j) {
						Real diff = data[i + j] - mean;
						sumSq += diff * diff;
						compensation += diff;
					}
					// Correction term improves numerical stability
					result[i] = (sumSq - compensation * compensation / windowSize) / (windowSize - 1);
				}

				return result;
			}

			/// @brief Rolling Standard Deviation
			/// @param data Input time series data
			/// @param windowSize Number of points in the rolling window (must be >= 2)
			/// @return Vector of size (n - windowSize + 1) containing rolling standard deviations
			static Vector<Real> RollingStdDev(const Vector<Real>& data, int windowSize) {
				Vector<Real> variance = RollingVariance(data, windowSize);

				for (int i = 0; i < variance.size(); ++i)
					variance[i] = std::sqrt(variance[i]);

				return variance;
			}

			/// @brief Rolling Minimum using deque-based O(n) algorithm
			/// @param data Input time series data
			/// @param windowSize Number of points in the rolling window
			/// @return Vector of size (n - windowSize + 1) containing rolling minimums
			static Vector<Real> RollingMin(const Vector<Real>& data, int windowSize) {
				int n = data.size();

				if (n == 0)
					throw StatisticsError("RollingMin: Input data cannot be empty");
				if (windowSize < 1)
					throw StatisticsError("RollingMin: Window size must be at least 1");
				if (windowSize > n)
					throw StatisticsError("RollingMin: Window size cannot exceed data length");

				int resultSize = n - windowSize + 1;
				Vector<Real> result(resultSize);

				// Monotonic deque storing indices of potential minimums
				std::deque<int> dq;

				for (int i = 0; i < n; ++i) {
					// Remove elements outside current window
					while (!dq.empty() && dq.front() <= i - windowSize)
						dq.pop_front();

					// Remove elements larger than current (they can never be minimum)
					while (!dq.empty() && data[dq.back()] >= data[i])
						dq.pop_back();

					dq.push_back(i);

					// Start recording results once we have a full window
					if (i >= windowSize - 1)
						result[i - windowSize + 1] = data[dq.front()];
				}

				return result;
			}

			/// @brief Rolling Maximum using deque-based O(n) algorithm
			/// @param data Input time series data
			/// @param windowSize Number of points in the rolling window
			/// @return Vector of size (n - windowSize + 1) containing rolling maximums
			static Vector<Real> RollingMax(const Vector<Real>& data, int windowSize) {
				int n = data.size();

				if (n == 0)
					throw StatisticsError("RollingMax: Input data cannot be empty");
				if (windowSize < 1)
					throw StatisticsError("RollingMax: Window size must be at least 1");
				if (windowSize > n)
					throw StatisticsError("RollingMax: Window size cannot exceed data length");

				int resultSize = n - windowSize + 1;
				Vector<Real> result(resultSize);

				// Monotonic deque storing indices of potential maximums
				std::deque<int> dq;

				for (int i = 0; i < n; ++i) {
					// Remove elements outside current window
					while (!dq.empty() && dq.front() <= i - windowSize)
						dq.pop_front();

					// Remove elements smaller than current (they can never be maximum)
					while (!dq.empty() && data[dq.back()] <= data[i])
						dq.pop_back();

					dq.push_back(i);

					// Start recording results once we have a full window
					if (i >= windowSize - 1)
						result[i - windowSize + 1] = data[dq.front()];
				}

				return result;
			}


			/////////////////////////////////////////////////////////////////////////////////////
			///                         AUTOCORRELATION                                       ///
			/////////////////////////////////////////////////////////////////////////////////////

			/// @brief Single lag autocorrelation coefficient
			/// @details Measures correlation between x_t and x_{t-lag}
			///          Formula: ACF(k) = Cov(x_t, x_{t-k}) / Var(x_t)
			/// @param data Input time series data
			/// @param lag The lag value (must be >= 0 and < data.size())
			/// @return Autocorrelation coefficient in [-1, 1], ACF(0) = 1.0
			static Real Autocorrelation(const Vector<Real>& data, int lag) {
				int n = data.size();

				if (n == 0)
					throw StatisticsError("Autocorrelation: Input data cannot be empty");
				if (lag < 0)
					throw StatisticsError("Autocorrelation: Lag must be non-negative");
				if (lag >= n)
					throw StatisticsError("Autocorrelation: Lag must be less than data length");

				if (lag == 0)
					return 1.0;

				// Compute mean
				Real mean = 0.0;
				for (int i = 0; i < n; ++i)
					mean += data[i];
				mean /= n;

				// Compute variance and covariance at lag
				Real variance = 0.0;
				Real covariance = 0.0;

				for (int i = 0; i < n; ++i) {
					Real diff = data[i] - mean;
					variance += diff * diff;

					if (i >= lag)
						covariance += diff * (data[i - lag] - mean);
				}

				if (std::abs(variance) < Constants::Eps)
					return 0.0; // Constant series has no autocorrelation

				return covariance / variance;
			}

			/// @brief Autocorrelation Function (ACF) for multiple lags
			/// @details Computes ACF for lags 0, 1, 2, ..., maxLag
			/// @param data Input time series data
			/// @param maxLag Maximum lag to compute (inclusive)
			/// @return Vector of size (maxLag + 1) with ACF values
			static Vector<Real> AutocorrelationFunction(const Vector<Real>& data, int maxLag) {
				int n = data.size();

				if (n == 0)
					throw StatisticsError("AutocorrelationFunction: Input data cannot be empty");
				if (maxLag < 0)
					throw StatisticsError("AutocorrelationFunction: maxLag must be non-negative");
				if (maxLag >= n)
					throw StatisticsError("AutocorrelationFunction: maxLag must be less than data length");

				Vector<Real> acf(maxLag + 1);

				// Compute mean once
				Real mean = 0.0;
				for (int i = 0; i < n; ++i)
					mean += data[i];
				mean /= n;

				// Compute variance once
				Real variance = 0.0;
				for (int i = 0; i < n; ++i) {
					Real diff = data[i] - mean;
					variance += diff * diff;
				}

				acf[0] = 1.0; // ACF(0) is always 1

				if (std::abs(variance) < Constants::Eps) {
					// Constant series: ACF is 0 for all lags > 0
					for (int k = 1; k <= maxLag; ++k)
						acf[k] = 0.0;
					return acf;
				}

				// Compute ACF for each lag
				for (int k = 1; k <= maxLag; ++k) {
					Real covariance = 0.0;
					for (int i = k; i < n; ++i)
						covariance += (data[i] - mean) * (data[i - k] - mean);

					acf[k] = covariance / variance;
				}

				return acf;
			}

			/// @brief Partial Autocorrelation using Durbin-Levinson recursion
			/// @details Measures direct correlation between x_t and x_{t-k}, removing
			///          indirect effects through intermediate lags.
			/// @param data Input time series data
			/// @param lag The lag value (must be >= 0 and < data.size())
			/// @return Partial autocorrelation coefficient at given lag
			static Real PartialAutocorrelation(const Vector<Real>& data, int lag) {
				int n = data.size();

				if (n == 0)
					throw StatisticsError("PartialAutocorrelation: Input data cannot be empty");
				if (lag < 0)
					throw StatisticsError("PartialAutocorrelation: Lag must be non-negative");
				if (lag >= n)
					throw StatisticsError("PartialAutocorrelation: Lag must be less than data length");

				if (lag == 0)
					return 1.0;

				// Get ACF values we need
				Vector<Real> acf = AutocorrelationFunction(data, lag);

				// Durbin-Levinson algorithm
				Vector<Real> phi(lag + 1);
				Vector<Real> phiPrev(lag + 1);

				phi[1] = acf[1];

				for (int k = 2; k <= lag; ++k) {
					// Save previous phi values
					for (int j = 1; j < k; ++j)
						phiPrev[j] = phi[j];

					// Compute numerator and denominator
					Real num = acf[k];
					Real den = 1.0;

					for (int j = 1; j < k; ++j) {
						num -= phiPrev[j] * acf[k - j];
						den -= phiPrev[j] * acf[j];
					}

					if (std::abs(den) < Constants::Eps)
						phi[k] = 0.0;
					else
						phi[k] = num / den;

					// Update phi values
					for (int j = 1; j < k; ++j)
						phi[j] = phiPrev[j] - phi[k] * phiPrev[k - j];
				}

				return phi[lag];
			}

			/// @brief Partial Autocorrelation Function (PACF) for multiple lags
			/// @details Computes PACF for lags 0, 1, 2, ..., maxLag using Durbin-Levinson
			/// @param data Input time series data
			/// @param maxLag Maximum lag to compute (inclusive)
			/// @return Vector of size (maxLag + 1) with PACF values
			static Vector<Real> PartialAutocorrelationFunction(const Vector<Real>& data, int maxLag) {
				int n = data.size();

				if (n == 0)
					throw StatisticsError("PartialAutocorrelationFunction: Input data cannot be empty");
				if (maxLag < 0)
					throw StatisticsError("PartialAutocorrelationFunction: maxLag must be non-negative");
				if (maxLag >= n)
					throw StatisticsError("PartialAutocorrelationFunction: maxLag must be less than data length");

				Vector<Real> pacf(maxLag + 1);
				pacf[0] = 1.0;

				if (maxLag == 0)
					return pacf;

				// Get ACF values
				Vector<Real> acf = AutocorrelationFunction(data, maxLag);

				// Durbin-Levinson algorithm - compute all PACF values efficiently
				Vector<Real> phi(maxLag + 1);
				Vector<Real> phiPrev(maxLag + 1);

				phi[1] = acf[1];
				pacf[1] = phi[1];

				for (int k = 2; k <= maxLag; ++k) {
					// Save previous phi values
					for (int j = 1; j < k; ++j)
						phiPrev[j] = phi[j];

					// Compute phi[k]
					Real num = acf[k];
					Real den = 1.0;

					for (int j = 1; j < k; ++j) {
						num -= phiPrev[j] * acf[k - j];
						den -= phiPrev[j] * acf[j];
					}

					if (std::abs(den) < Constants::Eps)
						phi[k] = 0.0;
					else
						phi[k] = num / den;

					pacf[k] = phi[k];

					// Update phi values for next iteration
					for (int j = 1; j < k; ++j)
						phi[j] = phiPrev[j] - phi[k] * phiPrev[k - j];
				}

				return pacf;
			}


			/////////////////////////////////////////////////////////////////////////////////////
			///                      DIFFERENCING AND LAG OPERATIONS                          ///
			/////////////////////////////////////////////////////////////////////////////////////

			/// @brief Create lagged version of time series
			/// @details Returns data shifted by k positions: {x_k, x_{k+1}, ..., x_{n-1}}
			/// @param data Input time series data
			/// @param k Lag amount (number of positions to shift)
			/// @return Vector of size (n - k) with lagged values
			static Vector<Real> Lag(const Vector<Real>& data, int k) {
				int n = data.size();

				if (n == 0)
					throw StatisticsError("Lag: Input data cannot be empty");
				if (k < 0)
					throw StatisticsError("Lag: Lag must be non-negative");
				if (k >= n)
					throw StatisticsError("Lag: Lag must be less than data length");

				Vector<Real> result(n - k);
				for (int i = 0; i < n - k; ++i)
					result[i] = data[i + k];

				return result;
			}

			/// @brief Create lag matrix with columns for lags 0, 1, ..., maxLag
			/// @details Each column i contains the series lagged by i positions.
			///          Row j, Column i = data[j + i] for valid indices.
			/// @param data Input time series data
			/// @param maxLag Maximum lag (number of columns - 1)
			/// @return Matrix of size (n - maxLag) x (maxLag + 1)
			static Matrix<Real> LagMatrix(const Vector<Real>& data, int maxLag) {
				int n = data.size();

				if (n == 0)
					throw StatisticsError("LagMatrix: Input data cannot be empty");
				if (maxLag < 0)
					throw StatisticsError("LagMatrix: maxLag must be non-negative");
				if (maxLag >= n)
					throw StatisticsError("LagMatrix: maxLag must be less than data length");

				int nRows = n - maxLag;
				int nCols = maxLag + 1;
				Matrix<Real> result(nRows, nCols);

				for (int i = 0; i < nRows; ++i)
					for (int j = 0; j <= maxLag; ++j)
						result(i, j) = data[maxLag + i - j];

				return result;
			}

			/// @brief First difference of time series
			/// @details Computes d_i = x_{i+1} - x_i
			/// @param data Input time series data
			/// @return Vector of size (n - 1) containing first differences
			static Vector<Real> FirstDifference(const Vector<Real>& data) {
				int n = data.size();

				if (n < 2)
					throw StatisticsError("FirstDifference: Data must have at least 2 elements");

				Vector<Real> result(n - 1);
				for (int i = 0; i < n - 1; ++i)
					result[i] = data[i + 1] - data[i];

				return result;
			}

			/// @brief Apply differencing multiple times
			/// @details order=1 is first difference, order=2 is difference of differences, etc.
			/// @param data Input time series data
			/// @param order Number of times to difference (must be >= 1)
			/// @return Vector of size (n - order) containing differenced values
			static Vector<Real> Difference(const Vector<Real>& data, int order) {
				int n = data.size();

				if (order < 1)
					throw StatisticsError("Difference: Order must be at least 1");
				if (order >= n)
					throw StatisticsError("Difference: Order must be less than data length");

				Vector<Real> result = FirstDifference(data);

				for (int d = 2; d <= order; ++d)
					result = FirstDifference(result);

				return result;
			}

			/// @brief Seasonal difference
			/// @details Computes d_i = x_i - x_{i-period}
			/// @param data Input time series data
			/// @param period Seasonal period (e.g., 12 for monthly, 4 for quarterly)
			/// @return Vector of size (n - period) containing seasonal differences
			static Vector<Real> SeasonalDifference(const Vector<Real>& data, int period) {
				int n = data.size();

				if (period < 1)
					throw StatisticsError("SeasonalDifference: Period must be at least 1");
				if (period >= n)
					throw StatisticsError("SeasonalDifference: Period must be less than data length");

				Vector<Real> result(n - period);
				for (int i = 0; i < n - period; ++i)
					result[i] = data[i + period] - data[i];

				return result;
			}

			/// @brief Cumulative sum (undoes differencing)
			/// @details cumsum[0] = initial, cumsum[i] = cumsum[i-1] + differences[i-1]
			/// @param differences The differenced series
			/// @param initial Starting value (typically the first value of original series)
			/// @return Vector of size (differences.size() + 1) with reconstructed values
			static Vector<Real> CumulativeSum(const Vector<Real>& differences, Real initial) {
				int n = differences.size();

				Vector<Real> result(n + 1);
				result[0] = initial;

				for (int i = 0; i < n; ++i)
					result[i + 1] = result[i] + differences[i];

				return result;
			}

			/// @brief Cumulative sum without initial value (starts at 0)
			/// @param differences The differenced series
			/// @return Vector of size differences.size() with cumulative sums
			static Vector<Real> CumulativeSum(const Vector<Real>& differences) {
				int n = differences.size();

				if (n == 0)
					return Vector<Real>();

				Vector<Real> result(n);
				result[0] = differences[0];

				for (int i = 1; i < n; ++i)
					result[i] = result[i - 1] + differences[i];

				return result;
			}

		} // namespace TimeSeries
	} // namespace Statistics
} // namespace MML

#endif // MML_TIME_SERIES_H

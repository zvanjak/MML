///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Histogram.h                                                         ///
///  Description: Histogram computation and frequency analysis functions including    ///
///               multiple binning methods, CDF, and quantile computation             ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_HISTOGRAM_H
#define MML_HISTOGRAM_H

#include "MMLBase.h"
#include "MMLExceptions.h"
#include "base/Vector.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>

namespace MML {
	namespace Statistics {
		namespace Histogram {

			/////////////////////////////////////////////////////////////////////////////////////
			///                         BINNING METHOD ENUM                                   ///
			/////////////////////////////////////////////////////////////////////////////////////

			/// @brief Methods for automatic bin width/count selection
			enum class BinningMethod {
				Sturges,           ///< Sturges' formula: k = ceil(log2(n) + 1)
				Scott,             ///< Scott's rule: h = 3.49 * std * n^(-1/3)
				FreedmanDiaconis,  ///< Freedman-Diaconis: h = 2 * IQR * n^(-1/3)
				SquareRoot,        ///< Square root choice: k = ceil(sqrt(n))
				Rice               ///< Rice's rule: k = ceil(2 * n^(1/3))
			};

			/////////////////////////////////////////////////////////////////////////////////////
			///                         HISTOGRAM RESULT STRUCT                               ///
			/////////////////////////////////////////////////////////////////////////////////////

			/// @brief Result of histogram computation
			struct HistogramResult {
				Vector<Real> binEdges;      ///< Bin edges (size = numBins + 1)
				Vector<int> counts;          ///< Raw counts per bin
				Vector<Real> frequencies;    ///< Relative frequencies (counts / total)
				Vector<Real> density;        ///< Probability density (frequency / binWidth)
				Real binWidth;               ///< Width of each bin (uniform)
				int numBins;                 ///< Number of bins
				int totalCount;              ///< Total number of data points

				/// @brief Get bin centers
				Vector<Real> GetBinCenters() const {
					Vector<Real> centers(numBins);
					for (int i = 0; i < numBins; ++i) {
						centers[i] = (binEdges[i] + binEdges[i + 1]) / 2.0;
					}
					return centers;
				}

				/// @brief Get cumulative counts
				Vector<int> GetCumulativeCounts() const {
					Vector<int> cumCounts(numBins);
					cumCounts[0] = counts[0];
					for (int i = 1; i < numBins; ++i) {
						cumCounts[i] = cumCounts[i - 1] + counts[i];
					}
					return cumCounts;
				}

				/// @brief Get cumulative frequencies (empirical CDF at bin edges)
				Vector<Real> GetCumulativeFrequencies() const {
					Vector<Real> cumFreq(numBins);
					cumFreq[0] = frequencies[0];
					for (int i = 1; i < numBins; ++i) {
						cumFreq[i] = cumFreq[i - 1] + frequencies[i];
					}
					return cumFreq;
				}
			};

			/////////////////////////////////////////////////////////////////////////////////////
			///                         BIN COUNT ESTIMATION                                  ///
			/////////////////////////////////////////////////////////////////////////////////////

			/// @brief Estimate optimal number of bins using Sturges' formula
			/// @details k = ceil(log2(n) + 1)
			///          Works well for normally distributed data
			/// @param n Number of data points
			/// @return Recommended number of bins
			static int SturgesBinCount(int n) {
				if (n <= 0) return 1;
				return static_cast<int>(std::ceil(std::log2(static_cast<double>(n)) + 1.0));
			}

			/// @brief Estimate optimal number of bins using Rice's rule
			/// @details k = ceil(2 * n^(1/3))
			/// @param n Number of data points
			/// @return Recommended number of bins
			static int RiceBinCount(int n) {
				if (n <= 0) return 1;
				double value = 2.0 * std::cbrt(static_cast<double>(n));
				// Round to nearest integer if very close (within 1e-9) to handle floating-point precision
				double rounded = std::round(value);
				if (std::abs(value - rounded) < 1e-9)
					return static_cast<int>(rounded);
				return static_cast<int>(std::ceil(value));
			}

			/// @brief Estimate optimal number of bins using square root rule
			/// @details k = ceil(sqrt(n))
			/// @param n Number of data points
			/// @return Recommended number of bins
			static int SquareRootBinCount(int n) {
				if (n <= 0) return 1;
				return static_cast<int>(std::ceil(std::sqrt(static_cast<double>(n))));
			}

			/// @brief Estimate optimal bin width using Scott's rule
			/// @details h = 3.49 * std * n^(-1/3)
			///          Optimal for normally distributed data
			/// @param data Input data
			/// @return Recommended bin width
			static Real ScottBinWidth(const Vector<Real>& data) {
				int n = data.size();
				if (n <= 1) return 1.0;

				// Compute standard deviation
				Real mean = 0.0;
				for (int i = 0; i < n; ++i) mean += data[i];
				mean /= n;

				Real variance = 0.0;
				for (int i = 0; i < n; ++i) {
					Real diff = data[i] - mean;
					variance += diff * diff;
				}
				variance /= (n - 1);
				Real stdDev = std::sqrt(variance);

				if (stdDev < Constants::Eps) return 1.0;

				return 3.49 * stdDev * std::pow(static_cast<double>(n), -1.0 / 3.0);
			}

			/// @brief Estimate optimal bin width using Freedman-Diaconis rule
			/// @details h = 2 * IQR * n^(-1/3)
			///          More robust to outliers than Scott's rule
			/// @param data Input data
			/// @return Recommended bin width
			static Real FreedmanDiaconisBinWidth(const Vector<Real>& data) {
				int n = data.size();
				if (n <= 1) return 1.0;

				// Sort data for percentile calculation
				std::vector<Real> sorted(n);
				for (int i = 0; i < n; ++i) sorted[i] = data[i];
				std::sort(sorted.begin(), sorted.end());

				// Compute Q1 (25th percentile) and Q3 (75th percentile)
				auto percentile = [&sorted, n](Real p) -> Real {
					Real rank = (p / 100.0) * (n - 1);
					int lower = static_cast<int>(std::floor(rank));
					int upper = static_cast<int>(std::ceil(rank));
					if (lower == upper) return sorted[lower];
					Real fraction = rank - lower;
					return sorted[lower] + fraction * (sorted[upper] - sorted[lower]);
				};

				Real iqr = percentile(75.0) - percentile(25.0);

				if (iqr < Constants::Eps) {
					// Fallback to Scott's rule if IQR is zero
					return ScottBinWidth(data);
				}

				return 2.0 * iqr * std::pow(static_cast<double>(n), -1.0 / 3.0);
			}

			/// @brief Get bin count for a given method
			/// @param data Input data
			/// @param method Binning method to use
			/// @return Recommended number of bins
			static int GetBinCount(const Vector<Real>& data, BinningMethod method) {
				int n = data.size();
				if (n <= 0) return 1;

				// Find data range
				Real minVal = data[0], maxVal = data[0];
				for (int i = 1; i < n; ++i) {
					if (data[i] < minVal) minVal = data[i];
					if (data[i] > maxVal) maxVal = data[i];
				}
				Real range = maxVal - minVal;

				switch (method) {
				case BinningMethod::Sturges:
					return SturgesBinCount(n);

				case BinningMethod::Rice:
					return RiceBinCount(n);

				case BinningMethod::SquareRoot:
					return SquareRootBinCount(n);

				case BinningMethod::Scott: {
					Real width = ScottBinWidth(data);
					if (width < Constants::Eps || range < Constants::Eps) return SturgesBinCount(n);
					return std::max(1, static_cast<int>(std::ceil(range / width)));
				}

				case BinningMethod::FreedmanDiaconis: {
					Real width = FreedmanDiaconisBinWidth(data);
					if (width < Constants::Eps || range < Constants::Eps) return SturgesBinCount(n);
					return std::max(1, static_cast<int>(std::ceil(range / width)));
				}

				default:
					return SturgesBinCount(n);
				}
			}

			/////////////////////////////////////////////////////////////////////////////////////
			///                         CORE HISTOGRAM FUNCTIONS                              ///
			/////////////////////////////////////////////////////////////////////////////////////

			/// @brief Compute histogram with specified number of bins
			/// @details Creates uniform bins spanning [min(data), max(data)]
			/// @param data Input data
			/// @param numBins Number of bins to use
			/// @return HistogramResult containing bin edges, counts, frequencies, and density
			/// @throws StatisticsError if data is empty or numBins < 1
			static HistogramResult ComputeHistogram(const Vector<Real>& data, int numBins) {
				int n = data.size();

				if (n == 0)
					throw StatisticsError("ComputeHistogram: Input data cannot be empty");
				if (numBins < 1)
					throw StatisticsError("ComputeHistogram: Number of bins must be at least 1");

				// Find data range
				Real minVal = data[0], maxVal = data[0];
				for (int i = 1; i < n; ++i) {
					if (data[i] < minVal) minVal = data[i];
					if (data[i] > maxVal) maxVal = data[i];
				}

				// Handle constant data case
				if (maxVal - minVal < Constants::Eps) {
					HistogramResult result;
					result.numBins = 1;
					result.totalCount = n;
					result.binWidth = 1.0;  // Arbitrary width for constant data
					result.binEdges = Vector<Real>(2);
					result.binEdges[0] = minVal - 0.5;
					result.binEdges[1] = minVal + 0.5;
					result.counts = Vector<int>(1);
					result.counts[0] = n;
					result.frequencies = Vector<Real>(1);
					result.frequencies[0] = 1.0;
					result.density = Vector<Real>(1);
					result.density[0] = 1.0;
					return result;
				}

				// Add small padding to include max value in last bin
				Real padding = (maxVal - minVal) * 1e-10;
				maxVal += padding;

				Real binWidth = (maxVal - minVal) / numBins;

				// Create bin edges
				Vector<Real> binEdges(numBins + 1);
				for (int i = 0; i <= numBins; ++i) {
					binEdges[i] = minVal + i * binWidth;
				}

				// Count values in each bin
				Vector<int> counts(numBins);
				for (int i = 0; i < numBins; ++i) counts[i] = 0;

				for (int i = 0; i < n; ++i) {
					int binIndex = static_cast<int>((data[i] - minVal) / binWidth);
					// Clamp to valid range (handles edge cases)
					if (binIndex < 0) binIndex = 0;
					if (binIndex >= numBins) binIndex = numBins - 1;
					counts[binIndex]++;
				}

				// Compute frequencies and density
				Vector<Real> frequencies(numBins);
				Vector<Real> density(numBins);
				for (int i = 0; i < numBins; ++i) {
					frequencies[i] = static_cast<Real>(counts[i]) / n;
					density[i] = frequencies[i] / binWidth;
				}

				HistogramResult result;
				result.binEdges = binEdges;
				result.counts = counts;
				result.frequencies = frequencies;
				result.density = density;
				result.binWidth = binWidth;
				result.numBins = numBins;
				result.totalCount = n;

				return result;
			}

			/// @brief Compute histogram with custom bin edges
			/// @details Bins are [edge[i], edge[i+1]) except last bin is [edge[n-1], edge[n]]
			/// @param data Input data
			/// @param binEdges Custom bin edges (must be sorted ascending)
			/// @return HistogramResult containing counts, frequencies, and density
			/// @throws StatisticsError if data is empty or binEdges has fewer than 2 elements
			static HistogramResult ComputeHistogram(const Vector<Real>& data, const Vector<Real>& binEdges) {
				int n = data.size();
				int numBins = binEdges.size() - 1;

				if (n == 0)
					throw StatisticsError("ComputeHistogram: Input data cannot be empty");
				if (numBins < 1)
					throw StatisticsError("ComputeHistogram: Need at least 2 bin edges");

				// Count values in each bin using binary search
				Vector<int> counts(numBins);
				for (int i = 0; i < numBins; ++i) counts[i] = 0;

				for (int i = 0; i < n; ++i) {
					Real val = data[i];
					// Find bin using linear search (for small numBins) or binary search
					for (int b = 0; b < numBins; ++b) {
						if (val >= binEdges[b] && (val < binEdges[b + 1] || (b == numBins - 1 && val <= binEdges[b + 1]))) {
							counts[b]++;
							break;
						}
					}
				}

				// Compute frequencies and density
				Vector<Real> frequencies(numBins);
				Vector<Real> density(numBins);
				for (int i = 0; i < numBins; ++i) {
					Real width = binEdges[i + 1] - binEdges[i];
					frequencies[i] = static_cast<Real>(counts[i]) / n;
					density[i] = (width > Constants::Eps) ? frequencies[i] / width : 0.0;
				}

				// Compute average bin width
				Real avgBinWidth = (binEdges[numBins] - binEdges[0]) / numBins;

				HistogramResult result;
				result.binEdges = binEdges;
				result.counts = counts;
				result.frequencies = frequencies;
				result.density = density;
				result.binWidth = avgBinWidth;
				result.numBins = numBins;
				result.totalCount = n;

				return result;
			}

			/// @brief Compute histogram with automatic bin selection
			/// @details Uses Sturges' rule by default
			/// @param data Input data
			/// @param method Binning method (default: Sturges)
			/// @return HistogramResult
			static HistogramResult ComputeHistogramAuto(const Vector<Real>& data,
				BinningMethod method = BinningMethod::Sturges) {
				int numBins = GetBinCount(data, method);
				return ComputeHistogram(data, numBins);
			}

			/////////////////////////////////////////////////////////////////////////////////////
			///                         CONVENIENCE METHODS                                   ///
			/////////////////////////////////////////////////////////////////////////////////////

			/// @brief Compute histogram using Sturges' rule
			static HistogramResult HistogramSturges(const Vector<Real>& data) {
				return ComputeHistogramAuto(data, BinningMethod::Sturges);
			}

			/// @brief Compute histogram using Scott's rule
			static HistogramResult HistogramScott(const Vector<Real>& data) {
				return ComputeHistogramAuto(data, BinningMethod::Scott);
			}

			/// @brief Compute histogram using Freedman-Diaconis rule
			static HistogramResult HistogramFD(const Vector<Real>& data) {
				return ComputeHistogramAuto(data, BinningMethod::FreedmanDiaconis);
			}

			/// @brief Compute histogram using square root rule
			static HistogramResult HistogramSqrt(const Vector<Real>& data) {
				return ComputeHistogramAuto(data, BinningMethod::SquareRoot);
			}

			/// @brief Compute histogram using Rice's rule
			static HistogramResult HistogramRice(const Vector<Real>& data) {
				return ComputeHistogramAuto(data, BinningMethod::Rice);
			}

			/////////////////////////////////////////////////////////////////////////////////////
			///                         FREQUENCY TABLE                                       ///
			/////////////////////////////////////////////////////////////////////////////////////

			/// @brief Result of frequency table computation for discrete values
			struct FrequencyTableResult {
				std::map<Real, int> counts;       ///< Value -> count mapping
				std::map<Real, Real> frequencies; ///< Value -> relative frequency mapping
				int totalCount;                   ///< Total number of data points
				int uniqueCount;                  ///< Number of unique values

				/// @brief Get sorted unique values
				Vector<Real> GetValues() const {
					Vector<Real> values(uniqueCount);
					int i = 0;
					for (const auto& pair : counts) {
						values[i++] = pair.first;
					}
					return values;
				}

				/// @brief Get counts in same order as GetValues()
				Vector<int> GetCounts() const {
					Vector<int> result(uniqueCount);
					int i = 0;
					for (const auto& pair : counts) {
						result[i++] = pair.second;
					}
					return result;
				}

				/// @brief Get frequencies in same order as GetValues()
				Vector<Real> GetFrequencies() const {
					Vector<Real> result(uniqueCount);
					int i = 0;
					for (const auto& pair : frequencies) {
						result[i++] = pair.second;
					}
					return result;
				}
			};

			/// @brief Compute frequency table for discrete values
			/// @details Counts occurrences of each unique value
			/// @param data Input data (typically integers or discrete values)
			/// @return FrequencyTableResult containing value counts and frequencies
			/// @throws StatisticsError if data is empty
			static FrequencyTableResult FrequencyTable(const Vector<Real>& data) {
				int n = data.size();

				if (n == 0)
					throw StatisticsError("FrequencyTable: Input data cannot be empty");

				FrequencyTableResult result;
				result.totalCount = n;

				// Count occurrences
				for (int i = 0; i < n; ++i) {
					result.counts[data[i]]++;
				}

				// Compute frequencies
				for (const auto& pair : result.counts) {
					result.frequencies[pair.first] = static_cast<Real>(pair.second) / n;
				}

				result.uniqueCount = static_cast<int>(result.counts.size());

				return result;
			}

			/////////////////////////////////////////////////////////////////////////////////////
			///                         CUMULATIVE DISTRIBUTION                               ///
			/////////////////////////////////////////////////////////////////////////////////////

			/// @brief Result of empirical CDF computation
			struct ECDFResult {
				Vector<Real> x;        ///< Sorted unique x values
				Vector<Real> cdf;      ///< CDF values at each x
				int n;                 ///< Total number of data points
			};

			/// @brief Compute the Empirical Cumulative Distribution Function (ECDF)
			/// @details F(x) = (number of observations <= x) / n
			///          Returns step function values at each unique data point
			/// @param data Input data
			/// @return ECDFResult containing sorted x values and corresponding CDF values
			/// @throws StatisticsError if data is empty
			static ECDFResult EmpiricalCDF(const Vector<Real>& data) {
				int n = data.size();

				if (n == 0)
					throw StatisticsError("EmpiricalCDF: Input data cannot be empty");

				// Sort data
				std::vector<Real> sorted(n);
				for (int i = 0; i < n; ++i) sorted[i] = data[i];
				std::sort(sorted.begin(), sorted.end());

				// Find unique values and their CDF values
				std::vector<Real> uniqueX;
				std::vector<Real> cdfValues;

				Real prevVal = sorted[0];
				int count = 1;

				for (int i = 1; i <= n; ++i) {
					if (i == n || sorted[i] != prevVal) {
						uniqueX.push_back(prevVal);
						cdfValues.push_back(static_cast<Real>(i) / n);
						if (i < n) {
							prevVal = sorted[i];
							count = 1;
						}
					}
					else {
						count++;
					}
				}

				ECDFResult result;
				result.n = n;
				result.x = Vector<Real>(static_cast<int>(uniqueX.size()));
				result.cdf = Vector<Real>(static_cast<int>(cdfValues.size()));
				for (size_t i = 0; i < uniqueX.size(); ++i) {
					result.x[static_cast<int>(i)] = uniqueX[i];
					result.cdf[static_cast<int>(i)] = cdfValues[i];
				}

				return result;
			}

			/// @brief Evaluate ECDF at a given point
			/// @details Returns proportion of data <= x
			/// @param data Input data (not necessarily sorted)
			/// @param x Point at which to evaluate CDF
			/// @return F(x) = P(X <= x)
			static Real EvaluateECDF(const Vector<Real>& data, Real x) {
				int n = data.size();
				if (n == 0) return 0.0;

				int count = 0;
				for (int i = 0; i < n; ++i) {
					if (data[i] <= x) count++;
				}
				return static_cast<Real>(count) / n;
			}

			/////////////////////////////////////////////////////////////////////////////////////
			///                         QUANTILE COMPUTATION                                  ///
			/////////////////////////////////////////////////////////////////////////////////////

			/// @brief Compute multiple quantiles efficiently
			/// @details Sorts data once and computes all requested quantiles
			///          Uses linear interpolation (similar to Excel PERCENTILE.INC)
			/// @param data Input data
			/// @param probabilities Vector of probabilities in [0, 1]
			/// @return Vector of quantile values corresponding to each probability
			/// @throws StatisticsError if data is empty or any probability is out of range
			static Vector<Real> Quantiles(const Vector<Real>& data, const Vector<Real>& probabilities) {
				int n = data.size();
				int numQuantiles = probabilities.size();

				if (n == 0)
					throw StatisticsError("Quantiles: Input data cannot be empty");
				if (numQuantiles == 0)
					throw StatisticsError("Quantiles: Probabilities vector cannot be empty");

				// Validate probabilities
				for (int i = 0; i < numQuantiles; ++i) {
					if (probabilities[i] < 0.0 || probabilities[i] > 1.0)
						throw StatisticsError("Quantiles: All probabilities must be in [0, 1]");
				}

				// Sort data once
				std::vector<Real> sorted(n);
				for (int i = 0; i < n; ++i) sorted[i] = data[i];
				std::sort(sorted.begin(), sorted.end());

				// Compute each quantile
				Vector<Real> result(numQuantiles);
				for (int i = 0; i < numQuantiles; ++i) {
					Real p = probabilities[i];

					if (p == 0.0) {
						result[i] = sorted[0];
					}
					else if (p == 1.0) {
						result[i] = sorted[n - 1];
					}
					else {
						Real rank = p * (n - 1);
						int lower = static_cast<int>(std::floor(rank));
						int upper = static_cast<int>(std::ceil(rank));
						Real fraction = rank - lower;
						result[i] = sorted[lower] + fraction * (sorted[upper] - sorted[lower]);
					}
				}

				return result;
			}

			/// @brief Compute a single quantile
			/// @param data Input data
			/// @param p Probability in [0, 1]
			/// @return The p-th quantile
			static Real Quantile(const Vector<Real>& data, Real p) {
				Vector<Real> probs(1);
				probs[0] = p;
				return Quantiles(data, probs)[0];
			}

			/////////////////////////////////////////////////////////////////////////////////////
			///                         BIN COUNTING AND DIGITIZATION                         ///
			/////////////////////////////////////////////////////////////////////////////////////

			/// @brief Count data points in each bin (without normalization)
			/// @details Simple bin counting, bins are [edge[i], edge[i+1])
			///          except last bin which is [edge[n-1], edge[n]]
			/// @param data Input data
			/// @param binEdges Bin edges (must be sorted ascending, size = numBins + 1)
			/// @return Vector of counts per bin
			/// @throws StatisticsError if data is empty or binEdges has fewer than 2 elements
			static Vector<int> BinCount(const Vector<Real>& data, const Vector<Real>& binEdges) {
				int n = data.size();
				int numBins = binEdges.size() - 1;

				if (n == 0)
					throw StatisticsError("BinCount: Input data cannot be empty");
				if (numBins < 1)
					throw StatisticsError("BinCount: Need at least 2 bin edges");

				Vector<int> counts(numBins);
				for (int i = 0; i < numBins; ++i) counts[i] = 0;

				for (int i = 0; i < n; ++i) {
					Real val = data[i];
					// Find appropriate bin
					for (int b = 0; b < numBins; ++b) {
						if (val >= binEdges[b] && (val < binEdges[b + 1] || (b == numBins - 1 && val <= binEdges[b + 1]))) {
							counts[b]++;
							break;
						}
					}
				}

				return counts;
			}

			/// @brief Determine bin index for each data point (digitization)
			/// @details Returns the index of the bin each value falls into
			///          Values outside range get -1 (below) or numBins (above)
			/// @param data Input data
			/// @param binEdges Bin edges (must be sorted ascending)
			/// @return Vector of bin indices (0 to numBins-1, or -1/numBins for out of range)
			/// @throws StatisticsError if data is empty or binEdges has fewer than 2 elements
			static Vector<int> Digitize(const Vector<Real>& data, const Vector<Real>& binEdges) {
				int n = data.size();
				int numBins = binEdges.size() - 1;

				if (n == 0)
					throw StatisticsError("Digitize: Input data cannot be empty");
				if (numBins < 1)
					throw StatisticsError("Digitize: Need at least 2 bin edges");

				Vector<int> indices(n);

				for (int i = 0; i < n; ++i) {
					Real val = data[i];

					// Below range
					if (val < binEdges[0]) {
						indices[i] = -1;
						continue;
					}

					// Above range
					if (val > binEdges[numBins]) {
						indices[i] = numBins;
						continue;
					}

					// Find bin using linear search
					bool found = false;
					for (int b = 0; b < numBins; ++b) {
						if (val >= binEdges[b] && (val < binEdges[b + 1] || (b == numBins - 1 && val <= binEdges[b + 1]))) {
							indices[i] = b;
							found = true;
							break;
						}
					}

					if (!found) {
						indices[i] = numBins - 1;  // Fallback to last bin
					}
				}

				return indices;
			}

			/////////////////////////////////////////////////////////////////////////////////////
			///                         UTILITY FUNCTIONS                                     ///
			/////////////////////////////////////////////////////////////////////////////////////

			/// @brief Create uniform bin edges from min to max
			/// @param minVal Minimum value (left edge of first bin)
			/// @param maxVal Maximum value (right edge of last bin)
			/// @param numBins Number of bins
			/// @return Vector of bin edges (size = numBins + 1)
			static Vector<Real> CreateUniformBinEdges(Real minVal, Real maxVal, int numBins) {
				if (numBins < 1)
					throw StatisticsError("CreateUniformBinEdges: Number of bins must be at least 1");

				Vector<Real> edges(numBins + 1);
				Real width = (maxVal - minVal) / numBins;
				for (int i = 0; i <= numBins; ++i) {
					edges[i] = minVal + i * width;
				}
				return edges;
			}

			/// @brief Create logarithmic bin edges (useful for log-scale histograms)
			/// @param minVal Minimum value (must be positive)
			/// @param maxVal Maximum value
			/// @param numBins Number of bins
			/// @return Vector of bin edges with logarithmic spacing
			static Vector<Real> CreateLogBinEdges(Real minVal, Real maxVal, int numBins) {
				if (numBins < 1)
					throw StatisticsError("CreateLogBinEdges: Number of bins must be at least 1");
				if (minVal <= 0)
					throw StatisticsError("CreateLogBinEdges: Minimum value must be positive");
				if (maxVal <= minVal)
					throw StatisticsError("CreateLogBinEdges: Maximum must be greater than minimum");

				Vector<Real> edges(numBins + 1);
				Real logMin = std::log10(minVal);
				Real logMax = std::log10(maxVal);
				Real logWidth = (logMax - logMin) / numBins;

				for (int i = 0; i <= numBins; ++i) {
					edges[i] = std::pow(10.0, logMin + i * logWidth);
				}
				return edges;
			}

		}  // namespace Histogram
	}  // namespace Statistics
}  // namespace MML

#endif  // MML_HISTOGRAM_H

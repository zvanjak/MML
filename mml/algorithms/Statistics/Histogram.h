///////////////////////////////////////////////////////////////////////////////////////////
// Histogram.h
//
// Core histogram result types and binning-method selection.
// Histogram computation, frequency tables, ECDF, and quantile helpers are declared
// here as their migration tasks land.
///////////////////////////////////////////////////////////////////////////////////////////

// Intentionally distinct from the package's MML_HISTOGRAM_H guard so the
// compatibility wrapper can include this core header after migration completes.
#if !defined MML_STATISTICS_HISTOGRAM_H
#define MML_STATISTICS_HISTOGRAM_H

#include <mml/algorithms/Statistics.h>

#include <map>

namespace MML
{
	namespace Statistics
	{
		namespace Histogram
		{
			namespace Detail
			{
				inline void ValidateBinningData(const Vector<Real>& data, const char* functionName)
				{
					if (data.size() == 0)
						throw StatisticsError(std::string(functionName) + ": Input data cannot be empty");

					for (int i = 0; i < data.size(); ++i) {
						if (!std::isfinite(data[i]))
							throw StatisticsError(std::string(functionName) + ": Input data must be finite");
					}
				}

				inline int StableCeilToBinCount(double value)
				{
					if (!std::isfinite(value) || value <= 0.0)
						throw StatisticsError("Histogram bin-count calculation produced an invalid value");
					if (value > static_cast<double>(std::numeric_limits<int>::max()))
						throw StatisticsError("Histogram bin count exceeds the supported integer range");

					double nearestInteger = std::round(value);
					double tolerance = 32.0 * std::numeric_limits<double>::epsilon()
						* std::max(1.0, std::abs(value));
					if (std::abs(value - nearestInteger) <= tolerance)
						value = nearestInteger;

					return std::max(1, static_cast<int>(std::ceil(value)));
				}
			} // namespace Detail

			enum class BinningMethod
			{
				Sturges,
				Scott,
				FreedmanDiaconis,
				SquareRoot,
				Rice
			};

			struct HistogramResult
			{
				Vector<Real> binEdges;
				Vector<int> counts;
				Vector<Real> frequencies;
				Vector<Real> density;
				Real binWidth = 0.0;
				int numBins = 0;
				int totalCount = 0;

				[[nodiscard]] Vector<Real> GetBinCenters() const
				{
					Vector<Real> centers(numBins);
					for (int i = 0; i < numBins; ++i)
						centers[i] = binEdges[i] + (binEdges[i + 1] - binEdges[i]) / REAL(2.0);
					return centers;
				}

				[[nodiscard]] Vector<int> GetCumulativeCounts() const
				{
					Vector<int> cumulative(numBins);
					int runningCount = 0;
					for (int i = 0; i < numBins; ++i) {
						runningCount += counts[i];
						cumulative[i] = runningCount;
					}
					return cumulative;
				}

				[[nodiscard]] Vector<Real> GetCumulativeFrequencies() const
				{
					Vector<Real> cumulative(numBins);
					Real runningFrequency = 0.0;
					for (int i = 0; i < numBins; ++i) {
						runningFrequency += frequencies[i];
						cumulative[i] = runningFrequency;
					}
					return cumulative;
				}
			};

			struct FrequencyTableResult
			{
				std::map<Real, int> counts;
				std::map<Real, Real> frequencies;
				int totalCount = 0;
				int uniqueCount = 0;

				[[nodiscard]] Vector<Real> GetValues() const
				{
					Vector<Real> values(uniqueCount);
					int i = 0;
					for (const auto& pair : counts)
						values[i++] = pair.first;
					return values;
				}

				[[nodiscard]] Vector<int> GetCounts() const
				{
					Vector<int> result(uniqueCount);
					int i = 0;
					for (const auto& pair : counts)
						result[i++] = pair.second;
					return result;
				}

				[[nodiscard]] Vector<Real> GetFrequencies() const
				{
					Vector<Real> result(uniqueCount);
					int i = 0;
					for (const auto& pair : frequencies)
						result[i++] = pair.second;
					return result;
				}
			};

			struct ECDFResult
			{
				Vector<Real> x;
				Vector<Real> cdf;
				int n = 0;
			};

			/// Sturges' rule: ceil(log2(n) + 1).
			inline int SturgesBinCount(int n)
			{
				if (n <= 0)
					throw StatisticsError("SturgesBinCount requires a positive sample count");
				return Detail::StableCeilToBinCount(std::log2(static_cast<double>(n)) + 1.0);
			}

			/// Rice's rule: ceil(2 * cbrt(n)).
			inline int RiceBinCount(int n)
			{
				if (n <= 0)
					throw StatisticsError("RiceBinCount requires a positive sample count");
				return Detail::StableCeilToBinCount(2.0 * std::cbrt(static_cast<double>(n)));
			}

			/// Square-root rule: ceil(sqrt(n)).
			inline int SquareRootBinCount(int n)
			{
				if (n <= 0)
					throw StatisticsError("SquareRootBinCount requires a positive sample count");
				return Detail::StableCeilToBinCount(std::sqrt(static_cast<double>(n)));
			}

			/// Scott's rule: 3.49 * sample_stddev * n^(-1/3).
			inline Real ScottBinWidth(const Vector<Real>& data)
			{
				Detail::ValidateBinningData(data, "ScottBinWidth");
				if (data.size() == 1)
					return REAL(1.0);

				Real stdDev = SampleStdDev(data);
				if (!std::isfinite(stdDev))
					throw StatisticsError("ScottBinWidth produced a non-finite standard deviation");
				if (stdDev <= Constants::Eps)
					return REAL(1.0);

				Real sampleScale = static_cast<Real>(std::pow(static_cast<double>(data.size()), -1.0 / 3.0));
				return REAL(3.49) * stdDev * sampleScale;
			}

			/// Freedman-Diaconis rule: 2 * IQR * n^(-1/3).
			inline Real FreedmanDiaconisBinWidth(const Vector<Real>& data)
			{
				Detail::ValidateBinningData(data, "FreedmanDiaconisBinWidth");
				if (data.size() == 1)
					return REAL(1.0);

				Real iqr = IQR(data);
				if (!std::isfinite(iqr))
					throw StatisticsError("FreedmanDiaconisBinWidth produced a non-finite IQR");
				if (iqr <= Constants::Eps)
					return ScottBinWidth(data);

				Real sampleScale = static_cast<Real>(std::pow(static_cast<double>(data.size()), -1.0 / 3.0));
				return REAL(2.0) * iqr * sampleScale;
			}

			inline int GetBinCount(const Vector<Real>& data, BinningMethod method)
			{
				Detail::ValidateBinningData(data, "GetBinCount");

				Real minValue = data[0];
				Real maxValue = data[0];
				for (int i = 1; i < data.size(); ++i) {
					minValue = std::min(minValue, data[i]);
					maxValue = std::max(maxValue, data[i]);
				}

				Real range = maxValue - minValue;
				if (!std::isfinite(range))
					throw StatisticsError("GetBinCount: Data range overflow");
				switch (method) {
				case BinningMethod::Sturges:
					return SturgesBinCount(data.size());
				case BinningMethod::Rice:
					return RiceBinCount(data.size());
				case BinningMethod::SquareRoot:
					return SquareRootBinCount(data.size());
				case BinningMethod::Scott: {
					Real width = ScottBinWidth(data);
					if (width <= Constants::Eps || range <= Constants::Eps)
						return SturgesBinCount(data.size());
					return Detail::StableCeilToBinCount(static_cast<double>(range / width));
				}
				case BinningMethod::FreedmanDiaconis: {
					Real width = FreedmanDiaconisBinWidth(data);
					if (width <= Constants::Eps || range <= Constants::Eps)
						return SturgesBinCount(data.size());
					return Detail::StableCeilToBinCount(static_cast<double>(range / width));
				}
				default:
					throw StatisticsError("GetBinCount: Unknown binning method");
				}
			}

			inline HistogramResult ComputeHistogram(const Vector<Real>& data, int numBins)
			{
				Detail::ValidateBinningData(data, "ComputeHistogram");
				if (numBins < 1)
					throw StatisticsError("ComputeHistogram: Number of bins must be at least 1");
				if (numBins == std::numeric_limits<int>::max())
					throw StatisticsError("ComputeHistogram: Number of bins is too large");

				Real minValue = data[0];
				Real maxValue = data[0];
				for (int i = 1; i < data.size(); ++i) {
					minValue = std::min(minValue, data[i]);
					maxValue = std::max(maxValue, data[i]);
				}

				Real range = maxValue - minValue;
				if (!std::isfinite(range))
					throw StatisticsError("ComputeHistogram: Data range overflow");

				if (range < Constants::Eps) {
					Real lowerEdge = minValue - REAL(0.5);
					Real upperEdge = minValue + REAL(0.5);
					if (!(lowerEdge < minValue))
						lowerEdge = std::nextafter(minValue, -std::numeric_limits<Real>::infinity());
					if (!(upperEdge > minValue))
						upperEdge = std::nextafter(minValue, std::numeric_limits<Real>::infinity());
					if (!std::isfinite(lowerEdge))
						lowerEdge = minValue;
					if (!std::isfinite(upperEdge))
						upperEdge = minValue;
					if (!(lowerEdge < upperEdge))
						throw StatisticsError("ComputeHistogram: Cannot construct finite edges for constant data");

					HistogramResult result;
					result.binEdges = Vector<Real>({lowerEdge, upperEdge});
					result.counts = Vector<int>({data.size()});
					result.frequencies = Vector<Real>({REAL(1.0)});
					result.density = Vector<Real>({REAL(1.0) / (upperEdge - lowerEdge)});
					result.binWidth = upperEdge - lowerEdge;
					result.numBins = 1;
					result.totalCount = data.size();
					return result;
				}

				Real binWidth = range / static_cast<Real>(numBins);
				if (!std::isfinite(binWidth) || binWidth <= 0.0)
					throw StatisticsError("ComputeHistogram: Bin width is not representable");

				HistogramResult result;
				result.binEdges = Vector<Real>(numBins + 1);
				result.counts = Vector<int>(numBins);
				result.frequencies = Vector<Real>(numBins);
				result.density = Vector<Real>(numBins);
				result.binWidth = binWidth;
				result.numBins = numBins;
				result.totalCount = data.size();

				for (int i = 0; i < numBins; ++i)
					result.binEdges[i] = minValue + static_cast<Real>(i) * binWidth;
				result.binEdges[numBins] = maxValue;
				for (int i = 1; i <= numBins; ++i) {
					if (!std::isfinite(result.binEdges[i]) || !(result.binEdges[i] > result.binEdges[i - 1]))
						throw StatisticsError("ComputeHistogram: Uniform bin edges are not representable");
				}

				for (int i = 0; i < data.size(); ++i) {
					int binIndex;
					if (data[i] == maxValue) {
						binIndex = numBins - 1;
					} else {
						auto first = &result.binEdges[0];
						auto upper = std::upper_bound(first, first + result.binEdges.size(), data[i]);
						binIndex = static_cast<int>(upper - first) - 1;
					}
					++result.counts[binIndex];
				}

				for (int i = 0; i < numBins; ++i) {
					result.frequencies[i] = static_cast<Real>(result.counts[i]) / static_cast<Real>(data.size());
					result.density[i] = result.frequencies[i] / binWidth;
					if (!std::isfinite(result.density[i]))
						throw StatisticsError("ComputeHistogram: Density overflow");
				}

				return result;
			}

			inline HistogramResult ComputeHistogram(const Vector<Real>& data, const Vector<Real>& binEdges)
			{
				Detail::ValidateBinningData(data, "ComputeHistogram");
				if (binEdges.size() < 2)
					throw StatisticsError("ComputeHistogram: Need at least 2 bin edges");

				for (int i = 0; i < binEdges.size(); ++i) {
					if (!std::isfinite(binEdges[i]))
						throw StatisticsError("ComputeHistogram: Bin edges must be finite");
					if (i > 0 && !(binEdges[i] > binEdges[i - 1]))
						throw StatisticsError("ComputeHistogram: Bin edges must be strictly increasing");
				}

				int numBins = binEdges.size() - 1;
				Real span = binEdges[numBins] - binEdges[0];
				if (!std::isfinite(span) || span <= 0.0)
					throw StatisticsError("ComputeHistogram: Bin-edge span is not representable");

				HistogramResult result;
				result.binEdges = binEdges;
				result.counts = Vector<int>(numBins);
				result.frequencies = Vector<Real>(numBins);
				result.density = Vector<Real>(numBins);
				result.binWidth = span / static_cast<Real>(numBins);
				if (!std::isfinite(result.binWidth) || result.binWidth <= 0.0)
					throw StatisticsError("ComputeHistogram: Average bin width is not representable");
				result.numBins = numBins;
				result.totalCount = data.size();

				for (int i = 0; i < data.size(); ++i) {
					Real value = data[i];
					if (value < binEdges[0] || value > binEdges[numBins])
						continue;

					int binIndex;
					if (value == binEdges[numBins]) {
						binIndex = numBins - 1;
					} else {
						auto first = &binEdges[0];
						auto upper = std::upper_bound(first, first + binEdges.size(), value);
						binIndex = static_cast<int>(upper - first) - 1;
					}
					++result.counts[binIndex];
				}

				for (int i = 0; i < numBins; ++i) {
					Real width = binEdges[i + 1] - binEdges[i];
					if (!std::isfinite(width) || width <= 0.0)
						throw StatisticsError("ComputeHistogram: Bin width is not representable");
					result.frequencies[i] = static_cast<Real>(result.counts[i]) / static_cast<Real>(data.size());
					result.density[i] = result.frequencies[i] / width;
					if (!std::isfinite(result.density[i]))
						throw StatisticsError("ComputeHistogram: Density overflow");
				}

				return result;
			}

			inline HistogramResult ComputeHistogramAuto(
				const Vector<Real>& data,
				BinningMethod method = BinningMethod::Sturges)
			{
				return ComputeHistogram(data, GetBinCount(data, method));
			}

			inline HistogramResult HistogramSturges(const Vector<Real>& data)
			{
				return ComputeHistogramAuto(data, BinningMethod::Sturges);
			}

			inline HistogramResult HistogramScott(const Vector<Real>& data)
			{
				return ComputeHistogramAuto(data, BinningMethod::Scott);
			}

			inline HistogramResult HistogramFD(const Vector<Real>& data)
			{
				return ComputeHistogramAuto(data, BinningMethod::FreedmanDiaconis);
			}

			inline HistogramResult HistogramSqrt(const Vector<Real>& data)
			{
				return ComputeHistogramAuto(data, BinningMethod::SquareRoot);
			}

			inline HistogramResult HistogramRice(const Vector<Real>& data)
			{
				return ComputeHistogramAuto(data, BinningMethod::Rice);
			}

			inline FrequencyTableResult FrequencyTable(const Vector<Real>& data)
			{
				Detail::ValidateBinningData(data, "FrequencyTable");

				FrequencyTableResult result;
				result.totalCount = data.size();

				for (int i = 0; i < data.size(); ++i)
					++result.counts[data[i]];

				for (const auto& pair : result.counts)
					result.frequencies[pair.first] = static_cast<Real>(pair.second) / static_cast<Real>(data.size());

				result.uniqueCount = static_cast<int>(result.counts.size());
				return result;
			}

			inline ECDFResult EmpiricalCDF(const Vector<Real>& data)
			{
				Detail::ValidateBinningData(data, "EmpiricalCDF");

				std::vector<Real> sorted(data.size());
				for (int i = 0; i < data.size(); ++i)
					sorted[i] = data[i];
				std::sort(sorted.begin(), sorted.end());

				std::vector<Real> uniqueValues;
				std::vector<Real> cumulativeValues;
				uniqueValues.reserve(sorted.size());
				cumulativeValues.reserve(sorted.size());

				for (int i = 1; i <= data.size(); ++i) {
					if (i == data.size() || sorted[i] != sorted[i - 1]) {
						uniqueValues.push_back(sorted[i - 1]);
						cumulativeValues.push_back(static_cast<Real>(i) / static_cast<Real>(data.size()));
					}
				}

				ECDFResult result;
				result.n = data.size();
				result.x = Vector<Real>(static_cast<int>(uniqueValues.size()));
				result.cdf = Vector<Real>(static_cast<int>(cumulativeValues.size()));
				for (int i = 0; i < result.x.size(); ++i) {
					result.x[i] = uniqueValues[i];
					result.cdf[i] = cumulativeValues[i];
				}
				return result;
			}

			inline Real EvaluateECDF(const Vector<Real>& data, Real x)
			{
				Detail::ValidateBinningData(data, "EvaluateECDF");
				if (std::isnan(x))
					throw StatisticsError("EvaluateECDF: Evaluation point cannot be NaN");

				int count = 0;
				for (int i = 0; i < data.size(); ++i) {
					if (data[i] <= x)
						++count;
				}
				return static_cast<Real>(count) / static_cast<Real>(data.size());
			}

			inline Vector<Real> Quantiles(const Vector<Real>& data, const Vector<Real>& probabilities)
			{
				Detail::ValidateBinningData(data, "Quantiles");
				if (probabilities.size() == 0)
					throw StatisticsError("Quantiles: Probabilities vector cannot be empty");

				for (int i = 0; i < probabilities.size(); ++i) {
					if (!std::isfinite(probabilities[i]) || probabilities[i] < 0.0 || probabilities[i] > 1.0)
						throw StatisticsError("Quantiles: All probabilities must be finite and in [0, 1]");
				}

				std::vector<Real> sorted(data.size());
				for (int i = 0; i < data.size(); ++i)
					sorted[i] = data[i];
				std::sort(sorted.begin(), sorted.end());

				Vector<Real> result(probabilities.size());
				for (int i = 0; i < probabilities.size(); ++i) {
					Real rank = probabilities[i] * static_cast<Real>(data.size() - 1);
					int lower = static_cast<int>(std::floor(rank));
					int upper = static_cast<int>(std::ceil(rank));
					Real fraction = rank - static_cast<Real>(lower);
					result[i] = std::lerp(sorted[lower], sorted[upper], fraction);
				}
				return result;
			}

			inline Real Quantile(const Vector<Real>& data, Real p)
			{
				Vector<Real> probabilities({p});
				return Quantiles(data, probabilities)[0];
			}

			/// @brief Count how many samples fall into each bin defined by binEdges (numBins = size-1).
			inline Vector<int> BinCount(const Vector<Real>& data, const Vector<Real>& binEdges)
			{
				int n = data.size();
				int numBins = binEdges.size() - 1;
				if (n == 0)
					throw StatisticsError("BinCount: Input data cannot be empty");
				if (numBins < 1)
					throw StatisticsError("BinCount: Need at least 2 bin edges");

				Vector<int> counts(numBins);
				for (int i = 0; i < numBins; ++i)
					counts[i] = 0;

				for (int i = 0; i < n; ++i) {
					Real val = data[i];
					for (int b = 0; b < numBins; ++b) {
						if (val >= binEdges[b] && (val < binEdges[b + 1] || (b == numBins - 1 && val <= binEdges[b + 1]))) {
							++counts[b];
							break;
						}
					}
				}
				return counts;
			}

			/// @brief Return the bin index for each sample (-1 below range, numBins above range).
			inline Vector<int> Digitize(const Vector<Real>& data, const Vector<Real>& binEdges)
			{
				int n = data.size();
				int numBins = binEdges.size() - 1;
				if (n == 0)
					throw StatisticsError("Digitize: Input data cannot be empty");
				if (numBins < 1)
					throw StatisticsError("Digitize: Need at least 2 bin edges");

				Vector<int> indices(n);
				for (int i = 0; i < n; ++i) {
					Real val = data[i];
					if (val < binEdges[0]) { indices[i] = -1; continue; }
					if (val > binEdges[numBins]) { indices[i] = numBins; continue; }

					bool found = false;
					for (int b = 0; b < numBins; ++b) {
						if (val >= binEdges[b] && (val < binEdges[b + 1] || (b == numBins - 1 && val <= binEdges[b + 1]))) {
							indices[i] = b;
							found = true;
							break;
						}
					}
					if (!found)
						indices[i] = numBins - 1;
				}
				return indices;
			}

			/// @brief Create numBins+1 uniformly spaced bin edges over [minVal, maxVal].
			inline Vector<Real> CreateUniformBinEdges(Real minVal, Real maxVal, int numBins)
			{
				if (numBins < 1)
					throw StatisticsError("CreateUniformBinEdges: Number of bins must be at least 1");

				Vector<Real> edges(numBins + 1);
				Real width = (maxVal - minVal) / numBins;
				for (int i = 0; i <= numBins; ++i)
					edges[i] = minVal + i * width;
				return edges;
			}

			/// @brief Create numBins+1 logarithmically spaced bin edges over [minVal, maxVal] (minVal>0).
			inline Vector<Real> CreateLogBinEdges(Real minVal, Real maxVal, int numBins)
			{
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
				for (int i = 0; i <= numBins; ++i)
					edges[i] = std::pow(REAL(10.0), logMin + i * logWidth);
				return edges;
			}

		} // namespace Histogram
	} // namespace Statistics
} // namespace MML

#endif // MML_STATISTICS_HISTOGRAM_H
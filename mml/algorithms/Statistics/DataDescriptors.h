///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DataDescriptors.h                                                   ///
///  Description: Typed data analysis classes beyond Real numbers                     ///
///               BoolDataStats, TimeDataStats, DateDataStats, GeoDataStats,         ///
///               CategoryStats                                                       ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DATA_DESCRIPTORS_H
#define MML_DATA_DESCRIPTORS_H

#include "MMLBase.h"
#include "MMLExceptions.h"
#include "base/Vector.h"

#include <cmath>
#include <map>
#include <string>
#include <vector>
#include <algorithm>
#include <numeric>
#include <ctime>
#include <sstream>
#include <iomanip>

namespace MML
{
namespace Statistics
{

/////////////////////////////////////////////////////////////////////////////////////
///                          BOOLEAN DATA STATISTICS                               ///
/////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Statistics for boolean data
 * 
 * Computes proportion-based statistics including confidence intervals
 * and chi-square test for expected proportion.
 */
struct BoolDataStats
{
	size_t count;           ///< Total count
	size_t countTrue;       ///< Count of true values
	size_t countFalse;      ///< Count of false values
	Real proportionTrue;    ///< Proportion of true values (p = countTrue/count)
	Real standardError;     ///< Standard error of proportion: sqrt(p(1-p)/n)
	Real ciLower95;         ///< 95% confidence interval lower bound
	Real ciUpper95;         ///< 95% confidence interval upper bound

	/**
	 * @brief Compute statistics from boolean data
	 * @param data Vector of boolean values
	 * @return BoolDataStats with computed values
	 * @throws StatisticsError if data is empty
	 */
	static BoolDataStats Compute(const std::vector<bool>& data)
	{
		if (data.empty())
			throw StatisticsError("BoolDataStats: Data cannot be empty");

		BoolDataStats stats;
		stats.count = data.size();
		stats.countTrue = std::count(data.begin(), data.end(), true);
		stats.countFalse = stats.count - stats.countTrue;
		
		stats.proportionTrue = static_cast<Real>(stats.countTrue) / stats.count;
		
		// Standard error of proportion: sqrt(p(1-p)/n)
		Real p = stats.proportionTrue;
		Real n = static_cast<Real>(stats.count);
		stats.standardError = std::sqrt(p * (1.0 - p) / n);
		
		// 95% CI using normal approximation (z = 1.96)
		Real z = 1.96;
		stats.ciLower95 = std::max(0.0, p - z * stats.standardError);
		stats.ciUpper95 = std::min(1.0, p + z * stats.standardError);
		
		return stats;
	}

	/**
	 * @brief Chi-square test statistic for expected proportion
	 * @param expectedProportion The hypothesized proportion of true values
	 * @return Chi-square statistic value
	 * 
	 * Formula: chi2 = sum((O-E)^2/E) for both categories
	 * Under H0, this follows chi-square distribution with df=1
	 */
	Real ChiSquare(Real expectedProportion) const
	{
		if (expectedProportion <= 0.0 || expectedProportion >= 1.0)
			throw StatisticsError("BoolDataStats::ChiSquare: Expected proportion must be in (0, 1)");

		Real expectedTrue = count * expectedProportion;
		Real expectedFalse = count * (1.0 - expectedProportion);
		
		Real chi2 = (countTrue - expectedTrue) * (countTrue - expectedTrue) / expectedTrue
		          + (countFalse - expectedFalse) * (countFalse - expectedFalse) / expectedFalse;
		
		return chi2;
	}

	/// @brief Print summary to stream
	void PrintSummary(std::ostream& os = std::cout) const
	{
		os << "Boolean Data Statistics:\n";
		os << "  Count: " << count << " (True: " << countTrue << ", False: " << countFalse << ")\n";
		os << "  Proportion True: " << std::fixed << std::setprecision(4) << proportionTrue << "\n";
		os << "  Standard Error: " << standardError << "\n";
		os << "  95% CI: [" << ciLower95 << ", " << ciUpper95 << "]\n";
	}
};

/////////////////////////////////////////////////////////////////////////////////////
///                           TIME DATA STATISTICS                                 ///
/////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Statistics for time/duration data (seconds since epoch or duration in seconds)
 * 
 * Supports both timestamps and durations. For timestamps, provides
 * distributions by hour of day and day of week.
 */
struct TimeDataStats
{
	size_t count;           ///< Number of time values
	Real minSeconds;        ///< Minimum time in seconds
	Real maxSeconds;        ///< Maximum time in seconds
	Real avgSeconds;        ///< Average time in seconds
	Real medianSeconds;     ///< Median time in seconds
	Real stdDevSeconds;     ///< Standard deviation in seconds
	Real rangeSeconds;      ///< Range (max - min) in seconds

	std::vector<size_t> hourDistribution;    ///< Count by hour (0-23)
	std::vector<size_t> dayOfWeekDistribution; ///< Count by day (0=Sun, 6=Sat)

	/**
	 * @brief Compute statistics from time data (seconds)
	 * @param data Vector of time values in seconds (epoch or duration)
	 * @param isTimestamp If true, compute hour/day distributions (default: true)
	 * @return TimeDataStats with computed values
	 */
	static TimeDataStats Compute(const Vector<Real>& data, bool isTimestamp = true)
	{
		if (data.size() == 0)
			throw StatisticsError("TimeDataStats: Data cannot be empty");

		TimeDataStats stats;
		stats.count = data.size();

		// Basic statistics
		stats.minSeconds = data[0];
		stats.maxSeconds = data[0];
		Real sum = 0.0;
		
		for (int i = 0; i < data.size(); ++i) {
			if (data[i] < stats.minSeconds) stats.minSeconds = data[i];
			if (data[i] > stats.maxSeconds) stats.maxSeconds = data[i];
			sum += data[i];
		}
		
		stats.avgSeconds = sum / stats.count;
		stats.rangeSeconds = stats.maxSeconds - stats.minSeconds;

		// Median
		std::vector<Real> sorted(data.size());
		for (int i = 0; i < data.size(); ++i)
			sorted[i] = data[i];
		std::sort(sorted.begin(), sorted.end());
		
		if (stats.count % 2 == 0)
			stats.medianSeconds = (sorted[stats.count / 2 - 1] + sorted[stats.count / 2]) / 2.0;
		else
			stats.medianSeconds = sorted[stats.count / 2];

		// Standard deviation
		Real sumSq = 0.0;
		for (int i = 0; i < data.size(); ++i) {
			Real diff = data[i] - stats.avgSeconds;
			sumSq += diff * diff;
		}
		stats.stdDevSeconds = std::sqrt(sumSq / (stats.count - 1));

		// Hour and day distributions (only for timestamps)
		stats.hourDistribution.resize(24, 0);
		stats.dayOfWeekDistribution.resize(7, 0);
		
		if (isTimestamp) {
			for (int i = 0; i < data.size(); ++i) {
				time_t t = static_cast<time_t>(data[i]);
				std::tm* tm = std::gmtime(&t);
				if (tm) {
					stats.hourDistribution[tm->tm_hour]++;
					stats.dayOfWeekDistribution[tm->tm_wday]++;
				}
			}
		}

		return stats;
	}

	/// @brief Convert seconds to human-readable HH:MM:SS format
	static std::string SecondsToHMS(Real seconds)
	{
		bool negative = seconds < 0;
		if (negative) seconds = -seconds;
		
		int h = static_cast<int>(seconds / 3600);
		int m = static_cast<int>(std::fmod(seconds, 3600) / 60);
		int s = static_cast<int>(std::fmod(seconds, 60));
		
		std::ostringstream oss;
		if (negative) oss << "-";
		oss << std::setfill('0') << std::setw(2) << h << ":"
		    << std::setw(2) << m << ":" << std::setw(2) << s;
		return oss.str();
	}

	/// @brief Convert seconds to days (with decimal)
	static Real SecondsToDays(Real seconds) { return seconds / 86400.0; }

	/// @brief Print summary to stream
	void PrintSummary(std::ostream& os = std::cout) const
	{
		os << "Time Data Statistics:\n";
		os << "  Count: " << count << "\n";
		os << "  Range: " << SecondsToHMS(rangeSeconds) << " (" << SecondsToDays(rangeSeconds) << " days)\n";
		os << "  Min: " << SecondsToHMS(minSeconds) << "\n";
		os << "  Max: " << SecondsToHMS(maxSeconds) << "\n";
		os << "  Avg: " << SecondsToHMS(avgSeconds) << "\n";
		os << "  Median: " << SecondsToHMS(medianSeconds) << "\n";
		os << "  Std Dev: " << SecondsToHMS(stdDevSeconds) << "\n";
	}
};

/////////////////////////////////////////////////////////////////////////////////////
///                           DATE DATA STATISTICS                                 ///
/////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Statistics for date data (days since epoch or structured dates)
 * 
 * Provides year/month/day-of-week frequency distributions and
 * weekday vs weekend analysis.
 */
struct DateDataStats
{
	size_t count;           ///< Number of date values
	int minYear;            ///< Minimum year
	int maxYear;            ///< Maximum year
	int spanDays;           ///< Span in days (max - min)
	
	std::map<int, size_t> yearDistribution;      ///< Count by year
	std::vector<size_t> monthDistribution;       ///< Count by month (0=Jan, 11=Dec)
	std::vector<size_t> dayOfWeekDistribution;   ///< Count by day (0=Sun, 6=Sat)
	size_t weekdayCount;    ///< Mon-Fri count
	size_t weekendCount;    ///< Sat-Sun count

	/**
	 * @brief Compute statistics from date data (epoch days or timestamps)
	 * @param epochDays Vector of days since epoch (or timestamps in seconds)
	 * @param isSeconds If true, treat input as seconds; otherwise days (default: false)
	 * @return DateDataStats with computed values
	 */
	static DateDataStats Compute(const Vector<Real>& epochDays, bool isSeconds = false)
	{
		if (epochDays.size() == 0)
			throw StatisticsError("DateDataStats: Data cannot be empty");

		DateDataStats stats;
		stats.count = epochDays.size();
		stats.monthDistribution.resize(12, 0);
		stats.dayOfWeekDistribution.resize(7, 0);
		stats.weekdayCount = 0;
		stats.weekendCount = 0;

		Real minVal = epochDays[0];
		Real maxVal = epochDays[0];

		for (int i = 0; i < epochDays.size(); ++i) {
			Real val = epochDays[i];
			if (val < minVal) minVal = val;
			if (val > maxVal) maxVal = val;

			// Convert to time_t
			time_t t;
			if (isSeconds)
				t = static_cast<time_t>(val);
			else
				t = static_cast<time_t>(val * 86400.0);

			std::tm* tm = std::gmtime(&t);
			if (tm) {
				int year = tm->tm_year + 1900;
				stats.yearDistribution[year]++;
				stats.monthDistribution[tm->tm_mon]++;
				stats.dayOfWeekDistribution[tm->tm_wday]++;
				
				if (tm->tm_wday == 0 || tm->tm_wday == 6)
					stats.weekendCount++;
				else
					stats.weekdayCount++;
			}
		}

		// Compute year range
		if (!stats.yearDistribution.empty()) {
			stats.minYear = stats.yearDistribution.begin()->first;
			stats.maxYear = stats.yearDistribution.rbegin()->first;
		}

		// Span in days
		if (isSeconds)
			stats.spanDays = static_cast<int>((maxVal - minVal) / 86400.0);
		else
			stats.spanDays = static_cast<int>(maxVal - minVal);

		return stats;
	}

	/// @brief Get most common month (0-11)
	int GetModeMonth() const
	{
		return static_cast<int>(std::max_element(monthDistribution.begin(), 
		                                          monthDistribution.end()) - monthDistribution.begin());
	}

	/// @brief Get weekday proportion (Mon-Fri / total)
	Real GetWeekdayProportion() const
	{
		return static_cast<Real>(weekdayCount) / count;
	}

	/// @brief Print summary to stream
	void PrintSummary(std::ostream& os = std::cout) const
	{
		os << "Date Data Statistics:\n";
		os << "  Count: " << count << "\n";
		os << "  Year Range: " << minYear << " - " << maxYear << "\n";
		os << "  Span: " << spanDays << " days\n";
		os << "  Weekday/Weekend: " << weekdayCount << "/" << weekendCount 
		   << " (" << std::fixed << std::setprecision(1) << (GetWeekdayProportion() * 100) << "% weekday)\n";
		
		os << "  Month Distribution: ";
		const char* months[] = {"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"};
		for (int i = 0; i < 12; ++i) {
			if (monthDistribution[i] > 0)
				os << months[i] << ":" << monthDistribution[i] << " ";
		}
		os << "\n";
	}
};

/////////////////////////////////////////////////////////////////////////////////////
///                      GEOGRAPHIC DATA STATISTICS                                ///
/////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Statistics for geographic (lat/lon) data
 * 
 * Computes centroid, bounding box, and average distance from centroid.
 * Uses WGS84 coordinates (standard GPS).
 */
struct GeoDataStats
{
	size_t count;           ///< Number of points
	Real centroidLat;       ///< Centroid latitude
	Real centroidLon;       ///< Centroid longitude
	Real minLat, maxLat;    ///< Latitude bounding box
	Real minLon, maxLon;    ///< Longitude bounding box
	Real avgDistanceFromCentroid;  ///< Average distance from centroid (km)
	Real maxDistanceFromCentroid;  ///< Maximum distance from centroid (km)

	/// @brief Earth radius in kilometers (WGS84 mean)
	static constexpr Real EARTH_RADIUS_KM = 6371.0;

	/**
	 * @brief Compute Haversine distance between two points
	 * @param lat1, lon1 First point (degrees)
	 * @param lat2, lon2 Second point (degrees)
	 * @return Distance in kilometers
	 */
	static Real HaversineDistance(Real lat1, Real lon1, Real lat2, Real lon2)
	{
		// Convert to radians
		Real lat1Rad = lat1 * Constants::PI / 180.0;
		Real lat2Rad = lat2 * Constants::PI / 180.0;
		Real dLat = (lat2 - lat1) * Constants::PI / 180.0;
		Real dLon = (lon2 - lon1) * Constants::PI / 180.0;

		Real a = std::sin(dLat / 2) * std::sin(dLat / 2) +
		         std::cos(lat1Rad) * std::cos(lat2Rad) *
		         std::sin(dLon / 2) * std::sin(dLon / 2);
		Real c = 2 * std::atan2(std::sqrt(a), std::sqrt(1 - a));

		return EARTH_RADIUS_KM * c;
	}

	/**
	 * @brief Compute statistics from geographic data
	 * @param latitudes Vector of latitudes (degrees, -90 to 90)
	 * @param longitudes Vector of longitudes (degrees, -180 to 180)
	 * @return GeoDataStats with computed values
	 * @throws StatisticsError if vectors have different sizes or are empty
	 */
	static GeoDataStats Compute(const Vector<Real>& latitudes, const Vector<Real>& longitudes)
	{
		if (latitudes.size() == 0 || longitudes.size() == 0)
			throw StatisticsError("GeoDataStats: Data cannot be empty");
		if (latitudes.size() != longitudes.size())
			throw StatisticsError("GeoDataStats: Latitude and longitude vectors must have same size");

		GeoDataStats stats;
		stats.count = latitudes.size();

		// Compute centroid (simple average - works for small regions)
		Real sumLat = 0, sumLon = 0;
		stats.minLat = stats.maxLat = latitudes[0];
		stats.minLon = stats.maxLon = longitudes[0];

		for (int i = 0; i < latitudes.size(); ++i) {
			sumLat += latitudes[i];
			sumLon += longitudes[i];
			
			if (latitudes[i] < stats.minLat) stats.minLat = latitudes[i];
			if (latitudes[i] > stats.maxLat) stats.maxLat = latitudes[i];
			if (longitudes[i] < stats.minLon) stats.minLon = longitudes[i];
			if (longitudes[i] > stats.maxLon) stats.maxLon = longitudes[i];
		}

		stats.centroidLat = sumLat / stats.count;
		stats.centroidLon = sumLon / stats.count;

		// Compute distances from centroid
		Real sumDist = 0;
		stats.maxDistanceFromCentroid = 0;

		for (int i = 0; i < latitudes.size(); ++i) {
			Real dist = HaversineDistance(latitudes[i], longitudes[i],
			                               stats.centroidLat, stats.centroidLon);
			sumDist += dist;
			if (dist > stats.maxDistanceFromCentroid)
				stats.maxDistanceFromCentroid = dist;
		}

		stats.avgDistanceFromCentroid = sumDist / stats.count;

		return stats;
	}

	/// @brief Get bounding box diagonal distance (km)
	Real GetBoundingBoxDiagonal() const
	{
		return HaversineDistance(minLat, minLon, maxLat, maxLon);
	}

	/// @brief Print summary to stream
	void PrintSummary(std::ostream& os = std::cout) const
	{
		os << "Geographic Data Statistics:\n";
		os << "  Count: " << count << " points\n";
		os << std::fixed << std::setprecision(4);
		os << "  Centroid: (" << centroidLat << ", " << centroidLon << ")\n";
		os << "  Bounding Box: [" << minLat << ", " << minLon << "] to [" 
		   << maxLat << ", " << maxLon << "]\n";
		os << std::setprecision(2);
		os << "  Avg Distance from Centroid: " << avgDistanceFromCentroid << " km\n";
		os << "  Max Distance from Centroid: " << maxDistanceFromCentroid << " km\n";
		os << "  Bounding Box Diagonal: " << GetBoundingBoxDiagonal() << " km\n";
	}
};

/////////////////////////////////////////////////////////////////////////////////////
///                        CATEGORICAL DATA STATISTICS                             ///
/////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Statistics for categorical (string) data
 * 
 * Computes frequency distribution, unique count, mode, and entropy.
 */
struct CategoryStats
{
	size_t count;           ///< Total number of values
	size_t uniqueCount;     ///< Number of unique categories
	std::string mode;       ///< Most frequent category
	size_t modeCount;       ///< Count of most frequent category
	Real entropy;           ///< Shannon entropy (bits)
	
	std::map<std::string, size_t> frequencies;  ///< Category -> count

	/**
	 * @brief Compute statistics from categorical data
	 * @param data Vector of string categories
	 * @return CategoryStats with computed values
	 */
	static CategoryStats Compute(const std::vector<std::string>& data)
	{
		if (data.empty())
			throw StatisticsError("CategoryStats: Data cannot be empty");

		CategoryStats stats;
		stats.count = data.size();
		stats.modeCount = 0;

		// Build frequency map
		for (const auto& val : data) {
			stats.frequencies[val]++;
		}

		stats.uniqueCount = stats.frequencies.size();

		// Find mode and compute entropy
		stats.entropy = 0.0;
		Real n = static_cast<Real>(stats.count);

		for (const auto& [category, freq] : stats.frequencies) {
			if (freq > stats.modeCount) {
				stats.modeCount = freq;
				stats.mode = category;
			}
			
			// Shannon entropy: -sum(p * log2(p))
			Real p = freq / n;
			if (p > 0)
				stats.entropy -= p * std::log2(p);
		}

		return stats;
	}

	/**
	 * @brief Get top N categories by frequency
	 * @param n Number of categories to return
	 * @return Vector of (category, count) pairs sorted by count descending
	 */
	std::vector<std::pair<std::string, size_t>> GetTopN(size_t n) const
	{
		std::vector<std::pair<std::string, size_t>> sorted(frequencies.begin(), frequencies.end());
		std::sort(sorted.begin(), sorted.end(),
		          [](const auto& a, const auto& b) { return a.second > b.second; });
		
		if (sorted.size() > n)
			sorted.resize(n);
		
		return sorted;
	}

	/// @brief Get proportion of mode (mode_count / total)
	Real GetModeProportion() const
	{
		return static_cast<Real>(modeCount) / count;
	}

	/// @brief Print summary to stream
	void PrintSummary(std::ostream& os = std::cout) const
	{
		os << "Categorical Data Statistics:\n";
		os << "  Count: " << count << "\n";
		os << "  Unique Categories: " << uniqueCount << "\n";
		os << "  Mode: \"" << mode << "\" (count: " << modeCount 
		   << ", " << std::fixed << std::setprecision(1) << (GetModeProportion() * 100) << "%)\n";
		os << "  Entropy: " << std::setprecision(3) << entropy << " bits\n";
		
		os << "  Top 5 Categories:\n";
		auto top5 = GetTopN(5);
		for (const auto& [cat, freq] : top5) {
			os << "    \"" << cat << "\": " << freq 
			   << " (" << std::setprecision(1) << (100.0 * freq / count) << "%)\n";
		}
	}
};

/////////////////////////////////////////////////////////////////////////////////////
///                         COMBINED DATA DESCRIPTOR                               ///
/////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Unified interface for computing type-appropriate statistics
 * 
 * Provides static methods to compute stats based on detected or specified data type.
 */
struct DataDescriptor
{
	/// @brief Data type enumeration
	enum class Type { REAL, BOOL, TIME, DATE, GEO, CATEGORY };

	/// @brief Detect likely type from string data sample
	static Type DetectType(const std::vector<std::string>& samples)
	{
		if (samples.empty()) return Type::CATEGORY;

		// Check patterns in first few non-empty samples
		int numericCount = 0;
		int boolCount = 0;
		int dateCount = 0;
		int timeCount = 0;
		int checked = 0;

		for (const auto& s : samples) {
			if (s.empty()) continue;
			if (++checked > 100) break;  // Sample first 100

			// Boolean check
			std::string lower = s;
			std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
			if (lower == "true" || lower == "false" || lower == "yes" || lower == "no") {
				boolCount++;
				continue;
			}

			// Date check (YYYY-MM-DD or similar)
			if (s.length() >= 8 && s.length() <= 10 && 
			    (s[4] == '-' || s[4] == '/') && std::isdigit(s[0])) {
				dateCount++;
				continue;
			}

			// Time check (HH:MM:SS or similar)
			if (s.length() >= 5 && s.length() <= 8 && s[2] == ':') {
				timeCount++;
				continue;
			}

			// Numeric check
			try {
				(void)std::stod(s);  // Discard result, only checking if valid
				numericCount++;
			} catch (...) {
				// Not numeric - likely category
			}
		}

		// Determine type based on majority
		int threshold = checked / 2;
		if (boolCount > threshold) return Type::BOOL;
		if (dateCount > threshold) return Type::DATE;
		if (timeCount > threshold) return Type::TIME;
		if (numericCount > threshold) return Type::REAL;
		
		return Type::CATEGORY;
	}
};

} // namespace Statistics
} // namespace MML

#endif // MML_DATA_DESCRIPTORS_H

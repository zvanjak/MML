///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  Tests for DataDescriptors.h - Typed data statistics                              ///
///////////////////////////////////////////////////////////////////////////////////////////

#include <catch2/catch_all.hpp>
#include <catch2/catch_approx.hpp>

#include "MMLBase.h"
#include "base/Vector.h"
#include "algorithms/Statistics/DataDescriptors.h"

#include <cmath>
#include <vector>
#include <string>
#include <sstream>

using namespace MML;
using namespace MML::Statistics;

// Test precision helper
#define TEST_PRECISION_INFO() INFO("Test precision: " << std::numeric_limits<Real>::epsilon() << " (double)")

/////////////////////////////////////////////////////////////////////////////////////
///                         BOOLEAN DATA STATS TESTS                               ///
/////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("BoolDataStats - Basic statistics", "[datadescriptors][bool]") {
	TEST_PRECISION_INFO();

	SECTION("All true") {
		std::vector<bool> data = {true, true, true, true, true};
		auto stats = BoolDataStats::Compute(data);
		
		REQUIRE(stats.count == 5);
		REQUIRE(stats.countTrue == 5);
		REQUIRE(stats.countFalse == 0);
		REQUIRE(stats.proportionTrue == Catch::Approx(1.0));
		REQUIRE(stats.standardError == Catch::Approx(0.0));
	}

	SECTION("All false") {
		std::vector<bool> data = {false, false, false, false};
		auto stats = BoolDataStats::Compute(data);
		
		REQUIRE(stats.count == 4);
		REQUIRE(stats.countTrue == 0);
		REQUIRE(stats.countFalse == 4);
		REQUIRE(stats.proportionTrue == Catch::Approx(0.0));
	}

	SECTION("Mixed data") {
		std::vector<bool> data = {true, false, true, true, false, false, true, false, true, false};
		auto stats = BoolDataStats::Compute(data);
		
		REQUIRE(stats.count == 10);
		REQUIRE(stats.countTrue == 5);
		REQUIRE(stats.countFalse == 5);
		REQUIRE(stats.proportionTrue == Catch::Approx(0.5));
		
		// SE = sqrt(0.5 * 0.5 / 10) = sqrt(0.025) ≈ 0.158
		REQUIRE(stats.standardError == Catch::Approx(0.158113883).epsilon(0.001));
	}

	SECTION("Confidence interval bounds") {
		std::vector<bool> data(100);
		std::fill(data.begin(), data.begin() + 70, true);  // 70% true
		
		auto stats = BoolDataStats::Compute(data);
		
		REQUIRE(stats.proportionTrue == Catch::Approx(0.7));
		REQUIRE(stats.ciLower95 < stats.proportionTrue);
		REQUIRE(stats.ciUpper95 > stats.proportionTrue);
		REQUIRE(stats.ciLower95 >= 0.0);
		REQUIRE(stats.ciUpper95 <= 1.0);
	}

	SECTION("Chi-square test") {
		std::vector<bool> data(100);
		std::fill(data.begin(), data.begin() + 60, true);  // 60 true, 40 false
		
		auto stats = BoolDataStats::Compute(data);
		
		// Test against 50% expectation
		// Expected: 50 true, 50 false
		// Chi2 = (60-50)^2/50 + (40-50)^2/50 = 100/50 + 100/50 = 4.0
		Real chi2 = stats.ChiSquare(0.5);
		REQUIRE(chi2 == Catch::Approx(4.0));
		
		// Test against actual proportion (should be 0)
		Real chi2_actual = stats.ChiSquare(0.6);
		REQUIRE(chi2_actual == Catch::Approx(0.0).margin(0.001));
	}

	SECTION("Empty data throws") {
		std::vector<bool> empty;
		REQUIRE_THROWS_AS(BoolDataStats::Compute(empty), StatisticsError);
	}
}

/////////////////////////////////////////////////////////////////////////////////////
///                          TIME DATA STATS TESTS                                 ///
/////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("TimeDataStats - Basic statistics", "[datadescriptors][time]") {
	TEST_PRECISION_INFO();

	SECTION("Duration statistics") {
		// Durations in seconds: 60, 120, 180, 240, 300 (1-5 minutes)
		Vector<Real> durations(5);
		durations[0] = 60;
		durations[1] = 120;
		durations[2] = 180;
		durations[3] = 240;
		durations[4] = 300;
		
		auto stats = TimeDataStats::Compute(durations, false);  // Not timestamps
		
		REQUIRE(stats.count == 5);
		REQUIRE(stats.minSeconds == Catch::Approx(60.0));
		REQUIRE(stats.maxSeconds == Catch::Approx(300.0));
		REQUIRE(stats.avgSeconds == Catch::Approx(180.0));  // (60+120+180+240+300)/5 = 180
		REQUIRE(stats.medianSeconds == Catch::Approx(180.0));
		REQUIRE(stats.rangeSeconds == Catch::Approx(240.0));  // 300 - 60
	}

	SECTION("SecondsToHMS conversion") {
		REQUIRE(TimeDataStats::SecondsToHMS(0) == "00:00:00");
		REQUIRE(TimeDataStats::SecondsToHMS(61) == "00:01:01");
		REQUIRE(TimeDataStats::SecondsToHMS(3661) == "01:01:01");
		REQUIRE(TimeDataStats::SecondsToHMS(86400) == "24:00:00");  // 1 day
		REQUIRE(TimeDataStats::SecondsToHMS(-3600) == "-01:00:00");
	}

	SECTION("SecondsToDays conversion") {
		REQUIRE(TimeDataStats::SecondsToDays(86400) == Catch::Approx(1.0));
		REQUIRE(TimeDataStats::SecondsToDays(43200) == Catch::Approx(0.5));
		REQUIRE(TimeDataStats::SecondsToDays(604800) == Catch::Approx(7.0));  // 1 week
	}

	SECTION("Timestamp hour distribution") {
		// Create timestamps at different hours
		// Using epoch: Jan 1, 1970 00:00:00 UTC
		Vector<Real> timestamps(3);
		timestamps[0] = 0;           // 00:00 UTC
		timestamps[1] = 3600;        // 01:00 UTC
		timestamps[2] = 7200;        // 02:00 UTC
		
		auto stats = TimeDataStats::Compute(timestamps, true);
		
		REQUIRE(stats.hourDistribution[0] == 1);
		REQUIRE(stats.hourDistribution[1] == 1);
		REQUIRE(stats.hourDistribution[2] == 1);
	}

	SECTION("Empty data throws") {
		Vector<Real> empty(0);
		REQUIRE_THROWS_AS(TimeDataStats::Compute(empty), StatisticsError);
	}
}

/////////////////////////////////////////////////////////////////////////////////////
///                          DATE DATA STATS TESTS                                 ///
/////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("DateDataStats - Basic statistics", "[datadescriptors][date]") {
	TEST_PRECISION_INFO();

	SECTION("Date range statistics") {
		// Days since epoch for year 2020 (approx 18262) and 2024 (approx 19724)
		Vector<Real> epochDays(3);
		epochDays[0] = 18262;  // ~2020
		epochDays[1] = 18627;  // ~2021
		epochDays[2] = 19724;  // ~2024
		
		auto stats = DateDataStats::Compute(epochDays, false);
		
		REQUIRE(stats.count == 3);
		REQUIRE(stats.spanDays == 19724 - 18262);
	}

	SECTION("From timestamps") {
		// Unix timestamps for 2024 dates
		Vector<Real> timestamps(2);
		timestamps[0] = 1704067200;  // Jan 1, 2024 00:00:00 UTC
		timestamps[1] = 1706745600;  // Feb 1, 2024 00:00:00 UTC
		
		auto stats = DateDataStats::Compute(timestamps, true);  // isSeconds = true
		
		REQUIRE(stats.count == 2);
		REQUIRE(stats.minYear == 2024);
		REQUIRE(stats.maxYear == 2024);
		REQUIRE(stats.monthDistribution[0] == 1);  // January
		REQUIRE(stats.monthDistribution[1] == 1);  // February
	}

	SECTION("Weekday vs weekend") {
		// Create dates that span a full week
		Vector<Real> timestamps(7);
		// Mon Jan 1, 2024 to Sun Jan 7, 2024
		for (int i = 0; i < 7; ++i) {
			timestamps[i] = 1704067200 + i * 86400;  // Add days
		}
		
		auto stats = DateDataStats::Compute(timestamps, true);
		
		// Jan 1, 2024 is Monday
		REQUIRE(stats.weekdayCount == 5);
		REQUIRE(stats.weekendCount == 2);
		REQUIRE(stats.GetWeekdayProportion() == Catch::Approx(5.0/7.0));
	}

	SECTION("Empty data throws") {
		Vector<Real> empty(0);
		REQUIRE_THROWS_AS(DateDataStats::Compute(empty), StatisticsError);
	}
}

/////////////////////////////////////////////////////////////////////////////////////
///                          GEO DATA STATS TESTS                                  ///
/////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("GeoDataStats - Basic statistics", "[datadescriptors][geo]") {
	TEST_PRECISION_INFO();

	SECTION("Haversine distance") {
		// New York to Los Angeles: ~3944 km
		Real nyLat = 40.7128, nyLon = -74.0060;
		Real laLat = 34.0522, laLon = -118.2437;
		
		Real dist = GeoDataStats::HaversineDistance(nyLat, nyLon, laLat, laLon);
		REQUIRE(dist == Catch::Approx(3935.75).epsilon(0.01));  // Allow 1% error
	}

	SECTION("Same point distance") {
		Real dist = GeoDataStats::HaversineDistance(45.0, 10.0, 45.0, 10.0);
		REQUIRE(dist == Catch::Approx(0.0).margin(0.001));
	}

	SECTION("Centroid calculation") {
		// Three points forming a triangle
		Vector<Real> lats(3), lons(3);
		lats[0] = 0.0;  lons[0] = 0.0;
		lats[1] = 0.0;  lons[1] = 2.0;
		lats[2] = 2.0;  lons[2] = 1.0;
		
		auto stats = GeoDataStats::Compute(lats, lons);
		
		REQUIRE(stats.count == 3);
		REQUIRE(stats.centroidLat == Catch::Approx(2.0/3.0));
		REQUIRE(stats.centroidLon == Catch::Approx(1.0));
	}

	SECTION("Bounding box") {
		Vector<Real> lats(4), lons(4);
		lats[0] = 10.0;  lons[0] = -5.0;
		lats[1] = 20.0;  lons[1] = 5.0;
		lats[2] = 15.0;  lons[2] = -10.0;
		lats[3] = 12.0;  lons[3] = 15.0;
		
		auto stats = GeoDataStats::Compute(lats, lons);
		
		REQUIRE(stats.minLat == Catch::Approx(10.0));
		REQUIRE(stats.maxLat == Catch::Approx(20.0));
		REQUIRE(stats.minLon == Catch::Approx(-10.0));
		REQUIRE(stats.maxLon == Catch::Approx(15.0));
	}

	SECTION("Average distance from centroid") {
		// Square around origin: corners at (±1, ±1)
		Vector<Real> lats(4), lons(4);
		lats[0] = 1.0;   lons[0] = 1.0;
		lats[1] = 1.0;   lons[1] = -1.0;
		lats[2] = -1.0;  lons[2] = 1.0;
		lats[3] = -1.0;  lons[3] = -1.0;
		
		auto stats = GeoDataStats::Compute(lats, lons);
		
		REQUIRE(stats.centroidLat == Catch::Approx(0.0));
		REQUIRE(stats.centroidLon == Catch::Approx(0.0));
		REQUIRE(stats.avgDistanceFromCentroid > 0);  // All points equidistant
	}

	SECTION("Mismatched sizes throws") {
		Vector<Real> lats(3), lons(2);
		lats[0] = 0; lats[1] = 1; lats[2] = 2;
		lons[0] = 0; lons[1] = 1;
		
		REQUIRE_THROWS_AS(GeoDataStats::Compute(lats, lons), StatisticsError);
	}

	SECTION("Empty data throws") {
		Vector<Real> empty(0);
		REQUIRE_THROWS_AS(GeoDataStats::Compute(empty, empty), StatisticsError);
	}
}

/////////////////////////////////////////////////////////////////////////////////////
///                       CATEGORY STATS TESTS                                     ///
/////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("CategoryStats - Basic statistics", "[datadescriptors][category]") {
	TEST_PRECISION_INFO();

	SECTION("Frequency counting") {
		std::vector<std::string> data = {"A", "B", "A", "C", "A", "B"};
		auto stats = CategoryStats::Compute(data);
		
		REQUIRE(stats.count == 6);
		REQUIRE(stats.uniqueCount == 3);
		REQUIRE(stats.frequencies["A"] == 3);
		REQUIRE(stats.frequencies["B"] == 2);
		REQUIRE(stats.frequencies["C"] == 1);
	}

	SECTION("Mode detection") {
		std::vector<std::string> data = {"red", "blue", "red", "green", "red", "blue"};
		auto stats = CategoryStats::Compute(data);
		
		REQUIRE(stats.mode == "red");
		REQUIRE(stats.modeCount == 3);
		REQUIRE(stats.GetModeProportion() == Catch::Approx(0.5));
	}

	SECTION("Entropy calculation") {
		// Uniform distribution: 2 categories, equal frequency
		// Entropy = -sum(0.5 * log2(0.5)) = -2 * (-0.5) = 1.0 bit
		std::vector<std::string> uniform = {"A", "B", "A", "B"};
		auto stats = CategoryStats::Compute(uniform);
		
		REQUIRE(stats.entropy == Catch::Approx(1.0));
	}

	SECTION("Entropy - single category") {
		// All same: entropy = 0
		std::vector<std::string> single = {"X", "X", "X", "X"};
		auto stats = CategoryStats::Compute(single);
		
		REQUIRE(stats.entropy == Catch::Approx(0.0));
		REQUIRE(stats.uniqueCount == 1);
	}

	SECTION("Top N categories") {
		std::vector<std::string> data = {"A", "A", "A", "B", "B", "C", "D", "D", "D", "D"};
		auto stats = CategoryStats::Compute(data);
		
		auto top2 = stats.GetTopN(2);
		REQUIRE(top2.size() == 2);
		REQUIRE(top2[0].first == "D");  // 4 occurrences
		REQUIRE(top2[0].second == 4);
		REQUIRE(top2[1].first == "A");  // 3 occurrences
		REQUIRE(top2[1].second == 3);
	}

	SECTION("Species example (Iris-like)") {
		std::vector<std::string> species = {
			"setosa", "setosa", "setosa",
			"versicolor", "versicolor", "versicolor",
			"virginica", "virginica", "virginica"
		};
		auto stats = CategoryStats::Compute(species);
		
		REQUIRE(stats.uniqueCount == 3);
		// Uniform distribution: entropy = log2(3) ≈ 1.585 bits
		REQUIRE(stats.entropy == Catch::Approx(std::log2(3.0)).epsilon(0.001));
	}

	SECTION("Empty data throws") {
		std::vector<std::string> empty;
		REQUIRE_THROWS_AS(CategoryStats::Compute(empty), StatisticsError);
	}
}

/////////////////////////////////////////////////////////////////////////////////////
///                       DATA DESCRIPTOR TYPE DETECTION                           ///
/////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("DataDescriptor - Type detection", "[datadescriptors][detection]") {
	TEST_PRECISION_INFO();

	SECTION("Detect numeric") {
		std::vector<std::string> samples = {"1.5", "2.7", "3.14", "42", "-10.5"};
		auto type = DataDescriptor::DetectType(samples);
		REQUIRE(type == DataDescriptor::Type::REAL);
	}

	SECTION("Detect boolean") {
		std::vector<std::string> samples = {"true", "false", "True", "FALSE", "yes"};
		auto type = DataDescriptor::DetectType(samples);
		REQUIRE(type == DataDescriptor::Type::BOOL);
	}

	SECTION("Detect date") {
		std::vector<std::string> samples = {"2024-01-15", "2024-02-20", "2023/12/31"};
		auto type = DataDescriptor::DetectType(samples);
		REQUIRE(type == DataDescriptor::Type::DATE);
	}

	SECTION("Detect time") {
		std::vector<std::string> samples = {"10:30:00", "14:45:30", "09:00:00"};
		auto type = DataDescriptor::DetectType(samples);
		REQUIRE(type == DataDescriptor::Type::TIME);
	}

	SECTION("Detect category") {
		std::vector<std::string> samples = {"red", "blue", "green", "yellow"};
		auto type = DataDescriptor::DetectType(samples);
		REQUIRE(type == DataDescriptor::Type::CATEGORY);
	}

	SECTION("Empty returns category") {
		std::vector<std::string> empty;
		auto type = DataDescriptor::DetectType(empty);
		REQUIRE(type == DataDescriptor::Type::CATEGORY);
	}
}

/////////////////////////////////////////////////////////////////////////////////////
///                         PRINT SUMMARY TESTS                                    ///
/////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("DataDescriptors - PrintSummary", "[datadescriptors][output]") {
	TEST_PRECISION_INFO();

	SECTION("BoolDataStats output") {
		std::vector<bool> data = {true, true, false, true, false};
		auto stats = BoolDataStats::Compute(data);
		
		std::ostringstream oss;
		stats.PrintSummary(oss);
		std::string output = oss.str();
		
		REQUIRE(output.find("Boolean Data Statistics") != std::string::npos);
		REQUIRE(output.find("True: 3") != std::string::npos);
		REQUIRE(output.find("False: 2") != std::string::npos);
	}

	SECTION("TimeDataStats output") {
		Vector<Real> durations(3);
		durations[0] = 60; durations[1] = 120; durations[2] = 180;
		auto stats = TimeDataStats::Compute(durations, false);
		
		std::ostringstream oss;
		stats.PrintSummary(oss);
		std::string output = oss.str();
		
		REQUIRE(output.find("Time Data Statistics") != std::string::npos);
		REQUIRE(output.find("Count: 3") != std::string::npos);
	}

	SECTION("GeoDataStats output") {
		Vector<Real> lats(2), lons(2);
		lats[0] = 40.7; lons[0] = -74.0;
		lats[1] = 34.0; lons[1] = -118.2;
		auto stats = GeoDataStats::Compute(lats, lons);
		
		std::ostringstream oss;
		stats.PrintSummary(oss);
		std::string output = oss.str();
		
		REQUIRE(output.find("Geographic Data Statistics") != std::string::npos);
		REQUIRE(output.find("Centroid") != std::string::npos);
	}

	SECTION("CategoryStats output") {
		std::vector<std::string> data = {"A", "B", "A", "A"};
		auto stats = CategoryStats::Compute(data);
		
		std::ostringstream oss;
		stats.PrintSummary(oss);
		std::string output = oss.str();
		
		REQUIRE(output.find("Categorical Data Statistics") != std::string::npos);
		REQUIRE(output.find("Mode: \"A\"") != std::string::npos);
	}
}

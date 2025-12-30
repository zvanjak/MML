///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  Tests for StatisticsTestbed.h - Reference dataset verification                   ///
///  Validates that computed statistics match known reference values                  ///
///////////////////////////////////////////////////////////////////////////////////////////

#include <catch2/catch_all.hpp>
#include <catch2/catch_approx.hpp>

#include "MMLBase.h"
#include "base/Vector.h"
#include "algorithms/Statistics.h"
#include "testbeds/StatisticsTestbed.h"

#include <cmath>
#include <numeric>
#include <algorithm>

using namespace MML;
using namespace MML::Statistics;
using namespace MML::Testbeds::Statistics;

// Test precision helper
#define TEST_PRECISION_INFO() INFO("Test precision: " << std::numeric_limits<Real>::epsilon() << " (double)")

/////////////////////////////////////////////////////////////////////////////////////
///                              IRIS DATASET TESTS                                ///
/////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("StatisticsTestbed - Iris Dataset", "[testbed][iris]") {
	TEST_PRECISION_INFO();

	SECTION("Sepal Length statistics") {
		auto data = IrisDataset::SepalLength();
		
		REQUIRE(data.size() == 150);
		REQUIRE(Mean(data) == Catch::Approx(IrisDataset::SepalLengthMean).epsilon(0.0001));
		REQUIRE(StdDev(data) == Catch::Approx(IrisDataset::SepalLengthStd).epsilon(0.001));
		REQUIRE(Median(data) == Catch::Approx(IrisDataset::SepalLengthMedian));
		
		// Min/Max
		Real minVal = *std::min_element(&data[0], &data[0] + data.size());
		Real maxVal = *std::max_element(&data[0], &data[0] + data.size());
		REQUIRE(minVal == Catch::Approx(IrisDataset::SepalLengthMin));
		REQUIRE(maxVal == Catch::Approx(IrisDataset::SepalLengthMax));
	}

	SECTION("Sepal Width statistics") {
		auto data = IrisDataset::SepalWidth();
		
		REQUIRE(data.size() == 150);
		REQUIRE(Mean(data) == Catch::Approx(IrisDataset::SepalWidthMean).epsilon(0.0001));
		REQUIRE(StdDev(data) == Catch::Approx(IrisDataset::SepalWidthStd).epsilon(0.001));
	}

	SECTION("Petal Length statistics") {
		auto data = IrisDataset::PetalLength();
		
		REQUIRE(data.size() == 150);
		REQUIRE(Mean(data) == Catch::Approx(IrisDataset::PetalLengthMean).epsilon(0.0001));
		REQUIRE(StdDev(data) == Catch::Approx(IrisDataset::PetalLengthStd).epsilon(0.001));
	}

	SECTION("Petal Width statistics") {
		auto data = IrisDataset::PetalWidth();
		
		REQUIRE(data.size() == 150);
		REQUIRE(Mean(data) == Catch::Approx(IrisDataset::PetalWidthMean).epsilon(0.0001));
		REQUIRE(StdDev(data) == Catch::Approx(IrisDataset::PetalWidthStd).epsilon(0.001));
	}

	SECTION("Species distribution") {
		auto species = IrisDataset::Species();
		
		REQUIRE(species.size() == 150);
		
		int setosaCount = std::count(species.begin(), species.end(), "setosa");
		int versicolorCount = std::count(species.begin(), species.end(), "versicolor");
		int virginicaCount = std::count(species.begin(), species.end(), "virginica");
		
		REQUIRE(setosaCount == 50);
		REQUIRE(versicolorCount == 50);
		REQUIRE(virginicaCount == 50);
	}

	SECTION("Per-species Sepal Length means") {
		auto data = IrisDataset::SepalLength();
		
		// Setosa (first 50)
		Vector<Real> setosa(50);
		for (int i = 0; i < 50; ++i) setosa[i] = data[i];
		REQUIRE(Mean(setosa) == Catch::Approx(IrisDataset::SetosaStats::SepalLengthMean).epsilon(0.0001));
		
		// Versicolor (50-99)
		Vector<Real> versicolor(50);
		for (int i = 0; i < 50; ++i) versicolor[i] = data[50 + i];
		REQUIRE(Mean(versicolor) == Catch::Approx(IrisDataset::VersicolorStats::SepalLengthMean).epsilon(0.0001));
		
		// Virginica (100-149)
		Vector<Real> virginica(50);
		for (int i = 0; i < 50; ++i) virginica[i] = data[100 + i];
		REQUIRE(Mean(virginica) == Catch::Approx(IrisDataset::VirginicaStats::SepalLengthMean).epsilon(0.0001));
	}

	SECTION("Correlation between features") {
		auto sepalL = IrisDataset::SepalLength();
		auto petalL = IrisDataset::PetalLength();
		
		// Compute Pearson correlation
		Real corr = PearsonCorrelation(sepalL, petalL);
		
		// Should be strongly positive
		REQUIRE(corr == Catch::Approx(IrisDataset::CorrSepalLengthPetalLength).epsilon(0.001));
	}
}

/////////////////////////////////////////////////////////////////////////////////////
///                          ANSCOMBE'S QUARTET TESTS                              ///
/////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("StatisticsTestbed - Anscombe's Quartet", "[testbed][anscombe]") {
	TEST_PRECISION_INFO();

	SECTION("All datasets have same X mean") {
		auto x1 = AnscombesQuartet::Dataset1X();
		auto x2 = AnscombesQuartet::Dataset2X();
		auto x3 = AnscombesQuartet::Dataset3X();
		auto x4 = AnscombesQuartet::Dataset4X();
		
		REQUIRE(Mean(x1) == Catch::Approx(AnscombesQuartet::MeanX));
		REQUIRE(Mean(x2) == Catch::Approx(AnscombesQuartet::MeanX));
		REQUIRE(Mean(x3) == Catch::Approx(AnscombesQuartet::MeanX));
		REQUIRE(Mean(x4) == Catch::Approx(AnscombesQuartet::MeanX));
	}

	SECTION("All datasets have approximately same Y mean") {
		auto y1 = AnscombesQuartet::Dataset1Y();
		auto y2 = AnscombesQuartet::Dataset2Y();
		auto y3 = AnscombesQuartet::Dataset3Y();
		auto y4 = AnscombesQuartet::Dataset4Y();
		
		Real mean1 = Mean(y1);
		Real mean2 = Mean(y2);
		Real mean3 = Mean(y3);
		Real mean4 = Mean(y4);
		
		// All should be approximately 7.5
		REQUIRE(mean1 == Catch::Approx(AnscombesQuartet::MeanY).epsilon(0.01));
		REQUIRE(mean2 == Catch::Approx(AnscombesQuartet::MeanY).epsilon(0.01));
		REQUIRE(mean3 == Catch::Approx(AnscombesQuartet::MeanY).epsilon(0.01));
		REQUIRE(mean4 == Catch::Approx(AnscombesQuartet::MeanY).epsilon(0.01));
		
		// Verify precise values
		REQUIRE(mean1 == Catch::Approx(AnscombesQuartet::Dataset1MeanY));
	}

	SECTION("All datasets have same X variance") {
		auto x1 = AnscombesQuartet::Dataset1X();
		auto x2 = AnscombesQuartet::Dataset2X();
		auto x3 = AnscombesQuartet::Dataset3X();
		auto x4 = AnscombesQuartet::Dataset4X();
		
		REQUIRE(Variance(x1) == Catch::Approx(AnscombesQuartet::VarianceX));
		REQUIRE(Variance(x2) == Catch::Approx(AnscombesQuartet::VarianceX));
		REQUIRE(Variance(x3) == Catch::Approx(AnscombesQuartet::VarianceX));
		REQUIRE(Variance(x4) == Catch::Approx(AnscombesQuartet::VarianceX));
	}

	SECTION("All datasets have approximately same Y variance") {
		auto y1 = AnscombesQuartet::Dataset1Y();
		auto y2 = AnscombesQuartet::Dataset2Y();
		auto y3 = AnscombesQuartet::Dataset3Y();
		auto y4 = AnscombesQuartet::Dataset4Y();
		
		// All should be approximately 4.127
		REQUIRE(Variance(y1) == Catch::Approx(AnscombesQuartet::VarianceY).epsilon(0.01));
		REQUIRE(Variance(y2) == Catch::Approx(AnscombesQuartet::VarianceY).epsilon(0.01));
		REQUIRE(Variance(y3) == Catch::Approx(AnscombesQuartet::VarianceY).epsilon(0.01));
		REQUIRE(Variance(y4) == Catch::Approx(AnscombesQuartet::VarianceY).epsilon(0.01));
	}

	SECTION("All datasets have approximately same correlation") {
		auto x1 = AnscombesQuartet::Dataset1X();
		auto y1 = AnscombesQuartet::Dataset1Y();
		auto x2 = AnscombesQuartet::Dataset2X();
		auto y2 = AnscombesQuartet::Dataset2Y();
		auto x3 = AnscombesQuartet::Dataset3X();
		auto y3 = AnscombesQuartet::Dataset3Y();
		auto x4 = AnscombesQuartet::Dataset4X();
		auto y4 = AnscombesQuartet::Dataset4Y();
		
		Real corr1 = PearsonCorrelation(x1, y1);
		Real corr2 = PearsonCorrelation(x2, y2);
		Real corr3 = PearsonCorrelation(x3, y3);
		Real corr4 = PearsonCorrelation(x4, y4);
		
		// All should be approximately 0.816
		REQUIRE(corr1 == Catch::Approx(AnscombesQuartet::Correlation).epsilon(0.01));
		REQUIRE(corr2 == Catch::Approx(AnscombesQuartet::Correlation).epsilon(0.01));
		REQUIRE(corr3 == Catch::Approx(AnscombesQuartet::Correlation).epsilon(0.01));
		REQUIRE(corr4 == Catch::Approx(AnscombesQuartet::Correlation).epsilon(0.01));
	}
}

/////////////////////////////////////////////////////////////////////////////////////
///                       STANDARD NORMAL SAMPLES TESTS                            ///
/////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("StatisticsTestbed - Standard Normal Samples", "[testbed][normal]") {
	TEST_PRECISION_INFO();

	SECTION("N100 reproducibility") {
		auto data1 = StandardNormalSamples::N100();
		auto data2 = StandardNormalSamples::N100();
		
		// Same seed should produce same values
		REQUIRE(data1.size() == 100);
		REQUIRE(data2.size() == 100);
		for (int i = 0; i < 100; ++i) {
			REQUIRE(data1[i] == data2[i]);
		}
	}

	SECTION("N100 statistics within reasonable bounds") {
		auto data = StandardNormalSamples::N100();
		
		// With 100 samples from N(0,1), mean should be within ~0.2 of 0
		// and std within ~0.2 of 1 (3 sigma rule for sample statistics)
		REQUIRE(std::abs(Mean(data)) < 0.3);
		REQUIRE(StdDev(data) > 0.7);
		REQUIRE(StdDev(data) < 1.3);
	}

	SECTION("CLT convergence - mean approaches 0") {
		auto n100 = StandardNormalSamples::N100();
		auto n1000 = StandardNormalSamples::N1000();
		auto n10000 = StandardNormalSamples::N10000();
		
		Real mean100 = std::abs(Mean(n100));
		Real mean1000 = std::abs(Mean(n1000));
		Real mean10000 = std::abs(Mean(n10000));
		
		// Larger samples should have means closer to 0 (generally)
		// Allow some statistical fluctuation
		INFO("mean100 = " << mean100 << ", mean1000 = " << mean1000 << ", mean10000 = " << mean10000);
		
		// All should be reasonably small
		REQUIRE(mean10000 < 0.1);  // Very close to 0 with 10k samples
	}

	SECTION("CLT convergence - std approaches 1") {
		auto n1000 = StandardNormalSamples::N1000();
		auto n10000 = StandardNormalSamples::N10000();
		
		Real std1000 = StdDev(n1000);
		Real std10000 = StdDev(n10000);
		
		// Both should be close to 1.0 (allow statistical fluctuation)
		REQUIRE(std1000 == Catch::Approx(1.0).epsilon(0.1));
		REQUIRE(std10000 == Catch::Approx(1.0).epsilon(0.05));
	}
}

/////////////////////////////////////////////////////////////////////////////////////
///                        TIME SERIES EXAMPLES TESTS                              ///
/////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("StatisticsTestbed - Time Series Examples", "[testbed][timeseries]") {
	TEST_PRECISION_INFO();

	SECTION("Random walk reproducibility") {
		auto rw1 = TimeSeriesExamples::RandomWalk(100, 42);
		auto rw2 = TimeSeriesExamples::RandomWalk(100, 42);
		
		REQUIRE(rw1.size() == 100);
		for (int i = 0; i < 100; ++i) {
			REQUIRE(rw1[i] == rw2[i]);
		}
	}

	SECTION("Random walk starts at 0") {
		auto rw = TimeSeriesExamples::RandomWalk(100, 42);
		REQUIRE(rw[0] == Catch::Approx(0.0));
	}

	SECTION("Seasonal pattern has correct period") {
		auto seasonal = TimeSeriesExamples::Seasonal12Month(3, 10.0, 0.0, 42);  // No noise
		
		REQUIRE(seasonal.size() == 36);  // 3 years * 12 months
		
		// Check periodicity: value at month 0 and month 12 should be similar
		// sin(0) = 0, sin(2*pi) = 0 (with floating point tolerance)
		REQUIRE(seasonal[0] == Catch::Approx(seasonal[12]).margin(1e-10));
		REQUIRE(seasonal[0] == Catch::Approx(seasonal[24]).margin(1e-10));
		
		// Maximum at month 3 (sin(pi/2) = 1)
		// Month 3: sin(2*pi*3/12) = sin(pi/2) = 1
		Real expected_max = 10.0;  // amplitude * sin(pi/2)
		REQUIRE(seasonal[3] == Catch::Approx(expected_max).epsilon(0.01));
	}

	SECTION("Trend with noise has correct slope") {
		auto trend = TimeSeriesExamples::TrendWithNoise(100, 0.5, 10.0, 0.0, 42);  // No noise
		
		// y = 0.5 * t + 10
		REQUIRE(trend[0] == Catch::Approx(10.0));
		REQUIRE(trend[10] == Catch::Approx(15.0));  // 0.5 * 10 + 10
		REQUIRE(trend[99] == Catch::Approx(59.5));  // 0.5 * 99 + 10
	}

	SECTION("AR(1) process is stationary") {
		auto ar1 = TimeSeriesExamples::AR1Process(1000, 0.7, 42);
		
		// Mean should be close to 0 for stationary AR(1)
		REQUIRE(Mean(ar1) == Catch::Approx(0.0).margin(0.2));
		
		// Variance should be approximately 1 (by construction)
		REQUIRE(Variance(ar1) == Catch::Approx(1.0).epsilon(0.2));
	}

	SECTION("AR(1) theoretical ACF") {
		// phi^k for k = 0, 1, 2, 3
		REQUIRE(TimeSeriesExamples::AR1TheoreticalACF(0, 0.7) == Catch::Approx(1.0));
		REQUIRE(TimeSeriesExamples::AR1TheoreticalACF(1, 0.7) == Catch::Approx(0.7));
		REQUIRE(TimeSeriesExamples::AR1TheoreticalACF(2, 0.7) == Catch::Approx(0.49));
		REQUIRE(TimeSeriesExamples::AR1TheoreticalACF(3, 0.7) == Catch::Approx(0.343));
	}
}

/////////////////////////////////////////////////////////////////////////////////////
///                       NUMERICAL EDGE CASES TESTS                               ///
/////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("StatisticsTestbed - Numerical Edge Cases", "[testbed][edgecases]") {
	TEST_PRECISION_INFO();

	SECTION("Large values - mean computation") {
		auto data = NumericalEdgeCases::LargeValues();
		
		Real mean = Mean(data);
		REQUIRE(mean == Catch::Approx(NumericalEdgeCases::LargeValues_Mean).epsilon(0.01));
	}

	SECTION("Large values - variance computation") {
		auto data = NumericalEdgeCases::LargeValues();
		
		// Should compute without overflow
		Real std = StdDev(data);
		
		// Verify approximate magnitude
		REQUIRE(std > 1e12);  // Should be in the trillions
		REQUIRE(std < 1e13);
	}

	SECTION("Small values - mean computation") {
		auto data = NumericalEdgeCases::SmallValues();
		
		Real mean = Mean(data);
		REQUIRE(mean == Catch::Approx(NumericalEdgeCases::SmallValues_Mean).epsilon(0.01));
	}

	SECTION("Small values - variance computation") {
		auto data = NumericalEdgeCases::SmallValues();
		
		// Should compute without underflow
		Real var = Variance(data);
		REQUIRE(var > 0);  // Should be positive
		REQUIRE(var < 1e-29);  // Should be very small
	}

	SECTION("Nearly constant data") {
		auto data = NumericalEdgeCases::NearlyConstant();
		
		Real mean = Mean(data);
		Real var = Variance(data);
		
		// Mean should be close to 1000
		REQUIRE(mean == Catch::Approx(NumericalEdgeCases::NearlyConstant_Mean).epsilon(1e-8));
		
		// Variance should be very small (near zero)
		REQUIRE(var < 1e-15);
	}

	SECTION("With outliers - mean vs median") {
		auto data = NumericalEdgeCases::WithOutliers();
		
		Real mean = Mean(data);
		Real median = Median(data);
		
		// Median should be robust (around 10)
		REQUIRE(median == Catch::Approx(NumericalEdgeCases::WithOutliers_MedianApprox).epsilon(0.5));
		
		// Mean should be pulled up by outliers
		REQUIRE(mean > median);
		REQUIRE(mean == Catch::Approx(NumericalEdgeCases::WithOutliers_MeanApprox).epsilon(0.5));
	}

	SECTION("Uniform distribution theoretical values") {
		auto data = NumericalEdgeCases::Uniform01(10000, 42);
		
		Real mean = Mean(data);
		Real var = Variance(data);
		
		// Should be close to theoretical values with large n
		REQUIRE(mean == Catch::Approx(NumericalEdgeCases::Uniform01_TheoreticalMean).epsilon(0.02));
		REQUIRE(var == Catch::Approx(NumericalEdgeCases::Uniform01_TheoreticalVariance).epsilon(0.01));
	}
}

/////////////////////////////////////////////////////////////////////////////////////
///                          MTCARS DATASET TESTS                                  ///
/////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("StatisticsTestbed - MTCars Dataset", "[testbed][mtcars]") {
	TEST_PRECISION_INFO();

	SECTION("MPG statistics") {
		auto mpg = MtcarsSubset::MPG();
		
		REQUIRE(mpg.size() == 32);
		REQUIRE(Mean(mpg) == Catch::Approx(MtcarsSubset::MPG_Mean).epsilon(0.0001));
		REQUIRE(StdDev(mpg) == Catch::Approx(MtcarsSubset::MPG_Std).epsilon(0.001));
		REQUIRE(Median(mpg) == Catch::Approx(MtcarsSubset::MPG_Median));
	}

	SECTION("HP statistics") {
		auto hp = MtcarsSubset::HP();
		
		REQUIRE(hp.size() == 32);
		REQUIRE(Mean(hp) == Catch::Approx(MtcarsSubset::HP_Mean).epsilon(0.0001));
		REQUIRE(StdDev(hp) == Catch::Approx(MtcarsSubset::HP_Std).epsilon(0.001));
	}

	SECTION("Weight statistics") {
		auto wt = MtcarsSubset::Weight();
		
		REQUIRE(wt.size() == 32);
		REQUIRE(Mean(wt) == Catch::Approx(MtcarsSubset::Weight_Mean).epsilon(0.0001));
		REQUIRE(StdDev(wt) == Catch::Approx(MtcarsSubset::Weight_Std).epsilon(0.001));
	}

	SECTION("Negative correlation: MPG vs HP") {
		auto mpg = MtcarsSubset::MPG();
		auto hp = MtcarsSubset::HP();
		
		Real corr = PearsonCorrelation(mpg, hp);
		
		// Should be negative - more HP = less fuel efficiency
		REQUIRE(corr < 0);
		REQUIRE(corr == Catch::Approx(MtcarsSubset::Corr_MPG_HP).epsilon(0.001));
	}

	SECTION("Negative correlation: MPG vs Weight") {
		auto mpg = MtcarsSubset::MPG();
		auto wt = MtcarsSubset::Weight();
		
		Real corr = PearsonCorrelation(mpg, wt);
		
		// Should be negative - heavier cars = worse fuel efficiency
		REQUIRE(corr < 0);
		REQUIRE(corr == Catch::Approx(MtcarsSubset::Corr_MPG_Weight).epsilon(0.001));
	}
}

/////////////////////////////////////////////////////////////////////////////////////
///                      COMBINED USAGE PATTERNS                                   ///
/////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("StatisticsTestbed - Usage patterns", "[testbed][usage]") {
	TEST_PRECISION_INFO();

	SECTION("Typical test pattern - verify mean") {
		// Example: How to use testbed in tests
		auto data = IrisDataset::SepalLength();
		
		REQUIRE(Mean(data) == Catch::Approx(IrisDataset::SepalLengthMean).epsilon(1e-5));
	}

	SECTION("Generate reproducible random data") {
		// For tests that need random data but must be reproducible
		auto normal = StandardNormalSamples::Generate(50, 123);  // Custom seed
		auto normal2 = StandardNormalSamples::Generate(50, 123);
		
		// Same seed = same data
		for (int i = 0; i < 50; ++i) {
			REQUIRE(normal[i] == normal2[i]);
		}
	}

	SECTION("Edge case validation pattern") {
		// Test algorithm robustness
		auto large = NumericalEdgeCases::LargeValues();
		auto small = NumericalEdgeCases::SmallValues();
		auto outliers = NumericalEdgeCases::WithOutliers();
		
		// All should compute without exceptions
		REQUIRE_NOTHROW(Mean(large));
		REQUIRE_NOTHROW(Variance(small));
		REQUIRE_NOTHROW(StdDev(outliers));
	}
}

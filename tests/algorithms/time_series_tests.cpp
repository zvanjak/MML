#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "algorithms/Statistics/TimeSeries.h"
#endif

using namespace MML;
using namespace MML::Testing;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace MML::Tests::Algorithms::TimeSeriesTests
{
    /////////////////////////////////////////////////////////////////////////////////////
    ///                        SIMPLE MOVING AVERAGE TESTS                            ///
    /////////////////////////////////////////////////////////////////////////////////////
    
    TEST_CASE("TimeSeries::SimpleMovingAverage_BasicCalculation", "[timeseries][sma]")
    {
        TEST_PRECISION_INFO();
        // SMA of [1, 2, 3, 4, 5] with window=3 should be [2, 3, 4]
        Vector<Real> data({1.0, 2.0, 3.0, 4.0, 5.0});
        
        auto sma = Statistics::TimeSeries::SimpleMovingAverage(data, 3);
        
        REQUIRE(sma.size() == 3);
        REQUIRE_THAT(sma[0], WithinAbs(2.0, 1e-10));  // (1+2+3)/3
        REQUIRE_THAT(sma[1], WithinAbs(3.0, 1e-10));  // (2+3+4)/3
        REQUIRE_THAT(sma[2], WithinAbs(4.0, 1e-10));  // (3+4+5)/3
    }
    
    TEST_CASE("TimeSeries::SimpleMovingAverage_WindowSize1", "[timeseries][sma]")
    {
        TEST_PRECISION_INFO();
        // SMA with window=1 should return the original data
        Vector<Real> data({1.5, 2.5, 3.5, 4.5});
        
        auto sma = Statistics::TimeSeries::SimpleMovingAverage(data, 1);
        
        REQUIRE(sma.size() == 4);
        for (int i = 0; i < 4; ++i)
            REQUIRE_THAT(sma[i], WithinAbs(data[i], 1e-10));
    }
    
    TEST_CASE("TimeSeries::SimpleMovingAverage_WindowEqualsDataLength", "[timeseries][sma]")
    {
        TEST_PRECISION_INFO();
        // When window equals data length, result is single average
        Vector<Real> data({2.0, 4.0, 6.0, 8.0});
        
        auto sma = Statistics::TimeSeries::SimpleMovingAverage(data, 4);
        
        REQUIRE(sma.size() == 1);
        REQUIRE_THAT(sma[0], WithinAbs(5.0, 1e-10));  // (2+4+6+8)/4
    }
    
    TEST_CASE("TimeSeries::SimpleMovingAverage_ConstantData", "[timeseries][sma]")
    {
        TEST_PRECISION_INFO();
        // SMA of constant data should be that constant
        Vector<Real> data({5.0, 5.0, 5.0, 5.0, 5.0});
        
        auto sma = Statistics::TimeSeries::SimpleMovingAverage(data, 3);
        
        REQUIRE(sma.size() == 3);
        for (int i = 0; i < 3; ++i)
            REQUIRE_THAT(sma[i], WithinAbs(5.0, 1e-10));
    }
    
    TEST_CASE("TimeSeries::SimpleMovingAverage_EmptyData_Throws", "[timeseries][sma][error]")
    {
        TEST_PRECISION_INFO();
        Vector<Real> data;
        
        REQUIRE_THROWS_AS(Statistics::TimeSeries::SimpleMovingAverage(data, 3), StatisticsError);
    }
    
    TEST_CASE("TimeSeries::SimpleMovingAverage_WindowTooLarge_Throws", "[timeseries][sma][error]")
    {
        TEST_PRECISION_INFO();
        Vector<Real> data({1.0, 2.0, 3.0});
        
        REQUIRE_THROWS_AS(Statistics::TimeSeries::SimpleMovingAverage(data, 5), StatisticsError);
    }
    
    TEST_CASE("TimeSeries::SimpleMovingAverage_InvalidWindowSize_Throws", "[timeseries][sma][error]")
    {
        TEST_PRECISION_INFO();
        Vector<Real> data({1.0, 2.0, 3.0});
        
        REQUIRE_THROWS_AS(Statistics::TimeSeries::SimpleMovingAverage(data, 0), StatisticsError);
        REQUIRE_THROWS_AS(Statistics::TimeSeries::SimpleMovingAverage(data, -1), StatisticsError);
    }
    
    
    /////////////////////////////////////////////////////////////////////////////////////
    ///                     EXPONENTIAL MOVING AVERAGE TESTS                          ///
    /////////////////////////////////////////////////////////////////////////////////////
    
    TEST_CASE("TimeSeries::ExponentialMovingAverage_BasicCalculation", "[timeseries][ema]")
    {
        TEST_PRECISION_INFO();
        // EMA with alpha=0.5: EMA_t = 0.5*x_t + 0.5*EMA_{t-1}
        Vector<Real> data({1.0, 2.0, 3.0, 4.0, 5.0});
        
        auto ema = Statistics::TimeSeries::ExponentialMovingAverage(data, 0.5);
        
        REQUIRE(ema.size() == 5);
        REQUIRE_THAT(ema[0], WithinAbs(1.0, 1e-10));      // Initial = x_0
        REQUIRE_THAT(ema[1], WithinAbs(1.5, 1e-10));      // 0.5*2 + 0.5*1 = 1.5
        REQUIRE_THAT(ema[2], WithinAbs(2.25, 1e-10));     // 0.5*3 + 0.5*1.5 = 2.25
        REQUIRE_THAT(ema[3], WithinAbs(3.125, 1e-10));    // 0.5*4 + 0.5*2.25 = 3.125
        REQUIRE_THAT(ema[4], WithinAbs(4.0625, 1e-10));   // 0.5*5 + 0.5*3.125 = 4.0625
    }
    
    TEST_CASE("TimeSeries::ExponentialMovingAverage_HighAlpha", "[timeseries][ema]")
    {
        TEST_PRECISION_INFO();
        // High alpha (0.9) = more weight on recent values
        Vector<Real> data({10.0, 20.0, 30.0});
        
        auto ema = Statistics::TimeSeries::ExponentialMovingAverage(data, 0.9);
        
        REQUIRE(ema.size() == 3);
        REQUIRE_THAT(ema[0], WithinAbs(10.0, 1e-10));
        REQUIRE_THAT(ema[1], WithinAbs(19.0, 1e-10));   // 0.9*20 + 0.1*10 = 19
        REQUIRE_THAT(ema[2], WithinAbs(28.9, 1e-10));   // 0.9*30 + 0.1*19 = 28.9
    }
    
    TEST_CASE("TimeSeries::ExponentialMovingAverage_LowAlpha", "[timeseries][ema]")
    {
        TEST_PRECISION_INFO();
        // Low alpha (0.1) = more weight on history
        Vector<Real> data({10.0, 20.0, 30.0});
        
        auto ema = Statistics::TimeSeries::ExponentialMovingAverage(data, 0.1);
        
        REQUIRE(ema.size() == 3);
        REQUIRE_THAT(ema[0], WithinAbs(10.0, 1e-10));
        REQUIRE_THAT(ema[1], WithinAbs(11.0, 1e-10));   // 0.1*20 + 0.9*10 = 11
        REQUIRE_THAT(ema[2], WithinAbs(12.9, 1e-10));   // 0.1*30 + 0.9*11 = 12.9
    }
    
    TEST_CASE("TimeSeries::ExponentialMovingAverage_ConstantData", "[timeseries][ema]")
    {
        TEST_PRECISION_INFO();
        // EMA of constant data should be that constant
        Vector<Real> data({7.0, 7.0, 7.0, 7.0});
        
        auto ema = Statistics::TimeSeries::ExponentialMovingAverage(data, 0.5);
        
        for (int i = 0; i < 4; ++i)
            REQUIRE_THAT(ema[i], WithinAbs(7.0, 1e-10));
    }
    
    TEST_CASE("TimeSeries::ExponentialMovingAverage_InvalidAlpha_Throws", "[timeseries][ema][error]")
    {
        TEST_PRECISION_INFO();
        Vector<Real> data({1.0, 2.0, 3.0});
        
        REQUIRE_THROWS_AS(Statistics::TimeSeries::ExponentialMovingAverage(data, 0.0), StatisticsError);
        REQUIRE_THROWS_AS(Statistics::TimeSeries::ExponentialMovingAverage(data, 1.0), StatisticsError);
        REQUIRE_THROWS_AS(Statistics::TimeSeries::ExponentialMovingAverage(data, -0.5), StatisticsError);
        REQUIRE_THROWS_AS(Statistics::TimeSeries::ExponentialMovingAverage(data, 1.5), StatisticsError);
    }
    
    TEST_CASE("TimeSeries::ExponentialMovingAverageSpan_CalculatesAlphaCorrectly", "[timeseries][ema]")
    {
        TEST_PRECISION_INFO();
        // span=9 gives alpha = 2/(9+1) = 0.2
        Vector<Real> data({10.0, 20.0, 30.0});
        
        auto emaSpan = Statistics::TimeSeries::ExponentialMovingAverageSpan(data, 9);
        auto emaAlpha = Statistics::TimeSeries::ExponentialMovingAverage(data, 0.2);
        
        REQUIRE(emaSpan.size() == emaAlpha.size());
        for (int i = 0; i < emaSpan.size(); ++i)
            REQUIRE_THAT(emaSpan[i], WithinAbs(emaAlpha[i], 1e-10));
    }
    
    
    /////////////////////////////////////////////////////////////////////////////////////
    ///                      WEIGHTED MOVING AVERAGE TESTS                            ///
    /////////////////////////////////////////////////////////////////////////////////////
    
    TEST_CASE("TimeSeries::WeightedMovingAverage_CustomWeights", "[timeseries][wma]")
    {
        TEST_PRECISION_INFO();
        // WMA with weights [1, 2, 3] (most recent gets weight 3)
        Vector<Real> data({10.0, 20.0, 30.0, 40.0});
        Vector<Real> weights({1.0, 2.0, 3.0});
        
        auto wma = Statistics::TimeSeries::WeightedMovingAverage(data, weights);
        
        REQUIRE(wma.size() == 2);
        // First window [10, 20, 30]: (1*10 + 2*20 + 3*30) / 6 = 140/6 = 23.333...
        REQUIRE_THAT(wma[0], WithinAbs(140.0/6.0, 1e-10));
        // Second window [20, 30, 40]: (1*20 + 2*30 + 3*40) / 6 = 200/6 = 33.333...
        REQUIRE_THAT(wma[1], WithinAbs(200.0/6.0, 1e-10));
    }
    
    TEST_CASE("TimeSeries::WeightedMovingAverage_EqualWeights", "[timeseries][wma]")
    {
        TEST_PRECISION_INFO();
        // Equal weights should give same result as SMA
        Vector<Real> data({1.0, 2.0, 3.0, 4.0, 5.0});
        Vector<Real> weights({1.0, 1.0, 1.0});
        
        auto wma = Statistics::TimeSeries::WeightedMovingAverage(data, weights);
        auto sma = Statistics::TimeSeries::SimpleMovingAverage(data, 3);
        
        REQUIRE(wma.size() == sma.size());
        for (int i = 0; i < wma.size(); ++i)
            REQUIRE_THAT(wma[i], WithinAbs(sma[i], 1e-10));
    }
    
    TEST_CASE("TimeSeries::WeightedMovingAverageLinear_BasicCalculation", "[timeseries][wma]")
    {
        TEST_PRECISION_INFO();
        // Linear weights for window=3: [1, 2, 3]
        Vector<Real> data({10.0, 20.0, 30.0, 40.0});
        
        auto wmaLinear = Statistics::TimeSeries::WeightedMovingAverageLinear(data, 3);
        Vector<Real> linearWeights({1.0, 2.0, 3.0});
        auto wmaCustom = Statistics::TimeSeries::WeightedMovingAverage(data, linearWeights);
        
        REQUIRE(wmaLinear.size() == wmaCustom.size());
        for (int i = 0; i < wmaLinear.size(); ++i)
            REQUIRE_THAT(wmaLinear[i], WithinAbs(wmaCustom[i], 1e-10));
    }
    
    TEST_CASE("TimeSeries::WeightedMovingAverage_ZeroWeights_Throws", "[timeseries][wma][error]")
    {
        TEST_PRECISION_INFO();
        Vector<Real> data({1.0, 2.0, 3.0});
        Vector<Real> zeroWeights({0.0, 0.0, 0.0});
        
        REQUIRE_THROWS_AS(Statistics::TimeSeries::WeightedMovingAverage(data, zeroWeights), StatisticsError);
    }
    
    
    /////////////////////////////////////////////////////////////////////////////////////
    ///                        ROLLING STATISTICS TESTS                               ///
    /////////////////////////////////////////////////////////////////////////////////////
    
    TEST_CASE("TimeSeries::RollingSum_BasicCalculation", "[timeseries][rolling]")
    {
        TEST_PRECISION_INFO();
        Vector<Real> data({1.0, 2.0, 3.0, 4.0, 5.0});
        
        auto rsum = Statistics::TimeSeries::RollingSum(data, 3);
        
        REQUIRE(rsum.size() == 3);
        REQUIRE_THAT(rsum[0], WithinAbs(6.0, 1e-10));   // 1+2+3
        REQUIRE_THAT(rsum[1], WithinAbs(9.0, 1e-10));   // 2+3+4
        REQUIRE_THAT(rsum[2], WithinAbs(12.0, 1e-10)); // 3+4+5
    }
    
    TEST_CASE("TimeSeries::RollingVariance_BasicCalculation", "[timeseries][rolling]")
    {
        TEST_PRECISION_INFO();
        // Variance of [1, 2, 3] = ((1-2)² + (2-2)² + (3-2)²) / 2 = 2/2 = 1
        Vector<Real> data({1.0, 2.0, 3.0, 4.0, 5.0});
        
        auto rvar = Statistics::TimeSeries::RollingVariance(data, 3);
        
        REQUIRE(rvar.size() == 3);
        REQUIRE_THAT(rvar[0], WithinAbs(1.0, 1e-10));   // Var([1,2,3])
        REQUIRE_THAT(rvar[1], WithinAbs(1.0, 1e-10));   // Var([2,3,4])
        REQUIRE_THAT(rvar[2], WithinAbs(1.0, 1e-10));   // Var([3,4,5])
    }
    
    TEST_CASE("TimeSeries::RollingVariance_ConstantData", "[timeseries][rolling]")
    {
        TEST_PRECISION_INFO();
        // Variance of constant data should be 0
        Vector<Real> data({5.0, 5.0, 5.0, 5.0});
        
        auto rvar = Statistics::TimeSeries::RollingVariance(data, 3);
        
        REQUIRE(rvar.size() == 2);
        REQUIRE_THAT(rvar[0], WithinAbs(0.0, 1e-10));
        REQUIRE_THAT(rvar[1], WithinAbs(0.0, 1e-10));
    }
    
    TEST_CASE("TimeSeries::RollingStdDev_BasicCalculation", "[timeseries][rolling]")
    {
        TEST_PRECISION_INFO();
        Vector<Real> data({1.0, 2.0, 3.0, 4.0, 5.0});
        
        auto rstd = Statistics::TimeSeries::RollingStdDev(data, 3);
        auto rvar = Statistics::TimeSeries::RollingVariance(data, 3);
        
        REQUIRE(rstd.size() == rvar.size());
        for (int i = 0; i < rstd.size(); ++i)
            REQUIRE_THAT(rstd[i], WithinAbs(std::sqrt(rvar[i]), 1e-10));
    }
    
    TEST_CASE("TimeSeries::RollingMin_BasicCalculation", "[timeseries][rolling]")
    {
        TEST_PRECISION_INFO();
        Vector<Real> data({3.0, 1.0, 4.0, 1.0, 5.0, 9.0, 2.0});
        
        auto rmin = Statistics::TimeSeries::RollingMin(data, 3);
        
        REQUIRE(rmin.size() == 5);
        REQUIRE_THAT(rmin[0], WithinAbs(1.0, 1e-10));  // min(3,1,4) = 1
        REQUIRE_THAT(rmin[1], WithinAbs(1.0, 1e-10));  // min(1,4,1) = 1
        REQUIRE_THAT(rmin[2], WithinAbs(1.0, 1e-10));  // min(4,1,5) = 1
        REQUIRE_THAT(rmin[3], WithinAbs(1.0, 1e-10));  // min(1,5,9) = 1
        REQUIRE_THAT(rmin[4], WithinAbs(2.0, 1e-10));  // min(5,9,2) = 2
    }
    
    TEST_CASE("TimeSeries::RollingMax_BasicCalculation", "[timeseries][rolling]")
    {
        TEST_PRECISION_INFO();
        Vector<Real> data({3.0, 1.0, 4.0, 1.0, 5.0, 9.0, 2.0});
        
        auto rmax = Statistics::TimeSeries::RollingMax(data, 3);
        
        REQUIRE(rmax.size() == 5);
        REQUIRE_THAT(rmax[0], WithinAbs(4.0, 1e-10));  // max(3,1,4) = 4
        REQUIRE_THAT(rmax[1], WithinAbs(4.0, 1e-10));  // max(1,4,1) = 4
        REQUIRE_THAT(rmax[2], WithinAbs(5.0, 1e-10));  // max(4,1,5) = 5
        REQUIRE_THAT(rmax[3], WithinAbs(9.0, 1e-10));  // max(1,5,9) = 9
        REQUIRE_THAT(rmax[4], WithinAbs(9.0, 1e-10));  // max(5,9,2) = 9
    }
    
    TEST_CASE("TimeSeries::RollingMinMax_MonotonicIncreasing", "[timeseries][rolling]")
    {
        TEST_PRECISION_INFO();
        // For monotonic increasing data, min is first element, max is last
        Vector<Real> data({1.0, 2.0, 3.0, 4.0, 5.0});
        
        auto rmin = Statistics::TimeSeries::RollingMin(data, 3);
        auto rmax = Statistics::TimeSeries::RollingMax(data, 3);
        
        // min should always be the leftmost element of window
        REQUIRE_THAT(rmin[0], WithinAbs(1.0, 1e-10));
        REQUIRE_THAT(rmin[1], WithinAbs(2.0, 1e-10));
        REQUIRE_THAT(rmin[2], WithinAbs(3.0, 1e-10));
        
        // max should always be the rightmost element of window
        REQUIRE_THAT(rmax[0], WithinAbs(3.0, 1e-10));
        REQUIRE_THAT(rmax[1], WithinAbs(4.0, 1e-10));
        REQUIRE_THAT(rmax[2], WithinAbs(5.0, 1e-10));
    }
    
    
    /////////////////////////////////////////////////////////////////////////////////////
    ///                         AUTOCORRELATION TESTS                                 ///
    /////////////////////////////////////////////////////////////////////////////////////
    
    TEST_CASE("TimeSeries::Autocorrelation_Lag0_IsOne", "[timeseries][acf]")
    {
        TEST_PRECISION_INFO();
        Vector<Real> data({1.0, 2.0, 3.0, 4.0, 5.0});
        
        Real acf0 = Statistics::TimeSeries::Autocorrelation(data, 0);
        
        REQUIRE_THAT(acf0, WithinAbs(1.0, 1e-10));
    }
    
    TEST_CASE("TimeSeries::Autocorrelation_ConstantData", "[timeseries][acf]")
    {
        TEST_PRECISION_INFO();
        // Constant data has no autocorrelation (variance = 0)
        Vector<Real> data({5.0, 5.0, 5.0, 5.0, 5.0});
        
        Real acf1 = Statistics::TimeSeries::Autocorrelation(data, 1);
        Real acf2 = Statistics::TimeSeries::Autocorrelation(data, 2);
        
        REQUIRE_THAT(acf1, WithinAbs(0.0, 1e-10));
        REQUIRE_THAT(acf2, WithinAbs(0.0, 1e-10));
    }
    
    TEST_CASE("TimeSeries::Autocorrelation_LinearTrend", "[timeseries][acf]")
    {
        TEST_PRECISION_INFO();
        // Linear trend has high autocorrelation
        Vector<Real> data({1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0});
        
        Real acf1 = Statistics::TimeSeries::Autocorrelation(data, 1);
        
        // Should be close to 1 for highly correlated trend
        // The exact value for this series is 0.7 (n=10, lag=1)
        REQUIRE(acf1 >= 0.7);
    }
    
    TEST_CASE("TimeSeries::AutocorrelationFunction_BasicCalculation", "[timeseries][acf]")
    {
        TEST_PRECISION_INFO();
        Vector<Real> data({1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0});
        
        auto acf = Statistics::TimeSeries::AutocorrelationFunction(data, 3);
        
        REQUIRE(acf.size() == 4);
        REQUIRE_THAT(acf[0], WithinAbs(1.0, 1e-10));  // ACF(0) always 1
        // Alternating pattern should have negative ACF at lag 1, positive at lag 2
        REQUIRE(acf[1] < 0);  // Negative correlation at lag 1
        REQUIRE(acf[2] > 0);  // Positive correlation at lag 2
    }
    
    TEST_CASE("TimeSeries::PartialAutocorrelation_Lag0_IsOne", "[timeseries][pacf]")
    {
        TEST_PRECISION_INFO();
        Vector<Real> data({1.0, 2.0, 3.0, 4.0, 5.0});
        
        Real pacf0 = Statistics::TimeSeries::PartialAutocorrelation(data, 0);
        
        REQUIRE_THAT(pacf0, WithinAbs(1.0, 1e-10));
    }
    
    TEST_CASE("TimeSeries::PartialAutocorrelationFunction_BasicCalculation", "[timeseries][pacf]")
    {
        TEST_PRECISION_INFO();
        Vector<Real> data({1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0});
        
        auto pacf = Statistics::TimeSeries::PartialAutocorrelationFunction(data, 3);
        
        REQUIRE(pacf.size() == 4);
        REQUIRE_THAT(pacf[0], WithinAbs(1.0, 1e-10));  // PACF(0) always 1
    }
    
    TEST_CASE("TimeSeries::Autocorrelation_InvalidLag_Throws", "[timeseries][acf][error]")
    {
        TEST_PRECISION_INFO();
        Vector<Real> data({1.0, 2.0, 3.0});
        
        REQUIRE_THROWS_AS(Statistics::TimeSeries::Autocorrelation(data, -1), StatisticsError);
        REQUIRE_THROWS_AS(Statistics::TimeSeries::Autocorrelation(data, 3), StatisticsError);
        REQUIRE_THROWS_AS(Statistics::TimeSeries::Autocorrelation(data, 10), StatisticsError);
    }
    
    
    /////////////////////////////////////////////////////////////////////////////////////
    ///                         DIFFERENCING TESTS                                    ///
    /////////////////////////////////////////////////////////////////////////////////////
    
    TEST_CASE("TimeSeries::FirstDifference_BasicCalculation", "[timeseries][diff]")
    {
        TEST_PRECISION_INFO();
        Vector<Real> data({1.0, 3.0, 6.0, 10.0, 15.0});
        
        auto diff = Statistics::TimeSeries::FirstDifference(data);
        
        REQUIRE(diff.size() == 4);
        REQUIRE_THAT(diff[0], WithinAbs(2.0, 1e-10));   // 3 - 1
        REQUIRE_THAT(diff[1], WithinAbs(3.0, 1e-10));   // 6 - 3
        REQUIRE_THAT(diff[2], WithinAbs(4.0, 1e-10));   // 10 - 6
        REQUIRE_THAT(diff[3], WithinAbs(5.0, 1e-10));   // 15 - 10
    }
    
    TEST_CASE("TimeSeries::FirstDifference_LinearTrendBecomesConstant", "[timeseries][diff]")
    {
        TEST_PRECISION_INFO();
        // Linear trend y = 2x: [0, 2, 4, 6, 8]
        // First difference should be constant [2, 2, 2, 2]
        Vector<Real> data({0.0, 2.0, 4.0, 6.0, 8.0});
        
        auto diff = Statistics::TimeSeries::FirstDifference(data);
        
        for (int i = 0; i < diff.size(); ++i)
            REQUIRE_THAT(diff[i], WithinAbs(2.0, 1e-10));
    }
    
    TEST_CASE("TimeSeries::Difference_Order2", "[timeseries][diff]")
    {
        TEST_PRECISION_INFO();
        // Quadratic: [1, 4, 9, 16, 25] (squares)
        // First diff: [3, 5, 7, 9]
        // Second diff: [2, 2, 2]
        Vector<Real> data({1.0, 4.0, 9.0, 16.0, 25.0});
        
        auto diff2 = Statistics::TimeSeries::Difference(data, 2);
        
        REQUIRE(diff2.size() == 3);
        for (int i = 0; i < diff2.size(); ++i)
            REQUIRE_THAT(diff2[i], WithinAbs(2.0, 1e-10));
    }
    
    TEST_CASE("TimeSeries::SeasonalDifference_BasicCalculation", "[timeseries][diff]")
    {
        TEST_PRECISION_INFO();
        // Quarterly data with trend: [10, 20, 30, 40, 15, 25, 35, 45]
        // Seasonal diff (period=4): [15-10, 25-20, 35-30, 45-40] = [5, 5, 5, 5]
        Vector<Real> data({10.0, 20.0, 30.0, 40.0, 15.0, 25.0, 35.0, 45.0});
        
        auto sdiff = Statistics::TimeSeries::SeasonalDifference(data, 4);
        
        REQUIRE(sdiff.size() == 4);
        for (int i = 0; i < sdiff.size(); ++i)
            REQUIRE_THAT(sdiff[i], WithinAbs(5.0, 1e-10));
    }
    
    TEST_CASE("TimeSeries::FirstDifference_TooShort_Throws", "[timeseries][diff][error]")
    {
        TEST_PRECISION_INFO();
        Vector<Real> data({1.0});
        
        REQUIRE_THROWS_AS(Statistics::TimeSeries::FirstDifference(data), StatisticsError);
    }
    
    
    /////////////////////////////////////////////////////////////////////////////////////
    ///                         LAG OPERATIONS TESTS                                  ///
    /////////////////////////////////////////////////////////////////////////////////////
    
    TEST_CASE("TimeSeries::Lag_BasicCalculation", "[timeseries][lag]")
    {
        TEST_PRECISION_INFO();
        Vector<Real> data({1.0, 2.0, 3.0, 4.0, 5.0});
        
        auto lag2 = Statistics::TimeSeries::Lag(data, 2);
        
        REQUIRE(lag2.size() == 3);
        REQUIRE_THAT(lag2[0], WithinAbs(3.0, 1e-10));
        REQUIRE_THAT(lag2[1], WithinAbs(4.0, 1e-10));
        REQUIRE_THAT(lag2[2], WithinAbs(5.0, 1e-10));
    }
    
    TEST_CASE("TimeSeries::Lag_Zero", "[timeseries][lag]")
    {
        TEST_PRECISION_INFO();
        Vector<Real> data({1.0, 2.0, 3.0});
        
        auto lag0 = Statistics::TimeSeries::Lag(data, 0);
        
        REQUIRE(lag0.size() == 3);
        for (int i = 0; i < 3; ++i)
            REQUIRE_THAT(lag0[i], WithinAbs(data[i], 1e-10));
    }
    
    TEST_CASE("TimeSeries::LagMatrix_BasicCalculation", "[timeseries][lag]")
    {
        TEST_PRECISION_INFO();
        Vector<Real> data({1.0, 2.0, 3.0, 4.0, 5.0});
        
        auto lagMat = Statistics::TimeSeries::LagMatrix(data, 2);
        
        // Matrix should be (5-2) x (2+1) = 3 x 3
        REQUIRE(lagMat.RowNum() == 3);
        REQUIRE(lagMat.ColNum() == 3);
        
        // Row 0: [3, 2, 1] (data at t=2, t=1, t=0)
        REQUIRE_THAT(lagMat(0, 0), WithinAbs(3.0, 1e-10));
        REQUIRE_THAT(lagMat(0, 1), WithinAbs(2.0, 1e-10));
        REQUIRE_THAT(lagMat(0, 2), WithinAbs(1.0, 1e-10));
        
        // Row 1: [4, 3, 2]
        REQUIRE_THAT(lagMat(1, 0), WithinAbs(4.0, 1e-10));
        REQUIRE_THAT(lagMat(1, 1), WithinAbs(3.0, 1e-10));
        REQUIRE_THAT(lagMat(1, 2), WithinAbs(2.0, 1e-10));
        
        // Row 2: [5, 4, 3]
        REQUIRE_THAT(lagMat(2, 0), WithinAbs(5.0, 1e-10));
        REQUIRE_THAT(lagMat(2, 1), WithinAbs(4.0, 1e-10));
        REQUIRE_THAT(lagMat(2, 2), WithinAbs(3.0, 1e-10));
    }
    
    TEST_CASE("TimeSeries::Lag_InvalidLag_Throws", "[timeseries][lag][error]")
    {
        TEST_PRECISION_INFO();
        Vector<Real> data({1.0, 2.0, 3.0});
        
        REQUIRE_THROWS_AS(Statistics::TimeSeries::Lag(data, -1), StatisticsError);
        REQUIRE_THROWS_AS(Statistics::TimeSeries::Lag(data, 3), StatisticsError);
    }
    
    
    /////////////////////////////////////////////////////////////////////////////////////
    ///                         CUMULATIVE SUM TESTS                                  ///
    /////////////////////////////////////////////////////////////////////////////////////
    
    TEST_CASE("TimeSeries::CumulativeSum_WithInitial", "[timeseries][cumsum]")
    {
        TEST_PRECISION_INFO();
        Vector<Real> diffs({2.0, 3.0, 4.0, 5.0});
        
        auto cumsum = Statistics::TimeSeries::CumulativeSum(diffs, 1.0);
        
        REQUIRE(cumsum.size() == 5);
        REQUIRE_THAT(cumsum[0], WithinAbs(1.0, 1e-10));   // initial
        REQUIRE_THAT(cumsum[1], WithinAbs(3.0, 1e-10));   // 1 + 2
        REQUIRE_THAT(cumsum[2], WithinAbs(6.0, 1e-10));   // 3 + 3
        REQUIRE_THAT(cumsum[3], WithinAbs(10.0, 1e-10)); // 6 + 4
        REQUIRE_THAT(cumsum[4], WithinAbs(15.0, 1e-10)); // 10 + 5
    }
    
    TEST_CASE("TimeSeries::CumulativeSum_WithoutInitial", "[timeseries][cumsum]")
    {
        TEST_PRECISION_INFO();
        Vector<Real> data({1.0, 2.0, 3.0, 4.0});
        
        auto cumsum = Statistics::TimeSeries::CumulativeSum(data);
        
        REQUIRE(cumsum.size() == 4);
        REQUIRE_THAT(cumsum[0], WithinAbs(1.0, 1e-10));
        REQUIRE_THAT(cumsum[1], WithinAbs(3.0, 1e-10));
        REQUIRE_THAT(cumsum[2], WithinAbs(6.0, 1e-10));
        REQUIRE_THAT(cumsum[3], WithinAbs(10.0, 1e-10));
    }
    
    TEST_CASE("TimeSeries::DifferenceAndCumulativeSum_Roundtrip", "[timeseries][integration]")
    {
        TEST_PRECISION_INFO();
        // Difference then cumsum should recover original (with first element)
        Vector<Real> original({1.0, 3.0, 6.0, 10.0, 15.0});
        
        auto diff = Statistics::TimeSeries::FirstDifference(original);
        auto recovered = Statistics::TimeSeries::CumulativeSum(diff, original[0]);
        
        REQUIRE(recovered.size() == original.size());
        for (int i = 0; i < original.size(); ++i)
            REQUIRE_THAT(recovered[i], WithinAbs(original[i], 1e-10));
    }
    
    TEST_CASE("TimeSeries::CumulativeSum_EmptyInput", "[timeseries][cumsum]")
    {
        TEST_PRECISION_INFO();
        Vector<Real> empty;
        
        auto result = Statistics::TimeSeries::CumulativeSum(empty);
        
        REQUIRE(result.size() == 0);
    }
    
}  // namespace MML::Tests::Algorithms::TimeSeriesTests

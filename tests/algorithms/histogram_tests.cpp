///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        histogram_tests.cpp                                                 ///
///  Description: Unit tests for Histogram.h - histogram computation and frequency    ///
///               analysis including binning methods, CDF, and quantiles              ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///////////////////////////////////////////////////////////////////////////////////////////

#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "algorithms/Statistics/Histogram.h"
#endif

using namespace MML;
using namespace MML::Testing;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace MML::Tests::Algorithms::HistogramTests
{
	using namespace Statistics::Histogram;

	///////////////////////////////////////////////////////////////////////////////////////////
	///                         BIN COUNT ESTIMATION TESTS                                  ///
	///////////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("Histogram::SturgesBinCount_SmallN", "[histogram][bincount]")
	{
		TEST_PRECISION_INFO();
		// Sturges: k = ceil(log2(n) + 1)
		REQUIRE(SturgesBinCount(1) == 1);    // ceil(0 + 1) = 1
		REQUIRE(SturgesBinCount(2) == 2);    // ceil(1 + 1) = 2
		REQUIRE(SturgesBinCount(4) == 3);    // ceil(2 + 1) = 3
		REQUIRE(SturgesBinCount(8) == 4);    // ceil(3 + 1) = 4
		REQUIRE(SturgesBinCount(16) == 5);   // ceil(4 + 1) = 5
		REQUIRE(SturgesBinCount(100) == 8);  // ceil(6.64 + 1) = 8
		REQUIRE(SturgesBinCount(1000) == 11); // ceil(9.97 + 1) = 11
	}

	TEST_CASE("Histogram::RiceBinCount", "[histogram][bincount]")
	{
		TEST_PRECISION_INFO();
		// Rice: k = ceil(2 * n^(1/3))
		REQUIRE(RiceBinCount(8) == 4);     // ceil(2 * 2) = 4
		REQUIRE(RiceBinCount(27) == 6);    // ceil(2 * 3) = 6
		REQUIRE(RiceBinCount(125) == 10);  // ceil(2 * 5) = 10
		REQUIRE(RiceBinCount(1000) == 20); // ceil(2 * 10) = 20
	}

	TEST_CASE("Histogram::SquareRootBinCount", "[histogram][bincount]")
	{
		TEST_PRECISION_INFO();
		// k = ceil(sqrt(n))
		REQUIRE(SquareRootBinCount(4) == 2);
		REQUIRE(SquareRootBinCount(9) == 3);
		REQUIRE(SquareRootBinCount(100) == 10);
		REQUIRE(SquareRootBinCount(1000) == 32); // ceil(31.62) = 32
	}

	TEST_CASE("Histogram::ScottBinWidth_NormalData", "[histogram][binwidth]")
	{
		TEST_PRECISION_INFO();
		// Generate data with known std dev
		Vector<Real> data(100);
		for (int i = 0; i < 100; ++i) {
			data[i] = static_cast<Real>(i);  // 0 to 99, std ~= 29.01
		}
		
		Real width = ScottBinWidth(data);
		// Scott: h = 3.49 * std * n^(-1/3)
		// For uniform 0-99: std ≈ 29.01, so h ≈ 3.49 * 29.01 * 100^(-1/3) ≈ 21.8
		REQUIRE(width > 10.0);
		REQUIRE(width < 30.0);
	}

	TEST_CASE("Histogram::FreedmanDiaconisBinWidth", "[histogram][binwidth]")
	{
		TEST_PRECISION_INFO();
		// Generate data with known IQR
		Vector<Real> data(100);
		for (int i = 0; i < 100; ++i) {
			data[i] = static_cast<Real>(i);  // 0 to 99
		}
		// IQR for 0-99 is approximately Q3 - Q1 = 74.25 - 24.75 = 49.5
		
		Real width = FreedmanDiaconisBinWidth(data);
		// FD: h = 2 * IQR * n^(-1/3) ≈ 2 * 49.5 * 0.2154 ≈ 21.3
		REQUIRE(width > 10.0);
		REQUIRE(width < 30.0);
	}

	TEST_CASE("Histogram::GetBinCount_AllMethods", "[histogram][bincount]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data(100);
		for (int i = 0; i < 100; ++i) {
			data[i] = static_cast<Real>(i);
		}
		
		int sturges = GetBinCount(data, BinningMethod::Sturges);
		int rice = GetBinCount(data, BinningMethod::Rice);
		int sqrt_rule = GetBinCount(data, BinningMethod::SquareRoot);
		int scott = GetBinCount(data, BinningMethod::Scott);
		int fd = GetBinCount(data, BinningMethod::FreedmanDiaconis);
		
		// All should be reasonable for n=100
		REQUIRE(sturges >= 5);
		REQUIRE(sturges <= 10);
		REQUIRE(rice >= 5);
		REQUIRE(rice <= 15);
		REQUIRE(sqrt_rule == 10);
		REQUIRE(scott >= 3);
		REQUIRE(scott <= 15);
		REQUIRE(fd >= 3);
		REQUIRE(fd <= 15);
	}

	///////////////////////////////////////////////////////////////////////////////////////////
	///                         CORE HISTOGRAM TESTS                                        ///
	///////////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("Histogram::ComputeHistogram_BasicUniform", "[histogram][core]")
	{
		TEST_PRECISION_INFO();
		// Uniform data 0-9, 10 bins should put ~1 in each
		Vector<Real> data(10);
		for (int i = 0; i < 10; ++i) {
			data[i] = static_cast<Real>(i) + 0.5;  // 0.5, 1.5, ..., 9.5
		}
		
		auto result = ComputeHistogram(data, 10);
		
		REQUIRE(result.numBins == 10);
		REQUIRE(result.totalCount == 10);
		
		// Total counts should sum to n
		int totalCounts = 0;
		for (int i = 0; i < 10; ++i) {
			totalCounts += result.counts[i];
		}
		REQUIRE(totalCounts == 10);
	}

	TEST_CASE("Histogram::ComputeHistogram_ConcentratedData", "[histogram][core]")
	{
		TEST_PRECISION_INFO();
		// All values in one bin
		Vector<Real> data(100);
		for (int i = 0; i < 100; ++i) {
			data[i] = 5.0 + static_cast<Real>(i % 10) * 0.01;  // 5.00 to 5.09
		}
		
		auto result = ComputeHistogram(data, 10);
		
		REQUIRE(result.numBins == 10);
		REQUIRE(result.totalCount == 100);
		
		// Sum of counts should equal total
		int totalCounts = 0;
		for (int i = 0; i < result.numBins; ++i) {
			totalCounts += result.counts[i];
		}
		REQUIRE(totalCounts == 100);
	}

	TEST_CASE("Histogram::ComputeHistogram_ConstantData", "[histogram][core]")
	{
		TEST_PRECISION_INFO();
		// All same value - edge case
		Vector<Real> data(50);
		for (int i = 0; i < 50; ++i) {
			data[i] = 42.0;
		}
		
		auto result = ComputeHistogram(data, 10);
		
		// Should handle gracefully - all in one bin
		REQUIRE(result.totalCount == 50);
		
		int totalCounts = 0;
		for (int i = 0; i < result.numBins; ++i) {
			totalCounts += result.counts[i];
		}
		REQUIRE(totalCounts == 50);
	}

	TEST_CASE("Histogram::ComputeHistogram_CustomBinEdges", "[histogram][core]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data(20);
		// Data: 0, 1, 2, ..., 19
		for (int i = 0; i < 20; ++i) {
			data[i] = static_cast<Real>(i);
		}
		
		// Custom bins: [0,5), [5,10), [10,15), [15,20]
		Vector<Real> edges(5);
		edges[0] = 0.0;
		edges[1] = 5.0;
		edges[2] = 10.0;
		edges[3] = 15.0;
		edges[4] = 20.0;
		
		auto result = ComputeHistogram(data, edges);
		
		REQUIRE(result.numBins == 4);
		REQUIRE(result.counts[0] == 5);  // 0,1,2,3,4
		REQUIRE(result.counts[1] == 5);  // 5,6,7,8,9
		REQUIRE(result.counts[2] == 5);  // 10,11,12,13,14
		REQUIRE(result.counts[3] == 5);  // 15,16,17,18,19
	}

	TEST_CASE("Histogram::ComputeHistogram_BinCenters", "[histogram][core]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data(10);
		for (int i = 0; i < 10; ++i) {
			data[i] = static_cast<Real>(i);
		}
		
		auto result = ComputeHistogram(data, 5);
		auto centers = result.GetBinCenters();
		
		REQUIRE(centers.size() == 5);
		// Centers should be evenly spaced
		Real spacing = centers[1] - centers[0];
		for (int i = 2; i < 5; ++i) {
			REQUIRE_THAT(centers[i] - centers[i-1], WithinAbs(spacing, 0.01));
		}
	}

	TEST_CASE("Histogram::ComputeHistogram_CumulativeCounts", "[histogram][core]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data(100);
		for (int i = 0; i < 100; ++i) {
			data[i] = static_cast<Real>(i);
		}
		
		auto result = ComputeHistogram(data, 10);
		auto cumCounts = result.GetCumulativeCounts();
		auto cumFreq = result.GetCumulativeFrequencies();
		
		// Last cumulative count should equal total
		REQUIRE(cumCounts[result.numBins - 1] == 100);
		
		// Last cumulative frequency should be ~1.0
		REQUIRE_THAT(cumFreq[result.numBins - 1], WithinAbs(1.0, 0.001));
		
		// Cumulative should be monotonically non-decreasing
		for (int i = 1; i < result.numBins; ++i) {
			REQUIRE(cumCounts[i] >= cumCounts[i-1]);
			REQUIRE(cumFreq[i] >= cumFreq[i-1]);
		}
	}

	TEST_CASE("Histogram::ComputeHistogram_DensityIntegratesToOne", "[histogram][core]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data(100);
		for (int i = 0; i < 100; ++i) {
			data[i] = static_cast<Real>(i);
		}
		
		auto result = ComputeHistogram(data, 10);
		
		// Integral of density should be approximately 1
		Real integral = 0.0;
		for (int i = 0; i < result.numBins; ++i) {
			integral += result.density[i] * result.binWidth;
		}
		
		REQUIRE_THAT(integral, WithinAbs(1.0, 0.01));
	}

	///////////////////////////////////////////////////////////////////////////////////////////
	///                         AUTO HISTOGRAM TESTS                                        ///
	///////////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("Histogram::HistogramSturges", "[histogram][auto]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data(100);
		for (int i = 0; i < 100; ++i) {
			data[i] = static_cast<Real>(i);
		}
		
		auto result = HistogramSturges(data);
		
		// Sturges for n=100 gives 8 bins
		REQUIRE(result.numBins == 8);
		REQUIRE(result.totalCount == 100);
	}

	TEST_CASE("Histogram::HistogramScott", "[histogram][auto]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data(100);
		for (int i = 0; i < 100; ++i) {
			data[i] = static_cast<Real>(i);
		}
		
		auto result = HistogramScott(data);
		
		REQUIRE(result.numBins > 0);
		REQUIRE(result.totalCount == 100);
	}

	TEST_CASE("Histogram::HistogramFD", "[histogram][auto]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data(100);
		for (int i = 0; i < 100; ++i) {
			data[i] = static_cast<Real>(i);
		}
		
		auto result = HistogramFD(data);
		
		REQUIRE(result.numBins > 0);
		REQUIRE(result.totalCount == 100);
	}

	TEST_CASE("Histogram::HistogramSqrt", "[histogram][auto]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data(100);
		for (int i = 0; i < 100; ++i) {
			data[i] = static_cast<Real>(i);
		}
		
		auto result = HistogramSqrt(data);
		
		// Sqrt(100) = 10
		REQUIRE(result.numBins == 10);
		REQUIRE(result.totalCount == 100);
	}

	TEST_CASE("Histogram::HistogramRice", "[histogram][auto]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data(1000);
		for (int i = 0; i < 1000; ++i) {
			data[i] = static_cast<Real>(i);
		}
		
		auto result = HistogramRice(data);
		
		// Rice for n=1000: 2 * 10 = 20
		REQUIRE(result.numBins == 20);
		REQUIRE(result.totalCount == 1000);
	}

	///////////////////////////////////////////////////////////////////////////////////////////
	///                         FREQUENCY TABLE TESTS                                       ///
	///////////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("Histogram::FrequencyTable_Basic", "[histogram][frequency]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data(10);
		data[0] = 1.0; data[1] = 2.0; data[2] = 2.0; data[3] = 3.0; data[4] = 3.0;
		data[5] = 3.0; data[6] = 4.0; data[7] = 4.0; data[8] = 4.0; data[9] = 4.0;
		
		auto result = FrequencyTable(data);
		
		REQUIRE(result.totalCount == 10);
		REQUIRE(result.uniqueCount == 4);
		
		REQUIRE(result.counts.at(1.0) == 1);
		REQUIRE(result.counts.at(2.0) == 2);
		REQUIRE(result.counts.at(3.0) == 3);
		REQUIRE(result.counts.at(4.0) == 4);
		
		REQUIRE_THAT(result.frequencies.at(1.0), WithinAbs(0.1, 1e-10));
		REQUIRE_THAT(result.frequencies.at(2.0), WithinAbs(0.2, 1e-10));
		REQUIRE_THAT(result.frequencies.at(3.0), WithinAbs(0.3, 1e-10));
		REQUIRE_THAT(result.frequencies.at(4.0), WithinAbs(0.4, 1e-10));
	}

	TEST_CASE("Histogram::FrequencyTable_GetMethods", "[histogram][frequency]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data(6);
		data[0] = 1.0; data[1] = 1.0; data[2] = 2.0;
		data[3] = 2.0; data[4] = 2.0; data[5] = 3.0;
		
		auto result = FrequencyTable(data);
		
		auto values = result.GetValues();
		auto counts = result.GetCounts();
		auto freqs = result.GetFrequencies();
		
		REQUIRE(values.size() == 3);
		REQUIRE(counts.size() == 3);
		REQUIRE(freqs.size() == 3);
		
		// Values should be sorted
		REQUIRE_THAT(values[0], WithinAbs(1.0, 1e-10));
		REQUIRE_THAT(values[1], WithinAbs(2.0, 1e-10));
		REQUIRE_THAT(values[2], WithinAbs(3.0, 1e-10));
		
		REQUIRE(counts[0] == 2);  // count of 1.0
		REQUIRE(counts[1] == 3);  // count of 2.0
		REQUIRE(counts[2] == 1);  // count of 3.0
	}

	///////////////////////////////////////////////////////////////////////////////////////////
	///                         EMPIRICAL CDF TESTS                                         ///
	///////////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("Histogram::EmpiricalCDF_Basic", "[histogram][ecdf]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data(5);
		data[0] = 1.0; data[1] = 2.0; data[2] = 3.0; data[3] = 4.0; data[4] = 5.0;
		
		auto result = EmpiricalCDF(data);
		
		REQUIRE(result.n == 5);
		REQUIRE(result.x.size() == 5);
		REQUIRE(result.cdf.size() == 5);
		
		// CDF values should be i/n for sorted unique values
		REQUIRE_THAT(result.cdf[0], WithinAbs(0.2, 1e-10));  // 1/5
		REQUIRE_THAT(result.cdf[1], WithinAbs(0.4, 1e-10));  // 2/5
		REQUIRE_THAT(result.cdf[2], WithinAbs(0.6, 1e-10));  // 3/5
		REQUIRE_THAT(result.cdf[3], WithinAbs(0.8, 1e-10));  // 4/5
		REQUIRE_THAT(result.cdf[4], WithinAbs(1.0, 1e-10));  // 5/5
	}

	TEST_CASE("Histogram::EmpiricalCDF_WithDuplicates", "[histogram][ecdf]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data(6);
		data[0] = 1.0; data[1] = 1.0; data[2] = 2.0;
		data[3] = 2.0; data[4] = 2.0; data[5] = 3.0;
		
		auto result = EmpiricalCDF(data);
		
		// Should have 3 unique values
		REQUIRE(result.x.size() == 3);
		
		REQUIRE_THAT(result.x[0], WithinAbs(1.0, 1e-10));
		REQUIRE_THAT(result.x[1], WithinAbs(2.0, 1e-10));
		REQUIRE_THAT(result.x[2], WithinAbs(3.0, 1e-10));
		
		REQUIRE_THAT(result.cdf[0], WithinAbs(2.0/6.0, 1e-10));  // 2 values <= 1.0
		REQUIRE_THAT(result.cdf[1], WithinAbs(5.0/6.0, 1e-10));  // 5 values <= 2.0
		REQUIRE_THAT(result.cdf[2], WithinAbs(1.0, 1e-10));      // all values <= 3.0
	}

	TEST_CASE("Histogram::EvaluateECDF_Basic", "[histogram][ecdf]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data(10);
		for (int i = 0; i < 10; ++i) {
			data[i] = static_cast<Real>(i);  // 0 to 9
		}
		
		REQUIRE_THAT(EvaluateECDF(data, -1.0), WithinAbs(0.0, 1e-10));   // Below all data
		REQUIRE_THAT(EvaluateECDF(data, 0.0), WithinAbs(0.1, 1e-10));    // 1 value <= 0
		REQUIRE_THAT(EvaluateECDF(data, 4.5), WithinAbs(0.5, 1e-10));    // 5 values <= 4.5
		REQUIRE_THAT(EvaluateECDF(data, 9.0), WithinAbs(1.0, 1e-10));    // All values <= 9
		REQUIRE_THAT(EvaluateECDF(data, 100.0), WithinAbs(1.0, 1e-10));  // Above all data
	}

	///////////////////////////////////////////////////////////////////////////////////////////
	///                         QUANTILE TESTS                                              ///
	///////////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("Histogram::Quantiles_Uniform", "[histogram][quantile]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data(100);
		for (int i = 0; i < 100; ++i) {
			data[i] = static_cast<Real>(i);  // 0 to 99
		}
		
		Vector<Real> probs(5);
		probs[0] = 0.0;
		probs[1] = 0.25;
		probs[2] = 0.5;
		probs[3] = 0.75;
		probs[4] = 1.0;
		
		auto quantiles = Quantiles(data, probs);
		
		REQUIRE_THAT(quantiles[0], WithinAbs(0.0, 1e-10));     // Min
		REQUIRE_THAT(quantiles[1], WithinAbs(24.75, 1e-10));   // Q1
		REQUIRE_THAT(quantiles[2], WithinAbs(49.5, 1e-10));    // Median
		REQUIRE_THAT(quantiles[3], WithinAbs(74.25, 1e-10));   // Q3
		REQUIRE_THAT(quantiles[4], WithinAbs(99.0, 1e-10));    // Max
	}

	TEST_CASE("Histogram::Quantile_Single", "[histogram][quantile]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data(100);
		for (int i = 0; i < 100; ++i) {
			data[i] = static_cast<Real>(i);
		}
		
		Real median = Quantile(data, 0.5);
		REQUIRE_THAT(median, WithinAbs(49.5, 1e-10));
	}

	TEST_CASE("Histogram::Quantiles_SmallData", "[histogram][quantile]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data(3);
		data[0] = 1.0; data[1] = 2.0; data[2] = 3.0;
		
		Vector<Real> probs(3);
		probs[0] = 0.0;
		probs[1] = 0.5;
		probs[2] = 1.0;
		
		auto quantiles = Quantiles(data, probs);
		
		REQUIRE_THAT(quantiles[0], WithinAbs(1.0, 1e-10));
		REQUIRE_THAT(quantiles[1], WithinAbs(2.0, 1e-10));
		REQUIRE_THAT(quantiles[2], WithinAbs(3.0, 1e-10));
	}

	///////////////////////////////////////////////////////////////////////////////////////////
	///                         BIN COUNT AND DIGITIZE TESTS                                ///
	///////////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("Histogram::BinCount_Basic", "[histogram][bincount]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data(10);
		for (int i = 0; i < 10; ++i) {
			data[i] = static_cast<Real>(i);  // 0 to 9
		}
		
		Vector<Real> edges(4);
		edges[0] = 0.0; edges[1] = 3.0; edges[2] = 6.0; edges[3] = 10.0;
		
		auto counts = BinCount(data, edges);
		
		REQUIRE(counts.size() == 3);
		REQUIRE(counts[0] == 3);  // 0, 1, 2
		REQUIRE(counts[1] == 3);  // 3, 4, 5
		REQUIRE(counts[2] == 4);  // 6, 7, 8, 9
	}

	TEST_CASE("Histogram::Digitize_Basic", "[histogram][digitize]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data(5);
		data[0] = -1.0;   // Below range
		data[1] = 0.5;    // Bin 0
		data[2] = 1.5;    // Bin 1
		data[3] = 2.5;    // Bin 2
		data[4] = 10.0;   // Above range
		
		Vector<Real> edges(4);
		edges[0] = 0.0; edges[1] = 1.0; edges[2] = 2.0; edges[3] = 3.0;
		
		auto indices = Digitize(data, edges);
		
		REQUIRE(indices.size() == 5);
		REQUIRE(indices[0] == -1);  // Below range
		REQUIRE(indices[1] == 0);   // Bin 0
		REQUIRE(indices[2] == 1);   // Bin 1
		REQUIRE(indices[3] == 2);   // Bin 2
		REQUIRE(indices[4] == 3);   // Above range (numBins)
	}

	TEST_CASE("Histogram::Digitize_BoundaryValues", "[histogram][digitize]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data(4);
		data[0] = 0.0;   // Left edge of first bin
		data[1] = 1.0;   // Edge between bins
		data[2] = 2.0;   // Edge between bins
		data[3] = 3.0;   // Right edge of last bin
		
		Vector<Real> edges(4);
		edges[0] = 0.0; edges[1] = 1.0; edges[2] = 2.0; edges[3] = 3.0;
		
		auto indices = Digitize(data, edges);
		
		REQUIRE(indices[0] == 0);   // 0.0 is in bin 0
		REQUIRE(indices[1] == 1);   // 1.0 is in bin 1
		REQUIRE(indices[2] == 2);   // 2.0 is in bin 2
		REQUIRE(indices[3] == 2);   // 3.0 is in last bin (inclusive)
	}

	///////////////////////////////////////////////////////////////////////////////////////////
	///                         UTILITY FUNCTION TESTS                                      ///
	///////////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("Histogram::CreateUniformBinEdges", "[histogram][utility]")
	{
		TEST_PRECISION_INFO();
		auto edges = CreateUniformBinEdges(0.0, 10.0, 5);
		
		REQUIRE(edges.size() == 6);
		REQUIRE_THAT(edges[0], WithinAbs(0.0, 1e-10));
		REQUIRE_THAT(edges[1], WithinAbs(2.0, 1e-10));
		REQUIRE_THAT(edges[2], WithinAbs(4.0, 1e-10));
		REQUIRE_THAT(edges[3], WithinAbs(6.0, 1e-10));
		REQUIRE_THAT(edges[4], WithinAbs(8.0, 1e-10));
		REQUIRE_THAT(edges[5], WithinAbs(10.0, 1e-10));
	}

	TEST_CASE("Histogram::CreateLogBinEdges", "[histogram][utility]")
	{
		TEST_PRECISION_INFO();
		auto edges = CreateLogBinEdges(1.0, 1000.0, 3);
		
		REQUIRE(edges.size() == 4);
		REQUIRE_THAT(edges[0], WithinAbs(1.0, 1e-10));
		REQUIRE_THAT(edges[1], WithinAbs(10.0, 1e-10));
		REQUIRE_THAT(edges[2], WithinAbs(100.0, 1e-10));
		REQUIRE_THAT(edges[3], WithinAbs(1000.0, 1e-10));
	}

	TEST_CASE("Histogram::CreateLogBinEdges_CustomRange", "[histogram][utility]")
	{
		TEST_PRECISION_INFO();
		auto edges = CreateLogBinEdges(0.01, 100.0, 4);
		
		REQUIRE(edges.size() == 5);
		REQUIRE_THAT(edges[0], WithinAbs(0.01, 1e-12));
		REQUIRE_THAT(edges[1], WithinAbs(0.1, 1e-12));
		REQUIRE_THAT(edges[2], WithinAbs(1.0, 1e-12));
		REQUIRE_THAT(edges[3], WithinAbs(10.0, 1e-10));
		REQUIRE_THAT(edges[4], WithinAbs(100.0, 1e-10));
	}

	///////////////////////////////////////////////////////////////////////////////////////////
	///                         ERROR HANDLING TESTS                                        ///
	///////////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("Histogram::ComputeHistogram_EmptyData_Throws", "[histogram][error]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> empty;
		
		REQUIRE_THROWS_AS(ComputeHistogram(empty, 10), StatisticsError);
	}

	TEST_CASE("Histogram::ComputeHistogram_ZeroBins_Throws", "[histogram][error]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data(10);
		for (int i = 0; i < 10; ++i) data[i] = static_cast<Real>(i);
		
		REQUIRE_THROWS_AS(ComputeHistogram(data, 0), StatisticsError);
	}

	TEST_CASE("Histogram::FrequencyTable_EmptyData_Throws", "[histogram][error]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> empty;
		
		REQUIRE_THROWS_AS(FrequencyTable(empty), StatisticsError);
	}

	TEST_CASE("Histogram::EmpiricalCDF_EmptyData_Throws", "[histogram][error]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> empty;
		
		REQUIRE_THROWS_AS(EmpiricalCDF(empty), StatisticsError);
	}

	TEST_CASE("Histogram::Quantiles_EmptyData_Throws", "[histogram][error]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> empty;
		Vector<Real> probs(1);
		probs[0] = 0.5;
		
		REQUIRE_THROWS_AS(Quantiles(empty, probs), StatisticsError);
	}

	TEST_CASE("Histogram::Quantiles_InvalidProbability_Throws", "[histogram][error]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data(10);
		for (int i = 0; i < 10; ++i) data[i] = static_cast<Real>(i);
		
		Vector<Real> badProbs(1);
		badProbs[0] = 1.5;  // > 1.0
		
		REQUIRE_THROWS_AS(Quantiles(data, badProbs), StatisticsError);
	}

	TEST_CASE("Histogram::CreateLogBinEdges_NonPositiveMin_Throws", "[histogram][error]")
	{
		TEST_PRECISION_INFO();
		REQUIRE_THROWS_AS(CreateLogBinEdges(0.0, 10.0, 5), StatisticsError);
		REQUIRE_THROWS_AS(CreateLogBinEdges(-1.0, 10.0, 5), StatisticsError);
	}

}  // namespace MML::Tests::Algorithms::HistogramTests

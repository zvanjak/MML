#include <catch2/catch_all.hpp>

#include "TestMatchers.h"
#include "TestPrecision.h"

#include <mml/algorithms/Statistics/Histogram.h>

using namespace MML;
using namespace MML::Testing;
using Catch::Matchers::WithinAbs;

namespace MML::Tests::Algorithms::HistogramTests
{
	using namespace Statistics::Histogram;

	TEST_CASE("HistogramResult defaults to an empty result", "[statistics][histogram]")
	{
		HistogramResult result;

		REQUIRE(result.numBins == 0);
		REQUIRE(result.totalCount == 0);
		REQUIRE(result.binWidth == REAL(0.0));
		REQUIRE(result.GetBinCenters().size() == 0);
		REQUIRE(result.GetCumulativeCounts().size() == 0);
		REQUIRE(result.GetCumulativeFrequencies().size() == 0);
	}

	TEST_CASE("HistogramResult derives centers and cumulative values", "[statistics][histogram]")
	{
		TEST_PRECISION_INFO();
		HistogramResult result;
		result.binEdges = Vector<Real>({REAL(0.0), REAL(1.0), REAL(2.0), REAL(4.0)});
		result.counts = Vector<int>({2, 1, 3});
		result.frequencies = Vector<Real>({REAL(2.0) / REAL(6.0), REAL(1.0) / REAL(6.0), REAL(3.0) / REAL(6.0)});
		result.numBins = 3;
		result.totalCount = 6;

		auto centers = result.GetBinCenters();
		auto cumulativeCounts = result.GetCumulativeCounts();
		auto cumulativeFrequencies = result.GetCumulativeFrequencies();

		REQUIRE_THAT(centers[0], WithinAbs(REAL(0.5), TOL(1e-12, 1e-5)));
		REQUIRE_THAT(centers[1], WithinAbs(REAL(1.5), TOL(1e-12, 1e-5)));
		REQUIRE_THAT(centers[2], WithinAbs(REAL(3.0), TOL(1e-12, 1e-5)));
		REQUIRE(cumulativeCounts == Vector<int>({2, 3, 6}));
		REQUIRE_THAT(cumulativeFrequencies[0], WithinAbs(REAL(2.0) / REAL(6.0), TOL(1e-12, 1e-5)));
		REQUIRE_THAT(cumulativeFrequencies[1], WithinAbs(REAL(3.0) / REAL(6.0), TOL(1e-12, 1e-5)));
		REQUIRE_THAT(cumulativeFrequencies[2], WithinAbs(REAL(1.0), TOL(1e-12, 1e-5)));
	}

	TEST_CASE("FrequencyTable counts discrete values and relative frequencies", "[statistics][histogram][frequency]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(4.0), REAL(2.0), REAL(4.0), REAL(1.0), REAL(2.0), REAL(4.0)});

		auto result = FrequencyTable(data);

		REQUIRE(result.totalCount == 6);
		REQUIRE(result.uniqueCount == 3);
		REQUIRE(result.counts.at(REAL(1.0)) == 1);
		REQUIRE(result.counts.at(REAL(2.0)) == 2);
		REQUIRE(result.counts.at(REAL(4.0)) == 3);
		REQUIRE_THAT(result.frequencies.at(REAL(1.0)), WithinAbs(REAL(1.0) / REAL(6.0), TOL(1e-12, 1e-5)));
		REQUIRE_THAT(result.frequencies.at(REAL(2.0)), WithinAbs(REAL(2.0) / REAL(6.0), TOL(1e-12, 1e-5)));
		REQUIRE_THAT(result.frequencies.at(REAL(4.0)), WithinAbs(REAL(3.0) / REAL(6.0), TOL(1e-12, 1e-5)));
	}

	TEST_CASE("FrequencyTable accessors preserve sorted value order", "[statistics][histogram][frequency]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(3.0), REAL(1.0), REAL(3.0), REAL(2.0), REAL(2.0), REAL(2.0)});

		auto result = FrequencyTable(data);
		auto values = result.GetValues();
		auto counts = result.GetCounts();
		auto frequencies = result.GetFrequencies();

		REQUIRE(values == Vector<Real>({REAL(1.0), REAL(2.0), REAL(3.0)}));
		REQUIRE(counts == Vector<int>({1, 3, 2}));
		REQUIRE_THAT(frequencies[0], WithinAbs(REAL(1.0) / REAL(6.0), TOL(1e-12, 1e-5)));
		REQUIRE_THAT(frequencies[1], WithinAbs(REAL(3.0) / REAL(6.0), TOL(1e-12, 1e-5)));
		REQUIRE_THAT(frequencies[2], WithinAbs(REAL(2.0) / REAL(6.0), TOL(1e-12, 1e-5)));
	}

	TEST_CASE("FrequencyTable rejects empty and non-finite data", "[statistics][histogram][frequency][error]")
	{
		Vector<Real> empty;
		Vector<Real> withNaN({REAL(1.0), std::numeric_limits<Real>::quiet_NaN()});
		Vector<Real> withInfinity({REAL(1.0), std::numeric_limits<Real>::infinity()});

		REQUIRE_THROWS_AS(FrequencyTable(empty), StatisticsError);
		REQUIRE_THROWS_AS(FrequencyTable(withNaN), StatisticsError);
		REQUIRE_THROWS_AS(FrequencyTable(withInfinity), StatisticsError);
	}

	TEST_CASE("BinningMethod exposes all migrated package choices", "[statistics][histogram]")
	{
		REQUIRE(BinningMethod::Sturges != BinningMethod::Scott);
		REQUIRE(BinningMethod::FreedmanDiaconis != BinningMethod::SquareRoot);
		REQUIRE(BinningMethod::SquareRoot != BinningMethod::Rice);
	}

	TEST_CASE("Histogram bin-count rules are exact at integer boundaries", "[statistics][histogram][binning]")
	{
		REQUIRE(SturgesBinCount(1) == 1);
		REQUIRE(SturgesBinCount(8) == 4);
		REQUIRE(SturgesBinCount(1024) == 11);

		REQUIRE(RiceBinCount(8) == 4);
		REQUIRE(RiceBinCount(27) == 6);
		REQUIRE(RiceBinCount(125) == 10);
		REQUIRE(RiceBinCount(1000) == 20);

		REQUIRE(SquareRootBinCount(1) == 1);
		REQUIRE(SquareRootBinCount(9) == 3);
		REQUIRE(SquareRootBinCount(1000) == 32);

		constexpr int firstIntegerBeyondFloatPrecision = 16777217;
		REQUIRE(SturgesBinCount(firstIntegerBeyondFloatPrecision) == 26);
		REQUIRE(RiceBinCount(firstIntegerBeyondFloatPrecision) == 513);
		REQUIRE(SquareRootBinCount(firstIntegerBeyondFloatPrecision) == 4097);
	}

	TEST_CASE("Histogram Scott and Freedman-Diaconis widths match formulas", "[statistics][histogram][binning]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(0.0), REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0)});
		Real sampleStdDev = std::sqrt(REAL(2.5));
		Real scale = std::pow(REAL(5.0), -REAL(1.0) / REAL(3.0));

		REQUIRE_THAT(ScottBinWidth(data), WithinAbs(REAL(3.49) * sampleStdDev * scale, TOL(1e-12, 1e-5)));
		REQUIRE_THAT(FreedmanDiaconisBinWidth(data), WithinAbs(REAL(4.0) * scale, TOL(1e-12, 1e-5)));
	}

	TEST_CASE("Histogram GetBinCount handles every method and constant data", "[statistics][histogram][binning]")
	{
		Vector<Real> data(100);
		for (int i = 0; i < data.size(); ++i)
			data[i] = static_cast<Real>(i);

		REQUIRE(GetBinCount(data, BinningMethod::Sturges) == 8);
		REQUIRE(GetBinCount(data, BinningMethod::Rice) == 10);
		REQUIRE(GetBinCount(data, BinningMethod::SquareRoot) == 10);
		REQUIRE(GetBinCount(data, BinningMethod::Scott) == 5);
		REQUIRE(GetBinCount(data, BinningMethod::FreedmanDiaconis) == 5);

		Vector<Real> constant({REAL(7.0), REAL(7.0), REAL(7.0)});
		REQUIRE(GetBinCount(constant, BinningMethod::Sturges) == 3);
		REQUIRE(GetBinCount(constant, BinningMethod::Rice) == 3);
		REQUIRE(GetBinCount(constant, BinningMethod::SquareRoot) == 2);
		REQUIRE(GetBinCount(constant, BinningMethod::Scott) == 3);
		REQUIRE(GetBinCount(constant, BinningMethod::FreedmanDiaconis) == 3);
	}

	TEST_CASE("Histogram binning rules reject invalid input", "[statistics][histogram][binning][error]")
	{
		Vector<Real> empty;
		Vector<Real> nonFinite({REAL(1.0), std::numeric_limits<Real>::quiet_NaN(), REAL(2.0)});

		REQUIRE_THROWS_AS(SturgesBinCount(0), StatisticsError);
		REQUIRE_THROWS_AS(RiceBinCount(-1), StatisticsError);
		REQUIRE_THROWS_AS(SquareRootBinCount(0), StatisticsError);
		REQUIRE_THROWS_AS(ScottBinWidth(empty), StatisticsError);
		REQUIRE_THROWS_AS(FreedmanDiaconisBinWidth(empty), StatisticsError);
		REQUIRE_THROWS_AS(GetBinCount(empty, BinningMethod::Sturges), StatisticsError);
		REQUIRE_THROWS_AS(ScottBinWidth(nonFinite), StatisticsError);
		REQUIRE_THROWS_AS(FreedmanDiaconisBinWidth(nonFinite), StatisticsError);
		REQUIRE_THROWS_AS(GetBinCount(nonFinite, BinningMethod::Rice), StatisticsError);
		REQUIRE(Detail::StableCeilToBinCount(static_cast<double>(std::numeric_limits<int>::max()))
		        == std::numeric_limits<int>::max());
		REQUIRE_THROWS_AS(
			Detail::StableCeilToBinCount(static_cast<double>(std::numeric_limits<int>::max()) + 0.25),
			StatisticsError);
	}

	TEST_CASE("ComputeHistogram uses exact uniform edges and includes the maximum", "[statistics][histogram][compute]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(0.0), REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0)});

		auto result = ComputeHistogram(data, 4);

		REQUIRE(result.numBins == 4);
		REQUIRE(result.totalCount == 5);
		REQUIRE(result.binEdges == Vector<Real>({REAL(0.0), REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0)}));
		REQUIRE(result.counts == Vector<int>({1, 1, 1, 2}));
		REQUIRE_THAT(result.binWidth, WithinAbs(REAL(1.0), TOL(1e-12, 1e-5)));

		Real densityIntegral = 0.0;
		for (int i = 0; i < result.numBins; ++i)
			densityIntegral += result.density[i] * result.binWidth;
		REQUIRE_THAT(densityIntegral, WithinAbs(REAL(1.0), TOL(1e-12, 1e-5)));
	}

	TEST_CASE("ComputeHistogram collapses constant data to one bin", "[statistics][histogram][compute]")
	{
		Vector<Real> data({REAL(42.0), REAL(42.0), REAL(42.0)});

		auto result = ComputeHistogram(data, 10);

		REQUIRE(result.numBins == 1);
		REQUIRE(result.totalCount == 3);
		REQUIRE(result.counts == Vector<int>({3}));
		REQUIRE(result.frequencies == Vector<Real>({REAL(1.0)}));
		REQUIRE(result.density == Vector<Real>({REAL(1.0)}));
		REQUIRE(result.binEdges == Vector<Real>({REAL(41.5), REAL(42.5)}));
	}

	TEST_CASE("ComputeHistogram supports constant data at finite numeric limits", "[statistics][histogram][compute]")
	{
		for (Real value : {std::numeric_limits<Real>::lowest(), std::numeric_limits<Real>::max()}) {
			Vector<Real> data({value, value});

			auto result = ComputeHistogram(data, 2);

			REQUIRE(result.numBins == 1);
			REQUIRE(result.counts == Vector<int>({2}));
			REQUIRE(std::isfinite(result.binEdges[0]));
			REQUIRE(std::isfinite(result.binEdges[1]));
			REQUIRE(result.binEdges[0] < result.binEdges[1]);
			REQUIRE(std::isfinite(result.GetBinCenters()[0]));
			REQUIRE(std::isfinite(result.density[0]));
		}
	}

	TEST_CASE("ComputeHistogram preserves the package exact-epsilon boundary", "[statistics][histogram][compute]")
	{
		Vector<Real> data({REAL(0.0), Constants::Eps});

		auto result = ComputeHistogram(data, 2);

		REQUIRE(result.numBins == 2);
		REQUIRE(result.counts == Vector<int>({1, 1}));
		REQUIRE(result.binEdges[0] == REAL(0.0));
		REQUIRE(result.binEdges[2] == Constants::Eps);
	}

	TEST_CASE("ComputeHistogram assigns an exact generated edge to the following bin", "[statistics][histogram][compute]")
	{
		Real minValue = REAL(0.1);
		Real maxValue = REAL(0.7);
		Real generatedEdge = minValue + (maxValue - minValue) / REAL(3.0);
		Vector<Real> data({minValue, generatedEdge, maxValue});

		auto result = ComputeHistogram(data, 3);

		REQUIRE(result.binEdges[1] == generatedEdge);
		REQUIRE(result.counts == Vector<int>({1, 1, 1}));
	}

	TEST_CASE("ComputeHistogram custom edges preserve package out-of-range semantics", "[statistics][histogram][compute]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data({-REAL(1.0), REAL(0.0), REAL(0.5), REAL(1.0), REAL(2.0), REAL(4.0), REAL(5.0)});
		Vector<Real> edges({REAL(0.0), REAL(1.0), REAL(2.0), REAL(4.0)});

		auto result = ComputeHistogram(data, edges);

		REQUIRE(result.totalCount == 7);
		REQUIRE(result.numBins == 3);
		REQUIRE(result.counts == Vector<int>({2, 1, 2}));
		REQUIRE_THAT(result.frequencies[0], WithinAbs(REAL(2.0) / REAL(7.0), TOL(1e-12, 1e-5)));
		REQUIRE_THAT(result.frequencies[1], WithinAbs(REAL(1.0) / REAL(7.0), TOL(1e-12, 1e-5)));
		REQUIRE_THAT(result.frequencies[2], WithinAbs(REAL(2.0) / REAL(7.0), TOL(1e-12, 1e-5)));
		REQUIRE_THAT(result.density[0], WithinAbs(REAL(2.0) / REAL(7.0), TOL(1e-12, 1e-5)));
		REQUIRE_THAT(result.density[1], WithinAbs(REAL(1.0) / REAL(7.0), TOL(1e-12, 1e-5)));
		REQUIRE_THAT(result.density[2], WithinAbs(REAL(1.0) / REAL(7.0), TOL(1e-12, 1e-5)));
		REQUIRE_THAT(result.binWidth, WithinAbs(REAL(4.0) / REAL(3.0), TOL(1e-12, 1e-5)));
	}

	TEST_CASE("ComputeHistogram custom variable-width density integrates to one", "[statistics][histogram][compute]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(0.0), REAL(0.5), REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0)});
		Vector<Real> edges({REAL(0.0), REAL(1.0), REAL(2.0), REAL(4.0)});

		auto result = ComputeHistogram(data, edges);
		Real densityIntegral = 0.0;
		for (int i = 0; i < result.numBins; ++i)
			densityIntegral += result.density[i] * (result.binEdges[i + 1] - result.binEdges[i]);

		REQUIRE_THAT(densityIntegral, WithinAbs(REAL(1.0), TOL(1e-12, 1e-5)));
		REQUIRE_THAT(result.GetCumulativeFrequencies()[result.numBins - 1],
		             WithinAbs(REAL(1.0), TOL(1e-12, 1e-5)));
	}

	TEST_CASE("ComputeHistogram automatic and convenience APIs use migrated binning rules", "[statistics][histogram][compute]")
	{
		Vector<Real> data(100);
		for (int i = 0; i < data.size(); ++i)
			data[i] = static_cast<Real>(i);

		REQUIRE(ComputeHistogramAuto(data).numBins == SturgesBinCount(data.size()));
		REQUIRE(HistogramSturges(data).numBins == GetBinCount(data, BinningMethod::Sturges));
		REQUIRE(HistogramScott(data).numBins == GetBinCount(data, BinningMethod::Scott));
		REQUIRE(HistogramFD(data).numBins == GetBinCount(data, BinningMethod::FreedmanDiaconis));
		REQUIRE(HistogramSqrt(data).numBins == GetBinCount(data, BinningMethod::SquareRoot));
		REQUIRE(HistogramRice(data).numBins == GetBinCount(data, BinningMethod::Rice));
	}

	TEST_CASE("ComputeHistogram rejects invalid data and edges", "[statistics][histogram][compute][error]")
	{
		Vector<Real> empty;
		Vector<Real> data({REAL(1.0), REAL(2.0), REAL(3.0)});
		Vector<Real> nonFiniteData({REAL(1.0), std::numeric_limits<Real>::infinity()});
		Vector<Real> tooFewEdges({REAL(0.0)});
		Vector<Real> repeatedEdges({REAL(0.0), REAL(1.0), REAL(1.0)});
		Vector<Real> descendingEdges({REAL(0.0), REAL(2.0), REAL(1.0)});
		Vector<Real> nonFiniteEdges({REAL(0.0), std::numeric_limits<Real>::quiet_NaN(), REAL(2.0)});
		Vector<Real> overflowingSpan({-std::numeric_limits<Real>::max(), std::numeric_limits<Real>::max()});

		REQUIRE_THROWS_AS(ComputeHistogram(empty, 2), StatisticsError);
		REQUIRE_THROWS_AS(ComputeHistogram(data, 0), StatisticsError);
		REQUIRE_THROWS_AS(ComputeHistogram(nonFiniteData, 2), StatisticsError);
		REQUIRE_THROWS_AS(ComputeHistogram(data, tooFewEdges), StatisticsError);
		REQUIRE_THROWS_AS(ComputeHistogram(data, repeatedEdges), StatisticsError);
		REQUIRE_THROWS_AS(ComputeHistogram(data, descendingEdges), StatisticsError);
		REQUIRE_THROWS_AS(ComputeHistogram(data, nonFiniteEdges), StatisticsError);
		REQUIRE_THROWS_AS(ComputeHistogram(data, overflowingSpan), StatisticsError);
		REQUIRE_THROWS_AS(ComputeHistogram(data, std::numeric_limits<int>::max()), StatisticsError);

		if (std::numeric_limits<Real>::denorm_min() > 0.0) {
			Vector<Real> tinyEdges({REAL(0.0), std::numeric_limits<Real>::denorm_min()});
			Vector<Real> zero({REAL(0.0)});
			REQUIRE_THROWS_AS(ComputeHistogram(zero, tinyEdges), StatisticsError);
		}
	}

	TEST_CASE("EmpiricalCDF is right-continuous at sorted unique values", "[statistics][histogram][ecdf]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(3.0), REAL(1.0), REAL(2.0), REAL(2.0), REAL(1.0), REAL(2.0)});

		auto result = EmpiricalCDF(data);

		REQUIRE(result.n == 6);
		REQUIRE(result.x == Vector<Real>({REAL(1.0), REAL(2.0), REAL(3.0)}));
		REQUIRE_THAT(result.cdf[0], WithinAbs(REAL(2.0) / REAL(6.0), TOL(1e-12, 1e-5)));
		REQUIRE_THAT(result.cdf[1], WithinAbs(REAL(5.0) / REAL(6.0), TOL(1e-12, 1e-5)));
		REQUIRE_THAT(result.cdf[2], WithinAbs(REAL(1.0), TOL(1e-12, 1e-5)));
	}

	TEST_CASE("EvaluateECDF handles boundaries and infinities", "[statistics][histogram][ecdf]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(0.0), REAL(1.0), REAL(2.0), REAL(3.0)});

		REQUIRE_THAT(EvaluateECDF(data, -std::numeric_limits<Real>::infinity()), WithinAbs(REAL(0.0), TOL(1e-12, 1e-5)));
		REQUIRE_THAT(EvaluateECDF(data, REAL(0.0)), WithinAbs(REAL(0.25), TOL(1e-12, 1e-5)));
		REQUIRE_THAT(EvaluateECDF(data, REAL(1.5)), WithinAbs(REAL(0.5), TOL(1e-12, 1e-5)));
		REQUIRE_THAT(EvaluateECDF(data, REAL(3.0)), WithinAbs(REAL(1.0), TOL(1e-12, 1e-5)));
		REQUIRE_THAT(EvaluateECDF(data, std::numeric_limits<Real>::infinity()), WithinAbs(REAL(1.0), TOL(1e-12, 1e-5)));
	}

	TEST_CASE("Quantiles preserve probability order and match Percentile", "[statistics][histogram][quantile]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(9.0), REAL(1.0), REAL(5.0), REAL(3.0), REAL(7.0)});
		Vector<Real> probabilities({REAL(1.0), REAL(0.25), REAL(0.5), REAL(0.0), REAL(0.75)});

		auto quantiles = Quantiles(data, probabilities);

		for (int i = 0; i < probabilities.size(); ++i) {
			REQUIRE_THAT(
				quantiles[i],
				WithinAbs(Statistics::Percentile(data, probabilities[i] * REAL(100.0)), TOL(1e-12, 1e-5)));
		}
		REQUIRE_THAT(Quantile(data, REAL(0.5)), WithinAbs(Statistics::Median(data), TOL(1e-12, 1e-5)));
	}

	TEST_CASE("Quantile interpolation remains finite at numeric extremes", "[statistics][histogram][quantile]")
	{
		Vector<Real> data({std::numeric_limits<Real>::lowest(), std::numeric_limits<Real>::max()});
		Vector<Real> probabilities({REAL(0.25), REAL(0.5), REAL(0.75)});

		auto quantiles = Quantiles(data, probabilities);

		for (int i = 0; i < probabilities.size(); ++i) {
			Real percentile = Statistics::Percentile(data, probabilities[i] * REAL(100.0));
			REQUIRE(std::isfinite(quantiles[i]));
			REQUIRE(std::isfinite(percentile));
			REQUIRE(quantiles[i] == percentile);
		}
		REQUIRE(std::abs(Quantile(data, REAL(0.5))) <= REAL(1.0));
	}

	TEST_CASE("ECDF and quantiles reject invalid inputs", "[statistics][histogram][ecdf][quantile][error]")
	{
		Vector<Real> empty;
		Vector<Real> data({REAL(1.0), REAL(2.0), REAL(3.0)});
		Vector<Real> nonFiniteData({REAL(1.0), std::numeric_limits<Real>::infinity()});
		Vector<Real> noProbabilities;
		Vector<Real> badProbabilities({-REAL(0.1), REAL(1.1)});
		Vector<Real> nanProbability({std::numeric_limits<Real>::quiet_NaN()});

		REQUIRE_THROWS_AS(EmpiricalCDF(empty), StatisticsError);
		REQUIRE_THROWS_AS(EmpiricalCDF(nonFiniteData), StatisticsError);
		REQUIRE_THROWS_AS(EvaluateECDF(empty, REAL(0.0)), StatisticsError);
		REQUIRE_THROWS_AS(EvaluateECDF(nonFiniteData, REAL(0.0)), StatisticsError);
		REQUIRE_THROWS_AS(EvaluateECDF(data, std::numeric_limits<Real>::quiet_NaN()), StatisticsError);
		REQUIRE_THROWS_AS(Quantiles(empty, Vector<Real>({REAL(0.5)})), StatisticsError);
		REQUIRE_THROWS_AS(Quantiles(nonFiniteData, Vector<Real>({REAL(0.5)})), StatisticsError);
		REQUIRE_THROWS_AS(Quantiles(data, noProbabilities), StatisticsError);
		REQUIRE_THROWS_AS(Quantiles(data, badProbabilities), StatisticsError);
		REQUIRE_THROWS_AS(Quantiles(data, nanProbability), StatisticsError);
	}

	TEST_CASE("CreateUniformBinEdges produces evenly spaced edges", "[statistics][histogram]")
	{
		Vector<Real> edges = CreateUniformBinEdges(REAL(0.0), REAL(10.0), 5);
		REQUIRE(edges.size() == 6);
		REQUIRE_THAT(edges[0], WithinAbs(REAL(0.0), REAL(1e-12)));
		REQUIRE_THAT(edges[3], WithinAbs(REAL(6.0), REAL(1e-12)));
		REQUIRE_THAT(edges[5], WithinAbs(REAL(10.0), REAL(1e-12)));
		REQUIRE_THROWS_AS(CreateUniformBinEdges(REAL(0.0), REAL(1.0), 0), StatisticsError);
	}

	TEST_CASE("CreateLogBinEdges produces logarithmically spaced edges", "[statistics][histogram]")
	{
		Vector<Real> edges = CreateLogBinEdges(REAL(1.0), REAL(1000.0), 3);
		REQUIRE(edges.size() == 4);
		REQUIRE_THAT(edges[0], WithinAbs(REAL(1.0), REAL(1e-9)));
		REQUIRE_THAT(edges[1], WithinAbs(REAL(10.0), REAL(1e-9)));
		REQUIRE_THAT(edges[2], WithinAbs(REAL(100.0), REAL(1e-9)));
		REQUIRE_THAT(edges[3], WithinAbs(REAL(1000.0), REAL(1e-9)));
		REQUIRE_THROWS_AS(CreateLogBinEdges(REAL(0.0), REAL(10.0), 3), StatisticsError);
		REQUIRE_THROWS_AS(CreateLogBinEdges(REAL(10.0), REAL(1.0), 3), StatisticsError);
	}

	TEST_CASE("BinCount counts samples per bin", "[statistics][histogram]")
	{
		Vector<Real> data({ REAL(0.5), REAL(1.5), REAL(2.5), REAL(2.6) });
		Vector<Real> edges({ REAL(0.0), REAL(1.0), REAL(2.0), REAL(3.0) });
		Vector<int> counts = BinCount(data, edges);
		REQUIRE(counts.size() == 3);
		REQUIRE(counts[0] == 1);
		REQUIRE(counts[1] == 1);
		REQUIRE(counts[2] == 2);
		REQUIRE_THROWS_AS(BinCount(Vector<Real>(), edges), StatisticsError);
	}

	TEST_CASE("Digitize assigns bin indices with out-of-range sentinels", "[statistics][histogram]")
	{
		Vector<Real> data({ REAL(-1.0), REAL(0.5), REAL(1.5), REAL(5.0) });
		Vector<Real> edges({ REAL(0.0), REAL(1.0), REAL(2.0), REAL(3.0) });
		Vector<int> idx = Digitize(data, edges);
		REQUIRE(idx.size() == 4);
		REQUIRE(idx[0] == -1);
		REQUIRE(idx[1] == 0);
		REQUIRE(idx[2] == 1);
		REQUIRE(idx[3] == 3);
	}
}
#include <catch2/catch_all.hpp>
#include "TestPrecision.h"
#include "TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "algorithms/Statistics.h"
#include "algorithms/Statistics/Distributions.h"
#include "algorithms/Statistics/StatisticsRank.h"
#endif

using namespace MML;
using namespace MML::Testing;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace MML::Tests::Algorithms::RankCorrelationTests
{

	/*********************************************************************/
	/*****              Rank Computation Tests                       *****/
	/*********************************************************************/
	TEST_CASE("Statistics::ComputeRanks_SimpleSequence", "[statistics][rank]")
	{
		TEST_PRECISION_INFO();
		// Already sorted: 1, 2, 3, 4, 5
		Vector<Real> data({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0)});
		
		Vector<Real> ranks = Statistics::ComputeRanks(data);
		
		REQUIRE(ranks.size() == 5);
		REQUIRE_THAT(ranks[0], WithinAbs(REAL(1.0), REAL(1e-10)));
		REQUIRE_THAT(ranks[1], WithinAbs(REAL(2.0), REAL(1e-10)));
		REQUIRE_THAT(ranks[2], WithinAbs(REAL(3.0), REAL(1e-10)));
		REQUIRE_THAT(ranks[3], WithinAbs(REAL(4.0), REAL(1e-10)));
		REQUIRE_THAT(ranks[4], WithinAbs(REAL(5.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::ComputeRanks_UnsortedData", "[statistics][rank]")
	{
		TEST_PRECISION_INFO();
		// Unsorted: 3, 1, 4, 1.5, 2  -> sorted order: 1, 1.5, 2, 3, 4
		// Ranks:    4, 1, 5, 2, 3
		Vector<Real> data({REAL(3.0), REAL(1.0), REAL(4.0), REAL(1.5), REAL(2.0)});
		
		Vector<Real> ranks = Statistics::ComputeRanks(data);
		
		REQUIRE_THAT(ranks[0], WithinAbs(REAL(4.0), REAL(1e-10)));  // 3.0 is 4th smallest
		REQUIRE_THAT(ranks[1], WithinAbs(REAL(1.0), REAL(1e-10)));  // 1.0 is smallest
		REQUIRE_THAT(ranks[2], WithinAbs(REAL(5.0), REAL(1e-10)));  // 4.0 is largest
		REQUIRE_THAT(ranks[3], WithinAbs(REAL(2.0), REAL(1e-10)));  // 1.5 is 2nd smallest
		REQUIRE_THAT(ranks[4], WithinAbs(REAL(3.0), REAL(1e-10)));  // 2.0 is 3rd smallest
	}

	TEST_CASE("Statistics::ComputeRanks_WithTies", "[statistics][rank]")
	{
		TEST_PRECISION_INFO();
		// Data with ties: 1, 2, 2, 3
		// Without ties: ranks would be 1, 2, 3, 4
		// With ties at positions 2,3: average rank = (2+3)/2 = 2.5
		Vector<Real> data({REAL(1.0), REAL(2.0), REAL(2.0), REAL(3.0)});
		
		Vector<Real> ranks = Statistics::ComputeRanks(data);
		
		REQUIRE_THAT(ranks[0], WithinAbs(REAL(1.0), REAL(1e-10)));
		REQUIRE_THAT(ranks[1], WithinAbs(REAL(2.5), REAL(1e-10)));  // Average of 2 and 3
		REQUIRE_THAT(ranks[2], WithinAbs(REAL(2.5), REAL(1e-10)));  // Average of 2 and 3
		REQUIRE_THAT(ranks[3], WithinAbs(REAL(4.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::ComputeRanks_AllTied", "[statistics][rank]")
	{
		TEST_PRECISION_INFO();
		// All equal: average rank for all
		Vector<Real> data({REAL(5.0), REAL(5.0), REAL(5.0), REAL(5.0)});
		
		Vector<Real> ranks = Statistics::ComputeRanks(data);
		
		// Average rank = (1+2+3+4)/4 = 2.5
		for (int i = 0; i < 4; i++) {
			REQUIRE_THAT(ranks[i], WithinAbs(REAL(2.5), REAL(1e-10)));
		}
	}

	TEST_CASE("Statistics::ComputeRanks_EmptyVector", "[statistics][rank]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data;
		
		Vector<Real> ranks = Statistics::ComputeRanks(data);
		
		REQUIRE(ranks.size() == 0);
	}

	/*********************************************************************/
	/*****              Spearman Correlation Tests                   *****/
	/*********************************************************************/
	TEST_CASE("Statistics::SpearmanCorrelation_PerfectMonotonic", "[statistics][rank]")
	{
		TEST_PRECISION_INFO();
		// Perfect positive monotonic relationship (but not linear)
		Vector<Real> x({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0)});
		Vector<Real> y({REAL(1.0), REAL(4.0), REAL(9.0), REAL(16.0), REAL(25.0)});  // y = xÂ²
		
		Real rho = Statistics::SpearmanCorrelation(x, y);
		
		// Perfect monotonic -> rho = 1
		REQUIRE_THAT(rho, WithinAbs(REAL(1.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::SpearmanCorrelation_PerfectNegativeMonotonic", "[statistics][rank]")
	{
		TEST_PRECISION_INFO();
		// Perfect negative monotonic
		Vector<Real> x({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0)});
		Vector<Real> y({REAL(25.0), REAL(16.0), REAL(9.0), REAL(4.0), REAL(1.0)});
		
		Real rho = Statistics::SpearmanCorrelation(x, y);
		
		REQUIRE_THAT(rho, WithinAbs(REAL(-1.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::SpearmanCorrelation_NoMonotonicRelationship", "[statistics][rank]")
	{
		TEST_PRECISION_INFO();
		// U-shaped relationship: not monotonic
		Vector<Real> x({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0)});
		Vector<Real> y({REAL(4.0), REAL(1.0), REAL(0.0), REAL(1.0), REAL(4.0)});  // y = (x-3)Â²
		
		Real rho = Statistics::SpearmanCorrelation(x, y);
		
		// Should be close to zero (not exactly due to discretization)
		REQUIRE(std::abs(rho) < REAL(0.5));
	}

	TEST_CASE("Statistics::SpearmanCorrelation_WithTies", "[statistics][rank]")
	{
		TEST_PRECISION_INFO();
		// Data with ties
		Vector<Real> x({REAL(1.0), REAL(2.0), REAL(2.0), REAL(3.0), REAL(4.0)});
		Vector<Real> y({REAL(1.0), REAL(2.0), REAL(3.0), REAL(3.0), REAL(4.0)});
		
		Real rho = Statistics::SpearmanCorrelation(x, y);
		
		// Should be positive (monotonically increasing overall)
		REQUIRE(rho > REAL(0.8));
		REQUIRE(rho <= REAL(1.0));
	}

	TEST_CASE("Statistics::SpearmanCorrelation_SelfCorrelation", "[statistics][rank]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> x({REAL(3.0), REAL(1.0), REAL(4.0), REAL(1.5), REAL(9.0)});
		
		Real rho = Statistics::SpearmanCorrelation(x, x);
		
		REQUIRE_THAT(rho, WithinAbs(REAL(1.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::SpearmanCorrelation_RangeIsMinusOneToOne", "[statistics][rank]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> x({REAL(2.3), REAL(5.1), REAL(1.2), REAL(8.7), REAL(4.4)});
		Vector<Real> y({REAL(7.8), REAL(3.2), REAL(9.1), REAL(2.5), REAL(6.6)});
		
		Real rho = Statistics::SpearmanCorrelation(x, y);
		
		REQUIRE(rho >= REAL(-1.0));
		REQUIRE(rho <= REAL(1.0));
	}

	TEST_CASE("Statistics::SpearmanCorrelation_DifferentSizes_ThrowsException", "[statistics][rank][error]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> x({REAL(1.0), REAL(2.0), REAL(3.0)});
		Vector<Real> y({REAL(1.0), REAL(2.0)});
		
		REQUIRE_THROWS_AS(Statistics::SpearmanCorrelation(x, y), StatisticsError);
	}

	TEST_CASE("Statistics::SpearmanCorrelation_SingleElement_ThrowsException", "[statistics][rank][error]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> x({REAL(1.0)});
		Vector<Real> y({REAL(2.0)});
		
		REQUIRE_THROWS_AS(Statistics::SpearmanCorrelation(x, y), StatisticsError);
	}

	TEST_CASE("Statistics::SpearmanCorrelationWithTest_HighCorrelation", "[statistics][rank]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> x({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0),
		                REAL(6.0), REAL(7.0), REAL(8.0), REAL(9.0), REAL(10.0)});
		Vector<Real> y({REAL(2.0), REAL(4.0), REAL(6.0), REAL(8.0), REAL(10.0),
		                REAL(12.0), REAL(14.0), REAL(16.0), REAL(18.0), REAL(20.0)});
		
		auto result = Statistics::SpearmanCorrelationWithTest(x, y);
		
		REQUIRE_THAT(result.rho, WithinAbs(REAL(1.0), REAL(1e-10)));
		REQUIRE(result.n == 10);
		REQUIRE(result.zScore > REAL(2.0));  // Significant
		REQUIRE(result.IsSignificant(0.05));
	}

	TEST_CASE("Statistics::SpearmanCorrelationWithTest_InsufficientData_ThrowsException", "[statistics][rank][error]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> x({REAL(1.0), REAL(2.0)});
		Vector<Real> y({REAL(2.0), REAL(4.0)});
		
		REQUIRE_THROWS_AS(Statistics::SpearmanCorrelationWithTest(x, y), StatisticsError);
	}

	/*********************************************************************/
	/*****              Kendall Correlation Tests                    *****/
	/*********************************************************************/
	TEST_CASE("Statistics::KendallCorrelation_PerfectConcordance", "[statistics][rank]")
	{
		TEST_PRECISION_INFO();
		// All pairs concordant
		Vector<Real> x({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0)});
		Vector<Real> y({REAL(10.0), REAL(20.0), REAL(30.0), REAL(40.0), REAL(50.0)});
		
		Real tau = Statistics::KendallCorrelation(x, y);
		
		REQUIRE_THAT(tau, WithinAbs(REAL(1.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::KendallCorrelation_PerfectDiscordance", "[statistics][rank]")
	{
		TEST_PRECISION_INFO();
		// All pairs discordant
		Vector<Real> x({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0)});
		Vector<Real> y({REAL(50.0), REAL(40.0), REAL(30.0), REAL(20.0), REAL(10.0)});
		
		Real tau = Statistics::KendallCorrelation(x, y);
		
		REQUIRE_THAT(tau, WithinAbs(REAL(-1.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::KendallCorrelation_NoRelationship", "[statistics][rank]")
	{
		TEST_PRECISION_INFO();
		// Designed to have balanced concordant/discordant pairs
		Vector<Real> x({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0)});
		Vector<Real> y({REAL(1.0), REAL(4.0), REAL(2.0), REAL(3.0)});
		
		Real tau = Statistics::KendallCorrelation(x, y);
		
		// tau should be between -1 and 1, close to 0 or moderate
		REQUIRE(tau >= REAL(-1.0));
		REQUIRE(tau <= REAL(1.0));
	}

	TEST_CASE("Statistics::KendallCorrelation_WithTiesInX", "[statistics][rank]")
	{
		TEST_PRECISION_INFO();
		// Ties in x only
		Vector<Real> x({REAL(1.0), REAL(2.0), REAL(2.0), REAL(3.0)});
		Vector<Real> y({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0)});
		
		Real tau = Statistics::KendallCorrelation(x, y);
		
		// Should be positive but less than 1 due to tie
		REQUIRE(tau > REAL(0.5));
		REQUIRE(tau < REAL(1.0));
	}

	TEST_CASE("Statistics::KendallCorrelation_WithTiesInY", "[statistics][rank]")
	{
		TEST_PRECISION_INFO();
		// Ties in y only
		Vector<Real> x({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0)});
		Vector<Real> y({REAL(1.0), REAL(2.0), REAL(2.0), REAL(3.0)});
		
		Real tau = Statistics::KendallCorrelation(x, y);
		
		REQUIRE(tau > REAL(0.5));
		REQUIRE(tau < REAL(1.0));
	}

	TEST_CASE("Statistics::KendallCorrelation_AllTied", "[statistics][rank]")
	{
		TEST_PRECISION_INFO();
		// All tied in both x and y
		Vector<Real> x({REAL(5.0), REAL(5.0), REAL(5.0), REAL(5.0)});
		Vector<Real> y({REAL(3.0), REAL(3.0), REAL(3.0), REAL(3.0)});
		
		Real tau = Statistics::KendallCorrelation(x, y);
		
		// No valid pairs -> tau = 0
		REQUIRE_THAT(tau, WithinAbs(REAL(0.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::KendallCorrelation_SelfCorrelation", "[statistics][rank]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> x({REAL(3.0), REAL(1.0), REAL(4.0), REAL(1.5), REAL(9.0)});
		
		Real tau = Statistics::KendallCorrelation(x, x);
		
		REQUIRE_THAT(tau, WithinAbs(REAL(1.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::KendallCorrelation_RangeIsMinusOneToOne", "[statistics][rank]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> x({REAL(2.3), REAL(5.1), REAL(1.2), REAL(8.7), REAL(4.4)});
		Vector<Real> y({REAL(7.8), REAL(3.2), REAL(9.1), REAL(2.5), REAL(6.6)});
		
		Real tau = Statistics::KendallCorrelation(x, y);
		
		REQUIRE(tau >= REAL(-1.0));
		REQUIRE(tau <= REAL(1.0));
	}

	TEST_CASE("Statistics::KendallCorrelation_DifferentSizes_ThrowsException", "[statistics][rank][error]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> x({REAL(1.0), REAL(2.0), REAL(3.0)});
		Vector<Real> y({REAL(1.0), REAL(2.0)});
		
		REQUIRE_THROWS_AS(Statistics::KendallCorrelation(x, y), StatisticsError);
	}

	TEST_CASE("Statistics::KendallCorrelation_SingleElement_ThrowsException", "[statistics][rank][error]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> x({REAL(1.0)});
		Vector<Real> y({REAL(2.0)});
		
		REQUIRE_THROWS_AS(Statistics::KendallCorrelation(x, y), StatisticsError);
	}

	TEST_CASE("Statistics::KendallCorrelationWithTest_SignificantCorrelation", "[statistics][rank]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> x({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0),
		                REAL(6.0), REAL(7.0), REAL(8.0), REAL(9.0), REAL(10.0)});
		Vector<Real> y({REAL(2.0), REAL(4.0), REAL(6.0), REAL(8.0), REAL(10.0),
		                REAL(12.0), REAL(14.0), REAL(16.0), REAL(18.0), REAL(20.0)});
		
		auto result = Statistics::KendallCorrelationWithTest(x, y);
		
		REQUIRE_THAT(result.rho, WithinAbs(REAL(1.0), REAL(1e-10)));  // tau stored in rho
		REQUIRE(result.n == 10);
		REQUIRE(result.zScore > REAL(2.0));
		REQUIRE(result.IsSignificant(0.05));
	}

	TEST_CASE("Statistics::KendallCorrelationWithTest_InsufficientData_ThrowsException", "[statistics][rank][error]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> x({REAL(1.0), REAL(2.0)});
		Vector<Real> y({REAL(2.0), REAL(4.0)});
		
		REQUIRE_THROWS_AS(Statistics::KendallCorrelationWithTest(x, y), StatisticsError);
	}

	/*********************************************************************/
	/*****          Spearman vs Kendall Comparison Tests             *****/
	/*********************************************************************/
	TEST_CASE("Statistics::SpearmanVsKendall_SameSign", "[statistics][rank]")
	{
		TEST_PRECISION_INFO();
		// Both should agree on sign (direction of relationship)
		Vector<Real> x({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0)});
		Vector<Real> y({REAL(2.0), REAL(3.0), REAL(5.0), REAL(4.0), REAL(6.0)});
		
		Real rho = Statistics::SpearmanCorrelation(x, y);
		Real tau = Statistics::KendallCorrelation(x, y);
		
		// Same sign
		REQUIRE(rho * tau > 0);
		
		// Both positive
		REQUIRE(rho > 0);
		REQUIRE(tau > 0);
		
		// Spearman typically has larger magnitude than Kendall
		// (for same data, |rho| >= |tau| is common but not guaranteed)
	}

	TEST_CASE("Statistics::RankCorrelation_PerfectMonotonic_AllAgree", "[statistics][rank]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> x({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0)});
		Vector<Real> y({REAL(1.0), REAL(8.0), REAL(27.0), REAL(64.0), REAL(125.0)});  // y = xÂ³
		
		Real pearson = Statistics::PearsonCorrelation(x, y);
		Real spearman = Statistics::SpearmanCorrelation(x, y);
		Real kendall = Statistics::KendallCorrelation(x, y);
		
		// All three should be 1 for monotonic relationship
		// (Pearson may be less than 1 due to non-linearity)
		REQUIRE_THAT(spearman, WithinAbs(REAL(1.0), REAL(1e-10)));
		REQUIRE_THAT(kendall, WithinAbs(REAL(1.0), REAL(1e-10)));
		
		// Pearson less than 1 for non-linear but monotonic
		REQUIRE(pearson > REAL(0.9));
		REQUIRE(pearson <= REAL(1.0));
	}

	/******************************************************************************/
	/*****             Detailed API Tests - Rank Correlations                 *****/
	/******************************************************************************/

	TEST_CASE("Statistics::SpearmanCorrelationDetailed - values match simple API", "[statistics][rank][Detailed]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> x({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0)});
		Vector<Real> y({REAL(1.0), REAL(4.0), REAL(9.0), REAL(16.0), REAL(25.0)});

		Real simple = Statistics::SpearmanCorrelation(x, y);
		auto detailed = Statistics::SpearmanCorrelationDetailed(x, y);

		REQUIRE(detailed.IsSuccess());
		REQUIRE(detailed.algorithm_name == "SpearmanCorrelation");
		REQUIRE_THAT(detailed.rho, WithinAbs(simple, REAL(1e-12)));
		REQUIRE(detailed.n == 5);
		REQUIRE(detailed.elapsed_time_ms >= 0.0);
	}

	TEST_CASE("Statistics::SpearmanCorrelationWithTestDetailed - values match simple API", "[statistics][rank][Detailed]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> x({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0), REAL(6.0), REAL(7.0), REAL(8.0), REAL(9.0), REAL(10.0)});
		Vector<Real> y({REAL(2.0), REAL(4.0), REAL(6.0), REAL(8.0), REAL(10.0), REAL(12.0), REAL(14.0), REAL(16.0), REAL(18.0), REAL(20.0)});

		auto simple = Statistics::SpearmanCorrelationWithTest(x, y);
		auto detailed = Statistics::SpearmanCorrelationWithTestDetailed(x, y);

		REQUIRE(detailed.IsSuccess());
		REQUIRE(detailed.algorithm_name == "SpearmanCorrelationWithTest");
		REQUIRE_THAT(detailed.rho, WithinAbs(simple.rho, REAL(1e-12)));
		REQUIRE_THAT(detailed.zScore, WithinAbs(simple.zScore, REAL(1e-12)));
		REQUIRE(detailed.n == simple.n);
	}

	TEST_CASE("Statistics::KendallCorrelationDetailed - values match simple API", "[statistics][rank][Detailed]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> x({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0)});
		Vector<Real> y({REAL(1.0), REAL(4.0), REAL(9.0), REAL(16.0), REAL(25.0)});

		Real simple = Statistics::KendallCorrelation(x, y);
		auto detailed = Statistics::KendallCorrelationDetailed(x, y);

		REQUIRE(detailed.IsSuccess());
		REQUIRE(detailed.algorithm_name == "KendallCorrelation");
		REQUIRE_THAT(detailed.rho, WithinAbs(simple, REAL(1e-12)));
		REQUIRE(detailed.n == 5);
	}

	TEST_CASE("Statistics::KendallCorrelationWithTestDetailed - values match simple API", "[statistics][rank][Detailed]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> x({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0), REAL(6.0), REAL(7.0), REAL(8.0), REAL(9.0), REAL(10.0)});
		Vector<Real> y({REAL(2.0), REAL(4.0), REAL(6.0), REAL(8.0), REAL(10.0), REAL(12.0), REAL(14.0), REAL(16.0), REAL(18.0), REAL(20.0)});

		auto simple = Statistics::KendallCorrelationWithTest(x, y);
		auto detailed = Statistics::KendallCorrelationWithTestDetailed(x, y);

		REQUIRE(detailed.IsSuccess());
		REQUIRE(detailed.algorithm_name == "KendallCorrelationWithTest");
		REQUIRE_THAT(detailed.rho, WithinAbs(simple.rho, REAL(1e-12)));
		REQUIRE_THAT(detailed.zScore, WithinAbs(simple.zScore, REAL(1e-12)));
		REQUIRE(detailed.n == simple.n);
	}

	TEST_CASE("Statistics::SpearmanCorrelationDetailed - ConvertToStatus on invalid input", "[statistics][rank][Detailed]")
	{
		Vector<Real> x({REAL(1.0)});
		Vector<Real> y({REAL(2.0)});
		Statistics::StatisticsConfig config;
		config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;

		auto result = Statistics::SpearmanCorrelationDetailed(x, y, config);

		REQUIRE_FALSE(result.IsSuccess());
		REQUIRE(result.status == AlgorithmStatus::InvalidInput);
		REQUIRE_FALSE(result.error_message.empty());
	}

	TEST_CASE("Statistics::SpearmanCorrelationDetailed - Propagate policy throws", "[statistics][rank][Detailed]")
	{
		Vector<Real> x({REAL(1.0)});
		Vector<Real> y({REAL(2.0)});
		Statistics::StatisticsConfig config;
		config.exception_policy = EvaluationExceptionPolicy::Propagate;

		REQUIRE_THROWS(Statistics::SpearmanCorrelationDetailed(x, y, config));
	}

	TEST_CASE("Statistics::KendallCorrelationDetailed - ConvertToStatus on invalid input", "[statistics][rank][Detailed]")
	{
		Vector<Real> x({REAL(1.0), REAL(2.0)});
		Vector<Real> y({REAL(3.0)});
		Statistics::StatisticsConfig config;
		config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;

		auto result = Statistics::KendallCorrelationDetailed(x, y, config);

		REQUIRE_FALSE(result.IsSuccess());
		REQUIRE(result.status == AlgorithmStatus::InvalidInput);
	}

}
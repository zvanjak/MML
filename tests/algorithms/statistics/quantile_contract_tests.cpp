#include "../../TestPrecision.h"
#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/algorithms/Statistics/Distributions.h>
#endif
#include <catch2/catch_all.hpp>

using namespace MML;
using namespace MML::Statistics;
using Catch::Matchers::WithinRel;

namespace {
	void CheckQuantile(Real actual, Real expected) {
		CHECK(std::isfinite(actual));
		CHECK_THAT(actual, WithinRel(expected,
			std::max(Precision::EqualityTolerance, Precision::DistributionQuantileRelativeTolerance)));
	}
	void CheckProbability(Real actual, Real expected) {
		CHECK(std::isfinite(actual));
		CHECK_THAT(actual, WithinRel(expected,
			std::max(Precision::EqualityTolerance, Precision::DistributionQuantileProbabilityTolerance)));
	}
}

TEST_CASE("Continuous quantiles repair the three independently derived audit witnesses", "[statistics][quantiles]") {
	CheckQuantile(FDistribution(1, 1).criticalValue(Real(0.05)),
		1 / std::pow(std::tan(Constants::PI * Real(0.05) / 2), MML::Real(2)));
	CheckQuantile(BetaDistribution(Real(0.5), 1).inverseCdf(Real(1e-6)), Real(1e-12));
	Real p = Real(0.999999);
	CheckQuantile(TDistribution(2).inverseCdf(p), (2 * p - 1) / std::sqrt(2 * p * (1 - p)));
}

TEST_CASE("Student quantiles retain symmetric small tails and central resolution", "[statistics][quantiles]") {
	for (Real p : {Real(1e-150), Real(1e-20), Real(1e-6), Real(0.1), Real(0.25)}) {
		if (p <= 0) continue;
		CAPTURE(p);
		const Real magnitude = (1 - 2 * p) / (std::sqrt(p) * std::sqrt(2 * (1 - p)));
		const Real t = TDistribution(2).inverseCdf(p);
		CheckQuantile(t, -magnitude);
		CheckProbability(TDistribution(2).cdf(t), p);
		CheckProbability(TDistribution(2).sf(-t), p);
		CheckQuantile(TDistribution(1).inverseCdf(p), -1 / std::tan(Constants::PI * p));
	}
	for (Real p : {std::nextafter(Real(0.5), Real(0)), std::nextafter(Real(0.5), Real(1))}) {
		CAPTURE(p);
		CheckQuantile(TDistribution(2).inverseCdf(p), (2 * p - 1) / std::sqrt(2 * p * (1 - p)));
	}
	CHECK(TDistribution(2).inverseCdf(Real(0.5)) == 0);
}

TEST_CASE("Student two-sided critical values never subtract a tiny significance from one", "[statistics][quantiles]") {
	for (Real alpha : {Real(1e-150), Real(1e-20), Real(0.05), Real(0.9), std::nextafter(Real(1), Real(0))}) {
		if (alpha <= 0) continue;
		CAPTURE(alpha);
		Real expected = (1 - alpha) / (std::sqrt(alpha) * std::sqrt(1 - alpha / 2));
		Real t = TDistribution(2).criticalValue(alpha);
		CheckQuantile(t, expected);
		CheckProbability(TDistribution(2).twoTailedPValue(t), alpha);
	}
}

TEST_CASE("F critical values expand their bracket over many orders of magnitude", "[statistics][quantiles]") {
	for (Real alpha : {Real(1e-150), Real(1e-20), Real(1e-6), Real(0.05), Real(0.9)}) {
		if (alpha <= 0) continue;
		CAPTURE(alpha);
		Real cotangent = 1 / std::tan(Constants::PI * alpha / 2);
		if (!std::isfinite(cotangent * cotangent)) {
			CHECK_THROWS_AS(FDistribution(1, 1).criticalValue(alpha), ConvergenceError);
			continue;
		}
		Real f = FDistribution(1, 1).criticalValue(alpha);
		CheckQuantile(f, cotangent * cotangent);
		CheckProbability(FDistribution(1, 1).rightTailPValue(f), alpha);
	}
	CheckProbability(FDistribution(2, 2).cdf(std::numeric_limits<Real>::max()), Real(1));
	CheckProbability(FDistribution(2, 2).rightTailPValue(std::numeric_limits<Real>::max()),
		1 / std::numeric_limits<Real>::max());
}

TEST_CASE("Chi-square critical values use the direct tail and a verified upper bracket", "[statistics][quantiles]") {
	for (Real alpha : {Real(1e-300), Real(1e-20), Real(0.05), Real(0.9), std::nextafter(Real(1), Real(0))}) {
		if (alpha <= 0) continue;
		CAPTURE(alpha);
		Real x = ChiSquareDistribution(2).criticalValue(alpha);
		CheckQuantile(x, -2 * std::log(alpha));
		CheckProbability(ChiSquareDistribution(2).rightTailPValue(x), alpha);
	}
}

TEST_CASE("Gamma quantiles preserve closed-form answers across scale and probability ranges", "[statistics][quantiles]") {
	for (Real scale : {TOL(1e-200, 1e-20), Real(1), TOL(1e200, 1e20)}) {
		if (!(std::isfinite(scale) && scale > 0)) continue;
		for (Real p : {Real(1e-100), Real(1e-6), Real(0.1), Real(0.9), std::nextafter(Real(1), Real(0))}) {
			if (p <= 0) continue;
			CAPTURE(scale, p);
			CheckQuantile(GammaDistribution(1, scale).inverseCdf(p), -scale * std::log1p(-p));
		}
	}
}

TEST_CASE("Beta quantiles have no arbitrary interior support clamps", "[statistics][quantiles]") {
	for (Real p : {Real(1e-150), Real(1e-20), Real(1e-6), Real(0.1), Real(0.9)}) {
		if (p <= 0) continue;
		CAPTURE(p);
		CheckQuantile(BetaDistribution(Real(0.5), 1).inverseCdf(p), p * p);
		CheckQuantile(BetaDistribution(2, 1).inverseCdf(p), std::sqrt(p));
	}
}

TEST_CASE("Quantile domain checks preserve Gamma Beta clipping but reject NaN", "[statistics][quantiles]") {
	for (Real p : {-std::numeric_limits<Real>::infinity(), Real(-1), Real(0)}) {
		CHECK(GammaDistribution(2).inverseCdf(p) == 0);
		CHECK(BetaDistribution(2, 3).inverseCdf(p) == 0);
	}
	for (Real p : {Real(1), Real(2), std::numeric_limits<Real>::infinity()}) {
		CHECK(GammaDistribution(2).inverseCdf(p) == std::numeric_limits<Real>::infinity());
		CHECK(BetaDistribution(2, 3).inverseCdf(p) == 1);
	}
	for (Real p : {-std::numeric_limits<Real>::infinity(), Real(-1), Real(0), Real(1), Real(2),
		std::numeric_limits<Real>::infinity(), std::numeric_limits<Real>::quiet_NaN()}) {
		CHECK_THROWS_AS(TDistribution(2).inverseCdf(p), StatisticsError);
		CHECK_THROWS_AS(TDistribution(2).criticalValue(p), StatisticsError);
		CHECK_THROWS_AS(ChiSquareDistribution(2).criticalValue(p), StatisticsError);
		CHECK_THROWS_AS(FDistribution(1, 1).criticalValue(p), StatisticsError);
	}
	CHECK_THROWS_AS(GammaDistribution(2).inverseCdf(std::numeric_limits<Real>::quiet_NaN()), StatisticsError);
	CHECK_THROWS_AS(BetaDistribution(2, 3).inverseCdf(std::numeric_limits<Real>::quiet_NaN()), StatisticsError);
}

TEST_CASE("Unresolvable interior quantiles report convergence failure instead of clamping", "[statistics][quantiles]") {
	CHECK_THROWS_AS(BetaDistribution(Real(0.5), 1).inverseCdf(std::numeric_limits<Real>::min()), ConvergenceError);
	CHECK_THROWS_AS(BetaDistribution(1, Real(0.5)).inverseCdf(std::nextafter(Real(1), Real(0))), ConvergenceError);
	CHECK_THROWS_AS(FDistribution(1, 1).criticalValue(std::numeric_limits<Real>::min()), ConvergenceError);
	CHECK_THROWS_AS(GammaDistribution(1, std::numeric_limits<Real>::max()).inverseCdf(Real(0.9)), ConvergenceError);
	CHECK_THROWS_AS(BetaDistribution(2, 1).inverseCdf(std::numeric_limits<Real>::denorm_min()), ConvergenceError);
}

TEST_CASE("Quantile refinement and bracketing budgets report exact exhaustion", "[statistics][quantiles]") {
	try {
		StatisticsDetail::PositiveQuantile(Real(0.25), [](Real x) { return x; }, true, 1, false, 1);
		FAIL("Expected exhausted refinement budget");
	} catch (const ConvergenceError& error) {
		CHECK(error.iterations() == 1);
		CHECK(error.residual() > Precision::DistributionQuantileProbabilityTolerance);
	}
	try {
		StatisticsDetail::PositiveQuantile(Real(0.1), [](Real x) { return 1 / (1 + x); }, false, 1, true,
			Precision::DistributionQuantileMaxIterations, 0);
		FAIL("Expected exhausted bracketing budget");
	} catch (const ConvergenceError& error) {
		CHECK(error.iterations() == 0);
		CHECK(error.residual() > Precision::DistributionQuantileProbabilityTolerance);
	}
	CHECK_THROWS_AS(StatisticsDetail::PositiveQuantile(Real(0.1),
		[](Real) { return std::numeric_limits<Real>::quiet_NaN(); }, true), NumericalMethodError);
}

TEST_CASE("Quantile beta kernel cannot silently accept an exhausted fraction", "[statistics][quantiles][special]") {
	try {
		SpecialFunctions::SpecialFunctionsDetail::BetaFraction(Real(0.5), Real(1e12), Real(1e12));
		FAIL("Expected exhausted incomplete beta budget");
	} catch (const ConvergenceError& error) {
		CHECK(error.iterations() == Precision::IncompleteBetaMaxIterations);
		CHECK(error.residual() > Precision::Eps);
	}
	CHECK_THAT(SpecialFunctions::RegularizedBetaI(Real(0.5), Real(1e12), Real(1e12)),
		WithinRel(Real(0.5), Precision::EqualityTolerance));
	CHECK_THAT(BetaDistribution(Real(1e12), Real(1e12)).inverseCdf(Real(0.5)),
		WithinRel(Real(0.5), Precision::EqualityTolerance));
	// Float's x grid at 1e6 cannot resolve the distribution's requested
	// relative tail accuracy. 1000 still exercises the large-shape kernel.
	const Real shape = TOL(Real(1e6), Real(1000));
	const Real median = GammaDistribution(shape).inverseCdf(Real(0.5));
	CHECK(median < shape);
	CHECK(median > shape - 1);
}

TEST_CASE("Continuous quantiles match independent high-precision nonintegral-shape references", "[statistics][quantiles]") {
	// mpmath, 100 decimal digits and 220 bisections, independent of the MML kernels.
	CheckQuantile(BetaDistribution(Real(0.25), Real(2.5)).inverseCdf(Real(0.01)), Real(3.14989955591680043771807673017e-9));
	CheckQuantile(BetaDistribution(Real(2.5), Real(5.5)).inverseCdf(Real(0.5)), Real(0.296243090737929216948475912164));
	CheckQuantile(BetaDistribution(20, 3).inverseCdf(Real(0.99)), Real(0.979432958161774813204450117484));
	CheckQuantile(GammaDistribution(Real(0.25)).inverseCdf(Real(0.01)), Real(6.74969792955846769128227522726e-9));
	CheckQuantile(GammaDistribution(Real(2.5)).inverseCdf(Real(0.5)), Real(2.17573009554776365857905388831));
	CheckQuantile(GammaDistribution(20).inverseCdf(Real(0.99)), Real(31.8453698757822315261912728947));
	CheckQuantile(TDistribution(3).inverseCdf(Real(1e-12)), Real(-10331.1082442924861343106912452));
	CheckQuantile(TDistribution(5).inverseCdf(Real(0.01)), Real(-3.36492999890721859276001278247));
	CheckQuantile(TDistribution(30).inverseCdf(Real(0.1)), Real(-1.31041502539139557822300917976));
}

TEST_CASE("Continuous quantiles remain monotone across shapes and reciprocal F tails", "[statistics][quantiles]") {
	for (int df : {1, 2, 5, 30}) {
		Real previousT = -std::numeric_limits<Real>::infinity();
		Real previousChi = std::numeric_limits<Real>::infinity();
		Real previousGamma = 0;
		Real previousBeta = 0;
		for (int index = 1; index < 32; ++index) {
			const Real p = Real(index) / 32;
			CAPTURE(df, p);
			const Real t = TDistribution(df).inverseCdf(p);
			const Real chi = ChiSquareDistribution(df).criticalValue(p);
			const Real gamma = GammaDistribution(Real(df) / 2).inverseCdf(p);
			const Real beta = BetaDistribution(Real(df) / 2, Real(2.5)).inverseCdf(p);
			CHECK(t > previousT);
			CHECK(chi < previousChi);
			CHECK(gamma > previousGamma);
			CHECK(beta > previousBeta);
			CheckProbability(TDistribution(df).cdf(t), p);
			CheckProbability(ChiSquareDistribution(df).rightTailPValue(chi), p);
			CheckProbability(GammaDistribution(Real(df) / 2).cdf(gamma), p);
			CheckProbability(BetaDistribution(Real(df) / 2, Real(2.5)).cdf(beta), p);
			previousT = t;
			previousChi = chi;
			previousGamma = gamma;
			previousBeta = beta;
		}
	}
	for (Real alpha : {Real(0.01), Real(0.1), Real(0.5), Real(0.9), Real(0.99)})
		CheckQuantile(FDistribution(3, 7).criticalValue(alpha) * FDistribution(7, 3).criticalValue(1 - alpha), 1);
}

TEST_CASE("Iterative quantiles reject invalid public shape and scale state explicitly", "[statistics][quantiles]") {
	for (Real invalid : {Real(-1), Real(0), std::numeric_limits<Real>::infinity(), std::numeric_limits<Real>::quiet_NaN()}) {
		GammaDistribution gamma(2);
		BetaDistribution beta(2, 3);
		gamma.shape = invalid;
		beta.alpha = invalid;
		CHECK_THROWS_AS(gamma.inverseCdf(Real(0.5)), StatisticsError);
		CHECK_THROWS_AS(beta.inverseCdf(Real(0.5)), StatisticsError);
		gamma.shape = 2;
		beta.alpha = 2;
		gamma.scale = invalid;
		beta.beta = invalid;
		CHECK_THROWS_AS(gamma.inverseCdf(Real(0.5)), StatisticsError);
		CHECK_THROWS_AS(beta.inverseCdf(Real(0.5)), StatisticsError);
	}
	for (Real invalid : {Real(-1), Real(0), std::numeric_limits<Real>::infinity(), std::numeric_limits<Real>::quiet_NaN()}) {
		CHECK_THROWS_AS(SpecialFunctions::RegularizedBetaI(Real(0.5), invalid, 2), DomainError);
		CHECK_THROWS_AS(SpecialFunctions::RegularizedBetaI(Real(0.5), 2, invalid), DomainError);
	}
	CHECK_THROWS_AS(SpecialFunctions::RegularizedBetaI(std::numeric_limits<Real>::quiet_NaN(), 2, 3), DomainError);
	CHECK_THROWS_AS(SpecialFunctions::RegularizedBetaI(Real(0.5),
		std::numeric_limits<Real>::max(), std::numeric_limits<Real>::max()), NumericalMethodError);
}

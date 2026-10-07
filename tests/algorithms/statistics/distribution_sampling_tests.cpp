///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML) Tests                            ///
///                                                                                   ///
///  File:        distribution_sampling_tests.cpp                                     ///
///  Description: Tests for distribution sample() variate generation                  ///
///               (continuous + discrete), plus RNG-engine flexibility                 ///
///////////////////////////////////////////////////////////////////////////////////////////

#include "../../TestPrecision.h"
#include "../../../mml/algorithms/Statistics/Distributions.h"
#include "../../../mml/algorithms/Statistics/DiscreteDistributions.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <random>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace MML;
using namespace MML::Statistics;
using Catch::Approx;

namespace MML::Tests::Algorithms::Statistics::SamplingTests {

// Draw n samples from a continuous distribution with a deterministic engine.
template<class Dist>
static std::vector<double> draw(const Dist& d, unsigned seed, int n)
{
	std::mt19937 gen(seed);
	std::vector<double> out;
	out.reserve(n);
	for (int i = 0; i < n; ++i) out.push_back(static_cast<double>(d.sample(gen)));
	return out;
}

static void meanVar(const std::vector<double>& v, double& m, double& var)
{
	double s = 0.0;
	for (double x : v) s += x;
	m = s / v.size();
	double s2 = 0.0;
	for (double x : v) s2 += (x - m) * (x - m);
	var = s2 / v.size();
}

static double median(std::vector<double> v)
{
	std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
	return v[v.size() / 2];
}

constexpr int N = 200000;

TEST_CASE("Sampling::Normal_moments_and_range", "[sampling][normal]") {
	TEST_PRECISION_INFO();
	NormalDistribution d(2.0, 3.0);
	auto s = draw(d, 101, N);
	double m, var; meanVar(s, m, var);
	REQUIRE(m == Approx(2.0).margin(0.05));
	REQUIRE(var == Approx(9.0).epsilon(0.05));
}

TEST_CASE("Sampling::Exponential_moments_and_positivity", "[sampling][exponential]") {
	TEST_PRECISION_INFO();
	ExponentialDistribution d(2.0);  // rate 2 -> mean 0.5
	auto s = draw(d, 102, N);
	double m, var; meanVar(s, m, var);
	REQUIRE(m == Approx(0.5).epsilon(0.03));
	REQUIRE(var == Approx(0.25).epsilon(0.05));
	REQUIRE(*std::min_element(s.begin(), s.end()) >= 0.0);
}

TEST_CASE("Sampling::Uniform_moments_and_bounds", "[sampling][uniform]") {
	TEST_PRECISION_INFO();
	UniformDistribution d(1.0, 5.0);
	auto s = draw(d, 103, N);
	double m, var; meanVar(s, m, var);
	REQUIRE(m == Approx(3.0).epsilon(0.01));
	REQUIRE(var == Approx(16.0 / 12.0).epsilon(0.05));
	REQUIRE(*std::min_element(s.begin(), s.end()) >= 1.0);
	REQUIRE(*std::max_element(s.begin(), s.end()) <= 5.0);
}

TEST_CASE("Sampling::Gamma_moments_and_positivity", "[sampling][gamma]") {
	TEST_PRECISION_INFO();
	GammaDistribution d(3.0, 2.0);  // mean 6, var 12
	auto s = draw(d, 104, N);
	double m, var; meanVar(s, m, var);
	REQUIRE(m == Approx(6.0).epsilon(0.03));
	REQUIRE(var == Approx(12.0).epsilon(0.06));
	REQUIRE(*std::min_element(s.begin(), s.end()) > 0.0);
}

TEST_CASE("Sampling::Gamma_shape_below_one", "[sampling][gamma]") {
	TEST_PRECISION_INFO();
	GammaDistribution d(0.5, 2.0);  // mean 1, var 2
	auto s = draw(d, 105, N);
	double m, var; meanVar(s, m, var);
	REQUIRE(m == Approx(1.0).epsilon(0.04));
	REQUIRE(*std::min_element(s.begin(), s.end()) > 0.0);
}

TEST_CASE("Sampling::Beta_moments_and_bounds", "[sampling][beta]") {
	TEST_PRECISION_INFO();
	BetaDistribution d(2.0, 5.0);
	auto s = draw(d, 106, N);
	double m, var; meanVar(s, m, var);
	REQUIRE(m == Approx(2.0 / 7.0).epsilon(0.02));
	REQUIRE(var == Approx(10.0 / (49.0 * 8.0)).epsilon(0.08));
	REQUIRE(*std::min_element(s.begin(), s.end()) > 0.0);
	REQUIRE(*std::max_element(s.begin(), s.end()) < 1.0);
}

TEST_CASE("Sampling::ChiSquare_moments", "[sampling][chisquare]") {
	TEST_PRECISION_INFO();
	ChiSquareDistribution d(4);
	auto s = draw(d, 107, N);
	double m, var; meanVar(s, m, var);
	REQUIRE(m == Approx(4.0).epsilon(0.03));
	REQUIRE(var == Approx(8.0).epsilon(0.06));
}

TEST_CASE("Sampling::TDistribution_mean_near_zero", "[sampling][tdist]") {
	TEST_PRECISION_INFO();
	TDistribution d(10);
	auto s = draw(d, 108, N);
	double m, var; meanVar(s, m, var);
	REQUIRE(m == Approx(0.0).margin(0.05));
	REQUIRE(var == Approx(10.0 / 8.0).epsilon(0.10));
}

TEST_CASE("Sampling::FDistribution_mean_and_positivity", "[sampling][fdist]") {
	TEST_PRECISION_INFO();
	FDistribution d(10, 10);  // mean df2/(df2-2) = 1.25
	auto s = draw(d, 109, N);
	double m, var; meanVar(s, m, var);
	REQUIRE(m == Approx(1.25).epsilon(0.05));
	REQUIRE(*std::min_element(s.begin(), s.end()) > 0.0);
}

TEST_CASE("Sampling::Cauchy_median_near_location", "[sampling][cauchy]") {
	TEST_PRECISION_INFO();
	CauchyDistribution d(0.0, 1.0);
	auto s = draw(d, 110, N);
	REQUIRE(median(s) == Approx(0.0).margin(0.03));  // mean/var undefined; median is robust
}

TEST_CASE("Sampling::Logistic_moments", "[sampling][logistic]") {
	TEST_PRECISION_INFO();
	LogisticDistribution d(0.0, 1.0);
	auto s = draw(d, 111, N);
	double m, var; meanVar(s, m, var);
	REQUIRE(m == Approx(0.0).margin(0.05));
	REQUIRE(var == Approx(Constants::PI * Constants::PI / 3.0).epsilon(0.06));
}

TEST_CASE("Sampling::Weibull_and_Pareto_and_LogNormal", "[sampling][weibull][pareto][lognormal]") {
	TEST_PRECISION_INFO();
	WeibullDistribution w(2.0, 1.0);
	auto sw = draw(w, 112, N);
	double mw, vw; meanVar(sw, mw, vw);
	REQUIRE(mw == Approx(std::tgamma(1.5)).epsilon(0.03));
	REQUIRE(*std::min_element(sw.begin(), sw.end()) > 0.0);

	ParetoDistribution p(3.0, 1.0);  // mean alpha*xm/(alpha-1) = 1.5
	auto sp = draw(p, 113, N);
	double mp, vp; meanVar(sp, mp, vp);
	REQUIRE(mp == Approx(1.5).epsilon(0.05));
	REQUIRE(*std::min_element(sp.begin(), sp.end()) >= 1.0);

	LogNormalDistribution ln(0.0, 1.0);  // mean exp(0.5)
	auto sl = draw(ln, 114, N);
	double ml, vl; meanVar(sl, ml, vl);
	REQUIRE(ml == Approx(std::exp(0.5)).epsilon(0.05));
	REQUIRE(*std::min_element(sl.begin(), sl.end()) > 0.0);
}

TEST_CASE("Sampling::Discrete_moments_and_support", "[sampling][discrete]") {
	TEST_PRECISION_INFO();
	const int M = 100000;

	{
		BernoulliDistribution d(0.3);
		auto s = draw(d, 201, M);
		double m, var; meanVar(s, m, var);
		REQUIRE(m == Approx(0.3).epsilon(0.03));
		for (double x : s) REQUIRE((x == 0.0 || x == 1.0));
	}
	{
		BinomialDistribution d(20, 0.3);
		auto s = draw(d, 202, M);
		double m, var; meanVar(s, m, var);
		REQUIRE(m == Approx(6.0).epsilon(0.02));
		REQUIRE(var == Approx(4.2).epsilon(0.06));
		REQUIRE(*std::min_element(s.begin(), s.end()) >= 0.0);
		REQUIRE(*std::max_element(s.begin(), s.end()) <= 20.0);
	}
	{
		PoissonDistribution d(4.0);
		auto s = draw(d, 203, M);
		double m, var; meanVar(s, m, var);
		REQUIRE(m == Approx(4.0).epsilon(0.03));
		REQUIRE(var == Approx(4.0).epsilon(0.06));
		REQUIRE(*std::min_element(s.begin(), s.end()) >= 0.0);
	}
	{
		GeometricDistribution d(0.25);  // mean 3, var 12
		auto s = draw(d, 204, M);
		double m, var; meanVar(s, m, var);
		REQUIRE(m == Approx(3.0).epsilon(0.04));
		REQUIRE(*std::min_element(s.begin(), s.end()) >= 0.0);
	}
	{
		NegativeBinomialDistribution d(5, 0.5);  // mean 5, var 10
		auto s = draw(d, 205, M);
		double m, var; meanVar(s, m, var);
		REQUIRE(m == Approx(5.0).epsilon(0.04));
		REQUIRE(*std::min_element(s.begin(), s.end()) >= 0.0);
	}
	{
		HypergeometricDistribution d(50, 15, 10);  // mean n*K/N = 3
		auto s = draw(d, 206, M);
		double m, var; meanVar(s, m, var);
		REQUIRE(m == Approx(3.0).epsilon(0.03));
		REQUIRE(*std::min_element(s.begin(), s.end()) >= 0.0);
		REQUIRE(*std::max_element(s.begin(), s.end()) <= 10.0);
	}
}

TEST_CASE("Sampling::default_engine_and_custom_engine", "[sampling][engine]") {
	TEST_PRECISION_INFO();
	// Default no-arg sample() uses MML's thread-local engine; just check it is in range.
	Random::SetSeed(7);
	UniformDistribution d(0.0, 1.0);
	for (int i = 0; i < 100; ++i) {
		Real x = d.sample();
		REQUIRE(x >= 0.0);
		REQUIRE(x <= 1.0);
	}
	// A different engine type also works via the template (expanded RNG engines).
	std::minstd_rand lcg(123);
	NormalDistribution nd(0.0, 1.0);
	Real y = nd.sample(lcg);
	REQUIRE(std::isfinite(y));
}

} // namespace MML::Tests::Algorithms::Statistics::SamplingTests

#include <catch2/catch_all.hpp>
#include "../../TestPrecision.h"

#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/algorithms/Statistics.h>
#endif

#include <limits>

using namespace MML;
using Catch::Matchers::WithinRel;

namespace {
	void CheckRelative(Real actual, Real expected) {
		CHECK(std::isfinite(actual));
		CHECK_THAT(actual, WithinRel(expected, Precision::EqualityTolerance));
	}

	Vector<Real> Transform(const Vector<Real>& data, Real scale, Real shift = 0) {
		Vector<Real> result(data.size());
		for (int i = 0; i < data.size(); ++i)
			result[i] = scale * data[i] + shift;
		return result;
	}
}

TEST_CASE("Statistics summaries avoid overflowing finite means and midpoints", "[statistics][statistics-scaling]") {
	const Real huge = std::numeric_limits<Real>::max() * Real(0.75);
	for (Real value : {huge, -huge, std::numeric_limits<Real>::denorm_min()}) {
		Vector<Real> data({value, value});
		CheckRelative(Statistics::Avg(data), value);
		CheckRelative(Statistics::Mean(data), value);
		CheckRelative(Statistics::Median(data), value);
		CheckRelative(Statistics::TrimmedMean(data, 0), value);
		Real avg, variance, stddev;
		Statistics::AvgVar(data, avg, variance);
		CheckRelative(avg, value);
		CHECK(variance == 0);
		Statistics::AvgStdDev(data, avg, stddev);
		CheckRelative(avg, value);
		CHECK(stddev == 0);
		CHECK(Statistics::PopulationVariance(data) == 0);
	}
}

TEST_CASE("Statistics recorded double-range summary witnesses", "[statistics][statistics-scaling]") {
	const Real value = std::pow(Real(10), MML::Real(std::numeric_limits<Real>::max_exponent10));
	Vector<Real> data({value, value});
	CheckRelative(Statistics::Avg(data), value);
	CheckRelative(Statistics::Median(data), value);
}

TEST_CASE("Statistics means retain cancellation residuals across extreme ranges", "[statistics][statistics-scaling]") {
	const Real huge = std::numeric_limits<Real>::max() * Real(0.75);
	const Real small = std::pow(Real(10), Real(-20));
	for (const auto& data : {Vector<Real>({huge, small, -huge}),
	                         Vector<Real>({huge, huge, small, -huge, -huge}),
	                         Vector<Real>({-huge, -huge, small, huge, huge})}) {
		CheckRelative(Statistics::Avg(data), small / data.size());
		CheckRelative(Statistics::TrimmedMean(data, 0), small / data.size());
	}
	CHECK(Statistics::Median(Vector<Real>({-huge, huge})) == 0);
}

TEST_CASE("Statistics sample and population summaries avoid overflowing square sums", "[statistics][statistics-scaling]") {
	const Real unit = std::sqrt(std::numeric_limits<Real>::max()) * Real(0.875);
	const Real squared = unit * unit;
	Vector<Real> data({-unit, 0, unit});
	Real avg, variance;
	Statistics::AvgVar(data, avg, variance);
	CHECK(avg == 0);
	CheckRelative(variance, squared);
	CheckRelative(Statistics::SampleVariance(data), squared);
	CheckRelative(Statistics::Variance(data), squared);
	CheckRelative(Statistics::PopulationVariance(data), (squared / 3) * 2);
	CheckRelative(Statistics::SampleStdDev(data), unit);
	CheckRelative(Statistics::PopulationStdDev(data), unit * std::sqrt(Real(2) / 3));
	CheckRelative(Statistics::Covariance(data, data), squared);
	CheckRelative(Statistics::Covariance(data, Transform(data, -1)), -squared);
}

TEST_CASE("Statistics standard deviations do not materialize unrepresentable variances", "[statistics][statistics-scaling]") {
	const int exponent = std::min(std::numeric_limits<Real>::max_exponent - 4,
		(std::numeric_limits<Real>::digits - std::numeric_limits<Real>::min_exponent) / 2 + 8);
	const Real large = std::scalbn(Real(1), exponent);
	for (Real unit : {large, Real(1) / large}) {
		Vector<Real> data({-unit, 0, unit});
		CheckRelative(Statistics::SampleStdDev(data), unit);
		CheckRelative(Statistics::StdDev(data), unit);
		CheckRelative(Statistics::PopulationStdDev(data), unit * std::sqrt(Real(2) / 3));
		Real avg, stddev;
		Statistics::AvgStdDev(data, avg, stddev);
		CHECK(avg == 0);
		CheckRelative(stddev, unit);
	}
}

TEST_CASE("Statistics centered summaries are translation invariant", "[statistics][statistics-scaling]") {
	const Real shift = std::scalbn(Real(1), std::numeric_limits<Real>::digits - 4);
	Vector<Real> data({-2, -1, 0, 1, 3});
	Vector<Real> y({3, 1, 0, 2, -1});
	const auto translated = Transform(data, 1, shift);
	const auto translatedY = Transform(y, 1, -shift);
	CheckRelative(Statistics::SampleVariance(translated), Statistics::SampleVariance(data));
	CheckRelative(Statistics::PopulationVariance(translated), Statistics::PopulationVariance(data));
	CheckRelative(Statistics::Covariance(translated, translatedY), Statistics::Covariance(data, y));
	CheckRelative(Statistics::PearsonCorrelation(translated, translatedY), Statistics::PearsonCorrelation(data, y));
}

TEST_CASE("Statistics Pearson correlation survives independent unit scales", "[statistics][statistics-scaling]") {
	const Real large = std::pow(Real(10), MML::Real(std::min(100, std::numeric_limits<Real>::max_exponent10 / 3)));
	Vector<Real> x({-2, -1, 0, 1, 3});
	Vector<Real> y({3, 1, 0, 2, -1});
	const Real expected = Statistics::PearsonCorrelation(x, y);
	for (Real xScale : {large, Real(1), Real(1) / large}) {
		for (Real yScale : {large, Real(1), Real(1) / large}) {
			CheckRelative(Statistics::PearsonCorrelation(Transform(x, xScale), Transform(y, yScale)), expected);
			CheckRelative(Statistics::PearsonCorrelation(Transform(x, xScale), Transform(y, -yScale)), -expected);
		}
		Vector<Real> witness({-xScale, 0, xScale});
		CheckRelative(Statistics::PearsonCorrelation(witness, witness), 1);
		CheckRelative(Statistics::PearsonCorrelation(witness, Transform(witness, -1)), -1);
	}
}

TEST_CASE("Statistics Pearson handles maximal and subnormal nonconstant data", "[statistics][statistics-scaling]") {
	for (Real unit : {std::numeric_limits<Real>::max(), std::numeric_limits<Real>::denorm_min()}) {
		Vector<Real> data({-unit, 0, unit});
		CheckRelative(Statistics::PearsonCorrelation(data, data), 1);
		CheckRelative(Statistics::RSquared(data, data), 1);
		const auto result = Statistics::PearsonCorrelationWithTest(data, Transform(data, -1));
		CheckRelative(result.r, -1);
		CHECK(result.degreesOfFreedom == 1);
		CHECK(result.tStatistic < 0);
	}
}

TEST_CASE("Statistics covariance restores independent binary units safely", "[statistics][statistics-scaling]") {
	const Real large = std::scalbn(Real(1), std::numeric_limits<Real>::max_exponent - 4);
	Vector<Real> base({-2, -1, 0, 1, 3});
	Vector<Real> other({3, 1, 0, 2, -1});
	CheckRelative(Statistics::Covariance(Transform(base, large), Transform(other, Real(1) / large)),
	              Statistics::Covariance(base, other));
}

TEST_CASE("Statistics moments keep shape statistics invariant under unit changes", "[statistics][statistics-scaling]") {
	Vector<Real> base({-2, -1, 0, 1, 4});
	Real mean, deviation, stddev, variance, skew, kurtosis;
	Statistics::Moments(base, mean, deviation, stddev, variance, skew, kurtosis);
	const Real large = std::pow(Real(10), MML::Real(std::min(110, std::numeric_limits<Real>::max_exponent10 / 3 + 3)));
	for (Real unit : {large, Real(1) / large, -large}) {
		Real m, d, s, v, sk, k;
		Statistics::Moments(Transform(base, unit), m, d, s, v, sk, k);
		CheckRelative(m, mean * unit);
		CheckRelative(d, deviation * std::abs(unit));
		CheckRelative(s, stddev * std::abs(unit));
		CheckRelative(v, (variance * unit) * unit);
		CheckRelative(sk, unit > 0 ? skew : -skew);
		CheckRelative(k, kurtosis);
	}
}

TEST_CASE("Statistics moments distinguish underflowed variance from constant data", "[statistics][statistics-scaling]") {
	const Real unit = std::scalbn(Real(1),
		(std::numeric_limits<Real>::min_exponent - std::numeric_limits<Real>::digits) / 2 - 8);
	Real mean, deviation, stddev, variance, skew, kurtosis;
	Statistics::Moments(Vector<Real>({-unit, 0, unit}), mean, deviation, stddev, variance, skew, kurtosis);
	CHECK(mean == 0);
	CHECK(variance == 0);
	CheckRelative(deviation, unit * (Real(2) / 3));
	CheckRelative(stddev, unit);
	CHECK(skew == 0);
	CheckRelative(kurtosis, -Real(7) / 3);
}

TEST_CASE("Statistics matrix wrappers share stable scalar arithmetic", "[statistics][statistics-scaling]") {
	const Real large = std::sqrt(std::numeric_limits<Real>::max()) * Real(0.875);
	Matrix<Real> data(3, 2, {-large, large, 0, 0, large, -large});
	const auto cov = Statistics::CovarianceMatrix(data);
	CheckRelative(cov(0, 0), large * large);
	CheckRelative(cov(0, 1), -(large * large));
	CHECK(cov(0, 1) == cov(1, 0));

	for (Real unit : {std::numeric_limits<Real>::max() / 2, std::numeric_limits<Real>::denorm_min()}) {
		Matrix<Real> observations(3, 2, {-unit, unit, 0, 0, unit, -unit});
		const auto correlation = Statistics::CorrelationMatrix(observations);
		CHECK(correlation(0, 0) == 1);
		CHECK(correlation(1, 1) == 1);
		CheckRelative(correlation(0, 1), -1);
		CHECK(correlation(0, 1) == correlation(1, 0));
	}
}

TEST_CASE("Statistics median deviation tolerates discarded overflowing deviations", "[statistics][statistics-scaling]") {
	const Real huge = std::numeric_limits<Real>::max();
	CHECK(Statistics::MAD(Vector<Real>({-huge, -huge, huge})) == 0);
	const Real small = std::pow(Real(10), Real(-20));
	CheckRelative(Statistics::MAD(Vector<Real>({-huge, 0, small, 2 * small, huge})), small);
}

TEST_CASE("Statistics near-maximal translated data keep resolvable variation", "[statistics][statistics-scaling]") {
	const Real origin = std::scalbn(Real(1), std::numeric_limits<Real>::max_exponent - 2);
	const Real step = std::scalbn(origin, 2 - std::numeric_limits<Real>::digits);
	Vector<Real> data({origin - step, origin, origin + step});
	CheckRelative(Statistics::SampleStdDev(data), step);
	CheckRelative(Statistics::PearsonCorrelation(data, data), 1);

	const Real maximum = std::numeric_limits<Real>::max();
	Vector<Real> constant(13);
	for (int i = 0; i < constant.size(); ++i)
		constant[i] = maximum;
	CheckRelative(Statistics::Avg(constant), maximum);
	CHECK(Statistics::SampleVariance(constant) == 0);
	Matrix<Real> observations(3, 1, {maximum, maximum, maximum});
	CHECK(Statistics::CovarianceMatrix(observations)(0, 0) == 0);
}

TEST_CASE("Statistics matrices preserve correlations between disparate units", "[statistics][statistics-scaling]") {
	const Real large = std::scalbn(Real(1), std::numeric_limits<Real>::max_exponent - 4);
	const Real small = Real(1) / large;
	Matrix<Real> data(3, 2, {-large, small, 0, 0, large, -small});
	const auto result = Statistics::CorrelationMatrix(data);
	CheckRelative(result(0, 1), -1);
	CHECK(result(0, 1) == result(1, 0));
	CHECK_THROWS_AS(Statistics::TrimmedMean(Vector<Real>({1, 2, 3}), std::numeric_limits<Real>::quiet_NaN()),
	                StatisticsError);
}

TEST_CASE("Statistics summaries reject nonfinite input consistently", "[statistics][statistics-scaling]") {
	for (Real invalid : {std::numeric_limits<Real>::infinity(), -std::numeric_limits<Real>::infinity(),
	                     std::numeric_limits<Real>::quiet_NaN()}) {
		Vector<Real> bad({1, invalid, 3});
		Vector<Real> good({1, 2, 3});
		Real a = 7, b = 8, c, d, e, f;
		CHECK_THROWS_AS(Statistics::Avg(bad), StatisticsError);
		CHECK_THROWS_AS(Statistics::AvgVar(bad, a, b), StatisticsError);
		CHECK_THROWS_AS(Statistics::AvgStdDev(bad, a, b), StatisticsError);
		CHECK_THROWS_AS(Statistics::PopulationVariance(bad), StatisticsError);
		CHECK_THROWS_AS(Statistics::SampleStdDev(bad), StatisticsError);
		CHECK_THROWS_AS(Statistics::PopulationStdDev(bad), StatisticsError);
		CHECK_THROWS_AS(Statistics::Moments(bad, a, b, c, d, e, f), StatisticsError);
		CHECK_THROWS_AS(Statistics::Median(bad), StatisticsError);
		CHECK_THROWS_AS(Statistics::TrimmedMean(bad, 0), StatisticsError);
		CHECK_THROWS_AS(Statistics::MAD(bad), StatisticsError);
		CHECK_THROWS_AS(Statistics::Covariance(bad, good), StatisticsError);
		CHECK_THROWS_AS(Statistics::Covariance(good, bad), StatisticsError);
		CHECK_THROWS_AS(Statistics::PearsonCorrelation(bad, good), StatisticsError);
		CHECK_THROWS_AS(Statistics::PearsonCorrelation(good, bad), StatisticsError);
	}
}

TEST_CASE("Statistics genuine output overflow is explicit", "[statistics][statistics-scaling]") {
	const Real unit = std::numeric_limits<Real>::max();
	Vector<Real> data({-unit, 0, unit});
	Real a, b, c, d, e, f;
	CHECK_THROWS_AS(Statistics::SampleVariance(data), StatisticsError);
	CHECK_THROWS_AS(Statistics::PopulationVariance(data), StatisticsError);
	CHECK_THROWS_AS(Statistics::AvgVar(data, a, b), StatisticsError);
	CHECK_THROWS_AS(Statistics::Moments(data, a, b, c, d, e, f), StatisticsError);
	CHECK_THROWS_AS(Statistics::Covariance(data, data), StatisticsError);
	Matrix<Real> observations(3, 1, {-unit, 0, unit});
	CHECK_THROWS_AS(Statistics::CovarianceMatrix(observations), StatisticsError);
	CheckRelative(Statistics::SampleStdDev(data), unit);
	CheckRelative(Statistics::PopulationStdDev(data), unit * std::sqrt(Real(2) / 3));
}

TEST_CASE("Statistics zero variance and shape contracts remain explicit", "[statistics][statistics-scaling]") {
	Vector<Real> constant({2, 2, 2});
	Vector<Real> nonconstant({1, 2, 3});
	CHECK_THROWS_AS(Statistics::PearsonCorrelation(constant, nonconstant), StatisticsError);
	CHECK_THROWS_AS(Statistics::PearsonCorrelation(nonconstant, constant), StatisticsError);
	CHECK_THROWS_AS(Statistics::PearsonCorrelation(nonconstant, Vector<Real>({1, 2})), StatisticsError);
	CHECK(Statistics::PopulationVariance(Vector<Real>({2})) == 0);
	CHECK(Statistics::PopulationStdDev(Vector<Real>({2})) == 0);
	CHECK_THROWS_AS(Statistics::SampleVariance(Vector<Real>({2})), StatisticsError);
	CHECK_THROWS_AS(Statistics::SampleStdDev(Vector<Real>({2})), StatisticsError);
	CHECK_THROWS_AS(Statistics::CorrelationMatrix(Matrix<Real>(3, 1, {2, 2, 2})), StatisticsError);
}

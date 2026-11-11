#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "../../TestPrecision.h"

#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/core/DifferentialGeometry/Metric.h>
#include <mml/base/DifferentialGeometry/Hodge.h>
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <type_traits>

using namespace MML;

namespace {
	struct HodgeFrame4 {};
	Real LargeMetricScale() { return std::ldexp(Real(1), 2 * std::numeric_limits<Real>::max_exponent / 3); }
	void CheckHodgeValue(Real actual, Real expected) {
		CHECK(std::isfinite(actual));
		if (expected == 0) CHECK(actual == 0);
		else CHECK_THAT(actual / expected, Catch::Matchers::WithinAbs(1, TOL(2e-12, 2e-4)));
	}
	auto TopForm4(Real value) {
		DifferentialForm<Real, 4, 4, HodgeFrame4> form;
		std::array<int, 4> indices{0, 1, 2, 3};
		do {
			int inversions = 0;
			for (int i = 0; i < 4; ++i)
				for (int j = i + 1; j < 4; ++j)
					inversions += indices[i] > indices[j];
			form.ComponentAt(indices) = inversions % 2 == 0 ? value : -value;
		} while (std::next_permutation(indices.begin(), indices.end()));
		return form;
	}
}

TEST_CASE("Hodge scaling preserves the recorded large conformal metric witness", "[Hodge][hodge-scaling]")
{
	const Real large = LargeMetricScale();
	const Real recorded = std::pow(Real(10), Real(std::min(200, std::numeric_limits<Real>::max_exponent / 4)));
	for (Real scale : {Real(1), large, 1 / large, recorded, 1 / recorded}) {
		INFO("metric scale " << scale);
		auto metric = Metric<2, Cartesian2>::Diagonal({scale, scale});
		auto dx = BasisOneForm<0, 2, Cartesian2>();
		auto result = HodgeStar(dx, metric);
		CheckHodgeValue(result.Component(0), 0);
		CheckHodgeValue(result.Component(1), 1);
		CheckHodgeValue(hodge_star(dx, metric, Orientation::Negative).Component(1), -1);
		CheckHodgeValue(HodgeStar(result, metric).Component(0), -1);
	}
}

TEST_CASE("Hodge scaling combines unrepresentable four dimensional middle degree factors", "[Hodge][hodge-scaling]")
{
	for (Real scale : {LargeMetricScale(), 1 / LargeMetricScale()}) {
		auto metric = Metric<4, HodgeFrame4>::Diagonal({scale, scale, scale, scale});
		DifferentialForm<Real, 4, 2, HodgeFrame4> form;
		form.Component(0, 1) = 1;
		form.Component(1, 0) = -1;
		const auto result = HodgeStar(form, metric);
		for (int i = 0; i < 4; ++i)
			for (int j = 0; j < 4; ++j)
				CheckHodgeValue(result.Component(i, j), Real(i == 2 && j == 3 ? 1 : i == 3 && j == 2 ? -1 : 0));
		const auto twice = HodgeStar(result, metric);
		CheckHodgeValue(twice.Component(0, 1), 1);
		CheckHodgeValue(twice.Component(1, 0), -1);
	}
}

TEST_CASE("Hodge scaling includes form coefficients before materializing metric factors", "[Hodge][hodge-scaling]")
{
	const Real scale = LargeMetricScale();
	auto metric = Metric<4, HodgeFrame4>::Diagonal({scale, scale, scale, scale});
	DifferentialForm<Real, 4, 0, HodgeFrame4> scalar;
	scalar.Component() = 1 / scale;
	const auto volume = HodgeStar(scalar, metric);
	CheckHodgeValue(volume.Component(0, 1, 2, 3), scale);
	CheckHodgeValue(volume.Component(1, 0, 2, 3), -scale);
	CheckHodgeValue(HodgeStar(TopForm4(scale), metric).Component(), 1 / scale);
	const auto threeForm = Wedge(Wedge(BasisOneForm<0, 4, HodgeFrame4>(),
		BasisOneForm<1, 4, HodgeFrame4>()), BasisOneForm<2, 4, HodgeFrame4>());
	for (Real factor : {scale, 1 / scale}) {
		const auto degreeMetric = Metric<4, HodgeFrame4>::Diagonal({factor, factor, factor, factor});
		const auto oneForm = HodgeStar(threeForm, degreeMetric);
		CheckHodgeValue(oneForm.Component(3), 1 / factor);
		CheckHodgeValue(HodgeStar(oneForm, degreeMetric).Component(0, 1, 2), -1);
	}
}

TEST_CASE("Hodge scaling handles anisotropic and indefinite diagonal metrics", "[Hodge][hodge-scaling]")
{
	const Real large = LargeMetricScale();
	for (Real sign : {Real(-1), Real(1)}) {
		auto metric = Metric<2, Cartesian2>::Diagonal({sign * large, 1 / large});
		const auto dx = BasisOneForm<0, 2, Cartesian2>();
		const auto dy = BasisOneForm<1, 2, Cartesian2>();
		CheckHodgeValue(HodgeStar(dx, metric).Component(1), sign / large);
		CheckHodgeValue(HodgeStar(dy, metric).Component(0), -large);
		CheckHodgeValue(HodgeStar(HodgeStar(dx, metric), metric).Component(0), -sign);
	}
}

TEST_CASE("Hodge scaling preserves non diagonal metric mixing", "[Hodge][hodge-scaling]")
{
	for (Real scale : {Real(1), LargeMetricScale(), 1 / LargeMetricScale()}) {
		MatrixNM<Real, 2, 2> components;
		components(0, 0) = components(1, 1) = 2 * scale;
		components(0, 1) = components(1, 0) = scale;
		Metric<2, Cartesian2> metric(components);
		const auto dx = BasisOneForm<0, 2, Cartesian2>();
		const auto dy = BasisOneForm<1, 2, Cartesian2>();
		const auto starDx = HodgeStar(dx, metric), starDy = HodgeStar(dy, metric);
		const Real root = std::sqrt(Real(3));
		CheckHodgeValue(starDx.Component(0), 1 / root);
		CheckHodgeValue(starDx.Component(1), 2 / root);
		CheckHodgeValue(starDy.Component(0), -2 / root);
		CheckHodgeValue(starDy.Component(1), -1 / root);
		CheckHodgeValue(HodgeStar(starDx, metric).Component(0), -1);
	}
}

TEST_CASE("Hodge scaling protects LU updates and pivot interchanges", "[Hodge][hodge-scaling]")
{
	const Real large = std::numeric_limits<Real>::max() / 2;
	MatrixNM<Real, 2, 2> components;
	components(0, 0) = components(0, 1) = components(1, 0) = large + large / 2;
	components(1, 1) = -components(0, 0);
	Metric<2, Cartesian2> extreme(components);
	const auto result = HodgeStar(BasisOneForm<0, 2, Cartesian2>(), extreme);
	CheckHodgeValue(result.Component(0), -1 / std::sqrt(Real(2)));
	CheckHodgeValue(result.Component(1), 1 / std::sqrt(Real(2)));

	components(0, 0) = components(1, 1) = 1;
	components(0, 1) = components(1, 0) = 2;
	Metric<2, Cartesian2> pivoted(components);
	const auto swapped = HodgeStar(BasisOneForm<0, 2, Cartesian2>(), pivoted);
	CheckHodgeValue(swapped.Component(0), -2 / std::sqrt(Real(3)));
	CheckHodgeValue(swapped.Component(1), -1 / std::sqrt(Real(3)));
}

TEST_CASE("Hodge scaling preserves real and complex scalar types and separated component ranges", "[Hodge][hodge-scaling]")
{
	auto metric = Metric<2, Cartesian2>::Diagonal({LargeMetricScale(), LargeMetricScale()});
	DifferentialForm<std::complex<Real>, 2, 1, Cartesian2> complexForm;
	complexForm.Component(0) = {LargeMetricScale(), 1 / LargeMetricScale()};
	const auto complexResult = HodgeStar(complexForm, metric);
	static_assert(std::is_same_v<typename decltype(complexResult)::value_type, std::complex<Real>>);
	CheckHodgeValue(complexResult.Component(1).real(), LargeMetricScale());
	CheckHodgeValue(complexResult.Component(1).imag(), 1 / LargeMetricScale());
	DifferentialForm<float, 2, 1, Cartesian2> floatForm;
	floatForm.Component(0) = 1;
	const auto floatResult = HodgeStar(floatForm, metric);
	static_assert(std::is_same_v<typename decltype(floatResult)::value_type, float>);
	CHECK(floatResult.Component(1) == 1);
	DifferentialForm<int, 2, 1, Cartesian2> integerForm;
	integerForm.Component(0) = 3;
	const auto integerResult = HodgeStar(integerForm, Metric<2, Cartesian2>::Euclidean());
	static_assert(std::is_same_v<typename decltype(integerResult)::value_type, int>);
	CHECK(integerResult.Component(1) == 3);
}

TEST_CASE("Hodge scaling postpones factorial division until after a scaled contraction", "[Hodge][hodge-scaling]")
{
	const Real value = std::numeric_limits<Real>::max() / 2;
	auto metric = Metric<4, HodgeFrame4>::Euclidean();
	CheckHodgeValue(HodgeStar(TopForm4(value), metric).Component(), value);
	DifferentialForm<Real, 2, 2, Cartesian2> raw;
	raw.Component(0, 1) = 2;
	raw.Component(1, 0) = 1;
	CheckHodgeValue(HodgeStar(raw, Metric<2, Cartesian2>::Euclidean()).Component(), Real(1) / 2);
}

TEST_CASE("Hodge scaling reports truly unrepresentable outputs but preserves exact zeros", "[Hodge][hodge-scaling]")
{
	for (Real scale : {LargeMetricScale(), 1 / LargeMetricScale()}) {
		auto metric = Metric<4, HodgeFrame4>::Diagonal({scale, scale, scale, scale});
		DifferentialForm<Real, 4, 0, HodgeFrame4> scalar;
		scalar.Component() = 1;
		CHECK_THROWS_AS(HodgeStar(scalar, metric), GeometryError);
		CHECK_THROWS_AS(HodgeStar(TopForm4(1), metric), GeometryError);
		scalar.Component() = 0;
		const auto zeroVolume = HodgeStar(scalar, metric);
		for (int i = 0; i < zeroVolume.ComponentCount; ++i)
			CHECK(zeroVolume[i] == 0);
		CHECK(HodgeStar(TopForm4(0), metric).Component() == 0);
	}
}

TEST_CASE("Hodge scaling rejects nonfinite form coefficients and invalid orientation", "[Hodge][hodge-scaling]")
{
	auto metric = Metric<2, Cartesian2>::Euclidean();
	for (Real invalid : {std::numeric_limits<Real>::infinity(), -std::numeric_limits<Real>::infinity(),
		std::numeric_limits<Real>::quiet_NaN()}) {
		DifferentialForm<Real, 2, 1, Cartesian2> form;
		form.Component(0) = invalid;
		CHECK_THROWS_AS(HodgeStar(form, metric), GeometryError);
	}
	DifferentialForm<std::complex<Real>, 2, 1, Cartesian2> complexForm;
	complexForm.Component(0) = {1, std::numeric_limits<Real>::quiet_NaN()};
	CHECK_THROWS_AS(HodgeStar(complexForm, metric), GeometryError);
	CHECK_THROWS_AS(HodgeStar(BasisOneForm<0, 2, Cartesian2>(), metric, static_cast<Orientation>(0)), GeometryError);
}

TEST_CASE("Hodge scaling accepts subnormal results and checks the actual scalar output range", "[Hodge][hodge-scaling]")
{
	DifferentialForm<Real, 2, 1, Cartesian2> tiny;
	tiny.Component(0) = std::numeric_limits<Real>::denorm_min();
	const auto result = HodgeStar(tiny, Metric<2, Cartesian2>::Euclidean());
	CHECK(result.Component(1) == tiny.Component(0));
	constexpr int extraExponentRange = std::numeric_limits<Real>::max_exponent - std::numeric_limits<float>::max_exponent;
	if constexpr (extraExponentRange > 0) {
		DifferentialForm<float, 2, 1, Cartesian2> floatForm;
		floatForm.Component(0) = 1;
		auto overflow = Metric<2, Cartesian2>::Diagonal({1 / LargeMetricScale(), LargeMetricScale()});
		auto underflow = Metric<2, Cartesian2>::Diagonal({LargeMetricScale(), 1 / LargeMetricScale()});
		CHECK_THROWS_AS(HodgeStar(floatForm, overflow), GeometryError);
		CHECK_THROWS_AS(HodgeStar(floatForm, underflow), GeometryError);
	}
}

TEST_CASE("Hodge scaling reaches the derived metric cross product", "[Hodge][hodge-scaling]")
{
	for (Real scale : {LargeMetricScale(), 1 / LargeMetricScale()}) {
		const Real length = 1 / std::sqrt(scale);
		auto metric = Metric<3, Cartesian3>::Diagonal({scale, scale, scale});
		TangentVector<3, Cartesian3> first{length, 0, 0}, second{0, length, 0};
		const auto result = Cross(first, second, metric);
		CheckHodgeValue(result[0], 0);
		CheckHodgeValue(result[1], 0);
		CheckHodgeValue(result[2], length);
	}
}

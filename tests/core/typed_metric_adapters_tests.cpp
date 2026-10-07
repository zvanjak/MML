#include <catch2/catch_all.hpp>

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/base/DifferentialGeometry/Metric.h>
#include <mml/core/MetricTensor.h>
#include <mml/core/DifferentialGeometry/MetricAdapters.h>
#endif

using namespace MML;

namespace MML::Tests::Core::TypedMetricAdaptersTests
{
	TEST_CASE("TypedMetricAdapters - Cartesian MetricTensorField samples to typed Metric", "[Metric][MetricTensor][DifferentialForms]")
	{
		MetricTensorCartesian3D metricField;
		VectorN<Real, 3> pos{ REAL(1.0), REAL(2.0), REAL(3.0) };

		Metric<3, Cartesian3> metric = MetricAt<3, Cartesian3>(metricField, pos);

		REQUIRE(metric.covariantComponents()(0, 0) == REAL(1.0));
		REQUIRE(metric.covariantComponents()(1, 1) == REAL(1.0));
		REQUIRE(metric.covariantComponents()(2, 2) == REAL(1.0));
		REQUIRE(metric.covariantComponents()(0, 1) == REAL(0.0));
		REQUIRE(metric.contravariantComponents()(2, 2) == REAL(1.0));
	}

	TEST_CASE("TypedMetricAdapters - Spherical MetricTensorField lowers typed tangent vectors", "[Metric][MetricTensor][DifferentialForms]")
	{
		MetricTensorSpherical metricField;
		VectorN<Real, 3> pos{ REAL(2.0), Constants::PI / REAL(2.0), REAL(0.0) };

		Metric<3, Spherical3> metric = MetricAt<3, Spherical3>(metricField, pos);
		TangentVector<3, Spherical3> v{ REAL(1.0), REAL(2.0), REAL(3.0) };

		Covector<3, Spherical3> alpha = Flat(metric, v);

		REQUIRE(metric.covariantComponents()(0, 0) == REAL(1.0));
		REQUIRE(metric.covariantComponents()(1, 1) == REAL(4.0));
		REQUIRE(metric.covariantComponents()(2, 2) == Catch::Approx(REAL(4.0)));
		REQUIRE(alpha[0] == REAL(1.0));
		REQUIRE(alpha[1] == REAL(8.0));
		REQUIRE(alpha[2] == Catch::Approx(REAL(12.0)));
	}
}
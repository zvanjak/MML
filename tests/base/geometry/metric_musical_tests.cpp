#include <catch2/catch_all.hpp>

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/base/DifferentialGeometry/Metric.h>
#endif

#include <type_traits>
#include <utility>

using namespace MML;

namespace MML::Tests::Base::MetricMusicalTests
{
	template<class, class, class = void>
	struct HasFlat : std::false_type { };

	template<class MetricType, class VectorType>
	struct HasFlat<MetricType, VectorType, std::void_t<decltype(Flat(std::declval<MetricType>(), std::declval<VectorType>()))>> : std::true_type { };

	template<class, class, class = void>
	struct HasSharp : std::false_type { };

	template<class MetricType, class VectorType>
	struct HasSharp<MetricType, VectorType, std::void_t<decltype(Sharp(std::declval<MetricType>(), std::declval<VectorType>()))>> : std::true_type { };

	TEST_CASE("Metric - Euclidean construction and component access", "[Metric][DifferentialForms]")
	{
		using G = Metric<3, Cartesian3>;

		G metric = G::Euclidean();
		static_assert(G::Dimension == 3, "Metric dimension metadata is fixed at compile time");
		static_assert(std::is_same<typename G::frame_type, Cartesian3>::value, "Frame tag is part of the metric type");

		REQUIRE(metric(0, 0) == REAL(1.0));
		REQUIRE(metric(1, 1) == REAL(1.0));
		REQUIRE(metric(2, 2) == REAL(1.0));
		REQUIRE(metric(0, 1) == REAL(0.0));
		REQUIRE(metric.contravariantComponents()(0, 0) == REAL(1.0));
		REQUIRE(metric.contravariantComponents()(1, 2) == REAL(0.0));
	}

	TEST_CASE("Metric - Diagonal construction computes inverse components", "[Metric][DifferentialForms]")
	{
		Metric<3, Cartesian3> metric = Metric<3, Cartesian3>::Diagonal({ REAL(2.0), REAL(4.0), REAL(5.0) });

		REQUIRE(metric.covariantComponents()(0, 0) == REAL(2.0));
		REQUIRE(metric.covariantComponents()(1, 1) == REAL(4.0));
		REQUIRE(metric.covariantComponents()(2, 2) == REAL(5.0));
		REQUIRE(metric.contravariantComponents()(0, 0) == Catch::Approx(REAL(0.5)));
		REQUIRE(metric.contravariantComponents()(1, 1) == Catch::Approx(REAL(0.25)));
		REQUIRE(metric.contravariantComponents()(2, 2) == Catch::Approx(REAL(0.2)));
	}

	TEST_CASE("Metric - Flat sharp and inner products use explicit metric", "[Metric][DifferentialForms]")
	{
		Metric<3, Cartesian3> metric = Metric<3, Cartesian3>::Diagonal({ REAL(2.0), REAL(4.0), REAL(5.0) });
		TangentVector<3, Cartesian3> a{ REAL(1.0), REAL(2.0), REAL(3.0) };
		TangentVector<3, Cartesian3> b{ -REAL(2.0), REAL(0.5), REAL(4.0) };

		Covector<3, Cartesian3> aFlat = Flat(metric, a);
		REQUIRE(aFlat[0] == REAL(2.0));
		REQUIRE(aFlat[1] == REAL(8.0));
		REQUIRE(aFlat[2] == REAL(15.0));

		TangentVector<3, Cartesian3> roundTrip = Sharp(metric, aFlat);
		REQUIRE(roundTrip[0] == Catch::Approx(a[0]));
		REQUIRE(roundTrip[1] == Catch::Approx(a[1]));
		REQUIRE(roundTrip[2] == Catch::Approx(a[2]));

		Real expectedInner = REAL(2.0) * REAL(1.0) * -REAL(2.0)
			+ REAL(4.0) * REAL(2.0) * REAL(0.5)
			+ REAL(5.0) * REAL(3.0) * REAL(4.0);
		REQUIRE(Inner(metric, a, b) == expectedInner);
		REQUIRE(Inner(metric, a, b) == Pair(aFlat, b));
		REQUIRE(flat(metric, a).components() == aFlat.components());
		REQUIRE(sharp(metric, aFlat).components() == roundTrip.components());
		REQUIRE(inner(metric, a, b) == Inner(metric, a, b));
	}

	TEST_CASE("Metric - Musical operations reject invalid categories and frames", "[Metric][DifferentialForms]")
	{
		using CartesianMetric = Metric<3, Cartesian3>;
		using CartesianVector = TangentVector<3, Cartesian3>;
		using CartesianCovector = Covector<3, Cartesian3>;
		using SphericalVector = TangentVector<3, Spherical3>;

		static_assert(HasFlat<CartesianMetric, CartesianVector>::value, "flat lowers same-frame tangent vectors");
		static_assert(!HasFlat<CartesianMetric, CartesianCovector>::value, "flat does not accept covectors");
		static_assert(!HasFlat<CartesianMetric, SphericalVector>::value, "flat rejects different frames");
		static_assert(HasSharp<CartesianMetric, CartesianCovector>::value, "sharp raises same-frame covectors");
		static_assert(!HasSharp<CartesianMetric, CartesianVector>::value, "sharp does not accept tangent vectors");

		SUCCEED("Invalid musical operations are rejected by the overload set.");
	}

	TEST_CASE("Metric - Singular covariant components throw during construction", "[Metric][DifferentialForms]")
	{
		MatrixNM<Real, 3, 3> singular;
		singular(0, 0) = REAL(1.0);
		singular(1, 1) = REAL(1.0);

		REQUIRE_THROWS_AS((Metric<3, Cartesian3>(singular)), SingularMatrixError);
	}
}
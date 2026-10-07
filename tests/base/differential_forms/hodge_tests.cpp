#include <catch2/catch_all.hpp>

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/base/DifferentialGeometry/Hodge.h>
#endif

#include <type_traits>
#include <utility>

using namespace MML;

namespace MML::Tests::Base::HodgeTests
{
	template<class, class, class, class = void>
	struct HasCross : std::false_type { };

	template<class A, class B, class G>
	struct HasCross<A, B, G, std::void_t<decltype(Cross(std::declval<A>(), std::declval<B>(), std::declval<G>()))>> : std::true_type { };

	TEST_CASE("HodgeStar - Euclidean 2D basis one-form signs", "[Hodge][DifferentialForms]")
	{
		Metric<2, Cartesian2> metric = Metric<2, Cartesian2>::Euclidean();
		auto dx = BasisOneForm<0, 2, Cartesian2>();
		auto dy = BasisOneForm<1, 2, Cartesian2>();

		auto starDx = HodgeStar(dx, metric);
		auto starDy = HodgeStar(dy, metric);
		auto negativeStarDx = HodgeStar(dx, metric, Orientation::Negative);

		REQUIRE(starDx.Component(0) == REAL(0.0));
		REQUIRE(starDx.Component(1) == REAL(1.0));
		REQUIRE(starDy.Component(0) == -REAL(1.0));
		REQUIRE(starDy.Component(1) == REAL(0.0));
		REQUIRE(negativeStarDx.Component(1) == -REAL(1.0));
	}

	TEST_CASE("HodgeStar - Double star identities in Euclidean dimensions", "[Hodge][DifferentialForms]")
	{
		Metric<2, Cartesian2> metric2 = Metric<2, Cartesian2>::Euclidean();
		auto dx2 = BasisOneForm<0, 2, Cartesian2>();
		auto doubleDx2 = HodgeStar(HodgeStar(dx2, metric2), metric2);
		REQUIRE(doubleDx2.Component(0) == -REAL(1.0));
		REQUIRE(doubleDx2.Component(1) == REAL(0.0));

		Metric<3, Cartesian3> metric3 = Metric<3, Cartesian3>::Euclidean();
		auto dx3 = BasisOneForm<0, 3, Cartesian3>();
		auto doubleDx3 = HodgeStar(HodgeStar(dx3, metric3), metric3);
		REQUIRE(doubleDx3.Component(0) == REAL(1.0));
		REQUIRE(doubleDx3.Component(1) == REAL(0.0));
		REQUIRE(doubleDx3.Component(2) == REAL(0.0));
	}

	TEST_CASE("HodgeStar - Diagonal metric changes one-form scale", "[Hodge][DifferentialForms]")
	{
		Metric<2, Cartesian2> metric = Metric<2, Cartesian2>::Diagonal({ REAL(4.0), REAL(9.0) });
		auto dx = BasisOneForm<0, 2, Cartesian2>();

		auto starDx = HodgeStar(dx, metric);

		REQUIRE(starDx.Component(0) == REAL(0.0));
		REQUIRE(starDx.Component(1) == Catch::Approx(REAL(1.5)));
	}

	TEST_CASE("Cross - Derived typed 3D cross product matches Cartesian formula", "[Hodge][CrossProduct][DifferentialForms]")
	{
		Metric<3, Cartesian3> metric = Metric<3, Cartesian3>::Euclidean();
		TangentVector<3, Cartesian3> a{ REAL(1.0), REAL(2.0), REAL(3.0) };
		TangentVector<3, Cartesian3> b{ REAL(4.0), -REAL(5.0), REAL(6.0) };

		TangentVector<3, Cartesian3> c = Cross(a, b, metric);
		TangentVector<3, Cartesian3> flipped = Cross(a, b, metric, Orientation::Negative);

		REQUIRE(c[0] == REAL(27.0));
		REQUIRE(c[1] == REAL(6.0));
		REQUIRE(c[2] == -REAL(13.0));
		REQUIRE(flipped[0] == -c[0]);
		REQUIRE(flipped[1] == -c[1]);
		REQUIRE(flipped[2] == -c[2]);
	}

	TEST_CASE("Cross - Only 3D typed vectors are supported", "[Hodge][CrossProduct][DifferentialForms]")
	{
		using Vector2 = TangentVector<2, Cartesian2>;
		using Metric2 = Metric<2, Cartesian2>;
		using Vector3 = TangentVector<3, Cartesian3>;
		using Metric3 = Metric<3, Cartesian3>;

		static_assert(!HasCross<Vector2, Vector2, Metric2>::value, "2D cross product is not available");
		static_assert(HasCross<Vector3, Vector3, Metric3>::value, "3D cross product is available");
		SUCCEED("Invalid cross product dimensions are rejected by the overload set.");
	}
}
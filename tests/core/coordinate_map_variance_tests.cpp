#include <catch2/catch_all.hpp>

#include "../TestPrecision.h"

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/base/DifferentialGeometry/DifferentialForm.h>
#include <mml/core/CoordTransf/CoordTransf2D.h>
#include <mml/core/CoordTransf/CoordTransfCylindrical.h>
#include <mml/core/DifferentialGeometry/CoordinateMap.h>
#endif

#include <type_traits>
#include <utility>

using namespace MML;

namespace MML::Tests::Core::CoordinateMapVarianceTests
{
	constexpr Real CoordinateTolerance = Testing::Tol(REAL(1e-12), REAL(1e-5));

	template<class, class, class = void>
	struct HasAddition : std::false_type { };

	template<class A, class B>
	struct HasAddition<A, B, std::void_t<decltype(std::declval<A>() + std::declval<B>())>> : std::true_type { };

	using PolarToCartesianMap = PolarToCartesian2DMap;

	TEST_CASE("CoordinateMap - identity map preserves typed points vectors covectors and forms", "[CoordinateMap][DifferentialForms]")
	{
		auto map = MakeIdentityCoordinateMap<Cartesian3, 3>();
		Point<Real, 3, Cartesian3> point{ REAL(1.0), REAL(2.0), REAL(3.0) };
		TangentVector<3, Cartesian3> vector{ REAL(4.0), REAL(5.0), REAL(6.0) };
		Covector<3, Cartesian3> covector{ REAL(7.0), REAL(8.0), REAL(9.0) };
		Form2<3, Cartesian3> form;
		form.SetAlternatingComponent(REAL(3.0), 0, 1);

		REQUIRE(map_point(map, point).coordinates() == point.coordinates());
		REQUIRE(push_forward(map, vector, point).components() == vector.components());
		REQUIRE(pull_back(map, covector, point).components() == covector.components());

		Form2<3, Cartesian3> pulledForm = pull_back(map, form, point);
		REQUIRE(pulledForm.Component(0, 1) == REAL(3.0));
		REQUIRE(pulledForm.Component(1, 0) == -REAL(3.0));
		REQUIRE(pulledForm.IsAlternating());
	}

	TEST_CASE("CoordinateMap - map_point adapts existing polar to Cartesian transform", "[CoordinateMap][DifferentialForms]")
	{
		CoordTransfPolarToCartesian2D transform;
		PolarToCartesianMap map = MakePolarToCartesian2DMap(transform);
		Point<Real, 2, Polar2> polarPoint{ REAL(2.0), Constants::PI / REAL(2.0) };

		Point<Real, 2, Cartesian2> cartPoint = map_point(map, polarPoint);

		REQUIRE(cartPoint[0] == Catch::Approx(REAL(0.0)).margin(CoordinateTolerance));
		REQUIRE(cartPoint[1] == Catch::Approx(REAL(2.0)).epsilon(CoordinateTolerance));
	}

	TEST_CASE("CoordinateMap - push_forward transforms tangent vectors by Jacobian", "[CoordinateMap][DifferentialForms]")
	{
		CoordTransfPolarToCartesian2D transform;
		PolarToCartesianMap map(transform);
		Point<Real, 2, Polar2> atPoint{ REAL(2.0), Constants::PI / REAL(2.0) };
		TangentVector<2, Polar2> radial{ REAL(1.0), REAL(0.0) };
		TangentVector<2, Polar2> angular{ REAL(0.0), REAL(1.0) };

		TangentVector<2, Cartesian2> radialPush = push_forward(map, radial, atPoint);
		TangentVector<2, Cartesian2> angularPush = map.push_forward(angular, atPoint);

		REQUIRE(radialPush[0] == Catch::Approx(REAL(0.0)).margin(CoordinateTolerance));
		REQUIRE(radialPush[1] == Catch::Approx(REAL(1.0)).epsilon(CoordinateTolerance));
		REQUIRE(angularPush[0] == Catch::Approx(-REAL(2.0)).epsilon(CoordinateTolerance));
		REQUIRE(angularPush[1] == Catch::Approx(REAL(0.0)).margin(CoordinateTolerance));

		static_assert(!HasAddition<TangentVector<2, Polar2>, TangentVector<2, Cartesian2>>::value,
			"Different-frame tangent vectors cannot be mixed before transformation");
	}

	TEST_CASE("CoordinateMap - pull_back transforms covectors by Jacobian transpose", "[CoordinateMap][DifferentialForms]")
	{
		CoordTransfPolarToCartesian2D transform;
		PolarToCartesianMap map(transform);
		Point<Real, 2, Polar2> atPoint{ REAL(2.0), Constants::PI / REAL(2.0) };
		Covector<2, Cartesian2> dx{ REAL(1.0), REAL(0.0) };
		Covector<2, Cartesian2> dy{ REAL(0.0), REAL(1.0) };

		Covector<2, Polar2> dxPull = pull_back(map, dx, atPoint);
		Covector<2, Polar2> dyPull = map.pull_back(dy, atPoint);

		REQUIRE(dxPull[0] == Catch::Approx(REAL(0.0)).margin(CoordinateTolerance));
		REQUIRE(dxPull[1] == Catch::Approx(-REAL(2.0)).epsilon(CoordinateTolerance));
		REQUIRE(dyPull[0] == Catch::Approx(REAL(1.0)).epsilon(CoordinateTolerance));
		REQUIRE(dyPull[1] == Catch::Approx(REAL(0.0)).margin(CoordinateTolerance));
	}

	TEST_CASE("CoordinateMap - pull_back supports one-forms and two-forms", "[CoordinateMap][DifferentialForms]")
	{
		CoordTransfPolarToCartesian2D transform;
		PolarToCartesianMap map(transform);
		Point<Real, 2, Polar2> atPoint{ REAL(2.0), Constants::PI / REAL(2.0) };
		Form1<2, Cartesian2> dx = BasisOneForm<0, 2, Cartesian2>();
		Form2<2, Cartesian2> area = Wedge(dx, BasisOneForm<1, 2, Cartesian2>());
		TangentVector<2, Polar2> radial{ REAL(1.0), REAL(0.0) };
		TangentVector<2, Polar2> angular{ REAL(0.0), REAL(1.0) };

		Form1<2, Polar2> dxPull = map.pull_back(dx, atPoint);
		Form2<2, Polar2> areaPull = pull_back(map, area, atPoint);

		REQUIRE(dxPull.Component(0) == Catch::Approx(REAL(0.0)).margin(CoordinateTolerance));
		REQUIRE(dxPull.Component(1) == Catch::Approx(-REAL(2.0)).epsilon(CoordinateTolerance));
		REQUIRE(areaPull.Component(0, 1) == Catch::Approx(REAL(2.0)).epsilon(CoordinateTolerance));
		REQUIRE(areaPull(radial, angular) == Catch::Approx(REAL(2.0)).epsilon(CoordinateTolerance));
	}

	TEST_CASE("CoordinateMap - cylindrical 3D map transforms points tangent vectors and covectors", "[CoordinateMap][DifferentialForms]")
	{
		CoordTransfCylindricalToCartesian transform;
		CylindricalToCartesian3DMap map = MakeCylindricalToCartesian3DMap(transform);
		Point<Real, 3, Cylindrical3> atPoint{ REAL(2.0), Constants::PI / REAL(2.0), REAL(5.0) };

		Point<Real, 3, Cartesian3> cart = map_point(map, atPoint);
		REQUIRE(cart[0] == Catch::Approx(REAL(0.0)).margin(CoordinateTolerance));
		REQUIRE(cart[1] == Catch::Approx(REAL(2.0)).epsilon(CoordinateTolerance));
		REQUIRE(cart[2] == Catch::Approx(REAL(5.0)).epsilon(CoordinateTolerance));

		TangentVector<3, Cylindrical3> radial{ REAL(1.0), REAL(0.0), REAL(0.0) };
		TangentVector<3, Cylindrical3> angular{ REAL(0.0), REAL(1.0), REAL(0.0) };
		TangentVector<3, Cylindrical3> vertical{ REAL(0.0), REAL(0.0), REAL(1.0) };

		auto radialPush = push_forward(map, radial, atPoint);
		auto angularPush = push_forward(map, angular, atPoint);
		auto verticalPush = push_forward(map, vertical, atPoint);

		REQUIRE(radialPush[0] == Catch::Approx(REAL(0.0)).margin(CoordinateTolerance));
		REQUIRE(radialPush[1] == Catch::Approx(REAL(1.0)).epsilon(CoordinateTolerance));
		REQUIRE(radialPush[2] == Catch::Approx(REAL(0.0)).margin(CoordinateTolerance));
		REQUIRE(angularPush[0] == Catch::Approx(-REAL(2.0)).epsilon(CoordinateTolerance));
		REQUIRE(angularPush[1] == Catch::Approx(REAL(0.0)).margin(CoordinateTolerance));
		REQUIRE(verticalPush[2] == Catch::Approx(REAL(1.0)).epsilon(CoordinateTolerance));

		Covector<3, Cartesian3> dx{ REAL(1.0), REAL(0.0), REAL(0.0) };
		Covector<3, Cartesian3> dy{ REAL(0.0), REAL(1.0), REAL(0.0) };
		Covector<3, Cartesian3> dz{ REAL(0.0), REAL(0.0), REAL(1.0) };

		auto dxPull = pull_back(map, dx, atPoint);
		auto dyPull = pull_back(map, dy, atPoint);
		auto dzPull = pull_back(map, dz, atPoint);

		REQUIRE(dxPull[0] == Catch::Approx(REAL(0.0)).margin(CoordinateTolerance));
		REQUIRE(dxPull[1] == Catch::Approx(-REAL(2.0)).epsilon(CoordinateTolerance));
		REQUIRE(dxPull[2] == Catch::Approx(REAL(0.0)).margin(CoordinateTolerance));
		REQUIRE(dyPull[0] == Catch::Approx(REAL(1.0)).epsilon(CoordinateTolerance));
		REQUIRE(dyPull[1] == Catch::Approx(REAL(0.0)).margin(CoordinateTolerance));
		REQUIRE(dzPull[2] == Catch::Approx(REAL(1.0)).epsilon(CoordinateTolerance));
	}

	TEST_CASE("CoordinateMap - cylindrical 3D pull_back supports higher-degree forms", "[CoordinateMap][DifferentialForms]")
	{
		CoordTransfCylindricalToCartesian transform;
		CylindricalToCartesian3DMap map = MakeCylindricalToCartesian3DMap(transform);
		Point<Real, 3, Cylindrical3> atPoint{ REAL(2.0), Constants::PI / REAL(2.0), REAL(5.0) };
		Form3<3, Cartesian3> volume;
		volume.SetAlternatingComponent(REAL(1.0), 0, 1, 2);

		Form3<3, Cylindrical3> pulledVolume = pull_back(map, volume, atPoint);

		REQUIRE(pulledVolume.IsAlternating());
		REQUIRE(pulledVolume.Component(0, 1, 2) == Catch::Approx(REAL(2.0)).epsilon(CoordinateTolerance));
		REQUIRE(pulledVolume.Component(1, 0, 2) == Catch::Approx(-REAL(2.0)).epsilon(CoordinateTolerance));
	}
}
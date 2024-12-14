#include "../catch/catch.hpp"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "core/CoordTransf.h"
#include "core/CoordTransf/CoordTransfSpherical.h"
#include "core/CoordTransf/CoordTransfCylindrical.h"
#include "core/FieldOperations.h"
#endif

#include "../test_data/Fields.h"

using namespace MML;

// TODO - HIGH!!! Test coord transf
// make sure it is clear what convention is used
namespace MML::Tests::Core::CoordTransfTests
{
	/*********************************************************************/
	/*****              CoordTransfSpherToCart tests                 *****/
	/*********************************************************************/


	TEST_CASE("Test_CoordTransf_Cartesian_to_Spherical", "[simple]")
	{
		Vector3Spherical posSpher;

		posSpher = CoordTransfCartToSpher.transf(Vector3Cartesian{ 1.0, 0.0, 0.0 });
		REQUIRE(posSpher == Vector3Spherical{ 1.0, Constants::PI / 2, 0.0 });

		posSpher = CoordTransfCartToSpher.transf(Vector3Cartesian{ -1.0, 0.0, 0.0 });
		REQUIRE(posSpher == Vector3Spherical{ 1.0, Constants::PI / 2, Constants::PI });

		posSpher = CoordTransfCartToSpher.transf(Vector3Cartesian{ 0.0, 1.0, 0.0 });
		REQUIRE(posSpher == Vector3Spherical{ 1.0, Constants::PI / 2, Constants::PI / 2 });

		REQUIRE(CoordTransfCartToSpher.transf(Vec3Cart{ 0.0, -1.0, 0.0 }) == Vec3Sph{ 1.0, Constants::PI / 2, -Constants::PI / 2 });

		posSpher = CoordTransfCartToSpher.transf(Vector3Cartesian{ 0.0, 0.0, 1.0 });
		REQUIRE(posSpher == Vector3Spherical{ 1.0, 0.0, 0.0 });

		posSpher = CoordTransfCartToSpher.transf(Vector3Cartesian{ 0.0, 0.0, -1.0 });
		REQUIRE(posSpher == Vector3Spherical{ 1.0, Constants::PI, 0.0 });
	}

	TEST_CASE("Test_CoordTransf_Cartesian_to_Cylindrical", "[simple]")
	{
		Vector3Cylindrical posCyl;

		posCyl = CoordTransfCartToCyl.transf(Vector3Cartesian{ 1.0, 0.0, 0.0 });
		REQUIRE(posCyl == Vector3Cylindrical{ 1.0, 0.0, 0.0 });

		posCyl = CoordTransfCartToCyl.transf(Vector3Cartesian{ -1.0, 0.0, 0.0 });
		REQUIRE(posCyl == Vector3Cylindrical{ 1.0, Constants::PI, 0.0 });

		posCyl = CoordTransfCartToCyl.transf(Vector3Cartesian{ 0.0, 1.0, 0.0 });
		REQUIRE(posCyl == Vector3Cylindrical{ 1.0, Constants::PI / 2, 0.0 });

		posCyl = CoordTransfCartToCyl.transf(Vector3Cartesian{ 0.0, -1.0, 0.0 });
		REQUIRE(posCyl == Vector3Cylindrical{ 1.0, -Constants::PI / 2, 0.0 });

		posCyl = CoordTransfCartToCyl.transf(Vector3Cartesian{ 0.0, 0.0, 1.0 });
		REQUIRE(posCyl == Vector3Cylindrical{ 0.0, 0.0, 1.0 });

		posCyl = CoordTransfCartToCyl.transf(Vector3Cartesian{ 0.0, 0.0, -1.0 });
		REQUIRE(posCyl == Vector3Cylindrical{ 0.0, 0.0, -1.0 });
	}

	TEST_CASE("Test_CoordTransf_Spherical_to_Cartesian", "[simple]")
	{
		Vector3Cartesian posCart;

		posCart = CoordTransfSpherToCart.transf(Vector3Spherical{ 1.0, 0.0, 0.0 });
		REQUIRE(posCart == Vector3Cartesian{ 0.0, 0.0, 1.0 });

		posCart = CoordTransfSpherToCart.transf(Vector3Spherical{ 1.0, Constants::PI / 2, 0.0 });
		REQUIRE(posCart.IsEqualTo(Vector3Cartesian{ 1.0, 0.0, 0.0 }, 1e-16));

		posCart = CoordTransfSpherToCart.transf(Vector3Spherical{ 1.0, Constants::PI / 2, Constants::PI / 2 });
		REQUIRE(posCart.IsEqualTo(Vector3Cartesian{ 0.0, 1.0, 0.0 }, 1e-16));
	}

	TEST_CASE("Test_GetUnitVector")
	{
		Vector3Cartesian p1{ 2.0, 1.0, -2.0 };
		auto p1Spher = CoordTransfCartToSpher.transf(p1);

		// Vector3Cartesian vec_i = CoordTransfSpherToCart.getUnitVector(0, p1Spher);
		// Vector3Cartesian vec_j = CoordTransfSpherToCart.getUnitVector(1, p1Spher);
		// Vector3Cartesian vec_k = CoordTransfSpherToCart.getUnitVector(2, p1Spher);
		// std::cout << "Unit vector i        : " << vec_i.GetUnitVector() << std::endl;
		// std::cout << "Unit vector j        : " << vec_j.GetUnitVector() << std::endl;
		// std::cout << "Unit vector k        : " << vec_k.GetUnitVector() << std::endl<< std::endl;

		// // Explicit formulas for unit vectors in spherical coordinates
		// double r = p1Spher[0];
		// double theta = p1Spher[1];
		// double phi = p1Spher[2];
		// Vector3Spherical vec_i2{ sin(theta) * cos(phi), cos(theta) * cos(phi), -sin(theta) };
		// Vector3Spherical vec_j2{ sin(theta) * sin(phi), cos(theta) * sin(phi), cos(theta) };
		// Vector3Spherical vec_k2{ cos(theta), -sin(theta), 0.0 };
		// std::cout << "Unit vector i (calc) : " << vec_i2 << std::endl;
		// std::cout << "Unit vector j (calc) : " << vec_j2 << std::endl;
		// std::cout << "Unit vector k (calc) : " << vec_k2 << std::endl << std::endl;

		// Vector3Cartesian vec_r{ sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta) };
		// Vector3Cartesian vec_theta{ cos(theta) * cos(phi), cos(theta) * sin(phi), -sin(theta) };
		// Vector3Cartesian vec_phi{ -sin(phi), cos(phi), 0.0 };
		// std::cout << "Unit vector r (calc) : " << vec_r << std::endl;
		// std::cout << "Unit vector theta    : " << vec_theta << std::endl;
		// std::cout << "Unit vector phi      : " << vec_phi << std::endl;
	}
}
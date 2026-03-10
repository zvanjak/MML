#include <catch2/catch_all.hpp>
#include "../../TestPrecision.h"
#include "../../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "core/CoordTransf.h"
#include "core/CoordTransf/CoordTransfSpherical.h"
#include "core/CoordTransf/CoordTransfCylindrical.h"
#include "core/FieldOperations.h"
#include "core/Fields.h"
#endif


using namespace MML;
using namespace MML::Testing;

// TODO - HIGH!!! Test coord transf
// make sure it is clear what convention is used
namespace MML::Tests::Core::CoordTransfTests
{
	/*********************************************************************/
	/*****              CoordTransfSpherToCart tests                 *****/
	/*********************************************************************/


	TEST_CASE("Test_CoordTransf_Cartesian_to_Spherical", "[simple]")
	{
			TEST_PRECISION_INFO();
		Vector3Spherical posSpher;

		posSpher = CoordTransfCartToSpher.transf(Vector3Cartesian{ REAL(1.0), REAL(0.0), REAL(0.0) });
		REQUIRE(posSpher == Vector3Spherical{ REAL(1.0), Constants::PI / 2, REAL(0.0) });

		posSpher = CoordTransfCartToSpher.transf(Vector3Cartesian{ -REAL(1.0), REAL(0.0), REAL(0.0) });
		REQUIRE(posSpher == Vector3Spherical{ REAL(1.0), Constants::PI / 2, Constants::PI });

		posSpher = CoordTransfCartToSpher.transf(Vector3Cartesian{ REAL(0.0), REAL(1.0), REAL(0.0) });
		REQUIRE(posSpher == Vector3Spherical{ REAL(1.0), Constants::PI / 2, Constants::PI / 2 });

		REQUIRE(CoordTransfCartToSpher.transf(Vec3Cart{ REAL(0.0), -REAL(1.0), REAL(0.0) }) == Vec3Sph{ REAL(1.0), Constants::PI / 2, -Constants::PI / 2 });

		posSpher = CoordTransfCartToSpher.transf(Vector3Cartesian{ REAL(0.0), REAL(0.0), REAL(1.0) });
		REQUIRE(posSpher == Vector3Spherical{ REAL(1.0), REAL(0.0), REAL(0.0) });

		posSpher = CoordTransfCartToSpher.transf(Vector3Cartesian{ REAL(0.0), REAL(0.0), -REAL(1.0) });
		REQUIRE(posSpher == Vector3Spherical{ REAL(1.0), Constants::PI, REAL(0.0) });
	}

	TEST_CASE("Test_CoordTransf_Cartesian_to_Cylindrical", "[simple]")
	{
			TEST_PRECISION_INFO();
		Vector3Cylindrical posCyl;

		posCyl = CoordTransfCartToCyl.transf(Vector3Cartesian{ REAL(1.0), REAL(0.0), REAL(0.0) });
		REQUIRE(posCyl == Vector3Cylindrical{ REAL(1.0), REAL(0.0), REAL(0.0) });

		posCyl = CoordTransfCartToCyl.transf(Vector3Cartesian{ -REAL(1.0), REAL(0.0), REAL(0.0) });
		REQUIRE(posCyl == Vector3Cylindrical{ REAL(1.0), Constants::PI, REAL(0.0) });

		posCyl = CoordTransfCartToCyl.transf(Vector3Cartesian{ REAL(0.0), REAL(1.0), REAL(0.0) });
		REQUIRE(posCyl == Vector3Cylindrical{ REAL(1.0), Constants::PI / 2, REAL(0.0) });

		posCyl = CoordTransfCartToCyl.transf(Vector3Cartesian{ REAL(0.0), -REAL(1.0), REAL(0.0) });
		REQUIRE(posCyl == Vector3Cylindrical{ REAL(1.0), -Constants::PI / 2, REAL(0.0) });

		posCyl = CoordTransfCartToCyl.transf(Vector3Cartesian{ REAL(0.0), REAL(0.0), REAL(1.0) });
		REQUIRE(posCyl == Vector3Cylindrical{ REAL(0.0), REAL(0.0), REAL(1.0) });

		posCyl = CoordTransfCartToCyl.transf(Vector3Cartesian{ REAL(0.0), REAL(0.0), -REAL(1.0) });
		REQUIRE(posCyl == Vector3Cylindrical{ REAL(0.0), REAL(0.0), -REAL(1.0) });
	}

	TEST_CASE("Test_CoordTransf_Spherical_to_Cartesian", "[simple]")
	{
			TEST_PRECISION_INFO();
		Vector3Cartesian posCart;

		posCart = CoordTransfSpherToCart.transf(Vector3Spherical{ REAL(1.0), REAL(0.0), REAL(0.0) });
		REQUIRE(posCart == Vector3Cartesian{ REAL(0.0), REAL(0.0), REAL(1.0) });

		posCart = CoordTransfSpherToCart.transf(Vector3Spherical{ REAL(1.0), Constants::PI / 2, REAL(0.0) });
		REQUIRE(posCart.IsEqualTo(Vector3Cartesian{ REAL(1.0), REAL(0.0), REAL(0.0) }, 1e-16));

		posCart = CoordTransfSpherToCart.transf(Vector3Spherical{ REAL(1.0), Constants::PI / 2, Constants::PI / 2 });
		REQUIRE(posCart.IsEqualTo(Vector3Cartesian{ REAL(0.0), REAL(1.0), REAL(0.0) }, 1e-16));
	}

	TEST_CASE("Test_GetUnitVector")
	{
			TEST_PRECISION_INFO();
		Vector3Cartesian p1{ REAL(2.0), REAL(1.0), -REAL(2.0) };
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
		// Vector3Spherical vec_k2{ cos(theta), -sin(theta), REAL(0.0) };
		// std::cout << "Unit vector i (calc) : " << vec_i2 << std::endl;
		// std::cout << "Unit vector j (calc) : " << vec_j2 << std::endl;
		// std::cout << "Unit vector k (calc) : " << vec_k2 << std::endl << std::endl;

		// Vector3Cartesian vec_r{ sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta) };
		// Vector3Cartesian vec_theta{ cos(theta) * cos(phi), cos(theta) * sin(phi), -sin(theta) };
		// Vector3Cartesian vec_phi{ -sin(phi), cos(phi), REAL(0.0) };
		// std::cout << "Unit vector r (calc) : " << vec_r << std::endl;
		// std::cout << "Unit vector theta    : " << vec_theta << std::endl;
		// std::cout << "Unit vector phi      : " << vec_phi << std::endl;
	}

	/*********************************************************************/
	/*****     VectorFrom/VectorTo type consistency tests            *****/
	/*****     (om8d.16: Ensure correct return types when distinct)  *****/
	/*********************************************************************/

	TEST_CASE("CoordTransf - transfVecContravariant returns VectorTo type", "[coordtransf][types]")
	{
		TEST_PRECISION_INFO();
		// Use spherical->cartesian where VectorFrom=Vector3Spherical, VectorTo=Vector3Cartesian
		// These are distinct types, so this test verifies the return type is correct
		
		Vector3Spherical posSpher{ REAL(2.0), Constants::PI / 4, Constants::PI / 6 };
		Vector3Spherical vecSpher{ REAL(1.0), REAL(0.5), REAL(0.3) };  // A contravariant vector in spherical coords
		
		// transfVecContravariant should return VectorTo (Vector3Cartesian)
		Vector3Cartesian result = CoordTransfSpherToCart.transfVecContravariant(vecSpher, posSpher);
		
		// Verify the result is a valid Vector3Cartesian with finite values
		REQUIRE(std::isfinite(result[0]));
		REQUIRE(std::isfinite(result[1]));
		REQUIRE(std::isfinite(result[2]));
		
		// The transformed vector should have non-trivial values
		REQUIRE(result.NormL2() > 0);
	}

	TEST_CASE("CoordTransf - getBasisVec returns VectorTo type", "[coordtransf][types]")
	{
		TEST_PRECISION_INFO();
		// Use spherical->cartesian where VectorFrom=Vector3Spherical, VectorTo=Vector3Cartesian
		
		Vector3Spherical posSpher{ REAL(2.0), Constants::PI / 4, Constants::PI / 6 };
		
		// getBasisVec should return VectorTo (Vector3Cartesian)
		Vector3Cartesian e_r = CoordTransfSpherToCart.getBasisVec(0, posSpher);
		Vector3Cartesian e_theta = CoordTransfSpherToCart.getBasisVec(1, posSpher);
		Vector3Cartesian e_phi = CoordTransfSpherToCart.getBasisVec(2, posSpher);
		
		// All basis vectors should be finite
		REQUIRE(std::isfinite(e_r.NormL2()));
		REQUIRE(std::isfinite(e_theta.NormL2()));
		REQUIRE(std::isfinite(e_phi.NormL2()));
		
		// Basis vectors should be non-zero (except potentially at singularities)
		REQUIRE(e_r.NormL2() > REAL(0.1));
	}

	TEST_CASE("CoordTransfWithInverse - getInverseContravarBasisVec returns VectorTo type", "[coordtransf][types]")
	{
		TEST_PRECISION_INFO();
		// Use spherical->cartesian where VectorFrom=Vector3Spherical, VectorTo=Vector3Cartesian
		
		Vector3Cartesian posCart{ REAL(1.0), REAL(1.0), REAL(1.0) };
		
		// getInverseContravarBasisVec should return VectorTo (Vector3Cartesian)
		Vector3Cartesian result = CoordTransfSpherToCart.getInverseContravarBasisVec(0, posCart);
		
		// Result should be finite
		REQUIRE(std::isfinite(result[0]));
		REQUIRE(std::isfinite(result[1]));
		REQUIRE(std::isfinite(result[2]));
	}

	TEST_CASE("CoordTransfWithInverse - transfVecCovariant returns VectorTo type", "[coordtransf][types]")
	{
		TEST_PRECISION_INFO();
		// Use spherical->cartesian where VectorFrom=Vector3Spherical, VectorTo=Vector3Cartesian
		
		Vector3Cartesian posCart{ REAL(1.0), REAL(1.0), REAL(1.0) };
		Vector3Spherical gradSpher{ REAL(1.0), REAL(0.5), REAL(0.3) };  // A covariant vector (gradient) in spherical
		
		// transfVecCovariant should return VectorTo (Vector3Cartesian)
		Vector3Cartesian result = CoordTransfSpherToCart.transfVecCovariant(gradSpher, posCart);
		
		// Result should be finite
		REQUIRE(std::isfinite(result[0]));
		REQUIRE(std::isfinite(result[1]));
		REQUIRE(std::isfinite(result[2]));
	}
}

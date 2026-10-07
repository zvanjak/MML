#include <catch2/catch_all.hpp>
#include "../../TestPrecision.h"
#include "../../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/core/CoordTransf/CoordTransfBase.h>
#include <mml/core/CoordTransf/CoordTransfSpherical.h>
#include <mml/core/CoordTransf/CoordTransfCylindrical.h>
#include <mml/core/Fields/FieldOperations.h>
#include <mml/core/Fields/Fields.h>
#endif


using namespace MML;
using namespace MML::Testing;

// TODO - HIGH!!! Test coord transf
// make sure it is clear what convention is used
namespace MML::Tests::Core::CoordTransfTests
{
	class TrackingCylindricalTransform : public CoordTransfCylindricalToCartesian
	{
	public:
		mutable int forwardJacobianCalls = 0;
		mutable int inverseJacobianCalls = 0;

		MatrixNM<Real, 3, 3> jacobian(const VectorN<Real, 3>&) const override
		{
			++forwardJacobianCalls;
			MatrixNM<Real, 3, 3> jac;
			for (int i = 0; i < 3; ++i)
				for (int j = 0; j < 3; ++j)
					jac(i, j) = i == j ? REAL(2.0) : REAL(0.0);
			return jac;
		}

		MatrixNM<Real, 3, 3> inverseJacobian(const VectorN<Real, 3>&) const override
		{
			++inverseJacobianCalls;
			MatrixNM<Real, 3, 3> jac;
			for (int i = 0; i < 3; ++i)
				for (int j = 0; j < 3; ++j)
					jac(i, j) = i == j ? REAL(3.0) : REAL(0.0);
			return jac;
		}
	};

	/*********************************************************************/
	/*****              CoordTransfSpherToCart tests                 *****/
	/*********************************************************************/


	TEST_CASE("Test_CoordTransf_Cartesian_to_Spherical", "[simple]")
	{
			TEST_PRECISION_INFO();
		Vector3Spherical posSpher;

		posSpher = CoordTransfCartToSpher.transf(Vector3Cartesian{ REAL(1.0), REAL(0.0), REAL(0.0) });
		REQUIRE(posSpher.IsEqualTo(Vector3Spherical{ REAL(1.0), Constants::PI / 2, REAL(0.0) }));

		posSpher = CoordTransfCartToSpher.transf(Vector3Cartesian{ -REAL(1.0), REAL(0.0), REAL(0.0) });
		REQUIRE(posSpher.IsEqualTo(Vector3Spherical{ REAL(1.0), Constants::PI / 2, Constants::PI }));

		posSpher = CoordTransfCartToSpher.transf(Vector3Cartesian{ REAL(0.0), REAL(1.0), REAL(0.0) });
		REQUIRE(posSpher.IsEqualTo(Vector3Spherical{ REAL(1.0), Constants::PI / 2, Constants::PI / 2 }));

		REQUIRE(CoordTransfCartToSpher.transf(Vec3Cart{ REAL(0.0), -REAL(1.0), REAL(0.0) }).IsEqualTo(Vec3Sph{ REAL(1.0), Constants::PI / 2, -Constants::PI / 2 }));

		posSpher = CoordTransfCartToSpher.transf(Vector3Cartesian{ REAL(0.0), REAL(0.0), REAL(1.0) });
		REQUIRE(posSpher.IsEqualTo(Vector3Spherical{ REAL(1.0), REAL(0.0), REAL(0.0) }));

		posSpher = CoordTransfCartToSpher.transf(Vector3Cartesian{ REAL(0.0), REAL(0.0), -REAL(1.0) });
		REQUIRE(posSpher.IsEqualTo(Vector3Spherical{ REAL(1.0), Constants::PI, REAL(0.0) }));
	}

	TEST_CASE("Test_CoordTransf_Cartesian_to_Cylindrical", "[simple]")
	{
			TEST_PRECISION_INFO();
		Vector3Cylindrical posCyl;

		posCyl = CoordTransfCartToCyl.transf(Vector3Cartesian{ REAL(1.0), REAL(0.0), REAL(0.0) });
		REQUIRE(posCyl == Vector3Cylindrical{ REAL(1.0), REAL(0.0), REAL(0.0) });

		posCyl = CoordTransfCartToCyl.transf(Vector3Cartesian{ -REAL(1.0), REAL(0.0), REAL(0.0) });
		REQUIRE(posCyl.IsEqualTo(Vector3Cylindrical{ REAL(1.0), Constants::PI, REAL(0.0) }));

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
		REQUIRE(posCart.IsEqualTo(Vector3Cartesian{ REAL(1.0), REAL(0.0), REAL(0.0) }, TOL(1e-16, 1e-6)));

		posCart = CoordTransfSpherToCart.transf(Vector3Spherical{ REAL(1.0), Constants::PI / 2, Constants::PI / 2 });
		REQUIRE(posCart.IsEqualTo(Vector3Cartesian{ REAL(0.0), REAL(1.0), REAL(0.0) }, TOL(1e-16, 1e-6)));
	}

	TEST_CASE("Standard coordinate transforms provide analytical Jacobians", "[coordtransf][jacobian]")
	{
		TEST_PRECISION_INFO();
		const Vector3Spherical sphericalPos{ REAL(2.0), Constants::PI / 2, REAL(0.0) };
		const CoordTransf<Vector3Spherical, Vector3Cartesian, 3>& sphericalBase = CoordTransfSpherToCart;
		const auto sphericalJac = sphericalBase.jacobian(sphericalPos);

		REQUIRE(sphericalJac(0, 0) == Catch::Approx(REAL(1.0)).margin(TOL(1e-15, 1e-5)));
		REQUIRE(sphericalJac(1, 2) == Catch::Approx(REAL(2.0)).margin(TOL(1e-15, 1e-5)));
		REQUIRE(sphericalJac(2, 1) == Catch::Approx(-REAL(2.0)).margin(TOL(1e-15, 1e-5)));

		const Vector3Cylindrical cylindricalPos{ REAL(2.0), Constants::PI / 2, REAL(3.0) };
		const CoordTransf<Vector3Cylindrical, Vector3Cartesian, 3>& cylindricalBase = CoordTransfCylToCart;
		const auto cylindricalJac = cylindricalBase.jacobian(cylindricalPos);

		REQUIRE(cylindricalJac(0, 1) == Catch::Approx(-REAL(2.0)).margin(TOL(1e-15, 1e-5)));
		REQUIRE(cylindricalJac(1, 0) == Catch::Approx(REAL(1.0)).margin(TOL(1e-15, 1e-5)));
		REQUIRE(cylindricalJac(2, 2) == REAL(1.0));
	}

	TEST_CASE("Analytical coordinate Jacobians compose with their inverses", "[coordtransf][jacobian]")
	{
		TEST_PRECISION_INFO();
		const Vector3Spherical sphericalPos{ REAL(2.5), REAL(1.1), REAL(0.7) };
		const auto cartesianPos = CoordTransfSpherToCart.transf(sphericalPos);
		const auto forwardJac = CoordTransfSpherToCart.jacobian(sphericalPos);
		const auto inverseJac = CoordTransfSpherToCart.inverseJacobian(cartesianPos);

		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j)
			{
				Real product = REAL(0.0);
				for (int k = 0; k < 3; ++k)
					product += inverseJac(i, k) * forwardJac(k, j);
				REQUIRE(product == Catch::Approx(i == j ? REAL(1.0) : REAL(0.0)).margin(TOL(1e-14, 1e-4)));
			}
	}

	TEST_CASE("Coordinate helpers dispatch through virtual Jacobian hooks once", "[coordtransf][jacobian]")
	{
		TEST_PRECISION_INFO();
		TrackingCylindricalTransform transform;
		const CoordTransf<Vector3Cylindrical, Vector3Cartesian, 3>& forwardBase = transform;
		const CoordTransfWithInverse<Vector3Cylindrical, Vector3Cartesian, 3>& inverseBase = transform;
		const Vector3Cylindrical sourceVec{ REAL(1.0), REAL(2.0), REAL(3.0) };
		const Vector3Cylindrical sourcePos{ REAL(2.0), REAL(0.4), REAL(1.0) };
		const Vector3Cartesian targetPos = transform.transf(sourcePos);

		const auto contravariant = forwardBase.transfVecContravariant(sourceVec, sourcePos);
		REQUIRE(contravariant == Vector3Cartesian{ REAL(2.0), REAL(4.0), REAL(6.0) });
		REQUIRE(transform.forwardJacobianCalls == 1);

		const auto covariant = inverseBase.transfVecCovariant(sourceVec, targetPos);
		REQUIRE(covariant == Vector3Cartesian{ REAL(3.0), REAL(6.0), REAL(9.0) });
		REQUIRE(transform.inverseJacobianCalls == 1);

		transform.forwardJacobianCalls = 0;
		transform.inverseJacobianCalls = 0;
		Tensor2<3> tensor(1, 1);
		tensor(0, 0) = REAL(1.0);
		inverseBase.transfTensor2(tensor, sourcePos);
		REQUIRE(transform.forwardJacobianCalls == 1);
		REQUIRE(transform.inverseJacobianCalls == 1);
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

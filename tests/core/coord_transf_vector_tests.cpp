#include "../catch/catch.hpp"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "core/CoordTransf.h"
#include "core/CoordTransf/CoordTransfSpherical.h"

#include "core/FieldOperations.h"
#endif

#include "../test_data/Fields.h"

using namespace MML;

namespace MML::Tests::Core::CoordTransfVectorTests
{

	TEST_CASE("Test_Contravariant_transf_cart_to_spher")
	{
		// TODO 0.9 - dodati još 5-6 test caseova
		Vector3Cartesian v_cart{ 1.0, 1.0, 0.0 };
		Vector3Cartesian posCart{ 1.0, 1.0, 0.0 };

		Vector3Spherical v_transf_to_spher = CoordTransfCartToSpher.transfVecContravariant(v_cart, posCart);
		REQUIRE(true == v_transf_to_spher.IsEqualTo(Vector3Spherical(sqrt(2), 0.0, 0.0), 1e-8));

		Vector3Spherical x1_spher{ CoordTransfCartToSpher.transf(posCart) };
		Vector3Cartesian v_back_transf_to_cart = CoordTransfSpherToCart.transfVecContravariant(v_transf_to_spher, x1_spher);
		REQUIRE(true == v_back_transf_to_cart.IsEqualTo(v_cart, 1e-7));

		double r = 2;
		double phi = Constants::PI / 3;
		Vec3Sph contravar_sph{1.0, r, 0.0};
		Vec3Sph covar_sph{0.0, -r*r, cos(phi) * cos(phi)};

		Vec3Cart res = CoordTransfSpherToCart.transfVecContravariant(contravar_sph, Vector3Spherical{ r, Constants::PI / 2, phi });

		Vector3Cartesian v_cart1{ 3.0, 4.0, 5.0 };
		Vector3Cartesian posCart1{ 1.0, 2.0, 3.0 };
		Vector3Spherical v_spher_real{6.9488, 0.95622, -0.89445};
		Vector3Spherical v_spher1{6.9487922897236354, 0.25555062599984346	, -0.4};

		 v_transf_to_spher = CoordTransfCartToSpher.transfVecContravariant(v_cart1, posCart1);
		 REQUIRE(true == v_transf_to_spher.IsEqualTo(v_spher1, 1e-8));

		 // TODO - pokazati kako se do real dolazi projekcijom na bazne vektore

	}

	TEST_CASE("Test_Covariant_transf_cart_to_spher")
	{
		Vector3Cartesian p_cart{ 1.0, 1.0, 1.0 };
		Vector3Spherical p_spher(CoordTransfSpherToCart.transfInverse(p_cart));

		ScalarFunction<3> fPotCart(Fields::InverseRadialPotentialFieldCart);
		Vector3Cartesian  grad_cart = ScalarFieldOperations::GradientCart<3>(fPotCart, p_cart);

		Vector3Spherical grad_transf_to_spher = CoordTransfCartToSpher.transfVecCovariant(grad_cart, p_spher);
		// TODO 0.9 - nekad je -1/3 bilo sqrt(2)???
		REQUIRE(true == grad_transf_to_spher.IsEqualTo(Vector3Spherical(-1.0 / 3, 0.0, 0.0), 1e-6));

		Vector3Cartesian back_transf_to_cart = CoordTransfSpherToCart.transfVecCovariant(grad_transf_to_spher, p_cart);
		REQUIRE(true == back_transf_to_cart.IsEqualTo(grad_cart, 1e-7));
	}

} // namespace MML::Tests::Core::CoordTransfVectorTests
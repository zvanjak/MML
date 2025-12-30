#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "core/CoordTransf.h"
#include "core/CoordTransf/CoordTransfSpherical.h"

#include "core/FieldOperations.h"
#include "core/Fields.h"
#endif


using namespace MML;
using namespace MML::Testing;

namespace MML::Tests::Core::CoordTransfVectorTests
{

	TEST_CASE("Test_Contravariant_transf_cart_to_spher")
	{
			TEST_PRECISION_INFO();
		// TODO REAL(0.9) - dodati jos 5-6 test caseova
		Vector3Cartesian v_cart{ REAL(1.0), REAL(1.0), REAL(0.0) };
		Vector3Cartesian posCart{ REAL(1.0), REAL(1.0), REAL(0.0) };

		Vector3Spherical v_transf_to_spher = CoordTransfCartToSpher.transfVecContravariant(v_cart, posCart);
		REQUIRE(true == v_transf_to_spher.IsEqualTo(Vector3Spherical(sqrt(2), REAL(0.0), REAL(0.0)), 1e-8));

		Vector3Spherical x1_spher{ CoordTransfCartToSpher.transf(posCart) };
		Vector3Cartesian v_back_transf_to_cart = CoordTransfSpherToCart.transfVecContravariant(v_transf_to_spher, x1_spher);
		REQUIRE(true == v_back_transf_to_cart.IsEqualTo(v_cart, 1e-7));

		double r = 2;
		double phi = Constants::PI / 3;
		Vec3Sph contravar_sph{REAL(1.0), REAL(r), REAL(0.0)};
		Vec3Sph covar_sph{REAL(0.0), REAL(-r*r), REAL(cos(phi) * cos(phi))};

		Vec3Cart res = CoordTransfSpherToCart.transfVecContravariant(contravar_sph, Vector3Spherical{ REAL(r), Constants::PI / REAL(2), REAL(phi) });

		Vector3Cartesian v_cart1{ REAL(3.0), REAL(4.0), REAL(5.0) };
		Vector3Cartesian posCart1{ REAL(1.0), REAL(2.0), REAL(3.0) };
		Vector3Spherical v_spher_real{REAL(6.9488), REAL(0.95622), -REAL(0.89445)};
		Vector3Spherical v_spher1{REAL(6.9487922897236354), REAL(0.25555062599984346)	, -REAL(0.4)};

		 v_transf_to_spher = CoordTransfCartToSpher.transfVecContravariant(v_cart1, posCart1);
		 REQUIRE(true == v_transf_to_spher.IsEqualTo(v_spher1, 1e-8));

		 // TODO - pokazati kako se do real dolazi projekcijom na bazne vektore

	}

	TEST_CASE("Test_Covariant_transf_cart_to_spher")
	{
			TEST_PRECISION_INFO();
		Vector3Cartesian p_cart{ REAL(1.0), REAL(1.0), REAL(1.0) };
		Vector3Spherical p_spher(CoordTransfSpherToCart.transfInverse(p_cart));

		ScalarFunction<3> fPotCart(Fields::InverseRadialPotentialFieldCart);
		Vector3Cartesian  grad_cart = ScalarFieldOperations::GradientCart<3>(fPotCart, p_cart);

		Vector3Spherical grad_transf_to_spher = CoordTransfCartToSpher.transfVecCovariant(grad_cart, p_spher);
		// TODO REAL(0.9) - nekad je -1/3 bilo sqrt(2)???
		REQUIRE(true == grad_transf_to_spher.IsEqualTo(Vector3Spherical(-REAL(1.0) / 3, REAL(0.0), REAL(0.0)), 1e-6));

		Vector3Cartesian back_transf_to_cart = CoordTransfSpherToCart.transfVecCovariant(grad_transf_to_spher, p_cart);
		REQUIRE(true == back_transf_to_cart.IsEqualTo(grad_cart, 1e-7));
	}

} // namespace MML::Tests::Core::CoordTransfVectorTests
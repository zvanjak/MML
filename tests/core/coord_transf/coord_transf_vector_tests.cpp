#include <catch2/catch_all.hpp>
#include "../../TestPrecision.h"
#include "../../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/core/CoordTransf/CoordTransfBase.h>
#include <mml/core/CoordTransf/CoordTransfSpherical.h>

#include <mml/core/Fields/FieldOperations.h>
#include <mml/core/Fields/Fields.h>
#endif


using namespace MML;
using namespace MML::Testing;

namespace MML::Tests::Core::CoordTransfVectorTests
{

	TEST_CASE("Test_Contravariant_transf_cart_to_spher")
	{
			TEST_PRECISION_INFO();
		// TODO static_cast<MML::Real>(0.9) - dodati jos 5-6 test caseova
		Vector3Cartesian v_cart{ static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0) };
		Vector3Cartesian posCart{ static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0) };

		Vector3Spherical v_transf_to_spher = CoordTransfCartToSpher.transfVecContravariant(v_cart, posCart);
		REQUIRE(true == v_transf_to_spher.IsEqualTo(Vector3Spherical(sqrt(Real(2)), static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0)), TOL(1e-8, 1e-4)));

		Vector3Spherical x1_spher{ CoordTransfCartToSpher.transf(posCart) };
		Vector3Cartesian v_back_transf_to_cart = CoordTransfSpherToCart.transfVecContravariant(v_transf_to_spher, x1_spher);
		REQUIRE(true == v_back_transf_to_cart.IsEqualTo(v_cart, TOL(1e-7, 1e-3)));

		double r = 2;
		double phi = Constants::PI / 3;
		Vec3Sph contravar_sph{static_cast<MML::Real>(1.0), static_cast<MML::Real>(r), static_cast<MML::Real>(0.0)};
		Vec3Sph covar_sph{static_cast<MML::Real>(0.0), static_cast<MML::Real>(-r*r), static_cast<MML::Real>(cos(phi) * cos(phi))};

		Vec3Cart res = CoordTransfSpherToCart.transfVecContravariant(contravar_sph, Vector3Spherical{ static_cast<MML::Real>(r), Constants::PI / static_cast<MML::Real>(2), static_cast<MML::Real>(phi) });
		for (int i = 0; i < 3; ++i)
			REQUIRE(std::isfinite(res[i]));

		Vector3Cartesian v_cart1{ static_cast<MML::Real>(3.0), static_cast<MML::Real>(4.0), static_cast<MML::Real>(5.0) };
		Vector3Cartesian posCart1{ static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0) };
		Vector3Spherical v_spher_real{static_cast<MML::Real>(6.9488), static_cast<MML::Real>(0.95622), -static_cast<MML::Real>(0.89445)};
		Vector3Spherical v_spher1{static_cast<MML::Real>(6.9487922897236354), static_cast<MML::Real>(0.25555062599984346)	, -static_cast<MML::Real>(0.4)};

		 v_transf_to_spher = CoordTransfCartToSpher.transfVecContravariant(v_cart1, posCart1);
		 REQUIRE(true == v_transf_to_spher.IsEqualTo(v_spher1, TOL(1e-8, 1e-4)));

		 // TODO - pokazati kako se do real dolazi projekcijom na bazne vektore

	}

	TEST_CASE("Test_Covariant_transf_cart_to_spher")
	{
			TEST_PRECISION_INFO();
		Vector3Cartesian p_cart{ static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0) };
		Vector3Spherical p_spher(CoordTransfSpherToCart.transfInverse(p_cart));

		ScalarFunction<3> fPotCart(Fields::InverseRadialPotentialFieldCart);
		Vector3Cartesian  grad_cart = ScalarFieldOperations::GradientCart<3>(fPotCart, p_cart);

		Vector3Spherical grad_transf_to_spher = CoordTransfCartToSpher.transfVecCovariant(grad_cart, p_spher);
		// TODO static_cast<MML::Real>(0.9) - nekad je -1/3 bilo sqrt(2)???
		REQUIRE(true == grad_transf_to_spher.IsEqualTo(Vector3Spherical(-static_cast<MML::Real>(1.0) / 3, static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0)), TOL(1e-6, 1e-3)));

		Vector3Cartesian back_transf_to_cart = CoordTransfSpherToCart.transfVecCovariant(grad_transf_to_spher, p_cart);
		REQUIRE(true == back_transf_to_cart.IsEqualTo(grad_cart, TOL(1e-7, 1e-3)));
	}

} // namespace MML::Tests::Core::CoordTransfVectorTests

#include "../catch/catch.hpp"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "core/CoordTransf.h"
#include "core/Fields.h"
#include "core/FieldOperations.h"
#endif

using namespace MML;

// TODO - HIGH!!! Test coord transf
// make sure it is clear what convention is used

TEST_CASE("Test_CoordTransf_Cartesian_to_Spherical", "[simple]")
{
    Vector3Spherical posSpher;

    posSpher = CoordTransfCartToSpher.transf(Vector3Cartesian{1.0, 0.0, 0.0});
    REQUIRE(posSpher == Vector3Spherical{1.0, Constants::PI / 2, 0.0});

    posSpher = CoordTransfCartToSpher.transf(Vector3Cartesian{-1.0, 0.0, 0.0});
    REQUIRE(posSpher == Vector3Spherical{1.0, Constants::PI / 2, Constants::PI});

    posSpher = CoordTransfCartToSpher.transf(Vector3Cartesian{0.0, 1.0, 0.0});
    REQUIRE(posSpher == Vector3Spherical{1.0, Constants::PI / 2, Constants::PI / 2});
    
    posSpher = CoordTransfCartToSpher.transf(Vector3Cartesian{0.0, -1.0, 0.0});
    REQUIRE(posSpher == Vector3Spherical{1.0, Constants::PI / 2, -Constants::PI / 2});
    
    posSpher = CoordTransfCartToSpher.transf(Vector3Cartesian{0.0, 0.0, 1.0});
    REQUIRE(posSpher == Vector3Spherical{1.0, 0.0, 0.0});
    
    posSpher = CoordTransfCartToSpher.transf(Vector3Cartesian{0.0, 0.0, -1.0});
    REQUIRE(posSpher == Vector3Spherical{1.0, Constants::PI, 0.0});
}

TEST_CASE("Test_CoordTransf_Cartesian_to_Cylindrical", "[simple]")
{
    Vector3Cylindrical posCyl;

    posCyl = CoordTransfCartToCyl.transf(Vector3Cartesian{1.0, 0.0, 0.0});
    REQUIRE(posCyl == Vector3Cylindrical{1.0, 0.0, 0.0});

    posCyl = CoordTransfCartToCyl.transf(Vector3Cartesian{-1.0, 0.0, 0.0});
    REQUIRE(posCyl == Vector3Cylindrical{1.0, Constants::PI, 0.0});

    posCyl = CoordTransfCartToCyl.transf(Vector3Cartesian{0.0, 1.0, 0.0});
    REQUIRE(posCyl == Vector3Cylindrical{1.0, Constants::PI / 2, 0.0});
    
    posCyl = CoordTransfCartToCyl.transf(Vector3Cartesian{0.0, -1.0, 0.0});
    REQUIRE(posCyl == Vector3Cylindrical{1.0, -Constants::PI / 2, 0.0});
    
    posCyl = CoordTransfCartToCyl.transf(Vector3Cartesian{0.0, 0.0, 1.0});
    REQUIRE(posCyl == Vector3Cylindrical{0.0, 0.0, 1.0});
    
    posCyl = CoordTransfCartToCyl.transf(Vector3Cartesian{0.0, 0.0, -1.0});
    REQUIRE(posCyl == Vector3Cylindrical{0.0, 0.0, -1.0});
}

TEST_CASE("Test_CoordTransf_Spherical_to_Cartesian", "[simple]")
{
    Vector3Cartesian posCart;

    posCart  = CoordTransfSpherToCart.transf(Vector3Spherical{1.0, 0.0, 0.0});
    REQUIRE(posCart == Vector3Cartesian{0.0, 0.0, 1.0});
    
    posCart  = CoordTransfSpherToCart.transf(Vector3Spherical{1.0, Constants::PI / 2, 0.0});
    
    REQUIRE(1.0 == posCart.X());
    REQUIRE(0.0 == posCart.Y());
    REQUIRE_THAT(0.0, Catch::Matchers::WithinAbs(posCart.Z(), 1e-16));

    posCart  = CoordTransfSpherToCart.transf(Vector3Spherical{1.0, Constants::PI / 2, Constants::PI/2});
    
    REQUIRE_THAT(0.0, Catch::Matchers::WithinAbs(posCart.X(), 1e-16));
    REQUIRE(1.0 == posCart.Y());
    REQUIRE_THAT(0.0, Catch::Matchers::WithinAbs(posCart.Z(), 1e-16));
}

TEST_CASE("Test_Contravariant_transf_cart_to_spher")
{
    // TODO 0.9 - dodati još 5-6 test caseova
    Vector3Cartesian v_cart{1.0, 1.0, 0.0};
    Vector3Cartesian posCart{1.0, 1.0, 0.0};

    Vector3Spherical v_transf_to_spher = CoordTransfCartToSpher.transfVecContravariant(v_cart, posCart);
    REQUIRE(true == v_transf_to_spher.IsEqual(Vector3Spherical(sqrt(2), 0.0, 0.0), 1e-8));

    Vector3Spherical x1_spher{ CoordTransfCartToSpher.transf(posCart) };
    Vector3Cartesian v_back_transf_to_cart = CoordTransfSpherToCart.transfVecContravariant(v_transf_to_spher, x1_spher);
    REQUIRE(true == v_back_transf_to_cart.IsEqual(v_cart, 1e-7));
}

TEST_CASE("Test_Covariant_transf_cart_to_spher")
{
    Vector3Cartesian p_cart{1.0, 1.0, 1.0};
    Vector3Spherical p_spher(CoordTransfSpherToCart.transfInverse(p_cart));

    ScalarFunction<3> fPotCart(Fields::InverseRadialPotentialFieldCart);    
    Vector3Cartesian  grad_cart = ScalarFieldOperations::GradientCart<3>(fPotCart, p_cart);

    Vector3Spherical grad_transf_to_spher = CoordTransfCartToSpher.transfVecCovariant(grad_cart, p_spher);
    // TODO 0.9 - nekad je -1/3 bilo sqrt(2)???
    REQUIRE(true == grad_transf_to_spher.IsEqual(Vector3Spherical(-1.0/3, 0.0, 0.0), 1e-6));

    Vector3Cartesian back_transf_to_cart = CoordTransfSpherToCart.transfVecCovariant(grad_transf_to_spher, p_cart);
    REQUIRE(true == back_transf_to_cart.IsEqual(grad_cart, 1e-7));
}

TEST_CASE("Test_GetUnitVector")
{
    Vector3Cartesian p1{2.0, 1.0, -2.0};
    auto p1Spher      = CoordTransfCartToSpher.transf(p1);

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
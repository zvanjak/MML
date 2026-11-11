///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        icoord_transf_tests.cpp                                             ///
///  Description: Contract tests for ICoordTransf / ICoordTransfWithInverse after the  ///
///               redesign (MinimalMathLibrary-wewi.6.3): dimensions are deduced from ///
///               the vector types, the base is IFunction<VectorTo, const VectorFrom&>,///
///               and non-square (From != To) transformations are first-class.        ///
///////////////////////////////////////////////////////////////////////////////////////////

#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"

#include <mml/MMLBase.h>
#include <mml/base/Function.h>
#include <mml/base/Vector/VectorTypes2D.h>
#include <mml/base/Vector/VectorTypes3D.h>
#include <mml/interfaces/ICoordTransf.h>
#include <mml/core/CoordTransf/CoordTransfBase.h>
#include <mml/core/CoordTransf/CoordTransfSpherical.h>

#include <type_traits>

using namespace MML;

namespace MML::Tests::Interfaces::ICoordTransfTests {

namespace {
    // Embedding of the unit cylinder: (phi, z) in R^2 -> (cos phi, sin phi, z) in R^3.
    // A genuine From != To transformation: Jacobian is 3 x 2.
    class CylinderEmbedding : public CoordTransf<Vector2Cartesian, Vector3Cartesian>
    {
        static Real fx(const VectorN<Real, 2>& q) { return std::cos(q[0]); }
        static Real fy(const VectorN<Real, 2>& q) { return std::sin(q[0]); }
        static Real fz(const VectorN<Real, 2>& q) { return q[1]; }
        inline static ScalarFunction<2> _f[3] = { ScalarFunction<2>{fx}, ScalarFunction<2>{fy}, ScalarFunction<2>{fz} };
    public:
        Vector3Cartesian transf(const Vector2Cartesian& q) const override { return Vector3Cartesian(fx(q), fy(q), fz(q)); }
        const IScalarFunction<2>& coordTransfFunc(int i) const override { return _f[i]; }
    };
}

TEST_CASE("ICoordTransf deduces FromDim/ToDim from the vector types and derives from IFunction<To, From>", "[interfaces][coord-transf][wewi.6.3]")
{
    using T = ICoordTransf<Vector2Cartesian, Vector3Cartesian>;
    STATIC_REQUIRE(T::FromDim == 2);
    STATIC_REQUIRE(T::ToDim == 3);
    STATIC_REQUIRE(std::is_base_of_v<IFunction<Vector3Cartesian, const Vector2Cartesian&>, T>);
    // the old base claimed the map was N -> N; it must be gone
    STATIC_REQUIRE_FALSE(std::is_base_of_v<IVectorFunction<2>, T>);
    STATIC_REQUIRE_FALSE(std::is_base_of_v<IVectorFunction<3>, T>);

    // the legacy explicit-N spelling still compiles when N == FromDim
    STATIC_REQUIRE(std::is_same_v<ICoordTransf<Vector3Spherical, Vector3Cartesian, 3>, ICoordTransf<Vector3Spherical, Vector3Cartesian>>);
    STATIC_REQUIRE(Detail::VectorDimension<Vector2Polar> == 2);
    STATIC_REQUIRE(Detail::VectorDimension<Vector3Spherical> == 3);
    STATIC_REQUIRE(Detail::VectorDimension<VectorN<Real, 7>> == 7);
}

TEST_CASE("Non-square transformation: Jacobian is ToDim x FromDim and vector transforms honour both dims", "[interfaces][coord-transf][wewi.6.3]")
{
    CylinderEmbedding emb;
    STATIC_REQUIRE(std::is_same_v<CylinderEmbedding::JacobianType, MatrixNM<Real, 3, 2>>);

    Vector2Cartesian q(Constants::PI / 3, static_cast<MML::Real>(2.0));
    Vector3Cartesian p = emb.transf(q);
    REQUIRE(p[0] == Catch::Approx(0.5));
    REQUIRE(p[1] == Catch::Approx(std::sqrt(3.0) / 2));
    REQUIRE(p[2] == Catch::Approx(2.0));

    // IFunction call form is transf(); a plain VectorN<2> is accepted via the converting ctor
    Vector3Cartesian p2 = emb(q);
    REQUIRE(p2.IsEqualTo(p));
    Vector3Cartesian p3 = emb(VectorN<Real, 2>{Constants::PI / 3, static_cast<MML::Real>(2.0)});
    REQUIRE(p3.IsEqualTo(p));
    const IFunction<Vector3Cartesian, const Vector2Cartesian&>& asFunction = emb;
    REQUIRE(asFunction(q).IsEqualTo(p));

    // J = [[-sin, 0], [cos, 0], [0, 1]]
    auto J = emb.jacobian(q);
    REQUIRE(J(0, 0) == Catch::Approx(-std::sin(Constants::PI / 3)).margin(1e-8));
    REQUIRE(J(1, 0) == Catch::Approx(std::cos(Constants::PI / 3)).margin(1e-8));
    REQUIRE(J(2, 0) == Catch::Approx(0.0).margin(1e-8));
    REQUIRE(J(0, 1) == Catch::Approx(0.0).margin(1e-8));
    REQUIRE(J(1, 1) == Catch::Approx(0.0).margin(1e-8));
    REQUIRE(J(2, 1) == Catch::Approx(1.0).margin(1e-8));

    // tangent (contravariant) push-forward of d/dphi is the circle tangent
    Vector3Cartesian t = emb.transfVecContravariant(Vector2Cartesian(1.0, 0.0), q);
    REQUIRE(t[0] == Catch::Approx(-std::sin(Constants::PI / 3)).margin(1e-8));
    REQUIRE(t[1] == Catch::Approx(std::cos(Constants::PI / 3)).margin(1e-8));
    REQUIRE(t[2] == Catch::Approx(0.0).margin(1e-8));

    // basis vectors: e_phi (column 0) and e_z (column 1) of the 3 x 2 Jacobian
    REQUIRE(emb.getBasisVec(0, q).IsEqualTo(t, 1e-8));
    Vector3Cartesian ez = emb.getBasisVec(1, q);
    REQUIRE(ez.IsEqualTo(Vector3Cartesian(0, 0, 1), 1e-8));

    // pull-back of a covector (gradient of z in R^3) onto the chart is (0, 1)
    Vector2Cartesian w = emb.transfInverseVecCovariant(Vector3Cartesian(0, 0, 1), q);
    REQUIRE(w[0] == Catch::Approx(0.0).margin(1e-8));
    REQUIRE(w[1] == Catch::Approx(1.0).margin(1e-8));
}

TEST_CASE("Square transformations keep their behaviour and types", "[interfaces][coord-transf][wewi.6.3]")
{
    CoordTransfSphericalToCartesian sph2cart;
    using S = CoordTransfSphericalToCartesian;
    STATIC_REQUIRE(S::FromDim == 3);
    STATIC_REQUIRE(S::ToDim == 3);
    STATIC_REQUIRE(std::is_same_v<S::JacobianType, MatrixNM<Real, 3, 3>>);
    STATIC_REQUIRE(std::is_same_v<S::InverseJacobianType, MatrixNM<Real, 3, 3>>);

    Vector3Spherical s(2.0, Constants::PI / 2, 0.0);
    Vector3Cartesian c = sph2cart.transf(s);
    REQUIRE(c[0] == Catch::Approx(2.0).margin(1e-12));
    REQUIRE(c[1] == Catch::Approx(0.0).margin(1e-12));
    REQUIRE(c[2] == Catch::Approx(0.0).margin(1e-12));
    REQUIRE(sph2cart(s).IsEqualTo(c, 1e-12));                       // IFunction form agrees with transf()
    REQUIRE(sph2cart.transfInverse(c).IsEqualTo(s, 1e-10));

    // inverse Jacobian is FromDim x ToDim and K * J = I at a regular point
    Vector3Spherical s2(1.5, 1.0, 0.7);
    auto J = sph2cart.jacobian(s2);
    auto K = sph2cart.inverseJacobian(sph2cart.transf(s2));
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            Real sum = 0;
            for (int k = 0; k < 3; ++k) sum += K(i, k) * J(k, j);
            REQUIRE(sum == Catch::Approx(i == j ? 1.0 : 0.0).margin(1e-6));
        }
}

} // namespace MML::Tests::Interfaces::ICoordTransfTests

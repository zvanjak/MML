///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML) Tests                            ///
///                                                                                   ///
///  File:        vector_types_tests.cpp                                              ///
///  Description: Comprehensive tests for dimensional vector types                    ///
///               Tests all vector type classes: Vec2Cart, Vec2Pol, Vec3Cart,         ///
///               Vec3Sph, Vec3Cyl, Vec4Mink                                          ///
///////////////////////////////////////////////////////////////////////////////////////////

#include <catch2/catch_all.hpp>
#include "../../TestPrecision.h"
#include "../../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/base/Vector/VectorTypes2D.h>
#include <mml/base/Vector/VectorTypes3D.h>
#include <mml/base/Vector/VectorTypes4D.h>
#include <mml/base/Geometry/Geometry.h>
#endif

#include <cmath>

using namespace MML;
using namespace MML::Testing;
using namespace MML::Testing::Matchers;

namespace MML::Tests::Base::VectorTypesTests
{

///////////////////////////////////////////////////////////////////////////////////////////
///                         Vector2Cartesian Tests                                     ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Vector2Cartesian - Constructors", "[VectorTypes][Vector2Cartesian]")
{
    TEST_PRECISION_INFO();

    SECTION("Default constructor creates zero vector")
    {
        Vector2Cartesian v;
        REQUIRE(v.X() == static_cast<MML::Real>(0.0));
        REQUIRE(v.Y() == static_cast<MML::Real>(0.0));
    }

    SECTION("Component constructor")
    {
        Vector2Cartesian v(static_cast<MML::Real>(3.0), static_cast<MML::Real>(4.0));
        REQUIRE(v.X() == static_cast<MML::Real>(3.0));
        REQUIRE(v.Y() == static_cast<MML::Real>(4.0));
    }

    SECTION("Construct from VectorN")
    {
        VectorN<Real, 2> base{static_cast<MML::Real>(5.0), static_cast<MML::Real>(6.0)};
        Vector2Cartesian v(base);
        REQUIRE(v.X() == static_cast<MML::Real>(5.0));
        REQUIRE(v.Y() == static_cast<MML::Real>(6.0));
    }

    SECTION("Construct from two points")
    {
        Point2Cartesian a(static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0));
        Point2Cartesian b(static_cast<MML::Real>(4.0), static_cast<MML::Real>(6.0));
        Vector2Cartesian v(a, b);
        REQUIRE(v.X() == static_cast<MML::Real>(3.0));
        REQUIRE(v.Y() == static_cast<MML::Real>(4.0));
    }

    SECTION("Initializer list constructor")
    {
        Vector2Cartesian v{static_cast<MML::Real>(7.0), static_cast<MML::Real>(8.0)};
        REQUIRE(v.X() == static_cast<MML::Real>(7.0));
        REQUIRE(v.Y() == static_cast<MML::Real>(8.0));
    }
}

TEST_CASE("Vector2Cartesian - Arithmetic operations", "[VectorTypes][Vector2Cartesian]")
{
    TEST_PRECISION_INFO();
    Vector2Cartesian a(static_cast<MML::Real>(3.0), static_cast<MML::Real>(4.0));
    Vector2Cartesian b(static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0));

    SECTION("Vector addition")
    {
        Vector2Cartesian c = a + b;
        REQUIRE(c.X() == static_cast<MML::Real>(4.0));
        REQUIRE(c.Y() == static_cast<MML::Real>(6.0));
    }

    SECTION("Vector subtraction")
    {
        Vector2Cartesian c = a - b;
        REQUIRE(c.X() == static_cast<MML::Real>(2.0));
        REQUIRE(c.Y() == static_cast<MML::Real>(2.0));
    }

    SECTION("Unary minus (negation)")
    {
        Vector2Cartesian c = -a;
        REQUIRE(c.X() == static_cast<MML::Real>(-3.0));
        REQUIRE(c.Y() == static_cast<MML::Real>(-4.0));
    }

    SECTION("Scalar multiplication (vector * scalar)")
    {
        Vector2Cartesian c = a * static_cast<MML::Real>(2.0);
        REQUIRE(c.X() == static_cast<MML::Real>(6.0));
        REQUIRE(c.Y() == static_cast<MML::Real>(8.0));
    }

    SECTION("Scalar multiplication (scalar * vector)")
    {
        Vector2Cartesian c = static_cast<MML::Real>(2.0) * a;
        REQUIRE(c.X() == static_cast<MML::Real>(6.0));
        REQUIRE(c.Y() == static_cast<MML::Real>(8.0));
    }

    SECTION("Scalar division")
    {
        Vector2Cartesian c = a / static_cast<MML::Real>(2.0);
        REQUIRE(c.X() == static_cast<MML::Real>(1.5));
        REQUIRE(c.Y() == static_cast<MML::Real>(2.0));
    }

    SECTION("Scalar product (dot product)")
    {
        Real dot = a * b;
        REQUIRE(dot == static_cast<MML::Real>(11.0));  // 3*1 + 4*2 = 11
    }

    SECTION("ScalarProduct friend function")
    {
        Real dot = ScalarProduct(a, b);
        REQUIRE(dot == static_cast<MML::Real>(11.0));
    }
}

TEST_CASE("Vector2Cartesian - Comparison operations", "[VectorTypes][Vector2Cartesian]")
{
    TEST_PRECISION_INFO();
    Vector2Cartesian a(static_cast<MML::Real>(3.0), static_cast<MML::Real>(4.0));
    Vector2Cartesian b(static_cast<MML::Real>(3.0), static_cast<MML::Real>(4.0));
    Vector2Cartesian c(static_cast<MML::Real>(3.1), static_cast<MML::Real>(4.0));

    SECTION("Equality operator")
    {
        REQUIRE(a == b);
        REQUIRE_FALSE(a == c);
    }

    SECTION("Inequality operator")
    {
        REQUIRE_FALSE(a != b);
        REQUIRE(a != c);
    }

    SECTION("IsEqualTo with tolerance")
    {
        REQUIRE(a.IsEqualTo(b, TOL(1e-10, 1e-5)));
        REQUIRE(a.IsEqualTo(c, static_cast<MML::Real>(0.2)));  // Within tolerance
        REQUIRE_FALSE(a.IsEqualTo(c, static_cast<MML::Real>(0.05)));  // Outside tolerance
    }
}

TEST_CASE("Vector2Cartesian - Geometric operations", "[VectorTypes][Vector2Cartesian]")
{
    TEST_PRECISION_INFO();

    SECTION("NormL2 (magnitude)")
    {
        Vector2Cartesian v(static_cast<MML::Real>(3.0), static_cast<MML::Real>(4.0));
        REQUIRE(v.NormL2() == static_cast<MML::Real>(5.0));
    }

    SECTION("Normalized returns unit vector")
    {
        Vector2Cartesian v(static_cast<MML::Real>(3.0), static_cast<MML::Real>(4.0));
        Vector2Cartesian unit = v.Normalized();
        REQUIRE_THAT(unit.X(), RealWithinRel(static_cast<MML::Real>(0.6), TOL(1e-10, 1e-5)));
        REQUIRE_THAT(unit.Y(), RealWithinRel(static_cast<MML::Real>(0.8), TOL(1e-10, 1e-5)));
        REQUIRE_THAT(unit.NormL2(), RealWithinRel(static_cast<MML::Real>(1.0), TOL(1e-10, 1e-5)));
    }

    SECTION("Normalized of zero vector throws")
    {
        Vector2Cartesian v(static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0));
        REQUIRE_THROWS_AS(v.Normalized(), VectorDimensionError);
    }

    SECTION("GetAsUnitVector")
    {
        Vector2Cartesian v(static_cast<MML::Real>(0.0), static_cast<MML::Real>(5.0));
        Vector2Cartesian unit = v.Normalized();
        REQUIRE(unit.X() == static_cast<MML::Real>(0.0));
        REQUIRE_THAT(unit.Y(), RealWithinRel(static_cast<MML::Real>(1.0), TOL(1e-10, 1e-5)));
    }

    SECTION("getPerpendicularVectors")
    {
        Vector2Cartesian v(static_cast<MML::Real>(3.0), static_cast<MML::Real>(4.0));
        Vector2Cartesian v1, v2;
        v.GetPerpendicularVectors(v1, v2);
        
        // v1 and v2 should be perpendicular to v
        REQUIRE_THAT((v * v1), RealApprox(static_cast<MML::Real>(0.0)).margin(Tolerance::Loose));
        REQUIRE_THAT((v * v2), RealApprox(static_cast<MML::Real>(0.0)).margin(Tolerance::Loose));
        
        // v1 and v2 should be opposite
        REQUIRE_THAT(v1.X(), RealWithinRel(-v2.X(), TOL(1e-10, 1e-5)));
        REQUIRE_THAT(v1.Y(), RealWithinRel(-v2.Y(), TOL(1e-10, 1e-5)));
    }
}

TEST_CASE("Vector2Cartesian - Point operations", "[VectorTypes][Vector2Cartesian]")
{
    TEST_PRECISION_INFO();
    Point2Cartesian p(static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0));
    Vector2Cartesian v(static_cast<MML::Real>(3.0), static_cast<MML::Real>(4.0));

    SECTION("Point + Vector")
    {
        Point2Cartesian result = p + v;
        REQUIRE(result.X() == static_cast<MML::Real>(4.0));
        REQUIRE(result.Y() == static_cast<MML::Real>(6.0));
    }

    SECTION("Point - Vector")
    {
        Point2Cartesian result = p - v;
        REQUIRE(result.X() == static_cast<MML::Real>(-2.0));
        REQUIRE(result.Y() == static_cast<MML::Real>(-2.0));
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         Vector2Polar Tests                                         ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Vector2Polar - Constructors", "[VectorTypes][Vector2Polar]")
{
    TEST_PRECISION_INFO();

    SECTION("Default constructor")
    {
        Vector2Polar v;
        REQUIRE(v.R() == static_cast<MML::Real>(0.0));
        REQUIRE(v.Phi() == static_cast<MML::Real>(0.0));
    }

    SECTION("Component constructor")
    {
        Vector2Polar v(static_cast<MML::Real>(5.0), Constants::PI / 4);
        REQUIRE(v.R() == static_cast<MML::Real>(5.0));
        REQUIRE_THAT(v.Phi(), RealWithinRel(Constants::PI / 4, TOL(1e-10, 1e-5)));
    }
}

TEST_CASE("Vector2Polar - Arithmetic operations", "[VectorTypes][Vector2Polar]")
{
    TEST_PRECISION_INFO();

    SECTION("Unary minus adds π to angle")
    {
        Vector2Polar v(static_cast<MML::Real>(5.0), static_cast<MML::Real>(0.0));
        Vector2Polar neg = -v;
        REQUIRE(neg.R() == static_cast<MML::Real>(5.0));
        REQUIRE_THAT(neg.Phi(), RealWithinRel(Constants::PI, TOL(1e-10, 1e-5)));
    }

    SECTION("Scalar multiplication scales radius only")
    {
        Vector2Polar v(static_cast<MML::Real>(5.0), Constants::PI / 4);
        Vector2Polar scaled = v * static_cast<MML::Real>(2.0);
        REQUIRE(scaled.R() == static_cast<MML::Real>(10.0));
        REQUIRE_THAT(scaled.Phi(), RealWithinRel(Constants::PI / 4, TOL(1e-10, 1e-5)));
    }

    SECTION("Scalar multiplication (scalar * vector)")
    {
        Vector2Polar v(static_cast<MML::Real>(5.0), Constants::PI / 4);
        Vector2Polar scaled = static_cast<MML::Real>(3.0) * v;
        REQUIRE(scaled.R() == static_cast<MML::Real>(15.0));
        REQUIRE_THAT(scaled.Phi(), RealWithinRel(Constants::PI / 4, TOL(1e-10, 1e-5)));
    }

    SECTION("Scalar division divides radius only")
    {
        Vector2Polar v(static_cast<MML::Real>(10.0), Constants::PI / 3);
        Vector2Polar divided = v / static_cast<MML::Real>(2.0);
        REQUIRE(divided.R() == static_cast<MML::Real>(5.0));
        REQUIRE_THAT(divided.Phi(), RealWithinRel(Constants::PI / 3, TOL(1e-10, 1e-5)));
    }

    SECTION("Vector addition - same direction")
    {
        Vector2Polar v1(static_cast<MML::Real>(3.0), static_cast<MML::Real>(0.0));
        Vector2Polar v2(static_cast<MML::Real>(4.0), static_cast<MML::Real>(0.0));
        Vector2Polar sum = v1 + v2;
        REQUIRE_THAT(sum.R(), RealWithinRel(static_cast<MML::Real>(7.0), TOL(1e-10, 1e-5)));
    }

    SECTION("Vector addition - opposite directions")
    {
        Vector2Polar v1(static_cast<MML::Real>(5.0), static_cast<MML::Real>(0.0));
        Vector2Polar v2(static_cast<MML::Real>(3.0), Constants::PI);
        Vector2Polar sum = v1 + v2;
        REQUIRE_THAT(sum.R(), RealWithinRel(static_cast<MML::Real>(2.0), TOL(1e-10, 1e-5)));
    }

    SECTION("Vector addition - perpendicular")
    {
        Vector2Polar v1(static_cast<MML::Real>(3.0), static_cast<MML::Real>(0.0));
        Vector2Polar v2(static_cast<MML::Real>(4.0), Constants::PI / 2);
        Vector2Polar sum = v1 + v2;
        REQUIRE_THAT(sum.R(), RealWithinRel(static_cast<MML::Real>(5.0), TOL(1e-10, 1e-5)));  // 3-4-5 triangle
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         Vector3Cartesian Tests                                     ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Vector3Cartesian - Constructors", "[VectorTypes][Vector3Cartesian]")
{
    TEST_PRECISION_INFO();

    SECTION("Default constructor creates zero vector")
    {
        Vector3Cartesian v;
        REQUIRE(v.X() == static_cast<MML::Real>(0.0));
        REQUIRE(v.Y() == static_cast<MML::Real>(0.0));
        REQUIRE(v.Z() == static_cast<MML::Real>(0.0));
    }

    SECTION("Component constructor")
    {
        Vector3Cartesian v(static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0));
        REQUIRE(v.X() == static_cast<MML::Real>(1.0));
        REQUIRE(v.Y() == static_cast<MML::Real>(2.0));
        REQUIRE(v.Z() == static_cast<MML::Real>(3.0));
    }

    SECTION("Construct from VectorN")
    {
        VectorN<Real, 3> base{static_cast<MML::Real>(4.0), static_cast<MML::Real>(5.0), static_cast<MML::Real>(6.0)};
        Vector3Cartesian v(base);
        REQUIRE(v.X() == static_cast<MML::Real>(4.0));
        REQUIRE(v.Y() == static_cast<MML::Real>(5.0));
        REQUIRE(v.Z() == static_cast<MML::Real>(6.0));
    }

    SECTION("Construct from two points")
    {
        Point3Cartesian a(static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0));
        Point3Cartesian b(static_cast<MML::Real>(4.0), static_cast<MML::Real>(6.0), static_cast<MML::Real>(8.0));
        Vector3Cartesian v(a, b);
        REQUIRE(v.X() == static_cast<MML::Real>(3.0));
        REQUIRE(v.Y() == static_cast<MML::Real>(4.0));
        REQUIRE(v.Z() == static_cast<MML::Real>(5.0));
    }

    SECTION("Construct from Point3Cartesian")
    {
        Point3Cartesian p(static_cast<MML::Real>(7.0), static_cast<MML::Real>(8.0), static_cast<MML::Real>(9.0));
        Vector3Cartesian v(p);
        REQUIRE(v.X() == static_cast<MML::Real>(7.0));
        REQUIRE(v.Y() == static_cast<MML::Real>(8.0));
        REQUIRE(v.Z() == static_cast<MML::Real>(9.0));
    }

    SECTION("Initializer list constructor")
    {
        Vector3Cartesian v{static_cast<MML::Real>(10.0), static_cast<MML::Real>(11.0), static_cast<MML::Real>(12.0)};
        REQUIRE(v.X() == static_cast<MML::Real>(10.0));
        REQUIRE(v.Y() == static_cast<MML::Real>(11.0));
        REQUIRE(v.Z() == static_cast<MML::Real>(12.0));
    }
}

TEST_CASE("Vector3Cartesian - Arithmetic operations", "[VectorTypes][Vector3Cartesian]")
{
    TEST_PRECISION_INFO();
    Vector3Cartesian a(static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0));
    Vector3Cartesian b(static_cast<MML::Real>(4.0), static_cast<MML::Real>(5.0), static_cast<MML::Real>(6.0));

    SECTION("Vector addition")
    {
        Vector3Cartesian c = a + b;
        REQUIRE(c.X() == static_cast<MML::Real>(5.0));
        REQUIRE(c.Y() == static_cast<MML::Real>(7.0));
        REQUIRE(c.Z() == static_cast<MML::Real>(9.0));
    }

    SECTION("Vector subtraction")
    {
        Vector3Cartesian c = b - a;
        REQUIRE(c.X() == static_cast<MML::Real>(3.0));
        REQUIRE(c.Y() == static_cast<MML::Real>(3.0));
        REQUIRE(c.Z() == static_cast<MML::Real>(3.0));
    }

    SECTION("Unary minus")
    {
        Vector3Cartesian c = -a;
        REQUIRE(c.X() == static_cast<MML::Real>(-1.0));
        REQUIRE(c.Y() == static_cast<MML::Real>(-2.0));
        REQUIRE(c.Z() == static_cast<MML::Real>(-3.0));
    }

    SECTION("Scalar multiplication")
    {
        Vector3Cartesian c = a * static_cast<MML::Real>(2.0);
        REQUIRE(c.X() == static_cast<MML::Real>(2.0));
        REQUIRE(c.Y() == static_cast<MML::Real>(4.0));
        REQUIRE(c.Z() == static_cast<MML::Real>(6.0));
    }

    SECTION("Scalar division")
    {
        Vector3Cartesian c = b / static_cast<MML::Real>(2.0);
        REQUIRE(c.X() == static_cast<MML::Real>(2.0));
        REQUIRE(c.Y() == static_cast<MML::Real>(2.5));
        REQUIRE(c.Z() == static_cast<MML::Real>(3.0));
    }

    SECTION("Scalar product (dot product)")
    {
        Real dot = a * b;
        REQUIRE(dot == static_cast<MML::Real>(32.0));  // 1*4 + 2*5 + 3*6 = 32
    }

    SECTION("ScalarProduct member function")
    {
        REQUIRE(a.ScalarProduct(b) == static_cast<MML::Real>(32.0));
    }

    SECTION("ScalarProduct friend function")
    {
        REQUIRE(ScalarProduct(a, b) == static_cast<MML::Real>(32.0));
    }
}

TEST_CASE("Vector3Cartesian - Cross product", "[VectorTypes][Vector3Cartesian]")
{
    TEST_PRECISION_INFO();

    SECTION("Cross product of unit vectors")
    {
        Vector3Cartesian i(static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0));
        Vector3Cartesian j(static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0));
        Vector3Cartesian k(static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0));

        // i × j = k
        Vector3Cartesian ixj = VectorProduct(i, j);
        REQUIRE(ixj.X() == static_cast<MML::Real>(0.0));
        REQUIRE(ixj.Y() == static_cast<MML::Real>(0.0));
        REQUIRE(ixj.Z() == static_cast<MML::Real>(1.0));

        // j × k = i
        Vector3Cartesian jxk = VectorProduct(j, k);
        REQUIRE(jxk.X() == static_cast<MML::Real>(1.0));
        REQUIRE(jxk.Y() == static_cast<MML::Real>(0.0));
        REQUIRE(jxk.Z() == static_cast<MML::Real>(0.0));

        // k × i = j
        Vector3Cartesian kxi = VectorProduct(k, i);
        REQUIRE(kxi.X() == static_cast<MML::Real>(0.0));
        REQUIRE(kxi.Y() == static_cast<MML::Real>(1.0));
        REQUIRE(kxi.Z() == static_cast<MML::Real>(0.0));
    }

    SECTION("Cross product anti-commutativity")
    {
        Vector3Cartesian a(static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0));
        Vector3Cartesian b(static_cast<MML::Real>(4.0), static_cast<MML::Real>(5.0), static_cast<MML::Real>(6.0));

        Vector3Cartesian axb = VectorProduct(a, b);
        Vector3Cartesian bxa = VectorProduct(b, a);

        REQUIRE_THAT(axb.X(), RealWithinRel(-bxa.X(), TOL(1e-10, 1e-5)));
        REQUIRE_THAT(axb.Y(), RealWithinRel(-bxa.Y(), TOL(1e-10, 1e-5)));
        REQUIRE_THAT(axb.Z(), RealWithinRel(-bxa.Z(), TOL(1e-10, 1e-5)));
    }

    SECTION("Cross product of parallel vectors is zero")
    {
        Vector3Cartesian a(static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0));
        Vector3Cartesian b(static_cast<MML::Real>(2.0), static_cast<MML::Real>(4.0), static_cast<MML::Real>(6.0));

        Vector3Cartesian cross = VectorProduct(a, b);
        REQUIRE_THAT(cross.NormL2(), RealApprox(static_cast<MML::Real>(0.0)).margin(Tolerance::Loose));
    }

    SECTION("Cross product perpendicular to both inputs")
    {
        Vector3Cartesian a(static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0));
        Vector3Cartesian b(static_cast<MML::Real>(4.0), static_cast<MML::Real>(5.0), static_cast<MML::Real>(6.0));

        Vector3Cartesian cross = VectorProduct(a, b);
        REQUIRE_THAT(ScalarProduct(cross, a), RealApprox(static_cast<MML::Real>(0.0)).margin(Tolerance::Loose));
        REQUIRE_THAT(ScalarProduct(cross, b), RealApprox(static_cast<MML::Real>(0.0)).margin(Tolerance::Loose));
    }
}

TEST_CASE("Vector3Cartesian - Geometric operations", "[VectorTypes][Vector3Cartesian]")
{
    TEST_PRECISION_INFO();

    SECTION("NormL2")
    {
        Vector3Cartesian v(static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0), static_cast<MML::Real>(6.0));
        if constexpr (std::is_same_v<Real, float>)
            REQUIRE_THAT(v.NormL2(), RealWithinRelNoScale(Real(7), Real(4) * Precision::Eps));
        else
            REQUIRE(v.NormL2() == Real(7));
    }

    SECTION("Normalized")
    {
        Vector3Cartesian v(static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(5.0));
        Vector3Cartesian unit = v.Normalized();
        REQUIRE(unit.X() == static_cast<MML::Real>(0.0));
        REQUIRE(unit.Y() == static_cast<MML::Real>(0.0));
        REQUIRE_THAT(unit.Z(), RealWithinRel(static_cast<MML::Real>(1.0), TOL(1e-10, 1e-5)));
    }

    SECTION("Normalized of zero vector throws")
    {
        Vector3Cartesian v(static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0));
        REQUIRE_THROWS_AS(v.Normalized(), VectorDimensionError);
    }

    SECTION("GetAsUnitVector")
    {
        Vector3Cartesian v(static_cast<MML::Real>(3.0), static_cast<MML::Real>(4.0), static_cast<MML::Real>(0.0));
        Vector3Cartesian unit = v.Normalized();
        REQUIRE_THAT(unit.NormL2(), RealWithinRel(static_cast<MML::Real>(1.0), TOL(1e-10, 1e-5)));
        REQUIRE_THAT(unit.X(), RealWithinRel(static_cast<MML::Real>(0.6), TOL(1e-10, 1e-5)));
        REQUIRE_THAT(unit.Y(), RealWithinRel(static_cast<MML::Real>(0.8), TOL(1e-10, 1e-5)));
    }

    SECTION("GetAsUnitVector of zero vector throws; OrZero variant returns zero")
    {
        Vector3Cartesian v(static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0));
        REQUIRE_THROWS_AS(v.Normalized(), VectorDimensionError);

        Vector3Cartesian unit = v.NormalizedOrZero();
        REQUIRE(unit.X() == static_cast<MML::Real>(0.0));
        REQUIRE(unit.Y() == static_cast<MML::Real>(0.0));
        REQUIRE(unit.Z() == static_cast<MML::Real>(0.0));
    }

    SECTION("IsParallelTo")
    {
        Vector3Cartesian a(static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0));
        Vector3Cartesian b(static_cast<MML::Real>(2.0), static_cast<MML::Real>(4.0), static_cast<MML::Real>(6.0));  // Parallel
        Vector3Cartesian c(static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0));  // Not parallel

        REQUIRE(a.IsParallelTo(b));
        REQUIRE_FALSE(a.IsParallelTo(c));
    }

    SECTION("IsPerpendicularTo")
    {
        Vector3Cartesian a(static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0));
        Vector3Cartesian b(static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0));  // Perpendicular
        Vector3Cartesian c(static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0));  // Not perpendicular

        REQUIRE(a.IsPerpendicularTo(b));
        REQUIRE_FALSE(a.IsPerpendicularTo(c));
    }

    SECTION("AngleToVector")
    {
        Vector3Cartesian a(static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0));
        Vector3Cartesian b(static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0));
        Vector3Cartesian c(static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0));
        Vector3Cartesian d(static_cast<MML::Real>(-1.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0));

        REQUIRE_THAT(a.AngleToVector(b), RealWithinRel(Constants::PI / 2, TOL(1e-10, 1e-5)));  // 90°
        REQUIRE_THAT(a.AngleToVector(c), RealApprox(static_cast<MML::Real>(0.0)).margin(Tolerance::Loose));  // 0°
        REQUIRE_THAT(a.AngleToVector(d), RealWithinRel(Constants::PI, TOL(1e-10, 1e-5)));  // 180°
    }

    SECTION("AngleToVector - zero vector throws")
    {
        Vector3Cartesian a(static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0));
        Vector3Cartesian zero(static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0));

        REQUIRE_THROWS_AS(a.AngleToVector(zero), DivisionByZeroError);
        REQUIRE_THROWS_AS(zero.AngleToVector(a), DivisionByZeroError);
        REQUIRE_THROWS_AS(zero.AngleToVector(zero), DivisionByZeroError);
    }

    SECTION("GetPerpendicularVectors")
    {
        Vector3Cartesian v(static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0));
        Vector3Cartesian v1, v2;

        bool success = v.GetPerpendicularVectors(v1, v2);
        REQUIRE(success);

        // Both should be perpendicular to v
        REQUIRE_THAT(v.ScalarProduct(v1), RealApprox(static_cast<MML::Real>(0.0)).margin(Tolerance::Loose));
        REQUIRE_THAT(v.ScalarProduct(v2), RealApprox(static_cast<MML::Real>(0.0)).margin(Tolerance::Loose));

        // v1 and v2 should be perpendicular to each other
        REQUIRE_THAT(v1.ScalarProduct(v2), RealApprox(static_cast<MML::Real>(0.0)).margin(Tolerance::Loose));

        // Both should be unit vectors
        REQUIRE_THAT(v1.NormL2(), RealWithinRel(static_cast<MML::Real>(1.0), TOL(1e-10, 1e-5)));
        REQUIRE_THAT(v2.NormL2(), RealWithinRel(static_cast<MML::Real>(1.0), TOL(1e-10, 1e-5)));
    }

    SECTION("GetPerpendicularVectors of zero vector fails")
    {
        Vector3Cartesian v(static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0));
        Vector3Cartesian v1, v2;

        bool success = v.GetPerpendicularVectors(v1, v2);
        REQUIRE_FALSE(success);
    }
}

TEST_CASE("Vector3Cartesian - Point operations", "[VectorTypes][Vector3Cartesian]")
{
    TEST_PRECISION_INFO();
    Point3Cartesian p(static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0));
    Vector3Cartesian v(static_cast<MML::Real>(4.0), static_cast<MML::Real>(5.0), static_cast<MML::Real>(6.0));

    SECTION("Point + Vector")
    {
        Point3Cartesian result = p + v;
        REQUIRE(result.X() == static_cast<MML::Real>(5.0));
        REQUIRE(result.Y() == static_cast<MML::Real>(7.0));
        REQUIRE(result.Z() == static_cast<MML::Real>(9.0));
    }

    SECTION("Point - Vector")
    {
        Point3Cartesian result = p - v;
        REQUIRE(result.X() == static_cast<MML::Real>(-3.0));
        REQUIRE(result.Y() == static_cast<MML::Real>(-3.0));
        REQUIRE(result.Z() == static_cast<MML::Real>(-3.0));
    }

    SECTION("getAsPoint")
    {
        Point3Cartesian pt = v.ToPoint();
        REQUIRE(pt.X() == static_cast<MML::Real>(4.0));
        REQUIRE(pt.Y() == static_cast<MML::Real>(5.0));
        REQUIRE(pt.Z() == static_cast<MML::Real>(6.0));
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         Vector3Spherical Tests                                     ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Vector3Spherical - Constructors", "[VectorTypes][Vector3Spherical]")
{
    TEST_PRECISION_INFO();

    SECTION("Default constructor creates zero vector")
    {
        Vector3Spherical v;
        REQUIRE(v.R() == static_cast<MML::Real>(0.0));
        REQUIRE(v.Theta() == static_cast<MML::Real>(0.0));
        REQUIRE(v.Phi() == static_cast<MML::Real>(0.0));
    }

    SECTION("Component constructor")
    {
        Vector3Spherical v(static_cast<MML::Real>(5.0), Constants::PI / 4, Constants::PI / 3);
        REQUIRE(v.R() == static_cast<MML::Real>(5.0));
        REQUIRE_THAT(v.Theta(), RealWithinRel(Constants::PI / 4, TOL(1e-10, 1e-5)));
        REQUIRE_THAT(v.Phi(), RealWithinRel(Constants::PI / 3, TOL(1e-10, 1e-5)));
    }
}

TEST_CASE("Vector3Spherical - Arithmetic operations", "[VectorTypes][Vector3Spherical]")
{
    TEST_PRECISION_INFO();

    SECTION("Scalar multiplication scales radius only")
    {
        Vector3Spherical v(static_cast<MML::Real>(5.0), Constants::PI / 4, Constants::PI / 3);
        Vector3Spherical scaled = v * static_cast<MML::Real>(2.0);
        REQUIRE(scaled.R() == static_cast<MML::Real>(10.0));
        REQUIRE_THAT(scaled.Theta(), RealWithinRel(Constants::PI / 4, TOL(1e-10, 1e-5)));
        REQUIRE_THAT(scaled.Phi(), RealWithinRel(Constants::PI / 3, TOL(1e-10, 1e-5)));
    }

    SECTION("Scalar division")
    {
        Vector3Spherical v(static_cast<MML::Real>(10.0), Constants::PI / 4, Constants::PI / 3);
        Vector3Spherical divided = v / static_cast<MML::Real>(2.0);
        REQUIRE(divided.R() == static_cast<MML::Real>(5.0));
        REQUIRE_THAT(divided.Theta(), RealWithinRel(Constants::PI / 4, TOL(1e-10, 1e-5)));
        REQUIRE_THAT(divided.Phi(), RealWithinRel(Constants::PI / 3, TOL(1e-10, 1e-5)));
    }

    SECTION("Division by zero throws")
    {
        Vector3Spherical v(static_cast<MML::Real>(10.0), Constants::PI / 4, Constants::PI / 3);
        REQUIRE_THROWS_AS(v / static_cast<MML::Real>(0.0), DivisionByZeroError);
    }

    SECTION("Vector addition - same direction")
    {
        Vector3Spherical v1(static_cast<MML::Real>(3.0), Constants::PI / 2, static_cast<MML::Real>(0.0));  // Along +x axis
        Vector3Spherical v2(static_cast<MML::Real>(4.0), Constants::PI / 2, static_cast<MML::Real>(0.0));
        Vector3Spherical sum = v1 + v2;
        REQUIRE_THAT(sum.R(), RealWithinRel(static_cast<MML::Real>(7.0), TOL(1e-10, 1e-5)));
    }
}

TEST_CASE("Vector3Spherical - Comparison", "[VectorTypes][Vector3Spherical]")
{
    TEST_PRECISION_INFO();

    SECTION("Equality operator")
    {
        Vector3Spherical a(static_cast<MML::Real>(5.0), Constants::PI / 4, Constants::PI / 3);
        Vector3Spherical b(static_cast<MML::Real>(5.0), Constants::PI / 4, Constants::PI / 3);
        REQUIRE(a == b);
    }

    SECTION("IsEqualTo handles angle wrap-around")
    {
        Vector3Spherical a(static_cast<MML::Real>(5.0), Constants::PI / 4, Constants::PI - static_cast<MML::Real>(0.001));
        Vector3Spherical b(static_cast<MML::Real>(5.0), Constants::PI / 4, -Constants::PI + static_cast<MML::Real>(0.001));
        // These represent nearly the same direction (phi near ±π)
        REQUIRE(a.IsEqualTo(b, static_cast<MML::Real>(0.01)));
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         Vector3Cylindrical Tests                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Vector3Cylindrical - Constructors", "[VectorTypes][Vector3Cylindrical]")
{
    TEST_PRECISION_INFO();

    SECTION("Default constructor creates zero vector")
    {
        Vector3Cylindrical v;
        REQUIRE(v.R() == static_cast<MML::Real>(0.0));
        REQUIRE(v.Phi() == static_cast<MML::Real>(0.0));
        REQUIRE(v.Z() == static_cast<MML::Real>(0.0));
    }

    SECTION("Component constructor")
    {
        Vector3Cylindrical v(static_cast<MML::Real>(5.0), Constants::PI / 4, static_cast<MML::Real>(3.0));
        REQUIRE(v.R() == static_cast<MML::Real>(5.0));
        REQUIRE_THAT(v.Phi(), RealWithinRel(Constants::PI / 4, TOL(1e-10, 1e-5)));
        REQUIRE(v.Z() == static_cast<MML::Real>(3.0));
    }
}

TEST_CASE("Vector3Cylindrical - Arithmetic operations", "[VectorTypes][Vector3Cylindrical]")
{
    TEST_PRECISION_INFO();

    SECTION("Unary minus adds π to phi and negates Z")
    {
        Vector3Cylindrical v(static_cast<MML::Real>(5.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(3.0));
        Vector3Cylindrical neg = -v;
        REQUIRE(neg.R() == static_cast<MML::Real>(5.0));
        REQUIRE_THAT(neg.Phi(), RealWithinRel(Constants::PI, TOL(1e-10, 1e-5)));
        REQUIRE(neg.Z() == static_cast<MML::Real>(-3.0));
    }

    SECTION("Scalar multiplication scales R and Z")
    {
        Vector3Cylindrical v(static_cast<MML::Real>(5.0), Constants::PI / 4, static_cast<MML::Real>(3.0));
        Vector3Cylindrical scaled = v * static_cast<MML::Real>(2.0);
        REQUIRE(scaled.R() == static_cast<MML::Real>(10.0));
        REQUIRE_THAT(scaled.Phi(), RealWithinRel(Constants::PI / 4, TOL(1e-10, 1e-5)));
        REQUIRE(scaled.Z() == static_cast<MML::Real>(6.0));
    }

    SECTION("Division by zero throws")
    {
        Vector3Cylindrical v(static_cast<MML::Real>(5.0), Constants::PI / 4, static_cast<MML::Real>(3.0));
        REQUIRE_THROWS_AS(v / static_cast<MML::Real>(0.0), DivisionByZeroError);
    }
}

TEST_CASE("Vector3Cylindrical - GetAsUnitVectorAtPos", "[VectorTypes][Vector3Cylindrical]")
{
    TEST_PRECISION_INFO();

    SECTION("Radial vector at same position returns unit radial")
    {
        // Vector along +R (phi=0), evaluated at pos with phi=0
        Vector3Cylindrical v(static_cast<MML::Real>(5.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0));
        Vector3Cylindrical pos(static_cast<MML::Real>(3.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0));
        auto unit = v.GetAsUnitVectorAtPos(pos);
        REQUIRE_THAT(unit.R(), RealWithinAbs(static_cast<MML::Real>(1.0), TOL(1e-12, 1e-5)));
        REQUIRE_THAT(unit.Phi(), RealWithinAbs(static_cast<MML::Real>(0.0), TOL(1e-12, 1e-5)));
        REQUIRE_THAT(unit.Z(), RealWithinAbs(static_cast<MML::Real>(0.0), TOL(1e-12, 1e-5)));
    }

    SECTION("Z-only vector returns unit Z regardless of position")
    {
        Vector3Cylindrical v(static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(7.0));
        Vector3Cylindrical pos(static_cast<MML::Real>(2.0), Constants::PI / 3, static_cast<MML::Real>(5.0));
        auto unit = v.GetAsUnitVectorAtPos(pos);
        REQUIRE_THAT(unit.R(), RealWithinAbs(static_cast<MML::Real>(0.0), TOL(1e-12, 1e-5)));
        REQUIRE_THAT(unit.Phi(), RealWithinAbs(static_cast<MML::Real>(0.0), TOL(1e-12, 1e-5)));
        REQUIRE_THAT(unit.Z(), RealWithinAbs(static_cast<MML::Real>(1.0), TOL(1e-12, 1e-5)));
    }

    SECTION("Result is unit length")
    {
        Vector3Cylindrical v(static_cast<MML::Real>(3.0), Constants::PI / 6, static_cast<MML::Real>(4.0));
        Vector3Cylindrical pos(static_cast<MML::Real>(1.0), Constants::PI / 4, static_cast<MML::Real>(0.0));
        auto unit = v.GetAsUnitVectorAtPos(pos);
        Real normSq = unit.R() * unit.R() + unit.Phi() * unit.Phi() + unit.Z() * unit.Z();
        REQUIRE_THAT(normSq, RealWithinAbs(static_cast<MML::Real>(1.0), TOL(1e-12, 1e-5)));
    }

    SECTION("Zero vector returns zero")
    {
        Vector3Cylindrical v(static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0));
        Vector3Cylindrical pos(static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0));
        auto unit = v.GetAsUnitVectorAtPos(pos);
        REQUIRE(unit.R() == static_cast<MML::Real>(0.0));
        REQUIRE(unit.Phi() == static_cast<MML::Real>(0.0));
        REQUIRE(unit.Z() == static_cast<MML::Real>(0.0));
    }

    SECTION("Cross-validate with Cartesian unit vector")
    {
        // A cylindrical vector at phi=PI/4: Cartesian = (R*cos(pi/4), R*sin(pi/4), Z)
        Vector3Cylindrical v(static_cast<MML::Real>(2.0), Constants::PI / 4, static_cast<MML::Real>(2.0));
        Vector3Cylindrical pos(static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0));
        auto unit = v.GetAsUnitVectorAtPos(pos);

        // Convert original to Cartesian and normalize manually
        Real x = static_cast<MML::Real>(2.0) * std::cos(Constants::PI / 4);
        Real y = static_cast<MML::Real>(2.0) * std::sin(Constants::PI / 4);
        Real z = static_cast<MML::Real>(2.0);
        Real norm = std::sqrt(x * x + y * y + z * z);

        // At pos phi=0: ê_R=(1,0,0), ê_Phi=(0,1,0), ê_Z=(0,0,1)
        // So components equal Cartesian components
        REQUIRE_THAT(unit.R(), RealWithinAbs(x / norm, TOL(1e-12, 1e-5)));
        REQUIRE_THAT(unit.Phi(), RealWithinAbs(y / norm, TOL(1e-12, 1e-5)));
        REQUIRE_THAT(unit.Z(), RealWithinAbs(z / norm, TOL(1e-12, 1e-5)));
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         Vector4Minkowski Tests                                     ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Vector4Minkowski - Constructors", "[VectorTypes][Vector4Minkowski]")
{
    TEST_PRECISION_INFO();

    SECTION("Default constructor creates zero 4-vector")
    {
        Vector4Minkowski v;
        REQUIRE(v.T() == static_cast<MML::Real>(0.0));
        REQUIRE(v.X() == static_cast<MML::Real>(0.0));
        REQUIRE(v.Y() == static_cast<MML::Real>(0.0));
        REQUIRE(v.Z() == static_cast<MML::Real>(0.0));
    }

    SECTION("Initializer list constructor")
    {
        Vector4Minkowski v{static_cast<MML::Real>(5.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0)};
        REQUIRE(v.T() == static_cast<MML::Real>(5.0));
        REQUIRE(v.X() == static_cast<MML::Real>(1.0));
        REQUIRE(v.Y() == static_cast<MML::Real>(2.0));
        REQUIRE(v.Z() == static_cast<MML::Real>(3.0));
    }
}

TEST_CASE("Vector4Minkowski - Minkowski scalar product", "[VectorTypes][Vector4Minkowski]")
{
    TEST_PRECISION_INFO();

    SECTION("Minkowski metric: (-,+,+,+) convention")
    {
        Vector4Minkowski a{static_cast<MML::Real>(2.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0)};
        Vector4Minkowski b{static_cast<MML::Real>(3.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0)};
        
        // Minkowski product: -T1*T2 + X1*X2 + Y1*Y2 + Z1*Z2 = -2*3 + 1*1 = -5
        REQUIRE(ScalarProduct(a, b) == -static_cast<MML::Real>(5.0));
    }

    SECTION("Self product of timelike vector is negative")
    {
        Vector4Minkowski v{static_cast<MML::Real>(5.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0)};  // T² > X² + Y² + Z²
        Real selfProd = ScalarProduct(v, v);
        REQUIRE(selfProd < static_cast<MML::Real>(0.0));  // -25 + 3 = -22 < 0
    }

    SECTION("Self product of spacelike vector is positive")
    {
        Vector4Minkowski v{static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(2.0)};  // T² < X² + Y² + Z²
        Real selfProd = ScalarProduct(v, v);
        REQUIRE(selfProd > static_cast<MML::Real>(0.0));  // -1 + 12 = 11 > 0
    }

    SECTION("Self product of lightlike vector is zero")
    {
        Vector4Minkowski v{static_cast<MML::Real>(3.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(2.0)};  // T² = X² + Y² + Z² = 9
        Real selfProd = ScalarProduct(v, v);
        REQUIRE_THAT(selfProd, RealApprox(static_cast<MML::Real>(0.0)).margin(Tolerance::Loose));
    }
}

TEST_CASE("Vector4Minkowski - Classification", "[VectorTypes][Vector4Minkowski]")
{
    TEST_PRECISION_INFO();

    SECTION("Timelike vector (inside light cone)")
    {
        Vector4Minkowski v{static_cast<MML::Real>(5.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0)};
        REQUIRE(v.IsTimelike());
        REQUIRE_FALSE(v.IsSpacelike());
        REQUIRE_FALSE(v.IsLightlike());
    }

    SECTION("Spacelike vector (outside light cone)")
    {
        Vector4Minkowski v{static_cast<MML::Real>(1.0), static_cast<MML::Real>(3.0), static_cast<MML::Real>(3.0), static_cast<MML::Real>(3.0)};
        REQUIRE_FALSE(v.IsTimelike());
        REQUIRE(v.IsSpacelike());
        REQUIRE_FALSE(v.IsLightlike());
    }

    SECTION("Lightlike vector (on light cone)")
    {
        // T² = X² + Y² + Z²: 9 = 1 + 4 + 4
        Vector4Minkowski v{static_cast<MML::Real>(3.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(2.0)};
        Real interval = ScalarProduct(v, v);
        REQUIRE_THAT(interval, RealApprox(static_cast<MML::Real>(0.0)).margin(Tolerance::Loose));
    }
}

TEST_CASE("Vector4Minkowski - Norm", "[VectorTypes][Vector4Minkowski]")
{
    TEST_PRECISION_INFO();

    SECTION("Norm of zero vector")
    {
        Vector4Minkowski v;
        REQUIRE(v.Norm() == static_cast<MML::Real>(0.0));
    }

    SECTION("Norm of timelike vector")
    {
        Vector4Minkowski v{static_cast<MML::Real>(5.0), static_cast<MML::Real>(3.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0)};  // interval = 25 - 9 = 16
        REQUIRE_THAT(v.Norm(), RealWithinRel(static_cast<MML::Real>(4.0), TOL(1e-10, 1e-5)));  // sqrt(16)
    }

    SECTION("Norm of spacelike vector is negative")
    {
        Vector4Minkowski v{static_cast<MML::Real>(3.0), static_cast<MML::Real>(5.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0)};  // interval = 9 - 25 = -16
        REQUIRE_THAT(v.Norm(), RealWithinRel(static_cast<MML::Real>(-4.0), TOL(1e-10, 1e-5)));  // -sqrt(16)
    }

    SECTION("Norm of lightlike vector is zero")
    {
        Vector4Minkowski v{static_cast<MML::Real>(5.0), static_cast<MML::Real>(3.0), static_cast<MML::Real>(4.0), static_cast<MML::Real>(0.0)};  // interval = 25 - 9 - 16 = 0
        REQUIRE_THAT(v.Norm(), RealApprox(static_cast<MML::Real>(0.0)).margin(Tolerance::Loose));
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         Type Aliases Tests                                         ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("VectorTypes - Type aliases", "[VectorTypes][Aliases]")
{
    TEST_PRECISION_INFO();

    SECTION("Vec2Cart is Vector2Cartesian")
    {
        Vec2Cart v(static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0));
        REQUIRE(v.X() == static_cast<MML::Real>(1.0));
        REQUIRE(v.Y() == static_cast<MML::Real>(2.0));
    }

    SECTION("Vec3Cart is Vector3Cartesian")
    {
        Vec3Cart v(static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0));
        REQUIRE(v.X() == static_cast<MML::Real>(1.0));
        REQUIRE(v.Y() == static_cast<MML::Real>(2.0));
        REQUIRE(v.Z() == static_cast<MML::Real>(3.0));
    }

    SECTION("Vec3Sph is Vector3Spherical")
    {
        Vec3Sph v(static_cast<MML::Real>(5.0), Constants::PI / 4, Constants::PI / 3);
        REQUIRE(v.R() == static_cast<MML::Real>(5.0));
    }

    SECTION("Vec3Cyl is Vector3Cylindrical")
    {
        Vec3Cyl v(static_cast<MML::Real>(5.0), Constants::PI / 4, static_cast<MML::Real>(3.0));
        REQUIRE(v.R() == static_cast<MML::Real>(5.0));
    }

    SECTION("Vec4Mink is Vector4Minkowski")
    {
        Vec4Mink v{static_cast<MML::Real>(5.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0)};
        REQUIRE(v.T() == static_cast<MML::Real>(5.0));
    }
}

} // namespace MML::Tests::Base::VectorTypesTests

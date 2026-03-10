///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML) Tests                            ///
///                                                                                   ///
///  File:        vector_types_tests.cpp                                              ///
///  Description: Comprehensive tests for VectorTypes.h                               ///
///               Tests all vector type classes: Vec2Cart, Vec2Pol, Vec3Cart,         ///
///               Vec3Sph, Vec3Cyl, Vec4Mink                                          ///
///////////////////////////////////////////////////////////////////////////////////////////

#include <catch2/catch_all.hpp>
#include "../../TestPrecision.h"
#include "../../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Vector/VectorTypes.h"
#include "mml/base/Geometry/Geometry.h"
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
        REQUIRE(v.X() == REAL(0.0));
        REQUIRE(v.Y() == REAL(0.0));
    }

    SECTION("Component constructor")
    {
        Vector2Cartesian v(REAL(3.0), REAL(4.0));
        REQUIRE(v.X() == REAL(3.0));
        REQUIRE(v.Y() == REAL(4.0));
    }

    SECTION("Construct from VectorN")
    {
        VectorN<Real, 2> base{REAL(5.0), REAL(6.0)};
        Vector2Cartesian v(base);
        REQUIRE(v.X() == REAL(5.0));
        REQUIRE(v.Y() == REAL(6.0));
    }

    SECTION("Construct from two points")
    {
        Point2Cartesian a(REAL(1.0), REAL(2.0));
        Point2Cartesian b(REAL(4.0), REAL(6.0));
        Vector2Cartesian v(a, b);
        REQUIRE(v.X() == REAL(3.0));
        REQUIRE(v.Y() == REAL(4.0));
    }

    SECTION("Initializer list constructor")
    {
        Vector2Cartesian v{REAL(7.0), REAL(8.0)};
        REQUIRE(v.X() == REAL(7.0));
        REQUIRE(v.Y() == REAL(8.0));
    }
}

TEST_CASE("Vector2Cartesian - Arithmetic operations", "[VectorTypes][Vector2Cartesian]")
{
    TEST_PRECISION_INFO();
    Vector2Cartesian a(REAL(3.0), REAL(4.0));
    Vector2Cartesian b(REAL(1.0), REAL(2.0));

    SECTION("Vector addition")
    {
        Vector2Cartesian c = a + b;
        REQUIRE(c.X() == REAL(4.0));
        REQUIRE(c.Y() == REAL(6.0));
    }

    SECTION("Vector subtraction")
    {
        Vector2Cartesian c = a - b;
        REQUIRE(c.X() == REAL(2.0));
        REQUIRE(c.Y() == REAL(2.0));
    }

    SECTION("Unary minus (negation)")
    {
        Vector2Cartesian c = -a;
        REQUIRE(c.X() == REAL(-3.0));
        REQUIRE(c.Y() == REAL(-4.0));
    }

    SECTION("Scalar multiplication (vector * scalar)")
    {
        Vector2Cartesian c = a * REAL(2.0);
        REQUIRE(c.X() == REAL(6.0));
        REQUIRE(c.Y() == REAL(8.0));
    }

    SECTION("Scalar multiplication (scalar * vector)")
    {
        Vector2Cartesian c = REAL(2.0) * a;
        REQUIRE(c.X() == REAL(6.0));
        REQUIRE(c.Y() == REAL(8.0));
    }

    SECTION("Scalar division")
    {
        Vector2Cartesian c = a / REAL(2.0);
        REQUIRE(c.X() == REAL(1.5));
        REQUIRE(c.Y() == REAL(2.0));
    }

    SECTION("Scalar product (dot product)")
    {
        Real dot = a * b;
        REQUIRE(dot == REAL(11.0));  // 3*1 + 4*2 = 11
    }

    SECTION("ScalarProduct friend function")
    {
        Real dot = ScalarProduct(a, b);
        REQUIRE(dot == REAL(11.0));
    }
}

TEST_CASE("Vector2Cartesian - Comparison operations", "[VectorTypes][Vector2Cartesian]")
{
    TEST_PRECISION_INFO();
    Vector2Cartesian a(REAL(3.0), REAL(4.0));
    Vector2Cartesian b(REAL(3.0), REAL(4.0));
    Vector2Cartesian c(REAL(3.1), REAL(4.0));

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
        REQUIRE(a.IsEqualTo(b, REAL(1e-10)));
        REQUIRE(a.IsEqualTo(c, REAL(0.2)));  // Within tolerance
        REQUIRE_FALSE(a.IsEqualTo(c, REAL(0.05)));  // Outside tolerance
    }
}

TEST_CASE("Vector2Cartesian - Geometric operations", "[VectorTypes][Vector2Cartesian]")
{
    TEST_PRECISION_INFO();

    SECTION("NormL2 (magnitude)")
    {
        Vector2Cartesian v(REAL(3.0), REAL(4.0));
        REQUIRE(v.NormL2() == REAL(5.0));
    }

    SECTION("Normalized returns unit vector")
    {
        Vector2Cartesian v(REAL(3.0), REAL(4.0));
        Vector2Cartesian unit = v.Normalized();
        REQUIRE_THAT(unit.X(), RealWithinRel(REAL(0.6), REAL(1e-10)));
        REQUIRE_THAT(unit.Y(), RealWithinRel(REAL(0.8), REAL(1e-10)));
        REQUIRE_THAT(unit.NormL2(), RealWithinRel(REAL(1.0), REAL(1e-10)));
    }

    SECTION("Normalized of zero vector returns zero")
    {
        Vector2Cartesian v(REAL(0.0), REAL(0.0));
        Vector2Cartesian unit = v.Normalized();
        REQUIRE(unit.X() == REAL(0.0));
        REQUIRE(unit.Y() == REAL(0.0));
    }

    SECTION("GetAsUnitVector")
    {
        Vector2Cartesian v(REAL(0.0), REAL(5.0));
        Vector2Cartesian unit = v.GetAsUnitVector();
        REQUIRE(unit.X() == REAL(0.0));
        REQUIRE_THAT(unit.Y(), RealWithinRel(REAL(1.0), REAL(1e-10)));
    }

    SECTION("getPerpendicularVectors")
    {
        Vector2Cartesian v(REAL(3.0), REAL(4.0));
        Vector2Cartesian v1, v2;
        v.getPerpendicularVectors(v1, v2);
        
        // v1 and v2 should be perpendicular to v
        REQUIRE_THAT((v * v1), RealApprox(REAL(0.0)).margin(Tolerance::Loose));
        REQUIRE_THAT((v * v2), RealApprox(REAL(0.0)).margin(Tolerance::Loose));
        
        // v1 and v2 should be opposite
        REQUIRE_THAT(v1.X(), RealWithinRel(-v2.X(), REAL(1e-10)));
        REQUIRE_THAT(v1.Y(), RealWithinRel(-v2.Y(), REAL(1e-10)));
    }
}

TEST_CASE("Vector2Cartesian - Point operations", "[VectorTypes][Vector2Cartesian]")
{
    TEST_PRECISION_INFO();
    Point2Cartesian p(REAL(1.0), REAL(2.0));
    Vector2Cartesian v(REAL(3.0), REAL(4.0));

    SECTION("Point + Vector")
    {
        Point2Cartesian result = p + v;
        REQUIRE(result.X() == REAL(4.0));
        REQUIRE(result.Y() == REAL(6.0));
    }

    SECTION("Point - Vector")
    {
        Point2Cartesian result = p - v;
        REQUIRE(result.X() == REAL(-2.0));
        REQUIRE(result.Y() == REAL(-2.0));
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
        REQUIRE(v.R() == REAL(0.0));
        REQUIRE(v.Phi() == REAL(0.0));
    }

    SECTION("Component constructor")
    {
        Vector2Polar v(REAL(5.0), Constants::PI / 4);
        REQUIRE(v.R() == REAL(5.0));
        REQUIRE_THAT(v.Phi(), RealWithinRel(Constants::PI / 4, REAL(1e-10)));
    }
}

TEST_CASE("Vector2Polar - Arithmetic operations", "[VectorTypes][Vector2Polar]")
{
    TEST_PRECISION_INFO();

    SECTION("Unary minus adds π to angle")
    {
        Vector2Polar v(REAL(5.0), REAL(0.0));
        Vector2Polar neg = -v;
        REQUIRE(neg.R() == REAL(5.0));
        REQUIRE_THAT(neg.Phi(), RealWithinRel(Constants::PI, REAL(1e-10)));
    }

    SECTION("Scalar multiplication scales radius only")
    {
        Vector2Polar v(REAL(5.0), Constants::PI / 4);
        Vector2Polar scaled = v * REAL(2.0);
        REQUIRE(scaled.R() == REAL(10.0));
        REQUIRE_THAT(scaled.Phi(), RealWithinRel(Constants::PI / 4, REAL(1e-10)));
    }

    SECTION("Scalar multiplication (scalar * vector)")
    {
        Vector2Polar v(REAL(5.0), Constants::PI / 4);
        Vector2Polar scaled = REAL(3.0) * v;
        REQUIRE(scaled.R() == REAL(15.0));
        REQUIRE_THAT(scaled.Phi(), RealWithinRel(Constants::PI / 4, REAL(1e-10)));
    }

    SECTION("Scalar division divides radius only")
    {
        Vector2Polar v(REAL(10.0), Constants::PI / 3);
        Vector2Polar divided = v / REAL(2.0);
        REQUIRE(divided.R() == REAL(5.0));
        REQUIRE_THAT(divided.Phi(), RealWithinRel(Constants::PI / 3, REAL(1e-10)));
    }

    SECTION("Vector addition - same direction")
    {
        Vector2Polar v1(REAL(3.0), REAL(0.0));
        Vector2Polar v2(REAL(4.0), REAL(0.0));
        Vector2Polar sum = v1 + v2;
        REQUIRE_THAT(sum.R(), RealWithinRel(REAL(7.0), REAL(1e-10)));
    }

    SECTION("Vector addition - opposite directions")
    {
        Vector2Polar v1(REAL(5.0), REAL(0.0));
        Vector2Polar v2(REAL(3.0), Constants::PI);
        Vector2Polar sum = v1 + v2;
        REQUIRE_THAT(sum.R(), RealWithinRel(REAL(2.0), REAL(1e-10)));
    }

    SECTION("Vector addition - perpendicular")
    {
        Vector2Polar v1(REAL(3.0), REAL(0.0));
        Vector2Polar v2(REAL(4.0), Constants::PI / 2);
        Vector2Polar sum = v1 + v2;
        REQUIRE_THAT(sum.R(), RealWithinRel(REAL(5.0), REAL(1e-10)));  // 3-4-5 triangle
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
        REQUIRE(v.X() == REAL(0.0));
        REQUIRE(v.Y() == REAL(0.0));
        REQUIRE(v.Z() == REAL(0.0));
    }

    SECTION("Component constructor")
    {
        Vector3Cartesian v(REAL(1.0), REAL(2.0), REAL(3.0));
        REQUIRE(v.X() == REAL(1.0));
        REQUIRE(v.Y() == REAL(2.0));
        REQUIRE(v.Z() == REAL(3.0));
    }

    SECTION("Construct from VectorN")
    {
        VectorN<Real, 3> base{REAL(4.0), REAL(5.0), REAL(6.0)};
        Vector3Cartesian v(base);
        REQUIRE(v.X() == REAL(4.0));
        REQUIRE(v.Y() == REAL(5.0));
        REQUIRE(v.Z() == REAL(6.0));
    }

    SECTION("Construct from two points")
    {
        Point3Cartesian a(REAL(1.0), REAL(2.0), REAL(3.0));
        Point3Cartesian b(REAL(4.0), REAL(6.0), REAL(8.0));
        Vector3Cartesian v(a, b);
        REQUIRE(v.X() == REAL(3.0));
        REQUIRE(v.Y() == REAL(4.0));
        REQUIRE(v.Z() == REAL(5.0));
    }

    SECTION("Construct from Point3Cartesian")
    {
        Point3Cartesian p(REAL(7.0), REAL(8.0), REAL(9.0));
        Vector3Cartesian v(p);
        REQUIRE(v.X() == REAL(7.0));
        REQUIRE(v.Y() == REAL(8.0));
        REQUIRE(v.Z() == REAL(9.0));
    }

    SECTION("Initializer list constructor")
    {
        Vector3Cartesian v{REAL(10.0), REAL(11.0), REAL(12.0)};
        REQUIRE(v.X() == REAL(10.0));
        REQUIRE(v.Y() == REAL(11.0));
        REQUIRE(v.Z() == REAL(12.0));
    }
}

TEST_CASE("Vector3Cartesian - Arithmetic operations", "[VectorTypes][Vector3Cartesian]")
{
    TEST_PRECISION_INFO();
    Vector3Cartesian a(REAL(1.0), REAL(2.0), REAL(3.0));
    Vector3Cartesian b(REAL(4.0), REAL(5.0), REAL(6.0));

    SECTION("Vector addition")
    {
        Vector3Cartesian c = a + b;
        REQUIRE(c.X() == REAL(5.0));
        REQUIRE(c.Y() == REAL(7.0));
        REQUIRE(c.Z() == REAL(9.0));
    }

    SECTION("Vector subtraction")
    {
        Vector3Cartesian c = b - a;
        REQUIRE(c.X() == REAL(3.0));
        REQUIRE(c.Y() == REAL(3.0));
        REQUIRE(c.Z() == REAL(3.0));
    }

    SECTION("Unary minus")
    {
        Vector3Cartesian c = -a;
        REQUIRE(c.X() == REAL(-1.0));
        REQUIRE(c.Y() == REAL(-2.0));
        REQUIRE(c.Z() == REAL(-3.0));
    }

    SECTION("Scalar multiplication")
    {
        Vector3Cartesian c = a * REAL(2.0);
        REQUIRE(c.X() == REAL(2.0));
        REQUIRE(c.Y() == REAL(4.0));
        REQUIRE(c.Z() == REAL(6.0));
    }

    SECTION("Scalar division")
    {
        Vector3Cartesian c = b / REAL(2.0);
        REQUIRE(c.X() == REAL(2.0));
        REQUIRE(c.Y() == REAL(2.5));
        REQUIRE(c.Z() == REAL(3.0));
    }

    SECTION("Scalar product (dot product)")
    {
        Real dot = a * b;
        REQUIRE(dot == REAL(32.0));  // 1*4 + 2*5 + 3*6 = 32
    }

    SECTION("ScalarProduct member function")
    {
        REQUIRE(a.ScalarProduct(b) == REAL(32.0));
    }

    SECTION("ScalarProduct friend function")
    {
        REQUIRE(ScalarProduct(a, b) == REAL(32.0));
    }
}

TEST_CASE("Vector3Cartesian - Cross product", "[VectorTypes][Vector3Cartesian]")
{
    TEST_PRECISION_INFO();

    SECTION("Cross product of unit vectors")
    {
        Vector3Cartesian i(REAL(1.0), REAL(0.0), REAL(0.0));
        Vector3Cartesian j(REAL(0.0), REAL(1.0), REAL(0.0));
        Vector3Cartesian k(REAL(0.0), REAL(0.0), REAL(1.0));

        // i × j = k
        Vector3Cartesian ixj = VectorProduct(i, j);
        REQUIRE(ixj.X() == REAL(0.0));
        REQUIRE(ixj.Y() == REAL(0.0));
        REQUIRE(ixj.Z() == REAL(1.0));

        // j × k = i
        Vector3Cartesian jxk = VectorProduct(j, k);
        REQUIRE(jxk.X() == REAL(1.0));
        REQUIRE(jxk.Y() == REAL(0.0));
        REQUIRE(jxk.Z() == REAL(0.0));

        // k × i = j
        Vector3Cartesian kxi = VectorProduct(k, i);
        REQUIRE(kxi.X() == REAL(0.0));
        REQUIRE(kxi.Y() == REAL(1.0));
        REQUIRE(kxi.Z() == REAL(0.0));
    }

    SECTION("Cross product anti-commutativity")
    {
        Vector3Cartesian a(REAL(1.0), REAL(2.0), REAL(3.0));
        Vector3Cartesian b(REAL(4.0), REAL(5.0), REAL(6.0));

        Vector3Cartesian axb = VectorProduct(a, b);
        Vector3Cartesian bxa = VectorProduct(b, a);

        REQUIRE_THAT(axb.X(), RealWithinRel(-bxa.X(), REAL(1e-10)));
        REQUIRE_THAT(axb.Y(), RealWithinRel(-bxa.Y(), REAL(1e-10)));
        REQUIRE_THAT(axb.Z(), RealWithinRel(-bxa.Z(), REAL(1e-10)));
    }

    SECTION("Cross product of parallel vectors is zero")
    {
        Vector3Cartesian a(REAL(1.0), REAL(2.0), REAL(3.0));
        Vector3Cartesian b(REAL(2.0), REAL(4.0), REAL(6.0));

        Vector3Cartesian cross = VectorProduct(a, b);
        REQUIRE_THAT(cross.NormL2(), RealApprox(REAL(0.0)).margin(Tolerance::Loose));
    }

    SECTION("Cross product perpendicular to both inputs")
    {
        Vector3Cartesian a(REAL(1.0), REAL(2.0), REAL(3.0));
        Vector3Cartesian b(REAL(4.0), REAL(5.0), REAL(6.0));

        Vector3Cartesian cross = VectorProduct(a, b);
        REQUIRE_THAT(ScalarProduct(cross, a), RealApprox(REAL(0.0)).margin(Tolerance::Loose));
        REQUIRE_THAT(ScalarProduct(cross, b), RealApprox(REAL(0.0)).margin(Tolerance::Loose));
    }
}

TEST_CASE("Vector3Cartesian - Geometric operations", "[VectorTypes][Vector3Cartesian]")
{
    TEST_PRECISION_INFO();

    SECTION("NormL2")
    {
        Vector3Cartesian v(REAL(2.0), REAL(3.0), REAL(6.0));
        REQUIRE(v.NormL2() == REAL(7.0));  // 2² + 3² + 6² = 49 = 7²
    }

    SECTION("Normalized")
    {
        Vector3Cartesian v(REAL(0.0), REAL(0.0), REAL(5.0));
        Vector3Cartesian unit = v.Normalized();
        REQUIRE(unit.X() == REAL(0.0));
        REQUIRE(unit.Y() == REAL(0.0));
        REQUIRE_THAT(unit.Z(), RealWithinRel(REAL(1.0), REAL(1e-10)));
    }

    SECTION("Normalized of zero vector")
    {
        Vector3Cartesian v(REAL(0.0), REAL(0.0), REAL(0.0));
        Vector3Cartesian unit = v.Normalized();
        REQUIRE(unit.X() == REAL(0.0));
        REQUIRE(unit.Y() == REAL(0.0));
        REQUIRE(unit.Z() == REAL(0.0));
    }

    SECTION("GetAsUnitVector")
    {
        Vector3Cartesian v(REAL(3.0), REAL(4.0), REAL(0.0));
        Vector3Cartesian unit = v.GetAsUnitVector();
        REQUIRE_THAT(unit.NormL2(), RealWithinRel(REAL(1.0), REAL(1e-10)));
        REQUIRE_THAT(unit.X(), RealWithinRel(REAL(0.6), REAL(1e-10)));
        REQUIRE_THAT(unit.Y(), RealWithinRel(REAL(0.8), REAL(1e-10)));
    }

    SECTION("GetAsUnitVector of zero vector returns zero")
    {
        Vector3Cartesian v(REAL(0.0), REAL(0.0), REAL(0.0));
        Vector3Cartesian unit = v.GetAsUnitVector();
        REQUIRE(unit.X() == REAL(0.0));
        REQUIRE(unit.Y() == REAL(0.0));
        REQUIRE(unit.Z() == REAL(0.0));
    }

    SECTION("IsParallelTo")
    {
        Vector3Cartesian a(REAL(1.0), REAL(2.0), REAL(3.0));
        Vector3Cartesian b(REAL(2.0), REAL(4.0), REAL(6.0));  // Parallel
        Vector3Cartesian c(REAL(1.0), REAL(0.0), REAL(0.0));  // Not parallel

        REQUIRE(a.IsParallelTo(b));
        REQUIRE_FALSE(a.IsParallelTo(c));
    }

    SECTION("IsPerpendicularTo")
    {
        Vector3Cartesian a(REAL(1.0), REAL(0.0), REAL(0.0));
        Vector3Cartesian b(REAL(0.0), REAL(1.0), REAL(0.0));  // Perpendicular
        Vector3Cartesian c(REAL(1.0), REAL(1.0), REAL(0.0));  // Not perpendicular

        REQUIRE(a.IsPerpendicularTo(b));
        REQUIRE_FALSE(a.IsPerpendicularTo(c));
    }

    SECTION("AngleToVector")
    {
        Vector3Cartesian a(REAL(1.0), REAL(0.0), REAL(0.0));
        Vector3Cartesian b(REAL(0.0), REAL(1.0), REAL(0.0));
        Vector3Cartesian c(REAL(1.0), REAL(0.0), REAL(0.0));
        Vector3Cartesian d(REAL(-1.0), REAL(0.0), REAL(0.0));

        REQUIRE_THAT(a.AngleToVector(b), RealWithinRel(Constants::PI / 2, REAL(1e-10)));  // 90°
        REQUIRE_THAT(a.AngleToVector(c), RealApprox(REAL(0.0)).margin(Tolerance::Loose));  // 0°
        REQUIRE_THAT(a.AngleToVector(d), RealWithinRel(Constants::PI, REAL(1e-10)));  // 180°
    }

    SECTION("GetPerpendicularVectors")
    {
        Vector3Cartesian v(REAL(1.0), REAL(2.0), REAL(3.0));
        Vector3Cartesian v1, v2;

        bool success = v.GetPerpendicularVectors(v1, v2);
        REQUIRE(success);

        // Both should be perpendicular to v
        REQUIRE_THAT(v.ScalarProduct(v1), RealApprox(REAL(0.0)).margin(Tolerance::Loose));
        REQUIRE_THAT(v.ScalarProduct(v2), RealApprox(REAL(0.0)).margin(Tolerance::Loose));

        // v1 and v2 should be perpendicular to each other
        REQUIRE_THAT(v1.ScalarProduct(v2), RealApprox(REAL(0.0)).margin(Tolerance::Loose));

        // Both should be unit vectors
        REQUIRE_THAT(v1.NormL2(), RealWithinRel(REAL(1.0), REAL(1e-10)));
        REQUIRE_THAT(v2.NormL2(), RealWithinRel(REAL(1.0), REAL(1e-10)));
    }

    SECTION("GetPerpendicularVectors of zero vector fails")
    {
        Vector3Cartesian v(REAL(0.0), REAL(0.0), REAL(0.0));
        Vector3Cartesian v1, v2;

        bool success = v.GetPerpendicularVectors(v1, v2);
        REQUIRE_FALSE(success);
    }
}

TEST_CASE("Vector3Cartesian - Point operations", "[VectorTypes][Vector3Cartesian]")
{
    TEST_PRECISION_INFO();
    Point3Cartesian p(REAL(1.0), REAL(2.0), REAL(3.0));
    Vector3Cartesian v(REAL(4.0), REAL(5.0), REAL(6.0));

    SECTION("Point + Vector")
    {
        Point3Cartesian result = p + v;
        REQUIRE(result.X() == REAL(5.0));
        REQUIRE(result.Y() == REAL(7.0));
        REQUIRE(result.Z() == REAL(9.0));
    }

    SECTION("Point - Vector")
    {
        Point3Cartesian result = p - v;
        REQUIRE(result.X() == REAL(-3.0));
        REQUIRE(result.Y() == REAL(-3.0));
        REQUIRE(result.Z() == REAL(-3.0));
    }

    SECTION("getAsPoint")
    {
        Point3Cartesian pt = v.getAsPoint();
        REQUIRE(pt.X() == REAL(4.0));
        REQUIRE(pt.Y() == REAL(5.0));
        REQUIRE(pt.Z() == REAL(6.0));
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
        REQUIRE(v.R() == REAL(0.0));
        REQUIRE(v.Theta() == REAL(0.0));
        REQUIRE(v.Phi() == REAL(0.0));
    }

    SECTION("Component constructor")
    {
        Vector3Spherical v(REAL(5.0), Constants::PI / 4, Constants::PI / 3);
        REQUIRE(v.R() == REAL(5.0));
        REQUIRE_THAT(v.Theta(), RealWithinRel(Constants::PI / 4, REAL(1e-10)));
        REQUIRE_THAT(v.Phi(), RealWithinRel(Constants::PI / 3, REAL(1e-10)));
    }
}

TEST_CASE("Vector3Spherical - Arithmetic operations", "[VectorTypes][Vector3Spherical]")
{
    TEST_PRECISION_INFO();

    SECTION("Scalar multiplication scales radius only")
    {
        Vector3Spherical v(REAL(5.0), Constants::PI / 4, Constants::PI / 3);
        Vector3Spherical scaled = v * REAL(2.0);
        REQUIRE(scaled.R() == REAL(10.0));
        REQUIRE_THAT(scaled.Theta(), RealWithinRel(Constants::PI / 4, REAL(1e-10)));
        REQUIRE_THAT(scaled.Phi(), RealWithinRel(Constants::PI / 3, REAL(1e-10)));
    }

    SECTION("Scalar division")
    {
        Vector3Spherical v(REAL(10.0), Constants::PI / 4, Constants::PI / 3);
        Vector3Spherical divided = v / REAL(2.0);
        REQUIRE(divided.R() == REAL(5.0));
        REQUIRE_THAT(divided.Theta(), RealWithinRel(Constants::PI / 4, REAL(1e-10)));
        REQUIRE_THAT(divided.Phi(), RealWithinRel(Constants::PI / 3, REAL(1e-10)));
    }

    SECTION("Division by zero throws")
    {
        Vector3Spherical v(REAL(10.0), Constants::PI / 4, Constants::PI / 3);
        REQUIRE_THROWS_AS(v / REAL(0.0), DivisionByZeroError);
    }

    SECTION("Vector addition - same direction")
    {
        Vector3Spherical v1(REAL(3.0), Constants::PI / 2, REAL(0.0));  // Along +x axis
        Vector3Spherical v2(REAL(4.0), Constants::PI / 2, REAL(0.0));
        Vector3Spherical sum = v1 + v2;
        REQUIRE_THAT(sum.R(), RealWithinRel(REAL(7.0), REAL(1e-10)));
    }
}

TEST_CASE("Vector3Spherical - Comparison", "[VectorTypes][Vector3Spherical]")
{
    TEST_PRECISION_INFO();

    SECTION("Equality operator")
    {
        Vector3Spherical a(REAL(5.0), Constants::PI / 4, Constants::PI / 3);
        Vector3Spherical b(REAL(5.0), Constants::PI / 4, Constants::PI / 3);
        REQUIRE(a == b);
    }

    SECTION("IsEqualTo handles angle wrap-around")
    {
        Vector3Spherical a(REAL(5.0), Constants::PI / 4, Constants::PI - REAL(0.001));
        Vector3Spherical b(REAL(5.0), Constants::PI / 4, -Constants::PI + REAL(0.001));
        // These represent nearly the same direction (phi near ±π)
        REQUIRE(a.IsEqualTo(b, REAL(0.01)));
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
        REQUIRE(v.R() == REAL(0.0));
        REQUIRE(v.Phi() == REAL(0.0));
        REQUIRE(v.Z() == REAL(0.0));
    }

    SECTION("Component constructor")
    {
        Vector3Cylindrical v(REAL(5.0), Constants::PI / 4, REAL(3.0));
        REQUIRE(v.R() == REAL(5.0));
        REQUIRE_THAT(v.Phi(), RealWithinRel(Constants::PI / 4, REAL(1e-10)));
        REQUIRE(v.Z() == REAL(3.0));
    }
}

TEST_CASE("Vector3Cylindrical - Arithmetic operations", "[VectorTypes][Vector3Cylindrical]")
{
    TEST_PRECISION_INFO();

    SECTION("Unary minus adds π to phi")
    {
        Vector3Cylindrical v(REAL(5.0), REAL(0.0), REAL(3.0));
        Vector3Cylindrical neg = -v;
        REQUIRE(neg.R() == REAL(5.0));
        REQUIRE_THAT(neg.Phi(), RealWithinRel(Constants::PI, REAL(1e-10)));
        REQUIRE(neg.Z() == REAL(3.0));
    }

    SECTION("Scalar multiplication scales R and Z")
    {
        Vector3Cylindrical v(REAL(5.0), Constants::PI / 4, REAL(3.0));
        Vector3Cylindrical scaled = v * REAL(2.0);
        REQUIRE(scaled.R() == REAL(10.0));
        REQUIRE_THAT(scaled.Phi(), RealWithinRel(Constants::PI / 4, REAL(1e-10)));
        REQUIRE(scaled.Z() == REAL(6.0));
    }

    SECTION("Division by zero throws")
    {
        Vector3Cylindrical v(REAL(5.0), Constants::PI / 4, REAL(3.0));
        REQUIRE_THROWS_AS(v / REAL(0.0), DivisionByZeroError);
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
        REQUIRE(v.T() == REAL(0.0));
        REQUIRE(v.X() == REAL(0.0));
        REQUIRE(v.Y() == REAL(0.0));
        REQUIRE(v.Z() == REAL(0.0));
    }

    SECTION("Initializer list constructor")
    {
        Vector4Minkowski v{REAL(5.0), REAL(1.0), REAL(2.0), REAL(3.0)};
        REQUIRE(v.T() == REAL(5.0));
        REQUIRE(v.X() == REAL(1.0));
        REQUIRE(v.Y() == REAL(2.0));
        REQUIRE(v.Z() == REAL(3.0));
    }
}

TEST_CASE("Vector4Minkowski - Minkowski scalar product", "[VectorTypes][Vector4Minkowski]")
{
    TEST_PRECISION_INFO();

    SECTION("Minkowski metric: (+,-,-,-) convention")
    {
        Vector4Minkowski a{REAL(2.0), REAL(1.0), REAL(0.0), REAL(0.0)};
        Vector4Minkowski b{REAL(3.0), REAL(1.0), REAL(0.0), REAL(0.0)};
        
        // Minkowski product: T₁T₂ - X₁X₂ - Y₁Y₂ - Z₁Z₂ = 2*3 - 1*1 = 5
        REQUIRE(ScalarProduct(a, b) == REAL(5.0));
    }

    SECTION("Self product of timelike vector is positive")
    {
        Vector4Minkowski v{REAL(5.0), REAL(1.0), REAL(1.0), REAL(1.0)};  // T² > X² + Y² + Z²
        Real selfProd = ScalarProduct(v, v);
        REQUIRE(selfProd > REAL(0.0));  // 25 - 3 = 22 > 0
    }

    SECTION("Self product of spacelike vector is negative")
    {
        Vector4Minkowski v{REAL(1.0), REAL(2.0), REAL(2.0), REAL(2.0)};  // T² < X² + Y² + Z²
        Real selfProd = ScalarProduct(v, v);
        REQUIRE(selfProd < REAL(0.0));  // 1 - 12 = -11 < 0
    }

    SECTION("Self product of lightlike vector is zero")
    {
        Vector4Minkowski v{REAL(3.0), REAL(1.0), REAL(2.0), REAL(2.0)};  // T² = X² + Y² + Z² = 9
        Real selfProd = ScalarProduct(v, v);
        REQUIRE_THAT(selfProd, RealApprox(REAL(0.0)).margin(Tolerance::Loose));
    }
}

TEST_CASE("Vector4Minkowski - Classification", "[VectorTypes][Vector4Minkowski]")
{
    TEST_PRECISION_INFO();

    SECTION("Timelike vector (inside light cone)")
    {
        Vector4Minkowski v{REAL(5.0), REAL(1.0), REAL(1.0), REAL(1.0)};
        REQUIRE(v.isTimelike());
        REQUIRE_FALSE(v.isSpacelike());
        REQUIRE_FALSE(v.isLightlike());
    }

    SECTION("Spacelike vector (outside light cone)")
    {
        Vector4Minkowski v{REAL(1.0), REAL(3.0), REAL(3.0), REAL(3.0)};
        REQUIRE_FALSE(v.isTimelike());
        REQUIRE(v.isSpacelike());
        REQUIRE_FALSE(v.isLightlike());
    }

    SECTION("Lightlike vector (on light cone)")
    {
        // T² = X² + Y² + Z²: 9 = 1 + 4 + 4
        Vector4Minkowski v{REAL(3.0), REAL(1.0), REAL(2.0), REAL(2.0)};
        Real interval = v.T() * v.T() - v.X() * v.X() - v.Y() * v.Y() - v.Z() * v.Z();
        REQUIRE_THAT(interval, RealApprox(REAL(0.0)).margin(Tolerance::Loose));
    }
}

TEST_CASE("Vector4Minkowski - Norm", "[VectorTypes][Vector4Minkowski]")
{
    TEST_PRECISION_INFO();

    SECTION("Norm of zero vector")
    {
        Vector4Minkowski v;
        REQUIRE(v.Norm() == REAL(0.0));
    }

    SECTION("Norm of timelike vector")
    {
        Vector4Minkowski v{REAL(5.0), REAL(3.0), REAL(0.0), REAL(0.0)};  // interval = 25 - 9 = 16
        REQUIRE_THAT(v.Norm(), RealWithinRel(REAL(4.0), REAL(1e-10)));  // sqrt(16)
    }

    SECTION("Norm of spacelike vector is negative")
    {
        Vector4Minkowski v{REAL(3.0), REAL(5.0), REAL(0.0), REAL(0.0)};  // interval = 9 - 25 = -16
        REQUIRE_THAT(v.Norm(), RealWithinRel(REAL(-4.0), REAL(1e-10)));  // -sqrt(16)
    }

    SECTION("Norm of lightlike vector is zero")
    {
        Vector4Minkowski v{REAL(5.0), REAL(3.0), REAL(4.0), REAL(0.0)};  // interval = 25 - 9 - 16 = 0
        REQUIRE_THAT(v.Norm(), RealApprox(REAL(0.0)).margin(Tolerance::Loose));
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
        Vec2Cart v(REAL(1.0), REAL(2.0));
        REQUIRE(v.X() == REAL(1.0));
        REQUIRE(v.Y() == REAL(2.0));
    }

    SECTION("Vec3Cart is Vector3Cartesian")
    {
        Vec3Cart v(REAL(1.0), REAL(2.0), REAL(3.0));
        REQUIRE(v.X() == REAL(1.0));
        REQUIRE(v.Y() == REAL(2.0));
        REQUIRE(v.Z() == REAL(3.0));
    }

    SECTION("Vec3Sph is Vector3Spherical")
    {
        Vec3Sph v(REAL(5.0), Constants::PI / 4, Constants::PI / 3);
        REQUIRE(v.R() == REAL(5.0));
    }

    SECTION("Vec3Cyl is Vector3Cylindrical")
    {
        Vec3Cyl v(REAL(5.0), Constants::PI / 4, REAL(3.0));
        REQUIRE(v.R() == REAL(5.0));
    }

    SECTION("Vec4Mink is Vector4Minkowski")
    {
        Vec4Mink v{REAL(5.0), REAL(1.0), REAL(2.0), REAL(3.0)};
        REQUIRE(v.T() == REAL(5.0));
    }
}

} // namespace MML::Tests::Base::VectorTypesTests

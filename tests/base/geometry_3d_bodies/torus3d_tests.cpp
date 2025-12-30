#include <catch2/catch_all.hpp>
#include "../../TestPrecision.h"
#include "../../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Geometry3DBodies.h"
#endif

using namespace MML;
using namespace MML::Testing;

TEST_CASE("Torus3D::Volume", "[geometry][torus][volume]")
{
    SECTION("Torus with R=3, r=1")
    {
        Torus3D torus(3.0, 1.0);
        // Volume = 2π²Rr² = 2π² * 3 * 1 = 6π²
        Real expected = 6.0 * Constants::PI * Constants::PI;
        REQUIRE_THAT(torus.Volume() , RealApprox(expected));
    }

    SECTION("Torus with R=5, r=2")
    {
        Torus3D torus(5.0, 2.0);
        // Volume = 2π² * 5 * 4 = 40π²
        Real expected = 40.0 * Constants::PI * Constants::PI;
        REQUIRE_THAT(torus.Volume() , RealApprox(expected));
    }

    SECTION("Torus centered at arbitrary point")
    {
        Pnt3Cart center(10.0, 20.0, 30.0);
        Torus3D torus(4.0, 1.5, center);
        // Volume = 2π² * 4 * 2.25 = 18π²
        Real expected = 18.0 * Constants::PI * Constants::PI;
        REQUIRE_THAT(torus.Volume() , RealApprox(expected));
    }

    SECTION("Thin torus - large R, small r")
    {
        Torus3D torus(10.0, 0.5);
        // Volume = 2π² * 10 * 0.25 = 5π²
        Real expected = 5.0 * Constants::PI * Constants::PI;
        REQUIRE_THAT(torus.Volume() , RealApprox(expected));
    }
}

TEST_CASE("Torus3D::SurfaceArea", "[geometry][torus][surface]")
{
    SECTION("Torus with R=3, r=1")
    {
        Torus3D torus(3.0, 1.0);
        // SA = 4π²Rr = 4π² * 3 * 1 = 12π²
        Real expected = 12.0 * Constants::PI * Constants::PI;
        REQUIRE_THAT(torus.SurfaceArea() , RealApprox(expected));
    }

    SECTION("Torus with R=5, r=2")
    {
        Torus3D torus(5.0, 2.0);
        // SA = 4π² * 5 * 2 = 40π²
        Real expected = 40.0 * Constants::PI * Constants::PI;
        REQUIRE_THAT(torus.SurfaceArea() , RealApprox(expected));
    }

    SECTION("Torus centered at arbitrary point")
    {
        Pnt3Cart center(10.0, 20.0, 30.0);
        Torus3D torus(4.0, 1.5, center);
        // SA = 4π² * 4 * 1.5 = 24π²
        Real expected = 24.0 * Constants::PI * Constants::PI;
        REQUIRE_THAT(torus.SurfaceArea() , RealApprox(expected));
    }

    SECTION("Large torus")
    {
        Torus3D torus(20.0, 3.0);
        // SA = 4π² * 20 * 3 = 240π²
        Real expected = 240.0 * Constants::PI * Constants::PI;
        REQUIRE_THAT(torus.SurfaceArea() , RealApprox(expected));
    }
}

TEST_CASE("Torus3D::GetCenter", "[geometry][torus][center]")
{
    SECTION("Torus centered at origin")
    {
        Torus3D torus(5.0, 2.0);
        Pnt3Cart center = torus.GetCenter();
        REQUIRE_THAT(center.X() , RealApprox(0.0));
        REQUIRE_THAT(center.Y() , RealApprox(0.0));
        REQUIRE_THAT(center.Z() , RealApprox(0.0));
    }

    SECTION("Torus centered at arbitrary point")
    {
        Pnt3Cart expected(10.0, 20.0, 30.0);
        Torus3D torus(5.0, 2.0, expected);
        Pnt3Cart center = torus.GetCenter();
        REQUIRE_THAT(center.X() , RealApprox(10.0));
        REQUIRE_THAT(center.Y() , RealApprox(20.0));
        REQUIRE_THAT(center.Z() , RealApprox(30.0));
    }
}

TEST_CASE("Torus3D::GetBoundingBox", "[geometry][torus][bounding]")
{
    SECTION("Torus centered at origin with R=3, r=1")
    {
        Torus3D torus(3.0, 1.0);
        BoundingBox3D bbox = torus.GetBoundingBox();
        
        // Outer radius = R + r = 4, extends ±4 in XY, ±r in Z
        REQUIRE_THAT(bbox.Min().X() , RealApprox(-4.0));
        REQUIRE_THAT(bbox.Min().Y() , RealApprox(-4.0));
        REQUIRE_THAT(bbox.Min().Z() , RealApprox(-1.0));
        
        REQUIRE_THAT(bbox.Max().X() , RealApprox(4.0));
        REQUIRE_THAT(bbox.Max().Y() , RealApprox(4.0));
        REQUIRE_THAT(bbox.Max().Z() , RealApprox(1.0));
    }

    SECTION("Torus centered at arbitrary point")
    {
        Pnt3Cart center(10.0, 20.0, 30.0);
        Torus3D torus(5.0, 2.0, center);
        BoundingBox3D bbox = torus.GetBoundingBox();
        
        // Outer radius = 7, so XY extends ±7 from center, Z extends ±2
        REQUIRE_THAT(bbox.Min().X() , RealApprox(3.0));
        REQUIRE_THAT(bbox.Min().Y() , RealApprox(13.0));
        REQUIRE_THAT(bbox.Min().Z() , RealApprox(28.0));
        
        REQUIRE_THAT(bbox.Max().X() , RealApprox(17.0));
        REQUIRE_THAT(bbox.Max().Y() , RealApprox(27.0));
        REQUIRE_THAT(bbox.Max().Z() , RealApprox(32.0));
    }

    SECTION("Bounding box dimensions")
    {
        Torus3D torus(6.0, 2.0);
        BoundingBox3D bbox = torus.GetBoundingBox();
        Real width = bbox.Max().X() - bbox.Min().X();
        Real depth = bbox.Max().Y() - bbox.Min().Y();
        Real height = bbox.Max().Z() - bbox.Min().Z();
        
        REQUIRE_THAT(width , RealApprox(16.0));   // 2 * (R + r) = 2 * 8
        REQUIRE_THAT(depth , RealApprox(16.0));   // 2 * (R + r) = 2 * 8
        REQUIRE_THAT(height , RealApprox(4.0));   // 2 * r = 2 * 2
    }
}

TEST_CASE("Torus3D::GetBoundingSphere", "[geometry][torus][bounding]")
{
    SECTION("Bounding sphere for torus at origin")
    {
        Torus3D torus(3.0, 1.0);
        BoundingSphere3D bsphere = torus.GetBoundingSphere();
        
        REQUIRE_THAT(bsphere.Center().X() , RealApprox(0.0));
        REQUIRE_THAT(bsphere.Center().Y() , RealApprox(0.0));
        REQUIRE_THAT(bsphere.Center().Z() , RealApprox(0.0));
        
        // Radius = R + r = 4
        REQUIRE_THAT(bsphere.Radius() , RealApprox(4.0));
    }

    SECTION("Bounding sphere for torus at arbitrary center")
    {
        Pnt3Cart center(10.0, 20.0, 30.0);
        Torus3D torus(5.0, 2.0, center);
        BoundingSphere3D bsphere = torus.GetBoundingSphere();
        
        REQUIRE_THAT(bsphere.Center().X() , RealApprox(10.0));
        REQUIRE_THAT(bsphere.Center().Y() , RealApprox(20.0));
        REQUIRE_THAT(bsphere.Center().Z() , RealApprox(30.0));
        
        // Radius = R + r = 7
        REQUIRE_THAT(bsphere.Radius() , RealApprox(7.0));
    }

    SECTION("Bounding sphere contains outermost points")
    {
        Torus3D torus(4.0, 1.5);
        BoundingSphere3D bsphere = torus.GetBoundingSphere();
        
        // Outermost point in XY plane at distance R+r from center
        Pnt3Cart outerPoint(5.5, 0.0, 0.0);
        REQUIRE(bsphere.Contains(outerPoint));
        
        // Point at top of tube
        Pnt3Cart topPoint(4.0, 0.0, 1.5);
        REQUIRE(bsphere.Contains(topPoint));
    }
}

TEST_CASE("Torus3D::ToString", "[geometry][torus][string]")
{
    SECTION("Torus at origin")
    {
        Torus3D torus(5.0, 2.0);
        std::string str = torus.ToString();
        
        REQUIRE(str.find("Torus3D") != std::string::npos);
        REQUIRE(str.find("Center") != std::string::npos);
        REQUIRE(str.find("MajorRadius") != std::string::npos);
        REQUIRE(str.find("MinorRadius") != std::string::npos);
        REQUIRE(str.find("Volume") != std::string::npos);
        REQUIRE(str.find("SurfaceArea") != std::string::npos);
    }

    SECTION("Torus at arbitrary center")
    {
        Pnt3Cart center(10.0, 20.0, 30.0);
        Torus3D torus(3.0, 1.5, center);
        std::string str = torus.ToString();
        
        REQUIRE(str.find("3") != std::string::npos);   // Major radius
        REQUIRE(str.find("1.5") != std::string::npos); // Minor radius
    }
}

TEST_CASE("Torus3D::Getters", "[geometry][torus][getters]")
{
    SECTION("GetMajorRadius and GetMinorRadius")
    {
        Torus3D torus(7.5, 2.5);
        REQUIRE_THAT(torus.GetMajorRadius() , RealApprox(7.5));
        REQUIRE_THAT(torus.GetMinorRadius() , RealApprox(2.5));
    }

    SECTION("GetNumU and GetNumV with custom tessellation")
    {
        Torus3D torus(5.0, 2.0, 32, 16);
        REQUIRE(torus.GetNumU() == 32);
        REQUIRE(torus.GetNumV() == 16);
    }

    SECTION("GetNumU and GetNumV with default tessellation")
    {
        Torus3D torus(5.0, 2.0);
        REQUIRE(torus.GetNumU() == 20);
        REQUIRE(torus.GetNumV() == 12);
    }
}

TEST_CASE("Torus3D::IsInside", "[geometry][torus][inside]")
{
    SECTION("Point at center is NOT inside (hollow torus)")
    {
        Torus3D torus(5.0, 2.0);
        // Center of torus is hollow - distance from center to tube is R=5, but tube radius is only 2
        REQUIRE(torus.IsInside(Pnt3Cart(0.0, 0.0, 0.0)) == false);
    }

    SECTION("Point on major circle at center height is inside")
    {
        Torus3D torus(5.0, 2.0);
        // Point on major circle (distance R from center, z=0)
        REQUIRE(torus.IsInside(Pnt3Cart(5.0, 0.0, 0.0)) == true);
        REQUIRE(torus.IsInside(Pnt3Cart(0.0, 5.0, 0.0)) == true);
        REQUIRE(torus.IsInside(Pnt3Cart(-5.0, 0.0, 0.0)) == true);
    }

    SECTION("Points on outer edge of torus")
    {
        Torus3D torus(5.0, 2.0);
        // Outermost point at distance R+r = 7 from center
        REQUIRE(torus.IsInside(Pnt3Cart(7.0, 0.0, 0.0)) == true);
        REQUIRE(torus.IsInside(Pnt3Cart(0.0, 7.0, 0.0)) == true);
    }

    SECTION("Points on inner edge of torus")
    {
        Torus3D torus(5.0, 2.0);
        // Innermost point at distance R-r = 3 from center
        REQUIRE(torus.IsInside(Pnt3Cart(3.0, 0.0, 0.0)) == true);
        REQUIRE(torus.IsInside(Pnt3Cart(0.0, 3.0, 0.0)) == true);
    }

    SECTION("Points at top and bottom of tube")
    {
        Torus3D torus(5.0, 2.0);
        // Point directly above major circle
        REQUIRE(torus.IsInside(Pnt3Cart(5.0, 0.0, 2.0)) == true);
        // Point directly below major circle
        REQUIRE(torus.IsInside(Pnt3Cart(5.0, 0.0, -2.0)) == true);
    }

    SECTION("Points outside radially are not inside")
    {
        Torus3D torus(5.0, 2.0);
        // Beyond outer radius
        REQUIRE(torus.IsInside(Pnt3Cart(7.1, 0.0, 0.0)) == false);
        REQUIRE(torus.IsInside(Pnt3Cart(0.0, 7.1, 0.0)) == false);
    }

    SECTION("Points inside inner hole are not inside")
    {
        Torus3D torus(5.0, 2.0);
        // Inside the hole (distance < R-r)
        REQUIRE(torus.IsInside(Pnt3Cart(2.9, 0.0, 0.0)) == false);
        REQUIRE(torus.IsInside(Pnt3Cart(2.0, 0.0, 0.0)) == false);
        REQUIRE(torus.IsInside(Pnt3Cart(1.0, 0.0, 0.0)) == false);
    }

    SECTION("Points above/below torus height are not inside")
    {
        Torus3D torus(5.0, 2.0);
        REQUIRE(torus.IsInside(Pnt3Cart(5.0, 0.0, 2.1)) == false);
        REQUIRE(torus.IsInside(Pnt3Cart(5.0, 0.0, -2.1)) == false);
        REQUIRE(torus.IsInside(Pnt3Cart(5.0, 0.0, 3.0)) == false);
    }

    SECTION("Torus at arbitrary position")
    {
        Pnt3Cart center(10.0, 20.0, 30.0);
        Torus3D torus(5.0, 2.0, center);
        
        // Center is hollow
        REQUIRE(torus.IsInside(center) == false);
        
        // Point on major circle
        REQUIRE(torus.IsInside(Pnt3Cart(15.0, 20.0, 30.0)) == true);
        
        // Point on outer edge
        REQUIRE(torus.IsInside(Pnt3Cart(17.0, 20.0, 30.0)) == true);
        
        // Point outside
        REQUIRE(torus.IsInside(Pnt3Cart(17.1, 20.0, 30.0)) == false);
        
        // Point inside hole
        REQUIRE(torus.IsInside(Pnt3Cart(12.9, 20.0, 30.0)) == false);
    }

    SECTION("Diagonal points in tube cross-section")
    {
        Torus3D torus(5.0, 2.0);
        // Point at 45° in tube cross-section: distance from major circle = √2 ≈ 1.414 < 2
        Real offset = std::sqrt(2.0) / 2.0;  // ≈ 0.707
        REQUIRE(torus.IsInside(Pnt3Cart(5.0 + offset, 0.0, offset)) == true);
    }
}

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

TEST_CASE("Cylinder3D::Volume", "[geometry][cylinder][volume]")
{
    SECTION("Unit cylinder - radius 1, height 1")
    {
        Cylinder3D cyl(1.0, 1.0);
        Real expected = Constants::PI;  // π * 1² * 1 = π
        REQUIRE_THAT(cyl.Volume() , RealApprox(expected));
    }

    SECTION("Cylinder with radius 2, height 5")
    {
        Cylinder3D cyl(2.0, 5.0);
        Real expected = Constants::PI * 4.0 * 5.0;  // π * 4 * 5 = 20π
        REQUIRE_THAT(cyl.Volume() , RealApprox(expected));
    }

    SECTION("Cylinder centered at arbitrary point")
    {
        Pnt3Cart base(10.0, 20.0, 30.0);
        Cylinder3D cyl(3.0, 10.0, base);
        Real expected = Constants::PI * 9.0 * 10.0;  // π * 9 * 10 = 90π
        REQUIRE_THAT(cyl.Volume() , RealApprox(expected));
    }

    SECTION("Tall thin cylinder")
    {
        Cylinder3D cyl(0.5, 20.0);
        Real expected = Constants::PI * 0.25 * 20.0;  // π * 0.25 * 20 = 5π
        REQUIRE_THAT(cyl.Volume() , RealApprox(expected));
    }
}

TEST_CASE("Cylinder3D::SurfaceArea", "[geometry][cylinder][surface]")
{
    SECTION("Unit cylinder - radius 1, height 1")
    {
        Cylinder3D cyl(1.0, 1.0);
        // SA = 2πR² + 2πRH = 2π(1) + 2π(1)(1) = 2π + 2π = 4π
        Real expected = 4.0 * Constants::PI;
        REQUIRE_THAT(cyl.SurfaceArea() , RealApprox(expected));
    }

    SECTION("Cylinder with radius 2, height 5")
    {
        Cylinder3D cyl(2.0, 5.0);
        // SA = 2π(4) + 2π(2)(5) = 8π + 20π = 28π
        Real expected = 28.0 * Constants::PI;
        REQUIRE_THAT(cyl.SurfaceArea() , RealApprox(expected));
    }

    SECTION("Cylinder centered at arbitrary point")
    {
        Pnt3Cart base(10.0, 20.0, 30.0);
        Cylinder3D cyl(3.0, 10.0, base);
        // SA = 2π(9) + 2π(3)(10) = 18π + 60π = 78π
        Real expected = 78.0 * Constants::PI;
        REQUIRE_THAT(cyl.SurfaceArea() , RealApprox(expected));
    }

    SECTION("Wide flat cylinder")
    {
        Cylinder3D cyl(10.0, 1.0);
        // SA = 2π(100) + 2π(10)(1) = 200π + 20π = 220π
        Real expected = 220.0 * Constants::PI;
        REQUIRE_THAT(cyl.SurfaceArea() , RealApprox(expected));
    }
}

TEST_CASE("Cylinder3D::GetCenter", "[geometry][cylinder][center]")
{
    SECTION("Cylinder with base at origin")
    {
        Cylinder3D cyl(5.0, 10.0);
        Pnt3Cart center = cyl.GetCenter();
        // Base at (0,0,-5), top at (0,0,5), center at (0,0,0)
        REQUIRE_THAT(center.X() , RealApprox(0.0));
        REQUIRE_THAT(center.Y() , RealApprox(0.0));
        REQUIRE_THAT(center.Z() , RealApprox(0.0));
    }

    SECTION("Cylinder with custom base position")
    {
        Pnt3Cart base(10.0, 20.0, 30.0);
        Cylinder3D cyl(5.0, 20.0, base);
        Pnt3Cart center = cyl.GetCenter();
        // Base at (10,20,30), top at (10,20,50), center at (10,20,40)
        REQUIRE_THAT(center.X() , RealApprox(10.0));
        REQUIRE_THAT(center.Y() , RealApprox(20.0));
        REQUIRE_THAT(center.Z() , RealApprox(40.0));
    }
}

TEST_CASE("Cylinder3D::GetBoundingBox", "[geometry][cylinder][bounding]")
{
    SECTION("Cylinder with base at origin")
    {
        Cylinder3D cyl(5.0, 10.0);
        BoundingBox3D bbox = cyl.GetBoundingBox();
        
        // Base at z=-5, so cylinder goes from z=-5 to z=5
        // But constructor sets _center to (0,0,-H/2), so actually z=0 to z=10
        REQUIRE_THAT(bbox.Min().X() , RealApprox(-5.0));
        REQUIRE_THAT(bbox.Min().Y() , RealApprox(-5.0));
        REQUIRE_THAT(bbox.Min().Z() , RealApprox(-5.0));
        
        REQUIRE_THAT(bbox.Max().X() , RealApprox(5.0));
        REQUIRE_THAT(bbox.Max().Y() , RealApprox(5.0));
        REQUIRE_THAT(bbox.Max().Z() , RealApprox(5.0));
    }

    SECTION("Cylinder with custom base position")
    {
        Pnt3Cart base(10.0, 20.0, 30.0);
        Cylinder3D cyl(3.0, 12.0, base);
        BoundingBox3D bbox = cyl.GetBoundingBox();
        
        REQUIRE_THAT(bbox.Min().X() , RealApprox(7.0));
        REQUIRE_THAT(bbox.Min().Y() , RealApprox(17.0));
        REQUIRE_THAT(bbox.Min().Z() , RealApprox(30.0));
        
        REQUIRE_THAT(bbox.Max().X() , RealApprox(13.0));
        REQUIRE_THAT(bbox.Max().Y() , RealApprox(23.0));
        REQUIRE_THAT(bbox.Max().Z() , RealApprox(42.0));
    }

    SECTION("Bounding box dimensions")
    {
        Cylinder3D cyl(4.0, 8.0);
        BoundingBox3D bbox = cyl.GetBoundingBox();
        Real width = bbox.Max().X() - bbox.Min().X();
        Real depth = bbox.Max().Y() - bbox.Min().Y();
        Real height = bbox.Max().Z() - bbox.Min().Z();
        
        REQUIRE_THAT(width , RealApprox(8.0));   // 2 * radius
        REQUIRE_THAT(depth , RealApprox(8.0));   // 2 * radius
        REQUIRE_THAT(height , RealApprox(8.0));  // height
    }
}

TEST_CASE("Cylinder3D::GetBoundingSphere", "[geometry][cylinder][bounding]")
{
    SECTION("Bounding sphere for symmetric cylinder")
    {
        Cylinder3D cyl(3.0, 8.0);
        BoundingSphere3D bsphere = cyl.GetBoundingSphere();
        
        // Geometric center at (0, 0, 0) since base is at (0,0,-4)
        REQUIRE_THAT(bsphere.Center().X() , RealApprox(0.0));
        REQUIRE_THAT(bsphere.Center().Y() , RealApprox(0.0));
        REQUIRE_THAT(bsphere.Center().Z() , RealApprox(0.0));
        
        // Radius = sqrt(R² + (H/2)²) = sqrt(9 + 16) = sqrt(25) = 5
        Real expected_radius = std::sqrt(9.0 + 16.0);
        REQUIRE_THAT(bsphere.Radius() , RealApprox(expected_radius));
    }

    SECTION("Bounding sphere for tall cylinder")
    {
        Pnt3Cart base(0.0, 0.0, 0.0);
        Cylinder3D cyl(2.0, 10.0, base);
        BoundingSphere3D bsphere = cyl.GetBoundingSphere();
        
        // Center at (0, 0, 5)
        REQUIRE_THAT(bsphere.Center().Z() , RealApprox(5.0));
        
        // Radius = sqrt(4 + 25) = sqrt(29) ≈ 5.385
        Real expected_radius = std::sqrt(4.0 + 25.0);
        REQUIRE_THAT(bsphere.Radius() , RealApprox(expected_radius));
    }

    SECTION("Bounding sphere contains cylinder corners")
    {
        Cylinder3D cyl(3.0, 6.0);
        BoundingSphere3D bsphere = cyl.GetBoundingSphere();
        
        // Check that top edge points are contained
        Pnt3Cart topEdge(3.0, 0.0, 3.0);  // Top of cylinder at max radius
        REQUIRE(bsphere.Contains(topEdge));
        
        // Check that bottom edge points are contained
        Pnt3Cart bottomEdge(3.0, 0.0, -3.0);  // Bottom of cylinder at max radius
        REQUIRE(bsphere.Contains(bottomEdge));
    }
}

TEST_CASE("Cylinder3D::ToString", "[geometry][cylinder][string]")
{
    SECTION("Cylinder at default position")
    {
        Cylinder3D cyl(5.0, 10.0);
        std::string str = cyl.ToString();
        
        REQUIRE(str.find("Cylinder3D") != std::string::npos);
        REQUIRE(str.find("Center") != std::string::npos);
        REQUIRE(str.find("Radius") != std::string::npos);
        REQUIRE(str.find("Height") != std::string::npos);
        REQUIRE(str.find("Volume") != std::string::npos);
        REQUIRE(str.find("SurfaceArea") != std::string::npos);
    }

    SECTION("Cylinder at custom position")
    {
        Pnt3Cart base(10.0, 20.0, 30.0);
        Cylinder3D cyl(3.0, 12.0, base);
        std::string str = cyl.ToString();
        
        REQUIRE(str.find("3") != std::string::npos);  // Radius
        REQUIRE(str.find("12") != std::string::npos); // Height
    }
}

TEST_CASE("Cylinder3D::Getters", "[geometry][cylinder][getters]")
{
    SECTION("GetRadius and GetHeight")
    {
        Cylinder3D cyl(7.5, 15.0);
        REQUIRE_THAT(cyl.GetRadius() , RealApprox(7.5));
        REQUIRE_THAT(cyl.GetHeight() , RealApprox(15.0));
    }

    SECTION("GetBaseCenter")
    {
        Pnt3Cart base(10.0, 20.0, 30.0);
        Cylinder3D cyl(5.0, 10.0, base);
        Pnt3Cart baseCenter = cyl.GetBaseCenter();
        
        REQUIRE_THAT(baseCenter.X() , RealApprox(10.0));
        REQUIRE_THAT(baseCenter.Y() , RealApprox(20.0));
        REQUIRE_THAT(baseCenter.Z() , RealApprox(30.0));
    }

    SECTION("GetNumSegments")
    {
        Cylinder3D cyl1(5.0, 10.0, 32);
        REQUIRE(cyl1.GetNumSegments() == 32);
        
        Cylinder3D cyl2(5.0, 10.0);  // Default segments
        REQUIRE(cyl2.GetNumSegments() == 20);
    }
}

TEST_CASE("Cylinder3D::IsInside", "[geometry][cylinder][inside]")
{
    SECTION("Point at geometric center is inside")
    {
        Cylinder3D cyl(5.0, 10.0);
        REQUIRE(cyl.IsInside(Pnt3Cart(0.0, 0.0, 0.0)) == true);
    }

    SECTION("Point at base center is inside")
    {
        Pnt3Cart base(0.0, 0.0, 0.0);
        Cylinder3D cyl(5.0, 10.0, base);
        REQUIRE(cyl.IsInside(base) == true);
    }

    SECTION("Point on axis at various heights")
    {
        Pnt3Cart base(0.0, 0.0, 0.0);
        Cylinder3D cyl(5.0, 10.0, base);
        
        REQUIRE(cyl.IsInside(Pnt3Cart(0.0, 0.0, 0.0)) == true);   // Bottom
        REQUIRE(cyl.IsInside(Pnt3Cart(0.0, 0.0, 5.0)) == true);   // Middle
        REQUIRE(cyl.IsInside(Pnt3Cart(0.0, 0.0, 10.0)) == true);  // Top
    }

    SECTION("Points on circular edge at various heights")
    {
        Pnt3Cart base(0.0, 0.0, 0.0);
        Cylinder3D cyl(5.0, 10.0, base);
        
        REQUIRE(cyl.IsInside(Pnt3Cart(5.0, 0.0, 0.0)) == true);   // Bottom edge
        REQUIRE(cyl.IsInside(Pnt3Cart(5.0, 0.0, 5.0)) == true);   // Middle edge
        REQUIRE(cyl.IsInside(Pnt3Cart(5.0, 0.0, 10.0)) == true);  // Top edge
        REQUIRE(cyl.IsInside(Pnt3Cart(0.0, 5.0, 5.0)) == true);   // Different angle
    }

    SECTION("Points outside radius are not inside")
    {
        Pnt3Cart base(0.0, 0.0, 0.0);
        Cylinder3D cyl(5.0, 10.0, base);
        
        REQUIRE(cyl.IsInside(Pnt3Cart(5.1, 0.0, 5.0)) == false);
        REQUIRE(cyl.IsInside(Pnt3Cart(0.0, 5.1, 5.0)) == false);
        REQUIRE(cyl.IsInside(Pnt3Cart(10.0, 0.0, 5.0)) == false);
    }

    SECTION("Points outside height range are not inside")
    {
        Pnt3Cart base(0.0, 0.0, 0.0);
        Cylinder3D cyl(5.0, 10.0, base);
        
        REQUIRE(cyl.IsInside(Pnt3Cart(0.0, 0.0, -0.1)) == false);   // Below base
        REQUIRE(cyl.IsInside(Pnt3Cart(0.0, 0.0, 10.1)) == false);   // Above top
        REQUIRE(cyl.IsInside(Pnt3Cart(2.0, 0.0, -1.0)) == false);   // Below and within radius
        REQUIRE(cyl.IsInside(Pnt3Cart(2.0, 0.0, 11.0)) == false);   // Above and within radius
    }

    SECTION("Cylinder at arbitrary position")
    {
        Pnt3Cart base(10.0, 20.0, 30.0);
        Cylinder3D cyl(3.0, 12.0, base);
        
        REQUIRE(cyl.IsInside(base) == true);                                  // Base center
        REQUIRE(cyl.IsInside(Pnt3Cart(10.0, 20.0, 36.0)) == true);           // Geometric center
        REQUIRE(cyl.IsInside(Pnt3Cart(13.0, 20.0, 36.0)) == true);           // On edge
        REQUIRE(cyl.IsInside(Pnt3Cart(10.0, 20.0, 42.0)) == true);           // Top center
        
        REQUIRE(cyl.IsInside(Pnt3Cart(14.0, 20.0, 36.0)) == false);          // Outside radius
        REQUIRE(cyl.IsInside(Pnt3Cart(10.0, 20.0, 29.9)) == false);          // Below base
        REQUIRE(cyl.IsInside(Pnt3Cart(10.0, 20.0, 42.1)) == false);          // Above top
    }

    SECTION("3-4-5 right triangle test")
    {
        Pnt3Cart base(0.0, 0.0, 0.0);
        Cylinder3D cyl(5.0, 10.0, base);
        
        // Point at (3, 4, 5): distance from axis = 5.0 (exactly on edge)
        REQUIRE(cyl.IsInside(Pnt3Cart(3.0, 4.0, 5.0)) == true);
        
        // Point at (3, 4.1, 5): distance from axis > 5.0
        REQUIRE(cyl.IsInside(Pnt3Cart(3.0, 4.1, 5.0)) == false);
    }
}

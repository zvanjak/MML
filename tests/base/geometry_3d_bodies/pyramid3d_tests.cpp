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

TEST_CASE("Pyramid3D::Volume", "[geometry][pyramid][volume]")
{
    SECTION("Pyramid with base 6, height 4")
    {
        Pyramid3D pyramid(6.0, 4.0);
        // Volume = (1/3) * a² * h = (1/3) * 36 * 4 = 48
        Real expected = 48.0;
        REQUIRE_THAT(pyramid.Volume() , RealApprox(expected));
    }

    SECTION("Pyramid with base 3, height 9")
    {
        Pyramid3D pyramid(3.0, 9.0);
        // Volume = (1/3) * 9 * 9 = 27
        Real expected = 27.0;
        REQUIRE_THAT(pyramid.Volume() , RealApprox(expected));
    }

    SECTION("Unit pyramid - base 1, height 1")
    {
        Pyramid3D pyramid(1.0, 1.0);
        // Volume = (1/3) * 1 * 1 = 1/3
        Real expected = 1.0 / 3.0;
        REQUIRE_THAT(pyramid.Volume() , RealApprox(expected));
    }

    SECTION("Pyramid at arbitrary position")
    {
        Vec3Cart base(10.0, 20.0, 30.0);
        Pyramid3D pyramid(4.0, 6.0, base);
        // Volume = (1/3) * 16 * 6 = 32
        Real expected = 32.0;
        REQUIRE_THAT(pyramid.Volume() , RealApprox(expected));
    }
}

TEST_CASE("Pyramid3D::SurfaceArea", "[geometry][pyramid][surface]")
{
    SECTION("Pyramid with base 6, height 4")
    {
        Pyramid3D pyramid(6.0, 4.0);
        // Base = 36
        // Slant height s = sqrt(16 + 9) = sqrt(25) = 5
        // SA = 36 + 2*6*5 = 36 + 60 = 96
        Real expected = 96.0;
        REQUIRE_THAT(pyramid.SurfaceArea() , RealApprox(expected));
    }

    SECTION("Pyramid with base 8, height 3")
    {
        Pyramid3D pyramid(8.0, 3.0);
        // Base = 64
        // Slant height s = sqrt(9 + 16) = 5
        // SA = 64 + 2*8*5 = 64 + 80 = 144
        Real expected = 144.0;
        REQUIRE_THAT(pyramid.SurfaceArea() , RealApprox(expected));
    }

    SECTION("Unit pyramid")
    {
        Pyramid3D pyramid(1.0, 1.0);
        // Base = 1
        // Slant height s = sqrt(1 + 0.25) = sqrt(1.25)
        // SA = 1 + 2*sqrt(1.25)
        Real slant = std::sqrt(1.25);
        Real expected = 1.0 + 2.0 * slant;
        REQUIRE_THAT(pyramid.SurfaceArea() , RealApprox(expected));
    }

    SECTION("3-4-5 right triangle pyramid")
    {
        Pyramid3D pyramid(6.0, 4.0);
        // Base side = 6, height = 4
        // Half base = 3, so slant height = sqrt(4² + 3²) = 5
        Real slant_height = std::sqrt(16.0 + 9.0);
        REQUIRE_THAT(slant_height , RealApprox(5.0));
    }
}

TEST_CASE("Pyramid3D::GetCenter", "[geometry][pyramid][center]")
{
    SECTION("Pyramid at origin")
    {
        Pyramid3D pyramid(6.0, 12.0);
        Pnt3Cart center = pyramid.GetCenter();
        // Centroid at h/4 = 3 from base
        REQUIRE_THAT(center.X() , RealApprox(0.0));
        REQUIRE_THAT(center.Y() , RealApprox(0.0));
        REQUIRE_THAT(center.Z() , RealApprox(3.0));
    }

    SECTION("Pyramid at arbitrary position")
    {
        Vec3Cart base(10.0, 20.0, 30.0);
        Pyramid3D pyramid(8.0, 20.0, base);
        Pnt3Cart center = pyramid.GetCenter();
        // Centroid at base_z + h/4 = 30 + 5 = 35
        REQUIRE_THAT(center.X() , RealApprox(10.0));
        REQUIRE_THAT(center.Y() , RealApprox(20.0));
        REQUIRE_THAT(center.Z() , RealApprox(35.0));
    }

    SECTION("Centroid is at 1/4 height")
    {
        Pyramid3D pyramid(10.0, 8.0);
        Pnt3Cart center = pyramid.GetCenter();
        REQUIRE_THAT(center.Z() , RealApprox(2.0));  // 8/4 = 2
    }
}

TEST_CASE("Pyramid3D::GetBoundingBox", "[geometry][pyramid][bounding]")
{
    SECTION("Pyramid at origin")
    {
        Pyramid3D pyramid(10.0, 8.0);
        BoundingBox3D bbox = pyramid.GetBoundingBox();
        
        REQUIRE_THAT(bbox.Min().X() , RealApprox(-5.0));
        REQUIRE_THAT(bbox.Min().Y() , RealApprox(-5.0));
        REQUIRE_THAT(bbox.Min().Z() , RealApprox(0.0));
        
        REQUIRE_THAT(bbox.Max().X() , RealApprox(5.0));
        REQUIRE_THAT(bbox.Max().Y() , RealApprox(5.0));
        REQUIRE_THAT(bbox.Max().Z() , RealApprox(8.0));
    }

    SECTION("Pyramid at arbitrary position")
    {
        Vec3Cart base(10.0, 20.0, 30.0);
        Pyramid3D pyramid(6.0, 12.0, base);
        BoundingBox3D bbox = pyramid.GetBoundingBox();
        
        REQUIRE_THAT(bbox.Min().X() , RealApprox(7.0));
        REQUIRE_THAT(bbox.Min().Y() , RealApprox(17.0));
        REQUIRE_THAT(bbox.Min().Z() , RealApprox(30.0));
        
        REQUIRE_THAT(bbox.Max().X() , RealApprox(13.0));
        REQUIRE_THAT(bbox.Max().Y() , RealApprox(23.0));
        REQUIRE_THAT(bbox.Max().Z() , RealApprox(42.0));
    }

    SECTION("Bounding box dimensions")
    {
        Pyramid3D pyramid(8.0, 10.0);
        BoundingBox3D bbox = pyramid.GetBoundingBox();
        Real width = bbox.Max().X() - bbox.Min().X();
        Real depth = bbox.Max().Y() - bbox.Min().Y();
        Real height = bbox.Max().Z() - bbox.Min().Z();
        
        REQUIRE_THAT(width , RealApprox(8.0));   // Base side
        REQUIRE_THAT(depth , RealApprox(8.0));   // Base side
        REQUIRE_THAT(height , RealApprox(10.0)); // Pyramid height
    }
}

TEST_CASE("Pyramid3D::GetBoundingSphere", "[geometry][pyramid][bounding]")
{
    SECTION("Bounding sphere centered at centroid")
    {
        Pyramid3D pyramid(6.0, 8.0);
        BoundingSphere3D bsphere = pyramid.GetBoundingSphere();
        
        // Centroid at (0, 0, 2) [h/4 from base]
        REQUIRE_THAT(bsphere.Center().X() , RealApprox(0.0));
        REQUIRE_THAT(bsphere.Center().Y() , RealApprox(0.0));
        REQUIRE_THAT(bsphere.Center().Z() , RealApprox(2.0));
        
        // Radius is max of:
        // - Distance to base corner (3,3,0): sqrt(9 + 9 + 4) = sqrt(22) ≈ 4.69
        // - Distance to apex (0,0,8): 8 - 2 = 6.0
        Real distToApex = 6.0;  // 3h/4 = 3*8/4 = 6
        REQUIRE_THAT(bsphere.Radius() , RealApprox(distToApex));
    }

    SECTION("Bounding sphere contains all vertices")
    {
        Pyramid3D pyramid(10.0, 12.0);
        BoundingSphere3D bsphere = pyramid.GetBoundingSphere();
        
        // Check base corners
        REQUIRE(bsphere.Contains(Pnt3Cart(5.0, 5.0, 0.0)));
        REQUIRE(bsphere.Contains(Pnt3Cart(-5.0, -5.0, 0.0)));
        
        // Check apex
        REQUIRE(bsphere.Contains(Pnt3Cart(0.0, 0.0, 12.0)));
    }

    SECTION("Bounding sphere at arbitrary position")
    {
        Vec3Cart base(10.0, 20.0, 30.0);
        Pyramid3D pyramid(8.0, 16.0, base);
        BoundingSphere3D bsphere = pyramid.GetBoundingSphere();
        
        // Centroid at (10, 20, 34)
        REQUIRE_THAT(bsphere.Center().X() , RealApprox(10.0));
        REQUIRE_THAT(bsphere.Center().Y() , RealApprox(20.0));
        REQUIRE_THAT(bsphere.Center().Z() , RealApprox(34.0));
    }
}

TEST_CASE("Pyramid3D::ToString", "[geometry][pyramid][string]")
{
    SECTION("Pyramid at origin")
    {
        Pyramid3D pyramid(6.0, 8.0);
        std::string str = pyramid.ToString();
        
        REQUIRE(str.find("Pyramid3D") != std::string::npos);
        REQUIRE(str.find("Center") != std::string::npos);
        REQUIRE(str.find("BaseSize") != std::string::npos);
        REQUIRE(str.find("Height") != std::string::npos);
        REQUIRE(str.find("Volume") != std::string::npos);
        REQUIRE(str.find("SurfaceArea") != std::string::npos);
    }

    SECTION("Pyramid at arbitrary position")
    {
        Vec3Cart base(10.0, 20.0, 30.0);
        Pyramid3D pyramid(4.0, 6.0, base);
        std::string str = pyramid.ToString();
        
        REQUIRE(str.find("4") != std::string::npos);  // Base size
        REQUIRE(str.find("6") != std::string::npos);  // Height
    }
}

TEST_CASE("Pyramid3D::Getters", "[geometry][pyramid][getters]")
{
    SECTION("GetBaseSize and GetHeight")
    {
        Pyramid3D pyramid(7.5, 15.0);
        REQUIRE_THAT(pyramid.GetBaseSize() , RealApprox(7.5));
        REQUIRE_THAT(pyramid.GetHeight() , RealApprox(15.0));
    }

    SECTION("GetBaseCenter")
    {
        Vec3Cart base(10.0, 20.0, 30.0);
        Pyramid3D pyramid(5.0, 10.0, base);
        Vec3Cart baseCenter = pyramid.GetBaseCenter();
        
        REQUIRE_THAT(baseCenter.X() , RealApprox(10.0));
        REQUIRE_THAT(baseCenter.Y() , RealApprox(20.0));
        REQUIRE_THAT(baseCenter.Z() , RealApprox(30.0));
    }
}

TEST_CASE("Pyramid3D::IsInside", "[geometry][pyramid][inside]")
{
    SECTION("Point at base center is inside")
    {
        Pyramid3D pyramid(10.0, 12.0);
        REQUIRE(pyramid.IsInside(Pnt3Cart(0.0, 0.0, 0.0)) == true);
    }

    SECTION("Point at geometric center is inside")
    {
        Pyramid3D pyramid(10.0, 12.0);
        REQUIRE(pyramid.IsInside(Pnt3Cart(0.0, 0.0, 3.0)) == true);  // At h/4
    }

    SECTION("Apex is inside")
    {
        Pyramid3D pyramid(10.0, 12.0);
        REQUIRE(pyramid.IsInside(Pnt3Cart(0.0, 0.0, 12.0)) == true);
    }

    SECTION("Points on base edges are inside")
    {
        Pyramid3D pyramid(10.0, 12.0);
        REQUIRE(pyramid.IsInside(Pnt3Cart(5.0, 0.0, 0.0)) == true);    // Base edge
        REQUIRE(pyramid.IsInside(Pnt3Cart(-5.0, 0.0, 0.0)) == true);
        REQUIRE(pyramid.IsInside(Pnt3Cart(0.0, 5.0, 0.0)) == true);
        REQUIRE(pyramid.IsInside(Pnt3Cart(0.0, -5.0, 0.0)) == true);
    }

    SECTION("Points on base corners are inside")
    {
        Pyramid3D pyramid(10.0, 12.0);
        REQUIRE(pyramid.IsInside(Pnt3Cart(5.0, 5.0, 0.0)) == true);
        REQUIRE(pyramid.IsInside(Pnt3Cart(-5.0, -5.0, 0.0)) == true);
        REQUIRE(pyramid.IsInside(Pnt3Cart(5.0, -5.0, 0.0)) == true);
        REQUIRE(pyramid.IsInside(Pnt3Cart(-5.0, 5.0, 0.0)) == true);
    }

    SECTION("Points on pyramid surface at mid-height")
    {
        Pyramid3D pyramid(10.0, 10.0);
        // At z=5 (mid-height), cross-section half-side = 5 * (1 - 5/10) = 2.5
        REQUIRE(pyramid.IsInside(Pnt3Cart(2.5, 0.0, 5.0)) == true);
        REQUIRE(pyramid.IsInside(Pnt3Cart(-2.5, 0.0, 5.0)) == true);
        REQUIRE(pyramid.IsInside(Pnt3Cart(0.0, 2.5, 5.0)) == true);
        REQUIRE(pyramid.IsInside(Pnt3Cart(2.5, 2.5, 5.0)) == true);
    }

    SECTION("Points below base are not inside")
    {
        Pyramid3D pyramid(10.0, 12.0);
        REQUIRE(pyramid.IsInside(Pnt3Cart(0.0, 0.0, -0.1)) == false);
        REQUIRE(pyramid.IsInside(Pnt3Cart(0.0, 0.0, -1.0)) == false);
    }

    SECTION("Points above apex are not inside")
    {
        Pyramid3D pyramid(10.0, 12.0);
        REQUIRE(pyramid.IsInside(Pnt3Cart(0.0, 0.0, 12.1)) == false);
        REQUIRE(pyramid.IsInside(Pnt3Cart(0.0, 0.0, 13.0)) == false);
    }

    SECTION("Points outside base are not inside")
    {
        Pyramid3D pyramid(10.0, 12.0);
        REQUIRE(pyramid.IsInside(Pnt3Cart(5.1, 0.0, 0.0)) == false);
        REQUIRE(pyramid.IsInside(Pnt3Cart(0.0, 5.1, 0.0)) == false);
        REQUIRE(pyramid.IsInside(Pnt3Cart(6.0, 6.0, 0.0)) == false);
    }

    SECTION("Points outside slanted sides are not inside")
    {
        Pyramid3D pyramid(10.0, 10.0);
        // At z=5, half-side = 2.5, so x or y > 2.5 should be outside
        REQUIRE(pyramid.IsInside(Pnt3Cart(3.0, 0.0, 5.0)) == false);
        REQUIRE(pyramid.IsInside(Pnt3Cart(0.0, 3.0, 5.0)) == false);
        REQUIRE(pyramid.IsInside(Pnt3Cart(2.6, 2.6, 5.0)) == false);
    }

    SECTION("Pyramid at arbitrary position")
    {
        Vec3Cart base(10.0, 20.0, 30.0);
        Pyramid3D pyramid(8.0, 12.0, base);
        
        // Base center
        REQUIRE(pyramid.IsInside(Pnt3Cart(10.0, 20.0, 30.0)) == true);
        
        // Apex
        REQUIRE(pyramid.IsInside(Pnt3Cart(10.0, 20.0, 42.0)) == true);
        
        // Base corner
        REQUIRE(pyramid.IsInside(Pnt3Cart(14.0, 24.0, 30.0)) == true);
        
        // Outside
        REQUIRE(pyramid.IsInside(Pnt3Cart(15.0, 20.0, 30.0)) == false);
        REQUIRE(pyramid.IsInside(Pnt3Cart(10.0, 20.0, 29.9)) == false);
        REQUIRE(pyramid.IsInside(Pnt3Cart(10.0, 20.0, 42.1)) == false);
    }

    SECTION("Linear cross-section reduction")
    {
        Pyramid3D pyramid(10.0, 10.0);
        // At z=0 (base): half-side = 5
        // At z=2.5: half-side = 5 * (1 - 0.25) = 3.75
        // At z=5: half-side = 5 * 0.5 = 2.5
        // At z=7.5: half-side = 5 * 0.25 = 1.25
        // At z=10: half-side = 0
        
        REQUIRE(pyramid.IsInside(Pnt3Cart(3.75, 0.0, 2.5)) == true);
        REQUIRE(pyramid.IsInside(Pnt3Cart(2.5, 0.0, 5.0)) == true);
        REQUIRE(pyramid.IsInside(Pnt3Cart(1.25, 0.0, 7.5)) == true);
        
        REQUIRE(pyramid.IsInside(Pnt3Cart(3.76, 0.0, 2.5)) == false);
        REQUIRE(pyramid.IsInside(Pnt3Cart(2.51, 0.0, 5.0)) == false);
    }
}

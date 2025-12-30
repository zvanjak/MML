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

TEST_CASE("Sphere3D::Volume", "[geometry][sphere][volume]")
{
    SECTION("Unit sphere centered at origin")
    {
        Sphere3D sphere(1.0);
        Real expected = (4.0 / 3.0) * Constants::PI;  // (4/3)π ≈ 4.1888
        REQUIRE_THAT(sphere.Volume() , RealApprox(expected));
    }

    SECTION("Sphere with radius 5 centered at origin")
    {
        Sphere3D sphere(5.0);
        Real expected = (4.0 / 3.0) * Constants::PI * 125.0;  // (4/3)π·125 ≈ 523.599
        REQUIRE_THAT(sphere.Volume() , RealApprox(expected));
    }

    SECTION("Sphere centered at arbitrary point")
    {
        Pnt3Cart center(10.0, 20.0, 30.0);
        Sphere3D sphere(3.0, center);
        Real expected = (4.0 / 3.0) * Constants::PI * 27.0;  // (4/3)π·27 ≈ 113.097
        REQUIRE_THAT(sphere.Volume() , RealApprox(expected));
    }

    SECTION("Small sphere")
    {
        Sphere3D sphere(0.1);
        Real expected = (4.0 / 3.0) * Constants::PI * 0.001;  // (4/3)π·0.001 ≈ 0.004189
        REQUIRE_THAT(sphere.Volume() , RealApprox(expected));
    }
}

TEST_CASE("Sphere3D::SurfaceArea", "[geometry][sphere][surface]")
{
    SECTION("Unit sphere centered at origin")
    {
        Sphere3D sphere(1.0);
        Real expected = 4.0 * Constants::PI;  // 4π ≈ 12.566
        REQUIRE_THAT(sphere.SurfaceArea() , RealApprox(expected));
    }

    SECTION("Sphere with radius 5 centered at origin")
    {
        Sphere3D sphere(5.0);
        Real expected = 4.0 * Constants::PI * 25.0;  // 4π·25 ≈ 314.159
        REQUIRE_THAT(sphere.SurfaceArea() , RealApprox(expected));
    }

    SECTION("Sphere centered at arbitrary point")
    {
        Pnt3Cart center(10.0, 20.0, 30.0);
        Sphere3D sphere(3.0, center);
        Real expected = 4.0 * Constants::PI * 9.0;  // 4π·9 ≈ 113.097
        REQUIRE_THAT(sphere.SurfaceArea() , RealApprox(expected));
    }

    SECTION("Large sphere")
    {
        Sphere3D sphere(100.0);
        Real expected = 4.0 * Constants::PI * 10000.0;  // 4π·10000 ≈ 125663.7
        REQUIRE_THAT(sphere.SurfaceArea() , RealApprox(expected));
    }
}

TEST_CASE("Sphere3D::GetCenter", "[geometry][sphere][center]")
{
    SECTION("Sphere centered at origin")
    {
        Sphere3D sphere(5.0);
        Pnt3Cart center = sphere.GetCenter();
        REQUIRE_THAT(center.X() , RealApprox(0.0));
        REQUIRE_THAT(center.Y() , RealApprox(0.0));
        REQUIRE_THAT(center.Z() , RealApprox(0.0));
    }

    SECTION("Sphere centered at arbitrary point")
    {
        Pnt3Cart expected(10.0, 20.0, 30.0);
        Sphere3D sphere(5.0, expected);
        Pnt3Cart center = sphere.GetCenter();
        REQUIRE_THAT(center.X() , RealApprox(10.0));
        REQUIRE_THAT(center.Y() , RealApprox(20.0));
        REQUIRE_THAT(center.Z() , RealApprox(30.0));
    }
}

TEST_CASE("Sphere3D::GetBoundingBox", "[geometry][sphere][bounding]")
{
    SECTION("Unit sphere centered at origin")
    {
        Sphere3D sphere(1.0);
        BoundingBox3D bbox = sphere.GetBoundingBox();
        
        REQUIRE_THAT(bbox.Min().X() , RealApprox(-1.0));
        REQUIRE_THAT(bbox.Min().Y() , RealApprox(-1.0));
        REQUIRE_THAT(bbox.Min().Z() , RealApprox(-1.0));
        
        REQUIRE_THAT(bbox.Max().X() , RealApprox(1.0));
        REQUIRE_THAT(bbox.Max().Y() , RealApprox(1.0));
        REQUIRE_THAT(bbox.Max().Z() , RealApprox(1.0));
    }

    SECTION("Sphere with radius 5 centered at origin")
    {
        Sphere3D sphere(5.0);
        BoundingBox3D bbox = sphere.GetBoundingBox();
        
        REQUIRE_THAT(bbox.Min().X() , RealApprox(-5.0));
        REQUIRE_THAT(bbox.Min().Y() , RealApprox(-5.0));
        REQUIRE_THAT(bbox.Min().Z() , RealApprox(-5.0));
        
        REQUIRE_THAT(bbox.Max().X() , RealApprox(5.0));
        REQUIRE_THAT(bbox.Max().Y() , RealApprox(5.0));
        REQUIRE_THAT(bbox.Max().Z() , RealApprox(5.0));
    }

    SECTION("Sphere centered at arbitrary point")
    {
        Pnt3Cart center(10.0, 20.0, 30.0);
        Sphere3D sphere(3.0, center);
        BoundingBox3D bbox = sphere.GetBoundingBox();
        
        REQUIRE_THAT(bbox.Min().X() , RealApprox(7.0));
        REQUIRE_THAT(bbox.Min().Y() , RealApprox(17.0));
        REQUIRE_THAT(bbox.Min().Z() , RealApprox(27.0));
        
        REQUIRE_THAT(bbox.Max().X() , RealApprox(13.0));
        REQUIRE_THAT(bbox.Max().Y() , RealApprox(23.0));
        REQUIRE_THAT(bbox.Max().Z() , RealApprox(33.0));
    }

    SECTION("Bounding box dimensions match diameter")
    {
        Sphere3D sphere(7.5);
        BoundingBox3D bbox = sphere.GetBoundingBox();
        Real width = bbox.Max().X() - bbox.Min().X();
        Real height = bbox.Max().Y() - bbox.Min().Y();
        Real depth = bbox.Max().Z() - bbox.Min().Z();
        
        REQUIRE_THAT(width , RealApprox(15.0));
        REQUIRE_THAT(height , RealApprox(15.0));
        REQUIRE_THAT(depth , RealApprox(15.0));
    }
}

TEST_CASE("Sphere3D::GetBoundingSphere", "[geometry][sphere][bounding]")
{
    SECTION("Bounding sphere equals sphere itself - origin")
    {
        Sphere3D sphere(5.0);
        BoundingSphere3D bsphere = sphere.GetBoundingSphere();
        
        REQUIRE_THAT(bsphere.Center().X() , RealApprox(0.0));
        REQUIRE_THAT(bsphere.Center().Y() , RealApprox(0.0));
        REQUIRE_THAT(bsphere.Center().Z() , RealApprox(0.0));
        REQUIRE_THAT(bsphere.Radius() , RealApprox(5.0));
    }

    SECTION("Bounding sphere equals sphere itself - arbitrary center")
    {
        Pnt3Cart center(10.0, 20.0, 30.0);
        Sphere3D sphere(3.0, center);
        BoundingSphere3D bsphere = sphere.GetBoundingSphere();
        
        REQUIRE_THAT(bsphere.Center().X() , RealApprox(10.0));
        REQUIRE_THAT(bsphere.Center().Y() , RealApprox(20.0));
        REQUIRE_THAT(bsphere.Center().Z() , RealApprox(30.0));
        REQUIRE_THAT(bsphere.Radius() , RealApprox(3.0));
    }

    SECTION("Bounding sphere volume matches sphere volume")
    {
        Sphere3D sphere(7.0);
        BoundingSphere3D bsphere = sphere.GetBoundingSphere();
        
        REQUIRE_THAT(bsphere.Volume() , RealApprox(sphere.Volume()).epsilon(1e-10));
    }
}

TEST_CASE("Sphere3D::ToString", "[geometry][sphere][string]")
{
    SECTION("Sphere at origin")
    {
        Sphere3D sphere(5.0);
        std::string str = sphere.ToString();
        
        // Verify string contains key components
        REQUIRE(str.find("Sphere3D") != std::string::npos);
        REQUIRE(str.find("Center") != std::string::npos);
        REQUIRE(str.find("Radius") != std::string::npos);
        REQUIRE(str.find("Volume") != std::string::npos);
        REQUIRE(str.find("SurfaceArea") != std::string::npos);
    }

    SECTION("Sphere at arbitrary center")
    {
        Pnt3Cart center(10.0, 20.0, 30.0);
        Sphere3D sphere(3.0, center);
        std::string str = sphere.ToString();
        
        // Verify string contains coordinates and radius
        REQUIRE(str.find("10") != std::string::npos);
        REQUIRE(str.find("20") != std::string::npos);
        REQUIRE(str.find("30") != std::string::npos);
        REQUIRE(str.find("3") != std::string::npos);
    }
}

TEST_CASE("Sphere3D::Getters", "[geometry][sphere][getters]")
{
    SECTION("GetRadius returns correct value")
    {
        Sphere3D sphere(7.5);
        REQUIRE_THAT(sphere.GetRadius() , RealApprox(7.5));
    }

    SECTION("GetNumLatitude and GetNumLongitude return tessellation parameters")
    {
        Sphere3D sphere(5.0, 32, 48);
        REQUIRE(sphere.GetNumLatitude() == 32);
        REQUIRE(sphere.GetNumLongitude() == 48);
    }

    SECTION("Default tessellation parameters")
    {
        Sphere3D sphere(5.0);
        REQUIRE(sphere.GetNumLatitude() == 16);
        REQUIRE(sphere.GetNumLongitude() == 20);
    }
}

TEST_CASE("Sphere3D::IsInside", "[geometry][sphere][inside]")
{
    SECTION("Point at center is inside")
    {
        Sphere3D sphere(5.0);
        REQUIRE(sphere.IsInside(Pnt3Cart(0.0, 0.0, 0.0)) == true);
    }

    SECTION("Point on surface is inside (boundary)")
    {
        Sphere3D sphere(5.0);
        REQUIRE(sphere.IsInside(Pnt3Cart(5.0, 0.0, 0.0)) == true);
        REQUIRE(sphere.IsInside(Pnt3Cart(0.0, 5.0, 0.0)) == true);
        REQUIRE(sphere.IsInside(Pnt3Cart(0.0, 0.0, 5.0)) == true);
    }

    SECTION("Point just outside is not inside")
    {
        Sphere3D sphere(5.0);
        REQUIRE(sphere.IsInside(Pnt3Cart(5.1, 0.0, 0.0)) == false);
        REQUIRE(sphere.IsInside(Pnt3Cart(0.0, 5.1, 0.0)) == false);
        REQUIRE(sphere.IsInside(Pnt3Cart(0.0, 0.0, 5.1)) == false);
    }

    SECTION("Point far outside is not inside")
    {
        Sphere3D sphere(5.0);
        REQUIRE(sphere.IsInside(Pnt3Cart(100.0, 0.0, 0.0)) == false);
        REQUIRE(sphere.IsInside(Pnt3Cart(0.0, -100.0, 0.0)) == false);
        REQUIRE(sphere.IsInside(Pnt3Cart(0.0, 0.0, 100.0)) == false);
    }

    SECTION("Points with sphere centered at arbitrary location")
    {
        Pnt3Cart center(10.0, 20.0, 30.0);
        Sphere3D sphere(3.0, center);
        
        REQUIRE(sphere.IsInside(center) == true);
        REQUIRE(sphere.IsInside(Pnt3Cart(13.0, 20.0, 30.0)) == true);  // On surface
        REQUIRE(sphere.IsInside(Pnt3Cart(10.0, 23.0, 30.0)) == true);  // On surface
        REQUIRE(sphere.IsInside(Pnt3Cart(10.0, 20.0, 33.0)) == true);  // On surface
        
        REQUIRE(sphere.IsInside(Pnt3Cart(14.0, 20.0, 30.0)) == false);  // Outside
        REQUIRE(sphere.IsInside(Pnt3Cart(10.0, 24.0, 30.0)) == false);  // Outside
    }

    SECTION("3-4-5 triangle test - point at distance sqrt(50)")
    {
        Sphere3D sphere(7.0);  // radius < sqrt(50) ≈ 7.071
        REQUIRE(sphere.IsInside(Pnt3Cart(3.0, 4.0, 5.0)) == false);  // distance = sqrt(50) ≈ 7.071 > 7.0
        
        Sphere3D sphere2(8.0);  // radius > sqrt(50)
        REQUIRE(sphere2.IsInside(Pnt3Cart(3.0, 4.0, 5.0)) == true);
        
        Sphere3D sphere3(6.0);  // radius < sqrt(50)
        REQUIRE(sphere3.IsInside(Pnt3Cart(3.0, 4.0, 5.0)) == false);
    }
}

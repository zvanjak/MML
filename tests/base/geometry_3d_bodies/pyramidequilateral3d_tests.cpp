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

TEST_CASE("PyramidEquilateral3D::Constructor", "[geometry][pyramidequilateral][constructor]")
{
    SECTION("Equilateral pyramid with base side 6")
    {
        PyramidEquilateral3D pyramid(6.0);
        
        // Height should be a/sqrt(3) = 6/sqrt(3) = 2*sqrt(3)
        Real expected_height = 6.0 / std::sqrt(3.0);
        REQUIRE_THAT(pyramid.GetHeight() , RealApprox(expected_height));
        REQUIRE_THAT(pyramid.GetBaseSize() , RealApprox(6.0));
    }

    SECTION("Equilateral pyramid at arbitrary position")
    {
        Vec3Cart baseCenter(10.0, 20.0, 30.0);
        PyramidEquilateral3D pyramid(9.0, baseCenter);
        
        Real expected_height = 9.0 / std::sqrt(3.0);
        REQUIRE_THAT(pyramid.GetHeight() , RealApprox(expected_height));
        REQUIRE_THAT(pyramid.GetBaseSize() , RealApprox(9.0));
        
        Vec3Cart center = pyramid.GetBaseCenter();
        REQUIRE_THAT(center.X() , RealApprox(10.0));
        REQUIRE_THAT(center.Y() , RealApprox(20.0));
        REQUIRE_THAT(center.Z() , RealApprox(30.0));
    }
}

TEST_CASE("PyramidEquilateral3D::Volume", "[geometry][pyramidequilateral][volume]")
{
    SECTION("Equilateral pyramid inherits volume calculation")
    {
        PyramidEquilateral3D pyramid(6.0);
        // h = 6/sqrt(3) = 2*sqrt(3)
        // Volume = (1/3) * 6² * 2*sqrt(3) = (1/3) * 36 * 2*sqrt(3) = 24*sqrt(3)
        Real expected = 24.0 * std::sqrt(3.0);
        REQUIRE_THAT(pyramid.Volume() , RealApprox(expected));
    }

    SECTION("Equilateral pyramid with base 3")
    {
        PyramidEquilateral3D pyramid(3.0);
        // h = 3/sqrt(3) = sqrt(3)
        // Volume = (1/3) * 9 * sqrt(3) = 3*sqrt(3)
        Real expected = 3.0 * std::sqrt(3.0);
        REQUIRE_THAT(pyramid.Volume() , RealApprox(expected));
    }
}

TEST_CASE("PyramidEquilateral3D::SurfaceArea", "[geometry][pyramidequilateral][surface]")
{
    SECTION("Surface area calculation")
    {
        PyramidEquilateral3D pyramid(6.0);
        // h = 2*sqrt(3)
        // Slant height s = sqrt(h² + (a/2)²) = sqrt(12 + 9) = sqrt(21)
        // SA = a² + 2*a*s = 36 + 12*sqrt(21)
        Real h = 6.0 / std::sqrt(3.0);
        Real slant = std::sqrt(h * h + 3.0 * 3.0);
        Real expected = 36.0 + 12.0 * slant;
        REQUIRE_THAT(pyramid.SurfaceArea() , RealApprox(expected));
    }
}

TEST_CASE("PyramidEquilateral3D::GetCenter", "[geometry][pyramidequilateral][center]")
{
    SECTION("Center at h/4 from base")
    {
        PyramidEquilateral3D pyramid(12.0);
        Pnt3Cart center = pyramid.GetCenter();
        
        // h = 12/sqrt(3) = 4*sqrt(3)
        // Centroid at h/4 = sqrt(3)
        Real expected_z = 12.0 / (4.0 * std::sqrt(3.0));
        REQUIRE_THAT(center.X() , RealApprox(0.0));
        REQUIRE_THAT(center.Y() , RealApprox(0.0));
        REQUIRE_THAT(center.Z() , RealApprox(expected_z));
    }

    SECTION("Center at arbitrary position")
    {
        Vec3Cart baseCenter(5.0, 10.0, 15.0);
        PyramidEquilateral3D pyramid(6.0, baseCenter);
        Pnt3Cart center = pyramid.GetCenter();
        
        // h = 6/sqrt(3) = 2*sqrt(3), centroid at h/4 = sqrt(3)/2
        Real offset_z = 6.0 / (4.0 * std::sqrt(3.0));
        REQUIRE_THAT(center.X() , RealApprox(5.0));
        REQUIRE_THAT(center.Y() , RealApprox(10.0));
        REQUIRE_THAT(center.Z() , RealApprox(15.0 + offset_z));
    }
}

TEST_CASE("PyramidEquilateral3D::GetBoundingBox", "[geometry][pyramidequilateral][bounding]")
{
    SECTION("Bounding box from base to apex")
    {
        PyramidEquilateral3D pyramid(8.0);
        BoundingBox3D bbox = pyramid.GetBoundingBox();
        
        Real h = 8.0 / std::sqrt(3.0);
        REQUIRE_THAT(bbox.Min().X() , RealApprox(-4.0));
        REQUIRE_THAT(bbox.Min().Y() , RealApprox(-4.0));
        REQUIRE_THAT(bbox.Min().Z() , RealApprox(0.0));
        
        REQUIRE_THAT(bbox.Max().X() , RealApprox(4.0));
        REQUIRE_THAT(bbox.Max().Y() , RealApprox(4.0));
        REQUIRE_THAT(bbox.Max().Z() , RealApprox(h));
    }
}

TEST_CASE("PyramidEquilateral3D::GetBoundingSphere", "[geometry][pyramidequilateral][bounding]")
{
    SECTION("Bounding sphere contains all vertices")
    {
        PyramidEquilateral3D pyramid(6.0);
        BoundingSphere3D bsphere = pyramid.GetBoundingSphere();
        
        Real h = 6.0 / std::sqrt(3.0);
        
        // Check base corners
        REQUIRE(bsphere.Contains(Pnt3Cart(3.0, 3.0, 0.0)));
        REQUIRE(bsphere.Contains(Pnt3Cart(-3.0, -3.0, 0.0)));
        
        // Check apex
        REQUIRE(bsphere.Contains(Pnt3Cart(0.0, 0.0, h)));
    }

    SECTION("Bounding sphere centered at centroid")
    {
        PyramidEquilateral3D pyramid(12.0);
        BoundingSphere3D bsphere = pyramid.GetBoundingSphere();
        
        Real h = 12.0 / std::sqrt(3.0);
        Real centroid_z = h / 4.0;
        
        REQUIRE_THAT(bsphere.Center().X() , RealApprox(0.0));
        REQUIRE_THAT(bsphere.Center().Y() , RealApprox(0.0));
        REQUIRE_THAT(bsphere.Center().Z() , RealApprox(centroid_z));
    }
}

TEST_CASE("PyramidEquilateral3D::ToString", "[geometry][pyramidequilateral][string]")
{
    SECTION("String contains pyramid info")
    {
        PyramidEquilateral3D pyramid(8.0);
        std::string str = pyramid.ToString();
        
        REQUIRE(str.find("PyramidEquilateral3D") != std::string::npos);
        REQUIRE(str.find("Equilateral") != std::string::npos);
        REQUIRE(str.find("BaseSize") != std::string::npos);
        REQUIRE(str.find("Height") != std::string::npos);
        REQUIRE(str.find("Volume") != std::string::npos);
        REQUIRE(str.find("SurfaceArea") != std::string::npos);
    }

    SECTION("ToString at arbitrary position")
    {
        Vec3Cart baseCenter(10.0, 20.0, 30.0);
        PyramidEquilateral3D pyramid(6.0, baseCenter);
        std::string str = pyramid.ToString();
        
        REQUIRE(str.find("6") != std::string::npos);  // BaseSize
    }
}

TEST_CASE("PyramidEquilateral3D::IsInside", "[geometry][pyramidequilateral][inside]")
{
    SECTION("Base center is inside")
    {
        PyramidEquilateral3D pyramid(10.0);
        REQUIRE(pyramid.IsInside(Pnt3Cart(0.0, 0.0, 0.0)) == true);
    }

    SECTION("Centroid is inside")
    {
        PyramidEquilateral3D pyramid(12.0);
        Real h = 12.0 / std::sqrt(3.0);
        Real centroid_z = h / 4.0;
        REQUIRE(pyramid.IsInside(Pnt3Cart(0.0, 0.0, centroid_z)) == true);
    }

    SECTION("Apex is inside")
    {
        PyramidEquilateral3D pyramid(9.0);
        Real h = 9.0 / std::sqrt(3.0);
        REQUIRE(pyramid.IsInside(Pnt3Cart(0.0, 0.0, h)) == true);
    }

    SECTION("Base corners are inside")
    {
        PyramidEquilateral3D pyramid(8.0);
        REQUIRE(pyramid.IsInside(Pnt3Cart(4.0, 4.0, 0.0)) == true);
        REQUIRE(pyramid.IsInside(Pnt3Cart(-4.0, -4.0, 0.0)) == true);
    }

    SECTION("Points outside are not inside")
    {
        PyramidEquilateral3D pyramid(6.0);
        Real h = 6.0 / std::sqrt(3.0);
        
        REQUIRE(pyramid.IsInside(Pnt3Cart(0.0, 0.0, -0.1)) == false);  // Below base
        REQUIRE(pyramid.IsInside(Pnt3Cart(0.0, 0.0, h + 0.1)) == false);  // Above apex
        REQUIRE(pyramid.IsInside(Pnt3Cart(3.1, 0.0, 0.0)) == false);  // Outside base
    }

    SECTION("Equilateral pyramid at arbitrary position")
    {
        Vec3Cart baseCenter(10.0, 20.0, 30.0);
        PyramidEquilateral3D pyramid(6.0, baseCenter);
        Real h = 6.0 / std::sqrt(3.0);
        
        // Base center
        REQUIRE(pyramid.IsInside(Pnt3Cart(10.0, 20.0, 30.0)) == true);
        
        // Near apex (slightly below to avoid boundary condition)
        REQUIRE(pyramid.IsInside(Pnt3Cart(10.0, 20.0, 30.0 + h - 0.01)) == true);
        
        // Outside
        REQUIRE(pyramid.IsInside(Pnt3Cart(10.0, 20.0, 29.9)) == false);
    }
}

TEST_CASE("PyramidEquilateral3D::InheritedMethods", "[geometry][pyramidequilateral][inheritance]")
{
    SECTION("Inherits all IBody methods from Pyramid3D")
    {
        PyramidEquilateral3D pyramid(6.0);
        
        // These should all work without error
        Real volume = pyramid.Volume();
        Real surfaceArea = pyramid.SurfaceArea();
        Pnt3Cart center = pyramid.GetCenter();
        BoundingBox3D bbox = pyramid.GetBoundingBox();
        BoundingSphere3D bsphere = pyramid.GetBoundingSphere();
        std::string str = pyramid.ToString();
        
        REQUIRE(volume > 0);
        REQUIRE(surfaceArea > 0);
        REQUIRE(str.length() > 0);
    }

    SECTION("Getters work correctly")
    {
        PyramidEquilateral3D pyramid(9.0);
        
        REQUIRE_THAT(pyramid.GetBaseSize() , RealApprox(9.0));
        REQUIRE_THAT(pyramid.GetHeight() , RealApprox(9.0 / std::sqrt(3.0)).epsilon(1e-10));
        
        Vec3Cart baseCenter = pyramid.GetBaseCenter();
        REQUIRE_THAT(baseCenter.X() , RealApprox(0.0));
        REQUIRE_THAT(baseCenter.Y() , RealApprox(0.0));
        REQUIRE_THAT(baseCenter.Z() , RealApprox(0.0));
    }
}

TEST_CASE("PyramidEquilateral3D::GeometricProperties", "[geometry][pyramidequilateral][properties]")
{
    SECTION("Height-to-base ratio is 1/sqrt(3)")
    {
        for (Real a : {3.0, 6.0, 9.0, 12.0}) {
            PyramidEquilateral3D pyramid(a);
            Real expected_ratio = 1.0 / std::sqrt(3.0);
            Real actual_ratio = pyramid.GetHeight() / pyramid.GetBaseSize();
            REQUIRE_THAT(actual_ratio , RealApprox(expected_ratio));
        }
    }

    SECTION("Volume scales as a³")
    {
        PyramidEquilateral3D pyramid1(2.0);
        PyramidEquilateral3D pyramid2(4.0);
        
        // Volume should scale as a³, so V2/V1 should be 8
        Real ratio = pyramid2.Volume() / pyramid1.Volume();
        REQUIRE_THAT(ratio , RealApprox(8.0));
    }

    SECTION("Surface area scales as a²")
    {
        PyramidEquilateral3D pyramid1(3.0);
        PyramidEquilateral3D pyramid2(6.0);
        
        // Surface area should scale as a², so SA2/SA1 should be 4
        Real ratio = pyramid2.SurfaceArea() / pyramid1.SurfaceArea();
        REQUIRE_THAT(ratio , RealApprox(4.0));
    }
}

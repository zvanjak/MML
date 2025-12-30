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

TEST_CASE("CubeWithTriangles3D::Volume", "[geometry][cubewithtriangles][volume]")
{
    SECTION("Unit cube")
    {
        CubeWithTriangles3D cube(1.0);
        REQUIRE_THAT(cube.Volume() , RealApprox(1.0));
    }

    SECTION("Cube with side 5")
    {
        CubeWithTriangles3D cube(5.0);
        Real expected = 125.0;  // 5³
        REQUIRE_THAT(cube.Volume() , RealApprox(expected));
    }

    SECTION("Cube at arbitrary position")
    {
        Pnt3Cart center(10.0, 20.0, 30.0);
        CubeWithTriangles3D cube(3.0, center);
        Real expected = 27.0;  // 3³
        REQUIRE_THAT(cube.Volume() , RealApprox(expected));
    }
}

TEST_CASE("CubeWithTriangles3D::SurfaceArea", "[geometry][cubewithtriangles][surface]")
{
    SECTION("Unit cube")
    {
        CubeWithTriangles3D cube(1.0);
        Real expected = 6.0;  // 6 * 1²
        REQUIRE_THAT(cube.SurfaceArea() , RealApprox(expected));
    }

    SECTION("Cube with side 4")
    {
        CubeWithTriangles3D cube(4.0);
        Real expected = 96.0;  // 6 * 4²
        REQUIRE_THAT(cube.SurfaceArea() , RealApprox(expected));
    }

    SECTION("Cube at arbitrary position")
    {
        Pnt3Cart center(5.0, 10.0, 15.0);
        CubeWithTriangles3D cube(6.0, center);
        Real expected = 216.0;  // 6 * 6²
        REQUIRE_THAT(cube.SurfaceArea() , RealApprox(expected));
    }
}

TEST_CASE("CubeWithTriangles3D::GetCenter", "[geometry][cubewithtriangles][center]")
{
    SECTION("Cube at origin")
    {
        CubeWithTriangles3D cube(5.0);
        Pnt3Cart center = cube.GetCenter();
        REQUIRE_THAT(center.X() , RealApprox(0.0));
        REQUIRE_THAT(center.Y() , RealApprox(0.0));
        REQUIRE_THAT(center.Z() , RealApprox(0.0));
    }

    SECTION("Cube at arbitrary position")
    {
        Pnt3Cart expectedCenter(10.0, 20.0, 30.0);
        CubeWithTriangles3D cube(7.0, expectedCenter);
        Pnt3Cart center = cube.GetCenter();
        REQUIRE_THAT(center.X() , RealApprox(10.0));
        REQUIRE_THAT(center.Y() , RealApprox(20.0));
        REQUIRE_THAT(center.Z() , RealApprox(30.0));
    }
}

TEST_CASE("CubeWithTriangles3D::GetBoundingBox", "[geometry][cubewithtriangles][bounding]")
{
    SECTION("Cube at origin")
    {
        CubeWithTriangles3D cube(10.0);
        BoundingBox3D bbox = cube.GetBoundingBox();
        
        REQUIRE_THAT(bbox.Min().X() , RealApprox(-5.0));
        REQUIRE_THAT(bbox.Min().Y() , RealApprox(-5.0));
        REQUIRE_THAT(bbox.Min().Z() , RealApprox(-5.0));
        
        REQUIRE_THAT(bbox.Max().X() , RealApprox(5.0));
        REQUIRE_THAT(bbox.Max().Y() , RealApprox(5.0));
        REQUIRE_THAT(bbox.Max().Z() , RealApprox(5.0));
    }

    SECTION("Cube at arbitrary position")
    {
        Pnt3Cart center(10.0, 20.0, 30.0);
        CubeWithTriangles3D cube(6.0, center);
        BoundingBox3D bbox = cube.GetBoundingBox();
        
        REQUIRE_THAT(bbox.Min().X() , RealApprox(7.0));
        REQUIRE_THAT(bbox.Min().Y() , RealApprox(17.0));
        REQUIRE_THAT(bbox.Min().Z() , RealApprox(27.0));
        
        REQUIRE_THAT(bbox.Max().X() , RealApprox(13.0));
        REQUIRE_THAT(bbox.Max().Y() , RealApprox(23.0));
        REQUIRE_THAT(bbox.Max().Z() , RealApprox(33.0));
    }
}

TEST_CASE("CubeWithTriangles3D::GetBoundingSphere", "[geometry][cubewithtriangles][bounding]")
{
    SECTION("Sphere centered at cube center")
    {
        CubeWithTriangles3D cube(6.0);
        BoundingSphere3D bsphere = cube.GetBoundingSphere();
        
        REQUIRE_THAT(bsphere.Center().X() , RealApprox(0.0));
        REQUIRE_THAT(bsphere.Center().Y() , RealApprox(0.0));
        REQUIRE_THAT(bsphere.Center().Z() , RealApprox(0.0));
        
        // Radius = (a/2)*sqrt(3) = 3*sqrt(3)
        Real expected_radius = 3.0 * std::sqrt(3.0);
        REQUIRE_THAT(bsphere.Radius() , RealApprox(expected_radius));
    }

    SECTION("Bounding sphere contains all corners")
    {
        CubeWithTriangles3D cube(8.0);
        BoundingSphere3D bsphere = cube.GetBoundingSphere();
        
        // Check all 8 corners
        REQUIRE(bsphere.Contains(Pnt3Cart(4.0, 4.0, 4.0)));
        REQUIRE(bsphere.Contains(Pnt3Cart(4.0, 4.0, -4.0)));
        REQUIRE(bsphere.Contains(Pnt3Cart(4.0, -4.0, 4.0)));
        REQUIRE(bsphere.Contains(Pnt3Cart(4.0, -4.0, -4.0)));
        REQUIRE(bsphere.Contains(Pnt3Cart(-4.0, 4.0, 4.0)));
        REQUIRE(bsphere.Contains(Pnt3Cart(-4.0, 4.0, -4.0)));
        REQUIRE(bsphere.Contains(Pnt3Cart(-4.0, -4.0, 4.0)));
        REQUIRE(bsphere.Contains(Pnt3Cart(-4.0, -4.0, -4.0)));
    }

    SECTION("Bounding sphere at arbitrary position")
    {
        Pnt3Cart center(10.0, 20.0, 30.0);
        CubeWithTriangles3D cube(4.0, center);
        BoundingSphere3D bsphere = cube.GetBoundingSphere();
        
        REQUIRE_THAT(bsphere.Center().X() , RealApprox(10.0));
        REQUIRE_THAT(bsphere.Center().Y() , RealApprox(20.0));
        REQUIRE_THAT(bsphere.Center().Z() , RealApprox(30.0));
        
        Real expected_radius = 2.0 * std::sqrt(3.0);
        REQUIRE_THAT(bsphere.Radius() , RealApprox(expected_radius));
    }
}

TEST_CASE("CubeWithTriangles3D::ToString", "[geometry][cubewithtriangles][string]")
{
    SECTION("Cube at origin")
    {
        CubeWithTriangles3D cube(5.0);
        std::string str = cube.ToString();
        
        REQUIRE(str.find("CubeWithTriangles3D") != std::string::npos);
        REQUIRE(str.find("Center") != std::string::npos);
        REQUIRE(str.find("Side") != std::string::npos);
        REQUIRE(str.find("Volume") != std::string::npos);
        REQUIRE(str.find("SurfaceArea") != std::string::npos);
        REQUIRE(str.find("Triangles") != std::string::npos);
    }

    SECTION("Cube at arbitrary position")
    {
        Pnt3Cart center(10.0, 20.0, 30.0);
        CubeWithTriangles3D cube(8.0, center);
        std::string str = cube.ToString();
        
        REQUIRE(str.find("8") != std::string::npos);  // Side
    }
}

TEST_CASE("CubeWithTriangles3D::Getters", "[geometry][cubewithtriangles][getters]")
{
    SECTION("GetSide")
    {
        CubeWithTriangles3D cube(7.5);
        REQUIRE_THAT(cube.GetSide() , RealApprox(7.5));
    }

    SECTION("GetSide at arbitrary position")
    {
        Pnt3Cart center(5.0, 10.0, 15.0);
        CubeWithTriangles3D cube(12.0, center);
        REQUIRE_THAT(cube.GetSide() , RealApprox(12.0));
    }
}

TEST_CASE("CubeWithTriangles3D::IsInside", "[geometry][cubewithtriangles][inside]")
{
    SECTION("Point at center is inside")
    {
        CubeWithTriangles3D cube(10.0);
        REQUIRE(cube.IsInside(Pnt3Cart(0.0, 0.0, 0.0)) == true);
    }

    SECTION("Points on faces are inside")
    {
        CubeWithTriangles3D cube(10.0);
        REQUIRE(cube.IsInside(Pnt3Cart(5.0, 0.0, 0.0)) == true);
        REQUIRE(cube.IsInside(Pnt3Cart(-5.0, 0.0, 0.0)) == true);
        REQUIRE(cube.IsInside(Pnt3Cart(0.0, 5.0, 0.0)) == true);
        REQUIRE(cube.IsInside(Pnt3Cart(0.0, -5.0, 0.0)) == true);
        REQUIRE(cube.IsInside(Pnt3Cart(0.0, 0.0, 5.0)) == true);
        REQUIRE(cube.IsInside(Pnt3Cart(0.0, 0.0, -5.0)) == true);
    }

    SECTION("Points on corners are inside")
    {
        CubeWithTriangles3D cube(10.0);
        REQUIRE(cube.IsInside(Pnt3Cart(5.0, 5.0, 5.0)) == true);
        REQUIRE(cube.IsInside(Pnt3Cart(-5.0, -5.0, -5.0)) == true);
        REQUIRE(cube.IsInside(Pnt3Cart(5.0, -5.0, 5.0)) == true);
        REQUIRE(cube.IsInside(Pnt3Cart(-5.0, 5.0, -5.0)) == true);
    }

    SECTION("Points outside are not inside")
    {
        CubeWithTriangles3D cube(10.0);
        REQUIRE(cube.IsInside(Pnt3Cart(5.1, 0.0, 0.0)) == false);
        REQUIRE(cube.IsInside(Pnt3Cart(0.0, 5.1, 0.0)) == false);
        REQUIRE(cube.IsInside(Pnt3Cart(0.0, 0.0, 5.1)) == false);
        REQUIRE(cube.IsInside(Pnt3Cart(6.0, 6.0, 6.0)) == false);
    }

    SECTION("Cube at arbitrary position")
    {
        Pnt3Cart center(10.0, 20.0, 30.0);
        CubeWithTriangles3D cube(8.0, center);
        
        // Center is inside
        REQUIRE(cube.IsInside(Pnt3Cart(10.0, 20.0, 30.0)) == true);
        
        // Face centers are inside
        REQUIRE(cube.IsInside(Pnt3Cart(14.0, 20.0, 30.0)) == true);
        REQUIRE(cube.IsInside(Pnt3Cart(6.0, 20.0, 30.0)) == true);
        
        // Corners are inside
        REQUIRE(cube.IsInside(Pnt3Cart(14.0, 24.0, 34.0)) == true);
        REQUIRE(cube.IsInside(Pnt3Cart(6.0, 16.0, 26.0)) == true);
        
        // Outside
        REQUIRE(cube.IsInside(Pnt3Cart(14.1, 20.0, 30.0)) == false);
        REQUIRE(cube.IsInside(Pnt3Cart(10.0, 24.1, 30.0)) == false);
    }
}

TEST_CASE("CubeWithTriangles3D::TriangleSurfaces", "[geometry][cubewithtriangles][triangles]")
{
    SECTION("Cube has 12 triangles (2 per face)")
    {
        CubeWithTriangles3D cube(5.0);
        // Access through BodyWithTriangleSurfaces interface
        // A cube has 6 faces, each decomposed into 2 triangles = 12 triangles
        REQUIRE(cube.ToString().find("Triangles=12") != std::string::npos);
    }

    SECTION("Cube at arbitrary position has 12 triangles")
    {
        Pnt3Cart center(10.0, 20.0, 30.0);
        CubeWithTriangles3D cube(8.0, center);
        REQUIRE(cube.ToString().find("Triangles=12") != std::string::npos);
    }
}

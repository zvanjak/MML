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

TEST_CASE("Cube3D::Volume", "[geometry][cube][volume]")
{
    SECTION("Unit cube centered at origin")
    {
        Cube3D cube(1.0);
        REQUIRE_THAT(cube.Volume() , RealApprox(1.0));
    }

    SECTION("Cube with side 5 centered at origin")
    {
        Cube3D cube(5.0);
        REQUIRE_THAT(cube.Volume() , RealApprox(125.0));  // 5³ = 125
    }

    SECTION("Cube centered at arbitrary point")
    {
        Pnt3Cart center(10.0, 20.0, 30.0);
        Cube3D cube(3.0, center);
        REQUIRE_THAT(cube.Volume() , RealApprox(27.0));  // 3³ = 27
    }

    SECTION("Small cube")
    {
        Cube3D cube(0.1);
        REQUIRE_THAT(cube.Volume() , RealApprox(0.001));  // 0.1³ = 0.001
    }
}

TEST_CASE("Cube3D::SurfaceArea", "[geometry][cube][surface]")
{
    SECTION("Unit cube centered at origin")
    {
        Cube3D cube(1.0);
        REQUIRE_THAT(cube.SurfaceArea() , RealApprox(6.0));  // 6 * 1² = 6
    }

    SECTION("Cube with side 5 centered at origin")
    {
        Cube3D cube(5.0);
        REQUIRE_THAT(cube.SurfaceArea() , RealApprox(150.0));  // 6 * 25 = 150
    }

    SECTION("Cube centered at arbitrary point")
    {
        Pnt3Cart center(10.0, 20.0, 30.0);
        Cube3D cube(3.0, center);
        REQUIRE_THAT(cube.SurfaceArea() , RealApprox(54.0));  // 6 * 9 = 54
    }

    SECTION("Large cube")
    {
        Cube3D cube(10.0);
        REQUIRE_THAT(cube.SurfaceArea() , RealApprox(600.0));  // 6 * 100 = 600
    }
}

TEST_CASE("Cube3D::GetCenter", "[geometry][cube][center]")
{
    SECTION("Cube centered at origin")
    {
        Cube3D cube(5.0);
        Pnt3Cart center = cube.GetCenter();
        REQUIRE_THAT(center.X() , RealApprox(0.0));
        REQUIRE_THAT(center.Y() , RealApprox(0.0));
        REQUIRE_THAT(center.Z() , RealApprox(0.0));
    }

    SECTION("Cube centered at arbitrary point")
    {
        Pnt3Cart expected(10.0, 20.0, 30.0);
        Cube3D cube(5.0, expected);
        Pnt3Cart center = cube.GetCenter();
        REQUIRE_THAT(center.X() , RealApprox(10.0));
        REQUIRE_THAT(center.Y() , RealApprox(20.0));
        REQUIRE_THAT(center.Z() , RealApprox(30.0));
    }
}

TEST_CASE("Cube3D::GetBoundingBox", "[geometry][cube][bounding]")
{
    SECTION("Unit cube centered at origin")
    {
        Cube3D cube(1.0);
        BoundingBox3D bbox = cube.GetBoundingBox();
        
        REQUIRE_THAT(bbox.Min().X() , RealApprox(-0.5));
        REQUIRE_THAT(bbox.Min().Y() , RealApprox(-0.5));
        REQUIRE_THAT(bbox.Min().Z() , RealApprox(-0.5));
        
        REQUIRE_THAT(bbox.Max().X() , RealApprox(0.5));
        REQUIRE_THAT(bbox.Max().Y() , RealApprox(0.5));
        REQUIRE_THAT(bbox.Max().Z() , RealApprox(0.5));
    }

    SECTION("Cube with side 10 centered at origin")
    {
        Cube3D cube(10.0);
        BoundingBox3D bbox = cube.GetBoundingBox();
        
        REQUIRE_THAT(bbox.Min().X() , RealApprox(-5.0));
        REQUIRE_THAT(bbox.Min().Y() , RealApprox(-5.0));
        REQUIRE_THAT(bbox.Min().Z() , RealApprox(-5.0));
        
        REQUIRE_THAT(bbox.Max().X() , RealApprox(5.0));
        REQUIRE_THAT(bbox.Max().Y() , RealApprox(5.0));
        REQUIRE_THAT(bbox.Max().Z() , RealApprox(5.0));
    }

    SECTION("Cube centered at arbitrary point")
    {
        Pnt3Cart center(10.0, 20.0, 30.0);
        Cube3D cube(6.0, center);
        BoundingBox3D bbox = cube.GetBoundingBox();
        
        REQUIRE_THAT(bbox.Min().X() , RealApprox(7.0));
        REQUIRE_THAT(bbox.Min().Y() , RealApprox(17.0));
        REQUIRE_THAT(bbox.Min().Z() , RealApprox(27.0));
        
        REQUIRE_THAT(bbox.Max().X() , RealApprox(13.0));
        REQUIRE_THAT(bbox.Max().Y() , RealApprox(23.0));
        REQUIRE_THAT(bbox.Max().Z() , RealApprox(33.0));
    }

    SECTION("Bounding box dimensions match cube side")
    {
        Cube3D cube(8.0);
        BoundingBox3D bbox = cube.GetBoundingBox();
        Real width = bbox.Max().X() - bbox.Min().X();
        Real height = bbox.Max().Y() - bbox.Min().Y();
        Real depth = bbox.Max().Z() - bbox.Min().Z();
        
        REQUIRE_THAT(width , RealApprox(8.0));
        REQUIRE_THAT(height , RealApprox(8.0));
        REQUIRE_THAT(depth , RealApprox(8.0));
    }
}

TEST_CASE("Cube3D::GetBoundingSphere", "[geometry][cube][bounding]")
{
    SECTION("Bounding sphere centered at origin")
    {
        Cube3D cube(2.0);
        BoundingSphere3D bsphere = cube.GetBoundingSphere();
        
        REQUIRE_THAT(bsphere.Center().X() , RealApprox(0.0));
        REQUIRE_THAT(bsphere.Center().Y() , RealApprox(0.0));
        REQUIRE_THAT(bsphere.Center().Z() , RealApprox(0.0));
        
        // Radius = (a√3)/2 = 2√3/2 = √3 ≈ 1.732
        Real expected_radius = 2.0 * std::sqrt(3.0) / 2.0;
        REQUIRE_THAT(bsphere.Radius() , RealApprox(expected_radius));
    }

    SECTION("Bounding sphere at arbitrary center")
    {
        Pnt3Cart center(10.0, 20.0, 30.0);
        Cube3D cube(4.0, center);
        BoundingSphere3D bsphere = cube.GetBoundingSphere();
        
        REQUIRE_THAT(bsphere.Center().X() , RealApprox(10.0));
        REQUIRE_THAT(bsphere.Center().Y() , RealApprox(20.0));
        REQUIRE_THAT(bsphere.Center().Z() , RealApprox(30.0));
        
        // Radius = (4√3)/2 = 2√3 ≈ 3.464
        Real expected_radius = 4.0 * std::sqrt(3.0) / 2.0;
        REQUIRE_THAT(bsphere.Radius() , RealApprox(expected_radius));
    }

    SECTION("Bounding sphere contains all corners")
    {
        Cube3D cube(6.0);
        BoundingSphere3D bsphere = cube.GetBoundingSphere();
        
        // Corner at maximum distance: (3, 3, 3) from origin
        Pnt3Cart corner(3.0, 3.0, 3.0);
        Real dist_sq = 3.0*3.0 + 3.0*3.0 + 3.0*3.0;  // = 27
        Real dist = std::sqrt(dist_sq);  // = 3√3 ≈ 5.196
        
        REQUIRE(dist <= bsphere.Radius());
        REQUIRE(bsphere.Contains(corner));
    }
}

TEST_CASE("Cube3D::ToString", "[geometry][cube][string]")
{
    SECTION("Cube at origin")
    {
        Cube3D cube(5.0);
        std::string str = cube.ToString();
        
        REQUIRE(str.find("Cube3D") != std::string::npos);
        REQUIRE(str.find("Center") != std::string::npos);
        REQUIRE(str.find("Side") != std::string::npos);
        REQUIRE(str.find("Volume") != std::string::npos);
        REQUIRE(str.find("SurfaceArea") != std::string::npos);
    }

    SECTION("Cube at arbitrary center")
    {
        Pnt3Cart center(10.0, 20.0, 30.0);
        Cube3D cube(3.0, center);
        std::string str = cube.ToString();
        
        REQUIRE(str.find("10") != std::string::npos);
        REQUIRE(str.find("20") != std::string::npos);
        REQUIRE(str.find("30") != std::string::npos);
        REQUIRE(str.find("3") != std::string::npos);
    }
}

TEST_CASE("Cube3D::GetSide", "[geometry][cube][getters]")
{
    SECTION("GetSide returns correct value")
    {
        Cube3D cube(7.5);
        REQUIRE_THAT(cube.GetSide() , RealApprox(7.5));
    }

    SECTION("GetSide for cube at arbitrary center")
    {
        Pnt3Cart center(1.0, 2.0, 3.0);
        Cube3D cube(12.0, center);
        REQUIRE_THAT(cube.GetSide() , RealApprox(12.0));
    }
}

TEST_CASE("Cube3D::IsInside", "[geometry][cube][inside]")
{
    SECTION("Point at center is inside")
    {
        Cube3D cube(10.0);
        REQUIRE(cube.IsInside(Pnt3Cart(0.0, 0.0, 0.0)) == true);
    }

    SECTION("Points on faces are inside (boundary)")
    {
        Cube3D cube(10.0);
        REQUIRE(cube.IsInside(Pnt3Cart(5.0, 0.0, 0.0)) == true);   // Right face
        REQUIRE(cube.IsInside(Pnt3Cart(-5.0, 0.0, 0.0)) == true);  // Left face
        REQUIRE(cube.IsInside(Pnt3Cart(0.0, 5.0, 0.0)) == true);   // Front face
        REQUIRE(cube.IsInside(Pnt3Cart(0.0, -5.0, 0.0)) == true);  // Back face
        REQUIRE(cube.IsInside(Pnt3Cart(0.0, 0.0, 5.0)) == true);   // Top face
        REQUIRE(cube.IsInside(Pnt3Cart(0.0, 0.0, -5.0)) == true);  // Bottom face
    }

    SECTION("Points on edges are inside")
    {
        Cube3D cube(10.0);
        REQUIRE(cube.IsInside(Pnt3Cart(5.0, 5.0, 0.0)) == true);
        REQUIRE(cube.IsInside(Pnt3Cart(5.0, -5.0, 0.0)) == true);
        REQUIRE(cube.IsInside(Pnt3Cart(-5.0, 5.0, 0.0)) == true);
        REQUIRE(cube.IsInside(Pnt3Cart(-5.0, -5.0, 0.0)) == true);
    }

    SECTION("Points at corners are inside")
    {
        Cube3D cube(10.0);
        REQUIRE(cube.IsInside(Pnt3Cart(5.0, 5.0, 5.0)) == true);
        REQUIRE(cube.IsInside(Pnt3Cart(-5.0, -5.0, -5.0)) == true);
        REQUIRE(cube.IsInside(Pnt3Cart(5.0, -5.0, 5.0)) == true);
        REQUIRE(cube.IsInside(Pnt3Cart(-5.0, 5.0, -5.0)) == true);
    }

    SECTION("Points just outside are not inside")
    {
        Cube3D cube(10.0);
        REQUIRE(cube.IsInside(Pnt3Cart(5.1, 0.0, 0.0)) == false);
        REQUIRE(cube.IsInside(Pnt3Cart(0.0, 5.1, 0.0)) == false);
        REQUIRE(cube.IsInside(Pnt3Cart(0.0, 0.0, 5.1)) == false);
    }

    SECTION("Points far outside are not inside")
    {
        Cube3D cube(10.0);
        REQUIRE(cube.IsInside(Pnt3Cart(100.0, 0.0, 0.0)) == false);
        REQUIRE(cube.IsInside(Pnt3Cart(0.0, -100.0, 0.0)) == false);
        REQUIRE(cube.IsInside(Pnt3Cart(0.0, 0.0, 100.0)) == false);
    }

    SECTION("Cube centered at arbitrary location")
    {
        Pnt3Cart center(10.0, 20.0, 30.0);
        Cube3D cube(6.0, center);
        
        REQUIRE(cube.IsInside(center) == true);
        REQUIRE(cube.IsInside(Pnt3Cart(13.0, 20.0, 30.0)) == true);  // On right face
        REQUIRE(cube.IsInside(Pnt3Cart(10.0, 23.0, 30.0)) == true);  // On front face
        REQUIRE(cube.IsInside(Pnt3Cart(10.0, 20.0, 33.0)) == true);  // On top face
        
        REQUIRE(cube.IsInside(Pnt3Cart(14.0, 20.0, 30.0)) == false);  // Outside right
        REQUIRE(cube.IsInside(Pnt3Cart(10.0, 24.0, 30.0)) == false);  // Outside front
        REQUIRE(cube.IsInside(Pnt3Cart(10.0, 20.0, 34.0)) == false);  // Outside top
    }

    SECTION("Interior points are inside")
    {
        Cube3D cube(20.0);
        REQUIRE(cube.IsInside(Pnt3Cart(1.0, 2.0, 3.0)) == true);
        REQUIRE(cube.IsInside(Pnt3Cart(-4.0, 5.0, -6.0)) == true);
    }
}

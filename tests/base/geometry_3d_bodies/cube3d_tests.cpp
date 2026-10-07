#include <catch2/catch_all.hpp>
#include "../../TestPrecision.h"
#include "../../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/base/Geometry/Geometry3DBodies.h>
#endif

using namespace MML;
using namespace MML::Testing;

namespace MML::Tests::Base::Geometry3DBodies::Cube3DTests {

namespace {
    Real BoundaryBoxY1(Real) { return -1.0; }
    Real BoundaryBoxY2(Real) { return 2.0; }
    Real BoundaryBoxZ1(Real, Real) { return 5.0; }
    Real BoundaryBoxZ2(Real, Real) { return 9.0; }
    Real UnitY1(Real) { return 0.0; }
    Real UnitY2(Real) { return 1.0; }
    Real UnitZ1(Real, Real) { return 0.0; }
    Real SlopedZ2(Real x, Real) { return 1.0 + x; }
}

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
        Box3D bbox = cube.GetBoundingBox();
        
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
        Box3D bbox = cube.GetBoundingBox();
        
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
        Box3D bbox = cube.GetBoundingBox();
        
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
        Box3D bbox = cube.GetBoundingBox();
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
        
        // Different computation paths for sqrt(27) may differ by ULPs
        REQUIRE_THAT(dist, RealApprox(bsphere.Radius()));
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

TEST_CASE("BodyWithRectSurfaces derives geometry from a closed quad mesh", "[geometry][mesh][rect]")
{
    Cube3D cube(2.0, Pnt3Cart(3.0, 4.0, 5.0));

    REQUIRE_THAT(cube.BodyWithRectSurfaces::Volume(), RealApprox(8.0));
    REQUIRE_THAT(cube.BodyWithRectSurfaces::SurfaceArea(), RealApprox(24.0));

    const Pnt3Cart center = cube.BodyWithRectSurfaces::GetCenter();
    REQUIRE_THAT(center.X(), RealApprox(3.0));
    REQUIRE_THAT(center.Y(), RealApprox(4.0));
    REQUIRE_THAT(center.Z(), RealApprox(5.0));

    const Box3D box = cube.BodyWithRectSurfaces::GetBoundingBox();
    REQUIRE_THAT(box.Min().X(), RealApprox(2.0));
    REQUIRE_THAT(box.Min().Y(), RealApprox(3.0));
    REQUIRE_THAT(box.Min().Z(), RealApprox(4.0));
    REQUIRE_THAT(box.Max().X(), RealApprox(4.0));
    REQUIRE_THAT(box.Max().Y(), RealApprox(5.0));
    REQUIRE_THAT(box.Max().Z(), RealApprox(6.0));

    REQUIRE(cube.BodyWithRectSurfaces::IsInside(Pnt3Cart(3.0, 4.0, 5.0)));
    REQUIRE(cube.BodyWithRectSurfaces::IsInside(Pnt3Cart(4.0, 4.0, 5.0)));
    REQUIRE_FALSE(cube.BodyWithRectSurfaces::IsInside(Pnt3Cart(4.1, 4.0, 5.0)));
}

TEST_CASE("ISolidBodyWithBoundary numerically derives body geometry", "[geometry][boundary-body]")
{
    SolidBodyWithBoundaryConstDensity box(
        2.0, 4.0, BoundaryBoxY1, BoundaryBoxY2, BoundaryBoxZ1, BoundaryBoxZ2, 1.0);

    REQUIRE_THAT(box.Volume(), RealApprox(24.0));
    REQUIRE_THAT(box.SurfaceArea(), RealApprox(52.0));

    const Pnt3Cart center = box.GetCenter();
    REQUIRE_THAT(center.X(), RealApprox(3.0));
    REQUIRE_THAT(center.Y(), RealApprox(0.5));
    REQUIRE_THAT(center.Z(), RealApprox(7.0));

    const Box3D bounds = box.GetBoundingBox();
    REQUIRE_THAT(bounds.Min().X(), RealApprox(2.0));
    REQUIRE_THAT(bounds.Min().Y(), RealApprox(-1.0));
    REQUIRE_THAT(bounds.Min().Z(), RealApprox(5.0));
    REQUIRE_THAT(bounds.Max().X(), RealApprox(4.0));
    REQUIRE_THAT(bounds.Max().Y(), RealApprox(2.0));
    REQUIRE_THAT(bounds.Max().Z(), RealApprox(9.0));
    REQUIRE_THAT(box.GetBoundingSphere().Radius(), RealApprox(std::sqrt(7.25)));

    REQUIRE(box.IsInside(Pnt3Cart(3.0, 0.0, 7.0)));
    REQUIRE_FALSE(box.IsInside(Pnt3Cart(4.1, 0.0, 7.0)));
}

TEST_CASE("ISolidBodyWithBoundary includes graph slopes in surface area", "[geometry][boundary-body]")
{
    SolidBodyWithBoundaryConstDensity wedge(
        0.0, 1.0, UnitY1, UnitY2, UnitZ1, SlopedZ2, 1.0);

    REQUIRE_THAT(wedge.Volume(), RealApprox(1.5));
    REQUIRE_THAT(wedge.SurfaceArea(), RealApprox(7.0 + std::sqrt(2.0)).epsilon(TOL(1e-9, 2e-5)).margin(TOL(1e-9, 2e-5)));
}

} // namespace MML::Tests::Base::Geometry3DBodies::Cube3DTests

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

namespace MML::Tests::Base::Geometry3DBodies::PyramidEquilateral3DTests {

TEST_CASE("PyramidEquilateral3D::Constructor", "[geometry][pyramidequilateral][constructor]")
{
    SECTION("Equilateral pyramid with base side 6")
    {
        PyramidEquilateral3D pyramid(6.0);
        
        Real expected_height = MML::Real(6.0) / std::sqrt(MML::Real(2.0));
        REQUIRE_THAT(pyramid.getHeight() , RealApprox(expected_height));
        REQUIRE_THAT(pyramid.getBaseSize() , RealApprox(6.0));
    }

    SECTION("Equilateral pyramid at arbitrary position")
    {
        Vec3Cart baseCenter(10.0, 20.0, 30.0);
        PyramidEquilateral3D pyramid(9.0, baseCenter);
        
        Real expected_height = MML::Real(9.0) / std::sqrt(MML::Real(2.0));
        REQUIRE_THAT(pyramid.getHeight() , RealApprox(expected_height));
        REQUIRE_THAT(pyramid.getBaseSize() , RealApprox(9.0));
        
        Vec3Cart center = pyramid.getBaseCenter();
        REQUIRE_THAT(center.X() , RealApprox(10.0));
        REQUIRE_THAT(center.Y() , RealApprox(20.0));
        REQUIRE_THAT(center.Z() , RealApprox(30.0));
    }
}

TEST_CASE("PyramidEquilateral3D::Volume", "[geometry][pyramidequilateral][volume]")
{
    const Real tolerance = TOL(Tolerance::Standard, Real(6) * Precision::Eps);
    SECTION("Equilateral pyramid inherits volume calculation")
    {
        PyramidEquilateral3D pyramid(6.0);
        Real expected = MML::Real(36.0) * std::sqrt(MML::Real(2.0));
        REQUIRE_THAT(pyramid.Volume(), RealApprox(expected).epsilon(tolerance));
    }

    SECTION("Equilateral pyramid with base 3")
    {
        PyramidEquilateral3D pyramid(3.0);
        Real expected = MML::Real(9.0) / std::sqrt(MML::Real(2.0));
        REQUIRE_THAT(pyramid.Volume(), RealApprox(expected).epsilon(tolerance));
    }
}

TEST_CASE("PyramidEquilateral3D::SurfaceArea", "[geometry][pyramidequilateral][surface]")
{
    SECTION("Surface area calculation")
    {
        PyramidEquilateral3D pyramid(6.0);
        // Square plus four equilateral triangles.
        Real expected = MML::Real(36.0) + MML::Real(36.0) * std::sqrt(MML::Real(3.0));
        REQUIRE_THAT(pyramid.SurfaceArea() , RealApprox(expected));
    }
}

TEST_CASE("PyramidEquilateral3D::GetCenter", "[geometry][pyramidequilateral][center]")
{
    SECTION("Center at h/4 from base")
    {
        PyramidEquilateral3D pyramid(12.0);
        Pnt3Cart center = pyramid.getCenter();
        
        Real expected_z = MML::Real(12.0) / (MML::Real(4.0) * std::sqrt(MML::Real(2.0)));
        REQUIRE_THAT(center.X() , RealApprox(0.0));
        REQUIRE_THAT(center.Y() , RealApprox(0.0));
        REQUIRE_THAT(center.Z() , RealApprox(expected_z));
    }

    SECTION("Center at arbitrary position")
    {
        Vec3Cart baseCenter(5.0, 10.0, 15.0);
        PyramidEquilateral3D pyramid(6.0, baseCenter);
        Pnt3Cart center = pyramid.getCenter();
        
        Real offset_z = MML::Real(6.0) / (MML::Real(4.0) * std::sqrt(MML::Real(2.0)));
        REQUIRE_THAT(center.X() , RealApprox(5.0));
        REQUIRE_THAT(center.Y() , RealApprox(10.0));
        REQUIRE_THAT(center.Z() , RealApprox(MML::Real(15.0) + offset_z));
    }
}

TEST_CASE("PyramidEquilateral3D::GetBoundingBox", "[geometry][pyramidequilateral][bounding]")
{
    SECTION("Bounding box from base to apex")
    {
        PyramidEquilateral3D pyramid(8.0);
        Box3D bbox = pyramid.getBoundingBox();
        
        Real h = MML::Real(8.0) / std::sqrt(MML::Real(2.0));
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
        BoundingSphere3D bsphere = pyramid.getBoundingSphere();
        
        Real h = MML::Real(6.0) / std::sqrt(MML::Real(2.0));
        
        // Check base corners
        REQUIRE(bsphere.Contains(Pnt3Cart(3.0, 3.0, 0.0)));
        REQUIRE(bsphere.Contains(Pnt3Cart(-3.0, -3.0, 0.0)));
        
        // Check apex
        REQUIRE(bsphere.Contains(Pnt3Cart(0.0, 0.0, h)));
    }

    SECTION("Bounding sphere centered at centroid")
    {
        PyramidEquilateral3D pyramid(12.0);
        BoundingSphere3D bsphere = pyramid.getBoundingSphere();
        
        Real h = MML::Real(12.0) / std::sqrt(MML::Real(2.0));
        Real centroid_z = h / MML::Real(4.0);
        
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
        std::string str = pyramid.to_string();
        
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
        std::string str = pyramid.to_string();
        
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
        Real h = MML::Real(12.0) / std::sqrt(MML::Real(2.0));
        Real centroid_z = h / MML::Real(4.0);
        REQUIRE(pyramid.IsInside(Pnt3Cart(0.0, 0.0, centroid_z)) == true);
    }

    SECTION("Apex is inside")
    {
        PyramidEquilateral3D pyramid(9.0);
        Real h = MML::Real(9.0) / std::sqrt(MML::Real(2.0));
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
        Real h = MML::Real(6.0) / std::sqrt(MML::Real(2.0));
        
        REQUIRE(pyramid.IsInside(Pnt3Cart(MML::Real(0.0), MML::Real(0.0), -MML::Real(0.1))) == false);  // Below base
        REQUIRE(pyramid.IsInside(Pnt3Cart(MML::Real(0.0), MML::Real(0.0), h + MML::Real(0.1))) == false);  // Above apex
        REQUIRE(pyramid.IsInside(Pnt3Cart(MML::Real(3.1), MML::Real(0.0), MML::Real(0.0))) == false);  // Outside base
    }

    SECTION("Equilateral pyramid at arbitrary position")
    {
        Vec3Cart baseCenter(10.0, 20.0, 30.0);
        PyramidEquilateral3D pyramid(6.0, baseCenter);
        Real h = MML::Real(6.0) / std::sqrt(MML::Real(2.0));
        
        // Base center
        REQUIRE(pyramid.IsInside(Pnt3Cart(10.0, 20.0, 30.0)) == true);
        
        // Near apex (slightly below to avoid boundary condition)
        REQUIRE(pyramid.IsInside(Pnt3Cart(MML::Real(10.0), MML::Real(20.0), MML::Real(30.0) + h - MML::Real(0.01))) == true);
        
        // Outside
        REQUIRE(pyramid.IsInside(Pnt3Cart(MML::Real(10.0), MML::Real(20.0), MML::Real(29.9))) == false);
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
        Pnt3Cart center = pyramid.getCenter();
        Box3D bbox = pyramid.getBoundingBox();
        BoundingSphere3D bsphere = pyramid.getBoundingSphere();
        std::string str = pyramid.to_string();
        REQUIRE(center.Z() >= bbox.Min().Z());
        REQUIRE(center.Z() <= bbox.Max().Z());
        REQUIRE(bsphere.Contains(center));
        
        REQUIRE(volume > 0);
        REQUIRE(surfaceArea > 0);
        REQUIRE(str.length() > 0);
    }

    SECTION("Getters work correctly")
    {
        PyramidEquilateral3D pyramid(9.0);
        
        REQUIRE_THAT(pyramid.getBaseSize() , RealApprox(9.0));
        REQUIRE_THAT(pyramid.getHeight() , RealApprox(MML::Real(9.0) / std::sqrt(MML::Real(2.0))).epsilon(TOL(1e-10, 1e-5)));
        
        Vec3Cart baseCenter = pyramid.getBaseCenter();
        REQUIRE_THAT(baseCenter.X() , RealApprox(0.0));
        REQUIRE_THAT(baseCenter.Y() , RealApprox(0.0));
        REQUIRE_THAT(baseCenter.Z() , RealApprox(0.0));
    }
}

TEST_CASE("PyramidEquilateral3D::GeometricProperties", "[geometry][pyramidequilateral][properties]")
{
    SECTION("Height-to-base ratio is 1/sqrt(2)")
    {
        for (Real a : {MML::Real(3.0), MML::Real(6.0), MML::Real(9.0), MML::Real(12.0)}) {
            PyramidEquilateral3D pyramid(a);
            Real expected_ratio = MML::Real(1.0) / std::sqrt(MML::Real(2.0));
            Real actual_ratio = pyramid.getHeight() / pyramid.getBaseSize();
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

} // namespace MML::Tests::Base::Geometry3DBodies::PyramidEquilateral3DTests

#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Geometry3DBodies.h"
#endif

using namespace MML;
using namespace MML::Testing;

namespace MML::Tests::Base::BoundingVolumesTests
{
	TEST_CASE("BoundingBox3D::Construction", "[BoundingBox3D]") {
			TEST_PRECISION_INFO();
		SECTION("From min/max points") {
			Pnt3Cart min(1, 2, 3);
			Pnt3Cart max(4, 6, 9);
			BoundingBox3D box(min, max);
			
			REQUIRE_THAT(box.Min().X() , RealApprox(REAL(1.0)));
			REQUIRE_THAT(box.Min().Y() , RealApprox(REAL(2.0)));
			REQUIRE_THAT(box.Min().Z() , RealApprox(REAL(3.0)));
			REQUIRE_THAT(box.Max().X() , RealApprox(REAL(4.0)));
			REQUIRE_THAT(box.Max().Y() , RealApprox(REAL(6.0)));
			REQUIRE_THAT(box.Max().Z() , RealApprox(REAL(9.0)));
		}

		SECTION("Validation - min must be <= max") {
			Pnt3Cart min(5, 5, 5);
			Pnt3Cart max(1, 1, 1);  // Invalid: max < min
			REQUIRE_THROWS_AS(BoundingBox3D(min, max), std::invalid_argument);
		}
	}

	TEST_CASE("BoundingBox3D::Dimensions", "[BoundingBox3D]") {
			TEST_PRECISION_INFO();
		Pnt3Cart min(1, 2, 3);
		Pnt3Cart max(4, 7, 13);
		BoundingBox3D box(min, max);

		REQUIRE_THAT(box.Width() , RealApprox(REAL(3.0)));   // 4 - 1
		REQUIRE_THAT(box.Height() , RealApprox(REAL(5.0)));  // 7 - 2
		REQUIRE_THAT(box.Depth() , RealApprox(REAL(10.0)));  // 13 - 3
		REQUIRE_THAT(box.Volume() , RealApprox(REAL(150.0))); // 3 * 5 * 10
	}

	TEST_CASE("BoundingBox3D::Center", "[BoundingBox3D]") {
			TEST_PRECISION_INFO();
		Pnt3Cart min(0, 0, 0);
		Pnt3Cart max(10, 20, 30);
		BoundingBox3D box(min, max);

		Pnt3Cart center = box.Center();
		REQUIRE_THAT(center.X() , RealApprox(REAL(5.0)));
		REQUIRE_THAT(center.Y() , RealApprox(REAL(10.0)));
		REQUIRE_THAT(center.Z() , RealApprox(REAL(15.0)));
	}

	TEST_CASE("BoundingBox3D::Contains", "[BoundingBox3D]") {
			TEST_PRECISION_INFO();
		Pnt3Cart min(0, 0, 0);
		Pnt3Cart max(10, 10, 10);
		BoundingBox3D box(min, max);

		SECTION("Point inside") {
			REQUIRE(box.Contains(Pnt3Cart(5, 5, 5)));
			REQUIRE(box.Contains(Pnt3Cart(0, 0, 0)));  // Min corner
			REQUIRE(box.Contains(Pnt3Cart(10, 10, 10)));  // Max corner
		}

		SECTION("Point outside") {
			REQUIRE_FALSE(box.Contains(Pnt3Cart(-1, 5, 5)));
			REQUIRE_FALSE(box.Contains(Pnt3Cart(5, 11, 5)));
			REQUIRE_FALSE(box.Contains(Pnt3Cart(5, 5, 15)));
		}
	}

	TEST_CASE("BoundingBox3D::Intersects", "[BoundingBox3D]") {
			TEST_PRECISION_INFO();
		BoundingBox3D box1(Pnt3Cart(0, 0, 0), Pnt3Cart(10, 10, 10));

		SECTION("Overlapping boxes") {
			BoundingBox3D box2(Pnt3Cart(5, 5, 5), Pnt3Cart(15, 15, 15));
			REQUIRE(box1.Intersects(box2));
			REQUIRE(box2.Intersects(box1));  // Symmetric
		}

		SECTION("Separated boxes") {
			BoundingBox3D box2(Pnt3Cart(15, 15, 15), Pnt3Cart(25, 25, 25));
			REQUIRE_FALSE(box1.Intersects(box2));
		}
	}

	// ============================================================================
	// BoundingSphere3D Tests
	// ============================================================================

	TEST_CASE("BoundingSphere3D::Construction", "[BoundingSphere3D]") {
			TEST_PRECISION_INFO();
		SECTION("From center and radius") {
			Pnt3Cart center(1, 2, 3);
			BoundingSphere3D sphere(center, REAL(5.0));
			
			REQUIRE_THAT(sphere.Center().X() , RealApprox(REAL(1.0)));
			REQUIRE_THAT(sphere.Center().Y() , RealApprox(REAL(2.0)));
			REQUIRE_THAT(sphere.Center().Z() , RealApprox(REAL(3.0)));
			REQUIRE_THAT(sphere.Radius() , RealApprox(REAL(5.0)));
		}

		SECTION("Validation - radius must be non-negative") {
			REQUIRE_THROWS_AS(BoundingSphere3D(Pnt3Cart(0, 0, 0), -REAL(1.0)), std::invalid_argument);
		}
	}

	TEST_CASE("BoundingSphere3D::Volume", "[BoundingSphere3D]") {
			TEST_PRECISION_INFO();
		BoundingSphere3D sphere(Pnt3Cart(0, 0, 0), REAL(5.0));
		
		// Volume = (4/3) * π * r³ = (4/3) * π * 125
		Real expected = (REAL(4.0) / REAL(3.0)) * Constants::PI * REAL(125.0);
		REQUIRE_THAT(sphere.Volume() , RealApprox(expected));
	}

	TEST_CASE("BoundingSphere3D::Contains", "[BoundingSphere3D]") {
			TEST_PRECISION_INFO();
		BoundingSphere3D sphere(Pnt3Cart(0, 0, 0), REAL(5.0));

		SECTION("Point inside") {
			REQUIRE(sphere.Contains(Pnt3Cart(0, 0, 0)));  // Center
			REQUIRE(sphere.Contains(Pnt3Cart(3, 0, 0)));  // Inside
			REQUIRE(sphere.Contains(Pnt3Cart(5, 0, 0)));  // On surface
		}

		SECTION("Point outside") {
			REQUIRE_FALSE(sphere.Contains(Pnt3Cart(6, 0, 0)));
			REQUIRE_FALSE(sphere.Contains(Pnt3Cart(10, 10, 10)));
		}
	}

	TEST_CASE("BoundingSphere3D::Intersects", "[BoundingSphere3D]") {
			TEST_PRECISION_INFO();
		BoundingSphere3D sphere1(Pnt3Cart(0, 0, 0), REAL(5.0));

		SECTION("Overlapping spheres") {
			BoundingSphere3D sphere2(Pnt3Cart(7, 0, 0), REAL(5.0));
			REQUIRE(sphere1.Intersects(sphere2));  // Distance = 7, sum of radii = 10
		}

		SECTION("Separated spheres") {
			BoundingSphere3D sphere2(Pnt3Cart(15, 0, 0), REAL(5.0));
			REQUIRE_FALSE(sphere1.Intersects(sphere2));  // Distance = 15, sum = 10
		}
	}
}

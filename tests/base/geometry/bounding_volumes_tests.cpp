#include <catch2/catch_all.hpp>
#include "../../TestPrecision.h"
#include "../../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "mml/base/Geometry/Geometry3DBodies.h"
#endif

using namespace MML;
using namespace MML::Testing;

namespace MML::Tests::Base::BoundingVolumesTests
{
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

	// ============================================================================
	// Box3D Tests
	// ============================================================================

	TEST_CASE("Box3D::Construction", "[Box3D]") {
		TEST_PRECISION_INFO();
		SECTION("Default constructor") {
			Box3D box;
			REQUIRE_THAT(box.Min().X(), RealApprox(REAL(0.0)));
			REQUIRE_THAT(box.Max().X(), RealApprox(REAL(1.0)));
			REQUIRE_THAT(box.Volume(), RealApprox(REAL(1.0)));
		}

		SECTION("From min/max points") {
			Pnt3Cart min(1, 2, 3);
			Pnt3Cart max(4, 6, 9);
			Box3D box(min, max);

			REQUIRE_THAT(box.Min().X(), RealApprox(REAL(1.0)));
			REQUIRE_THAT(box.Min().Y(), RealApprox(REAL(2.0)));
			REQUIRE_THAT(box.Min().Z(), RealApprox(REAL(3.0)));
			REQUIRE_THAT(box.Max().X(), RealApprox(REAL(4.0)));
			REQUIRE_THAT(box.Max().Y(), RealApprox(REAL(6.0)));
			REQUIRE_THAT(box.Max().Z(), RealApprox(REAL(9.0)));
		}

		SECTION("From coordinates - auto-normalizes") {
			Box3D box(5, 6, 7, 1, 2, 3);  // max < min, should swap
			REQUIRE_THAT(box.MinX(), RealApprox(REAL(1.0)));
			REQUIRE_THAT(box.MaxX(), RealApprox(REAL(5.0)));
		}
	}

	TEST_CASE("Box3D::FactoryMethods", "[Box3D]") {
		TEST_PRECISION_INFO();
		SECTION("FromCenterAndSize") {
			Box3D box = Box3D::FromCenterAndSize(Pnt3Cart(5, 5, 5), 4, 6, 8);
			REQUIRE_THAT(box.MinX(), RealApprox(REAL(3.0)));  // 5 - 2
			REQUIRE_THAT(box.MaxX(), RealApprox(REAL(7.0)));  // 5 + 2
			REQUIRE_THAT(box.Width(), RealApprox(REAL(4.0)));
			REQUIRE_THAT(box.Height(), RealApprox(REAL(6.0)));
			REQUIRE_THAT(box.Depth(), RealApprox(REAL(8.0)));
		}

		SECTION("FromPoints") {
			std::vector<Pnt3Cart> points = {
				Pnt3Cart(1, 5, 3), Pnt3Cart(4, 2, 7), Pnt3Cart(0, 8, 1)
			};
			Box3D box = Box3D::FromPoints(points);
			REQUIRE_THAT(box.MinX(), RealApprox(REAL(0.0)));
			REQUIRE_THAT(box.MinY(), RealApprox(REAL(2.0)));
			REQUIRE_THAT(box.MinZ(), RealApprox(REAL(1.0)));
			REQUIRE_THAT(box.MaxX(), RealApprox(REAL(4.0)));
			REQUIRE_THAT(box.MaxY(), RealApprox(REAL(8.0)));
			REQUIRE_THAT(box.MaxZ(), RealApprox(REAL(7.0)));
		}

		SECTION("UnitCube") {
			Box3D box = Box3D::UnitCube();
			REQUIRE_THAT(box.Volume(), RealApprox(REAL(1.0)));
		}

		SECTION("CenteredUnitCube") {
			Box3D box = Box3D::CenteredUnitCube();
			REQUIRE_THAT(box.Center().X(), RealApprox(REAL(0.0)));
			REQUIRE_THAT(box.Center().Y(), RealApprox(REAL(0.0)));
			REQUIRE_THAT(box.Center().Z(), RealApprox(REAL(0.0)));
			REQUIRE_THAT(box.Volume(), RealApprox(REAL(1.0)));
		}
	}

	TEST_CASE("Box3D::Dimensions", "[Box3D]") {
		TEST_PRECISION_INFO();
		Box3D box(1, 2, 3, 4, 7, 13);

		REQUIRE_THAT(box.Width(), RealApprox(REAL(3.0)));   // 4 - 1
		REQUIRE_THAT(box.Height(), RealApprox(REAL(5.0)));  // 7 - 2
		REQUIRE_THAT(box.Depth(), RealApprox(REAL(10.0)));  // 13 - 3

		Vec3Cart size = box.Size();
		REQUIRE_THAT(size.X(), RealApprox(REAL(3.0)));
		REQUIRE_THAT(size.Y(), RealApprox(REAL(5.0)));
		REQUIRE_THAT(size.Z(), RealApprox(REAL(10.0)));

		Vec3Cart half = box.HalfExtents();
		REQUIRE_THAT(half.X(), RealApprox(REAL(1.5)));
		REQUIRE_THAT(half.Y(), RealApprox(REAL(2.5)));
		REQUIRE_THAT(half.Z(), RealApprox(REAL(5.0)));

		// Diagonal = sqrt(9 + 25 + 100) = sqrt(134)
		REQUIRE_THAT(box.Diagonal(), RealApprox(std::sqrt(REAL(134.0))));
	}

	TEST_CASE("Box3D::GeometricProperties", "[Box3D]") {
		TEST_PRECISION_INFO();
		Box3D box(0, 0, 0, 2, 3, 5);  // 2 × 3 × 5 box

		REQUIRE_THAT(box.Volume(), RealApprox(REAL(30.0)));  // 2 * 3 * 5 = 30
		// Surface = 2*(2*3 + 3*5 + 5*2) = 2*(6 + 15 + 10) = 62
		REQUIRE_THAT(box.SurfaceArea(), RealApprox(REAL(62.0)));

		Pnt3Cart center = box.Center();
		REQUIRE_THAT(center.X(), RealApprox(REAL(1.0)));
		REQUIRE_THAT(center.Y(), RealApprox(REAL(1.5)));
		REQUIRE_THAT(center.Z(), RealApprox(REAL(2.5)));
	}

	TEST_CASE("Box3D::Corners", "[Box3D]") {
		TEST_PRECISION_INFO();
		Box3D box(0, 0, 0, 1, 2, 3);

		REQUIRE_THAT(box.Corner000().X(), RealApprox(REAL(0.0)));
		REQUIRE_THAT(box.Corner000().Y(), RealApprox(REAL(0.0)));
		REQUIRE_THAT(box.Corner000().Z(), RealApprox(REAL(0.0)));

		REQUIRE_THAT(box.Corner111().X(), RealApprox(REAL(1.0)));
		REQUIRE_THAT(box.Corner111().Y(), RealApprox(REAL(2.0)));
		REQUIRE_THAT(box.Corner111().Z(), RealApprox(REAL(3.0)));

		std::vector<Pnt3Cart> corners = box.Corners();
		REQUIRE(corners.size() == 8);
	}

	TEST_CASE("Box3D::Contains", "[Box3D]") {
		TEST_PRECISION_INFO();
		Box3D box(0, 0, 0, 10, 10, 10);

		SECTION("Point inside") {
			REQUIRE(box.Contains(Pnt3Cart(5, 5, 5)));
			REQUIRE(box.Contains(Pnt3Cart(0, 0, 0)));   // Min corner
			REQUIRE(box.Contains(Pnt3Cart(10, 10, 10))); // Max corner
		}

		SECTION("Point outside") {
			REQUIRE_FALSE(box.Contains(Pnt3Cart(-1, 5, 5)));
			REQUIRE_FALSE(box.Contains(Pnt3Cart(5, 11, 5)));
			REQUIRE_FALSE(box.Contains(Pnt3Cart(5, 5, 15)));
		}

		SECTION("Box contains another box") {
			Box3D inner(2, 2, 2, 8, 8, 8);
			REQUIRE(box.Contains(inner));

			Box3D outer(0, 0, 0, 20, 20, 20);
			REQUIRE_FALSE(box.Contains(outer));
		}
	}

	TEST_CASE("Box3D::DistanceAndClosestPoint", "[Box3D]") {
		TEST_PRECISION_INFO();
		Box3D box(0, 0, 0, 10, 10, 10);

		SECTION("Point inside - distance is 0") {
			REQUIRE_THAT(box.DistanceToPoint(Pnt3Cart(5, 5, 5)), RealApprox(REAL(0.0)));
		}

		SECTION("Point outside on axis") {
			REQUIRE_THAT(box.DistanceToPoint(Pnt3Cart(15, 5, 5)), RealApprox(REAL(5.0)));
		}

		SECTION("Closest point") {
			Pnt3Cart closest = box.ClosestPoint(Pnt3Cart(15, 5, 5));
			REQUIRE_THAT(closest.X(), RealApprox(REAL(10.0)));
			REQUIRE_THAT(closest.Y(), RealApprox(REAL(5.0)));
			REQUIRE_THAT(closest.Z(), RealApprox(REAL(5.0)));
		}
	}

	TEST_CASE("Box3D::Intersects", "[Box3D]") {
		TEST_PRECISION_INFO();
		Box3D box1(0, 0, 0, 10, 10, 10);

		SECTION("Overlapping boxes") {
			Box3D box2(5, 5, 5, 15, 15, 15);
			REQUIRE(box1.Intersects(box2));
			REQUIRE(box2.Intersects(box1));  // Symmetric
		}

		SECTION("Separated boxes") {
			Box3D box2(15, 15, 15, 25, 25, 25);
			REQUIRE_FALSE(box1.Intersects(box2));
		}

		SECTION("Touching boxes") {
			Box3D box2(10, 0, 0, 20, 10, 10);  // Touching at X=10
			REQUIRE(box1.Intersects(box2));
		}
	}

	TEST_CASE("Box3D::IntersectionAndUnion", "[Box3D]") {
		TEST_PRECISION_INFO();
		Box3D box1(0, 0, 0, 10, 10, 10);
		Box3D box2(5, 5, 5, 15, 15, 15);

		SECTION("Intersection") {
			Box3D inter = box1.Intersection(box2);
			REQUIRE_THAT(inter.MinX(), RealApprox(REAL(5.0)));
			REQUIRE_THAT(inter.MinY(), RealApprox(REAL(5.0)));
			REQUIRE_THAT(inter.MinZ(), RealApprox(REAL(5.0)));
			REQUIRE_THAT(inter.MaxX(), RealApprox(REAL(10.0)));
			REQUIRE_THAT(inter.MaxY(), RealApprox(REAL(10.0)));
			REQUIRE_THAT(inter.MaxZ(), RealApprox(REAL(10.0)));
			REQUIRE_THAT(inter.Volume(), RealApprox(REAL(125.0)));  // 5×5×5
		}

		SECTION("Union") {
			Box3D uni = box1.Union(box2);
			REQUIRE_THAT(uni.MinX(), RealApprox(REAL(0.0)));
			REQUIRE_THAT(uni.MinY(), RealApprox(REAL(0.0)));
			REQUIRE_THAT(uni.MinZ(), RealApprox(REAL(0.0)));
			REQUIRE_THAT(uni.MaxX(), RealApprox(REAL(15.0)));
			REQUIRE_THAT(uni.MaxY(), RealApprox(REAL(15.0)));
			REQUIRE_THAT(uni.MaxZ(), RealApprox(REAL(15.0)));
		}
	}

	TEST_CASE("Box3D::Transformations", "[Box3D]") {
		TEST_PRECISION_INFO();
		Box3D box(0, 0, 0, 2, 2, 2);

		SECTION("Translated") {
			Box3D moved = box.Translated(Vec3Cart(5, 5, 5));
			REQUIRE_THAT(moved.MinX(), RealApprox(REAL(5.0)));
			REQUIRE_THAT(moved.MaxX(), RealApprox(REAL(7.0)));
		}

		SECTION("Scaled") {
			Box3D scaled = box.Scaled(2.0);
			REQUIRE_THAT(scaled.Width(), RealApprox(REAL(4.0)));
			REQUIRE_THAT(scaled.Volume(), RealApprox(REAL(64.0)));  // 4×4×4
			// Center should remain the same
			REQUIRE_THAT(scaled.Center().X(), RealApprox(box.Center().X()));
		}

		SECTION("Expanded") {
			Box3D expanded = box.Expanded(1.0);
			REQUIRE_THAT(expanded.MinX(), RealApprox(REAL(-1.0)));
			REQUIRE_THAT(expanded.MaxX(), RealApprox(REAL(3.0)));
			REQUIRE_THAT(expanded.Width(), RealApprox(REAL(4.0)));
		}
	}

	TEST_CASE("Box3D::IsDegenerate", "[Box3D]") {
		TEST_PRECISION_INFO();
		SECTION("Normal box") {
			Box3D box(0, 0, 0, 1, 1, 1);
			REQUIRE_FALSE(box.IsDegenerate());
		}

		SECTION("Zero width") {
			Box3D box(0, 0, 0, 0, 1, 1);
			REQUIRE(box.IsDegenerate());
		}

		SECTION("All zero") {
			Box3D box(5, 5, 5, 5, 5, 5);
			REQUIRE(box.IsDegenerate());
		}
	}}

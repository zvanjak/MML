#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "../../TestPrecision.h"

#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/base/Geometry/Geometry.h>
#include <mml/base/Geometry/Geometry3DBodies.h>
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <utility>
#include <vector>

using namespace MML;

namespace {
	void CheckClose(Real actual, Real expected) {
		CHECK_THAT(actual, Catch::Matchers::WithinAbs(
			expected, Precision::EqualityTolerance * std::max(Real(1), std::abs(expected))));
	}

	void CheckPoint(const Pnt3Cart& actual, const Pnt3Cart& expected) {
		CheckClose(actual.X(), expected.X());
		CheckClose(actual.Y(), expected.Y());
		CheckClose(actual.Z(), expected.Z());
	}

	Pnt3Cart SurfacePoint(const RectSurface3D& surface, Real u, Real w) {
		const auto point = surface(u, w);
		return Pnt3Cart(point[0], point[1], point[2]);
	}

	void CheckClosedOutwardMesh(const BodyWithTriangleSurfaces& body, const Pnt3Cart& interior) {
		std::vector<Pnt3Cart> vertices;
		std::map<std::pair<std::size_t, std::size_t>, std::pair<int, int>> edges;
		auto vertexIndex = [&](const Pnt3Cart& point) {
			for (std::size_t index = 0; index < vertices.size(); ++index)
				if (vertices[index].Dist(point) <= Precision::GeometryEpsilon)
					return index;
			vertices.push_back(point);
			return vertices.size() - 1;
		};
		for (int index = 0; index < body.getSurfaceCount(); ++index) {
			const auto& face = body.getSurface(index);
			const std::array<Pnt3Cart, 3> points{face.Pnt1(), face.Pnt2(), face.Pnt3()};
			const Vec3Cart normal = VectorProduct(Vec3Cart(points[0], points[1]), Vec3Cart(points[0], points[2]));
			CHECK(ScalarProduct(normal, Vec3Cart(interior, face.Centroid())) > 0);
			for (int edge = 0; edge < 3; ++edge) {
				const auto from = vertexIndex(points[edge]);
				const auto to = vertexIndex(points[(edge + 1) % 3]);
				REQUIRE(from != to);
				auto& counts = edges[{std::min(from, to), std::max(from, to)}];
				++counts.first;
				counts.second += from < to ? 1 : -1;
			}
		}
		for (const auto& entry : edges) {
			CHECK(entry.second.first == 2);
			CHECK(entry.second.second == 0);
		}
		CHECK(vertices.size() + body.getSurfaceCount() == edges.size() + 2);
	}
}

TEST_CASE("Plane intersections satisfy both equations away from the origin", "[geometry-formulas][Plane3D]") {
	for (const Pnt3Cart point : {Pnt3Cart(1, 2, 0), Pnt3Cart(-3, 4, 5)}) {
		for (const auto& normals : {
				 std::array<Vec3Cart, 2>{Vec3Cart(1, 0, 0), Vec3Cart(0, 1, 0)},
				 std::array<Vec3Cart, 2>{Vec3Cart(2, -1, 3), Vec3Cart(-1, 4, 2)}}) {
			const Plane3D first(point, normals[0]);
			const Plane3D second(point, normals[1]);
			for (bool reverse : {false, true}) {
				Line3D line;
				REQUIRE((reverse ? second.IntersectionWithPlane(first, line)
								 : first.IntersectionWithPlane(second, line)));
				for (Real offset : {Real(-7), Real(0), Real(11)}) {
					const Pnt3Cart sample = line.StartPoint() + offset * line.Direction();
					CheckClose(ScalarProduct(first.Normal(), Vec3Cart(point, sample)), 0);
					CheckClose(ScalarProduct(second.Normal(), Vec3Cart(point, sample)), 0);
				}
				CheckClose(ScalarProduct(first.Normal(), line.Direction()), 0);
				CheckClose(ScalarProduct(second.Normal(), line.Direction()), 0);
			}
		}
	}
	Line3D line;
	REQUIRE(Plane3D(2, 3, 4).IntersectionWithPlane(Plane3D(3, 5, 7), line));
	CheckClose(line.StartPoint().X() / 2 + line.StartPoint().Y() / 3 + line.StartPoint().Z() / 4, 1);
	CheckClose(line.StartPoint().X() / 3 + line.StartPoint().Y() / 5 + line.StartPoint().Z() / 7, 1);
	CHECK_FALSE(Plane3D::XYPlane().IntersectionWithPlane(
		Plane3D(Pnt3Cart(0, 0, 2), Vec3Cart(0, 0, -1)), line));
}

TEST_CASE("Tetrahedron edge volumes agree with vertex determinants under relabeling", "[geometry-formulas][Tetrahedron]") {
	for (const auto& vertices : {
			 std::array<Pnt3Cart, 4>{Pnt3Cart(0, 0, 0), Pnt3Cart(2, 0, 0), Pnt3Cart(0, 3, 0), Pnt3Cart(0, 0, 4)},
			 std::array<Pnt3Cart, 4>{Pnt3Cart(1, -2, 3), Pnt3Cart(4, -1, 2), Pnt3Cart(0, 2, 4), Pnt3Cart(2, 1, 8)}}) {
		const Real expected = std::abs(ScalarProduct(Vec3Cart(vertices[0], vertices[1]),
			VectorProduct(Vec3Cart(vertices[0], vertices[2]), Vec3Cart(vertices[0], vertices[3])))) / 6;
		const int scaleExponent = std::min(100, (std::numeric_limits<Real>::max_exponent - 16) / 3);
		for (Real scale : {std::ldexp(Real(1), -scaleExponent), Real(0.5), Real(1), Real(8),
				 std::ldexp(Real(1), scaleExponent)}) {
			REQUIRE(std::isnormal(expected * scale * scale * scale));
			std::array<int, 4> order{0, 1, 2, 3};
			do {
				const auto& a = vertices[order[0]];
				const auto& b = vertices[order[1]];
				const auto& c = vertices[order[2]];
				const auto& d = vertices[order[3]];
				const Tetrahedron shape(scale * a.Dist(b), scale * a.Dist(c), scale * a.Dist(d),
					scale * b.Dist(c), scale * b.Dist(d), scale * c.Dist(d));
				CheckClose(shape.Volume() / scale / scale / scale, expected);
			} while (std::next_permutation(order.begin(), order.end()));
		}
	}
	CheckClose(Tetrahedron().Volume(), 0);
	CheckClose(Tetrahedron(3, 4, 0, 5, 3, 4).Volume(), 0);
}

TEST_CASE("Tetrahedron volume rejects nonrealizable and nonfinite edges", "[geometry-formulas][Tetrahedron]") {
	const Real nan = std::numeric_limits<Real>::quiet_NaN();
	const Real inf = std::numeric_limits<Real>::infinity();
	for (Real invalid : {Real(-1), nan, inf})
		for (int edge = 0; edge < 6; ++edge) {
			std::array<Real, 6> lengths{1, 1, 1, 1, 1, 1};
			lengths[edge] = invalid;
			CHECK_THROWS_AS(Tetrahedron(lengths[0], lengths[1], lengths[2], lengths[3], lengths[4], lengths[5]).Volume(),
				GeometryError);
		}

	CHECK_THROWS_AS(Tetrahedron(1, 1, 1, 3, 1, 1).Volume(), GeometryError);
	CHECK_THROWS_AS(Tetrahedron::Regular(std::numeric_limits<Real>::max()).Volume(), GeometryError);
	// Every face is a valid triangle, but these three unit edges cannot all subtend 120 degrees.
	CHECK_THROWS_AS(Tetrahedron(1, 1, 1, std::sqrt(Real(3)), std::sqrt(Real(3)), std::sqrt(Real(3)) + Real(0.1)).Volume(),
		GeometryError);
}

TEST_CASE("Parallelogram diagonals reconstruct vertex lengths and shoelace area", "[geometry-formulas][Parallelogram]") {
	for (Real p : {Real(1), Real(3), Real(7)})
		for (Real q : {Real(1), Real(2), Real(5)})
			for (Real theta : {Constants::PI / 6, Constants::PI / 2, 5 * Constants::PI / 6}) {
				const auto shape = Parallelogram::FromDiagonals(p, q, theta);
				const Pnt2Cart first(p / 2, 0);
				const Pnt2Cart second(q * std::cos(theta) / 2, q * std::sin(theta) / 2);
				const Pnt2Cart third(-p / 2, 0);
				const Pnt2Cart fourth(-second.X(), -second.Y());
				const std::array<Pnt2Cart, 4> points{first, second, third, fourth};
				Real twiceArea = 0;
				for (int index = 0; index < 4; ++index) {
					const auto& next = points[(index + 1) % 4];
					twiceArea += points[index].X() * next.Y() - points[index].Y() * next.X();
				}
				CheckClose(shape.SideA(), first.Dist(fourth));
				CheckClose(shape.SideB(), first.Dist(second));
				CheckClose(shape.DiagonalP(), p);
				CheckClose(shape.DiagonalQ(), q);
				CheckClose(shape.Area(), std::abs(twiceArea) / 2);
			}
	const auto square = Parallelogram::FromDiagonals(std::sqrt(Real(2)), std::sqrt(Real(2)), Constants::PI / 2);
	CHECK(square.IsSquare());
	CheckClose(square.Area(), 1);
}

TEST_CASE("Equilateral square pyramids have equal base and apex edge lengths", "[geometry-formulas][pyramidequilateral]") {
	for (Real side : {Real(1), Real(3), Real(6)}) {
		for (const Vec3Cart center : {Vec3Cart(0, 0, 0), Vec3Cart(4, -5, 6)}) {
			const PyramidEquilateral3D pyramid(side, center);
			const Pnt3Cart apex(center.X(), center.Y(), center.Z() + pyramid.getHeight());
			const Pnt3Cart corner(center.X() + side / 2, center.Y() + side / 2, center.Z());
			CheckClose(apex.Dist(corner), side);
		}

	}
}

TEST_CASE("Pyramid faces match analytic bounds moments and translation", "[geometry-mesh][pyramid]") {
	for (const auto dimensions : {std::array<Real, 2>{2, 3}, std::array<Real, 2>{4, 1}}) {
		const Pyramid3D origin(dimensions[0], dimensions[1]);
		for (const Vec3Cart offset : {Vec3Cart(0, 0, 0), Vec3Cart(4, -5, 6)}) {
			const Pyramid3D body(dimensions[0], dimensions[1], offset);
			CheckClosedOutwardMesh(body, body.getCenter());
			CheckClose(body.BodyWithTriangleSurfaces::Volume(), body.Volume());
			CheckClose(body.BodyWithTriangleSurfaces::SurfaceArea(), body.SurfaceArea());
			CheckPoint(body.BodyWithTriangleSurfaces::getCenter(), body.getCenter());
			CheckPoint(body.BodyWithTriangleSurfaces::getBoundingBox().Min(), body.getBoundingBox().Min());
			CheckPoint(body.BodyWithTriangleSurfaces::getBoundingBox().Max(), body.getBoundingBox().Max());
			for (int index = 0; index < body.getSurfaceCount(); ++index) {
				const auto& face = body.getSurface(index);
				const auto& original = origin.getSurface(index);
				CheckPoint(face.Pnt1(), original.Pnt1() + offset);
				CheckPoint(face.Pnt2(), original.Pnt2() + offset);
				CheckPoint(face.Pnt3(), original.Pnt3() + offset);
			}
		}
	}
	const PyramidEquilateral3D equal(2);
	CheckClosedOutwardMesh(equal, equal.getCenter());
	for (int index = 2; index < equal.getSurfaceCount(); ++index) {
		const auto& face = equal.getSurface(index);
		CheckClose(face.A(), 2);
		CheckClose(face.B(), 2);
		CheckClose(face.C(), 2);
	}
}

TEST_CASE("Sphere triangles are outward closed and converge to analytic geometry", "[geometry-mesh][sphere]") {
	Real previousVolume = 0;
	Real previousArea = 0;
	for (int latitude : {2, 4, 8}) {
		const int longitude = 2 * latitude;
		const Pnt3Cart center(3, -4, 5);
		const Sphere3D body(2, center, latitude, longitude);
		CHECK(body.getSurfaceCount() == 2 * (latitude - 1) * longitude);
		CheckClosedOutwardMesh(body, center);
		const Real volume = body.BodyWithTriangleSurfaces::Volume();
		const Real area = body.BodyWithTriangleSurfaces::SurfaceArea();
		CHECK(volume > previousVolume);
		CHECK(volume < body.Volume());
		CHECK(area > previousArea);
		CHECK(area < body.SurfaceArea());
		CheckPoint(body.BodyWithTriangleSurfaces::getCenter(), center);
		for (int index = 0; index < body.getSurfaceCount(); ++index) {
			const auto& face = body.getSurface(index);
			for (const auto point : {face.Pnt1(), face.Pnt2(), face.Pnt3()})
				CheckClose(center.Dist(point), 2);
		}
		previousVolume = volume;
		previousArea = area;
	}
}

TEST_CASE("Coplanar quad maps honor every corner and integrated area", "[geometry-mesh][RectSurface3D]") {
	for (const Pnt3Cart third : {Pnt3Cart(4, 3, 0), Pnt3Cart(2, 3, 0), Pnt3Cart(5, 3, 0)}) {
		const RectSurface3D quad(Pnt3Cart(0, 0, 0), Pnt3Cart(4, 0, 0), third, Pnt3Cart(1, 3, 0));
		CheckPoint(SurfacePoint(quad, quad.getMinU(), quad.getMinW()), quad._pnt1);
		CheckPoint(SurfacePoint(quad, quad.getMaxU(), quad.getMinW()), quad._pnt2);
		CheckPoint(SurfacePoint(quad, quad.getMaxU(), quad.getMaxW()), quad._pnt3);
		CheckPoint(SurfacePoint(quad, quad.getMinU(), quad.getMaxW()), quad._pnt4);
		CheckPoint(SurfacePoint(quad, 0, 0), quad.getCenter());
		const Vec3Cart acrossU(SurfacePoint(quad, quad.getMinU(), 0), SurfacePoint(quad, quad.getMaxU(), 0));
		const Vec3Cart acrossW(SurfacePoint(quad, 0, quad.getMinW()), SurfacePoint(quad, 0, quad.getMaxW()));
		CheckClose(quad.getArea(), VectorProduct(acrossU, acrossW).NormL2());
		CheckClose(quad.getArea(), 3 * (4 + third.X() - 1) / 2);
	}
}

TEST_CASE("Torus patches agree with their stored corners and triangulation", "[geometry-mesh][torus]") {
	for (const Pnt3Cart center : {Pnt3Cart(0, 0, 0), Pnt3Cart(3, -4, 5)}) {
		const Torus3D body(3, 1, center, 8, 6);
		Real patchArea = 0;
		for (int index = 0; index < body.getSurfaceCount(); ++index) {
			const auto& patch = body.getSurface(index);
			CheckPoint(SurfacePoint(patch, patch.getMinU(), patch.getMinW()), patch._pnt1);
			CheckPoint(SurfacePoint(patch, patch.getMaxU(), patch.getMinW()), patch._pnt2);
			CheckPoint(SurfacePoint(patch, patch.getMaxU(), patch.getMaxW()), patch._pnt3);
			CheckPoint(SurfacePoint(patch, patch.getMinU(), patch.getMaxW()), patch._pnt4);
			const Vec3Cart diagonal(patch._pnt1, patch._pnt3);
			const Real triangleArea = (VectorProduct(Vec3Cart(patch._pnt1, patch._pnt2), diagonal).NormL2()
				+ VectorProduct(diagonal, Vec3Cart(patch._pnt1, patch._pnt4)).NormL2()) / 2;
			CheckClose(patch.getArea(), triangleArea);
			patchArea += patch.getArea();
		}
		CheckClose(patchArea, body.BodyWithRectSurfaces::SurfaceArea());
	}
}

TEST_CASE("Tetrahedron range fixtures stress sixth powers without losing the volume", "[geometry-formulas][Tetrahedron]") {
	const int exponent = std::min(100, (std::numeric_limits<Real>::max_exponent - 16) / 3);
	for (int sign : {-1, 1}) {
		const Real scale = std::ldexp(Real(1), sign * exponent);
		const Real cube = scale * scale * scale;
		const Real sixth = cube * cube;
		REQUIRE(std::isnormal(cube));
		if (sign < 0) {
			if constexpr (std::is_same_v<Real, float>) REQUIRE(sixth == 0);
		} else {
			if constexpr (std::is_same_v<Real, float>) REQUIRE(std::isinf(sixth));
		}
		const Real volume = Tetrahedron(scale, scale, scale, std::sqrt(Real(2)) * scale,
			std::sqrt(Real(2)) * scale, std::sqrt(Real(2)) * scale).Volume();
		REQUIRE(std::isnormal(volume));
		CheckClose(volume / cube, Real(1) / 6);
	}
}

TEST_CASE("Parallelogram short diagonals survive rounded equal side squares", "[geometry-formulas][Parallelogram]") {
	const Real longDiagonal = std::ldexp(Real(1), std::numeric_limits<Real>::digits / 2 + 1);
	REQUIRE(longDiagonal * longDiagonal + 1 == longDiagonal * longDiagonal);
	const auto shape = Parallelogram::FromDiagonals(1, longDiagonal, Constants::PI / 2);
	REQUIRE(shape.IsValid());
	CHECK(shape.IsRhombus());
	CheckClose(shape.DiagonalP(), 1);
	CheckClose(shape.DiagonalQ(), longDiagonal);
	CheckClose(shape.Area(), longDiagonal / 2);
}

TEST_CASE("Mesh centroids retain local coordinates under translation and face reordering", "[geometry-mesh][pyramid]") {
	struct Mesh : BodyWithTriangleSurfaces {
		explicit Mesh(std::vector<Triangle3D> faces) { _surfaces = std::move(faces); }
	};
	const Real translation = std::ldexp(Real(1), std::numeric_limits<Real>::digits - 5);
	for (const Vec3Cart offset : {Vec3Cart(0, 0, 0), Vec3Cart(translation, -translation, translation),
			 Vec3Cart(-translation, translation, -translation)}) {
		const Pyramid3D source(8, 4, offset);
		for (int first = 0; first < source.getSurfaceCount(); ++first) {
			for (bool reverse : {false, true}) {
				std::vector<Triangle3D> faces;
				for (int index = 0; index < source.getSurfaceCount(); ++index) {
					const auto& face = source.getSurface((first + index) % source.getSurfaceCount());
					faces.emplace_back(face.Pnt1(), reverse ? face.Pnt3() : face.Pnt2(),
						reverse ? face.Pnt2() : face.Pnt3());
				}
				const Mesh mesh(std::move(faces));
				const auto center = mesh.getCenter();
				CAPTURE(offset, first, reverse);
				CheckPoint(Pnt3Cart(center.X() - offset.X(), center.Y() - offset.Y(), center.Z() - offset.Z()),
					Pnt3Cart(0, 0, 1));
				CheckClose(mesh.Volume(), Real(256) / 3);
			}
		}
	}
}

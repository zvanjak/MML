#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../TestPrecision.h"
#include "../TestMatchers.h"

#include "base/Geometry2D.h"

using namespace MML;
using namespace MML::Testing;
using Catch::Matchers::WithinAbs;

TEST_CASE("Polygon2D - Construction and basic properties", "[Polygon2D][Construction]")
{
    // Empty polygon
    Polygon2D empty;
    REQUIRE(empty.NumVertices() == 0);
    
    // From vector
    std::vector<Point2Cartesian> vertices = {
        Point2Cartesian(0, 0),
        Point2Cartesian(1, 0),
        Point2Cartesian(1, 1),
        Point2Cartesian(0, 1)
    };
    Polygon2D square(vertices);
    REQUIRE(square.NumVertices() == 4);
    
    // From initializer list
    Polygon2D triangle{
        Point2Cartesian(0, 0),
        Point2Cartesian(1, 0),
        Point2Cartesian(REAL(0.5), 1)
    };
    REQUIRE(triangle.NumVertices() == 3);
}

TEST_CASE("Polygon2D - Vertex access", "[Polygon2D][Access]")
{
    Polygon2D poly{
        Point2Cartesian(0, 0),
        Point2Cartesian(2, 0),
        Point2Cartesian(2, 2),
        Point2Cartesian(0, 2)
    };
    
    // Operator[] access
    REQUIRE_THAT(poly[0].X(), WithinAbs(REAL(0.0), REAL(1e-10)));
    REQUIRE_THAT(poly[1].X(), WithinAbs(REAL(2.0), REAL(1e-10)));
    
    // Vertex() access
    REQUIRE_THAT(poly.Vertex(2).Y(), WithinAbs(REAL(2.0), REAL(1e-10)));
    
    // Out of range should throw
    REQUIRE_THROWS_AS(poly.Vertex(-1), std::out_of_range);
    REQUIRE_THROWS_AS(poly.Vertex(4), std::out_of_range);
}

TEST_CASE("Polygon2D - AddVertex and Clear", "[Polygon2D][Modification]")
{
    Polygon2D poly;
    REQUIRE(poly.NumVertices() == 0);
    
    poly.AddVertex(Point2Cartesian(0, 0));
    poly.AddVertex(Point2Cartesian(1, 0));
    poly.AddVertex(Point2Cartesian(1, 1));
    REQUIRE(poly.NumVertices() == 3);
    
    poly.Clear();
    REQUIRE(poly.NumVertices() == 0);
}

TEST_CASE("Polygon2D - Edge access", "[Polygon2D][Edges]")
{
    Polygon2D triangle{
        Point2Cartesian(0, 0),
        Point2Cartesian(3, 0),
        Point2Cartesian(0, 4)
    };
    
    // Edge 0: (0,0) -> (3,0)
    auto edge0 = triangle.Edge(0);
    REQUIRE_THAT(edge0.StartPoint().X(), WithinAbs(REAL(0.0), REAL(1e-10)));
    REQUIRE_THAT(edge0.EndPoint().X(), WithinAbs(REAL(3.0), REAL(1e-10)));
    REQUIRE_THAT(edge0.Length(), WithinAbs(REAL(3.0), REAL(1e-10)));
    
    // Edge 2: (0,4) -> (0,0) (wraps around)
    auto edge2 = triangle.Edge(2);
    REQUIRE_THAT(edge2.StartPoint().Y(), WithinAbs(REAL(4.0), REAL(1e-10)));
    REQUIRE_THAT(edge2.EndPoint().Y(), WithinAbs(REAL(0.0), REAL(1e-10)));
}

TEST_CASE("Polygon2D - Perimeter", "[Polygon2D][Perimeter]")
{
    // Unit square: perimeter = 4
    Polygon2D square{
        Point2Cartesian(0, 0),
        Point2Cartesian(1, 0),
        Point2Cartesian(1, 1),
        Point2Cartesian(0, 1)
    };
    REQUIRE_THAT(square.Perimeter(), WithinAbs(REAL(4.0), REAL(1e-10)));
    
    // 3-4-5 right triangle: perimeter = 12
    Polygon2D triangle{
        Point2Cartesian(0, 0),
        Point2Cartesian(3, 0),
        Point2Cartesian(0, 4)
    };
    REQUIRE_THAT(triangle.Perimeter(), WithinAbs(REAL(12.0), REAL(1e-10)));
}

TEST_CASE("Polygon2D - Area and SignedArea", "[Polygon2D][Area]")
{
    // Unit square CCW: area = 1, signed area = +1
    Polygon2D squareCCW{
        Point2Cartesian(0, 0),
        Point2Cartesian(1, 0),
        Point2Cartesian(1, 1),
        Point2Cartesian(0, 1)
    };
    REQUIRE_THAT(squareCCW.Area(), WithinAbs(REAL(1.0), REAL(1e-10)));
    REQUIRE_THAT(squareCCW.SignedArea(), WithinAbs(REAL(1.0), REAL(1e-10)));
    REQUIRE(squareCCW.IsCounterClockwise());
    
    // Unit square CW: area = 1, signed area = -1
    Polygon2D squareCW{
        Point2Cartesian(0, 0),
        Point2Cartesian(0, 1),
        Point2Cartesian(1, 1),
        Point2Cartesian(1, 0)
    };
    REQUIRE_THAT(squareCW.Area(), WithinAbs(REAL(1.0), REAL(1e-10)));
    REQUIRE_THAT(squareCW.SignedArea(), WithinAbs(-REAL(1.0), REAL(1e-10)));
    REQUIRE(squareCW.IsClockwise());
    
    // 3-4-5 right triangle: area = 6
    Polygon2D triangle{
        Point2Cartesian(0, 0),
        Point2Cartesian(3, 0),
        Point2Cartesian(0, 4)
    };
    REQUIRE_THAT(triangle.Area(), WithinAbs(REAL(6.0), REAL(1e-10)));
}

TEST_CASE("Polygon2D - Centroid", "[Polygon2D][Centroid]")
{
    // Unit square centered at (REAL(0.5), REAL(0.5))
    Polygon2D square{
        Point2Cartesian(0, 0),
        Point2Cartesian(1, 0),
        Point2Cartesian(1, 1),
        Point2Cartesian(0, 1)
    };
    auto centroid = square.Centroid();
    REQUIRE_THAT(centroid.X(), WithinAbs(REAL(0.5), REAL(1e-10)));
    REQUIRE_THAT(centroid.Y(), WithinAbs(REAL(0.5), REAL(1e-10)));
    
    // Triangle
    Polygon2D triangle{
        Point2Cartesian(0, 0),
        Point2Cartesian(3, 0),
        Point2Cartesian(0, 3)
    };
    auto triCentroid = triangle.Centroid();
    REQUIRE_THAT(triCentroid.X(), WithinAbs(REAL(1.0), REAL(1e-10)));
    REQUIRE_THAT(triCentroid.Y(), WithinAbs(REAL(1.0), REAL(1e-10)));
}

TEST_CASE("Polygon2D - BoundingBox", "[Polygon2D][BoundingBox]")
{
    Polygon2D poly{
        Point2Cartesian(-1, 2),
        Point2Cartesian(3, 5),
        Point2Cartesian(2, -1),
        Point2Cartesian(-2, 0)
    };
    
    auto bbox = poly.GetBoundingBox();
    REQUIRE_THAT(bbox.minX, WithinAbs(-REAL(2.0), REAL(1e-10)));
    REQUIRE_THAT(bbox.maxX, WithinAbs(REAL(3.0), REAL(1e-10)));
    REQUIRE_THAT(bbox.minY, WithinAbs(-REAL(1.0), REAL(1e-10)));
    REQUIRE_THAT(bbox.maxY, WithinAbs(REAL(5.0), REAL(1e-10)));
    REQUIRE_THAT(bbox.Width(), WithinAbs(REAL(5.0), REAL(1e-10)));
    REQUIRE_THAT(bbox.Height(), WithinAbs(REAL(6.0), REAL(1e-10)));
    
    auto center = bbox.Center();
    REQUIRE_THAT(center.X(), WithinAbs(REAL(0.5), REAL(1e-10)));
    REQUIRE_THAT(center.Y(), WithinAbs(REAL(2.0), REAL(1e-10)));
}

TEST_CASE("Polygon2D - Reverse orientation", "[Polygon2D][Reverse]")
{
    Polygon2D poly{
        Point2Cartesian(0, 0),
        Point2Cartesian(1, 0),
        Point2Cartesian(1, 1),
        Point2Cartesian(0, 1)
    };
    
    REQUIRE(poly.IsCounterClockwise());
    
    poly.Reverse();
    
    REQUIRE(poly.IsClockwise());
    REQUIRE_THAT(poly[0].X(), WithinAbs(REAL(0.0), REAL(1e-10)));
    REQUIRE_THAT(poly[0].Y(), WithinAbs(REAL(1.0), REAL(1e-10)));
}

TEST_CASE("Polygon2D - IsSimple (simple polygons)", "[Polygon2D][IsSimple]")
{
    // Simple square
    Polygon2D square{
        Point2Cartesian(0, 0),
        Point2Cartesian(1, 0),
        Point2Cartesian(1, 1),
        Point2Cartesian(0, 1)
    };
    REQUIRE(square.IsSimple());
    
    // Simple pentagon
    Polygon2D pentagon{
        Point2Cartesian(0, 0),
        Point2Cartesian(2, 0),
        Point2Cartesian(REAL(2.5), REAL(1.5)),
        Point2Cartesian(1, 2),
        Point2Cartesian(-REAL(0.5), REAL(1.5))
    };
    REQUIRE(pentagon.IsSimple());
}

TEST_CASE("Polygon2D - IsSimple (self-intersecting)", "[Polygon2D][IsSimple]")
{
    // Bowtie / figure-eight (self-intersecting)
    Polygon2D bowtie{
        Point2Cartesian(0, 0),
        Point2Cartesian(2, 2),
        Point2Cartesian(2, 0),
        Point2Cartesian(0, 2)
    };
    REQUIRE_FALSE(bowtie.IsSimple());
}

TEST_CASE("Polygon2D - IsConvex", "[Polygon2D][IsConvex]")
{
    // Convex square
    Polygon2D square{
        Point2Cartesian(0, 0),
        Point2Cartesian(1, 0),
        Point2Cartesian(1, 1),
        Point2Cartesian(0, 1)
    };
    REQUIRE(square.IsConvex());
    
    // Convex regular pentagon
    Polygon2D pentagon{
        Point2Cartesian(0, 0),
        Point2Cartesian(2, 0),
        Point2Cartesian(REAL(2.5), REAL(1.5)),
        Point2Cartesian(1, REAL(2.5)),
        Point2Cartesian(-REAL(0.5), REAL(1.5))
    };
    REQUIRE(pentagon.IsConvex());
    
    // Non-convex (L-shape)
    Polygon2D lShape{
        Point2Cartesian(0, 0),
        Point2Cartesian(2, 0),
        Point2Cartesian(2, 1),
        Point2Cartesian(1, 1),
        Point2Cartesian(1, 2),
        Point2Cartesian(0, 2)
    };
    REQUIRE_FALSE(lShape.IsConvex());
}

TEST_CASE("Polygon2D - Point containment (Contains)", "[Polygon2D][Contains]")
{
    Polygon2D square{
        Point2Cartesian(0, 0),
        Point2Cartesian(2, 0),
        Point2Cartesian(2, 2),
        Point2Cartesian(0, 2)
    };
    
    // Inside
    REQUIRE(square.Contains(Point2Cartesian(1, 1)));
    REQUIRE(square.Contains(Point2Cartesian(REAL(0.5), REAL(0.5))));
    REQUIRE(square.Contains(Point2Cartesian(REAL(1.9), REAL(1.9))));
    
    // Outside
    REQUIRE_FALSE(square.Contains(Point2Cartesian(-1, 1)));
    REQUIRE_FALSE(square.Contains(Point2Cartesian(3, 1)));
    REQUIRE_FALSE(square.Contains(Point2Cartesian(1, -1)));
    REQUIRE_FALSE(square.Contains(Point2Cartesian(1, 3)));
    
    // Corners and edges (boundary - may vary by implementation)
    // Most ray-casting implementations treat boundary as outside
}

TEST_CASE("Polygon2D - Point containment (non-convex)", "[Polygon2D][Contains]")
{
    // L-shaped polygon
    Polygon2D lShape{
        Point2Cartesian(0, 0),
        Point2Cartesian(2, 0),
        Point2Cartesian(2, 1),
        Point2Cartesian(1, 1),
        Point2Cartesian(1, 2),
        Point2Cartesian(0, 2)
    };
    
    // Inside the L
    REQUIRE(lShape.Contains(Point2Cartesian(REAL(0.5), REAL(0.5))));
    REQUIRE(lShape.Contains(Point2Cartesian(REAL(1.5), REAL(0.5))));
    REQUIRE(lShape.Contains(Point2Cartesian(REAL(0.5), REAL(1.5))));
    
    // Outside (in the "missing" corner)
    REQUIRE_FALSE(lShape.Contains(Point2Cartesian(REAL(1.5), REAL(1.5))));
    
    // Clearly outside
    REQUIRE_FALSE(lShape.Contains(Point2Cartesian(3, 3)));
    REQUIRE_FALSE(lShape.Contains(Point2Cartesian(-1, -1)));
}

TEST_CASE("Polygon2D - WindingNumber", "[Polygon2D][WindingNumber]")
{
    Polygon2D square{
        Point2Cartesian(0, 0),
        Point2Cartesian(2, 0),
        Point2Cartesian(2, 2),
        Point2Cartesian(0, 2)
    };
    
    // Inside: winding number = ±1 (depends on orientation)
    REQUIRE(square.WindingNumber(Point2Cartesian(1, 1)) != 0);
    
    // Outside: winding number = 0
    REQUIRE(square.WindingNumber(Point2Cartesian(-1, 1)) == 0);
    REQUIRE(square.WindingNumber(Point2Cartesian(3, 1)) == 0);
}

TEST_CASE("Polygon2D - Triangularization", "[Polygon2D][Triangularization]")
{
    // Triangle should return 1 triangle
    Polygon2D triangle{
        Point2Cartesian(0, 0),
        Point2Cartesian(1, 0),
        Point2Cartesian(REAL(0.5), 1)
    };
    auto triangles1 = triangle.Triangularization();
    REQUIRE(triangles1.size() == 1);
    
    // Square should return 2 triangles
    Polygon2D square{
        Point2Cartesian(0, 0),
        Point2Cartesian(1, 0),
        Point2Cartesian(1, 1),
        Point2Cartesian(0, 1)
    };
    auto triangles2 = square.Triangularization();
    REQUIRE(triangles2.size() == 2);
    
    // Pentagon should return 3 triangles
    Polygon2D pentagon{
        Point2Cartesian(0, 0),
        Point2Cartesian(2, 0),
        Point2Cartesian(REAL(2.5), REAL(1.5)),
        Point2Cartesian(1, 2),
        Point2Cartesian(-REAL(0.5), REAL(1.5))
    };
    auto triangles3 = pentagon.Triangularization();
    REQUIRE(triangles3.size() == 3);
}

TEST_CASE("Polygon2D - Legacy compatibility", "[Polygon2D][Legacy]")
{
    Polygon2D poly{
        Point2Cartesian(0, 0),
        Point2Cartesian(1, 0),
        Point2Cartesian(1, 1)
    };
    
    // Points() should return vertices
    auto points = poly.Points();
    REQUIRE(points.size() == 3);
    
    // IsInside() should work like Contains()
    REQUIRE(poly.IsInside(Point2Cartesian(REAL(0.5), REAL(0.25))));
    REQUIRE_FALSE(poly.IsInside(Point2Cartesian(2, 2)));
}

TEST_CASE("Polygon2D - Edge cases", "[Polygon2D][EdgeCases]")
{
    // Empty polygon
    Polygon2D empty;
    REQUIRE_THROWS_AS(empty.GetBoundingBox(), GeometryError);
    REQUIRE_THROWS_AS(empty.Centroid(), GeometryError);
    REQUIRE_THAT(empty.Area(), WithinAbs(REAL(0.0), REAL(1e-10)));
    REQUIRE_THAT(empty.Perimeter(), WithinAbs(REAL(0.0), REAL(1e-10)));
    
    // Single vertex
    Polygon2D singleVertex{Point2Cartesian(1, 2)};
    REQUIRE(singleVertex.NumVertices() == 1);
    auto centroid1 = singleVertex.Centroid();
    REQUIRE_THAT(centroid1.X(), WithinAbs(REAL(1.0), REAL(1e-10)));
    REQUIRE_THAT(centroid1.Y(), WithinAbs(REAL(2.0), REAL(1e-10)));
    
    // Two vertices (degenerate)
    Polygon2D twoVertices{
        Point2Cartesian(0, 0),
        Point2Cartesian(2, 0)
    };
    REQUIRE(twoVertices.NumVertices() == 2);
    auto centroid2 = twoVertices.Centroid();
    REQUIRE_THAT(centroid2.X(), WithinAbs(REAL(1.0), REAL(1e-10)));
    REQUIRE_THAT(centroid2.Y(), WithinAbs(REAL(0.0), REAL(1e-10)));
}

TEST_CASE("Polygon2D - Complex polygon containment", "[Polygon2D][Contains][Complex]")
{
    // Star-shaped polygon (non-convex but simple)
    Polygon2D star{
        Point2Cartesian(0, 2),
        Point2Cartesian(REAL(0.5), REAL(0.5)),
        Point2Cartesian(2, 0),
        Point2Cartesian(REAL(0.5), -REAL(0.5)),
        Point2Cartesian(0, -2),
        Point2Cartesian(-REAL(0.5), -REAL(0.5)),
        Point2Cartesian(-2, 0),
        Point2Cartesian(-REAL(0.5), REAL(0.5))
    };
    
    // Center should be inside
    REQUIRE(star.Contains(Point2Cartesian(0, 0)));
    
    // Points in the "arms" should be inside
    REQUIRE(star.Contains(Point2Cartesian(0, 1)));
    REQUIRE(star.Contains(Point2Cartesian(1, 0)));
    
    // Points far outside
    REQUIRE_FALSE(star.Contains(Point2Cartesian(0, 3)));
    REQUIRE_FALSE(star.Contains(Point2Cartesian(3, 0)));
}

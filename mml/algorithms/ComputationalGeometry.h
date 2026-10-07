///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        ComputationalGeometry.h                                             ///
///  Description: Aggregate header for computational geometry algorithms              ///
///               Includes all modular CompGeometry headers                           ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

/// @file ComputationalGeometry.h
/// @brief Aggregate header for computational geometry algorithms.
/// This header includes robust predicates, convex hulls, intersections, circle
/// geometry, triangulation, polygon operations, Voronoi diagrams, and selected
/// convenience aliases in the MML namespace.
///
/// Main types and algorithms included here:
/// - RobustPredicates - robust orientation and in-circle predicates
/// - ConvexHull2D, ConvexHull3D, and ConvexHull3DComputer
/// - Intersections - 2D/3D line, segment, and ray-triangle intersection helpers
/// - Circles - circle construction, containment, and closest-pair helpers
/// - Triangulation and DelaunayTriangulation - polygon and Delaunay triangulation
/// - PolygonOps - clipping, approximate boolean operations, and polygon measurements
/// - Voronoi, VoronoiDiagram, and VoronoiEdge - Voronoi diagram construction and storage
///
/// For more focused includes, use the individual headers under CompGeometry/.

#ifndef MML_COMPUTATIONAL_GEOMETRY_H
#define MML_COMPUTATIONAL_GEOMETRY_H

#include <mml/algorithms/CompGeometry/CompGeometryBase.h>
#include <mml/algorithms/CompGeometry/RobustPredicates.h>
#include <mml/algorithms/CompGeometry/ConvexHull.h>
#include <mml/algorithms/CompGeometry/Intersections.h>
#include <mml/algorithms/CompGeometry/Circles.h>
#include <mml/algorithms/CompGeometry/Triangulation.h>
#include <mml/algorithms/CompGeometry/PolygonOps.h>
#include <mml/algorithms/CompGeometry/VoronoiDiagram.h>
#include <mml/algorithms/CompGeometry/ConvexHull3D.h>

namespace MML {

using DelaunayTriangulation = CompGeometry::DelaunayTriangulation;
using VoronoiDiagram = CompGeometry::VoronoiDiagram;
using VoronoiEdge = CompGeometry::VoronoiEdge;
using ConvexHull3D = CompGeometry::ConvexHull3D;

} // namespace MML

#endif // MML_COMPUTATIONAL_GEOMETRY_H
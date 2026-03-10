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
#ifndef MML_COMPUTATIONAL_GEOMETRY_H
#define MML_COMPUTATIONAL_GEOMETRY_H

// Include all computational geometry headers
#include "mml/algorithms/CompGeometry/CompGeometryBase.h"
#include "mml/algorithms/CompGeometry/ConvexHull.h"
#include "mml/algorithms/CompGeometry/Intersections.h"
#include "mml/algorithms/CompGeometry/Circles.h"
#include "mml/algorithms/CompGeometry/Triangulation.h"
#include "mml/algorithms/CompGeometry/PolygonOps.h"
#include "mml/algorithms/CompGeometry/VoronoiDiagram.h"
#include "mml/algorithms/CompGeometry/ConvexHull3D.h"

namespace MML {

// ============================================================================
// Type Aliases for convenience
// ============================================================================
// These bring key types from CompGeometry namespace into MML namespace

using DelaunayTriangulation = CompGeometry::DelaunayTriangulation;
using VoronoiDiagram = CompGeometry::VoronoiDiagram;
using VoronoiEdge = CompGeometry::VoronoiEdge;
using ConvexHull3D = CompGeometry::ConvexHull3D;

} // namespace MML

#endif // MML_COMPUTATIONAL_GEOMETRY_H
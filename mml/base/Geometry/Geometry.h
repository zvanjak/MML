///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Geometry.h                                                          ///
///  Description: Aggregate header for core geometry components                      ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

/// @file Geometry.h
/// @brief Aggregate header for core geometry components.
/// This header includes geometry point types, coordinate representations, and
/// foundational two- and three-dimensional shape types.
///
/// Main types included here:
/// - Point2Cartesian, Point2Polar, Point3Cartesian, Point3Spherical, and Point3Cylindrical
/// - Triangle, Circle, Ellipse, CircularSector, CircularSegment, and Annulus
/// - RegularPolygon, Rectangle, Parallelogram, Rhombus, and Trapezoid
/// - SphereGeom, CylinderGeom, ConeGeom, Frustum, Tetrahedron, Spheroid, and TorusGeom
///
/// For more focused includes, use the individual headers under GeometryBase/.

#if !defined MML_GEOMETRY_H
#define MML_GEOMETRY_H

#include <mml/base/Geometry/GeometryBase/GeometryPoints.h>
#include <mml/base/Geometry/GeometryBase/Geometry2DShapes.h>
#include <mml/base/Geometry/GeometryBase/Geometry3DShapes.h>

#endif // MML_GEOMETRY_H

///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Geometry3D.h                                                        ///
///  Description: Aggregate header for 3D geometry classes                            ///
///               Includes lines, planes, triangles, and parametric surfaces          ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

/// @file Geometry3D.h
/// @brief Aggregate header for three-dimensional geometry classes.
/// This header includes lines, finite line segments, planes, triangles, and
/// parametric surface primitives for spatial geometry.
///
/// Main types included here:
/// - LineIntersectionType3D and LineIntersection3D - line intersection classification
/// - Line3D - infinite line with point projection, distance, and intersection operations
/// - SegmentLine3D - finite line segment
/// - Plane3D - plane representation with point, line, and plane operations
/// - Triangle3D - 3D triangle with geometric predicates and measurements
/// - TriangleSurface3D and RectSurface3D - parametric surface primitives
///
/// For more focused includes, use the individual headers under Geometry3D/.

#if !defined MML_GEOMETRY_3D_H
#define MML_GEOMETRY_3D_H

#include <mml/base/Geometry/Geometry3D/Geometry3DLines.h>
#include <mml/base/Geometry/Geometry3D/Geometry3DPlane.h>
#include <mml/base/Geometry/Geometry3D/Geometry3DSurfaces.h>

#endif // MML_GEOMETRY_3D_H

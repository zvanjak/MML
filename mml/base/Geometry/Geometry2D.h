///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Geometry2D.h                                                        ///
///  Description: Aggregate header for 2D geometry classes                           ///
///               Includes Line2D, SegmentLine2D, Triangle2D, Polygon2D,             ///
///               Circle2D, and Box2D                                                ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                        ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

/// @file Geometry2D.h
/// @brief Aggregate header for two-dimensional geometry classes.
/// This header includes line, segment, triangle, polygon, circle, and box types
/// for planar geometry.
///
/// Main types included here:
/// - LineIntersectionType2D and LineIntersection2D - line intersection classification and result data
/// - SegmentIntersectionType and SegmentIntersection - segment intersection classification and result data
/// - Line2D and SegmentLine2D - infinite and finite planar line geometry
/// - Triangle2D - triangle measurements, containment, and geometric centers
/// - Polygon2D - polygon geometry and computational-geometry operations
/// - Circle2D and Box2D - circle and axis-aligned box primitives
///
/// For more focused includes, use the individual headers under Geometry2D/.

#if !defined MML_GEOMETRY_2D_H
#define MML_GEOMETRY_2D_H

#include <mml/base/Geometry/Geometry2D/Geometry2DLines.h>
#include <mml/base/Geometry/Geometry2D/Geometry2DTriangle.h>
#include <mml/base/Geometry/Geometry2D/Geometry2DPolygon.h>
#include <mml/base/Geometry/Geometry2D/Geometry2DCircleBox.h>

#endif // MML_GEOMETRY_2D_H

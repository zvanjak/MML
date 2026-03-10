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
#if !defined MML_GEOMETRY_3D_H
#define MML_GEOMETRY_3D_H

// ============================================================================
// 3D Geometry Primitives (aggregate header)
// ============================================================================
//
// This header provides all 3D geometry classes:
//
// Lines (Geometry3DLines.h):
//   - LineIntersectionType3D  - Enum for line intersection classification
//   - LineIntersection3D      - Struct for detailed intersection results
//   - Line3D                  - Infinite line with distance/intersection ops
//   - SegmentLine3D          - Finite line segment
//
// Planes (Geometry3DPlane.h):
//   - Plane3D                - Plane with point/line/plane operations
//
// Surfaces (Geometry3DSurfaces.h):
//   - Triangle3D             - 3D triangle with area/containment tests
//   - TriangleSurface3D      - Parametric triangular surface
//   - RectSurface3D          - Parametric rectangular surface
//
// ============================================================================

#include "mml/base/Geometry/Geometry3DCore/Geometry3DLines.h"
#include "mml/base/Geometry/Geometry3DCore/Geometry3DPlane.h"
#include "mml/base/Geometry/Geometry3DCore/Geometry3DSurfaces.h"

#endif // MML_GEOMETRY_3D_H

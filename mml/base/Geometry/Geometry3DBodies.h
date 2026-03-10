///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Geometry3DBodies.h                                                  ///
///  Description: 3D solid bodies - aggregate header                                  ///
///               Includes bounding volumes, body interfaces, and solid primitives   ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                        ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

/// @file Geometry3DBodies.h
/// @brief Aggregate header for 3D solid body representations.
/// 
/// This header includes all 3D body modules:
/// - **Bounding volumes**: BoundingSphere3D, Box3D (AABB)
/// - **Body interfaces**: IBody, ISolidBodyWithBoundary, mesh body base classes
/// - **Simple solids**: Cube3D, CubeWithTriangles3D
/// - **Parametric shapes**: Torus3D, Cylinder3D
/// - **Other primitives**: Sphere3D, Pyramid3D, PyramidEquilateral3D
/// 
/// @see Geometry3DBodiesCore/ for individual module headers

#if !defined MML_GEOMETRY_3D_BODIES_H
#define MML_GEOMETRY_3D_BODIES_H

// Bounding volumes: BoundingSphere3D, Box3D (AABB)
#include "mml/base/Geometry/Geometry3DBodiesCore/Geometry3DBounding.h"

// Body interfaces and base classes: IBody, ISolidBodyWithBoundary,
// BodyWithTriangleSurfaces, BodyWithRectSurfaces, ComposedSolidSurfaces3D
#include "mml/base/Geometry/Geometry3DBodiesCore/Geometry3DBodyBase.h"

// Cube classes: Cube3D, CubeWithTriangles3D
#include "mml/base/Geometry/Geometry3DBodiesCore/Geometry3DCubes.h"

// Torus and Cylinder: Torus3D, Cylinder3D
#include "mml/base/Geometry/Geometry3DBodiesCore/Geometry3DTorusCylinder.h"

// Sphere and Pyramid: Sphere3D, Pyramid3D, PyramidEquilateral3D
#include "mml/base/Geometry/Geometry3DBodiesCore/Geometry3DSpherePyramid.h"

#endif

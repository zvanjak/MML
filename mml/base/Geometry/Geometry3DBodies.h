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
/// - Bounding volumes
/// - Body interfaces and mesh-backed body base classes
/// - Cube, torus, cylinder, sphere, and pyramid primitives
///
/// Main types included here:
/// - BoundingSphere3D and Box3D - bounding volume primitives and factories
/// - IBody and ISolidBodyWithBoundary - common body interfaces
/// - SolidBodyWithBoundary and SolidBodyWithBoundaryConstDensity - function-bounded solids
/// - BodyWithTriangleSurfaces, BodyWithRectSurfaces, and ComposedSolidSurfaces3D
/// - Cube3D, CubeWithTriangles3D, Torus3D, Cylinder3D, Sphere3D, Pyramid3D, and PyramidEquilateral3D
///
/// For more focused includes, use the individual headers under Geometry3DBodies/.

#if !defined MML_GEOMETRY_3D_BODIES_H
#define MML_GEOMETRY_3D_BODIES_H

#include <mml/base/Geometry/Geometry3DBodies/Geometry3DBounding.h>
#include <mml/base/Geometry/Geometry3DBodies/Geometry3DBodyBase.h>
#include <mml/base/Geometry/Geometry3DBodies/Geometry3DCubes.h>
#include <mml/base/Geometry/Geometry3DBodies/Geometry3DTorusCylinder.h>
#include <mml/base/Geometry/Geometry3DBodies/Geometry3DSpherePyramid.h>

#endif

///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Geometry.h                                                          ///
///  Description: AGGREGATE HEADER - includes all core geometry components            ///
///               Points (2D/3D), Coordinate systems, Pure 2D/3D shapes               ///
///                                                                                   ///
///  This file provides backward compatibility by including all split headers.        ///
///  For new code, prefer including specific headers:                                 ///
///    - GeometryCore/GeometryPoints.h   - Point classes (Cartesian, Polar, etc.)     ///
///    - GeometryCore/Geometry2DShapes.h - 2D shapes (Triangle, Circle, etc.)         ///
///    - GeometryCore/Geometry3DShapes.h - 3D shapes (Sphere, Cylinder, etc.)         ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_GEOMETRY_H
#define MML_GEOMETRY_H

// Core geometry components - split for maintainability
#include "mml/base/Geometry/GeometryCore/GeometryPoints.h"
#include "mml/base/Geometry/GeometryCore/Geometry2DShapes.h"
#include "mml/base/Geometry/GeometryCore/Geometry3DShapes.h"

#endif // MML_GEOMETRY_H

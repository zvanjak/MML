///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DifferentialGeometry_base.h                                         ///
///  Description: Aggregate header for pointwise differential-geometry primitives     ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

/// @file DifferentialGeometry_base.h
/// @brief Aggregate header for pointwise differential-geometry primitives.
/// This header includes frame tags, points, differential forms,
/// Hodge operations, and typed field wrappers.
///
/// Main types included here:
/// - Cartesian2, Cartesian3, Polar2, Cylindrical3, Spherical3, and Parameter2 frame tags
/// - Point - typed point coordinates bound to a scalar type, dimension, and frame tag
/// - DifferentialForm and DifferentialFormComponent - exterior-form storage and operations
/// - Orientation and HodgeStar - oriented Hodge dual operations
/// - IScalarField, IVectorField, IFormField, ScalarFunctionFieldAdapter, and VectorFunctionFieldAdapter
///
/// Constructing metrics and algorithmic geometry operations require DifferentialGeometry_core.h.

#if !defined MML_BASE_DIFFERENTIAL_GEOMETRY_H
#define MML_BASE_DIFFERENTIAL_GEOMETRY_H

// Constructing typed metrics requires core/DifferentialGeometry/Metric.h.
#include <mml/base/DifferentialGeometry/FrameTags.h>
#include <mml/base/DifferentialGeometry/Point.h>
#include <mml/base/DifferentialGeometry/DifferentialForm.h>
#include <mml/base/DifferentialGeometry/Hodge.h>
#include <mml/base/DifferentialGeometry/TypedFields.h>

#endif // MML_BASE_DIFFERENTIAL_GEOMETRY_H
///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DifferentialGeometry_core.h                                         ///
///  Description: Aggregate header for core differential-geometry algorithms          ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

/// @file DifferentialGeometry_core.h
/// @brief Aggregate header for core differential-geometry algorithms.
/// This header includes charts, atlases, coordinate maps, field operations,
/// induced metrics, metric adapters, form integration, and frame helpers.
///
/// Main types and functions included here:
/// - Chart and FunctionChart - chart embeddings with push-forward and pull-back helpers
/// - Atlas, UnitSphereStereographicChart, and MakeUnitSphereStereographicAtlas
/// - CoordinateMap, IdentityCoordinateMap, and CoordTransfCoordinateMap
/// - ExteriorDerivative, exterior_derivative, Gradient, and GradientCart
/// - InducedMetric2D, GaussianCurvatureIntrinsic, and MetricAt
/// - IntegrateOneForm, IntegrateTwoForm, rectangle boundary/integral helpers, and Stokes checks
/// - FrenetFrame, DarbouxFrame, ComputeFrenetFrame, and ComputeDarbouxFrame
///
/// For pointwise geometry primitives, include DifferentialGeometry_base.h.

#if !defined MML_CORE_DIFFERENTIAL_GEOMETRY_H
#define MML_CORE_DIFFERENTIAL_GEOMETRY_H

#include <mml/core/DifferentialGeometry/Chart.h>
#include <mml/core/DifferentialGeometry/Atlas.h>
#include <mml/core/DifferentialGeometry/CoordinateMap.h>
#include <mml/core/DifferentialGeometry/FieldOperations.h>
#include <mml/core/DifferentialGeometry/InducedMetric.h>
#include <mml/core/DifferentialGeometry/MetricAdapters.h>
#include <mml/core/DifferentialGeometry/FormIntegration.h>
#include <mml/core/DifferentialGeometry/Frames.h>

#endif // MML_CORE_DIFFERENTIAL_GEOMETRY_H
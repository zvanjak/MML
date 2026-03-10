///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        InterpolatedFunction.h                                              ///
///  Description: Aggregate header for all interpolation classes                      ///
///               Includes 1D, 2D, and parametric curve interpolation                 ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

/// @file InterpolatedFunction.h
/// @brief Aggregate header for all interpolation classes.
/// This header includes all interpolation functionality:
/// - 1D interpolation (linear, polynomial, spline, rational)
/// - 2D grid interpolation (bilinear, bicubic spline)
/// - Parametric curve interpolation (linear, spline)
///
/// For more focused includes, use the individual headers:
/// - InterpolatedFunctions/InterpolatedRealFunction.h - 1D interpolation
/// - InterpolatedFunctions/Interpolation2DFunction.h - 2D grid interpolation
/// - InterpolatedFunctions/InterpolationParametricCurve.h - Parametric curves
///
/// @see InterpolatedRealFunction.h for detailed 1D interpolation documentation
/// @ingroup Interpolation

#if !defined MML_INTERPOLATEDFUNCTION_H
#define MML_INTERPOLATEDFUNCTION_H

// Include all interpolation headers
#include "base/InterpolatedFunctions/InterpolatedRealFunction.h"
#include "base/InterpolatedFunctions/Interpolation2DFunction.h"
#include "base/InterpolatedFunctions/InterpolationParametricCurve.h"

// Include Function.h for backward compatibility (RealFunction typedef, etc.)
#include "base/Function.h"

#endif // MML_INTERPOLATEDFUNCTION_H

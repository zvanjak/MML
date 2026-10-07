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
/// Main types included here:
/// - InterpolatedRealFunctionLinear - shared base and linear interpolation
/// - InterpolatedRealFunctionPolynomial, InterpolatedRealFunctionRational, and InterpolatedRealFunctionBarycentric
/// - InterpolatedFunctionSpline and monotone cubic spline support
/// - InterpolatedRealFunctionAdvanced - Hermite and Akima interpolation
/// - Interpolation2DFunction - two-dimensional grid interpolation
/// - InterpolationParametricCurve - parametric curve interpolation
///
/// @see InterpolatedRealFunctionLinear.h for the shared 1D interpolation base
/// @ingroup Interpolation

#if !defined MML_INTERPOLATEDFUNCTION_H
#define MML_INTERPOLATEDFUNCTION_H

#include <mml/base/InterpolatedFunctions/InterpolatedRealFunctionLinear.h>
#include <mml/base/InterpolatedFunctions/InterpolatedRealFunctionPolynomial.h>
#include <mml/base/InterpolatedFunctions/InterpolatedRealFunctionRational.h>
#include <mml/base/InterpolatedFunctions/InterpolatedRealFunctionBarycentric.h>
#include <mml/base/InterpolatedFunctions/InterpolatedFunctionSpline.h>
#include <mml/base/InterpolatedFunctions/InterpolatedRealFunctionAdvanced.h>
#include <mml/base/InterpolatedFunctions/Interpolation2DFunction.h>
#include <mml/base/InterpolatedFunctions/InterpolationParametricCurve.h>

#include <mml/base/Function.h>

#endif // MML_INTERPOLATEDFUNCTION_H

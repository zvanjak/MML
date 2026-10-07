///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        BaseUtils.h                                                         ///
///  Description: Aggregate header for utility functions                              ///
///               Includes all focused utility headers                                ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

/// @file BaseUtils.h
/// @brief Aggregate header for common utility functions.
/// This header includes symbolic utilities, angle and coordinate helpers,
/// comparisons, vector operations, matrix operations, and mixed Real/Complex helpers.
///
/// Main utility groups included here:
/// - LeviCivita and KroneckerDelta - symbolic tensor helpers
/// - DegToRad, RadToDeg, angle normalization, and explicit degree/minute/second conversion
/// - Cartesian, polar, spherical, and cylindrical coordinate conversions
/// - AreEqual and AreEqualAbs overloads for Complex, Vector<Real>, and Vector<Complex>
/// - ScalarProduct, VectorsAngle, projections, and OuterProduct
/// - Matrix construction, comparison, transformation, orthogonalization, and matrix functions
/// - Mixed Real/Complex vector and matrix arithmetic helpers
///
/// For more focused includes, use the individual headers under BaseUtils/.

#ifndef MML_BASEUTILS_H
#define MML_BASEUTILS_H

#include <vector>

#include <mml/MMLBase.h>
#include <mml/base/Vector/Vector.h>
#include <mml/base/Vector/VectorN.h>
#include <mml/base/Matrix/Matrix.h>
#include <mml/base/Matrix/MatrixNM.h>
#include <mml/base/BaseUtils/SymbolUtils.h>
#include <mml/base/BaseUtils/AngleCoordUtils.h>
#include <mml/base/BaseUtils/ComparisonUtils.h>
#include <mml/base/BaseUtils/VectorOps.h>
#include <mml/base/BaseUtils/MatrixOps.h>
#include <mml/base/BaseUtils/MixedTypeOps.h>

#endif // MML_BASEUTILS_H

///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        BaseUtils.h                                                         ///
///  Description: Aggregate header for utility functions                              ///
///               Includes all split utility headers for backward compatibility       ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#ifndef MML_BASEUTILS_H
#define MML_BASEUTILS_H

// Standard headers
#include <vector>

#include "MMLBase.h"

// Include MML base types needed by utilities
#include "base/Vector/Vector.h"
#include "base/Vector/VectorN.h"
#include "base/Matrix/Matrix.h"
#include "base/Matrix/MatrixNM.h"

// ============================================================================
// Split Utility Headers
// ============================================================================
// The functionality previously in this monolithic file has been split into
// focused utility modules for better organization and faster compile times.
// ============================================================================

#include "mml/base/BaseUtils/SymbolUtils.h"       // LeviCivita, KroneckerDelta
#include "mml/base/BaseUtils/AngleCoordUtils.h"   // Angle conversions, coordinate transforms
#include "mml/base/BaseUtils/ComparisonUtils.h"   // AreEqual for Complex, Vector
#include "mml/base/BaseUtils/VectorOps.h"         // ScalarProduct, VectorsAngle, projections
#include "mml/base/BaseUtils/MatrixOps.h"         // Matrix functions, properties, creation
#include "mml/base/BaseUtils/MixedTypeOps.h"      // Complex+Real mixed operations

// ============================================================================
// Backward Compatibility
// ============================================================================
// All functions are now available in MML::Utils namespace through the includes.
// Existing code using MML::Utils::FunctionName() will continue to work.
// ============================================================================

#endif // MML_BASEUTILS_H

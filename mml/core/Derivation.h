///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Derivation.h                                                        ///
///  Description: Numerical differentiation umbrella header                           ///
///               Includes all derivative computation modules                         ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

/// @file Derivation.h
/// @brief Aggregate header for numerical differentiation.
/// This header includes finite-difference derivatives for real, scalar, vector,
/// complex, parametric curve, parametric surface, and tensor-field functions,
/// plus complex-step differentiation and Jacobian helpers.
///
/// Main functions included here:
/// - NDer1, NDer2, NDer4, NDer6, NDer8, and left/right one-sided real derivatives
/// - NDer1Partial, NDer2Partial, gradients, Hessians, directional derivatives, and Laplacians for scalar functions
/// - Vector-function partial derivatives, full Jacobian matrices, divergence, curl, and vector Laplacian helpers
/// - Parametric-curve and parametric-surface first/second derivative helpers
/// - Tensor-field component derivative helpers
/// - ComplexStep, NDer1Complex, NDer2Complex, and NDer4Complex
/// - CalcJacobian and CalcJacobian2 for static and dynamic vector functions
///
/// For more focused includes, use the individual headers under Derivation/.

#if !defined MML_DERIVATION_H
#define MML_DERIVATION_H

#include <mml/MMLBase.h>
#include <mml/core/Derivation/DerivationBase.h>
#include <mml/core/Derivation/DerivationRealFunction.h>
#include <mml/core/Derivation/DerivationScalarFunction.h>
#include <mml/core/Derivation/DerivationVectorFunction.h>
#include <mml/core/Derivation/DerivationParametricCurve.h>
#include <mml/core/Derivation/DerivationParametricSurface.h>
#include <mml/core/Derivation/DerivationTensorField.h>
#include <mml/core/Derivation/DerivationComplexStep.h>
#include <mml/core/Derivation/DerivationComplex.h>
#include <mml/core/Derivation/Jacobians.h>
#include <mml/core/Derivation/Hessians.h>

#endif // MML_DERIVATION_H
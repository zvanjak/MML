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
#if !defined MML_DERIVATION_H
#define MML_DERIVATION_H

#include "MMLBase.h"

/// @file DerivationBase.h
/// @brief Core numerical differentiation algorithms and utilities
#include "core/Derivation/DerivationBase.h"

/// @file DerivationRealFunction.h
/// @brief Derivatives of real-valued functions f: ℝ → ℝ
#include "core/Derivation/DerivationRealFunction.h"
/// @file DerivationScalarFunction.h
/// @brief Derivatives of scalar functions f: ℝⁿ → ℝ (gradients)
#include "core/Derivation/DerivationScalarFunction.h"
/// @file DerivationVectorFunction.h
/// @brief Derivatives of vector functions f: ℝⁿ → ℝᵐ
#include "core/Derivation/DerivationVectorFunction.h"
/// @file DerivationParametricCurve.h
/// @brief Derivatives of parametric curves (tangent vectors)
#include "core/Derivation/DerivationParametricCurve.h"
/// @file DerivationParametricSurface.h
/// @brief Derivatives of parametric surfaces (tangent vectors, normal vectors)
#include "core/Derivation/DerivationParametricSurface.h"
/// @file DerivationTensorField.h
/// @brief Derivatives of tensor fields (covariant derivatives)
#include "core/Derivation/DerivationTensorField.h"

/// @file Jacobians.h
/// @brief Jacobian matrix computations for multivariable functions
#include "core/Derivation/Jacobians.h"

#endif // MML_DERIVATION_H
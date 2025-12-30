///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Derivation.h                                                        ///
///  Description: Numerical differentiation umbrella header                           ///
///               Includes all derivative computation modules                         ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DERIVATION_H
#define MML_DERIVATION_H

#include "MMLBase.h"

#include "core/Derivation/DerivationBase.h"

#include "core/Derivation/DerivationRealFunction.h"
#include "core/Derivation/DerivationScalarFunction.h"
#include "core/Derivation/DerivationVectorFunction.h"
#include "core/Derivation/DerivationParametricCurve.h"
#include "core/Derivation/DerivationParametricSurface.h"
#include "core/Derivation/DerivationTensorField.h"

#include "core/Derivation/Jacobians.h"

#endif // MML_DERIVATION_H
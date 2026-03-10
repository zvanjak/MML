///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        RootFinding.h                                                       ///
///  Description: Root finding algorithms - Umbrella header                           ///
///               Includes: Bisection, Newton, Secant, Brent, Ridders methods         ///
///               Bracketing, multi-root finding, and convergence control             ///
///                                                                                   ///
///  REFACTORED: Split into modular headers for better organization:                  ///
///    - RootFindingBase.h:      Config and Result types                              ///
///    - RootFindingBracketing.h: Bracketing utilities + Bisection + False Position   ///
///    - RootFindingMethods.h:   Newton, Secant, Ridders, Brent methods               ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_ROOTFINDING_H
#define MML_ROOTFINDING_H

// Modular root-finding headers (in root_finding subfolder)
#include "mml/algorithms/RootFinding/RootFindingBase.h"
#include "mml/algorithms/RootFinding/RootFindingBracketing.h"
#include "mml/algorithms/RootFinding/RootFindingMethods.h"
#include "mml/algorithms/RootFinding/RootFindingPolynoms.h"

#endif // MML_ROOTFINDING_H

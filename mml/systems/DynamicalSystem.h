///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DynamicalSystem.h                                                   ///
///  Description: Umbrella header for dynamical systems framework                     ///
///               Includes all dynamical system components                            ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                        ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DYNAMICAL_SYSTEM_H
#define MML_DYNAMICAL_SYSTEM_H

// Types and result structures
#include "systems/DynamicalSystemTypes.h"

// Base class for dynamical systems
#include "systems/DynamicalSystemBase.h"

// Analysis tools (fixed points, Lyapunov, bifurcations, phase space)
#include "systems/DynamicalSystemAnalyzers.h"

// Classic continuous systems (Lorenz, Rössler, Van der Pol, etc.)
#include "systems/ContinuousSystems.h"

// Discrete maps (Logistic, Hénon, Standard, Tent)
#include "systems/DiscreteMaps.h"

// Unified analyzer facade
#include "systems/DynamicalSystemAnalyzer.h"

#endif // MML_DYNAMICAL_SYSTEM_H

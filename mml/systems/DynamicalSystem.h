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
#include <mml/systems/DynamicalSystem/DynamicalSystemTypes.h>

// Base class for dynamical systems
#include <mml/systems/DynamicalSystem/DynamicalSystemBase.h>

// Analysis tools (fixed points, Lyapunov, bifurcations, phase space)
#include <mml/systems/DynamicalSystem/FixedPointAnalysis.h>
#include <mml/systems/DynamicalSystem/LyapunovAnalysis.h>
#include <mml/systems/DynamicalSystem/BifurcationAnalysis.h>
#include <mml/systems/DynamicalSystem/PhaseSpaceAnalysis.h>

// Classic continuous systems (Lorenz, Rössler, Van der Pol, etc.)
#include <mml/systems/ContinuousSystems.h>

// Discrete maps (Logistic, Hénon, Standard, Tent)
#include <mml/systems/DiscreteMaps.h>

// Unified analyzer facade
#include <mml/systems/DynamicalSystemAnalyzer.h>

#endif // MML_DYNAMICAL_SYSTEM_H

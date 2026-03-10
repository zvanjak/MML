///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Serializer.h                                                        ///
///  Description: Data serialization utilities for MML objects (umbrella include)     ///
///               Includes all serialization modules for functions, curves, etc.      ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
///
/// @file Serializer.h
/// @brief Umbrella include for all serialization utilities
/// 
/// This file provides backward compatibility by including all serialization modules.
/// For more granular includes, you can use the individual headers:
///
/// - tools/serializer/SerializerBase.h       - Base types (SerializeError, SerializeResult) and header writers
/// - tools/serializer/SerializerFunctions.h  - Real function serialization (SaveRealFunc, etc.)
/// - tools/serializer/SerializerCurves.h     - Parametric curve serialization
/// - tools/serializer/SerializerSurfaces.h   - Surface and scalar function serialization
/// - tools/serializer/SerializerVectors.h    - Vector field serialization
/// - tools/serializer/SerializerODE.h        - ODE solution serialization
/// - tools/serializer/SerializerSimulation.h - Particle simulation serialization
/// - tools/serializer/SerializerFieldLines.h - Field line serialization
///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_SERIALIZER_H
#define MML_SERIALIZER_H

// Include all serialization modules
#include "mml/tools/serializer/SerializerBase.h"
#include "mml/tools/serializer/SerializerFunctions.h"
#include "mml/tools/serializer/SerializerCurves.h"
#include "mml/tools/serializer/SerializerSurfaces.h"
#include "mml/tools/serializer/SerializerVectors.h"
#include "mml/tools/serializer/SerializerODE.h"
#include "mml/tools/serializer/SerializerSimulation.h"
#include "mml/tools/serializer/SerializerFieldLines.h"

#endif // MML_SERIALIZER_H

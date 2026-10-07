///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Serializer.h                                                        ///
///  Description: Presentation-data exporters for MML visualizers                     ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
///
/// @file Serializer.h
/// @brief Umbrella include for visualizer presentation-data exporters
/// 
/// Durable object save/load APIs are provided separately by tools/Persistence.h.
/// This aggregate header includes only visualizer export modules.
/// For more granular includes, you can use the individual headers:
///
/// - tools/serializer/SerializerBase.h       - Export format constants and header writers
/// - tools/serializer/SerializerFunctions.h  - Real function serialization (SaveRealFunc, etc.)
/// - tools/serializer/SerializerCurves.h     - Parametric curve serialization
/// - tools/serializer/SerializerSurfaces.h   - Surface and scalar function serialization
/// - tools/serializer/SerializerVectorFields.h - Vector field serialization
/// - tools/serializer/SerializerODE.h        - ODE solution serialization
/// - tools/serializer/SerializerSimulation.h - Particle simulation serialization
///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_SERIALIZER_H
#define MML_SERIALIZER_H

// Include all presentation-data export modules
#include <mml/tools/serializer/SerializerBase.h>
#include <mml/tools/serializer/SerializerFunctions.h>
#include <mml/tools/serializer/SerializerCurves.h>
#include <mml/tools/serializer/SerializerSurfaces.h>
#include <mml/tools/serializer/SerializerVectorFields.h>
#include <mml/tools/serializer/SerializerODE.h>
#include <mml/tools/serializer/SerializerSimulation.h>


#endif // MML_SERIALIZER_H

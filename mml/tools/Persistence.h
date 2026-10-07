///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Persistence.h                                                       ///
///  Description: Durable object persistence umbrella include                        ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_PERSISTENCE_H
#define MML_PERSISTENCE_H

#include <mml/tools/persistence/PersistenceBase.h>
#include <mml/tools/persistence/JSON.h>
#include <mml/tools/persistence/BinaryBase.h>
#include <mml/tools/persistence/VectorBinary.h>
#include <mml/tools/persistence/MatrixBinary.h>
#include <mml/tools/persistence/VectorJSON.h>
#include <mml/tools/persistence/MatrixJSON.h>
#include <mml/tools/persistence/FunctionJSON.h>
#include <mml/tools/persistence/MatrixIO.h>

#endif // MML_PERSISTENCE_H
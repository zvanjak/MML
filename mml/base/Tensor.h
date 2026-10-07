///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Tensor.h                                                            ///
///  Description: Public umbrella for concrete tensors of ranks 1 through 5           ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

/// @file Tensor.h
/// @brief Aggregate header for concrete tensor ranks.
/// This header includes the rank-1 tensor alias and concrete rank-2 through rank-5
/// tensor classes.
///
/// Main types included here:
/// - Rank1Tensor - rank-1 tensor alias for vector-like tensor quantities
/// - Tensor2, Tensor3, Tensor4, and Tensor5 - fixed-rank tensor containers
/// - Rank-specific indexing, arithmetic, and component access operations
///
/// For more focused includes, use the individual headers under Tensor/.

#if !defined MML_TENSOR_H
#define MML_TENSOR_H

#include <mml/base/Tensor/Rank1Tensor.h>
#include <mml/base/Tensor/Tensor2.h>
#include <mml/base/Tensor/Tensor3.h>
#include <mml/base/Tensor/Tensor4.h>
#include <mml/base/Tensor/Tensor5.h>

#endif

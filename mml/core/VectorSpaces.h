///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        VectorSpaces.h                                                      ///
///  Description: Finite-dimensional vector-space abstractions                        ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

/// @file VectorSpaces.h
/// @brief Aggregate header for finite-dimensional vector-space abstractions.
/// This header includes affine spaces, bases, dual spaces, inner-product spaces,
/// linear maps, and subspaces.
///
/// Main types included here:
/// - PointInSpace, AffineFrame, AffineSpace, and AffineMap
/// - Basis - ordered basis storage and coordinate conversion helpers
/// - VectorInSpace, CovectorInSpace, and DualSpace
/// - InnerProductSpace - finite-dimensional inner-product structure and Euclidean factory
/// - LinearMap - fixed-size linear maps with composition and identity helpers
/// - Subspace - null space, column space, row space, left null space, sums, and intersections
///
/// For more focused includes, use the individual headers under VectorSpaces/.

#if !defined MML_VECTOR_SPACES_H
#define MML_VECTOR_SPACES_H

#include <mml/core/VectorSpaces/AffineSpace.h>
#include <mml/core/VectorSpaces/Basis.h>
#include <mml/core/VectorSpaces/DualSpace.h>
#include <mml/core/VectorSpaces/InnerProductSpace.h>
#include <mml/core/VectorSpaces/LinearMap.h>
#include <mml/core/VectorSpaces/Subspace.h>

#endif // MML_VECTOR_SPACES_H
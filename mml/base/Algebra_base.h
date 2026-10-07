///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Algebra_base.h                                                      ///
///  Description: Aggregate header for foundational algebraic types                   ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

/// @file Algebra_base.h
/// @brief Aggregate header for algebraic types and algorithms.
/// This header includes finite groups, permutations, modular arithmetic,
/// finite fields, polynomial types, representations, Lie group primitives,
/// law checks, and algorithms over these structures.
///
/// Main types included here:
/// - AlgebraTraits and AlgebraEquality - common algebraic structure traits
/// - Permutation and PermutationParity - finite permutation operations
/// - FiniteGroup, CyclicGroup, DihedralGroup, and GroupAction
/// - ModInt and ModularRing - modular arithmetic over integer residue classes
/// - Polynomial - generic polynomial type over arbitrary coefficient domains
/// - FiniteFieldElement - finite field element representation
/// - Representation - finite-dimensional matrix representation of group elements
/// - SO2, SO3, SE2, and SE3 - Lie groups for rotations and rigid motions
/// - Group law checks, finite-group and group-action algorithms
/// - Modular, finite-field matrix, polynomial, and representation algorithms

#if !defined MML_ALGEBRA_H
#define MML_ALGEBRA_H

#include <mml/base/Algebra/AlgebraTraits.h>
#include <mml/base/Algebra/Permutation.h>
#include <mml/base/Algebra/FiniteGroup.h>
#include <mml/base/Algebra/CyclicGroup.h>
#include <mml/base/Algebra/DihedralGroup.h>
#include <mml/base/Algebra/GroupAction.h>
#include <mml/base/Algebra/ModInt.h>
#include <mml/base/Algebra/ModularStructures.h>
#include <mml/base/Algebra/Polynomial.h>
#include <mml/base/Algebra/FiniteField.h>
#include <mml/base/Algebra/Representation.h>
#include <mml/base/Algebra/LieGroups.h>

#include <mml/base/Algebra/Algorithms/AlgebraLawChecks.h>
#include <mml/base/Algebra/Algorithms/FiniteGroupAlgorithms.h>
#include <mml/base/Algebra/Algorithms/GroupActionAlgorithms.h>
#include <mml/base/Algebra/Algorithms/ModularArithmeticAlgorithms.h>
#include <mml/base/Algebra/Algorithms/FieldMatrixAlgorithms.h>
#include <mml/base/Algebra/Algorithms/PolynomialAlgorithms.h>
#include <mml/base/Algebra/Algorithms/RepresentationAlgorithms.h>

#endif // MML_ALGEBRA_H
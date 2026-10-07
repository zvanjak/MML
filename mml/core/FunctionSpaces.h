///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        FunctionSpaces.h                                                    ///
///  Description: Finite function-space and operator-discretization aggregate header  ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

/// @file FunctionSpaces.h
/// @brief Aggregate header for finite function spaces and discretized operators.
/// This header includes one-dimensional function spaces, trial spaces, expansions,
/// projections, interpolation, collocation spaces, boundary conditions, assembled
/// operators, matrix-free operators, and dense BVP helpers.
///
/// Main types included here:
/// - DiscretizationMethod, FunctionSpaceError, and FunctionSpaceInputError
/// - FunctionSpace1D, L2IntervalSpace, WeightedL2IntervalSpace, and TrialSpace1D
/// - FunctionExpansion1D and FunctionSpaceOperationResult1D
/// - OrthogonalBasisFunctionSpace1D, OrthogonalBasisTrialSpace1D, and ChebyshevCollocationSpace1D
/// - LinearDifferentialOperator1D, LinearOperator, and DirichletSecondDerivativeOperator1D
/// - BoundaryConditionKind, BoundaryCondition1D, and BoundaryConditions1D
/// - DenseBVPSolveResult1D and dense one-dimensional boundary-value solver helpers
///
/// For more focused includes, use the individual headers under FunctionSpaces/.

#if !defined MML_FUNCTION_SPACES_H
#define MML_FUNCTION_SPACES_H

#include <mml/core/FunctionSpaces/FunctionSpacesBase.h>
#include <mml/core/FunctionSpaces/FunctionSpace1D.h>
#include <mml/core/FunctionSpaces/TrialSpace1D.h>
#include <mml/core/FunctionSpaces/FunctionExpansion1D.h>
#include <mml/core/FunctionSpaces/FunctionSpaceResult.h>
#include <mml/core/FunctionSpaces/OrthogonalBasisTrialSpace1D.h>
#include <mml/core/FunctionSpaces/Projection.h>
#include <mml/core/FunctionSpaces/Interpolation.h>
#include <mml/core/FunctionSpaces/ChebyshevCollocationSpace1D.h>
#include <mml/core/FunctionSpaces/LinearDifferentialOperator1D.h>
#include <mml/core/FunctionSpaces/BoundaryCondition1D.h>
#include <mml/core/FunctionSpaces/OperatorAssembly.h>
#include <mml/core/FunctionSpaces/DenseBVPSolver1D.h>
#include <mml/core/FunctionSpaces/BVPDiagnostics1D.h>
#include <mml/core/FunctionSpaces/LinearOperator.h>
#include <mml/core/FunctionSpaces/MatrixFreeOperators1D.h>

#endif // MML_FUNCTION_SPACES_H

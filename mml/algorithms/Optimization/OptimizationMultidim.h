///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Optimization/OptimizationMultidim.h                                 ///
///  Description: Multidimensional optimization aggregate header                      ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_OPTIMIZATION_MULTIDIM_H
#define MML_OPTIMIZATION_MULTIDIM_H

#include <mml/algorithms/Optimization/Multidim/MultidimTypes.h>
#include <mml/algorithms/Optimization/Multidim/NelderMead.h>
#include <mml/algorithms/Optimization/Multidim/LineSearch.h>
#include <mml/algorithms/Optimization/Multidim/Powell.h>
#include <mml/algorithms/Optimization/Multidim/QuasiNewton.h>
#include <mml/algorithms/Optimization/Multidim/MultidimSolvers.h>
#include <mml/algorithms/Optimization/Constraints/BoundConstraints.h>
#include <mml/algorithms/Optimization/Constraints/ProjectedGradient.h>
#include <mml/algorithms/Optimization/Constraints/BoxConstrainedNelderMead.h>
#include <mml/algorithms/Optimization/Constraints/BoxConstrainedPowell.h>

#endif // MML_OPTIMIZATION_MULTIDIM_H

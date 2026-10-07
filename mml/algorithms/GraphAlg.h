///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        GraphAlg.h                                                          ///
///  Description: Aggregate header for graph algorithms                              ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

/// @file GraphAlg.h
/// @brief Aggregate header for graph algorithms.
/// This header includes graph traversal, shortest paths, structural analysis,
/// network flow, spanning tree, and spectral graph algorithms.
///
/// Main algorithms included here:
/// - Breadth-first and depth-first graph traversal helpers
/// - Dijkstra, Bellman-Ford, Floyd-Warshall, and A* shortest-path routines
/// - Connectivity, cycle, topological-sort, and graph-structure analysis
/// - Max-flow and minimum-cut algorithms
/// - Kruskal and Prim spanning-tree algorithms
/// - Spectral graph helpers built on matrix representations
///
/// For more focused includes, use the individual headers under Graphs/.

#if !defined MML_ALGORITHMS_GRAPHALG_H
#define MML_ALGORITHMS_GRAPHALG_H

#include <mml/algorithms/Graphs/GraphTraversal.h>
#include <mml/algorithms/Graphs/GraphShortestPaths.h>
#include <mml/algorithms/Graphs/GraphStructure.h>
#include <mml/algorithms/Graphs/GraphFlow.h>
#include <mml/algorithms/Graphs/GraphSpanningTree.h>

#endif

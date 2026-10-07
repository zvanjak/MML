# Graph Algorithms

MML provides an adjacency-list graph type plus traversal, shortest-path, connectivity, DAG, spanning-tree, flow, cut, matching, and matrix-conversion helpers. The algorithms are header-only and operate on `MML::Graph<V, E>`.

Main headers:

- `mml/base/Graph.h`
- `mml/algorithms/GraphAlg.h`
- `mml/algorithms/Graphs/GraphTraversal.h`
- `mml/algorithms/Graphs/GraphShortestPaths.h`
- `mml/algorithms/Graphs/GraphStructure.h`
- `mml/algorithms/Graphs/GraphSpanningTree.h`
- `mml/algorithms/Graphs/GraphFlow.h`

Include `mml/algorithms/GraphAlg.h` for the aggregate graph-algorithm header, or include individual headers when you only need one family.

## Graph Type

`Graph<V, E>` stores vertex data of type `V` and edge weights of type `E`. Defaults are `Graph<int, Real>`.

| API | Purpose |
|-----|---------|
| `Graph<>(type, allowSelfLoops)` | Create an empty graph; type defaults to undirected |
| `Graph<>(numVertices, type, allowSelfLoops)` | Create a graph with indexed vertices `0..n-1` |
| `addVertex(data)` | Add a vertex and return its index |
| `addEdge(from, to, weight)` | Add an edge; undirected graphs add the reverse adjacency entry automatically |
| `hasEdge(from, to)` / `edgeWeight(from, to)` | Query edge existence and weight |
| `neighbors(v)` | Access outgoing adjacency list |
| `edges()` | Return edge list, with undirected edges represented once |
| `directedEdges()` | Return directed edge entries |
| `toAdjacencyMatrix()` | Convert to adjacency matrix |
| `toDegreeMatrix()` | Convert to weighted degree matrix |
| `toLaplacianMatrix()` | Convert to graph Laplacian matrix |
| `toNormalizedLaplacian()` | Convert to normalized Laplacian matrix |

Factory helpers in `mml/base/Graph.h` create common undirected graphs:

- `createPathGraph<V, E>(n)`
- `createCycleGraph<V, E>(n)`
- `createCompleteGraph<V, E>(n)`
- `createStarGraph<V, E>(n)`
- `createGridGraph<V, E>(rows, cols)`

```cpp
#include <mml/algorithms/GraphAlg.h>
using namespace MML;

Graph<int, Real> graph(4);
graph.addEdge(0, 1, 2.0);
graph.addEdge(1, 2, 3.0);
graph.addEdge(2, 3, 1.0);

auto path = DijkstraPath(graph, 0, 3);
```

## Result Objects and Status

Most algorithms return a structured result rather than throwing for ordinary graph-algorithm outcomes. Result types include:

| Result type | Used by |
|-------------|---------|
| `TraversalResult` | BFS, DFS, Dijkstra, Bellman-Ford |
| `PathResult` | Targeted shortest-path routines |
| `ShortestPathTreeResult` | Single-source path-tree wrappers |
| `AllPairsShortestPathsResult` | All-pairs shortest paths |
| `ComponentsResult` | Undirected connected components |
| `MSTResult` | Kruskal and Prim |
| `TopologicalSortResult` | Topological sorting |
| `StronglyConnectedComponentsResult` | SCC and condensation DAG analysis |
| `UndirectedConnectivityResult` | Articulation points, bridges, and biconnected components |
| `CycleResult` | Cycle detection |
| `BipartiteResult` | Bipartite check and two-coloring |
| `TransitiveResult` | Transitive closure/reduction |
| `LongestPathDAGResult` | Longest paths in DAGs |
| `MaxFlowMinCutResult` | Max-flow/min-cut algorithms |
| `BipartiteMatchingResult` | Hopcroft-Karp matching |

Graph-specific outcomes are reported with `GraphAlgorithmStatus`:

- `Success`
- `InvalidVertex`
- `UnreachableTarget`
- `NegativeCycle`
- `CycleDetected`
- `DisconnectedGraph`
- `GraphTypeMismatch`
- `UnsupportedNegativeWeights`

Result objects provide `succeeded()` when the operation uses both general `AlgorithmStatus` and graph-specific status.

## Traversal and Components

| Function | Description |
|----------|-------------|
| `BFS(graph, start)` | Breadth-first traversal from `start`; returns visit order, parents, hop distances, and node count |
| `BFSVisit(graph, start, visitor)` | BFS with a callback; returning `false` stops traversal |
| `DFS(graph, start)` | Iterative depth-first traversal |
| `DFSVisit(graph, start, visitor)` | DFS with a callback; returning `false` stops traversal |
| `ConnectedComponents(graph)` | Undirected connected components |
| `IsConnected(graph)` | Boolean connectivity check for undirected graphs |

```cpp
auto graph = createPathGraph<int, Real>(5);
auto bfs = BFS(graph, 0);
// bfs.distance[4] == 4 on this path graph.
```

## Shortest Paths

| Function | Description |
|----------|-------------|
| `Dijkstra(graph, start)` | Single-source shortest paths for non-negative weights |
| `DijkstraShortestPathTree(graph, start)` | Dijkstra result packaged as a `ShortestPathTreeResult` |
| `DijkstraPath(graph, start, target)` | Targeted Dijkstra path and total weight |
| `ShortestPathUnweighted(graph, start, target)` | BFS-based shortest path in unweighted graphs |
| `AStarPath(graph, start, target, heuristic)` | A* search with caller-supplied heuristic |
| `BellmanFord(graph, start)` | Single-source shortest paths with negative edges and negative-cycle detection |
| `FloydWarshall(graph)` | All-pairs shortest paths with negative-cycle detection |
| `HasNegativeCycle(graph, start)` | Negative-cycle check reachable from a start vertex |

`Dijkstra` and `AStarPath` reject negative weights with `UnsupportedNegativeWeights`. `BellmanFord` and `FloydWarshall` report `NegativeCycle` when applicable.

## Structural Algorithms

| Function | Description |
|----------|-------------|
| `TopologicalSort(graph)` | Kahn topological sort for directed acyclic graphs |
| `IsDAG(graph)` | Boolean DAG check using topological sort |
| `StronglyConnectedComponents(graph)` | Tarjan SCC decomposition for directed graphs |
| `CondensationDAG(graph)` | SCC result with condensation edge list |
| `UndirectedConnectivity(graph)` | Articulation points, bridges, and biconnected edge components |
| `FindCycle(graph)` | Cycle detection with one returned cycle when found |
| `BipartiteColoring(graph)` | Bipartite check and two-coloring |
| `TransitiveClosure(graph)` | Directed reachability closure |
| `TransitiveReductionDAG(graph)` | Transitive reduction edge list for DAGs |
| `LongestPathDAG(graph, source)` | Longest paths from a source in a DAG |

Directed-only routines return `GraphTypeMismatch` when called on an undirected graph; undirected-only routines do the same when called on a directed graph.

## Spanning Trees, Flow, and Matching

| Function | Description |
|----------|-------------|
| `Kruskal(graph)` | Minimum spanning tree using Kruskal's algorithm |
| `Prim(graph, start)` | Minimum spanning tree using Prim's algorithm |
| `EdmondsKarp(graph, source, sink)` | Directed max-flow/min-cut using Edmonds-Karp |
| `Dinic(graph, source, sink)` | Directed max-flow/min-cut using Dinic's algorithm |
| `HopcroftKarp(graph, leftPartition)` | Maximum cardinality bipartite matching with explicit left partition |
| `HopcroftKarp(graph)` | Maximum cardinality bipartite matching after deriving a two-coloring |

MST routines report `DisconnectedGraph` if the result cannot span all vertices. Flow routines require directed graphs and non-negative capacities. Matching routines require a bipartite graph or a valid explicit partition.

## Matrix Conversions

The graph class can produce matrix forms used by spectral or linear-algebra workflows:

```cpp
auto graph = createCycleGraph<int, Real>(4);
Matrix<Real> adjacency = graph.toAdjacencyMatrix();
Matrix<Real> laplacian = graph.toLaplacianMatrix();
Matrix<Real> normalized = graph.toNormalizedLaplacian();
```

The current graph-algorithm aggregate includes matrix conversion support through `Graph`, but it does not expose a separate named spectral-analysis algorithm header.

## See Also

- [Algebra & discrete mathematics](../base/Algebra.md)
- [Matrices](../base/Matrices.md)
- [Sparse matrices](../base/SparseMatrices.md)
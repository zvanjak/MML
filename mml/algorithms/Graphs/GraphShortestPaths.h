///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        GraphShortestPaths.h                                                ///
///  Description: Single-source and all-pairs shortest-path algorithms               ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_GRAPH_ALGORITHMS_GRAPHS_GRAPHSHORTESTPATHS_H
#define MML_GRAPH_ALGORITHMS_GRAPHS_GRAPHSHORTESTPATHS_H

#include <mml/algorithms/Graphs/GraphTraversal.h>
#include <algorithm>
#include <functional>
#include <limits>
#include <queue>
#include <tuple>
#include <utility>
#include <vector>

namespace MML
{

	///////////////////////////////////////////////////////////////////////////
	///                    DIJKSTRA'S ALGORITHM                             ///
	///////////////////////////////////////////////////////////////////////////

	/// Find shortest path from start to all vertices using Dijkstra's algorithm
	/// Requires non-negative edge weights
	/// Complexity: O((V + E) log V) with priority queue
	///
	/// @param graph  The graph (must have non-negative weights)
	/// @param start  Starting vertex index
	/// @return TraversalResult with distances and parent pointers for path reconstruction
	template<typename V, typename E>
	TraversalResult Dijkstra(const Graph<V, E>& graph, size_t start)
	{
		auto startTime = std::chrono::steady_clock::now();
		TraversalResult result;
		result.algorithm_name = "Dijkstra";
		size_t n = graph.numVertices();

		if (start >= n)
		{
			result.status = AlgorithmStatus::InvalidInput;
			result.graphStatus = GraphAlgorithmStatus::InvalidVertex;
			detail::SetGraphMessage(result.message, result.diagnostics, "Invalid start vertex");
			result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
			return result;
		}

		if (detail::HasNegativeWeight(graph))
		{
			result.status = AlgorithmStatus::InvalidInput;
			result.graphStatus = GraphAlgorithmStatus::UnsupportedNegativeWeights;
			detail::SetGraphMessage(result.message, result.diagnostics, "Dijkstra requires non-negative edge weights");
			result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
			return result;
		}

		result.parent.resize(n, GRAPH_NOT_FOUND);
		result.distance.resize(n, std::numeric_limits<Real>::infinity());
		result.distance[start] = 0;

		// Priority queue: (distance, vertex)
		using PQEntry = std::pair<Real, size_t>;
		std::priority_queue<PQEntry, std::vector<PQEntry>, std::greater<PQEntry>> pq;

		pq.push({REAL(0.0), start});

		while (!pq.empty())
		{
			auto [dist, u] = pq.top();
			pq.pop();

			// Skip if we've already found a better path
			if (dist > result.distance[u])
				continue;

			result.nodesVisited++;

			for (const auto& edge : graph.neighbors(u))
			{
				size_t v = edge.to;
				Real newDist = result.distance[u] + static_cast<Real>(edge.weight);

				if (newDist < result.distance[v])
				{
					result.distance[v] = newDist;
					result.parent[v] = u;
					pq.push({newDist, v});
				}
			}
		}

		result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
		return result;
	}

	/// Find a shortest-path tree from start using Dijkstra's algorithm.
	template<typename V, typename E>
	ShortestPathTreeResult DijkstraShortestPathTree(const Graph<V, E>& graph, size_t start)
	{
		auto traversal = Dijkstra(graph, start);
		return detail::MakeShortestPathTree(traversal, start, "DijkstraShortestPathTree");
	}

	/// Find shortest path from start to specific target using Dijkstra
	/// Early termination when target is reached
	template<typename V, typename E>
	PathResult DijkstraPath(const Graph<V, E>& graph, size_t start, size_t target)
	{
		auto startTime = std::chrono::steady_clock::now();
		PathResult result;
		result.algorithm_name = "DijkstraPath";
		size_t n = graph.numVertices();

		if (start >= n || target >= n)
		{
			result.status = AlgorithmStatus::InvalidInput;
			result.graphStatus = GraphAlgorithmStatus::InvalidVertex;
			detail::SetGraphMessage(result.message, result.diagnostics, "Invalid start or target vertex");
			result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
			return result;
		}

		if (detail::HasNegativeWeight(graph))
		{
			result.status = AlgorithmStatus::InvalidInput;
			result.graphStatus = GraphAlgorithmStatus::UnsupportedNegativeWeights;
			detail::SetGraphMessage(result.message, result.diagnostics, "Dijkstra requires non-negative edge weights");
			result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
			return result;
		}

		if (start == target)
		{
			result.found = true;
			result.path.push_back(start);
			result.totalWeight = 0;
			result.nodesExplored = 1;
			result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
			return result;
		}

		std::vector<Real> dist(n, std::numeric_limits<Real>::infinity());
		std::vector<size_t> parent(n, GRAPH_NOT_FOUND);
		dist[start] = 0;

		using PQEntry = std::pair<Real, size_t>;
		std::priority_queue<PQEntry, std::vector<PQEntry>, std::greater<PQEntry>> pq;
		pq.push({REAL(0.0), start});

		while (!pq.empty())
		{
			auto [d, u] = pq.top();
			pq.pop();

			result.nodesExplored++;

			if (u == target)
			{
				// Reconstruct path
				result.found = true;
				result.totalWeight = dist[target];

				size_t curr = target;
				while (curr != GRAPH_NOT_FOUND)
				{
					result.path.push_back(curr);
					curr = parent[curr];
				}
				std::reverse(result.path.begin(), result.path.end());
				result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
				return result;
			}

			if (d > dist[u])
				continue;

			for (const auto& edge : graph.neighbors(u))
			{
				size_t v = edge.to;
				Real newDist = dist[u] + static_cast<Real>(edge.weight);

				if (newDist < dist[v])
				{
					dist[v] = newDist;
					parent[v] = u;
					pq.push({newDist, v});
				}
			}
		}

		result.found = false;
		result.status = AlgorithmStatus::AlgorithmSpecificFailure;
		result.graphStatus = GraphAlgorithmStatus::UnreachableTarget;
		detail::SetGraphMessage(result.message, result.diagnostics, "No path exists from start to target");
		result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
		return result;
	}

	///////////////////////////////////////////////////////////////////////////
	///                    PATH UTILITIES                                   ///
	///////////////////////////////////////////////////////////////////////////

	/// Reconstruct path from parent array (used by BFS/DFS/Dijkstra)
	inline std::vector<size_t> ReconstructPath(
		const std::vector<size_t>& parent,
		size_t start,
		size_t target)
	{
		std::vector<size_t> path;

		if (target >= parent.size())
			return path;

		if (parent[target] == GRAPH_NOT_FOUND && target != start)
			return path;  // No path exists

		size_t curr = target;
		while (curr != GRAPH_NOT_FOUND)
		{
			path.push_back(curr);
			if (curr == start)
				break;
			curr = parent[curr];
		}

		std::reverse(path.begin(), path.end());
		return path;
	}

	/// Find shortest path (unweighted) using BFS
	template<typename V, typename E>
	PathResult ShortestPathUnweighted(const Graph<V, E>& graph, size_t start, size_t target)
	{
		auto startTime = std::chrono::steady_clock::now();
		PathResult result;
		result.algorithm_name = "ShortestPathUnweighted";
		size_t n = graph.numVertices();

		if (start >= n || target >= n)
		{
			result.status = AlgorithmStatus::InvalidInput;
			result.graphStatus = GraphAlgorithmStatus::InvalidVertex;
			detail::SetGraphMessage(result.message, result.diagnostics, "Invalid start or target vertex");
			result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
			return result;
		}

		if (start == target)
		{
			result.found = true;
			result.path.push_back(start);
			result.totalWeight = 0;
			result.nodesExplored = 1;
			result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
			return result;
		}

		auto bfsResult = BFS(graph, start);
		result.nodesExplored = bfsResult.nodesVisited;

		if (bfsResult.distance[target] == std::numeric_limits<Real>::infinity())
		{
			result.found = false;
			result.status = AlgorithmStatus::AlgorithmSpecificFailure;
			result.graphStatus = GraphAlgorithmStatus::UnreachableTarget;
			detail::SetGraphMessage(result.message, result.diagnostics, "No path exists from start to target");
			result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
			return result;
		}

		result.found = true;
		result.path = ReconstructPath(bfsResult.parent, start, target);
		result.totalWeight = bfsResult.distance[target];
		result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);

		return result;
	}

	/// Find shortest path from start to target using A* with a user-provided heuristic.
	/// The heuristic is called as heuristic(vertex, target) and should be admissible for optimality.
	template<typename V, typename E, typename Heuristic>
	PathResult AStarPath(const Graph<V, E>& graph, size_t start, size_t target, Heuristic heuristic)
	{
		auto startTime = std::chrono::steady_clock::now();
		PathResult result;
		result.algorithm_name = "AStarPath";
		size_t n = graph.numVertices();

		if (start >= n || target >= n)
		{
			result.status = AlgorithmStatus::InvalidInput;
			result.graphStatus = GraphAlgorithmStatus::InvalidVertex;
			detail::SetGraphMessage(result.message, result.diagnostics, "Invalid start or target vertex");
			result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
			return result;
		}

		if (detail::HasNegativeWeight(graph))
		{
			result.status = AlgorithmStatus::InvalidInput;
			result.graphStatus = GraphAlgorithmStatus::UnsupportedNegativeWeights;
			detail::SetGraphMessage(result.message, result.diagnostics, "A* requires non-negative edge weights");
			result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
			return result;
		}

		if (start == target)
		{
			result.found = true;
			result.path.push_back(start);
			result.totalWeight = 0;
			result.nodesExplored = 1;
			result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
			return result;
		}

		std::vector<Real> gScore(n, std::numeric_limits<Real>::infinity());
		std::vector<size_t> parent(n, GRAPH_NOT_FOUND);
		std::vector<bool> closed(n, false);
		gScore[start] = 0;

		using PQEntry = std::pair<Real, size_t>;
		std::priority_queue<PQEntry, std::vector<PQEntry>, std::greater<PQEntry>> open;
		open.push({ static_cast<Real>(heuristic(start, target)), start });

		while (!open.empty())
		{
			auto [fScore, current] = open.top();
			(void)fScore;
			open.pop();

			if (closed[current])
				continue;

			closed[current] = true;
			++result.nodesExplored;

			if (current == target)
			{
				result.found = true;
				result.totalWeight = gScore[target];
				result.path = ReconstructPath(parent, start, target);
				result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
				return result;
			}

			for (const auto& edge : graph.neighbors(current))
			{
				size_t neighbor = edge.to;
				Real tentative = gScore[current] + static_cast<Real>(edge.weight);
				if (tentative < gScore[neighbor])
				{
					gScore[neighbor] = tentative;
					parent[neighbor] = current;
					Real estimated = tentative + static_cast<Real>(heuristic(neighbor, target));
					open.push({ estimated, neighbor });
				}
			}
		}

		result.status = AlgorithmStatus::AlgorithmSpecificFailure;
		result.graphStatus = GraphAlgorithmStatus::UnreachableTarget;
		detail::SetGraphMessage(result.message, result.diagnostics, "No path exists from start to target");
		result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
		return result;
	}

	/// Compute all-pairs shortest paths using Floyd-Warshall.
	/// Supports negative edge weights and reports a negative cycle if one exists.
	template<typename V, typename E>
	AllPairsShortestPathsResult FloydWarshall(const Graph<V, E>& graph)
	{
		auto startTime = std::chrono::steady_clock::now();
		AllPairsShortestPathsResult result;
		result.algorithm_name = "FloydWarshall";
		size_t n = graph.numVertices();
		Real inf = std::numeric_limits<Real>::infinity();

		result.distance.assign(n, std::vector<Real>(n, inf));
		result.next.assign(n, std::vector<size_t>(n, GRAPH_NOT_FOUND));

		for (size_t i = 0; i < n; ++i)
		{
			result.distance[i][i] = 0;
			result.next[i][i] = i;
		}

		for (const auto& [from, to, weight] : graph.directedEdges())
		{
			Real w = static_cast<Real>(weight);
			if (w < result.distance[from][to])
			{
				result.distance[from][to] = w;
				result.next[from][to] = to;
			}
		}

		for (size_t k = 0; k < n; ++k)
		{
			for (size_t i = 0; i < n; ++i)
			{
				if (result.distance[i][k] == inf)
					continue;
				for (size_t j = 0; j < n; ++j)
				{
					if (result.distance[k][j] == inf)
						continue;
					Real candidate = result.distance[i][k] + result.distance[k][j];
					if (candidate < result.distance[i][j])
					{
						result.distance[i][j] = candidate;
						result.next[i][j] = result.next[i][k];
					}
				}
			}
		}

		for (size_t i = 0; i < n; ++i)
		{
			if (result.distance[i][i] < 0)
			{
				result.status = AlgorithmStatus::AlgorithmSpecificFailure;
				result.graphStatus = GraphAlgorithmStatus::NegativeCycle;
				detail::SetGraphMessage(result.message, result.diagnostics, "Negative cycle detected");
				result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
				return result;
			}
		}

		result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
		return result;
	}

	///////////////////////////////////////////////////////////////////////////
	///                    BELLMAN-FORD ALGORITHM                           ///
	///////////////////////////////////////////////////////////////////////////

	/// Find shortest paths from start vertex using Bellman-Ford algorithm
	/// Handles negative edge weights, detects negative cycles
	/// Complexity: O(V * E)
	///
	/// @param graph  The graph (can have negative weights)
	/// @param start  Starting vertex index
	/// @return TraversalResult with distances and parent pointers
	///         If negative cycle detected, returns empty result with diagnostics
	template<typename V, typename E>
	TraversalResult BellmanFord(const Graph<V, E>& graph, size_t start)
	{
		auto startTime = std::chrono::steady_clock::now();
		TraversalResult result;
		result.algorithm_name = "BellmanFord";
		size_t n = graph.numVertices();

		if (start >= n)
		{
			result.status = AlgorithmStatus::InvalidInput;
			result.graphStatus = GraphAlgorithmStatus::InvalidVertex;
			detail::SetGraphMessage(result.message, result.diagnostics, "Invalid start vertex");
			result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
			return result;
		}

		result.parent.resize(n, GRAPH_NOT_FOUND);
		result.distance.resize(n, std::numeric_limits<Real>::infinity());
		result.distance[start] = 0;

		std::vector<std::tuple<size_t, size_t, E>> edges = graph.directedEdges();

		// Relax all edges (V-1) times
		for (size_t i = 0; i < n - 1; ++i)
		{
			bool changed = false;
			for (const auto& [u, v, w] : edges)
			{
				if (result.distance[u] != std::numeric_limits<Real>::infinity())
				{
					Real newDist = result.distance[u] + static_cast<Real>(w);
					if (newDist < result.distance[v])
					{
						result.distance[v] = newDist;
						result.parent[v] = u;
						changed = true;
					}
				}
			}
			// Early termination if no changes
			if (!changed)
				break;
		}

		// Check for negative cycles
		for (const auto& [u, v, w] : edges)
		{
			if (result.distance[u] != std::numeric_limits<Real>::infinity())
			{
				if (result.distance[u] + static_cast<Real>(w) < result.distance[v])
				{
					// Negative cycle detected
					result.status = AlgorithmStatus::AlgorithmSpecificFailure;
					result.graphStatus = GraphAlgorithmStatus::NegativeCycle;
					detail::SetGraphMessage(result.message, result.diagnostics, "Negative cycle detected");
					result.parent.clear();
					result.distance.clear();
					result.visitOrder.clear();
					result.nodesVisited = 0;
					result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
					return result;
				}
			}
		}

		// Count reachable nodes
		for (size_t i = 0; i < n; ++i)
		{
			if (result.distance[i] != std::numeric_limits<Real>::infinity())
				result.nodesVisited++;
		}

		result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
		return result;
	}

	/// Check if graph contains a negative cycle reachable from start
	template<typename V, typename E>
	bool HasNegativeCycle(const Graph<V, E>& graph, size_t start)
	{
		auto result = BellmanFord(graph, start);
		return result.graphStatus == GraphAlgorithmStatus::NegativeCycle;
	}

	/// Compute all-pairs shortest paths using Johnson's algorithm.
	/// Supports negative edge weights but reports negative cycles explicitly.
	template<typename V, typename E>
	AllPairsShortestPathsResult Johnson(const Graph<V, E>& graph)
	{
		auto startTime = std::chrono::steady_clock::now();
		AllPairsShortestPathsResult result;
		result.algorithm_name = "Johnson";
		size_t n = graph.numVertices();
		Real inf = std::numeric_limits<Real>::infinity();

		result.distance.assign(n, std::vector<Real>(n, inf));
		result.next.assign(n, std::vector<size_t>(n, GRAPH_NOT_FOUND));
		if (n == 0)
		{
			result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
			return result;
		}

		std::vector<Real> h(n, 0);
		auto edges = graph.directedEdges();

		for (size_t i = 0; i < n - 1; ++i)
		{
			bool changed = false;
			for (const auto& [u, v, weight] : edges)
			{
				Real candidate = h[u] + static_cast<Real>(weight);
				if (candidate < h[v])
				{
					h[v] = candidate;
					changed = true;
				}
			}
			if (!changed)
				break;
		}

		for (const auto& [u, v, weight] : edges)
		{
			if (h[u] + static_cast<Real>(weight) < h[v])
			{
				result.status = AlgorithmStatus::AlgorithmSpecificFailure;
				result.graphStatus = GraphAlgorithmStatus::NegativeCycle;
				detail::SetGraphMessage(result.message, result.diagnostics, "Negative cycle detected");
				result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
				return result;
			}
		}

		for (size_t source = 0; source < n; ++source)
		{
			std::vector<Real> dist(n, inf);
			std::vector<size_t> parent(n, GRAPH_NOT_FOUND);
			dist[source] = 0;

			using PQEntry = std::pair<Real, size_t>;
			std::priority_queue<PQEntry, std::vector<PQEntry>, std::greater<PQEntry>> pq;
			pq.push({REAL(0.0), source});

			while (!pq.empty())
			{
				auto [distU, u] = pq.top();
				pq.pop();
				if (distU > dist[u])
					continue;

				for (const auto& edge : graph.neighbors(u))
				{
					size_t v = edge.to;
					Real reweighted = static_cast<Real>(edge.weight) + h[u] - h[v];
					Real candidate = dist[u] + reweighted;
					if (candidate < dist[v])
					{
						dist[v] = candidate;
						parent[v] = u;
						pq.push({ candidate, v });
					}
				}
			}

			for (size_t target = 0; target < n; ++target)
			{
				if (source == target)
				{
					result.distance[source][target] = 0;
					result.next[source][target] = target;
				}
				else if (dist[target] != inf)
				{
					result.distance[source][target] = dist[target] - h[source] + h[target];
					auto path = ReconstructPath(parent, source, target);
					if (path.size() >= 2)
						result.next[source][target] = path[1];
				}
			}
		}

		result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
		return result;
	}

}

#endif

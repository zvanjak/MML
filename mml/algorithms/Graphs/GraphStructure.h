///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        GraphStructure.h                                                    ///
///  Description: Connectivity, cycle, and DAG structure algorithms                  ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_GRAPH_ALGORITHMS_GRAPHS_GRAPHSTRUCTURE_H
#define MML_GRAPH_ALGORITHMS_GRAPHS_GRAPHSTRUCTURE_H

#include <mml/algorithms/Graphs/detail/GraphAlgorithmUtils.h>
#include <algorithm>
#include <functional>
#include <limits>
#include <queue>
#include <set>
#include <utility>
#include <vector>

namespace MML
{

	///////////////////////////////////////////////////////////////////////////
	///                    TOPOLOGICAL SORT                                 ///
	///////////////////////////////////////////////////////////////////////////

	/// Perform topological sort on a directed acyclic graph (DAG)
	/// Uses Kahn's algorithm (BFS-based)
	/// Complexity: O(V + E)
	///
	/// @param graph  The graph (must be directed and acyclic)
	/// @return TopologicalSortResult with sorted order and cycle detection
	template<typename V, typename E>
	TopologicalSortResult TopologicalSort(const Graph<V, E>& graph)
	{
		auto startTime = std::chrono::steady_clock::now();
		TopologicalSortResult result;
		result.algorithm_name = "TopologicalSort";
		size_t n = graph.numVertices();

		if (n == 0)
		{
			result.isDAG = true;
			result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
			return result;
		}

		if (!graph.isDirected())
		{
			result.status = AlgorithmStatus::InvalidInput;
			result.graphStatus = GraphAlgorithmStatus::GraphTypeMismatch;
			result.isDAG = false;
			detail::SetGraphMessage(result.message, result.diagnostics, "Topological sort requires a directed graph");
			result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
			return result;
		}

		// Calculate in-degrees
		std::vector<size_t> inDegree(n, 0);
		for (size_t u = 0; u < n; ++u)
		{
			for (const auto& edge : graph.neighbors(u))
			{
				inDegree[edge.to]++;
			}
		}

		// Initialize queue with zero in-degree vertices
		std::queue<size_t> queue;
		for (size_t v = 0; v < n; ++v)
		{
			if (inDegree[v] == 0)
				queue.push(v);
		}

		// Process vertices
		while (!queue.empty())
		{
			size_t u = queue.front();
			queue.pop();
			result.order.push_back(u);

			for (const auto& edge : graph.neighbors(u))
			{
				if (--inDegree[edge.to] == 0)
					queue.push(edge.to);
			}
		}

		// Check if all vertices were processed
		if (result.order.size() != n)
		{
			result.status = AlgorithmStatus::AlgorithmSpecificFailure;
			result.graphStatus = GraphAlgorithmStatus::CycleDetected;
			result.isDAG = false;
			result.order.clear();
			detail::SetGraphMessage(result.message, result.diagnostics, "Graph contains a cycle");
			result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
			return result;
		}

		result.isDAG = true;
		result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
		return result;
	}

	/// Check if a directed graph is a DAG (has no cycles)
	template<typename V, typename E>
	bool IsDAG(const Graph<V, E>& graph)
	{
		auto result = TopologicalSort(graph);
		return result.isDAG;
	}

	///////////////////////////////////////////////////////////////////////////
	///                    STRUCTURAL GRAPH ALGORITHMS                      ///
	///////////////////////////////////////////////////////////////////////////

	/// Find strongly connected components in a directed graph using Tarjan's algorithm.
	template<typename V, typename E>
	StronglyConnectedComponentsResult StronglyConnectedComponents(const Graph<V, E>& graph)
	{
		auto startTime = std::chrono::steady_clock::now();
		StronglyConnectedComponentsResult result;
		result.algorithm_name = "StronglyConnectedComponents";
		size_t n = graph.numVertices();

		if (!graph.isDirected())
		{
			result.status = AlgorithmStatus::InvalidInput;
			result.graphStatus = GraphAlgorithmStatus::GraphTypeMismatch;
			detail::SetGraphMessage(result.message, result.diagnostics, "Strongly connected components require a directed graph");
			result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
			return result;
		}

		result.componentId.assign(n, GRAPH_NOT_FOUND);
		std::vector<size_t> index(n, GRAPH_NOT_FOUND), lowlink(n, 0), stack;
		std::vector<bool> onStack(n, false);
		size_t nextIndex = 0;

		std::function<void(size_t)> dfs = [&](size_t v)
		{
			index[v] = lowlink[v] = nextIndex++;
			stack.push_back(v);
			onStack[v] = true;

			for (const auto& edge : graph.neighbors(v))
			{
				size_t w = edge.to;
				if (index[w] == GRAPH_NOT_FOUND)
				{
					dfs(w);
					lowlink[v] = std::min(lowlink[v], lowlink[w]);
				}
				else if (onStack[w])
				{
					lowlink[v] = std::min(lowlink[v], index[w]);
				}
			}

			if (lowlink[v] == index[v])
			{
				std::vector<size_t> component;
				while (true)
				{
					size_t w = stack.back();
					stack.pop_back();
					onStack[w] = false;
					result.componentId[w] = result.components.size();
					component.push_back(w);
					if (w == v)
						break;
				}
				std::sort(component.begin(), component.end());
				result.components.push_back(std::move(component));
			}
		};

		for (size_t v = 0; v < n; ++v)
		{
			if (index[v] == GRAPH_NOT_FOUND)
				dfs(v);
		}

		result.numComponents = result.components.size();
		std::set<std::pair<size_t, size_t>> condensation;
		for (const auto& [u, v, weight] : graph.directedEdges())
		{
			(void)weight;
			size_t cu = result.componentId[u];
			size_t cv = result.componentId[v];
			if (cu != cv)
				condensation.insert({ cu, cv });
		}
		result.condensationEdges.assign(condensation.begin(), condensation.end());
		result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
		return result;
	}

	/// Build condensation DAG edge list for a directed graph.
	template<typename V, typename E>
	StronglyConnectedComponentsResult CondensationDAG(const Graph<V, E>& graph)
	{
		auto result = StronglyConnectedComponents(graph);
		result.algorithm_name = "CondensationDAG";
		return result;
	}

	/// Find articulation points, bridges, and biconnected edge components in an undirected graph.
	template<typename V, typename E>
	UndirectedConnectivityResult UndirectedConnectivity(const Graph<V, E>& graph)
	{
		auto startTime = std::chrono::steady_clock::now();
		UndirectedConnectivityResult result;
		result.algorithm_name = "UndirectedConnectivity";
		size_t n = graph.numVertices();

		if (!graph.isUndirected())
		{
			result.status = AlgorithmStatus::InvalidInput;
			result.graphStatus = GraphAlgorithmStatus::GraphTypeMismatch;
			detail::SetGraphMessage(result.message, result.diagnostics, "Undirected connectivity analysis requires an undirected graph");
			result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
			return result;
		}

		std::vector<size_t> disc(n, GRAPH_NOT_FOUND), low(n, 0), parent(n, GRAPH_NOT_FOUND);
		std::vector<bool> articulation(n, false);
		std::vector<std::pair<size_t, size_t>> edgeStack;
		size_t time = 0;

		auto normalized = [](size_t a, size_t b) {
			return std::pair<size_t, size_t>{std::min(a, b), std::max(a, b)};
		};
		std::function<void(size_t)> flushComponentTo = [&](size_t u)
		{
			std::vector<std::pair<size_t, size_t>> component;
			while (!edgeStack.empty())
			{
				auto edge = edgeStack.back();
				edgeStack.pop_back();
				component.push_back(edge);
				if (edge.first == std::min(u, edge.second) && edge.second == std::max(u, edge.first))
					break;
			}
			if (!component.empty())
			{
				std::sort(component.begin(), component.end());
				result.biconnectedComponents.push_back(std::move(component));
			}
		};

		std::function<void(size_t)> dfs = [&](size_t u)
		{
			disc[u] = low[u] = time++;
			size_t children = 0;

			for (const auto& graphEdge : graph.neighbors(u))
			{
				size_t v = graphEdge.to;
				if (disc[v] == GRAPH_NOT_FOUND)
				{
					++children;
					parent[v] = u;
					edgeStack.push_back(normalized(u, v));
					dfs(v);
					low[u] = std::min(low[u], low[v]);

					if ((parent[u] == GRAPH_NOT_FOUND && children > 1) ||
						(parent[u] != GRAPH_NOT_FOUND && low[v] >= disc[u]))
					{
						articulation[u] = true;
						std::vector<std::pair<size_t, size_t>> component;
						auto stop = normalized(u, v);
						while (!edgeStack.empty())
						{
							auto edge = edgeStack.back();
							edgeStack.pop_back();
							component.push_back(edge);
							if (edge == stop)
								break;
						}
						std::sort(component.begin(), component.end());
						result.biconnectedComponents.push_back(std::move(component));
					}

					if (low[v] > disc[u])
						result.bridges.push_back(normalized(u, v));
				}
				else if (v != parent[u] && disc[v] < disc[u])
				{
					low[u] = std::min(low[u], disc[v]);
					edgeStack.push_back(normalized(u, v));
				}
			}
		};

		for (size_t v = 0; v < n; ++v)
		{
			if (disc[v] == GRAPH_NOT_FOUND)
			{
				dfs(v);
				if (!edgeStack.empty())
				{
					std::vector<std::pair<size_t, size_t>> component(edgeStack.begin(), edgeStack.end());
					edgeStack.clear();
					std::sort(component.begin(), component.end());
					result.biconnectedComponents.push_back(std::move(component));
				}
			}
		}

		for (size_t v = 0; v < n; ++v)
		{
			if (articulation[v])
				result.articulationPoints.push_back(v);
		}
		std::sort(result.bridges.begin(), result.bridges.end());
		result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
		return result;
	}

	/// Detect a cycle and return one cycle when present.
	template<typename V, typename E>
	CycleResult FindCycle(const Graph<V, E>& graph)
	{
		auto startTime = std::chrono::steady_clock::now();
		CycleResult result;
		result.algorithm_name = "FindCycle";
		size_t n = graph.numVertices();
		std::vector<size_t> parent(n, GRAPH_NOT_FOUND);

		if (graph.isDirected())
		{
			std::vector<int> color(n, 0);
			std::function<bool(size_t)> dfs = [&](size_t u)
			{
				color[u] = 1;
				for (const auto& edge : graph.neighbors(u))
				{
					size_t v = edge.to;
					if (color[v] == 0)
					{
						parent[v] = u;
						if (dfs(v))
							return true;
					}
					else if (color[v] == 1)
					{
						result.hasCycle = true;
						result.cycle.push_back(v);
						for (size_t curr = u; curr != v && curr != GRAPH_NOT_FOUND; curr = parent[curr])
							result.cycle.push_back(curr);
						result.cycle.push_back(v);
						std::reverse(result.cycle.begin(), result.cycle.end());
						return true;
					}
				}
				color[u] = 2;
				return false;
			};

			for (size_t v = 0; v < n && !result.hasCycle; ++v)
				if (color[v] == 0)
					dfs(v);
		}
		else
		{
			std::vector<bool> visited(n, false);
			std::function<bool(size_t)> dfs = [&](size_t u)
			{
				visited[u] = true;
				for (const auto& edge : graph.neighbors(u))
				{
					size_t v = edge.to;
					if (!visited[v])
					{
						parent[v] = u;
						if (dfs(v))
							return true;
					}
					else if (v != parent[u])
					{
						result.hasCycle = true;
						result.cycle.push_back(v);
						for (size_t curr = u; curr != v && curr != GRAPH_NOT_FOUND; curr = parent[curr])
							result.cycle.push_back(curr);
						result.cycle.push_back(v);
						std::reverse(result.cycle.begin(), result.cycle.end());
						return true;
					}
				}
				return false;
			};

			for (size_t v = 0; v < n && !result.hasCycle; ++v)
				if (!visited[v])
					dfs(v);
		}

		if (result.hasCycle)
			result.graphStatus = GraphAlgorithmStatus::CycleDetected;
		result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
		return result;
	}

	/// Check whether a graph is bipartite and return a two-coloring when possible.
	template<typename V, typename E>
	BipartiteResult BipartiteColoring(const Graph<V, E>& graph)
	{
		auto startTime = std::chrono::steady_clock::now();
		BipartiteResult result;
		result.algorithm_name = "BipartiteColoring";
		size_t n = graph.numVertices();
		result.color.assign(n, -1);

		for (size_t start = 0; start < n; ++start)
		{
			if (result.color[start] != -1)
				continue;
			std::queue<size_t> queue;
			result.color[start] = 0;
			queue.push(start);

			while (!queue.empty())
			{
				size_t u = queue.front();
				queue.pop();
				for (const auto& edge : graph.neighbors(u))
				{
					size_t v = edge.to;
					if (result.color[v] == -1)
					{
						result.color[v] = 1 - result.color[u];
						queue.push(v);
					}
					else if (result.color[v] == result.color[u])
					{
						result.isBipartite = false;
						result.status = AlgorithmStatus::AlgorithmSpecificFailure;
						result.graphStatus = GraphAlgorithmStatus::CycleDetected;
						detail::SetGraphMessage(result.message, result.diagnostics, "Odd cycle prevents bipartite coloring");
						result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
						return result;
					}
				}
			}
		}

		result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
		return result;
	}

	/// Compute transitive closure for a directed graph.
	template<typename V, typename E>
	TransitiveResult TransitiveClosure(const Graph<V, E>& graph)
	{
		auto startTime = std::chrono::steady_clock::now();
		TransitiveResult result;
		result.algorithm_name = "TransitiveClosure";
		size_t n = graph.numVertices();

		if (!graph.isDirected())
		{
			result.status = AlgorithmStatus::InvalidInput;
			result.graphStatus = GraphAlgorithmStatus::GraphTypeMismatch;
			detail::SetGraphMessage(result.message, result.diagnostics, "Transitive closure requires a directed graph");
			result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
			return result;
		}

		result.reachable.assign(n, std::vector<bool>(n, false));
		for (size_t i = 0; i < n; ++i)
			result.reachable[i][i] = true;
		for (const auto& [u, v, weight] : graph.directedEdges())
		{
			(void)weight;
			result.reachable[u][v] = true;
		}

		for (size_t k = 0; k < n; ++k)
			for (size_t i = 0; i < n; ++i)
				if (result.reachable[i][k])
					for (size_t j = 0; j < n; ++j)
						result.reachable[i][j] = result.reachable[i][j] || result.reachable[k][j];

		for (size_t i = 0; i < n; ++i)
			for (size_t j = 0; j < n; ++j)
				if (i != j && result.reachable[i][j])
					result.edges.push_back({ i, j });

		result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
		return result;
	}

	/// Compute transitive reduction edge list for a DAG.
	template<typename V, typename E>
	TransitiveResult TransitiveReductionDAG(const Graph<V, E>& graph)
	{
		auto startTime = std::chrono::steady_clock::now();
		TransitiveResult result;
		result.algorithm_name = "TransitiveReductionDAG";

		if (!graph.isDirected())
		{
			result.status = AlgorithmStatus::InvalidInput;
			result.graphStatus = GraphAlgorithmStatus::GraphTypeMismatch;
			detail::SetGraphMessage(result.message, result.diagnostics, "Transitive reduction requires a directed graph");
			result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
			return result;
		}

		auto topo = TopologicalSort(graph);
		if (!topo.isDAG)
		{
			result.status = topo.status;
			result.graphStatus = topo.graphStatus;
			detail::SetGraphMessage(result.message, result.diagnostics, "Transitive reduction requires a DAG");
			result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
			return result;
		}

		auto closure = TransitiveClosure(graph);
		result.reachable = closure.reachable;
		for (const auto& [u, v, weight] : graph.directedEdges())
		{
			(void)weight;
			bool redundant = false;
			for (size_t mid = 0; mid < graph.numVertices(); ++mid)
			{
				if (mid != u && mid != v && closure.reachable[u][mid] && closure.reachable[mid][v])
				{
					redundant = true;
					break;
				}
			}
			if (!redundant)
				result.edges.push_back({ u, v });
		}

		std::sort(result.edges.begin(), result.edges.end());
		result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
		return result;
	}

	/// Compute longest paths from a source in a DAG.
	template<typename V, typename E>
	LongestPathDAGResult LongestPathDAG(const Graph<V, E>& graph, size_t source)
	{
		auto startTime = std::chrono::steady_clock::now();
		LongestPathDAGResult result;
		result.algorithm_name = "LongestPathDAG";
		size_t n = graph.numVertices();

		if (source >= n)
		{
			result.status = AlgorithmStatus::InvalidInput;
			result.graphStatus = GraphAlgorithmStatus::InvalidVertex;
			detail::SetGraphMessage(result.message, result.diagnostics, "Invalid source vertex");
			result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
			return result;
		}

		auto topo = TopologicalSort(graph);
		if (!topo.isDAG)
		{
			result.status = topo.status;
			result.graphStatus = topo.graphStatus;
			detail::SetGraphMessage(result.message, result.diagnostics, "Longest path requires a DAG");
			result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
			return result;
		}

		result.source = source;
		result.distance.assign(n, -std::numeric_limits<Real>::infinity());
		result.parent.assign(n, GRAPH_NOT_FOUND);
		result.distance[source] = 0;

		for (size_t u : topo.order)
		{
			if (result.distance[u] == -std::numeric_limits<Real>::infinity())
				continue;
			for (const auto& edge : graph.neighbors(u))
			{
				size_t v = edge.to;
				Real candidate = result.distance[u] + static_cast<Real>(edge.weight);
				if (candidate > result.distance[v])
				{
					result.distance[v] = candidate;
					result.parent[v] = u;
				}
			}
		}

		result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
		return result;
	}

}

#endif

///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        GraphFlow.h                                                         ///
///  Description: Maximum-flow, minimum-cut, and matching algorithms                 ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_GRAPH_ALGORITHMS_GRAPHS_GRAPHFLOW_H
#define MML_GRAPH_ALGORITHMS_GRAPHS_GRAPHFLOW_H

#include <mml/algorithms/Graphs/GraphStructure.h>
#include <algorithm>
#include <functional>
#include <limits>
#include <queue>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace MML
{

	///////////////////////////////////////////////////////////////////////////
	///                    FLOWS, CUTS, AND MATCHING                        ///
	///////////////////////////////////////////////////////////////////////////

	namespace detail
	{
		template<typename V, typename E>
		bool InitializeCapacityMatrix(
			const Graph<V, E>& graph,
			std::vector<std::vector<Real>>& capacity,
			std::string& message,
			std::string& diagnostics)
		{
			size_t n = graph.numVertices();
			capacity.assign(n, std::vector<Real>(n, 0));
			for (const auto& [u, v, weight] : graph.directedEdges())
			{
				Real c = static_cast<Real>(weight);
				if (c < 0)
				{
					SetGraphMessage(message, diagnostics, "Flow capacities must be non-negative");
					return false;
				}
				capacity[u][v] += c;
			}
			return true;
		}

		template<typename V, typename E>
		void PopulateMinCutAndFlows(
			const Graph<V, E>& graph,
			size_t source,
			const std::vector<std::vector<Real>>& capacity,
			const std::vector<std::vector<Real>>& residual,
			MaxFlowMinCutResult& result)
		{
			size_t n = graph.numVertices();
			std::vector<bool> reachable(n, false);
			std::queue<size_t> queue;
			reachable[source] = true;
			queue.push(source);

			while (!queue.empty())
			{
				size_t u = queue.front();
				queue.pop();
				for (size_t v = 0; v < n; ++v)
				{
					if (!reachable[v] && residual[u][v] > 0)
					{
						reachable[v] = true;
						queue.push(v);
					}
				}
			}

			for (size_t v = 0; v < n; ++v)
			{
				if (reachable[v])
					result.sourceSide.push_back(v);
				else
					result.sinkSide.push_back(v);
			}

			for (const auto& [u, v, weight] : graph.directedEdges())
			{
				(void)weight;
				Real flow = capacity[u][v] - residual[u][v];
				if (capacity[u][v] > 0)
					result.flowEdges.push_back({ u, v, flow, capacity[u][v] });
				if (reachable[u] && !reachable[v] && capacity[u][v] > 0)
					result.cutEdges.push_back({ u, v, capacity[u][v] });
			}
		}
	}  // namespace detail

	/// Compute max-flow/min-cut using Edmonds-Karp.
	template<typename V, typename E>
	MaxFlowMinCutResult EdmondsKarp(const Graph<V, E>& graph, size_t source, size_t sink)
	{
		auto startTime = std::chrono::steady_clock::now();
		MaxFlowMinCutResult result;
		result.algorithm_name = "EdmondsKarp";
		size_t n = graph.numVertices();

		if (!graph.isDirected())
		{
			result.status = AlgorithmStatus::InvalidInput;
			result.graphStatus = GraphAlgorithmStatus::GraphTypeMismatch;
			detail::SetGraphMessage(result.message, result.diagnostics, "Max-flow requires a directed graph");
			result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
			return result;
		}
		if (source >= n || sink >= n || source == sink)
		{
			result.status = AlgorithmStatus::InvalidInput;
			result.graphStatus = GraphAlgorithmStatus::InvalidVertex;
			detail::SetGraphMessage(result.message, result.diagnostics, "Invalid source or sink vertex");
			result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
			return result;
		}

		std::vector<std::vector<Real>> capacity;
		if (!detail::InitializeCapacityMatrix(graph, capacity, result.message, result.diagnostics))
		{
			result.status = AlgorithmStatus::InvalidInput;
			result.graphStatus = GraphAlgorithmStatus::UnsupportedNegativeWeights;
			result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
			return result;
		}

		auto residual = capacity;
		std::vector<size_t> parent(n, GRAPH_NOT_FOUND);

		auto bfs = [&]() -> Real
		{
			std::fill(parent.begin(), parent.end(), GRAPH_NOT_FOUND);
			parent[source] = source;
			std::queue<std::pair<size_t, Real>> queue;
			queue.push({ source, std::numeric_limits<Real>::infinity() });

			while (!queue.empty())
			{
				auto [u, flow] = queue.front();
				queue.pop();
				for (size_t v = 0; v < n; ++v)
				{
					if (parent[v] == GRAPH_NOT_FOUND && residual[u][v] > 0)
					{
						parent[v] = u;
						Real newFlow = std::min(flow, residual[u][v]);
						if (v == sink)
							return newFlow;
						queue.push({ v, newFlow });
					}
				}
			}
			return Real{0};
		};

		Real augment = 0;
		while ((augment = bfs()) > 0)
		{
			result.maxFlow += augment;
			size_t v = sink;
			while (v != source)
			{
				size_t u = parent[v];
				residual[u][v] -= augment;
				residual[v][u] += augment;
				v = u;
			}
		}

		detail::PopulateMinCutAndFlows(graph, source, capacity, residual, result);
		result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
		return result;
	}

	/// Compute max-flow/min-cut using Dinic's algorithm.
	template<typename V, typename E>
	MaxFlowMinCutResult Dinic(const Graph<V, E>& graph, size_t source, size_t sink)
	{
		auto startTime = std::chrono::steady_clock::now();
		MaxFlowMinCutResult result;
		result.algorithm_name = "Dinic";
		size_t n = graph.numVertices();

		if (!graph.isDirected())
		{
			result.status = AlgorithmStatus::InvalidInput;
			result.graphStatus = GraphAlgorithmStatus::GraphTypeMismatch;
			detail::SetGraphMessage(result.message, result.diagnostics, "Max-flow requires a directed graph");
			result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
			return result;
		}
		if (source >= n || sink >= n || source == sink)
		{
			result.status = AlgorithmStatus::InvalidInput;
			result.graphStatus = GraphAlgorithmStatus::InvalidVertex;
			detail::SetGraphMessage(result.message, result.diagnostics, "Invalid source or sink vertex");
			result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
			return result;
		}

		std::vector<std::vector<Real>> capacity;
		if (!detail::InitializeCapacityMatrix(graph, capacity, result.message, result.diagnostics))
		{
			result.status = AlgorithmStatus::InvalidInput;
			result.graphStatus = GraphAlgorithmStatus::UnsupportedNegativeWeights;
			result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
			return result;
		}

		auto residual = capacity;
		std::vector<std::vector<size_t>> adjacency(n);
		for (size_t u = 0; u < n; ++u)
		{
			for (size_t v = 0; v < n; ++v)
			{
				if (capacity[u][v] > 0 || capacity[v][u] > 0)
					adjacency[u].push_back(v);
			}
		}

		std::vector<int> level(n, -1);
		auto buildLevelGraph = [&]()
		{
			std::fill(level.begin(), level.end(), -1);
			std::queue<size_t> queue;
			level[source] = 0;
			queue.push(source);
			while (!queue.empty())
			{
				size_t u = queue.front();
				queue.pop();
				for (size_t v : adjacency[u])
				{
					if (level[v] < 0 && residual[u][v] > 0)
					{
						level[v] = level[u] + 1;
						queue.push(v);
					}
				}
			}
			return level[sink] >= 0;
		};

		std::vector<size_t> next(n, 0);
		std::function<Real(size_t, Real)> sendFlow = [&](size_t u, Real flow) -> Real
		{
			if (u == sink)
				return flow;
			for (; next[u] < adjacency[u].size(); ++next[u])
			{
				size_t v = adjacency[u][next[u]];
				if (level[v] == level[u] + 1 && residual[u][v] > 0)
				{
					Real pushed = sendFlow(v, std::min(flow, residual[u][v]));
					if (pushed > 0)
					{
						residual[u][v] -= pushed;
						residual[v][u] += pushed;
						return pushed;
					}
				}
			}
			return Real{0};
		};

		while (buildLevelGraph())
		{
			std::fill(next.begin(), next.end(), 0);
			while (true)
			{
				Real pushed = sendFlow(source, std::numeric_limits<Real>::infinity());
				if (pushed <= 0)
					break;
				result.maxFlow += pushed;
			}
		}

		detail::PopulateMinCutAndFlows(graph, source, capacity, residual, result);
		result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
		return result;
	}

	/// Compute maximum cardinality bipartite matching using Hopcroft-Karp and an explicit left partition.
	template<typename V, typename E>
	BipartiteMatchingResult HopcroftKarp(const Graph<V, E>& graph, const std::vector<size_t>& leftPartition)
	{
		auto startTime = std::chrono::steady_clock::now();
		BipartiteMatchingResult result;
		result.algorithm_name = "HopcroftKarp";
		size_t n = graph.numVertices();
		result.mate.assign(n, GRAPH_NOT_FOUND);

		std::vector<bool> isLeft(n, false);
		for (size_t u : leftPartition)
		{
			if (u >= n || isLeft[u])
			{
				result.status = AlgorithmStatus::InvalidInput;
				result.graphStatus = GraphAlgorithmStatus::InvalidVertex;
				detail::SetGraphMessage(result.message, result.diagnostics, "Invalid left partition vertex");
				result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
				return result;
			}
			isLeft[u] = true;
		}

		std::vector<std::vector<size_t>> adjacency(n);
		std::vector<std::set<size_t>> uniqueAdjacency(n);
		for (const auto& [u, v, weight] : graph.directedEdges())
		{
			(void)weight;
			if (isLeft[u] == isLeft[v])
			{
				result.status = AlgorithmStatus::InvalidInput;
				result.graphStatus = GraphAlgorithmStatus::GraphTypeMismatch;
				detail::SetGraphMessage(result.message, result.diagnostics, "Edges must cross the bipartition");
				result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
				return result;
			}
			if (isLeft[u])
				uniqueAdjacency[u].insert(v);
			else
				uniqueAdjacency[v].insert(u);
		}
		for (size_t u = 0; u < n; ++u)
			adjacency[u].assign(uniqueAdjacency[u].begin(), uniqueAdjacency[u].end());

		std::vector<int> dist(n, -1);
		auto bfs = [&]()
		{
			std::queue<size_t> queue;
			bool foundFreeRight = false;
			for (size_t u : leftPartition)
			{
				if (result.mate[u] == GRAPH_NOT_FOUND)
				{
					dist[u] = 0;
					queue.push(u);
				}
				else
					dist[u] = -1;
			}

			while (!queue.empty())
			{
				size_t u = queue.front();
				queue.pop();
				for (size_t v : adjacency[u])
				{
					size_t matchedLeft = result.mate[v];
					if (matchedLeft == GRAPH_NOT_FOUND)
						foundFreeRight = true;
					else if (dist[matchedLeft] < 0)
					{
						dist[matchedLeft] = dist[u] + 1;
						queue.push(matchedLeft);
					}
				}
			}
			return foundFreeRight;
		};

		std::function<bool(size_t)> dfs = [&](size_t u)
		{
			for (size_t v : adjacency[u])
			{
				size_t matchedLeft = result.mate[v];
				if (matchedLeft == GRAPH_NOT_FOUND || (dist[matchedLeft] == dist[u] + 1 && dfs(matchedLeft)))
				{
					result.mate[u] = v;
					result.mate[v] = u;
					return true;
				}
			}
			dist[u] = -1;
			return false;
		};

		while (bfs())
		{
			for (size_t u : leftPartition)
				if (result.mate[u] == GRAPH_NOT_FOUND && dfs(u))
					++result.cardinality;
		}

		for (size_t u : leftPartition)
			if (result.mate[u] != GRAPH_NOT_FOUND)
				result.matching.push_back({ u, result.mate[u] });
		std::sort(result.matching.begin(), result.matching.end());
		result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
		return result;
	}

	/// Compute maximum cardinality bipartite matching by first deriving a two-coloring.
	template<typename V, typename E>
	BipartiteMatchingResult HopcroftKarp(const Graph<V, E>& graph)
	{
		auto coloring = BipartiteColoring(graph);
		if (!coloring.isBipartite)
		{
			BipartiteMatchingResult result;
			result.algorithm_name = "HopcroftKarp";
			result.status = coloring.status;
			result.graphStatus = coloring.graphStatus;
			result.message = coloring.message;
			result.diagnostics = coloring.diagnostics;
			result.elapsed_time_ms = coloring.elapsed_time_ms;
			return result;
		}

		std::vector<size_t> leftPartition;
		for (size_t v = 0; v < coloring.color.size(); ++v)
			if (coloring.color[v] == 0)
				leftPartition.push_back(v);
		return HopcroftKarp(graph, leftPartition);
	}

}

#endif

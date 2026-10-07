///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        GraphSpanningTree.h                                                 ///
///  Description: Minimum-spanning-tree algorithms                                   ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_GRAPH_ALGORITHMS_GRAPHS_GRAPHSPANNINGTREE_H
#define MML_GRAPH_ALGORITHMS_GRAPHS_GRAPHSPANNINGTREE_H

#include <mml/algorithms/Graphs/detail/GraphAlgorithmUtils.h>
#include <algorithm>
#include <functional>
#include <queue>
#include <tuple>
#include <utility>
#include <vector>

namespace MML
{

	///////////////////////////////////////////////////////////////////////////
	///                    MINIMUM SPANNING TREE (MST)                      ///
	///////////////////////////////////////////////////////////////////////////

	namespace detail
	{
		/// Union-Find (Disjoint Set Union) data structure for Kruskal's algorithm
		class UnionFind
		{
		public:
			explicit UnionFind(size_t n) : parent(n), rank(n, 0)
			{
				for (size_t i = 0; i < n; ++i)
					parent[i] = i;
			}

			size_t find(size_t x)
			{
				if (parent[x] != x)
					parent[x] = find(parent[x]);  // Path compression
				return parent[x];
			}

			bool unite(size_t x, size_t y)
			{
				size_t px = find(x);
				size_t py = find(y);
				if (px == py)
					return false;  // Already in same set

				// Union by rank
				if (rank[px] < rank[py])
					std::swap(px, py);
				parent[py] = px;
				if (rank[px] == rank[py])
					++rank[px];
				return true;
			}

		private:
			std::vector<size_t> parent;
			std::vector<size_t> rank;
		};
	}  // namespace detail

	/// Find Minimum Spanning Tree using Kruskal's algorithm
	/// Complexity: O(E log E)
	///
	/// @param graph  The graph (must be connected for a valid MST)
	/// @return MSTResult with MST edges and total weight
	template<typename V, typename E>
	MSTResult Kruskal(const Graph<V, E>& graph)
	{
		auto startTime = std::chrono::steady_clock::now();
		MSTResult result;
		result.algorithm_name = "Kruskal";
		size_t n = graph.numVertices();

		if (n == 0)
		{
			result.isComplete = true;
			result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
			return result;
		}

		std::vector<std::tuple<E, size_t, size_t>> edges;
		for (const auto& [u, v, weight] : graph.edges())
		{
			edges.push_back({ weight, u, v });
		}

		// Sort edges by weight
		std::sort(edges.begin(), edges.end());

		// Kruskal's algorithm
		detail::UnionFind uf(n);

		for (const auto& [weight, u, v] : edges)
		{
			if (uf.unite(u, v))
			{
				result.edges.push_back({u, v, static_cast<Real>(weight)});
				result.totalWeight += static_cast<Real>(weight);

				// MST has exactly V-1 edges
				if (result.edges.size() == n - 1)
					break;
			}
		}

		result.isComplete = (result.edges.size() == n - 1);
		if (!result.isComplete)
		{
			result.status = AlgorithmStatus::AlgorithmSpecificFailure;
			result.graphStatus = GraphAlgorithmStatus::DisconnectedGraph;
			detail::SetGraphMessage(result.message, result.diagnostics, "Graph is not connected - no spanning tree exists");
		}

		result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
		return result;
	}

	/// Alternative MST using Prim's algorithm
	/// Better for dense graphs
	/// Complexity: O((V + E) log V) with priority queue
	template<typename V, typename E>
	MSTResult Prim(const Graph<V, E>& graph, size_t start = 0)
	{
		auto startTime = std::chrono::steady_clock::now();
		MSTResult result;
		result.algorithm_name = "Prim";
		size_t n = graph.numVertices();

		if (n == 0)
		{
			result.isComplete = true;
			result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
			return result;
		}

		if (start >= n)
		{
			result.status = AlgorithmStatus::InvalidInput;
			result.graphStatus = GraphAlgorithmStatus::InvalidVertex;
			detail::SetGraphMessage(result.message, result.diagnostics, "Invalid start vertex");
			result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
			return result;
		}

		std::vector<bool> inMST(n, false);

		// Priority queue: (weight, from, to)
		using PQEntry = std::tuple<E, size_t, size_t>;
		std::priority_queue<PQEntry, std::vector<PQEntry>, std::greater<PQEntry>> pq;

		// Start from the given vertex
		inMST[start] = true;
		for (const auto& edge : graph.neighbors(start))
		{
			pq.push({edge.weight, start, edge.to});
		}

		size_t edgesAdded = 0;

		while (!pq.empty() && edgesAdded < n - 1)
		{
			auto [weight, u, v] = pq.top();
			pq.pop();

			if (inMST[v])
				continue;

			inMST[v] = true;
			result.edges.push_back({u, v, static_cast<Real>(weight)});
			result.totalWeight += static_cast<Real>(weight);
			edgesAdded++;

			for (const auto& edge : graph.neighbors(v))
			{
				if (!inMST[edge.to])
					pq.push({edge.weight, v, edge.to});
			}
		}

		result.isComplete = (result.edges.size() == n - 1);
		if (!result.isComplete)
		{
			result.status = AlgorithmStatus::AlgorithmSpecificFailure;
			result.graphStatus = GraphAlgorithmStatus::DisconnectedGraph;
			detail::SetGraphMessage(result.message, result.diagnostics, "Graph is not connected - no spanning tree exists");
		}

		result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
		return result;
	}

}

#endif

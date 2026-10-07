///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        GraphTraversal.h                                                    ///
///  Description: Graph traversal and connected-components algorithms               ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_GRAPH_ALGORITHMS_GRAPHS_GRAPHTRAVERSAL_H
#define MML_GRAPH_ALGORITHMS_GRAPHS_GRAPHTRAVERSAL_H

#include <mml/algorithms/Graphs/detail/GraphAlgorithmUtils.h>
#include <limits>
#include <queue>
#include <stack>
#include <utility>
#include <vector>

namespace MML
{

	///////////////////////////////////////////////////////////////////////////
	///                    BREADTH-FIRST SEARCH (BFS)                       ///
	///////////////////////////////////////////////////////////////////////////

	/// Perform BFS from start vertex, visiting all reachable vertices
	/// Complexity: O(V + E)
	///
	/// @param graph  The graph to traverse
	/// @param start  Starting vertex index
	/// @return TraversalResult with visit order, parent pointers, and distances (in hops)
	template<typename V, typename E>
	TraversalResult BFS(const Graph<V, E>& graph, size_t start)
	{
		auto startTime = std::chrono::steady_clock::now();
		TraversalResult result;
		result.algorithm_name = "BFS";
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

		std::vector<bool> visited(n, false);
		std::queue<size_t> queue;

		visited[start] = true;
		result.distance[start] = 0;
		queue.push(start);

		while (!queue.empty())
		{
			size_t u = queue.front();
			queue.pop();

			result.visitOrder.push_back(u);
			result.nodesVisited++;

			for (const auto& edge : graph.neighbors(u))
			{
				size_t v = edge.to;
				if (!visited[v])
				{
					visited[v] = true;
					result.parent[v] = u;
					result.distance[v] = result.distance[u] + 1;
					queue.push(v);
				}
			}
		}

		result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
		return result;
	}

	/// BFS with visitor callback - called for each visited vertex
	template<typename V, typename E, typename Visitor>
	void BFSVisit(const Graph<V, E>& graph, size_t start, Visitor visitor)
	{
		size_t n = graph.numVertices();
		if (start >= n) return;

		std::vector<bool> visited(n, false);
		std::queue<size_t> queue;

		visited[start] = true;
		queue.push(start);

		while (!queue.empty())
		{
			size_t u = queue.front();
			queue.pop();

			if (!visitor(u))  // Visitor returns false to stop
				return;

			for (const auto& edge : graph.neighbors(u))
			{
				if (!visited[edge.to])
				{
					visited[edge.to] = true;
					queue.push(edge.to);
				}
			}
		}
	}

	///////////////////////////////////////////////////////////////////////////
	///                    DEPTH-FIRST SEARCH (DFS)                         ///
	///////////////////////////////////////////////////////////////////////////

	/// Perform DFS from start vertex (iterative)
	/// Complexity: O(V + E)
	///
	/// @param graph  The graph to traverse
	/// @param start  Starting vertex index
	/// @return TraversalResult with visit order and parent pointers
	template<typename V, typename E>
	TraversalResult DFS(const Graph<V, E>& graph, size_t start)
	{
		auto startTime = std::chrono::steady_clock::now();
		TraversalResult result;
		result.algorithm_name = "DFS";
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

		std::vector<bool> visited(n, false);
		std::stack<size_t> stack;

		stack.push(start);
		result.distance[start] = 0;

		while (!stack.empty())
		{
			size_t u = stack.top();
			stack.pop();

			if (visited[u])
				continue;

			visited[u] = true;
			result.visitOrder.push_back(u);
			result.nodesVisited++;

			// Push neighbors in reverse order for consistent traversal
			const auto& neighbors = graph.neighbors(u);
			for (auto it = neighbors.rbegin(); it != neighbors.rend(); ++it)
			{
				size_t v = it->to;
				if (!visited[v])
				{
					if (result.parent[v] == GRAPH_NOT_FOUND)
					{
						result.parent[v] = u;
						result.distance[v] = result.distance[u] + 1;
					}
					stack.push(v);
				}
			}
		}

		result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
		return result;
	}

	/// DFS with visitor callback
	template<typename V, typename E, typename Visitor>
	void DFSVisit(const Graph<V, E>& graph, size_t start, Visitor visitor)
	{
		size_t n = graph.numVertices();
		if (start >= n) return;

		std::vector<bool> visited(n, false);
		std::stack<size_t> stack;

		stack.push(start);

		while (!stack.empty())
		{
			size_t u = stack.top();
			stack.pop();

			if (visited[u])
				continue;

			visited[u] = true;

			if (!visitor(u))
				return;

			const auto& neighbors = graph.neighbors(u);
			for (auto it = neighbors.rbegin(); it != neighbors.rend(); ++it)
			{
				if (!visited[it->to])
					stack.push(it->to);
			}
		}
	}

	///////////////////////////////////////////////////////////////////////////
	///                    CONNECTED COMPONENTS                             ///
	///////////////////////////////////////////////////////////////////////////

	/// Find all connected components in an undirected graph
	/// Complexity: O(V + E)
	///
	/// @param graph  The graph (should be undirected for meaningful results)
	/// @return ComponentsResult with component IDs and lists
	template<typename V, typename E>
	ComponentsResult ConnectedComponents(const Graph<V, E>& graph)
	{
		auto startTime = std::chrono::steady_clock::now();
		ComponentsResult result;
		result.algorithm_name = "ConnectedComponents";
		size_t n = graph.numVertices();

		if (n == 0)
		{
			result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
			return result;
		}

		result.componentId.resize(n, GRAPH_NOT_FOUND);

		size_t componentNum = 0;

		for (size_t v = 0; v < n; ++v)
		{
			if (result.componentId[v] != GRAPH_NOT_FOUND)
				continue;  // Already assigned

			// BFS to find all vertices in this component
			std::vector<size_t> component;
			std::queue<size_t> queue;

			queue.push(v);
			result.componentId[v] = componentNum;

			while (!queue.empty())
			{
				size_t u = queue.front();
				queue.pop();
				component.push_back(u);

				for (const auto& edge : graph.neighbors(u))
				{
					if (result.componentId[edge.to] == GRAPH_NOT_FOUND)
					{
						result.componentId[edge.to] = componentNum;
						queue.push(edge.to);
					}
				}
			}

			result.components.push_back(std::move(component));
			++componentNum;
		}

		result.numComponents = componentNum;
		result.elapsed_time_ms = detail::ElapsedMilliseconds(startTime);
		return result;
	}

	/// Check if the graph is connected (single component)
	template<typename V, typename E>
	bool IsConnected(const Graph<V, E>& graph)
	{
		if (graph.numVertices() == 0)
			return true;

		auto result = BFS(graph, 0);
		return result.nodesVisited == graph.numVertices();
	}

}

#endif

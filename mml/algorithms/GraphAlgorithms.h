///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        GraphAlgorithms.h                                                   ///
///  Description: Graph traversal and path-finding algorithms                         ///
///               BFS, DFS, Dijkstra, Connected Components                            ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_GRAPH_ALGORITHMS_H
#define MML_GRAPH_ALGORITHMS_H

#include "base/Graph.h"

#include <queue>
#include <stack>
#include <limits>
#include <functional>

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
		TraversalResult result;
		size_t n = graph.numVertices();

		if (start >= n)
			return result;

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
		TraversalResult result;
		size_t n = graph.numVertices();

		if (start >= n)
			return result;

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
		TraversalResult result;
		size_t n = graph.numVertices();

		if (start >= n)
			return result;

		result.parent.resize(n, GRAPH_NOT_FOUND);
		result.distance.resize(n, std::numeric_limits<Real>::infinity());
		result.distance[start] = 0;

		// Priority queue: (distance, vertex)
		using PQEntry = std::pair<Real, size_t>;
		std::priority_queue<PQEntry, std::vector<PQEntry>, std::greater<PQEntry>> pq;

		pq.push({0, start});

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

		return result;
	}

	/// Find shortest path from start to specific target using Dijkstra
	/// Early termination when target is reached
	template<typename V, typename E>
	PathResult DijkstraPath(const Graph<V, E>& graph, size_t start, size_t target)
	{
		PathResult result;
		size_t n = graph.numVertices();

		if (start >= n || target >= n)
		{
			result.diagnostics = "Invalid start or target vertex";
			return result;
		}

		if (start == target)
		{
			result.found = true;
			result.path.push_back(start);
			result.totalWeight = 0;
			result.nodesExplored = 1;
			return result;
		}

		std::vector<Real> dist(n, std::numeric_limits<Real>::infinity());
		std::vector<size_t> parent(n, GRAPH_NOT_FOUND);
		dist[start] = 0;

		using PQEntry = std::pair<Real, size_t>;
		std::priority_queue<PQEntry, std::vector<PQEntry>, std::greater<PQEntry>> pq;
		pq.push({0, start});

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
		result.diagnostics = "No path exists from start to target";
		return result;
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
		ComponentsResult result;
		size_t n = graph.numVertices();

		if (n == 0)
			return result;

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
		PathResult result;
		size_t n = graph.numVertices();

		if (start >= n || target >= n)
		{
			result.diagnostics = "Invalid start or target vertex";
			return result;
		}

		if (start == target)
		{
			result.found = true;
			result.path.push_back(start);
			result.totalWeight = 0;
			result.nodesExplored = 1;
			return result;
		}

		auto bfsResult = BFS(graph, start);
		result.nodesExplored = bfsResult.nodesVisited;

		if (bfsResult.distance[target] == std::numeric_limits<Real>::infinity())
		{
			result.found = false;
			result.diagnostics = "No path exists from start to target";
			return result;
		}

		result.found = true;
		result.path = ReconstructPath(bfsResult.parent, start, target);
		result.totalWeight = bfsResult.distance[target];

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
		TraversalResult result;
		size_t n = graph.numVertices();

		if (start >= n)
			return result;

		result.parent.resize(n, GRAPH_NOT_FOUND);
		result.distance.resize(n, std::numeric_limits<Real>::infinity());
		result.distance[start] = 0;

		// Collect all edges
		std::vector<std::tuple<size_t, size_t, E>> edges;
		for (size_t u = 0; u < n; ++u)
		{
			for (const auto& edge : graph.neighbors(u))
			{
				edges.push_back({u, edge.to, edge.weight});
			}
		}

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
					result.parent.clear();
					result.distance.clear();
					result.visitOrder.clear();
					result.nodesVisited = 0;
					return result;  // Empty result indicates negative cycle
				}
			}
		}

		// Count reachable nodes
		for (size_t i = 0; i < n; ++i)
		{
			if (result.distance[i] != std::numeric_limits<Real>::infinity())
				result.nodesVisited++;
		}

		return result;
	}

	/// Check if graph contains a negative cycle reachable from start
	template<typename V, typename E>
	bool HasNegativeCycle(const Graph<V, E>& graph, size_t start)
	{
		auto result = BellmanFord(graph, start);
		return result.parent.empty();  // Empty if negative cycle found
	}

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
		TopologicalSortResult result;
		size_t n = graph.numVertices();

		if (n == 0)
		{
			result.isDAG = true;
			return result;
		}

		if (!graph.isDirected())
		{
			result.isDAG = false;
			result.diagnostics = "Topological sort requires a directed graph";
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
			result.isDAG = false;
			result.order.clear();
			result.diagnostics = "Graph contains a cycle";
			return result;
		}

		result.isDAG = true;
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
		MSTResult result;
		size_t n = graph.numVertices();

		if (n == 0)
		{
			result.isComplete = true;
			return result;
		}

		// Collect all edges (for undirected graph, avoid duplicates)
		std::vector<std::tuple<E, size_t, size_t>> edges;
		for (size_t u = 0; u < n; ++u)
		{
			for (const auto& edge : graph.neighbors(u))
			{
				// For undirected graphs, only add edge once (u < v)
				if (graph.isDirected() || u < edge.to)
					edges.push_back({edge.weight, u, edge.to});
			}
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
			result.diagnostics = "Graph is not connected - no spanning tree exists";
		}

		return result;
	}

	/// Alternative MST using Prim's algorithm
	/// Better for dense graphs
	/// Complexity: O((V + E) log V) with priority queue
	template<typename V, typename E>
	MSTResult Prim(const Graph<V, E>& graph, size_t start = 0)
	{
		MSTResult result;
		size_t n = graph.numVertices();

		if (n == 0)
		{
			result.isComplete = true;
			return result;
		}

		if (start >= n)
			start = 0;

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
			result.diagnostics = "Graph is not connected - no spanning tree exists";
		}

		return result;
	}

}  // namespace MML

#endif // MML_GRAPH_ALGORITHMS_H

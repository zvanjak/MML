///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        GraphSpectral.h                                                     ///
///  Description: Spectral graph analysis algorithms                                  ///
///               Algebraic connectivity, spectral partitioning, centrality measures  ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_GRAPH_SPECTRAL_H
#define MML_GRAPH_SPECTRAL_H

#include "base/Graph.h"
#include "GraphAlgorithms.h"
#include "algorithms/EigenSystemSolvers.h"

#include <cmath>
#include <algorithm>
#include <limits>
#include <vector>

namespace MML
{
	///////////////////////////////////////////////////////////////////////////
	///                    SPECTRAL PROPERTIES                              ///
	///////////////////////////////////////////////////////////////////////////

	/// Compute all eigenvalues of the Laplacian matrix using Jacobi iteration.
	/// For connected graphs: lambda_1 = 0, lambda_2 > 0 (algebraic connectivity)
	/// Complexity: O(V^3) for full eigenvalue decomposition
	///
	/// @param graph  The graph
	/// @param tol    Convergence tolerance (default: 1e-10)
	/// @param maxIter Maximum Jacobi iterations (default: 100)
	/// @return Vector of eigenvalues in ascending order
	template<typename V, typename E>
	std::vector<Real> LaplacianEigenvalues(const Graph<V, E>& graph, Real tol = 1e-10, int maxIter = 100)
	{
		size_t n = graph.numVertices();
		if (n == 0)
			return {};

		Matrix<Real> L = graph.toLaplacianMatrix();
		auto result = SymmMatEigenSolverJacobi::Solve(L, tol, maxIter);

		// Convert Vector<Real> to std::vector<Real>
		std::vector<Real> eigenvalues;
		eigenvalues.reserve(n);
		for (size_t i = 0; i < n; ++i)
			eigenvalues.push_back(result.eigenvalues[static_cast<int>(i)]);

		return eigenvalues;
	}

	/// Compute algebraic connectivity (Fiedler value) exactly.
	/// This is the second-smallest eigenvalue of the Laplacian matrix.
	/// Measures how well-connected the graph is.
	/// Complexity: O(V^3) for full eigenvalue decomposition
	///
	/// @param graph  The graph (should be connected for meaningful results)
	/// @param tol    Convergence tolerance (default: 1e-10)
	/// @param maxIter Maximum Jacobi iterations (default: 100)
	/// @return Algebraic connectivity (0 if graph is disconnected or has < 2 vertices)
	template<typename V, typename E>
	Real AlgebraicConnectivity(const Graph<V, E>& graph, Real tol = 1e-10, int maxIter = 100)
	{
		size_t n = graph.numVertices();
		if (n <= 1)
			return 0;

		auto eigenvalues = LaplacianEigenvalues(graph, tol, maxIter);
		
		// Return second smallest eigenvalue (eigenvalues are sorted ascending)
		// The smallest is always ~0 for the constant eigenvector
		return (eigenvalues.size() >= 2) ? eigenvalues[1] : 0;
	}

	/// Compute degree sequence from the Laplacian matrix diagonal
	/// @note Returns vertex degrees (diagonal of Laplacian), NOT eigenvalue approximations.
	///       For actual eigenvalues, use LaplacianEigenvalues() instead.
	/// Complexity: O(V)
	/// @see LaplacianEigenvalues() for exact eigenvalue computation using Jacobi iteration
	///
	/// @param graph  The graph
	/// @return Vector of vertex degrees in ascending order
	template<typename V, typename E>
	std::vector<Real> LaplacianDegreeSequence(const Graph<V, E>& graph)
	{
		size_t n = graph.numVertices();
		if (n == 0)
			return {};

		Matrix<Real> L = graph.toLaplacianMatrix();

		std::vector<Real> eigenvalues;
		eigenvalues.reserve(n);

		// Return diagonal elements as approximation (works for special cases like star graphs)
		// The diagonal elements are the vertex degrees for unweighted graphs
		for (size_t i = 0; i < n; ++i)
		{
			eigenvalues.push_back(L(i, i));
		}

		std::sort(eigenvalues.begin(), eigenvalues.end());
		return eigenvalues;
	}

	/// Compute lower bound on algebraic connectivity using Cheeger-like inequality
	/// @note Returns a lower bound, not the exact Fiedler value.
	///       For exact algebraic connectivity, use AlgebraicConnectivity() instead.
	/// Complexity: O(V)
	/// @see AlgebraicConnectivity() for exact computation using Jacobi iteration
	///
	/// @param graph  The graph (should be connected for meaningful results)
	/// @return Lower bound on algebraic connectivity (0 if graph is disconnected)
	template<typename V, typename E>
	Real AlgebraicConnectivityBound(const Graph<V, E>& graph)
	{
		size_t n = graph.numVertices();
		if (n <= 1)
			return 0;

		// Check if disconnected (return 0)
		auto components = ConnectedComponents(graph);
		if (components.numComponents > 1)
			return 0;

		// Use Cheeger-like bound: lambda_2 >= 2 * (1 - cos(pi/n)) * min_degree
		Real minDegree = std::numeric_limits<Real>::max();
		for (size_t v = 0; v < n; ++v)
		{
			Real d = static_cast<Real>(graph.degree(v));
			if (d < minDegree)
				minDegree = d;
		}

		return minDegree > 0 ? 2.0 * (1.0 - std::cos(Constants::PI / n)) * minDegree : 0;
	}

	/// Check if graph is connected using spectral method
	/// Graph is connected iff algebraic connectivity > 0
	template<typename V, typename E>
	bool IsConnectedSpectral(const Graph<V, E>& graph)
	{
		// Count multiplicity of eigenvalue 0 in Laplacian
		// = number of connected components
		auto components = ConnectedComponents(graph);
		return components.numComponents <= 1;
	}

	///////////////////////////////////////////////////////////////////////////
	///                    CENTRALITY MEASURES                              ///
	///////////////////////////////////////////////////////////////////////////

	/// Compute degree centrality for all vertices
	/// Normalized by (n-1) where n is number of vertices
	///
	/// @param graph  The graph
	/// @return Vector of degree centralities (normalized to [0,1])
	template<typename V, typename E>
	std::vector<Real> DegreeCentrality(const Graph<V, E>& graph)
	{
		size_t n = graph.numVertices();
		if (n <= 1)
			return std::vector<Real>(n, 0);

		std::vector<Real> centrality(n);
		Real maxDegree = static_cast<Real>(n - 1);

		for (size_t v = 0; v < n; ++v)
		{
			centrality[v] = static_cast<Real>(graph.degree(v)) / maxDegree;
		}

		return centrality;
	}

	/// Compute closeness centrality for all vertices
	/// Closeness(v) = (n-1) / sum of shortest distances from v to all other vertices
	/// Uses BFS for unweighted graphs
	///
	/// @param graph  The graph
	/// @return Vector of closeness centralities
	template<typename V, typename E>
	std::vector<Real> ClosenessCentrality(const Graph<V, E>& graph)
	{
		size_t n = graph.numVertices();
		if (n <= 1)
			return std::vector<Real>(n, 0);

		std::vector<Real> centrality(n);

		for (size_t v = 0; v < n; ++v)
		{
			auto bfsResult = BFS(graph, v);
			
			Real totalDist = 0;
			size_t reachable = 0;

			for (size_t u = 0; u < n; ++u)
			{
				if (u != v && bfsResult.distance[u] != std::numeric_limits<Real>::infinity())
				{
					totalDist += bfsResult.distance[u];
					++reachable;
				}
			}

			if (reachable > 0 && totalDist > 0)
			{
				centrality[v] = static_cast<Real>(reachable) / totalDist;
			}
			else
			{
				centrality[v] = 0;
			}
		}

		return centrality;
	}

	/// Compute betweenness centrality for all vertices
	/// Measures how often a vertex lies on shortest paths between other vertices
	/// Uses Brandes' algorithm: O(VE) for unweighted graphs
	///
	/// @param graph  The graph
	/// @param normalized  If true, normalize by 2/((n-1)(n-2))
	/// @return Vector of betweenness centralities
	template<typename V, typename E>
	std::vector<Real> BetweennessCentrality(const Graph<V, E>& graph, bool normalized = true)
	{
		size_t n = graph.numVertices();
		if (n <= 2)
			return std::vector<Real>(n, 0);

		std::vector<Real> betweenness(n, 0);

		// Brandes' algorithm
		for (size_t s = 0; s < n; ++s)
		{
			// BFS from s
			std::vector<size_t> dist(n, std::numeric_limits<size_t>::max());
			std::vector<Real> sigma(n, 0);  // Number of shortest paths through v
			std::vector<std::vector<size_t>> pred(n);  // Predecessors on shortest paths

			dist[s] = 0;
			sigma[s] = 1;

			std::queue<size_t> queue;
			std::stack<size_t> stack;
			queue.push(s);

			while (!queue.empty())
			{
				size_t v = queue.front();
				queue.pop();
				stack.push(v);

				for (const auto& edge : graph.neighbors(v))
				{
					size_t w = edge.to;

					// First time visiting w
					if (dist[w] == std::numeric_limits<size_t>::max())
					{
						dist[w] = dist[v] + 1;
						queue.push(w);
					}

					// Shortest path to w via v
					if (dist[w] == dist[v] + 1)
					{
						sigma[w] += sigma[v];
						pred[w].push_back(v);
					}
				}
			}

			// Back-propagation of dependencies
			std::vector<Real> delta(n, 0);

			while (!stack.empty())
			{
				size_t w = stack.top();
				stack.pop();

				for (size_t v : pred[w])
				{
					delta[v] += (sigma[v] / sigma[w]) * (1.0 + delta[w]);
				}

				if (w != s)
				{
					betweenness[w] += delta[w];
				}
			}
		}

		// For undirected graphs, divide by 2 (each path counted twice)
		if (graph.isUndirected())
		{
			for (size_t v = 0; v < n; ++v)
			{
				betweenness[v] /= 2.0;
			}
		}

		// Normalize if requested
		if (normalized && n > 2)
		{
			Real normFactor = 2.0 / ((n - 1) * (n - 2));
			for (size_t v = 0; v < n; ++v)
			{
				betweenness[v] *= normFactor;
			}
		}

		return betweenness;
	}

	/// Compute eigenvector centrality using power iteration
	/// Vertices with high eigenvector centrality are connected to many
	/// vertices that themselves have high eigenvector centrality
	///
	/// @param graph  The graph
	/// @param maxIterations  Maximum iterations for power method
	/// @param tolerance  Convergence tolerance
	/// @return Vector of eigenvector centralities (L2-normalized)
	template<typename V, typename E>
	std::vector<Real> EigenvectorCentrality(
		const Graph<V, E>& graph,
		size_t maxIterations = 100,
		Real tolerance = 1e-6)
	{
		size_t n = graph.numVertices();
		if (n == 0)
			return {};

		// Initialize with uniform values
		std::vector<Real> centrality(n, 1.0 / std::sqrt(static_cast<Real>(n)));
		std::vector<Real> next(n);

		Matrix<Real> A = graph.toAdjacencyMatrix();

		for (size_t iter = 0; iter < maxIterations; ++iter)
		{
			// Matrix-vector multiply: next = A * centrality
			for (size_t i = 0; i < n; ++i)
			{
				next[i] = 0;
				for (size_t j = 0; j < n; ++j)
				{
					next[i] += A(i, j) * centrality[j];
				}
			}

			// Normalize
			Real norm = 0;
			for (size_t i = 0; i < n; ++i)
			{
				norm += next[i] * next[i];
			}
			norm = std::sqrt(norm);

			if (norm > 0)
			{
				for (size_t i = 0; i < n; ++i)
				{
					next[i] /= norm;
				}
			}

			// Check convergence
			Real diff = 0;
			for (size_t i = 0; i < n; ++i)
			{
				diff += std::abs(next[i] - centrality[i]);
			}

			centrality = next;

			if (diff < tolerance)
				break;
		}

		return centrality;
	}

	/// Compute PageRank centrality
	/// Models a random surfer on the graph
	///
	/// @param graph  The graph
	/// @param dampingFactor  Probability of following a link (typically 0.85)
	/// @param maxIterations  Maximum iterations
	/// @param tolerance  Convergence tolerance
	/// @return Vector of PageRank scores (sum to 1)
	template<typename V, typename E>
	std::vector<Real> PageRank(
		const Graph<V, E>& graph,
		Real dampingFactor = 0.85,
		size_t maxIterations = 100,
		Real tolerance = 1e-6)
	{
		size_t n = graph.numVertices();
		if (n == 0)
			return {};

		std::vector<Real> rank(n, 1.0 / n);
		std::vector<Real> next(n);
		Real teleport = (1.0 - dampingFactor) / n;

		for (size_t iter = 0; iter < maxIterations; ++iter)
		{
			// Initialize with teleport probability
			std::fill(next.begin(), next.end(), teleport);

			// Add contributions from neighbors
			for (size_t v = 0; v < n; ++v)
			{
				size_t outDegree = graph.neighbors(v).size();
				if (outDegree > 0)
				{
					Real contribution = dampingFactor * rank[v] / outDegree;
					for (const auto& edge : graph.neighbors(v))
					{
						next[edge.to] += contribution;
					}
				}
				else
				{
					// Dangling node: distribute evenly
					Real contribution = dampingFactor * rank[v] / n;
					for (size_t u = 0; u < n; ++u)
					{
						next[u] += contribution;
					}
				}
			}

			// Normalize (should sum to 1, but ensure numerical stability)
			Real sum = 0;
			for (size_t i = 0; i < n; ++i)
				sum += next[i];
			for (size_t i = 0; i < n; ++i)
				next[i] /= sum;

			// Check convergence
			Real diff = 0;
			for (size_t i = 0; i < n; ++i)
			{
				diff += std::abs(next[i] - rank[i]);
			}

			rank = next;

			if (diff < tolerance)
				break;
		}

		return rank;
	}

	///////////////////////////////////////////////////////////////////////////
	///                    GRAPH DENSITY AND STRUCTURE                      ///
	///////////////////////////////////////////////////////////////////////////

	/// Compute graph density: ratio of edges to maximum possible edges
	/// Density = 2E / (V(V-1)) for undirected, E / (V(V-1)) for directed
	template<typename V, typename E>
	Real GraphDensity(const Graph<V, E>& graph)
	{
		size_t n = graph.numVertices();
		if (n <= 1)
			return 0;

		size_t maxEdges = n * (n - 1);
		if (graph.isUndirected())
			maxEdges /= 2;

		return static_cast<Real>(graph.numEdges()) / maxEdges;
	}

	/// Compute average clustering coefficient
	/// Measures how much vertices tend to cluster together
	template<typename V, typename E>
	Real AverageClusteringCoefficient(const Graph<V, E>& graph)
	{
		size_t n = graph.numVertices();
		if (n == 0)
			return 0;

		Real totalCC = 0;
		size_t validVertices = 0;

		for (size_t v = 0; v < n; ++v)
		{
			const auto& neighbors = graph.neighbors(v);
			size_t k = neighbors.size();

			if (k < 2)
				continue;  // Cannot form triangles

			// Count edges between neighbors
			size_t triangles = 0;
			for (size_t i = 0; i < k; ++i)
			{
				for (size_t j = i + 1; j < k; ++j)
				{
					if (graph.hasEdge(neighbors[i].to, neighbors[j].to))
						++triangles;
				}
			}

			// Local clustering coefficient
			Real maxTriangles = static_cast<Real>(k * (k - 1)) / 2;
			totalCC += static_cast<Real>(triangles) / maxTriangles;
			++validVertices;
		}

		return validVertices > 0 ? totalCC / validVertices : 0;
	}

	/// Compute graph diameter (longest shortest path)
	/// Returns infinity for disconnected graphs
	template<typename V, typename E>
	Real GraphDiameter(const Graph<V, E>& graph)
	{
		size_t n = graph.numVertices();
		if (n <= 1)
			return 0;

		Real diameter = 0;

		for (size_t v = 0; v < n; ++v)
		{
			auto bfsResult = BFS(graph, v);

			for (size_t u = 0; u < n; ++u)
			{
				if (bfsResult.distance[u] > diameter)
					diameter = bfsResult.distance[u];
			}
		}

		return diameter;
	}

	/// Compute average path length
	/// Average of all shortest paths between all pairs of vertices
	template<typename V, typename E>
	Real AveragePathLength(const Graph<V, E>& graph)
	{
		size_t n = graph.numVertices();
		if (n <= 1)
			return 0;

		Real totalLength = 0;
		size_t pathCount = 0;

		for (size_t v = 0; v < n; ++v)
		{
			auto bfsResult = BFS(graph, v);

			for (size_t u = v + 1; u < n; ++u)
			{
				if (bfsResult.distance[u] != std::numeric_limits<Real>::infinity())
				{
					totalLength += bfsResult.distance[u];
					++pathCount;
				}
			}
		}

		return pathCount > 0 ? totalLength / pathCount : 0;
	}

}  // namespace MML

#endif // MML_GRAPH_SPECTRAL_H

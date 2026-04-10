///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Graph.h                                                             ///
///  Description: Graph data structure with adjacency list representation             ///
///               Supports directed/undirected, weighted/unweighted graphs            ///
///               Provides matrix conversions for spectral analysis                   ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_GRAPH_H
#define MML_GRAPH_H

#include "MMLBase.h"
#include "MMLExceptions.h"

#include "base/Vector/Vector.h"
#include "base/Matrix/Matrix.h"

#include <vector>
#include <string>
#include <limits>
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace MML
{
	/// @brief Sentinel value for "no parent" / "not found" in graph algorithms
	static constexpr size_t GRAPH_NOT_FOUND = static_cast<size_t>(-1);

	///////////////////////////////////////////////////////////////////////////
	///                         RESULT STRUCTS                              ///
	///////////////////////////////////////////////////////////////////////////

	/// Result of a path-finding algorithm (Dijkstra, BellmanFord, etc.)
	struct PathResult
	{
		bool        found;           ///< Whether a path was found
		std::vector<size_t> path;    ///< Sequence of vertex indices from start to end
		Real        totalWeight;     ///< Total weight/distance of the path
		size_t      nodesExplored;   ///< Number of nodes visited during search
		std::string diagnostics;     ///< Additional diagnostic information

		PathResult() : found(false), totalWeight(0), nodesExplored(0) {}
	};

	/// Result of a traversal algorithm (BFS, DFS)
	struct TraversalResult
	{
		std::vector<size_t> visitOrder;   ///< Order in which vertices were visited
		std::vector<size_t> parent;       ///< Parent of each vertex in traversal tree (GRAPH_NOT_FOUND if root/unvisited)
		std::vector<Real>   distance;     ///< Distance from start (BFS: hops, Dijkstra: weighted)
		size_t              nodesVisited; ///< Total nodes visited

		TraversalResult() : nodesVisited(0) {}
	};

	/// Result of connected components analysis
	struct ComponentsResult
	{
		size_t                          numComponents;  ///< Number of connected components
		std::vector<size_t>             componentId;    ///< componentId[v] = which component vertex v belongs to
		std::vector<std::vector<size_t>> components;    ///< List of vertices in each component

		ComponentsResult() : numComponents(0) {}
	};

	/// Result of minimum spanning tree algorithm
	struct MSTResult
	{
		bool                                         isComplete;   ///< True if MST spans all vertices (graph is connected)
		std::vector<std::tuple<size_t, size_t, Real>> edges;       ///< MST edges: (from, to, weight)
		Real                                         totalWeight;  ///< Sum of edge weights in MST
		std::string                                  diagnostics;

		MSTResult() : isComplete(false), totalWeight(0) {}
	};

	/// Result of topological sort
	struct TopologicalSortResult
	{
		bool                 isDAG;       ///< True if graph is a DAG (no cycles)
		std::vector<size_t>  order;       ///< Topological order of vertices
		std::string          diagnostics; ///< Error message if cycle detected

		TopologicalSortResult() : isDAG(false) {}
	};

	///////////////////////////////////////////////////////////////////////////
	///                         EDGE STRUCTURE                              ///
	///////////////////////////////////////////////////////////////////////////

	/// Represents an edge with destination vertex and weight
	template<typename E = Real>
	struct Edge
	{
		size_t to;      ///< Destination vertex index
		E      weight;  ///< Edge weight

		Edge(size_t dest, E w = E{1}) : to(dest), weight(w) {}

		bool operator<(const Edge& other) const { return weight < other.weight; }
	};

	///////////////////////////////////////////////////////////////////////////
	///                         GRAPH CLASS                                 ///
	///////////////////////////////////////////////////////////////////////////

	/// Graph data structure using adjacency list representation
	/// 
	/// Template parameters:
	/// - V: Vertex data type (default: int for simple indexed graphs)
	/// - E: Edge weight type (default: Real for weighted graphs)
	///
	/// Supports:
	/// - Directed and undirected graphs
	/// - Weighted and unweighted edges
	/// - Self-loops (optional)
	/// - Conversion to adjacency/Laplacian matrices for spectral analysis
	///
	/// Example usage:
	/// @code
	///   Graph<> g(Graph<>::Type::Undirected);
	///   g.addVertex(0);
	///   g.addVertex(1);
	///   g.addVertex(2);
	///   g.addEdge(0, 1, 1.5);
	///   g.addEdge(1, 2, 2.0);
	///   
	///   Matrix<Real> adj = g.toAdjacencyMatrix();
	///   Matrix<Real> lap = g.toLaplacianMatrix();
	/// @endcode
	template<typename V = int, typename E = Real>
	class Graph
	{
	public:
		/// Graph type enumeration
		enum class Type { Directed, Undirected };

	private:
		Type                                   _type;
		std::vector<V>                         _vertexData;
		std::vector<std::vector<Edge<E>>>      _adjList;
		size_t                                 _numEdges;
		bool                                   _allowSelfLoops;

	public:
		//////////////////////////////////////////////////////////////////////////
		///                      CONSTRUCTORS                                  ///
		//////////////////////////////////////////////////////////////////////////

		/// Create empty graph
		explicit Graph(Type type = Type::Undirected, bool allowSelfLoops = false)
			: _type(type), _numEdges(0), _allowSelfLoops(allowSelfLoops)
		{}

		/// Create graph with pre-allocated vertices
		Graph(size_t numVertices, Type type = Type::Undirected, bool allowSelfLoops = false)
			: _type(type), _numEdges(0), _allowSelfLoops(allowSelfLoops)
		{
			_vertexData.resize(numVertices);
			_adjList.resize(numVertices);
			for (size_t i = 0; i < numVertices; ++i)
				_vertexData[i] = static_cast<V>(i);
		}

		//////////////////////////////////////////////////////////////////////////
		///                      VERTEX OPERATIONS                             ///
		//////////////////////////////////////////////////////////////////////////

		/// Add a vertex with optional data, returns vertex index
		size_t addVertex(const V& data = V{})
		{
			size_t idx = _vertexData.size();
			_vertexData.push_back(data);
			_adjList.push_back({});
			return idx;
		}

		/// Get vertex data (const)
		const V& vertexData(size_t v) const
		{
			if (v >= _vertexData.size())
				throw std::out_of_range("Graph::vertexData - vertex index out of range");
			return _vertexData[v];
		}

		/// Get vertex data (mutable)
		V& vertexData(size_t v)
		{
			if (v >= _vertexData.size())
				throw std::out_of_range("Graph::vertexData - vertex index out of range");
			return _vertexData[v];
		}

		/// Number of vertices
		size_t numVertices() const noexcept { return _vertexData.size(); }

		//////////////////////////////////////////////////////////////////////////
		///                      EDGE OPERATIONS                               ///
		//////////////////////////////////////////////////////////////////////////

		/// Add an edge from vertex 'from' to vertex 'to' with given weight
		/// For undirected graphs, adds edge in both directions
		void addEdge(size_t from, size_t to, const E& weight = E{1})
		{
			if (from >= _vertexData.size() || to >= _vertexData.size())
				throw std::out_of_range("Graph::addEdge - vertex index out of range");

			if (!_allowSelfLoops && from == to)
				throw std::invalid_argument("Graph::addEdge - self-loops not allowed");

			// Check if edge already exists
			for (const auto& edge : _adjList[from])
			{
				if (edge.to == to)
					throw std::invalid_argument("Graph::addEdge - edge already exists");
			}

			_adjList[from].push_back(Edge<E>(to, weight));
			++_numEdges;

			if (_type == Type::Undirected && from != to)
			{
				_adjList[to].push_back(Edge<E>(from, weight));
				// Note: we count undirected edge as 1 edge, not 2
			}
		}

		/// Check if edge exists from 'from' to 'to'
		bool hasEdge(size_t from, size_t to) const
		{
			if (from >= _vertexData.size() || to >= _vertexData.size())
				return false;

			for (const auto& edge : _adjList[from])
			{
				if (edge.to == to)
					return true;
			}
			return false;
		}

		/// Get edge weight (throws if edge doesn't exist)
		E edgeWeight(size_t from, size_t to) const
		{
			if (from >= _vertexData.size() || to >= _vertexData.size())
				throw std::out_of_range("Graph::edgeWeight - vertex index out of range");

			for (const auto& edge : _adjList[from])
			{
				if (edge.to == to)
					return edge.weight;
			}
			throw std::invalid_argument("Graph::edgeWeight - edge does not exist");
		}

		/// Number of edges (for undirected graphs, each edge counted once)
		size_t numEdges() const { return _numEdges; }

		/// Remove edge from 'from' to 'to'
		bool removeEdge(size_t from, size_t to)
		{
			if (from >= _vertexData.size() || to >= _vertexData.size())
				return false;

			auto& edges = _adjList[from];
			auto it = std::find_if(edges.begin(), edges.end(),
				[to](const Edge<E>& e) { return e.to == to; });

			if (it == edges.end())
				return false;

			edges.erase(it);

			if (_type == Type::Undirected && from != to)
			{
				auto& reverseEdges = _adjList[to];
				auto rit = std::find_if(reverseEdges.begin(), reverseEdges.end(),
					[from](const Edge<E>& e) { return e.to == from; });
				if (rit != reverseEdges.end())
					reverseEdges.erase(rit);
			}

			--_numEdges;

			return true;
		}

		//////////////////////////////////////////////////////////////////////////
		///                      ADJACENCY ACCESS                              ///
		//////////////////////////////////////////////////////////////////////////

		/// Get all edges from vertex v
		const std::vector<Edge<E>>& neighbors(size_t v) const
		{
			if (v >= _adjList.size())
				throw std::out_of_range("Graph::neighbors - vertex index out of range");
			return _adjList[v];
		}

		/// Degree of vertex (number of edges)
		size_t degree(size_t v) const
		{
			if (v >= _adjList.size())
				throw std::out_of_range("Graph::degree - vertex index out of range");
			return _adjList[v].size();
		}

		/// Out-degree for directed graphs
		size_t outDegree(size_t v) const { return degree(v); }

		/// In-degree for directed graphs (O(V+E) - expensive!)
		size_t inDegree(size_t v) const
		{
			if (v >= _adjList.size())
				throw std::out_of_range("Graph::inDegree - vertex index out of range");

			if (_type == Type::Undirected)
				return degree(v);

			size_t count = 0;
			for (size_t i = 0; i < _adjList.size(); ++i)
			{
				for (const auto& edge : _adjList[i])
				{
					if (edge.to == v)
						++count;
				}
			}
			return count;
		}

		//////////////////////////////////////////////////////////////////////////
		///                      GRAPH PROPERTIES                              ///
		//////////////////////////////////////////////////////////////////////////

		/// Get graph type
		Type type() const { return _type; }

		/// Is this a directed graph?
		bool isDirected() const { return _type == Type::Directed; }

		/// Is this an undirected graph?
		bool isUndirected() const { return _type == Type::Undirected; }

		/// Check if graph has any edges with non-unit weights
		bool isWeighted() const
		{
			for (const auto& edges : _adjList)
			{
				for (const auto& edge : edges)
				{
					if (std::abs(static_cast<Real>(edge.weight) - Real{1}) > Real{1e-10})
						return true;
				}
			}
			return false;
		}

		/// Check if graph is empty (no vertices)
		bool isEmpty() const noexcept { return _vertexData.empty(); }

		/// Are self-loops allowed?
		bool allowsSelfLoops() const noexcept { return _allowSelfLoops; }

		//////////////////////////////////////////////////////////////////////////
		///                      MATRIX REPRESENTATIONS                        ///
		//////////////////////////////////////////////////////////////////////////

		/// Convert to adjacency matrix A where A[i][j] = weight of edge i->j
		Matrix<E> toAdjacencyMatrix() const
		{
			size_t n = _vertexData.size();
			Matrix<E> A(n, n);

			// Initialize to zero
			for (size_t i = 0; i < n; ++i)
				for (size_t j = 0; j < n; ++j)
					A(i, j) = E{0};

			// Fill in edges
			for (size_t i = 0; i < n; ++i)
			{
				for (const auto& edge : _adjList[i])
				{
					A(i, edge.to) = edge.weight;
				}
			}

			return A;
		}

		/// Convert to degree matrix D where D[i][i] = degree of vertex i
		Matrix<E> toDegreeMatrix() const
		{
			size_t n = _vertexData.size();
			Matrix<E> D(n, n);

			// Initialize to zero
			for (size_t i = 0; i < n; ++i)
				for (size_t j = 0; j < n; ++j)
					D(i, j) = E{0};

			// Set diagonal to degrees
			for (size_t i = 0; i < n; ++i)
			{
				E deg = E{0};
				for (const auto& edge : _adjList[i])
					deg += edge.weight;  // Weighted degree
				D(i, i) = deg;
			}

			return D;
		}

		/// Convert to Laplacian matrix L = D - A
		/// For spectral graph theory applications
		Matrix<E> toLaplacianMatrix() const
		{
			size_t n = _vertexData.size();
			Matrix<E> L(n, n);

			// Initialize to zero
			for (size_t i = 0; i < n; ++i)
				for (size_t j = 0; j < n; ++j)
					L(i, j) = E{0};

			for (size_t i = 0; i < n; ++i)
			{
				E deg = E{0};
				for (const auto& edge : _adjList[i])
				{
					L(i, edge.to) = -edge.weight;  // Off-diagonal: -weight
					deg += edge.weight;
				}
				L(i, i) = deg;  // Diagonal: sum of weights
			}

			return L;
		}

		/// Convert to normalized Laplacian L_norm = I - D^(-1/2) A D^(-1/2)
		/// Useful for spectral clustering
		Matrix<Real> toNormalizedLaplacian() const
		{
			size_t n = _vertexData.size();
			Matrix<Real> L(n, n);

			// Compute degrees
			std::vector<Real> deg(n, 0);
			for (size_t i = 0; i < n; ++i)
			{
				for (const auto& edge : _adjList[i])
					deg[i] += static_cast<Real>(edge.weight);
			}

			// Compute D^(-1/2)
			std::vector<Real> invSqrtDeg(n, 0);
			for (size_t i = 0; i < n; ++i)
			{
				if (deg[i] > 0)
					invSqrtDeg[i] = 1.0 / std::sqrt(deg[i]);
			}

			// Initialize to identity
			for (size_t i = 0; i < n; ++i)
			{
				for (size_t j = 0; j < n; ++j)
					L(i, j) = (i == j) ? Real{1} : Real{0};
			}

			// Subtract D^(-1/2) A D^(-1/2)
			for (size_t i = 0; i < n; ++i)
			{
				for (const auto& edge : _adjList[i])
				{
					size_t j = edge.to;
					Real val = invSqrtDeg[i] * static_cast<Real>(edge.weight) * invSqrtDeg[j];
					L(i, j) -= val;
				}
			}

			return L;
		}

		//////////////////////////////////////////////////////////////////////////
		///                      UTILITY METHODS                               ///
		//////////////////////////////////////////////////////////////////////////

		/// Clear all vertices and edges
		void clear()
		{
			_vertexData.clear();
			_adjList.clear();
			_numEdges = 0;
		}

		/// Reserve space for vertices
		void reserve(size_t numVertices)
		{
			_vertexData.reserve(numVertices);
			_adjList.reserve(numVertices);
		}

		/// Get string representation for debugging
		std::string toString() const
		{
			std::string s = isDirected() ? "Directed Graph" : "Undirected Graph";
			s += " (V=" + std::to_string(numVertices()) + ", E=" + std::to_string(numEdges()) + ")\n";

			for (size_t i = 0; i < _adjList.size(); ++i)
			{
				s += "  " + std::to_string(i) + " -> ";
				for (size_t j = 0; j < _adjList[i].size(); ++j)
				{
					if (j > 0) s += ", ";
					s += std::to_string(_adjList[i][j].to);
					s += "(" + std::to_string(_adjList[i][j].weight) + ")";
				}
				s += "\n";
			}
			return s;
		}
	};

	///////////////////////////////////////////////////////////////////////////
	///                      GRAPH FACTORY FUNCTIONS                        ///
	///////////////////////////////////////////////////////////////////////////

	/// Create a path graph P_n: 0 -- 1 -- 2 -- ... -- (n-1)
	template<typename V = int, typename E = Real>
	Graph<V, E> createPathGraph(size_t n)
	{
		Graph<V, E> g(n, Graph<V, E>::Type::Undirected);
		for (size_t i = 0; i < n - 1; ++i)
			g.addEdge(i, i + 1, E{1});
		return g;
	}

	/// Create a cycle graph C_n: 0 -- 1 -- 2 -- ... -- (n-1) -- 0
	template<typename V = int, typename E = Real>
	Graph<V, E> createCycleGraph(size_t n)
	{
		Graph<V, E> g(n, Graph<V, E>::Type::Undirected);
		for (size_t i = 0; i < n - 1; ++i)
			g.addEdge(i, i + 1, E{1});
		if (n > 2)
			g.addEdge(n - 1, 0, E{1});
		return g;
	}

	/// Create a complete graph K_n (all vertices connected)
	template<typename V = int, typename E = Real>
	Graph<V, E> createCompleteGraph(size_t n)
	{
		Graph<V, E> g(n, Graph<V, E>::Type::Undirected);
		for (size_t i = 0; i < n; ++i)
			for (size_t j = i + 1; j < n; ++j)
				g.addEdge(i, j, E{1});
		return g;
	}

	/// Create a star graph S_n (vertex 0 connected to all others)
	template<typename V = int, typename E = Real>
	Graph<V, E> createStarGraph(size_t n)
	{
		Graph<V, E> g(n, Graph<V, E>::Type::Undirected);
		for (size_t i = 1; i < n; ++i)
			g.addEdge(0, i, E{1});
		return g;
	}

	/// Create a grid graph (rows x cols)
	template<typename V = int, typename E = Real>
	Graph<V, E> createGridGraph(size_t rows, size_t cols)
	{
		Graph<V, E> g(rows * cols, Graph<V, E>::Type::Undirected);

		auto idx = [cols](size_t r, size_t c) { return r * cols + c; };

		for (size_t r = 0; r < rows; ++r)
		{
			for (size_t c = 0; c < cols; ++c)
			{
				if (c + 1 < cols)
					g.addEdge(idx(r, c), idx(r, c + 1), E{1});
				if (r + 1 < rows)
					g.addEdge(idx(r, c), idx(r + 1, c), E{1});
			}
		}
		return g;
	}

}  // namespace MML

#endif // MML_GRAPH_H

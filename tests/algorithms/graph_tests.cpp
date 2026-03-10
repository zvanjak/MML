///////////////////////////////////////////////////////////////////////////////////////////
// MML Graph Tests
// Tests for Graph data structure (Graph.h)
///////////////////////////////////////////////////////////////////////////////////////////

#include "../TestPrecision.h"

#include "base/Graph.h"
#include "algorithms/GraphAlgorithms.h"
#include "algorithms/GraphSpectral.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace MML;

namespace MML::Tests::Algorithms::GraphTests
{

///////////////////////////////////////////////////////////////////////////////////////////
///                         GRAPH CONSTRUCTION TESTS                                    ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Graph - default construction creates empty undirected graph", "[Graph][construction]")
{
	Graph<> g;

	REQUIRE(g.numVertices() == 0);
	REQUIRE(g.numEdges() == 0);
	REQUIRE(g.isUndirected());
	REQUIRE_FALSE(g.isDirected());
	REQUIRE(g.isEmpty());
}

TEST_CASE("Graph - construction with size creates vertices", "[Graph][construction]")
{
	Graph<> g(5);

	REQUIRE(g.numVertices() == 5);
	REQUIRE(g.numEdges() == 0);
	REQUIRE(g.isUndirected());

	// Vertices should be initialized with their indices
	for (size_t i = 0; i < 5; ++i)
		REQUIRE(g.vertexData(i) == static_cast<int>(i));
}

TEST_CASE("Graph - directed graph construction", "[Graph][construction]")
{
	Graph<> g(3, Graph<>::Type::Directed);

	REQUIRE(g.numVertices() == 3);
	REQUIRE(g.isDirected());
	REQUIRE_FALSE(g.isUndirected());
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         VERTEX OPERATIONS TESTS                                     ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Graph - addVertex adds vertices and returns index", "[Graph][vertex]")
{
	Graph<std::string> g;

	size_t idx0 = g.addVertex("A");
	size_t idx1 = g.addVertex("B");
	size_t idx2 = g.addVertex("C");

	REQUIRE(idx0 == 0);
	REQUIRE(idx1 == 1);
	REQUIRE(idx2 == 2);
	REQUIRE(g.numVertices() == 3);

	REQUIRE(g.vertexData(0) == "A");
	REQUIRE(g.vertexData(1) == "B");
	REQUIRE(g.vertexData(2) == "C");
}

TEST_CASE("Graph - vertex data can be modified", "[Graph][vertex]")
{
	Graph<int> g(3);

	g.vertexData(0) = 100;
	g.vertexData(1) = 200;
	g.vertexData(2) = 300;

	REQUIRE(g.vertexData(0) == 100);
	REQUIRE(g.vertexData(1) == 200);
	REQUIRE(g.vertexData(2) == 300);
}

TEST_CASE("Graph - vertex access out of range throws", "[Graph][vertex]")
{
	Graph<> g(3);

	REQUIRE_THROWS_AS(g.vertexData(3), std::out_of_range);
	REQUIRE_THROWS_AS(g.vertexData(100), std::out_of_range);
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         EDGE OPERATIONS TESTS                                       ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Graph - addEdge for undirected graph", "[Graph][edge]")
{
	Graph<> g(4);

	g.addEdge(0, 1, 1.0);
	g.addEdge(1, 2, 2.0);
	g.addEdge(2, 3, 3.0);

	REQUIRE(g.numEdges() == 3);

	// Check edges exist in both directions
	REQUIRE(g.hasEdge(0, 1));
	REQUIRE(g.hasEdge(1, 0));
	REQUIRE(g.hasEdge(1, 2));
	REQUIRE(g.hasEdge(2, 1));
	REQUIRE(g.hasEdge(2, 3));
	REQUIRE(g.hasEdge(3, 2));

	// Check non-existent edges
	REQUIRE_FALSE(g.hasEdge(0, 2));
	REQUIRE_FALSE(g.hasEdge(0, 3));
}

TEST_CASE("Graph - addEdge for directed graph", "[Graph][edge]")
{
	Graph<> g(4, Graph<>::Type::Directed);

	g.addEdge(0, 1, 1.0);
	g.addEdge(1, 2, 2.0);

	REQUIRE(g.numEdges() == 2);

	// Check directed edges
	REQUIRE(g.hasEdge(0, 1));
	REQUIRE_FALSE(g.hasEdge(1, 0));  // No reverse edge
	REQUIRE(g.hasEdge(1, 2));
	REQUIRE_FALSE(g.hasEdge(2, 1));
}

TEST_CASE("Graph - edge weights are stored correctly", "[Graph][edge]")
{
	Graph<> g(3);

	g.addEdge(0, 1, 1.5);
	g.addEdge(1, 2, 2.7);

	REQUIRE_THAT(g.edgeWeight(0, 1), Catch::Matchers::WithinRel(1.5, 1e-10));
	REQUIRE_THAT(g.edgeWeight(1, 0), Catch::Matchers::WithinRel(1.5, 1e-10));
	REQUIRE_THAT(g.edgeWeight(1, 2), Catch::Matchers::WithinRel(2.7, 1e-10));
}

TEST_CASE("Graph - addEdge throws for duplicate edge", "[Graph][edge]")
{
	Graph<> g(3);

	g.addEdge(0, 1);
	REQUIRE_THROWS_AS(g.addEdge(0, 1), std::invalid_argument);
}

TEST_CASE("Graph - addEdge throws for self-loop when not allowed", "[Graph][edge]")
{
	Graph<> g(3);  // Self-loops not allowed by default

	REQUIRE_THROWS_AS(g.addEdge(0, 0), std::invalid_argument);
}

TEST_CASE("Graph - self-loops allowed when enabled", "[Graph][edge]")
{
	Graph<> g(3, Graph<>::Type::Undirected, true);  // Allow self-loops

	g.addEdge(0, 0, 5.0);

	REQUIRE(g.hasEdge(0, 0));
	REQUIRE_THAT(g.edgeWeight(0, 0), Catch::Matchers::WithinRel(5.0, 1e-10));
}

TEST_CASE("Graph - removeEdge removes edges", "[Graph][edge]")
{
	Graph<> g(3);

	g.addEdge(0, 1);
	g.addEdge(1, 2);

	REQUIRE(g.numEdges() == 2);
	REQUIRE(g.hasEdge(0, 1));

	bool removed = g.removeEdge(0, 1);

	REQUIRE(removed);
	REQUIRE(g.numEdges() == 1);
	REQUIRE_FALSE(g.hasEdge(0, 1));
	REQUIRE_FALSE(g.hasEdge(1, 0));  // Both directions removed for undirected
}

TEST_CASE("Graph - removeEdge returns false for non-existent edge", "[Graph][edge]")
{
	Graph<> g(3);

	bool removed = g.removeEdge(0, 1);

	REQUIRE_FALSE(removed);
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         DEGREE TESTS                                                ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Graph - degree returns correct values for undirected graph", "[Graph][degree]")
{
	Graph<> g(4);

	g.addEdge(0, 1);
	g.addEdge(0, 2);
	g.addEdge(0, 3);
	g.addEdge(1, 2);

	REQUIRE(g.degree(0) == 3);  // Connected to 1, 2, 3
	REQUIRE(g.degree(1) == 2);  // Connected to 0, 2
	REQUIRE(g.degree(2) == 2);  // Connected to 0, 1
	REQUIRE(g.degree(3) == 1);  // Connected to 0
}

TEST_CASE("Graph - in/out degree for directed graph", "[Graph][degree]")
{
	Graph<> g(4, Graph<>::Type::Directed);

	g.addEdge(0, 1);
	g.addEdge(0, 2);
	g.addEdge(1, 2);
	g.addEdge(3, 2);

	REQUIRE(g.outDegree(0) == 2);
	REQUIRE(g.outDegree(1) == 1);
	REQUIRE(g.outDegree(2) == 0);
	REQUIRE(g.outDegree(3) == 1);

	REQUIRE(g.inDegree(0) == 0);
	REQUIRE(g.inDegree(1) == 1);
	REQUIRE(g.inDegree(2) == 3);
	REQUIRE(g.inDegree(3) == 0);
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         NEIGHBORS ACCESS TESTS                                      ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Graph - neighbors returns adjacent vertices", "[Graph][neighbors]")
{
	Graph<> g(4);

	g.addEdge(0, 1);
	g.addEdge(0, 2);
	g.addEdge(0, 3);

	const auto& neighbors = g.neighbors(0);

	REQUIRE(neighbors.size() == 3);

	std::vector<size_t> neighborIds;
	for (const auto& edge : neighbors)
		neighborIds.push_back(edge.to);

	std::sort(neighborIds.begin(), neighborIds.end());
	REQUIRE(neighborIds == std::vector<size_t>{1, 2, 3});
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         MATRIX CONVERSION TESTS                                     ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Graph - toAdjacencyMatrix for undirected graph", "[Graph][matrix]")
{
	Graph<> g(3);

	g.addEdge(0, 1, 1.0);
	g.addEdge(1, 2, 2.0);

	Matrix<Real> A = g.toAdjacencyMatrix();

	REQUIRE(A.rows() == 3);
	REQUIRE(A.cols() == 3);

	// Check symmetric structure
	REQUIRE_THAT(A(0, 1), Catch::Matchers::WithinRel(1.0, 1e-10));
	REQUIRE_THAT(A(1, 0), Catch::Matchers::WithinRel(1.0, 1e-10));
	REQUIRE_THAT(A(1, 2), Catch::Matchers::WithinRel(2.0, 1e-10));
	REQUIRE_THAT(A(2, 1), Catch::Matchers::WithinRel(2.0, 1e-10));

	// Check zeros
	REQUIRE_THAT(A(0, 0), Catch::Matchers::WithinAbs(0.0, 1e-10));
	REQUIRE_THAT(A(0, 2), Catch::Matchers::WithinAbs(0.0, 1e-10));
	REQUIRE_THAT(A(2, 0), Catch::Matchers::WithinAbs(0.0, 1e-10));
}

TEST_CASE("Graph - toLaplacianMatrix", "[Graph][matrix]")
{
	// Path graph: 0 -- 1 -- 2
	Graph<> g(3);
	g.addEdge(0, 1);
	g.addEdge(1, 2);

	Matrix<Real> L = g.toLaplacianMatrix();

	// Laplacian of path graph P_3:
	// [ 1  -1   0 ]
	// [-1   2  -1 ]
	// [ 0  -1   1 ]

	REQUIRE_THAT(L(0, 0), Catch::Matchers::WithinRel(1.0, 1e-10));
	REQUIRE_THAT(L(0, 1), Catch::Matchers::WithinRel(-1.0, 1e-10));
	REQUIRE_THAT(L(0, 2), Catch::Matchers::WithinAbs(0.0, 1e-10));

	REQUIRE_THAT(L(1, 0), Catch::Matchers::WithinRel(-1.0, 1e-10));
	REQUIRE_THAT(L(1, 1), Catch::Matchers::WithinRel(2.0, 1e-10));
	REQUIRE_THAT(L(1, 2), Catch::Matchers::WithinRel(-1.0, 1e-10));

	REQUIRE_THAT(L(2, 0), Catch::Matchers::WithinAbs(0.0, 1e-10));
	REQUIRE_THAT(L(2, 1), Catch::Matchers::WithinRel(-1.0, 1e-10));
	REQUIRE_THAT(L(2, 2), Catch::Matchers::WithinRel(1.0, 1e-10));
}

TEST_CASE("Graph - Laplacian row sums are zero", "[Graph][matrix]")
{
	// For any graph, Laplacian row sums should be zero
	Graph<> g(5);
	g.addEdge(0, 1);
	g.addEdge(0, 2);
	g.addEdge(1, 2);
	g.addEdge(2, 3);
	g.addEdge(3, 4);

	Matrix<Real> L = g.toLaplacianMatrix();

	for (int i = 0; i < 5; ++i)
	{
		Real rowSum = 0;
		for (int j = 0; j < 5; ++j)
			rowSum += L(i, j);
		REQUIRE_THAT(rowSum, Catch::Matchers::WithinAbs(0.0, 1e-10));
	}
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         GRAPH FACTORY TESTS                                         ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Graph - createPathGraph creates correct structure", "[Graph][factory]")
{
	auto g = createPathGraph<int, Real>(5);

	REQUIRE(g.numVertices() == 5);
	REQUIRE(g.numEdges() == 4);

	// Check path structure
	for (size_t i = 0; i < 4; ++i)
		REQUIRE(g.hasEdge(i, i + 1));

	// End vertices have degree 1, middle have degree 2
	REQUIRE(g.degree(0) == 1);
	REQUIRE(g.degree(1) == 2);
	REQUIRE(g.degree(2) == 2);
	REQUIRE(g.degree(3) == 2);
	REQUIRE(g.degree(4) == 1);
}

TEST_CASE("Graph - createCycleGraph creates correct structure", "[Graph][factory]")
{
	auto g = createCycleGraph<int, Real>(5);

	REQUIRE(g.numVertices() == 5);
	REQUIRE(g.numEdges() == 5);

	// All vertices should have degree 2
	for (size_t i = 0; i < 5; ++i)
		REQUIRE(g.degree(i) == 2);

	// Check cycle structure
	REQUIRE(g.hasEdge(4, 0));  // Closing edge
}

TEST_CASE("Graph - createCompleteGraph creates correct structure", "[Graph][factory]")
{
	auto g = createCompleteGraph<int, Real>(4);

	REQUIRE(g.numVertices() == 4);
	REQUIRE(g.numEdges() == 6);  // n*(n-1)/2 = 4*3/2 = 6

	// All vertices should have degree n-1 = 3
	for (size_t i = 0; i < 4; ++i)
		REQUIRE(g.degree(i) == 3);

	// All pairs should be connected
	for (size_t i = 0; i < 4; ++i)
		for (size_t j = i + 1; j < 4; ++j)
			REQUIRE(g.hasEdge(i, j));
}

TEST_CASE("Graph - createStarGraph creates correct structure", "[Graph][factory]")
{
	auto g = createStarGraph<int, Real>(5);

	REQUIRE(g.numVertices() == 5);
	REQUIRE(g.numEdges() == 4);  // n-1 edges

	// Center vertex has degree n-1
	REQUIRE(g.degree(0) == 4);

	// Leaf vertices have degree 1
	for (size_t i = 1; i < 5; ++i)
		REQUIRE(g.degree(i) == 1);
}

TEST_CASE("Graph - createGridGraph creates correct structure", "[Graph][factory]")
{
	auto g = createGridGraph<int, Real>(3, 4);  // 3x4 grid

	REQUIRE(g.numVertices() == 12);
	// Edges: (rows-1)*cols + rows*(cols-1) = 2*4 + 3*3 = 8 + 9 = 17
	REQUIRE(g.numEdges() == 17);

	// Corner vertices have degree 2
	REQUIRE(g.degree(0) == 2);
	REQUIRE(g.degree(3) == 2);
	REQUIRE(g.degree(8) == 2);
	REQUIRE(g.degree(11) == 2);
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         PROPERTY TESTS                                              ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Graph - isWeighted detects weighted graphs", "[Graph][property]")
{
	Graph<> g1(3);
	g1.addEdge(0, 1);  // Default weight 1.0
	g1.addEdge(1, 2);

	REQUIRE_FALSE(g1.isWeighted());

	Graph<> g2(3);
	g2.addEdge(0, 1, 2.5);
	g2.addEdge(1, 2);

	REQUIRE(g2.isWeighted());
}

TEST_CASE("Graph - clear removes all vertices and edges", "[Graph][utility]")
{
	Graph<> g(5);
	g.addEdge(0, 1);
	g.addEdge(1, 2);

	REQUIRE(g.numVertices() == 5);
	REQUIRE(g.numEdges() == 2);

	g.clear();

	REQUIRE(g.numVertices() == 0);
	REQUIRE(g.numEdges() == 0);
	REQUIRE(g.isEmpty());
}

TEST_CASE("Graph - toString produces non-empty string", "[Graph][utility]")
{
	Graph<> g(3);
	g.addEdge(0, 1);
	g.addEdge(1, 2);

	std::string s = g.toString();

	REQUIRE_FALSE(s.empty());
	REQUIRE(s.find("Undirected") != std::string::npos);
	REQUIRE(s.find("V=3") != std::string::npos);
	REQUIRE(s.find("E=2") != std::string::npos);
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         GRAPH ALGORITHMS TESTS                                      ///
///                         BFS, DFS, Dijkstra, Connected Components                    ///
///////////////////////////////////////////////////////////////////////////////////////////

// GraphAlgorithms.h included at top of file

///////////////////////////////////////////////////////////////////////////////////////////
///                         BFS TESTS                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("BFS - on simple path graph", "[Graph][BFS]")
{
	// 0 -- 1 -- 2 -- 3 -- 4
	auto g = createPathGraph<int, Real>(5);

	auto result = BFS(g, 0);

	REQUIRE(result.nodesVisited == 5);
	REQUIRE(result.visitOrder.size() == 5);
	REQUIRE(result.visitOrder[0] == 0);  // Start vertex first

	// Check distances (hops)
	REQUIRE(result.distance[0] == 0);
	REQUIRE(result.distance[1] == 1);
	REQUIRE(result.distance[2] == 2);
	REQUIRE(result.distance[3] == 3);
	REQUIRE(result.distance[4] == 4);

	// Check parent pointers
	REQUIRE(result.parent[0] == static_cast<size_t>(-1));  // No parent for start
	REQUIRE(result.parent[1] == 0);
	REQUIRE(result.parent[2] == 1);
	REQUIRE(result.parent[3] == 2);
	REQUIRE(result.parent[4] == 3);
}

TEST_CASE("BFS - on cycle graph visits all vertices", "[Graph][BFS]")
{
	// 0 -- 1 -- 2 -- 3 -- 0 (cycle)
	auto g = createCycleGraph<int, Real>(4);

	auto result = BFS(g, 0);

	REQUIRE(result.nodesVisited == 4);
	REQUIRE(result.visitOrder.size() == 4);

	// Distance to opposite vertex (2) should be 2 (two paths of length 2)
	REQUIRE(result.distance[2] == 2);
}

TEST_CASE("BFS - on disconnected graph only visits reachable vertices", "[Graph][BFS]")
{
	Graph<> g(6);
	// Component 1: 0-1-2
	g.addEdge(0, 1);
	g.addEdge(1, 2);
	// Component 2: 3-4-5
	g.addEdge(3, 4);
	g.addEdge(4, 5);

	auto result = BFS(g, 0);

	REQUIRE(result.nodesVisited == 3);
	REQUIRE(result.distance[3] == std::numeric_limits<Real>::infinity());
	REQUIRE(result.distance[4] == std::numeric_limits<Real>::infinity());
	REQUIRE(result.distance[5] == std::numeric_limits<Real>::infinity());
}

TEST_CASE("BFS - on complete graph", "[Graph][BFS]")
{
	auto g = createCompleteGraph<int, Real>(5);

	auto result = BFS(g, 0);

	REQUIRE(result.nodesVisited == 5);

	// All vertices reachable in 1 hop from start
	for (size_t i = 1; i < 5; ++i)
		REQUIRE(result.distance[i] == 1);
}

TEST_CASE("BFS - with invalid start returns empty result", "[Graph][BFS]")
{
	Graph<> g(3);

	auto result = BFS(g, 10);  // Invalid vertex

	REQUIRE(result.nodesVisited == 0);
	REQUIRE(result.visitOrder.empty());
}

TEST_CASE("BFS - BFSVisit with callback", "[Graph][BFS]")
{
	auto g = createPathGraph<int, Real>(5);
	std::vector<size_t> visited;

	BFSVisit(g, 0, [&](size_t v) {
		visited.push_back(v);
		return true;  // Continue
	});

	REQUIRE(visited.size() == 5);
	REQUIRE(visited[0] == 0);
}

TEST_CASE("BFS - BFSVisit early termination", "[Graph][BFS]")
{
	auto g = createPathGraph<int, Real>(10);
	std::vector<size_t> visited;

	BFSVisit(g, 0, [&](size_t v) {
		visited.push_back(v);
		return v < 3;  // Stop after visiting vertex 3
	});

	REQUIRE(visited.size() == 4);  // Visited 0, 1, 2, 3 then stopped
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         DFS TESTS                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("DFS - on simple path graph", "[Graph][DFS]")
{
	auto g = createPathGraph<int, Real>(5);

	auto result = DFS(g, 0);

	REQUIRE(result.nodesVisited == 5);
	REQUIRE(result.visitOrder.size() == 5);
	REQUIRE(result.visitOrder[0] == 0);

	// All vertices should be visited
	std::vector<bool> seen(5, false);
	for (size_t v : result.visitOrder)
		seen[v] = true;
	for (bool s : seen)
		REQUIRE(s);
}

TEST_CASE("DFS - on tree structure", "[Graph][DFS]")
{
	//       0
	//      / \
	//     1   2
	//    / \
	//   3   4
	Graph<> g(5);
	g.addEdge(0, 1);
	g.addEdge(0, 2);
	g.addEdge(1, 3);
	g.addEdge(1, 4);

	auto result = DFS(g, 0);

	REQUIRE(result.nodesVisited == 5);
	REQUIRE(result.visitOrder[0] == 0);

	// Parent structure should form a tree
	REQUIRE(result.parent[0] == static_cast<size_t>(-1));
}

TEST_CASE("DFS - on disconnected graph", "[Graph][DFS]")
{
	Graph<> g(6);
	g.addEdge(0, 1);
	g.addEdge(1, 2);
	g.addEdge(3, 4);
	g.addEdge(4, 5);

	auto result = DFS(g, 0);

	REQUIRE(result.nodesVisited == 3);
}

TEST_CASE("DFS - DFSVisit with callback", "[Graph][DFS]")
{
	auto g = createCompleteGraph<int, Real>(4);
	std::vector<size_t> visited;

	DFSVisit(g, 0, [&](size_t v) {
		visited.push_back(v);
		return true;
	});

	REQUIRE(visited.size() == 4);
}

TEST_CASE("DFS - DFSVisit early termination", "[Graph][DFS]")
{
	auto g = createPathGraph<int, Real>(10);
	size_t count = 0;

	DFSVisit(g, 0, [&](size_t) {
		++count;
		return count < 5;  // Stop after 5 visits
	});

	REQUIRE(count == 5);
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         DIJKSTRA TESTS                                              ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Dijkstra - on unweighted path graph", "[Graph][Dijkstra]")
{
	auto g = createPathGraph<int, Real>(5);

	auto result = Dijkstra(g, 0);

	REQUIRE(result.distance[0] == 0);
	REQUIRE(result.distance[1] == 1);
	REQUIRE(result.distance[2] == 2);
	REQUIRE(result.distance[3] == 3);
	REQUIRE(result.distance[4] == 4);
}

TEST_CASE("Dijkstra - on weighted graph", "[Graph][Dijkstra]")
{
	//     1
	//   0---1
	//   |   |
	// 5 |   | 1
	//   |   |
	//   2---3
	//     1
	// In undirected: 0->2 can go 0->1->3->2 = 1+1+1 = 3 (shorter than direct 5)
	Graph<int, Real> g(4);
	g.addEdge(0, 1, 1.0);
	g.addEdge(0, 2, 5.0);
	g.addEdge(1, 3, 1.0);
	g.addEdge(2, 3, 1.0);

	auto result = Dijkstra(g, 0);

	REQUIRE(result.distance[0] == 0);
	REQUIRE(result.distance[1] == 1);
	REQUIRE(result.distance[2] == 3);  // 0->1->3->2 = 3 (shorter than direct 5)
	REQUIRE(result.distance[3] == 2);  // 0->1->3
}

TEST_CASE("Dijkstra - finds shorter path through intermediate vertices", "[Graph][Dijkstra]")
{
	// Classic example: direct path is longer than indirect
	//   0 ----10---- 2
	//    \          /
	//     1---1---1
	Graph<int, Real> g(3);
	g.addEdge(0, 2, 10.0);  // Direct: cost 10
	g.addEdge(0, 1, 1.0);   // Via 1: cost 1+1=2
	g.addEdge(1, 2, 1.0);

	auto result = Dijkstra(g, 0);

	REQUIRE(result.distance[2] == 2.0);  // Should find shorter path via 1
	REQUIRE(result.parent[2] == 1);
	REQUIRE(result.parent[1] == 0);
}

TEST_CASE("DijkstraPath - finds path and weight", "[Graph][Dijkstra]")
{
	Graph<int, Real> g(4);
	g.addEdge(0, 1, 2.0);
	g.addEdge(1, 2, 3.0);
	g.addEdge(2, 3, 1.0);
	g.addEdge(0, 3, 10.0);  // Direct but longer

	auto result = DijkstraPath(g, 0, 3);

	REQUIRE(result.found);
	REQUIRE(result.totalWeight == 6.0);  // 2+3+1
	REQUIRE(result.path.size() == 4);
	REQUIRE(result.path[0] == 0);
	REQUIRE(result.path[3] == 3);
}

TEST_CASE("DijkstraPath - start equals target", "[Graph][Dijkstra]")
{
	Graph<> g(3);
	g.addEdge(0, 1);
	g.addEdge(1, 2);

	auto result = DijkstraPath(g, 1, 1);

	REQUIRE(result.found);
	REQUIRE(result.totalWeight == 0);
	REQUIRE(result.path.size() == 1);
	REQUIRE(result.path[0] == 1);
}

TEST_CASE("DijkstraPath - no path exists", "[Graph][Dijkstra]")
{
	Graph<> g(4);
	g.addEdge(0, 1);
	// 2, 3 are disconnected

	auto result = DijkstraPath(g, 0, 3);

	REQUIRE_FALSE(result.found);
	REQUIRE_FALSE(result.diagnostics.empty());
}

TEST_CASE("DijkstraPath - invalid vertices", "[Graph][Dijkstra]")
{
	Graph<> g(3);

	auto result1 = DijkstraPath(g, 10, 0);
	REQUIRE_FALSE(result1.found);

	auto result2 = DijkstraPath(g, 0, 10);
	REQUIRE_FALSE(result2.found);
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         CONNECTED COMPONENTS TESTS                                  ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("ConnectedComponents - single component", "[Graph][Components]")
{
	auto g = createCompleteGraph<int, Real>(5);

	auto result = ConnectedComponents(g);

	REQUIRE(result.numComponents == 1);
	REQUIRE(result.components.size() == 1);
	REQUIRE(result.components[0].size() == 5);

	// All vertices should have same component ID
	for (size_t i = 0; i < 5; ++i)
		REQUIRE(result.componentId[i] == 0);
}

TEST_CASE("ConnectedComponents - two components", "[Graph][Components]")
{
	Graph<> g(6);
	// Component 0: 0-1-2
	g.addEdge(0, 1);
	g.addEdge(1, 2);
	// Component 1: 3-4-5
	g.addEdge(3, 4);
	g.addEdge(4, 5);

	auto result = ConnectedComponents(g);

	REQUIRE(result.numComponents == 2);
	REQUIRE(result.components.size() == 2);

	// Check component sizes
	REQUIRE(result.components[0].size() == 3);
	REQUIRE(result.components[1].size() == 3);

	// Vertices 0,1,2 should have same component ID (different from 3,4,5)
	REQUIRE(result.componentId[0] == result.componentId[1]);
	REQUIRE(result.componentId[1] == result.componentId[2]);
	REQUIRE(result.componentId[3] == result.componentId[4]);
	REQUIRE(result.componentId[4] == result.componentId[5]);
	REQUIRE(result.componentId[0] != result.componentId[3]);
}

TEST_CASE("ConnectedComponents - all isolated vertices", "[Graph][Components]")
{
	Graph<> g(4);
	// No edges - each vertex is its own component

	auto result = ConnectedComponents(g);

	REQUIRE(result.numComponents == 4);
	REQUIRE(result.components.size() == 4);

	for (size_t i = 0; i < 4; ++i)
	{
		REQUIRE(result.components[i].size() == 1);
		REQUIRE(result.componentId[i] == i);
	}
}

TEST_CASE("ConnectedComponents - empty graph", "[Graph][Components]")
{
	Graph<> g;

	auto result = ConnectedComponents(g);

	REQUIRE(result.numComponents == 0);
	REQUIRE(result.components.empty());
}

TEST_CASE("IsConnected - connected graph", "[Graph][Components]")
{
	auto g = createCycleGraph<int, Real>(5);

	REQUIRE(IsConnected(g));
}

TEST_CASE("IsConnected - disconnected graph", "[Graph][Components]")
{
	Graph<> g(4);
	g.addEdge(0, 1);
	g.addEdge(2, 3);
	// Two separate components

	REQUIRE_FALSE(IsConnected(g));
}

TEST_CASE("IsConnected - empty graph is connected", "[Graph][Components]")
{
	Graph<> g;

	REQUIRE(IsConnected(g));
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         PATH UTILITIES TESTS                                        ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("ReconstructPath - valid path", "[Graph][PathUtils]")
{
	// Parent array representing: 0 <- 1 <- 2 <- 3
	std::vector<size_t> parent = {static_cast<size_t>(-1), 0, 1, 2};

	auto path = ReconstructPath(parent, 0, 3);

	REQUIRE(path.size() == 4);
	REQUIRE(path[0] == 0);
	REQUIRE(path[1] == 1);
	REQUIRE(path[2] == 2);
	REQUIRE(path[3] == 3);
}

TEST_CASE("ReconstructPath - no path (disconnected)", "[Graph][PathUtils]")
{
	std::vector<size_t> parent = {static_cast<size_t>(-1), 0, static_cast<size_t>(-1)};

	auto path = ReconstructPath(parent, 0, 2);

	REQUIRE(path.empty());  // No path exists
}

TEST_CASE("ShortestPathUnweighted - finds shortest path", "[Graph][PathUtils]")
{
	//   0---1---2
	//   |       |
	//   3-------4
	Graph<> g(5);
	g.addEdge(0, 1);
	g.addEdge(1, 2);
	g.addEdge(0, 3);
	g.addEdge(3, 4);
	g.addEdge(2, 4);

	auto result = ShortestPathUnweighted(g, 0, 4);

	REQUIRE(result.found);
	REQUIRE(result.totalWeight == 2);  // 0->3->4 or 0->1->2->4 would be 3, but 0->3->4 is 2
	REQUIRE(result.path.size() == 3);
}

TEST_CASE("ShortestPathUnweighted - same start and target", "[Graph][PathUtils]")
{
	Graph<> g(3);
	g.addEdge(0, 1);

	auto result = ShortestPathUnweighted(g, 0, 0);

	REQUIRE(result.found);
	REQUIRE(result.totalWeight == 0);
	REQUIRE(result.path.size() == 1);
}

TEST_CASE("ShortestPathUnweighted - no path exists", "[Graph][PathUtils]")
{
	Graph<> g(4);
	g.addEdge(0, 1);
	// 2, 3 disconnected

	auto result = ShortestPathUnweighted(g, 0, 2);

	REQUIRE_FALSE(result.found);
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         BELLMAN-FORD TESTS                                          ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("BellmanFord - on simple path graph", "[Graph][BellmanFord]")
{
	auto g = createPathGraph<int, Real>(5);

	auto result = BellmanFord(g, 0);

	REQUIRE(result.nodesVisited == 5);
	REQUIRE(result.distance[0] == 0);
	REQUIRE(result.distance[1] == 1);
	REQUIRE(result.distance[2] == 2);
	REQUIRE(result.distance[3] == 3);
	REQUIRE(result.distance[4] == 4);
}

TEST_CASE("BellmanFord - with negative weights (no cycle)", "[Graph][BellmanFord]")
{
	// Directed graph with negative weights but no negative cycle
	// 0 --2--> 1 --(-1)--> 2
	Graph<int, Real> g(3, Graph<int, Real>::Type::Directed);
	g.addEdge(0, 1, 2.0);
	g.addEdge(1, 2, -1.0);

	auto result = BellmanFord(g, 0);

	REQUIRE(result.distance[0] == 0);
	REQUIRE(result.distance[1] == 2);
	REQUIRE(result.distance[2] == 1);  // 2 + (-1) = 1
}

TEST_CASE("BellmanFord - detects negative cycle", "[Graph][BellmanFord]")
{
	// Create a negative cycle: 0 -> 1 -> 2 -> 0 with total weight < 0
	Graph<int, Real> g(3, Graph<int, Real>::Type::Directed);
	g.addEdge(0, 1, 1.0);
	g.addEdge(1, 2, 1.0);
	g.addEdge(2, 0, -5.0);  // Creates negative cycle: 1 + 1 + (-5) = -3

	auto result = BellmanFord(g, 0);

	REQUIRE(result.parent.empty());  // Empty indicates negative cycle
	REQUIRE(result.distance.empty());
}

TEST_CASE("BellmanFord - finds shorter path through negative edge", "[Graph][BellmanFord]")
{
	// 0 --10--> 2
	// 0 --2--> 1 --(-5)--> 2
	Graph<int, Real> g(3, Graph<int, Real>::Type::Directed);
	g.addEdge(0, 2, 10.0);
	g.addEdge(0, 1, 2.0);
	g.addEdge(1, 2, -5.0);

	auto result = BellmanFord(g, 0);

	REQUIRE(result.distance[2] == -3.0);  // 2 + (-5) = -3 < 10
	REQUIRE(result.parent[2] == 1);
}

TEST_CASE("HasNegativeCycle - no cycle", "[Graph][BellmanFord]")
{
	Graph<int, Real> g(3, Graph<int, Real>::Type::Directed);
	g.addEdge(0, 1, 1.0);
	g.addEdge(1, 2, 1.0);

	REQUIRE_FALSE(HasNegativeCycle(g, 0));
}

TEST_CASE("HasNegativeCycle - with cycle", "[Graph][BellmanFord]")
{
	Graph<int, Real> g(3, Graph<int, Real>::Type::Directed);
	g.addEdge(0, 1, 1.0);
	g.addEdge(1, 2, 1.0);
	g.addEdge(2, 0, -5.0);

	REQUIRE(HasNegativeCycle(g, 0));
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         TOPOLOGICAL SORT TESTS                                      ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("TopologicalSort - simple DAG", "[Graph][TopologicalSort]")
{
	// 0 -> 1 -> 2 -> 3
	Graph<> g(4, Graph<>::Type::Directed);
	g.addEdge(0, 1);
	g.addEdge(1, 2);
	g.addEdge(2, 3);

	auto result = TopologicalSort(g);

	REQUIRE(result.isDAG);
	REQUIRE(result.order.size() == 4);

	// In topological order, earlier vertices come before later ones
	// Find positions
	std::vector<size_t> position(4);
	for (size_t i = 0; i < 4; ++i)
		position[result.order[i]] = i;

	REQUIRE(position[0] < position[1]);
	REQUIRE(position[1] < position[2]);
	REQUIRE(position[2] < position[3]);
}

TEST_CASE("TopologicalSort - diamond DAG", "[Graph][TopologicalSort]")
{
	//     1
	//    / \
	//   0   3
	//    \ /
	//     2
	Graph<> g(4, Graph<>::Type::Directed);
	g.addEdge(0, 1);
	g.addEdge(0, 2);
	g.addEdge(1, 3);
	g.addEdge(2, 3);

	auto result = TopologicalSort(g);

	REQUIRE(result.isDAG);
	REQUIRE(result.order.size() == 4);

	// Verify ordering constraints
	std::vector<size_t> position(4);
	for (size_t i = 0; i < 4; ++i)
		position[result.order[i]] = i;

	REQUIRE(position[0] < position[1]);
	REQUIRE(position[0] < position[2]);
	REQUIRE(position[1] < position[3]);
	REQUIRE(position[2] < position[3]);
}

TEST_CASE("TopologicalSort - detects cycle", "[Graph][TopologicalSort]")
{
	// 0 -> 1 -> 2 -> 0 (cycle)
	Graph<> g(3, Graph<>::Type::Directed);
	g.addEdge(0, 1);
	g.addEdge(1, 2);
	g.addEdge(2, 0);

	auto result = TopologicalSort(g);

	REQUIRE_FALSE(result.isDAG);
	REQUIRE(result.order.empty());
	REQUIRE_FALSE(result.diagnostics.empty());
}

TEST_CASE("TopologicalSort - fails on undirected graph", "[Graph][TopologicalSort]")
{
	Graph<> g(3);  // Undirected by default
	g.addEdge(0, 1);
	g.addEdge(1, 2);

	auto result = TopologicalSort(g);

	REQUIRE_FALSE(result.isDAG);
	REQUIRE(result.diagnostics.find("directed") != std::string::npos);
}

TEST_CASE("TopologicalSort - empty graph", "[Graph][TopologicalSort]")
{
	Graph<> g(0, Graph<>::Type::Directed);

	auto result = TopologicalSort(g);

	REQUIRE(result.isDAG);
	REQUIRE(result.order.empty());
}

TEST_CASE("TopologicalSort - single vertex", "[Graph][TopologicalSort]")
{
	Graph<> g(1, Graph<>::Type::Directed);

	auto result = TopologicalSort(g);

	REQUIRE(result.isDAG);
	REQUIRE(result.order.size() == 1);
	REQUIRE(result.order[0] == 0);
}

TEST_CASE("IsDAG - true for DAG", "[Graph][TopologicalSort]")
{
	Graph<> g(3, Graph<>::Type::Directed);
	g.addEdge(0, 1);
	g.addEdge(1, 2);

	REQUIRE(IsDAG(g));
}

TEST_CASE("IsDAG - false for cyclic graph", "[Graph][TopologicalSort]")
{
	Graph<> g(3, Graph<>::Type::Directed);
	g.addEdge(0, 1);
	g.addEdge(1, 2);
	g.addEdge(2, 0);

	REQUIRE_FALSE(IsDAG(g));
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         MINIMUM SPANNING TREE TESTS                                 ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Kruskal - simple triangle", "[Graph][MST]")
{
	//     2
	//   0---1
	//    \ /
	//   1 X 3
	//    / \
	//     2
	Graph<int, Real> g(3);
	g.addEdge(0, 1, 2.0);
	g.addEdge(1, 2, 3.0);
	g.addEdge(0, 2, 1.0);

	auto result = Kruskal(g);

	REQUIRE(result.isComplete);
	REQUIRE(result.edges.size() == 2);
	REQUIRE(result.totalWeight == 3.0);  // 1 + 2 (smallest edges)
}

TEST_CASE("Kruskal - path graph (already a tree)", "[Graph][MST]")
{
	auto g = createPathGraph<int, Real>(5);

	auto result = Kruskal(g);

	REQUIRE(result.isComplete);
	REQUIRE(result.edges.size() == 4);
	REQUIRE(result.totalWeight == 4.0);  // 4 edges of weight 1
}

TEST_CASE("Kruskal - complete graph", "[Graph][MST]")
{
	auto g = createCompleteGraph<int, Real>(4);

	auto result = Kruskal(g);

	REQUIRE(result.isComplete);
	REQUIRE(result.edges.size() == 3);  // V-1 edges
	REQUIRE(result.totalWeight == 3.0); // All edges have weight 1
}

TEST_CASE("Kruskal - weighted graph finds minimum", "[Graph][MST]")
{
	// Square graph with diagonals
	//  0--1--1
	//  |\/|
	//  |/\|
	//  2--3--2
	Graph<int, Real> g(4);
	g.addEdge(0, 1, 1.0);
	g.addEdge(1, 3, 2.0);
	g.addEdge(3, 2, 2.0);
	g.addEdge(2, 0, 1.0);
	g.addEdge(0, 3, 5.0);  // Diagonal (expensive)
	g.addEdge(1, 2, 5.0);  // Diagonal (expensive)

	auto result = Kruskal(g);

	REQUIRE(result.isComplete);
	REQUIRE(result.edges.size() == 3);
	REQUIRE(result.totalWeight == 4.0);  // Should pick 1+1+2 = 4
}

TEST_CASE("Kruskal - disconnected graph fails", "[Graph][MST]")
{
	Graph<int, Real> g(4);
	g.addEdge(0, 1, 1.0);
	// 2, 3 disconnected

	auto result = Kruskal(g);

	REQUIRE_FALSE(result.isComplete);
	REQUIRE_FALSE(result.diagnostics.empty());
}

TEST_CASE("Kruskal - empty graph", "[Graph][MST]")
{
	Graph<int, Real> g;

	auto result = Kruskal(g);

	REQUIRE(result.isComplete);
	REQUIRE(result.edges.empty());
	REQUIRE(result.totalWeight == 0);
}

TEST_CASE("Prim - simple triangle", "[Graph][MST]")
{
	Graph<int, Real> g(3);
	g.addEdge(0, 1, 2.0);
	g.addEdge(1, 2, 3.0);
	g.addEdge(0, 2, 1.0);

	auto result = Prim(g, 0);

	REQUIRE(result.isComplete);
	REQUIRE(result.edges.size() == 2);
	REQUIRE(result.totalWeight == 3.0);
}

TEST_CASE("Prim - complete graph", "[Graph][MST]")
{
	auto g = createCompleteGraph<int, Real>(4);

	auto result = Prim(g);

	REQUIRE(result.isComplete);
	REQUIRE(result.edges.size() == 3);
	REQUIRE(result.totalWeight == 3.0);
}

TEST_CASE("Prim - disconnected graph fails", "[Graph][MST]")
{
	Graph<int, Real> g(4);
	g.addEdge(0, 1, 1.0);

	auto result = Prim(g, 0);

	REQUIRE_FALSE(result.isComplete);
}

TEST_CASE("Prim and Kruskal produce same total weight", "[Graph][MST]")
{
	// Create a random-ish weighted graph
	Graph<int, Real> g(5);
	g.addEdge(0, 1, 4.0);
	g.addEdge(0, 2, 1.0);
	g.addEdge(1, 2, 2.0);
	g.addEdge(1, 3, 5.0);
	g.addEdge(2, 3, 3.0);
	g.addEdge(2, 4, 6.0);
	g.addEdge(3, 4, 2.0);

	auto kruskalResult = Kruskal(g);
	auto primResult = Prim(g);

	REQUIRE(kruskalResult.isComplete);
	REQUIRE(primResult.isComplete);
	REQUIRE(kruskalResult.totalWeight == primResult.totalWeight);
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         SPECTRAL GRAPH ANALYSIS TESTS                               ///
///////////////////////////////////////////////////////////////////////////////////////////

// GraphSpectral.h included at top of file

///////////////////////////////////////////////////////////////////////////////////////////
///                         CENTRALITY TESTS                                            ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("DegreeCentrality - complete graph has uniform centrality", "[Graph][Spectral]")
{
	auto g = createCompleteGraph<int, Real>(5);

	auto centrality = DegreeCentrality(g);

	REQUIRE(centrality.size() == 5);

	// All vertices have degree 4 out of max 4, so centrality = 1.0
	for (size_t i = 0; i < 5; ++i)
	{
		REQUIRE_THAT(centrality[i], Catch::Matchers::WithinAbs(1.0, 1e-10));
	}
}

TEST_CASE("DegreeCentrality - star graph", "[Graph][Spectral]")
{
	auto g = createStarGraph<int, Real>(5);  // Center + 4 leaves

	auto centrality = DegreeCentrality(g);

	REQUIRE(centrality.size() == 5);

	// Center (vertex 0) has degree 4, max is 4
	REQUIRE_THAT(centrality[0], Catch::Matchers::WithinAbs(1.0, 1e-10));

	// Leaves have degree 1, max is 4
	for (size_t i = 1; i < 5; ++i)
	{
		REQUIRE_THAT(centrality[i], Catch::Matchers::WithinAbs(0.25, 1e-10));
	}
}

TEST_CASE("ClosenessCentrality - path graph", "[Graph][Spectral]")
{
	auto g = createPathGraph<int, Real>(5);

	auto centrality = ClosenessCentrality(g);

	REQUIRE(centrality.size() == 5);

	// Middle vertex (2) should have highest closeness
	REQUIRE(centrality[2] > centrality[0]);
	REQUIRE(centrality[2] > centrality[4]);
	// Symmetry
	REQUIRE_THAT(centrality[0], Catch::Matchers::WithinAbs(centrality[4], 1e-10));
	REQUIRE_THAT(centrality[1], Catch::Matchers::WithinAbs(centrality[3], 1e-10));
}

TEST_CASE("ClosenessCentrality - star graph center is most central", "[Graph][Spectral]")
{
	auto g = createStarGraph<int, Real>(5);

	auto centrality = ClosenessCentrality(g);

	// Center is closest to all others (distance 1 to each)
	for (size_t i = 1; i < 5; ++i)
	{
		REQUIRE(centrality[0] > centrality[i]);
	}
}

TEST_CASE("BetweennessCentrality - path graph", "[Graph][Spectral]")
{
	auto g = createPathGraph<int, Real>(5);

	auto centrality = BetweennessCentrality(g, true);

	REQUIRE(centrality.size() == 5);

	// Middle vertex lies on all shortest paths
	REQUIRE(centrality[2] > centrality[1]);
	REQUIRE(centrality[2] > centrality[3]);

	// End vertices have 0 betweenness (not on any shortest path between others)
	REQUIRE_THAT(centrality[0], Catch::Matchers::WithinAbs(0.0, 1e-10));
	REQUIRE_THAT(centrality[4], Catch::Matchers::WithinAbs(0.0, 1e-10));
}

TEST_CASE("BetweennessCentrality - complete graph has zero betweenness", "[Graph][Spectral]")
{
	auto g = createCompleteGraph<int, Real>(5);

	auto centrality = BetweennessCentrality(g);

	// In complete graph, no vertex is on any shortest path between others
	// (all pairs are directly connected)
	for (size_t i = 0; i < 5; ++i)
	{
		REQUIRE_THAT(centrality[i], Catch::Matchers::WithinAbs(0.0, 1e-10));
	}
}

TEST_CASE("EigenvectorCentrality - regular graph has uniform centrality", "[Graph][Spectral]")
{
	auto g = createCycleGraph<int, Real>(5);

	auto centrality = EigenvectorCentrality(g);

	REQUIRE(centrality.size() == 5);

	// All vertices should have same eigenvector centrality in regular graph
	for (size_t i = 1; i < 5; ++i)
	{
		REQUIRE_THAT(centrality[i], Catch::Matchers::WithinAbs(centrality[0], 1e-5));
	}
}

TEST_CASE("PageRank - sums to 1", "[Graph][Spectral]")
{
	auto g = createCompleteGraph<int, Real>(5);

	auto ranks = PageRank(g);

	Real sum = 0;
	for (Real r : ranks)
		sum += r;

	REQUIRE_THAT(sum, Catch::Matchers::WithinAbs(1.0, 1e-10));
}

TEST_CASE("PageRank - regular graph has uniform ranks", "[Graph][Spectral]")
{
	auto g = createCycleGraph<int, Real>(5);

	auto ranks = PageRank(g);

	// All should be equal (0.2)
	for (size_t i = 0; i < 5; ++i)
	{
		REQUIRE_THAT(ranks[i], Catch::Matchers::WithinAbs(0.2, 1e-5));
	}
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         GRAPH STRUCTURE TESTS                                       ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("GraphDensity - complete graph has density 1", "[Graph][Spectral]")
{
	auto g = createCompleteGraph<int, Real>(5);

	Real density = GraphDensity(g);

	REQUIRE_THAT(density, Catch::Matchers::WithinAbs(1.0, 1e-10));
}

TEST_CASE("GraphDensity - path graph", "[Graph][Spectral]")
{
	auto g = createPathGraph<int, Real>(5);

	Real density = GraphDensity(g);

	// Path has 4 edges, max is 10 (5*4/2)
	REQUIRE_THAT(density, Catch::Matchers::WithinAbs(0.4, 1e-10));
}

TEST_CASE("GraphDensity - empty graph", "[Graph][Spectral]")
{
	Graph<> g(5);  // No edges

	Real density = GraphDensity(g);

	REQUIRE_THAT(density, Catch::Matchers::WithinAbs(0.0, 1e-10));
}

TEST_CASE("AverageClusteringCoefficient - complete graph is 1", "[Graph][Spectral]")
{
	auto g = createCompleteGraph<int, Real>(5);

	Real cc = AverageClusteringCoefficient(g);

	REQUIRE_THAT(cc, Catch::Matchers::WithinAbs(1.0, 1e-10));
}

TEST_CASE("AverageClusteringCoefficient - star graph is 0", "[Graph][Spectral]")
{
	auto g = createStarGraph<int, Real>(5);

	Real cc = AverageClusteringCoefficient(g);

	// Leaves have degree 1 (no triangles possible)
	// Center has neighbors that are not connected to each other
	REQUIRE_THAT(cc, Catch::Matchers::WithinAbs(0.0, 1e-10));
}

TEST_CASE("GraphDiameter - path graph", "[Graph][Spectral]")
{
	auto g = createPathGraph<int, Real>(5);

	Real diameter = GraphDiameter(g);

	REQUIRE_THAT(diameter, Catch::Matchers::WithinAbs(4.0, 1e-10));
}

TEST_CASE("GraphDiameter - complete graph is 1", "[Graph][Spectral]")
{
	auto g = createCompleteGraph<int, Real>(5);

	Real diameter = GraphDiameter(g);

	REQUIRE_THAT(diameter, Catch::Matchers::WithinAbs(1.0, 1e-10));
}

TEST_CASE("GraphDiameter - cycle graph", "[Graph][Spectral]")
{
	auto g = createCycleGraph<int, Real>(6);

	Real diameter = GraphDiameter(g);

	// Diameter is n/2 for even cycle
	REQUIRE_THAT(diameter, Catch::Matchers::WithinAbs(3.0, 1e-10));
}

TEST_CASE("AveragePathLength - complete graph is 1", "[Graph][Spectral]")
{
	auto g = createCompleteGraph<int, Real>(5);

	Real apl = AveragePathLength(g);

	REQUIRE_THAT(apl, Catch::Matchers::WithinAbs(1.0, 1e-10));
}

TEST_CASE("AveragePathLength - path graph", "[Graph][Spectral]")
{
	auto g = createPathGraph<int, Real>(5);

	Real apl = AveragePathLength(g);

	// Sum of all shortest paths / number of pairs
	// Pairs: (0,1)=1, (0,2)=2, (0,3)=3, (0,4)=4, (1,2)=1, (1,3)=2, (1,4)=3, (2,3)=1, (2,4)=2, (3,4)=1
	// Total = 20, pairs = 10
	REQUIRE_THAT(apl, Catch::Matchers::WithinAbs(2.0, 1e-10));
}

TEST_CASE("AlgebraicConnectivity - connected graph is positive", "[Graph][Spectral]")
{
	auto g = createCompleteGraph<int, Real>(5);

	Real lambda2 = AlgebraicConnectivityApprox(g);

	REQUIRE(lambda2 > 0);
}

TEST_CASE("AlgebraicConnectivity - disconnected graph is zero", "[Graph][Spectral]")
{
	Graph<> g(4);
	g.addEdge(0, 1);
	// 2, 3 isolated

	Real lambda2 = AlgebraicConnectivityApprox(g);

	REQUIRE_THAT(lambda2, Catch::Matchers::WithinAbs(0.0, 1e-10));
}

TEST_CASE("IsConnectedSpectral - matches IsConnected", "[Graph][Spectral]")
{
	auto g1 = createCompleteGraph<int, Real>(5);
	auto g2 = Graph<>(4);
	g2.addEdge(0, 1);

	REQUIRE(IsConnectedSpectral(g1) == IsConnected(g1));
	REQUIRE(IsConnectedSpectral(g2) == IsConnected(g2));
}

} // namespace MML::Tests::Algorithms::GraphTests

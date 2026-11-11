#include <catch2/catch_test_macros.hpp>

#ifdef MML_SINGLE_HEADER_INCLUDE
#include MML_SINGLE_HEADER_INCLUDE
#else
#include <mml/algorithms/GraphAlg.h>
#endif

#include <algorithm>
#include <limits>
#include <numeric>
#include <tuple>
#include <vector>

using namespace MML;

TEST_CASE("MST algorithms reject directed graphs before empty or start shortcuts", "[Graph][GraphContracts][MST]")
{
	for (const size_t count : {size_t{0}, size_t{1}, size_t{3}}) {
		Graph<> graph(count, Graph<>::Type::Directed);
		if (count == 3) {
			graph.addEdge(1, 0, 1);
			graph.addEdge(2, 0, 2);
		}
		for (const auto& result : {Kruskal(graph), Prim(graph), Prim(graph, GRAPH_NOT_FOUND)}) {
			CHECK(result.status == AlgorithmStatus::InvalidInput);
			CHECK(result.graphStatus == GraphAlgorithmStatus::GraphTypeMismatch);
			CHECK_FALSE(result.succeeded());
			CHECK_FALSE(result.isComplete);
			CHECK(result.edges.empty());
			CHECK(result.totalWeight == 0);
			CHECK_FALSE(result.message.empty());
			CHECK(result.diagnostics == result.message);
			CHECK(result.elapsed_time_ms >= 0);
		}
	}
}

TEST_CASE("Undirected MST controls retain weights and disconnected diagnostics", "[Graph][GraphContracts][MST]")
{
	Graph<> graph(3);
	graph.addEdge(1, 0, 1);
	graph.addEdge(2, 0, 2);
	for (const auto& result : {Kruskal(graph), Prim(graph), Prim(graph, 2)}) {
		CHECK(result.succeeded());
		CHECK(result.isComplete);
		CHECK(result.totalWeight == 3);
		CHECK(result.edges.size() == 2);
	}
	graph.addVertex();
	for (const auto& result : {Kruskal(graph), Prim(graph)}) {
		CHECK_FALSE(result.succeeded());
		CHECK_FALSE(result.isComplete);
		CHECK(result.graphStatus == GraphAlgorithmStatus::DisconnectedGraph);
	}
	for (const auto& result : {Kruskal(Graph<>{}), Prim(Graph<>{})}) {
		CHECK(result.succeeded());
		CHECK(result.isComplete);
		CHECK(result.edges.empty());
	}
}

TEST_CASE("Bipartite coloring and automatic matching accept an incoming directed edge", "[Graph][GraphContracts][Matching]")
{
	Graph<> graph(2, Graph<>::Type::Directed);
	graph.addEdge(1, 0);
	const auto coloring = BipartiteColoring(graph);
	CHECK(coloring.succeeded());
	CHECK(coloring.isBipartite);
	CHECK(coloring.color[0] != coloring.color[1]);
	for (const auto& matching : {HopcroftKarp(graph), HopcroftKarp(graph, std::vector<size_t>{1}),
		HopcroftKarp(graph, std::vector<size_t>{0})}) {
		CHECK(matching.succeeded());
		CHECK(matching.cardinality == 1);
		CHECK(matching.matching.size() == 1);
	}
}

TEST_CASE("Parent paths must arrive at the requested start", "[Graph][GraphContracts][PathUtils]")
{
	const std::vector<size_t> parent{GRAPH_NOT_FOUND, 0, 1, GRAPH_NOT_FOUND};
	CHECK(ReconstructPath(parent, 3, 2).empty());
	CHECK(ReconstructPath(parent, 0, 2) == std::vector<size_t>{0, 1, 2});
	CHECK(ReconstructPath(parent, 1, 2) == std::vector<size_t>{1, 2});
	CHECK(ReconstructPath(parent, 2, 2) == std::vector<size_t>{2});
	CHECK(ReconstructPath(parent, 0, 3).empty());
	CHECK_THROWS_AS(ReconstructPath(parent, parent.size(), 2), IndexError);
	CHECK_THROWS_AS(ReconstructPath(parent, 0, parent.size()), IndexError);
	CHECK_THROWS_AS(ReconstructPath({}, 0, 0), IndexError);
}

TEST_CASE("Parent path walks reject malformed links and terminate on cycles", "[Graph][GraphContracts][PathUtils]")
{
	CHECK_THROWS_AS(ReconstructPath({GRAPH_NOT_FOUND, 9}, 0, 1), IndexError);
	CHECK_THROWS_AS(ReconstructPath({GRAPH_NOT_FOUND, 1}, 0, 1), ArgumentError);
	CHECK_THROWS_AS(ReconstructPath({GRAPH_NOT_FOUND, 2, 1}, 0, 2), ArgumentError);
	CHECK(ReconstructPath({GRAPH_NOT_FOUND, 9}, 1, 1) == std::vector<size_t>{1});

	std::vector<size_t> parent(4096);
	std::iota(parent.begin() + 1, parent.end(), size_t{0});
	parent[0] = GRAPH_NOT_FOUND;
	const auto path = ReconstructPath(parent, 0, parent.size() - 1);
	REQUIRE(path.size() == parent.size());
	for (size_t i = 0; i < path.size(); ++i)
		CHECK(path[i] == i);
	parent[1] = parent.size() - 1;
	CHECK_THROWS_AS(ReconstructPath(parent, 0, parent.size() - 1), ArgumentError);
}

TEST_CASE("Parent path result accessors share bounded validation and unreachable behavior", "[Graph][GraphContracts][PathUtils]")
{
	const auto verify = [](auto result, Real unreachable) {
		CHECK(result.pathTo(0).empty());
		result.source = 0;
		result.parent = {GRAPH_NOT_FOUND, 0, 1, GRAPH_NOT_FOUND};
		result.distance = {0, 1, 2, unreachable};
		CHECK(result.pathTo(2) == std::vector<size_t>{0, 1, 2});
		CHECK(result.pathTo(0) == std::vector<size_t>{0});
		CHECK(result.pathTo(3).empty());
		CHECK_THROWS_AS(result.pathTo(4), IndexError);
		result.source = 3;
		CHECK(result.pathTo(2).empty());
		result.source = 8;
		CHECK_THROWS_AS(result.pathTo(2), IndexError);
		result.source = 0;
		result.parent[1] = 2;
		CHECK_THROWS_AS(result.pathTo(2), ArgumentError);
		result.parent[1] = 9;
		CHECK_THROWS_AS(result.pathTo(2), IndexError);
		result.distance.pop_back();
		CHECK_THROWS_AS(result.pathTo(2), ArgumentError);
	};
	verify(ShortestPathTreeResult{}, std::numeric_limits<Real>::infinity());
	verify(LongestPathDAGResult{}, -std::numeric_limits<Real>::infinity());
}

TEST_CASE("Directed bipartite matching is independent of every two-by-two edge orientation", "[Graph][GraphContracts][Matching]")
{
	for (unsigned pattern = 0; pattern < 256; ++pattern) {
		CAPTURE(pattern);
		Graph<> graph(6, Graph<>::Type::Directed);
		bool present[4]{};
		for (unsigned edge = 0; edge < 4; ++edge) {
			const unsigned orientation = (pattern >> (2 * edge)) & 3;
			const size_t u = edge / 2;
			const size_t v = 2 + edge % 2;
			present[edge] = orientation != 0;
			if (orientation & 1) graph.addEdge(u, v);
			if (orientation & 2) graph.addEdge(v, u);
		}
		const size_t expected = (present[0] && present[3]) || (present[1] && present[2]) ? 2 :
			(present[0] || present[1] || present[2] || present[3] ? 1 : 0);
		const auto coloring = BipartiteColoring(graph);
		REQUIRE(coloring.succeeded());
		CHECK(coloring.isBipartite);
		REQUIRE(coloring.color.size() == graph.numVertices());
		for (const auto& [u, v, weight] : graph.directedEdges()) {
			(void)weight;
			CHECK(coloring.color[u] != coloring.color[v]);
		}
		for (const auto& matching : {HopcroftKarp(graph), HopcroftKarp(graph, std::vector<size_t>{0, 1}),
			HopcroftKarp(graph, std::vector<size_t>{2, 3})}) {
			CHECK(matching.succeeded());
			CHECK(matching.cardinality == expected);
			CHECK(matching.matching.size() == expected);
			REQUIRE(matching.mate.size() == graph.numVertices());
			std::vector<bool> used(graph.numVertices(), false);
			for (const auto& [u, v] : matching.matching) {
				REQUIRE(u < used.size());
				REQUIRE(v < used.size());
				CHECK_FALSE(used[u]);
				CHECK_FALSE(used[v]);
				CHECK(matching.mate[u] == v);
				CHECK(matching.mate[v] == u);
				CHECK((graph.hasEdge(u, v) || graph.hasEdge(v, u)));
				used[u] = used[v] = true;
			}
			for (size_t v = 0; v < used.size(); ++v)
				if (!used[v]) CHECK(matching.mate[v] == GRAPH_NOT_FOUND);
		}
	}
}

TEST_CASE("Directed odd cycles and self loops cannot be bipartite in any orientation", "[Graph][GraphContracts][Matching]")
{
	for (unsigned pattern = 0; pattern < 27; ++pattern) {
		Graph<> graph(3, Graph<>::Type::Directed);
		unsigned orientations = pattern;
		for (unsigned edge = 0; edge < 3; ++edge) {
			const unsigned orientation = orientations % 3;
			orientations /= 3;
			const size_t u = edge;
			const size_t v = (edge + 1) % 3;
			if (orientation != 1) graph.addEdge(u, v);
			if (orientation != 0) graph.addEdge(v, u);
		}
		const auto coloring = BipartiteColoring(graph);
		CHECK_FALSE(coloring.succeeded());
		CHECK_FALSE(coloring.isBipartite);
		CHECK(coloring.graphStatus == GraphAlgorithmStatus::CycleDetected);
		const auto matching = HopcroftKarp(graph);
		CHECK_FALSE(matching.succeeded());
		CHECK(matching.graphStatus == GraphAlgorithmStatus::CycleDetected);
	}
	for (const auto type : {Graph<>::Type::Directed, Graph<>::Type::Undirected}) {
		Graph<> graph(1, type, true);
		graph.addEdge(0, 0);
		CHECK_FALSE(BipartiteColoring(graph).isBipartite);
		CHECK_FALSE(HopcroftKarp(graph).succeeded());
	}
}

TEST_CASE("Antiparallel flow entries are signed net flow with a valid max flow and cut", "[Graph][GraphContracts][Flow]")
{
	Graph<> graph(3, Graph<>::Type::Directed);
	graph.addEdge(0, 1, 5);
	graph.addEdge(1, 0, 3);
	graph.addEdge(1, 2, 5);
	graph.addEdge(2, 1, 2);
	graph.addEdge(0, 2, 0);
	for (const auto& result : {EdmondsKarp(graph, 0, 2), Dinic(graph, 0, 2)}) {
		CHECK(result.succeeded());
		CHECK(result.maxFlow == 5);
		REQUIRE(result.flowEdges.size() == 4);
		std::vector<std::vector<Real>> net(3, std::vector<Real>(3, 0));
		bool negative = false;
		for (const auto& [u, v, flow, capacity] : result.flowEdges) {
			CHECK(capacity == graph.edgeWeight(u, v));
			CHECK(flow <= capacity);
			CHECK(flow >= -graph.edgeWeight(v, u));
			net[u][v] = flow;
			negative = negative || flow < 0;
		}
		CHECK(negative);
		CHECK(net[0][1] == 5);
		CHECK(net[1][0] == -5);
		CHECK(net[1][2] == 5);
		CHECK(net[2][1] == -5);
		CHECK(net[0][1] == -net[1][0]);
		CHECK(net[1][2] == -net[2][1]);
		CHECK(net[1][0] + net[1][2] == 0);
		CHECK(net[0][1] + net[0][2] == result.maxFlow);
		CHECK(net[2][0] + net[2][1] == -result.maxFlow);
		Real cutCapacity = 0;
		for (const auto& [u, v, capacity] : result.cutEdges) {
			CHECK(std::find(result.sourceSide.begin(), result.sourceSide.end(), u) != result.sourceSide.end());
			CHECK(std::find(result.sinkSide.begin(), result.sinkSide.end(), v) != result.sinkSide.end());
			cutCapacity += capacity;
		}
		CHECK(cutCapacity == result.maxFlow);
	}
}

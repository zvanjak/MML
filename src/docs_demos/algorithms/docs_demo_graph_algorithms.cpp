///////////////////////////////////////////////////////////////////////////////////////////
///  File:        docs_demo_graph_algorithms.cpp                                      ///
///  Description: Brief demonstration of GraphAlg.h - graph traversal                 ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/MMLBase.h>
#include <mml/algorithms/GraphAlg.h>
#endif

#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <vector>

using namespace MML;

void Docs_Demo_GraphContracts()
{
    const std::vector<size_t> parent{GRAPH_NOT_FOUND, 0, 1, GRAPH_NOT_FOUND};
    const auto connected = ReconstructPath(parent, 0, 2);
    const auto differentRoot = ReconstructPath(parent, 3, 2);
    if (connected != std::vector<size_t>{0, 1, 2} || !differentRoot.empty())
        throw std::runtime_error("Graph parent-root contract failed");
    bool cycleRejected = false;
    try {
        (void)ReconstructPath({GRAPH_NOT_FOUND, 2, 1}, 0, 2);
    } catch (const ArgumentError& error) {
        cycleRejected = true;
        std::cout << "Malformed parent cycle: " << error.what() << '\n';
    }
    if (!cycleRejected)
        throw std::runtime_error("Graph parent-cycle contract failed");

    Graph<> incoming(2, Graph<>::Type::Directed);
    incoming.addEdge(1, 0);
    const auto coloring = BipartiteColoring(incoming);
    const auto matching = HopcroftKarp(incoming);
    if (!coloring.succeeded() || !coloring.isBipartite || !matching.succeeded() || matching.cardinality != 1)
        throw std::runtime_error("Directed bipartite contract failed");
    for (const auto& mst : {Kruskal(incoming), Prim(incoming)}) {
        if (mst.status != AlgorithmStatus::InvalidInput || mst.graphStatus != GraphAlgorithmStatus::GraphTypeMismatch ||
            mst.isComplete || !mst.edges.empty())
            throw std::runtime_error("Directed MST rejection contract failed");
        std::cout << "Directed MST: " << mst.message << '\n';
    }

    Graph<> network(3, Graph<>::Type::Directed);
    network.addEdge(0, 1, 5);
    network.addEdge(1, 0, 3);
    network.addEdge(1, 2, 5);
    network.addEdge(2, 1, 2);
    const auto flow = Dinic(network, 0, 2);
    const auto comparison = EdmondsKarp(network, 0, 2);
    if (!flow.succeeded() || !comparison.succeeded() || flow.maxFlow != 5 || comparison.maxFlow != 5)
        throw std::runtime_error("Antiparallel max-flow contract failed");
    bool negative = false;
    for (const auto& [from, to, signedNetFlow, capacity] : flow.flowEdges) {
        const Real oneWayFlow = std::max(signedNetFlow, Real{0});
        negative = negative || signedNetFlow < 0;
        if (signedNetFlow > capacity || signedNetFlow < -network.edgeWeight(to, from))
            throw std::runtime_error("Signed net-flow capacity contract failed");
        std::cout << from << " -> " << to << ": net=" << signedNetFlow
                  << ", one-way=" << oneWayFlow << ", capacity=" << capacity << '\n';
    }
    if (!negative)
        throw std::runtime_error("Antiparallel signed net-flow witness is missing");
}

void Docs_Demo_GraphAlgorithms()
{
    std::cout << "=== GraphAlgorithms Demo ===" << std::endl;
    
    // Create a weighted graph: 0--1--2
    //                          |     |
    //                          3-----4
    Graph<int, Real> graph;
    for (int i = 0; i < 5; i++) graph.addVertex(i);
    
    graph.addEdge(0, 1, 1.0);
    graph.addEdge(1, 2, 2.0);
    graph.addEdge(0, 3, 1.5);
    graph.addEdge(2, 4, 1.0);
    graph.addEdge(3, 4, 2.5);
    
    std::cout << "Graph: 5 vertices, " << graph.numEdges() << " edges" << std::endl;
    
    // BFS traversal
    std::cout << "\nBFS from vertex 0:" << std::endl;
    TraversalResult bfsResult = BFS(graph, 0);
    std::cout << "  Visit order: ";
    for (size_t v : bfsResult.visitOrder) std::cout << v << " ";
    std::cout << std::endl;
    std::cout << "  Nodes visited: " << bfsResult.nodesVisited << std::endl;
    
    // DFS traversal
    std::cout << "\nDFS from vertex 0:" << std::endl;
    TraversalResult dfsResult = DFS(graph, 0);
    std::cout << "  Visit order: ";
    for (size_t v : dfsResult.visitOrder) std::cout << v << " ";
    std::cout << std::endl;
    
    // Dijkstra shortest paths
    std::cout << "\nDijkstra from vertex 0:" << std::endl;
    TraversalResult dijkResult = Dijkstra(graph, 0);
    std::cout << std::fixed << std::setprecision(1);
    std::cout << "  Distances: ";
    for (size_t i = 0; i < 5; i++) 
        std::cout << i << ":" << dijkResult.distance[i] << " ";
    std::cout << std::endl;
    
    // Find shortest path 0 -> 4 using DijkstraPath
    PathResult pathResult = DijkstraPath(graph, 0, 4);
    std::cout << "  Path 0->4: ";
    for (size_t v : pathResult.path) std::cout << v << " ";
    std::cout << "(length: " << pathResult.totalWeight << ")" << std::endl;
    
    // Connected components
    ComponentsResult compResult = ConnectedComponents(graph);
    std::cout << "\nConnected components: " << compResult.numComponents << std::endl;
    
    std::cout << "=== Demo Complete ===" << std::endl;
}

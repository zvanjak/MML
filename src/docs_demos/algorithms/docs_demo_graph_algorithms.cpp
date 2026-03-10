///////////////////////////////////////////////////////////////////////////////////////////
///  File:        docs_demo_graph_algorithms.cpp                                      ///
///  Description: Brief demonstration of GraphAlgorithms.h - graph traversal          ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#include "algorithms/GraphAlgorithms.h"
#endif

#include <iostream>
#include <iomanip>

using namespace MML;

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

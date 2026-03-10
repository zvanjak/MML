///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        docs_demo_graph.cpp                                                 ///
///  Description: Documentation examples for Graph.h                                  ///
///               Graph data structure with adjacency list representation             ///
///                                                                                   ///
///  Usage:       Run MML_DocsDemo application to execute these examples              ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Graph.h"
#include "base/Matrix/Matrix.h"
#endif

#include <iostream>
#include <iomanip>

using namespace MML;

///////////////////////////////////////////////////////////////////////////////////////////
///                              GRAPH CREATION                                         ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_GraphCreation()
{
    std::cout << "\n=== Graph Creation ===\n\n";

    // Create an undirected graph
    Graph<> g1(Graph<>::Type::Undirected);
    
    // Add vertices (returns index)
    size_t v0 = g1.addVertex(0);
    size_t v1 = g1.addVertex(1);
    size_t v2 = g1.addVertex(2);
    size_t v3 = g1.addVertex(3);
    
    std::cout << "Created undirected graph with " << g1.numVertices() << " vertices.\n";
    std::cout << "Vertex indices: " << v0 << ", " << v1 << ", " << v2 << ", " << v3 << "\n\n";

    // Add edges with weights
    g1.addEdge(0, 1, 1.0);
    g1.addEdge(0, 2, 2.0);
    g1.addEdge(1, 2, 1.5);
    g1.addEdge(2, 3, 3.0);
    
    std::cout << "Added " << g1.numEdges() << " edges.\n";
    std::cout << "Graph properties:\n";
    std::cout << "  isDirected: " << (g1.isDirected() ? "yes" : "no") << "\n";
    std::cout << "  isWeighted: " << (g1.isWeighted() ? "yes" : "no") << "\n\n";

    // Create with pre-allocated vertices
    Graph<> g2(5, Graph<>::Type::Directed);
    std::cout << "Created directed graph with " << g2.numVertices() << " pre-allocated vertices.\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                              EDGE OPERATIONS                                        ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_EdgeOperations()
{
    std::cout << "\n=== Edge Operations ===\n\n";

    Graph<> g(4, Graph<>::Type::Undirected);
    
    // Add edges
    g.addEdge(0, 1, 1.5);
    g.addEdge(1, 2, 2.0);
    g.addEdge(2, 3, 2.5);
    g.addEdge(0, 3, 4.0);

    // Check edge existence
    std::cout << "Edge queries:\n";
    std::cout << "  hasEdge(0, 1): " << (g.hasEdge(0, 1) ? "yes" : "no") << "\n";
    std::cout << "  hasEdge(0, 2): " << (g.hasEdge(0, 2) ? "yes" : "no") << "\n";
    std::cout << "  hasEdge(1, 0): " << (g.hasEdge(1, 0) ? "yes" : "no") << " (undirected)\n\n";

    // Get edge weights
    std::cout << "Edge weights:\n";
    std::cout << "  edgeWeight(0, 1) = " << g.edgeWeight(0, 1) << "\n";
    std::cout << "  edgeWeight(2, 3) = " << g.edgeWeight(2, 3) << "\n\n";

    // Get neighbors
    std::cout << "Neighbors of vertex 1:\n";
    for (const auto& edge : g.neighbors(1)) {
        std::cout << "  -> vertex " << edge.to << " (weight " << edge.weight << ")\n";
    }

    // Vertex degrees
    std::cout << "\nVertex degrees:\n";
    for (size_t i = 0; i < g.numVertices(); ++i) {
        std::cout << "  degree(" << i << ") = " << g.degree(i) << "\n";
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
///                              DIRECTED GRAPHS                                        ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_DirectedGraph()
{
    std::cout << "\n=== Directed Graph ===\n\n";

    Graph<> g(4, Graph<>::Type::Directed);
    
    // Add directed edges
    g.addEdge(0, 1, 1.0);  // 0 -> 1
    g.addEdge(0, 2, 1.0);  // 0 -> 2
    g.addEdge(1, 2, 1.0);  // 1 -> 2
    g.addEdge(2, 3, 1.0);  // 2 -> 3
    g.addEdge(3, 0, 1.0);  // 3 -> 0 (creates cycle)

    std::cout << "Directed graph edges:\n";
    for (size_t i = 0; i < g.numVertices(); ++i) {
        std::cout << "  From vertex " << i << ": ";
        for (const auto& edge : g.neighbors(i)) {
            std::cout << i << "->" << edge.to << " ";
        }
        std::cout << "\n";
    }

    std::cout << "\nIn-degrees and Out-degrees:\n";
    for (size_t i = 0; i < g.numVertices(); ++i) {
        std::cout << "  Vertex " << i << ": out=" << g.outDegree(i) 
                  << ", in=" << g.inDegree(i) << "\n";
    }

    // Check directed edge asymmetry
    std::cout << "\nDirected edge asymmetry:\n";
    std::cout << "  hasEdge(0, 1): " << (g.hasEdge(0, 1) ? "yes" : "no") << "\n";
    std::cout << "  hasEdge(1, 0): " << (g.hasEdge(1, 0) ? "yes" : "no") << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                              MATRIX REPRESENTATIONS                                 ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_MatrixRepresentations()
{
    std::cout << "\n=== Matrix Representations ===\n\n";

    Graph<> g(4, Graph<>::Type::Undirected);
    g.addEdge(0, 1, 1.0);
    g.addEdge(0, 2, 1.0);
    g.addEdge(1, 2, 1.0);
    g.addEdge(2, 3, 1.0);

    std::cout << "Graph: 0--1, 0--2, 1--2, 2--3 (unit weights)\n\n";

    // Adjacency matrix
    Matrix<Real> A = g.toAdjacencyMatrix();
    std::cout << "Adjacency Matrix A:\n";
    A.Print(std::cout, 6, 1);

    // Degree matrix
    Matrix<Real> D = g.toDegreeMatrix();
    std::cout << "\nDegree Matrix D (diagonal):\n";
    D.Print(std::cout, 6, 1);

    // Laplacian matrix
    Matrix<Real> L = g.toLaplacianMatrix();
    std::cout << "\nLaplacian Matrix L = D - A:\n";
    L.Print(std::cout, 6, 1);

    std::cout << "\nLaplacian properties:\n";
    std::cout << "  - Row sums are zero\n";
    std::cout << "  - Positive semi-definite\n";
    std::cout << "  - Smallest eigenvalue is 0 (multiplicity = # components)\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                              GRAPH WITH DATA                                        ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_GraphWithData()
{
    std::cout << "\n=== Graph with Vertex Data ===\n\n";

    // Graph with string vertex labels
    Graph<std::string, Real> g(Graph<std::string, Real>::Type::Undirected);
    
    size_t a = g.addVertex("Alice");
    size_t b = g.addVertex("Bob");
    size_t c = g.addVertex("Carol");
    size_t d = g.addVertex("Dave");

    g.addEdge(a, b, 1.0);  // Alice -- Bob
    g.addEdge(a, c, 1.0);  // Alice -- Carol
    g.addEdge(b, c, 1.0);  // Bob -- Carol
    g.addEdge(c, d, 1.0);  // Carol -- Dave

    std::cout << "Social network graph:\n";
    for (size_t i = 0; i < g.numVertices(); ++i) {
        std::cout << "  " << g.vertexData(i) << " is connected to: ";
        for (const auto& edge : g.neighbors(i)) {
            std::cout << g.vertexData(edge.to) << " ";
        }
        std::cout << "\n";
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
///                              MAIN DEMO FUNCTION                                     ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_Graph()
{
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "               GRAPH DEMO\n";
    std::cout << "               Graph.h\n";
    std::cout << std::string(70, '=') << "\n";
    std::cout << "This demo covers:\n";
    std::cout << "  - Creating undirected and directed graphs\n";
    std::cout << "  - Edge operations (add, query, remove)\n";
    std::cout << "  - Vertex degrees and neighbors\n";
    std::cout << "  - Matrix representations (adjacency, Laplacian)\n";
    std::cout << "  - Graphs with custom vertex data\n";
    std::cout << std::string(70, '=') << "\n";

    Demo_GraphCreation();
    Demo_EdgeOperations();
    Demo_DirectedGraph();
    Demo_MatrixRepresentations();
    Demo_GraphWithData();

    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "               END OF GRAPH DEMO\n";
    std::cout << std::string(70, '=') << "\n";
}

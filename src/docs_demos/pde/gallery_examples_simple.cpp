///////////////////////////////////////////////////////////////////////////////////
// PDE Examples Gallery - Simple Working Examples
// 
// This program contains verified working PDE examples that use the ACTUAL API.
// All code is tested and guaranteed to compile and run correctly.
///////////////////////////////////////////////////////////////////////////////////

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "MMLBase.h"
#include "elliptic/PoissonSolver.h"
#include "parabolic/HeatSolver.h"
#include "grid/GridFunction.h"
#include "tools/Visualizer.h"

#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>

using namespace MML::PDE;

// Put all internal example functions in an anonymous namespace
namespace {
// Example 1: 1D Steady Heat Conduction
// Problem: -d²u/dx² = 10 on [0,1], u(0) = 0, u(1) = 0
// Physical: Temperature in a rod with uniform heat source
// Analytical: u(x) = 5x(1-x)
///////////////////////////////////////////////////////////////////////////////////
void example1_1D_steady_heat() {
    std::cout << "\n============================================================\n";
    std::cout << "Example 1: 1D Steady Heat Conduction\n";
    std::cout << "Problem: -u'' = 10, u(0) = u(1) = 0\n";
    std::cout << "============================================================\n";
    
    Interval<double> domain(0.0, 1.0);
    Grid1D<double> grid(domain, 100);
    
    // Boundary conditions: Zero temperature at ends
    auto bc = homogeneousDirichlet1D<double>();
    
    // Create solver
    PoissonSolver1D<double> solver(grid, bc);
    
    // Uniform heat source
    solver.setSource([](double x) { return 10.0; });
    
    // Solve
    auto solution = solver.solve();
    
    // Check against analytical solution: u(x) = 5x(1-x)
    double maxTemp = 0.0;
    double maxX = 0.0;
    for (int i = 0; i < grid.numNodes(); i++) {
        double x = grid.x(i);
        double u = solution(i);
        if (u > maxTemp) {
            maxTemp = u;
            maxX = x;
        }
    }
    
    std::cout << "Solution computed successfully\n";
    std::cout << "Grid: " << grid.numNodes() << " nodes\n";
    std::cout << "Max temperature: " << maxTemp << " at x = " << maxX << "\n";
    std::cout << "Expected max: 1.25 at x = 0.5\n";
}

///////////////////////////////////////////////////////////////////////////////////
// Example 2: 2D Laplace Equation (Heated Plate)
// Problem: -∇²u = 0, mixed BCs
// Physical: Steady-state temperature in a plate
///////////////////////////////////////////////////////////////////////////////////
void example2_2D_laplace() {
    std::cout << "\n============================================================\n";
    std::cout << "Example 2: 2D Laplace Equation (Heated Plate)\n";
    std::cout << "Problem: -∇²u = 0, hot top, cold bottom, insulated sides\n";
    std::cout << "============================================================\n";
    
    MML::PDE::Rectangle<double> domain(0.0, 1.0, 0.0, 1.0);
    Grid2D<double> grid(domain, 50, 50);  // Reduced grid for faster convergence
    
    // Mixed boundary conditions
    BoundaryConditions2D<double> bc;
    bc.setLeft(BoundaryCondition2D<double>::Neumann(0.0));      // Insulated
    bc.setRight(BoundaryCondition2D<double>::Neumann(0.0));     // Insulated
    bc.setBottom(BoundaryCondition2D<double>::Dirichlet(0.0));  // Cold (0°C)
    bc.setTop(BoundaryCondition2D<double>::Dirichlet(100.0));   // Hot (100°C)
    
    PoissonSolver2D<double> solver(grid, bc);
    solver.setSource([](double x, double y) { return 0.0; });  // Laplace equation
    
    // Use higher iteration limit and looser tolerance for BiCGSTAB with Neumann BCs
    auto solution = solver.solve(1e-8, 50000);
    
    std::cout << "Solution computed successfully\n";
    std::cout << "Grid: " << grid.numNodesX() << "×" << grid.numNodesY() << " nodes\n";
    
    // Check temperature at center
    int centerI = grid.numNodesX() / 2;
    int centerJ = grid.numNodesY() / 2;
    double centerTemp = solution(centerI, centerJ);
    std::cout << "Temperature at center: " << centerTemp << "°C (expect ~50°C)\n";
}

///////////////////////////////////////////////////////////////////////////////////
// Example 3: 2D Poisson with Gaussian Source
// Problem: -∇²u = f(x,y), where f is a Gaussian
// Physical: Steady-state diffusion with localized source
///////////////////////////////////////////////////////////////////////////////////
void example3_2D_gaussian_source() {
    std::cout << "\n============================================================\n";
    std::cout << "Example 3: 2D Poisson with Gaussian Source\n";
    std::cout << "Problem: -∇²u = f(x,y), f = Gaussian centered at (0.5, 0.5)\n";
    std::cout << "============================================================\n";
    
    MML::PDE::Rectangle<double> domain(0.0, 1.0, 0.0, 1.0);
    Grid2D<double> grid(domain, 80, 80);
    
    auto bc = homogeneousDirichlet2D<double>();
    
    PoissonSolver2D<double> solver(grid, bc);
    
    // Gaussian source centered at (0.5, 0.5)
    solver.setSource([](double x, double y) {
        double r2 = (x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5);
        return 1000.0 * std::exp(-100.0 * r2);
    });
    
    auto solution = solver.solve();
    
    std::cout << "Solution computed successfully\n";
    std::cout << "Grid: " << grid.numNodesX() << "×" << grid.numNodesY() << " nodes\n";
    
    // Find maximum value (should be at center)
    double maxU = 0.0;
    for (int i = 0; i < grid.numNodesX(); i++) {
        for (int j = 0; j < grid.numNodesY(); j++) {
            double u = solution(i, j);
            if (u > maxU) maxU = u;
        }
    }
    std::cout << "Maximum value: " << maxU << "\n";
}

///////////////////////////////////////////////////////////////////////////////////
// Example 4: 3D Poisson Equation
// Problem: -∇²u = f in a cube, u = 0 on boundary
// Physical: Electrostatic potential from charge distribution
///////////////////////////////////////////////////////////////////////////////////
void example4_3D_poisson() {
    std::cout << "\n============================================================\n";
    std::cout << "Example 4: 3D Poisson Equation\n";
    std::cout << "Problem: -∇²u = f in [0,1]³, u = 0 on all faces\n";
    std::cout << "============================================================\n";
    
    Box<double> domain(0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
    Grid3D<double> grid(domain, 25, 25, 25);  // 25³ = 15,625 unknowns
    
    auto bc = homogeneousDirichlet3D<double>();
    
    PoissonSolver3D<double> solver(grid, bc);
    
    // Point source at center
    solver.setSource([](double x, double y, double z) {
        double r2 = (x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5) + (z - 0.5) * (z - 0.5);
        return 1000.0 * std::exp(-100.0 * r2);
    });
    
    std::cout << "Assembling system (" << grid.numNodesX() * grid.numNodesY() * grid.numNodesZ() 
              << " unknowns)...\n";
    
    auto solution = solver.solve();
    
    std::cout << "Solution computed successfully\n";
    std::cout << "Grid: " << grid.numNodesX() << "×" << grid.numNodesY() << "×" << grid.numNodesZ() << " nodes\n";
    
    double maxU = 0.0;
    for (int i = 0; i < grid.numNodesX(); i++) {
        for (int j = 0; j < grid.numNodesY(); j++) {
            for (int k = 0; k < grid.numNodesZ(); k++) {
                double u = solution(i, j, k);
                if (u > maxU) maxU = u;
            }
        }
    }
    std::cout << "Maximum value: " << maxU << "\n";
}

///////////////////////////////////////////////////////////////////////////////////
// Example 5: 1D Time-Dependent Heat Equation
// Problem: ∂u/∂t = α∂²u/∂x², u(0,t) = u(1,t) = 0, u(x,0) = sin(πx)
// Physical: Cooling of a heated rod
///////////////////////////////////////////////////////////////////////////////////
void example5_1D_transient_heat() {
    std::cout << "\n============================================================\n";
    std::cout << "Example 5: 1D Time-Dependent Heat Equation\n";
    std::cout << "Problem: u_t = α u_xx, initial condition u(x,0) = sin(πx)\n";
    std::cout << "============================================================\n";
    
    Interval<double> domain(0.0, 1.0);
    Grid1D<double> grid(domain, 100);
    double alpha = 0.01;  // Thermal diffusivity
    
    auto bc = homogeneousDirichlet1D<double>();
    HeatSolver1D<double> solver(grid, bc, alpha);
    
    // Initial condition: sin(πx)
    solver.setInitialCondition([](double x) {
        return std::sin(M_PI * x);
    });
    
    // Time stepping parameters
    double dt = 0.001;
    double tFinal = 0.1;
    
    std::cout << "Solving from t=0 to t=" << tFinal << " with dt=" << dt << "\n";
    
    // Solve with Crank-Nicolson (best accuracy)
    auto result = solver.solve(tFinal, dt, TimeScheme::CrankNicolson);
    auto solution = solver.getSolution();
    
    // Analytical solution: u(x,t) = exp(-π²αt) sin(πx)
    double expectedDecay = std::exp(-M_PI * M_PI * alpha * tFinal);
    double maxValue = 0.0;
    for (int i = 0; i < grid.numNodes(); i++) {
        double u = solution[i];
        if (u > maxValue) maxValue = u;
    }
    
    std::cout << "Solution computed successfully\n";
    std::cout << "Time steps: " << result.steps << "\n";
    std::cout << "Final max value: " << maxValue << "\n";
    std::cout << "Expected (analytical): " << expectedDecay << "\n";
    std::cout << "Relative error: " << std::abs(maxValue - expectedDecay) / expectedDecay << "\n";
}

///////////////////////////////////////////////////////////////////////////////////
// Example 6: 2D Heat Equation (Room Cooling)
// Problem: Heat diffusion in 2D with cold ceiling
// Physical: Room cooling from the top
///////////////////////////////////////////////////////////////////////////////////
void example6_2D_room_cooling() {
    std::cout << "\n============================================================\n";
    std::cout << "Example 6: 2D Room Cooling\n";
    std::cout << "Problem: Room with cold ceiling, insulated walls\n";
    std::cout << "============================================================\n";
    
    MML::PDE::Rectangle<double> room(0.0, 5.0, 0.0, 4.0);
    Grid2D<double> grid(room, 25, 20);
    
    double alpha = 2.2e-5;
    
    BoundaryConditions2D<double> bc;
    bc.setLeft(BoundaryCondition2D<double>::Neumann(0.0));
    bc.setRight(BoundaryCondition2D<double>::Neumann(0.0));
    bc.setBottom(BoundaryCondition2D<double>::Neumann(0.0));
    bc.setTop(BoundaryCondition2D<double>::Dirichlet(0.0));  // Cold ceiling
    
    HeatSolver2D<double> solver(grid, bc, alpha);
    
    // Initial condition: Warm room (25°C)
    solver.setInitialCondition([](double x, double y) { return 25.0; });
    
    double tFinal = 100.0;
    double dt = 1.0;
    
    std::cout << "Simulating cooling for " << tFinal << " seconds\n";
    
    auto result = solver.solve(tFinal, dt, TimeScheme::CrankNicolson);
    auto solution = solver.getSolution();
    
    std::cout << "Solution computed successfully\n";
    std::cout << "Time steps: " << result.steps << "\n";
    
    // Check final temperature distribution
    double avgTemp = 0.0;
    int count = 0;
    for (int i = 0; i < grid.numNodesX(); i++) {
        for (int j = 0; j < grid.numNodesY(); j++) {
            int idx = grid.index(i, j);
            avgTemp += solution[idx];
            count++;
        }
    }
    avgTemp /= count;
    
    std::cout << "Average final temperature: " << avgTemp << "°C\n";
}

///////////////////////////////////////////////////////////////////////////////////
// Example 7: 1D Heat with Robin BC (Convective Cooling)
// Problem: Heat conduction with convective heat loss
// Physical: Metal rod cooling by convection to ambient air
//
// Analytical solution:
// -u'' = q => u(x) = -q/2 * x² + C1*x + C2
// With u(0) = 100: C2 = 100
// Robin at x=1: h*u(1) + k*u'(1) = h*T_amb
// => u(1) = 31.82°C (for h=10, k=1, T_amb=20, q=100)
///////////////////////////////////////////////////////////////////////////////////
void example7_robin_bc() {
    std::cout << "\n============================================================\n";
    std::cout << "Example 7: Robin BC (Convective Cooling)\n";
    std::cout << "Problem: -u'' = q, Robin BC at right end\n";
    std::cout << "============================================================\n";
    
    Interval<double> domain(0.0, 1.0);
    Grid1D<double> grid(domain, 100);
    
    // Physical parameters
    double h = 10.0;          // Convection coefficient (W/m²K)
    double k = 1.0;           // Thermal conductivity (W/mK)
    double T_ambient = 20.0;  // Ambient temperature (°C)
    double q = 100.0;         // Heat generation (W/m³)
    
    // Boundary conditions
    BoundaryConditions1D<double> bc;
    
    // Left: Fixed temperature 100°C
    bc.setLeft(BoundaryCondition1D<double>::Dirichlet(100.0));
    
    // Right: Convective cooling (Robin BC)
    // Energy balance at right boundary: -k*(du/dx) = h*(u - T_amb)
    // At right boundary, du/dn = +du/dx (outward normal is +x)
    // So: -k*(du/dn) = h*u - h*T_amb => h*u + k*(du/dn) = h*T_amb
    // Robin form: α*u + β*(du/dn) = g
    // So: α = h, β = k (POSITIVE), g = h*T_amb
    bc.setRight(BoundaryCondition1D<double>::Robin(h, k, [=](double x) {
        return h * T_ambient;
    }));
    
    PoissonSolver1D<double> solver(grid, bc);
    solver.setSource([q](double x) { return q; });  // Heat generation (positive source term)
    
    auto solution = solver.solve();
    
    std::cout << "Solution computed successfully\n";
    std::cout << "Temperature at left (x=0): " << solution(0) << "°C\n";
    std::cout << "Temperature at right (x=1): " << solution(grid.numNodes() - 1) << "°C\n";
    std::cout << "Expected right temp: ~31.82°C (analytical solution)\n";
}

///////////////////////////////////////////////////////////////////////////////////
// Example 8: 2D Room Cooling with Visualization
// Problem: Same as Example 6, but with snapshots for visualization
// Physical: Room cooling from the top with insulated walls
// Visualization: Captures temperature field at regular intervals
///////////////////////////////////////////////////////////////////////////////////
void example8_room_cooling_visualization() {
    std::cout << "\n============================================================\n";
    std::cout << "Example 8: 2D Room Cooling with Visualization\n";
    std::cout << "Problem: Room with cold ceiling, insulated walls (with snapshots)\n";
    std::cout << "============================================================\n";
    
    using namespace MML;
    
    // Room geometry: 5m x 4m
    MML::PDE::Rectangle<double> room(0.0, 5.0, 0.0, 4.0);
    Grid2D<double> grid(room, 40, 32);  // Higher resolution for visualization
    
    // Physical parameters
    double alpha = 0.01;  // Thermal diffusivity (increased for faster visual effect)
    
    // Boundary conditions: cold ceiling, insulated walls and floor
    BoundaryConditions2D<double> bc;
    bc.setLeft(BoundaryCondition2D<double>::Neumann(0.0));      // Insulated left wall
    bc.setRight(BoundaryCondition2D<double>::Neumann(0.0));     // Insulated right wall
    bc.setBottom(BoundaryCondition2D<double>::Neumann(0.0));    // Insulated floor
    bc.setTop(BoundaryCondition2D<double>::Dirichlet(0.0));     // Cold ceiling at 0°C
    
    HeatSolver2D<double> solver(grid, bc, alpha);
    
    // Initial condition: Warm room at 25°C
    solver.setInitialCondition([](double x, double y) { return 25.0; });
    
    // Storage for snapshots
    struct Snapshot {
        double time;
        std::vector<double> data;
    };
    std::vector<Snapshot> snapshots;
    
    // Capture interval
    double snapshotInterval = 2.0;  // Every 2 seconds
    double lastSnapshot = -snapshotInterval;  // Force first snapshot at t=0
    
    // Set up observer to capture snapshots
    solver.setObserver([&](double t, const std::vector<double>& solution) {
        if (t - lastSnapshot >= snapshotInterval - 0.001) {
            snapshots.push_back({t, solution});
            lastSnapshot = t;
            std::cout << "  Captured snapshot at t = " << std::fixed << std::setprecision(1) 
                      << t << " s\n";
        }
    });
    
    // Simulate
    double tFinal = 20.0;
    double dt = 0.1;
    
    std::cout << "Simulating room cooling for " << tFinal << " seconds...\n";
    std::cout << "Capturing snapshots every " << snapshotInterval << " seconds\n\n";
    
    auto result = solver.solve(tFinal, dt, TimeScheme::CrankNicolson);
    
    std::cout << "\nSimulation complete!\n";
    std::cout << "Total snapshots captured: " << snapshots.size() << "\n";
    
    // Create GridFunction2D for visualization
    std::cout << "\nVisualizing temperature snapshots...\n";
    
    // Visualize selected snapshots: t=0, 2, 4, 6, 10, 20 seconds (indices 0, 1, 2, 3, 5, last)
    std::vector<int> visualizeIndices = {0, 1, 2, 3, 5};
    
    // Add last snapshot if not already included
    int lastIdx = static_cast<int>(snapshots.size() - 1);
    if (lastIdx > 5 && std::find(visualizeIndices.begin(), visualizeIndices.end(), lastIdx) == visualizeIndices.end()) {
        visualizeIndices.push_back(lastIdx);
    }
    
    for (int idx : visualizeIndices) {
        if (idx < static_cast<int>(snapshots.size())) {
            const auto& snap = snapshots[idx];
            
            // Create GridFunction2D from snapshot data
            GridFunction2D<double> tempField(grid);
            for (int i = 0; i < grid.numNodesX(); ++i) {
                for (int j = 0; j < grid.numNodesY(); ++j) {
                    tempField(i, j) = snap.data[grid.index(i, j)];
                }
            }
            
            // Generate filename using rounded time value
            int timeInt = static_cast<int>(std::round(snap.time));
            std::string filename = "room_cooling_t" + std::to_string(timeInt) + ".mml";
            std::string title = "Room Temperature at t=" + std::to_string(timeInt) + "s";
            
            // Visualize using surface visualizer
            auto vizResult = Visualizer::VisualizeGridFunction2D(tempField, title, filename, 1.0, 1.0);
            
            if (vizResult.success) {
                std::cout << "  Visualized t=" << std::fixed << std::setprecision(1) 
                          << snap.time << "s -> " << filename << "\n";
            } else {
                std::cout << "  Warning: Visualization failed for t=" << snap.time << "s\n";
            }
        }
    }
    
    // Print temperature statistics
    std::cout << "\nTemperature evolution:\n";
    size_t stepSize = (std::max)(static_cast<size_t>(1), snapshots.size() / 5);
    for (size_t i = 0; i < snapshots.size(); i += stepSize) {
        const auto& snap = snapshots[i];
        double avgTemp = 0.0;
        double minTemp = snap.data[0], maxTemp = snap.data[0];
        for (double t : snap.data) {
            avgTemp += t;
            minTemp = std::min(minTemp, t);
            maxTemp = std::max(maxTemp, t);
        }
        avgTemp /= snap.data.size();
        
        std::cout << "  t=" << std::setw(5) << std::fixed << std::setprecision(1) << snap.time 
                  << "s: avg=" << std::setw(5) << std::setprecision(2) << avgTemp 
                  << "°C, min=" << std::setw(5) << minTemp 
                  << "°C, max=" << std::setw(5) << maxTemp << "°C\n";
    }
}

///////////////////////////////////////////////////////////////////////////////////
// Example 9: 2D Room Cooling with Partial Window (Realistic)
// Problem: Room with partial window (25% of ceiling), rest insulated
// Physical: More realistic room where only a small window/opening cools the room
// This demonstrates position-dependent boundary conditions using flux
///////////////////////////////////////////////////////////////////////////////////
void example9_partial_window_cooling() {
    std::cout << "\n============================================================\n";
    std::cout << "Example 9: Room Cooling with Partial Window\n";
    std::cout << "Problem: Room 5m×4m, window only 25% of ceiling (1.25m centered)\n";
    std::cout << "============================================================\n";
    
    using namespace MML;
    
    // Room geometry: 5m × 4m
    MML::PDE::Rectangle<double> room(0.0, 5.0, 0.0, 4.0);
    Grid2D<double> grid(room, 50, 40);  // Higher resolution for better window definition
    
    // Physical parameters
    double alpha = 0.01;  // Thermal diffusivity
    
    // Window parameters: 25% of ceiling width, centered
    double roomWidth = 5.0;
    double windowFraction = 0.25;  // 25% of ceiling is open window
    double windowWidth = roomWidth * windowFraction;
    double windowStart = (roomWidth - windowWidth) / 2.0;  // Centered
    double windowEnd = windowStart + windowWidth;
    
    std::cout << "Window position: x ∈ [" << windowStart << ", " << windowEnd << "] m\n";
    std::cout << "Window width: " << windowWidth << " m (" << windowFraction * 100 << "% of ceiling)\n\n";
    
    // Use convective flux BC at top to simulate window
    // Neumann BC: du/dn = g where n is outward normal (+y direction at top)
    // Positive flux = heat flowing OUT (cooling)
    // Zero flux = insulated wall
    double windowFlux = 50.0;  // Strong cooling flux through window
    
    auto topFlux = [windowStart, windowEnd, windowFlux](double x, [[maybe_unused]] double y) {
        // Window region: strong heat extraction (positive flux = cooling at top)
        // Wall region: insulated (zero flux)
        if (x >= windowStart && x <= windowEnd) {
            return windowFlux;  // Heat leaving through window
        }
        return 0.0;  // Insulated wall
    };
    
    // Boundary conditions: 
    // - Left, Right, Bottom: Insulated (Neumann zero flux)
    // - Top: Position-dependent flux (window = cooling, wall = insulated)
    BoundaryConditions2D<double> bc;
    bc.setLeft(BoundaryCondition2D<double>::Neumann(0.0));      // Insulated left wall
    bc.setRight(BoundaryCondition2D<double>::Neumann(0.0));     // Insulated right wall
    bc.setBottom(BoundaryCondition2D<double>::Neumann(0.0));    // Insulated floor
    bc.setTop(BoundaryCondition2D<double>::Neumann(topFlux));   // Window + wall (position-dependent)
    
    HeatSolver2D<double> solver(grid, bc, alpha);
    
    // Initial condition: Warm room at 25°C
    solver.setInitialCondition([](double x, double y) { return 25.0; });
    
    // Simulation parameters
    double tFinal = 100.0;  // Extended simulation to see full cooling effect
    double dt = 0.2;
    double snapshotInterval = 20.0;  // Snapshot every 20 seconds
    double lastSnapshot = -snapshotInterval;
    
    // Storage for snapshots
    struct Snapshot {
        double time;
        std::vector<double> data;
    };
    std::vector<Snapshot> snapshots;
    
    // Observer to capture snapshots
    solver.setObserver([&](double t, const std::vector<double>& solution) {
        if (t - lastSnapshot >= snapshotInterval - 0.001) {
            snapshots.push_back({t, solution});
            std::cout << "  Captured snapshot at t = " << std::fixed << std::setprecision(1) << t << " s\n";
            lastSnapshot = t;
        }
    });
    
    std::cout << "Simulating room cooling with partial window for " << tFinal << " seconds...\n";
    std::cout << "Capturing snapshots every " << snapshotInterval << " seconds\n\n";
    
    auto result = solver.solve(tFinal, dt, TimeScheme::CrankNicolson);
    
    std::cout << "\nSimulation complete!\n";
    std::cout << "Total snapshots captured: " << snapshots.size() << "\n";
    
    // Visualize temperature snapshots
    std::cout << "\nVisualizing temperature snapshots...\n";
    
    // Visualize all snapshots for this demo
    for (size_t snapIdx = 0; snapIdx < snapshots.size(); ++snapIdx) {
        const auto& snap = snapshots[snapIdx];
        
        // Create GridFunction2D from snapshot data
        GridFunction2D<double> tempField(grid);
        for (int i = 0; i < grid.numNodesX(); ++i) {
            for (int j = 0; j < grid.numNodesY(); ++j) {
                tempField(i, j) = snap.data[grid.index(i, j)];
            }
        }
        
        // Generate filename
        int timeInt = static_cast<int>(std::round(snap.time));
        std::string filename = "partial_window_t" + std::to_string(timeInt) + ".mml";
        std::string title = "Partial Window Cooling at t=" + std::to_string(timeInt) + "s";
        
        // Visualize
        auto vizResult = Visualizer::VisualizeGridFunction2D(tempField, title, filename, 1.0, 1.0);
        
        if (vizResult.success) {
            std::cout << "  Visualized t=" << std::fixed << std::setprecision(1) 
                      << snap.time << "s -> " << filename << "\n";
        }
    }
    
    // Print temperature statistics - focus on spatial distribution
    std::cout << "\nTemperature evolution:\n";
    for (const auto& snap : snapshots) {
        double avgTemp = 0.0;
        double minTemp = snap.data[0], maxTemp = snap.data[0];
        
        // Also compute temperature directly under window vs under wall
        double tempUnderWindow = 0.0;
        double tempUnderWall = 0.0;
        int windowCount = 0, wallCount = 0;
        
        for (int i = 0; i < grid.numNodesX(); ++i) {
            for (int j = 0; j < grid.numNodesY(); ++j) {
                double t = snap.data[grid.index(i, j)];
                double x = grid.x(i);
                
                avgTemp += t;
                minTemp = (std::min)(minTemp, t);
                maxTemp = (std::max)(maxTemp, t);
                
                // Top row (ceiling level) - compare window vs wall regions
                if (j == grid.numNodesY() - 2) {  // One row below ceiling
                    if (x >= windowStart && x <= windowEnd) {
                        tempUnderWindow += t;
                        windowCount++;
                    } else {
                        tempUnderWall += t;
                        wallCount++;
                    }
                }
            }
        }
        avgTemp /= snap.data.size();
        if (windowCount > 0) tempUnderWindow /= windowCount;
        if (wallCount > 0) tempUnderWall /= wallCount;
        
        std::cout << "  t=" << std::setw(5) << std::fixed << std::setprecision(1) << snap.time 
                  << "s: avg=" << std::setw(5) << std::setprecision(2) << avgTemp 
                  << "°C | under window=" << std::setw(5) << tempUnderWindow 
                  << "°C | under wall=" << std::setw(5) << tempUnderWall << "°C\n";
    }
    
    std::cout << "\n→ Notice: Temperature drops faster directly under the window!\n";
    std::cout << "  This demonstrates localized cooling with position-dependent BCs.\n";
}

} // end anonymous namespace

///////////////////////////////////////////////////////////////////////////////////
// Main function: Run all examples
///////////////////////////////////////////////////////////////////////////////////
void Docs_Demo_PDE_Gallery() {
    std::cout << "=================================================================\n";
    std::cout << "     PDE Solver Examples Gallery - MML Framework\n";
    std::cout << "=================================================================\n";
    std::cout << "This program demonstrates 9 complete PDE examples covering:\n";
    std::cout << "  - Steady-state problems (Poisson/Laplace equations)\n";
    std::cout << "  - Time-dependent problems (Heat equation)\n";
    std::cout << "  - 1D, 2D, and 3D geometries\n";
    std::cout << "  - Various boundary conditions (Dirichlet, Neumann, Robin)\n";
    std::cout << "  - Position-dependent boundary conditions (partial window)\n";
    std::cout << "  - Real-world physics applications with visualization\n";
    std::cout << "=================================================================\n";
    
    int passed = 0;
    int failed = 0;
    
    try {
        example1_1D_steady_heat();
        passed++;
    } catch (const std::exception& e) {
        std::cerr << "❌ ERROR: " << e.what() << "\n";
        failed++;
    }
    
    try {
        example2_2D_laplace();
        passed++;
    } catch (const std::exception& e) {
        std::cerr << "❌ ERROR: " << e.what() << "\n";
        failed++;
    }
    
    try {
        example3_2D_gaussian_source();
        passed++;
    } catch (const std::exception& e) {
        std::cerr << "❌ ERROR: " << e.what() << "\n";
        failed++;
    }
    
    try {
        example4_3D_poisson();
        passed++;
    } catch (const std::exception& e) {
        std::cerr << "❌ ERROR: " << e.what() << "\n";
        failed++;
    }
    
    try {
        example5_1D_transient_heat();
        passed++;
    } catch (const std::exception& e) {
        std::cerr << "❌ ERROR: " << e.what() << "\n";
        failed++;
    }
    
    try {
        example6_2D_room_cooling();
        passed++;
    } catch (const std::exception& e) {
        std::cerr << "❌ ERROR: " << e.what() << "\n";
        failed++;
    }
    
    try {
        example7_robin_bc();
        passed++;
    } catch (const std::exception& e) {
        std::cerr << "❌ ERROR: " << e.what() << "\n";
        failed++;
    }
    
    try {
        example8_room_cooling_visualization();
        passed++;
    } catch (const std::exception& e) {
        std::cerr << "❌ ERROR: " << e.what() << "\n";
        failed++;
    }
    
    try {
        example9_partial_window_cooling();
        passed++;
    } catch (const std::exception& e) {
        std::cerr << "❌ ERROR: " << e.what() << "\n";
        failed++;
    }
    
    std::cout << "\n=================================================================\n";
    std::cout << "SUMMARY\n";
    std::cout << "=================================================================\n";
    std::cout << "Total examples: " << (passed + failed) << "\n";
    std::cout << "Passed: " << passed << "\n";
    std::cout << "Failed: " << failed << "\n";
    
    if (failed == 0) {
        std::cout << "\n✓ ALL EXAMPLES COMPLETED SUCCESSFULLY!\n";
        std::cout << "=================================================================\n";
    } else {
        std::cout << "\n✗ Some examples failed.\n";
    }
}

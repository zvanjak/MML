/**
 * Example 4: Room Temperature Simulation
 * 
 * Simulates heat diffusion in a 2D room with:
 * - Initial temperature: 20°C everywhere
 * - Cold window on left wall (0°C)
 * - Warm heater on right wall (25°C)
 * - Fixed temperature walls (15°C on top/bottom)
 * 
 * NOTE: Currently using Dirichlet BCs on all walls because
 * Neumann BC support in HeatSolver2D is incomplete.
 * 
 * Physical interpretation:
 * - Room dimensions: 5m × 4m
 * - Thermal diffusivity: 0.02 m²/s
 * - Simulate 30 seconds to see temperature gradient
 * 
 * From: docs/pde/Quick_Start_Guide.md
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

// Define M_PI for MSVC
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "MMLBase.h"
#include "parabolic/HeatSolver.h"
#include "tools/Visualizer.h"

using namespace MML;
using namespace MML::PDE;

// Helper functions and classes for Example 4
namespace {

// Initial condition: 20°C everywhere
double ex4_initial_temperature(double x, double y) {
    return 20.0;
}

// Observer to monitor and save temperature evolution
class Ex4_TemperatureMonitor {
private:
    std::ofstream outfile;
    double room_width;
    double room_height;
    int save_interval;
    int step_count;
    const Grid2D<double>* grid_ptr;
    
public:
    Ex4_TemperatureMonitor(const std::string& filename, double width, double height, int interval = 100)
        : room_width(width), room_height(height), save_interval(interval), step_count(0), grid_ptr(nullptr) {
        outfile.open(filename);
        if (!outfile) {
            throw std::runtime_error("Cannot open file: " + filename);
        }
        // CSV header
        outfile << "time,x,y,temperature\n";
    }
    
    void setGrid(const Grid2D<double>* grid) {
        grid_ptr = grid;
    }
    
    void operator()(double t, const std::vector<double>& solution) {
        if (!grid_ptr) return;
        
        if (step_count % save_interval == 0) {
            int nx = grid_ptr->numNodesX();
            int ny = grid_ptr->numNodesY();
            
            // Save full grid
            for (int i = 0; i < nx; ++i) {
                for (int j = 0; j < ny; ++j) {
                    double x = grid_ptr->x(i);
                    double y = grid_ptr->y(j);
                    int idx = grid_ptr->index(i, j);
                    double temp = solution[idx];
                    
                    outfile << std::fixed << std::setprecision(4)
                            << t << "," << x << "," << y << "," << temp << "\n";
                }
            }
            
            // Console progress
            int centerIdx = grid_ptr->index(nx/2, ny/2);
            double center_temp = solution[centerIdx];
            std::cout << "t = " << std::setw(6) << std::fixed << std::setprecision(2) 
                      << t << " s: Center temperature = " 
                      << std::setw(6) << std::setprecision(2) 
                      << center_temp << " °C\n";
        }
        step_count++;
    }
    
    ~Ex4_TemperatureMonitor() {
        if (outfile.is_open()) {
            outfile.close();
        }
    }
};

} // end anonymous namespace

void Docs_Demo_PDE_RoomTemperature() {
    // Physical parameters
    const double room_width = 5.0;   // meters
    const double room_height = 4.0;  // meters
    const double alpha = 0.02;       // Thermal diffusivity - reduced for slower cooling
    const double T_final = 30.0;     // Simulation time (seconds) - shows interesting gradient
    const double dt = 0.1;           // Time step (seconds)
    
    // Spatial discretization
    MML::PDE::Rectangle<double> domain(0.0, room_width, 0.0, room_height);
    Grid2D<double> grid(domain, 50, 40);  // 50x40 cells = 51x41 nodes
    
    std::cout << "Room Temperature Simulation\n";
    std::cout << "===========================\n\n";
    
    std::cout << "Physical parameters:\n";
    std::cout << "  Room size: " << room_width << " m × " << room_height << " m\n";
    std::cout << "  Thermal diffusivity: " << alpha << " m²/s (reduced for slower cooling)\n";
    std::cout << "  Initial temperature: 20.0 °C\n";
    std::cout << "  Window (left): 0.0 °C\n";
    std::cout << "  Heater (right): 25.0 °C\n";
    std::cout << "  Top/Bottom walls: 15.0 °C\n";
    std::cout << "  Simulation time: " << T_final << " s (to capture interesting gradient)\n";
    std::cout << "  Time step: " << dt << " s\n\n";
    
    std::cout << "Grid parameters:\n";
    std::cout << "  Grid: " << grid.numNodesX() << " × " << grid.numNodesY() << " nodes\n";
    std::cout << "  dx = " << grid.dx() << " m, dy = " << grid.dy() << " m\n\n";
    
    // Boundary conditions (all Dirichlet - Neumann not fully supported in 2D yet):
    BoundaryConditions2D<double> bc;
    bc.setLeft(BoundaryCondition2D<double>::Dirichlet(0.0));    // Cold window at 0°C
    bc.setRight(BoundaryCondition2D<double>::Dirichlet(25.0));  // Warm heater at 25°C
    bc.setBottom(BoundaryCondition2D<double>::Dirichlet(15.0)); // Wall at 15°C
    bc.setTop(BoundaryCondition2D<double>::Dirichlet(15.0));    // Wall at 15°C
    
    // Create solver
    HeatSolver2D<double> heat_solver(grid, bc, alpha);
    heat_solver.setInitialCondition(ex4_initial_temperature);
    
    // Create observer to monitor and save results
    Ex4_TemperatureMonitor monitor("room_cooling.csv", room_width, room_height, 100);
    monitor.setGrid(&grid);
    
    // Set observer before solving
    heat_solver.setObserver([&](double t, const std::vector<double>& sol) {
        monitor(t, sol);
    });
    
    std::cout << "Starting simulation...\n\n";
    
    // Solve (Crank-Nicolson: unconditionally stable, 2nd order accurate)
    auto result = heat_solver.solve(T_final, dt, TimeScheme::CrankNicolson);
    
    // Get final solution from solver
    const auto& solution = heat_solver.getSolution();
    
    std::cout << "\n======================\n";
    std::cout << "Simulation Complete!\n";
    std::cout << "======================\n\n";
    
    std::cout << "Solver completed " << result.steps << " time steps\n";
    std::cout << "Solution stable: " << (result.stable ? "Yes" : "No") << "\n\n";
    
    // Final temperature distribution at key points
    std::cout << "Final temperature distribution:\n";
    std::cout << "Location\t\t\tTemperature (°C)\n";
    std::cout << "------------------------------------------------\n";
    
    int nx = grid.numNodesX();
    int ny = grid.numNodesY();
    
    auto temp_at_point = [&](int i, int j, const std::string& name) {
        int idx = grid.index(i, j);
        double temp = solution[idx];
        std::cout << std::left << std::setw(30) << name 
                  << std::fixed << std::setprecision(2) << temp << "\n";
    };
    
    temp_at_point(0, ny/2, "Window (left, cold)");
    temp_at_point(nx/4, ny/2, "Near window");
    temp_at_point(nx/2, ny/2, "Room center");
    temp_at_point(3*nx/4, ny/2, "Near heater");
    temp_at_point(nx-1, ny/2, "Heater (right, warm)");
    temp_at_point(nx/2, 0, "Bottom wall");
    temp_at_point(nx/2, ny-1, "Top wall");
    
    std::cout << "\nResults saved to: room_cooling.csv\n";
    
    // ====================================================================
    // VISUALIZATION: Room cooling temperature distribution
    // Domain is already in meters (5m x 4m), scale by 20 for nice display
    // Temperature in °C, scale by 4 for better height visibility
    // ====================================================================
    std::cout << "\n========================================\n";
    std::cout << "Launching temperature visualization...\n";
    std::cout << "(Domain scaled 20x, temperature scaled 4x)\n";
    std::cout << "========================================\n";
    
    // Create GridFunction2D from solution vector for visualization
    GridFunction2D<double> solutionGrid(grid, solution);
    
    auto vizResult = Visualizer::VisualizeGridFunction2D(
        solutionGrid, 
        "Room Temperature - Cold Window (Left) vs Warm Heater (Right)",
        "pde_room_cooling_final.mml",
        20.0,   // Scale coordinates by 20 (5m → 100 units)
        4.0     // Scale temperature by 4 for better height visibility
    );
    
    if (!vizResult.success) {
        std::cerr << "Visualization failed: " << vizResult.errorMessage << "\n";
        return;
    }
    
    std::cout << "Visualization complete!\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
/// @file main.cpp
/// @brief Example 05: Rigid Body Collision Simulator
/// @details Simulates two parallelepipeds in a cubic container with elastic collisions.
///
/// SCENARIO:
/// - Container: 10m × 10m × 10m cubic box (walls at ±5m)
/// - Box 1: 2m × 1m × 0.6m, mass 10kg, starting at (-2, 0, 0)
/// - Box 2: 1.5m × 1m × 0.8m, mass 8kg, starting at (+2, 0, 0)
/// - Both boxes have initial translational and rotational velocities
/// - Perfectly elastic collisions (coefficient of restitution = 1.0)
///
/// OUTPUT:
/// - Console: Energy/momentum conservation tracking
/// - CSV file: Trajectory data for visualization
///
/// @author Generated for MinimalMathLibrary
/// @date January 2026
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#include "mml/tools/Visualizer.h"
#endif

// Self-contained rigid body simulation (no MPL dependency)
#include "RigidBodyCore.h"
#include "RigidBodySimulator.h"
#include "RigidBodySerializer.h"

#include <filesystem>

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace MML;
using namespace RigidBodySim;

/// @brief Print simulation header
void PrintHeader()
{
    std::cout << "\n";
    std::cout << "╔═══════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║     RIGID BODY COLLISION SIMULATOR - Two Boxes + Sphere           ║\n";
    std::cout << "╠═══════════════════════════════════════════════════════════════════╣\n";
    std::cout << "║  Container: 10m × 10m × 10m cubic box                             ║\n";
    std::cout << "║  Box 1: 2m × 1m × 0.6m, 10kg (red)                                ║\n";
    std::cout << "║  Box 2: 1.5m × 1m × 0.8m, 8kg (blue)                              ║\n";
    std::cout << "║  Sphere: radius 0.5m, 5kg (green)                                 ║\n";
    std::cout << "║  Perfectly elastic collisions (e = 1.0)                           ║\n";
    std::cout << "╚═══════════════════════════════════════════════════════════════════╝\n\n";
}

int main()
{
    PrintHeader();
    
    try {
    // ===== Create Box 1 =====
    Real mass1 = 10.0;  // kg
    Real halfA1 = 1.0, halfB1 = 0.5, halfC1 = 0.3;  // 2m × 1m × 0.6m
    RigidBodyBox box1(mass1, halfA1, halfB1, halfC1);
    
    box1.Position() = Vec3Cart(-2.0, 0.5, 0.0);
    box1.Velocity() = Vec3Cart(3.0, 0.5, 0.2);
    box1.Orientation() = Quaternion::FromAxisAngle(Vec3Cart(0, 0, 1), 0.1);  // Slight tilt
    box1.AngularVel() = Vec3Cart(0.3, 0.5, 0.2);
    
    // ===== Create Box 2 =====
    Real mass2 = 8.0;  // kg
    Real halfA2 = 0.75, halfB2 = 0.5, halfC2 = 0.4;  // 1.5m × 1m × 0.8m
    RigidBodyBox box2(mass2, halfA2, halfB2, halfC2);
    
    box2.Position() = Vec3Cart(2.0, -0.5, 0.3);
    box2.Velocity() = Vec3Cart(-2.0, 0.3, -0.1);
    box2.Orientation() = Quaternion::FromAxisAngle(Vec3Cart(1, 0, 0), -0.15);
    box2.AngularVel() = Vec3Cart(-0.2, 0.4, -0.3);
    
    // ===== Create Sphere =====
    Real mass3 = 5.0;  // kg
    Real radius = 0.35; // Smaller, faster sphere (0.7m diameter)
    RigidBodySphere sphere(mass3, radius);
    
    sphere.Position() = Vec3Cart(-3.0, 3.5, -2.0);  // Far corner, away from boxes
    sphere.Velocity() = Vec3Cart(3.0, -6.0, 1.5); // 3x speed - HAVOC MODE!
    sphere.AngularVel() = Vec3Cart(3.0, 1.5, 1.0); // Fast spin!
    
    // ===== Print Initial Conditions =====
    std::cout << "Initial Conditions:\n";
    std::cout << "───────────────────────────────────────────────────────────────────\n";
    std::cout << "Box 1: pos = (" << box1.Position().X() << ", " << box1.Position().Y() << ", " << box1.Position().Z() << ")\n";
    std::cout << "       vel = (" << box1.Velocity().X() << ", " << box1.Velocity().Y() << ", " << box1.Velocity().Z() << ")\n";
    std::cout << "       ω   = (" << box1.AngularVel().X() << ", " << box1.AngularVel().Y() << ", " << box1.AngularVel().Z() << ") rad/s\n";
    std::cout << "       KE  = " << box1.KineticEnergy() << " J (trans: " << box1.TranslationalKineticEnergy() 
              << ", rot: " << box1.RotationalKineticEnergy() << ")\n\n";
    
    std::cout << "Box 2: pos = (" << box2.Position().X() << ", " << box2.Position().Y() << ", " << box2.Position().Z() << ")\n";
    std::cout << "       vel = (" << box2.Velocity().X() << ", " << box2.Velocity().Y() << ", " << box2.Velocity().Z() << ")\n";
    std::cout << "       ω   = (" << box2.AngularVel().X() << ", " << box2.AngularVel().Y() << ", " << box2.AngularVel().Z() << ") rad/s\n";
    std::cout << "       KE  = " << box2.KineticEnergy() << " J (trans: " << box2.TranslationalKineticEnergy() 
              << ", rot: " << box2.RotationalKineticEnergy() << ")\n\n";
    
    std::cout << "Sphere: pos = (" << sphere.Position().X() << ", " << sphere.Position().Y() << ", " << sphere.Position().Z() << ")\n";
    std::cout << "        vel = (" << sphere.Velocity().X() << ", " << sphere.Velocity().Y() << ", " << sphere.Velocity().Z() << ")\n";
    std::cout << "        ω   = (" << sphere.AngularVel().X() << ", " << sphere.AngularVel().Y() << ", " << sphere.AngularVel().Z() << ") rad/s\n";
    std::cout << "        KE  = " << sphere.KineticEnergy() << " J (trans: " << sphere.TranslationalKineticEnergy() 
              << ", rot: " << sphere.RotationalKineticEnergy() << ")\n\n";
    
    // ===== Configure Simulator =====
    SimulationConfig config;
    config.timeStep = 0.001;          // 1ms timestep
    config.totalTime = 40.0;          // 40 second simulation
    config.containerHalfSize = 5.0;   // 10m × 10m × 10m container
    config.coeffRestitution = 1.0;    // Perfectly elastic
    
    RigidBodySimulator simulator(config);
    simulator.AddBody(box1);
    simulator.AddBody(box2);
    simulator.AddBody(sphere);
    
    // Store initial conservation values
    Real initialEnergy = simulator.TotalKineticEnergy();
    Real initialTransKE = box1.TranslationalKineticEnergy() + box2.TranslationalKineticEnergy() + sphere.TranslationalKineticEnergy();
    Real initialRotKE = box1.RotationalKineticEnergy() + box2.RotationalKineticEnergy() + sphere.RotationalKineticEnergy();
    Vec3Cart initialMomentum = simulator.TotalLinearMomentum();
    Vec3Cart initialAngMomentum = simulator.TotalAngularMomentum();
    
    std::cout << "System Totals:\n";
    std::cout << "  Initial kinetic energy:  " << std::fixed << std::setprecision(6) << initialEnergy << " J\n";
    std::cout << "  Initial lin. momentum:   (" << initialMomentum.X() << ", " 
              << initialMomentum.Y() << ", " << initialMomentum.Z() << ") kg·m/s\n";
    std::cout << "  Initial ang. momentum:   (" << initialAngMomentum.X() << ", " 
              << initialAngMomentum.Y() << ", " << initialAngMomentum.Z() << ") kg·m²/s\n\n";
    
    // ===== Run Simulation =====
    std::cout << "Running simulation for " << config.totalTime << " seconds...\n";
    std::cout << "───────────────────────────────────────────────────────────────────\n";
    std::cout << std::flush;  // Flush output before heavy computation
    
    Real lastReportTime = 0.0;
    
    simulator.SetProgressCallback([&](Real t, const RigidBodySimulator::BodyContainer& bodies)
    {
        // Print progress every second
        if (t - lastReportTime >= 1.0)
        {
            Real energy = 0;
            for (const auto& b : bodies)
                energy += b->KineticEnergy();
            
            Real energyError = std::abs(energy - initialEnergy) / initialEnergy * 100.0;
            
            std::cout << "t = " << std::fixed << std::setprecision(1) << std::setw(5) << t 
                      << "s  |  E = " << std::setprecision(4) << energy 
                      << " J  |  ΔE = " << std::setprecision(4) << energyError << "%\n" << std::flush;
            
            lastReportTime = t;
        }
    });
    
    simulator.Run();
    
    std::cout << "\n───────────────────────────────────────────────────────────────────\n";
    std::cout << "Simulation complete!\n\n";
    
    // ===== Final Results =====
    const auto& fb1 = simulator.GetBody(0);
    const auto& fb2 = simulator.GetBody(1);
    const auto& fs = simulator.GetBody(2);
    Real finalEnergy = simulator.TotalKineticEnergy();
    Real finalTransKE = fb1.TranslationalKineticEnergy() + fb2.TranslationalKineticEnergy() + fs.TranslationalKineticEnergy();
    Real finalRotKE = fb1.RotationalKineticEnergy() + fb2.RotationalKineticEnergy() + fs.RotationalKineticEnergy();
    Vec3Cart finalMomentum = simulator.TotalLinearMomentum();
    Vec3Cart finalAngMomentum = simulator.TotalAngularMomentum();
    
    Real energyError = std::abs(finalEnergy - initialEnergy) / initialEnergy * 100.0;
    Real initialMomMag = initialMomentum.NormL2();
    Real finalMomMag = finalMomentum.NormL2();
    Real initialAngMomMag = initialAngMomentum.NormL2();
    Real finalAngMomMag = finalAngMomentum.NormL2();
    
    std::cout << "╔═══════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║                        CONSERVATION CHECK                         ║\n";
    std::cout << "╠═══════════════════════════════════════════════════════════════════╣\n";
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "║  KINETIC ENERGY:                                                  ║\n";
    std::cout << "║  Initial Total:  " << std::setw(12) << initialEnergy << " J  (trans: " 
              << std::setw(9) << initialTransKE << ", rot: " << std::setw(9) << initialRotKE << ") ║\n";
    std::cout << "║  Final Total:    " << std::setw(12) << finalEnergy << " J  (trans: "
              << std::setw(9) << finalTransKE << ", rot: " << std::setw(9) << finalRotKE << ") ║\n";
    std::cout << "║  Energy Error:   " << std::setw(12) << energyError << " %";
    if (energyError < 1.0) std::cout << "  ✓ CONSERVED"; else std::cout << "  ✗ ERROR";
    std::cout << "                ║\n";
    std::cout << "╠═══════════════════════════════════════════════════════════════════╣\n";
    std::cout << "║  LINEAR MOMENTUM (vector sum changes on wall bounces):            ║\n";
    std::cout << "║  Initial: (" << std::setw(9) << initialMomentum.X() << ", " 
              << std::setw(9) << initialMomentum.Y() << ", " << std::setw(9) << initialMomentum.Z() << ")   ║\n";
    std::cout << "║  Final:   (" << std::setw(9) << finalMomentum.X() << ", " 
              << std::setw(9) << finalMomentum.Y() << ", " << std::setw(9) << finalMomentum.Z() << ")   ║\n";
    std::cout << "║  |p| Initial: " << std::setw(12) << initialMomMag 
              << "   |p| Final: " << std::setw(12) << finalMomMag << "      ║\n";
    std::cout << "╠═══════════════════════════════════════════════════════════════════╣\n";
    std::cout << "║  ANGULAR MOMENTUM (vector sum changes on wall bounces):           ║\n";
    std::cout << "║  Initial: (" << std::setw(9) << initialAngMomentum.X() << ", " 
              << std::setw(9) << initialAngMomentum.Y() << ", " << std::setw(9) << initialAngMomentum.Z() << ")   ║\n";
    std::cout << "║  Final:   (" << std::setw(9) << finalAngMomentum.X() << ", " 
              << std::setw(9) << finalAngMomentum.Y() << ", " << std::setw(9) << finalAngMomentum.Z() << ")   ║\n";
    std::cout << "║  |L| Initial: " << std::setw(12) << initialAngMomMag 
              << "   |L| Final: " << std::setw(12) << finalAngMomMag << "      ║\n";
    std::cout << "╚═══════════════════════════════════════════════════════════════════╝\n\n";
    
    // ===== Write Trajectory Data =====
    RigidBodySerializer serializer;
    serializer.SetContainerSize(config.containerHalfSize);
    serializer.AddBody(RigidBodyInfo::Box(mass1, halfA1, halfB1, halfC1, "Red", "Box1"));
    serializer.AddBody(RigidBodyInfo::Box(mass2, halfA2, halfB2, halfC2, "Blue", "Box2"));
    serializer.AddBody(RigidBodyInfo::Sphere(mass3, radius, "LimeGreen", "Sphere"));
    
    // Save combined simulation (every 10th frame to reduce file size)
    auto result = serializer.SaveSimulation("rigid_body_simulation.mml", 
                                            simulator.GetHistory(), 
                                            config.timeStep, 10);
    if (result)
        std::cout << result.message << "\n";
    else
        std::cerr << "Error: " << result.message << "\n";
    
    // Also save individual body trajectories
    serializer.SaveSingleBody("rigid_body_box1.mml", simulator.GetHistory(), 0, config.timeStep, 10);
    serializer.SaveSingleBody("rigid_body_box2.mml", simulator.GetHistory(), 1, config.timeStep, 10);
    serializer.SaveSingleBody("rigid_body_sphere.mml", simulator.GetHistory(), 2, config.timeStep, 10);
    
    // ===== Launch Visualizer with All Bodies =====
    std::cout << "\nLaunching visualizer for all bodies...\n";
    auto cwd = std::filesystem::current_path();
    std::vector<std::string> trajectoryFiles = {
        (cwd / "rigid_body_box1.mml").string(),
        (cwd / "rigid_body_box2.mml").string(),
        (cwd / "rigid_body_sphere.mml").string()
    };
    auto vizResult = Visualizer::VisualizeRigidBodyTrajectory(trajectoryFiles);
    if (!vizResult) {
        std::cerr << "Visualizer failed: " << vizResult.errorMessage << "\n";
        std::cerr << "You can manually run:\n";
        std::cerr << "  tools\\visualizers\\win\\WPF\\MML_RigidBodyMovement_Visualizer\\MML_RigidBodyMovement_Visualizer.exe ";
        for (const auto& f : trajectoryFiles) std::cerr << "\"" << f << "\" ";
        std::cerr << "\n";
    }
    
    // ===== Final Positions =====
    std::cout << "\nFinal State:\n";
    std::cout << "───────────────────────────────────────────────────────────────────\n";
    
    std::cout << "Box 1: pos = (" << fb1.Position().X() << ", " 
              << fb1.Position().Y() << ", " << fb1.Position().Z() << ")\n";
    std::cout << "       vel = (" << fb1.Velocity().X() << ", " 
              << fb1.Velocity().Y() << ", " << fb1.Velocity().Z() << ")\n";
    
    std::cout << "Box 2: pos = (" << fb2.Position().X() << ", " 
              << fb2.Position().Y() << ", " << fb2.Position().Z() << ")\n";
    std::cout << "       vel = (" << fb2.Velocity().X() << ", " 
              << fb2.Velocity().Y() << ", " << fb2.Velocity().Z() << ")\n";
    
    std::cout << "Sphere: pos = (" << fs.Position().X() << ", " 
              << fs.Position().Y() << ", " << fs.Position().Z() << ")\n";
    std::cout << "        vel = (" << fs.Velocity().X() << ", " 
              << fs.Velocity().Y() << ", " << fs.Velocity().Z() << ")\n";
    
    } catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}

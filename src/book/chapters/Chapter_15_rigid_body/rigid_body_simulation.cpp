#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "mml/algorithms/ODESolvers.h"
#include "mml/tools/Visualizer.h"

#include "mpl/RigidBody/RigidBodies.h"
#include "mpl/RigidBody/RigidBodySimulator.h"
#include "mpl/RigidBody/RigidBodySerializer.h"
#endif

#include <filesystem>
#include <iomanip>

using namespace MML;
using namespace MPL;

///////////////////////////////////////////////////////////////////////////////////////////
/// @brief Demo 1: Single tumbling box in free space (no gravity, no container)
/// @details Demonstrates Euler equations for rigid body rotation.
///          Validates conservation of:
///          - Total kinetic energy
///          - Angular momentum magnitude
///          - Quaternion normalization
///
/// Uses RigidBodySimulator in single-body mode with no container collisions.
///////////////////////////////////////////////////////////////////////////////////////////
void Chapter15_Demo_TumblingBox()
{
    std::cout << "\n";
    std::cout << "╔═══════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║       CHAPTER 15 DEMO 1: Tumbling Box (Euler Equations)           ║\n";
    std::cout << "╠═══════════════════════════════════════════════════════════════════╣\n";
    std::cout << "║  Box: 2m × 1m × 0.6m, mass = 10 kg                                ║\n";
    std::cout << "║  Initial ω = (2.0, 1.0, 0.5) rad/s                                ║\n";
    std::cout << "║  Simulation: 10 seconds, dt = 0.001s                              ║\n";
    std::cout << "║  Free space (no gravity, no container)                            ║\n";
    std::cout << "╚═══════════════════════════════════════════════════════════════════╝\n\n";

    // Create asymmetric box (different dimensions = different principal moments)
    Real mass = 10.0;  // kg
    Real halfA = 1.0;  // half-width X (full width = 2m)
    Real halfB = 0.5;  // half-width Y (full width = 1m)
    Real halfC = 0.3;  // half-width Z (full width = 0.6m)
    
    RigidBodyBox box(mass, halfA, halfB, halfC);
    
    // Set initial state - asymmetric spin for interesting tumbling motion
    box.Position() = Vec3Cart(0, 0, 0);
    box.Velocity() = Vec3Cart(0.5, 0.2, 0);      // Slow drift
    box.Orientation() = Quaternion::Identity();
    box.AngularVel() = Vec3Cart(2.0, 1.0, 0.5);  // Asymmetric spin
    
    // Print initial inertia tensor
    std::cout << "Inertia tensor (body frame, diagonal):\n";
    auto I = box.InertiaTensorBody();
    std::cout << "  I_xx = " << std::fixed << std::setprecision(4) << I(0,0) << " kg·m²\n";
    std::cout << "  I_yy = " << I(1,1) << " kg·m²\n";
    std::cout << "  I_zz = " << I(2,2) << " kg·m²\n\n";
    
    // Configure simulator - no container (huge bounds so no wall collisions)
    SimulationConfig config;
    config.timeStep = 0.001;          // 1ms timestep for accuracy
    config.totalTime = 10.0;          // 10 second simulation
    config.containerHalfSize = 1000.0; // Effectively infinite (no wall collisions)
    config.coeffRestitution = 1.0;
    
    RigidBodySimulator simulator(config);
    simulator.AddBody(box);
    
    // Store initial conservation values
    Real initialEnergy = simulator.TotalKineticEnergy();
    Vec3Cart initialAngMom = simulator.TotalAngularMomentum();
    Real initialAngMomMag = initialAngMom.NormL2();
    
    std::cout << "Initial conditions:\n";
    std::cout << "  Kinetic energy:     " << std::setprecision(6) << initialEnergy << " J\n";
    std::cout << "  Angular momentum:   |L| = " << initialAngMomMag << " kg·m²/s\n";
    std::cout << "  Quaternion norm:    " << box.Orientation().Norm() << "\n\n";
    
    // Track max errors during simulation
    Real maxEnergyError = 0.0;
    Real maxAngMomError = 0.0;
    Real lastReportTime = 0.0;
    
    std::cout << "Time (s) | Kinetic E (J) | ΔE/E₀ (%)  | |L| (kg·m²/s) | ΔL/L₀ (%)\n";
    std::cout << "---------|---------------|------------|---------------|----------\n";
    
    // Progress callback to track conservation
    simulator.SetProgressCallback([&](Real t, const RigidBodySimulator::BodyContainer& bodies)
    {
        Real energy = bodies[0]->KineticEnergy();
        Vec3Cart angMom = bodies[0]->AngularMomentum();
        Real angMomMag = angMom.NormL2();
        
        Real energyError = std::abs(energy - initialEnergy) / initialEnergy * 100.0;
        Real angMomError = std::abs(angMomMag - initialAngMomMag) / initialAngMomMag * 100.0;
        
        maxEnergyError = std::max(maxEnergyError, energyError);
        maxAngMomError = std::max(maxAngMomError, angMomError);
        
        // Print every second
        if (t - lastReportTime >= 1.0)
        {
            std::cout << std::fixed << std::setprecision(2) << std::setw(8) << t << " | "
                      << std::setprecision(6) << std::setw(13) << energy << " | "
                      << std::setprecision(6) << std::setw(10) << energyError << " | "
                      << std::setprecision(6) << std::setw(13) << angMomMag << " | "
                      << std::setprecision(6) << std::setw(9) << angMomError << "\n";
            lastReportTime = t;
        }
    });
    
    // Run simulation
    std::cout << std::flush;
    simulator.Run();
    
    // Final validation
    const auto& finalBody = simulator.GetBody(0);
    Real finalEnergy = finalBody.KineticEnergy();
    Vec3Cart finalAngMom = finalBody.AngularMomentum();
    Real finalAngMomMag = finalAngMom.NormL2();
    Real qNorm = finalBody.Orientation().Norm();
    
    std::cout << "\n";
    std::cout << "╔═══════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║                        VALIDATION RESULTS                         ║\n";
    std::cout << "╠═══════════════════════════════════════════════════════════════════╣\n";
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "║  Max energy error:    " << std::setw(12) << maxEnergyError << " %";
    if (maxEnergyError < 0.01) std::cout << "  ✓ PASS"; else std::cout << "  ✗ FAIL";
    std::cout << "               ║\n";
    std::cout << "║  Max ang.mom. error:  " << std::setw(12) << maxAngMomError << " %";
    if (maxAngMomError < 0.01) std::cout << "  ✓ PASS"; else std::cout << "  ✗ FAIL";
    std::cout << "               ║\n";
    std::cout << "║  Final quaternion |q|: " << std::setw(11) << qNorm;
    if (std::abs(qNorm - 1.0) < 0.0001) std::cout << "   ✓ NORMALIZED"; else std::cout << "   ✗ DRIFT";
    std::cout << "          ║\n";
    std::cout << "╚═══════════════════════════════════════════════════════════════════╝\n";
    
    // ===== Save trajectory and visualize =====
    RigidBodySerializer serializer;
    serializer.SetContainerSize(config.containerHalfSize);
    serializer.AddBody(RigidBodyInfo(mass, halfA, halfB, halfC, "DodgerBlue", "TumblingBox"));
    
    auto cwd = std::filesystem::current_path();
    std::string filename = (cwd / "chapter15_tumbling_box.mml").string();
    
    // Save every 10th frame for smoother visualization
    auto result = serializer.SaveSingleBody(filename, simulator.GetHistory(), 0, config.timeStep, 10);
    
    if (result)
    {
        std::cout << "\nTrajectory saved to: " << filename << "\n";
        std::cout << "Launching visualizer...\n";
        
        std::vector<std::string> files = { filename };
        Visualizer::VisualizeRigidBodyTrajectory(files);
    }
    else
    {
        std::cerr << "Error saving trajectory: " << result.message << "\n";
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
/// @brief Demo 2: Symmetric top precession and nutation
/// @details Demonstrates gyroscopic effects with a symmetric top.
///          The top has I_xx = I_yy ≠ I_zz (axial symmetry).
///////////////////////////////////////////////////////////////////////////////////////////
void Chapter15_Demo_SymmetricTop()
{
    std::cout << "\n";
    std::cout << "╔═══════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║       CHAPTER 15 DEMO 2: Symmetric Top (Precession)               ║\n";
    std::cout << "╠═══════════════════════════════════════════════════════════════════╣\n";
    std::cout << "║  Top: 0.6m × 0.6m × 2m (tall cylinder-like shape)                 ║\n";
    std::cout << "║  Initial spin about Z-axis with slight tilt                       ║\n";
    std::cout << "║  Simulation: 15 seconds                                           ║\n";
    std::cout << "╚═══════════════════════════════════════════════════════════════════╝\n\n";
    
    // Symmetric top - square cross-section, elongated along Z
    Real mass = 5.0;
    Real halfA = 0.3;  // X half-width
    Real halfB = 0.3;  // Y half-width (same as X for symmetry)
    Real halfC = 1.0;  // Z half-height (elongated)
    
    RigidBodyBox top(mass, halfA, halfB, halfC);
    
    // Initial orientation: slight tilt from vertical
    Quaternion tilt = Quaternion::FromAxisAngle(Vec3Cart(1, 0, 0), 0.2);  // 0.2 rad tilt
    top.Position() = Vec3Cart(0, 0, 0);
    top.Velocity() = Vec3Cart(0, 0, 0);
    top.Orientation() = tilt;
    top.AngularVel() = Vec3Cart(0.3, 0.1, 5.0);  // Fast spin about Z, slow wobble
    
    // Print inertia
    auto I = top.InertiaTensorBody();
    std::cout << "Inertia tensor (symmetric: I_xx ≈ I_yy):\n";
    std::cout << "  I_xx = " << std::fixed << std::setprecision(4) << I(0,0) << " kg·m²\n";
    std::cout << "  I_yy = " << I(1,1) << " kg·m²\n";
    std::cout << "  I_zz = " << I(2,2) << " kg·m² (spin axis)\n\n";
    
    SimulationConfig config;
    config.timeStep = 0.001;
    config.totalTime = 15.0;
    config.containerHalfSize = 1000.0;  // No container
    
    RigidBodySimulator simulator(config);
    simulator.AddBody(top);
    
    Real initialEnergy = simulator.TotalKineticEnergy();
    std::cout << "Initial kinetic energy: " << initialEnergy << " J\n";
    std::cout << "Running simulation...\n" << std::flush;
    
    simulator.Run();
    
    Real finalEnergy = simulator.TotalKineticEnergy();
    Real energyError = std::abs(finalEnergy - initialEnergy) / initialEnergy * 100.0;
    
    std::cout << "Final kinetic energy:   " << finalEnergy << " J\n";
    std::cout << "Energy conservation:    " << std::setprecision(6) << energyError << " %";
    if (energyError < 0.01) std::cout << " ✓\n"; else std::cout << " ✗\n";
    
    // Save and visualize
    RigidBodySerializer serializer;
    serializer.SetContainerSize(config.containerHalfSize);
    serializer.AddBody(RigidBodyInfo(mass, halfA, halfB, halfC, "Gold", "SymmetricTop"));
    
    auto cwd = std::filesystem::current_path();
    std::string filename = (cwd / "chapter15_symmetric_top.mml").string();
    
    auto result = serializer.SaveSingleBody(filename, simulator.GetHistory(), 0, config.timeStep, 10);
    
    if (result)
    {
        std::cout << "\nTrajectory saved to: " << filename << "\n";
        std::cout << "Launching visualizer...\n";
        
        std::vector<std::string> files = { filename };
        Visualizer::VisualizeRigidBodyTrajectory(files);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
/// @brief Demo 3: Box bouncing in a container
/// @details Demonstrates rigid body dynamics with wall collisions.
///          Validates energy conservation through elastic bounces.
///////////////////////////////////////////////////////////////////////////////////////////
void Chapter15_Demo_BoxInContainer()
{
    std::cout << "\n";
    std::cout << "╔═══════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║       CHAPTER 15 DEMO 3: Box in Container (Collisions)            ║\n";
    std::cout << "╠═══════════════════════════════════════════════════════════════════╣\n";
    std::cout << "║  Box: 2m × 1.4m × 1m, mass = 8 kg                                 ║\n";
    std::cout << "║  Container: 8m × 8m × 8m cubic box                                ║\n";
    std::cout << "║  Elastic wall collisions (e = 1.0)                                ║\n";
    std::cout << "║  Simulation: 20 seconds                                           ║\n";
    std::cout << "╚═══════════════════════════════════════════════════════════════════╝\n\n";
    
    Real mass = 8.0;
    Real halfA = 1.0, halfB = 0.7, halfC = 0.5;
    
    RigidBodyBox box(mass, halfA, halfB, halfC);
    
    // Start off-center with velocity toward walls
    box.Position() = Vec3Cart(1.0, 0.5, -0.5);
    box.Velocity() = Vec3Cart(2.5, 1.5, -1.0);
    box.Orientation() = Quaternion::FromAxisAngle(Vec3Cart(1, 1, 0).GetAsUnitVector(), 0.3);
    box.AngularVel() = Vec3Cart(1.0, 0.5, 0.8);
    
    SimulationConfig config;
    config.timeStep = 0.001;
    config.totalTime = 20.0;
    config.containerHalfSize = 4.0;   // 8m cube container
    config.coeffRestitution = 1.0;    // Perfectly elastic
    
    RigidBodySimulator simulator(config);
    simulator.AddBody(box);
    
    Real initialEnergy = simulator.TotalKineticEnergy();
    Real lastReportTime = 0.0;
    int collisionCount = 0;
    
    std::cout << "Initial kinetic energy: " << std::fixed << std::setprecision(4) 
              << initialEnergy << " J\n";
    std::cout << "Running simulation with wall collisions...\n\n";
    
    simulator.SetProgressCallback([&](Real t, const RigidBodySimulator::BodyContainer& bodies)
    {
        if (t - lastReportTime >= 2.0)
        {
            Real energy = bodies[0]->KineticEnergy();
            Real error = std::abs(energy - initialEnergy) / initialEnergy * 100.0;
            std::cout << "t = " << std::setw(5) << std::setprecision(1) << t 
                      << "s  |  E = " << std::setprecision(4) << energy 
                      << " J  |  ΔE = " << std::setprecision(4) << error << "%\n";
            lastReportTime = t;
        }
    });
    
    simulator.Run();
    
    Real finalEnergy = simulator.TotalKineticEnergy();
    Real energyError = std::abs(finalEnergy - initialEnergy) / initialEnergy * 100.0;
    
    std::cout << "\n";
    std::cout << "╔═══════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║                           RESULTS                                 ║\n";
    std::cout << "╠═══════════════════════════════════════════════════════════════════╣\n";
    std::cout << "║  Initial energy:  " << std::setw(12) << std::setprecision(6) << initialEnergy << " J                            ║\n";
    std::cout << "║  Final energy:    " << std::setw(12) << finalEnergy << " J                            ║\n";
    std::cout << "║  Energy error:    " << std::setw(12) << energyError << " %";
    if (energyError < 1.0) std::cout << "  ✓ CONSERVED"; else std::cout << "  ✗ DRIFT";
    std::cout << "              ║\n";
    std::cout << "╚═══════════════════════════════════════════════════════════════════╝\n";
    
    // Save and visualize
    RigidBodySerializer serializer;
    serializer.SetContainerSize(config.containerHalfSize);
    serializer.AddBody(RigidBodyInfo(mass, halfA, halfB, halfC, "Crimson", "BouncingBox"));
    
    auto cwd = std::filesystem::current_path();
    std::string filename = (cwd / "chapter15_box_in_container.mml").string();
    
    auto result = serializer.SaveSingleBody(filename, simulator.GetHistory(), 0, config.timeStep, 10);
    
    if (result)
    {
        std::cout << "\nTrajectory saved to: " << filename << "\n";
        std::cout << "Launching visualizer...\n";
        
        std::vector<std::string> files = { filename };
        Visualizer::VisualizeRigidBodyTrajectory(files);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
/// @brief Demo 4: Two Boxes and a Sphere - Mixed Shape Collision Chaos
/// @details Demonstrates the new multi-shape rigid body system with:
///          - Box-Box collisions (SAT algorithm)
///          - Sphere-Box collisions (closest-point-on-OBB)
///          - Sphere-Wall collisions (O(1) distance-based)
///          - Box-Wall collisions (vertex-based)
///          All shapes use the unified collision dispatcher.
///
/// The dream scenario: "Two tumbling boxes with one sphere to wreak havoc!"
///////////////////////////////////////////////////////////////////////////////////////////
void Chapter15_Demo_MixedShapes()
{
    std::cout << "\n";
    std::cout << "╔═══════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║    CHAPTER 15 DEMO 4: Two Boxes + Sphere (Multi-Shape Chaos!)     ║\n";
    std::cout << "╠═══════════════════════════════════════════════════════════════════╣\n";
    std::cout << "║  Box 1: 2m × 1m × 0.6m, 10 kg (red)                               ║\n";
    std::cout << "║  Box 2: 1.5m × 1m × 0.8m, 8 kg (blue)                             ║\n";
    std::cout << "║  Sphere: radius 0.6m, 5 kg (green) - the havoc maker!             ║\n";
    std::cout << "║  Container: 10m × 10m × 10m cubic box                             ║\n";
    std::cout << "║  Elastic collisions (e = 1.0)                                     ║\n";
    std::cout << "║  Simulation: 20 seconds                                           ║\n";
    std::cout << "╚═══════════════════════════════════════════════════════════════════╝\n\n";
    
    // ===== Create Box 1 (red) =====
    Real mass1 = 10.0;
    Real half1a = 1.0, half1b = 0.5, half1c = 0.3;  // 2m × 1m × 0.6m
    RigidBodyBox box1(mass1, half1a, half1b, half1c);
    
    box1.Position() = Vec3Cart(-2.0, 0.5, 0.0);
    box1.Velocity() = Vec3Cart(3.0, 0.5, 0.2);
    box1.Orientation() = Quaternion::FromAxisAngle(Vec3Cart(1, 0, 0), 0.3);
    box1.AngularVel() = Vec3Cart(0.3, 0.5, 0.2);
    
    // ===== Create Box 2 (blue) =====
    Real mass2 = 8.0;
    Real half2a = 0.75, half2b = 0.5, half2c = 0.4;  // 1.5m × 1m × 0.8m
    RigidBodyBox box2(mass2, half2a, half2b, half2c);
    
    box2.Position() = Vec3Cart(2.0, -0.5, 0.3);
    box2.Velocity() = Vec3Cart(-2.0, 0.3, -0.1);
    box2.Orientation() = Quaternion::FromAxisAngle(Vec3Cart(0, 1, 0), 0.5);
    box2.AngularVel() = Vec3Cart(-0.2, 0.4, -0.3);
    
    // ===== Create Sphere (green) - THE HAVOC MAKER! =====
    Real massSphere = 5.0;
    Real radius = 0.156;  // 30% bigger than 0.12 - chaos ball!
    RigidBodySphere sphere(massSphere, radius);  // Sphere constructor!
    
    // Start sphere from above, coming down FAST to cause maximum chaos
    sphere.Position() = Vec3Cart(0.0, 0.0, 2.0);
    sphere.Velocity() = Vec3Cart(3.0, -2.0, -6.0);  // DOUBLE velocity - maximum havoc!
    sphere.Orientation() = Quaternion::Identity();
    sphere.AngularVel() = Vec3Cart(2.0, 1.5, 1.0);  // Fast spin
    
    // Print shape information
    std::cout << "Shape Configuration:\n";
    std::cout << "────────────────────\n";
    std::cout << "  Box 1:   " << (IsBox(box1) ? "BOX" : "???") 
              << " - half-extents (" << half1a << ", " << half1b << ", " << half1c << ")\n";
    std::cout << "  Box 2:   " << (IsBox(box2) ? "BOX" : "???") 
              << " - half-extents (" << half2a << ", " << half2b << ", " << half2c << ")\n";
    std::cout << "  Sphere:  " << (IsSphere(sphere) ? "SPHERE" : "???") 
              << " - radius " << radius << "\n\n";
    
    // Print inertia tensors
    std::cout << "Inertia Tensors (body frame, diagonal):\n";
    std::cout << "────────────────────────────────────────\n";
    auto I1 = box1.InertiaTensorBody();
    auto I2 = box2.InertiaTensorBody();
    auto ISph = sphere.InertiaTensorBody();
    
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "  Box 1:   I = (" << I1(0,0) << ", " << I1(1,1) << ", " << I1(2,2) << ") kg·m²\n";
    std::cout << "  Box 2:   I = (" << I2(0,0) << ", " << I2(1,1) << ", " << I2(2,2) << ") kg·m²\n";
    std::cout << "  Sphere:  I = " << ISph(0,0) << " kg·m² (isotropic)\n\n";
    
    // ===== Configure Simulator =====
    SimulationConfig config;
    config.timeStep = 0.001;
    config.totalTime = 20.0;
    config.containerHalfSize = 5.0;   // 10m cube container
    config.coeffRestitution = 1.0;    // Perfectly elastic
    config.maxCollisionIterations = 10;
    
    RigidBodySimulator simulator(config);
    simulator.AddBody(box1);
    simulator.AddBody(box2);
    simulator.AddBody(sphere);
    
    Real initialEnergy = simulator.TotalKineticEnergy();
    Vec3Cart initialMomentum = simulator.TotalLinearMomentum();
    
    std::cout << "Initial System State:\n";
    std::cout << "─────────────────────\n";
    std::cout << "  Total kinetic energy: " << initialEnergy << " J\n";
    std::cout << "  Total linear momentum: (" << initialMomentum.X() << ", " 
              << initialMomentum.Y() << ", " << initialMomentum.Z() << ") kg·m/s\n\n";
    
    Real lastReportTime = 0.0;
    Real maxEnergyError = 0.0;
    
    std::cout << "Running simulation with mixed-shape collisions...\n";
    std::cout << "═══════════════════════════════════════════════\n";
    
    simulator.SetProgressCallback([&](Real t, const RigidBodySimulator::BodyContainer& bodies)
    {
        Real energy = 0.0;
        for (const auto& body : bodies)
            energy += body->KineticEnergy();
        
        Real error = std::abs(energy - initialEnergy) / initialEnergy * 100.0;
        maxEnergyError = std::max(maxEnergyError, error);
        
        if (t - lastReportTime >= 1.0)
        {
            std::cout << "t = " << std::setw(5) << std::setprecision(1) << t 
                      << "s  |  E = " << std::setprecision(4) << energy 
                      << " J  |  ΔE = " << std::setprecision(6) << error << "%\n";
            lastReportTime = t;
        }
    });
    
    simulator.Run();
    
    // ===== Final Validation =====
    Real finalEnergy = simulator.TotalKineticEnergy();
    Real energyError = std::abs(finalEnergy - initialEnergy) / initialEnergy * 100.0;
    
    std::cout << "\n";
    std::cout << "╔═══════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║                    CONSERVATION CHECK                             ║\n";
    std::cout << "╠═══════════════════════════════════════════════════════════════════╣\n";
    std::cout << "║  Initial energy:    " << std::setw(12) << std::setprecision(6) << initialEnergy << " J                          ║\n";
    std::cout << "║  Final energy:      " << std::setw(12) << finalEnergy << " J                          ║\n";
    std::cout << "║  Max energy error:  " << std::setw(12) << maxEnergyError << " %";
    if (maxEnergyError < 1.0) std::cout << "  ✓ CONSERVED"; else std::cout << "  ✗ DRIFT";
    std::cout << "            ║\n";
    std::cout << "╚═══════════════════════════════════════════════════════════════════╝\n";
    
    // ===== Save and Visualize =====
    // The WPF visualizer expects RIGID_BODY_TRAJECTORY_3D format (single body per file).
    // Save each body separately, then load all files into the visualizer.
    RigidBodySerializer serializer;
    serializer.SetContainerSize(config.containerHalfSize);
    serializer.AddBody(RigidBodyInfo::Box(mass1, half1a, half1b, half1c, "Crimson", "Box1"));
    serializer.AddBody(RigidBodyInfo::Box(mass2, half2a, half2b, half2c, "RoyalBlue", "Box2"));
    serializer.AddBody(RigidBodyInfo::Sphere(massSphere, radius, "LimeGreen", "Sphere"));
    
    auto cwd = std::filesystem::current_path();
    std::string file1 = (cwd / "chapter15_mixed_box1.mml").string();
    std::string file2 = (cwd / "chapter15_mixed_box2.mml").string();
    std::string file3 = (cwd / "chapter15_mixed_sphere.mml").string();
    
    // Save each body to separate file (RIGID_BODY_TRAJECTORY_3D format)
    auto result1 = serializer.SaveSingleBody(file1, simulator.GetHistory(), 0, config.timeStep, 10);
    auto result2 = serializer.SaveSingleBody(file2, simulator.GetHistory(), 1, config.timeStep, 10);
    auto result3 = serializer.SaveSingleBody(file3, simulator.GetHistory(), 2, config.timeStep, 10);
    
    if (result1 && result2 && result3)
    {
        std::cout << "\nTrajectories saved to:\n";
        std::cout << "  - " << file1 << "\n";
        std::cout << "  - " << file2 << "\n";
        std::cout << "  - " << file3 << " (sphere as cube)\n";
        std::cout << "Launching visualizer...\n";
        
        // Load all three files into visualizer
        std::vector<std::string> files = { file1, file2, file3 };
        Visualizer::VisualizeRigidBodyTrajectory(files);
    }
    else
    {
        if (!result1) std::cerr << "Error saving Box1: " << result1.message << "\n";
        if (!result2) std::cerr << "Error saving Box2: " << result2.message << "\n";
        if (!result3) std::cerr << "Error saving Sphere: " << result3.message << "\n";
    }
}


// ===== Entry points for chapter coordinator =====

void Chapter15_Parallelopid_rotation_no_gravity()
{
    Chapter15_Demo_TumblingBox();
}

void Chapter15_Symmerical_top_simulation()
{
    Chapter15_Demo_SymmetricTop();
}

void Chapter15_Box_in_container()
{
    Chapter15_Demo_BoxInContainer();
}

void Chapter15_Mixed_shapes_chaos()
{
    Chapter15_Demo_MixedShapes();
}

void Chapter15_Rigid_body()
{
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                   CHAPTER 15 - Rigid body                       ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    Chapter15_Parallelopid_rotation_no_gravity();
    Chapter15_Symmerical_top_simulation();
    Chapter15_Box_in_container();
    Chapter15_Mixed_shapes_chaos();
}


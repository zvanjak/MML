/**
 * @file show_readme_files.cpp
 * @brief Opens all README visualization data files using appropriate visualizers
 * 
 * Contains one sample file for each visualization type supported by MML.
 * Files are loaded from src/visualization_examples/visualization_data/readme/
 * and opened sequentially with the matching visualizer.
 * 
 * All FromFile/file-based Visualizer methods accept full paths directly,
 * so no copying to the results folder is needed.
 */

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#include "tools/Visualizer.h"
#endif

#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

using namespace MML;

namespace {

    /// @brief Find the readme data folder, searching common locations
    /// @return Path to the readme folder, or empty if not found
    std::filesystem::path FindReadmeDataFolder()
    {
        std::string dataFolder = "src/visualization_examples/visualization_data/readme/";

        // Try relative to current working directory (running from project root)
        if (std::filesystem::exists(dataFolder))
            return std::filesystem::absolute(dataFolder);

        // Try relative to parent directories (running from build folder)
        std::vector<std::string> prefixes = {"../", "../../", "../../../"};
        for (const auto& prefix : prefixes) {
            std::string path = prefix + dataFolder;
            if (std::filesystem::exists(path))
                return std::filesystem::absolute(path);
        }

        return {};
    }

    /// @brief Get the full (absolute) path to a file in the readme folder
    std::string FullPath(const std::filesystem::path& folder, const std::string& fileName)
    {
        return (folder / fileName).string();
    }

    /// @brief Show a single readme file using its visualizer, with status output
    void ShowFile(const std::filesystem::path& folder, const std::string& fileName,
                  const std::string& description,
                  std::function<VisualizerResult(const std::string&)> visualize)
    {
        std::cout << "  " << description << " (" << fileName << ")...\n";

        std::string path = FullPath(folder, fileName);
        if (!std::filesystem::exists(path)) {
            std::cerr << "    File not found: " << path << "\n";
            return;
        }

        auto result = visualize(path);
        if (!result.success)
            std::cerr << "    Visualizer error: " << result.errorMessage << "\n";
    }

} // anonymous namespace

void Show_Readme_Files()
{
    std::cout << "\n=== README Visualization Files ===\n\n";
    std::cout << "Opening one sample file for each visualization type.\n";
    std::cout << "Close each visualizer window to proceed to the next.\n\n";

    auto folder = FindReadmeDataFolder();
    if (folder.empty()) {
        std::cerr << "ERROR: Could not find readme data folder.\n";
        std::cerr << "Expected: src/visualization_examples/visualization_data/readme/\n";
        return;
    }

    std::cout << "Data folder: " << folder << "\n\n";

    // --- Real functions (MULTI_REAL_FUNCTION) ---
    std::cout << "[1/11] Real Functions\n";
    ShowFile(folder, "multireal_func_damped.mml",
             "Damped oscillation (multi real function)",
             [](const std::string& f) { return Visualizer::VisualizeRealFunctionFromFile(f); });

    ShowFile(folder, "multireal_func_lorenz_time_series.mml",
             "Lorenz time series (multi real function)",
             [](const std::string& f) { return Visualizer::VisualizeRealFunctionFromFile(f); });

    // --- Scalar function 2D ---
    std::cout << "\n[2/11] Scalar Function 2D\n";
    ShowFile(folder, "scalar_func_2d_ripple.mml",
             "Ripple pattern (2D scalar field as surface)",
             [](const std::string& f) { return Visualizer::VisualizeScalarFunc2DFromFile(f); });

    // --- Scalar function 3D ---
    std::cout << "\n[3/11] Scalar Function 3D\n";
    ShowFile(folder, "scalar_func_3d_gyroid_3d.mml",
             "Gyroid (3D scalar field / isosurface)",
             [](const std::string& f) { return Visualizer::VisualizeScalarFunc3DFromFile(f); });

    // --- Parametric curve 2D ---
    std::cout << "\n[4/11] Parametric Curve 2D\n";
    ShowFile(folder, "param_curve_2d_butterfly.mml",
             "Butterfly curve (2D parametric)",
             [](const std::string& f) { return Visualizer::VisualizeMultiParamCurve2D({f}); });

    // --- Parametric curve 3D ---
    std::cout << "\n[5/11] Parametric Curve 3D\n";
    ShowFile(folder, "param_curve_3d_trefoil.mml",
             "Trefoil knot (3D parametric)",
             [](const std::string& f) { return Visualizer::VisualizeMultiParamCurve3D({f}); });

    // --- Parametric surface ---
    std::cout << "\n[6/11] Parametric Surface\n";
    ShowFile(folder, "parametric_surface_klein_bottle.mml",
             "Klein bottle (parametric surface)",
             [](const std::string& f) { return Visualizer::VisualizeParametricSurfaceFromFile(f); });

    // --- Vector field 2D ---
    std::cout << "\n[7/11] Vector Field 2D\n";
    ShowFile(folder, "vector_field_2d.mml",
             "2D vector field",
             [](const std::string& f) { return Visualizer::VisualizeVectorField2DFromFile(f); });

    // --- Vector field 3D ---
    std::cout << "\n[8/11] Vector Field 3D\n";
    ShowFile(folder, "vector_field_3d.mml",
             "3D vector field",
             [](const std::string& f) { return Visualizer::VisualizeVectorField3DFromFile(f); });

    // --- Particle simulation 2D ---
    std::cout << "\n[9/11] Particle Simulation 2D\n";
    ShowFile(folder, "particle_visualizer_2d_1.mml",
             "2D particle collision simulation",
             [](const std::string& f) { return Visualizer::VisualizeParticleSimulation2D(f); });

    // --- Particle simulation 3D ---
    std::cout << "\n[10/11] Particle Simulation 3D\n";
    ShowFile(folder, "particle_visualizer_3d_1.mml",
             "3D particle simulation",
             [](const std::string& f) { return Visualizer::VisualizeParticleSimulation3D(f); });

    // --- Rigid body trajectories (3 files together) ---
    std::cout << "\n[11/11] Rigid Body Trajectories\n";
    std::cout << "  Mixed rigid body simulation (box1 + box2 + sphere)...\n";
    {
        std::vector<std::string> rigidFiles = {
            FullPath(folder, "rigid_body_mixed_box1.mml"),
            FullPath(folder, "rigid_body_mixed_box2.mml"),
            FullPath(folder, "rigid_body_mixed_sphere.mml")
        };

        // Check all files exist
        bool allExist = true;
        for (const auto& rf : rigidFiles) {
            if (!std::filesystem::exists(rf)) {
                std::cerr << "    File not found: " << rf << "\n";
                allExist = false;
            }
        }

        if (allExist) {
            auto result = Visualizer::VisualizeRigidBodyTrajectory(rigidFiles);
            if (!result.success)
                std::cerr << "    Visualizer error: " << result.errorMessage << "\n";
        }
    }

    std::cout << "\nAll README visualization files shown.\n";
}
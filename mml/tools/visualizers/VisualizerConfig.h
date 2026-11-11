///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        tools/visualizers/VisualizerConfig.h                                                  ///
///  Description: Visualization utilities for functions, curves, fields               ///
///               Cross-platform support: Windows (WPF/Qt/FLTK),                      ///
///               Linux (Qt/FLTK), macOS (Qt/FLTK)                                    ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_VISUALIZER_CONFIG_H
#define MML_VISUALIZER_CONFIG_H

// Standard headers - include what we use
#include <algorithm>
#include <array>
#include <cctype>
#include <cstddef>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>

// Cross-platform executable extension - defined inside MML namespace (see below)

// Helper to get environment variable without compiler warnings
// Note: C++20 still doesn't provide a standard alternative to std::getenv()
namespace MML::detail::visualizer {
inline const char *GetEnv(const char *name) noexcept {
#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable : 4996)
#elif defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#elif defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif

  const char *result = std::getenv(name);

#if defined(_MSC_VER)
#pragma warning(pop)
#elif defined(__clang__)
#pragma clang diagnostic pop
#elif defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

  return result;
}
}

#include <mml/MMLExceptions.h>

namespace MML
{

#ifdef _WIN32
  static constexpr const char *EXECUTABLE_EXT = ".exe";
#else
  static constexpr const char *EXECUTABLE_EXT = "";
#endif

  ///////////////////////////////////////////////////////////////////////////////
  //                    LINUX VISUALIZER BACKEND SELECTION
  ///////////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////////
  //          PROJECT PATH RESOLUTION (REQUIRED BY BACKEND SYSTEM)
  ///////////////////////////////////////////////////////////////////////////////

  // Get project root path from environment or use current directory as fallback
  inline std::string GetProjectPath() {
    // 1. Try environment variable
    if (const char *envPath = detail::visualizer::GetEnv("MML_PROJECT_PATH")) {
      return std::string(envPath);
    }

    // Result-file placement is independent of visualizer bundle configuration.
    return std::filesystem::current_path().string();
  }

  namespace detail {
    inline std::filesystem::path& GetVisualizerRootState() {
      static std::filesystem::path root;
      return root;
    }
  }

  /// @brief Set the bundle directory containing win/, linux/, and mac/.
  /// @details An empty path clears the runtime override. Relative paths are
  ///          resolved when set. Configure before launching from other threads.
  inline void SetVisualizerRoot(const std::filesystem::path& root) {
    if (root.native().find(std::filesystem::path::value_type{}) != std::filesystem::path::string_type::npos)
      throw VisualizerError("Visualizer root contains an embedded null character");
    detail::GetVisualizerRootState() = root.empty() ? root : std::filesystem::absolute(root).lexically_normal();
  }

  /// @brief Resolve runtime override, MML_VISUALIZER_ROOT, then cwd/visualizers.
  inline std::string GetVisualizerRoot() {
    const auto& root = detail::GetVisualizerRootState();
    if (!root.empty())
      return root.string();
    if (const char* envRoot = detail::visualizer::GetEnv("MML_VISUALIZER_ROOT"); envRoot && *envRoot)
      return std::filesystem::absolute(envRoot).lexically_normal().string();
    return (std::filesystem::current_path() / "visualizers").string();
  }

  ///////////////////////////////////////////////////////////////////////////////
  //          UNIFIED CROSS-PLATFORM VISUALIZER BACKEND SYSTEM
  ///////////////////////////////////////////////////////////////////////////////

  /// Visualizer backend options (Windows and Linux)
  /// Enables runtime selection of visualization backend on both platforms
  enum class VisualizerBackend {
    WPF,  ///< Windows Presentation Foundation (Windows only, highest quality)
    Qt,   ///< Qt5/Qt6 (Cross-platform, feature-rich)
    FLTK, ///< Fast Light Toolkit (Cross-platform, lightweight)
    Auto  ///< Platform-specific default (WPF on Windows, Qt on Linux)
  };

  namespace detail::visualizer {
    enum class Kind {
      RealFunction, ScalarFunction2D, Curve2D, Curve3D, VectorField2D, VectorField3D,
      Particles2D, Particles3D, ParametricSurface, ScalarFunction3D, RigidBody, World
    };

    inline constexpr std::array<const char*, 12> ExecutableNames{
      "MML_RealFunctionVisualizer", "MML_ScalarFunction2D_Visualizer",
      "MML_ParametricCurve2D_Visualizer", "MML_ParametricCurve3D_Visualizer",
      "MML_VectorField2D_Visualizer", "MML_VectorField3D_Visualizer",
      "MML_ParticleVisualizer2D", "MML_ParticleVisualizer3D",
      "MML_ParametricSurface_Visualizer", "MML_ScalarFunction3D_Visualizer",
      "MML_RigidBodyMovement_Visualizer", "MML_WorldVisualizer"
    };

    inline const char* ExecutableName(Kind kind) {
      const auto index = static_cast<std::size_t>(kind);
      if (index >= ExecutableNames.size())
        throw VisualizerError("Invalid visualizer kind");
      return ExecutableNames[index];
    }
  }

  namespace detail {
    /// Global state for unified visualizer backend selection
    inline VisualizerBackend& GetUnifiedBackendState() {
      static VisualizerBackend backend = VisualizerBackend::Auto;
      return backend;
    }
  }

  inline std::string GetVisualizerPath(const std::string& visualizerName, VisualizerBackend backend);
  inline VisualizerBackend GetVisualizerBackend();

  /// Check if at least one supported executable is installed for a backend.
  /// @param backend The backend to check (WPF, Qt, FLTK, or Auto)
  /// @return true if the backend is available, false otherwise
  inline bool IsBackendAvailable(VisualizerBackend backend) {
    if (backend == VisualizerBackend::Auto)
      backend = GetVisualizerBackend();
#ifndef _WIN32
    if (backend == VisualizerBackend::WPF)
      return false;
#endif
    if (backend != VisualizerBackend::WPF && backend != VisualizerBackend::Qt && backend != VisualizerBackend::FLTK)
      return false;
    for (const auto* name : detail::visualizer::ExecutableNames)
      if (std::filesystem::is_regular_file(GetVisualizerPath(name, backend)))
        return true;
    return false;
  }

  /// Set the visualizer backend at runtime
  /// @param backend The backend to use (WPF, Qt, FLTK, or Auto)
  /// @note If the selected backend is not available, defaults to Auto
  inline void SetVisualizerBackend(VisualizerBackend backend) {
    if (backend != VisualizerBackend::Auto && !IsBackendAvailable(backend)) {
      std::cerr << "Warning: Selected visualizer backend not available. "
                << "Using Auto selection.\n";
      backend = VisualizerBackend::Auto;
    }
    detail::GetUnifiedBackendState() = backend;
  }

  /// Get the currently selected visualizer backend
  /// Checks (in order):
  ///   1. Environment variable MML_VISUALIZER_BACKEND
  ///   2. Runtime API setting via SetVisualizerBackend()
  ///   3. Platform-specific default (WPF on Windows, Qt on Linux)
  ///
  /// Environment variable examples:
  ///   Windows: set MML_VISUALIZER_BACKEND=Qt
  ///   Linux:   export MML_VISUALIZER_BACKEND=FLTK
  ///
  /// @return Current backend selection (never returns Auto, always resolves to concrete backend)
  inline VisualizerBackend GetVisualizerBackend() {
    bool environmentAuto = false;
    // 1. Check environment variable for override
    if (const char* envBackend = detail::visualizer::GetEnv("MML_VISUALIZER_BACKEND")) {
      std::string backend(envBackend);
      // Convert to uppercase for case-insensitive comparison
      std::transform(backend.begin(), backend.end(), backend.begin(),
                     [](unsigned char c) { return static_cast<char>(std::toupper(c)); });

      if (backend == "WPF")  return VisualizerBackend::WPF;
      if (backend == "QT" || backend == "QT5" || backend == "QT6")
        return VisualizerBackend::Qt;
      if (backend == "FLTK") return VisualizerBackend::FLTK;
      environmentAuto = backend == "AUTO";
    }

    // 2. Get runtime setting
    VisualizerBackend backend = environmentAuto ? VisualizerBackend::Auto : detail::GetUnifiedBackendState();

    // 3. Resolve Auto to platform default
    if (backend == VisualizerBackend::Auto) {
#ifdef _WIN32
      return VisualizerBackend::WPF;  // Default to WPF on Windows
#elif defined(__linux__)
      return VisualizerBackend::Qt;   // Default to Qt on Linux
#elif defined(__APPLE__)
      return VisualizerBackend::Qt;   // Default to Qt on macOS
#else
      return VisualizerBackend::FLTK; // FLTK for other platforms
#endif
    }

    return backend;
  }

  // WPF bundle layout: win/WPF/<directory>/<executable>.exe.
  inline std::string GetWPFVisualizerPath(const std::string &name,
                                          const std::string &wpfDir) {
  #ifdef _WIN32
    return (std::filesystem::path(GetVisualizerRoot()) / "win" / "WPF" / wpfDir / (name + ".exe")).string();
  #else
    return ""; // Not available on non-Windows
  #endif
  }

  // Global paths with lazy evaluation
  inline std::string GetGlobalPath() {
    static std::string path = GetProjectPath();
    return path;
  }

  inline std::string GetResultFilesPath() {
    return (std::filesystem::path(GetGlobalPath()) / "results").string() + std::string(1, std::filesystem::path::preferred_separator);
  }

  ///////////////////////////////////////////////////////////////////////////////
  //          GENERIC VISUALIZER PATH DISPATCHER
  ///////////////////////////////////////////////////////////////////////////////

  /// Helper: Get WPF directory name for a visualizer
  /// New WPF structure uses MML_ prefixed folder names matching the executable
  /// Example: "MML_RealFunctionVisualizer" -> "MML_RealFunctionVisualizer"
  inline std::string ConvertToWPFDirName(const std::string& visualizerName) {
    // New WPF folder structure: folder name matches executable name (MML_ prefix)
    return visualizerName;
  }

  /// Generic visualizer path resolution with backend dispatch
  /// Constructs the full path to a visualizer executable based on:
  ///   - Visualizer name (e.g., "MML_RealFunctionVisualizer")
  ///   - Selected backend (from GetVisualizerBackend())
  ///   - Current platform (Windows/Linux)
  ///
  /// @param visualizerName Base name of the visualizer executable
  /// @param backend Backend to use (if not specified, uses GetVisualizerBackend())
  /// @return Full path to the visualizer executable
  inline std::string GetVisualizerPath(const std::string& visualizerName,
                                       VisualizerBackend backend = VisualizerBackend::Auto) {
    std::filesystem::path root(GetVisualizerRoot());

    // Resolve Auto to actual backend
    if (backend == VisualizerBackend::Auto) {
      backend = GetVisualizerBackend();
    }

#ifdef _WIN32
    if (backend == VisualizerBackend::WPF) {
      // WPF has one subdirectory per application.
      std::string wpfDir = ConvertToWPFDirName(visualizerName);
      return (root / "win" / "WPF" / wpfDir / (visualizerName + ".exe")).string();
    } else if (backend == VisualizerBackend::Qt) {
      // Qt flat structure: all visualizers share Qt DLLs in one folder
      return (root / "win" / "Qt" / (visualizerName + ".exe")).string();
    } else if (backend == VisualizerBackend::FLTK) {
      // FLTK uses a flat directory with suffixed executable names.
      return (root / "win" / "FLTK" / (visualizerName + "_FLTK.exe")).string();
    }
#elif defined(__linux__)
    if (backend == VisualizerBackend::Qt) {
      return (root / "linux" / "Qt" / visualizerName).string();
    } else if (backend == VisualizerBackend::FLTK) {
      return (root / "linux" / "FLTK" / (visualizerName + "_FLTK")).string();
    }
#elif defined(__APPLE__)
    if (backend == VisualizerBackend::Qt) {
      // Qt apps are .app bundles on macOS - executable is inside Contents/MacOS/
      return (root / "mac" / "Qt" / (visualizerName + ".app") / "Contents" / "MacOS" / visualizerName).string();
    } else if (backend == VisualizerBackend::FLTK) {
      return (root / "mac" / "FLTK" / (visualizerName + "_FLTK")).string();
    }
#endif

    // Should never reach here if backend selection is working correctly
    throw VisualizerError("Visualizer backend not supported on this platform: " +
                             visualizerName);
  }

  ///////////////////////////////////////////////////////////////////////////////
  //          SPECIFIC VISUALIZER PATH FUNCTIONS (USING DISPATCHER)
  ///////////////////////////////////////////////////////////////////////////////

  // Visualizer paths - cross-platform support
  namespace detail::visualizer {
    inline std::string ExecutablePath(Kind kind) {
      return GetVisualizerPath(ExecutableName(kind));
    }
  }

  inline std::string GetRealFuncVisualizerPath() {
    return detail::visualizer::ExecutablePath(detail::visualizer::Kind::RealFunction);
  }

  inline std::string GetSurfaceVisualizerPath() {
    return detail::visualizer::ExecutablePath(detail::visualizer::Kind::ScalarFunction2D);
  }

  inline std::string GetParametricCurve2DVisualizerPath() {
    return detail::visualizer::ExecutablePath(detail::visualizer::Kind::Curve2D);
  }

  inline std::string GetParametricCurve3DVisualizerPath() {
    return detail::visualizer::ExecutablePath(detail::visualizer::Kind::Curve3D);
  }

  inline std::string GetVectorField2DVisualizerPath() {
    return detail::visualizer::ExecutablePath(detail::visualizer::Kind::VectorField2D);
  }

  inline std::string GetVectorField3DVisualizerPath() {
    return detail::visualizer::ExecutablePath(detail::visualizer::Kind::VectorField3D);
  }

  inline std::string GetParticle2DVisualizerPath() {
    return detail::visualizer::ExecutablePath(detail::visualizer::Kind::Particles2D);
  }

  inline std::string GetParticle3DVisualizerPath() {
    return detail::visualizer::ExecutablePath(detail::visualizer::Kind::Particles3D);
  }

  inline std::string GetParametricSurfaceVisualizerPath() {
    return detail::visualizer::ExecutablePath(detail::visualizer::Kind::ParametricSurface);
  }

  inline std::string GetScalarFunction3DVisualizerPath() {
    return detail::visualizer::ExecutablePath(detail::visualizer::Kind::ScalarFunction3D);
  }

  inline std::string GetRigidBodyVisualizerPath() {
    return detail::visualizer::ExecutablePath(detail::visualizer::Kind::RigidBody);
  }

  inline std::string GetWorldVisualizerPath() {
    return detail::visualizer::ExecutablePath(detail::visualizer::Kind::World);
  }

} // namespace MML

#endif // MML_VISUALIZER_CONFIG_H

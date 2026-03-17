///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        MMLVisualizators.h                                                  ///
///  Description: Visualization utilities for functions, curves, fields               ///
///               Cross-platform support: Windows (WPF/Qt/FLTK),                      ///
///               Linux (Qt/FLTK), macOS (Qt/FLTK)                                    ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_VISUALIZATORS_H
#define MML_VISUALIZATORS_H

// Standard headers - include what we use
#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>

// Cross-platform executable extension
#ifdef _WIN32
static constexpr const char *EXECUTABLE_EXT = ".exe";
#else
static constexpr const char *EXECUTABLE_EXT = "";
#endif

// Helper to get environment variable without compiler warnings
// Note: C++20 still doesn't provide a standard alternative to std::getenv()
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

#include "MMLExceptions.h"

namespace MML
{

  ///////////////////////////////////////////////////////////////////////////////
  //                    LINUX VISUALIZER BACKEND SELECTION
  ///////////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////////
  //          PROJECT PATH RESOLUTION (REQUIRED BY BACKEND SYSTEM)
  ///////////////////////////////////////////////////////////////////////////////

  // Get project root path from environment or use current directory as fallback
  inline std::string GetProjectPath() {
    // 1. Try environment variable
    if (const char *envPath = GetEnv("MML_PROJECT_PATH")) {
      return std::string(envPath);
    }

    // 2. Try to find project root by looking for characteristic directories
    std::filesystem::path current = std::filesystem::current_path();
    std::filesystem::path prev;
    while (current != prev && current.has_parent_path()) {
      if (std::filesystem::exists(current / "mml") &&
          std::filesystem::exists(current / "tools")) {
        return current.string();
      }
      prev = current;
      current = current.parent_path();
    }

    // 3. Fallback to current directory
    return std::filesystem::current_path().string();
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

  namespace detail {
    /// Global state for unified visualizer backend selection
    inline VisualizerBackend& GetUnifiedBackendState() {
      static VisualizerBackend backend = VisualizerBackend::Auto;
      return backend;
    }
  }

  /// Check if a visualizer backend is available on the current platform
  /// @param backend The backend to check (WPF, Qt, FLTK, or Auto)
  /// @return true if the backend is available, false otherwise
  inline bool IsBackendAvailable(VisualizerBackend backend) {
    std::string projectPath = GetProjectPath();
    
#ifdef _WIN32
    if (backend == VisualizerBackend::WPF) {
      return std::filesystem::exists(std::filesystem::path(projectPath) / "tools" / "visualizers" / "win" / "WPF");
    } else if (backend == VisualizerBackend::Qt) {
      return std::filesystem::exists(std::filesystem::path(projectPath) / "tools" / "visualizers" / "win" / "Qt");
    } else if (backend == VisualizerBackend::FLTK) {
      return std::filesystem::exists(std::filesystem::path(projectPath) / "tools" / "visualizers" / "win" / "FLTK");
    } else if (backend == VisualizerBackend::Auto) {
      return true;  // Auto always available (selects default)
    }
#elif defined(__linux__)
    if (backend == VisualizerBackend::WPF) {
      return false;  // WPF is Windows-only
    } else if (backend == VisualizerBackend::Qt) {
      return std::filesystem::exists(std::filesystem::path(projectPath) / "tools" / "visualizers" / "linux" / "Qt");
    } else if (backend == VisualizerBackend::FLTK) {
      return std::filesystem::exists(std::filesystem::path(projectPath) / "tools" / "visualizers" / "linux" / "FLTK");
    } else if (backend == VisualizerBackend::Auto) {
      return true;  // Auto always available (selects default)
    }
#elif defined(__APPLE__)
    if (backend == VisualizerBackend::WPF) {
      return false;  // WPF is Windows-only
    } else if (backend == VisualizerBackend::Qt) {
      return std::filesystem::exists(std::filesystem::path(projectPath) / "tools" / "visualizers" / "mac" / "Qt");
    } else if (backend == VisualizerBackend::FLTK) {
      return std::filesystem::exists(std::filesystem::path(projectPath) / "tools" / "visualizers" / "mac" / "FLTK");
    } else if (backend == VisualizerBackend::Auto) {
      return true;  // Auto always available (selects default)
    }
#endif
    
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
    // 1. Check environment variable for override
    if (const char* envBackend = GetEnv("MML_VISUALIZER_BACKEND")) {
      std::string backend(envBackend);
      // Convert to uppercase for case-insensitive comparison
      std::transform(backend.begin(), backend.end(), backend.begin(), ::toupper);
      
      if (backend == "WPF")  return VisualizerBackend::WPF;
      if (backend == "QT" || backend == "QT5" || backend == "QT6") 
        return VisualizerBackend::Qt;
      if (backend == "FLTK") return VisualizerBackend::FLTK;
      if (backend == "AUTO") return VisualizerBackend::Auto;
    }
    
    // 2. Get runtime setting
    VisualizerBackend backend = detail::GetUnifiedBackendState();
    
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

  // Get Windows WPF visualizer path (hardcoded structure)
  // tools/visualizers/win/WPF/<lowercase_dir>/<ExeName>.exe
  inline std::string GetWPFVisualizerPath(const std::string &name,
                                          const std::string &wpfDir) {
    std::string projectPath = GetProjectPath();
  #ifdef _WIN32
    return (std::filesystem::path(projectPath) / "tools" / "visualizers" / "win" / "WPF" / wpfDir / (name + ".exe")).string();
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
    std::string projectPath = GetProjectPath();
    
    // Resolve Auto to actual backend
    if (backend == VisualizerBackend::Auto) {
      backend = GetVisualizerBackend();
    }
    
#ifdef _WIN32
    if (backend == VisualizerBackend::WPF) {
      // WPF has subdirectory structure: tools/visualizers/win/WPF/<dir>/<name>.exe
      std::string wpfDir = ConvertToWPFDirName(visualizerName);
      return (std::filesystem::path(projectPath) / "tools" / "visualizers" / "win" / "WPF" / wpfDir / (visualizerName + ".exe")).string();
    } else if (backend == VisualizerBackend::Qt) {
      // Qt flat structure: all visualizers share Qt DLLs in one folder
      // tools/visualizers/win/Qt/<name>.exe
      return (std::filesystem::path(projectPath) / "tools" / "visualizers" / "win" / "Qt" / (visualizerName + ".exe")).string();
    } else if (backend == VisualizerBackend::FLTK) {
      // FLTK flat structure: tools/visualizers/win/FLTK/<name>_FLTK.exe
      return (std::filesystem::path(projectPath) / "tools" / "visualizers" / "win" / "FLTK" / (visualizerName + "_FLTK.exe")).string();
    }
#elif defined(__linux__)
    if (backend == VisualizerBackend::Qt) {
      // Qt structure: tools/visualizers/linux/Qt/<name>
      return (std::filesystem::path(projectPath) / "tools" / "visualizers" / "linux" / "Qt" / visualizerName).string();
    } else if (backend == VisualizerBackend::FLTK) {
      // FLTK structure: tools/visualizers/linux/FLTK/<name>_FLTK
      return (std::filesystem::path(projectPath) / "tools" / "visualizers" / "linux" / "FLTK" / (visualizerName + "_FLTK")).string();
    }
#elif defined(__APPLE__)
    if (backend == VisualizerBackend::Qt) {
      // Qt apps are .app bundles on macOS - executable is inside Contents/MacOS/
      // Structure: tools/visualizers/mac/Qt/<name>.app/Contents/MacOS/<name>
      return (std::filesystem::path(projectPath) / "tools" / "visualizers" / "mac" / "Qt" / (visualizerName + ".app") / "Contents" / "MacOS" / visualizerName).string();
    } else if (backend == VisualizerBackend::FLTK) {
      // FLTK structure: tools/visualizers/mac/FLTK/<name>_FLTK
      return (std::filesystem::path(projectPath) / "tools" / "visualizers" / "mac" / "FLTK" / (visualizerName + "_FLTK")).string();
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
  inline std::string GetRealFuncVisualizerPath() {
    return GetVisualizerPath("MML_RealFunctionVisualizer");
  }

  inline std::string GetSurfaceVisualizerPath() {
    return GetVisualizerPath("MML_ScalarFunction2D_Visualizer");
  }

  inline std::string GetParametricCurve2DVisualizerPath() {
    return GetVisualizerPath("MML_ParametricCurve2D_Visualizer");
  }

  inline std::string GetParametricCurve3DVisualizerPath() {
    return GetVisualizerPath("MML_ParametricCurve3D_Visualizer");
  }

  inline std::string GetVectorField2DVisualizerPath() {
    return GetVisualizerPath("MML_VectorField2D_Visualizer");
  }

  inline std::string GetVectorField3DVisualizerPath() {
    return GetVisualizerPath("MML_VectorField3D_Visualizer");
  }

  inline std::string GetParticle2DVisualizerPath() {
    return GetVisualizerPath("MML_ParticleVisualizer2D");
  }

  inline std::string GetParticle3DVisualizerPath() {
    return GetVisualizerPath("MML_ParticleVisualizer3D");
  }

  inline std::string GetParametricSurfaceVisualizerPath() {
    return GetVisualizerPath("MML_ParametricSurface_Visualizer");
  }

  inline std::string GetScalarFunction3DVisualizerPath() {
    return GetVisualizerPath("MML_ScalarFunction3D_Visualizer");
  }

  inline std::string GetRigidBodyVisualizerPath() {
    return GetVisualizerPath("MML_RigidBodyMovement_Visualizer");
  }

} // namespace MML

#endif // MML_VISUALIZATORS_H
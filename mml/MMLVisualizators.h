///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        MMLVisualizators.h                                                  ///
///  Description: Visualization utilities for functions, curves, fields               ///
///               Cross-platform support via FLTK/Qt backend                          ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
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

namespace MML
{

  // Matrix printing format options
  struct MatrixPrintFormat {
    int width = 10;               // Field width for each element
    int precision = 3;            // Number of decimal places
    bool scientific = false;      // Use scientific notation (e.g., 1.23e+05)
    bool fixed = true;            // Use fixed-point notation
    bool showHeader = true;       // Show "Rows: N Cols: M" header
    bool showBrackets = true;     // Show [ ] brackets around rows
    bool compactMode = false;     // Single line for small matrices
    std::string delimiter = ", "; // Delimiter between elements

    // Predefined formats
    static MatrixPrintFormat Default() {
      return MatrixPrintFormat{10, 3, false, true, true, true, false, ", "};
    }

    static MatrixPrintFormat Compact() {
      return MatrixPrintFormat{8, 2, false, true, false, true, true, ", "};
    }

    static MatrixPrintFormat Scientific() {
      return MatrixPrintFormat{12, 6, true, false, true, true, false, ", "};
    }

    static MatrixPrintFormat HighPrecision() {
      return MatrixPrintFormat{15, 10, false, true, true, true, false, ", "};
    }

    static MatrixPrintFormat NoDelimiter() {
      return MatrixPrintFormat{10, 3, false, true, true, true, false, " "};
    }
  };

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
    while (current.has_parent_path()) {
      if (std::filesystem::exists(current / "mml") &&
          std::filesystem::exists(current / "results") &&
          std::filesystem::exists(current / "tools")) {
        return current.string();
      }
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
      return std::filesystem::exists(projectPath + "/tools/visualizers/win/WPF");
    } else if (backend == VisualizerBackend::Qt) {
      return std::filesystem::exists(projectPath + "/tools/visualizers/win/Qt");
    } else if (backend == VisualizerBackend::FLTK) {
      return std::filesystem::exists(projectPath + "/tools/visualizers/win/FLTK");
    } else if (backend == VisualizerBackend::Auto) {
      return true;  // Auto always available (selects default)
    }
#elif defined(__linux__)
    if (backend == VisualizerBackend::WPF) {
      return false;  // WPF is Windows-only
    } else if (backend == VisualizerBackend::Qt) {
      return std::filesystem::exists(projectPath + "/tools/visualizers/linux/Qt");
    } else if (backend == VisualizerBackend::FLTK) {
      return std::filesystem::exists(projectPath + "/tools/visualizers/linux/FLTK");
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
    return projectPath + "/tools/visualizers/win/WPF/" + wpfDir + "/" + name +
          ".exe";
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
    return GetGlobalPath() + "/results/";
  }

  ///////////////////////////////////////////////////////////////////////////////
  //          GENERIC VISUALIZER PATH DISPATCHER
  ///////////////////////////////////////////////////////////////////////////////

  /// Helper: Convert visualizer name to WPF directory name
  /// WPF uses lowercase directory names with underscores
  /// Example: "MML_RealFunctionVisualizer" -> "real_function_visualizer"
  inline std::string ConvertToWPFDirName(const std::string& visualizerName) {
    // Special-case mappings to existing WPF directory names
    if (visualizerName == "MML_ParticleVisualizer2D") return "particle_2d_visualizer";
    if (visualizerName == "MML_ParticleVisualizer3D") return "particle_3d_visualizer";
    if (visualizerName == "MML_ParametricCurve2D_Visualizer") return "parametric_curve_2d_visualizer";
    if (visualizerName == "MML_ParametricCurve3D_Visualizer") return "parametric_curve_3d_visualizer";
    if (visualizerName == "MML_VectorField2D_Visualizer") return "vector_field_2d_visualizer";
    if (visualizerName == "MML_VectorField3D_Visualizer") return "vector_field_3d_visualizer";
    if (visualizerName == "MML_ScalarFunction2D_Visualizer") return "scalar_function_2d_visualizer";
    if (visualizerName == "MML_RealFunctionVisualizer") return "real_function_visualizer";

    std::string result = visualizerName;
    
    // Remove "MML_" prefix if present
    if (result.find("MML_") == 0) {
      result = result.substr(4);
    }
    
    // Convert to lowercase and replace capital letters with underscore + lowercase
    std::string converted;
    for (size_t i = 0; i < result.length(); ++i) {
      if (std::isupper(result[i])) {
        if (i > 0 && result[i-1] != '_') {
          converted += '_';
        }
        converted += std::tolower(result[i]);
      } else {
        converted += result[i];
      }
    }
    
    // Remove duplicate underscores
    std::string final;
    char prev = '\0';
    for (char c : converted) {
      if (c != '_' || prev != '_') {
        final += c;
      }
      prev = c;
    }
    
    return final;
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
      return projectPath + "/tools/visualizers/win/WPF/" + wpfDir + "/" + 
             visualizerName + ".exe";
    } else if (backend == VisualizerBackend::Qt) {
      // Qt structure: tools/visualizers/win/Qt/<name>/<name>.exe
      return projectPath + "/tools/visualizers/win/Qt/" + visualizerName + "/" +
             visualizerName + ".exe";
    } else if (backend == VisualizerBackend::FLTK) {
      // FLTK flat structure: tools/visualizers/win/FLTK/<name>_FLTK.exe
      return projectPath + "/tools/visualizers/win/FLTK/" + 
             visualizerName + "_FLTK.exe";
    }
#elif defined(__linux__)
    if (backend == VisualizerBackend::Qt) {
      // Qt structure: tools/visualizers/linux/Qt/<name>_Qt
      return projectPath + "/tools/visualizers/linux/Qt/" + 
             visualizerName + "_Qt";
    } else if (backend == VisualizerBackend::FLTK) {
      // FLTK structure: tools/visualizers/linux/FLTK/<name>_FLTK
      return projectPath + "/tools/visualizers/linux/FLTK/" + 
             visualizerName + "_FLTK";
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

} // namespace MML

#endif // MML_VISUALIZATORS_H
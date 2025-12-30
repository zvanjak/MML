# Qt Visualizers for Windows

This directory contains self-contained Qt visualizers with all necessary DLLs for Windows deployment.

## Contents

### MML_ParametricCurve3D_Visualizer
- **Executable**: `MML_ParametricCurve3D_Visualizer.exe` (111 KB)
- **Purpose**: 3D parametric curve visualization with OpenGL
- **Features**: Interactive 3D camera controls, multiple curves, orbit/pan/zoom
- **Data Format**: PARAMETRIC_CURVE_CARTESIAN_3D

### MML_RealFunctionVisualizer
- **Executable**: `MML_RealFunctionVisualizer.exe` (101 KB)
- **Purpose**: 2D real function visualization with OpenGL
- **Features**: Multi-function support, pan/zoom, grid, axes
- **Data Format**: REAL_FUNCTION

### MML_ParticleVisualizer3D
- **Executable**: `MML_ParticleVisualizer3D.exe` (118 KB)
- **Purpose**: 3D particle animation with OpenGL sphere rendering
- **Features**: Play/Pause/Restart, animation timeline, lighting
- **Data Format**: PARTICLE_SIMULATION_DATA_3D

### MML_VectorField3D_Visualizer
- **Executable**: `MML_VectorField3D_Visualizer.exe` (122 KB)
- **Purpose**: 3D vector field visualization with arrow rendering
- **Features**: 27,000+ vectors, color by magnitude, adjustable scale, interactive camera
- **Data Format**: VECTOR_FIELD_3D_CARTESIAN

### MML_ScalarFunction2D_Visualizer
- **Executable**: `MML_ScalarFunction2D_Visualizer.exe` (105 KB)
- **Purpose**: 2D scalar function f(x,y) as 3D surface visualization
- **Features**: Quad mesh surface, color height map, grid points, X/Y scaling
- **Data Format**: SCALAR_FUNCTION_CARTESIAN_2D

## Deployment

Each visualizer directory contains:
- **Main executable** (*.exe)
- **Qt DLLs** (Qt6Core, Qt6Gui, Qt6Widgets, Qt6OpenGL, Qt6OpenGLWidgets, Qt6Network, Qt6Svg)
- **System DLLs** (D3Dcompiler_47.dll, opengl32sw.dll, icuuc.dll)
- **Plugin directories**:
  - `platforms/` - Platform integration (qwindows.dll)
  - `styles/` - UI styles (qmodernwindowsstyle.dll)
  - `iconengines/` - Icon rendering (qsvgicon.dll)
  - `imageformats/` - Image loading (qgif, qico, qjpeg, qsvg)
  - `generic/`, `networkinformation/`, `tls/` - Additional plugins
  - `translations/` - Qt translations (31 languages)

**Total Size per Visualizer**: ~56 MB (with all Qt dependencies)

## Usage

### Command Line

```bash
# 3D Parametric Curves
.\MML_ParametricCurve3D_Visualizer\MML_ParametricCurve3D_Visualizer.exe <curve-file1> <curve-file2> ...

# 2D Real Functions
.\MML_RealFunctionVisualizer\MML_RealFunctionVisualizer.exe <function-file1> <function-file2> ...

# 3D Particle Animation
.\MML_ParticleVisualizer3D\MML_ParticleVisualizer3D.exe <simulation-file>

# 3D Vector Field
.\MML_VectorField3D_Visualizer\MML_VectorField3D_Visualizer.exe <vector-field-file>

# 2D Scalar Function (3D Surface)
.\MML_ScalarFunction2D_Visualizer\MML_ScalarFunction2D_Visualizer.exe <surface-file>
```

### Examples

```bash
# 3D: Lorenz attractor
.\MML_ParametricCurve3D_Visualizer\MML_ParametricCurve3D_Visualizer.exe ^
  ..\..\WPF\MML_ParametricCurve3D_Visualizer\data\example3_wpf_lorentz_system.txt

# 3D: 4-body trajectories
.\MML_ParametricCurve3D_Visualizer\MML_ParametricCurve3D_Visualizer.exe ^
  ..\..\WPF\MML_ParametricCurve3D_Visualizer\data\body0.txt ^
  ..\..\WPF\MML_ParametricCurve3D_Visualizer\data\body1.txt ^
  ..\..\WPF\MML_ParametricCurve3D_Visualizer\data\body2.txt ^
  ..\..\WPF\MML_ParametricCurve3D_Visualizer\data\body3.txt

# 2D: Real functions
.\MML_RealFunctionVisualizer\MML_RealFunctionVisualizer.exe ^
  ..\..\WPF\MML_RealFunctionVisualizer\data\real_func1.txt ^
  ..\..\WPF\MML_RealFunctionVisualizer\data\real_func2.txt

# 3D: Particle animation
.\MML_ParticleVisualizer3D\MML_ParticleVisualizer3D.exe ^
  ..\..\WPF\MML_ParticleVisualizer3D\data\SimData3D.txt

# 3D: Vector field (gravity)
.\MML_VectorField3D_Visualizer\MML_VectorField3D_Visualizer.exe ^
  ..\..\WPF\MML_VectorField3D_Visualizer\data\vector_field.txt

# 2D Scalar function (3D surface)
.\MML_ScalarFunction2D_Visualizer\MML_ScalarFunction2D_Visualizer.exe ^
  ..\..\WPF\MML_ScalarFunction2D_Visualizer\data\example3_wpf_surface1.txt
```

## Requirements

- **Windows 10/11** (64-bit)
- **OpenGL support** (typically available with graphics drivers)
- **No Qt installation required** - all DLLs included

## Controls

### MML_ParametricCurve3D_Visualizer (3D)
- **Left Mouse + Drag**: Rotate (orbit camera)
- **Right/Middle Mouse + Drag**: Pan view
- **Mouse Wheel**: Zoom in/out
- **Reset View Button**: Restore default camera

### MML_RealFunctionVisualizer (2D)
- **Left/Right Mouse + Drag**: Pan view
- **Mouse Wheel**: Zoom in/out
- **Load Function Button**: Add more functions
- **Reset View Button**: Auto-fit all functions

### MML_ParticleVisualizer3D (3D Animation)
- **Left Mouse + Drag**: Rotate (orbit camera)
- **Mouse Wheel**: Zoom in/out
- **Start Button**: Begin animation
- **Pause Button**: Pause/Resume animation
- **Restart Button**: Reset to first frame

### MML_VectorField3D_Visualizer (3D)
- **Left Mouse + Drag**: Rotate (orbit camera)
- **Mouse Wheel**: Zoom in/out
- **Vector Scale Slider**: Adjust arrow sizes (0.1x to 10.0x)
- **Color by Magnitude Checkbox**: Enable/disable magnitude-based coloring

### MML_ScalarFunction2D_Visualizer (3D Surface)
- **Left Mouse + Drag**: Rotate (orbit camera)
- **Mouse Wheel**: Zoom in/out
- **Scale X/Y Sliders**: Adjust domain scaling (0.1x to 20.0x)
- **Show Grid Points Checkbox**: Display/hide surface grid points
- **Color by Height Checkbox**: Enable/disable height-based coloring

## Known Issues

- HTML color format warnings in console (harmless, Qt QTextHtmlParser quirk)
- These are cosmetic and don't affect functionality

## Build Information

- **Qt Version**: 6.10.0
- **Compiler**: MSVC 19.44.35221.0 (Visual Studio 2022 Enterprise)
- **Build Type**: Release
- **C++ Standard**: C++17
- **Deployment Tool**: windeployqt

## Architecture

Both visualizers follow the same pattern:
- **MainWindow**: Qt Widgets UI with controls and info display
- **GLWidget**: QOpenGLWidget for OpenGL rendering
- **MMLData**: Data structures (curves, functions, colors)
- **MMLFileParser**: Data file parsing
- **CMake build system** with Qt6 + OpenGL

## Directory Structure

```
Qt/
├── MML_ParametricCurve3D_Visualizer/
│   ├── MML_ParametricCurve3D_Visualizer.exe
│   └── [Qt DLLs and plugins]
├── MML_RealFunctionVisualizer/
│   ├── MML_RealFunctionVisualizer.exe
│   └── [Qt DLLs and plugins]
├── MML_ParticleVisualizer3D/
│   ├── MML_ParticleVisualizer3D.exe
│   └── [Qt DLLs and plugins]
└── MML_VectorField3D_Visualizer/
    ├── MML_VectorField3D_Visualizer.exe
    └── [Qt DLLs and plugins]

Each visualizer includes:
- Qt6Core.dll, Qt6Gui.dll, Qt6Widgets.dll
- Qt6OpenGL.dll, Qt6OpenGLWidgets.dll
- D3Dcompiler_47.dll, opengl32sw.dll
- platforms/, styles/, iconengines/, etc.
```

## Distribution

These visualizers are fully self-contained and can be:
- Run directly from this directory
- Copied to other Windows machines
- Distributed as standalone applications
- No Qt installation needed on target machine

## Performance

- **Hardware Accelerated**: Uses OpenGL for smooth rendering
- **Responsive**: Interactive pan/zoom at 60+ FPS
- **Scalable**: Handles thousands of points efficiently
- **Multi-threaded**: Qt event loop + OpenGL rendering

## See Also

- Source code: `../../Qt/MML_ParametricCurve3D_Visualizer/`
- Source code: `../../Qt/MML_RealFunctionVisualizer/`
- Source code: `../../Qt/MML_ParticleVisualizer3D/`
- Source code: `../../Qt/MML_VectorField3D_Visualizer/`
- Source code: `../../Qt/MML_ScalarFunction2D_Visualizer/`
- FLTK versions: `../../visualizers/win/FLTK/` (coming soon)
- WPF versions: `../../WPF/`

## Version History

- **v1.0** (2025-12-12): Initial deployment
  - MML_ParametricCurve3D_Visualizer: 3D curve visualization
  - MML_RealFunctionVisualizer: 2D function visualization
  - MML_ParticleVisualizer3D: 3D particle animation with spheres
  - MML_VectorField3D_Visualizer: 3D vector field with 27,000+ arrows
  - MML_ScalarFunction2D_Visualizer: 2D scalar function as 3D surface
  - Qt 6.10.0 with OpenGL acceleration
  - Full Qt DLL deployment via windeployqt

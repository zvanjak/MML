# MML Visualizers - Windows Distribution

Pre-built Windows executables for all MML visualizers.

## FLTK Visualizers (4 apps)

Self-contained single executables for 2D visualization:

| Visualizer | Description |
|------------|-------------|
| `MML_RealFunctionVisualizer_FLTK.exe` | Plot real functions y=f(x) |
| `MML_ParametricCurve2D_Visualizer_FLTK.exe` | 2D parametric curves (x(t), y(t)) |
| `MML_VectorField2D_Visualizer_FLTK.exe` | 2D vector field visualization |
| `MML_ParticleVisualizer2D_FLTK.exe` | 2D particle simulation animation |

**Usage:** Simply run the `.exe` file - no installation required.

## Qt Visualizers (8 apps)

Modern themed visualizers with dark/light mode support:

| Folder | Description |
|--------|-------------|
| `MML_RealFunctionVisualizer/` | Plot real functions y=f(x) |
| `MML_ParametricCurve2D_Visualizer/` | 2D parametric curves |
| `MML_ParametricCurve3D_Visualizer/` | 3D parametric curves with OpenGL |
| `MML_VectorField2D_Visualizer/` | 2D vector field visualization |
| `MML_VectorField3D_Visualizer/` | 3D vector field with OpenGL |
| `MML_ParticleVisualizer2D/` | 2D particle simulation |
| `MML_ParticleVisualizer3D/` | 3D particle simulation with OpenGL |
| `MML_ScalarFunction2D_Visualizer/` | 2D scalar function heatmaps |

**Usage:** Each folder contains the `.exe` plus required Qt DLLs. Run the `.exe` directly from its folder.

## Features

- **Dark/Light Theme Toggle**: All Qt visualizers support theme switching via title bar icon
- **Interactive Controls**: Zoom, pan, and animation controls
- **File Loading**: Load data files from WPF/data folders
- **OpenGL Acceleration**: 3D visualizers use hardware-accelerated rendering

## Data Files

Sample data files are located in:
- `WPF/MML_RealFunctionVisualizer/data/`
- `WPF/MML_ParametricCurve2D_Visualizer/data/`
- (and similarly for other visualizers)

---

*Built with FLTK 1.3.x and Qt 6.10.0*

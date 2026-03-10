# MML Linux Visualizers

Complete deployment of MML mathematical visualizers for Linux.

## FLTK Visualizers (2D, Cross-platform)

Located in `FLTK/`

**Features:**
- Modern themed UI (Dark/Light theme toggle)
- Tokyo Night dark theme
- Teal/Emerald light theme
- Click 🌓 icon to toggle themes
- Standalone executables (no external dependencies except FLTK)

**Visualizers:**
- `MML_RealFunctionVisualizer_FLTK` - Plot 1D real functions
- `MML_ParametricCurve2D_Visualizer_FLTK` - 2D parametric curves
- `MML_ParticleVisualizer2D_FLTK` - 2D particle simulations
- `MML_VectorField2D_Visualizer_FLTK` - 2D vector fields

**Usage:**
```bash
cd FLTK
./MML_RealFunctionVisualizer_FLTK
```

## Qt Visualizers (3D + 2D, OpenGL)

Located in `Qt/`

**Features:**
- Full 3D OpenGL rendering
- Interactive camera controls
- Self-contained deployment (bundled Qt 6.9.3 libraries)
- No Qt installation required
- ICU libraries included

**Visualizers (3D):**
- `MML_ParametricCurve3D_Visualizer` - 3D parametric curves
- `MML_ParticleVisualizer3D` - 3D particle simulations
- `MML_VectorField3D_Visualizer` - 3D vector fields

**Visualizers (2D):**
- `MML_RealFunctionVisualizer` - Plot 1D real functions
- `MML_ParametricCurve2D_Visualizer` - 2D parametric curves  
- `MML_ParticleVisualizer2D` - 2D particle simulations
- `MML_ScalarFunction2D_Visualizer` - 2D scalar function heatmaps
- `MML_VectorField2D_Visualizer` - 2D vector fields

**Dependencies (one-time install):**
```bash
sudo apt install libxkbcommon-x11-0 libxcb-xkb1
```

**Usage:**
```bash
cd Qt
./run_MML_ParametricCurve3D_Visualizer.sh <datafile>
./run_MML_VectorField3D_Visualizer.sh <datafile>
```

## Example Data Files

Data files are located in the WPF project directories:
- `../../WPF/MML_ParametricCurve3D_Visualizer/data/`
- `../../WPF/MML_VectorField3D_Visualizer/data/`
- `../../WPF/MML_ScalarFunction2D_Visualizer/data/`

## Build Information

- FLTK visualizers: Built Jan 27, 2026 with theming support
- Qt visualizers: Built Jan 27, 2026 with Qt 6.9.3
- Total size: ~70MB (including all Qt libraries)


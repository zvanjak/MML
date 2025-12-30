# WPF Visualizers - Windows Deployment

This directory contains self-contained deployments of all WPF-based MML Visualizers for Windows.

## Available Visualizers

| Visualizer | Description | Executable |
|------------|-------------|------------|
| **real_function_visualizer** | 2D real function plotting | MML_RealFunctionVisualizer.exe |
| **parametric_curve_2d_visualizer** | 2D parametric curve visualization | MML_ParametricCurve2D_Visualizer.exe |
| **parametric_curve_3d_visualizer** | 3D parametric curve visualization | MML_ParametricCurve3D_Visualizer.exe |
| **particle_2d_visualizer** | 2D particle simulation animation | MML_ParticleVisualizer2D.exe |
| **particle_3d_visualizer** | 3D particle simulation animation | MML_ParticleVisualizer3D.exe |
| **vector_field_2d_visualizer** | 2D vector field visualization | MML_VectorField2D_Visualizer.exe |
| **vector_field_3d_visualizer** | 3D vector field visualization | MML_VectorField3D_Visualizer.exe |
| **scalar_function_2d_visualizer** | 2D scalar function heatmap | MML_ScalarFunction2D_Visualizer.exe |

## Requirements

- **Operating System**: Windows 10/11
- **.NET Runtime**: .NET 8.0 Desktop Runtime (Windows x64)
  - Download from: https://dotnet.microsoft.com/download/dotnet/8.0

## Installation

No installation required! Each visualizer directory contains:
- Main executable (`.exe`)
- Required .NET assemblies (`.dll`)
- Runtime configuration (`.json`)
- Debug symbols (`.pdb`) - optional, can be deleted to save space

## Usage

### Running Visualizers

Each visualizer accepts data files as command-line arguments:

```powershell
# Real function visualizer
.\real_function_visualizer\MML_RealFunctionVisualizer.exe path\to\data.txt

# Parametric curve 3D
.\parametric_curve_3d_visualizer\MML_ParametricCurve3D_Visualizer.exe curve1.txt curve2.txt

# Particle simulation
.\particle_3d_visualizer\MML_ParticleVisualizer3D.exe simulation.txt
```

### Data File Formats

Each visualizer expects specific data formats:

- **REAL_FUNCTION**: Real function data (x1, x2, NumPoints, function values)
- **PARAMETRIC_CURVE_CARTESIAN_2D/3D**: Parametric curve points
- **PARTICLE_SIMULATION_DATA_2D/3D**: Particle positions over time
- **VECTOR_FIELD_2D/3D_CARTESIAN**: Vector field data at grid points
- **SCALAR_FUNCTION_2D**: Scalar values on 2D grid

See the WPF project data folders for example files.

## Deployment Structure

```
visualizers/win/WPF/
├── real_function_visualizer/
│   ├── MML_RealFunctionVisualizer.exe
│   ├── MML_RealFunctionVisualizer.dll
│   ├── MML_VisualizersBase.dll
│   ├── MathLibCSharp.dll
│   ├── WPF3DHelperLib.dll (for 3D visualizers)
│   └── *.json, *.pdb
├── parametric_curve_2d_visualizer/
│   └── [same structure]
├── parametric_curve_3d_visualizer/
│   └── [same structure]
... [and so on]
```

## Updating Deployments

To update deployed visualizers after rebuilding in Visual Studio:

```powershell
# From project root
.\WPF\update_apps.ps1

# Force update all files
.\WPF\update_apps.ps1 -Force

# Show help
.\WPF\update_apps.ps1 -Help
```

The update script:
- ✅ Copies only files that are newer (efficient)
- ✅ Creates deployment folders automatically
- ✅ Validates that builds exist before copying
- ✅ Reports detailed statistics

## File Sizes

Each WPF visualizer deployment is approximately **200-260 KB**:
- Executable: ~150 KB
- .NET assemblies: ~50-110 KB
- Configuration/symbols: ~10 KB

Total for all 8 visualizers: **~1.7 MB**

## Technology

- **Framework**: .NET 8.0 (Windows Desktop)
- **UI**: WPF (Windows Presentation Foundation)
- **3D Rendering**: WPF 3D with XAML
- **Language**: C# 12

## Dependencies

All required assemblies are included in each visualizer's folder:
- `MML_VisualizersBase.dll` - Common visualizer base classes
- `MathLibCSharp.dll` - Mathematical operations
- `WPF3DHelperLib.dll` - 3D rendering helpers (3D visualizers only)

No additional dependencies needed!

## Troubleshooting

### "Application cannot start"
- Install .NET 8.0 Desktop Runtime from Microsoft
- Ensure you're running on Windows 10 or later

### "File not found" errors
- Check that all DLL files are in the same folder as the EXE
- Verify data file paths are correct

### Performance Issues
- WPF 3D visualizers use software rendering by default
- Enable GPU acceleration in Windows graphics settings for better performance

## Development

### Building from Source
1. Open `WPF\MML_Visualizers.sln` in Visual Studio 2022
2. Select **Release** configuration
3. Build solution (Ctrl+Shift+B)
4. Run `.\WPF\update_apps.ps1` to deploy

### Source Locations
- **Source Code**: `WPF\[VisualizerName]\`
- **Build Output**: `WPF\[VisualizerName]\bin\Release\net8.0-windows\`
- **Deployment**: `visualizers\win\WPF\[visualizer_name]\`

## Version History

| Version | Date | Changes |
|---------|------|---------|
| 1.0 | 2025-12-12 | Initial deployment with automated update script |

## License

Part of the MML_Visualizers project.

## Support

For issues or questions:
- Check example data files in `WPF\[VisualizerName]\data\`
- Review project README at repository root
- Examine source code in `WPF\` directory

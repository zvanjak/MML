# 📊 MML Visualization Suite

MML includes a powerful **cross-platform visualization suite** for functions, fields, curves, surfaces, and particle systems. All visualizers work on Windows (WPF), Linux (Qt), and macOS (Qt).

> Back to the [main README](../README.md) • See also [Code examples](README_Code_examples.md) and [Usage examples](README_Usage_examples.md).

---

## Available Visualizers

All visualizers have **complete cross-platform support** with multiple backend options:
- **Windows**: WPF-based visualizers (primary), Qt and FLTK also available
- **Linux**: Qt-based visualizers (primary), FLTK for lightweight 2D
- **macOS**: Qt-based visualizers (primary), FLTK for lightweight 2D

| Visualizer | Purpose | Output | Platform |
|------------|---------|--------|----------|
| **RealFunctionVisualizer** | Plot 1D functions | 2D graphs | Windows (WPF), Linux/macOS (Qt/FLTK) |
| **MultiRealFunctionVisualizer** | Compare multiple 1D functions | 2D overlaid graphs | Windows (WPF), Linux/macOS (Qt/FLTK) |
| **ScalarFunction2DVisualizer** | 3D surface plots | 3D surfaces | Windows (WPF), Linux/macOS (Qt) |
| **ParametricCurve2DVisualizer** | 2D parametric curves | 2D curves | Windows (WPF), Linux/macOS (Qt/FLTK) |
| **ParametricCurve3DVisualizer** | 3D parametric curves | 3D curves | Windows (WPF), Linux/macOS (Qt) |
| **VectorField2DVisualizer** | 2D vector field arrows | 2D arrow plots | Windows (WPF), Linux/macOS (Qt/FLTK) |
| **VectorField3DVisualizer** | 3D vector field arrows | 3D arrow plots | Windows (WPF), Linux/macOS (Qt) |
| **ParticleVisualizer2D** | Animated 2D particle systems | 2D animations with playback | Windows (WPF), Linux/macOS (Qt) |
| **ParticleVisualizer3D** | Animated 3D particle systems | 3D animations with playback | Windows (WPF), Linux/macOS (Qt) |

### Where to find examples
- `../src/visualization_examples/` — Ready-to-run demos for every visualizer type
- `../src/book/chapters/Chapter_02_vizualization/` — Comprehensive examples from the companion book
- Examples of custom styling, multi-plot layouts, and animation controls

---

## Gallery

### Windows — WPF Visualizers

<table>
<tr>
<td align="center" width="25%">

**Real Functions**

![WPF Real](images/readme/visualization_suite/win/wpf_real_func_multi_damped_oscillations.png)

</td>
<td align="center" width="25%">

**Real Functions (Lorentz)**

![WPF Lorentz](images/readme/visualization_suite/win/wpf_real_func_multi_Lorentz.png)

</td>
<td align="center" width="25%">

**Scalar Function 2D**

![WPF Scalar 2D](images/readme/visualization_suite/win/wpf_scalar_func_2d.png)

</td>
<td align="center" width="25%">

**Scalar Function 3D**

![WPF Scalar 3D](images/readme/visualization_suite/win/wpf_scalar_func_3d.png)

</td>
</tr>
<tr>
<td align="center">

**Parametric Curve 2D**

![WPF Curve 2D](images/readme/visualization_suite/win/wpf_param_curve_2d_butterfly.png)

</td>
<td align="center">

**Parametric Curve 3D**

![WPF Curve 3D](images/readme/visualization_suite/win/wpf_param_curve_3d.png)

</td>
<td align="center">

**Parametric Surface**

![WPF Surface](images/readme/visualization_suite/win/wpf_param_surface.png)

</td>
<td align="center">

**Vector Field 3D**

![WPF VecField](images/readme/visualization_suite/win/wpf_vector_field_3d.png)

</td>
</tr>
<tr>
<td align="center">

**Particle Visualizer 2D**

![WPF Particle 2D](images/readme/visualization_suite/win/wpf_particle_vis_2d.png)

</td>
<td align="center">

**Particle Visualizer 3D**

![WPF Particle 3D](images/readme/visualization_suite/win/wpf_particle_vis_3d.png)

</td>
<td align="center">

**Rigid Body Simulation**

![WPF Rigid](images/readme/visualization_suite/win/wpf_rigid_body.png)

</td>
<td align="center">

</td>
</tr>
</table>

### Windows — Qt Visualizers

<table>
<tr>
<td align="center" width="25%">

**Real Functions**

![Qt Real](images/readme/visualization_suite/win/win_qt_real_func_multi.png)

</td>
<td align="center" width="25%">

**Scalar Function 2D**

![Qt Scalar 2D](images/readme/visualization_suite/win/win_qt_scalar_func_2d.png)

</td>
<td align="center" width="25%">

**Scalar Function 3D**

![Qt Scalar 3D](images/readme/visualization_suite/win/win_qt_scalar_func_3d.png)

</td>
<td align="center" width="25%">

**Parametric Surface**

![Qt Surface](images/readme/visualization_suite/win/win_qt_param_surface.png)

</td>
</tr>
<tr>
<td align="center">

**Parametric Curve 2D**

![Qt Curve 2D](images/readme/visualization_suite/win/win_qt_param_curve_2d.png)

</td>
<td align="center">

**Parametric Curve 3D**

![Qt Curve 3D](images/readme/visualization_suite/win/win_qt_param_curve_3d.png)

</td>
<td align="center">

**Vector Field 3D**

![Qt VecField](images/readme/visualization_suite/win/win_qt_vector_field_3d.png)

</td>
<td align="center">

**Rigid Body Simulation**

![Qt Rigid](images/readme/visualization_suite/win/win_qt_rigid_body.png)

</td>
</tr>
<tr>
<td align="center">

**Particle Visualizer 2D**

![Qt Particle 2D](images/readme/visualization_suite/win/win_qt_particle_vis_2d.png)

</td>
<td align="center">

**Particle Visualizer 3D**

![Qt Particle 3D](images/readme/visualization_suite/win/win_qt_particle_vis_3d.png)

</td>
<td align="center">

</td>
<td align="center">

</td>
</tr>
</table>

### Linux — Qt Visualizers

<table>
<tr>
<td align="center" width="25%">

**Real Functions**

![Linux Real](images/readme/visualization_suite/linux/linux_qt_real_func_multri.png)

</td>
<td align="center" width="25%">

**Scalar Function 2D**

![Linux Scalar 2D](images/readme/visualization_suite/linux/linux_qt_scalar_func_2d.png)

</td>
<td align="center" width="25%">

**Scalar Function 3D**

![Linux Scalar 3D](images/readme/visualization_suite/linux/linux_qt_scalar_func_3d.png)

</td>
<td align="center" width="25%">

**Parametric Surface**

![Linux Surface](images/readme/visualization_suite/linux/linux_qt_param_surface.png)

</td>
</tr>
<tr>
<td align="center">

**Parametric Curve 2D**

![Linux Curve 2D](images/readme/visualization_suite/linux/linux_qt_param_curve_2d.png)

</td>
<td align="center">

**Parametric Curve 3D**

![Linux Curve 3D](images/readme/visualization_suite/linux/linux_qt_param_curve_3d.png)

</td>
<td align="center">

**Vector Field 3D**

![Linux VecField](images/readme/visualization_suite/linux/linux_qt_vector_field_3d.png)

</td>
<td align="center">

**Rigid Body Simulation**

![Linux Rigid](images/readme/visualization_suite/linux/linux_qt_rigid_body_vis.png)

</td>
</tr>
<tr>
<td align="center">

**Particle Visualizer 2D**

![Linux Particle 2D](images/readme/visualization_suite/linux/linux_qt_particle_vis_2d.png)

</td>
<td align="center">

**Particle Visualizer 3D**

![Linux Particle 3D](images/readme/visualization_suite/linux/linux_qt_particle_vis_3d.png)

</td>
<td align="center">

**Scalar 2D (Dark Theme)**

![Linux Scalar Dark](images/readme/visualization_suite/linux/linux_qt_scalar_func_2d_dark.png)

</td>
<td align="center">

**Parametric Surface (Dark)**

![Linux Surface Dark](images/readme/visualization_suite/linux/linux_qt_param_surface_dark.png)

</td>
</tr>
</table>

### macOS — Qt Visualizers

<table>
<tr>
<td align="center" width="25%">

**Real Functions**

![Mac Real](images/readme/visualization_suite/mac/mac_qt_real_func_multi.png)

</td>
<td align="center" width="25%">

**Real Functions (Lorentz)**

![Mac Lorentz](images/readme/visualization_suite/mac/mac_qt_real_func_multi_Lorentz_system.png)

</td>
<td align="center" width="25%">

**Scalar Function 2D**

![Mac Scalar 2D](images/readme/visualization_suite/mac/mac_qt_scalar_func_2d_monkey_saddle.png)

</td>
<td align="center" width="25%">

**Scalar Function 3D**

![Mac Scalar 3D](images/readme/visualization_suite/mac/mac_qt_scalar_function_3d_gyroid.png)

</td>
</tr>
<tr>
<td align="center">

**Parametric Curve 2D**

![Mac Curve 2D](images/readme/visualization_suite/mac/mac_qt_param_curve_2de_butterfly.png)

</td>
<td align="center">

**Parametric Curve 3D**

![Mac Curve 3D](images/readme/visualization_suite/mac/mac_qt_param_curve_3d_trefoil.png)

</td>
<td align="center">

**Parametric Surface**

![Mac Surface](images/readme/visualization_suite/mac/mac_qt_param_surface_klein.png)

</td>
<td align="center">

**Vector Field 3D**

![Mac VecField](images/readme/visualization_suite/mac/mac_qt_vector_field_3d_gravity.png)

</td>
</tr>
<tr>
<td align="center">

**Particle Visualizer 2D**

![Mac Particle 2D](images/readme/visualization_suite/mac/mac_qt_particle_visualizer_2d.png)

</td>
<td align="center">

**Particle Visualizer 3D**

![Mac Particle 3D](images/readme/visualization_suite/mac/mac_qt_partice_visualizer_3d.png)

</td>
<td align="center">

**Scalar Function 2D (Ripple)**

![Mac Ripple](images/readme/visualization_suite/mac/mac_qt_scalar_func_2d_ripple.png)

</td>
<td align="center">

**Rigid Body Simulation**

![Mac Rigid](images/readme/visualization_suite/mac/mac_qt_rigid_body.png)

</td>
</tr>
</table>

---

## Visualization Code Examples

**Real Function Plotting — Lorenz System Time Series:**

```cpp
// Lorenz attractor — chaotic time evolution of x(t), y(t), z(t)
ODESystem lorenz_system(3, [](Real t, const Vector<Real>& x, Vector<Real>& dxdt) {
    const Real sigma = 10.0, rho = 28.0, beta = 8.0 / 3.0;
    dxdt[0] = sigma * (x[1] - x[0]);
    dxdt[1] = x[0] * (rho - x[2]) - x[1];
    dxdt[2] = x[0] * x[1] - beta * x[2];
});

Vector<Real> initial_state({ 1.0, 1.0, 1.0 });
CashKarpIntegrator solver(lorenz_system);
ODESystemSolution sol = solver.integrate(initial_state, 0.0, 50.0, 0.001, 1e-10, 0.001);

// Extract solution components as real functions via spline interpolation
Vector<Real> t_vals = sol.getTValues();
SplineInterpRealFunc x_t(t_vals, sol.getXValues(0));
SplineInterpRealFunc y_t(t_vals, sol.getXValues(1));
SplineInterpRealFunc z_t(t_vals, sol.getXValues(2));

std::vector<IRealFunction*> time_series = { &x_t, &y_t, &z_t };
Visualizer::VisualizeMultiRealFunction(
    time_series, "Lorenz System: Chaotic Time Evolution",
    { "x(t)", "y(t)", "z(t)" }, 0.0, 50.0, 1000, "lorenz_time_series.mml");
```

| Windows (WPF) | Linux (Qt) | macOS (Qt) |
|:-------------:|:----------:|:----------:|
| ![Lorenz Win](images/readme/visualization_suite/win/wpf_real_func_multi_Lorentz.png) | ![Lorenz Linux](images/readme/visualization_suite/linux/linux_qt_real_func_multri.png) | ![Lorenz Mac](images/readme/visualization_suite/mac/mac_qt_real_func_multi_Lorentz_system.png) |

**Scalar Function 2D — Ripple (2D Sinc / Sombrero):**

```cpp
// 2D Sinc function — concentric ripples with central peak
ScalarFunction<2> sinc_2d{ [](const VectorN<Real, 2>& v) {
    Real r = std::sqrt(v[0]*v[0] + v[1]*v[1]);
    if (r < 1e-10) return Real(50.0);
    return 50.0 * std::sin(r) / r;
}};

Visualizer::VisualizeScalarFunc2DCartesian(
    sinc_2d, "2D Sinc (Sombrero): z = 50·sin(r)/r",
    -15.0, 15.0, 80, -15.0, 15.0, 80, "sinc_2d.mml");
```

| Windows (WPF) | Linux (Qt) | macOS (Qt) |
|:-------------:|:----------:|:----------:|
| ![Ripple Win](images/readme/visualization_suite/win/wpf_scalar_func_2d.png) | ![Ripple Linux](images/readme/visualization_suite/linux/linux_qt_scalar_func_2d.png) | ![Ripple Mac](images/readme/visualization_suite/mac/mac_qt_scalar_func_2d_ripple.png) |

**Parametric Curves — Trefoil Knot:**

```cpp
// Trefoil knot: x = sin(t) + 2sin(2t), y = cos(t) - 2cos(2t), z = -sin(3t)
auto trefoil = [](Real t) {
    Real x = 50.0 * (std::sin(t) + 2.0*std::sin(2.0*t));
    Real y = 50.0 * (std::cos(t) - 2.0*std::cos(2.0*t));
    Real z = 50.0 * (-std::sin(3.0*t));
    return VectorN<Real, 3>{x, y, z};
};

ParametricCurve<3> knot(trefoil);
Visualizer::VisualizeParamCurve3D(
    knot, "Trefoil Knot", 0.0, 2.0*Constants::PI, 500, "trefoil.mml");
```

| Windows (WPF) | Linux (Qt) | macOS (Qt) |
|:-------------:|:----------:|:----------:|
| ![Trefoil Win](images/readme/visualization_suite/win/wpf_param_curve_3d.png) | ![Trefoil Linux](images/readme/visualization_suite/linux/linux_qt_param_curve_3d.png) | ![Trefoil Mac](images/readme/visualization_suite/mac/mac_qt_param_curve_3d_trefoil.png) |

**Vector Fields — Two-Body Gravity:**

```cpp
// Gravity field of two masses
VectorFunction<3> gravity{[](const VectorN<Real, 3> &x) {
    const VectorN<Real, 3> x1{100, 0, 0}, x2{-100, 0, 0};
    const Real m1 = 1000, m2 = 1000, G = 10;
    return -G * m1 * (x - x1) / pow((x - x1).NormL2(), 3)
           -G * m2 * (x - x2) / pow((x - x2).NormL2(), 3);
}};

Visualizer::VisualizeVectorField3DCartesian(
    gravity, "Gravity Field",
    -200, 200, 15, -200, 200, 15, -200, 200, 15, "gravity.mml");
```

| Windows (WPF) | Linux (Qt) | macOS (Qt) |
|:-------------:|:----------:|:----------:|
| ![Field Win](images/readme/visualization_suite/win/wpf_vector_field_3d.png) | ![Field Linux](images/readme/visualization_suite/linux/linux_qt_vector_field_3d.png) | ![Field Mac](images/readme/visualization_suite/mac/mac_qt_vector_field_3d_gravity.png) |

---

## Data Export

All visualizers use serialized data files that can also be loaded by external tools (Python, MATLAB, etc.):

```cpp
// Serialize function data
Serializer::SaveRealFunc(f, "Function", 0, 10, 100, "data.txt");

// Serialize ODE solution
Serializer::SaveODESolutionAsMultiFunc(
    solution, "ODE Solution", {"x", "v"}, "ode_data.txt");

// Serialize 2D vector field
Serializer::SaveVectorFunc2DCartesian(
    field, "Vector Field", -10, 10, 20, -10, 10, 20, "field_data.txt");
```

See [Visualizers](tools/Visualizers.md) and [Serializer](tools/Serializer.md) for the complete API.
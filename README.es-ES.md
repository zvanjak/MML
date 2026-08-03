<div align="center">

# 🔢 MML — Biblioteca Matemática Minimalista

### **El Kit de Herramientas Completo de Cómputo Numérico en C++**

*Un solo encabezado • Multiplataforma • Visualización incluida*

[![Ubuntu](https://github.com/zvanjak/MML/workflows/Ubuntu/badge.svg)](https://github.com/zvanjak/MML/actions?query=workflow%3AUbuntu)
[![Windows](https://github.com/zvanjak/MML/workflows/Windows/badge.svg)](https://github.com/zvanjak/MML/actions?query=workflow%3AWindows)
[![macOS](https://github.com/zvanjak/MML/workflows/macOS/badge.svg)](https://github.com/zvanjak/MML/actions?query=workflow%3AmacOS)
[![C++17](https://img.shields.io/badge/C%2B%2B-17-blue.svg)](https://isocpp.org/std/the-standard)
[![Tests](https://img.shields.io/badge/tests-4540%20passing-brightgreen.svg)](tests/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE.md)

---

##  Misión 

**🚀 Simplemente `#include "MML.h"` y calcula** — vectores, matrices, tensores, solvers de EDO, valores propios, integración y más!

[Inicio Rápido](#-inicio-rápido) • [Características](#-características-clave) • [Documentación](#-documentación) • [Ejemplos](#-ejemplos-reales) • [Visualización](#-suite-de-visualización)

</div>

---

## 🎯 ¿Qué es MML?

MML es una **biblioteca matemática C++ completa de un solo encabezado** para cómputo numérico. Con solo una directiva `#include "MML.h"`, obtienes acceso a un kit de herramientas completo de objetos matemáticos y algoritmos — desde vectores y matrices básicos hasta solvers de EDO y operaciones de campos.
Y como bonificación, ¡también obtienes herramientas de visualización multiplataforma para trazar funciones, curvas, superficies y campos vectoriales!

<table>
<tr>
<td width="50%">

### El Problema

La mayoría de las bibliotecas matemáticas C++ requieren:
- Sistemas de compilación y dependencias complejos
- Vinculación contra múltiples bibliotecas
- Configuración específica de la plataforma
- Curvas de aprendizaje pronunciadas

**Resultado:** Posiblemente horas de configuración antes de escribir código real.

</td>
<td width="50%">

### La Solución MML

**Integración sin fricción:**

```cpp
#include "MML.h"
// Eso es todo. Comienza a calcular.
```

- ✅ Un solo archivo de encabezado
- ✅ Puro C++17, sin dependencias
- ✅ Funciona en Windows, Linux, Mac
- ✅ 4,540 pruebas garantizan la corrección

</td>
</tr>
</table>

Si lo necesitas, también puedes usarlo por partes incluyendo solo los encabezados seleccionados del directorio `mml/`.

---

## 🏛️ Filosofía de Diseño

MML se basa en tres principios fundamentales:

<table>
<tr>
<td align="center" width="33%">

### 🎯 Completitud y Simplicidad

Un **kit de herramientas completo** que cubre vectores, matrices, tensores, EDOs, integración, operaciones de campos y más. Una **sintaxis intuitiva** hace que los objetos matemáticos sean ciudadanos de primera clase en C++.

</td>
<td align="center" width="33%">

### 🔬 Correctitud y Precisión

**4,540 pruebas unitarias** validan cada algoritmo contra soluciones analíticas conocidas. Los métodos numéricos **reportan la precisión lograda**, y bancos de pruebas dedicados demuestran la exactitud.

</td>
<td align="center" width="33%">

### ⚡ Integración y Visualización Triviales

**Un solo encabezado, cero dependencias** — simplemente `#include "MML.h"` y comienza a calcular inmediatamente. **Visualizadores multiplataforma** para funciones, superficies, campos vectoriales y sistemas de partículas en Windows, Linux y Mac.

</td>
</tr>
</table>

### Ejemplo Destacado: Verificar el Teorema de la Divergencia de Gauss

📄 *[Ver código fuente completo](src/readme_examples/readme00_fundamental_theorems.cpp)*

```cpp
// Verificar ∫∫∫(∇·F)dV = ∮∮(F·n̂)dS sobre un cubo unitario
// Definir campo vectorial F(x,y,z) = (x², y², z²)
VectorFunction<3> F([](const VectorN<Real, 3>& p) {
    return VectorN<Real, 3>{ p[0]*p[0], p[1]*p[1], p[2]*p[2] };
});

// Calcular la divergencia NUMÉRICAMENTE - ¡MML calcula ∇·F automáticamente!
ScalarFunctionFromStdFunc<3> divF([&F](const VectorN<Real, 3>& p) {
    return VectorFieldOperations::DivCart<3>(F, p);
});

// Integral de volumen: ∫∫∫(∇·F)dV
auto y_lo = [](Real) { return 0.0; };  auto y_hi = [](Real) { return 1.0; };
auto z_lo = [](Real,Real) { return 0.0; };  auto z_hi = [](Real,Real) { return 1.0; };
Real volIntegral = Integrate3D(divF, GAUSS10, 0, 1, y_lo, y_hi, z_lo, z_hi).value;

// Integral de superficie: ∮∮(F·n̂)dS a través de las 6 caras
Cube3D unitCube(1.0, Point3Cartesian(0.5, 0.5, 0.5));
Real surfIntegral = SurfaceIntegration::SurfaceIntegral(F, unitCube, 1e-8);

std::cout << "Integral de volumen:  " << volIntegral << "\n";   // 3.0000000000
std::cout << "Integral de superficie: " << surfIntegral << "\n";  // 3.0000000000
std::cout << "Error: " << std::abs(volIntegral - surfIntegral) << "\n";  // 9.77e-15 ✓
```

---

## 🚀 Inicio Rápido

### Instalación

**Opción 1: Simplemente obtén el encabezado único (Recomendado)**
```bash
# Descarga MML.h directamente desde el repositorio
curl -O https://raw.githubusercontent.com/zvanjak/MML/master/mml/single_header/MML.h

# Inclúyelo en tu proyecto
#include "MML.h"
```

**Opción 2: Obtén todo el repositorio (si necesitas módulos individuales)**
```bash
git clone https://github.com/zvanjak/MML.git
cd MML
cmake -B build
cmake --build build
```


**Opción 3: Visual Studio Code**

1. Abre VS Code y presiona `Ctrl+Shift+P` (o `Cmd+Shift+P` en Mac)
2. Escribe **"Git: Clone"** y presiona Enter
3. Pega la URL del repositorio:
   ```
   https://github.com/zvanjak/MML.git
   ```
4. Selecciona una carpeta y abre el repositorio clonado
5. Cuando se te solicite, instala las extensiones recomendadas (C/C++, CMake Tools)
6. Presiona `Ctrl+Shift+P` → **"CMake: Configure"** para configurar la compilación
7. Presiona `F7` para compilar o usa la barra de estado de CMake

> 💡 **Consejo:** La [extensión CMake Tools](https://marketplace.visualstudio.com/items?itemName=ms-vscode.cmake-tools) proporciona IntelliSense, botones de compilación e integración de pruebas de forma nativa.

### Primer Programa

```cpp
#include "MML.h"
using namespace MML;

int main() {
    // Crear una matriz y resolver un sistema lineal
    Matrix<Real> A{3, 3, {4, 1, 2, 
                          1, 5, 1, 
                          2, 1, 6}};
    Vector<Real> b{1, 2, 3};
    
    LUSolver<Real> solver(A);
    Vector<Real> x = solver.Solve(b);
    
    std::cout << "Solución: " << x << std::endl;
    std::cout << "Residuo: " << (A * x - b).NormL2() << std::endl;
    
    return 0;
}
```

**Compilar:**
```bash
g++ -std=c++17 -O3 myprogram.cpp -o myprogram
```

---

## ✨ Características Clave

### 🏗️ Arquitectura de la Biblioteca

```
┌─────────────────────────────────────────────────────────────────────┐
│                          MML.h (encabezado único)                      │
├─────────────────────────────────────────────────────────────────────┤
│                                                                     │
│   mml/                                                              │
│   ├── base/        Vectores, Matrices, Tensores, Funciones, Geometría  │
│   ├── core/        Derivación, Integración, Solvers Lineales, Campos  │
│   ├── algorithms/  EDO, Búsqueda de Raíces, Interpolación, Solvers de Valores Propios  │
│   ├── systems/     Sistemas Dinámicos, Sistemas Lineales, Atractores    │
│   ├── interfaces/  Interfaces abstractas                              │
│   └── tools/       Visualización, Serialización, Impresión en Consola   │
│                                                                     │
└─────────────────────────────────────────────────────────────────────┘
```

MML está organizado en cuatro capas principales, cada una construida sobre la anterior:

- **[Base](docs/base/README_Base.md)** — El fundamento matemático. Vectores, matrices, tensores, funciones, polinomios, cuaterniones y primitivas geométricas. Estos son los objetos con los que computes — diseñados con una sintaxis intuitiva para que las expresiones matemáticas en el código se lean de forma natural.
- **[Core](docs/core/README_Core.md)** — Operaciones numéricas sobre objetos base. Diferenciación numérica (hasta precisión de orden 8), integración (1D/2D/3D con cuadratura adaptativa), solvers de sistemas lineales (LU, QR, SVD, Cholesky), operaciones de campos vectoriales (gradiente, divergencia, rotacional, Laplaciano), transformaciones de coordenadas y tensores métricos.
- **[Algorithms](docs/algorithms/README_Algorithms.md)** — Métodos de resolución de problemas de nivel superior. Solvers de EDO (paso fijo y adaptativo), búsqueda de raíces (Bisección, Newton, Brent), descomposición de valores propios, interpolación, ajuste de curvas, integración de trayectorias y superficies, geometría diferencial y análisis de funciones.
- **[Tools](docs/tools/)** — Puente entre el cómputo y la presentación. Impresión en consola con formato de calidad de publicación (6 formatos de exportación), serialización de archivos para todos los objetos matemáticos y visualizadores multiplataforma para funciones, superficies, campos vectoriales y sistemas de partículas.


### 📐 Objetos Matemáticos

| Categoría | Tipos | Descripción |
|----------|-------|-------------|
| [**Vectores**](docs/base/Vectors.md) | `Vector<T>`, `VectorN<T,N>`, Cartesianas 2D/3D, Polares, Esféricas | Vectores de tamaño dinámico y fijo en múltiples sistemas de coordenadas |
| [**Matrices**](docs/base/Matrices.md) | `Matrix<T>`, `MatrixNM<N,M>`, Simétrica, Tridiagonal, Banda | Álgebra matricial completa con almacenamiento especializado |
| [**Tensores**](docs/base/Tensors.md) | `Tensor2<N>` hasta `Tensor5<N>` | Tensores de rango 2-5 para cálculos avanzados |
| [**Funciones**](docs/base/Functions.md) | `IRealFunction`, `IScalarFunction<N>`, `IVectorFunction<N>` | Objetos de función de primera clase con soporte de cálculo |
| [**Curvas y Superficies**](docs/core/Curves_and_surfaces.md) | `ParametricCurve<N>`, `ParametricSurface<N>`, 2D/3D predefinidas | Curvas paramétricas, superficies, longitud de arco, marcos de Frenet |
| [**Geometría**](docs/base/Geometry.md) | Puntos, Líneas, Planos, Triángulos, Cuerpos | Primitivas geométricas 2D y 3D |
| [**Polinomios**](docs/base/Polynoms.md) | `Polynom<T>` | Polinomios genéricos sobre cualquier cuerpo |
| [**Cuaterniones**](docs/base/Quaternions.md) | Álgebra completa de cuaterniones | Rotaciones 3D, interpolación SLERP |

### 🔢 Algoritmos Numéricos

| Categoría | Algoritmos | Descripción |
|----------|------------|-------------|
| [**Álgebra Lineal**](docs/core/Linear_equations_solvers.md) | LUSolver, QRSolver, SVD, Cholesky | Descomposiciones y solvers matriciales |
| [**Solvers de Valores Propios**](docs/algorithms/Eigen_solvers.md) | EigenSolver (simétrico y general) | Cálculo de valores propios de matrices reales |
| [**Derivación**](docs/core/Derivation.md) | NDer1-8, NSecDer, NThirdDer, Gradient, Jacobian | Derivadas 1ª/2ª/3ª, precisión de O(h) a O(h⁸) |
| [**Integración 1D**](docs/core/Integration.md) | Trap, Simpson, Romberg, Gauss-Kronrod (G7K15, G10K21) | Cuadratura adaptativa con estimaciones de error |
| [**Integración 2D/3D**](docs/core/Multidim_integration.md) | Integrate2D, Integrate3D, Monte Carlo | Integración multidimensional |
| [**Integrales Impropias**](docs/core/Integration.md) | IntegrateUpperInf, IntegrateLowerInf, IntegrateInfInf | Límites seminfinitos e infinitos |
| [**Trayectoria y Superficie**](docs/algorithms/Path_integration.md) | PathIntegration, SurfaceIntegration | Integrales de curvas y superficies, flujo |
| [**Solvers de EDO**](docs/algorithms/Differential_equations_solvers.md) | ODESystemFixedStepSolver, Euler, RK4, Adaptativo | Ecuaciones diferenciales ordinarias |
| [**Sistemas Dinámicos**](docs/systems/) | DynamicalSystem, puntos fijos, exponentes de Lyapunov | Análisis de espacio de fases y estabilidad |
| [**Búsqueda de Raíces**](docs/algorithms/Root_finding.md) | Bisección, Newton, Secante, Brent | Resolución de ecuaciones |
| [**Interpolación**](docs/base/Interpolated_functions.md) | LinearInterpRealFunc, SplineInterpRealFunc, PolynomInterpRealFunc | Aproximación de funciones |

### 🌀 Capacidades Avanzadas

| Característica | Descripción |
|---------|-------------|
| [**Operaciones de Campos**](docs/core/Field_operations.md) | `GradientCart`, `DivCart`, `CurlCart`, `LaplacianCart` en Cartesianas, cilíndricas, esféricas |
| [**Transformaciones de Coordenadas**](docs/core/Coordinate_transformations.md) | `CoordTransfSphericalToCartesian`, `CoordTransfLorentzXAxis`, rotaciones |
| [**Tensores Métricos**](docs/core/Metric_tensor.md) | Soporte de sistemas de coordenadas generales |
| [**Geometría Diferencial**](docs/algorithms/Differential_geometry.md) | Curvatura, torsión, marcos de Frenet |
| [**Sistemas Dinámicos**](docs/systems/) | Clasificación de puntos fijos, exponentes de Lyapunov, análisis de bifurcación |
| [**Análisis de Funciones**](docs/algorithms/Function_analyzer.md) | Encontrar raíces, extremos, puntos de inflexión |

## 📚 Herramientas y Utilidades

**Herramientas prácticas** para trabajar con objetos matemáticos — la cuarta capa de MML.

| Herramienta | Descripción |
|------|-------------|
| [**ConsolePrinter**](docs/tools/ConsolePrinter.md) | Formato de tabla hermoso, 6 formatos de exportación (TXT, CSV, JSON, HTML, LaTeX, Markdown), 5 estilos de borde |
| [**Serializer**](docs/tools/Serializer.md) | Guardar funciones, soluciones de EDO, simulaciones de partículas, campos vectoriales a archivos |
| [**Visualizadores**](docs/tools/Visualizers.md) | Gráficos multiplataforma para funciones, campos, curvas, superficies, partículas |
| **Random** | Generadores de números aleatorios para simulaciones y Monte Carlo |

**Por qué las Herramientas son Importantes:**

La capa de Herramientas es el puente entre el cómputo y la presentación. Ya sea que estés depurando algoritmos numéricos, preparando resultados para publicación o creando animaciones para enseñar — estas utilidades lo hacen sin esfuerzo:

- **ConsolePrinter**: Formatea cualquier vector, matriz o tabla en una salida de calidad de publicación con una sola llamada
- **Serializer**: Persiste los resultados de simulación para análisis posterior o visualización en herramientas externas
- **Visualizadores**: Lanza visores en tiempo real directamente desde el código — sin necesidad de exportación manual de datos



Consulta la [Galería de Visualización](#-galería-de-visualización) completa a continuación para todos los tipos de visualizador.

---

## 🧪 Ejemplos Reales

> **Ejemplos autónomos y listos para producción** que demuestran las capacidades de MML en simulaciones físicas del mundo real. Cada ejemplo incluye todo el código de física dentro del directorio — sin dependencias externas.

### 🌌 [Ejemplo 00: Gravedad N-Cuerpos](docs/examples/Example_00_N_body_gravity.md) — *Ejemplo Destacado*

**Simulaciones del Sistema Solar y Cúmulos Estelares** — Ley de gravitación universal de Newton con **7 integradores** (Euler, RK4, Verlet, Leapfrog, RK5, DP5, DP8). Integradores simplecticos para estabilidad orbital a largo plazo, métodos adaptativos para máxima exactitud. Motor de física autónomo (~870 líneas).

<table>
<tr>
<td align="center" width="33%">

![Sistema Solar](docs/images/readme/examples/00_N_body_gravity/_Example00_solar_system.png)

*Mecánica orbital del sistema solar*

</td>
<td align="center" width="33%">

![Simulación de Partículas](docs/images/readme/examples/00_N_body_gravity/_Example00_solar_system_particle_sim.png)

*Visualización de partículas en tiempo real*

</td>
<td align="center" width="33%">

![Vista General del Cúmulo](docs/images/readme/examples/00_N_body_gravity/_Example00_star_clusters_particle_overview.png)

*Vista general de colisión de cúmulos estelares*

</td>
</tr>
<tr>
<td align="center">

![Paso del Cúmulo 1](docs/images/readme/examples/00_N_body_gravity/_Example00_star_clusters_particle_step%201.png)

*Aproximación del cúmulo*

</td>
<td align="center">

![Paso del Cúmulo 3](docs/images/readme/examples/00_N_body_gravity/_Example00_star_clusters_particle_step%203.png)

*Interacción gravitacional*

</td>
<td align="center">

![Trayectorias del Cúmulo](docs/images/readme/examples/00_N_body_gravity/_Example00_star_clusters_trajectories_visualization.png)

*Visualización completa de trayectorias*

</td>
</tr>
</table>

```cpp
NBodyGravitySimConfig config = NBodyGravityConfigGenerator::Config1_Solar_system();
NBodyGravitySimulator solver(config);

auto results_verlet = solver.SolveVerlet(0.01, 365000);       // 10 años, simplectico
auto results_dp8 = solver.SolveDP8(10.0, 1e-12, 0.01, 0.001); // 10 años, DP8 adaptativo
```

---

### 🎾 [Ejemplo 02: Péndulo Doble](docs/examples/Example_02_double_pendulum.md) — *Teoría del Caos*

Caos determinista en acción — diferencias infinitamente pequeñas en las condiciones iniciales llevan a resultados completamente diferentes.

<table>
<tr>
<td align="center" width="33%">

![Trayectoria](docs/images/readme/examples/02_double_pendulum/01_trajectory_for%20both%20angles.png)

*Trayectorias angulares*

</td>
<td align="center" width="33%">

![Espacio de Fases](docs/images/readme/examples/02_double_pendulum/02_phase_space.png)

*Retrato del espacio de fases*

</td>
<td align="center" width="33%">

![Efecto Mariposa](docs/images/readme/examples/02_double_pendulum/03_butterfly_effect.png)

*Divergencia del efecto mariposa*

</td>
</tr>
</table>

---

### 🏎️ [Ejemplo 03: Análisis de G-Force de Fórmula 1](docs/examples/Example_03_F1_GForce_analysis.md) — *Curvas Paramétricas*

Datos de telemetría real de F1 (Silverstone, Monza) analizado usando cálculos de curvas paramétricas y curvatura de MML. G lateral = v²κ/g, G longitudinal = (1/g)·dv/dt.

<table>
<tr>
<td align="center" width="33%">

![Trayectoria del Circuito](docs/images/readme/examples/03_formula_1_sim/01_track_path.png)

*Disposición del circuito desde telemetría*

</td>
<td align="center" width="33%">

![G-Forces](docs/images/readme/examples/03_formula_1_sim/02_g_forces.png)

*Perfil de G-force alrededor de la vuelta*

</td>
<td align="center" width="33%">

![Perfil de Velocidad](docs/images/readme/examples/03_formula_1_sim/03_speed_profile.png)

*Análisis del perfil de velocidad*

</td>
</tr>
</table>

---

### 💥 [Ejemplo 04: Simulador de Colisiones 2D](docs/examples/Example_04_collision_simulator_2d.md) — *Teoría Cinética*

**30,000+ partículas** con física de colisión elástica exacta, particionamiento espacial para rendimiento O(N) y ejecución multi-hilo. ¡Observa cómo se propagan las ondas de choque!

<table>
<tr>
<td align="center" width="33%">

![Choque 1](docs/images/readme/examples/04_collision_simulator_2d/01_shock_wave.png)

*Frente de choque inicial*

</td>
<td align="center" width="33%">

![Choque 2](docs/images/readme/examples/04_collision_simulator_2d/02_shock_wave.png)

*Propagación de la onda*

</td>
<td align="center" width="33%">

![Choque 3](docs/images/readme/examples/04_collision_simulator_2d/03_shock_wave.png)

*Dispersión de la onda de choque*

</td>
</tr>
</table>

---

### 📦 [Ejemplo 05: Colisiones de Cuerpos Rígidos](docs/examples/Example_05_rigid_body.md) — *Dinámica 3D*

Dos paralelepípedos y una esfera en un contenedor cúbico con colisiones elásticas, dinámica rotacional completa usando cuaterniones y tensores de inercia.

<table>
<tr>
<td align="center" width="50%">

![Inicio](docs/images/readme/examples/05_rigid_body/01_rigid_body_start.png)

*Configuración inicial*

</td>
<td align="center" width="50%">

![Colisión](docs/images/readme/examples/05_rigid_body/02_rigid_body.png)

*Dinámica a mitad de la colisión*

</td>
</tr>
</table>

---

### Más Ejemplos

| # | Ejemplo | Descripción |
|---|---------|-------------|
| 01 | [**Lanzamiento de Proyectil**](docs/examples/Example_01_projectile_launch.md) | Trayectoria balística con resistencia del aire — modelos de arrastre, comparación vacío vs aire |
| 06 | [**Transformaciones de Lorentz**](docs/examples/Example_06_Lorentz_transformations.md) | Relatividad especial — dilatación del tiempo, contracción de la longitud, Paradoja del Gemelo |

### 🚀 Pruébalo Ahora

```bash
cmake -B build && cmake --build build

# Ejecuta la simulación destacada de N-Cuerpos
./build/src/examples/Release/Example00_NBodyGravity      # Windows
./build/src/examples/Example00_NBodyGravity              # Linux
```

Todos los ejemplos producen salida de visualización visible con los visores incluidos basados en Qt.

---

## 📊 Suite de Visualización


MML incluye una poderosa **suite de visualización multiplataforma** para funciones, campos, curvas, superficies y sistemas de partículas. Todos los visualizadores funcionan en Windows (WPF), Linux (Qt) y macOS (Qt).

### Galería de Visualización

#### Windows — Visualizadores WPF

<table>
<tr>
<td align="center" width="25%">

**Funciones Reales**

![WPF Real](docs/images/readme/visualization_suite/win/wpf_real_func_multi_damped_oscillations.png)

</td>
<td align="center" width="25%">

**Funciones Reales (Lorentz)**

![WPF Lorentz](docs/images/readme/visualization_suite/win/wpf_real_func_multi_Lorentz.png)

</td>
<td align="center" width="25%">

**Función Escalar 2D**

![WPF Scalar 2D](docs/images/readme/visualization_suite/win/wpf_scalar_func_2d.png)

</td>
<td align="center" width="25%">

**Función Escalar 3D**

![WPF Scalar 3D](docs/images/readme/visualization_suite/win/wpf_scalar_func_3d.png)

</td>
</tr>
<tr>
<td align="center">

**Curva Paramétrica 2D**

![WPF Curve 2D](docs/images/readme/visualization_suite/win/wpf_param_curve_2d_butterfly.png)

</td>
<td align="center">

**Curva Paramétrica 3D**

![WPF Curve 3D](docs/images/readme/visualization_suite/win/wpf_param_curve_3d.png)

</td>
<td align="center">

**Superficie Paramétrica**

![WPF Surface](docs/images/readme/visualization_suite/win/wpf_param_surface.png)

</td>
<td align="center">

**Campo Vectorial 3D**

![WPF VecField](docs/images/readme/visualization_suite/win/wpf_vector_field_3d.png)

</td>
</tr>
<tr>
<td align="center">

**Visualizador de Partículas 2D**

![WPF Particle 2D](docs/images/readme/visualization_suite/win/wpf_particle_vis_2d.png)

</td>
<td align="center">

**Visualizador de Partículas 3D**

![WPF Particle 3D](docs/images/readme/visualization_suite/win/wpf_particle_vis_3d.png)

</td>
<td align="center">

**Simulación de Cuerpo Rígido**

![WPF Rigid](docs/images/readme/visualization_suite/win/wpf_rigid_body.png)

</td>
<td align="center">

</td>
</tr>
</table>

#### Windows — Visualizadores Qt

<table>
<tr>
<td align="center" width="25%">

**Funciones Reales**

![Qt Real](docs/images/readme/visualization_suite/win/win_qt_real_func_multi.png)

</td>
<td align="center" width="25%">

**Función Escalar 2D**

![Qt Scalar 2D](docs/images/readme/visualization_suite/win/win_qt_scalar_func_2d.png)

</td>
<td align="center" width="25%">

**Función Escalar 3D**

![Qt Scalar 3D](docs/images/readme/visualization_suite/win/win_qt_scalar_func_3d.png)

</td>
<td align="center" width="25%">

**Superficie Paramétrica**

![Qt Surface](docs/images/readme/visualization_suite/win/win_qt_param_surface.png)

</td>
</tr>
<tr>
<td align="center">

**Curva Paramétrica 2D**

![Qt Curve 2D](docs/images/readme/visualization_suite/win/win_qt_param_curve_2d.png)

</td>
<td align="center">

**Curva Paramétrica 3D**

![Qt Curve 3D](docs/images/readme/visualization_suite/win/win_qt_param_curve_3d.png)

</td>
<td align="center">

**Campo Vectorial 3D**

![Qt VecField](docs/images/readme/visualization_suite/win/win_qt_vector_field_3d.png)

</td>
<td align="center">

**Simulación de Cuerpo Rígido**

![Qt Rigid](docs/images/readme/visualization_suite/win/win_qt_rigid_body.png)

</td>
</tr>
<tr>
<td align="center">

**Visualizador de Partículas 2D**

![Qt Particle 2D](docs/images/readme/visualization_suite/win/win_qt_particle_vis_2d.png)

</td>
<td align="center">

**Visualizador de Partículas 3D**

![Qt Particle 3D](docs/images/readme/visualization_suite/win/win_qt_particle_vis_3d.png)

</td>
<td align="center">

</td>
<td align="center">

</td>
</tr>
</table>

#### Linux — Visualizadores Qt

<table>
<tr>
<td align="center" width="25%">

**Funciones Reales**

![Linux Real](docs/images/readme/visualization_suite/linux/linux_qt_real_func_multri.png)

</td>
<td align="center" width="25%">

**Función Escalar 2D**

![Linux Scalar 2D](docs/images/readme/visualization_suite/linux/linux_qt_scalar_func_2d.png)

</td>
<td align="center" width="25%">

**Función Escalar 3D**

![Linux Scalar 3D](docs/images/readme/visualization_suite/linux/linux_qt_scalar_func_3d.png)

</td>
<td align="center" width="25%">

**Superficie Paramétrica**

![Linux Surface](docs/images/readme/visualization_suite/linux/linux_qt_param_surface.png)

</td>
</tr>
<tr>
<td align="center">

**Curva Paramétrica 2D**

![Linux Curve 2D](docs/images/readme/visualization_suite/linux/linux_qt_param_curve_2d.png)

</td>
<td align="center">

**Curva Paramétrica 3D**

![Linux Curve 3D](docs/images/readme/visualization_suite/linux/linux_qt_param_curve_3d.png)

</td>
<td align="center">

**Campo Vectorial 3D**

![Linux VecField](docs/images/readme/visualization_suite/linux/linux_qt_vector_field_3d.png)

</td>
<td align="center">

**Simulación de Cuerpo Rígido**

![Linux Rigid](docs/images/readme/visualization_suite/linux/linux_qt_rigid_body_vis.png)

</td>
</tr>
<tr>
<td align="center">

**Visualizador de Partículas 2D**

![Linux Particle 2D](docs/images/readme/visualization_suite/linux/linux_qt_particle_vis_2d.png)

</td>
<td align="center">

**Visualizador de Partículas 3D**

![Linux Particle 3D](docs/images/readme/visualization_suite/linux/linux_qt_particle_vis_3d.png)

</td>
<td align="center">

**Escalar 2D (Tema Oscuro)**

![Linux Scalar Dark](docs/images/readme/visualization_suite/linux/linux_qt_scalar_func_2d_dark.png)

</td>
<td align="center">

**Superficie Paramétrica (Tema Oscuro)**

![Linux Surface Dark](docs/images/readme/visualization_suite/linux/linux_qt_param_surface_dark.png)

</td>
</tr>
</table>

#### macOS — Visualizadores Qt

<table>
<tr>
<td align="center" width="25%">

**Funciones Reales**

![Mac Real](docs/images/readme/visualization_suite/mac/mac_qt_real_func_multi.png)

</td>
<td align="center" width="25%">

**Funciones Reales (Lorentz)**

![Mac Lorentz](docs/images/readme/visualization_suite/mac/mac_qt_real_func_multi_Lorentz_system.png)

</td>
<td align="center" width="25%">

**Función Escalar 2D**

![Mac Scalar 2D](docs/images/readme/visualization_suite/mac/mac_qt_scalar_func_2d_monkey_saddle.png)

</td>
<td align="center" width="25%">

**Función Escalar 3D**

![Mac Scalar 3D](docs/images/readme/visualization_suite/mac/mac_qt_scalar_function_3d_gyroid.png)

</td>
</tr>
<tr>
<td align="center">

**Curva Paramétrica 2D**

![Mac Curve 2D](docs/images/readme/visualization_suite/mac/mac_qt_param_curve_2de_butterfly.png)

</td>
<td align="center">

**Curva Paramétrica 3D**

![Mac Curve 3D](docs/images/readme/visualization_suite/mac/mac_qt_param_curve_3d_trefoil.png)

</td>
<td align="center">

**Superficie Paramétrica**

![Mac Surface](docs/images/readme/visualization_suite/mac/mac_qt_param_surface_klein.png)

</td>
<td align="center">

**Campo Vectorial 3D**

![Mac VecField](docs/images/readme/visualization_suite/mac/mac_qt_vector_field_3d_gravity.png)

</td>
</tr>
<tr>
<td align="center">

**Visualizador de Partículas 2D**

![Mac Particle 2D](docs/images/readme/visualization_suite/mac/mac_qt_particle_visualizer_2d.png)

</td>
<td align="center">

**Visualizador de Partículas 3D**

![Mac Particle 3D](docs/images/readme/visualization_suite/mac/mac_qt_partice_visualizer_3d.png)

</td>
<td align="center">

**Función Escalar 2D (Ripple)**

![Mac Ripple](docs/images/readme/visualization_suite/mac/mac_qt_scalar_func_2d_ripple.png)

</td>
<td align="center">

**Simulación de Cuerpo Rígido**

![Mac Rigid](docs/images/readme/visualization_suite/mac/mac_qt_rigid_body.png)

</td>
</tr>
</table>

### Visualizadores Disponibles

Todos los visualizadores tienen **soporte multiplataforma completo** con múltiples opciones de backend:
- **Windows**: Visualizadores basados en WPF (principal), Qt y FLTK también disponibles
- **Linux**: Visualizadores basados en Qt (principal), FLTK para 2D ligero
- **macOS**: Visualizadores basados en Qt (principal), FLTK para 2D ligero

| Visualizador | Propósito | Salida | Plataforma |
|------------|---------|--------|----------|
| **RealFunctionVisualizer** | Trazar funciones 1D | Gráficos 2D | Windows (WPF), Linux/macOS (Qt/FLTK) |
| **MultiRealFunctionVisualizer** | Comparar múltiples funciones 1D | Gráficos 2D superpuestos | Windows (WPF), Linux/macOS (Qt/FLTK) |
| **ScalarFunction2DVisualizer** | Gráficos de superficie 3D | Superficies 3D | Windows (WPF), Linux/macOS (Qt) |
| **ParametricCurve2DVisualizer** | Curvas paramétricas 2D | Curvas 2D | Windows (WPF), Linux/macOS (Qt/FLTK) |
| **ParametricCurve3DVisualizer** | Curvas paramétricas 3D | Curvas 3D | Windows (WPF), Linux/macOS (Qt) |
| **VectorField2DVisualizer** | Flechas de campo vectorial 2D | Gráficos de flechas 2D | Windows (WPF), Linux/macOS (Qt/FLTK) |
| **VectorField3DVisualizer** | Flechas de campo vectorial 3D | Gráficos de flechas 3D | Windows (WPF), Linux/macOS (Qt) |
| **ParticleVisualizer2D** | Sistemas de partículas 2D animados | Animaciones 2D con reproducción | Windows (WPF), Linux/macOS (Qt) |
| **ParticleVisualizer3D** | Sistemas de partículas 3D animados | Animaciones 3D con reproducción | Windows (WPF), Linux/macOS (Qt) |


### Ejemplos de Visualización
- `src/visualization_examples/` — Demos listas para ejecutar para cada tipo de visualizador
- `src/book/chapters/Chapter_02_vizualization/` — Ejemplos completos del libro complementario
- Ejemplos de estilos personalizados, diseños de múltiples gráficos y controles de animación

**Trazado de Funciones Reales — Serie Temporal del Sistema de Lorenz:**

```cpp
// Atractor de Lorenz — evolución temporal caótica de x(t), y(t), z(t)
ODESystem lorenz_system(3, [](Real t, const Vector<Real>& x, Vector<Real>& dxdt) {
    const Real sigma = 10.0, rho = 28.0, beta = 8.0 / 3.0;
    dxdt[0] = sigma * (x[1] - x[0]);
    dxdt[1] = x[0] * (rho - x[2]) - x[1];
    dxdt[2] = x[0] * x[1] - beta * x[2];
});

Vector<Real> initial_state({ 1.0, 1.0, 1.0 });
CashKarpIntegrator solver(lorenz_system);
ODESystemSolution sol = solver.integrate(initial_state, 0.0, 50.0, 0.001, 1e-10, 0.001);

// Extraer componentes de la solución como funciones reales mediante interpolación spline
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
| ![Lorenz Win](docs/images/readme/visualization_suite/win/wpf_real_func_multi_Lorentz.png) | ![Lorenz Linux](docs/images/readme/visualization_suite/linux/linux_qt_real_func_multri.png) | ![Lorenz Mac](docs/images/readme/visualization_suite/mac/mac_qt_real_func_multi_Lorentz_system.png) |

**Función Escalar 2D — Ripple (2D Sinc / Sombrero):**

```cpp
// Función 2D Sinc — ondulaciones concéntricas con pico central
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
| ![Ripple Win](docs/images/readme/visualization_suite/win/wpf_scalar_func_2d.png) | ![Ripple Linux](docs/images/readme/visualization_suite/linux/linux_qt_scalar_func_2d.png) | ![Ripple Mac](docs/images/readme/visualization_suite/mac/mac_qt_scalar_func_2d_ripple.png) |

**Curvas Paramétricas — Nudo Trefoil:**

```cpp
// Nudo trefoil: x = sin(t) + 2sin(2t), y = cos(t) - 2cos(2t), z = -sin(3t)
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
| ![Trefoil Win](docs/images/readme/visualization_suite/win/wpf_param_curve_3d.png) | ![Trefoil Linux](docs/images/readme/visualization_suite/linux/linux_qt_param_curve_3d.png) | ![Trefoil Mac](docs/images/readme/visualization_suite/mac/mac_qt_param_curve_3d_trefoil.png) |

**Campos Vectoriales — Gravedad de Dos Cuerpos:**

```cpp
// Campo gravitatorio de dos masas
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
| ![Field Win](docs/images/readme/visualization_suite/win/wpf_vector_field_3d.png) | ![Field Linux](docs/images/readme/visualization_suite/linux/linux_qt_vector_field_3d.png) | ![Field Mac](docs/images/readme/visualization_suite/mac/mac_qt_vector_field_3d_gravity.png) |

### Exportación de Datos

Todos los visualizadores usan archivos de datos serializados que también pueden ser cargados por herramientas externas (Python, MATLAB, etc.):

```cpp
// Serializar datos de función
Serializer::SaveRealFunc(f, "Function", 0, 10, 100, "data.txt");

// Serializar solución de EDO
Serializer::SaveODESolutionAsMultiFunc(
    solution, "ODE Solution", {"x", "v"}, "ode_data.txt");

// Serializar campo vectorial 2D
Serializer::SaveVectorFunc2DCartesian(
    field, "Vector Field", -10, 10, 20, -10, 10, 20, "field_data.txt");
```


---

## 📝 Ejemplos de Código

### Vectores y Matrices

📄 *[Ver código fuente completo](src/readme_examples/readme02_vectors_matrices.cpp)*

```cpp
#include "MML.h"
using namespace MML;

// Vectores y matrices reales
Vector<Real> vec1{1.5, -2.0, 0.5}, vec2{1.0, 1.0, -3.0};
Matrix<Real> mat_3x3{3, 3, {1.0,  2.0, -1.0,
                           -1.0,  5.0,  6.0,
                            3.0, -2.0,  1.0}};

// Aritmética de vectores y matrices
Vector<Real> result = Real{2.0} * (vec1 + vec2) * mat_3x3 / vec1.NormL2();
std::cout << "Resultado: " << result << std::endl;

// Vectores y matrices complejos
VectorComplex vec_cmplx{{Complex(1,1), Complex(-1,2)}};
MatrixComplex mat_cmplx{2, 2, {Complex(0.5,1),  Complex(-1,2),
                               Complex(-1,-2), Complex(-2,2)}};

// Propiedades de la matriz
std::cout << "IsOrthogonal:  " << Utils::IsOrthogonal(mat_3x3) << std::endl;
std::cout << "IsHermitian:   " << Utils::IsHermitian(mat_cmplx) << std::endl;

/* SALIDA ESPERADA:
    Resultado: [   -3.137858162,     3.922322703,    -8.629109946]
    IsOrthogonal:  0
    IsHermitian:   1
*/
```

---

### Sistemas Lineales y Valores Propios

📄 *[Ver código fuente completo](src/readme_examples/readme03_linear_systems.cpp)*

```cpp
// Definir un sistema lineal Ax = b
Matrix<Real> A{5, 5, {0.2,  4.1, -2.1, -7.4,  1.6,
                      1.6,  1.5, -1.1,  0.7,  5.0,
                     -3.8, -8.0,  9.6, -5.4, -7.8,
                      4.6, -8.2,  8.4,  0.4,  8.0,
                     -2.6,  2.9,  0.1, -9.6, -2.7}};
Vector<Real> b{1.1, 4.7, 0.1, 9.3, 0.4};

// Resolver usando descomposición LU
LUSolver<Real> luSolver(A);
Vector<Real> x = luSolver.Solve(b);

std::cout << "Solución: " << x << std::endl;
std::cout << "Verificación A*x: " << (A * x) << std::endl;
std::cout << "Norma del residuo: " << (A*x - b).NormL2() << std::endl;

// Calcular valores y vectores propios
auto eigenResult = EigenSolver::Solve(A);

std::cout << "Valores propios:" << std::endl;
for (const auto& ev : eigenResult.eigenvalues)
    std::cout << "  λ = " << ev.real << " + " << ev.imag << "i" << std::endl;

// Descomposición QR
QRSolver<Real> qr(A);
Matrix<Real> Q = qr.GetQ();
Matrix<Real> R = qr.GetR();
std::cout << "QR válido: " << Matrix<Real>::AreEqual(A, Q * R, 1e-10) << std::endl;
```

---

### Polinomios y Álgebra

📄 *[Ver código fuente completo](src/readme_examples/readme09_polynomials.cpp)*

```cpp
// Crear polinomios: p(x) = 2x³ - 3x² + x - 5
PolynomReal p{-5, 1, -3, 2};  // Coeficientes: [constante, x, x², x³]

// Evaluar en un punto
Real val = p(2.0);  // p(2) = 2(8) - 3(4) + 2 - 5 = 16 - 12 + 2 - 5 = 1
std::cout << "p(2) = " << val << std::endl;

// Aritmética de polinomios
PolynomReal q{1, 2};           // q(x) = 2x + 1
PolynomReal sum = p + q;       // Suma
PolynomReal prod = p * q;      // Multiplicación
std::cout << "Grado de p*q: " << prod.degree() << std::endl;

// Cálculo con polinomios
PolynomReal dp = p.derivative();        // p'(x) = 6x² - 6x + 1
PolynomReal ip = p.integral();          // ∫p(x)dx (constante = 0)
std::cout << "p'(x) en x=1: " << dp(1.0) << std::endl;

// Resolver ecuación cuadrática: x² - 5x + 6 = 0 → x = 2, 3
Complex r1, r2;
int numReal = SolveQuadratic(1.0, -5.0, 6.0, r1, r2);
std::cout << "Raíces: " << r1.real() << ", " << r2.real() << std::endl;

// Resolver ecuación cúbica: x³ - 6x² + 11x - 6 = 0 → x = 1, 2, 3
Complex c1, c2, c3;
SolveCubic(1.0, -6.0, 11.0, -6.0, c1, c2, c3);
std::cout << "Raíces cúbicas: " << c1.real() << ", " << c2.real() 
          << ", " << c3.real() << std::endl;
```

---

### Definiendo Funciones

📄 *[Ver código fuente completo](src/readme_examples/readme04_functions_interpolation.cpp)*

MML proporciona múltiples formas de crear objetos de función que pueden ser derivados, integrados y analizados:

```cpp
// =====================================================================
// CASO 1: Desde función independiente
// =====================================================================
double MyFunction(double x) { return sin(x) * (1.0 + 0.5*x*x); }
RealFunction f1(MyFunction);

// =====================================================================
// CASO 2: Creación directa con lambda (más común y recomendado)
// =====================================================================
RealFunction f2{[](Real x) { return sin(x) * (1.0 + 0.5*x*x); }};

// Diferentes tipos de funciones — todos creados con lambdas
ScalarFunction<3> scalar([](const VectorN<Real,3>& x) { 
    return x[0]*x[0] + x[1]*x[1] + x[2]*x[2];   // campo escalar: r²
});

VectorFunction<3> vector([](const VectorN<Real,3>& x) -> VectorN<Real,3> { 
    return {x[1]*x[2], x[0]*x[2], x[0]*x[1]};   // campo vectorial
});

ParametricCurve<3> helix([](Real t) -> VectorN<Real,3> { 
    return {cos(t), sin(t), 0.2*t};             // curva paramétrica 3D
});

ParametricSurface<3> torus([](Real u, Real v) -> VectorN<Real,3> { 
    Real R = 3.0, r = 1.0;                      // superficie paramétrica
    return {(R + r*cos(v))*cos(u), (R + r*cos(v))*sin(u), r*sin(v)}; 
});

// =====================================================================
// CASO 3a: Desde clase con operator() — útil para funciones con estado
// =====================================================================
class FunctionFromClassOperator {
    double _amplitude;
public:
    FunctionFromClassOperator(double amp) : _amplitude(amp) {}
    double operator()(double x) const { return _amplitude * std::sin(x); }
};

FunctionFromClassOperator obj(2.5);  // amplitud = 2.5
RealFunctionFromStdFunc f3(std::function<double(double)>{obj});

// =====================================================================
// CASO 3b: Clase heredando IRealFunction directamente
// =====================================================================
class MyDerivedFunction : public IRealFunction {
    double _frequency;
public:
    MyDerivedFunction(double freq) : _frequency(freq) {}
    double operator()(double x) const override { return std::cos(_frequency * x); }
};
MyDerivedFunction f4(2.0);  // frecuencia = 2.0 → cos(2x)

// =====================================================================
// CASO 4: Adaptador para clase externa/legada que no puedes modificar
// =====================================================================
class ExternalComplexClass {
public:
    double ComputeValue(double x) const { return std::exp(-x*x); }
};

class ExternalClassWrapper : public IRealFunction {
    const ExternalComplexClass& _ref;
public:
    ExternalClassWrapper(const ExternalComplexClass& obj) : _ref(obj) {}
    double operator()(double x) const override { return _ref.ComputeValue(x); }
};

ExternalComplexClass externalObj;
ExternalClassWrapper f5(externalObj);  // ¡Ahora usable con todos los algoritmos de MML!

// Todas estas funciones ahora pueden ser derivadas, integradas, analizadas...
std::cout << "f2(1) = " << f2(1.0) << std::endl;           // 1.26147
std::cout << "f2'(1) = " << Derivation::NDer4(f2, 1.0) << std::endl;  // derivada
std::cout << "∫f2 = " << IntegrateSimpson(f2, 0, 1).value << std::endl;  // integral
```

---

### Interpolación

📄 *[Ver código fuente completo](src/readme_examples/readme04_functions_interpolation.cpp)*

Crea funciones suaves a partir de puntos de datos discretos:

```cpp
// Puntos de datos (podrían ser de un experimento, simulación, archivo...)
Vector<Real> x_data{0, 2.5, 5.0, 7.5, 10.0};
Vector<Real> y_data{0, 1.5, 4.0, 9.5, 16.0};

// Tres métodos de interpolación — cada uno crea una función invocable
LinearInterpRealFunc   linear_interp(x_data, y_data);      // Lineal por tramos
SplineInterpRealFunc   spline_interp(x_data, y_data);      // Spline cúbico (suave)
PolynomInterpRealFunc  poly_interp(x_data, y_data, 3);     // Polinomio grado 3

// Evaluar en cualquier punto dentro del rango de interpolación
Real x = 3.7;
std::cout << "Lineal:     " << linear_interp(x) << std::endl;  // 2.70
std::cout << "Spline:     " << spline_interp(x) << std::endl;  // 2.41 (más suave)
std::cout << "Polinomio: " << poly_interp(x) << std::endl;    // 2.43

// Comparación en múltiples puntos:
//     x     Lineal   Spline   Polinomio
//   -----  -------  -------  -------
//    1.0    0.600    0.576    0.480
//    3.0    2.000    1.842    1.760
//    5.0    4.000    4.000    4.000
//    7.0    7.400    6.954    7.680
//    9.0   12.200   12.654   14.240
```

---

### Derivadas Numéricas

📄 *[Ver código fuente completo](src/readme_examples/readme05_numerical_calculus.cpp)*

```cpp
// Comparar órdenes de derivada en f(x) = sin(x) en x = 1.0
// Derivada analítica: f'(1) = cos(1) ≈ 0.5403023058681398

RealFunction f{[](Real x) { return std::sin(x); }};
Real x = 1.0;
Real analytical = std::cos(1.0);

Real der1 = Derivation::NDer1(f, x);   // O(h) - diferencia hacia adelante
Real der2 = Derivation::NDer2(f, x);   // O(h²) - diferencia central
Real der4 = Derivation::NDer4(f, x);   // O(h⁴) - stencil de 5 puntos
Real der6 = Derivation::NDer6(f, x);   // O(h⁶) - stencil de 7 puntos
Real der8 = Derivation::NDer8(f, x);   // O(h⁸) - stencil de 9 puntos

std::cout << std::scientific << std::setprecision(10);
std::cout << "Error NDer1: " << std::abs(der1 - analytical) << std::endl;  // ~1e-8
std::cout << "Error NDer2: " << std::abs(der2 - analytical) << std::endl;  // ~1e-11
std::cout << "Error NDer8: " << std::abs(der8 - analytical) << std::endl;  // ~1e-14
```

---

### Integración Numérica

📄 *[Ver código fuente completo](src/readme_examples/readme05_numerical_calculus.cpp)*

```cpp
// Integrar f(x) = sin(x) de 0 a π — Resultado analítico: 2.0
RealFunction f{[](Real x) { return std::sin(x); }};
Real a = 0.0, b = Constants::PI;

auto trap = IntegrateTrap(f, a, b, 1e-8);
auto simp = IntegrateSimpson(f, a, b, 1e-8);
auto romb = IntegrateRomberg(f, a, b, 1e-10);

std::cout << "Trapecio: " << trap.value << " (error: " << trap.error_estimate << ")\n";
std::cout << "Simpson:   " << simp.value << " (error: " << simp.error_estimate << ")\n";
std::cout << "Romberg:   " << romb.value << " (error: " << romb.error_estimate << ")\n";

// Monte Carlo para integrales de alta dimensión
ScalarFunction<3> volume_func([](const VectorN<Real, 3>& v) {
    return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
});
MonteCarloIntegrator<3> mc_integrator;
VectorN<Real, 3> lower{0, 0, 0}, upper{1, 1, 1};
auto mc_result = mc_integrator.integrate(volume_func, lower, upper, 
                                          MonteCarloConfig().samples(100000));
std::cout << "Estimación MC: " << mc_result.value 
          << " +/- " << mc_result.error_estimate << std::endl;
```

---

### Algoritmos de Búsqueda de Raíces

📄 *[Ver código fuente completo](src/readme_examples/readme10_root_finding.cpp)*

```cpp
// Encontrar raíz de f(x) = x³ - 2x - 5 (tiene raíz cerca de x ≈ 2.0945)
RealFunction f{[](Real x) { return x*x*x - 2*x - 5; }};

// Comparar diferentes métodos
Real root_bisect = RootFinding::FindRootBisection(f, 2.0, 3.0, 1e-12);
Real root_brent  = RootFinding::FindRootBrent(f, 2.0, 3.0, 1e-12);
Real root_newton = RootFinding::FindRootNewton(f, 2.0, 3.0, 1e-12);
Real root_ridder = RootFinding::FindRootRidders(f, 2.0, 3.0, 1e-12);

std::cout << std::setprecision(15);
std::cout << "Bisección: " << root_bisect << std::endl;
std::cout << "Brent:     " << root_brent << std::endl;
std::cout << "Newton:    " << root_newton << std::endl;
std::cout << "Ridders:   " << root_ridder << std::endl;

// Obtener información detallada de convergencia usando configuración
RootFinding::RootFindingConfig config;
config.tolerance = 1e-14;
config.max_iterations = 100;

auto result = RootFinding::FindRootBrent(f, 2.0, 3.0, config);
std::cout << "Iteraciones: " << result.iterations_used << std::endl;
std::cout << "f(raíz) =   " << result.function_value << std::endl;

// Encontrar múltiples raíces mediante búsqueda de intervalos
Vector<Real> brackets_lo, brackets_hi;
int numRoots = RootFinding::FindRootBrackets(f, -10.0, 10.0, 100, 
                                              brackets_lo, brackets_hi);
std::cout << "Encontrados " << numRoots << " intervalo(s) con raíz" << std::endl;
```

---

### Operaciones de Campos

📄 *[Ver código fuente completo](src/readme_examples/readme06_field_operations.cpp)*

```cpp
// Campo escalar: potencial gravitacional φ(x,y,z) = -1/r
ScalarFunction<3> potential([](const VectorN<Real, 3>& x) {
    return -1.0 / x.NormL2();
});
VectorN<Real, 3> pos{1.0, 2.0, 2.0};

// Gradiente ∇φ (da la dirección de la fuerza)
auto grad = ScalarFieldOperations::GradientCart<3>(potential, pos);
std::cout << "Gradiente en (1,2,2): " << grad << std::endl;

// Laplaciano ∇²φ (cero fuera de la masa para la gravedad!)
Real laplacian = ScalarFieldOperations::LaplacianCart<3>(potential, pos);
std::cout << "Laplaciano en (1,2,2): " << laplacian << std::endl;

// Campo vectorial: campo de velocidad rotacional v = (y, -x, z)
VectorFunction<3> velocity([](const VectorN<Real, 3>& x) -> VectorN<Real, 3> {
    return {x[1], -x[0], x[2]};
});

// Divergencia ∇·v (tasa de compresión/expansión)
Real div = VectorFieldOperations::DivCart<3>(velocity, pos);
std::cout << "Divergencia en (1,2,2): " << div << std::endl;

// Rotacional ∇×v (rotación/vorticidad)
auto curl = VectorFieldOperations::CurlCart(velocity, pos);
std::cout << "Rotacional en (1,2,2): " << curl << std::endl;
```

---

### Transformaciones de Coordenadas

📄 *[Ver código fuente completo](src/readme_examples/readme11_coord_transforms.cpp)*

```cpp
// Punto en coordenadas cartesianas
Vector3Cartesian cart_pos{1.0, 1.0, 1.0};

// Convertir a Esférica (r, θ, φ) - convención Matemática/ISO
// θ = ángulo polar desde el eje z, φ = ángulo azimutal en el plano xy
Vector3Spherical sph_pos = CoordTransfCartToSpher.transf(cart_pos);
std::cout << "Esférica: r=" << sph_pos[0] << ", θ=" << sph_pos[1] 
          << ", φ=" << sph_pos[2] << std::endl;

// Convertir a Cilíndrica (r, φ, z)
Vector3Cylindrical cyl_pos = CoordTransfCartToCyl.transf(cart_pos);
std::cout << "Cilíndrica: r=" << cyl_pos[0] << ", φ=" << cyl_pos[1] 
          << ", z=" << cyl_pos[2] << std::endl;

// Convertir de vuelta a Cartesianas
Vector3Cartesian back = CoordTransfSpherToCart.transf(sph_pos);
std::cout << "Vuelta a Cartesianas: " << back << std::endl;

// Operaciones de campo en coordenadas esféricas
// Potencial de inverso del cuadrado: φ = -1/r
ScalarFunction<3> pot_spher([](const VectorN<Real, 3>& x) { return -1.0/x[0]; });

// Usar un punto alejado del origen para el gradiente
Vector3Spherical test_sph{2.0, Constants::PI/4, Constants::PI/4};
auto grad_spher = ScalarFieldOperations::GradientSpher(pot_spher, test_sph);
std::cout << "Gradiente de -1/r en esférico en r=2:" << std::endl;
std::cout << "  ∂φ/∂r = " << grad_spher[0] << " (analítico: 1/r² = 0.25)" << std::endl;

// Transformación de vector covariante (transformar gradiente a coordenadas cartesianas)
Vector3Cartesian test_cart = CoordTransfSpherToCart.transf(test_sph);
auto force_cart = CoordTransfSpherToCart.transfVecCovariant(grad_spher, test_cart);
std::cout << "Vector fuerza en Cartesianas: " << force_cart << std::endl;
```

---

### Curvas Paramétricas

📄 *[Ver código fuente completo](src/readme_examples/readme12_parametric_curves.cpp)*

```cpp
// Definir una hélice 3D: r(t) = (cos(t), sin(t), 0.2t)
Curves::CurveCartesian3D helix([](Real t) -> VectorN<Real, 3> { 
    return {cos(t), sin(t), 0.2*t}; 
});

Real t = Constants::PI / 4;

// Propiedades de la curva en el parámetro t
auto pos = helix(t);                    // Posición sobre la curva
auto tangent = helix.getTangent(t);     // Vector tangente dr/dt
auto unit_tan = helix.getTangentUnit(t);// Tangente unitaria T
auto normal = helix.getNormal(t);       // Vector normal (aceleración)
auto binormal = helix.getBinormal(t);   // Binormal B = T × N

std::cout << std::setprecision(6);
std::cout << "Posición:       " << pos << std::endl;
std::cout << "Tangente dr/dt:  " << tangent << std::endl;
std::cout << "Tangente unitaria:   " << unit_tan << std::endl;
std::cout << "Normal (acel): " << normal << std::endl;
std::cout << "Binormal:       " << binormal << std::endl;

// Curvatura κ (aparato de Frenet-Serret)
Real curvature = helix.getCurvature(t);
std::cout << "Curvatura κ:    " << curvature << std::endl;

// Curvas predefinidas
Curves::LemniscateCurve lemniscate;       // Curva en forma de ocho
Curves::ToroidalSpiralCurve torus(5, 2);  // Espiral sobre un toro
Curves::Circle3DXZCurve circle(3.0);      // Círculo en el plano XZ
```

---

### Integrales de Trayectoria y Línea

📄 *[Ver código fuente completo](src/readme_examples/readme13_path_integrals.cpp)*

```cpp
// Integral de línea de un campo vectorial a lo largo de una curva
// ∫ F·dr donde F = (y, -x, 0) a lo largo del círculo unitario

VectorFunction<3> F([](const VectorN<Real, 3>& p) -> VectorN<Real, 3> {
    return {p[1], -p[0], 0.0};
});

// Círculo unitario en el plano XY: r(t) = (cos(t), sin(t), 0), t ∈ [0, 2π]
ParametricCurve<3> circle([](Real t) -> VectorN<Real, 3> {
    return {cos(t), sin(t), 0.0};
});

// Integral de trabajo: W = ∫₀^{2π} F(r(t)) · r'(t) dt
Real work = PathIntegration::LineIntegral(F, circle, 0.0, 2*Constants::PI, 1e-8);
std::cout << "Integral de trabajo ∮ F·dr:" << std::endl;
std::cout << "  Numérica:  " << work << std::endl;
std::cout << "  Analítica: -2π = " << -2*Constants::PI << std::endl;

// Cálculo de longitud de arco: ∫ ds
Real arc_length = PathIntegration::ParametricCurveLength<3>(circle, 0.0, 2*Constants::PI);
std::cout << "Longitud de arco del círculo unitario:" << std::endl;
std::cout << "  Numérica:  " << arc_length << std::endl;
std::cout << "  Analítica: 2π = " << 2*Constants::PI << std::endl;
```

---

### Análisis de Funciones

📄 *[Ver código fuente completo](src/readme_examples/readme14_function_analysis.cpp)*

```cpp
// Analizar comportamiento de función en un intervalo
RealFunction f{[](Real x) { return x*x*x - 3*x + 1; }};

RealFunctionAnalyzer analyzer(f, "x³ - 3x + 1");
analyzer.PrintIntervalAnalysis(-3.0, 3.0, 100, 1e-6);
/*  Salida:
    f(x) = x³ - 3x + 1 - Análisis en intervalo [-3.00, 3.00]:
      Definida    : sí
      Continua : sí
      Monótona  : no
      Mín        : -0.999...
      Máx        : 2.999...
*/

// Encontrar raíces usando búsqueda de intervalos + bisección
Vector<Real> xb1, xb2;
int numBrackets = RootFinding::FindRootBrackets(f, -3.0, 3.0, 100, xb1, xb2);
std::cout << "Encontrados " << numBrackets << " intervalo(s) con raíz:" << std::endl;
std::cout << std::setprecision(10);
for (int i = 0; i < numBrackets; i++) {
    Real root = RootFinding::FindRootBisection(f, xb1[i], xb2[i], 1e-10);
    std::cout << "  x = " << root << ", f(x) = " << f(root) << std::endl;
}

// Analizar una función con discontinuidad
RealFunctionFromStdFunc step([](Real x) -> Real { 
    if (x < 0) return 0.0;
    else if (x > 0) return 1.0;
    else return 0.5;
});
RealFunctionAnalyzer step_analyzer(step, "step(x)");
step_analyzer.PrintIntervalAnalysis(-2.0, 2.0, 100, 1e-6);
// Detecta discontinuidad en x = 0
```

---

### Ecuaciones Diferenciales

📄 *[Ver código fuente completo](src/readme_examples/readme07_ode_solvers.cpp)*

```cpp
// Oscilador armónico simple: d²x/dt² = -ω²x
// Como sistema: dx/dt = v, dv/dt = -ω²x
const Real omega = 2.0;

ODESystem system(2, [](Real t, const Vector<Real>& y, Vector<Real>& dydt) {
    const Real omega = 2.0;
    dydt[0] = y[1];                    // dx/dt = v
    dydt[1] = -omega*omega * y[0];     // dv/dt = -ω²x
});

Vector<Real> initial_cond{1.0, 0.0};  // x(0) = 1, v(0) = 0

// Resolver con paso RK4
RungeKutta4_StepCalculator stepper;
ODESystemFixedStepSolver solver(system, stepper);
ODESystemSolution solution = solver.integrate(initial_cond, 0.0, 10.0, 1000);

// Acceder a valores finales (índice 999 para 1000 puntos)
int last = solution.size() - 1;
std::cout << "Posición final: " << solution.getXValue(last, 0) << std::endl;
std::cout << "Velocidad final: " << solution.getXValue(last, 1) << std::endl;

// Solución analítica en t=10: x(t) = cos(ωt), v(t) = -ω sin(ωt)
Real t_final = 10.0;
std::cout << "x analítica:   " << std::cos(omega * t_final) << std::endl;
std::cout << "v analítica:   " << -omega * std::sin(omega * t_final) << std::endl;
```

---

### 🦋 Análisis de Sistemas Dinámicos — *La Joya de la Corona*

📄 *[Ver código fuente completo](src/readme_examples/readme15_dynamical_systems.cpp)*

```cpp
// EL SISTEMA DE LORENZ - El efecto mariposa en acción
// dx/dt = σ(y - x), dy/dt = x(ρ - z) - y, dz/dt = xy - βz
using namespace MML::Systems;

LorenzSystem lorenz(10.0, 28.0, 8.0/3.0);  // Parámetros caóticos clásicos

std::cout << "Disipativo: " << (lorenz.isDissipative() ? "sí" : "no") << std::endl;
std::cout << "Divergencia del flujo: " << lorenz.getDivergence() << " (siempre negativa -> existe atractor)" << std::endl;

// ANÁLISIS DE PUNTOS FIJOS - Encontrar equilibrios y clasificar estabilidad
std::vector<Vector<Real>> guesses = {
    Vector<Real>{0.0, 0.0, 0.0},      // Origen
    Vector<Real>{8.0, 8.0, 27.0},     // Cerca de C+
    Vector<Real>{-8.0, -8.0, 27.0}    // Cerca de C-
};

auto fixedPoints = FixedPointFinder::FindMultiple(lorenz, guesses);

std::cout << "Encontrados " << fixedPoints.size() << " puntos fijos:" << std::endl;
for (const auto& fp : fixedPoints) {
    std::cout << "  Punto Fijo: (" << fp.location[0] << ", " 
              << fp.location[1] << ", " << fp.location[2] << ")" << std::endl;
    std::cout << "    Tipo: " << ToString(fp.type) << std::endl;
    std::cout << "    Estable: " << (fp.isStable ? "sí" : "NO (inestable)") << std::endl;
    std::cout << "    Valores propios: ";
    for (const auto& ev : fp.eigenvalues)
        std::cout << "(" << ev.real() << "+" << ev.imag() << "i) ";
    std::cout << std::endl;
}

// EXPONENTES DE LYAPUNOV - Cuantificar el caos mediante sensibilidad a condiciones iniciales
Vector<Real> x0 = lorenz.getDefaultInitialCondition();
auto lyapResult = LyapunovAnalyzer::Compute(lorenz, x0, 
                                             500.0,   // Tiempo total de integración
                                             1.0,     // Intervalo de ortonormalización
                                             0.01);   // Tamaño de paso

std::cout << "Espectro de Lyapunov: [" << lyapResult.exponents[0] << ", "
          << lyapResult.exponents[1] << ", " << lyapResult.exponents[2] << "]" << std::endl;
std::cout << "Exponente máximo: " << lyapResult.maxExponent;
if (lyapResult.maxExponent > 0) std::cout << " (POSITIVO -> ¡CAOS!)";
std::cout << std::endl;
std::cout << "Dimensión de Kaplan-Yorke: " << lyapResult.kaplanYorkeDimension 
          << " (atractor fractal!)" << std::endl;
std::cout << "El sistema es " << (lyapResult.isChaotic ? "CAÓTICO" : "regular") << std::endl;

// COMPARAR SISTEMAS CAÓTICOS CLÁSICOS
RosslerSystem rossler(0.2, 0.2, 5.7);       // Caos en espiral
VanDerPolSystem vanderpol(1.0);             // Ciclo límite (NO caótico)
DoublePendulumSystem pendulum(1.0, 1.0, 9.81);  // Caos mecánico

// Calcular sus espectros de Lyapunov...
// Lorenz:       λ₁ ≈ +0.91  (CAÓTICO)
// Rössler:      λ₁ ≈ +0.07  (CAÓTICO)
// Van der Pol:  λ₁ ≈  0.00  (periódico)
// Péndulo Doble: λ₁ ≈ +1.5   (CAÓTICO)

// ANÁLISIS DE BIFURCACIÓN - Ruta al caos mediante barridos de parámetros
LorenzSystem lorenzSweep;
Vector<Real> sweepIC{1.0, 1.0, 1.0};

auto bifurcation = BifurcationAnalyzer::Sweep(
    lorenzSweep,
    1,              // Índice del parámetro (rho)
    20.0, 30.0,     // Rango del parámetro
    6,              // Número de valores del parámetro
    sweepIC,
    2,              // Registrar máximos de la componente z
    50.0,           // Tiempo transitorio
    20.0,           // Tiempo de registro
    0.01            // Tamaño de paso
);

// Revela: rho < 24: periódico → cascada de duplicación de período → rho ≈ 24.74: ¡inicio del caos!
```

---

## 🔬 Precisión y Pruebas

MML toma la exactitud numérica en serio. Cada algoritmo es rigurosamente validado contra soluciones analíticas, y las características de precisión están documentadas.

### Resumen del Suite de Pruebas

**4,540 pruebas unitarias** en **111 archivos de prueba**, organizadas por dominio:

| Dominio | Archivos | Validaciones Clave |
|--------|-------|-----------------|
| **Álgebra Lineal** | 18 | Descomposiciones LU/QR/SVD, solvers de valores propios, números de condición hasta 10¹⁵ |
| **Cálculo** | 14 | Derivación (órdenes 1-8), integración (1D/2D/3D), Gauss-Kronrod |
| **Solvers de EDO** | 6 | Todos los pasos, detección de eventos, sistemas rígidos (λ = -10⁶) |
| **Geometría** | 22 | Primitivas 2D/3D, envolvente convexa, Voronoi, KD-tree, triangulación |
| **Funciones** | 12 | Real, escalar, vectorial, curvas/superficies paramétricas, interpolación |
| **Operaciones de Campos** | 8 | Gradiente, divergencia, rotacional en Cartesianas/esféricas/cilíndricas |
| **Búsqueda de Raíces** | 4 | Bisección, Newton, Brent, Ridders con verificación de convergencia |

### Bancos de Pruebas Predefinidos

El espacio de nombres `TestBeds::` proporciona **entradas probadas y confiables** para la validación de algoritmos:

| Banco de Pruebas | Contenidos | Propósito |
|----------|----------|---------|
| **Sistemas Lineales** | Hilbert, Vandermonde, Kahan | Probar solvers (κ = 10³ a 10¹⁵) |
| **Sistemas de EDO** | Lorenz, Van der Pol, Robertson problema rígido | Validar paso adaptativo |
| **Integración** | Oscilatoria, singular (1/√x), discontinua | Probar robustez de cuadratura |
| **Valores Propios** | Simétrico, no simétrico, raíces repetidas | Verificar exactitud de descomposición |
| **Funciones** | sin, exp, polinomios con derivadas conocidas | Bancos de precisión |
| **Curvas Paramétricas** | Helix, Lemniscata, Espiral de Toro | Geometría diferencial (curvatura analítica) |

### Bancos de Precisión

**Derivación Numérica** — Error vs derivada analítica de sin(x) en x = 1.0:

| Método |Stencil | Error |
|--------|---------|-------|
| `NDer1` | 2 puntos hacia adelante | ~10⁻⁸ |
| `NDer2` | 3 puntos central | ~10⁻¹¹ |
| `NDer4` | 5 puntos | ~10⁻¹³ |
| `NDer6` | 7 puntos | ~10⁻¹⁴ |
| `NDer8` | 9 puntos | ~10⁻¹⁵ (ε de máquina!) |

**Conservación de Energía de Solvers de EDO** — Órbita de Kepler después de 1000 periodos:

| Solver | Deriva Relativa de Energía |
|--------|----------------------|
| RK4 (paso fijo) | ~10⁻⁶ |
| RKF45 (adaptativo) | ~10⁻¹² |
| Dormand-Prince 8(5,3) | ~10⁻¹⁴ |

**Problemas Mal Condicionados** — MML maneja casos extremos:

```cpp
// Matriz de Hilbert — notoriamente mal condicionada (κ ≈ 10¹⁵ para n=10)
Matrix<Real> H = TestBeds::hilbert_10x10();   // número de condición ~10¹³
LUSolver<Real> solver(H);
Vector<Real> x = solver.Solve(b);
// Logra residuo ||Ax - b|| / ||b|| < 10⁻³ a pesar de la extrema mala condición

// EDO rígida — Cinética química de Robertson (λ ratio = 10⁸)
TestBeds::RobertsonStiffODE stiff;  // λ = {-0.04, -3×10⁴, -10⁸}
// Los solvers implícitos manejan esto; los explícitos necesitarían dt < 10⁻⁸
```

### Ejecutando el Suite de Pruebas

```bash
# Compilar y ejecutar todas las pruebas
cd build && ctest -j8 --output-on-failure

# Ejecutar categoría de prueba específica
ctest -R "integration" -V

# Salida: 4540 pruebas pasan en ~30 segundos
```

📚 **Informes de Análisis de Precisión:** [Resumen](docs/testing_precision/README.md) • [Derivación](docs/testing_precision/DERIVATION_ANALYSIS.md) • [Integración](docs/testing_precision/INTEGRATION_ANALYSIS.md) • [Solvers de EDO](docs/testing_precision/ODE_SOLVER_ANALYSIS.md)

📁 **Documentación de Bancos de Pruebas:** [Funciones](docs/testbeds/Functions_testbed.md) • [Sistemas de EDO](docs/testbeds/ODESystems_testbed.md) • [Sistemas Lineales](docs/testbeds/LinAlgSystems_testbed.md) • [Curvas y Superficies](docs/testbeds/ParametricCurvesSurfaces_testbed.md)

---

## 📚 Documentación

| Recurso | Descripción |
|----------|-------------|
| [Tipos Base](docs/base/README_Base.md) | Vectores, Matrices, Tensores, Funciones |
| [Operaciones Core](docs/core/README_Core.md) | Derivación, Integración, Operaciones de Campos |
| [Algoritmos](docs/algorithms/README_Algorithms.md) | Solvers de EDO, Búsqueda de Raíces, Interpolación |
| [Sistemas](docs/systems/) | Sistemas Dinámicos, Espacio de Fases, Análisis de Estabilidad |
| [Herramientas](docs/tools/) | Visualización, Serialización, Salida en Consola |
| [Ejemplos](examples/) | Ejemplos completos funcionales |
| [Suites de Pruebas](docs/testbeds/) | Simulaciones de física y validación |
| [📖 Referencias de Libros](references/book_references.md) | Libros fundamentales que informaron MML |
| [📄 Referencias de Artículos](references/paperes_references.md) | Artículos académicos detrás de los algoritmos |

---

## 🛠️ Compilación y Pruebas

```bash
# Configurar y compilar
cmake -B build
cmake --build build

# Ejecutar pruebas
cd build && ctest -j8 --output-on-failure

# Compilar ejemplos
cmake --build build --target examples
```

---

## 📜 Licencia

MML se distribuye bajo la **Licencia MIT** — libre para uso personal, académico y comercial.

Ver [LICENSE.md](LICENSE.md) para detalles.

---

## ☕ Apoyar MML

Si MML te ha sido útil, considera apoyar su desarrollo continuo:

<p align="center">
  <a href="https://github.com/sponsors/zvanjak">
    <img src="https://img.shields.io/badge/Sponsor-❤️-pink?logo=github&style=for-the-badge" alt="GitHub Sponsors">
  </a>

</p>

¡Tu apoyo ayuda a mantener y mejorar MML! 🚀

---

<div align="center">

**Hecho con ❤️ para la comunidad de cómputo científico en C++**

[⬆ Volver al Inicio](#-mml--biblioteca-matemática-minimalista)

</div>

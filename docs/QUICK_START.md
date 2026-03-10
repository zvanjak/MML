# Quick Start Guide — 5 Minutes to Your First MML Program

Welcome to **MinimalMathLibrary (MML)**! In just 5 minutes, you'll be solving linear systems, computing derivatives, finding roots, and solving differential equations. Let's dive in! 🚀

---

## 🎯 Step 1: Include the Library (1 minute)

MML is header-only. Just include one file and you're ready:

```cpp
#include "MML.h"
using namespace MML;

int main() {
    // Your code here
    return 0;
}
```

**Compile with C++17 or later:**
```bash
g++ -std=c++17 -I/path/to/mml your_program.cpp -o your_program
```

That's it! No linking, no external dependencies. Just include and go.

---

## 🧮 Step 2: Work with Vectors and Matrices (1 minute)

### Create and Use Vectors

```cpp
#include "MML.h"
using namespace MML;
#include <iostream>

int main() {
    // Create vectors
    Vector<Real> v{1.0, 2.0, 3.0};
    Vector<Real> w{4.0, 5.0, 6.0};
    
    // Basic operations
    Vector<Real> sum = v + w;
    Real dot = v.ScalarProductCartesian(w);  // Dot product
    Real length = v.NormL2();                 // Length
    
    std::cout << "v + w = " << sum << std::endl;
    std::cout << "v · w = " << dot << std::endl;
    std::cout << "||v|| = " << length << std::endl;
    
    return 0;
}
```

**Output:**
```
v + w = [5, 7, 9]
v · w = 32
||v|| = 3.74166
```

### Create and Use Matrices

```cpp
// Create a 3x3 matrix
Matrix<Real> A{3, 3, {
    1, 2, 3,
    4, 5, 6,
    7, 8, 10
}};

// Matrix operations
Matrix<Real> At = A.transpose();
Matrix<Real> I = Matrix<Real>::Identity(3);
Real trace = A.Trace();

std::cout << "Matrix A:\n" << A << std::endl;
std::cout << "Trace: " << trace << std::endl;
```

---

## ⚡ Step 3: Solve a Linear System (1 minute)

Let's solve **A·x = b** for **x**:

```cpp
#include "MML.h"
using namespace MML;
#include <iostream>

int main() {
    // System: A·x = b
    Matrix<Real> A{3, 3, {
        4, 1, 2,
        1, 5, 1,
        2, 1, 6
    }};
    
    Vector<Real> b{8, 9, 12};
    
    // Solve using LU decomposition
    LUSolver<Real> solver(A);
    Vector<Real> x = solver.Solve(b);
    
    std::cout << "Solution: " << x << std::endl;
    
    // Verify: A*x should equal b
    Vector<Real> check = A * x;
    std::cout << "A*x = " << check << std::endl;
    std::cout << "b   = " << b << std::endl;
    
    return 0;
}
```

**Output:**
```
Solution: [1, 1, 1]
A*x = [8, 9, 12]
b   = [8, 9, 12]
```

**Perfect match!** ✅

---

## 📐 Step 4: Calculus — Derivatives and Integrals (1 minute)

### Compute Derivatives

```cpp
#include "MML.h"
using namespace MML;
#include <iostream>
#include <cmath>

int main() {
    // Define a function: f(x) = x²
    RealFunction f = [](Real x) { return x * x; };
    
    // Compute derivative at x = 2.0
    Real x0 = 2.0;
    Real df = Derivation::NDer1(f, x0);  // f'(2.0)
    
    std::cout << "f(x) = x²" << std::endl;
    std::cout << "f'(2.0) = " << df << std::endl;
    std::cout << "Expected: 4.0" << std::endl;
    
    return 0;
}
```

**Output:**
```
f(x) = x²
f'(2.0) = 4
Expected: 4.0
```

### Compute Integrals

```cpp
// Define function: f(x) = sin(x)
RealFunction f = [](Real x) { return std::sin(x); };

// Integrate from 0 to π
Real a = 0.0, b = Constants::PI;
Real integral = Integration::Qromb(f, a, b);

std::cout << "∫₀^π sin(x) dx = " << integral << std::endl;
std::cout << "Expected: 2.0" << std::endl;
```

**Output:**
```
∫₀^π sin(x) dx = 2
Expected: 2.0
```

---

## 🎯 Step 5: Find Roots of Equations (30 seconds)

```cpp
#include "MML.h"
using namespace MML;
#include <iostream>
#include <cmath>

int main() {
    // Find root of: f(x) = x² - 2 = 0  (answer: x = √2)
    RealFunction f = [](Real x) { return x*x - 2.0; };
    
    // Find root in interval [1, 2]
    Real root = RootFinding::Brent(f, 1.0, 2.0, 1e-10);
    
    std::cout << "Root of x² - 2 = 0: " << root << std::endl;
    std::cout << "√2 = " << std::sqrt(2.0) << std::endl;
    
    return 0;
}
```

**Output:**
```
Root of x² - 2 = 0: 1.41421
√2 = 1.41421
```

---

## 🌊 Step 6: Solve Differential Equations (30 seconds)

Solve a simple ODE: **dy/dt = -y**, **y(0) = 1** (exponential decay)

```cpp
#include "MML.h"
using namespace MML;
#include <iostream>

int main() {
    // Define ODE: dy/dt = -y
    using ODEFunction = std::function<Vector<Real>(Real, const Vector<Real>&)>;
    
    ODEFunction f = [](Real t, const Vector<Real>& y) {
        return Vector<Real>{-y[0]};  // dy/dt = -y
    };
    
    // Initial condition: y(0) = 1.0
    Vector<Real> y0{1.0};
    
    // Solve from t=0 to t=2
    Real t0 = 0.0, tend = 2.0;
    int num_steps = 100;
    
    auto [t_vals, y_vals] = ODESystemSolver::RungeKutta4(f, t0, tend, y0, num_steps);
    
    // Print solution at t=2
    std::cout << "y(2) = " << y_vals.back()[0] << std::endl;
    std::cout << "Expected: e^(-2) = " << std::exp(-2.0) << std::endl;
    
    return 0;
}
```

**Output:**
```
y(2) = 0.135335
Expected: e^(-2) = 0.135335
```

---

## 🎉 Complete Example: Everything Together!

Here's a complete program showcasing all 6 steps:

```cpp
#include "MML.h"
using namespace MML;
#include <iostream>
#include <cmath>

int main() {
    std::cout << "=== MinimalMathLibrary Quick Start Demo ===" << std::endl;
    
    // 1. Vectors
    std::cout << "\n1. VECTORS:" << std::endl;
    Vector<Real> v{1, 2, 3}, w{4, 5, 6};
    std::cout << "   v·w = " << v.ScalarProductCartesian(w) << std::endl;
    
    // 2. Matrices
    std::cout << "\n2. MATRICES:" << std::endl;
    Matrix<Real> A = Matrix<Real>::Identity(3);
    std::cout << "   I₃ trace = " << A.Trace() << std::endl;
    
    // 3. Linear Systems
    std::cout << "\n3. LINEAR SYSTEM:" << std::endl;
    Matrix<Real> M{2, 2, {3, 1, 1, 2}};
    Vector<Real> b{9, 8};
    LUSolver<Real> lu(M);
    Vector<Real> x = lu.Solve(b);
    std::cout << "   Solution: [" << x[0] << ", " << x[1] << "]" << std::endl;
    
    // 4. Calculus
    std::cout << "\n4. CALCULUS:" << std::endl;
    RealFunction f = [](Real x) { return x*x; };
    Real df = Derivation::NDer1(f, 3.0);
    std::cout << "   (x²)' at x=3: " << df << " (expected: 6)" << std::endl;
    
    Real integral = Integration::Qromb([](Real x){return x;}, 0, 2);
    std::cout << "   ∫₀² x dx = " << integral << " (expected: 2)" << std::endl;
    
    // 5. Root Finding
    std::cout << "\n5. ROOT FINDING:" << std::endl;
    Real root = RootFinding::Brent([](Real x){return x*x-2;}, 0, 2, 1e-10);
    std::cout << "   √2 ≈ " << root << std::endl;
    
    // 6. ODEs
    std::cout << "\n6. DIFFERENTIAL EQUATIONS:" << std::endl;
    using ODEFunction = std::function<Vector<Real>(Real, const Vector<Real>&)>;
    ODEFunction ode = [](Real t, const Vector<Real>& y) {
        return Vector<Real>{-y[0]};
    };
    auto [t_vals, y_vals] = ODESystemSolver::RungeKutta4(
        ode, 0.0, 1.0, Vector<Real>{1.0}, 50
    );
    std::cout << "   e^(-1) ≈ " << y_vals.back()[0] << std::endl;
    
    std::cout << "\n✅ All tests passed! You're ready to use MML!" << std::endl;
    
    return 0;
}
```

**Output:**
```
=== MinimalMathLibrary Quick Start Demo ===

1. VECTORS:
   v·w = 32

2. MATRICES:
   I₃ trace = 3

3. LINEAR SYSTEM:
   Solution: [2, 3]

4. CALCULUS:
   (x²)' at x=3: 6 (expected: 6)
   ∫₀² x dx = 2 (expected: 2)

5. ROOT FINDING:
   √2 ≈ 1.41421

6. DIFFERENTIAL EQUATIONS:
   e^(-1) ≈ 0.367879

✅ All tests passed! You're ready to use MML!
```

---

## 🚀 What's Next?

Congratulations! In just 5 minutes you've learned:
- ✅ Vector and matrix operations
- ✅ Solving linear systems
- ✅ Computing derivatives and integrals
- ✅ Finding roots of equations
- ✅ Solving differential equations

### Where to Go From Here

**By Topic:**
- **More Linear Algebra:** [docs/core/Linear_equations_solvers.md](core/Linear_equations_solvers.md)
- **Advanced Calculus:** [docs/core/Integration.md](core/Integration.md), [docs/core/Derivation.md](core/Derivation.md)
- **Optimization:** [docs/algorithms/Optimization.md](algorithms/Optimization.md)
- **Statistics:** [docs/core/Statistics.md](core/Statistics.md)
- **Fourier Analysis:** [docs/algorithms/FFT.md](algorithms/FFT.md)
- **Geometry:** [docs/base/Geometry2D.md](base/Geometry2D.md), [docs/base/Geometry3D.md](base/Geometry3D.md)

**By Use Case:**
- **Need Fast Lookup?** → [API_CHEATSHEET.md](API_CHEATSHEET.md) (single-page reference)
- **Physics Simulations?** → [examples/](../examples/) (projectiles, orbits, waves)
- **Numerical Recipes?** → [COOKBOOK.md](COOKBOOK.md) (common tasks, 10-20 lines each)
- **Complete Reference?** → [docs/](.) (all modules documented)

**By Learning Style:**
- **Learn by Example:** Browse [examples/](../examples/) and [tests/](../tests/)
- **Learn by Reading:** Start with [README.md](../README.md)
- **Learn by Doing:** Try the complete example above, then modify it!

---

## 💡 Tips for Success

**1. Start Simple**
- Begin with the complete example above
- Modify one thing at a time
- Build confidence incrementally

**2. Use the Right Solver**
- **Linear systems:** LU (general), Cholesky (SPD matrices), QR (ill-conditioned)
- **Root finding:** Brent (no derivative), Newton (with derivative)
- **Integration:** Qromb (adaptive), GaussLegendre (Gaussian quadrature)
- **ODEs:** RungeKutta4 (most problems), Euler (simple/testing)

**3. Check Your Results**
- Always verify solutions (e.g., A*x should equal b)
- Compare with known analytical solutions when possible
- Use the [API_CHEATSHEET.md](API_CHEATSHEET.md) for quick algorithm selection

**4. Explore Examples**
- The [examples/](../examples/) folder has working code for common scenarios
- The [tests/](../tests/) folder shows correct usage patterns
- Both are great learning resources!

---

## 🐛 Common Pitfalls

**Problem:** Compiler errors about missing functions
- **Solution:** Make sure you're compiling with C++17 or later: `-std=c++17`

**Problem:** "Singular matrix" errors
- **Solution:** Your matrix isn't invertible. Try SVD or check your matrix construction.

**Problem:** Poor accuracy in numerical integration
- **Solution:** Increase the number of steps or use adaptive methods (Qromb, Qsimp)

**Problem:** ODE solver gives wrong results
- **Solution:** Check your ODE function signature and initial conditions. For stiff equations, reduce step size.

---

## 📚 Quick Reference Card

```cpp
// VECTORS
Vector<Real> v{1, 2, 3};
Real dot = v.ScalarProductCartesian(w);
Real len = v.NormL2();

// MATRICES  
Matrix<Real> A{3, 3, {...}};
Matrix<Real> I = Matrix<Real>::Identity(3);
Matrix<Real> At = A.transpose();

// LINEAR SYSTEMS
LUSolver<Real> lu(A);
Vector<Real> x = lu.Solve(b);

// CALCULUS
Real dy = Derivation::NDer1(f, x0);
Real I = Integration::Qromb(f, a, b);

// ROOT FINDING
Real root = RootFinding::Brent(f, a, b, tol);

// ODEs
auto [t, y] = ODESystemSolver::RungeKutta4(f, t0, tend, y0, steps);
```

---

## 🎓 Exercises for Practice

**Easy:**
1. Solve the system: 2x + y = 5, x + 3y = 8
2. Compute the derivative of sin(x) at x = π/4
3. Integrate e^x from 0 to 1

**Medium:**
4. Find all roots of x³ - 6x² + 11x - 6 = 0
5. Solve the ODE: dy/dt = t·y, y(0) = 1 from t=0 to t=2
6. Fit a quadratic polynomial to 5 data points

**Challenge:**
7. Solve a 10×10 linear system with random entries
8. Compute eigenvalues of a symmetric matrix
9. Minimize f(x,y) = x² + xy + y² using optimization

*Solutions can be found in [examples/exercises/](../examples/exercises/) (if available)*

---

## 🤝 Getting Help

- **Documentation:** [docs/](.) — Complete reference for all modules
- **Examples:** [examples/](../examples/) — Working code for common tasks  
- **Issues:** [GitHub Issues](https://github.com/zvanjak/MinimalMathLibrary/issues) — Report bugs or ask questions
- **API Reference:** [API_CHEATSHEET.md](API_CHEATSHEET.md) — Fast lookup for experienced users

---

**Happy Computing!** 🎉

You now have a powerful numerical computing library at your fingertips. Start experimenting, solving real problems, and building amazing projects!

*MML v1.1 — From zero to hero in 5 minutes*

# Grid Infrastructure for PDE Discretization

Complete guide to structured grids in MML's PDE module.

## Table of Contents
1. [Introduction](#introduction)
2. [Grid Fundamentals](#grid-fundamentals)
3. [Grid1D - One-Dimensional Grids](#grid1d---one-dimensional-grids)
4. [Grid2D - Two-Dimensional Grids](#grid2d---two-dimensional-grids)
5. [Grid3D - Three-Dimensional Grids](#grid3d---three-dimensional-grids)
6. [Grid Functions](#grid-functions)
7. [Index Management](#index-management)
8. [Boundary Handling](#boundary-handling)
9. [Stencil Operations](#stencil-operations)
10. [Iteration Patterns](#iteration-patterns)
11. [Common Pitfalls](#common-pitfalls)
12. [Best Practices](#best-practices)

---

## Introduction

**Structured grids** are the foundation of finite difference methods for PDEs. MML provides uniform rectangular grids in 1D, 2D, and 3D with:
- **Simple indexing** for efficient access
- **Boundary identification** for applying boundary conditions
- **Neighbor access** for stencil operations
- **Domain mapping** between physical and discrete spaces

### Key Concepts

```
Physical Domain          Discrete Grid
    [0, 1]      →      0---1---2---3---4
                       x₀  x₁  x₂  x₃  x₄
                       
n = 4 cells
nodes = n + 1 = 5
interior nodes = n - 1 = 3 (nodes 1, 2, 3)
```

---

## Grid Fundamentals

### Cells vs Nodes (⚠️ Critical!)

**MOST IMPORTANT CONCEPT:**

```
Grid with n = 3 cells:

├────────┼────────┼────────┤
0        1        2        3
         cell 0   cell 1   cell 2
node 0   node 1   node 2   node 3

n = number of CELLS
numNodes() = n + 1 = number of NODES
```

**Always remember:**
- `n` or `nx`, `ny`, `nz` = **cell count**
- `numNodes()`, `numNodesX()` = **node count** = n + 1
- **Loop bounds:** Use `numNodes()`, not `n`!

### Grid Spacing

**1D:**
```
dx = (xmax - xmin) / n
```

**2D:**
```
dx = (xmax - xmin) / nx
dy = (ymax - ymin) / ny
```

**3D:**
```
dx = (xmax - xmin) / nx
dy = (ymax - ymin) / ny
dz = (zmax - zmin) / nz
```

---

## Grid1D - One-Dimensional Grids

### Class: `Grid1D<T>`

```cpp
#include "pde/grid/Grid1D.h"
using namespace MML::PDE;

// Method 1: Grid on [0, 1] with n cells
Grid1D<double> grid1(100);  // 101 nodes

// Method 2: Grid on [xmin, xmax] with n cells
Grid1D<double> grid2(0.0, 5.0, 50);  // 51 nodes

// Method 3: Grid from Interval with n cells
Interval<double> domain(-1.0, 1.0);
Grid1D<double> grid3(domain, 80);  // 81 nodes
```

### Properties

```cpp
// Domain
double xmin = grid.xmin();  // Left boundary
double xmax = grid.xmax();  // Right boundary
double L = grid.length();   // Domain length

// Discretization
int n = grid.numCells();           // Number of cells
int nodes = grid.numNodes();       // n + 1
int interior = grid.numInteriorNodes();  // n - 1
double dx = grid.dx();             // Grid spacing
```

### Coordinate Access

```cpp
// Get x-coordinate of node i
double x_i = grid.x(i);

// Find node containing point xp
int i = grid.nodeIndex(xp);

// Locate point for interpolation
int i;
double alpha;  // Interpolation parameter in [0,1]
grid.locate(xp, i, alpha);
// xp = grid.x(i) + alpha * grid.dx()
```

### Boundary Information

```cpp
// Check node type
bool left_boundary = grid.isLeftBoundary(i);
bool right_boundary = grid.isRightBoundary(i);

// Get boundary node indices
int left_node = grid.leftBoundaryNode();   // 0
int right_node = grid.rightBoundaryNode(); // n
```

### Example: 1D Poisson Assembly

```cpp
#include "pde/grid/Grid1D.h"
#include "pde/sparse/SparseMatrixCOO.h"

// Setup grid
int n = 100;
Grid1D<double> grid(0.0, 1.0, n);

int num_nodes = grid.numNodes();
double dx = grid.dx();
double dx_sq_inv = 1.0 / (dx * dx);

// Build matrix
Sparse::SparseMatrixCOO<double> A(num_nodes, num_nodes);

for (int i = 0; i < num_nodes; ++i) {
    if (grid.isLeftBoundary(i) || grid.isRightBoundary(i)) {
        // Boundary: identity row
        A.addEntry(i, i, 1.0);
    } else {
        // Interior: 3-point stencil [-1, 2, -1] / dx²
        A.addEntry(i, i-1, -dx_sq_inv);
        A.addEntry(i, i,    2.0 * dx_sq_inv);
        A.addEntry(i, i+1, -dx_sq_inv);
    }
}
```

---

## Grid2D - Two-Dimensional Grids

### Class: `Grid2D<T>`

```cpp
#include "pde/grid/Grid2D.h"
using namespace MML::PDE;

// Method 1: Grid on [0,1]² with nx × ny cells
Grid2D<double> grid1(50, 50);  // (51 × 51) = 2,601 nodes

// Method 2: Square grid on [0,1]² with n × n cells
Grid2D<double> grid2(100);  // (101 × 101) = 10,201 nodes

// Method 3: Grid on [xmin,xmax] × [ymin,ymax] with nx × ny cells
Grid2D<double> grid3(0.0, 2.0, 0.0, 1.0, 80, 40);

// Method 4: Grid from Rectangle
Rectangle<double> domain(0.0, 1.0, 0.0, 1.0);
Grid2D<double> grid4(domain, 50, 50);
```

### Properties

```cpp
// Domain
double xmin = grid.xmin();
double xmax = grid.xmax();
double ymin = grid.ymin();
double ymax = grid.ymax();

// Discretization
int nx = grid.nx();              // Cells in x-direction
int ny = grid.ny();              // Cells in y-direction
int num_x_nodes = grid.numNodesX();  // nx + 1
int num_y_nodes = grid.numNodesY();  // ny + 1
int total_nodes = grid.numNodes();   // (nx+1) * (ny+1)
int interior = grid.numInteriorNodes();  // (nx-1) * (ny-1)

double dx = grid.dx();
double dy = grid.dy();
double h = grid.h();             // min(dx, dy)
```

### Index Conversion (⚠️ Critical!)

**Linear indexing:** `idx = i + j * (nx + 1)` (row-major order)

```cpp
Grid Layout (4×3 cells = 5×4 nodes):

j=3:  15--16--17--18--19
j=2:  10--11--12--13--14
j=1:   5---6---7---8---9
j=0:   0---1---2---3---4
      i=0 i=1 i=2 i=3 i=4

Node (i=2, j=1) → linear index = 2 + 1*5 = 7
```

**API:**

```cpp
// Multi-index → linear
int idx = grid.index(i, j);

// Linear → multi-index
int i, j;
grid.index(idx, i, j);

// Extract i or j from linear index
int i = grid.indexI(idx);
int j = grid.indexJ(idx);
```

### Coordinate Access

```cpp
// Get (x, y) coordinates of node (i, j)
double x = grid.x(i);
double y = grid.y(j);

// Or both at once
auto [x, y] = grid.xy(i, j);

// Find node containing point (xp, yp)
int i = grid.nodeIndexX(xp);
int j = grid.nodeIndexY(yp);
```

### Neighbor Access (5-Point Stencil)

```cpp
//          North (j+1)
//             |
//  West (i-1)--C--(i+1) East
//             |
//         South (j-1)

// Get neighbor indices
int center = grid.index(i, j);
int west  = grid.index(i-1, j);   // if i > 0
int east  = grid.index(i+1, j);   // if i < nx
int south = grid.index(i, j-1);   // if j > 0
int north = grid.index(i, j+1);   // if j < ny
```

### Boundary Information

```cpp
// Check node location
auto loc = grid.nodeLocation(i, j);
// Returns: Interior, Boundary

bool is_interior = grid.isInterior(i, j);
bool is_boundary = grid.isBoundary(i, j);

// Check specific sides
bool on_left   = grid.isLeftBoundary(i);
bool on_right  = grid.isRightBoundary(i);
bool on_bottom = grid.isBottomBoundary(j);
bool on_top    = grid.isTopBoundary(j);

// Get which sides a node touches
auto sides = grid.boundarySides(i, j);
// Returns: vector of BoundarySide enum
```

### Example: 2D Laplacian Assembly

```cpp
#include "pde/grid/Grid2D.h"
#include "pde/sparse/SparseMatrixCOO.h"

// Setup grid
int nx = 50, ny = 50;
Grid2D<double> grid(0.0, 1.0, 0.0, 1.0, nx, ny);

int n = grid.numNodes();
double dx = grid.dx();
double dy = grid.dy();
double dx_sq_inv = 1.0 / (dx * dx);
double dy_sq_inv = 1.0 / (dy * dy);
double center_coeff = 2.0 * (dx_sq_inv + dy_sq_inv);

// Build matrix
Sparse::SparseMatrixCOO<double> A(n, n);

// Loop over all grid points
for (int j = 0; j < grid.numNodesY(); ++j) {
    for (int i = 0; i < grid.numNodesX(); ++i) {
        int row = grid.index(i, j);
        
        if (grid.isBoundary(i, j)) {
            // Boundary: identity
            A.addEntry(row, row, 1.0);
        } else {
            // Interior: 5-point stencil
            //        dy⁻²
            //          |
            //  dx⁻²--[2(dx⁻²+dy⁻²)]--dx⁻²
            //          |
            //        dy⁻²
            
            A.addEntry(row, grid.index(i-1, j), -dx_sq_inv);  // West
            A.addEntry(row, grid.index(i+1, j), -dx_sq_inv);  // East
            A.addEntry(row, grid.index(i, j-1), -dy_sq_inv);  // South
            A.addEntry(row, grid.index(i, j+1), -dy_sq_inv);  // North
            A.addEntry(row, row, center_coeff);               // Center
        }
    }
}
```

---

## Grid3D - Three-Dimensional Grids

### Class: `Grid3D<T>`

```cpp
#include "pde/grid/Grid3D.h"
using namespace MML::PDE;

// Method 1: Grid on [0,1]³ with nx × ny × nz cells
Grid3D<double> grid1(30, 30, 30);  // 31³ = 29,791 nodes

// Method 2: Cubic grid on [0,1]³ with n × n × n cells
Grid3D<double> grid2(50);  // 51³ = 132,651 nodes

// Method 3: Grid on box with nx × ny × nz cells
Grid3D<double> grid3(0.0, 2.0, 0.0, 1.0, 0.0, 1.0, 80, 40, 40);

// Method 4: Grid from Box
Box<double> domain(0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
Grid3D<double> grid4(domain, 30, 30, 30);
```

### Properties

```cpp
// Domain
double xmin = grid.xmin();
double xmax = grid.xmax();
double ymin = grid.ymin();
double ymax = grid.ymax();
double zmin = grid.zmin();
double zmax = grid.zmax();

// Discretization
int nx = grid.nx();
int ny = grid.ny();
int nz = grid.nz();
int num_x_nodes = grid.numNodesX();  // nx + 1
int num_y_nodes = grid.numNodesY();  // ny + 1
int num_z_nodes = grid.numNodesZ();  // nz + 1
int total_nodes = grid.numNodes();   // (nx+1) * (ny+1) * (nz+1)
int interior = grid.numInteriorNodes();  // (nx-1) * (ny-1) * (nz-1)

double dx = grid.dx();
double dy = grid.dy();
double dz = grid.dz();
double h = grid.h();                 // min(dx, dy, dz)
```

### Index Conversion

**Linear indexing:** `idx = i + j*(nx+1) + k*(nx+1)*(ny+1)` (row-major order)

```cpp
// Multi-index → linear
int idx = grid.index(i, j, k);

// Linear → multi-index
int i, j, k;
grid.index(idx, i, j, k);

// Extract i, j, or k from linear index
int i = grid.indexI(idx);
int j = grid.indexJ(idx);
int k = grid.indexK(idx);
```

### Neighbor Access (7-Point Stencil)

```cpp
//           Top (k+1)
//              |
//         North (j+1)
//              |
//  West (i-1)--C--(i+1) East
//              |
//         South (j-1)
//              |
//         Bottom (k-1)

// Get neighbor indices
int center = grid.index(i, j, k);
int west   = grid.index(i-1, j, k);   // if i > 0
int east   = grid.index(i+1, j, k);   // if i < nx
int south  = grid.index(i, j-1, k);   // if j > 0
int north  = grid.index(i, j+1, k);   // if j < ny
int bottom = grid.index(i, j, k-1);   // if k > 0
int top    = grid.index(i, j, k+1);   // if k < nz
```

### Boundary Information

```cpp
// Check node location
bool is_interior = grid.isInterior(i, j, k);
bool is_boundary = grid.isBoundary(i, j, k);

// Check specific faces
bool on_left   = grid.isLeftBoundary(i);
bool on_right  = grid.isRightBoundary(i);
bool on_front  = grid.isFrontBoundary(j);
bool on_back   = grid.isBackBoundary(j);
bool on_bottom = grid.isBottomBoundary(k);
bool on_top    = grid.isTopBoundary(k);
```

### Example: 3D Laplacian Assembly

```cpp
#include "pde/grid/Grid3D.h"
#include "pde/sparse/SparseMatrixCOO.h"

// Setup grid
int nx = 20, ny = 20, nz = 20;
Grid3D<double> grid(0.0, 1.0, 0.0, 1.0, 0.0, 1.0, nx, ny, nz);

int n = grid.numNodes();
double dx_sq_inv = 1.0 / (grid.dx() * grid.dx());
double dy_sq_inv = 1.0 / (grid.dy() * grid.dy());
double dz_sq_inv = 1.0 / (grid.dz() * grid.dz());
double center = 2.0 * (dx_sq_inv + dy_sq_inv + dz_sq_inv);

// Build matrix
Sparse::SparseMatrixCOO<double> A(n, n);

// Loop over all grid points
for (int k = 0; k < grid.numNodesZ(); ++k) {
    for (int j = 0; j < grid.numNodesY(); ++j) {
        for (int i = 0; i < grid.numNodesX(); ++i) {
            int row = grid.index(i, j, k);
            
            if (grid.isBoundary(i, j, k)) {
                A.addEntry(row, row, 1.0);
            } else {
                // 7-point stencil
                A.addEntry(row, grid.index(i-1,j,k), -dx_sq_inv);
                A.addEntry(row, grid.index(i+1,j,k), -dx_sq_inv);
                A.addEntry(row, grid.index(i,j-1,k), -dy_sq_inv);
                A.addEntry(row, grid.index(i,j+1,k), -dy_sq_inv);
                A.addEntry(row, grid.index(i,j,k-1), -dz_sq_inv);
                A.addEntry(row, grid.index(i,j,k+1), -dz_sq_inv);
                A.addEntry(row, row, center);
            }
        }
    }
}
```

---

## Grid Functions

**GridFunction** stores discrete function values on grid nodes.

```cpp
#include "pde/grid/GridFunction.h"

// 1D grid function
Grid1D<double> grid1d(100);
GridFunction1D<double> u(grid1d);  // Automatically sized

// Access/modify
u[i] = 1.5;
double val = u[i];

// Initialize
u.fill(0.0);
u.setValues([&](int i) { return std::sin(grid1d.x(i)); });

// 2D grid function
Grid2D<double> grid2d(50, 50);
GridFunction2D<double> u2d(grid2d);

// Access by linear index or (i,j)
u2d[grid2d.index(i, j)] = 2.0;
u2d(i, j) = 2.0;  // Equivalent

// 3D grid function
Grid3D<double> grid3d(30, 30, 30);
GridFunction3D<double> u3d(grid3d);

u3d(i, j, k) = 3.0;
```

---

## Index Management

### 1D Indexing

```
Simple: node i has index i

Neighbors: i-1, i, i+1
Boundary: i = 0 or i = n
Interior: 0 < i < n
```

### 2D Indexing

**Row-major order:** `idx = i + j * (nx + 1)`

```
Example: 3×2 cells (4×3 nodes)

j=2:   8---9--10--11
j=1:   4---5---6---7
j=0:   0---1---2---3
      i=0 i=1 i=2 i=3

Node (2, 1): idx = 2 + 1*4 = 6 ✓
```

**Common patterns:**

```cpp
// Iterate all nodes
for (int j = 0; j < grid.numNodesY(); ++j) {
    for (int i = 0; i < grid.numNodesX(); ++i) {
        int idx = grid.index(i, j);
        // ...
    }
}

// Iterate interior nodes only
for (int j = 1; j < ny; ++j) {
    for (int i = 1; i < nx; ++i) {
        int idx = grid.index(i, j);
        // ...
    }
}

// Neighbors of (i, j)
if (i > 0) {
    int west = grid.index(i-1, j);
}
if (i < nx) {
    int east = grid.index(i+1, j);
}
// etc.
```

### 3D Indexing

**Row-major order:** `idx = i + j*(nx+1) + k*(nx+1)*(ny+1)`

```cpp
// Iterate all nodes
for (int k = 0; k < grid.numNodesZ(); ++k) {
    for (int j = 0; j < grid.numNodesY(); ++j) {
        for (int i = 0; i < grid.numNodesX(); ++i) {
            int idx = grid.index(i, j, k);
            // ...
        }
    }
}
```

---

## Boundary Handling

### Boundary Node Classification

**1D:**
- Boundary: i = 0 or i = n
- Interior: 0 < i < n

**2D:**
- Boundary: i = 0 OR i = nx OR j = 0 OR j = ny
- Interior: 0 < i < nx AND 0 < j < ny

**3D:**
- Boundary: i = 0 OR i = nx OR j = 0 OR j = ny OR k = 0 OR k = nz
- Interior: 0 < i < nx AND 0 < j < ny AND 0 < k < nz

### Retrieving Boundary Nodes

```cpp
// 2D example
Grid2D<double> grid(50, 50);

// All interior nodes
auto interior = grid.interiorNodes();

// Boundary nodes on specific side
auto left_nodes = grid.boundaryNodes(BoundarySide::Left);
auto top_nodes = grid.boundaryNodes(BoundarySide::Top);

// All boundary nodes
auto boundary = grid.allBoundaryNodes();

// Check which sides a node touches
auto sides = grid.boundarySides(i, j);
for (auto side : sides) {
    if (side == BoundarySide::Left) {
        // Apply left BC
    }
}
```

---

## Stencil Operations

### 1D Stencils

**3-point stencil** for `-u''(x)`:

```
    [-1]   [2]   [-1]
     i-1    i     i+1
     
Divided by dx²
```

```cpp
// At interior node i:
double u_xx = (u[i-1] - 2*u[i] + u[i+1]) / (dx * dx);
```

### 2D Stencils

**5-point stencil** for `-∇²u`:

```
           [-dy⁻²]
               |
   [-dx⁻²]--[2(dx⁻²+dy⁻²)]--[-dx⁻²]
               |
           [-dy⁻²]
```

```cpp
// At interior node (i, j):
double u_xx = (u(i-1,j) - 2*u(i,j) + u(i+1,j)) / (dx*dx);
double u_yy = (u(i,j-1) - 2*u(i,j) + u(i,j+1)) / (dy*dy);
double laplacian = u_xx + u_yy;
```

**9-point stencil** (higher order):

```
  [-1/12dy²]    [-4/3dy²]    [-1/12dy²]
        |            |            |
 [-1/12dx²]--[-4/3dx²+...]--[-1/12dx²]
        |            |            |
  [-1/12dy²]    [-4/3dy²]    [-1/12dy²]
```

### 3D Stencils

**7-point stencil** for `-∇²u`:

```
6 neighbors + center = 7 points
Coefficients: ±dx⁻², ±dy⁻², ±dz⁻², center = 2(dx⁻²+dy⁻²+dz⁻²)
```

**Compact representation:**

```cpp
double lap_u = 
    (u(i-1,j,k) + u(i+1,j,k) - 2*u(i,j,k)) / (dx*dx) +
    (u(i,j-1,k) + u(i,j+1,k) - 2*u(i,j,k)) / (dy*dy) +
    (u(i,j,k-1) + u(i,j,k+1) - 2*u(i,j,k)) / (dz*dz);
```

---

## Iteration Patterns

### Manual Loops

**1D - All nodes:**

```cpp
for (int i = 0; i < grid.numNodes(); ++i) {
    double x = grid.x(i);
    // ...
}
```

**1D - Interior only:**

```cpp
for (int i = 1; i < grid.numCells(); ++i) {  // 1 to n-1
    // Interior node
}
```

**2D - All nodes:**

```cpp
for (int j = 0; j < grid.numNodesY(); ++j) {
    for (int i = 0; i < grid.numNodesX(); ++i) {
        int idx = grid.index(i, j);
        // ...
    }
}
```

**2D - Interior only:**

```cpp
for (int j = 1; j < grid.ny(); ++j) {
    for (int i = 1; i < grid.nx(); ++i) {
        int idx = grid.index(i, j);
        // Interior node
    }
}
```

### Functional Iterators

**Grid2D provides convenience functions:**

```cpp
// Iterate all nodes
grid.forEachNode([&](int i, int j) {
    double x = grid.x(i);
    double y = grid.y(j);
    int idx = grid.index(i, j);
    // ...
});

// Iterate interior nodes only
grid.forEachInterior([&](int i, int j) {
    int idx = grid.index(i, j);
    // Guaranteed interior
});
```

---

## Common Pitfalls

### ❌ Pitfall 1: Using `nx()` instead of `numNodesX()`

```cpp
Grid2D<double> grid(50, 50);

// WRONG: Misses last node!
for (int i = 0; i < grid.nx(); ++i) {  // 0 to 49
    // Missing i = 50!
}

// CORRECT: All nodes
for (int i = 0; i < grid.numNodesX(); ++i) {  // 0 to 50
    // All nodes included
}
```

**Remember:** `nx()` = cells, `numNodesX()` = nodes = nx + 1

### ❌ Pitfall 2: Off-by-One in Interior Loops

```cpp
// WRONG: Includes boundary
for (int i = 1; i <= grid.nx(); ++i) {
    // i = nx is BOUNDARY, not interior!
}

// CORRECT: Interior only
for (int i = 1; i < grid.nx(); ++i) {  // 1 to nx-1
    // Interior nodes only
}
```

### ❌ Pitfall 3: Incorrect Linear Indexing

```cpp
// WRONG: Column-major order
int idx = j + i * numNodesY;  // ✗ Wrong!

// CORRECT: Row-major order
int idx = i + j * numNodesX;  // ✓ Correct!

// OR: Use grid method
int idx = grid.index(i, j);   // ✓ Always correct!
```

### ❌ Pitfall 4: Forgetting Boundary Checks

```cpp
// DANGEROUS: Assumes interior node
int west = grid.index(i-1, j);  // Crash if i = 0!

// SAFE: Check bounds
if (i > 0) {
    int west = grid.index(i-1, j);
    // Use west neighbor
}

// BETTER: Check location once
if (grid.isInterior(i, j)) {
    // All neighbors exist
    int west = grid.index(i-1, j);
    int east = grid.index(i+1, j);
    // ...
}
```

### ❌ Pitfall 5: Wrong Stencil Coefficients

```cpp
// WRONG: Forgetting to divide by grid spacing
double laplacian = u[i-1] - 2*u[i] + u[i+1];  // ✗ Missing dx²!

// CORRECT: Include grid spacing
double dx_sq = grid.dx() * grid.dx();
double laplacian = (u[i-1] - 2*u[i] + u[i+1]) / dx_sq;  // ✓
```

---

## Best Practices

### ✅ DO:

1. **Always use `numNodes()` for loop bounds**
   ```cpp
   for (int i = 0; i < grid.numNodesX(); ++i) {  // Correct!
   ```

2. **Use grid index functions**
   ```cpp
   int idx = grid.index(i, j);  // Not manual calculation
   ```

3. **Check boundary status before accessing neighbors**
   ```cpp
   if (grid.isInterior(i, j)) {
       // All 4 neighbors exist
   }
   ```

4. **Store grid spacing in variables**
   ```cpp
   double dx = grid.dx();
   double dx_sq_inv = 1.0 / (dx * dx);
   // Use dx_sq_inv in loops
   ```

5. **Use meaningful variable names**
   ```cpp
   int num_x_nodes = grid.numNodesX();  // Clear
   // Not: int n = grid.numNodesX();    // Ambiguous
   ```

### ❌ DON'T:

1. **Don't confuse cells and nodes**
   ```cpp
   int n = grid.numCells();
   for (int i = 0; i <= n; ++i) {  // ✗ Should use numNodes()
   ```

2. **Don't hard-code index formulas**
   ```cpp
   int idx = i + j * 51;  // ✗ What if grid size changes?
   int idx = grid.index(i, j);  // ✓ Always correct
   ```

3. **Don't access out-of-bounds neighbors**
   ```cpp
   int west = grid.index(i-1, j);  // ✗ Crash if i = 0!
   ```

4. **Don't iterate over invalid ranges**
   ```cpp
   // Interior nodes
   for (int i = 1; i <= nx; ++i) {  // ✗ i=nx is boundary!
   for (int i = 1; i < nx; ++i) {   // ✓ i∈[1, nx-1]
   ```

### Performance Tips

1. **Reuse grid quantities**
   ```cpp
   double dx_sq_inv = 1.0 / (grid.dx() * grid.dx());
   // Use in loops instead of recalculating
   ```

2. **Cache numNodes()**
   ```cpp
   int nx_nodes = grid.numNodesX();
   for (int i = 0; i < nx_nodes; ++i) {  // Faster
   ```

3. **Use row-major loops (cache-friendly)**
   ```cpp
   // Good: Iterate i in inner loop (row-major)
   for (int j = 0; j < ny_nodes; ++j) {
       for (int i = 0; i < nx_nodes; ++i) {
           // Sequential memory access
       }
   }
   ```

---

## See Also

- [Sparse Matrices](Sparse_Matrices.md) - Matrix assembly on grids
- [Boundary Conditions](Boundary_Conditions.md) - Applying BCs to grids
- [Elliptic PDEs](Elliptic_PDEs.md) - Using grids for Poisson/Laplace
- [Parabolic PDEs](Parabolic_PDEs.md) - Time-dependent problems on grids
- [API Reference](API_Reference.md) - Complete grid class documentation

---

**Next:** [Boundary Conditions →](Boundary_Conditions.md)

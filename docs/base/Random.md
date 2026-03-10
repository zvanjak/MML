# Random - Random Number Generation

**File**: `mml/base/Random.h`

Thread-safe random number generation utilities using Mersenne Twister (MT19937).

## Table of Contents
- [Overview](#overview)
- [Uniform Real Distribution](#uniform-real-distribution)
- [Uniform Integer Distribution](#uniform-integer-distribution)
- [Random Directions](#random-directions)
- [Thread Safety](#thread-safety)
- [Examples](#examples)

---

## Overview

The `Random` class provides static methods for generating random numbers with:
- ✅ **Thread safety** - Each thread has its own RNG (thread_local)
- ✅ **High quality** - Uses Mersenne Twister (MT19937)
- ✅ **High performance** - No locking, RNG reused across calls
- ✅ **Uniform distributions** - Real and integer ranges
- ✅ **Direction vectors** - Uniform on circle/sphere

**Note**: Results are non-reproducible by design. The RNG is seeded from `std::random_device` on first use per thread.

---

## Uniform Real Distribution

### `UniformReal(min, max)`

Generate a random real number uniformly distributed in `[min, max)`.

```cpp
// Random in [0, 1)
Real x = Random::UniformReal(0.0, 1.0);

// Random in [-5, 5)
Real y = Random::UniformReal(-5.0, 5.0);

// Random angle in [0, 2π)
Real angle = Random::UniformReal(0.0, 2.0 * Constants::PI);
```

**Parameters:**
- `min` - Lower bound (inclusive)
- `max` - Upper bound (exclusive)

**Returns:** Random value in `[min, max)`

---

## Uniform Integer Distribution

### `UniformInt(min, max)`

Generate a random integer uniformly distributed in `[min, max]`.

```cpp
// Dice roll (1-6)
int dice = Random::UniformInt(1, 6);

// Random array index (0 to n-1)
int idx = Random::UniformInt(0, array.size() - 1);

// Random boolean
bool coin = (Random::UniformInt(0, 1) == 1);
```

**Parameters:**
- `min` - Lower bound (inclusive)
- `max` - Upper bound (inclusive)

**Returns:** Random integer in `[min, max]`

---

## Random Directions

### `UniformVecDirection2(vx, vy, abs)`

Generate a random 2D unit direction uniformly distributed on the circle.

```cpp
Real vx, vy;

// Unit direction (|v| = 1)
Random::UniformVecDirection2(vx, vy, 1.0);

// Direction with magnitude 5
Random::UniformVecDirection2(vx, vy, 5.0);
```

**Parameters:**
- `vx` - Output: x-component
- `vy` - Output: y-component  
- `abs` - Desired vector magnitude

**Returns:** The magnitude (same as `abs` parameter)

**Algorithm:** Generates a uniform random angle θ in `[0, 2π)`, then computes:
- `vx = abs * cos(θ)`
- `vy = abs * sin(θ)`

### `UniformVecDirection3(vx, vy, vz, abs)`

Generate a random 3D unit direction uniformly distributed on the sphere.

```cpp
Real vx, vy, vz;

// Unit direction (|v| = 1)
Random::UniformVecDirection3(vx, vy, vz, 1.0);

// Direction with magnitude 10
Random::UniformVecDirection3(vx, vy, vz, 10.0);
```

**Parameters:**
- `vx` - Output: x-component
- `vy` - Output: y-component
- `vz` - Output: z-component
- `abs` - Desired vector magnitude

**Returns:** The magnitude (same as `abs` parameter)

**Algorithm:** Uses the correct uniform sphere sampling method:
1. Sample `u` uniformly in `[-1, 1]` (this is `cos(φ)`)
2. Sample `θ` uniformly in `[0, 2π)`
3. Compute `sin(φ) = sqrt(1 - u²)`
4. Return: `(abs * sin(φ) * cos(θ), abs * sin(φ) * sin(θ), abs * u)`

> **Note:** This algorithm ensures truly uniform distribution on the sphere surface, unlike naive methods that would cluster points near the poles.

---

## Thread Safety

The `Random` class uses a `thread_local` Mersenne Twister generator:

```cpp
inline static thread_local std::mt19937 gen{std::random_device{}()};
```

**Benefits:**
- **Thread safe** - No data races between threads
- **No locking** - Each thread has independent RNG
- **Fast** - RNG is created once per thread, reused for all calls
- **Auto-seeded** - Uses `std::random_device` for initial seed

**Performance Note:** Creating a new RNG on every call would be 100-1000x slower than reusing a thread-local RNG.

---

## Examples

> 📁 **Runnable Examples**: See [`src/docs_demos/base/docs_demo_random.cpp`](../../src/docs_demos/base/docs_demo_random.cpp)

### Monte Carlo Integration

Estimate π using random sampling:

```cpp
int insideCircle = 0;
const int samples = 10000;

for (int i = 0; i < samples; ++i) {
    Real x = Random::UniformReal(-1.0, 1.0);
    Real y = Random::UniformReal(-1.0, 1.0);
    if (x*x + y*y <= 1.0) {
        insideCircle++;
    }
}

Real piEstimate = 4.0 * insideCircle / samples;
// Typical result: ~3.14
```

### Random Walk

Simulate a 2D random walk:

```cpp
Real x = 0.0, y = 0.0;

for (int step = 0; step < 100; ++step) {
    Real dx, dy;
    Random::UniformVecDirection2(dx, dy, 1.0);  // Step size = 1
    x += dx;
    y += dy;
}

Real distance = std::sqrt(x*x + y*y);
// Expected RMS distance: sqrt(N) = sqrt(100) = 10
```

### Random Element Selection

Select random elements from a container:

```cpp
std::vector<double> data = {10, 20, 30, 40, 50};
int idx = Random::UniformInt(0, data.size() - 1);
double selected = data[idx];
```

### Particle Simulation Initial Conditions

Initialize particles with random velocities:

```cpp
for (auto& particle : particles) {
    Real vx, vy, vz;
    Random::UniformVecDirection3(vx, vy, vz, initialSpeed);
    particle.velocity = Vector3Cartesian(vx, vy, vz);
}
```

---

## Demo Functions

| Function | Description |
|----------|-------------|
| `Demo_Random_UniformReal()` | Real number generation and statistics |
| `Demo_Random_UniformInt()` | Integer generation and histogram |
| `Demo_Random_Direction2D()` | 2D direction vectors on circle |
| `Demo_Random_Direction3D()` | 3D direction vectors on sphere |
| `Demo_Random_UseCases()` | Monte Carlo, random walk, sampling |
| `Demo_Random_ThreadSafety()` | Thread safety notes |

---

## See Also
- [Vector.md](Vector.md) - Vector types for storing random data
- [BaseUtils.md](BaseUtils.md) - Additional utility functions

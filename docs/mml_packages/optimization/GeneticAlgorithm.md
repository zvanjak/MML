# Genetic Algorithm (GA)

## Overview

The Genetic Algorithm implementation in MML provides a **real-coded, mixed-integer** optimization framework for single-objective problems. It's designed for high-dimensional problems where variables can be continuous, integer, or categorical.

**Location:** `mml/algorithms/Optimization/GeneticAlgorithm.h`

## Features

- **Mixed Variable Types**: Continuous, Integer, and Categorical
- **Modular Operators**: Pluggable selection, crossover, and mutation strategies
- **Elitism**: Preserve best solutions across generations
- **Stagnation Detection**: Early termination when no improvement
- **Type-Aware Operations**: Different operators for different variable types

## Quick Start

### Basic Continuous Optimization

```cpp
#include "algorithms/Optimization/GeneticAlgorithm.h"

using namespace MML;

// Define objective function
auto sphere = [](const Vector<Real>& x) {
    Real sum = 0;
    for (int i = 0; i < x.size(); ++i)
        sum += x[i] * x[i];
    return sum;
};

// Create problem: 10 variables in [-5.12, 5.12]
auto spec = CreateContinuousProblem(10, -5.12, 5.12);

// Run GA
GeneticAlgorithm ga;
ga.SetProblem(spec);
GAResult result = ga.Minimize(sphere);

std::cout << "Best fitness: " << result.bestFitness << "\n";
std::cout << "Generations: " << result.generations << "\n";
```

### Mixed-Integer Optimization

```cpp
// Neural network hyperparameter tuning example
GAProblemSpec spec;
spec.AddContinuous(0.0001, 0.1, "learning_rate")
    .AddInteger(1, 5, "num_layers")
    .AddInteger(16, 256, "neurons_per_layer")
    .AddContinuous(0.0, 0.5, "dropout")
    .AddCategorical(3, "optimizer");  // 0=SGD, 1=Adam, 2=RMSprop

auto evaluate = [](const Vector<Real>& x) {
    Real lr = x[0];
    int layers = static_cast<int>(x[1]);
    int neurons = static_cast<int>(x[2]);
    Real dropout = x[3];
    int optimizer = static_cast<int>(x[4]);
    
    // Train model and return validation loss
    return trainAndEvaluate(lr, layers, neurons, dropout, optimizer);
};

GAConfig config;
config.populationSize = 50;
config.maxGenerations = 100;

GeneticAlgorithm ga(config);
ga.SetProblem(spec);
GAResult result = ga.Minimize(evaluate);
```

## Variable Types

### Continuous

Real-valued variables with any value in the specified range.

```cpp
// Variable x in [-10, 10]
spec.AddContinuous(-10.0, 10.0, "x");

// Or using factory
auto varSpec = VariableSpec::Continuous(-10.0, 10.0, "x");
```

### Integer

Integer-valued variables. Values are rounded during crossover and mutation.

```cpp
// Integer n in [1, 100]
spec.AddInteger(1, 100, "n");
```

### Categorical

Discrete set of options, encoded as integers 0, 1, 2, ...

```cpp
// 5 categories (values 0-4)
spec.AddCategorical(5, "category");
```

## Configuration

### GAConfig Structure

| Parameter | Default | Description |
|-----------|---------|-------------|
| `populationSize` | 100 | Number of individuals |
| `numElites` | 2 | Best individuals preserved unchanged |
| `crossoverRate` | 0.9 | Probability of crossover per pair |
| `mutationRate` | 0.1 | Per-gene mutation probability |
| `maxGenerations` | 1000 | Maximum generations |
| `maxFuncEvals` | 100000 | Max function evaluations (0=unlimited) |
| `stagnationLimit` | 100 | Generations without improvement before stop |
| `targetFitness` | -∞ | Stop if this fitness is reached |
| `fitnessThreshold` | 1e-10 | Convergence threshold |
| `minimize` | true | Minimize (true) or maximize (false) |
| `seed` | 0 | Random seed (0=time-based) |
| `verbose` | false | Print progress to stdout |
| `trackHistory` | false | Track fitness history per generation |

### Example Configuration

```cpp
GAConfig config;
config.populationSize = 200;
config.maxGenerations = 500;
config.crossoverRate = 0.8;
config.mutationRate = 0.05;
config.numElites = 5;
config.stagnationLimit = 50;
config.seed = 42;  // For reproducibility
config.trackHistory = true;

GeneticAlgorithm ga(config);
```

## Selection Strategies

### Tournament Selection (Default)

Picks k random individuals and selects the best.

```cpp
// Binary tournament (k=2) - moderate selection pressure
ga.SetSelection(std::make_unique<TournamentSelection<Real>>(2));

// Larger tournament (k=5) - higher selection pressure
ga.SetSelection(std::make_unique<TournamentSelection<Real>>(5));
```

### Rank Selection

Selection probability based on rank, not raw fitness. More stable when fitness values vary widely.

```cpp
// Selective pressure 1.5 (moderate)
ga.SetSelection(std::make_unique<RankSelection<Real>>(1.5));

// Selective pressure 2.0 (high)
ga.SetSelection(std::make_unique<RankSelection<Real>>(2.0));
```

## Crossover Operators

### BLX-α Crossover (Default)

Blend crossover samples offspring from an extended range around parents.

```cpp
// BLX-0.5 (most common)
ga.SetCrossover(std::make_unique<BLXCrossover>(0.5));

// BLX-0.0 (children strictly between parents)
ga.SetCrossover(std::make_unique<BLXCrossover>(0.0));
```

**How it works:**
- For gene i: `range = |parent1[i] - parent2[i]|`
- Child sampled from `[min - α*range, max + α*range]`

### SBX Crossover (Simulated Binary)

Mimics single-point binary crossover behavior. Widely used in NSGA-II.

```cpp
// η=15 (moderate spread)
ga.SetCrossover(std::make_unique<SBXCrossover>(15.0));

// η=2 (high spread, more exploration)
ga.SetCrossover(std::make_unique<SBXCrossover>(2.0));

// η=20 (low spread, children close to parents)
ga.SetCrossover(std::make_unique<SBXCrossover>(20.0));
```

### Uniform Crossover

Each gene independently comes from either parent.

```cpp
// 50% swap probability (default)
ga.SetCrossover(std::make_unique<UniformCrossover>(0.5));
```

## Mutation Operators

### Gaussian Mutation (Default)

Adds Gaussian noise proportional to variable range.

```cpp
// σ = 10% of variable range (default)
ga.SetMutation(std::make_unique<GaussianMutation>(0.1));

// σ = 5% of variable range (smaller perturbations)
ga.SetMutation(std::make_unique<GaussianMutation>(0.05));
```

**Type-aware behavior:**
- Continuous: Gaussian perturbation
- Integer/Categorical: Step mutation (±1) or resample

### Polynomial Mutation

Bounded polynomial distribution (NSGA-II style).

```cpp
// η=20 (small perturbations, common)
ga.SetMutation(std::make_unique<PolynomialMutation>(20.0));

// η=100 (very small perturbations)
ga.SetMutation(std::make_unique<PolynomialMutation>(100.0));
```

## Result Structure

### GAResult Fields

| Field | Type | Description |
|-------|------|-------------|
| `bestSolution` | `Vector<Real>` | Best solution found |
| `bestFitness` | `Real` | Best fitness value |
| `generations` | `int` | Generations completed |
| `funcEvals` | `int` | Total function evaluations |
| `stagnationCount` | `int` | Consecutive generations without improvement |
| `converged` | `bool` | True if converged/stagnated |
| `terminationReason` | `string` | Human-readable reason |
| `fitnessHistory` | `vector<Real>` | Best fitness per generation (if tracked) |
| `avgFitnessHistory` | `vector<Real>` | Average fitness per generation (if tracked) |

### Example Result Processing

```cpp
GAResult result = ga.Minimize(func);

std::cout << "=== GA Results ===\n";
std::cout << "Best fitness: " << result.bestFitness << "\n";
std::cout << "Generations: " << result.generations << "\n";
std::cout << "Function evals: " << result.funcEvals << "\n";
std::cout << "Termination: " << result.terminationReason << "\n";

std::cout << "Best solution: [";
for (int i = 0; i < result.bestSolution.size(); ++i) {
    if (i > 0) std::cout << ", ";
    std::cout << result.bestSolution[i];
}
std::cout << "]\n";

// If history was tracked
if (!result.fitnessHistory.empty()) {
    std::cout << "Fitness progression:\n";
    for (size_t i = 0; i < result.fitnessHistory.size(); i += 10) {
        std::cout << "  Gen " << i << ": " << result.fitnessHistory[i] << "\n";
    }
}
```

## Convenience Functions

### Quick Minimization

```cpp
// With problem spec
auto result = GAMinimize(func, spec, 100, 500);  // popSize=100, maxGen=500

// With bounds (all continuous)
Vector<Real> lb(5), ub(5);
// ... set bounds ...
auto result = GAMinimize(func, lb, ub, 100, 500);
```

### Problem Spec Factories

```cpp
// Uniform bounds
auto spec = CreateContinuousProblem(10, -5.0, 5.0);

// Different bounds per variable
Vector<Real> lb = {-1, -2, -3};
Vector<Real> ub = {1, 2, 3};
auto spec = CreateContinuousProblem(lb, ub);
```

## Benchmark Functions

The test suite includes several standard benchmark functions:

| Function | Dimension | Global Minimum | Characteristics |
|----------|-----------|----------------|-----------------|
| Sphere | n | 0 at origin | Unimodal, easy |
| Rastrigin | n | 0 at origin | Highly multimodal |
| Rosenbrock | 2 | 0 at (1,1) | Valley-shaped |
| Ackley | n | 0 at origin | Many local minima |

## Best Practices

### Population Size
- Small problems (n < 10): 50-100 individuals
- Medium problems (n < 50): 100-200 individuals
- Large problems (n > 50): 200-500 individuals

### Selection Pressure
- Start with tournament size k=3 (balanced)
- Increase to k=5 if premature convergence
- Use rank selection if fitness scaling is extreme

### Crossover Rate
- Keep high (0.8-0.95) for most problems
- Lower if population diversity is lost too quickly

### Mutation Rate
- Per-gene rate: 1/n to 0.1 (where n = dimension)
- Higher for multimodal problems
- Lower for fine-tuning near optimum

### Elitism
- Always use some elitism (2-5 individuals)
- Prevents loss of best solutions
- Ensures monotonic improvement in best fitness

## Comparison with Simulated Annealing

| Aspect | Genetic Algorithm | Simulated Annealing |
|--------|-------------------|---------------------|
| Population | Yes (parallel search) | No (single solution) |
| Memory | Higher (stores population) | Lower |
| Exploration | Crossover + mutation | Random neighbors |
| Mixed-integer | Native support | Requires custom neighbor |
| Parallelizable | Yes | Limited |
| Best for | Multimodal, mixed-type | Smooth, continuous |

## See Also

- [SimulatedAnnealing](SimulatedAnnealing.md) - Alternative metaheuristic
- [OptimizationFramework](OptimizationFramework.md) - Common optimization infrastructure
- [OptimizationConfig](OptimizationConfig.md) - Configurable termination and observers

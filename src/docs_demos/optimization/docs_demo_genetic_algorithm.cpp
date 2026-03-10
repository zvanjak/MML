///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        docs_demo_genetic_algorithm.cpp                                     ///
///  Description: Comprehensive demonstration of GeneticAlgorithm.h capabilities      ///
///               for single-objective optimization with mixed-integer variables      ///
///                                                                                   ///
///  Topics Covered:                                                                  ///
///    1. Basic GA optimization (continuous variables)                                ///
///    2. Classic benchmark functions (Rosenbrock, Rastrigin, Ackley, Schwefel)       ///
///    3. Selection strategies (Tournament, Rank)                                     ///
///    4. Crossover operators (BLX-α, SBX, Uniform)                                   ///
///    5. Mutation operators (Gaussian, Polynomial)                                   ///
///    6. Mixed-integer optimization                                                  ///
///    7. Configuration tuning (population size, rates, elitism)                      ///
///    8. Convergence monitoring and history tracking                                 ///
///    9. Practical engineering optimization example                                  ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                    ///
///////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <cmath>

#include "MMLBase.h"
#include "base/Vector/Vector.h"
#include "GeneticAlgorithm.h"
#include "test_problems/SingleObjective_TestProblems.h"

using namespace MML;
using namespace MML::Optimization;
using namespace MML::TestProblems;

///////////////////////////////////////////////////////////////////////////////////////////
// Helper function to print GA results
///////////////////////////////////////////////////////////////////////////////////////////
void PrintGAResult(const std::string& title, const GAResult& result, 
                   const Vector<Real>& trueOpt, Real trueVal)
{
    std::cout << "\n=== " << title << " ===" << std::endl;
    std::cout << std::fixed << std::setprecision(6);
    
    // Print best solution
    std::cout << "Best solution: (";
    for (int i = 0; i < result.bestSolution.size(); ++i) {
        std::cout << result.bestSolution[i];
        if (i < result.bestSolution.size() - 1) std::cout << ", ";
    }
    std::cout << ")" << std::endl;
    
    std::cout << "Best fitness: " << result.bestFitness << std::endl;
    
    // Compare to true optimum
    if (trueOpt.size() > 0) {
        Real error = 0;
        for (int i = 0; i < std::min(result.bestSolution.size(), trueOpt.size()); ++i) {
            Real diff = result.bestSolution[i] - trueOpt[i];
            error += diff * diff;
        }
        error = std::sqrt(error);
        std::cout << "Distance to optimum: " << error << std::endl;
        std::cout << "True optimal value: " << trueVal << std::endl;
    }
    
    std::cout << "Generations: " << result.generations << std::endl;
    std::cout << "Function evaluations: " << result.funcEvals << std::endl;
    std::cout << "Termination: " << result.terminationReason << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
// Example 1: Basic GA - Quadratic Function
///////////////////////////////////////////////////////////////////////////////////////////
void Example_BasicGA()
{
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Example 1: Basic Genetic Algorithm" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    std::cout << "\nMinimize f(x,y) = (x-2)^2 + (y-3)^2" << std::endl;
    std::cout << "Global minimum at (2, 3) with f = 0" << std::endl;
    std::cout << std::endl;
    
    // Simple quadratic function
    auto quadratic = [](const Vector<Real>& x) -> Real {
        Real dx = x[0] - 2.0;
        Real dy = x[1] - 3.0;
        return dx * dx + dy * dy;
    };
    
    // Create problem spec: 2D continuous in [-10, 10]
    GAProblemSpec spec = CreateContinuousProblem(2, -10.0, 10.0);
    
    // Configure GA
    GAConfig config;
    config.populationSize = 50;
    config.maxGenerations = 100;
    config.seed = 42;  // Reproducibility
    
    // Create and run GA
    GeneticAlgorithm ga(config);
    ga.SetProblem(spec);
    
    GAResult result = ga.Minimize(quadratic);
    
    Vector<Real> trueMin({2.0, 3.0});
    PrintGAResult("Basic GA Result", result, trueMin, 0.0);
}

///////////////////////////////////////////////////////////////////////////////////////////
// Example 2: Rosenbrock's Banana Function
///////////////////////////////////////////////////////////////////////////////////////////
void Example_GA_Rosenbrock()
{
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Example 2: Rosenbrock's Banana Function" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    std::cout << "\nMinimize f(x,y) = (1-x)^2 + 100*(y-x^2)^2" << std::endl;
    std::cout << "Global minimum at (1, 1) with f = 0" << std::endl;
    std::cout << "\nRosenbrock is notoriously difficult due to its narrow" << std::endl;
    std::cout << "curved valley - a classic test for global optimizers." << std::endl;
    std::cout << std::endl;
    
    // Use centralized test function
    auto rosenbrock = Functions::Rosenbrock();
    
    GAProblemSpec spec = CreateContinuousProblem(2, -5.0, 5.0);
    
    GAConfig config;
    config.populationSize = 100;
    config.maxGenerations = 500;
    config.seed = 42;
    
    GeneticAlgorithm ga(config);
    ga.SetProblem(spec);
    
    GAResult result = ga.Minimize(rosenbrock);
    
    Vector<Real> trueMin({1.0, 1.0});
    PrintGAResult("Rosenbrock Result", result, trueMin, 0.0);
}

///////////////////////////////////////////////////////////////////////////////////////////
// Example 3: Rastrigin - Highly Multimodal
///////////////////////////////////////////////////////////////////////////////////////////
void Example_GA_Rastrigin()
{
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Example 3: Rastrigin Function (Highly Multimodal)" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    std::cout << "\nMinimize f(x) = 10n + sum(x_i^2 - 10*cos(2*pi*x_i))" << std::endl;
    std::cout << "Global minimum at origin with f = 0" << std::endl;
    std::cout << "\nRastrigin has many local minima in a regular pattern." << std::endl;
    std::cout << "Tests the GA's ability to escape local optima." << std::endl;
    std::cout << std::endl;
    
    auto rastrigin = Functions::Rastrigin();
    
    // 5D Rastrigin
    GAProblemSpec spec = CreateContinuousProblem(5, -5.12, 5.12);
    
    GAConfig config;
    config.populationSize = 150;
    config.maxGenerations = 500;
    config.seed = 42;
    
    GeneticAlgorithm ga(config);
    ga.SetProblem(spec);
    
    GAResult result = ga.Minimize(rastrigin);
    
    Vector<Real> trueMin({0.0, 0.0, 0.0, 0.0, 0.0});
    PrintGAResult("Rastrigin 5D Result", result, trueMin, 0.0);
}

///////////////////////////////////////////////////////////////////////////////////////////
// Example 4: Comparing Selection Strategies
///////////////////////////////////////////////////////////////////////////////////////////
void Example_SelectionStrategies()
{
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Example 4: Comparing Selection Strategies" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    std::cout << "\nComparing Tournament vs Rank selection on Sphere function" << std::endl;
    std::cout << "f(x) = sum(x_i^2), minimum at origin" << std::endl;
    std::cout << std::endl;
    
    auto sphere = Functions::Sphere();
    GAProblemSpec spec = CreateContinuousProblem(5, -5.12, 5.12);
    
    GAConfig config;
    config.populationSize = 50;
    config.maxGenerations = 200;
    config.seed = 42;
    
    // 1. Tournament Selection (k=2, binary)
    std::cout << "--- Tournament Selection (k=2) ---" << std::endl;
    {
        GeneticAlgorithm ga(config);
        ga.SetProblem(spec);
        ga.SetSelection(std::make_unique<TournamentSelection<Real>>(2, 42));
        
        GAResult result = ga.Minimize(sphere);
        std::cout << "Best fitness: " << std::fixed << std::setprecision(8) 
                  << result.bestFitness << std::endl;
        std::cout << "Generations: " << result.generations << std::endl;
    }
    
    // 2. Tournament Selection (k=5, higher pressure)
    std::cout << "\n--- Tournament Selection (k=5) ---" << std::endl;
    {
        GeneticAlgorithm ga(config);
        ga.SetProblem(spec);
        ga.SetSelection(std::make_unique<TournamentSelection<Real>>(5, 42));
        
        GAResult result = ga.Minimize(sphere);
        std::cout << "Best fitness: " << std::fixed << std::setprecision(8) 
                  << result.bestFitness << std::endl;
        std::cout << "Generations: " << result.generations << std::endl;
    }
    
    // 3. Rank Selection
    std::cout << "\n--- Rank Selection (pressure=1.5) ---" << std::endl;
    {
        GeneticAlgorithm ga(config);
        ga.SetProblem(spec);
        ga.SetSelection(std::make_unique<RankSelection<Real>>(1.5, 42));
        
        GAResult result = ga.Minimize(sphere);
        std::cout << "Best fitness: " << std::fixed << std::setprecision(8) 
                  << result.bestFitness << std::endl;
        std::cout << "Generations: " << result.generations << std::endl;
    }
    
    // 4. Rank Selection (higher pressure)
    std::cout << "\n--- Rank Selection (pressure=2.0) ---" << std::endl;
    {
        GeneticAlgorithm ga(config);
        ga.SetProblem(spec);
        ga.SetSelection(std::make_unique<RankSelection<Real>>(2.0, 42));
        
        GAResult result = ga.Minimize(sphere);
        std::cout << "Best fitness: " << std::fixed << std::setprecision(8) 
                  << result.bestFitness << std::endl;
        std::cout << "Generations: " << result.generations << std::endl;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// Example 5: Comparing Crossover Operators
///////////////////////////////////////////////////////////////////////////////////////////
void Example_CrossoverOperators()
{
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Example 5: Comparing Crossover Operators" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    std::cout << "\nComparing BLX-α, SBX, and Uniform crossover on Ackley function" << std::endl;
    std::cout << "f(x) = -20*exp(-0.2*sqrt(mean(x^2))) - exp(mean(cos(2*pi*x))) + 20 + e" << std::endl;
    std::cout << "Global minimum at origin with f = 0" << std::endl;
    std::cout << std::endl;
    
    auto ackley = Functions::Ackley();
    GAProblemSpec spec = CreateContinuousProblem(5, -5.0, 5.0);
    
    GAConfig config;
    config.populationSize = 80;
    config.maxGenerations = 300;
    config.seed = 42;
    
    // 1. BLX-α Crossover (α=0.5)
    std::cout << "--- BLX-α Crossover (α=0.5) ---" << std::endl;
    {
        GeneticAlgorithm ga(config);
        ga.SetProblem(spec);
        ga.SetCrossover(std::make_unique<BLXCrossover>(0.5, 42));
        
        GAResult result = ga.Minimize(ackley);
        std::cout << "Best fitness: " << std::fixed << std::setprecision(8) 
                  << result.bestFitness << std::endl;
        std::cout << "Generations: " << result.generations << std::endl;
    }
    
    // 2. SBX Crossover (η=15)
    std::cout << "\n--- SBX Crossover (η=15) ---" << std::endl;
    {
        GeneticAlgorithm ga(config);
        ga.SetProblem(spec);
        ga.SetCrossover(std::make_unique<SBXCrossover>(15.0, 42));
        
        GAResult result = ga.Minimize(ackley);
        std::cout << "Best fitness: " << std::fixed << std::setprecision(8) 
                  << result.bestFitness << std::endl;
        std::cout << "Generations: " << result.generations << std::endl;
    }
    
    // 3. SBX Crossover (η=2, more exploration)
    std::cout << "\n--- SBX Crossover (η=2, more exploration) ---" << std::endl;
    {
        GeneticAlgorithm ga(config);
        ga.SetProblem(spec);
        ga.SetCrossover(std::make_unique<SBXCrossover>(2.0, 42));
        
        GAResult result = ga.Minimize(ackley);
        std::cout << "Best fitness: " << std::fixed << std::setprecision(8) 
                  << result.bestFitness << std::endl;
        std::cout << "Generations: " << result.generations << std::endl;
    }
    
    // 4. Uniform Crossover
    std::cout << "\n--- Uniform Crossover ---" << std::endl;
    {
        GeneticAlgorithm ga(config);
        ga.SetProblem(spec);
        ga.SetCrossover(std::make_unique<UniformCrossover>(0.5, 42));
        
        GAResult result = ga.Minimize(ackley);
        std::cout << "Best fitness: " << std::fixed << std::setprecision(8) 
                  << result.bestFitness << std::endl;
        std::cout << "Generations: " << result.generations << std::endl;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// Example 6: Comparing Mutation Operators
///////////////////////////////////////////////////////////////////////////////////////////
void Example_MutationOperators()
{
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Example 6: Comparing Mutation Operators" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    std::cout << "\nComparing Gaussian vs Polynomial mutation on Griewank function" << std::endl;
    std::cout << "f(x) = sum(x_i^2)/4000 - prod(cos(x_i/sqrt(i+1))) + 1" << std::endl;
    std::cout << "Global minimum at origin with f = 0" << std::endl;
    std::cout << std::endl;
    
    auto griewank = Functions::Griewank();
    GAProblemSpec spec = CreateContinuousProblem(5, -10.0, 10.0);
    
    GAConfig config;
    config.populationSize = 80;
    config.maxGenerations = 300;
    config.seed = 42;
    
    // 1. Gaussian Mutation (σ=0.1)
    std::cout << "--- Gaussian Mutation (σ=0.1) ---" << std::endl;
    {
        GeneticAlgorithm ga(config);
        ga.SetProblem(spec);
        ga.SetMutation(std::make_unique<GaussianMutation>(0.1, 42));
        
        GAResult result = ga.Minimize(griewank);
        std::cout << "Best fitness: " << std::fixed << std::setprecision(8) 
                  << result.bestFitness << std::endl;
        std::cout << "Generations: " << result.generations << std::endl;
    }
    
    // 2. Gaussian Mutation (σ=0.2, larger steps)
    std::cout << "\n--- Gaussian Mutation (σ=0.2) ---" << std::endl;
    {
        GeneticAlgorithm ga(config);
        ga.SetProblem(spec);
        ga.SetMutation(std::make_unique<GaussianMutation>(0.2, 42));
        
        GAResult result = ga.Minimize(griewank);
        std::cout << "Best fitness: " << std::fixed << std::setprecision(8) 
                  << result.bestFitness << std::endl;
        std::cout << "Generations: " << result.generations << std::endl;
    }
    
    // 3. Polynomial Mutation (η=20)
    std::cout << "\n--- Polynomial Mutation (η=20) ---" << std::endl;
    {
        GeneticAlgorithm ga(config);
        ga.SetProblem(spec);
        ga.SetMutation(std::make_unique<PolynomialMutation>(20.0, 42));
        
        GAResult result = ga.Minimize(griewank);
        std::cout << "Best fitness: " << std::fixed << std::setprecision(8) 
                  << result.bestFitness << std::endl;
        std::cout << "Generations: " << result.generations << std::endl;
    }
    
    // 4. Polynomial Mutation (η=5, larger perturbations)
    std::cout << "\n--- Polynomial Mutation (η=5) ---" << std::endl;
    {
        GeneticAlgorithm ga(config);
        ga.SetProblem(spec);
        ga.SetMutation(std::make_unique<PolynomialMutation>(5.0, 42));
        
        GAResult result = ga.Minimize(griewank);
        std::cout << "Best fitness: " << std::fixed << std::setprecision(8) 
                  << result.bestFitness << std::endl;
        std::cout << "Generations: " << result.generations << std::endl;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// Example 7: Mixed-Integer Optimization
///////////////////////////////////////////////////////////////////////////////////////////
void Example_MixedInteger()
{
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Example 7: Mixed-Integer Optimization" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    std::cout << "\nMinimize f(x, n, c) = (x - 2.5)^2 + (n - 3)^2 + category_penalty(c)" << std::endl;
    std::cout << "where x is continuous, n is integer, c is categorical" << std::endl;
    std::cout << "Optimal: x=2.5, n=3, c=2 (third category)" << std::endl;
    std::cout << std::endl;
    
    // Mixed-integer objective:
    // x: continuous in [0, 5], optimal at 2.5
    // n: integer in [1, 5], optimal at 3
    // c: categorical (3 options), optimal at category 2 (index 2)
    auto mixedObj = [](const Vector<Real>& x) -> Real {
        Real continuous = (x[0] - 2.5) * (x[0] - 2.5);
        Real integer = (x[1] - 3.0) * (x[1] - 3.0);
        
        // Category penalty: category 2 is best (penalty = 0)
        int category = static_cast<int>(std::round(x[2]));
        Real categoryPenalty = (category == 2) ? 0.0 : 1.0 + std::abs(category - 2);
        
        return continuous + integer + categoryPenalty;
    };
    
    // Build problem specification
    GAProblemSpec spec;
    spec.AddContinuous(0.0, 5.0, "x");      // Continuous variable
    spec.AddInteger(1, 5, "n");              // Integer variable
    spec.AddCategorical(4, "category");      // 4 categories (0, 1, 2, 3)
    
    GAConfig config;
    config.populationSize = 50;
    config.maxGenerations = 100;
    config.seed = 42;
    
    GeneticAlgorithm ga(config);
    ga.SetProblem(spec);
    
    GAResult result = ga.Minimize(mixedObj);
    
    std::cout << "=== Mixed-Integer Result ===" << std::endl;
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "x (continuous): " << result.bestSolution[0] << " (optimal: 2.5)" << std::endl;
    std::cout << "n (integer):    " << result.bestSolution[1] << " (optimal: 3)" << std::endl;
    std::cout << "c (category):   " << static_cast<int>(result.bestSolution[2]) << " (optimal: 2)" << std::endl;
    std::cout << "Best fitness:   " << result.bestFitness << std::endl;
    std::cout << "Generations:    " << result.generations << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
// Example 8: Configuration Tuning - Population Size and Rates
///////////////////////////////////////////////////////////////////////////////////////////
void Example_ConfigurationTuning()
{
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Example 8: Configuration Tuning" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    std::cout << "\nStudying effect of population size and mutation rate on Schwefel function" << std::endl;
    std::cout << "f(x) = 418.9829*n - sum(x_i * sin(sqrt(|x_i|)))" << std::endl;
    std::cout << "Global minimum at x_i = 420.9687 with f ≈ 0" << std::endl;
    std::cout << std::endl;
    
    auto schwefel = Functions::Schwefel();
    GAProblemSpec spec = CreateContinuousProblem(3, -500.0, 500.0);
    
    // Test different configurations
    std::cout << "--- Effect of Population Size (mutation rate = 0.1) ---" << std::endl;
    for (int popSize : {30, 50, 100, 200}) {
        GAConfig config;
        config.populationSize = popSize;
        config.maxGenerations = 300;
        config.mutationRate = 0.1;
        config.seed = 42;
        
        GeneticAlgorithm ga(config);
        ga.SetProblem(spec);
        
        GAResult result = ga.Minimize(schwefel);
        std::cout << "Pop=" << std::setw(3) << popSize 
                  << ": best=" << std::fixed << std::setprecision(4) << result.bestFitness
                  << ", evals=" << result.funcEvals << std::endl;
    }
    
    std::cout << "\n--- Effect of Mutation Rate (population = 100) ---" << std::endl;
    for (Real mutRate : {0.01, 0.05, 0.1, 0.2, 0.3}) {
        GAConfig config;
        config.populationSize = 100;
        config.maxGenerations = 300;
        config.mutationRate = mutRate;
        config.seed = 42;
        
        GeneticAlgorithm ga(config);
        ga.SetProblem(spec);
        
        GAResult result = ga.Minimize(schwefel);
        std::cout << "Mutation=" << std::fixed << std::setprecision(2) << mutRate 
                  << ": best=" << std::setprecision(4) << result.bestFitness << std::endl;
    }
    
    std::cout << "\n--- Effect of Elitism ---" << std::endl;
    for (int elites : {0, 1, 2, 5, 10}) {
        GAConfig config;
        config.populationSize = 100;
        config.maxGenerations = 300;
        config.numElites = elites;
        config.seed = 42;
        
        GeneticAlgorithm ga(config);
        ga.SetProblem(spec);
        
        GAResult result = ga.Minimize(schwefel);
        std::cout << "Elites=" << std::setw(2) << elites 
                  << ": best=" << std::fixed << std::setprecision(4) << result.bestFitness << std::endl;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// Example 9: Convergence History Tracking
///////////////////////////////////////////////////////////////////////////////////////////
void Example_ConvergenceHistory()
{
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Example 9: Convergence History Tracking" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    std::cout << "\nTracking fitness evolution over generations on Styblinski-Tang function" << std::endl;
    std::cout << "f(x) = 0.5 * sum(x_i^4 - 16*x_i^2 + 5*x_i)" << std::endl;
    std::cout << "Global minimum at x_i ≈ -2.903534 with f/n ≈ -39.166" << std::endl;
    std::cout << std::endl;
    
    auto styblinski = Functions::StyblinskiTang();
    GAProblemSpec spec = CreateContinuousProblem(3, -5.0, 5.0);
    
    GAConfig config;
    config.populationSize = 80;
    config.maxGenerations = 150;
    config.trackHistory = true;  // Enable history tracking
    config.seed = 42;
    
    GeneticAlgorithm ga(config);
    ga.SetProblem(spec);
    
    GAResult result = ga.Minimize(styblinski);
    
    std::cout << "=== Convergence History (every 10 generations) ===" << std::endl;
    std::cout << std::setw(6) << "Gen" << std::setw(15) << "Best" << std::setw(15) << "Average" << std::endl;
    std::cout << std::string(36, '-') << std::endl;
    
    for (size_t i = 0; i < result.fitnessHistory.size(); i += 10) {
        std::cout << std::setw(6) << i 
                  << std::setw(15) << std::fixed << std::setprecision(4) << result.fitnessHistory[i]
                  << std::setw(15) << result.avgFitnessHistory[i] << std::endl;
    }
    
    // Print final
    size_t last = result.fitnessHistory.size() - 1;
    std::cout << std::setw(6) << last 
              << std::setw(15) << std::fixed << std::setprecision(4) << result.fitnessHistory[last]
              << std::setw(15) << result.avgFitnessHistory[last] << std::endl;
    
    std::cout << "\nFinal solution: (";
    for (int i = 0; i < result.bestSolution.size(); ++i) {
        std::cout << std::setprecision(6) << result.bestSolution[i];
        if (i < result.bestSolution.size() - 1) std::cout << ", ";
    }
    std::cout << ")" << std::endl;
    std::cout << "True optimum: x_i ≈ -2.903534" << std::endl;
    std::cout << "Final fitness: " << result.bestFitness << std::endl;
    std::cout << "True minimum (3D): " << -39.16616570377 * 3 << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
// Example 10: Stagnation Detection and Early Termination
///////////////////////////////////////////////////////////////////////////////////////////
void Example_StagnationDetection()
{
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Example 10: Stagnation Detection" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    std::cout << "\nDemonstrating early termination when population converges" << std::endl;
    std::cout << "Using Sphere function (easy to converge)" << std::endl;
    std::cout << std::endl;
    
    auto sphere = Functions::Sphere();
    GAProblemSpec spec = CreateContinuousProblem(3, -5.12, 5.12);
    
    // Without stagnation limit
    std::cout << "--- Without stagnation limit (maxGen=500) ---" << std::endl;
    {
        GAConfig config;
        config.populationSize = 50;
        config.maxGenerations = 500;
        config.stagnationLimit = 1000;  // Effectively disabled
        config.seed = 42;
        
        GeneticAlgorithm ga(config);
        ga.SetProblem(spec);
        
        GAResult result = ga.Minimize(sphere);
        std::cout << "Generations: " << result.generations << std::endl;
        std::cout << "Best fitness: " << std::scientific << result.bestFitness << std::endl;
        std::cout << "Termination: " << result.terminationReason << std::endl;
    }
    
    // With stagnation limit
    std::cout << "\n--- With stagnation limit = 30 ---" << std::endl;
    {
        GAConfig config;
        config.populationSize = 50;
        config.maxGenerations = 500;
        config.stagnationLimit = 30;
        config.seed = 42;
        
        GeneticAlgorithm ga(config);
        ga.SetProblem(spec);
        
        GAResult result = ga.Minimize(sphere);
        std::cout << "Generations: " << result.generations << std::endl;
        std::cout << "Best fitness: " << std::scientific << result.bestFitness << std::endl;
        std::cout << "Stagnation count: " << result.stagnationCount << std::endl;
        std::cout << "Termination: " << result.terminationReason << std::endl;
    }
    
    // With target fitness
    std::cout << "\n--- With target fitness = 0.001 ---" << std::endl;
    {
        GAConfig config;
        config.populationSize = 50;
        config.maxGenerations = 500;
        config.targetFitness = 0.001;
        config.seed = 42;
        
        GeneticAlgorithm ga(config);
        ga.SetProblem(spec);
        
        GAResult result = ga.Minimize(sphere);
        std::cout << "Generations: " << result.generations << std::endl;
        std::cout << "Best fitness: " << std::scientific << result.bestFitness << std::endl;
        std::cout << "Termination: " << result.terminationReason << std::endl;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// Example 11: Higher Dimensional Problem
///////////////////////////////////////////////////////////////////////////////////////////
void Example_GA_HigherDimensional()
{
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Example 11: Higher Dimensional Optimization (10D)" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    std::cout << "\nOptimizing 10D Levy function" << std::endl;
    std::cout << "f(x) involves sin^2 and weighted sum terms" << std::endl;
    std::cout << "Global minimum at x_i = 1 with f = 0" << std::endl;
    std::cout << std::endl;
    
    auto levy = Functions::Levy();
    GAProblemSpec spec = CreateContinuousProblem(10, -10.0, 10.0);
    
    GAConfig config;
    config.populationSize = 200;
    config.maxGenerations = 800;
    config.numElites = 5;
    config.seed = 42;
    
    // Use NSGA-II style operators for better performance
    GeneticAlgorithm ga(config);
    ga.SetProblem(spec);
    ga.SetCrossover(std::make_unique<SBXCrossover>(15.0, 42));
    ga.SetMutation(std::make_unique<PolynomialMutation>(20.0, 42));
    
    GAResult result = ga.Minimize(levy);
    
    std::cout << "=== 10D Levy Result ===" << std::endl;
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Best fitness: " << result.bestFitness << " (optimal: 0)" << std::endl;
    std::cout << "Generations: " << result.generations << std::endl;
    std::cout << "Function evaluations: " << result.funcEvals << std::endl;
    
    std::cout << "\nSolution (first 5 components): ";
    for (int i = 0; i < 5; ++i) {
        std::cout << result.bestSolution[i] << " ";
    }
    std::cout << "..." << std::endl;
    std::cout << "True optimum: all x_i = 1" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
// Example 12: Practical Application - Beam Design Optimization
///////////////////////////////////////////////////////////////////////////////////////////
void Example_BeamDesignOptimization()
{
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Example 12: Practical Application - Beam Design" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    std::cout << "\nMinimize beam weight subject to stress and deflection constraints" << std::endl;
    std::cout << "Variables: width (b), height (h)" << std::endl;
    std::cout << "Weight = ρ * b * h * L" << std::endl;
    std::cout << "Stress constraint: σ = 6*P*L/(b*h^2) <= σ_max" << std::endl;
    std::cout << "Deflection constraint: δ = 4*P*L^3/(E*b*h^3) <= δ_max" << std::endl;
    std::cout << std::endl;
    
    // Problem parameters
    const Real L = 1.0;           // Length (m)
    const Real P = 10000.0;       // Load (N)
    const Real E = 200e9;         // Young's modulus (Pa) - steel
    const Real rho = 7850.0;      // Density (kg/m^3) - steel
    const Real sigma_max = 250e6; // Max stress (Pa)
    const Real delta_max = 0.005; // Max deflection (m)
    
    // Objective: minimize weight + penalty for constraint violations
    auto beamDesign = [&](const Vector<Real>& x) -> Real {
        Real b = x[0];  // width (m)
        Real h = x[1];  // height (m)
        
        // Weight objective
        Real weight = rho * b * h * L;
        
        // Stress constraint: σ <= σ_max
        Real sigma = 6.0 * P * L / (b * h * h);
        Real stressViolation = std::max(0.0, sigma - sigma_max);
        
        // Deflection constraint: δ <= δ_max
        Real delta = 4.0 * P * L * L * L / (E * b * h * h * h);
        Real deflectionViolation = std::max(0.0, delta - delta_max);
        
        // Penalty method for constraints
        Real penalty = 1e6 * (stressViolation / sigma_max + deflectionViolation / delta_max);
        
        return weight + penalty;
    };
    
    // Variables: b in [0.01, 0.2] m, h in [0.02, 0.3] m
    GAProblemSpec spec;
    spec.AddContinuous(0.01, 0.2, "width (b)");
    spec.AddContinuous(0.02, 0.3, "height (h)");
    
    GAConfig config;
    config.populationSize = 100;
    config.maxGenerations = 200;
    config.seed = 42;
    
    GeneticAlgorithm ga(config);
    ga.SetProblem(spec);
    
    GAResult result = ga.Minimize(beamDesign);
    
    Real b = result.bestSolution[0];
    Real h = result.bestSolution[1];
    Real weight = rho * b * h * L;
    Real sigma = 6.0 * P * L / (b * h * h);
    Real delta = 4.0 * P * L * L * L / (E * b * h * h * h);
    
    std::cout << "=== Beam Design Result ===" << std::endl;
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Optimal width (b):  " << b * 1000 << " mm" << std::endl;
    std::cout << "Optimal height (h): " << h * 1000 << " mm" << std::endl;
    std::cout << "Weight:             " << weight << " kg" << std::endl;
    std::cout << "Stress:             " << sigma / 1e6 << " MPa (max: " << sigma_max / 1e6 << " MPa)" << std::endl;
    std::cout << "Deflection:         " << delta * 1000 << " mm (max: " << delta_max * 1000 << " mm)" << std::endl;
    std::cout << "Stress margin:      " << (1 - sigma / sigma_max) * 100 << "%" << std::endl;
    std::cout << "Deflection margin:  " << (1 - delta / delta_max) * 100 << "%" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
// Example 13: Convenience Function Usage
///////////////////////////////////////////////////////////////////////////////////////////
void Example_ConvenienceFunctions()
{
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Example 13: Convenience Functions" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    std::cout << "\nUsing GAMinimize() convenience function for quick optimization" << std::endl;
    std::cout << std::endl;
    
    auto booth = Functions::Booth();
    
    // Method 1: With explicit bounds vectors
    std::cout << "--- Method 1: With bound vectors ---" << std::endl;
    {
        Vector<Real> lb(2), ub(2);
        lb[0] = -10.0; lb[1] = -10.0;
        ub[0] = 10.0;  ub[1] = 10.0;
        
        auto result = GAMinimize(booth, lb, ub, 50, 100);
        
        std::cout << "Best solution: (" << result.bestSolution[0] << ", " 
                  << result.bestSolution[1] << ")" << std::endl;
        std::cout << "Best fitness: " << result.bestFitness << std::endl;
        std::cout << "True optimum: (1, 3) with f = 0" << std::endl;
    }
    
    // Method 2: With problem specification
    std::cout << "\n--- Method 2: With problem specification ---" << std::endl;
    {
        auto spec = CreateContinuousProblem(2, -10.0, 10.0);
        auto result = GAMinimize(booth, spec, 50, 100);
        
        std::cout << "Best solution: (" << result.bestSolution[0] << ", " 
                  << result.bestSolution[1] << ")" << std::endl;
        std::cout << "Best fitness: " << result.bestFitness << std::endl;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// Main - Run all demos
///////////////////////////////////////////////////////////////////////////////////////////
void Docs_Demo_GeneticAlgorithm()
{
    std::cout << "\n";
    std::cout << "######################################################################" << std::endl;
    std::cout << "###                                                                ###" << std::endl;
    std::cout << "###       GENETIC ALGORITHM COMPREHENSIVE DEMONSTRATION            ###" << std::endl;
    std::cout << "###                                                                ###" << std::endl;
    std::cout << "###   Single-objective optimization with mixed-integer variables   ###" << std::endl;
    std::cout << "###   Features: Selection, Crossover, Mutation strategies          ###" << std::endl;
    std::cout << "###                                                                ###" << std::endl;
    std::cout << "######################################################################" << std::endl;
    
    Example_BasicGA();
    Example_GA_Rosenbrock();
    Example_GA_Rastrigin();
    Example_SelectionStrategies();
    Example_CrossoverOperators();
    Example_MutationOperators();
    Example_MixedInteger();
    Example_ConfigurationTuning();
    Example_ConvergenceHistory();
    Example_StagnationDetection();
    Example_GA_HigherDimensional();
    Example_BeamDesignOptimization();
    Example_ConvenienceFunctions();
    
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Genetic Algorithm demonstration complete!" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
}

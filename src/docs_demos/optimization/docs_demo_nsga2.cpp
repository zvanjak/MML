///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        docs_demo_nsga2.cpp                                                 ///
///  Description: Comprehensive demonstration of NSGA-II multi-objective optimizer    ///
///               based on Deb et al. (2002) "A fast and elitist multiobjective       ///
///               genetic algorithm: NSGA-II"                                         ///
///                                                                                   ///
///  Topics Covered:                                                                  ///
///    1. Basic NSGA-II optimization with ZDT test problems                           ///
///    2. Pareto front visualization and metrics                                      ///
///    3. Different Pareto front shapes (convex, concave, disconnected)               ///
///    4. Hypervolume quality indicator                                               ///
///    5. Configuration tuning (population size, generations, operator params)        ///
///    6. Convergence analysis and ideal/nadir points                                 ///
///    7. Practical bi-objective engineering example                                  ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                    ///
///////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>

#include "MMLBase.h"
#include "base/Vector/VectorN.h"
#include "NSGA2.h"
#include "MultiObjective.h"
#include "test_problems/MultiObjectiveConstrained_TestProblems.h"

using namespace MML;
using namespace MML::Optimization;
using namespace MML::TestProblems;

///////////////////////////////////////////////////////////////////////////////////////////
// Helper function to print NSGA-II results
///////////////////////////////////////////////////////////////////////////////////////////
template<int N, int M>
void PrintNSGA2Result(const std::string& title, const NSGA2Result<N, M>& result,
                      const VectorN<Real, M>& referencePoint)
{
    std::cout << "\n=== " << title << " ===" << std::endl;
    std::cout << std::fixed << std::setprecision(6);
    
    std::cout << "Pareto front size: " << result.paretoFront.size() << std::endl;
    std::cout << "Generations: " << result.generations << std::endl;
    std::cout << "Function evaluations: " << result.functionEvaluations << std::endl;
    
    // Calculate and print ideal/nadir points
    auto ideal = result.GetIdealPoint();
    auto nadir = result.GetNadirPoint();
    
    std::cout << "Ideal point: (";
    for (int m = 0; m < M; ++m) {
        std::cout << ideal[m];
        if (m < M - 1) std::cout << ", ";
    }
    std::cout << ")" << std::endl;
    
    std::cout << "Nadir point: (";
    for (int m = 0; m < M; ++m) {
        std::cout << nadir[m];
        if (m < M - 1) std::cout << ", ";
    }
    std::cout << ")" << std::endl;
    
    // Hypervolume for 2D problems
    if constexpr (M == 2) {
        Real hv = result.GetHypervolume2D(referencePoint);
        std::cout << "Hypervolume: " << hv << " (reference: " 
                  << referencePoint[0] << ", " << referencePoint[1] << ")" << std::endl;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// Helper to print selected Pareto front points
///////////////////////////////////////////////////////////////////////////////////////////
template<int N, int M>
void PrintParetoFrontSample(const NSGA2Result<N, M>& result, int numSamples = 10)
{
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "\nSample Pareto front points (f1, f2):" << std::endl;
    std::cout << std::string(35, '-') << std::endl;
    
    const auto& front = result.paretoFront;
    int step = std::max(1, static_cast<int>(front.size()) / numSamples);
    
    // Sort by first objective for display
    std::vector<const ObjectivePoint<N, M>*> sorted;
    for (const auto& pt : front) {
        sorted.push_back(&pt);
    }
    std::sort(sorted.begin(), sorted.end(),
        [](const auto* a, const auto* b) {
            return a->GetObjective(0) < b->GetObjective(0);
        });
    
    int count = 0;
    for (size_t i = 0; i < sorted.size() && count < numSamples; i += step) {
        std::cout << "  (" << std::setw(8) << sorted[i]->GetObjective(0) 
                  << ", " << std::setw(8) << sorted[i]->GetObjective(1) << ")" << std::endl;
        ++count;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// Example 1: Basic NSGA-II with ZDT1 (Convex Pareto Front)
///////////////////////////////////////////////////////////////////////////////////////////
void Example_ZDT1_Basic()
{
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Example 1: ZDT1 - Convex Pareto Front" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    std::cout << "\nZDT1 is a classic bi-objective test problem with:" << std::endl;
    std::cout << "  - 30 decision variables (N=30)" << std::endl;
    std::cout << "  - 2 objectives to minimize" << std::endl;
    std::cout << "  - Convex Pareto front: f2 = 1 - sqrt(f1)" << std::endl;
    std::cout << "  - Bounds: all x_i in [0, 1]" << std::endl;
    std::cout << std::endl;
    
    // Create problem
    ZDT1<30> problem;
    
    // Create problem specification for variable bounds
    GAProblemSpec spec;
    for (int i = 0; i < 30; ++i) {
        auto [lb, ub] = problem.GetVariableBounds(i);
        spec.AddContinuous(lb, ub, "x" + std::to_string(i));
    }
    
    // Configure NSGA-II
    NSGA2Config config;
    config.populationSize = 100;
    config.maxGenerations = 250;
    config.seed = 42;
    
    // Create and run optimizer
    NSGA2<30, 2> nsga(config);
    nsga.SetProblem(spec);
    
    auto result = nsga.Optimize(problem);
    
    // Print results
    VectorN<Real, 2> refPoint{1.1, 1.1};
    PrintNSGA2Result("ZDT1 Result", result, refPoint);
    PrintParetoFrontSample(result);
    
    // Verify convergence to true front
    std::cout << "\nCompare to true Pareto front f2 = 1 - sqrt(f1):" << std::endl;
    std::cout << "  f1=0.0: found=" << result.GetIdealPoint()[0] << ", true=0.0" << std::endl;
    std::cout << "  At f1=0.5, f2 should be ~0.293 (1 - sqrt(0.5))" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
// Example 2: ZDT2 - Concave Pareto Front
///////////////////////////////////////////////////////////////////////////////////////////
void Example_ZDT2_Concave()
{
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Example 2: ZDT2 - Concave (Non-Convex) Pareto Front" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    std::cout << "\nZDT2 tests algorithm behavior on non-convex fronts:" << std::endl;
    std::cout << "  - Concave Pareto front: f2 = 1 - f1^2" << std::endl;
    std::cout << "  - More challenging than ZDT1 for some algorithms" << std::endl;
    std::cout << std::endl;
    
    ZDT2<30> problem;
    
    GAProblemSpec spec;
    for (int i = 0; i < 30; ++i) {
        auto [lb, ub] = problem.GetVariableBounds(i);
        spec.AddContinuous(lb, ub, "x" + std::to_string(i));
    }
    
    NSGA2Config config;
    config.populationSize = 100;
    config.maxGenerations = 250;
    config.seed = 42;
    
    NSGA2<30, 2> nsga(config);
    nsga.SetProblem(spec);
    
    auto result = nsga.Optimize(problem);
    
    VectorN<Real, 2> refPoint{1.1, 1.1};
    PrintNSGA2Result("ZDT2 Result", result, refPoint);
    PrintParetoFrontSample(result);
    
    std::cout << "\nTrue Pareto front: f2 = 1 - f1^2 (concave curve)" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
// Example 3: ZDT3 - Disconnected Pareto Front
///////////////////////////////////////////////////////////////////////////////////////////
void Example_ZDT3_Disconnected()
{
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Example 3: ZDT3 - Disconnected Pareto Front" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    std::cout << "\nZDT3 has a discontinuous Pareto front with multiple segments:" << std::endl;
    std::cout << "  - f2 = 1 - sqrt(f1) - f1*sin(10*pi*f1)" << std::endl;
    std::cout << "  - Tests the algorithm's ability to find all front segments" << std::endl;
    std::cout << std::endl;
    
    ZDT3<30> problem;
    
    GAProblemSpec spec;
    for (int i = 0; i < 30; ++i) {
        auto [lb, ub] = problem.GetVariableBounds(i);
        spec.AddContinuous(lb, ub, "x" + std::to_string(i));
    }
    
    NSGA2Config config;
    config.populationSize = 100;
    config.maxGenerations = 250;
    config.seed = 42;
    
    NSGA2<30, 2> nsga(config);
    nsga.SetProblem(spec);
    
    auto result = nsga.Optimize(problem);
    
    VectorN<Real, 2> refPoint{1.0, 2.0};  // Different ref point due to negative f2 values
    PrintNSGA2Result("ZDT3 Result", result, refPoint);
    PrintParetoFrontSample(result, 15);  // More samples to show discontinuity
    
    std::cout << "\nNote: ZDT3's Pareto front has 5 disconnected segments" << std::endl;
    std::cout << "with f2 values ranging from about -0.77 to 1.0" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
// Example 4: ZDT4 - Multimodal (Many Local Fronts)
///////////////////////////////////////////////////////////////////////////////////////////
void Example_ZDT4_Multimodal()
{
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Example 4: ZDT4 - Multimodal (Many Local Pareto Fronts)" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    std::cout << "\nZDT4 is the most difficult ZDT problem:" << std::endl;
    std::cout << "  - Has 21^9 local Pareto fronts (for N=10)" << std::endl;
    std::cout << "  - x_i in [-5, 5] for i > 1 (larger search space)" << std::endl;
    std::cout << "  - Tests global vs local optima finding" << std::endl;
    std::cout << std::endl;
    
    ZDT4<10> problem;  // Use N=10 (standard for ZDT4)
    
    GAProblemSpec spec;
    for (int i = 0; i < 10; ++i) {
        auto [lb, ub] = problem.GetVariableBounds(i);
        spec.AddContinuous(lb, ub, "x" + std::to_string(i));
    }
    
    NSGA2Config config;
    config.populationSize = 150;  // Larger population for multimodal
    config.maxGenerations = 400;  // More generations needed
    config.seed = 42;
    
    NSGA2<10, 2> nsga(config);
    nsga.SetProblem(spec);
    
    auto result = nsga.Optimize(problem);
    
    VectorN<Real, 2> refPoint{1.1, 1.1};
    PrintNSGA2Result("ZDT4 Result", result, refPoint);
    PrintParetoFrontSample(result);
    
    // Check if we found the global front (g should be close to 1)
    std::cout << "\nGlobal Pareto front has f2 = 1 - sqrt(f1) (same as ZDT1)" << std::endl;
    std::cout << "Local fronts have higher f2 values (dominated)" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
// Example 5: ZDT6 - Non-uniform Distribution
///////////////////////////////////////////////////////////////////////////////////////////
void Example_ZDT6_NonUniform()
{
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Example 5: ZDT6 - Non-Uniform Pareto Front Distribution" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    std::cout << "\nZDT6 has two challenges:" << std::endl;
    std::cout << "  - Non-uniform distribution (solutions cluster near f1=1)" << std::endl;
    std::cout << "  - Low density of solutions near the true Pareto front" << std::endl;
    std::cout << "  - f1 range is approximately [0.28, 1.0]" << std::endl;
    std::cout << std::endl;
    
    ZDT6<10> problem;
    
    GAProblemSpec spec;
    for (int i = 0; i < 10; ++i) {
        auto [lb, ub] = problem.GetVariableBounds(i);
        spec.AddContinuous(lb, ub, "x" + std::to_string(i));
    }
    
    NSGA2Config config;
    config.populationSize = 100;
    config.maxGenerations = 300;
    config.seed = 42;
    
    NSGA2<10, 2> nsga(config);
    nsga.SetProblem(spec);
    
    auto result = nsga.Optimize(problem);
    
    VectorN<Real, 2> refPoint{1.1, 1.1};
    PrintNSGA2Result("ZDT6 Result", result, refPoint);
    PrintParetoFrontSample(result);
    
    std::cout << "\nNote: f1 starts around 0.28 (not 0) due to the exponential term" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
// Example 6: Schaffer's SCH1 - Simple 1D Problem
///////////////////////////////////////////////////////////////////////////////////////////
void Example_SCH1_Simple()
{
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Example 6: Schaffer's SCH1 - Simple 1D Bi-objective" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    std::cout << "\nSCH1 is the simplest bi-objective test problem:" << std::endl;
    std::cout << "  - Single decision variable x" << std::endl;
    std::cout << "  - f1(x) = x^2" << std::endl;
    std::cout << "  - f2(x) = (x-2)^2" << std::endl;
    std::cout << "  - Pareto front: x in [0, 2]" << std::endl;
    std::cout << std::endl;
    
    SCH1 problem;
    
    GAProblemSpec spec;
    auto [lb, ub] = problem.GetVariableBounds(0);
    spec.AddContinuous(lb, ub, "x");
    
    NSGA2Config config;
    config.populationSize = 50;
    config.maxGenerations = 100;
    config.seed = 42;
    
    NSGA2<1, 2> nsga(config);
    nsga.SetProblem(spec);
    
    auto result = nsga.Optimize(problem);
    
    VectorN<Real, 2> refPoint{5.0, 5.0};
    PrintNSGA2Result("SCH1 Result", result, refPoint);
    
    std::cout << "\nSample solutions with decision variable x:" << std::endl;
    std::cout << std::string(45, '-') << std::endl;
    
    // Show x values along with objectives
    auto& front = result.paretoFront;
    std::vector<const ObjectivePoint<1, 2>*> sorted;
    for (const auto& pt : front) sorted.push_back(&pt);
    std::sort(sorted.begin(), sorted.end(),
        [](const auto* a, const auto* b) {
            return a->GetObjective(0) < b->GetObjective(0);
        });
    
    int step = std::max(1, static_cast<int>(sorted.size()) / 8);
    for (size_t i = 0; i < sorted.size(); i += step) {
        Real x = sorted[i]->GetX()[0];
        std::cout << "  x=" << std::setw(6) << std::setprecision(3) << x
                  << " -> f1=" << std::setw(6) << sorted[i]->GetObjective(0)
                  << ", f2=" << std::setw(6) << sorted[i]->GetObjective(1) << std::endl;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// Example 7: Fonseca-Fleming (FON) - Concave Front
///////////////////////////////////////////////////////////////////////////////////////////
void Example_FON_Concave()
{
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Example 7: Fonseca-Fleming (FON) - Concave Pareto Front" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    std::cout << "\nFON problem with strongly concave Pareto front:" << std::endl;
    std::cout << "  - N=3 decision variables" << std::endl;
    std::cout << "  - f1 = 1 - exp(-sum((x_i - 1/sqrt(n))^2))" << std::endl;
    std::cout << "  - f2 = 1 - exp(-sum((x_i + 1/sqrt(n))^2))" << std::endl;
    std::cout << "  - Pareto optimal: all x_i equal" << std::endl;
    std::cout << std::endl;
    
    FON<3> problem;
    
    GAProblemSpec spec;
    for (int i = 0; i < 3; ++i) {
        auto [lb, ub] = problem.GetVariableBounds(i);
        spec.AddContinuous(lb, ub, "x" + std::to_string(i));
    }
    
    NSGA2Config config;
    config.populationSize = 100;
    config.maxGenerations = 200;
    config.seed = 42;
    
    NSGA2<3, 2> nsga(config);
    nsga.SetProblem(spec);
    
    auto result = nsga.Optimize(problem);
    
    VectorN<Real, 2> refPoint{1.0, 1.0};
    PrintNSGA2Result("FON Result", result, refPoint);
    PrintParetoFrontSample(result);
}

///////////////////////////////////////////////////////////////////////////////////////////
// Example 8: Kursawe (KUR) - Non-convex Discontinuous
///////////////////////////////////////////////////////////////////////////////////////////
void Example_KUR_Complex()
{
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Example 8: Kursawe (KUR) - Complex Non-Convex Front" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    std::cout << "\nKursawe problem with complex Pareto front shape:" << std::endl;
    std::cout << "  - N=3 decision variables" << std::endl;
    std::cout << "  - Non-convex, possibly discontinuous front" << std::endl;
    std::cout << "  - f1 involves exponential and square root terms" << std::endl;
    std::cout << "  - f2 involves sine function (non-convex)" << std::endl;
    std::cout << std::endl;
    
    KUR<3> problem;
    
    GAProblemSpec spec;
    for (int i = 0; i < 3; ++i) {
        auto [lb, ub] = problem.GetVariableBounds(i);
        spec.AddContinuous(lb, ub, "x" + std::to_string(i));
    }
    
    NSGA2Config config;
    config.populationSize = 100;
    config.maxGenerations = 250;
    config.seed = 42;
    
    NSGA2<3, 2> nsga(config);
    nsga.SetProblem(spec);
    
    auto result = nsga.Optimize(problem);
    
    VectorN<Real, 2> refPoint{0.0, 20.0};  // KUR has f1 < 0
    PrintNSGA2Result("KUR Result", result, refPoint);
    PrintParetoFrontSample(result);
    
    std::cout << "\nNote: KUR typically has f1 in [-20, -14] and f2 in [-12, 2]" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
// Example 9: Configuration Tuning - Effect of Population Size
///////////////////////////////////////////////////////////////////////////////////////////
void Example_NSGA2_ConfigTuning()
{
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Example 9: Configuration Tuning - Population Size Effect" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    std::cout << "\nStudying effect of population size on ZDT1 convergence:" << std::endl;
    std::cout << std::endl;
    
    ZDT1<30> problem;
    
    GAProblemSpec spec;
    for (int i = 0; i < 30; ++i) {
        auto [lb, ub] = problem.GetVariableBounds(i);
        spec.AddContinuous(lb, ub, "x" + std::to_string(i));
    }
    
    VectorN<Real, 2> refPoint{1.1, 1.1};
    
    std::cout << std::setw(12) << "PopSize" << std::setw(15) << "Front Size" 
              << std::setw(15) << "Hypervolume" << std::setw(15) << "Func Evals" << std::endl;
    std::cout << std::string(57, '-') << std::endl;
    
    for (int popSize : {50, 100, 150, 200}) {
        NSGA2Config config;
        config.populationSize = popSize;
        config.maxGenerations = 200;
        config.seed = 42;
        
        NSGA2<30, 2> nsga(config);
        nsga.SetProblem(spec);
        
        auto result = nsga.Optimize(problem);
        Real hv = result.GetHypervolume2D(refPoint);
        
        std::cout << std::setw(12) << popSize 
                  << std::setw(15) << result.paretoFront.size()
                  << std::setw(15) << std::fixed << std::setprecision(6) << hv
                  << std::setw(15) << result.functionEvaluations << std::endl;
    }
    
    std::cout << "\nLarger population = better coverage but more function evaluations" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
// Example 10: Effect of Generation Count
///////////////////////////////////////////////////////////////////////////////////////////
void Example_NSGA2_GenerationEffect()
{
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Example 10: Convergence Analysis - Generation Count Effect" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    std::cout << "\nTracking convergence of hypervolume over generations on ZDT1:" << std::endl;
    std::cout << std::endl;
    
    ZDT1<30> problem;
    
    GAProblemSpec spec;
    for (int i = 0; i < 30; ++i) {
        auto [lb, ub] = problem.GetVariableBounds(i);
        spec.AddContinuous(lb, ub, "x" + std::to_string(i));
    }
    
    VectorN<Real, 2> refPoint{1.1, 1.1};
    // Theoretical hypervolume for perfect ZDT1 front with ref (1.1, 1.1) ≈ 0.853
    
    std::cout << std::setw(12) << "Generations" << std::setw(15) << "Hypervolume" 
              << std::setw(15) << "% of Optimal" << std::endl;
    std::cout << std::string(42, '-') << std::endl;
    
    const Real optimalHV = 0.853;  // Approximate theoretical optimal
    
    for (int gens : {25, 50, 100, 150, 200, 250, 300}) {
        NSGA2Config config;
        config.populationSize = 100;
        config.maxGenerations = gens;
        config.seed = 42;
        
        NSGA2<30, 2> nsga(config);
        nsga.SetProblem(spec);
        
        auto result = nsga.Optimize(problem);
        Real hv = result.GetHypervolume2D(refPoint);
        Real pctOptimal = (hv / optimalHV) * 100;
        
        std::cout << std::setw(12) << gens 
                  << std::setw(15) << std::fixed << std::setprecision(6) << hv
                  << std::setw(14) << std::setprecision(1) << pctOptimal << "%" << std::endl;
    }
    
    std::cout << "\nConvergence typically plateaus after 200-250 generations" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
// Example 11: Operator Parameter Tuning
///////////////////////////////////////////////////////////////////////////////////////////
void Example_NSGA2_OperatorTuning()
{
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Example 11: Operator Parameter Tuning" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    std::cout << "\nEffect of SBX crossover eta (distribution index) on ZDT1:" << std::endl;
    std::cout << "  - Lower eta = more exploration (children far from parents)" << std::endl;
    std::cout << "  - Higher eta = more exploitation (children near parents)" << std::endl;
    std::cout << std::endl;
    
    ZDT1<30> problem;
    
    GAProblemSpec spec;
    for (int i = 0; i < 30; ++i) {
        auto [lb, ub] = problem.GetVariableBounds(i);
        spec.AddContinuous(lb, ub, "x" + std::to_string(i));
    }
    
    VectorN<Real, 2> refPoint{1.1, 1.1};
    
    std::cout << "--- Crossover eta (mutation eta = 20) ---" << std::endl;
    std::cout << std::setw(10) << "Eta" << std::setw(15) << "Hypervolume" << std::endl;
    std::cout << std::string(25, '-') << std::endl;
    
    for (Real eta : {2.0, 5.0, 10.0, 15.0, 20.0, 30.0}) {
        NSGA2Config config;
        config.populationSize = 100;
        config.maxGenerations = 200;
        config.crossoverEta = eta;
        config.mutationEta = 20.0;
        config.seed = 42;
        
        NSGA2<30, 2> nsga(config);
        nsga.SetProblem(spec);
        
        auto result = nsga.Optimize(problem);
        Real hv = result.GetHypervolume2D(refPoint);
        
        std::cout << std::setw(10) << std::fixed << std::setprecision(1) << eta 
                  << std::setw(15) << std::setprecision(6) << hv << std::endl;
    }
    
    std::cout << "\n--- Mutation eta (crossover eta = 15) ---" << std::endl;
    std::cout << std::setw(10) << "Eta" << std::setw(15) << "Hypervolume" << std::endl;
    std::cout << std::string(25, '-') << std::endl;
    
    for (Real eta : {5.0, 10.0, 20.0, 30.0, 50.0}) {
        NSGA2Config config;
        config.populationSize = 100;
        config.maxGenerations = 200;
        config.crossoverEta = 15.0;
        config.mutationEta = eta;
        config.seed = 42;
        
        NSGA2<30, 2> nsga(config);
        nsga.SetProblem(spec);
        
        auto result = nsga.Optimize(problem);
        Real hv = result.GetHypervolume2D(refPoint);
        
        std::cout << std::setw(10) << std::fixed << std::setprecision(1) << eta 
                  << std::setw(15) << std::setprecision(6) << hv << std::endl;
    }
    
    std::cout << "\nStandard values: crossoverEta=15, mutationEta=20 (Deb et al.)" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
// Example 12: Practical Application - Weight vs Cost Optimization
///////////////////////////////////////////////////////////////////////////////////////////
void Example_NSGA2_BeamDesign()
{
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Example 12: Practical Application - Beam Design Trade-off" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    std::cout << "\nBi-objective beam design:" << std::endl;
    std::cout << "  - Objective 1: Minimize weight (ρ * b * h * L)" << std::endl;
    std::cout << "  - Objective 2: Minimize deflection (P * L^3 / (E * b * h^3))" << std::endl;
    std::cout << "  - Variables: width b [10, 200] mm, height h [20, 300] mm" << std::endl;
    std::cout << "  - Constraint: Stress <= 250 MPa (handled via penalty)" << std::endl;
    std::cout << std::endl;
    
    // Problem parameters
    const Real L = 1.0;           // Length (m)
    const Real P = 10000.0;       // Load (N)
    const Real E = 200e9;         // Young's modulus (Pa) - steel
    const Real rho = 7850.0;      // Density (kg/m^3) - steel
    const Real sigma_max = 250e6; // Max stress (Pa)
    
    // Custom bi-objective problem class
    class BeamDesignProblem
    {
    public:
        Real L_, P_, E_, rho_, sigma_max_;
        
        BeamDesignProblem(Real L, Real P, Real E, Real rho, Real sigma_max)
            : L_(L), P_(P), E_(E), rho_(rho), sigma_max_(sigma_max) {}
        
        VectorN<Real, 2> Evaluate(const VectorN<Real, 2>& x) const
        {
            Real b = x[0] / 1000.0;  // Convert mm to m
            Real h = x[1] / 1000.0;
            
            // Weight objective
            Real weight = rho_ * b * h * L_;
            
            // Deflection objective (scaled for similar magnitude)
            Real deflection = 4.0 * P_ * L_ * L_ * L_ / (E_ * b * h * h * h);
            Real deflection_mm = deflection * 1000;  // Convert to mm
            
            // Stress constraint penalty
            Real sigma = 6.0 * P_ * L_ / (b * h * h);
            Real penalty = 0;
            if (sigma > sigma_max_) {
                penalty = 1000 * (sigma - sigma_max_) / sigma_max_;
            }
            
            return VectorN<Real, 2>{weight + penalty, deflection_mm + penalty};
        }
    };
    
    BeamDesignProblem problem(L, P, E, rho, sigma_max);
    
    // Variable bounds: b [10, 200] mm, h [20, 300] mm
    GAProblemSpec spec;
    spec.AddContinuous(10.0, 200.0, "width_mm");
    spec.AddContinuous(20.0, 300.0, "height_mm");
    
    NSGA2Config config;
    config.populationSize = 100;
    config.maxGenerations = 200;
    config.seed = 42;
    
    NSGA2<2, 2> nsga(config);
    nsga.SetProblem(spec);
    
    auto result = nsga.Optimize(problem);
    
    std::cout << "=== Beam Design Pareto Front ===" << std::endl;
    std::cout << "Pareto front size: " << result.paretoFront.size() << std::endl;
    std::cout << "\nSample trade-off solutions:" << std::endl;
    std::cout << std::string(60, '-') << std::endl;
    std::cout << std::setw(12) << "Width(mm)" << std::setw(12) << "Height(mm)"
              << std::setw(12) << "Weight(kg)" << std::setw(15) << "Deflect(mm)" << std::endl;
    std::cout << std::string(60, '-') << std::endl;
    
    // Sort by weight
    auto& front = result.paretoFront;
    std::vector<const ObjectivePoint<2, 2>*> sorted;
    for (const auto& pt : front) sorted.push_back(&pt);
    std::sort(sorted.begin(), sorted.end(),
        [](const auto* a, const auto* b) {
            return a->GetObjective(0) < b->GetObjective(0);
        });
    
    int step = std::max(1, static_cast<int>(sorted.size()) / 8);
    for (size_t i = 0; i < sorted.size(); i += step) {
        Real b_mm = sorted[i]->GetX()[0];
        Real h_mm = sorted[i]->GetX()[1];
        Real weight = sorted[i]->GetObjective(0);
        Real deflect = sorted[i]->GetObjective(1);
        
        std::cout << std::fixed << std::setprecision(2)
                  << std::setw(12) << b_mm << std::setw(12) << h_mm
                  << std::setw(12) << weight << std::setw(15) << deflect << std::endl;
    }
    
    std::cout << "\nThis Pareto front shows the trade-off between:" << std::endl;
    std::cout << "  - Light beams with high deflection" << std::endl;
    std::cout << "  - Heavy beams with low deflection" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
// Main - Run all demos
///////////////////////////////////////////////////////////////////////////////////////////
void Docs_Demo_NSGA2()
{
    std::cout << "\n";
    std::cout << "######################################################################" << std::endl;
    std::cout << "###                                                                ###" << std::endl;
    std::cout << "###        NSGA-II MULTI-OBJECTIVE OPTIMIZATION DEMONSTRATION      ###" << std::endl;
    std::cout << "###                                                                ###" << std::endl;
    std::cout << "###   Non-dominated Sorting Genetic Algorithm II (Deb et al. 2002) ###" << std::endl;
    std::cout << "###   Features: Pareto ranking, crowding distance, elitism         ###" << std::endl;
    std::cout << "###                                                                ###" << std::endl;
    std::cout << "######################################################################" << std::endl;
    
    Example_ZDT1_Basic();
    Example_ZDT2_Concave();
    Example_ZDT3_Disconnected();
    Example_ZDT4_Multimodal();
    Example_ZDT6_NonUniform();
    Example_SCH1_Simple();
    Example_FON_Concave();
    Example_KUR_Complex();
    Example_NSGA2_ConfigTuning();
    Example_NSGA2_GenerationEffect();
    Example_NSGA2_OperatorTuning();
    Example_NSGA2_BeamDesign();
    
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "NSGA-II demonstration complete!" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
}

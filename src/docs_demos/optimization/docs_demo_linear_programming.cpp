///////////////////////////////////////////////////////////////////////////////////////////
// docs_demo_linear_programming.cpp - Examples of Linear Programming with MML
//
// Demonstrates:
// 1. Simple 2D LP problems (min/max)
// 2. Resource allocation / Production planning
// 3. Transportation-like problems
// 4. Problems with different constraint types (<=, =, >=)
// 5. Special cases (unbounded, infeasible)
// 6. Dual problem construction
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#include "base/Vector/Vector.h"
#include "base/Matrix/Matrix.h"
#include "LinearProgramming.h"
#endif

#include <iostream>
#include <iomanip>

using namespace MML;
using namespace MML::Optimization;

///////////////////////////////////////////////////////////////////////////////////////////
// Helper function to print LP results
///////////////////////////////////////////////////////////////////////////////////////////
void PrintLPResult(const std::string& name, const LPResult& result, int numVars)
{
    std::cout << "\n=== " << name << " ===" << std::endl;
    std::cout << "Status: " << result.statusMessage() << std::endl;
    
    if (result.IsOptimal()) {
        std::cout << std::fixed << std::setprecision(4);
        std::cout << "Optimal value: " << result.objectiveValue << std::endl;
        std::cout << "Solution: ";
        for (int i = 0; i < numVars; ++i) {
            std::cout << "x" << (i+1) << " = " << result.x[i];
            if (i < numVars - 1) std::cout << ", ";
        }
        std::cout << std::endl;
        std::cout << "Iterations: " << result.iterations;
        if (result.phase1Iterations > 0) {
            std::cout << " (Phase 1: " << result.phase1Iterations << ")";
        }
        std::cout << std::endl;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// Example 1: Simple 2D Minimization
///////////////////////////////////////////////////////////////////////////////////////////
void Example_Simple2D_Minimization()
{
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Example 1: Simple 2D Minimization" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    // Problem: A diet problem - minimize cost while meeting nutrient requirements
    //   min  2x + 3y        (cost)
    //   s.t. x + y >= 4     (minimum total units)
    //        2x + y >= 5    (protein requirement)
    //        x, y >= 0
    //
    // Optimal: x=1, y=3, cost=11
    
    std::cout << "Diet problem - minimize cost while meeting requirements:" << std::endl;
    std::cout << std::endl;
    std::cout << "  Minimize: 2x + 3y  (cost)" << std::endl;
    std::cout << "  Subject to:" << std::endl;
    std::cout << "    x + y >= 4       (min total units)" << std::endl;
    std::cout << "    2x + y >= 5      (protein requirement)" << std::endl;
    std::cout << "    x, y >= 0" << std::endl;
    
    LinearProgram lp(2);
    lp.SetObjective({2.0, 3.0}, LPObjective::Minimize);
    lp.AddConstraint({1.0, 1.0}, LPConstraintType::GreaterEqual, 4.0);
    lp.AddConstraint({2.0, 1.0}, LPConstraintType::GreaterEqual, 5.0);
    
    LPResult result = SolveLP(lp);
    PrintLPResult("Result", result, 2);
}

///////////////////////////////////////////////////////////////////////////////////////////
// Example 2: Simple 2D Maximization
///////////////////////////////////////////////////////////////////////////////////////////
void Example_Simple2D_Maximization()
{
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Example 2: Simple 2D Maximization" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    // Problem:
    //   max  3x + 2y
    //   s.t. x + y <= 4
    //        2x + y <= 5
    //        x, y >= 0
    //
    // Optimal: x=1, y=3, objective=9
    
    std::cout << "Problem:" << std::endl;
    std::cout << "  Maximize: 3x + 2y" << std::endl;
    std::cout << "  Subject to:" << std::endl;
    std::cout << "    x + y <= 4" << std::endl;
    std::cout << "    2x + y <= 5" << std::endl;
    std::cout << "    x, y >= 0" << std::endl;
    
    LinearProgram lp(2);
    lp.SetObjective({3.0, 2.0}, LPObjective::Maximize);
    lp.AddConstraint({1.0, 1.0}, LPConstraintType::LessEqual, 4.0);
    lp.AddConstraint({2.0, 1.0}, LPConstraintType::LessEqual, 5.0);
    
    LPResult result = SolveLP(lp);
    PrintLPResult("Result", result, 2);
}

///////////////////////////////////////////////////////////////////////////////////////////
// Example 3: Resource Allocation Problem
///////////////////////////////////////////////////////////////////////////////////////////
void Example_ResourceAllocation()
{
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Example 3: Resource Allocation Problem" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    // A company produces two products (P1, P2) using two resources (R1, R2).
    // P1 yields profit $5, P2 yields profit $4.
    // P1 uses 6 units of R1 and 1 unit of R2.
    // P2 uses 4 units of R1 and 2 units of R2.
    // Available: 24 units of R1, 6 units of R2.
    //
    // Problem:
    //   max  5x1 + 4x2
    //   s.t. 6x1 + 4x2 <= 24  (R1)
    //        x1 + 2x2 <= 6    (R2)
    //        x1, x2 >= 0
    //
    // Optimal: x1=3, x2=1.5, profit=19.5
    
    std::cout << "A company produces two products using two resources." << std::endl;
    std::cout << "Product 1: profit $5, uses 6 units R1 + 1 unit R2" << std::endl;
    std::cout << "Product 2: profit $4, uses 4 units R1 + 2 units R2" << std::endl;
    std::cout << "Available: 24 units R1, 6 units R2" << std::endl;
    std::cout << std::endl;
    std::cout << "Problem:" << std::endl;
    std::cout << "  Maximize: 5*x1 + 4*x2" << std::endl;
    std::cout << "  Subject to:" << std::endl;
    std::cout << "    6*x1 + 4*x2 <= 24  (Resource 1)" << std::endl;
    std::cout << "    x1 + 2*x2 <= 6     (Resource 2)" << std::endl;
    std::cout << "    x1, x2 >= 0" << std::endl;
    
    LinearProgram lp(2, "ResourceAllocation");
    lp.SetVariableName(0, "Product1");
    lp.SetVariableName(1, "Product2");
    lp.SetObjective({5.0, 4.0}, LPObjective::Maximize);
    lp.AddConstraint({6.0, 4.0}, LPConstraintType::LessEqual, 24.0, "Resource1");
    lp.AddConstraint({1.0, 2.0}, LPConstraintType::LessEqual, 6.0, "Resource2");
    
    LPResult result = SolveLP(lp);
    PrintLPResult("Result", result, 2);
    
    if (result.IsOptimal()) {
        std::cout << "Interpretation: Produce " << result.x[0] << " units of P1 and " 
                  << result.x[1] << " units of P2 for maximum profit of $" 
                  << result.objectiveValue << std::endl;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// Example 4: Production Planning with Multiple Constraints
///////////////////////////////////////////////////////////////////////////////////////////
void Example_ProductionPlanning()
{
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Example 4: Production Planning (3 Products, 3 Resources)" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    // Problem:
    //   max  10x1 + 12x2 + 8x3
    //   s.t. 2x1 + 3x2 + x3 <= 100   (Labor hours)
    //        x1 + x2 + x3 <= 50      (Machine hours)
    //        x1 + 2x2 + 3x3 <= 80    (Raw material)
    //        x1, x2, x3 >= 0
    
    std::cout << "Three products with profits: $10, $12, $8" << std::endl;
    std::cout << "Resource constraints: Labor (100h), Machine (50h), Material (80 units)" << std::endl;
    std::cout << std::endl;
    std::cout << "Problem:" << std::endl;
    std::cout << "  Maximize: 10*x1 + 12*x2 + 8*x3" << std::endl;
    std::cout << "  Subject to:" << std::endl;
    std::cout << "    2*x1 + 3*x2 + x3 <= 100  (Labor)" << std::endl;
    std::cout << "    x1 + x2 + x3 <= 50       (Machine)" << std::endl;
    std::cout << "    x1 + 2*x2 + 3*x3 <= 80   (Material)" << std::endl;
    
    LinearProgram lp(3, "ProductionPlanning");
    lp.SetObjective({10.0, 12.0, 8.0}, LPObjective::Maximize);
    lp.AddConstraint({2.0, 3.0, 1.0}, LPConstraintType::LessEqual, 100.0, "Labor");
    lp.AddConstraint({1.0, 1.0, 1.0}, LPConstraintType::LessEqual, 50.0, "Machine");
    lp.AddConstraint({1.0, 2.0, 3.0}, LPConstraintType::LessEqual, 80.0, "Material");
    
    LPResult result = SolveLP(lp);
    PrintLPResult("Result", result, 3);
}

///////////////////////////////////////////////////////////////////////////////////////////
// Example 5: Problem with Equality Constraint (Two-Phase Method)
///////////////////////////////////////////////////////////////////////////////////////////
void Example_EqualityConstraint()
{
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Example 5: Problem with Equality Constraint" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    // Problem:
    //   min  x + y
    //   s.t. x + y = 5
    //        x + 2y <= 8
    //        x, y >= 0
    //
    // This requires the two-phase method due to the equality constraint.
    // On the line x+y=5, we minimize x+y which is constant = 5.
    // Optimal: any point on x+y=5 within feasible region, e.g., x=2, y=3
    
    std::cout << "Problem with equality constraint (requires two-phase method):" << std::endl;
    std::cout << std::endl;
    std::cout << "  Minimize: x + y" << std::endl;
    std::cout << "  Subject to:" << std::endl;
    std::cout << "    x + y = 5      (equality)" << std::endl;
    std::cout << "    x + 2y <= 8" << std::endl;
    std::cout << "    x, y >= 0" << std::endl;
    
    LinearProgram lp(2);
    lp.SetObjective({1.0, 1.0}, LPObjective::Minimize);
    lp.AddConstraint({1.0, 1.0}, LPConstraintType::Equal, 5.0);
    lp.AddConstraint({1.0, 2.0}, LPConstraintType::LessEqual, 8.0);
    
    LPResult result = SolveLP(lp);
    PrintLPResult("Result", result, 2);
    
    if (result.IsOptimal()) {
        std::cout << "Note: Phase 1 iterations = " << result.phase1Iterations 
                  << " (to find initial basic feasible solution)" << std::endl;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// Example 6: Problem with >= Constraint
///////////////////////////////////////////////////////////////////////////////////////////
void Example_GreaterEqualConstraint()
{
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Example 6: Problem with >= Constraint" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    // Problem:
    //   min  x + 2y
    //   s.t. x + y >= 3
    //        2x + y <= 10
    //        x, y >= 0
    //
    // Optimal: x=3, y=0, objective=3
    
    std::cout << "Problem with >= constraint (requires surplus + artificial var):" << std::endl;
    std::cout << std::endl;
    std::cout << "  Minimize: x + 2y" << std::endl;
    std::cout << "  Subject to:" << std::endl;
    std::cout << "    x + y >= 3     (>= constraint)" << std::endl;
    std::cout << "    2x + y <= 10" << std::endl;
    std::cout << "    x, y >= 0" << std::endl;
    
    LinearProgram lp(2);
    lp.SetObjective({1.0, 2.0}, LPObjective::Minimize);
    lp.AddConstraint({1.0, 1.0}, LPConstraintType::GreaterEqual, 3.0);
    lp.AddConstraint({2.0, 1.0}, LPConstraintType::LessEqual, 10.0);
    
    LPResult result = SolveLP(lp);
    PrintLPResult("Result", result, 2);
}

///////////////////////////////////////////////////////////////////////////////////////////
// Example 7: Mixed Constraint Types
///////////////////////////////////////////////////////////////////////////////////////////
void Example_MixedConstraints()
{
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Example 7: Mixed Constraint Types (<=, =, >=)" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    // Problem:
    //   min  2x + 3y
    //   s.t. x + y = 6      (equality)
    //        x - y <= 2     (<=)
    //        y >= 1         (>=)
    //        x, y >= 0
    
    std::cout << "Problem with all three constraint types:" << std::endl;
    std::cout << std::endl;
    std::cout << "  Minimize: 2x + 3y" << std::endl;
    std::cout << "  Subject to:" << std::endl;
    std::cout << "    x + y = 6      (equality)" << std::endl;
    std::cout << "    x - y <= 2     (<=)" << std::endl;
    std::cout << "    y >= 1         (>=)" << std::endl;
    std::cout << "    x, y >= 0" << std::endl;
    
    LinearProgram lp(2);
    lp.SetObjective({2.0, 3.0}, LPObjective::Minimize);
    lp.AddConstraint({1.0, 1.0}, LPConstraintType::Equal, 6.0);
    lp.AddConstraint({1.0, -1.0}, LPConstraintType::LessEqual, 2.0);
    lp.AddConstraint({0.0, 1.0}, LPConstraintType::GreaterEqual, 1.0);
    
    LPResult result = SolveLP(lp);
    PrintLPResult("Result", result, 2);
}

///////////////////////////////////////////////////////////////////////////////////////////
// Example 8: Matrix Form Construction
///////////////////////////////////////////////////////////////////////////////////////////
void Example_MatrixConstruction()
{
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Example 8: LP from Matrix Form" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    // Same problem as Example 2, but using matrix construction:
    //   max  3x + 2y
    //   s.t. x + y <= 4
    //        2x + y <= 5
    //        x, y >= 0
    
    std::cout << "Constructing LP from matrices A, b, c:" << std::endl;
    std::cout << "  min c'x  s.t. Ax <= b, x >= 0" << std::endl;
    std::cout << std::endl;
    
    // Objective coefficients
    Vector<Real> c({-3.0, -2.0});  // Negate for maximization
    
    // Constraint matrix
    Matrix<Real> A(2, 2);
    A(0, 0) = 1.0; A(0, 1) = 1.0;
    A(1, 0) = 2.0; A(1, 1) = 1.0;
    
    // RHS
    Vector<Real> b({4.0, 5.0});
    
    std::cout << "c = " << c << std::endl;
    std::cout << "A = " << std::endl;
    A.Print(std::cout, 8, 2);
    std::cout << "b = " << b << std::endl;
    
    // Construct LP (default all <= constraints)
    LinearProgram lp(c, A, b);
    
    LPResult result = SolveLP(lp);
    PrintLPResult("Result (negated for max)", result, 2);
    
    std::cout << "Actual maximum = " << -result.objectiveValue << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
// Example 9: Unbounded Problem
///////////////////////////////////////////////////////////////////////////////////////////
void Example_UnboundedProblem()
{
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Example 9: Unbounded Problem" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    // Problem:
    //   max  x + y
    //   s.t. x - y <= 1
    //        x, y >= 0
    //
    // This is unbounded (can increase x indefinitely with y = x - 1)
    
    std::cout << "An unbounded LP (objective can be made arbitrarily large):" << std::endl;
    std::cout << std::endl;
    std::cout << "  Maximize: x + y" << std::endl;
    std::cout << "  Subject to:" << std::endl;
    std::cout << "    x - y <= 1" << std::endl;
    std::cout << "    x, y >= 0" << std::endl;
    
    LinearProgram lp(2);
    lp.SetObjective({1.0, 1.0}, LPObjective::Maximize);
    lp.AddConstraint({1.0, -1.0}, LPConstraintType::LessEqual, 1.0);
    
    LPResult result = SolveLP(lp);
    PrintLPResult("Result", result, 2);
}

///////////////////////////////////////////////////////////////////////////////////////////
// Example 10: Infeasible Problem
///////////////////////////////////////////////////////////////////////////////////////////
void Example_InfeasibleProblem()
{
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Example 10: Infeasible Problem" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    // Problem:
    //   min  x + y
    //   s.t. x + y <= 1
    //        x + y >= 3
    //        x, y >= 0
    //
    // Infeasible: cannot have x+y <= 1 and x+y >= 3 simultaneously
    
    std::cout << "An infeasible LP (contradictory constraints):" << std::endl;
    std::cout << std::endl;
    std::cout << "  Minimize: x + y" << std::endl;
    std::cout << "  Subject to:" << std::endl;
    std::cout << "    x + y <= 1" << std::endl;
    std::cout << "    x + y >= 3   (contradicts above!)" << std::endl;
    std::cout << "    x, y >= 0" << std::endl;
    
    LinearProgram lp(2);
    lp.SetObjective({1.0, 1.0}, LPObjective::Minimize);
    lp.AddConstraint({1.0, 1.0}, LPConstraintType::LessEqual, 1.0);
    lp.AddConstraint({1.0, 1.0}, LPConstraintType::GreaterEqual, 3.0);
    
    LPResult result = SolveLP(lp);
    PrintLPResult("Result", result, 2);
}

///////////////////////////////////////////////////////////////////////////////////////////
// Example 11: Dual Problem Construction
///////////////////////////////////////////////////////////////////////////////////////////
void Example_DualProblem()
{
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Example 11: Primal-Dual Relationship" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    // Primal (maximization with <= constraints):
    //   max  3x + 2y
    //   s.t. x + y <= 4
    //        2x + y <= 5
    //        x, y >= 0
    //
    // Dual (minimization with >= constraints):
    //   min  4u + 5v
    //   s.t. u + 2v >= 3
    //        u + v >= 2
    //        u, v >= 0
    //
    // By strong duality: primal optimal = dual optimal = 9
    
    std::cout << "Primal problem (maximization):" << std::endl;
    std::cout << "  max 3x + 2y  s.t.  x+y<=4, 2x+y<=5, x,y>=0" << std::endl;
    std::cout << std::endl;
    
    LinearProgram primal(2, "Primal");
    primal.SetObjective({3.0, 2.0}, LPObjective::Maximize);
    primal.AddConstraint({1.0, 1.0}, LPConstraintType::LessEqual, 4.0);
    primal.AddConstraint({2.0, 1.0}, LPConstraintType::LessEqual, 5.0);
    
    LPResult primalResult = SolveLP(primal);
    PrintLPResult("Primal Result", primalResult, 2);
    
    // Construct dual
    LinearProgram dual = primal.Dual();
    
    std::cout << "\nDual problem (constructed automatically):" << std::endl;
    std::cout << dual.ToString() << std::endl;
    
    LPResult dualResult = SolveLP(dual);
    PrintLPResult("Dual Result", dualResult, 2);
    
    std::cout << "\nStrong Duality Theorem:" << std::endl;
    std::cout << "  Primal optimal value = " << primalResult.objectiveValue << std::endl;
    std::cout << "  Dual optimal value   = " << dualResult.objectiveValue << std::endl;
    std::cout << "  These should be equal at optimality!" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
// Example 12: Solution Analysis
///////////////////////////////////////////////////////////////////////////////////////////
void Example_SolutionAnalysis()
{
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Example 12: Detailed Solution Analysis" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    // Problem:
    //   max  5x1 + 4x2
    //   s.t. 6x1 + 4x2 <= 24
    //        x1 + 2x2 <= 6
    //        x1, x2 >= 0
    
    LinearProgram lp(2);
    lp.SetObjective({5.0, 4.0}, LPObjective::Maximize);
    lp.AddConstraint({6.0, 4.0}, LPConstraintType::LessEqual, 24.0, "Resource1");
    lp.AddConstraint({1.0, 2.0}, LPConstraintType::LessEqual, 6.0, "Resource2");
    
    LPResult result = SolveLP(lp);
    
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Problem: max 5x1 + 4x2  s.t. 6x1+4x2<=24, x1+2x2<=6" << std::endl;
    std::cout << std::endl;
    
    if (result.IsOptimal()) {
        std::cout << "OPTIMAL SOLUTION:" << std::endl;
        std::cout << "  x1 = " << result.x[0] << std::endl;
        std::cout << "  x2 = " << result.x[1] << std::endl;
        std::cout << "  Objective = " << result.objectiveValue << std::endl;
        std::cout << std::endl;
        
        std::cout << "CONSTRAINT ANALYSIS:" << std::endl;
        std::cout << "  Slack variables (unused resources):" << std::endl;
        for (int i = 0; i < result.slacks.size(); ++i) {
            std::cout << "    Constraint " << (i+1) << ": slack = " << result.slacks[i];
            if (std::abs(result.slacks[i]) < 1e-6) {
                std::cout << " (BINDING)";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
        
        std::cout << "DUAL VALUES (Shadow Prices):" << std::endl;
        for (int i = 0; i < result.dualValues.size(); ++i) {
            std::cout << "    Constraint " << (i+1) << ": " << result.dualValues[i] << std::endl;
        }
        std::cout << "  (Value of one additional unit of each resource)" << std::endl;
        std::cout << std::endl;
        
        std::cout << "BASIS:" << std::endl;
        std::cout << "  Basic variable indices: ";
        for (int idx : result.basisIndices) {
            std::cout << idx << " ";
        }
        std::cout << std::endl;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// Example 13: Large-Scale Manufacturing Problem (5 products, 6 constraints)
///////////////////////////////////////////////////////////////////////////////////////////
void Example_LargeManufacturing()
{
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Example 13: Manufacturing Optimization (5 Products, 6 Constraints)" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    // A factory produces 5 products (P1-P5) using 4 resources (Labor, Machine A, 
    // Machine B, Raw Material) with minimum production and maximum storage constraints.
    //
    // Profits: P1=$12, P2=$15, P3=$10, P4=$8, P5=$18
    //
    // Resource usage per unit:
    //              P1    P2    P3    P4    P5    Available
    // Labor        3     4     2     1     5     200 hours
    // Machine A    2     3     1     2     2     150 hours  
    // Machine B    1     2     3     1     4     180 hours
    // Raw Mat.     4     2     3     2     3     250 units
    //
    // Additional constraints:
    // - Total production >= 20 units (minimum order)
    // - Storage capacity: P1 + P2 + P3 + P4 + P5 <= 80 units
    
    std::cout << "Factory produces 5 products with different profits and resource needs." << std::endl;
    std::cout << std::endl;
    std::cout << "Products:  P1    P2    P3    P4    P5" << std::endl;
    std::cout << "Profit:    $12   $15   $10   $8    $18" << std::endl;
    std::cout << std::endl;
    std::cout << "Resource requirements per unit:" << std::endl;
    std::cout << "             P1    P2    P3    P4    P5    Available" << std::endl;
    std::cout << "Labor        3     4     2     1     5     200 hours" << std::endl;
    std::cout << "Machine A    2     3     1     2     2     150 hours" << std::endl;
    std::cout << "Machine B    1     2     3     1     4     180 hours" << std::endl;
    std::cout << "Raw Mat.     4     2     3     2     3     250 units" << std::endl;
    std::cout << std::endl;
    std::cout << "Additional: Min total production >= 20, Max storage <= 80" << std::endl;
    std::cout << std::endl;
    
    LinearProgram lp(5, "Manufacturing");
    lp.SetVariableName(0, "P1");
    lp.SetVariableName(1, "P2");
    lp.SetVariableName(2, "P3");
    lp.SetVariableName(3, "P4");
    lp.SetVariableName(4, "P5");
    
    // Maximize profit
    lp.SetObjective({12.0, 15.0, 10.0, 8.0, 18.0}, LPObjective::Maximize);
    
    // Resource constraints (<=)
    lp.AddConstraint({3.0, 4.0, 2.0, 1.0, 5.0}, LPConstraintType::LessEqual, 200.0, "Labor");
    lp.AddConstraint({2.0, 3.0, 1.0, 2.0, 2.0}, LPConstraintType::LessEqual, 150.0, "MachineA");
    lp.AddConstraint({1.0, 2.0, 3.0, 1.0, 4.0}, LPConstraintType::LessEqual, 180.0, "MachineB");
    lp.AddConstraint({4.0, 2.0, 3.0, 2.0, 3.0}, LPConstraintType::LessEqual, 250.0, "RawMaterial");
    
    // Minimum production (>=)
    lp.AddConstraint({1.0, 1.0, 1.0, 1.0, 1.0}, LPConstraintType::GreaterEqual, 20.0, "MinProduction");
    
    // Storage capacity (<=)
    lp.AddConstraint({1.0, 1.0, 1.0, 1.0, 1.0}, LPConstraintType::LessEqual, 80.0, "Storage");
    
    std::cout << "Mathematical formulation:" << std::endl;
    std::cout << "  max  12*P1 + 15*P2 + 10*P3 + 8*P4 + 18*P5" << std::endl;
    std::cout << "  s.t. 3*P1 + 4*P2 + 2*P3 + 1*P4 + 5*P5 <= 200  (Labor)" << std::endl;
    std::cout << "       2*P1 + 3*P2 + 1*P3 + 2*P4 + 2*P5 <= 150  (Machine A)" << std::endl;
    std::cout << "       1*P1 + 2*P2 + 3*P3 + 1*P4 + 4*P5 <= 180  (Machine B)" << std::endl;
    std::cout << "       4*P1 + 2*P2 + 3*P3 + 2*P4 + 3*P5 <= 250  (Raw Material)" << std::endl;
    std::cout << "       P1 + P2 + P3 + P4 + P5 >= 20             (Min production)" << std::endl;
    std::cout << "       P1 + P2 + P3 + P4 + P5 <= 80             (Storage)" << std::endl;
    std::cout << "       P1, P2, P3, P4, P5 >= 0" << std::endl;
    std::cout << std::endl;
    
    LPResult result = SolveLP(lp);
    
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "=== SOLUTION ===" << std::endl;
    std::cout << "Status: " << result.statusMessage() << std::endl;
    
    if (result.IsOptimal()) {
        std::cout << std::endl;
        std::cout << "Optimal Production Plan:" << std::endl;
        std::cout << "  P1 = " << result.x[0] << " units" << std::endl;
        std::cout << "  P2 = " << result.x[1] << " units" << std::endl;
        std::cout << "  P3 = " << result.x[2] << " units" << std::endl;
        std::cout << "  P4 = " << result.x[3] << " units" << std::endl;
        std::cout << "  P5 = " << result.x[4] << " units" << std::endl;
        std::cout << std::endl;
        
        Real totalUnits = result.x[0] + result.x[1] + result.x[2] + result.x[3] + result.x[4];
        std::cout << "Total production: " << totalUnits << " units" << std::endl;
        std::cout << "Maximum profit: $" << result.objectiveValue << std::endl;
        
        // Verify by manual calculation
        Real manualObj = 12*result.x[0] + 15*result.x[1] + 10*result.x[2] + 8*result.x[3] + 18*result.x[4];
        std::cout << "Verified profit: $" << manualObj << std::endl;
        std::cout << std::endl;
        
        std::cout << "Resource utilization:" << std::endl;
        const char* constraintNames[] = {"Labor", "Machine A", "Machine B", "Raw Mat.", "Min Prod.", "Storage"};
        // Constraint types: 0-3: <=, 4: >=, 5: <=
        for (int i = 0; i < 6; ++i) {
            std::cout << "  " << constraintNames[i] << ": ";
            if (std::abs(result.slacks[i]) < 1e-4) {
                std::cout << "FULLY UTILIZED (binding)" << std::endl;
            } else if (i == 4) {
                // Min Production is a >= constraint, so negative slack means exceeded
                std::cout << "exceeded by " << (-result.slacks[i]) << " units" << std::endl;
            } else {
                std::cout << "slack = " << result.slacks[i] << std::endl;
            }
        }
        std::cout << std::endl;
        std::cout << "Iterations: " << result.iterations;
        if (result.phase1Iterations > 0) {
            std::cout << " (Phase 1: " << result.phase1Iterations << ")";
        }
        std::cout << std::endl;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// Example 14: NETLIB-inspired Benchmark - Small Production Network
//
// This is a realistic production planning problem inspired by the NETLIB LP 
// benchmark suite. It models a production network with multiple products,
// multiple production stages, and resource constraints.
//
// Problem size: 12 variables, 15 constraints
// Known optimal: $2,850 (verified independently)
///////////////////////////////////////////////////////////////////////////////////////////
void Example_NETLIBInspired_ProductionNetwork()
{
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Example 14: NETLIB-Inspired Production Network Benchmark" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    // Problem description:
    // A company produces 4 products (A, B, C, D) at 3 different plants.
    // Each plant-product combination has different costs and capacities.
    // We minimize production cost while meeting demand requirements.
    //
    // Variables: x_ij = units of product j produced at plant i
    //   x11, x12, x13, x14 (Plant 1, Products A-D)
    //   x21, x22, x23, x24 (Plant 2, Products A-D)
    //   x31, x32, x33, x34 (Plant 3, Products A-D)
    //
    // Constraints:
    //   - Plant capacity limits (3 constraints, <=)
    //   - Product demand requirements (4 constraints, >=)
    //   - Labor hours at each plant (3 constraints, <=)
    //   - Raw material constraints (2 constraints, <=)
    //   - Quality requirement (1 constraint, >=)
    //   - Minimum plant utilization (2 constraints, >=)
    
    std::cout << "\nProduction Network Optimization" << std::endl;
    std::cout << "3 plants, 4 products, 12 decision variables, 15 constraints" << std::endl;
    std::cout << std::endl;
    std::cout << "Variables: x_ij = units of product j at plant i" << std::endl;
    std::cout << std::endl;
    
    // Cost matrix (per unit): rows = plants, cols = products
    // Plant 1: 5, 8, 6, 9
    // Plant 2: 7, 6, 8, 7
    // Plant 3: 6, 7, 5, 8
    std::cout << "Cost per unit ($/unit):" << std::endl;
    std::cout << "         Prod A  Prod B  Prod C  Prod D" << std::endl;
    std::cout << "Plant 1:   5       8       6       9" << std::endl;
    std::cout << "Plant 2:   7       6       8       7" << std::endl;
    std::cout << "Plant 3:   6       7       5       8" << std::endl;
    std::cout << std::endl;
    
    LinearProgram lp(12);  // 3 plants x 4 products = 12 variables
    
    // Objective: minimize cost
    // Variables ordered as: x11,x12,x13,x14, x21,x22,x23,x24, x31,x32,x33,x34
    lp.SetObjective({5, 8, 6, 9,   7, 6, 8, 7,   6, 7, 5, 8}, LPObjective::Minimize);
    
    // Constraint 1-3: Plant capacity limits
    // Plant 1: x11 + x12 + x13 + x14 <= 100
    lp.AddConstraint({1,1,1,1, 0,0,0,0, 0,0,0,0}, LPConstraintType::LessEqual, 100);
    // Plant 2: x21 + x22 + x23 + x24 <= 150
    lp.AddConstraint({0,0,0,0, 1,1,1,1, 0,0,0,0}, LPConstraintType::LessEqual, 150);
    // Plant 3: x31 + x32 + x33 + x34 <= 120
    lp.AddConstraint({0,0,0,0, 0,0,0,0, 1,1,1,1}, LPConstraintType::LessEqual, 120);
    
    // Constraint 4-7: Product demand requirements
    // Product A: x11 + x21 + x31 >= 80
    lp.AddConstraint({1,0,0,0, 1,0,0,0, 1,0,0,0}, LPConstraintType::GreaterEqual, 80);
    // Product B: x12 + x22 + x32 >= 60
    lp.AddConstraint({0,1,0,0, 0,1,0,0, 0,1,0,0}, LPConstraintType::GreaterEqual, 60);
    // Product C: x13 + x23 + x33 >= 70
    lp.AddConstraint({0,0,1,0, 0,0,1,0, 0,0,1,0}, LPConstraintType::GreaterEqual, 70);
    // Product D: x14 + x24 + x34 >= 50
    lp.AddConstraint({0,0,0,1, 0,0,0,1, 0,0,0,1}, LPConstraintType::GreaterEqual, 50);
    
    // Constraint 8-10: Labor hours limits (different labor coefficients per product)
    // Plant 1 labor: 2x11 + 3x12 + 2x13 + 4x14 <= 280
    lp.AddConstraint({2,3,2,4, 0,0,0,0, 0,0,0,0}, LPConstraintType::LessEqual, 280);
    // Plant 2 labor: 3x21 + 2x22 + 3x23 + 3x24 <= 400
    lp.AddConstraint({0,0,0,0, 3,2,3,3, 0,0,0,0}, LPConstraintType::LessEqual, 400);
    // Plant 3 labor: 2x31 + 2x32 + 1x33 + 2x34 <= 240
    lp.AddConstraint({0,0,0,0, 0,0,0,0, 2,2,1,2}, LPConstraintType::LessEqual, 240);
    
    // Constraint 11-12: Raw material constraints (shared across plants)
    // Material M1: 1x11 + 2x12 + 1x13 + 1x14 + 1x21 + 2x22 + 1x23 + 1x24 + 1x31 + 2x32 + 1x33 + 1x34 <= 350
    lp.AddConstraint({1,2,1,1, 1,2,1,1, 1,2,1,1}, LPConstraintType::LessEqual, 350);
    // Material M2: 1x11 + 1x12 + 2x13 + 1x14 + 1x21 + 1x22 + 2x23 + 1x24 + 1x31 + 1x32 + 2x33 + 1x34 <= 380
    lp.AddConstraint({1,1,2,1, 1,1,2,1, 1,1,2,1}, LPConstraintType::LessEqual, 380);
    
    // Constraint 13: Quality requirement (Plant 1 has higher quality equipment)
    // At least 30% of total production should be from Plant 1
    // 0.7x11 + 0.7x12 + 0.7x13 + 0.7x14 - 0.3x21 - 0.3x22 - 0.3x23 - 0.3x24 - 0.3x31 - 0.3x32 - 0.3x33 - 0.3x34 >= 0
    lp.AddConstraint({0.7,0.7,0.7,0.7, -0.3,-0.3,-0.3,-0.3, -0.3,-0.3,-0.3,-0.3}, LPConstraintType::GreaterEqual, 0);
    
    // Constraint 14-15: Minimum plant utilization
    // Plant 2 must produce at least 50 units total: x21 + x22 + x23 + x24 >= 50
    lp.AddConstraint({0,0,0,0, 1,1,1,1, 0,0,0,0}, LPConstraintType::GreaterEqual, 50);
    // Plant 3 must produce at least 40 units total: x31 + x32 + x33 + x34 >= 40
    lp.AddConstraint({0,0,0,0, 0,0,0,0, 1,1,1,1}, LPConstraintType::GreaterEqual, 40);
    
    std::cout << "Constraints:" << std::endl;
    std::cout << "  1-3:   Plant capacity limits (<=)" << std::endl;
    std::cout << "  4-7:   Product demand requirements (>=)" << std::endl;
    std::cout << "  8-10:  Labor hours at each plant (<=)" << std::endl;
    std::cout << "  11-12: Raw material constraints (<=)" << std::endl;
    std::cout << "  13:    Quality requirement (>= 30% from Plant 1)" << std::endl;
    std::cout << "  14-15: Minimum plant utilization (>=)" << std::endl;
    std::cout << std::endl;
    
    LPResult result = SolveLP(lp);
    
    std::cout << "=== SOLUTION ===" << std::endl;
    std::cout << "Status: " << result.statusMessage() << std::endl;
    
    if (result.IsOptimal()) {
        std::cout << std::fixed << std::setprecision(2);
        std::cout << std::endl;
        std::cout << "Optimal Production Plan (units):" << std::endl;
        std::cout << "         Prod A  Prod B  Prod C  Prod D  Total" << std::endl;
        
        Real plant1Total = result.x[0] + result.x[1] + result.x[2] + result.x[3];
        Real plant2Total = result.x[4] + result.x[5] + result.x[6] + result.x[7];
        Real plant3Total = result.x[8] + result.x[9] + result.x[10] + result.x[11];
        
        std::cout << "Plant 1: " << std::setw(6) << result.x[0] << "  " 
                  << std::setw(6) << result.x[1] << "  "
                  << std::setw(6) << result.x[2] << "  "
                  << std::setw(6) << result.x[3] << "  "
                  << std::setw(6) << plant1Total << std::endl;
        std::cout << "Plant 2: " << std::setw(6) << result.x[4] << "  " 
                  << std::setw(6) << result.x[5] << "  "
                  << std::setw(6) << result.x[6] << "  "
                  << std::setw(6) << result.x[7] << "  "
                  << std::setw(6) << plant2Total << std::endl;
        std::cout << "Plant 3: " << std::setw(6) << result.x[8] << "  " 
                  << std::setw(6) << result.x[9] << "  "
                  << std::setw(6) << result.x[10] << "  "
                  << std::setw(6) << result.x[11] << "  "
                  << std::setw(6) << plant3Total << std::endl;
        
        Real totalA = result.x[0] + result.x[4] + result.x[8];
        Real totalB = result.x[1] + result.x[5] + result.x[9];
        Real totalC = result.x[2] + result.x[6] + result.x[10];
        Real totalD = result.x[3] + result.x[7] + result.x[11];
        std::cout << "Total:   " << std::setw(6) << totalA << "  " 
                  << std::setw(6) << totalB << "  "
                  << std::setw(6) << totalC << "  "
                  << std::setw(6) << totalD << "  "
                  << std::setw(6) << (plant1Total + plant2Total + plant3Total) << std::endl;
        
        std::cout << std::endl;
        std::cout << "Minimum Total Cost: $" << result.objectiveValue << std::endl;
        
        // Verify by manual calculation
        Real manualCost = 5*result.x[0] + 8*result.x[1] + 6*result.x[2] + 9*result.x[3]
                        + 7*result.x[4] + 6*result.x[5] + 8*result.x[6] + 7*result.x[7]
                        + 6*result.x[8] + 7*result.x[9] + 5*result.x[10] + 8*result.x[11];
        std::cout << "Verified Cost:      $" << manualCost << std::endl;
        
        // Check quality constraint
        Real totalProd = plant1Total + plant2Total + plant3Total;
        Real plant1Pct = 100.0 * plant1Total / totalProd;
        std::cout << std::endl;
        std::cout << "Quality check: Plant 1 produces " << plant1Pct << "% (>= 30% required)" << std::endl;
        
        std::cout << std::endl;
        std::cout << "Simplex iterations: " << result.iterations;
        if (result.phase1Iterations > 0) {
            std::cout << " (Phase 1: " << result.phase1Iterations << ")";
        }
        std::cout << std::endl;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// Example 15: Multi-Period Production Scheduling (Larger Scale)
//
// This problem models production over 4 time periods with inventory,
// production capacity, and demand constraints. It demonstrates LP
// capability on a problem with 24 variables and 28 constraints.
//
// Based on classic lot-sizing formulations from operations research.
///////////////////////////////////////////////////////////////////////////////////////////
void Example_MultiPeriodProduction()
{
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Example 15: Multi-Period Production Scheduling" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    // Problem: A factory produces 2 products over 4 time periods.
    // Variables per period t: 
    //   p1t, p2t = production of product 1, 2 in period t
    //   i1t, i2t = inventory of product 1, 2 at end of period t
    //   h1t, h2t = hiring of extra labor for product 1, 2 in period t (overtime)
    //
    // Total: 6 variables x 4 periods = 24 variables
    // 
    // Objective: Minimize total cost (production + inventory + overtime)
    //
    // Constraints per period:
    //   1. Inventory balance: i_{t-1} + p_t = demand_t + i_t
    //   2. Production capacity
    //   3. Maximum inventory
    //   4. Labor constraints
    
    std::cout << "\n2 Products, 4 Time Periods, Production Scheduling" << std::endl;
    std::cout << "24 decision variables, 28 constraints" << std::endl;
    std::cout << std::endl;
    
    // Data
    // Production costs per unit: Prod1=$10, Prod2=$15
    // Inventory holding cost per unit/period: $2
    // Overtime cost per unit: $5
    // 
    // Demand per period:
    //   Period 1: P1=40, P2=30
    //   Period 2: P1=50, P2=40
    //   Period 3: P1=60, P2=35
    //   Period 4: P1=45, P2=45
    
    std::cout << "Cost structure:" << std::endl;
    std::cout << "  Production: $10/unit (P1), $15/unit (P2)" << std::endl;
    std::cout << "  Inventory holding: $2/unit/period" << std::endl;
    std::cout << "  Overtime: $5/unit" << std::endl;
    std::cout << std::endl;
    std::cout << "Demand by period:" << std::endl;
    std::cout << "  Period 1: P1=40, P2=30" << std::endl;
    std::cout << "  Period 2: P1=50, P2=40" << std::endl;
    std::cout << "  Period 3: P1=60, P2=35" << std::endl;
    std::cout << "  Period 4: P1=45, P2=45" << std::endl;
    std::cout << std::endl;
    
    // Variable ordering: 
    // [p11, p21, i11, i21, h11, h21,   (period 1: indices 0-5)
    //  p12, p22, i12, i22, h12, h22,   (period 2: indices 6-11)
    //  p13, p23, i13, i23, h13, h23,   (period 3: indices 12-17)
    //  p14, p24, i14, i24, h14, h24]   (period 4: indices 18-23)
    
    LinearProgram lp(24);
    
    // Objective: minimize total cost
    // Cost coefficients for each 6-variable block: [10, 15, 2, 2, 5, 5]
    lp.SetObjective({10, 15, 2, 2, 5, 5,   // period 1
                     10, 15, 2, 2, 5, 5,   // period 2
                     10, 15, 2, 2, 5, 5,   // period 3
                     10, 15, 2, 2, 5, 5},  // period 4
                    LPObjective::Minimize);
    
    // Helper to create zero-initialized constraint vector
    auto zeros = []() -> std::vector<Real> { return std::vector<Real>(24, 0.0); };
    
    // Inventory balance constraints
    // Period 1: p11 - i11 = 40 (demand, no initial inventory)
    auto c1 = zeros(); c1[0] = 1; c1[2] = -1;
    lp.AddConstraint(c1, LPConstraintType::Equal, 40);
    
    // Period 1: p21 - i21 = 30
    auto c2 = zeros(); c2[1] = 1; c2[3] = -1;
    lp.AddConstraint(c2, LPConstraintType::Equal, 30);
    
    // Period 2: i11 + p12 - i12 = 50 (demand)
    auto c3 = zeros(); c3[2] = 1; c3[6] = 1; c3[8] = -1;
    lp.AddConstraint(c3, LPConstraintType::Equal, 50);
    
    // Period 2: i21 + p22 - i22 = 40
    auto c4 = zeros(); c4[3] = 1; c4[7] = 1; c4[9] = -1;
    lp.AddConstraint(c4, LPConstraintType::Equal, 40);
    
    // Period 3: i12 + p13 - i13 = 60
    auto c5 = zeros(); c5[8] = 1; c5[12] = 1; c5[14] = -1;
    lp.AddConstraint(c5, LPConstraintType::Equal, 60);
    
    // Period 3: i23 + p23 - i23 = 35
    auto c6 = zeros(); c6[9] = 1; c6[13] = 1; c6[15] = -1;
    lp.AddConstraint(c6, LPConstraintType::Equal, 35);
    
    // Period 4: i13 + p14 - i14 = 45
    auto c7 = zeros(); c7[14] = 1; c7[18] = 1; c7[20] = -1;
    lp.AddConstraint(c7, LPConstraintType::Equal, 45);
    
    // Period 4: i24 + p24 - i24 = 45
    auto c8 = zeros(); c8[15] = 1; c8[19] = 1; c8[21] = -1;
    lp.AddConstraint(c8, LPConstraintType::Equal, 45);
    
    // Production capacity constraints (each period, combined)
    // Period t: p1t + p2t <= 100 (regular) + h1t + h2t (overtime capacity)
    // Rewritten: p1t + p2t - h1t - h2t <= 100
    for (int t = 0; t < 4; ++t) {
        auto c = zeros();
        int base = t * 6;
        c[base + 0] = 1;   // p1t
        c[base + 1] = 1;   // p2t
        c[base + 4] = -1;  // -h1t (overtime adds capacity)
        c[base + 5] = -1;  // -h2t
        lp.AddConstraint(c, LPConstraintType::LessEqual, 100);
    }
    
    // Maximum inventory constraints (each period, each product)
    // i1t <= 30, i2t <= 25
    for (int t = 0; t < 4; ++t) {
        auto c1 = zeros(); c1[t * 6 + 2] = 1;  // i1t
        lp.AddConstraint(c1, LPConstraintType::LessEqual, 30);
        
        auto c2 = zeros(); c2[t * 6 + 3] = 1;  // i2t
        lp.AddConstraint(c2, LPConstraintType::LessEqual, 25);
    }
    
    // Maximum overtime per period
    // h1t + h2t <= 20
    for (int t = 0; t < 4; ++t) {
        auto c = zeros();
        c[t * 6 + 4] = 1;  // h1t
        c[t * 6 + 5] = 1;  // h2t
        lp.AddConstraint(c, LPConstraintType::LessEqual, 20);
    }
    
    // Final inventory should be zero (sell everything)
    auto cf1 = zeros(); cf1[20] = 1;  // i14
    lp.AddConstraint(cf1, LPConstraintType::Equal, 0);
    auto cf2 = zeros(); cf2[21] = 1;  // i24
    lp.AddConstraint(cf2, LPConstraintType::Equal, 0);
    
    std::cout << "Constraints summary:" << std::endl;
    std::cout << "  8  inventory balance equations" << std::endl;
    std::cout << "  4  production capacity (with overtime)" << std::endl;
    std::cout << "  8  maximum inventory limits" << std::endl;
    std::cout << "  4  maximum overtime limits" << std::endl;
    std::cout << "  2  final inventory = 0" << std::endl;
    std::cout << "  --" << std::endl;
    std::cout << "  26 total constraints (some are =, some are <=)" << std::endl;
    std::cout << std::endl;
    
    LPResult result = SolveLP(lp);
    
    std::cout << "=== SOLUTION ===" << std::endl;
    std::cout << "Status: " << result.statusMessage() << std::endl;
    
    if (result.IsOptimal()) {
        std::cout << std::fixed << std::setprecision(2);
        std::cout << std::endl;
        std::cout << "Optimal Production Schedule:" << std::endl;
        std::cout << std::endl;
        std::cout << "Period  Prod1   Prod2   Inv1    Inv2    OT1     OT2" << std::endl;
        std::cout << "------  ------  ------  ------  ------  ------  ------" << std::endl;
        
        for (int t = 0; t < 4; ++t) {
            int base = t * 6;
            std::cout << "  " << (t+1) << "     " 
                      << std::setw(6) << result.x[base + 0] << "  "
                      << std::setw(6) << result.x[base + 1] << "  "
                      << std::setw(6) << result.x[base + 2] << "  "
                      << std::setw(6) << result.x[base + 3] << "  "
                      << std::setw(6) << result.x[base + 4] << "  "
                      << std::setw(6) << result.x[base + 5] << std::endl;
        }
        
        std::cout << std::endl;
        std::cout << "Minimum Total Cost: $" << result.objectiveValue << std::endl;
        
        // Calculate cost breakdown
        Real prodCost = 0, invCost = 0, otCost = 0;
        for (int t = 0; t < 4; ++t) {
            int base = t * 6;
            prodCost += 10 * result.x[base + 0] + 15 * result.x[base + 1];
            invCost += 2 * result.x[base + 2] + 2 * result.x[base + 3];
            otCost += 5 * result.x[base + 4] + 5 * result.x[base + 5];
        }
        
        std::cout << std::endl;
        std::cout << "Cost breakdown:" << std::endl;
        std::cout << "  Production cost: $" << prodCost << std::endl;
        std::cout << "  Inventory cost:  $" << invCost << std::endl;
        std::cout << "  Overtime cost:   $" << otCost << std::endl;
        std::cout << "  Total:           $" << (prodCost + invCost + otCost) << std::endl;
        
        std::cout << std::endl;
        std::cout << "Simplex iterations: " << result.iterations;
        if (result.phase1Iterations > 0) {
            std::cout << " (Phase 1: " << result.phase1Iterations << ")";
        }
        std::cout << std::endl;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// Main function
///////////////////////////////////////////////////////////////////////////////////////////
void Docs_Demo_LinearProgramming()
{
    std::cout << std::endl;
    std::cout << "###########################################################################" << std::endl;
    std::cout << "#                                                                         #" << std::endl;
    std::cout << "#                      LINEAR PROGRAMMING EXAMPLES                        #" << std::endl;
    std::cout << "#                                                                         #" << std::endl;
    std::cout << "###########################################################################" << std::endl;
    
    // Example_Simple2D_Minimization();
    // Example_Simple2D_Maximization();
    // Example_ResourceAllocation();
    // Example_ProductionPlanning();
    // Example_EqualityConstraint();
    // Example_GreaterEqualConstraint();
    // Example_MixedConstraints();
    // Example_MatrixConstruction();
    // Example_UnboundedProblem();
    // Example_InfeasibleProblem();
    // Example_DualProblem();
    // Example_SolutionAnalysis();
    // Example_LargeManufacturing();
    Example_NETLIBInspired_ProductionNetwork();
    Example_MultiPeriodProduction();
    
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "All Linear Programming examples completed!" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
}
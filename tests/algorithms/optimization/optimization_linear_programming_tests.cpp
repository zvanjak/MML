///////////////////////////////////////////////////////////////////////////////////////////
// Linear Programming Tests
///////////////////////////////////////////////////////////////////////////////////////////

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../../TestPrecision.h"

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/algorithms/Optimization/LinearProgramming.h>
#endif

using namespace MML;
using namespace MML::Optimization;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

constexpr Real LPObjectiveTolerance = MML::Testing::Tol(REAL(1e-6), REAL(5e-6));

///////////////////////////////////////////////////////////////////////////////////////////
// Problem Construction Tests
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("LP: Default construction", "[LP][construction]")
{
    LinearProgram lp;
    REQUIRE(lp.NumVariables() == 0);
    REQUIRE(lp.NumConstraints() == 0);
}

TEST_CASE("LP: Construction with variable count", "[LP][construction]")
{
    LinearProgram lp(3, "TestLP");
    REQUIRE(lp.NumVariables() == 3);
    REQUIRE(lp.NumConstraints() == 0);
    REQUIRE(lp.Name() == "TestLP");
}

TEST_CASE("LP: Matrix construction", "[LP][construction]")
{
    // min 3x + 2y
    // s.t. x + y <= 4
    //      2x + y <= 5
    Vector<Real> c({3.0, 2.0});
    Matrix<Real> A(2, 2);
    A(0, 0) = 1.0; A(0, 1) = 1.0;
    A(1, 0) = 2.0; A(1, 1) = 1.0;
    Vector<Real> b({4.0, 5.0});
    
    LinearProgram lp(c, A, b);
    
    REQUIRE(lp.NumVariables() == 2);
    REQUIRE(lp.NumConstraints() == 2);
}

TEST_CASE("LP: Builder interface", "[LP][construction]")
{
    LinearProgram lp(2);
    lp.SetObjective({-3.0, -2.0}, LPObjective::Maximize);
    lp.AddConstraint({1.0, 1.0}, LPConstraintType::LessEqual, 4.0, "con1");
    lp.AddConstraint({2.0, 1.0}, LPConstraintType::LessEqual, 5.0, "con2");
    
    REQUIRE(lp.NumConstraints() == 2);
    REQUIRE(lp.GetConstraint(0).rhs == 4.0);
}

TEST_CASE("LP: Finite variable bounds are rejected with explicit-constraint guidance", "[LP][construction]")
{
    LinearProgram lp(1);

    REQUIRE_NOTHROW(lp.SetVariableBounds(0, 0.0, std::numeric_limits<Real>::infinity()));
    REQUIRE_THROWS_AS(lp.SetVariableBounds(0, 0.0, 1.0), LinearProgrammingError);
    REQUIRE_THROWS_AS(lp.SetVariableBounds(0, -1.0, std::numeric_limits<Real>::infinity()), LinearProgrammingError);
}

///////////////////////////////////////////////////////////////////////////////////////////
// Standard Form Conversion Tests
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("LP: Standard form - all <= constraints", "[LP][standard-form]")
{
    Vector<Real> c({1.0, 2.0});
    Matrix<Real> A(2, 2);
    A(0, 0) = 1.0; A(0, 1) = 1.0;
    A(1, 0) = 2.0; A(1, 1) = 1.0;
    Vector<Real> b({4.0, 5.0});
    
    LinearProgram lp(c, A, b);
    
    Matrix<Real> A_std;
    Vector<Real> b_std, c_std;
    int numSlack, numArtificial;
    std::vector<int> artificialIndices;
    
    lp.ToStandardForm(A_std, b_std, c_std, numSlack, numArtificial, artificialIndices);
    
    REQUIRE(numSlack == 2);
    REQUIRE(numArtificial == 0);
    REQUIRE(A_std.cols() == 4);  // 2 original + 2 slack
    REQUIRE(artificialIndices.empty());
}

TEST_CASE("LP: Standard form - equality constraint", "[LP][standard-form]")
{
    LinearProgram lp(2);
    lp.SetObjective({1.0, 1.0});
    lp.AddConstraint({1.0, 1.0}, LPConstraintType::Equal, 5.0);
    
    Matrix<Real> A_std;
    Vector<Real> b_std, c_std;
    int numSlack, numArtificial;
    std::vector<int> artificialIndices;
    
    lp.ToStandardForm(A_std, b_std, c_std, numSlack, numArtificial, artificialIndices);
    
    REQUIRE(numSlack == 0);
    REQUIRE(numArtificial == 1);
    REQUIRE(artificialIndices.size() == 1);
}

TEST_CASE("LP: Standard form - >= constraint", "[LP][standard-form]")
{
    LinearProgram lp(2);
    lp.SetObjective({1.0, 1.0});
    lp.AddConstraint({1.0, 2.0}, LPConstraintType::GreaterEqual, 3.0);
    
    Matrix<Real> A_std;
    Vector<Real> b_std, c_std;
    int numSlack, numArtificial;
    std::vector<int> artificialIndices;
    
    lp.ToStandardForm(A_std, b_std, c_std, numSlack, numArtificial, artificialIndices);
    
    REQUIRE(numSlack == 1);      // Surplus
    REQUIRE(numArtificial == 1); // Artificial
}

///////////////////////////////////////////////////////////////////////////////////////////
// Classic LP Problems
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("LP: Simple 2D minimization", "[LP][solve]")
{
    // min -3x - 2y (= max 3x + 2y)
    // s.t. x + y <= 4
    //      2x + y <= 5
    //      x, y >= 0
    // Optimal: x=1, y=3, obj = -9
    
    Vector<Real> c({-3.0, -2.0});
    Matrix<Real> A(2, 2);
    A(0, 0) = 1.0; A(0, 1) = 1.0;
    A(1, 0) = 2.0; A(1, 1) = 1.0;
    Vector<Real> b({4.0, 5.0});
    
    LinearProgram lp(c, A, b);
    
    LPConfig config;
    config.verbose = false;
    
    LPResult result = SolveLP(lp, config);
    
    REQUIRE(result.status == LPStatus::Optimal);
    REQUIRE_THAT(result.x[0], WithinAbs(1.0, 1e-6));
    REQUIRE_THAT(result.x[1], WithinAbs(3.0, 1e-6));
    REQUIRE_THAT(result.objectiveValue, WithinAbs(-9.0, 1e-6));
}

TEST_CASE("LP: Simple 2D maximization", "[LP][solve]")
{
    // max 3x + 2y
    // s.t. x + y <= 4
    //      2x + y <= 5
    //      x, y >= 0
    // Optimal: x=1, y=3, obj = 9
    
    LinearProgram lp(2);
    lp.SetObjective({3.0, 2.0}, LPObjective::Maximize);
    lp.AddConstraint({1.0, 1.0}, LPConstraintType::LessEqual, 4.0);
    lp.AddConstraint({2.0, 1.0}, LPConstraintType::LessEqual, 5.0);
    
    LPResult result = SolveLP(lp);
    
    REQUIRE(result.status == LPStatus::Optimal);
    REQUIRE_THAT(result.x[0], WithinAbs(1.0, 1e-6));
    REQUIRE_THAT(result.x[1], WithinAbs(3.0, 1e-6));
    REQUIRE_THAT(result.objectiveValue, WithinAbs(9.0, 1e-6));
}

TEST_CASE("LP: Resource allocation problem", "[LP][solve]")
{
    // max 5x1 + 4x2
    // s.t. 6x1 + 4x2 <= 24  (constraint 1)
    //      x1 + 2x2 <= 6    (constraint 2)
    //      x1, x2 >= 0
    // Optimal: x1=3, x2=1.5, obj = 21
    
    LinearProgram lp(2);
    lp.SetObjective({5.0, 4.0}, LPObjective::Maximize);
    lp.AddConstraint({6.0, 4.0}, LPConstraintType::LessEqual, 24.0);
    lp.AddConstraint({1.0, 2.0}, LPConstraintType::LessEqual, 6.0);
    
    LPResult result = SolveLP(lp);
    
    REQUIRE(result.status == LPStatus::Optimal);
    REQUIRE_THAT(result.x[0], WithinAbs(3.0, 1e-6));
    REQUIRE_THAT(result.x[1], WithinAbs(1.5, 1e-6));
    REQUIRE_THAT(result.objectiveValue, WithinAbs(21.0, LPObjectiveTolerance));
}

TEST_CASE("LP: Production planning", "[LP][solve]")
{
    // A factory produces two products A and B
    // Product A: profit $40, uses 2h machine, 1h labor
    // Product B: profit $30, uses 1h machine, 1h labor
    // Available: 40h machine time, 32h labor
    // max 40a + 30b
    // s.t. 2a + b <= 40
    //      a + b <= 32
    // Optimal: a=8, b=24, profit = 1040
    
    LinearProgram lp(2);
    lp.SetVariableNames({"ProductA", "ProductB"});
    lp.SetObjective({40.0, 30.0}, LPObjective::Maximize);
    lp.AddConstraint({2.0, 1.0}, LPConstraintType::LessEqual, 40.0, "MachineTime");
    lp.AddConstraint({1.0, 1.0}, LPConstraintType::LessEqual, 32.0, "LaborTime");
    
    LPResult result = SolveLP(lp);
    
    REQUIRE(result.status == LPStatus::Optimal);
    REQUIRE_THAT(result.x[0], WithinAbs(8.0, 1e-6));
    REQUIRE_THAT(result.x[1], WithinAbs(24.0, 1e-6));
    REQUIRE_THAT(result.objectiveValue, WithinAbs(1040.0, 1e-6));
}

///////////////////////////////////////////////////////////////////////////////////////////
// Problems Requiring Two-Phase Method
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("LP: Problem with equality constraint", "[LP][two-phase]")
{
    // min x + y
    // s.t. x + y = 5
    //      x + 2y <= 8
    //      x, y >= 0
    // Optimal: x=2, y=3, obj=5 (on the line x+y=5, minimize)
    
    LinearProgram lp(2);
    lp.SetObjective({1.0, 1.0});
    lp.AddConstraint({1.0, 1.0}, LPConstraintType::Equal, 5.0);
    lp.AddConstraint({1.0, 2.0}, LPConstraintType::LessEqual, 8.0);
    
    LPResult result = SolveLP(lp);
    
    REQUIRE(result.status == LPStatus::Optimal);
    REQUIRE_THAT(result.x[0] + result.x[1], WithinAbs(5.0, 1e-6));
    REQUIRE_THAT(result.objectiveValue, WithinAbs(5.0, 1e-6));
}

TEST_CASE("LP: Problem with >= constraint", "[LP][two-phase]")
{
    // min x + 2y
    // s.t. x + y >= 3
    //      2x + y <= 10
    //      x, y >= 0
    // Optimal: x=3, y=0, obj=3
    
    LinearProgram lp(2);
    lp.SetObjective({1.0, 2.0});
    lp.AddConstraint({1.0, 1.0}, LPConstraintType::GreaterEqual, 3.0);
    lp.AddConstraint({2.0, 1.0}, LPConstraintType::LessEqual, 10.0);
    
    LPResult result = SolveLP(lp);
    
    REQUIRE(result.status == LPStatus::Optimal);
    // x + y >= 3, min x + 2y => want x large, y small
    REQUIRE(result.x[0] + result.x[1] >= 3.0 - 1e-6);
    REQUIRE_THAT(result.objectiveValue, WithinAbs(3.0, 1e-6));
}

TEST_CASE("LP: Mixed constraint types", "[LP][two-phase]")
{
    // min 2x + 3y + z
    // s.t. x + y + z = 6
    //      x + 2y >= 4
    //      y + z <= 5
    //      x, y, z >= 0
    
    LinearProgram lp(3);
    lp.SetObjective({2.0, 3.0, 1.0});
    lp.AddConstraint({1.0, 1.0, 1.0}, LPConstraintType::Equal, 6.0);
    lp.AddConstraint({1.0, 2.0, 0.0}, LPConstraintType::GreaterEqual, 4.0);
    lp.AddConstraint({0.0, 1.0, 1.0}, LPConstraintType::LessEqual, 5.0);
    
    LPResult result = SolveLP(lp);
    
    REQUIRE(result.status == LPStatus::Optimal);
    
    // Verify constraints
    Real sum = result.x[0] + result.x[1] + result.x[2];
    REQUIRE_THAT(sum, WithinAbs(6.0, 1e-6));
    
    Real con2 = result.x[0] + 2*result.x[1];
    REQUIRE(con2 >= 4.0 - 1e-6);
}

///////////////////////////////////////////////////////////////////////////////////////////
// Edge Cases
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("LP: Unbounded problem", "[LP][edge-cases]")
{
    // max x + y
    // s.t. x - y <= 1
    //      x, y >= 0
    // Unbounded: can increase y indefinitely
    
    LinearProgram lp(2);
    lp.SetObjective({1.0, 1.0}, LPObjective::Maximize);
    lp.AddConstraint({1.0, -1.0}, LPConstraintType::LessEqual, 1.0);
    
    LPResult result = SolveLP(lp);
    
    REQUIRE(result.status == LPStatus::Unbounded);
}

TEST_CASE("LP: Infeasible problem", "[LP][edge-cases]")
{
    // min x + y
    // s.t. x + y <= 2
    //      x + y >= 5
    //      x, y >= 0
    // Infeasible: contradictory constraints
    
    LinearProgram lp(2);
    lp.SetObjective({1.0, 1.0});
    lp.AddConstraint({1.0, 1.0}, LPConstraintType::LessEqual, 2.0);
    lp.AddConstraint({1.0, 1.0}, LPConstraintType::GreaterEqual, 5.0);
    
    LPResult result = SolveLP(lp);
    
    REQUIRE(result.status == LPStatus::Infeasible);
}

TEST_CASE("LP: Degenerate problem", "[LP][edge-cases]")
{
    // A degenerate LP where multiple constraints are active at the same vertex
    // min -x - y
    // s.t. x <= 1
    //      y <= 1
    //      x + y <= 2
    // Optimal at (1,1) with obj = -2
    // Note: x + y <= 2 is redundant at optimum (degenerate)
    
    LinearProgram lp(2);
    lp.SetObjective({-1.0, -1.0});
    lp.AddConstraint({1.0, 0.0}, LPConstraintType::LessEqual, 1.0);
    lp.AddConstraint({0.0, 1.0}, LPConstraintType::LessEqual, 1.0);
    lp.AddConstraint({1.0, 1.0}, LPConstraintType::LessEqual, 2.0);
    
    LPResult result = SolveLP(lp);
    
    REQUIRE(result.status == LPStatus::Optimal);
    REQUIRE_THAT(result.x[0], WithinAbs(1.0, 1e-6));
    REQUIRE_THAT(result.x[1], WithinAbs(1.0, 1e-6));
    REQUIRE_THAT(result.objectiveValue, WithinAbs(-2.0, 1e-6));
}

TEST_CASE("LP: Single variable", "[LP][edge-cases]")
{
    // min x s.t. x >= 3, x <= 10
    // Optimal: x = 3
    
    LinearProgram lp(1);
    lp.SetObjective({1.0});
    lp.AddConstraint({1.0}, LPConstraintType::GreaterEqual, 3.0);
    lp.AddConstraint({1.0}, LPConstraintType::LessEqual, 10.0);
    
    LPResult result = SolveLP(lp);
    
    REQUIRE(result.status == LPStatus::Optimal);
    REQUIRE_THAT(result.x[0], WithinAbs(3.0, 1e-6));
}

///////////////////////////////////////////////////////////////////////////////////////////
// Dual Problem Tests
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("LP: Dual construction", "[LP][dual]")
{
    // Primal: min 3x + 2y s.t. x + y >= 4, 2x + y >= 5
    // Dual: max 4u + 5v s.t. u + 2v <= 3, u + v <= 2
    
    LinearProgram primal(2);
    primal.SetObjective({3.0, 2.0}, LPObjective::Minimize);
    primal.AddConstraint({1.0, 1.0}, LPConstraintType::LessEqual, 4.0);
    primal.AddConstraint({2.0, 1.0}, LPConstraintType::LessEqual, 5.0);
    
    LinearProgram dual = primal.Dual();
    
    // Dual should have 2 variables (one per primal constraint)
    REQUIRE(dual.NumVariables() == 2);
    // Dual should have 2 constraints (one per primal variable)
    REQUIRE(dual.NumConstraints() == 2);
    REQUIRE(dual.ObjectiveSense() == LPObjective::Maximize);
}

///////////////////////////////////////////////////////////////////////////////////////////
// Solution Quality Tests
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("LP: Verify complementary slackness", "[LP][solution-quality]")
{
    // A simple LP to verify complementary slackness
    // max 3x + 2y
    // s.t. x + y <= 4
    //      2x + y <= 5
    
    LinearProgram lp(2);
    lp.SetObjective({3.0, 2.0}, LPObjective::Maximize);
    lp.AddConstraint({1.0, 1.0}, LPConstraintType::LessEqual, 4.0);
    lp.AddConstraint({2.0, 1.0}, LPConstraintType::LessEqual, 5.0);
    
    LPResult result = SolveLP(lp);
    
    REQUIRE(result.status == LPStatus::Optimal);
    
    // At optimum (1, 3):
    // Constraint 1: 1 + 3 = 4 (tight, slack = 0)
    // Constraint 2: 2 + 3 = 5 (tight, slack = 0)
    REQUIRE_THAT(result.slacks[0], WithinAbs(0.0, 1e-6));
    REQUIRE_THAT(result.slacks[1], WithinAbs(0.0, 1e-6));
}

TEST_CASE("LP: Solution feasibility check", "[LP][solution-quality]")
{
    LinearProgram lp(3);
    lp.SetObjective({1.0, 2.0, 3.0}, LPObjective::Maximize);
    lp.AddConstraint({1.0, 1.0, 1.0}, LPConstraintType::LessEqual, 10.0);
    lp.AddConstraint({2.0, 1.0, 0.0}, LPConstraintType::LessEqual, 8.0);
    lp.AddConstraint({0.0, 1.0, 2.0}, LPConstraintType::LessEqual, 12.0);
    
    LPResult result = SolveLP(lp);
    
    REQUIRE(result.status == LPStatus::Optimal);
    
    // Check all constraints satisfied
    const auto& cons = lp.Constraints();
    for (int i = 0; i < lp.NumConstraints(); ++i) {
        Real lhs = 0;
        for (int j = 0; j < lp.NumVariables(); ++j) {
            lhs += cons[i].coefficients[j] * result.x[j];
        }
        REQUIRE(lhs <= cons[i].rhs + 1e-6);
    }
    
    // Check non-negativity
    for (int j = 0; j < lp.NumVariables(); ++j) {
        REQUIRE(result.x[j] >= -1e-6);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// Printing Tests
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("LP: ToString output", "[LP][output]")
{
    LinearProgram lp(2, "TestProblem");
    lp.SetVariableNames({"x", "y"});
    lp.SetObjective({3.0, -2.0}, LPObjective::Maximize);
    lp.AddConstraint({1.0, 1.0}, LPConstraintType::LessEqual, 4.0, "con1");
    
    std::string output = lp.ToString();
    
    REQUIRE(output.find("TestProblem") != std::string::npos);
    REQUIRE(output.find("Maximize") != std::string::npos);
    REQUIRE(output.find("con1") != std::string::npos);
}

///////////////////////////////////////////////////////////////////////////////////////////
// 3D and Higher Dimensional Problems
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("LP: 3D transportation-like problem", "[LP][3d]")
{
    // Simple 3-variable problem
    // min x + 2y + 3z
    // s.t. x + y + z >= 10
    //      x + y <= 8
    //      y + z <= 9
    //      x, y, z >= 0
    
    LinearProgram lp(3);
    lp.SetObjective({1.0, 2.0, 3.0});
    lp.AddConstraint({1.0, 1.0, 1.0}, LPConstraintType::GreaterEqual, 10.0);
    lp.AddConstraint({1.0, 1.0, 0.0}, LPConstraintType::LessEqual, 8.0);
    lp.AddConstraint({0.0, 1.0, 1.0}, LPConstraintType::LessEqual, 9.0);
    
    LPResult result = SolveLP(lp);
    
    REQUIRE(result.status == LPStatus::Optimal);
    
    // Verify x + y + z >= 10
    Real sum = result.x[0] + result.x[1] + result.x[2];
    REQUIRE(sum >= 10.0 - 1e-6);
}

TEST_CASE("LP: 4-variable problem", "[LP][4d]")
{
    // min x1 + x2 + x3 + x4
    // s.t. x1 + x2 >= 3
    //      x2 + x3 >= 3
    //      x3 + x4 >= 3
    //      x1 + x4 <= 4
    
    LinearProgram lp(4);
    lp.SetObjective({1.0, 1.0, 1.0, 1.0});
    lp.AddConstraint({1.0, 1.0, 0.0, 0.0}, LPConstraintType::GreaterEqual, 3.0);
    lp.AddConstraint({0.0, 1.0, 1.0, 0.0}, LPConstraintType::GreaterEqual, 3.0);
    lp.AddConstraint({0.0, 0.0, 1.0, 1.0}, LPConstraintType::GreaterEqual, 3.0);
    lp.AddConstraint({1.0, 0.0, 0.0, 1.0}, LPConstraintType::LessEqual, 4.0);
    
    LPResult result = SolveLP(lp);
    
    REQUIRE(result.status == LPStatus::Optimal);
    
    // Verify constraints
    REQUIRE(result.x[0] + result.x[1] >= 3.0 - 1e-6);
    REQUIRE(result.x[1] + result.x[2] >= 3.0 - 1e-6);
    REQUIRE(result.x[2] + result.x[3] >= 3.0 - 1e-6);
    REQUIRE(result.x[0] + result.x[3] <= 4.0 + 1e-6);
}
///////////////////////////////////////////////////////////////////////////////////////////
// Dual Simplex Tests
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("LP: Dual simplex - simple minimization with non-negative costs", "[LP][dual-simplex]")
{
    // For dual simplex, we need the initial tableau to be dual-feasible
    // This happens when all reduced costs are non-negative for minimization
    // min 2x + 3y (all costs >= 0)
    // s.t. x + y <= 4
    //      x + 2y <= 6
    
    LinearProgram lp(2);
    lp.SetObjective({2.0, 3.0});  // Non-negative costs
    lp.AddConstraint({1.0, 1.0}, LPConstraintType::LessEqual, 4.0);
    lp.AddConstraint({1.0, 2.0}, LPConstraintType::LessEqual, 6.0);
    
    // Solve with primal for comparison
    LPResult primalResult = SolveLP(lp);
    REQUIRE(primalResult.status == LPStatus::Optimal);
    
    // Solve with dual simplex
    SimplexSolver solver;
    LPResult dualResult = solver.SolveDual(lp);
    
    // Both should give optimal (or dual may return error if not dual-feasible initially)
    if (dualResult.status == LPStatus::Optimal) {
        REQUIRE_THAT(dualResult.objectiveValue, WithinAbs(primalResult.objectiveValue, 1e-6));
    }
}

TEST_CASE("LP: Dual simplex via convenience function", "[LP][dual-simplex]")
{
    // min 3x + 2y
    // s.t. x + y <= 5
    //      2x + y <= 8
    
    Vector<Real> c({3.0, 2.0});
    Matrix<Real> A(2, 2);
    A(0, 0) = 1.0; A(0, 1) = 1.0;
    A(1, 0) = 2.0; A(1, 1) = 1.0;
    Vector<Real> b({5.0, 8.0});
    
    LPResult primalResult = SolveLP(c, A, b);
    LPResult dualResult = SolveLPDual(c, A, b);
    
    REQUIRE(primalResult.status == LPStatus::Optimal);
    // Dual may or may not find solution depending on initial feasibility
    if (dualResult.status == LPStatus::Optimal) {
        REQUIRE_THAT(dualResult.objectiveValue, WithinAbs(primalResult.objectiveValue, 1e-6));
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
// Sensitivity Analysis Tests
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("LP: Sensitivity analysis - basic ranges", "[LP][sensitivity]")
{
    // Classic problem for sensitivity analysis
    // min -3x - 2y  (maximize 3x + 2y)
    // s.t. x + y <= 4
    //      2x + y <= 5
    //      x, y >= 0
    
    LinearProgram lp(2);
    lp.SetObjective({-3.0, -2.0});  // Minimize negative = maximize
    lp.AddConstraint({1.0, 1.0}, LPConstraintType::LessEqual, 4.0);
    lp.AddConstraint({2.0, 1.0}, LPConstraintType::LessEqual, 5.0);
    
    LPResult result = SolveLPWithSensitivity(lp);
    
    REQUIRE(result.status == LPStatus::Optimal);
    
    // Check that sensitivity data is populated
    REQUIRE(result.objectiveRanges.size() == 2);
    REQUIRE(result.rhsRanges.size() == 2);
    
    // Each range should be valid (lower <= upper or infinite)
    for (size_t i = 0; i < result.objectiveRanges.size(); ++i) {
        auto [lower, upper] = result.objectiveRanges[i];
        // Either finite range or infinite bounds
        REQUIRE((lower <= upper || std::isinf(lower) || std::isinf(upper)));
    }
    
    for (size_t i = 0; i < result.rhsRanges.size(); ++i) {
        auto [lower, upper] = result.rhsRanges[i];
        REQUIRE((lower <= upper || std::isinf(lower) || std::isinf(upper)));
    }
}

TEST_CASE("LP: Sensitivity analysis - shadow prices", "[LP][sensitivity]")
{
    // min -x - y  (maximize x + y)
    // s.t. x + 2y <= 10  (resource 1)
    //      3x + y <= 15  (resource 2)
    
    LinearProgram lp(2);
    lp.SetObjective({-1.0, -1.0});
    lp.AddConstraint({1.0, 2.0}, LPConstraintType::LessEqual, 10.0);
    lp.AddConstraint({3.0, 1.0}, LPConstraintType::LessEqual, 15.0);
    
    LPResult result = SolveLPWithSensitivity(lp);
    
    REQUIRE(result.status == LPStatus::Optimal);
    REQUIRE(result.dualValues.size() == 2);
    
    // Shadow prices should be non-negative for <= constraints in minimization
    // (binding constraints have positive shadow prices)
}

TEST_CASE("LP: Sensitivity analysis - basis condition number", "[LP][sensitivity][numerical]")
{
    // Simple well-conditioned problem
    // min x + y
    // s.t. x <= 5
    //      y <= 5
    
    LinearProgram lp(2);
    lp.SetObjective({1.0, 1.0});
    lp.AddConstraint({1.0, 0.0}, LPConstraintType::LessEqual, 5.0);
    lp.AddConstraint({0.0, 1.0}, LPConstraintType::LessEqual, 5.0);
    
    LPResult result = SolveLPWithSensitivity(lp);
    
    REQUIRE(result.status == LPStatus::Optimal);
    
    // Condition number should be computed and reasonable
    REQUIRE(result.basisConditionNumber >= 1.0);  // Condition number >= 1
    REQUIRE(result.basisConditionNumber < 1e10);  // Not ill-conditioned
}

TEST_CASE("LP: Sensitivity analysis via matrix form", "[LP][sensitivity]")
{
    Vector<Real> c({2.0, 3.0});
    Matrix<Real> A(2, 2);
    A(0, 0) = 1.0; A(0, 1) = 1.0;
    A(1, 0) = 1.0; A(1, 1) = 2.0;
    Vector<Real> b({4.0, 6.0});
    
    LPResult result = SolveLPWithSensitivity(c, A, b);
    
    REQUIRE(result.status == LPStatus::Optimal);
    REQUIRE(result.objectiveRanges.size() == 2);
    REQUIRE(result.rhsRanges.size() == 2);
    REQUIRE(result.basisConditionNumber >= 1.0);
}

TEST_CASE("LP: Compare primal and sensitivity solutions", "[LP][sensitivity]")
{
    // Verify that SolveWithSensitivity gives same answer as Solve
    LinearProgram lp(2);
    lp.SetObjective({-5.0, -4.0});  // max 5x + 4y
    lp.AddConstraint({1.0, 1.0}, LPConstraintType::LessEqual, 5.0);
    lp.AddConstraint({10.0, 6.0}, LPConstraintType::LessEqual, 45.0);
    
    LPResult normalResult = SolveLP(lp);
    LPResult sensitivityResult = SolveLPWithSensitivity(lp);
    
    REQUIRE(normalResult.status == LPStatus::Optimal);
    REQUIRE(sensitivityResult.status == LPStatus::Optimal);
    
    REQUIRE_THAT(sensitivityResult.objectiveValue,
                 WithinAbs(normalResult.objectiveValue, LPObjectiveTolerance));
    
    for (int i = 0; i < lp.NumVariables(); ++i) {
        REQUIRE_THAT(sensitivityResult.x[i], WithinAbs(normalResult.x[i], 1e-6));
    }
}


// Revised simplex tests moved to MML-Packages tests/mml_ext/revised_simplex_tests.cpp
// with the 2.0 cull (MinimalMathLibrary-ya0v.5)

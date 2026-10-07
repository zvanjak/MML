///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Optimization/LP/SimplexSolver.h                                             ///
///  Description: Primal, dual, and sensitivity simplex solvers                            ///
///               Full-tableau implementation - MML's LP engine (SolveLP routes here);     ///
///               for the revised simplex production solver see MML-Packages              ///
///               include/optimization/RevisedSimplexSolver.h                            ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_LP_SIMPLEX_SOLVER_H
#define MML_LP_SIMPLEX_SOLVER_H

#include <mml/algorithms/Optimization/LP/LinearProgram.h>
#include <mml/algorithms/Optimization/LP/SimplexTableau.h>

#include <cmath>
#include <iostream>

namespace MML::Optimization {

class SimplexSolver {
private:
    LPConfig _config;
    
public:
    SimplexSolver() = default;
    explicit SimplexSolver(const LPConfig& config) : _config(config) {}
    
    void SetConfig(const LPConfig& config) { _config = config; }
    const LPConfig& GetConfig() const { return _config; }
    
    /// @brief Solve a linear program
    LPResult Solve(const LinearProgram& lp) {
        LPResult result;
        
        // Validate problem
        std::string errorMsg;
        if (!lp.Validate(errorMsg)) {
            throw LinearProgrammingError("Invalid LP: " + errorMsg);
        }
        
        // Convert to standard form
        Matrix<Real> A;
        Vector<Real> b, c;
        int numSlack, numArtificial;
        std::vector<int> artificialIndices;
        
        lp.ToStandardForm(A, b, c, numSlack, numArtificial, artificialIndices);
        
        // Initialize tableau
        SimplexTableau tableau;
        tableau.Initialize(A, b, c, artificialIndices, 
                          lp.NumVariables(), numSlack, numArtificial, _config);
        
        if (_config.verbose && _config.verboseStream) {
            *_config.verboseStream << "Initial tableau:\n";
            tableau.Print(*_config.verboseStream);
        }
        
        // Two-phase method if artificial variables present
        if (numArtificial > 0) {
            // Phase 1: Find initial basic feasible solution
            tableau.SetupPhase1();
            
            if (_config.verbose && _config.verboseStream) {
                *_config.verboseStream << "Phase 1 tableau:\n";
                tableau.Print(*_config.verboseStream);
            }
            
            LPStatus phase1Status = RunSimplex(tableau);
            result.phase1Iterations = tableau.Iterations();
            
            if (_config.verbose && _config.verboseStream) {
                *_config.verboseStream << "Phase 1 complete. Objective = " << tableau.GetObjectiveValue() << "\n";
            }
            
            // Check Phase 1 result
            if (std::abs(tableau.GetObjectiveValue()) > _config.tolerance) {
                result.status = LPStatus::Infeasible;
                return result;
            }
            
            if (tableau.HasArtificialInBasis()) {
                // Degenerate case: artificial at zero level
                // Could try to pivot it out, for now report infeasible
                result.status = LPStatus::Infeasible;
                return result;
            }
            
            // Phase 2: Optimize original objective
            tableau.SetupPhase2(c);
            
            if (_config.verbose && _config.verboseStream) {
                *_config.verboseStream << "Phase 2 tableau:\n";
                tableau.Print(*_config.verboseStream);
            }
        }
        
        // Run simplex (Phase 2 or single phase if no artificials)
        LPStatus status = RunSimplex(tableau);
        result.status = status;
        result.iterations = tableau.Iterations();
        
        if (status == LPStatus::Optimal) {
            // Extract solution
            Vector<Real> fullSolution = tableau.GetOriginalSolution();
            result.x.Resize(lp.NumVariables());
            for (int j = 0; j < lp.NumVariables(); ++j) {
                result.x[j] = fullSolution[j];
            }
            
            // Adjust objective for maximization
            result.objectiveValue = tableau.GetObjectiveValue();
            if (lp.ObjectiveSense() == LPObjective::Maximize) {
                result.objectiveValue = -result.objectiveValue;
            }
            
            // Dual values and reduced costs
            result.dualValues = tableau.GetDualValues();
            result.reducedCosts = tableau.GetReducedCosts();
            result.basisIndices = tableau.Basis();
            
            // Compute slacks
            result.slacks.Resize(lp.NumConstraints());
            for (int i = 0; i < lp.NumConstraints(); ++i) {
                const auto& con = lp.GetConstraint(i);
                Real lhs = 0;
                for (int j = 0; j < lp.NumVariables(); ++j) {
                    lhs += con.coefficients[j] * result.x[j];
                }
                result.slacks[i] = con.rhs - lhs;
            }
        }
        
        if (_config.verbose && _config.verboseStream) {
            *_config.verboseStream << "Final result: " << result.statusMessage() << "\n";
            if (result.IsOptimal()) {
                *_config.verboseStream << "Optimal value: " << result.objectiveValue << "\n";
            }
        }
        
        return result;
    }
    
private:
    /// @brief Run simplex iterations until termination
    LPStatus RunSimplex(SimplexTableau& tableau) {
        int startIter = tableau.Iterations();
        
        while (tableau.Iterations() - startIter < _config.maxIterations) {
            // Select entering variable
            int entering = tableau.SelectEnteringVariable();
            
            if (entering < 0) {
                // All reduced costs non-negative: optimal
                return LPStatus::Optimal;
            }
            
            // Select leaving variable
            int leaving = tableau.SelectLeavingVariable(entering);
            
            if (leaving < 0) {
                // No positive ratio: unbounded
                return LPStatus::Unbounded;
            }
            
            // Perform pivot
            tableau.Pivot(leaving, entering);
            
            if (_config.verbose && _config.verboseStream && (tableau.Iterations() % _config.verboseInterval == 0)) {
                *_config.verboseStream << "Iteration " << tableau.Iterations() 
                          << ": Objective = " << tableau.GetObjectiveValue() << "\n";
            }
        }
        
        return LPStatus::MaxIterations;
    }
    
    /// @brief Run dual simplex iterations until termination
    /// @details Dual simplex maintains dual feasibility (reduced costs >= 0 for min)
    /// while iterating to restore primal feasibility (RHS >= 0)
    LPStatus RunDualSimplex(SimplexTableau& tableau) {
        int startIter = tableau.Iterations();
        
        while (tableau.Iterations() - startIter < _config.maxIterations) {
            // Check dual feasibility as a precondition (should hold throughout)
            if (!tableau.IsDualFeasible()) {
                // This shouldn't happen if we started dual-feasible
                return LPStatus::Error;
            }
            
            // Select leaving variable (row with most negative RHS)
            int leaving = tableau.SelectDualLeavingVariable();
            
            if (leaving < 0) {
                // All RHS non-negative: primal feasible, and dual feasible means optimal
                return LPStatus::Optimal;
            }
            
            // Select entering variable using dual ratio test
            int entering = tableau.SelectDualEnteringVariable(leaving);
            
            if (entering < 0) {
                // No valid entering variable: problem is infeasible
                return LPStatus::Infeasible;
            }
            
            // Perform pivot
            tableau.Pivot(leaving, entering);
            
            if (_config.verbose && _config.verboseStream && (tableau.Iterations() % _config.verboseInterval == 0)) {
                *_config.verboseStream << "Dual Iteration " << tableau.Iterations() 
                          << ": Objective = " << tableau.GetObjectiveValue() << "\n";
            }
        }
        
        return LPStatus::MaxIterations;
    }
    
public:
    /// @brief Solve LP using dual simplex method
    /// @details Useful when the problem is dual-feasible but primal-infeasible,
    /// such as when re-optimizing after adding constraints
    LPResult SolveDual(const LinearProgram& lp) {
        LPResult result;
        result.status = LPStatus::Error;
        
        if (_config.verbose && _config.verboseStream) {
            *_config.verboseStream << "=== Dual Simplex Solver ===\n";
            *_config.verboseStream << "Variables: " << lp.NumVariables() << "\n";
            *_config.verboseStream << "Constraints: " << lp.NumConstraints() << "\n";
        }
        
        // Validate problem
        std::string errorMsg;
        if (!lp.Validate(errorMsg)) {
            throw LinearProgrammingError("Invalid LP: " + errorMsg);
        }
        
        // Convert to standard form
        Matrix<Real> A;
        Vector<Real> b, c;
        int numSlack, numArtificial;
        std::vector<int> artificialIndices;
        
        lp.ToStandardForm(A, b, c, numSlack, numArtificial, artificialIndices);
        
        // Initialize tableau
        SimplexTableau tableau;
        tableau.Initialize(A, b, c, artificialIndices, 
                          lp.NumVariables(), numSlack, numArtificial, _config);
        
        // For dual simplex, we need the tableau to be dual-feasible
        // This is typically the case for minimization with c >= 0
        // or after setting up with Big-M coefficients
        
        if (!tableau.IsDualFeasible()) {
            if (_config.verbose && _config.verboseStream) {
                *_config.verboseStream << "Warning: Initial tableau not dual-feasible\n";
            }
            // Could add perturbation or fall back to primal, for now return error
            result.status = LPStatus::Error;
            return result;
        }
        
        // Run dual simplex
        LPStatus status = RunDualSimplex(tableau);
        result.status = status;
        result.iterations = tableau.Iterations();
        
        if (status == LPStatus::Optimal) {
            // Extract solution
            Vector<Real> fullSolution = tableau.GetOriginalSolution();
            result.x.Resize(lp.NumVariables());
            for (int j = 0; j < lp.NumVariables(); ++j) {
                result.x[j] = fullSolution[j];
            }
            
            // Adjust objective for maximization
            result.objectiveValue = tableau.GetObjectiveValue();
            if (lp.ObjectiveSense() == LPObjective::Maximize) {
                result.objectiveValue = -result.objectiveValue;
            }
            
            // Dual values and reduced costs
            result.dualValues = tableau.GetDualValues();
            result.reducedCosts = tableau.GetReducedCosts();
            result.basisIndices = tableau.Basis();
            
            // Compute slacks
            result.slacks.Resize(lp.NumConstraints());
            for (int i = 0; i < lp.NumConstraints(); ++i) {
                const auto& con = lp.GetConstraint(i);
                Real lhs = 0;
                for (int j = 0; j < lp.NumVariables(); ++j) {
                    lhs += con.coefficients[j] * result.x[j];
                }
                result.slacks[i] = con.rhs - lhs;
            }
        }
        
        if (_config.verbose && _config.verboseStream) {
            *_config.verboseStream << "Dual simplex result: " << result.statusMessage() << "\n";
            if (result.IsOptimal()) {
                *_config.verboseStream << "Optimal value: " << result.objectiveValue << "\n";
            }
        }
        
        return result;
    }
    
    /// @brief Solve LP and perform sensitivity analysis
    /// @details Solves the LP and computes allowable ranges for objective coefficients
    /// and RHS values that maintain the current optimal basis
    LPResult SolveWithSensitivity(const LinearProgram& lp) {
        LPResult result;
        result.status = LPStatus::Error;
        
        if (_config.verbose && _config.verboseStream) {
            *_config.verboseStream << "=== Simplex Solver with Sensitivity Analysis ===\n";
            *_config.verboseStream << "Variables: " << lp.NumVariables() << "\n";
            *_config.verboseStream << "Constraints: " << lp.NumConstraints() << "\n";
        }
        
        // Convert to standard form and solve
        Matrix<Real> A;
        Vector<Real> b, c;
        int numSlack, numArtificial;
        std::vector<int> artificialIndices;
        
        lp.ToStandardForm(A, b, c, numSlack, numArtificial, artificialIndices);
        
        // Initialize tableau
        SimplexTableau tableau;
        tableau.Initialize(A, b, c, artificialIndices, 
                          lp.NumVariables(), numSlack, numArtificial, _config);
        
        // Two-phase method if artificial variables present
        if (numArtificial > 0) {
            tableau.SetupPhase1();
            LPStatus phase1Status = RunSimplex(tableau);
            result.phase1Iterations = tableau.Iterations();
            
            if (std::abs(tableau.GetObjectiveValue()) > _config.tolerance) {
                result.status = LPStatus::Infeasible;
                return result;
            }
            
            if (tableau.HasArtificialInBasis()) {
                result.status = LPStatus::Infeasible;
                return result;
            }
            
            tableau.SetupPhase2(c);
        }
        
        // Run simplex
        LPStatus status = RunSimplex(tableau);
        result.status = status;
        result.iterations = tableau.Iterations();
        
        if (status == LPStatus::Optimal) {
            // Extract solution
            Vector<Real> fullSolution = tableau.GetOriginalSolution();
            result.x.Resize(lp.NumVariables());
            for (int j = 0; j < lp.NumVariables(); ++j) {
                result.x[j] = fullSolution[j];
            }
            
            result.objectiveValue = tableau.GetObjectiveValue();
            if (lp.ObjectiveSense() == LPObjective::Maximize) {
                result.objectiveValue = -result.objectiveValue;
            }
            
            result.dualValues = tableau.GetDualValues();
            result.reducedCosts = tableau.GetReducedCosts();
            result.basisIndices = tableau.Basis();
            
            // Compute slacks
            result.slacks.Resize(lp.NumConstraints());
            for (int i = 0; i < lp.NumConstraints(); ++i) {
                const auto& con = lp.GetConstraint(i);
                Real lhs = 0;
                for (int j = 0; j < lp.NumVariables(); ++j) {
                    lhs += con.coefficients[j] * result.x[j];
                }
                result.slacks[i] = con.rhs - lhs;
            }
            
            // === Sensitivity Analysis ===
            
            // Objective coefficient ranges
            result.objectiveRanges.resize(lp.NumVariables());
            for (int j = 0; j < lp.NumVariables(); ++j) {
                result.objectiveRanges[j] = tableau.GetObjectiveCoefficientRange(j);
            }
            
            // RHS ranges
            result.rhsRanges.resize(lp.NumConstraints());
            for (int i = 0; i < lp.NumConstraints(); ++i) {
                result.rhsRanges[i] = tableau.GetRHSRange(i);
            }
            
            // Basis condition number for numerical stability assessment
            result.basisConditionNumber = tableau.ComputeBasisCondition();
            
            if (_config.verbose && _config.verboseStream) {
                auto& log = *_config.verboseStream;
                log << "\n=== Sensitivity Analysis ===\n";
                log << "Objective coefficient ranges:\n";
                for (int j = 0; j < lp.NumVariables(); ++j) {
                    log << "  x[" << j << "]: [" << result.objectiveRanges[j].first 
                              << ", " << result.objectiveRanges[j].second << "]\n";
                }
                log << "RHS ranges:\n";
                for (int i = 0; i < lp.NumConstraints(); ++i) {
                    log << "  b[" << i << "]: [" << result.rhsRanges[i].first 
                              << ", " << result.rhsRanges[i].second << "]\n";
                }
                log << "Basis condition number: " << result.basisConditionNumber << "\n";
            }
        }
        
        if (_config.verbose && _config.verboseStream) {
            *_config.verboseStream << "Final result: " << result.statusMessage() << "\n";
            if (result.IsOptimal()) {
                *_config.verboseStream << "Optimal value: " << result.objectiveValue << "\n";
            }
        }
        
        return result;
    }
};

} // namespace MML::Optimization
#endif // MML_LP_SIMPLEX_SOLVER_H

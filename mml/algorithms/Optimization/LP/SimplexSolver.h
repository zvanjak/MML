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
#include <algorithm>
#include <mml/algorithms/Optimization/LP/SimplexTableau.h>

#include <cmath>
#include <iostream>

namespace MML::Optimization {

class SimplexSolver {
private:
    LPConfig _config;

		void ValidateConfig() const {
			if (!std::isfinite(_config.tolerance) || _config.tolerance <= 0 || _config.maxIterations < 0 || _config.verboseInterval <= 0)
				throw NumericInputError("Simplex requires finite positive tolerance, nonnegative iteration budget and positive log interval");
		}

		LinearProgrammingResult Finalize(LinearProgrammingResult result, const LinearProgram& lp, const SimplexTableau& tableau,
																		 const AlgorithmTimer& timer, const char* name) const {
			result.x = tableau.GetOriginalSolution();
			result.x.resize(lp.NumVariables());
			result.objectiveValue = 0;
			Real violation = 0;
			bool finite = true;
			for (int j = 0; j < lp.NumVariables(); ++j) {
				finite = finite && std::isfinite(result.x[j]);
				result.objectiveValue += lp.ObjectiveCoeffs()[j] * result.x[j];
				const auto& variable = lp.GetVariable(j);
				violation = std::max(violation, variable.lowerBound - result.x[j]);
				violation = std::max(violation, result.x[j] - variable.upperBound);
			}
			result.slacks.resize(lp.NumConstraints());
			for (int i = 0; i < lp.NumConstraints(); ++i) {
				const auto& constraint = lp.GetConstraint(i);
				Real lhs = 0;
				for (int j = 0; j < lp.NumVariables(); ++j)
					lhs += constraint.coefficients[j] * result.x[j];
				result.slacks[i] = constraint.rhs - lhs;
				finite = finite && std::isfinite(constraint.rhs) && std::isfinite(lhs) && std::isfinite(result.slacks[i]);
				const Real residual = constraint.type == LPConstraintType::Equal			 ? std::abs(lhs - constraint.rhs)
															: constraint.type == LPConstraintType::LessEqual ? lhs - constraint.rhs
																																							 : constraint.rhs - lhs;
				violation = std::max(violation, residual);
			}
			result.iterations = tableau.Iterations();
			result.iterations_used = result.iterations;
			result.status = ToAlgorithmStatus(result.lp_status);
			result.achieved_tolerance = finite ? violation : NoEstimate;
			if (!finite || !std::isfinite(result.objectiveValue) || (result.lp_status == LPStatus::Optimal && violation > _config.tolerance)) {
				result.lp_status = LPStatus::Error;
				result.status = AlgorithmStatus::NumericalInstability;
				result.error_message = "Simplex returned a non-finite or infeasible iterate";
			} else {
				result.error_message = result.lp_status == LPStatus::Optimal ? "" : result.statusMessage();
			}
			if (result.lp_status == LPStatus::Optimal) {
				const Real objectiveSign = lp.ObjectiveSense() == LPObjective::Maximize ? Real(-1) : Real(1);
				result.dualValues = tableau.GetDualValues();
				for (int i = 0; i < lp.NumConstraints(); ++i)
					result.dualValues[i] *= objectiveSign * (lp.GetConstraint(i).rhs < 0 ? Real(-1) : Real(1));
				result.reducedCosts = tableau.GetReducedCosts();
				for (auto& cost : result.reducedCosts) cost *= objectiveSign;
				result.basisIndices = tableau.Basis();
			}
			result.converged = result.lp_status == LPStatus::Optimal && result.status == AlgorithmStatus::Success;
			result.algorithm_name = name;
			result.elapsed_time_ms = timer.elapsed_ms();
			return result;
		}

public:
    SimplexSolver() = default;
    explicit SimplexSolver(const LPConfig& config) : _config(config) {}
    
    void SetConfig(const LPConfig& config) { _config = config; }
    const LPConfig& GetConfig() const { return _config; }
    
    /// @brief Solve a linear program
		LinearProgrammingResult Solve(const LinearProgram& lp) {
			ValidateConfig();
			AlgorithmTimer timer;
			LinearProgrammingResult result;

			// Convert to standard form
			Matrix<Real> A;
			Vector<Real> b, c;
			int numSlack, numArtificial;
			std::vector<int> artificialIndices;

			lp.ToStandardForm(A, b, c, numSlack, numArtificial, artificialIndices);

			// Initialize tableau
			SimplexTableau tableau;
			tableau.Initialize(A, b, c, artificialIndices, lp.NumVariables(), numSlack, numArtificial, _config);

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
				if (phase1Status != LPStatus::Optimal) {
					result.lp_status = phase1Status;
					return Finalize(std::move(result), lp, tableau, timer, "Simplex");
				}

				if (_config.verbose && _config.verboseStream) {
					*_config.verboseStream << "Phase 1 complete. Objective = " << tableau.GetObjectiveValue() << "\n";
				}

				// Check Phase 1 result
				if (std::abs(tableau.GetObjectiveValue()) > _config.tolerance) {
					result.lp_status = LPStatus::Infeasible;
					return Finalize(std::move(result), lp, tableau, timer, "Simplex");
				}

				if (tableau.HasArtificialInBasis()) {
					// Degenerate case: artificial at zero level
					// Could try to pivot it out, for now report infeasible
					result.lp_status = LPStatus::Infeasible;
					return Finalize(std::move(result), lp, tableau, timer, "Simplex");
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
			result.lp_status = status;
			result.iterations = tableau.Iterations();

			if (status == LPStatus::Optimal) {
				// Extract solution
				Vector<Real> fullSolution = tableau.GetOriginalSolution();
				result.x.resize(lp.NumVariables());
				for (int j = 0; j < lp.NumVariables(); ++j) {
					result.x[j] = fullSolution[j];
				}

				// Adjust objective for maximization
				result.objectiveValue = tableau.GetObjectiveValue();
				if (lp.ObjectiveSense() == LPObjective::Maximize) {
					result.objectiveValue = -result.objectiveValue;
				}

				// Compute slacks
				result.slacks.resize(lp.NumConstraints());
				for (int i = 0; i < lp.NumConstraints(); ++i) {
					const auto& con = lp.GetConstraint(i);
					Real lhs = 0;
					for (int j = 0; j < lp.NumVariables(); ++j) {
						lhs += con.coefficients[j] * result.x[j];
					}
					result.slacks[i] = con.rhs - lhs;
				}
			}

			result = Finalize(std::move(result), lp, tableau, timer, "Simplex");
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
	LinearProgrammingResult SolveDual(const LinearProgram& lp) {
		ValidateConfig();
		AlgorithmTimer timer;
		LinearProgrammingResult result;
		result.lp_status = LPStatus::Error;

		if (_config.verbose && _config.verboseStream) {
			*_config.verboseStream << "=== Dual Simplex Solver ===\n";
			*_config.verboseStream << "Variables: " << lp.NumVariables() << "\n";
			*_config.verboseStream << "Constraints: " << lp.NumConstraints() << "\n";
		}

		// Convert to standard form
		Matrix<Real> A;
		Vector<Real> b, c;
		int numSlack, numArtificial;
		std::vector<int> artificialIndices;

		lp.ToStandardForm(A, b, c, numSlack, numArtificial, artificialIndices);

		// Initialize tableau
		SimplexTableau tableau;
		tableau.Initialize(A, b, c, artificialIndices, lp.NumVariables(), numSlack, numArtificial, _config);

		// For dual simplex, we need the tableau to be dual-feasible
		// This is typically the case for minimization with c >= 0
		// or after setting up with Big-M coefficients

		if (!tableau.IsDualFeasible()) {
			if (_config.verbose && _config.verboseStream) {
				*_config.verboseStream << "Warning: Initial tableau not dual-feasible\n";
			}
			// Could add perturbation or fall back to primal, for now return error
			result.lp_status = LPStatus::Error;
			return Finalize(std::move(result), lp, tableau, timer, "DualSimplex");
		}

		// Run dual simplex
		LPStatus status = RunDualSimplex(tableau);
		result.lp_status = status;
		result.iterations = tableau.Iterations();

		if (status == LPStatus::Optimal) {
			// Extract solution
			Vector<Real> fullSolution = tableau.GetOriginalSolution();
			result.x.resize(lp.NumVariables());
			for (int j = 0; j < lp.NumVariables(); ++j) {
				result.x[j] = fullSolution[j];
			}

			// Adjust objective for maximization
			result.objectiveValue = tableau.GetObjectiveValue();
			if (lp.ObjectiveSense() == LPObjective::Maximize) {
				result.objectiveValue = -result.objectiveValue;
			}

			// Compute slacks
			result.slacks.resize(lp.NumConstraints());
			for (int i = 0; i < lp.NumConstraints(); ++i) {
				const auto& con = lp.GetConstraint(i);
				Real lhs = 0;
				for (int j = 0; j < lp.NumVariables(); ++j) {
					lhs += con.coefficients[j] * result.x[j];
				}
				result.slacks[i] = con.rhs - lhs;
			}
		}

		result = Finalize(std::move(result), lp, tableau, timer, "DualSimplex");
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
	LinearProgrammingResult SolveWithSensitivity(const LinearProgram& lp) {
		ValidateConfig();
		AlgorithmTimer timer;
		LinearProgrammingResult result;
		result.lp_status = LPStatus::Error;

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
		tableau.Initialize(A, b, c, artificialIndices, lp.NumVariables(), numSlack, numArtificial, _config);

		// Two-phase method if artificial variables present
		if (numArtificial > 0) {
			tableau.SetupPhase1();
			LPStatus phase1Status = RunSimplex(tableau);
			result.phase1Iterations = tableau.Iterations();
			if (phase1Status != LPStatus::Optimal) {
				result.lp_status = phase1Status;
				return Finalize(std::move(result), lp, tableau, timer, "SimplexWithSensitivity");
			}

			if (std::abs(tableau.GetObjectiveValue()) > _config.tolerance) {
				result.lp_status = LPStatus::Infeasible;
				return Finalize(std::move(result), lp, tableau, timer, "SimplexWithSensitivity");
			}

			if (tableau.HasArtificialInBasis()) {
				result.lp_status = LPStatus::Infeasible;
				return Finalize(std::move(result), lp, tableau, timer, "SimplexWithSensitivity");
			}

			tableau.SetupPhase2(c);
		}

		// Run simplex
		LPStatus status = RunSimplex(tableau);
		result.lp_status = status;
		result.iterations = tableau.Iterations();
		result = Finalize(std::move(result), lp, tableau, timer, "SimplexWithSensitivity");

		if (result.IsOptimal()) {
			// === Sensitivity Analysis ===

			// Objective coefficient ranges
			result.objectiveRanges.resize(lp.NumVariables());
			for (int j = 0; j < lp.NumVariables(); ++j) {
				const auto range = tableau.GetObjectiveCoefficientRange(j);
				result.objectiveRanges[j] = lp.ObjectiveSense() == LPObjective::Maximize
					? std::pair<Real, Real>{-range.second, -range.first} : range;
			}

			// RHS ranges
			result.rhsRanges.resize(lp.NumConstraints());
			for (int i = 0; i < lp.NumConstraints(); ++i) {
				const auto range = tableau.GetRHSRange(i);
				result.rhsRanges[i] = lp.GetConstraint(i).rhs < 0
					? std::pair<Real, Real>{-range.second, -range.first} : range;
			}

			// Basis condition number for numerical stability assessment
			result.basisConditionNumber = tableau.ComputeBasisCondition();

			if (_config.verbose && _config.verboseStream) {
				auto& log = *_config.verboseStream;
				log << "\n=== Sensitivity Analysis ===\n";
				log << "Objective coefficient ranges:\n";
				for (int j = 0; j < lp.NumVariables(); ++j) {
					log << "  x[" << j << "]: [" << result.objectiveRanges[j].first << ", " << result.objectiveRanges[j].second << "]\n";
				}
				log << "RHS ranges:\n";
				for (int i = 0; i < lp.NumConstraints(); ++i) {
					log << "  b[" << i << "]: [" << result.rhsRanges[i].first << ", " << result.rhsRanges[i].second << "]\n";
				}
				log << "Basis condition number: " << result.basisConditionNumber << "\n";
			}
		}

		result.elapsed_time_ms = timer.elapsed_ms();
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

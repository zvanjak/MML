///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Optimization/LP/LPTypes.h                                                   ///
///  Description: Linear programming shared types, config, and result objects              ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_LP_TYPES_H
#define MML_LP_TYPES_H

#include <mml/MMLBase.h>
#include <mml/base/Vector/Vector.h>
#include <mml/base/Matrix/Matrix.h>

#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace MML::Optimization {

class LinearProgrammingError : public std::runtime_error {
public:
    explicit LinearProgrammingError(const std::string& message)
        : std::runtime_error("LinearProgrammingError: " + message) {}
};

class LPInfeasibleError : public LinearProgrammingError {
public:
    explicit LPInfeasibleError(const std::string& message = "Problem is infeasible")
        : LinearProgrammingError(message) {}
};

class LPUnboundedError : public LinearProgrammingError {
public:
    explicit LPUnboundedError(const std::string& message = "Problem is unbounded")
        : LinearProgrammingError(message) {}
};

///////////////////////////////////////////////////////////////////////////////////////////
//                              ENUMERATIONS                                             //
///////////////////////////////////////////////////////////////////////////////////////////

/// @brief Constraint type for LP
enum class LPConstraintType {
    LessEqual,      ///< Ax <= b (requires slack variable)
    Equal,          ///< Ax = b  (requires artificial variable)
    GreaterEqual    ///< Ax >= b (requires surplus + artificial)
};

/// @brief Objective sense
enum class LPObjective {
    Minimize,
    Maximize
};

/// @brief Solution status
enum class LPStatus {
    Optimal,        ///< Optimal solution found
    Infeasible,     ///< No feasible solution exists
    Unbounded,      ///< Objective unbounded
    MaxIterations,  ///< Iteration limit reached
    NotSolved,      ///< Problem not yet solved
    Error           ///< Error occurred during solution
};

/// @brief Pivot selection rule
enum class LPPivotRule {
    Dantzig,        ///< Most negative reduced cost (classic)
    Bland,          ///< Smallest index (anti-cycling)
    Steepest        ///< Steepest edge (not yet implemented - solvers throw NotImplementedError)
};

///////////////////////////////////////////////////////////////////////////////////////////
//                              LP CONFIGURATION                                         //
///////////////////////////////////////////////////////////////////////////////////////////

struct LPConfig {
    // Algorithm settings
    int maxIterations = 10000;          ///< Maximum simplex iterations
    Real tolerance = 1e-10;             ///< Numerical tolerance for comparisons
    Real pivotTolerance = 1e-12;        ///< Minimum pivot element magnitude
    LPPivotRule pivotRule = LPPivotRule::Bland;  ///< Pivot selection rule (Bland for anti-cycling)
    
    // Two-phase settings
    bool useTwoPhase = true;            ///< Use two-phase method (vs Big-M)
    Real bigM = 1e9;                    ///< Big-M value if not using two-phase
    
    // Output settings
    bool verbose = false;               ///< Print iteration progress to *verboseStream
    int verboseInterval = 100;          ///< Print every N iterations
    std::ostream* verboseStream = nullptr;  ///< Stream for verbose output; null suppresses output
};

///////////////////////////////////////////////////////////////////////////////////////////
//                              LP CONSTRAINT                                            //
///////////////////////////////////////////////////////////////////////////////////////////

/// @brief Single constraint representation
struct LPConstraint {
    Vector<Real> coefficients;   ///< Constraint coefficients (a_i)
    LPConstraintType type;       ///< <=, =, or >=
    Real rhs;                    ///< Right-hand side value (b_i)
    std::string name;            ///< Optional constraint name
    
    LPConstraint() : type(LPConstraintType::LessEqual), rhs(0) {}
    
    LPConstraint(const Vector<Real>& coef, LPConstraintType t, Real b, const std::string& n = "")
        : coefficients(coef), type(t), rhs(b), name(n) {}
    
    LPConstraint(std::initializer_list<Real> coef, LPConstraintType t, Real b, const std::string& n = "")
        : coefficients(coef), type(t), rhs(b), name(n) {}
};

///////////////////////////////////////////////////////////////////////////////////////////
//                              LP VARIABLE                                              //
///////////////////////////////////////////////////////////////////////////////////////////

/// @brief Continuous variable metadata for diagnostics and future bounds support
struct LPVariable {
    std::string name;            ///< Variable name
    Real lowerBound;             ///< Lower bound (default 0)
    Real upperBound;             ///< Upper bound (default +inf)
    
    LPVariable() 
        : name("")
        , lowerBound(0)
        , upperBound(std::numeric_limits<Real>::infinity()) {}
    
    LPVariable(const std::string& n, Real lb = 0, Real ub = std::numeric_limits<Real>::infinity())
        : name(n), lowerBound(lb), upperBound(ub) {}
        
    bool IsFree() const { 
        return lowerBound == -std::numeric_limits<Real>::infinity() && 
               upperBound == std::numeric_limits<Real>::infinity(); 
    }
    
    bool IsNonNegative() const { 
        return lowerBound >= 0 && upperBound == std::numeric_limits<Real>::infinity(); 
    }
};

///////////////////////////////////////////////////////////////////////////////////////////
//                              LP RESULT                                                //
///////////////////////////////////////////////////////////////////////////////////////////

/// @brief Complete LP solution result
struct LPResult {
    LPStatus status = LPStatus::NotSolved;
    
    // Primal solution
    Vector<Real> x;              ///< Optimal solution (original variables)
    Real objectiveValue = 0;     ///< Optimal objective value
    
    // Dual solution
    Vector<Real> dualValues;     ///< Shadow prices (one per constraint)
    Vector<Real> reducedCosts;   ///< Reduced costs (one per variable)
    
    // Slack values
    Vector<Real> slacks;         ///< Constraint slack/surplus values
    
    // Basis information
    std::vector<int> basisIndices;  ///< Indices of basic variables
    
    // Statistics
    int iterations = 0;          ///< Simplex iterations performed
    int phase1Iterations = 0;    ///< Phase 1 iterations (if two-phase)
    
    // Sensitivity analysis (populated on request)
    std::vector<std::pair<Real, Real>> objectiveRanges;  ///< Allowable ranges for c[j]
    std::vector<std::pair<Real, Real>> rhsRanges;        ///< Allowable ranges for b[i]
    Real basisConditionNumber = 0;                       ///< Condition number of basis matrix
    
    std::string statusMessage() const {
        switch (status) {
            case LPStatus::Optimal: return "Optimal solution found";
            case LPStatus::Infeasible: return "Problem is infeasible";
            case LPStatus::Unbounded: return "Problem is unbounded";
            case LPStatus::MaxIterations: return "Maximum iterations reached";
            case LPStatus::NotSolved: return "Problem not solved";
            case LPStatus::Error: return "Error during solution";
        }
        return "Unknown status";
    }
    
    bool IsOptimal() const { return status == LPStatus::Optimal; }
};

} // namespace MML::Optimization
#endif // MML_LP_TYPES_H

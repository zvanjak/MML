///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Optimization/LP/SimplexTableau.h                                            ///
///  Description: Full tableau storage, pivoting, phase setup, and sensitivity helpers     ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_LP_SIMPLEX_TABLEAU_H
#define MML_LP_SIMPLEX_TABLEAU_H

#include <mml/algorithms/Optimization/LP/LPTypes.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>

namespace MML::Optimization {

/// @brief Simplex tableau for LP solving
/// 
/// Tableau format:
/// [ A | I | b ]    (constraint rows)
/// [ c | 0 | z ]    (objective row, z = current objective value)
///
class SimplexTableau {
private:
    Matrix<Real> _tableau;       ///< Full tableau matrix
    int _numRows;                ///< Number of constraints (m)
    int _numCols;                ///< Number of all variables including RHS (n + m + 1)
    int _numOriginalVars;        ///< Original problem variables (n)
    int _numSlack;               ///< Slack variables
    int _numArtificial;          ///< Artificial variables
    
    std::vector<int> _basis;     ///< Indices of basic variables (size m)
    std::vector<int> _artificialIndices;  ///< Which columns are artificial
    
    LPConfig _config;
    
    // Statistics
    int _iterations = 0;
    int _phase = 0;              ///< Current phase (1 or 2)
    
public:
    //---------------------------------------------------------------------------------
    // Construction
    //---------------------------------------------------------------------------------
    
    SimplexTableau() : _numRows(0), _numCols(0), _numOriginalVars(0), 
                       _numSlack(0), _numArtificial(0) {}
    
    /// @brief Initialize tableau from standard form
    void Initialize(const Matrix<Real>& A, const Vector<Real>& b, const Vector<Real>& c,
                    const std::vector<int>& artificialIndices,
                    int numOriginalVars, int numSlack, int numArtificial,
                    const LPConfig& config = LPConfig()) {
        if (config.pivotRule == LPPivotRule::Steepest)
            throw NotImplementedError("SimplexTableau: LPPivotRule::Steepest is not implemented");
        _config = config;
        _numRows = A.rows();
        _numOriginalVars = numOriginalVars;
        _numSlack = numSlack;
        _numArtificial = numArtificial;
        _numCols = A.cols() + 1;  // +1 for RHS column
        _artificialIndices = artificialIndices;
        
        // Tableau: [A | b] with objective row at bottom
        _tableau = Matrix<Real>(_numRows + 1, _numCols);
        
        // Fill constraint rows
        for (int i = 0; i < _numRows; ++i) {
            for (int j = 0; j < A.cols(); ++j) {
                _tableau(i, j) = A(i, j);
            }
            _tableau(i, _numCols - 1) = b[i];  // RHS
        }
        
        // Fill objective row (bottom)
        for (int j = 0; j < c.size(); ++j) {
            _tableau(_numRows, j) = c[j];
        }
        _tableau(_numRows, _numCols - 1) = 0;  // Initial objective value
        
        // Initial basis: Find the unit column for each row
        // The basic variable for row i is the column j where A(i,j)=1 and A(k,j)=0 for k!=i
        _basis.resize(_numRows);
        
        for (int i = 0; i < _numRows; ++i) {
            _basis[i] = -1;  // Not found yet
            
            // First check artificial variables (they take priority for two-phase)
            for (int artIdx : artificialIndices) {
                if (std::abs(_tableau(i, artIdx) - 1.0) < 1e-10) {
                    // Verify it's a unit column (1 in row i, 0 elsewhere)
                    bool isUnit = true;
                    for (int k = 0; k < _numRows && isUnit; ++k) {
                        if (k != i && std::abs(_tableau(k, artIdx)) > 1e-10) {
                            isUnit = false;
                        }
                    }
                    if (isUnit) {
                        _basis[i] = artIdx;
                        break;
                    }
                }
            }
            
            if (_basis[i] >= 0) continue;  // Found artificial
            
            // Look for slack variable (unit column with +1, not -1 like surplus)
            int slackStart = numOriginalVars;
            int slackEnd = numOriginalVars + numSlack;
            for (int j = slackStart; j < slackEnd; ++j) {
                if (std::abs(_tableau(i, j) - 1.0) < 1e-10) {
                    // Verify it's a unit column
                    bool isUnit = true;
                    for (int k = 0; k < _numRows && isUnit; ++k) {
                        if (k != i && std::abs(_tableau(k, j)) > 1e-10) {
                            isUnit = false;
                        }
                    }
                    if (isUnit) {
                        _basis[i] = j;
                        break;
                    }
                }
            }
        }
    }
    
    //---------------------------------------------------------------------------------
    // Accessors
    //---------------------------------------------------------------------------------
    
    int NumRows() const { return _numRows; }
    int NumCols() const { return _numCols; }
    int NumVariables() const { return _numCols - 1; }  // Excluding RHS
    const std::vector<int>& Basis() const { return _basis; }
    int Iterations() const { return _iterations; }
    
    Real GetElement(int i, int j) const { return _tableau(i, j); }
    Real GetRHS(int i) const { return _tableau(i, _numCols - 1); }
    Real GetObjectiveValue() const { return -_tableau(_numRows, _numCols - 1); }
    Real GetReducedCost(int j) const { return _tableau(_numRows, j); }
    
    bool IsBasic(int j) const {
        return std::find(_basis.begin(), _basis.end(), j) != _basis.end();
    }
    
    bool IsArtificial(int j) const {
        return std::find(_artificialIndices.begin(), _artificialIndices.end(), j) 
               != _artificialIndices.end();
    }
    
    //---------------------------------------------------------------------------------
    // Pivot Operations
    //---------------------------------------------------------------------------------
    
    /// @brief Select entering variable (column to enter basis)
    /// @return Column index, or -1 if optimal
    int SelectEnteringVariable() const {
        int entering = -1;
        Real minCost = -_config.tolerance;  // Must be negative
        
        if (_config.pivotRule == LPPivotRule::Bland) {
            // Bland's rule: smallest index with negative reduced cost
            for (int j = 0; j < _numCols - 1; ++j) {
                // Skip artificial variables in Phase 2
                if (_phase == 2 && IsArtificial(j)) continue;
                
                if (_tableau(_numRows, j) < -_config.tolerance) {
                    return j;
                }
            }
        } else {
            // Dantzig's rule: most negative reduced cost
            for (int j = 0; j < _numCols - 1; ++j) {
                // Skip artificial variables in Phase 2
                if (_phase == 2 && IsArtificial(j)) continue;
                
                if (_tableau(_numRows, j) < minCost) {
                    minCost = _tableau(_numRows, j);
                    entering = j;
                }
            }
        }
        
        return entering;
    }
    
    /// @brief Select leaving variable (row to leave basis)
    /// @param entering Entering column index
    /// @return Row index, or -1 if unbounded
    int SelectLeavingVariable(int entering) const {
        int leaving = -1;
        Real minRatio = std::numeric_limits<Real>::infinity();
        
        for (int i = 0; i < _numRows; ++i) {
            Real aij = _tableau(i, entering);
            if (aij > _config.pivotTolerance) {
                Real ratio = _tableau(i, _numCols - 1) / aij;
                if (ratio >= 0 && ratio < minRatio) {
                    minRatio = ratio;
                    leaving = i;
                } else if (std::abs(ratio - minRatio) < _config.tolerance && 
                           _config.pivotRule == LPPivotRule::Bland) {
                    // Bland's rule: smallest index for tie-breaking
                    if (_basis[i] < _basis[leaving]) {
                        leaving = i;
                    }
                }
            }
        }
        
        return leaving;
    }
    
    /// @brief Perform pivot operation
    /// @param pivotRow Leaving variable row
    /// @param pivotCol Entering variable column
    void Pivot(int pivotRow, int pivotCol) {
        Real pivot = _tableau(pivotRow, pivotCol);
        
        if (std::abs(pivot) < _config.pivotTolerance) {
            throw LinearProgrammingError("Pivot element too small: " + std::to_string(pivot));
        }
        
        // Normalize pivot row
        for (int j = 0; j < _numCols; ++j) {
            _tableau(pivotRow, j) /= pivot;
        }
        
        // Eliminate column in other rows
        for (int i = 0; i <= _numRows; ++i) {  // Include objective row
            if (i != pivotRow) {
                Real factor = _tableau(i, pivotCol);
                for (int j = 0; j < _numCols; ++j) {
                    _tableau(i, j) -= factor * _tableau(pivotRow, j);
                }
            }
        }
        
        // Update basis
        _basis[pivotRow] = pivotCol;
        _iterations++;
    }
    
    //---------------------------------------------------------------------------------
    // Dual Simplex Operations
    //---------------------------------------------------------------------------------
    
    /// @brief Check if current solution is primal feasible (all RHS >= 0)
    bool IsPrimalFeasible() const {
        for (int i = 0; i < _numRows; ++i) {
            if (_tableau(i, _numCols - 1) < -_config.tolerance) {
                return false;
            }
        }
        return true;
    }
    
    /// @brief Check if current solution is dual feasible (all reduced costs >= 0 for minimization)
    bool IsDualFeasible() const {
        for (int j = 0; j < _numCols - 1; ++j) {
            // Skip artificial variables in Phase 2
            if (_phase == 2 && IsArtificial(j)) continue;
            
            if (_tableau(_numRows, j) < -_config.tolerance) {
                return false;
            }
        }
        return true;
    }
    
    /// @brief Select leaving variable for dual simplex (row with most negative RHS)
    /// @return Row index, or -1 if optimal (all RHS non-negative)
    int SelectDualLeavingVariable() const {
        int leaving = -1;
        Real minRHS = -_config.tolerance;  // Must be negative
        
        if (_config.pivotRule == LPPivotRule::Bland) {
            // Bland's rule: smallest index with negative RHS
            for (int i = 0; i < _numRows; ++i) {
                if (_tableau(i, _numCols - 1) < -_config.tolerance) {
                    return i;
                }
            }
        } else {
            // Most negative RHS (dual Dantzig's rule)
            for (int i = 0; i < _numRows; ++i) {
                if (_tableau(i, _numCols - 1) < minRHS) {
                    minRHS = _tableau(i, _numCols - 1);
                    leaving = i;
                }
            }
        }
        
        return leaving;
    }
    
    /// @brief Select entering variable for dual simplex
    /// @param leaving Leaving row index
    /// @return Column index, or -1 if dual unbounded (primal infeasible)
    int SelectDualEnteringVariable(int leaving) const {
        int entering = -1;
        Real minRatio = std::numeric_limits<Real>::infinity();
        
        for (int j = 0; j < _numCols - 1; ++j) {
            // Skip artificial variables in Phase 2
            if (_phase == 2 && IsArtificial(j)) continue;
            
            Real aij = _tableau(leaving, j);
            
            // Need negative pivot element for dual simplex
            if (aij < -_config.pivotTolerance) {
                // Dual ratio test: c_j / |a_{ij}| = -c_j / a_{ij}
                Real cj = _tableau(_numRows, j);
                Real ratio = -cj / aij;
                
                if (ratio >= 0 && ratio < minRatio) {
                    minRatio = ratio;
                    entering = j;
                } else if (std::abs(ratio - minRatio) < _config.tolerance &&
                           _config.pivotRule == LPPivotRule::Bland) {
                    // Bland's rule: smallest index for tie-breaking
                    if (j < entering) {
                        entering = j;
                    }
                }
            }
        }
        
        return entering;
    }
    
    //---------------------------------------------------------------------------------
    // Solution Extraction
    //---------------------------------------------------------------------------------
    
    /// @brief Extract primal solution
    Vector<Real> GetPrimalSolution() const {
        Vector<Real> x(_numCols - 1);  // Exclude RHS
        for (int i = 0; i < _numRows; ++i) {
            x[_basis[i]] = _tableau(i, _numCols - 1);
        }
        return x;
    }
    
    /// @brief Extract solution for original variables only
    Vector<Real> GetOriginalSolution() const {
        Vector<Real> x(_numOriginalVars);
        for (int i = 0; i < _numRows; ++i) {
            if (_basis[i] < _numOriginalVars) {
                x[_basis[i]] = _tableau(i, _numCols - 1);
            }
        }
        return x;
    }
    
    /// @brief Get dual values (shadow prices)
    Vector<Real> GetDualValues() const {
        // Dual values are the reduced costs of slack variables
        // (with sign adjustment)
        Vector<Real> dual(_numRows);
        for (int i = 0; i < _numRows; ++i) {
            int slackCol = _numOriginalVars + i;
            if (slackCol < _numCols - 1) {
                dual[i] = -_tableau(_numRows, slackCol);
            }
        }
        return dual;
    }
    
    /// @brief Get reduced costs for all variables
    Vector<Real> GetReducedCosts() const {
        Vector<Real> rc(_numCols - 1);
        for (int j = 0; j < _numCols - 1; ++j) {
            rc[j] = _tableau(_numRows, j);
        }
        return rc;
    }
    
    /// @brief Check if any artificial variable is in basis with positive value
    bool HasArtificialInBasis() const {
        for (int i = 0; i < _numRows; ++i) {
            if (IsArtificial(_basis[i])) {
                if (_tableau(i, _numCols - 1) > _config.tolerance) {
                    return true;
                }
            }
        }
        return false;
    }
    
    //---------------------------------------------------------------------------------
    // Phase 1 Setup
    //---------------------------------------------------------------------------------
    
    /// @brief Set up Phase 1 objective (minimize sum of artificials)
    void SetupPhase1() {
        _phase = 1;
        
        // Save original objective
        Vector<Real> originalObj(_numCols);
        for (int j = 0; j < _numCols; ++j) {
            originalObj[j] = _tableau(_numRows, j);
        }
        
        // Phase 1 objective: minimize sum of artificial variables
        for (int j = 0; j < _numCols; ++j) {
            _tableau(_numRows, j) = 0;
        }
        for (int artIdx : _artificialIndices) {
            _tableau(_numRows, artIdx) = 1.0;
        }
        
        // Make artificial variables non-basic in objective row
        // (row operations to create proper reduced costs)
        for (int i = 0; i < _numRows; ++i) {
            if (IsArtificial(_basis[i])) {
                // Subtract this row from objective row
                for (int j = 0; j < _numCols; ++j) {
                    _tableau(_numRows, j) -= _tableau(i, j);
                }
            }
        }
    }
    
    /// @brief Restore original objective for Phase 2
    void SetupPhase2(const Vector<Real>& originalC) {
        _phase = 2;
        
        // Reset objective row
        for (int j = 0; j < _numCols; ++j) {
            _tableau(_numRows, j) = 0;
        }
        for (int j = 0; j < originalC.size() && j < _numCols - 1; ++j) {
            _tableau(_numRows, j) = originalC[j];
        }
        
        // Make basic variables non-basic in objective row
        for (int i = 0; i < _numRows; ++i) {
            int basisCol = _basis[i];
            if (std::abs(_tableau(_numRows, basisCol)) > _config.tolerance) {
                Real factor = _tableau(_numRows, basisCol);
                for (int j = 0; j < _numCols; ++j) {
                    _tableau(_numRows, j) -= factor * _tableau(i, j);
                }
            }
        }
    }
    
    //---------------------------------------------------------------------------------
    // Sensitivity Analysis
    //---------------------------------------------------------------------------------
    
    /// @brief Get allowable range for an objective coefficient
    /// @param varIndex Variable index (original variable)
    /// @return Pair (lower_bound, upper_bound) for the coefficient
    std::pair<Real, Real> GetObjectiveCoefficientRange(int varIndex) const {
        Real lower = -std::numeric_limits<Real>::infinity();
        Real upper = std::numeric_limits<Real>::infinity();
        
        // Check if variable is basic
        int basicRow = -1;
        for (int i = 0; i < _numRows; ++i) {
            if (_basis[i] == varIndex) {
                basicRow = i;
                break;
            }
        }
        
        if (basicRow >= 0) {
            // Basic variable: coefficient change affects reduced costs of non-basic variables
            for (int j = 0; j < _numCols - 1; ++j) {
                if (!IsBasic(j) && !IsArtificial(j)) {
                    Real aij = _tableau(basicRow, j);
                    Real cj = _tableau(_numRows, j);
                    
                    if (std::abs(aij) > _config.tolerance) {
                        // Change delta must keep c_j - delta * a_{ij} >= 0
                        if (aij > 0) {
                            upper = std::min(upper, cj / aij);
                        } else {
                            lower = std::max(lower, cj / aij);
                        }
                    }
                }
            }
        } else {
            // Non-basic variable: reduced cost must remain non-negative
            Real cj = _tableau(_numRows, varIndex);
            lower = -cj;  // Current reduced cost is the margin
        }
        
        return {lower, upper};
    }
    
    /// @brief Get allowable range for a RHS value (constraint bound)
    /// @param constraintIndex Constraint index
    /// @return Pair (lower_bound, upper_bound) for the RHS
    std::pair<Real, Real> GetRHSRange(int constraintIndex) const {
        Real lower = -std::numeric_limits<Real>::infinity();
        Real upper = std::numeric_limits<Real>::infinity();
        
        Real currentRHS = _tableau(constraintIndex, _numCols - 1);
        
        // Change in RHS affects the solution through B^(-1)
        // The basic solution changes as: x_B = B^(-1)b, so x_B_i += delta * B^(-1)_{i,row}
        for (int i = 0; i < _numRows; ++i) {
            // The i-th element of B^(-1) column for this constraint
            // In the tableau, this is stored in the slack column for this constraint
            int slackCol = _numOriginalVars + constraintIndex;
            if (slackCol < _numCols - 1) {
                Real coef = _tableau(i, slackCol);
                Real xi = _tableau(i, _numCols - 1);
                
                if (std::abs(coef) > _config.tolerance) {
                    // Need x_B_i + delta * coef >= 0
                    if (coef > 0) {
                        lower = std::max(lower, -xi / coef);
                    } else {
                        upper = std::min(upper, -xi / coef);
                    }
                }
            }
        }
        
        // Convert to absolute bounds
        lower = currentRHS + lower;
        upper = currentRHS + upper;
        
        return {lower, upper};
    }
    
    /// @brief Get shadow price for a constraint
    Real GetShadowPrice(int constraintIndex) const {
        // Shadow price is the dual value (reduced cost of slack variable)
        int slackCol = _numOriginalVars + constraintIndex;
        if (slackCol < _numCols - 1) {
            return -_tableau(_numRows, slackCol);
        }
        return 0;
    }
    
    /// @brief Compute basis condition number (measure of numerical stability)
    Real ComputeBasisCondition() const {
        // Approximate condition using max/min absolute values in basis columns
        Real maxVal = 0;
        Real minVal = std::numeric_limits<Real>::infinity();
        
        for (int i = 0; i < _numRows; ++i) {
            int bCol = _basis[i];
            for (int k = 0; k < _numRows; ++k) {
                Real absVal = std::abs(_tableau(k, bCol));
                if (absVal > _config.tolerance) {
                    maxVal = std::max(maxVal, absVal);
                    minVal = std::min(minVal, absVal);
                }
            }
        }
        
        if (minVal > _config.tolerance && minVal != std::numeric_limits<Real>::infinity()) {
            return maxVal / minVal;
        }
        return std::numeric_limits<Real>::infinity();
    }
    
    //---------------------------------------------------------------------------------
    // Debug Output
    //---------------------------------------------------------------------------------
    
    void Print(std::ostream& os = std::cout) const {
        os << std::fixed << std::setprecision(4);
        os << "Simplex Tableau (Iteration " << _iterations << "):\n";
        os << "Basis: ";
        for (int b : _basis) os << b << " ";
        os << "\n\n";
        
        // Header
        os << std::setw(6) << "Basis";
        for (int j = 0; j < _numCols - 1; ++j) {
            os << std::setw(10) << ("x" + std::to_string(j));
        }
        os << std::setw(10) << "RHS" << "\n";
        os << std::string(6 + 10 * _numCols, '-') << "\n";
        
        // Constraint rows
        for (int i = 0; i < _numRows; ++i) {
            os << std::setw(6) << ("x" + std::to_string(_basis[i]));
            for (int j = 0; j < _numCols; ++j) {
                os << std::setw(10) << _tableau(i, j);
            }
            os << "\n";
        }
        
        // Objective row
        os << std::string(6 + 10 * _numCols, '-') << "\n";
        os << std::setw(6) << "z";
        for (int j = 0; j < _numCols; ++j) {
            os << std::setw(10) << _tableau(_numRows, j);
        }
        os << "\n\n";
    }
};

} // namespace MML::Optimization
#endif // MML_LP_SIMPLEX_TABLEAU_H

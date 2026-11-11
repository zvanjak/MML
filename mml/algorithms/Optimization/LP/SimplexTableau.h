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
#include <mml/core/LinAlgEqSolvers/LinAlgSVD.h>

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
    Matrix<Real> _originalA;     ///< Standard-form columns before pivoting
    Matrix<Real> _rowTransform;  ///< Constraint row operations, B^(-1)
    Vector<Real> _originalRHS;
    Vector<Real> _objective;
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
        const int rows = A.rows();
        const int cols = A.cols();
        if (rows <= 0 || cols <= 0 || b.size() != rows || c.size() != cols ||
            numOriginalVars <= 0 || numSlack < 0 || numArtificial < 0 ||
            numOriginalVars > cols || numSlack > cols - numOriginalVars ||
            numArtificial != cols - numOriginalVars - numSlack)
            throw LinearProgrammingError("Tableau dimensions or variable counts are inconsistent");
        if (rows == Matrix<Real>::MaxDimension() || cols == Matrix<Real>::MaxDimension() ||
            static_cast<size_t>(rows + 1) > Matrix<Real>::MaxElementCount() / static_cast<size_t>(cols + 1))
            throw LinearProgrammingError("Tableau dimensions exceed matrix storage limits");
        if (artificialIndices.size() != static_cast<size_t>(numArtificial))
            throw LinearProgrammingError("Artificial index count does not match the tableau");

        const int artificialStart = numOriginalVars + numSlack;
        std::vector<bool> seen(numArtificial, false);
        for (int index : artificialIndices) {
            if (index < artificialStart || index >= cols)
                throw LinearProgrammingError("Artificial index is outside the artificial column range");
            if (seen[index - artificialStart])
                throw LinearProgrammingError("Artificial indices must be unique");
            seen[index - artificialStart] = true;
        }
        for (int i = 0; i < rows; ++i) {
            if (!std::isfinite(b[i]))
                throw LinearProgrammingError("Tableau RHS must be finite");
            for (int j = 0; j < cols; ++j) {
                if (!std::isfinite(A(i, j)))
                    throw LinearProgrammingError("Tableau coefficients must be finite");
            }
        }
        for (int j = 0; j < cols; ++j) {
            if (!std::isfinite(c[j]))
                throw LinearProgrammingError("Tableau objective must be finite");
        }

        const auto isUnitColumn = [&](int row, int col) {
            if (std::abs(A(row, col) - Real(1)) >= Precision::EqualityTolerance)
                return false;
            for (int k = 0; k < rows; ++k) {
                if (k != row && std::abs(A(k, col)) > Precision::EqualityTolerance)
                    return false;
            }
            return true;
        };
        std::vector<int> basis(rows, -1);
        for (int i = 0; i < rows; ++i) {
            // Artificial unit columns take priority for two-phase simplex.
            for (int artIdx : artificialIndices) {
                if (isUnitColumn(i, artIdx)) {
                    basis[i] = artIdx;
                    break;
                }
            }
            if (basis[i] < 0) {
                for (int j = numOriginalVars; j < artificialStart; ++j) {
                    if (isUnitColumn(i, j)) {
                        basis[i] = j;
                        break;
                    }
                }
            }
            if (basis[i] < 0)
                throw LinearProgrammingError("Tableau row has no auxiliary unit column for its initial basis");
        }

        auto indices = artificialIndices;
        Matrix<Real> originalA(A);
        Matrix<Real> rowTransform = Matrix<Real>::Identity(rows);
        Vector<Real> originalRHS(b), objective(c);
        Matrix<Real> tableau(rows + 1, cols + 1);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j)
                tableau(i, j) = A(i, j);
            tableau(i, cols) = b[i];
        }
        for (int j = 0; j < cols; ++j)
            tableau(rows, j) = c[j];

        _tableau = std::move(tableau);
        _originalA = std::move(originalA);
        _rowTransform = std::move(rowTransform);
        _originalRHS = std::move(originalRHS);
        _objective = std::move(objective);
        _basis = std::move(basis);
        _artificialIndices = std::move(indices);
        _config = config;
        _numRows = rows;
        _numCols = cols + 1;
        _numOriginalVars = numOriginalVars;
        _numSlack = numSlack;
        _numArtificial = numArtificial;
        _iterations = 0;
        _phase = 0;
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
        for (int j = 0; j < _numRows; ++j)
            _rowTransform(pivotRow, j) /= pivot;
        
        // Eliminate column in other rows
        for (int i = 0; i <= _numRows; ++i) {  // Include objective row
            if (i != pivotRow) {
                Real factor = _tableau(i, pivotCol);
                for (int j = 0; j < _numCols; ++j) {
                    _tableau(i, j) -= factor * _tableau(pivotRow, j);
                }
                if (i < _numRows)
                    for (int j = 0; j < _numRows; ++j)
                        _rowTransform(i, j) -= factor * _rowTransform(pivotRow, j);
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
        Vector<Real> dual(_numRows);
        for (int i = 0; i < _numRows; ++i)
            dual[i] = GetShadowPrice(i);
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
        Vector<Real> objective(_numCols - 1);
        for (int j = 0; j < originalC.size() && j < _numCols - 1; ++j)
            objective[j] = originalC[j];
        _objective = std::move(objective);
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
        if (varIndex < 0 || varIndex >= _numOriginalVars)
            throw LinearProgrammingError("Objective sensitivity index is outside the original variables");
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
                    
                    if (aij != 0) {
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
        
        return {_objective[varIndex] + lower, _objective[varIndex] + upper};
    }
    
    /// @brief Get allowable range for a RHS value (constraint bound)
    /// @param constraintIndex Constraint index
    /// @return Pair (lower_bound, upper_bound) for the RHS
    std::pair<Real, Real> GetRHSRange(int constraintIndex) const {
        if (constraintIndex < 0 || constraintIndex >= _numRows)
            throw LinearProgrammingError("RHS sensitivity index is outside the original rows");
        Real lower = -std::numeric_limits<Real>::infinity();
        Real upper = std::numeric_limits<Real>::infinity();

        // x_B(delta) = B^(-1)b + delta B^(-1)e_j, with b in the input row coordinates.
        for (int i = 0; i < _numRows; ++i) {
            const Real coef = _rowTransform(i, constraintIndex);
            const Real xi = _tableau(i, _numCols - 1);
            if (!std::isfinite(coef) || !std::isfinite(xi))
                throw LinearProgrammingError("Cannot compute RHS sensitivity from a non-finite basis transform");
            if (IsArtificial(_basis[i])) {
                // A zero artificial left by a redundant row must remain zero, not merely nonnegative.
                if (coef != 0) {
                    lower = std::max(lower, -xi / coef);
                    upper = std::min(upper, -xi / coef);
                }
            } else {
                if (coef > 0) lower = std::max(lower, -xi / coef);
                if (coef < 0) upper = std::min(upper, -xi / coef);
            }
        }

        return {_originalRHS[constraintIndex] + lower, _originalRHS[constraintIndex] + upper};
    }
    
    /// @brief Get shadow price for a constraint
    Real GetShadowPrice(int constraintIndex) const {
        if (constraintIndex < 0 || constraintIndex >= _numRows)
            throw LinearProgrammingError("Shadow price index is outside the original rows");
        Real value = 0;
        for (int i = 0; i < _numRows; ++i)
            value += _objective[_basis[i]] * _rowTransform(i, constraintIndex);
        if (!std::isfinite(value))
            throw LinearProgrammingError("Shadow price is not representable");
        return value;
    }
    
    /// @brief Compute the spectral condition number of the input basis, not its canonical columns
    Real ComputeBasisCondition() const {
        if (_numRows == 0)
            throw LinearProgrammingError("Basis conditioning requires an initialized tableau");
        Matrix<Real> basis(_numRows, _numRows);
        Real scale = 0;
        for (int i = 0; i < _numRows; ++i)
            for (int j = 0; j < _numRows; ++j) {
                basis(i, j) = _originalA(i, _basis[j]);
                scale = std::max(scale, std::abs(basis(i, j)));
            }
        if (scale == 0) return std::numeric_limits<Real>::infinity();
        // Conditioning is scale invariant; the unscaled singular values may overflow.
        for (int i = 0; i < _numRows; ++i)
            for (int j = 0; j < _numRows; ++j)
                basis(i, j) /= scale;
        const auto values = SVDecompositionSolver<Real>(basis).getW();
        if (values[_numRows - 1] <= 0) return std::numeric_limits<Real>::infinity();
        return values[0] / values[_numRows - 1];
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

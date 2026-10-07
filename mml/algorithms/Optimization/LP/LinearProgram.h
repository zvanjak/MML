///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Optimization/LP/LinearProgram.h                                             ///
///  Description: Dense continuous linear programming model and standard-form conversion   ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_LP_LINEAR_PROGRAM_H
#define MML_LP_LINEAR_PROGRAM_H

#include <mml/algorithms/Optimization/LP/LPTypes.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>

namespace MML::Optimization {

class LinearProgram {
private:
    // Problem data
    Vector<Real> _c;                     ///< Objective coefficients
    std::vector<LPConstraint> _constraints;
    std::vector<LPVariable> _variables;
    LPObjective _objective = LPObjective::Minimize;
    std::string _name;
    
    // Derived dimensions
    int _numVars = 0;                    ///< Number of original variables
    int _numConstraints = 0;             ///< Number of constraints
    
public:
    //---------------------------------------------------------------------------------
    // Constructors
    //---------------------------------------------------------------------------------
    
    LinearProgram() = default;
    
    /// @brief Create LP with given number of variables
    explicit LinearProgram(int numVars, const std::string& name = "") 
        : _numVars(numVars)
        , _name(name) {
        _c.Resize(numVars);
        _variables.resize(numVars);
        for (int i = 0; i < numVars; ++i) {
            _variables[i] = LPVariable("x" + std::to_string(i + 1));
        }
    }
    
    /// @brief Create LP from objective and constraint matrix
    /// @param c Objective coefficients
    /// @param A Constraint matrix
    /// @param b Right-hand side
    /// @param types Constraint types (default all <=)
    LinearProgram(const Vector<Real>& c, 
                  const Matrix<Real>& A, 
                  const Vector<Real>& b,
                  const std::vector<LPConstraintType>& types = {})
        : _c(c)
        , _numVars(c.size())
        , _numConstraints(b.size()) {
        
        if (A.rows() != b.size() || A.cols() != c.size()) {
            throw LinearProgrammingError("Dimension mismatch in LP construction");
        }
        
        // Create variables
        _variables.resize(_numVars);
        for (int i = 0; i < _numVars; ++i) {
            _variables[i] = LPVariable("x" + std::to_string(i + 1));
        }
        
        // Create constraints
        _constraints.reserve(_numConstraints);
        for (int i = 0; i < _numConstraints; ++i) {
            Vector<Real> row(_numVars);
            for (int j = 0; j < _numVars; ++j) {
                row[j] = A(i, j);
            }
            LPConstraintType type = (i < static_cast<int>(types.size())) 
                                    ? types[i] 
                                    : LPConstraintType::LessEqual;
            _constraints.emplace_back(row, type, b[i]);
        }
    }
    
    //---------------------------------------------------------------------------------
    // Problem Building Interface
    //---------------------------------------------------------------------------------
    
    /// @brief Set objective function
    void SetObjective(const Vector<Real>& c, LPObjective sense = LPObjective::Minimize) {
        _c = c;
        _numVars = c.size();
        _objective = sense;
        
        // Ensure variables exist
        while (_variables.size() < static_cast<size_t>(_numVars)) {
            _variables.emplace_back("x" + std::to_string(_variables.size() + 1));
        }
    }
    
    /// @brief Set objective function from initializer list
    void SetObjective(std::initializer_list<Real> c, LPObjective sense = LPObjective::Minimize) {
        SetObjective(Vector<Real>(c), sense);
    }
    
    /// @brief Set objective coefficient for variable i
    void SetObjectiveCoeff(int i, Real value) {
        if (i < 0 || i >= _numVars) {
            throw LinearProgrammingError("Variable index out of range");
        }
        _c[i] = value;
    }
    
    /// @brief Set objective sense (min/max)
    void SetObjectiveSense(LPObjective sense) { _objective = sense; }
    
    /// @brief Add a constraint
    void AddConstraint(const Vector<Real>& coeffs, LPConstraintType type, Real rhs, 
                       const std::string& name = "") {
        if (coeffs.size() != _numVars && _numVars > 0) {
            throw LinearProgrammingError("Constraint dimension mismatch");
        }
        if (_numVars == 0) {
            _numVars = coeffs.size();
            _c.Resize(_numVars);
            _variables.resize(_numVars);
            for (int i = 0; i < _numVars; ++i) {
                _variables[i] = LPVariable("x" + std::to_string(i + 1));
            }
        }
        _constraints.emplace_back(coeffs, type, rhs, name);
        _numConstraints++;
    }
    
    /// @brief Add constraint using initializer list
    void AddConstraint(std::initializer_list<Real> coeffs, LPConstraintType type, Real rhs,
                       const std::string& name = "") {
        AddConstraint(Vector<Real>(coeffs), type, rhs, name);
    }
    
    /// @brief Add constraint using std::vector (for programmatic construction)
    void AddConstraint(const std::vector<Real>& coeffs, LPConstraintType type, Real rhs,
                       const std::string& name = "") {
        AddConstraint(Vector<Real>(coeffs), type, rhs, name);
    }
    
    /// @brief Set variable bounds metadata.
    /// @note The dense core simplex solver currently supports x >= 0 directly. Model finite upper/lower bounds as explicit constraints.
    void SetVariableBounds(int i, Real lb, Real ub) {
        if (i < 0 || i >= _numVars) {
            throw LinearProgrammingError("Variable index out of range");
        }
        if (lb != REAL(0.0) || ub != std::numeric_limits<Real>::infinity()) {
            throw LinearProgrammingError("Finite/free variable bounds are not transformed by the core LP solver; add them as explicit constraints");
        }
        _variables[i].lowerBound = lb;
        _variables[i].upperBound = ub;
    }
    
    /// @brief Set variable name
    void SetVariableName(int i, const std::string& name) {
        if (i < 0 || i >= _numVars) {
            throw LinearProgrammingError("Variable index out of range");
        }
        _variables[i].name = name;
    }
    
    /// @brief Set all variable names
    void SetVariableNames(const std::vector<std::string>& names) {
        for (size_t i = 0; i < names.size() && i < _variables.size(); ++i) {
            _variables[i].name = names[i];
        }
    }
    
    //---------------------------------------------------------------------------------
    // Accessors
    //---------------------------------------------------------------------------------
    
    int NumVariables() const { return _numVars; }
    int NumConstraints() const { return _numConstraints; }
    const Vector<Real>& ObjectiveCoeffs() const { return _c; }
    LPObjective ObjectiveSense() const { return _objective; }
    const std::vector<LPConstraint>& Constraints() const { return _constraints; }
    const std::vector<LPVariable>& Variables() const { return _variables; }
    const LPConstraint& GetConstraint(int i) const { return _constraints[i]; }
    const LPVariable& GetVariable(int i) const { return _variables[i]; }
    const std::string& Name() const { return _name; }
    void SetName(const std::string& name) { _name = name; }
    
    //---------------------------------------------------------------------------------
    // Standard Form Conversion
    //---------------------------------------------------------------------------------
    
    /// @brief Get problem in standard form (for simplex)
    /// 
    /// Standard form: min c'x, Ax = b, x >= 0
    /// - <= constraints: add slack variable
    /// - >= constraints: subtract surplus, add artificial
    /// - = constraints: add artificial
    /// - Variables with lb != 0: shift
    /// - Free variables: split into x+ - x-
    ///
    /// @param[out] A_std Standard form constraint matrix
    /// @param[out] b_std Standard form RHS
    /// @param[out] c_std Standard form objective
    /// @param[out] numSlack Number of slack variables added
    /// @param[out] numArtificial Number of artificial variables added
    /// @param[out] artificialIndices Indices of artificial variables
    void ToStandardForm(Matrix<Real>& A_std, Vector<Real>& b_std, Vector<Real>& c_std,
                        int& numSlack, int& numArtificial,
                        std::vector<int>& artificialIndices) const {
        
        // Count slack and artificial variables needed
        numSlack = 0;
        numArtificial = 0;
        
        for (const auto& con : _constraints) {
            switch (con.type) {
                case LPConstraintType::LessEqual:
                    numSlack++;
                    break;
                case LPConstraintType::Equal:
                    numArtificial++;
                    break;
                case LPConstraintType::GreaterEqual:
                    numSlack++;        // Surplus variable
                    numArtificial++;   // Artificial variable
                    break;
            }
        }
        
        int totalVars = _numVars + numSlack + numArtificial;
        
        // Build standard form matrix
        A_std = Matrix<Real>(_numConstraints, totalVars);
        b_std.Resize(_numConstraints);
        c_std.Resize(totalVars);
        artificialIndices.clear();
        
        // Copy original objective (negate if maximizing)
        Real objSign = (_objective == LPObjective::Maximize) ? -1.0 : 1.0;
        for (int j = 0; j < _numVars; ++j) {
            c_std[j] = objSign * _c[j];
        }
        
        // Fill constraint matrix
        int slackIdx = _numVars;
        int artificialIdx = _numVars + numSlack;
        
        for (int i = 0; i < _numConstraints; ++i) {
            const auto& con = _constraints[i];
            
            // Copy original coefficients
            for (int j = 0; j < _numVars; ++j) {
                A_std(i, j) = con.coefficients[j];
            }
            
            // Handle RHS sign (standard form requires b >= 0)
            Real rhsSign = 1.0;
            if (con.rhs < 0) {
                rhsSign = -1.0;
                // Flip constraint: multiply row by -1, flip type
                for (int j = 0; j < _numVars; ++j) {
                    A_std(i, j) = -A_std(i, j);
                }
            }
            b_std[i] = std::abs(con.rhs);
            
            // Adjust constraint type after potential flip
            LPConstraintType effectiveType = con.type;
            if (rhsSign < 0) {
                if (con.type == LPConstraintType::LessEqual) {
                    effectiveType = LPConstraintType::GreaterEqual;
                } else if (con.type == LPConstraintType::GreaterEqual) {
                    effectiveType = LPConstraintType::LessEqual;
                }
            }
            
            // Add slack/surplus/artificial variables
            switch (effectiveType) {
                case LPConstraintType::LessEqual:
                    A_std(i, slackIdx) = 1.0;  // Slack
                    slackIdx++;
                    break;
                    
                case LPConstraintType::Equal:
                    A_std(i, artificialIdx) = 1.0;  // Artificial
                    artificialIndices.push_back(artificialIdx);
                    artificialIdx++;
                    break;
                    
                case LPConstraintType::GreaterEqual:
                    A_std(i, slackIdx) = -1.0;  // Surplus
                    slackIdx++;
                    A_std(i, artificialIdx) = 1.0;  // Artificial
                    artificialIndices.push_back(artificialIdx);
                    artificialIdx++;
                    break;
            }
        }
    }
    
    //---------------------------------------------------------------------------------
    // Dual Problem
    //---------------------------------------------------------------------------------
    
    /// @brief Construct the dual LP
    /// 
    /// Primal: min c'x, Ax >= b, x >= 0
    /// Dual:   max b'y, A'y <= c, y >= 0
    ///
    /// Note: General duality rules apply based on constraint types
    LinearProgram Dual() const {
        // For now, handle standard case: min c'x, Ax <= b, x >= 0
        // Dual: max b'y, A'y <= c, y >= 0
        
        LinearProgram dual(_numConstraints);
        dual._name = _name + "_dual";
        
        // Dual objective is primal RHS
        Vector<Real> dualC(_numConstraints);
        for (int i = 0; i < _numConstraints; ++i) {
            dualC[i] = _constraints[i].rhs;
        }
        dual.SetObjective(dualC, _objective == LPObjective::Minimize 
                                  ? LPObjective::Maximize 
                                  : LPObjective::Minimize);
        
        // Dual constraints from primal columns
        for (int j = 0; j < _numVars; ++j) {
            Vector<Real> dualRow(_numConstraints);
            for (int i = 0; i < _numConstraints; ++i) {
                dualRow[i] = _constraints[i].coefficients[j];
            }
            
            // Constraint type depends on primal variable bounds and objective
            LPConstraintType dualType = LPConstraintType::LessEqual;
            if (_objective == LPObjective::Maximize) {
                dualType = LPConstraintType::GreaterEqual;
            }
            
            dual.AddConstraint(dualRow, dualType, _c[j]);
        }
        
        return dual;
    }
    
    //---------------------------------------------------------------------------------
    // Utility Functions
    //---------------------------------------------------------------------------------
    
    /// @brief Print problem in readable format
    std::string ToString() const {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(4);
        
        if (!_name.empty()) {
            oss << "Problem: " << _name << "\n";
        }
        
        // Objective
        oss << (_objective == LPObjective::Minimize ? "Minimize: " : "Maximize: ");
        bool first = true;
        for (int j = 0; j < _numVars; ++j) {
            if (std::abs(_c[j]) > 1e-10) {
                if (!first && _c[j] > 0) oss << "+ ";
                if (_c[j] < 0) oss << "- ";
                if (std::abs(std::abs(_c[j]) - 1.0) > 1e-10) {
                    oss << std::abs(_c[j]) << "*";
                }
                oss << _variables[j].name << " ";
                first = false;
            }
        }
        oss << "\n\nSubject to:\n";
        
        // Constraints
        for (int i = 0; i < _numConstraints; ++i) {
            const auto& con = _constraints[i];
            oss << "  ";
            if (!con.name.empty()) {
                oss << con.name << ": ";
            }
            
            first = true;
            for (int j = 0; j < _numVars; ++j) {
                if (std::abs(con.coefficients[j]) > 1e-10) {
                    if (!first && con.coefficients[j] > 0) oss << "+ ";
                    if (con.coefficients[j] < 0) oss << "- ";
                    if (std::abs(std::abs(con.coefficients[j]) - 1.0) > 1e-10) {
                        oss << std::abs(con.coefficients[j]) << "*";
                    }
                    oss << _variables[j].name << " ";
                    first = false;
                }
            }
            
            switch (con.type) {
                case LPConstraintType::LessEqual: oss << "<= "; break;
                case LPConstraintType::Equal: oss << "= "; break;
                case LPConstraintType::GreaterEqual: oss << ">= "; break;
            }
            oss << con.rhs << "\n";
        }
        
        // Variable bounds
        oss << "\nBounds:\n";
        for (int j = 0; j < _numVars; ++j) {
            const auto& var = _variables[j];
            if (var.lowerBound == 0 && var.upperBound == std::numeric_limits<Real>::infinity()) {
                oss << "  " << var.name << " >= 0\n";
            } else if (var.IsFree()) {
                oss << "  " << var.name << " free\n";
            } else {
                if (var.lowerBound > -std::numeric_limits<Real>::infinity()) {
                    oss << "  " << var.lowerBound << " <= ";
                }
                oss << var.name;
                if (var.upperBound < std::numeric_limits<Real>::infinity()) {
                    oss << " <= " << var.upperBound;
                }
                oss << "\n";
            }
        }
        
        return oss.str();
    }
    
    /// @brief Validate problem consistency
    bool Validate(std::string& errorMsg) const {
        if (_numVars <= 0) {
            errorMsg = "No variables defined";
            return false;
        }
        
        if (_numConstraints <= 0) {
            errorMsg = "No constraints defined";
            return false;
        }
        
        for (int i = 0; i < _numConstraints; ++i) {
            if (_constraints[i].coefficients.size() != _numVars) {
                errorMsg = "Constraint " + std::to_string(i) + " has wrong dimension";
                return false;
            }
        }
        
        return true;
    }
};


} // namespace MML::Optimization
#endif // MML_LP_LINEAR_PROGRAM_H

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
        : _name(name)
        , _numVars(numVars) {
        _c.resize(numVars);
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
            _c.resize(_numVars);
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
        if (lb != static_cast<MML::Real>(0.0) || ub != std::numeric_limits<Real>::infinity()) {
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
    /// - Negative RHS: negate the row and reverse its inequality before sizing
    /// Variables must have the supported default bounds x >= 0.
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
        std::string error;
        if (!Validate(error))
            throw LinearProgrammingError("Invalid LP: " + error);
        if (&b_std == &c_std || &numSlack == &numArtificial)
            throw LinearProgrammingError("Standard-form outputs must not alias");

        const auto effectiveType = [](const LPConstraint& con) {
            if (con.rhs < 0) {
                if (con.type == LPConstraintType::LessEqual)
                    return LPConstraintType::GreaterEqual;
                if (con.type == LPConstraintType::GreaterEqual)
                    return LPConstraintType::LessEqual;
            }
            return con.type;
        };

        // Size the same normalized rows that will be written below.
        numSlack = 0;
        numArtificial = 0;
        
        for (const auto& con : _constraints) {
            switch (effectiveType(con)) {
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
        
        const int maxDimension = Matrix<Real>::MaxDimension();
        if (numSlack > maxDimension - _numVars || numArtificial > maxDimension - _numVars - numSlack)
            throw LinearProgrammingError("Standard-form variable count exceeds the matrix index range");
        int totalVars = _numVars + numSlack + numArtificial;
        
        // Build standard form matrix
        A_std = Matrix<Real>(_numConstraints, totalVars);
        b_std.resize(_numConstraints);
        c_std.resize(totalVars);
        std::fill(c_std.begin(), c_std.end(), Real(0));
        artificialIndices.clear();
        
        // Copy original objective (negate if maximizing)
        Real objSign = (_objective == LPObjective::Maximize) ? -MML::Real(1.0) : MML::Real(1.0);
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
            if (con.rhs < 0) {
                // Flip constraint: multiply row by -1, flip type
                for (int j = 0; j < _numVars; ++j) {
                    A_std(i, j) = -A_std(i, j);
                }
            }
            b_std[i] = std::abs(con.rhs);
            
            // Add slack/surplus/artificial variables
            switch (effectiveType(con)) {
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
    /// Primal variables must have default bounds x >= 0. For minimization,
    /// <= rows have y <= 0, >= rows have y >= 0, and equality rows have free y.
    /// Maximization reverses the inequality multiplier signs.
    /// The returned LP uses nonnegative coordinates: its first m variables are
    /// the signed inequality multipliers or equality positive parts, followed
    /// by equality negative parts in original row order.
    LinearProgram Dual() const {
        std::string error;
        if (!Validate(error))
            throw LinearProgrammingError("Cannot construct dual: " + error);

        const Real senseSign = _objective == LPObjective::Minimize ? Real(1) : Real(-1);
        const auto multiplierSign = [senseSign](LPConstraintType type) -> Real {
            switch (type) {
                case LPConstraintType::LessEqual: return -senseSign;
                case LPConstraintType::GreaterEqual: return senseSign;
                case LPConstraintType::Equal: return Real(1);
            }
            throw LinearProgrammingError("Cannot construct dual: invalid constraint type");
        };

        int numDualVars = _numConstraints;
        for (const auto& con : _constraints) {
            if (con.type == LPConstraintType::Equal) {
                if (numDualVars == Matrix<Real>::MaxDimension())
                    throw LinearProgrammingError("Dual variable count exceeds the matrix index range");
                ++numDualVars;
            }
        }

        LinearProgram dual(numDualVars);
        dual._name = _name + "_dual";
        
        Vector<Real> dualC(numDualVars);
        int negativePart = _numConstraints;
        for (int i = 0; i < _numConstraints; ++i) {
            const auto& con = _constraints[i];
            const Real sign = multiplierSign(con.type);
            const std::string name = "y" + std::to_string(i + 1);
            dualC[i] = sign * con.rhs;
            dual.SetVariableName(i, name + (con.type == LPConstraintType::Equal ? "_plus" : sign < 0 ? "_minus" : ""));
            if (con.type == LPConstraintType::Equal) {
                dualC[negativePart] = -con.rhs;
                dual.SetVariableName(negativePart++, name + "_minus");
            }
        }
        dual.SetObjective(dualC, _objective == LPObjective::Minimize 
                                  ? LPObjective::Maximize 
                                  : LPObjective::Minimize);
        
        // Dual constraints from primal columns
        for (int j = 0; j < _numVars; ++j) {
            Vector<Real> dualRow(numDualVars);
            negativePart = _numConstraints;
            for (int i = 0; i < _numConstraints; ++i) {
                const auto& con = _constraints[i];
                dualRow[i] = multiplierSign(con.type) * con.coefficients[j];
                if (con.type == LPConstraintType::Equal)
                    dualRow[negativePart++] = -con.coefficients[j];
            }
            
            // Primal variables are nonnegative, so all transposed rows share a sense.
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
    std::string to_string() const {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(4);
        
        if (!_name.empty()) {
            oss << "Problem: " << _name << "\n";
        }
        
        // Objective
        oss << (_objective == LPObjective::Minimize ? "Minimize: " : "Maximize: ");
        bool first = true;
        for (int j = 0; j < _numVars; ++j) {
            if (std::abs(_c[j]) > Precision::EqualityTolerance) {
                if (!first && _c[j] > 0) oss << "+ ";
                if (_c[j] < 0) oss << "- ";
                if (std::abs(std::abs(_c[j]) - 1.0) > Precision::EqualityTolerance) {
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
                if (std::abs(con.coefficients[j]) > Precision::EqualityTolerance) {
                    if (!first && con.coefficients[j] > 0) oss << "+ ";
                    if (con.coefficients[j] < 0) oss << "- ";
                    if (std::abs(std::abs(con.coefficients[j]) - 1.0) > Precision::EqualityTolerance) {
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
    /// @deprecated Use to_string() (CONTRACT_POLICY section 11: STL spelling).
    [[deprecated("use to_string()")]]
    std::string ToString() const { return to_string(); }
    
    /// @brief Validate problem consistency
    bool Validate(std::string& errorMsg) const {
        errorMsg.clear();
        if (_numVars <= 0) {
            errorMsg = "No variables defined";
            return false;
        }
        
        if (_numConstraints <= 0) {
            errorMsg = "No constraints defined";
            return false;
        }
        if (_c.size() != _numVars || _variables.size() < static_cast<size_t>(_numVars) ||
            _constraints.size() != static_cast<size_t>(_numConstraints)) {
            errorMsg = "Inconsistent LP storage dimensions";
            return false;
        }
        if (_objective != LPObjective::Minimize && _objective != LPObjective::Maximize) {
            errorMsg = "Invalid objective sense";
            return false;
        }
        for (int j = 0; j < _numVars; ++j) {
            if (!std::isfinite(_c[j])) {
                errorMsg = "Objective coefficient " + std::to_string(j) + " is non-finite";
                return false;
            }
            if (_variables[j].lowerBound != 0 || _variables[j].upperBound != std::numeric_limits<Real>::infinity()) {
                errorMsg = "Variable " + std::to_string(j) + " must have default bounds x >= 0";
                return false;
            }
        }

        for (int i = 0; i < _numConstraints; ++i) {
            const auto& con = _constraints[i];
            if (con.coefficients.size() != _numVars) {
                errorMsg = "Constraint " + std::to_string(i) + " has wrong dimension";
                return false;
            }
            if (con.type != LPConstraintType::LessEqual && con.type != LPConstraintType::GreaterEqual &&
                con.type != LPConstraintType::Equal) {
                errorMsg = "Constraint " + std::to_string(i) + " has invalid type";
                return false;
            }
            if (!std::isfinite(con.rhs)) {
                errorMsg = "Constraint " + std::to_string(i) + " RHS is non-finite";
                return false;
            }
            for (int j = 0; j < _numVars; ++j) {
                if (!std::isfinite(con.coefficients[j])) {
                    errorMsg = "Constraint " + std::to_string(i) + " coefficient " + std::to_string(j) + " is non-finite";
                    return false;
                }
            }
        }
        
        return true;
    }
};


} // namespace MML::Optimization
#endif // MML_LP_LINEAR_PROGRAM_H

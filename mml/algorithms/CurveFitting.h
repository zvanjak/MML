///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        CurveFitting.h                                                      ///
///  Description: Curve fitting algorithms (linear, polynomial, exponential, spline)  ///
///               Least squares fitting and regression utilities                      ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_CURVE_FITTING_H
#define MML_CURVE_FITTING_H

#include "MMLBase.h"

#include "interfaces/IFunction.h"
#include "base/Vector.h"
#include "base/Matrix.h"
#include "core/LinAlgEqSolvers.h"

#include <functional>

namespace MML
{
    ///////////////////////////           LINEAR LEAST SQUARES           ///////////////////////////
    
    // Fits a linear function y = a*x + b to a set of data points using least squares method
    // 
    // Parameters:
    //   x_data - Vector of x-coordinates (independent variable)
    //   y_data - Vector of y-coordinates (dependent variable)
    //   a      - Output: slope coefficient
    //   b      - Output: intercept coefficient
    //
    // Returns:
    //   residual_norm - The norm of the residual vector (measure of fit quality)
    //
    // Mathematical formulation:
    //   Minimize: sum_i (y_i - (a*x_i + b))^2
    //   
    //   Normal equations: [sum(x_i^2)   sum(x_i)  ] [a]   [sum(x_i*y_i)]
    //                     [sum(x_i)     n         ] [b] = [sum(y_i)    ]
    //
    // Throws:
    //   std::invalid_argument if x_data and y_data have different sizes or are empty
    //   std::runtime_error if the normal equations system is singular (e.g., all x values identical)
    template<typename Real>
    Real LinearLeastSquares(const Vector<Real>& x_data, const Vector<Real>& y_data, 
                           Real& a, Real& b)
    {
        if (x_data.size() != y_data.size())
            throw CurveFittingError("LinearLeastSquares: x_data and y_data must have the same size");
        
        if (x_data.size() == 0)
            throw CurveFittingError("LinearLeastSquares: data vectors cannot be empty");
        
        int n = x_data.size();
        
        // Special case: single point - infinite solutions, use horizontal line through the point
        if (n == 1) {
            a = 0.0;
            b = y_data[0];
            return 0.0;
        }
        
        // Special case: two points - exact fit
        if (n == 2) {
            Real dx = x_data[1] - x_data[0];
            if (std::abs(dx) < std::numeric_limits<Real>::epsilon() * 10) {
                // Vertical line - can't fit y = ax + b, use horizontal line at average y
                a = 0.0;
                b = (y_data[0] + y_data[1]) / 2.0;
                // Compute residual norm for the horizontal line
                Real r0 = y_data[0] - b;
                Real r1 = y_data[1] - b;
                return std::sqrt(r0*r0 + r1*r1);
            }
            a = (y_data[1] - y_data[0]) / dx;
            b = y_data[0] - a * x_data[0];
            return 0.0;
        }
        
        // General case: n >= 3 points - least squares
        // Compute sums for normal equations
        Real sum_x = 0.0, sum_y = 0.0, sum_xx = 0.0, sum_xy = 0.0;
        
        for (int i = 0; i < n; i++) {
            sum_x  += x_data[i];
            sum_y  += y_data[i];
            sum_xx += x_data[i] * x_data[i];
            sum_xy += x_data[i] * y_data[i];
        }
        
        // Normal equations matrix: A = [sum_xx  sum_x ]
        //                              [sum_x   n     ]
        Real det = n * sum_xx - sum_x * sum_x;
        
        // Check for singularity (all x values essentially identical)
        if (std::abs(det) < std::numeric_limits<Real>::epsilon() * std::max(std::abs(n * sum_xx), std::abs(sum_x * sum_x))) {
            throw SingularMatrixError("LinearLeastSquares: singular system - all x values are nearly identical", det);
        }
        
        // Solve using Cramer's rule (explicit for 2x2 system)
        a = (n * sum_xy - sum_x * sum_y) / det;
        b = (sum_xx * sum_y - sum_x * sum_xy) / det;
        
        // Compute residual norm
        Real residual_sum = 0.0;
        for (int i = 0; i < n; i++) {
            Real residual = y_data[i] - (a * x_data[i] + b);
            residual_sum += residual * residual;
        }
        
        return std::sqrt(residual_sum);
    }
    
    
    // Fits a linear function y = a*x + b and returns detailed statistics
    template<typename Real>
    struct LinearFitResult
    {
        Real a;                  // Slope
        Real b;                  // Intercept
        Real residual_norm;      // ||y - (ax + b)||_2
        Real r_squared;          // Coefficient of determination (0 to 1, 1 = perfect fit)
        Real mean_squared_error; // Average squared residual
        
        LinearFitResult() : a(0), b(0), residual_norm(0), r_squared(0), mean_squared_error(0) {}
    };
    
    template<typename Real>
    LinearFitResult<Real> LinearLeastSquaresDetailed(const Vector<Real>& x_data, const Vector<Real>& y_data)
    {
        LinearFitResult<Real> result;
        
        int n = x_data.size();
        result.residual_norm = LinearLeastSquares(x_data, y_data, result.a, result.b);
        
        // Compute R^2 (coefficient of determination)
        Real y_mean = 0.0;
        for (int i = 0; i < n; i++)
            y_mean += y_data[i];
        y_mean /= n;
        
        Real ss_tot = 0.0;  // Total sum of squares
        Real ss_res = 0.0;  // Residual sum of squares
        
        for (int i = 0; i < n; i++) {
            Real y_pred = result.a * x_data[i] + result.b;
            Real residual = y_data[i] - y_pred;
            ss_res += residual * residual;
            ss_tot += (y_data[i] - y_mean) * (y_data[i] - y_mean);
        }
        
        // R^2 = 1 - (SS_res / SS_tot)
        if (ss_tot > std::numeric_limits<Real>::epsilon()) {
            result.r_squared = 1.0 - (ss_res / ss_tot);
        } else {
            // All y values are identical
            result.r_squared = (ss_res < std::numeric_limits<Real>::epsilon()) ? 1.0 : 0.0;
        }
        
        result.mean_squared_error = ss_res / n;
        
        return result;
    }
    
    
    ///////////////////////////    GENERAL LINEAR LEAST SQUARES    ///////////////////////////
    
    // Result structure for general linear least squares fitting
    // Contains fitted coefficients and comprehensive statistics
    template<typename Real>
    struct GeneralLinearFitResult
    {
        Vector<Real> coefficients;     // Fitted coefficients for each basis function
        Real residual_norm;            // ||y - A*c||_2 where A is the design matrix
        Real r_squared;                // Coefficient of determination (0 to 1, 1 = perfect fit)
        Real mean_squared_error;       // Average squared residual
        Real adjusted_r_squared;       // R² adjusted for number of parameters
        int num_data_points;           // Number of data points (m)
        int num_basis_functions;       // Number of basis functions (n)
        int effective_rank;            // Effective rank of the design matrix (from SVD)
        Real condition_number;         // Condition number of the design matrix (ratio of largest/smallest singular values)
        
        GeneralLinearFitResult() 
            : residual_norm(0), r_squared(0), mean_squared_error(0), adjusted_r_squared(0),
              num_data_points(0), num_basis_functions(0), effective_rank(0), condition_number(0) {}
        
        // Evaluate the fitted function at a point x using the basis functions
        Real evaluate(Real x, const Vector<const IRealFunction*>& basis_functions) const
        {
            if (coefficients.size() != basis_functions.size())
                throw CurveFittingError("GeneralLinearFitResult::evaluate: mismatched basis functions count");
            
            Real result = 0.0;
            for (int i = 0; i < coefficients.size(); i++)
                result += coefficients[i] * (*basis_functions[i])(x);
            return result;
        }
        
        // Evaluate using a function wrapper (std::function version)
        Real evaluate(Real x, const Vector<std::function<Real(Real)>>& basis_functions) const
        {
            if (coefficients.size() != basis_functions.size())
                throw CurveFittingError("GeneralLinearFitResult::evaluate: mismatched basis functions count");
            
            Real result = 0.0;
            for (int i = 0; i < coefficients.size(); i++)
                result += coefficients[i] * basis_functions[i](x);
            return result;
        }
    };
    
    // Fits a linear combination of arbitrary basis functions to data using least squares
    //
    // Mathematical formulation:
    //   Find coefficients c_0, c_1, ..., c_{n-1} that minimize:
    //     sum_i (y_i - sum_j c_j * f_j(x_i))^2
    //
    //   This is equivalent to solving the overdetermined system:
    //     A * c = y
    //   where A[i][j] = f_j(x_i) is the design matrix
    //
    //   The solution uses SVD: A = U * W * V^T
    //   Least squares solution: c = V * W^{-1} * U^T * y (pseudoinverse)
    //
    // Parameters:
    //   x_data          - Vector of x-coordinates (independent variable)
    //   y_data          - Vector of y-coordinates (dependent variable)
    //   basis_functions - Vector of pointers to basis functions f_j(x)
    //
    // Returns:
    //   GeneralLinearFitResult containing coefficients and statistics
    //
    // Example basis functions:
    //   - Polynomial: {1, x, x^2, x^3, ...}
    //   - Fourier:    {1, sin(x), cos(x), sin(2x), cos(2x), ...}
    //   - Custom:     {exp(-x^2), log(x+1), sqrt(x), ...}
    //
    // Throws:
    //   std::invalid_argument if vectors have mismatched sizes, are empty, 
    //                         or basis_functions is empty
    template<typename Real>
    GeneralLinearFitResult<Real> GeneralLinearLeastSquares(
        const Vector<Real>& x_data, 
        const Vector<Real>& y_data,
        const Vector<const IRealFunction*>& basis_functions)
    {
        if (x_data.size() != y_data.size())
            throw CurveFittingError("GeneralLinearLeastSquares: x_data and y_data must have the same size");
        
        if (x_data.size() == 0)
            throw CurveFittingError("GeneralLinearLeastSquares: data vectors cannot be empty");
        
        if (basis_functions.size() == 0)
            throw CurveFittingError("GeneralLinearLeastSquares: basis_functions cannot be empty");
        
        int m = x_data.size();                   // Number of data points
        int n = basis_functions.size();          // Number of basis functions
        
        if (m < n)
            throw CurveFittingError("GeneralLinearLeastSquares: need at least as many data points as basis functions");
        
        // Build the design matrix A[i][j] = f_j(x_i)
        Matrix<Real> A(m, n);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                A[i][j] = (*basis_functions[j])(x_data[i]);
            }
        }
        
        // Use SVD to solve the overdetermined system
        SVDecompositionSolver svd(A);
        
        // Solve A*c = y for coefficients c
        Vector<Real> c = svd.Solve(y_data);
        
        // Compute statistics
        GeneralLinearFitResult<Real> result;
        result.coefficients = c;
        result.num_data_points = m;
        result.num_basis_functions = n;
        
        // Compute residuals and fitted values
        Real ss_res = 0.0;  // Residual sum of squares
        for (int i = 0; i < m; i++) {
            Real y_pred = 0.0;
            for (int j = 0; j < n; j++)
                y_pred += c[j] * A[i][j];
            Real residual = y_data[i] - y_pred;
            ss_res += residual * residual;
        }
        
        result.residual_norm = std::sqrt(ss_res);
        result.mean_squared_error = ss_res / m;
        
        // Compute R^2
        Real y_mean = 0.0;
        for (int i = 0; i < m; i++)
            y_mean += y_data[i];
        y_mean /= m;
        
        Real ss_tot = 0.0;  // Total sum of squares
        for (int i = 0; i < m; i++) {
            Real diff = y_data[i] - y_mean;
            ss_tot += diff * diff;
        }
        
        if (ss_tot > std::numeric_limits<Real>::epsilon()) {
            result.r_squared = 1.0 - (ss_res / ss_tot);
            // Adjusted R^2: R^2_adj = 1 - (1 - R^2) * (m - 1) / (m - n)
            if (m > n)
                result.adjusted_r_squared = 1.0 - (1.0 - result.r_squared) * (m - 1.0) / (m - n);
            else
                result.adjusted_r_squared = result.r_squared;
        } else {
            result.r_squared = (ss_res < std::numeric_limits<Real>::epsilon()) ? 1.0 : 0.0;
            result.adjusted_r_squared = result.r_squared;
        }
        
        // Get condition number and rank from SVD
        Vector<Real> singular_values = svd.getW();
        result.effective_rank = 0;
        Real thresh = 0.5 * std::sqrt(m + n + 1.0) * singular_values[0] * std::numeric_limits<Real>::epsilon();
        for (int i = 0; i < n; i++) {
            if (singular_values[i] > thresh)
                result.effective_rank++;
        }
        
        // Condition number = largest singular value / smallest non-zero singular value
        if (result.effective_rank > 0 && singular_values[result.effective_rank - 1] > 0)
            result.condition_number = singular_values[0] / singular_values[result.effective_rank - 1];
        else
            result.condition_number = std::numeric_limits<Real>::infinity();
        
        return result;
    }
    
    // Overload using std::function for more flexible basis function specification
    // This allows lambda functions, function pointers, and functors
    template<typename Real>
    GeneralLinearFitResult<Real> GeneralLinearLeastSquares(
        const Vector<Real>& x_data, 
        const Vector<Real>& y_data,
        const Vector<std::function<Real(Real)>>& basis_functions)
    {
        if (x_data.size() != y_data.size())
            throw CurveFittingError("GeneralLinearLeastSquares: x_data and y_data must have the same size");
        
        if (x_data.size() == 0)
            throw CurveFittingError("GeneralLinearLeastSquares: data vectors cannot be empty");
        
        if (basis_functions.size() == 0)
            throw CurveFittingError("GeneralLinearLeastSquares: basis_functions cannot be empty");
        
        int m = x_data.size();
        int n = basis_functions.size();
        
        if (m < n)
            throw CurveFittingError("GeneralLinearLeastSquares: need at least as many data points as basis functions");
        
        // Build the design matrix
        Matrix<Real> A(m, n);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                A[i][j] = basis_functions[j](x_data[i]);
            }
        }
        
        // Use SVD to solve
        SVDecompositionSolver svd(A);
        Vector<Real> c = svd.Solve(y_data);
        
        // Compute statistics (same as above)
        GeneralLinearFitResult<Real> result;
        result.coefficients = c;
        result.num_data_points = m;
        result.num_basis_functions = n;
        
        Real ss_res = 0.0;
        for (int i = 0; i < m; i++) {
            Real y_pred = 0.0;
            for (int j = 0; j < n; j++)
                y_pred += c[j] * A[i][j];
            Real residual = y_data[i] - y_pred;
            ss_res += residual * residual;
        }
        
        result.residual_norm = std::sqrt(ss_res);
        result.mean_squared_error = ss_res / m;
        
        Real y_mean = 0.0;
        for (int i = 0; i < m; i++)
            y_mean += y_data[i];
        y_mean /= m;
        
        Real ss_tot = 0.0;
        for (int i = 0; i < m; i++) {
            Real diff = y_data[i] - y_mean;
            ss_tot += diff * diff;
        }
        
        if (ss_tot > std::numeric_limits<Real>::epsilon()) {
            result.r_squared = 1.0 - (ss_res / ss_tot);
            if (m > n)
                result.adjusted_r_squared = 1.0 - (1.0 - result.r_squared) * (m - 1.0) / (m - n);
            else
                result.adjusted_r_squared = result.r_squared;
        } else {
            result.r_squared = (ss_res < std::numeric_limits<Real>::epsilon()) ? 1.0 : 0.0;
            result.adjusted_r_squared = result.r_squared;
        }
        
        Vector<Real> singular_values = svd.getW();
        result.effective_rank = 0;
        Real thresh = 0.5 * std::sqrt(m + n + 1.0) * singular_values[0] * std::numeric_limits<Real>::epsilon();
        for (int i = 0; i < n; i++) {
            if (singular_values[i] > thresh)
                result.effective_rank++;
        }
        
        if (result.effective_rank > 0 && singular_values[result.effective_rank - 1] > 0)
            result.condition_number = singular_values[0] / singular_values[result.effective_rank - 1];
        else
            result.condition_number = std::numeric_limits<Real>::infinity();
        
        return result;
    }
    
    
    ///////////////////////////    POLYNOMIAL FITTING    ///////////////////////////
    
    // Fits a polynomial of specified degree to data: y = c_0 + c_1*x + c_2*x^2 + ... + c_n*x^n
    //
    // This is a convenience wrapper around GeneralLinearLeastSquares with polynomial basis functions
    //
    // Parameters:
    //   x_data - Vector of x-coordinates
    //   y_data - Vector of y-coordinates
    //   degree - Degree of the polynomial (0 = constant, 1 = linear, 2 = quadratic, etc.)
    //
    // Returns:
    //   GeneralLinearFitResult where coefficients[i] is the coefficient of x^i
    //
    // Example: degree=2 fits y = c[0] + c[1]*x + c[2]*x^2
    template<typename Real>
    GeneralLinearFitResult<Real> PolynomialFit(
        const Vector<Real>& x_data, 
        const Vector<Real>& y_data,
        int degree)
    {
        if (degree < 0)
            throw CurveFittingError("PolynomialFit: degree must be non-negative");
        
        // Create polynomial basis functions: {1, x, x^2, ..., x^degree}
        Vector<std::function<Real(Real)>> basis(degree + 1);
        
        for (int i = 0; i <= degree; i++) {
            int power = i;  // Capture by value for lambda
            basis[i] = [power](Real x) -> Real {
                if (power == 0) return static_cast<Real>(1.0);
                return std::pow(x, power);
            };
        }
        
        return GeneralLinearLeastSquares(x_data, y_data, basis);
    }
    
    // Evaluate a fitted polynomial at a point
    template<typename Real>
    Real EvaluatePolynomial(Real x, const Vector<Real>& coefficients)
    {
        // Use Horner's method for numerical stability: c_0 + x*(c_1 + x*(c_2 + ...))
        int n = coefficients.size();
        if (n == 0) return 0.0;
        
        Real result = coefficients[n - 1];
        for (int i = n - 2; i >= 0; i--) {
            result = coefficients[i] + x * result;
        }
        return result;
    }
    
} // namespace MML

#endif // MML_CURVE_FITTING_H
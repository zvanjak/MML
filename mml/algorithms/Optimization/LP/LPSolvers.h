///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Optimization/LP/LPSolvers.h                                                  ///
///  Description: Linear programming convenience solve functions                           ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_LP_SOLVERS_H
#define MML_LP_SOLVERS_H

#include <mml/algorithms/Optimization/LP/SimplexSolver.h>

#include <vector>

namespace MML::Optimization {

// Solver roles (MinimalMathLibrary-ya0v.5, Release 2.0):
//  - SimplexSolver / SimplexTableau: MML's LP engine (full tableau); also home of the
//    dual simplex (SolveLPDual) and sensitivity analysis (SolveLPWithSensitivity)
//  - RevisedSimplexSolver (production solver with explicit basis inverse) moved to
//    MML-Packages: include/optimization/RevisedSimplexSolver.h

/// @brief Solve LP using the full-tableau simplex solver
inline LPResult SolveLP(const LinearProgram& lp, const LPConfig& config = LPConfig()) {
    SimplexSolver solver(config);
    return solver.Solve(lp);
}

/// @brief Solve LP from matrix form: min c'x s.t. Ax <= b, x >= 0
inline LPResult SolveLP(const Vector<Real>& c, const Matrix<Real>& A, const Vector<Real>& b,
                        const LPConfig& config = LPConfig()) {
    std::vector<LPConstraintType> types(b.size(), LPConstraintType::LessEqual);
    LinearProgram lp(c, A, b, types);
    return SolveLP(lp, config);
}

/// @brief Solve LP using the full-tableau reference implementation
inline LPResult SolveLPTableau(const LinearProgram& lp, const LPConfig& config = LPConfig()) {
    SimplexSolver solver(config);
    return solver.Solve(lp);
}

/// @brief Solve LP using dual simplex method
inline LPResult SolveLPDual(const LinearProgram& lp, const LPConfig& config = LPConfig()) {
    SimplexSolver solver(config);
    return solver.SolveDual(lp);
}

/// @brief Solve LP from matrix form using dual simplex: min c'x s.t. Ax <= b, x >= 0
inline LPResult SolveLPDual(const Vector<Real>& c, const Matrix<Real>& A, const Vector<Real>& b,
                            const LPConfig& config = LPConfig()) {
    std::vector<LPConstraintType> types(b.size(), LPConstraintType::LessEqual);
    LinearProgram lp(c, A, b, types);
    return SolveLPDual(lp, config);
}

/// @brief Solve LP with sensitivity analysis
inline LPResult SolveLPWithSensitivity(const LinearProgram& lp, const LPConfig& config = LPConfig()) {
    SimplexSolver solver(config);
    return solver.SolveWithSensitivity(lp);
}

/// @brief Solve LP from matrix form with sensitivity analysis: min c'x s.t. Ax <= b, x >= 0
inline LPResult SolveLPWithSensitivity(const Vector<Real>& c, const Matrix<Real>& A, const Vector<Real>& b,
                                       const LPConfig& config = LPConfig()) {
    std::vector<LPConstraintType> types(b.size(), LPConstraintType::LessEqual);
    LinearProgram lp(c, A, b, types);
    return SolveLPWithSensitivity(lp, config);
}

} // namespace MML::Optimization
#endif // MML_LP_SOLVERS_H

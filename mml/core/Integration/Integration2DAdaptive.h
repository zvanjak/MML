///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Integration2DAdaptive.h                                             ///
///  Description: Adaptive 2D numerical integration with quadtree subdivision         ///
///               and error-driven cell refinement                                    ///
///                                                                                   ///
///  Features:    - Quadtree-based 2D integration with local error estimation         ///
///               - Automatic cell subdivision where error exceeds tolerance          ///
///               - Gauss-Kronrod embedded rules for reliable error bounds            ///
///                                                                                   ///
///  Usage:       auto r = IntegrateAdaptive2D(f, x1, x2, y1, y2, 1e-8);              ///
///               if (r.converged) std::cout << r.value << " ± " << r.error_estimate; ///
///                                                                                   ///
///  Algorithm:   1. Apply GK rule to entire domain → estimate + error                ///
///               2. If error < tolerance×(cell_area/total_area) → accept             ///
///               3. Else subdivide into 4 children (quadtree), recurse               ///
///               4. Accumulate results from all converged cells                      ///
///                                                                                   ///
///  References:  [PSWDK83] Piessens et al., QUADPACK, Springer 1983                  ///
///               [NR3] Press et al., Numerical Recipes 3rd ed., Ch. 4                ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                        ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#ifndef MML_INTEGRATION_2D_ADAPTIVE_H
#define MML_INTEGRATION_2D_ADAPTIVE_H

#include <functional>
#include <algorithm>
#include <cmath>

#include "MMLBase.h"
#include "core/Integration/GaussKronrod.h"
#include "core/Integration/IntegrationBase.h"

namespace MML
{
namespace Integration
{

///////////////////////////////////////////////////////////////////////////
///                    RESULT STRUCTURE                                 ///
///////////////////////////////////////////////////////////////////////////

/// @brief Result of adaptive 2D integration
/// 
/// Contains the integral estimate, error bound, and diagnostic information
/// about the quadtree subdivision process.
///
/// @par Example:
/// @code
/// auto result = IntegrateAdaptive2D(f, 0, 1, 0, 1, 1e-8);
/// if (result.converged) {
///     std::cout << "Integral = " << result.value 
///               << " ± " << result.error_estimate << "\n";
///     std::cout << "Used " << result.function_evaluations << " evaluations\n";
/// }
/// @endcode
struct AdaptiveResult2D
{
    Real value;               ///< Integral estimate (sum over all cells)
    Real error_estimate;      ///< Accumulated error bound from all cells
    int function_evaluations; ///< Total number of f(x,y) evaluations
    int cells_subdivided;     ///< Number of quadtree cell subdivisions performed
    int max_depth_reached;    ///< Maximum subdivision depth reached
    bool converged;           ///< True if total error < requested tolerance

    /// @brief Default constructor
    AdaptiveResult2D(Real val = 0.0, Real err = 0.0, int evals = 0, 
                     int cells = 0, int depth = 0, bool conv = false)
        : value(val)
        , error_estimate(err)
        , function_evaluations(evals)
        , cells_subdivided(cells)
        , max_depth_reached(depth)
        , converged(conv) 
    {}

    /// @brief Combine two results (for accumulating from child cells)
    AdaptiveResult2D& operator+=(const AdaptiveResult2D& other) {
        value += other.value;
        error_estimate += other.error_estimate;
        function_evaluations += other.function_evaluations;
        cells_subdivided += other.cells_subdivided;
        max_depth_reached = std::max(max_depth_reached, other.max_depth_reached);
        converged = converged && other.converged;
        return *this;
    }

    /// @brief Add two results
    friend AdaptiveResult2D operator+(AdaptiveResult2D lhs, const AdaptiveResult2D& rhs) {
        lhs += rhs;
        return lhs;
    }
};

///////////////////////////////////////////////////////////////////////////
///                    CONFIGURATION STRUCTURE                          ///
///////////////////////////////////////////////////////////////////////////

/// @brief Configuration for adaptive 2D integration
struct AdaptiveConfig2D
{
    Real tol_abs = 1e-10;        ///< Absolute error tolerance
    Real tol_rel = 1e-8;         ///< Relative error tolerance
    int max_depth = 20;          ///< Maximum quadtree subdivision depth
    int max_evaluations = 1000000; ///< Maximum function evaluations budget
    GKRule rule = GKRule::GK15;  ///< Gauss-Kronrod rule to use

    AdaptiveConfig2D& absoluteTolerance(Real tol) { tol_abs = tol; return *this; }
    AdaptiveConfig2D& relativeTolerance(Real tol) { tol_rel = tol; return *this; }
    AdaptiveConfig2D& maxDepth(int d) { max_depth = d; return *this; }
    AdaptiveConfig2D& maxEvaluations(int n) { max_evaluations = n; return *this; }
    AdaptiveConfig2D& gkRule(GKRule r) { rule = r; return *this; }
};

///////////////////////////////////////////////////////////////////////////
///                    2D ADAPTIVE INTEGRATION                          ///
///////////////////////////////////////////////////////////////////////////

namespace Detail {

/// @brief Tensor product 2D integration using GK15 with proper embedded error estimation
/// @details Evaluates f at all 15×15 = 225 GK15 node pairs, computing both:
///          - Kronrod estimate (15×15): high-order approximation
///          - Gauss estimate (7×7): low-order approximation (using embedded G7 nodes)
///          Error = |Kronrod - Gauss| is a reliable error bound.
///
/// @note This is more expensive than nested 1D but gives proper 2D error estimation.
///       The embedded G7 nodes are at GK15 indices {0, 2, 4, 6} (symmetric around 0).
template<typename Func>
static GKResult IntegrateCell2D_TensorGK15(Func f, Real x1, Real x2, Real y1, Real y2) {
    const Real cx = 0.5 * (x1 + x2);  // Center x
    const Real hx = 0.5 * (x2 - x1);  // Half-width x
    const Real cy = 0.5 * (y1 + y2);  // Center y
    const Real hy = 0.5 * (y2 - y1);  // Half-width y

    Real result_kronrod = 0.0;
    Real result_gauss = 0.0;
    int evals = 0;

    // Gauss node indices within GK15: the G7 nodes are at positions 0, 2, 4, 6
    // (GK15 node 0 = center, nodes 1-7 are positive, symmetric with negatives)
    // G7 uses nodes at positions (relative to GK15): 0 (center), ±2, ±4, ±6
    constexpr int gauss_indices[] = {0, 2, 4, 6};  // 4 indices (representing 7 points with symmetry)

    // Helper to check if a GK15 index maps to a Gauss node
    auto is_gauss_index = [](int i) -> bool {
        return (i == 0 || i == 2 || i == 4 || i == 6);
    };

    // Map GK15 index to G7 weight index
    auto gauss_weight_index = [](int i) -> int {
        return (i == 0) ? 0 : (i == 2) ? 1 : (i == 4) ? 2 : 3;
    };

    // Full tensor product: evaluate at all node pairs
    // GK15 has nodes at xGK15[0], ±xGK15[1], ±xGK15[2], ..., ±xGK15[7]
    // wGK15 gives weights for each node (symmetric)
    
    // Loop over all GK15 x-nodes (8 unique positions: center + 7 positive)
    for (int i = 0; i < GK15_N; ++i) {
        Real xi_pos = cx + hx * xGK15[i];
        Real xi_neg = cx - hx * xGK15[i];
        Real wi = wGK15[i];
        bool i_gauss = is_gauss_index(i);
        Real wgi = i_gauss ? wG7[gauss_weight_index(i)] : 0.0;
        
        // Loop over all GK15 y-nodes
        for (int j = 0; j < GK15_N; ++j) {
            Real yj_pos = cy + hy * xGK15[j];
            Real yj_neg = cy - hy * xGK15[j];
            Real wj = wGK15[j];
            bool j_gauss = is_gauss_index(j);
            Real wgj = j_gauss ? wG7[gauss_weight_index(j)] : 0.0;
            
            // Evaluate function at all quadrant combinations
            // For i=0 (center), xi_pos = xi_neg; same for j=0
            Real f_pp = f(xi_pos, yj_pos);
            evals++;
            
            Real sum_f = f_pp;
            int multiplicity = 1;
            
            if (i > 0 && j > 0) {
                // All 4 quadrants are distinct
                Real f_pn = f(xi_pos, yj_neg);
                Real f_np = f(xi_neg, yj_pos);
                Real f_nn = f(xi_neg, yj_neg);
                evals += 3;
                sum_f = f_pp + f_pn + f_np + f_nn;
                multiplicity = 4;
            } else if (i > 0) {  // j == 0: yj_pos == yj_neg
                Real f_np = f(xi_neg, yj_pos);
                evals++;
                sum_f = f_pp + f_np;
                multiplicity = 2;
            } else if (j > 0) {  // i == 0: xi_pos == xi_neg
                Real f_pn = f(xi_pos, yj_neg);
                evals++;
                sum_f = f_pp + f_pn;
                multiplicity = 2;
            }
            // else i == 0 && j == 0: center point, multiplicity = 1
            
            // Kronrod contribution (using GK15 weights for both dimensions)
            // Note: We sum f values and multiply by weight product once
            // This handles the symmetric pairs correctly
            result_kronrod += wi * wj * sum_f;
            
            // Gauss contribution (only if both indices are Gauss nodes)
            if (i_gauss && j_gauss) {
                result_gauss += wgi * wgj * sum_f;
            }
        }
    }

    // Scale by Jacobian (half-widths product)
    Real jacobian = hx * hy;
    result_kronrod *= jacobian;
    result_gauss *= jacobian;

    Real error = std::abs(result_kronrod - result_gauss);

    return GKResult(result_kronrod, error, evals, true);
}

} // namespace Detail

/// @brief Recursive helper for adaptive 2D integration
/// @param f Function f(x,y) -> Real
/// @param x1,x2,y1,y2 Cell bounds
/// @param total_area Total integration domain area (for tolerance scaling)
/// @param tol_abs Absolute tolerance
/// @param tol_rel Relative tolerance
/// @param depth Current recursion depth
/// @param max_depth Maximum allowed depth
/// @param evals_remaining Budget for function evaluations
/// @param rule Which GK rule to use
template<typename Func>
static AdaptiveResult2D IntegrateAdaptive2D_Recursive(
    Func f, 
    Real x1, Real x2, Real y1, Real y2,
    Real total_area,
    Real tol_abs, Real tol_rel,
    int depth, int max_depth,
    int& evals_remaining,
    GKRule rule)
{
    // Integrate this cell using tensor product GK15 with proper error estimation
    auto cell_result = Detail::IntegrateCell2D_TensorGK15(f, x1, x2, y1, y2);
    evals_remaining -= cell_result.function_evals;

    Real cell_area = (x2 - x1) * (y2 - y1);
    Real area_fraction = cell_area / total_area;
    
    // Tolerance for this cell: scale by area fraction, with epsilon floor to prevent underflow
    Real cell_tol = std::max({tol_abs * area_fraction, 
                              tol_rel * std::abs(cell_result.value),
                              std::numeric_limits<Real>::min()});

    // Check convergence
    if (cell_result.error_estimate <= cell_tol || depth >= max_depth || evals_remaining <= 0) {
        bool conv = (cell_result.error_estimate <= cell_tol);
        return AdaptiveResult2D(cell_result.value, cell_result.error_estimate,
                                 cell_result.function_evals, 0, depth, conv);
    }

    // Subdivide into 4 quadrants
    Real mx = 0.5 * (x1 + x2);
    Real my = 0.5 * (y1 + y2);

    AdaptiveResult2D result;
    result.cells_subdivided = 1;  // This cell was subdivided

    // Process all 4 children
    result += IntegrateAdaptive2D_Recursive(f, x1, mx, y1, my, total_area, tol_abs, tol_rel, depth+1, max_depth, evals_remaining, rule);
    result += IntegrateAdaptive2D_Recursive(f, mx, x2, y1, my, total_area, tol_abs, tol_rel, depth+1, max_depth, evals_remaining, rule);
    result += IntegrateAdaptive2D_Recursive(f, x1, mx, my, y2, total_area, tol_abs, tol_rel, depth+1, max_depth, evals_remaining, rule);
    result += IntegrateAdaptive2D_Recursive(f, mx, x2, my, y2, total_area, tol_abs, tol_rel, depth+1, max_depth, evals_remaining, rule);

    return result;
}

/// @brief Adaptive 2D integration over rectangular domain using quadtree subdivision
/// 
/// Computes ∬_R f(x,y) dx dy over rectangle R = [x1,x2] × [y1,y2] using
/// adaptive quadtree refinement. Cells with error exceeding tolerance are
/// subdivided until convergence or resource limits are reached.
///
/// @tparam Func Callable with signature Real(Real x, Real y)
/// @param f Function to integrate
/// @param x1,x2 X-integration bounds
/// @param y1,y2 Y-integration bounds  
/// @param tolerance Requested absolute accuracy
/// @param max_depth Maximum quadtree subdivision depth (default: 20)
/// @param max_evals Maximum function evaluations budget (default: 1000000)
/// @return AdaptiveResult2D with integral value, error estimate, and diagnostics
///
/// @par Algorithm:
/// 1. Apply tensor-product Gauss-Kronrod rule to the cell
/// 2. If error ≤ tolerance × (cell_area / total_area) → accept
/// 3. Otherwise subdivide into 4 quadrants and recurse
/// 4. Accumulate results from all leaf cells
///
/// @par Example:
/// @code
/// // Integrate x*y over unit square
/// auto result = IntegrateAdaptive2D(
///     [](Real x, Real y) { return x * y; },
///     0.0, 1.0, 0.0, 1.0, 1e-10
/// );
/// std::cout << result.value << " (should be 0.25)\n";
/// 
/// // Peak function - adaptive shines here
/// auto peak = IntegrateAdaptive2D(
///     [](Real x, Real y) { return std::exp(-100*((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5))); },
///     0.0, 1.0, 0.0, 1.0, 1e-8
/// );
/// std::cout << "Peak integral: " << peak.value << ", evals: " << peak.function_evaluations << "\n";
/// @endcode
template<typename Func>
static AdaptiveResult2D IntegrateAdaptive2D(
    Func f,
    Real x1, Real x2,
    Real y1, Real y2,
    Real tolerance = 1e-8,
    int max_depth = 20,
    int max_evals = 1000000)
{
    if (x2 <= x1 || y2 <= y1) {
        return AdaptiveResult2D(0.0, 0.0, 0, 0, 0, true);
    }

    Real total_area = (x2 - x1) * (y2 - y1);
    int evals_remaining = max_evals;

    return IntegrateAdaptive2D_Recursive(
        f, x1, x2, y1, y2,
        total_area,
        tolerance, tolerance,  // Use same tolerance for both abs and rel
        0, max_depth,
        evals_remaining,
        GKRule::GK15
    );
}

/// @brief Adaptive 2D integration with configuration object
/// 
/// Provides full control over integration parameters including
/// separate absolute/relative tolerances and choice of GK rule.
///
/// @param f Function to integrate
/// @param x1,x2,y1,y2 Integration bounds
/// @param config Configuration object with tolerances, limits, and rule selection
/// @return AdaptiveResult2D with integral value, error estimate, and diagnostics
template<typename Func>
static AdaptiveResult2D IntegrateAdaptive2D(
    Func f,
    Real x1, Real x2,
    Real y1, Real y2,
    const AdaptiveConfig2D& config)
{
    if (x2 <= x1 || y2 <= y1) {
        return AdaptiveResult2D(0.0, 0.0, 0, 0, 0, true);
    }

    Real total_area = (x2 - x1) * (y2 - y1);
    int evals_remaining = config.max_evaluations;

    return IntegrateAdaptive2D_Recursive(
        f, x1, x2, y1, y2,
        total_area,
        config.tol_abs, config.tol_rel,
        0, config.max_depth,
        evals_remaining,
        config.rule
    );
}

/**************************************************************************/
/*****       Detailed API - Adaptive 2D Integration                   *****/
/**************************************************************************/

/// Adaptive 2D integration with simple parameters and full diagnostics
template<typename Func>
static IntegrationDetailedResult IntegrateAdaptive2DDetailed(
    Func f,
    Real x1, Real x2,
    Real y1, Real y2,
    const IntegrationConfig& config = {},
    Real tolerance = 1e-8,
    int max_depth = 20,
    int max_evals = 1000000)
{
    return IntegrationDetail::ExecuteIntegrationDetailed<IntegrationDetailedResult>(
        "IntegrateAdaptive2D", config,
        [&](IntegrationDetailedResult& result) {
            auto r = IntegrateAdaptive2D(f, x1, x2, y1, y2, tolerance, max_depth, max_evals);
            result.value = r.value;
            result.error_estimate = r.error_estimate;
            result.function_evaluations = r.function_evaluations;
            result.converged = r.converged;
            if (!r.converged) {
                result.status = AlgorithmStatus::MaxIterationsExceeded;
                result.error_message = "IntegrateAdaptive2D did not converge";
            }
        });
}

/// Adaptive 2D integration with config object and full diagnostics
template<typename Func>
static IntegrationDetailedResult IntegrateAdaptive2DDetailed(
    Func f,
    Real x1, Real x2,
    Real y1, Real y2,
    const AdaptiveConfig2D& adaptive_config,
    const IntegrationConfig& config = {})
{
    return IntegrationDetail::ExecuteIntegrationDetailed<IntegrationDetailedResult>(
        "IntegrateAdaptive2D", config,
        [&](IntegrationDetailedResult& result) {
            auto r = IntegrateAdaptive2D(f, x1, x2, y1, y2, adaptive_config);
            result.value = r.value;
            result.error_estimate = r.error_estimate;
            result.function_evaluations = r.function_evaluations;
            result.converged = r.converged;
            if (!r.converged) {
                result.status = AlgorithmStatus::MaxIterationsExceeded;
                result.error_message = "IntegrateAdaptive2D did not converge";
            }
        });
}

} // namespace Integration
} // namespace MML

#endif // MML_INTEGRATION_2D_ADAPTIVE_H

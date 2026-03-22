///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Integration3DAdaptive.h                                             ///
///  Description: Adaptive 3D numerical integration with octree subdivision           ///
///               and error-driven cell refinement                                    ///
///                                                                                   ///
///  Features:    - Octree-based 3D integration with local error estimation           ///
///               - Automatic cell subdivision where error exceeds tolerance          ///
///               - Gauss-Kronrod embedded rules for reliable error bounds            ///
///                                                                                   ///
///  Usage:       auto r = IntegrateAdaptive3D(f, x1, x2, y1, y2, z1, z2, 1e-8);      ///
///               if (r.converged) std::cout << r.value << " ± " << r.error_estimate; ///
///                                                                                   ///
///  Algorithm:   1. Apply GK rule to entire domain → estimate + error                ///
///               2. If error < tolerance×(cell_volume/total_volume) → accept         ///
///               3. Else subdivide into 8 children (octree), recurse                 ///
///               4. Accumulate results from all converged cells                      ///
///                                                                                   ///
///  References:  [PSWDK83] Piessens et al., QUADPACK, Springer 1983                  ///
///               [NR3] Press et al., Numerical Recipes 3rd ed., Ch. 4                ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                        ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#ifndef MML_INTEGRATION_3D_ADAPTIVE_H
#define MML_INTEGRATION_3D_ADAPTIVE_H

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

/// @brief Result of adaptive 3D integration
/// 
/// Contains the integral estimate, error bound, and diagnostic information
/// about the octree subdivision process.
///
/// @par Example:
/// @code
/// auto result = IntegrateAdaptive3D(f, 0, 1, 0, 1, 0, 1, 1e-8);
/// if (result.converged) {
///     std::cout << "Integral = " << result.value 
///               << " ± " << result.error_estimate << "\n";
///     std::cout << "Subdivided " << result.cells_subdivided << " cells\n";
/// }
/// @endcode
struct AdaptiveResult3D
{
    Real value;               ///< Integral estimate (sum over all cells)
    Real error_estimate;      ///< Accumulated error bound from all cells
    int function_evaluations; ///< Total number of f(x,y,z) evaluations
    int cells_subdivided;     ///< Number of octree cell subdivisions performed
    int max_depth_reached;    ///< Maximum subdivision depth reached
    bool converged;           ///< True if total error < requested tolerance

    /// @brief Default constructor
    AdaptiveResult3D(Real val = 0.0, Real err = 0.0, int evals = 0, 
                     int cells = 0, int depth = 0, bool conv = false)
        : value(val)
        , error_estimate(err)
        , function_evaluations(evals)
        , cells_subdivided(cells)
        , max_depth_reached(depth)
        , converged(conv) 
    {}

    /// @brief Combine two results (for accumulating from child cells)
    AdaptiveResult3D& operator+=(const AdaptiveResult3D& other) {
        value += other.value;
        error_estimate += other.error_estimate;
        function_evaluations += other.function_evaluations;
        cells_subdivided += other.cells_subdivided;
        max_depth_reached = std::max(max_depth_reached, other.max_depth_reached);
        converged = converged && other.converged;
        return *this;
    }

    /// @brief Add two results
    friend AdaptiveResult3D operator+(AdaptiveResult3D lhs, const AdaptiveResult3D& rhs) {
        lhs += rhs;
        return lhs;
    }
};

///////////////////////////////////////////////////////////////////////////
///                    CONFIGURATION STRUCTURE                          ///
///////////////////////////////////////////////////////////////////////////

/// @brief Configuration for adaptive 3D integration
struct AdaptiveConfig3D
{
    Real tol_abs = 1e-10;        ///< Absolute error tolerance
    Real tol_rel = 1e-8;         ///< Relative error tolerance
    int max_depth = 15;          ///< Maximum octree subdivision depth (smaller than 2D!)
    int max_evaluations = 10000000; ///< Maximum function evaluations budget
    GKRule rule = GKRule::GK15;  ///< Gauss-Kronrod rule to use

    AdaptiveConfig3D& absoluteTolerance(Real tol) { tol_abs = tol; return *this; }
    AdaptiveConfig3D& relativeTolerance(Real tol) { tol_rel = tol; return *this; }
    AdaptiveConfig3D& maxDepth(int d) { max_depth = d; return *this; }
    AdaptiveConfig3D& maxEvaluations(int n) { max_evaluations = n; return *this; }
    AdaptiveConfig3D& gkRule(GKRule r) { rule = r; return *this; }
};

///////////////////////////////////////////////////////////////////////////
///                    3D ADAPTIVE INTEGRATION                          ///
///////////////////////////////////////////////////////////////////////////

namespace Detail {

/// @brief Tensor product 3D integration using GK15 with proper embedded error estimation
/// @details Evaluates f at all 15×15×15 = 3375 GK15 node triplets, computing both:
///          - Kronrod estimate (15×15×15): high-order approximation
///          - Gauss estimate (7×7×7): low-order approximation (using embedded G7 nodes)
///          Error = |Kronrod - Gauss| is a reliable error bound.
///
/// @note This is expensive (3375 evals vs 225 for 2D) but gives proper 3D error estimation.
template<typename Func>
static GKResult IntegrateCell3D_TensorGK15(Func f, Real x1, Real x2, Real y1, Real y2, Real z1, Real z2) {
    const Real cx = 0.5 * (x1 + x2);  // Center x
    const Real hx = 0.5 * (x2 - x1);  // Half-width x
    const Real cy = 0.5 * (y1 + y2);  // Center y
    const Real hy = 0.5 * (y2 - y1);  // Half-width y
    const Real cz = 0.5 * (z1 + z2);  // Center z
    const Real hz = 0.5 * (z2 - z1);  // Half-width z

    Real result_kronrod = 0.0;
    Real result_gauss = 0.0;
    int evals = 0;

    // Helper to check if a GK15 index maps to a Gauss node
    auto is_gauss_index = [](int i) -> bool {
        return (i == 0 || i == 2 || i == 4 || i == 6);
    };

    // Map GK15 index to G7 weight index
    auto gauss_weight_index = [](int i) -> int {
        return (i == 0) ? 0 : (i == 2) ? 1 : (i == 4) ? 2 : 3;
    };

    // Full tensor product: loop over all node triplets
    for (int i = 0; i < GK15_N; ++i) {
        Real xi_pos = cx + hx * xGK15[i];
        Real xi_neg = cx - hx * xGK15[i];
        Real wi = wGK15[i];
        bool i_gauss = is_gauss_index(i);
        Real wgi = i_gauss ? wG7[gauss_weight_index(i)] : 0.0;
        
        for (int j = 0; j < GK15_N; ++j) {
            Real yj_pos = cy + hy * xGK15[j];
            Real yj_neg = cy - hy * xGK15[j];
            Real wj = wGK15[j];
            bool j_gauss = is_gauss_index(j);
            Real wgj = j_gauss ? wG7[gauss_weight_index(j)] : 0.0;
            
            for (int k = 0; k < GK15_N; ++k) {
                Real zk_pos = cz + hz * xGK15[k];
                Real zk_neg = cz - hz * xGK15[k];
                Real wk = wGK15[k];
                bool k_gauss = is_gauss_index(k);
                Real wgk = k_gauss ? wG7[gauss_weight_index(k)] : 0.0;
                
                // Evaluate function at all octant combinations
                Real f_ppp = f(xi_pos, yj_pos, zk_pos);
                evals++;
                
                Real sum_f = f_ppp;
                
                // Handle symmetry: 2^n octants where n = number of non-zero indices
                if (i > 0 && j > 0 && k > 0) {
                    // All 8 octants are distinct
                    sum_f += f(xi_pos, yj_pos, zk_neg);
                    sum_f += f(xi_pos, yj_neg, zk_pos);
                    sum_f += f(xi_pos, yj_neg, zk_neg);
                    sum_f += f(xi_neg, yj_pos, zk_pos);
                    sum_f += f(xi_neg, yj_pos, zk_neg);
                    sum_f += f(xi_neg, yj_neg, zk_pos);
                    sum_f += f(xi_neg, yj_neg, zk_neg);
                    evals += 7;
                } else if (i > 0 && j > 0) {  // k == 0
                    sum_f += f(xi_pos, yj_neg, zk_pos);
                    sum_f += f(xi_neg, yj_pos, zk_pos);
                    sum_f += f(xi_neg, yj_neg, zk_pos);
                    evals += 3;
                } else if (i > 0 && k > 0) {  // j == 0
                    sum_f += f(xi_pos, yj_pos, zk_neg);
                    sum_f += f(xi_neg, yj_pos, zk_pos);
                    sum_f += f(xi_neg, yj_pos, zk_neg);
                    evals += 3;
                } else if (j > 0 && k > 0) {  // i == 0
                    sum_f += f(xi_pos, yj_pos, zk_neg);
                    sum_f += f(xi_pos, yj_neg, zk_pos);
                    sum_f += f(xi_pos, yj_neg, zk_neg);
                    evals += 3;
                } else if (i > 0) {  // j == 0 && k == 0
                    sum_f += f(xi_neg, yj_pos, zk_pos);
                    evals += 1;
                } else if (j > 0) {  // i == 0 && k == 0
                    sum_f += f(xi_pos, yj_neg, zk_pos);
                    evals += 1;
                } else if (k > 0) {  // i == 0 && j == 0
                    sum_f += f(xi_pos, yj_pos, zk_neg);
                    evals += 1;
                }
                // else all center: multiplicity = 1
                
                // Kronrod contribution
                result_kronrod += wi * wj * wk * sum_f;
                
                // Gauss contribution (only if all indices are Gauss nodes)
                if (i_gauss && j_gauss && k_gauss) {
                    result_gauss += wgi * wgj * wgk * sum_f;
                }
            }
        }
    }

    // Scale by Jacobian (product of half-widths)
    Real jacobian = hx * hy * hz;
    result_kronrod *= jacobian;
    result_gauss *= jacobian;

    Real error = std::abs(result_kronrod - result_gauss);

    return GKResult(result_kronrod, error, evals, true);
}

} // namespace Detail

/// @brief Recursive helper for adaptive 3D integration
/// @param f Function f(x,y,z) -> Real
/// @param x1,x2,y1,y2,z1,z2 Cell bounds
/// @param total_volume Total integration domain volume (for tolerance scaling)
/// @param tol_abs Absolute tolerance
/// @param tol_rel Relative tolerance
/// @param depth Current recursion depth
/// @param max_depth Maximum allowed depth
/// @param evals_remaining Budget for function evaluations
/// @param rule Which GK rule to use
template<typename Func>
static AdaptiveResult3D IntegrateAdaptive3D_Recursive(
    Func f, 
    Real x1, Real x2, Real y1, Real y2, Real z1, Real z2,
    Real total_volume,
    Real tol_abs, Real tol_rel,
    int depth, int max_depth,
    int& evals_remaining,
    GKRule rule)
{
    // Integrate this cell using tensor product GK15 with proper error estimation
    auto cell_result = Detail::IntegrateCell3D_TensorGK15(f, x1, x2, y1, y2, z1, z2);
    evals_remaining -= cell_result.function_evals;

    Real cell_volume = (x2 - x1) * (y2 - y1) * (z2 - z1);
    Real volume_fraction = cell_volume / total_volume;
    
    // Tolerance for this cell: scale by volume fraction
    Real cell_tol = std::max(tol_abs * volume_fraction, tol_rel * std::abs(cell_result.value));

    // Check convergence
    if (cell_result.error_estimate <= cell_tol || depth >= max_depth || evals_remaining <= 0) {
        bool conv = (cell_result.error_estimate <= cell_tol);
        return AdaptiveResult3D(cell_result.value, cell_result.error_estimate,
                                 cell_result.function_evals, 0, depth, conv);
    }

    // Subdivide into 8 octants
    Real mx = 0.5 * (x1 + x2);
    Real my = 0.5 * (y1 + y2);
    Real mz = 0.5 * (z1 + z2);

    AdaptiveResult3D result;
    result.cells_subdivided = 1;  // This cell was subdivided

    // Process all 8 children (2×2×2 = 8 octants)
    // Lower z half (z1 to mz)
    result += IntegrateAdaptive3D_Recursive(f, x1, mx, y1, my, z1, mz, total_volume, tol_abs, tol_rel, depth+1, max_depth, evals_remaining, rule);
    result += IntegrateAdaptive3D_Recursive(f, mx, x2, y1, my, z1, mz, total_volume, tol_abs, tol_rel, depth+1, max_depth, evals_remaining, rule);
    result += IntegrateAdaptive3D_Recursive(f, x1, mx, my, y2, z1, mz, total_volume, tol_abs, tol_rel, depth+1, max_depth, evals_remaining, rule);
    result += IntegrateAdaptive3D_Recursive(f, mx, x2, my, y2, z1, mz, total_volume, tol_abs, tol_rel, depth+1, max_depth, evals_remaining, rule);
    // Upper z half (mz to z2)
    result += IntegrateAdaptive3D_Recursive(f, x1, mx, y1, my, mz, z2, total_volume, tol_abs, tol_rel, depth+1, max_depth, evals_remaining, rule);
    result += IntegrateAdaptive3D_Recursive(f, mx, x2, y1, my, mz, z2, total_volume, tol_abs, tol_rel, depth+1, max_depth, evals_remaining, rule);
    result += IntegrateAdaptive3D_Recursive(f, x1, mx, my, y2, mz, z2, total_volume, tol_abs, tol_rel, depth+1, max_depth, evals_remaining, rule);
    result += IntegrateAdaptive3D_Recursive(f, mx, x2, my, y2, mz, z2, total_volume, tol_abs, tol_rel, depth+1, max_depth, evals_remaining, rule);

    return result;
}

/// @brief Adaptive 3D integration over rectangular box using octree subdivision
/// 
/// Computes ∭_B f(x,y,z) dx dy dz over box B = [x1,x2] × [y1,y2] × [z1,z2] using
/// adaptive octree refinement. Cells with error exceeding tolerance are
/// subdivided until convergence or resource limits are reached.
///
/// @tparam Func Callable with signature Real(Real x, Real y, Real z)
/// @param f Function to integrate
/// @param x1,x2 X-integration bounds
/// @param y1,y2 Y-integration bounds  
/// @param z1,z2 Z-integration bounds
/// @param tolerance Requested absolute accuracy
/// @param max_depth Maximum octree subdivision depth (default: 15, smaller than 2D!)
/// @param max_evals Maximum function evaluations budget (default: 10000000)
/// @return AdaptiveResult3D with integral value, error estimate, and diagnostics
///
/// @par Algorithm:
/// 1. Apply tensor-product Gauss-Kronrod integration to the cell
/// 2. If error ≤ tolerance × (cell_volume / total_volume) → accept
/// 3. Otherwise subdivide into 8 octants and recurse
/// 4. Accumulate results from all leaf cells
///
/// @par Example:
/// @code
/// // Integrate 1 over unit cube = 1.0
/// auto result = IntegrateAdaptive3D(
///     [](Real x, Real y, Real z) { return 1.0; },
///     0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1e-10
/// );
/// std::cout << result.value << " (should be 1.0)\n";
/// 
/// // Sphere indicator function: volume of unit sphere = 4π/3
/// auto sphere = IntegrateAdaptive3D(
///     [](Real x, Real y, Real z) { 
///         return (x*x + y*y + z*z <= 1.0) ? 1.0 : 0.0; 
///     },
///     -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, 1e-4
/// );
/// std::cout << "Sphere volume: " << sphere.value << " (4π/3 ≈ 4.189)\n";
/// @endcode
template<typename Func>
static AdaptiveResult3D IntegrateAdaptive3D(
    Func f,
    Real x1, Real x2,
    Real y1, Real y2,
    Real z1, Real z2,
    Real tolerance = 1e-8,
    int max_depth = 15,
    int max_evals = 10000000)
{
    if (x2 <= x1 || y2 <= y1 || z2 <= z1) {
        return AdaptiveResult3D(0.0, 0.0, 0, 0, 0, true);
    }

    Real total_volume = (x2 - x1) * (y2 - y1) * (z2 - z1);
    int evals_remaining = max_evals;

    return IntegrateAdaptive3D_Recursive(
        f, x1, x2, y1, y2, z1, z2,
        total_volume,
        tolerance, tolerance,  // Use same tolerance for both abs and rel
        0, max_depth,
        evals_remaining,
        GKRule::GK15
    );
}

/// @brief Adaptive 3D integration with configuration object
/// 
/// Provides full control over integration parameters including
/// separate absolute/relative tolerances and choice of GK rule.
///
/// @param f Function to integrate
/// @param x1,x2,y1,y2,z1,z2 Integration bounds
/// @param config Configuration object with tolerances, limits, and rule selection
/// @return AdaptiveResult3D with integral value, error estimate, and diagnostics
///
/// @par Example:
/// @code
/// AdaptiveConfig3D config;
/// config.absoluteTolerance(1e-6)
///       .relativeTolerance(1e-4)
///       .maxDepth(10)
///       .maxEvaluations(5000000);
///
/// auto result = IntegrateAdaptive3D(f, 0, 1, 0, 1, 0, 1, config);
/// @endcode
template<typename Func>
static AdaptiveResult3D IntegrateAdaptive3D(
    Func f,
    Real x1, Real x2,
    Real y1, Real y2,
    Real z1, Real z2,
    const AdaptiveConfig3D& config)
{
    if (x2 <= x1 || y2 <= y1 || z2 <= z1) {
        return AdaptiveResult3D(0.0, 0.0, 0, 0, 0, true);
    }

    Real total_volume = (x2 - x1) * (y2 - y1) * (z2 - z1);
    int evals_remaining = config.max_evaluations;

    return IntegrateAdaptive3D_Recursive(
        f, x1, x2, y1, y2, z1, z2,
        total_volume,
        config.tol_abs, config.tol_rel,
        0, config.max_depth,
        evals_remaining,
        config.rule
    );
}

/**************************************************************************/
/*****       Detailed API - Adaptive 3D Integration                   *****/
/**************************************************************************/

/// Adaptive 3D integration with simple parameters and full diagnostics
template<typename Func>
static IntegrationDetailedResult IntegrateAdaptive3DDetailed(
    Func f,
    Real x1, Real x2,
    Real y1, Real y2,
    Real z1, Real z2,
    const IntegrationConfig& config = {},
    Real tolerance = 1e-6,
    int max_depth = 15,
    int max_evals = 10000000)
{
    return IntegrationDetail::ExecuteIntegrationDetailed<IntegrationDetailedResult>(
        "IntegrateAdaptive3D", config,
        [&](IntegrationDetailedResult& result) {
            auto r = IntegrateAdaptive3D(f, x1, x2, y1, y2, z1, z2, tolerance, max_depth, max_evals);
            result.value = r.value;
            result.error_estimate = r.error_estimate;
            result.function_evaluations = r.function_evaluations;
            result.converged = r.converged;
            if (!r.converged) {
                result.status = AlgorithmStatus::MaxIterationsExceeded;
                result.error_message = "IntegrateAdaptive3D did not converge";
            }
        });
}

/// Adaptive 3D integration with config object and full diagnostics
template<typename Func>
static IntegrationDetailedResult IntegrateAdaptive3DDetailed(
    Func f,
    Real x1, Real x2,
    Real y1, Real y2,
    Real z1, Real z2,
    const AdaptiveConfig3D& adaptive_config,
    const IntegrationConfig& config = {})
{
    return IntegrationDetail::ExecuteIntegrationDetailed<IntegrationDetailedResult>(
        "IntegrateAdaptive3D", config,
        [&](IntegrationDetailedResult& result) {
            auto r = IntegrateAdaptive3D(f, x1, x2, y1, y2, z1, z2, adaptive_config);
            result.value = r.value;
            result.error_estimate = r.error_estimate;
            result.function_evaluations = r.function_evaluations;
            result.converged = r.converged;
            if (!r.converged) {
                result.status = AlgorithmStatus::MaxIterationsExceeded;
                result.error_message = "IntegrateAdaptive3D did not converge";
            }
        });
}

} // namespace Integration
} // namespace MML

#endif // MML_INTEGRATION_3D_ADAPTIVE_H

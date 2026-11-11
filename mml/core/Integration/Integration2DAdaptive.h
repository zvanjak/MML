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
///  Usage:       auto r = IntegrateAdaptive2DCore(f, x1, x2, y1, y2, 1e-8);              ///
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

#include <mml/MMLBase.h>
#include <mml/core/Integration/GaussKronrod.h>
#include <mml/core/Integration/IntegrationBase.h>

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
/// auto result = IntegrateAdaptive2DCore(f, 0, 1, 0, 1, 1e-8);
/// if (result.converged) {
///     std::cout << "Integral = " << result.value 
///               << " ± " << result.error_estimate << "\n";
///     std::cout << "Used " << result.function_evaluations << " evaluations\n";
/// }
/// @endcode
struct AdaptiveResult2D : public IntegrationResult
{
    // value, error_estimate, function_evaluations, converged, status come from IntegrationResult
    int cells_subdivided = 0;     ///< Number of quadtree cell subdivisions performed
    int max_depth_reached = 0;    ///< Maximum subdivision depth reached

    /// @brief Default constructor
    AdaptiveResult2D(Real val = 0.0, Real err = 0.0, int evals = 0, 
                     int cells = 0, int depth = 0, bool conv = false)
        : IntegrationResult(val, err, depth, conv, /*error_estimated*/ true)
        , cells_subdivided(cells)
        , max_depth_reached(depth)
    {
        function_evaluations = evals;
        algorithm_name = "IntegrateAdaptive2D";
    }

    /// @brief Combine two results (for accumulating from child cells)
    AdaptiveResult2D& operator+=(const AdaptiveResult2D& other) {
        const bool identity = function_evaluations == 0 && cells_subdivided == 0 && value == 0 && error_estimate == 0;
        value += other.value;
        error_estimate += other.error_estimate;
        function_evaluations += other.function_evaluations;
        cells_subdivided += other.cells_subdivided;
        max_depth_reached = std::max(max_depth_reached, other.max_depth_reached);
        iterations = max_depth_reached;
        converged = (identity || converged) && other.converged;
        status = converged ? AlgorithmStatus::Success : AlgorithmStatus::MaxIterationsExceeded;
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
    Real tol_abs = PrecisionValues<Real>::IntegrationTolerance;  ///< Absolute error tolerance
    Real tol_rel = PrecisionValues<Real>::IntegrationTolerance;  ///< Relative error tolerance
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
///          Error = |Kronrod - Gauss| is an embedded estimate, not a rigorous bound.
///
/// @note This is more expensive than nested 1D but gives proper 2D error estimation.
///       The embedded G7 nodes are at GK15 indices {0, 2, 4, 6} (symmetric around 0).
template<typename Func>
static GKResult IntegrateCell2D_TensorGK15(Func f, Real x1, Real x2, Real y1, Real y2) {
    const Real cx = std::midpoint(x1, x2);
    const Real hx = IntegrationDetail::HalfWidth(x1, x2);
    const Real cy = std::midpoint(y1, y2);
    const Real hy = IntegrationDetail::HalfWidth(y1, y2);

    Real result_kronrod = 0.0;
    Real result_gauss = 0.0;
    int evals = 0;

    // Gauss node indices within GK15: the G7 nodes are at positions 0, 2, 4, 6
    // (GK15 node 0 = center, nodes 1-7 are positive, symmetric with negatives)
    // G7 uses nodes at positions (relative to GK15): 0 (center), ±2, ±4, ±6

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
        Real wgi = i_gauss ? wG7[gauss_weight_index(i)] : MML::Real(0.0);
        
        // Loop over all GK15 y-nodes
        for (int j = 0; j < GK15_N; ++j) {
            Real yj_pos = cy + hy * xGK15[j];
            Real yj_neg = cy - hy * xGK15[j];
            Real wj = wGK15[j];
            bool j_gauss = is_gauss_index(j);
            Real wgj = j_gauss ? wG7[gauss_weight_index(j)] : MML::Real(0.0);
            
            // Evaluate function at all quadrant combinations
            // For i=0 (center), xi_pos = xi_neg; same for j=0
            Real f_pp = f(xi_pos, yj_pos);
            evals++;
            
            Real sum_f = f_pp;
            
            if (i > 0 && j > 0) {
                // All 4 quadrants are distinct
                Real f_pn = f(xi_pos, yj_neg);
                Real f_np = f(xi_neg, yj_pos);
                Real f_nn = f(xi_neg, yj_neg);
                evals += 3;
                sum_f = f_pp + f_pn + f_np + f_nn;
            } else if (i > 0) {  // j == 0: yj_pos == yj_neg
                Real f_np = f(xi_neg, yj_pos);
                evals++;
                sum_f = f_pp + f_np;
            } else if (j > 0) {  // i == 0: xi_pos == xi_neg
                Real f_pn = f(xi_pos, yj_neg);
                evals++;
                sum_f = f_pp + f_pn;
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
    result_kronrod = IntegrationDetail::Product({hx, hy, result_kronrod});
    result_gauss = IntegrationDetail::Product({hx, hy, result_gauss});

    Real error = std::abs(result_kronrod - result_gauss);

    return GKResult::Fixed(result_kronrod, error, evals);
}

} // namespace Detail

template<class Func>
static void RunAdaptive2D(Func& function, Real x1, Real x2, Real y1, Real y2,
    Real absolute, Real relative, int depth, int budget, GKRule rule, AdaptiveResult2D& result, bool checkFinite = true) {
    IntegrationDetail::HalfWidth(x1, x2);
    IntegrationDetail::HalfWidth(y1, y2);
    IntegrationDetail::ValidateAdaptive(absolute, relative, depth, budget);
    if (rule != GKRule::GK15)
        throw NotImplementedError("IntegrateAdaptive2D: only GKRule::GK15 is implemented");
    if (x2 <= x1 || y2 <= y1) {
        result = AdaptiveResult2D(0, 0, 0, 0, 0, true);
        return;
    }
    using Box = std::array<std::array<Real, 2>, 2>;
    auto sample = [&](Real x, Real y) {
        return IntegrationDetail::Sample(function, result.function_evaluations, checkFinite, x, y);
    };
    auto evaluate = [&](const Box& box) {
        return Detail::IntegrateCell2D_TensorGK15(sample, box[0][0], box[0][1], box[1][0], box[1][1]);
    };
    const Box domain{{{x1, x2}, {y1, y2}}};
    IntegrationDetail::RefineCells(result, std::vector<Box>{domain}, evaluate,
        IntegrationDetail::SplitBox<2>, absolute, relative, depth, budget, 225);
}

/// @brief Compatibility entry for adaptive 2D integration using global frontier refinement
/// @param f Function f(x,y) -> Real
/// @param x1,x2,y1,y2 Cell bounds
/// @param total_area Legacy parameter; the requested domain now has one global error target
/// @param tol_abs Absolute tolerance
/// @param tol_rel Relative tolerance
/// @param depth Current recursion depth
/// @param max_depth Maximum allowed depth
/// @param evals_remaining Budget for function evaluations
template<typename Func>
static AdaptiveResult2D IntegrateAdaptive2D_Recursive(
    Func f, 
    Real x1, Real x2, Real y1, Real y2,
    Real total_area,
    Real tol_abs, Real tol_rel,
    int depth, int max_depth,
    int& evals_remaining)
{
    AdaptiveResult2D result;
    RunAdaptive2D(f, x1, x2, y1, y2, tol_abs, tol_rel, max_depth - depth, evals_remaining, GKRule::GK15, result);
    evals_remaining -= result.function_evaluations;
    result.max_depth_reached += depth;
    result.iterations = result.max_depth_reached;
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
/// auto result = IntegrateAdaptive2DCore(
///     [](Real x, Real y) { return x * y; },
///     0.0, 1.0, 0.0, 1.0, 1e-10
/// );
/// std::cout << result.value << " (should be 0.25)\n";
/// 
/// // Peak function - adaptive shines here
/// auto peak = IntegrateAdaptive2DCore(
///     [](Real x, Real y) { return std::exp(-100*((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5))); },
///     0.0, 1.0, 0.0, 1.0, 1e-8
/// );
/// std::cout << "Peak integral: " << peak.value << ", evals: " << peak.function_evaluations << "\n";
/// @endcode
template<typename Func>
static AdaptiveResult2D IntegrateAdaptive2DCore(
    Func f,
    Real x1, Real x2,
    Real y1, Real y2,
    Real tolerance = PrecisionValues<Real>::IntegrationTolerance,
    int max_depth = 20,
    int max_evals = 1000000)
{
    AdaptiveResult2D result;
    RunAdaptive2D(f, x1, x2, y1, y2, tolerance, tolerance, max_depth, max_evals, GKRule::GK15, result);
    return result;
}

/// Result-returning variant with the plain argument list (RESULT_POLICY §4).
template<typename Func>
static auto IntegrateAdaptive2DDetailed(Func f, Real x1, Real x2, Real y1, Real y2, Real tolerance, int max_depth = 20, int max_evals = 1000000)
{
	auto r = IntegrateAdaptive2DCore(f, x1, x2, y1, y2, tolerance, max_depth, max_evals);
	r.algorithm_name = "IntegrateAdaptive2D";
	return r;
}

template<typename Func>
static Real IntegrateAdaptive2D(
    Func f,
    Real x1, Real x2,
    Real y1, Real y2,
    Real tolerance = PrecisionValues<Real>::IntegrationTolerance,
    int max_depth = 20,
    int max_evals = 1000000)
{
	return IntegrateAdaptive2DCore(f, x1, x2, y1, y2, tolerance, max_depth, max_evals).ValueOrThrow("IntegrateAdaptive2D");
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
static AdaptiveResult2D IntegrateAdaptive2DCore(
    Func f,
    Real x1, Real x2,
    Real y1, Real y2,
    const AdaptiveConfig2D& config)
{
    AdaptiveResult2D result;
    RunAdaptive2D(f, x1, x2, y1, y2, config.tol_abs, config.tol_rel, config.max_depth, config.max_evaluations, config.rule, result);
    return result;
}

template<typename Func>
static Real IntegrateAdaptive2D(
    Func f,
    Real x1, Real x2,
    Real y1, Real y2,
    const AdaptiveConfig2D& config)
{
	return IntegrateAdaptive2DCore(f, x1, x2, y1, y2, config).ValueOrThrow("IntegrateAdaptive2D");
}

/**************************************************************************/
/*****       Detailed API - Adaptive 2D Integration                   *****/
/**************************************************************************/

/// Adaptive 2D integration with simple parameters and full diagnostics
template<typename Func>
static AdaptiveResult2D IntegrateAdaptive2DDetailed(
    Func f,
    Real x1, Real x2,
    Real y1, Real y2,
    const IntegrationConfig& config = {},
    Real tolerance = PrecisionValues<Real>::IntegrationTolerance,
    int max_depth = 20,
    int max_evals = 1000000)
{
    return IntegrationDetail::ExecuteIntegrationDetailed<AdaptiveResult2D>(
        "IntegrateAdaptive2D", config,
        [&](AdaptiveResult2D& result) {
            RunAdaptive2D(f, x1, x2, y1, y2, tolerance, tolerance, max_depth, max_evals, GKRule::GK15, result, config.check_finite);
        });
}

/// Adaptive 2D integration with config object and full diagnostics
template<typename Func>
static AdaptiveResult2D IntegrateAdaptive2DDetailed(
    Func f,
    Real x1, Real x2,
    Real y1, Real y2,
    const AdaptiveConfig2D& adaptive_config,
    const IntegrationConfig& config = {})
{
    return IntegrationDetail::ExecuteIntegrationDetailed<AdaptiveResult2D>(
        "IntegrateAdaptive2D", config,
        [&](AdaptiveResult2D& result) {
            RunAdaptive2D(f, x1, x2, y1, y2, adaptive_config.tol_abs, adaptive_config.tol_rel, adaptive_config.max_depth,
                adaptive_config.max_evaluations, adaptive_config.rule, result, config.check_finite);
        });
}

} // namespace Integration
} // namespace MML

#endif // MML_INTEGRATION_2D_ADAPTIVE_H

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
///  Usage:       auto r = IntegrateAdaptive3DCore(f, x1, x2, y1, y2, z1, z2, 1e-8);      ///
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

/// @brief Result of adaptive 3D integration
/// 
/// Contains the integral estimate, error bound, and diagnostic information
/// about the octree subdivision process.
///
/// @par Example:
/// @code
/// auto result = IntegrateAdaptive3DCore(f, 0, 1, 0, 1, 0, 1, 1e-8);
/// if (result.converged) {
///     std::cout << "Integral = " << result.value 
///               << " ± " << result.error_estimate << "\n";
///     std::cout << "Subdivided " << result.cells_subdivided << " cells\n";
/// }
/// @endcode
struct AdaptiveResult3D : public IntegrationResult
{
    // value, error_estimate, function_evaluations, converged, status come from IntegrationResult
    int cells_subdivided = 0;     ///< Number of octree cell subdivisions performed
    int max_depth_reached = 0;    ///< Maximum subdivision depth reached

    /// @brief Default constructor
    AdaptiveResult3D(Real val = 0.0, Real err = 0.0, int evals = 0, 
                     int cells = 0, int depth = 0, bool conv = false)
        : IntegrationResult(val, err, depth, conv, /*error_estimated*/ true)
        , cells_subdivided(cells)
        , max_depth_reached(depth)
    {
        function_evaluations = evals;
        algorithm_name = "IntegrateAdaptive3D";
    }

    /// @brief Combine two results (for accumulating from child cells)
    AdaptiveResult3D& operator+=(const AdaptiveResult3D& other) {
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
    Real tol_abs = PrecisionValues<Real>::IntegrationTolerance;  ///< Absolute error tolerance
    Real tol_rel = PrecisionValues<Real>::IntegrationTolerance;  ///< Relative error tolerance
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
///          Error = |Kronrod - Gauss| is an embedded estimate, not a rigorous bound.
///
/// @note This is expensive (3375 evals vs 225 for 2D) but gives proper 3D error estimation.
template<typename Func>
static GKResult IntegrateCell3D_TensorGK15(Func f, Real x1, Real x2, Real y1, Real y2, Real z1, Real z2) {
    const Real cx = std::midpoint(x1, x2);
    const Real hx = IntegrationDetail::HalfWidth(x1, x2);
    const Real cy = std::midpoint(y1, y2);
    const Real hy = IntegrationDetail::HalfWidth(y1, y2);
    const Real cz = std::midpoint(z1, z2);
    const Real hz = IntegrationDetail::HalfWidth(z1, z2);

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
        Real wgi = i_gauss ? wG7[gauss_weight_index(i)] : MML::Real(0.0);
        
        for (int j = 0; j < GK15_N; ++j) {
            Real yj_pos = cy + hy * xGK15[j];
            Real yj_neg = cy - hy * xGK15[j];
            Real wj = wGK15[j];
            bool j_gauss = is_gauss_index(j);
            Real wgj = j_gauss ? wG7[gauss_weight_index(j)] : MML::Real(0.0);
            
            for (int k = 0; k < GK15_N; ++k) {
                Real zk_pos = cz + hz * xGK15[k];
                Real zk_neg = cz - hz * xGK15[k];
                Real wk = wGK15[k];
                bool k_gauss = is_gauss_index(k);
                Real wgk = k_gauss ? wG7[gauss_weight_index(k)] : MML::Real(0.0);
                
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
    result_kronrod = IntegrationDetail::Product({hx, hy, hz, result_kronrod});
    result_gauss = IntegrationDetail::Product({hx, hy, hz, result_gauss});

    Real error = std::abs(result_kronrod - result_gauss);

    return GKResult::Fixed(result_kronrod, error, evals);
}

} // namespace Detail

template<class Func>
static void RunAdaptive3D(Func& function, Real x1, Real x2, Real y1, Real y2, Real z1, Real z2,
    Real absolute, Real relative, int depth, int budget, GKRule rule, AdaptiveResult3D& result, bool checkFinite = true) {
    IntegrationDetail::HalfWidth(x1, x2);
    IntegrationDetail::HalfWidth(y1, y2);
    IntegrationDetail::HalfWidth(z1, z2);
    IntegrationDetail::ValidateAdaptive(absolute, relative, depth, budget);
    if (rule != GKRule::GK15)
        throw NotImplementedError("IntegrateAdaptive3D: only GKRule::GK15 is implemented");
    if (x2 <= x1 || y2 <= y1 || z2 <= z1) {
        result = AdaptiveResult3D(0, 0, 0, 0, 0, true);
        return;
    }
    using Box = std::array<std::array<Real, 2>, 3>;
    auto sample = [&](Real x, Real y, Real z) {
        return IntegrationDetail::Sample(function, result.function_evaluations, checkFinite, x, y, z);
    };
    auto evaluate = [&](const Box& box) {
        return Detail::IntegrateCell3D_TensorGK15(sample, box[0][0], box[0][1], box[1][0], box[1][1], box[2][0], box[2][1]);
    };
    const Box domain{{{x1, x2}, {y1, y2}, {z1, z2}}};
    IntegrationDetail::RefineCells(result, std::vector<Box>{domain}, evaluate,
        IntegrationDetail::SplitBox<3>, absolute, relative, depth, budget, 3375);
}

/// @brief Compatibility entry for adaptive 3D integration using global frontier refinement
/// @param f Function f(x,y,z) -> Real
/// @param x1,x2,y1,y2,z1,z2 Cell bounds
/// @param total_volume Legacy parameter; the requested domain now has one global error target
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
    AdaptiveResult3D result;
    RunAdaptive3D(f, x1, x2, y1, y2, z1, z2, tol_abs, tol_rel, max_depth - depth, evals_remaining, rule, result);
    evals_remaining -= result.function_evaluations;
    result.max_depth_reached += depth;
    result.iterations = result.max_depth_reached;
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
/// auto result = IntegrateAdaptive3DCore(
///     [](Real x, Real y, Real z) { return 1.0; },
///     0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1e-10
/// );
/// std::cout << result.value << " (should be 1.0)\n";
/// 
/// // Sphere indicator function: volume of unit sphere = 4π/3
/// auto sphere = IntegrateAdaptive3DCore(
///     [](Real x, Real y, Real z) { 
///         return (x*x + y*y + z*z <= 1.0) ? 1.0 : 0.0; 
///     },
///     -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, 1e-4
/// );
/// std::cout << "Sphere volume: " << sphere.value << " (4π/3 ≈ 4.189)\n";
/// @endcode
template<typename Func>
static AdaptiveResult3D IntegrateAdaptive3DCore(
    Func f,
    Real x1, Real x2,
    Real y1, Real y2,
    Real z1, Real z2,
    Real tolerance = PrecisionValues<Real>::IntegrationTolerance,
    int max_depth = 15,
    int max_evals = 10000000)
{
    AdaptiveResult3D result;
    RunAdaptive3D(f, x1, x2, y1, y2, z1, z2, tolerance, tolerance, max_depth, max_evals, GKRule::GK15, result);
    return result;
}

/// Result-returning variant with the plain argument list (RESULT_POLICY §4).
template<typename Func>
static auto IntegrateAdaptive3DDetailed(Func f, Real x1, Real x2, Real y1, Real y2, Real z1, Real z2, Real tolerance, int max_depth = 15, int max_evals = 10000000)
{
	auto r = IntegrateAdaptive3DCore(f, x1, x2, y1, y2, z1, z2, tolerance, max_depth, max_evals);
	r.algorithm_name = "IntegrateAdaptive3D";
	return r;
}

template<typename Func>
static Real IntegrateAdaptive3D(
    Func f,
    Real x1, Real x2,
    Real y1, Real y2,
    Real z1, Real z2,
    Real tolerance = PrecisionValues<Real>::IntegrationTolerance,
    int max_depth = 15,
    int max_evals = 10000000)
{
	return IntegrateAdaptive3DCore(f, x1, x2, y1, y2, z1, z2, tolerance, max_depth, max_evals).ValueOrThrow("IntegrateAdaptive3D");
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
/// auto result = IntegrateAdaptive3DCore(f, 0, 1, 0, 1, 0, 1, config);
/// @endcode
template<typename Func>
static AdaptiveResult3D IntegrateAdaptive3DCore(
    Func f,
    Real x1, Real x2,
    Real y1, Real y2,
    Real z1, Real z2,
    const AdaptiveConfig3D& config)
{
    AdaptiveResult3D result;
    RunAdaptive3D(f, x1, x2, y1, y2, z1, z2, config.tol_abs, config.tol_rel, config.max_depth, config.max_evaluations, config.rule, result);
    return result;
}

template<typename Func>
static Real IntegrateAdaptive3D(
    Func f,
    Real x1, Real x2,
    Real y1, Real y2,
    Real z1, Real z2,
    const AdaptiveConfig3D& config)
{
	return IntegrateAdaptive3DCore(f, x1, x2, y1, y2, z1, z2, config).ValueOrThrow("IntegrateAdaptive3D");
}

/**************************************************************************/
/*****       Detailed API - Adaptive 3D Integration                   *****/
/**************************************************************************/

/// Adaptive 3D integration with simple parameters and full diagnostics
template<typename Func>
static AdaptiveResult3D IntegrateAdaptive3DDetailed(
    Func f,
    Real x1, Real x2,
    Real y1, Real y2,
    Real z1, Real z2,
    const IntegrationConfig& config = {},
    Real tolerance = PrecisionValues<Real>::IntegrationTolerance,
    int max_depth = 15,
    int max_evals = 10000000)
{
    return IntegrationDetail::ExecuteIntegrationDetailed<AdaptiveResult3D>(
        "IntegrateAdaptive3D", config,
        [&](AdaptiveResult3D& result) {
            RunAdaptive3D(f, x1, x2, y1, y2, z1, z2, tolerance, tolerance, max_depth, max_evals, GKRule::GK15, result, config.check_finite);
        });
}

/// Adaptive 3D integration with config object and full diagnostics
template<typename Func>
static AdaptiveResult3D IntegrateAdaptive3DDetailed(
    Func f,
    Real x1, Real x2,
    Real y1, Real y2,
    Real z1, Real z2,
    const AdaptiveConfig3D& adaptive_config,
    const IntegrationConfig& config = {})
{
    return IntegrationDetail::ExecuteIntegrationDetailed<AdaptiveResult3D>(
        "IntegrateAdaptive3D", config,
        [&](AdaptiveResult3D& result) {
            RunAdaptive3D(f, x1, x2, y1, y2, z1, z2, adaptive_config.tol_abs, adaptive_config.tol_rel, adaptive_config.max_depth,
                adaptive_config.max_evaluations, adaptive_config.rule, result, config.check_finite);
        });
}

} // namespace Integration
} // namespace MML

#endif // MML_INTEGRATION_3D_ADAPTIVE_H

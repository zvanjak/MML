#ifndef MML_TEST_PRECISION_H
#define MML_TEST_PRECISION_H

#include "MMLBase.h"
#include <type_traits>
#include <cmath>
#include <limits>

namespace MML {
namespace Testing {

// ============================================================================
// Precision Detection and Type Traits
// ============================================================================

/// @brief Detect which Real type is currently active
enum class RealPrecision {
    Float,       // 32-bit, ~7 decimal digits
    Double,      // 64-bit, ~15 decimal digits  
    LongDouble   // 80/128-bit, ~18-33 decimal digits
};

/// @brief Get the current Real precision level
constexpr RealPrecision GetRealPrecision() {
    if constexpr (std::is_same_v<Real, float>) {
        return RealPrecision::Float;
    } else if constexpr (std::is_same_v<Real, long double>) {
        return RealPrecision::LongDouble;
    } else {
        return RealPrecision::Double;
    }
}

/// @brief Get precision name as string (for diagnostics)
constexpr const char* GetRealPrecisionName() {
    if constexpr (std::is_same_v<Real, float>) {
        return "float";
    } else if constexpr (std::is_same_v<Real, long double>) {
        return "long double";
    } else {
        return "double";
    }
}

/// @brief Get number of decimal digits of precision
constexpr int GetRealDigits() {
    return std::numeric_limits<Real>::digits10;
}

// ============================================================================
// Tolerance Scaling Utilities
// ============================================================================

/// @brief True when long double has genuinely more precision than double.
/// On GCC/Linux: long double is 80-bit (64 mantissa bits) vs double's 53 → true
/// On MSVC/Windows: long double == double (both 53 mantissa bits) → false
constexpr bool LongDoubleHasExtraPrecision =
    std::numeric_limits<long double>::digits > std::numeric_limits<double>::digits;

/// @brief Scale tolerance based on Real type precision
/// @param baseline_tol Tolerance value calibrated for double precision
/// @return Adjusted tolerance appropriate for current Real type
constexpr Real ScaleTolerance(Real baseline_tol) {
    if constexpr (std::is_same_v<Real, float>) {
        // Float has ~7 digits vs double's ~15 digits
        // Scale up tolerance by ~100x for numerical stability
        return baseline_tol * Real(100.0);
    } else if constexpr (std::is_same_v<Real, long double> && LongDoubleHasExtraPrecision) {
        // Long double has ~18-33 digits vs double's ~15 digits
        // Can use tighter tolerance (scale down by ~100x)
        return baseline_tol * Real(0.01);
    } else {
        // Double (or long double with no extra precision on MSVC) - use baseline as-is
        return baseline_tol;
    }
}

/// @brief Scale absolute tolerance based on value magnitude and precision
/// @param value The value being compared
/// @param baseline_tol Baseline absolute tolerance for double
/// @return Scaled absolute tolerance
inline Real ScaleAbsTolerance(Real value, Real baseline_tol) {
    Real scaled_tol = ScaleTolerance(baseline_tol);
    
    // For very small values, use direct tolerance
    if (std::abs(value) < Real(1e-10)) {
        return scaled_tol;
    }
    
    // For larger values, scale with magnitude
    return scaled_tol * std::max(Real(1.0), std::abs(value));
}

/// @brief Scale relative tolerance based on precision
/// @param baseline_rel_tol Baseline relative tolerance for double
/// @return Scaled relative tolerance  
constexpr Real ScaleRelTolerance(Real baseline_rel_tol) {
    return ScaleTolerance(baseline_rel_tol);
}

// ============================================================================
// Explicit Per-Precision Tolerances (NO automatic scaling)
// ============================================================================

/// @brief Select an explicit tolerance value based on current Real precision type.
/// Each tolerance is manually chosen per precision - no automatic scaling formula.
/// @param double_tol  Tolerance for double precision (baseline, well-tested)
/// @param float_tol   Tolerance for float precision (explicitly chosen per test)
/// @param ldouble_tol Tolerance for long double precision
/// @note When long double has no extra precision (MSVC), falls back to double_tol
constexpr Real Tol(Real double_tol, Real float_tol, Real ldouble_tol) {
    if constexpr (std::is_same_v<Real, float>) return float_tol;
    else if constexpr (std::is_same_v<Real, long double> && LongDoubleHasExtraPrecision)
        return ldouble_tol;
    else return double_tol;
}

/// @brief 2-arg version: long double defaults to same as double (refine later)
constexpr Real Tol(Real double_tol, Real float_tol) {
    if constexpr (std::is_same_v<Real, float>) return float_tol;
    else return double_tol;
}

// ============================================================================
// Common Test Tolerances (Pre-scaled for Current Precision)
// ============================================================================

namespace Tolerance {
    // General comparison tolerances
    constexpr Real Loose      = ScaleTolerance(Real(1e-6));   // Loose comparison
    constexpr Real Standard   = ScaleTolerance(Real(1e-9));   // Standard tests
    constexpr Real Strict     = ScaleTolerance(Real(1e-12));  // Strict requirements
    constexpr Real VeryStrict = ScaleTolerance(Real(1e-14));  // Very strict (near machine precision)
    
    // Operation-specific tolerances
    namespace Integration {
        constexpr Real Standard = ScaleTolerance(Real(1e-9));
        constexpr Real Adaptive = ScaleTolerance(Real(1e-11));
        constexpr Real Strict   = ScaleTolerance(Real(1e-12));
    }
    
    namespace Differentiation {
        constexpr Real Standard = ScaleTolerance(Real(1e-11));
        constexpr Real Strict   = ScaleTolerance(Real(1e-13));
    }
    
    namespace LinearAlgebra {
        constexpr Real Standard = ScaleTolerance(Real(1e-10));
        constexpr Real Strict   = ScaleTolerance(Real(1e-12));
        constexpr Real SVD      = ScaleTolerance(Real(1e-9));  // SVD is less precise
    }
    
    namespace RootFinding {
        constexpr Real Standard = ScaleTolerance(Real(1e-12));
        constexpr Real Strict   = ScaleTolerance(Real(1e-13));
    }
    
    namespace Geometry {
        constexpr Real Standard = ScaleTolerance(Real(1e-10));
        constexpr Real Strict   = ScaleTolerance(Real(1e-12));
    }
    
    namespace Statistics {
        constexpr Real Standard = ScaleTolerance(Real(1e-8));  // Statistics often less precise
        constexpr Real MonteCarlo = ScaleTolerance(Real(1e-3)); // Monte Carlo has inherent error
    }
}

// ============================================================================
// Helper Functions for Test Assertions
// ============================================================================

/// @brief Check if two Real values are approximately equal
/// @param a First value
/// @param b Second value  
/// @param rel_tol Relative tolerance (will be scaled for precision)
/// @param abs_tol Absolute tolerance (will be scaled for precision)
inline bool ApproxEqual(Real a, Real b, 
                       Real rel_tol = Tolerance::Standard,
                       Real abs_tol = Tolerance::Standard) {
    // Handle exact equality (including infinities)
    if (a == b) return true;
    
    // Handle NaN
    if (std::isnan(a) || std::isnan(b)) return false;
    
    // Absolute difference
    Real diff = std::abs(a - b);
    
    // Check absolute tolerance
    if (diff <= abs_tol) return true;
    
    // Check relative tolerance
    Real larger = std::max(std::abs(a), std::abs(b));
    return diff <= rel_tol * larger;
}

/// @brief Check if value is approximately zero
/// @param value Value to check
/// @param abs_tol Absolute tolerance
inline bool ApproxZero(Real value, Real abs_tol = Tolerance::Standard) {
    return std::abs(value) <= abs_tol;
}

/// @brief Get machine epsilon for current Real type
constexpr Real MachineEpsilon() {
    return std::numeric_limits<Real>::epsilon();
}

/// @brief Get smallest normalized positive value for current Real type
constexpr Real MinPositive() {
    return std::numeric_limits<Real>::min();
}

/// @brief Get maximum finite value for current Real type  
constexpr Real MaxFinite() {
    return std::numeric_limits<Real>::max();
}

// ============================================================================
// Diagnostic Helpers
// ============================================================================

/// @brief Get diagnostic string about current precision configuration
inline std::string GetPrecisionInfo() {
    std::ostringstream oss;
    oss << "Real type: " << GetRealPrecisionName() << "\n";
    oss << "Size: " << sizeof(Real) << " bytes\n";
    oss << "Digits: " << GetRealDigits() << " decimal\n";
    oss << "Epsilon: " << MachineEpsilon() << "\n";
    oss << "Min: " << MinPositive() << "\n";
    oss << "Max: " << MaxFinite() << "\n";
    return oss.str();
}

} // namespace Testing
} // namespace MML

// ============================================================================
// Convenience Macros for Common Test Patterns
// ============================================================================

/// @brief Per-precision tolerance: explicit double and float values, no auto-scaling
/// Usage: TOL(1e-10, 1e-4) returns 1e-10 for double, 1e-4 for float
#define TOL(DOUBLE_VAL, FLOAT_VAL) MML::Testing::Tol(REAL(DOUBLE_VAL), REAL(FLOAT_VAL))

/// @brief Per-precision tolerance with all three precision levels specified
#define TOL3(DOUBLE_VAL, FLOAT_VAL, LDOUBLE_VAL) \
    MML::Testing::Tol(REAL(DOUBLE_VAL), REAL(FLOAT_VAL), REAL(LDOUBLE_VAL))

/// @brief Print precision info once per test suite (useful for debugging)
#define TEST_PRECISION_INFO() \
    INFO("Running with " << MML::Testing::GetRealPrecisionName() \
         << " precision (" << MML::Testing::GetRealDigits() << " digits)")

// ============================================================================
// Long Double Support for Catch2 WithinRel Matcher
// ============================================================================
// Catch2 only provides WithinRel overloads for double and float.
// When Real = long double, calls like WithinRel(REAL(1.0), REAL(1e-5)) are
// ambiguous. We add explicit long double overloads that cast to double.
// On MSVC, long double == double so this is lossless. On GCC/Linux (80-bit
// extended), the cast loses some precision, but since we're comparing with
// tolerances (typically >= 1e-15), this is acceptable.

#ifdef MML_USE_LONG_DOUBLE
#include <catch2/matchers/catch_matchers_floating_point.hpp>
namespace Catch { namespace Matchers {
    inline WithinRelMatcher WithinRel(long double target, long double eps) {
        return WithinRel(static_cast<double>(target), static_cast<double>(eps));
    }
    inline WithinRelMatcher WithinRel(long double target) {
        return WithinRel(static_cast<double>(target));
    }
}} // namespace Catch::Matchers
#endif

#endif // MML_TEST_PRECISION_H

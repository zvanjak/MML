///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        MMLBase.h                                                           ///
///  Description: Core definitions, constants, type aliases, and precision settings   ///
///               Foundation header included by all MML components                    ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_BASE_H
#define MML_BASE_H

#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1

// MML headers first (catches missing includes in them)
#include "MMLExceptions.h"
#include "MMLPrecision.h"
#include "MMLVisualizators.h"

// Standard headers - only what MMLBase.h actually uses
#include <cmath>
#include <complex>
#include <limits>
#include <type_traits>

///////////////////////////////////////////////////////////////////////////////////////////
// Real Type Configuration
///////////////////////////////////////////////////////////////////////////////////////////
//
// CONFIGURATION MODEL: Per-Library Build (NOT Per-Translation-Unit)
//
// The Real typedef defines the floating-point precision for the entire library.
// This is a **build-time configuration** that must be consistent across ALL translation
// units that use MML together in the same program.
//
// USAGE:
//   1. Choose ONE Real type for your build (double, float, long double, or __float128)
//   2. Ensure ALL source files in your project see the SAME definition
//   3. Do NOT mix different Real types in the same executable
//
// CURRENT CONFIGURATION:
typedef double						 Real; // default real type
//
// OTHER SUPPORTED OPTIONS (uncomment ONE, comment others):
// typedef float						 Real;    // Lower precision, faster, smaller memory
// typedef long double       Real;    // Extended precision (80-bit on x86, 128-bit on some platforms)
// typedef __float128				 Real;    // Quad precision (GCC only, 128-bit IEEE 754)
//
// ABI AND COMPATIBILITY GUARANTEES:
//
// 1. BINARY COMPATIBILITY:
//    - Changing Real breaks ABI compatibility
//    - All libraries/object files must be recompiled with the same Real type
//    - Linking objects built with different Real types causes undefined behavior
//
// 2. SERIALIZATION:
//    - Data serialized with one Real type may NOT be readable with another
//    - Binary format depends on Real's size and representation
//    - Always document which Real type was used when saving data
//    - Consider text-based formats for cross-precision compatibility
//
// 3. PRECISION CONSTANTS:
//    - Constants in MML::Constants are stored as 'double' literals
//    - For float builds: implicit conversion (acceptable precision loss)
//    - For long double/__float128: constants may have less precision than Real
//    - Use REAL() macro for literals that should match Real precision
//
// 4. API CONTRACTS:
//    - All MML APIs use Real consistently
//    - Template parameters like Vector<Real>, Matrix<Real> must match the global Real
//    - Mixing Real with explicit float/double types is supported but requires care
//
// 5. THREAD SAFETY:
//    - Real type is compile-time constant (no runtime changes)
//    - Safe to use across multiple threads
//    - Thread-local contexts (Defaults) are independent of Real choice
//
// RECOMMENDATIONS:
//   - Use 'double' (default) for most applications (good balance of speed/precision)
//   - Use 'float' when memory or performance is critical and precision allows
//   - Use 'long double' for extended precision needs (note: compiler/platform dependent)
//   - Use '__float128' for maximum precision (GCC only, slower, requires libquadmath)
//
// PORTABILITY NOTES:
//   - 'double' and 'float' are fully portable (C++ standard)
//   - 'long double' precision varies by platform (80-bit x86, 128-bit ARM/PowerPC, 64-bit MSVC)
//   - '__float128' requires GCC and libquadmath linkage
//
///////////////////////////////////////////////////////////////////////////////////////////

// Macro for type-safe numeric literals that match Real type
// This ensures literals like 0.0, 1.0 match the current Real precision
// Usage: Real x = REAL(3.14159265358979323846);
#ifndef REAL
#define REAL(x) static_cast<Real>(x)
#endif

// Complex must have the same underlying type as Real
// Changing Real automatically changes Complex precision
typedef std::complex<Real> Complex; // default complex type

namespace MML 
{

  // Generic absolute value functions
  template <class Type> static Real Abs(const Type &a) { return std::abs(a); }
  template <class Type> static Real Abs(const std::complex<Type> &a) {
    return hypot(a.real(), a.imag());
  }

  ////////////         Floating-Point Comparison Functions         ////////////
  /// @brief Check if two values are equal within absolute tolerance.
  /// @details Use for values expected to be near zero.
  inline bool isWithinAbsPrec(Real a, Real b, Real eps) {
    return std::abs(a - b) < eps;
  }

  /// @brief Check if two values are equal within relative tolerance.
  /// @details Use for values of similar magnitude away from zero.
  inline bool isWithinRelPrec(Real a, Real b, Real eps) {
    return std::abs(a - b) < eps * std::max(Abs(a), Abs(b));
  }

  /// @brief Check if two values are nearly equal using combined absolute and relative tolerance.
  /// @details This is the recommended comparison for general floating-point equality.
  ///          Handles both small values (where absolute tolerance dominates) and
  ///          large values (where relative tolerance dominates) correctly.
  /// @param a First value
  /// @param b Second value  
  /// @param absEps Absolute tolerance (dominates for values near zero)
  /// @param relEps Relative tolerance (dominates for large values)
  /// @return true if |a - b| < absEps + relEps * max(|a|, |b|)
  inline bool isNearlyEqual(Real a, Real b, Real absEps, Real relEps) {
    return std::abs(a - b) < absEps + relEps * std::max(Abs(a), Abs(b));
  }

  /// @brief Check if two values are nearly equal using a single tolerance for both abs and rel.
  /// @details Convenience overload that uses the same tolerance for absolute and relative comparison.
  inline bool isNearlyEqual(Real a, Real b, Real eps) {
    return isNearlyEqual(a, b, eps, eps);
  }

  /// @brief Check if a value is nearly zero (within tolerance of zero).
  /// @details Use for checking if a value is effectively zero in numerical computations.
  /// @param a Value to check
  /// @param eps Tolerance (defaults to NumericalZeroThreshold)
  /// @return true if |a| < eps
  inline bool isNearlyZero(Real a, Real eps = PrecisionValues<Real>::NumericalZeroThreshold) {
    return std::abs(a) < eps;
  }

  template <class T> inline T POW2(const T &a) {
    const T &t = a;
    return t * t;
  }
  template <class T> inline T POW3(const T &a) {
    const T &t = a;
    return t * t * t;
  }
  template <class T> inline T POW4(const T &a) {
    const T &t = a * a;
    return t * t;
  }
  template <class T> inline T POW5(const T &a) {
    const T &t = a;
    return t * t * t * t * t;
  }

  ////////////                  Constants                ////////////////
  namespace Constants {
    // Mathematical constants with full long double precision
    // When Real=double, implicit narrowing is lossless
    // When Real=long double, these preserve maximum available precision
    static inline constexpr long double PI         = 3.141592653589793238462643383279502884L;  // pi
    static inline constexpr long double INV_PI     = 0.318309886183790671537767526745028724L;  // 1/pi
    static inline constexpr long double INV_SQRTPI = 0.564189583547756286948079451560772586L;  // 1/sqrt(pi)

    static inline constexpr long double E          = 2.718281828459045235360287471352662498L;  // e
    static inline constexpr long double LN2        = 0.693147180559945309417232121458176568L;  // ln(2)
    static inline constexpr long double LN10       = 2.302585092994045684017991454684364208L;  // ln(10)

    static inline constexpr long double SQRT2      = 1.414213562373095048801688724209698079L;  // sqrt(2)
    static inline constexpr long double SQRT3      = 1.732050807568877293527446341505872367L;  // sqrt(3)

    static inline constexpr long double GoldenRatio = 1.618033988749894848204586834365638118L; // (1 + sqrt(5)) / 2

    // Geometry epsilon for floating-point comparisons in geometric algorithms
    static inline constexpr double GEOMETRY_EPSILON = 1e-10;

    // Precision constants - use Real type for consistency with library's floating-point type
    static inline const Real Eps = std::numeric_limits<Real>::epsilon();
    static inline const Real PosInf = std::numeric_limits<Real>::infinity();
    static inline const Real NegInf = -std::numeric_limits<Real>::infinity();
  } // namespace Constants

  ////////////       Binary File Format Constants        ////////////
  /// @brief Magic numbers and version constants for MML binary file formats.
  /// @details These constants identify MML binary files and their format versions.
  ///          Magic numbers are 4-byte ASCII identifiers stored as uint32_t.
  ///          All binary files use little-endian byte order.
  namespace BinaryFormat {
    // Magic numbers (4-byte ASCII identifiers)
    static inline constexpr uint32_t MAGIC_MATRIX         = 0x4D4D4C4D;  // "MMLM" - MML Matrix
    static inline constexpr uint32_t MAGIC_VECTOR         = 0x4D4D4C56;  // "MMLV" - MML Vector (real)
    static inline constexpr uint32_t MAGIC_VECTOR_COMPLEX = 0x4D4D4C43;  // "MMLC" - MML Vector Complex
    static inline constexpr uint32_t MAGIC_SPARSE         = 0x4D4D4C53;  // "MMLS" - MML Sparse Matrix (reserved)
    static inline constexpr uint32_t MAGIC_TENSOR         = 0x4D4D4C54;  // "MMLT" - MML Tensor (reserved)
    
    // Current format versions
    static inline constexpr uint32_t VERSION_MATRIX         = 1;
    static inline constexpr uint32_t VERSION_VECTOR         = 1;
    static inline constexpr uint32_t VERSION_VECTOR_COMPLEX = 1;
    static inline constexpr uint32_t VERSION_SPARSE         = 1;  // reserved
    static inline constexpr uint32_t VERSION_TENSOR         = 1;  // reserved
    
    // File extensions (without dot)
    static inline constexpr const char* EXT_MATRIX         = "mmlm";
    static inline constexpr const char* EXT_VECTOR         = "mmlv";
    static inline constexpr const char* EXT_VECTOR_COMPLEX = "mmlc";  // complex vector
    static inline constexpr const char* EXT_SPARSE         = "mmls";  // reserved
    static inline constexpr const char* EXT_TENSOR         = "mmlt";  // reserved
  } // namespace BinaryFormat

  ////////////         Angle Comparison Functions         ////////////
  /// @brief Normalize angle to [-π, π) range for comparison.
  /// @param rad Angle in radians
  /// @return Normalized angle in [-π, π)
  inline Real normalizeAngle(Real rad) {
    rad = std::fmod(rad + Constants::PI, 2 * Constants::PI);
    if (rad < 0)
      rad += 2 * Constants::PI;
    return rad - Constants::PI;
  }

  /// @brief Check if two angles are equal, accounting for wrap-around at ±π.
  /// @details Normalizes the difference to [-π, π) before comparing.
  ///          Correctly handles cases like comparing -π and π (which are equal).
  /// @param a First angle in radians
  /// @param b Second angle in radians
  /// @param eps Tolerance for comparison
  /// @return true if angles are equivalent within tolerance
  inline bool AnglesAreEqual(Real a, Real b, Real eps) {
    Real diff = normalizeAngle(a - b);
    return std::abs(diff) < eps;
  }

  // is_simple_numeric helper
  template <typename T> struct is_simple_numeric : std::is_arithmetic<T> {};

  template <typename T>
  struct is_simple_numeric<std::complex<T>> : std::is_arithmetic<T> {};

  // Helper variable template (C++14 and later)
  template <typename T>
  inline constexpr bool is_MML_simple_numeric = is_simple_numeric<T>::value;


  struct AlgorithmContext {
    // Integration parameters (using PrecisionValues for consistency)
    Real trapezoidIntegrationEPS = REAL(1.0e-4);
    Real simpsonIntegrationEPS = REAL(1.0e-5);
    Real rombergIntegrationEPS = PrecisionValues<Real>::DefaultTolerance;
    Real workIntegralPrecision = REAL(1e-05);
    Real lineIntegralPrecision = REAL(1e-05);

    int trapezoidIntegrationMaxSteps = 20;
    int simpsonIntegrationMaxSteps = 20;
    int rombergIntegrationMaxSteps = 20;
    int rombergIntegrationUsedPnts = 5;

    // Root finding parameters
    int bisectionMaxSteps = 50;
    int newtonRaphsonMaxSteps = 20;
    int brentMaxSteps = 100;  // Brent typically needs more steps but converges reliably

    // ODE solver parameters
    int odeSolverMaxSteps = 100000;

    static AlgorithmContext &Get() {
      thread_local static AlgorithmContext ctx;
      return ctx;
    }
  };

  // Thread-safe configuration contexts
  struct PrintContext {
    int vectorWidth = 15;
    int vectorPrecision = 10;
    int vectorNWidth = 15;
    int vectorNPrecision = 10;

    static PrintContext &Get() {
      thread_local static PrintContext ctx;
      return ctx;
    }
  };  

  // Backward compatible Defaults namespace (now thread-safe via thread_local
  // contexts)
  namespace Defaults {
    // Output defaults (thread-safe - changed from static globals to thread_local)
    // Usage: Defaults::VectorPrintWidth = 20; (now thread-safe!)
    static inline int &VectorPrintWidth = PrintContext::Get().vectorWidth;
    static inline int &VectorPrintPrecision = PrintContext::Get().vectorPrecision;
    static inline int &VectorNPrintWidth = PrintContext::Get().vectorNWidth;
    static inline int &VectorNPrintPrecision = PrintContext::Get().vectorNPrecision;

    //////////               Default precisions             ///////////
    // Use the precision values based on the Real type
    static inline const double ComplexAreEqualTolerance =
        PrecisionValues<Real>::ComplexAreEqualTolerance;
    static inline const double ComplexAreEqualAbsTolerance =
        PrecisionValues<Real>::ComplexAreEqualAbsTolerance;
    static inline const double VectorIsEqualTolerance =
        PrecisionValues<Real>::VectorIsEqualTolerance;
    static inline const double MatrixIsEqualTolerance =
        PrecisionValues<Real>::MatrixIsEqualTolerance;

    static inline const double Pnt2CartIsEqualTolerance =
        PrecisionValues<Real>::Pnt2CartIsEqualTolerance;
    static inline const double Pnt2PolarIsEqualTolerance =
        PrecisionValues<Real>::Pnt2PolarIsEqualTolerance;
    static inline const double Pnt3CartIsEqualTolerance =
        PrecisionValues<Real>::Pnt3CartIsEqualTolerance;
    static inline const double Pnt3SphIsEqualTolerance =
        PrecisionValues<Real>::Pnt3SphIsEqualTolerance;
    static inline const double Pnt3CylIsEqualTolerance =
        PrecisionValues<Real>::Pnt3CylIsEqualTolerance;

    static inline const double Vec2CartIsEqualTolerance =
        PrecisionValues<Real>::Vec2CartIsEqualTolerance;
    static inline const double Vec3CartIsEqualTolerance =
        PrecisionValues<Real>::Vec3CartIsEqualTolerance;
    static inline const double Vec3CartIsParallelTolerance =
        PrecisionValues<Real>::Vec3CartIsParallelTolerance;

    static inline const double Vec3SphIsEqualTolerance =
        PrecisionValues<Real>::Vec3SphIsEqualTolerance;

    // Angle comparison tolerance (for wrap-aware angle equality)
    static inline const double AngleIsEqualTolerance =
        PrecisionValues<Real>::AngleIsEqualTolerance;

    // Shape property tolerance (for geometric shape classification)
    static inline const double ShapePropertyTolerance =
        PrecisionValues<Real>::ShapePropertyTolerance;

    static inline const double Line3DAreEqualTolerance =
        PrecisionValues<Real>::Line3DAreEqualTolerance;
    static inline const double Line3DIsPointOnLineTolerance =
        PrecisionValues<Real>::Line3DIsPointOnLineTolerance;
    static inline const double Line3DIsPerpendicularTolerance =
        PrecisionValues<Real>::Line3DIsPerpendicularTolerance;
    static inline const double Line3DIsParallelTolerance =
        PrecisionValues<Real>::Line3DIsParallelTolerance;
    static inline const double Line3DIntersectionTolerance =
        PrecisionValues<Real>::Line3DIntersectionTolerance;

    static inline const double Plane3DIsPointOnPlaneTolerance =
        PrecisionValues<Real>::Plane3DIsPointOnPlaneTolerance;

    static inline const double Triangle3DIsPointInsideTolerance =
        PrecisionValues<Real>::Triangle3DIsPointInsideTolerance;
    static inline const double Triangle3DIsRightTolerance =
        PrecisionValues<Real>::Triangle3DIsRightTolerance;
    static inline const double Triangle3DIsIsoscelesTolerance =
        PrecisionValues<Real>::Triangle3DIsIsoscelesTolerance;
    static inline const double Triangle3DIsEquilateralTolerance =
        PrecisionValues<Real>::Triangle3DIsEquilateralTolerance;

    static inline const double IsMatrixSymmetricTolerance =
        PrecisionValues<Real>::IsMatrixSymmetricTolerance;
    static inline const double IsMatrixDiagonalTolerance =
        PrecisionValues<Real>::IsMatrixDiagonalTolerance;
    static inline const double IsMatrixUnitTolerance =
        PrecisionValues<Real>::IsMatrixUnitTolerance;
    static inline const double IsMatrixZeroTolerance =
        PrecisionValues<Real>::IsMatrixZeroTolerance;
    static inline const double IsMatrixOrthogonalTolerance =
        PrecisionValues<Real>::IsMatrixOrthogonalTolerance;

    static inline const double RankAlgEPS = PrecisionValues<Real>::RankAlgEPS;

    // Numerical thresholds from PrecisionValues
    static inline const double DefaultTolerance =
        PrecisionValues<Real>::DefaultTolerance;
    static inline const double OrthogonalityTolerance =
        PrecisionValues<Real>::OrthogonalityTolerance;

    // Algorithm parameters (thread-safe - changed from static constants to
    // thread_local) Usage: Defaults::TrapezoidIntegrationEPS = 1e-6; (now
    // thread-safe and mutable!)
    static inline Real &TrapezoidIntegrationEPS =
        AlgorithmContext::Get().trapezoidIntegrationEPS;
    static inline Real &SimpsonIntegrationEPS =
        AlgorithmContext::Get().simpsonIntegrationEPS;
    static inline Real &RombergIntegrationEPS =
        AlgorithmContext::Get().rombergIntegrationEPS;
    static inline Real &WorkIntegralPrecision =
        AlgorithmContext::Get().workIntegralPrecision;
    static inline Real &LineIntegralPrecision =
        AlgorithmContext::Get().lineIntegralPrecision;

    static inline int &BisectionMaxSteps =
        AlgorithmContext::Get().bisectionMaxSteps;
    static inline int &NewtonRaphsonMaxSteps =
        AlgorithmContext::Get().newtonRaphsonMaxSteps;
    static inline int &BrentMaxSteps =
        AlgorithmContext::Get().brentMaxSteps;

    static inline int &TrapezoidIntegrationMaxSteps =
        AlgorithmContext::Get().trapezoidIntegrationMaxSteps;
    static inline int &SimpsonIntegrationMaxSteps =
        AlgorithmContext::Get().simpsonIntegrationMaxSteps;
    static inline int &RombergIntegrationMaxSteps =
        AlgorithmContext::Get().rombergIntegrationMaxSteps;
    static inline int &RombergIntegrationUsedPnts =
        AlgorithmContext::Get().rombergIntegrationUsedPnts;

    static inline int &ODESolverMaxSteps =
        AlgorithmContext::Get().odeSolverMaxSteps;
  } // namespace Defaults
} // namespace MML
#endif

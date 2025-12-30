///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        MMLBase.h                                                           ///
///  Description: Core definitions, constants, type aliases, and precision settings   ///
///               Foundation header included by all MML components                    ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
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

// HAJDUK ZIVI VJECNO!!!

typedef double						 Real; // default real type

// other possibilites:
// typedef float						 Real;
// typedef long double       Real;
// typedef __float128				 Real;    // only for GCC!

// Macro for type-safe numeric literals that match Real type
// This ensures literals like 0.0, 1.0 match the current Real precision
#ifndef REAL
#define REAL(x) static_cast<Real>(x)
#endif

// Complex must have the same underlaying type as Real
typedef std::complex<Real> Complex; // default complex type

namespace MML 
{

  // Generic absolute value functions
  template <class Type> static Real Abs(const Type &a) { return std::abs(a); }
  template <class Type> static Real Abs(const std::complex<Type> &a) {
    return hypot(a.real(), a.imag());
  }


  inline bool isWithinAbsPrec(Real a, Real b, Real eps) {
    return std::abs(a - b) < eps;
  }
  inline bool isWithinRelPrec(Real a, Real b, Real eps) {
    return std::abs(a - b) < eps * std::max(Abs(a), Abs(b));
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
    // Mathematical constants with maximum double precision (C++17 compatible)
    // These values match std::numbers from C++20 to full double precision
    static inline constexpr double PI         = 3.14159265358979323846;  // pi
    static inline constexpr double INV_PI     = 0.31830988618379067154;  // 1/pi
    static inline constexpr double INV_SQRTPI = 0.56418958354775628695;  // 1/sqrt(pi)

    static inline constexpr double E          = 2.71828182845904523536;  // e
    static inline constexpr double LN2        = 0.69314718055994530942;  // ln(2)
    static inline constexpr double LN10       = 2.30258509299404568402;  // ln(10)

    static inline constexpr double SQRT2      = 1.41421356237309504880;  // sqrt(2)
    static inline constexpr double SQRT3      = 1.73205080756887729353;  // sqrt(3)

    // Precision constants - use Real type for consistency with library's floating-point type
    static inline const Real Eps = std::numeric_limits<Real>::epsilon();
    static inline const Real PosInf = std::numeric_limits<Real>::infinity();
    static inline const Real NegInf = -std::numeric_limits<Real>::infinity();
  } // namespace Constants

  // a * x^2 + b * x + c = 0
  inline int SolveQuadratic(Real a, Real b, Real c, Complex &x1, Complex &x2) {
    Real D = b * b - 4 * a * c;
    if (D >= 0) {
      Real sqrtD = sqrt(D);
      x1 = (-b + sqrtD) / (2 * a);
      x2 = (-b - sqrtD) / (2 * a);
      return 2;
    } else {
      Complex sqrtD = std::sqrt(Complex(D));
      x1 = (-b + sqrtD) / (2 * a);
      x2 = (-b - sqrtD) / (2 * a);
      return 0;
    }
  }
  inline void SolveQuadratic(const Complex &a, const Complex &b, const Complex &c,
                            Complex &x1, Complex &x2) {
    Complex D = b * b - Real(4.0) * a * c;
    Complex sqrtD = std::sqrt(D);
    x1 = (-b + sqrtD) / (Real(2.0) * a);
    x2 = (-b - sqrtD) / (Real(2.0) * a);
  }
  // Solving cubic equation a * x^3 + b * x^2 + c * x + d = 0
  inline int SolveCubic(Real a, Real b, Real c, Real d, Complex &x1, Complex &x2,
                        Complex &x3) {
    // Normalize the coefficients
    Real A = b / a;
    Real B = c / a;
    Real C = d / a;

    // Calculate the discriminant
    Real Q = (Real(3.0) * B - POW2(A)) / Real(9.0);
    Real R =
        (Real(9.0) * A * B - Real(27.0) * C - Real(2.0) * POW3(A)) / Real(54.0);
    Real D = POW3(Q) + POW2(R); // Discriminant

    if (D >= 0) // Complex or duplicate roots
    {
      Real S = std::cbrt(R + std::sqrt(D));
      Real T = std::cbrt(R - std::sqrt(D));

      x1 = -A / Real(3.0) + (S + T); // Real root
      x2 = -A / Real(3.0) - (S + T) / Real(2.0) +
          Complex(0, std::sqrt(Real(3.0)) * (S - T) / Real(2.0)); // Complex root
      x3 = -A / Real(3.0) - (S + T) / Real(2.0) -
          Complex(0, std::sqrt(Real(3.0)) * (S - T) / Real(2.0)); // Complex root

      return 1;
    } else // Three real roots
    {
      Real theta = std::acos(R / std::sqrt(-POW3(Q)));
      x1 =
          Real(2.0) * std::sqrt(-Q) * std::cos(theta / Real(3.0)) - A / Real(3.0);
      x2 = Real(2.0) * std::sqrt(-Q) *
              std::cos((theta + Real(2.0) * Constants::PI) / Real(3.0)) -
          A / Real(3.0);
      x3 = Real(2.0) * std::sqrt(-Q) *
              std::cos((theta + Real(4.0) * Constants::PI) / Real(3.0)) -
          A / Real(3.0);

      return 3;
    }
  }
  // Solving quartic equation a * x^4 + b * x^3 + c * x^2 + d * x + e = 0
  inline void SolveQuartic(Real a, Real b, Real c, Real d, Real e, Complex &x1,
                          Complex &x2, Complex &x3, Complex &x4) {
    // Degenerate: reduce to cubic
    if (std::abs(a) < Constants::Eps) {
      SolveCubic(b, c, d, e, x1, x2, x3);
      x4 = Complex(0);
      return;
    }

    // Normalize coefficients
    Real A = b / a;
    Real B = c / a;
    Real C = d / a;
    Real D = e / a;

    // Depressed quartic y = x + A/4: y^4 + P y^2 + Q y + R = 0
    Real AA = A * A;
    Real P = B - Real(3.0) * AA / Real(8.0);
    Real Q = C - Real(0.5) * A * B + AA * A / Real(8.0);
    Real R = D - Real(0.25) * A * C + AA * B / Real(16.0) -
            Real(3.0) * AA * AA / Real(256.0);

    // Special case: biquadratic (Q ≈ 0) -> solve t^2 + P t + R = 0 where t = y^2
    if (std::abs(Q) <= Constants::Eps) {
      Complex t1, t2;
      SolveQuadratic(Complex(Real(1.0)), Complex(P), Complex(R), t1, t2);

      x1 = std::sqrt(t1) - A / Real(4.0);
      x2 = -std::sqrt(t1) - A / Real(4.0);
      x3 = std::sqrt(t2) - A / Real(4.0);
      x4 = -std::sqrt(t2) - A / Real(4.0);
      return;
    }

    // General case (Ferrari)
    Complex z1, z2, z3;
    SolveCubic(Real(1.0), -P / Real(2.0), -R,
              R * P / Real(2.0) - Q * Q / Real(8.0), z1, z2, z3);

    auto U_from = [P](const Complex &z) {
      return std::sqrt(Complex(Real(2.0)) * z - P);
    };

    // Choose z to maximize |U| to avoid division by small numbers
    Complex candidates[3] = {z1, z2, z3};
    Complex z = candidates[0];
    Complex U = U_from(z);
    for (int i = 1; i < 3; ++i) {
      Complex Ui = U_from(candidates[i]);
      if (std::abs(Ui) > std::abs(U)) {
        z = candidates[i];
        U = Ui;
      }
    }

    Complex W = Q / (Complex(Real(2.0)) * U);

    // Solve two quadratics in y
    Complex y1, y2, y3, y4;
    SolveQuadratic(Complex(Real(1.0)), U, z - W, y1, y2);
    SolveQuadratic(Complex(Real(1.0)), -U, z + W, y3, y4);

    // Back-substitute x = y - A/4
    x1 = y1 - A / Real(4.0);
    x2 = y2 - A / Real(4.0);
    x3 = y3 - A / Real(4.0);
    x4 = y4 - A / Real(4.0);
  }

  // is_simple_numeric helper
  template <typename T> struct is_simple_numeric : std::is_arithmetic<T> {};

  template <typename T>
  struct is_simple_numeric<std::complex<T>> : std::is_arithmetic<T> {};

  // Helper variable template (C++14 and later)
  template <typename T>
  inline constexpr bool is_MML_simple_numeric = is_simple_numeric<T>::value;

  ///////////////////////////////////////////////////////////////////////////////
  //                      ALGORITHM RESULT STRUCTURES
  ///////////////////////////////////////////////////////////////////////////////

  /// Result structure for numerical integration algorithms
  /// Provides convergence status and diagnostics
  /// @note For production code, check error_estimate and converged fields!
  struct IntegrationResult {
    Real value;          ///< Computed integral value
    Real error_estimate; ///< Estimated absolute error
    int iterations;      ///< Number of iterations/refinements performed
    bool converged;      ///< True if convergence criteria met

    /// Implicit conversion to Real for backward compatibility
    /// @warning Silently discards error_estimate, iterations, and converged!
    /// @deprecated Prefer explicit .value access in new code
    operator Real() const { return value; }

    /// Constructor for easy initialization
    IntegrationResult(Real val = 0.0, Real err = 0.0, int iter = 0,
                      bool conv = true)
        : value(val), error_estimate(err), iterations(iter), converged(conv) {}
  };

  /// Result structure for root finding algorithms
  /// Provides convergence status and diagnostics
  /// @note For production code, check function_value and converged fields!
  struct RootFindingResult {
    Real root;           ///< Found root value
    Real function_value; ///< f(root) - should be near zero
    int iterations;      ///< Number of iterations performed
    bool converged;      ///< True if convergence criteria met

    /// Implicit conversion to Real for backward compatibility
    /// @warning Silently discards function_value, iterations, and converged!
    /// @deprecated Prefer explicit .root access in new code
    operator Real() const { return root; }

    RootFindingResult(Real r = 0.0, Real fval = 0.0, int iter = 0,
                      bool conv = true)
        : root(r), function_value(fval), iterations(iter), converged(conv) {}
  };

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
    static inline const double IsMatrixOrthogonalTolerance =
        PrecisionValues<Real>::IsMatrixOrthogonalTolerance;

    static inline const double RankAlgEPS = PrecisionValues<Real>::RankAlgEPS;

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

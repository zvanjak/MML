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
#include <mml/MMLTypeDefs.h>
#include <mml/MMLExceptions.h>
#include <mml/MMLPrecision.h>
#include <mml/MMLConcepts.h>

// Standard headers - only what MMLBase.h actually uses
#include <cmath>
#include <complex>
#include <cstdint>
#include <limits>

// Macro for type-safe numeric literals that match Real type
// This ensures literals like 0.0, 1.0 match the current Real precision
// Usage: Real x = REAL(3.14159265358979323846);
#ifndef REAL
#define REAL(x) static_cast<Real>(x)
#endif

// Complex must have the same underlying type as Real
// Changing Real automatically changes Complex precision
typedef std::complex<Real> Complex; // default complex type

namespace MML {

	// Pull std math overloads into the MML namespace so unqualified calls
	// (sin, cos, sqrt, …) resolve to the correct overload for Real type.
	// Without this, on GCC/Linux, unqualified sin(long double) silently
	// calls the C sin(double), losing precision.
	using std::abs;
	using std::acos;
	using std::asin;
	using std::atan;
	using std::atan2;
	using std::cos;
	using std::exp;
	using std::fabs;
	using std::log;
	using std::log10;
	using std::pow;
	using std::sin;
	using std::sqrt;
	using std::tan;
	using std::hypot;

	// Generic absolute value functions
	template<class Type>
	static Real Abs(const Type& a) {
		return std::abs(a);
	}
	template<class Type>
	static Real Abs(const std::complex<Type>& a) {
		return hypot(a.real(), a.imag());
	}

	////////////         Floating-Point Comparison Functions         ////////////
	/// @brief Check if two values are equal within absolute tolerance.
	/// @details Use for values expected to be near zero.
	inline bool isWithinAbsPrec(Real a, Real b, Real eps) { return std::abs(a - b) < eps; }

	/// @brief Check if two values are equal within relative tolerance.
	/// @details Use for values of similar magnitude away from zero.
	inline bool isWithinRelPrec(Real a, Real b, Real eps) { return std::abs(a - b) < eps * std::max(Abs(a), Abs(b)); }

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
	inline bool isNearlyEqual(Real a, Real b, Real eps) { return isNearlyEqual(a, b, eps, eps); }

	/// @brief Check if a value is nearly zero (within tolerance of zero).
	/// @details Use for checking if a value is effectively zero in numerical computations.
	/// @param a Value to check
	/// @param eps Tolerance (defaults to NumericalZeroThreshold)
	/// @return true if |a| < eps
	inline bool isNearlyZero(Real a, Real eps = PrecisionValues<Real>::NumericalZeroThreshold) { return std::abs(a) < eps; }

	template<class T>
	inline T POW2(const T& a) {
		const T& t = a;
		return t * t;
	}
	template<class T>
	inline T POW3(const T& a) {
		const T& t = a;
		return t * t * t;
	}
	template<class T>
	inline T POW4(const T& a) {
		const T& t = a * a;
		return t * t;
	}
	template<class T>
	inline T POW5(const T& a) {
		return POW4(a) * a;
	}

	////////////                  Constants                ////////////////
	namespace Constants {
		// Mathematical constants — stored as Real for zero-friction use in expressions.
		// Long-double literals preserve maximum precision during compile-time narrowing.
		// When Real=double: lossless compile-time narrowing from long double literal
		// When Real=long double: full extended precision preserved
		static inline constexpr Real PI = Real(3.141592653589793238462643383279502884L);				 // pi
		static inline constexpr Real INV_PI = Real(0.318309886183790671537767526745028724L);		 // 1/pi
		static inline constexpr Real INV_SQRTPI = Real(0.564189583547756286948079451560772586L); // 1/sqrt(pi)

		static inline constexpr Real E = Real(2.718281828459045235360287471352662498L);		 // e
		static inline constexpr Real LN2 = Real(0.693147180559945309417232121458176568L);	 // ln(2)
		static inline constexpr Real LN10 = Real(2.302585092994045684017991454684364208L); // ln(10)

		static inline constexpr Real SQRT2 = Real(1.414213562373095048801688724209698079L); // sqrt(2)
		static inline constexpr Real SQRT3 = Real(1.732050807568877293527446341505872367L); // sqrt(3)

		static inline constexpr Real GoldenRatio = Real(1.618033988749894848204586834365638118L); // (1 + sqrt(5)) / 2

		// Geometry epsilon for floating-point comparisons in geometric algorithms
		static inline constexpr Real GEOMETRY_EPSILON = Real(1e-10L);

		// Precision constants - use Real type for consistency with library's floating-point type
		static inline constexpr Real Eps = std::numeric_limits<Real>::epsilon();
		static inline constexpr Real PosInf = std::numeric_limits<Real>::infinity();
		static inline constexpr Real NegInf = -std::numeric_limits<Real>::infinity();
	} // namespace Constants

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
		int brentMaxSteps = 100; // Brent typically needs more steps but converges reliably

		// ODE solver parameters
		int odeSolverMaxSteps = 100000;

		static AlgorithmContext& Get() {
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

		static PrintContext& Get() {
			thread_local static PrintContext ctx;
			return ctx;
		}
	};

	// Backward compatible Defaults namespace (thread-safe via thread_local contexts)
	namespace Defaults {
		namespace Detail {
			// Proxy resolving to the CALLING thread's context on every access.
			// (A plain `static inline T&` alias binds ONE thread's thread_local instance
			// at static initialization, silently breaking per-thread semantics.)
			template<typename T, typename Context, T Context::* Member>
			struct ThreadLocalAlias {
				operator T() const { return Context::Get().*Member; }
				ThreadLocalAlias& operator=(T value) {
					Context::Get().*Member = value;
					return *this;
				}
			};
		}

		// Output defaults (per-thread; reads/writes affect the calling thread only)
		// Usage: Defaults::VectorPrintWidth = 20;
		static inline Detail::ThreadLocalAlias<int, PrintContext, &PrintContext::vectorWidth> VectorPrintWidth;
		static inline Detail::ThreadLocalAlias<int, PrintContext, &PrintContext::vectorPrecision> VectorPrintPrecision;
		static inline Detail::ThreadLocalAlias<int, PrintContext, &PrintContext::vectorNWidth> VectorNPrintWidth;
		static inline Detail::ThreadLocalAlias<int, PrintContext, &PrintContext::vectorNPrecision> VectorNPrintPrecision;

		//////////               Default precisions             ///////////
		// Use the precision values based on the Real type
		static inline const Real ComplexAreEqualTolerance = PrecisionValues<Real>::ComplexAreEqualTolerance;
		static inline const Real ComplexAreEqualAbsTolerance = PrecisionValues<Real>::ComplexAreEqualAbsTolerance;
		static inline const Real VectorIsEqualTolerance = PrecisionValues<Real>::VectorIsEqualTolerance;
		static inline const Real MatrixIsEqualTolerance = PrecisionValues<Real>::MatrixIsEqualTolerance;

		static inline const Real Pnt2CartIsEqualTolerance = PrecisionValues<Real>::Pnt2CartIsEqualTolerance;
		static inline const Real Pnt2PolarIsEqualTolerance = PrecisionValues<Real>::Pnt2PolarIsEqualTolerance;
		static inline const Real Pnt3CartIsEqualTolerance = PrecisionValues<Real>::Pnt3CartIsEqualTolerance;
		static inline const Real Pnt3SphIsEqualTolerance = PrecisionValues<Real>::Pnt3SphIsEqualTolerance;
		static inline const Real Pnt3CylIsEqualTolerance = PrecisionValues<Real>::Pnt3CylIsEqualTolerance;

		static inline const Real Vec2CartIsEqualTolerance = PrecisionValues<Real>::Vec2CartIsEqualTolerance;
		static inline const Real Vec3CartIsEqualTolerance = PrecisionValues<Real>::Vec3CartIsEqualTolerance;
		static inline const Real Vec3CartIsParallelTolerance = PrecisionValues<Real>::Vec3CartIsParallelTolerance;

		static inline const Real Vec3SphIsEqualTolerance = PrecisionValues<Real>::Vec3SphIsEqualTolerance;

		// Angle comparison tolerance (for wrap-aware angle equality)
		static inline const Real AngleIsEqualTolerance = PrecisionValues<Real>::AngleIsEqualTolerance;

		// Shape property tolerance (for geometric shape classification)
		static inline const Real ShapePropertyTolerance = PrecisionValues<Real>::ShapePropertyTolerance;

		static inline const Real Line3DAreEqualTolerance = PrecisionValues<Real>::Line3DAreEqualTolerance;
		static inline const Real Line3DIsPointOnLineTolerance = PrecisionValues<Real>::Line3DIsPointOnLineTolerance;
		static inline const Real Line3DIsPerpendicularTolerance = PrecisionValues<Real>::Line3DIsPerpendicularTolerance;
		static inline const Real Line3DIsParallelTolerance = PrecisionValues<Real>::Line3DIsParallelTolerance;
		static inline const Real Line3DIntersectionTolerance = PrecisionValues<Real>::Line3DIntersectionTolerance;

		static inline const Real Plane3DIsPointOnPlaneTolerance = PrecisionValues<Real>::Plane3DIsPointOnPlaneTolerance;

		static inline const Real Triangle3DIsPointInsideTolerance = PrecisionValues<Real>::Triangle3DIsPointInsideTolerance;
		static inline const Real Triangle3DIsRightTolerance = PrecisionValues<Real>::Triangle3DIsRightTolerance;
		static inline const Real Triangle3DIsIsoscelesTolerance = PrecisionValues<Real>::Triangle3DIsIsoscelesTolerance;
		static inline const Real Triangle3DIsEquilateralTolerance = PrecisionValues<Real>::Triangle3DIsEquilateralTolerance;

		static inline const Real IsMatrixSymmetricTolerance = PrecisionValues<Real>::IsMatrixSymmetricTolerance;
		static inline const Real IsMatrixDiagonalTolerance = PrecisionValues<Real>::IsMatrixDiagonalTolerance;
		static inline const Real IsMatrixUnitTolerance = PrecisionValues<Real>::IsMatrixUnitTolerance;
		static inline const Real IsMatrixZeroTolerance = PrecisionValues<Real>::IsMatrixZeroTolerance;
		static inline const Real IsMatrixOrthogonalTolerance = PrecisionValues<Real>::IsMatrixOrthogonalTolerance;

		static inline const Real RankAlgEPS = PrecisionValues<Real>::RankAlgEPS;

		// Numerical thresholds from PrecisionValues
		static inline const Real DefaultTolerance = PrecisionValues<Real>::DefaultTolerance;
		static inline const Real OrthogonalityTolerance = PrecisionValues<Real>::OrthogonalityTolerance;

		// Algorithm parameters (per-thread; reads/writes affect the calling thread only)
		// Usage: Defaults::TrapezoidIntegrationEPS = 1e-6;
		static inline Detail::ThreadLocalAlias<Real, AlgorithmContext, &AlgorithmContext::trapezoidIntegrationEPS> TrapezoidIntegrationEPS;
		static inline Detail::ThreadLocalAlias<Real, AlgorithmContext, &AlgorithmContext::simpsonIntegrationEPS> SimpsonIntegrationEPS;
		static inline Detail::ThreadLocalAlias<Real, AlgorithmContext, &AlgorithmContext::rombergIntegrationEPS> RombergIntegrationEPS;
		static inline Detail::ThreadLocalAlias<Real, AlgorithmContext, &AlgorithmContext::workIntegralPrecision> WorkIntegralPrecision;
		static inline Detail::ThreadLocalAlias<Real, AlgorithmContext, &AlgorithmContext::lineIntegralPrecision> LineIntegralPrecision;

		static inline Detail::ThreadLocalAlias<int, AlgorithmContext, &AlgorithmContext::bisectionMaxSteps> BisectionMaxSteps;
		static inline Detail::ThreadLocalAlias<int, AlgorithmContext, &AlgorithmContext::newtonRaphsonMaxSteps> NewtonRaphsonMaxSteps;
		static inline Detail::ThreadLocalAlias<int, AlgorithmContext, &AlgorithmContext::brentMaxSteps> BrentMaxSteps;

		static inline Detail::ThreadLocalAlias<int, AlgorithmContext, &AlgorithmContext::trapezoidIntegrationMaxSteps> TrapezoidIntegrationMaxSteps;
		static inline Detail::ThreadLocalAlias<int, AlgorithmContext, &AlgorithmContext::simpsonIntegrationMaxSteps> SimpsonIntegrationMaxSteps;
		static inline Detail::ThreadLocalAlias<int, AlgorithmContext, &AlgorithmContext::rombergIntegrationMaxSteps> RombergIntegrationMaxSteps;
		static inline Detail::ThreadLocalAlias<int, AlgorithmContext, &AlgorithmContext::rombergIntegrationUsedPnts> RombergIntegrationUsedPnts;

		static inline Detail::ThreadLocalAlias<int, AlgorithmContext, &AlgorithmContext::odeSolverMaxSteps> ODESolverMaxSteps;
	} // namespace Defaults
} // namespace MML
#endif

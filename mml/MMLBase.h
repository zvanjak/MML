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

// MML headers first (catches missing includes in them)
#include <mml/MMLTypeDefs.h>
#include <mml/MMLExceptions.h>
#include <mml/MMLPrecision.h>
#include <mml/MMLConcepts.h>

// Standard headers - only what MMLBase.h actually uses
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <limits>
#include <type_traits>

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
		// Unsigned magnitudes are already nonnegative; std::abs has no unsigned overload.
		if constexpr (std::is_unsigned_v<Type>)
			return static_cast<Real>(a);
		else
			return std::abs(a);
	}
	template<class Type>
	static Real Abs(const std::complex<Type>& a) {
		return hypot(a.real(), a.imag());
	}

	/// @brief Overflow/underflow-safe Euclidean norm accumulator (LAPACK dnrm2 scale/ssq scheme).
	/// @details sqrt(sum x_i^2) overflows for |x_i| > ~1e154 and underflows to 0 for |x_i| < ~1e-154
	///          although the norm itself is representable. Feed |x_i| (a real magnitude) one at a
	///          time via add(), read the norm with value(). Used by Vector/VectorN::NormL2
	///          (analysis v4, roadmap item 18: "naive norms").
	struct StableNormL2Accumulator {
		Real scale = 0.0;   ///< Largest magnitude seen so far
		Real ssq = 1.0;     ///< Sum of (x_i / scale)^2

		void add(Real magnitude) {
			if (magnitude == 0.0)
				return;
			if (scale < magnitude) {
				const Real ratio = scale / magnitude;
				ssq = 1.0 + ssq * ratio * ratio;
				scale = magnitude;
			} else {
				const Real ratio = magnitude / scale;
				ssq += ratio * ratio;
			}
		}
		Real value() const { return scale * std::sqrt(ssq); }
	};

	////////////         Floating-Point Comparison Functions         ////////////
	// The canonical comparison vocabulary (docs/TOLERANCE_POLICY.md section 3). Inline tolerance
	// arithmetic elsewhere in the library should go through these so the semantics are in one place.

	/// @brief |a - b| < absEps. Use when the values are expected to be near zero (or the scale is known).
	inline bool IsWithinAbs(Real a, Real b, Real absEps) { return std::abs(a - b) < absEps; }

	/// @brief |a - b| < relEps * max(|a|, |b|), with a == b always true (so (0, 0) is "equal").
	/// Use when both values are away from zero and only their magnitude matters.
	inline bool IsWithinRel(Real a, Real b, Real relEps) {
		return a == b || std::abs(a - b) < relEps * std::max(Abs(a), Abs(b));
	}

	/// @brief Combined test: |a - b| < absEps + relEps * max(|a|, |b|). The recommended general
	/// floating-point equality: absEps dominates near zero, relEps dominates for large values.
	inline bool IsNearlyEqual(Real a, Real b, Real absEps, Real relEps) {
		return std::abs(a - b) < absEps + relEps * std::max(Abs(a), Abs(b));
	}

	/// @brief Combined test with one tolerance used for both parts.
	inline bool IsNearlyEqual(Real a, Real b, Real eps = Precision::EqualityTolerance) { return IsNearlyEqual(a, b, eps, eps); }

	/// @brief |a| < eps. "This computed quantity is negligible." Exact `== 0` is reserved for
	/// structural/literal inputs and IEEE special cases (TOLERANCE_POLICY section 4).
	inline bool IsNearlyZero(Real a, Real eps = Precision::ZeroTolerance) { return std::abs(a) < eps; }

	/// @brief |x| <= eps * scale: a scale-aware negligibility test, e.g. a pivot against the matrix norm
	/// (`IsNegligibleRelativeTo(pivot, norm * n, Precision::Eps)`), a residual against ‖b‖, a step against |x|.
	inline bool IsNegligibleRelativeTo(Real x, Real scale, Real eps = Precision::Eps) { return std::abs(x) <= eps * std::abs(scale); }

	/// @brief Guard for division: true when |denominator| is so small that the quotient is meaningless.
	inline bool IsUnsafeDenominator(Real d) { return std::abs(d) < Precision::DivisionSafetyThreshold; }

	// ---- 2.0 spellings (CONTRACT_POLICY section 11: free predicates are PascalCase) ----
	[[deprecated("use IsWithinAbs()")]]
	inline bool isWithinAbsPrec(Real a, Real b, Real eps) { return IsWithinAbs(a, b, eps); }
	[[deprecated("use IsWithinRel()")]]
	inline bool isWithinRelPrec(Real a, Real b, Real eps) { return IsWithinRel(a, b, eps); }
	[[deprecated("use IsNearlyEqual()")]]
	inline bool isNearlyEqual(Real a, Real b, Real absEps, Real relEps) { return IsNearlyEqual(a, b, absEps, relEps); }
	[[deprecated("use IsNearlyEqual()")]]
	inline bool isNearlyEqual(Real a, Real b, Real eps) { return IsNearlyEqual(a, b, eps, eps); }
	[[deprecated("use IsNearlyZero()")]]
	inline bool isNearlyZero(Real a, Real eps = Precision::ZeroTolerance) { return IsNearlyZero(a, eps); }

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
		// Integration parameters — runtime-adjustable per thread; the values come from Precision
		Real trapezoidIntegrationEPS = Precision::TrapezoidIntegrationTolerance;
		Real simpsonIntegrationEPS = Precision::SimpsonIntegrationTolerance;
		Real rombergIntegrationEPS = Precision::RombergIntegrationTolerance;
		Real workIntegralPrecision = Precision::LineIntegralTolerance;
		Real lineIntegralPrecision = Precision::LineIntegralTolerance;

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

		//////////   Deprecated tolerance mirrors (TOLERANCE_POLICY section 7)   ///////////
		// Precision::X is the canonical spelling; these copies are removed after one minor release.
		[[deprecated("use Precision::ComplexAreEqualTolerance")]] static inline const Real ComplexAreEqualTolerance = Precision::ComplexAreEqualTolerance;
		[[deprecated("use Precision::ComplexAreEqualAbsTolerance")]] static inline const Real ComplexAreEqualAbsTolerance = Precision::ComplexAreEqualAbsTolerance;
		[[deprecated("use Precision::VectorIsEqualTolerance")]] static inline const Real VectorIsEqualTolerance = Precision::VectorIsEqualTolerance;
		[[deprecated("use Precision::MatrixIsEqualTolerance")]] static inline const Real MatrixIsEqualTolerance = Precision::MatrixIsEqualTolerance;

		[[deprecated("use Precision::Pnt2CartIsEqualTolerance")]] static inline const Real Pnt2CartIsEqualTolerance = Precision::Pnt2CartIsEqualTolerance;
		[[deprecated("use Precision::Pnt2PolarIsEqualTolerance")]] static inline const Real Pnt2PolarIsEqualTolerance = Precision::Pnt2PolarIsEqualTolerance;
		[[deprecated("use Precision::Pnt3CartIsEqualTolerance")]] static inline const Real Pnt3CartIsEqualTolerance = Precision::Pnt3CartIsEqualTolerance;
		[[deprecated("use Precision::Pnt3SphIsEqualTolerance")]] static inline const Real Pnt3SphIsEqualTolerance = Precision::Pnt3SphIsEqualTolerance;
		[[deprecated("use Precision::Pnt3CylIsEqualTolerance")]] static inline const Real Pnt3CylIsEqualTolerance = Precision::Pnt3CylIsEqualTolerance;

		[[deprecated("use Precision::Vec2CartIsEqualTolerance")]] static inline const Real Vec2CartIsEqualTolerance = Precision::Vec2CartIsEqualTolerance;
		[[deprecated("use Precision::Vec3CartIsEqualTolerance")]] static inline const Real Vec3CartIsEqualTolerance = Precision::Vec3CartIsEqualTolerance;
		[[deprecated("use Precision::Vec3CartIsParallelTolerance")]] static inline const Real Vec3CartIsParallelTolerance = Precision::Vec3CartIsParallelTolerance;

		[[deprecated("use Precision::Vec3SphIsEqualTolerance")]] static inline const Real Vec3SphIsEqualTolerance = Precision::Vec3SphIsEqualTolerance;

		// Angle comparison tolerance (for wrap-aware angle equality)
		[[deprecated("use Precision::AngleIsEqualTolerance")]] static inline const Real AngleIsEqualTolerance = Precision::AngleIsEqualTolerance;

		// Shape property tolerance (for geometric shape classification)
		[[deprecated("use Precision::ShapePropertyTolerance")]] static inline const Real ShapePropertyTolerance = Precision::ShapePropertyTolerance;

		[[deprecated("use Precision::Line3DAreEqualTolerance")]] static inline const Real Line3DAreEqualTolerance = Precision::Line3DAreEqualTolerance;
		[[deprecated("use Precision::Line3DIsPointOnLineTolerance")]] static inline const Real Line3DIsPointOnLineTolerance = Precision::Line3DIsPointOnLineTolerance;
		[[deprecated("use Precision::Line3DIsPerpendicularTolerance")]] static inline const Real Line3DIsPerpendicularTolerance = Precision::Line3DIsPerpendicularTolerance;
		[[deprecated("use Precision::Line3DIsParallelTolerance")]] static inline const Real Line3DIsParallelTolerance = Precision::Line3DIsParallelTolerance;
		[[deprecated("use Precision::Line3DIntersectionTolerance")]] static inline const Real Line3DIntersectionTolerance = Precision::Line3DIntersectionTolerance;

		[[deprecated("use Precision::Plane3DIsPointOnPlaneTolerance")]] static inline const Real Plane3DIsPointOnPlaneTolerance = Precision::Plane3DIsPointOnPlaneTolerance;

		[[deprecated("use Precision::Triangle3DIsPointInsideTolerance")]] static inline const Real Triangle3DIsPointInsideTolerance = Precision::Triangle3DIsPointInsideTolerance;
		[[deprecated("use Precision::Triangle3DIsRightTolerance")]] static inline const Real Triangle3DIsRightTolerance = Precision::Triangle3DIsRightTolerance;
		[[deprecated("use Precision::Triangle3DIsIsoscelesTolerance")]] static inline const Real Triangle3DIsIsoscelesTolerance = Precision::Triangle3DIsIsoscelesTolerance;
		[[deprecated("use Precision::Triangle3DIsEquilateralTolerance")]] static inline const Real Triangle3DIsEquilateralTolerance = Precision::Triangle3DIsEquilateralTolerance;

		[[deprecated("use Precision::IsMatrixSymmetricTolerance")]] static inline const Real IsMatrixSymmetricTolerance = Precision::IsMatrixSymmetricTolerance;
		[[deprecated("use Precision::IsMatrixDiagonalTolerance")]] static inline const Real IsMatrixDiagonalTolerance = Precision::IsMatrixDiagonalTolerance;
		[[deprecated("use Precision::IsMatrixUnitTolerance")]] static inline const Real IsMatrixUnitTolerance = Precision::IsMatrixUnitTolerance;
		[[deprecated("use Precision::IsMatrixZeroTolerance")]] static inline const Real IsMatrixZeroTolerance = Precision::IsMatrixZeroTolerance;
		[[deprecated("use Precision::IsMatrixOrthogonalTolerance")]] static inline const Real IsMatrixOrthogonalTolerance = Precision::IsMatrixOrthogonalTolerance;

		[[deprecated("use Precision::RankAlgEPS")]] static inline const Real RankAlgEPS = Precision::RankAlgEPS;

		// Numerical thresholds from PrecisionValues
		[[deprecated("use Precision::DefaultTolerance")]] static inline const Real DefaultTolerance = Precision::DefaultTolerance;
		[[deprecated("use Precision::OrthogonalityTolerance")]] static inline const Real OrthogonalityTolerance = Precision::OrthogonalityTolerance;

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

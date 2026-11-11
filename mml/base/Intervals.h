///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Intervals.h                                                         ///
///  Description: Real and complex interval classes for interval arithmetic           ///
///               Bounds checking and range operations                                ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
/// @file Intervals.h
/// @brief Real number interval classes for interval arithmetic and set operations.
/// @ingroup Base
/// @section intervals_overview Overview
/// This file provides a comprehensive hierarchy of real number interval classes
/// supporting various endpoint combinations (open, closed, infinite) and set operations.
/// The classes are useful for:
/// - Domain specification for functions
/// - Numerical integration bounds
/// - Generating equidistant point distributions for sampling
/// - Interval arithmetic and containment checking
/// @section intervals_hierarchy Class Hierarchy
/// @code
/// BaseInterval (abstract base)
/// ├── CompleteRInterval              (-∞, +∞)
/// ├── CompleteRWithReccuringPointHoles  ℝ \ {periodic points}
/// ├── OpenInterval                    (a, b)
/// ├── ClosedInterval                  [a, b]
/// ├── OpenClosedInterval              (a, b]
/// ├── ClosedOpenInterval              [a, b)
/// ├── ClosedIntervalWithReccuringPointHoles  [a,b] \ {periodic points}
/// ├── NegInfToOpenInterval            (-∞, a)
/// ├── NegInfToClosedInterval          (-∞, a]
/// ├── OpenToInfInterval               (a, +∞)
/// └── ClosedToInfInterval             [a, +∞)
/// CompositeInterval (union of disjoint BaseInterval pieces; alias: Interval)
/// - Closed set algebra via operators: | (union), & (intersection), - (difference), ~ (complement)
/// - Static helpers: Union(), Intersection(), Difference(), Complement()
/// @endcode
/// @section intervals_notation Mathematical Notation
/// | Notation | Meaning | Class |
/// |----------|---------|-------|
/// | (a, b)   | Open interval | OpenInterval |
/// | [a, b]   | Closed interval | ClosedInterval |
/// | (a, b]   | Half-open (left open) | OpenClosedInterval |
/// | [a, b)   | Half-open (right open) | ClosedOpenInterval |
/// | (-∞, a)  | Left-infinite open | NegInfToOpenInterval |
/// | (-∞, a]  | Left-infinite closed | NegInfToClosedInterval |
/// | (a, +∞)  | Right-infinite open | OpenToInfInterval |
/// | [a, +∞)  | Right-infinite closed | ClosedToInfInterval |
/// @section intervals_usage Usage Examples
/// @code{.cpp}
/// // Create basic intervals
/// ClosedInterval unitInterval(0.0, 1.0);      // [0, 1]
/// OpenInterval openUnit(0.0, 1.0);            // (0, 1)
/// ClosedToInfInterval positive(0.0);          // [0, +∞)
/// // Check containment
/// bool hasPoint = unitInterval.contains(0.5); // true
/// bool hasEndpt = openUnit.contains(0.0);     // false (open at 0)
/// // Generate sampling points
/// std::vector<Real> points;
/// unitInterval.getEquidistantCovering(10, points);
/// // points = {0.0, 0.111..., 0.222..., ..., 1.0}
/// // Set operations with composite intervals
/// Interval result = Interval::Intersection(
/// ClosedInterval(-1.0, 1.0),
/// ClosedInterval(0.0, 2.0)
/// );  // [0, 1]
/// @endcode
/// @see IntegrationMethod Classes that use intervals for integration bounds
/// @see IRealFunction Functions that may have domain restrictions

#if !defined MML_INTERVALS_H
#define MML_INTERVALS_H

#include <mml/MMLBase.h>
#include <mml/MMLNumericValidation.h>

#include <mml/interfaces/IInterval.h>

// Standard headers - include what we use
#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <limits>
#include <memory>
#include <numeric>
#include <utility>
#include <vector>

namespace MML {
	/// /** @name Interval Type Definitions
	/// @{ */


	/// @brief Enumeration of interval endpoint types.
	/// Specifies whether an interval endpoint is open (excluded), closed (included),
	/// or extends to positive/negative infinity.

	enum class EndpointType {
		OPEN,	 ///< Endpoint is excluded from the interval
		CLOSED,	 ///< Endpoint is included in the interval
		NEG_INF, ///< Endpoint extends to negative infinity (-∞)
		POS_INF	 ///< Endpoint extends to positive infinity (+∞)
	};

	namespace IntervalDetail {
		struct SamplingRange { Real lower, upper, center; };

		inline void ValidateSamplingWindow(Real lower, Real upper) {
			ValidateFinite(lower, "Interval sampling window lower bound");
			ValidateFinite(upper, "Interval sampling window upper bound");
			if (lower > upper) throw ArgumentError("Interval sampling window is reversed");
		}

		inline bool ClipSamplingRange(Real lower, EndpointType lowerType, Real upper, EndpointType upperType,
			Real windowLower, Real windowUpper, SamplingRange& range) {
			const Real lo = std::max(lower, windowLower), hi = std::min(upper, windowUpper);
			if (lo > hi) return false;
			const bool openLower = lo == lower && lowerType == EndpointType::OPEN;
			const bool openUpper = hi == upper && upperType == EndpointType::OPEN;
			if (lo == hi && (openLower || openUpper)) return false;
			range.center = std::midpoint(lo, hi);
			range.lower = openLower ? std::nextafter(lo, hi) : lo;
			range.upper = openUpper ? std::nextafter(hi, lo) : hi;
			return range.lower <= range.upper;
		}

		inline bool FindSamplingMember(const IInterval& domain, Real lower, Real upper, Real target, Real& value) {
			if (domain.contains(target)) { value = target; return true; }
			if (lower == upper) return false;
			Real left = lower, right = upper;
			// Dyadic alternatives avoid periodic holes without iterating through their roundoff bands one ULP at a time.
			for (int bit = 0; bit < std::numeric_limits<Real>::digits; ++bit) {
				left = std::midpoint(left, target);
				right = std::midpoint(target, right);
				if (domain.contains(left)) { value = left; return true; }
				if (domain.contains(right)) { value = right; return true; }
			}
			if (domain.contains(lower)) { value = lower; return true; }
			if (domain.contains(upper)) { value = upper; return true; }
			if (std::nextafter(lower, upper) == upper) return false;
			throw NumericalMethodError("Interval sampling: no member resolved in a grid cell");
		}

		inline std::vector<Real> SampleRange(const IInterval& domain, int count, const SamplingRange& range) {
			std::vector<Real> result;
			if (count <= 0) return result;
			if (range.lower == range.upper) {
				if (domain.contains(range.lower)) result.push_back(range.lower);
				return result;
			}
			if (count > 1 && std::nextafter(range.lower, range.upper) == range.upper) {
				if (domain.contains(range.lower)) result.push_back(range.lower);
				if (domain.contains(range.upper)) result.push_back(range.upper);
				return result;
			}
			const auto grid = [&](int i) {
				if (count == 1) return std::clamp(range.center, range.lower, range.upper);
				if (i == 0) return range.lower;
				if (i == count - 1) return range.upper;
				return std::lerp(range.lower, range.upper, Real(i) / Real(count - 1));
			};
			for (int i = 0; i < count; ++i) {
				const Real target = grid(i);
				const Real left = i == 0 ? range.lower : std::midpoint(grid(i - 1), target);
				const Real right = i == count - 1 ? range.upper : std::midpoint(target, grid(i + 1));
				Real value;
				if (FindSamplingMember(domain, left, right, target, value)) result.push_back(value);
			}
			std::sort(result.begin(), result.end());
			result.erase(std::unique(result.begin(), result.end()), result.end());
			return result;
		}
	}
	/// /** @} */


	class CompositeInterval; // Forward declaration

	/// /** @name Base Interval Class
	/// @{ */


	/// @brief Abstract base class for all interval types.
	/// Provides common interface and storage for interval bounds and endpoint types.
	/// Derived classes implement specific endpoint semantics (open, closed, infinite).
	/// Key features:
	/// - Lower and upper bounds with configurable endpoint types
	/// - Point containment testing via contains()
	/// - Equidistant point generation for numerical sampling
	/// - Length and bound queries
	/// @note The Interval class (composite) is declared as friend to access
	/// protected members for set operations.

	class BaseInterval : public IInterval {
		friend class CompositeInterval; // Allow CompositeInterval to access protected members

	protected:
		Real _lower, _upper;				 ///< Lower and upper bounds
		EndpointType _lowerType, _upperType; ///< Endpoint types (open/closed/infinite)

		BaseInterval(Real lower, EndpointType lowerType, Real upper, EndpointType upperType)
			: _lower(lower)
			, _upper(upper)
			, _lowerType(lowerType)
			, _upperType(upperType) {
			ValidateFinite(lower, "Interval lower bound");
			ValidateFinite(upper, "Interval upper bound");
			if (lower > upper)
				throw ArgumentError("Interval: lower bound must not exceed upper bound");
		}

	public:
		~BaseInterval() override {}

		/// @brief Gets the lower bound of the interval.
		Real getLowerBound() const override { return _lower; }

		/// @brief Gets the upper bound of the interval.
		Real getUpperBound() const override { return _upper; }

		/// @brief Gets the length of the interval (upper - lower).
		Real getLength() const override { return _upper - _lower; }

		/// @brief Returns true if the interval is continuous (no holes).
		bool isContinuous() const override { return true; } // we suppose continuous intervals by default

		/// @brief Generates equidistant points covering the interval.
		/// Creates a uniform distribution of points across the interval, useful for:
		/// - Numerical integration sampling
		/// - Function plotting
		/// - Data point generation for curve fitting
		/// Open endpoints move inward by one representable step. Periodic holes
		/// use checked alternatives within grid cells, not a globally uniform grid.
		/// Unbounded sets are clipped to +/-Precision::UnboundedIntervalSamplingLimit.
		/// @param numPoints Maximum number of distinct points; nonpositive requests clear output.
		/// @param[out] points Vector to receive the generated points.
		/// @code{.cpp}
		/// ClosedInterval interval(0.0, 1.0);
		/// std::vector<Real> pts;
		/// interval.getEquidistantCovering(5, pts);
		/// // pts = {0.0, 0.25, 0.5, 0.75, 1.0}
		/// @endcode

		void getEquidistantCovering(int numPoints, std::vector<Real>& points) const override {
			const bool unbounded = _lowerType == EndpointType::NEG_INF || _upperType == EndpointType::POS_INF;
			getEquidistantCovering(numPoints, points,
				unbounded ? -Precision::UnboundedIntervalSamplingLimit : _lower,
				unbounded ? Precision::UnboundedIntervalSamplingLimit : _upper);
		}

		/// @brief Sample only members inside [windowLower,windowUpper]; output is unchanged on failure.
		void getEquidistantCovering(int numPoints, std::vector<Real>& points, Real windowLower, Real windowUpper) const override {
			if (numPoints <= 0) { points.clear(); return; }
			IntervalDetail::ValidateSamplingWindow(windowLower, windowUpper);
			IntervalDetail::SamplingRange range;
			std::vector<Real> result;
			if (IntervalDetail::ClipSamplingRange(_lower, _lowerType, _upper, _upperType, windowLower, windowUpper, range))
				result = IntervalDetail::SampleRange(*this, numPoints, range);
			points = std::move(result);
		}
	};
	/// /** @} */
	// End Base Interval Class

	/// /** @name Complete Real Line Intervals
	/// @brief Intervals spanning all of ℝ (with optional periodic holes)
	/// @{ */


	/// @brief The complete real line interval (-∞, +∞).
	/// Represents all real numbers. Every point is contained.
	/// Useful as a default domain for functions defined everywhere.

	class CompleteRInterval : public BaseInterval {
	public:
		/// @brief Constructs the interval (-∞, +∞).
		CompleteRInterval()
			: BaseInterval(-std::numeric_limits<Real>::max(), EndpointType::NEG_INF, std::numeric_limits<Real>::max(),
						   EndpointType::POS_INF) {}

		/// @brief Returns true for any real number x.
		bool contains(Real x) const { return std::isfinite(x); }
	};

	/// @brief The real line with periodic point exclusions.
	/// Represents ℝ \ {hole₀ + n·Δhole : n ∈ ℤ}, useful for domains of
	/// functions with periodic singularities such as tan(x) or csc(x).
	/// @code{.cpp}
	/// // Domain of tan(x): exclude x = π/2 + nπ
	/// CompleteRWithReccuringPointHoles tanDomain(Constants::PI/2, Constants::PI);
	/// tanDomain.contains(0.0);                    // true
	/// tanDomain.contains(Constants::PI / 2);     // false
	/// @endcode

	class CompleteRWithReccuringPointHoles : public BaseInterval {
		Real _hole0, _holeDelta; ///< First hole position and period
	public:
		/// @brief Constructs ℝ with periodic holes.
		/// @param hole0 Position of the first hole.
		/// @param holeDelta Period between consecutive holes.

		CompleteRWithReccuringPointHoles(Real hole0, Real holeDelta)
			: BaseInterval(-std::numeric_limits<Real>::max(), EndpointType::NEG_INF, std::numeric_limits<Real>::max(),
						   EndpointType::POS_INF) {
			_hole0 = hole0;
			_holeDelta = holeDelta;
			ValidateFinite(hole0, "Interval hole origin");
			ValidateFinite(holeDelta, "Interval hole period");
			if (holeDelta == 0) throw ArgumentError("Interval hole period must be nonzero");
		}

		/// @brief Returns true if x is not at a hole position.
		bool contains(Real x) const {
			if (!std::isfinite(x)) return false;
			Real diff = (x - _hole0) / _holeDelta;
			if (!std::isfinite(diff)) throw NumericalMethodError("Interval periodic membership is unrepresentable");
			if (IsNearlyZero(diff - std::round(diff), Precision::RoundoffZeroThreshold))
				return false;

			return true;
		}
		/// @brief Returns false since this interval has discontinuities.
		bool isContinuous() const { return false; }
	};
	/// /** @} */
	// End Complete Real Line Intervals

	/// /** @name Bounded Intervals
	/// @brief Intervals with finite bounds and various endpoint combinations
	/// @{ */


	/// @brief Open interval (a, b) - both endpoints excluded.
	/// Contains all x such that a < x < b.

	class OpenInterval : public BaseInterval {
	public:
		/// @brief Constructs the open interval (lower, upper).
		/// @param lower Left endpoint (excluded).
		/// @param upper Right endpoint (excluded).

		OpenInterval(Real lower, Real upper)
			: BaseInterval(lower, EndpointType::OPEN, upper, EndpointType::OPEN) {}

		/// @brief Returns true if lower < x < upper.
		bool contains(Real x) const { return (x > _lower) && (x < _upper); }
	};

	/// @brief Half-open interval (a, b] - left open, right closed.
	/// Contains all x such that a < x ≤ b.

	class OpenClosedInterval : public BaseInterval {
	public:
		/// @brief Constructs the interval (lower, upper].
		/// @param lower Left endpoint (excluded).
		/// @param upper Right endpoint (included).

		OpenClosedInterval(Real lower, Real upper)
			: BaseInterval(lower, EndpointType::OPEN, upper, EndpointType::CLOSED) {}

		/// @brief Returns true if lower < x ≤ upper.
		bool contains(Real x) const { return (x > _lower) && (x <= _upper); }
	};

	/// @brief Closed interval [a, b] - both endpoints included.
	/// Contains all x such that a ≤ x ≤ b. The most commonly used interval type.

	class ClosedInterval : public BaseInterval {
	public:
		/// @brief Constructs the closed interval [lower, upper].
		/// @param lower Left endpoint (included).
		/// @param upper Right endpoint (included).

		ClosedInterval(Real lower, Real upper)
			: BaseInterval(lower, EndpointType::CLOSED, upper, EndpointType::CLOSED) {}

		/// @brief Returns true if lower ≤ x ≤ upper.
		bool contains(Real x) const { return (x >= _lower) && (x <= _upper); }
	};

	/// @brief Closed interval with periodic point exclusions.
	/// Represents [a, b] \ {hole₀ + n·Δhole : n ∈ ℤ, hole₀ + n·Δhole ∈ [a,b]}.
	/// Useful for bounded domains of functions with internal singularities.

	class ClosedIntervalWithReccuringPointHoles : public BaseInterval {
		Real _hole0, _holeDelta; ///< First hole position and period
	public:
		/// @brief Constructs [lower, upper] with periodic holes.
		/// @param lower Left endpoint (included).
		/// @param upper Right endpoint (included).
		/// @param hole0 Position of the first hole.
		/// @param holeDelta Period between consecutive holes.

		ClosedIntervalWithReccuringPointHoles(Real lower, Real upper, Real hole0, Real holeDelta)
			: BaseInterval(lower, EndpointType::CLOSED, upper, EndpointType::CLOSED) {
			_hole0 = hole0;
			_holeDelta = holeDelta;
			ValidateFinite(hole0, "Interval hole origin");
			ValidateFinite(holeDelta, "Interval hole period");
			if (holeDelta == 0) throw ArgumentError("Interval hole period must be nonzero");
		}

		/// @brief Returns true if x is in [lower, upper] and not at a hole.
		bool contains(Real x) const {
			if (!std::isfinite(x) || x < _lower || x > _upper)
				return false;

			Real diff = (x - _hole0) / _holeDelta;
			if (!std::isfinite(diff)) throw NumericalMethodError("Interval periodic membership is unrepresentable");
			if (IsNearlyZero(diff - std::round(diff), Precision::RoundoffZeroThreshold))
				return false;

			return true;
		}
		/// @brief Returns false since this interval has discontinuities.
		bool isContinuous() const { return false; }
	};

	/// @brief Half-open interval [a, b) - left closed, right open.
	/// Contains all x such that a ≤ x < b.

	class ClosedOpenInterval : public BaseInterval {
	public:
		/// @brief Constructs the interval [lower, upper).
		/// @param lower Left endpoint (included).
		/// @param upper Right endpoint (excluded).

		ClosedOpenInterval(Real lower, Real upper)
			: BaseInterval(lower, EndpointType::CLOSED, upper, EndpointType::OPEN) {}

		/// @brief Returns true if lower ≤ x < upper.
		bool contains(Real x) const { return (x >= _lower) && (x < _upper); }
	};
	/// /** @} */
	// End Bounded Intervals

	/// /** @name Semi-Infinite Intervals
	/// @brief Intervals extending to positive or negative infinity
	/// @{ */


	/// @brief Left-infinite open interval (-∞, a).
	/// Contains all x such that x < a.

	class NegInfToOpenInterval : public BaseInterval {
	public:
		/// @brief Constructs the interval (-∞, upper).
		/// @param upper Right endpoint (excluded).

		NegInfToOpenInterval(Real upper)
			: BaseInterval(-std::numeric_limits<Real>::max(), EndpointType::NEG_INF, upper, EndpointType::OPEN) {}

		/// @brief Returns true if x < upper.
		bool contains(Real x) const { return x < _upper; }
	};

	/// @brief Left-infinite closed interval (-∞, a].
	/// Contains all x such that x ≤ a.

	class NegInfToClosedInterval : public BaseInterval {
	public:
		/// @brief Constructs the interval (-∞, upper].
		/// @param upper Right endpoint (included).

		NegInfToClosedInterval(Real upper)
			: BaseInterval(-std::numeric_limits<Real>::max(), EndpointType::NEG_INF, upper, EndpointType::CLOSED) {}

		/// @brief Returns true if x ≤ upper.
		bool contains(Real x) const { return x <= _upper; }
	};

	/// @brief Right-infinite open interval (a, +∞).
	/// Contains all x such that x > a.

	class OpenToInfInterval : public BaseInterval {
	public:
		/// @brief Constructs the interval (lower, +∞).
		/// @param lower Left endpoint (excluded).

		OpenToInfInterval(Real lower)
			: BaseInterval(lower, EndpointType::OPEN, std::numeric_limits<Real>::max(), EndpointType::POS_INF) {}

		/// @brief Returns true if x > lower.
		bool contains(Real x) const { return x > _lower; }
	};

	/// @brief Right-infinite closed interval [a, +∞).
	/// Contains all x such that x ≥ a.

	class ClosedToInfInterval : public BaseInterval {
	public:
		/// @brief Constructs the interval [lower, +∞).
		/// @param lower Left endpoint (included).

		ClosedToInfInterval(Real lower)
			: BaseInterval(lower, EndpointType::CLOSED, std::numeric_limits<Real>::max(), EndpointType::POS_INF) {}

		/// @brief Returns true if x ≥ lower.
		bool contains(Real x) const { return x >= _lower; }
	};
	/// /** @} */
	// End Semi-Infinite Intervals

	/// /** @name Composite Interval
	/// @brief Union of multiple intervals with set operations
	/// @{ */


	/// @brief Composite interval representing a union of disjoint subintervals.
	/// Stores a collection of BaseInterval objects and provides:
	/// - Union: AddInterval() to accumulate subintervals
	/// - Intersection: Static method for intersecting two intervals
	/// - Difference: Static method for set difference A \ B
	/// - Complement: Static method for ℝ \ A
	/// @section interval_set_ops Set Operations
	/// | Operation | Description | Result |
	/// |-----------|-------------|--------|
	/// | A ∩ B | Intersection | Points in both A and B |
	/// | A \ B | Difference | Points in A but not in B |
	/// | ℝ \ A | Complement | All points not in A |
	/// @code{.cpp}
	/// // Intersection of [0, 2] and [1, 3] = [1, 2]
	/// Interval result = Interval::Intersection(
	/// ClosedInterval(0.0, 2.0),
	/// ClosedInterval(1.0, 3.0)
	/// );
	/// // Difference [0, 3] \ [1, 2] = [0, 1) ∪ (2, 3]
	/// Interval diff = Interval::Difference(
	/// ClosedInterval(0.0, 3.0),
	/// ClosedInterval(1.0, 2.0)
	/// );
	/// @endcode

	class CompositeInterval : public IInterval {
		std::vector<std::shared_ptr<BaseInterval>> _intervals; ///< Disjoint component pieces

		void requirePlainAlgebra() const {
			for (const auto& piece : _intervals)
				if (!piece->isContinuous())
					throw NotImplementedError("CompositeInterval: set algebra for hole-bearing pieces is not supported");
		}

		/// @brief Endpoint ordering used when sorting pieces by starting edge.
		static int endpointRank(EndpointType t) {
			switch (t) {
			case EndpointType::NEG_INF: return 0;
			case EndpointType::CLOSED:  return 1;
			case EndpointType::OPEN:    return 2;
			case EndpointType::POS_INF: return 3;
			}
			return 4;
		}

		/// @brief Combines upper-endpoint types when two pieces share the same upper bound.
		static EndpointType mergeUpperType(EndpointType a, EndpointType b) {
			if (a == EndpointType::POS_INF || b == EndpointType::POS_INF) return EndpointType::POS_INF;
			if (a == EndpointType::CLOSED  || b == EndpointType::CLOSED)  return EndpointType::CLOSED;
			return EndpointType::OPEN;
		}

		/// @brief Appends a plain endpoint-typed piece, choosing the concrete leaf type.
		/// @details Empty pieces (lo > hi, or a degenerate point with an open end) are skipped.
		CompositeInterval& addPiece(Real lo, EndpointType loType, Real hi, EndpointType hiType) {
			if (lo > hi) return *this;
			if (lo == hi && (loType == EndpointType::OPEN || hiType == EndpointType::OPEN)) return *this;

			if (loType == EndpointType::NEG_INF) {
				if (hiType == EndpointType::POS_INF) _intervals.emplace_back(std::make_shared<CompleteRInterval>());
				else if (hiType == EndpointType::CLOSED) _intervals.emplace_back(std::make_shared<NegInfToClosedInterval>(hi));
				else _intervals.emplace_back(std::make_shared<NegInfToOpenInterval>(hi));
			}
			else if (hiType == EndpointType::POS_INF) {
				if (loType == EndpointType::CLOSED) _intervals.emplace_back(std::make_shared<ClosedToInfInterval>(lo));
				else _intervals.emplace_back(std::make_shared<OpenToInfInterval>(lo));
			}
			else if (loType == EndpointType::CLOSED && hiType == EndpointType::CLOSED) _intervals.emplace_back(std::make_shared<ClosedInterval>(lo, hi));
			else if (loType == EndpointType::OPEN   && hiType == EndpointType::CLOSED) _intervals.emplace_back(std::make_shared<OpenClosedInterval>(lo, hi));
			else if (loType == EndpointType::CLOSED && hiType == EndpointType::OPEN)   _intervals.emplace_back(std::make_shared<ClosedOpenInterval>(lo, hi));
			else _intervals.emplace_back(std::make_shared<OpenInterval>(lo, hi));
			return *this;
		}
	public:
		/// @brief Constructs an empty composite interval.
		CompositeInterval() {}

		/// @brief Wraps a single interval piece as a one-piece composite (implicit).
		/// @details Reconstructs a plain endpoint-typed piece from the source bounds, enabling
		///          any BaseInterval to be passed where a CompositeInterval is expected.
		/// @throws NotImplementedError for non-continuous pieces, which cannot be
		/// copied faithfully through this abstract view. Use AddInterval(concretePiece)
		/// to retain holes for membership and sampling.
		CompositeInterval(const BaseInterval& piece) {
			if (!piece.isContinuous())
				throw NotImplementedError("CompositeInterval: cannot reconstruct a hole-bearing piece from bounds");
			addPiece(piece._lower, piece._lowerType, piece._upper, piece._upperType);
		}

		/// @brief Adds a subinterval to the composite.
		/// @tparam IntervalType Type derived from BaseInterval.
		/// @param interval The interval to add.
		/// @return Reference to this for chaining.

		template<class IntervalType>
		CompositeInterval& AddInterval(const IntervalType& interval) {
			_intervals.emplace_back(std::make_shared<IntervalType>(interval));
			return *this;
		}

		/// @brief Canonicalizes the piece list: sorts pieces and merges overlapping or
		///        adjacent ones (respecting endpoint types), yielding disjoint pieces.
		/// @details No-op when any piece is non-continuous (periodic holes), to avoid
		///          collapsing hole information into a hull.
		void Normalize() {
			if (_intervals.empty())
				return;
			for (const auto& p : _intervals)
				if (!p->isContinuous())
					return;

			struct Seg { Real lo, hi; EndpointType loT, hiT; };
			std::vector<Seg> segs;
			segs.reserve(_intervals.size());
			for (const auto& p : _intervals) {
				if (p->_lower == p->_upper && (p->_lowerType == EndpointType::OPEN || p->_upperType == EndpointType::OPEN))
					continue;
				segs.push_back({ p->_lower, p->_upper, p->_lowerType, p->_upperType });
			}

			std::sort(segs.begin(), segs.end(), [](const Seg& x, const Seg& y) {
				if (x.lo != y.lo) return x.lo < y.lo;
				return endpointRank(x.loT) < endpointRank(y.loT);
			});

			std::vector<Seg> merged;
			for (const Seg& s : segs) {
				if (merged.empty()) { merged.push_back(s); continue; }
				Seg& cur = merged.back();
				bool overlapsOrAdjacent = (s.lo < cur.hi) ||
					(s.lo == cur.hi && (cur.hiT == EndpointType::CLOSED || s.loT == EndpointType::CLOSED));
				if (overlapsOrAdjacent) {
					if (s.hi > cur.hi) { cur.hi = s.hi; cur.hiT = s.hiT; }
					else if (s.hi == cur.hi) cur.hiT = mergeUpperType(cur.hiT, s.hiT);
				} else {
					merged.push_back(s);
				}
			}

			_intervals.clear();
			for (const Seg& s : merged)
				addPiece(s.lo, s.loT, s.hi, s.hiT);
		}

	private:
		/// @brief Flips CLOSED<->OPEN; infinite endpoint types are returned unchanged.
		static EndpointType flip(EndpointType t) {
			if (t == EndpointType::CLOSED) return EndpointType::OPEN;
			if (t == EndpointType::OPEN)   return EndpointType::CLOSED;
			return t;
		}

		/// @brief Intersects two plain pieces; returns false when the result is empty.
		static bool intersectPieces(Real lo1, EndpointType lo1T, Real hi1, EndpointType hi1T,
		                            Real lo2, EndpointType lo2T, Real hi2, EndpointType hi2T,
		                            Real& lo, EndpointType& loT, Real& hi, EndpointType& hiT) {
			lo = std::max(lo1, lo2);
			hi = std::min(hi1, hi2);
			if (lo > hi) return false;

			if (lo1T == EndpointType::NEG_INF && lo2T != EndpointType::NEG_INF) { lo = lo2; loT = lo2T; }
			else if (lo2T == EndpointType::NEG_INF && lo1T != EndpointType::NEG_INF) { lo = lo1; loT = lo1T; }
			else if (lo1 == lo2) {
				if (lo1T == EndpointType::NEG_INF && lo2T == EndpointType::NEG_INF) loT = EndpointType::NEG_INF;
				else loT = (lo1T == EndpointType::CLOSED && lo2T == EndpointType::CLOSED) ? EndpointType::CLOSED : EndpointType::OPEN;
			}
			else loT = (lo == lo1) ? lo1T : lo2T;

			if (hi1T == EndpointType::POS_INF && hi2T != EndpointType::POS_INF) { hi = hi2; hiT = hi2T; }
			else if (hi2T == EndpointType::POS_INF && hi1T != EndpointType::POS_INF) { hi = hi1; hiT = hi1T; }
			else if (hi1 == hi2) {
				if (hi1T == EndpointType::POS_INF && hi2T == EndpointType::POS_INF) hiT = EndpointType::POS_INF;
				else hiT = (hi1T == EndpointType::CLOSED && hi2T == EndpointType::CLOSED) ? EndpointType::CLOSED : EndpointType::OPEN;
			}
			else hiT = (hi == hi1) ? hi1T : hi2T;

			if (lo == hi && (loT == EndpointType::OPEN || hiT == EndpointType::OPEN)) return false;
			return true;
		}
	public:
		///////////////////////////////////////////////////////////////////////
		///                    Closed set algebra (member)                  ///
		///////////////////////////////////////////////////////////////////////

		/// @brief Union with another interval set: this ∪ other (normalized).
		CompositeInterval unionWith(const CompositeInterval& other) const {
			requirePlainAlgebra();
			other.requirePlainAlgebra();
			CompositeInterval ret;
			ret._intervals = _intervals;
			ret._intervals.insert(ret._intervals.end(), other._intervals.begin(), other._intervals.end());
			ret.Normalize();
			return ret;
		}

		/// @brief Intersection with another interval set: this ∩ other (normalized).
		CompositeInterval intersectWith(const CompositeInterval& other) const {
			requirePlainAlgebra();
			other.requirePlainAlgebra();
			CompositeInterval ret;
			for (const auto& p : _intervals)
				for (const auto& q : other._intervals) {
					Real lo, hi; EndpointType loT, hiT;
					if (intersectPieces(p->_lower, p->_lowerType, p->_upper, p->_upperType,
					                    q->_lower, q->_lowerType, q->_upper, q->_upperType,
					                    lo, loT, hi, hiT))
						ret.addPiece(lo, loT, hi, hiT);
				}
			ret.Normalize();
			return ret;
		}

		/// @brief Complement over the whole real line: ℝ \ this (normalized).
		CompositeInterval complement() const {
			requirePlainAlgebra();
			CompositeInterval norm = *this;
			norm.Normalize();

			CompositeInterval ret;
			const Real negInf = -std::numeric_limits<Real>::max();
			const Real posInf =  std::numeric_limits<Real>::max();

			if (norm._intervals.empty()) {
				ret.addPiece(negInf, EndpointType::NEG_INF, posInf, EndpointType::POS_INF);
				return ret;
			}

			struct Seg { Real lo, hi; EndpointType loT, hiT; };
			std::vector<Seg> segs;
			segs.reserve(norm._intervals.size());
			for (const auto& p : norm._intervals)
				segs.push_back({ p->_lower, p->_upper, p->_lowerType, p->_upperType });

			if (segs.front().loT != EndpointType::NEG_INF)
				ret.addPiece(negInf, EndpointType::NEG_INF, segs.front().lo, flip(segs.front().loT));

			for (size_t i = 0; i + 1 < segs.size(); ++i)
				ret.addPiece(segs[i].hi, flip(segs[i].hiT), segs[i + 1].lo, flip(segs[i + 1].loT));

			if (segs.back().hiT != EndpointType::POS_INF)
				ret.addPiece(segs.back().hi, flip(segs.back().hiT), posInf, EndpointType::POS_INF);

			return ret;
		}

		/// @brief Set difference: this \ other = this ∩ (ℝ \ other).
		CompositeInterval differenceWith(const CompositeInterval& other) const {
			requirePlainAlgebra();
			return intersectWith(other.complement());
		}

		CompositeInterval operator|(const CompositeInterval& other) const { return unionWith(other); }
		CompositeInterval operator&(const CompositeInterval& other) const { return intersectWith(other); }
		CompositeInterval operator-(const CompositeInterval& other) const { return differenceWith(other); }
		CompositeInterval operator~() const { return complement(); }

		/// @brief Closed-algebra static overloads accepting interval sets.
		/// @details A single BaseInterval converts implicitly, so these compose freely.
		///          The BaseInterval,BaseInterval overloads below win for two plain pieces.
		static CompositeInterval Union(const CompositeInterval& a, const CompositeInterval& b)        { return a.unionWith(b); }
		static CompositeInterval Intersection(const CompositeInterval& a, const CompositeInterval& b) { return a.intersectWith(b); }
		static CompositeInterval Difference(const CompositeInterval& a, const CompositeInterval& b)   { return a.differenceWith(b); }
		static CompositeInterval Complement(const CompositeInterval& a)                               { return a.complement(); }

		/// @brief Gets the minimum lower bound across all subintervals.
		Real getLowerBound() const override {
			if (_intervals.empty())
				return 0;
			Real minLower = _intervals[0]->getLowerBound();
			for (const auto& interval : _intervals)
				minLower = std::min(minLower, interval->getLowerBound());
			return minLower;
		}

		/// @brief Gets the maximum upper bound across all subintervals.
		Real getUpperBound() const override {
			if (_intervals.empty())
				return 0;
			Real maxUpper = _intervals[0]->getUpperBound();
			for (const auto& interval : _intervals)
				maxUpper = std::max(maxUpper, interval->getUpperBound());
			return maxUpper;
		}

		/// @brief Gets the hull length (span from min to max bound).
		/// @return getUpperBound() - getLowerBound()
		/// @see getMeasure() for total length of all subintervals
		Real getLength() const override {
			return getUpperBound() - getLowerBound();
		}

		/// @brief Gets the total measure (sum of all subinterval lengths).
		/// @return Sum of individual subinterval lengths
		Real getMeasure() const override {
			Real totalMeasure = 0;
			for (const auto& interval : _intervals)
				totalMeasure += interval->getLength();
			return totalMeasure;
		}

		/// @brief Returns false (composite intervals are not continuous by definition).
		bool isContinuous() const override { return false; }

		/// @brief Checks if x is contained in any subinterval.
		/// @param x The point to test.
		/// @return True if x is in at least one subinterval.

		bool contains(Real x) const override {
			// check for each interval if it contains x
			for (auto& interval : _intervals) {
				if (interval->contains(x))
					return true;
			}
			return false;
		}

		/// @brief Checks if another interval is fully contained.
		/// @param other The interval to test.
		/// @return True if other minus this set is empty (including joint coverage).
		/// @throws NotImplementedError for hole-bearing set predicates.

		bool contains(const BaseInterval& other) const {
			requirePlainAlgebra();
			return CompositeInterval(other).differenceWith(*this)._intervals.empty();
		}

		/// @brief Checks if this composite intersects with another interval.
		/// @param other The interval to test for intersection.
		/// @return True if any subinterval overlaps with other.

		bool intersects(const BaseInterval& other) const {
			return !intersectWith(CompositeInterval(other))._intervals.empty();
		}

		/// @brief Samples at most numPoints distinct ordered members.
		/// Budgets use normalized clipped occupied lengths, never hull gaps.
		/// If any component is unbounded the default window is
		/// +/-Precision::UnboundedIntervalSamplingLimit; otherwise the hull is used.
		/// @param numPoints Maximum number of points to generate.
		/// @param[out] points Vector to receive the generated points.

		void getEquidistantCovering(int numPoints, std::vector<Real>& points) const override {
			const bool unbounded = std::any_of(_intervals.begin(), _intervals.end(), [](const auto& piece) {
				return piece->_lowerType == EndpointType::NEG_INF || piece->_upperType == EndpointType::POS_INF;
			});
			getEquidistantCovering(numPoints, points,
				unbounded ? -Precision::UnboundedIntervalSamplingLimit : getLowerBound(),
				unbounded ? Precision::UnboundedIntervalSamplingLimit : getUpperBound());
		}

		/// @brief Sample inside a finite closed window, preserving holes and the caller output on failure.
		void getEquidistantCovering(int numPoints, std::vector<Real>& points, Real windowLower, Real windowUpper) const override {
			if (numPoints <= 0) { points.clear(); return; }
			IntervalDetail::ValidateSamplingWindow(windowLower, windowUpper);
			// Plain hulls are used only for grid allocation; every sample is checked against the original set.
			CompositeInterval grid;
			for (const auto& piece : _intervals) {
				const Real lo = std::max(piece->_lower, windowLower), hi = std::min(piece->_upper, windowUpper);
				const auto loType = lo > piece->_lower || piece->_lowerType == EndpointType::NEG_INF ? EndpointType::CLOSED : piece->_lowerType;
				const auto hiType = hi < piece->_upper || piece->_upperType == EndpointType::POS_INF ? EndpointType::CLOSED : piece->_upperType;
				grid.addPiece(lo, loType, hi, hiType);
			}
			grid.Normalize();
			std::vector<IntervalDetail::SamplingRange> ranges;
			for (const auto& piece : grid._intervals) {
				IntervalDetail::SamplingRange range;
				if (IntervalDetail::ClipSamplingRange(piece->_lower, piece->_lowerType, piece->_upper, piece->_upperType,
					windowLower, windowUpper, range))
					ranges.push_back(range);
			}

			std::vector<int> counts(ranges.size(), 0);
			int remaining = numPoints;
			if (ranges.size() <= static_cast<size_t>(remaining)) {
				for (auto& count : counts) count = 1;
				remaining -= static_cast<int>(ranges.size());
			}
			const bool halfWidths = std::any_of(ranges.begin(), ranges.end(), [](const auto& range) {
				return !std::isfinite(range.upper - range.lower);
			});
			std::vector<Real> weights;
			Real maxWeight = 0;
			for (const auto& range : ranges) {
				const Real width = halfWidths ? range.upper / 2 - range.lower / 2 : range.upper - range.lower;
				weights.push_back(width);
				maxWeight = std::max(maxWeight, width);
			}
			if (maxWeight > 0 && remaining > 0) {
				Real totalWeight = 0;
				for (auto& weight : weights) { weight /= maxWeight; totalWeight += weight; }
				std::vector<std::pair<Real, size_t>> remainders;
				const int budget = remaining;
				for (size_t i = 0; i < ranges.size(); ++i) {
					if (weights[i] == 0) continue;
					const Real quota = Real(budget) * (weights[i] / totalWeight);
					// Bound before conversion: even Real(INT_MAX) can round beyond int in float builds.
					const int count = quota >= Real(remaining) ? remaining : static_cast<int>(quota);
					counts[i] += count;
					remaining -= count;
					remainders.emplace_back(quota - Real(count), i);
				}
				std::stable_sort(remainders.begin(), remainders.end(), [](const auto& a, const auto& b) { return a.first > b.first; });
				size_t next = 0;
				while (remaining > 0) {
					++counts[remainders[next].second];
					--remaining;
					next = (next + 1) % remainders.size();
				}
			}
			else if (maxWeight == 0) {
				for (size_t i = 0; i < ranges.size() && remaining > 0; ++i)
					if (counts[i] == 0) { counts[i] = 1; --remaining; }
			}
			std::vector<Real> result;
			for (size_t i = 0; i < ranges.size(); ++i) {
				const auto samples = IntervalDetail::SampleRange(*this, counts[i], ranges[i]);
				result.insert(result.end(), samples.begin(), samples.end());
			}
			std::sort(result.begin(), result.end());
			result.erase(std::unique(result.begin(), result.end()), result.end());
			points = std::move(result);
		}
	};
	/// /** @} */
	// End Composite Interval

	/// @brief Back-compatibility alias. `CompositeInterval` is the composite (union-of-pieces) set.
	using Interval = CompositeInterval;
} // namespace MML

#endif
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
/// Interval (composite - stores multiple BaseInterval subintervals)
/// - Static methods: Intersection(), Difference(), Complement()
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
/// unitInterval.GetEquidistantCovering(10, points);
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

#include "MMLBase.h"

#include "interfaces/IInterval.h"

// Standard headers - include what we use
#include <algorithm>
#include <initializer_list>
#include <limits>
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
	/// /** @} */


	class Interval; // Forward declaration

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
		friend class Interval; // Allow Interval to access protected members

	protected:
		Real _lower, _upper;				 ///< Lower and upper bounds
		EndpointType _lowerType, _upperType; ///< Endpoint types (open/closed/infinite)

		BaseInterval(Real lower, EndpointType lowerType, Real upper, EndpointType upperType)
			: _lower(lower)
			, _lowerType(lowerType)
			, _upper(upper)
			, _upperType(upperType) {}

	public:
		virtual ~BaseInterval() {}

		/// @brief Gets the lower bound of the interval.
		Real getLowerBound() const { return _lower; }

		/// @brief Gets the upper bound of the interval.
		Real getUpperBound() const { return _upper; }

		/// @brief Gets the length of the interval (upper - lower).
		Real getLength() const { return _upper - _lower; }

		/// @brief Returns true if the interval is continuous (no holes).
		virtual bool isContinuous() const { return true; } // we suppose continuous intervals by default

		/// @brief Generates equidistant points covering the interval.
		/// Creates a uniform distribution of points across the interval, useful for:
		/// - Numerical integration sampling
		/// - Function plotting
		/// - Data point generation for curve fitting
		/// For open endpoints, points are slightly offset inward (by 0.1% of spacing)
		/// to ensure they remain strictly inside the interval.
		/// For infinite endpoints, practical limits (±10^10) are used.
		/// @param numPoints Number of points to generate (minimum 1).
		/// @param[out] points Vector to receive the generated points.
		/// @code{.cpp}
		/// ClosedInterval interval(0.0, 1.0);
		/// std::vector<Real> pts;
		/// interval.GetEquidistantCovering(5, pts);
		/// // pts = {0.0, 0.25, 0.5, 0.75, 1.0}
		/// @endcode

		void GetEquidistantCovering(int numPoints, std::vector<Real>& points) const {
			points.clear();
			if (numPoints <= 0)
				return;

			// Handle infinite endpoints
			Real lower = _lower;
			Real upper = _upper;

			// For infinite bounds, use practical limits (precision-dependent)
			if (_lowerType == EndpointType::NEG_INF)
				lower = std::is_same_v<Real, float> ? Real(-1e6) : Real(-1e10);
			if (_upperType == EndpointType::POS_INF)
				upper = std::is_same_v<Real, float> ? Real(1e6) : Real(1e10);

			if (numPoints == 1) {
				points.push_back((lower + upper) / 2.0);
				return;
			}

			Real delta = (upper - lower) / (numPoints - 1);

			for (int i = 0; i < numPoints; i++) {
				Real x = lower + i * delta;

				// Adjust for open endpoints
				if (i == 0 && _lowerType == EndpointType::OPEN)
					x += delta * 0.001; // Small epsilon inside
				if (i == numPoints - 1 && _upperType == EndpointType::OPEN)
					x -= delta * 0.001; // Small epsilon inside

				points.push_back(x);
			}
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
			: BaseInterval(-std::numeric_limits<double>::max(), EndpointType::NEG_INF, std::numeric_limits<double>::max(),
						   EndpointType::POS_INF) {}

		/// @brief Returns true for any real number x.
		bool contains(Real x) const { return true; }
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
			: BaseInterval(-std::numeric_limits<double>::max(), EndpointType::NEG_INF, std::numeric_limits<double>::max(),
						   EndpointType::POS_INF) {
			_hole0 = hole0;
			_holeDelta = holeDelta;
		}

		/// @brief Returns true if x is not at a hole position.
		bool contains(Real x) const {
			Real diff = (x - _hole0) / _holeDelta;
			if (std::abs(diff - std::round(diff)) < std::numeric_limits<Real>::epsilon() * 100)
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
		Real _lowerRealDif = REAL(0.0000001); ///< Small offset for internal calculations
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
		}

		/// @brief Returns true if x is in [lower, upper] and not at a hole.
		bool contains(Real x) const {
			if (x < _lower || x > _upper)
				return false;

			Real diff = (x - _hole0) / _holeDelta;
			if (std::abs(diff - std::round(diff)) < std::numeric_limits<Real>::epsilon() * 100)
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
			: BaseInterval(-std::numeric_limits<double>::max(), EndpointType::NEG_INF, upper, EndpointType::OPEN) {}

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
			: BaseInterval(-std::numeric_limits<double>::max(), EndpointType::NEG_INF, upper, EndpointType::CLOSED) {}

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
			: BaseInterval(lower, EndpointType::OPEN, std::numeric_limits<double>::max(), EndpointType::POS_INF) {}

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
			: BaseInterval(lower, EndpointType::CLOSED, std::numeric_limits<double>::max(), EndpointType::POS_INF) {}

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

	class Interval : public IInterval {
		Real _lower, _upper;
		std::vector<std::shared_ptr<BaseInterval>> _intervals; ///< Component intervals
	public:
		/// @brief Constructs an empty composite interval.
		Interval() {}

		/// @brief Constructs from a list of interval pointers.
		/// @param intervals Initializer list of BaseInterval pointers (takes ownership).

		Interval(std::initializer_list<BaseInterval*> intervals) {
			for (BaseInterval* interval : intervals) {
				_intervals.emplace_back(std::shared_ptr<BaseInterval>(interval));
			}
		}

		/// @brief Adds a subinterval to the composite.
		/// @tparam _IntervalType Type derived from BaseInterval.
		/// @param interval The interval to add.
		/// @return Reference to this for chaining.

		template<class _IntervalType>
		Interval& AddInterval(const _IntervalType& interval) {
			_intervals.emplace_back(std::make_shared<_IntervalType>(interval));
			return *this;
		}

		/// @brief Computes the intersection of two intervals.
		/// Returns A ∩ B, the set of points contained in both intervals.
		/// Handles endpoint types correctly (intersection of open and closed
		/// at the same point yields open).
		/// @param a First interval.
		/// @param b Second interval.
		/// @return Composite interval representing A ∩ B (may be empty).

		static Interval Intersection(const BaseInterval& a, const BaseInterval& b) {
			Interval ret;

			// Find intersection bounds
			Real lower = std::max(a._lower, b._lower);
			Real upper = std::min(a._upper, b._upper);

			// Check if intervals actually intersect
			if (lower > upper)
				return ret; // Empty intersection

			// Determine endpoint types
			EndpointType lowerType, upperType;

			if (lower == a._lower && lower == b._lower)
				lowerType = (a._lowerType == EndpointType::CLOSED && b._lowerType == EndpointType::CLOSED) ? EndpointType::CLOSED
																										   : EndpointType::OPEN;
			else if (lower == a._lower)
				lowerType = a._lowerType;
			else
				lowerType = b._lowerType;

			if (upper == a._upper && upper == b._upper)
				upperType = (a._upperType == EndpointType::CLOSED && b._upperType == EndpointType::CLOSED) ? EndpointType::CLOSED
																										   : EndpointType::OPEN;
			else if (upper == a._upper)
				upperType = a._upperType;
			else
				upperType = b._upperType;

			// Create appropriate interval based on endpoint types
			if (lowerType == EndpointType::CLOSED && upperType == EndpointType::CLOSED)
				ret.AddInterval(ClosedInterval(lower, upper));
			else if (lowerType == EndpointType::OPEN && upperType == EndpointType::CLOSED)
				ret.AddInterval(OpenClosedInterval(lower, upper));
			else if (lowerType == EndpointType::CLOSED && upperType == EndpointType::OPEN)
				ret.AddInterval(ClosedOpenInterval(lower, upper));
			else
				ret.AddInterval(OpenInterval(lower, upper));

			return ret;
		}

		/// @brief Computes the union of two intervals.
		/// Returns A ∪ B, the set of points contained in either interval.
		/// If the intervals overlap or are adjacent, returns a single merged interval.
		/// If disjoint, returns a composite interval with two parts.
		/// @param a First interval.
		/// @param b Second interval.
		/// @return Composite interval representing A ∪ B.

		static Interval Union(const BaseInterval& a, const BaseInterval& b) {
			Interval ret;

			// Check if intervals overlap or are adjacent
			// They overlap if a.upper >= b.lower AND b.upper >= a.lower
			// They are adjacent if one ends exactly where the other begins (considering endpoints)
			bool overlapsOrAdjacent = false;

			if (a._upper > b._lower && b._upper > a._lower) {
				// Strict overlap
				overlapsOrAdjacent = true;
			} else if (a._upper == b._lower) {
				// Adjacent: a ends where b begins
				// They merge if at least one endpoint is closed
				overlapsOrAdjacent = (a._upperType == EndpointType::CLOSED || b._lowerType == EndpointType::CLOSED);
			} else if (b._upper == a._lower) {
				// Adjacent: b ends where a begins
				overlapsOrAdjacent = (b._upperType == EndpointType::CLOSED || a._lowerType == EndpointType::CLOSED);
			}

			if (overlapsOrAdjacent) {
				// Merge into single interval
				Real lower = std::min(a._lower, b._lower);
				Real upper = std::max(a._upper, b._upper);

				// Determine endpoint types
				EndpointType lowerType, upperType;

				if (lower == a._lower && lower == b._lower)
					lowerType = (a._lowerType == EndpointType::CLOSED || b._lowerType == EndpointType::CLOSED) 
					            ? EndpointType::CLOSED : EndpointType::OPEN;
				else if (lower == a._lower)
					lowerType = a._lowerType;
				else
					lowerType = b._lowerType;

				if (upper == a._upper && upper == b._upper)
					upperType = (a._upperType == EndpointType::CLOSED || b._upperType == EndpointType::CLOSED) 
					            ? EndpointType::CLOSED : EndpointType::OPEN;
				else if (upper == a._upper)
					upperType = a._upperType;
				else
					upperType = b._upperType;

				// Create appropriate interval
				if (lowerType == EndpointType::CLOSED && upperType == EndpointType::CLOSED)
					ret.AddInterval(ClosedInterval(lower, upper));
				else if (lowerType == EndpointType::OPEN && upperType == EndpointType::CLOSED)
					ret.AddInterval(OpenClosedInterval(lower, upper));
				else if (lowerType == EndpointType::CLOSED && upperType == EndpointType::OPEN)
					ret.AddInterval(ClosedOpenInterval(lower, upper));
				else
					ret.AddInterval(OpenInterval(lower, upper));
			} else {
				// Disjoint intervals - add both (in order)
				const BaseInterval* first = (a._lower < b._lower) ? &a : &b;
				const BaseInterval* second = (a._lower < b._lower) ? &b : &a;

				// Add first interval
				if (first->_lowerType == EndpointType::CLOSED && first->_upperType == EndpointType::CLOSED)
					ret.AddInterval(ClosedInterval(first->_lower, first->_upper));
				else if (first->_lowerType == EndpointType::OPEN && first->_upperType == EndpointType::CLOSED)
					ret.AddInterval(OpenClosedInterval(first->_lower, first->_upper));
				else if (first->_lowerType == EndpointType::CLOSED && first->_upperType == EndpointType::OPEN)
					ret.AddInterval(ClosedOpenInterval(first->_lower, first->_upper));
				else
					ret.AddInterval(OpenInterval(first->_lower, first->_upper));

				// Add second interval
				if (second->_lowerType == EndpointType::CLOSED && second->_upperType == EndpointType::CLOSED)
					ret.AddInterval(ClosedInterval(second->_lower, second->_upper));
				else if (second->_lowerType == EndpointType::OPEN && second->_upperType == EndpointType::CLOSED)
					ret.AddInterval(OpenClosedInterval(second->_lower, second->_upper));
				else if (second->_lowerType == EndpointType::CLOSED && second->_upperType == EndpointType::OPEN)
					ret.AddInterval(ClosedOpenInterval(second->_lower, second->_upper));
				else
					ret.AddInterval(OpenInterval(second->_lower, second->_upper));
			}

			return ret;
		}

		/// @brief Computes the set difference of two intervals.
		/// Returns A \ B, the set of points in A but not in B.
		/// May result in zero, one, or two disjoint intervals depending on
		/// how B overlaps with A.
		/// @param a The interval to subtract from.
		/// @param b The interval to subtract.
		/// @return Composite interval representing A \ B.
		/// @note If B is entirely inside A, the result is two disjoint intervals.

		static Interval Difference(const BaseInterval& a, const BaseInterval& b) {
			Interval ret;

			// If b doesn't intersect a, return a
			if (b._upper <= a._lower || b._lower >= a._upper) {
				// Return copy of a based on its endpoint types
				if (a._lowerType == EndpointType::CLOSED && a._upperType == EndpointType::CLOSED)
					ret.AddInterval(ClosedInterval(a._lower, a._upper));
				else if (a._lowerType == EndpointType::OPEN && a._upperType == EndpointType::CLOSED)
					ret.AddInterval(OpenClosedInterval(a._lower, a._upper));
				else if (a._lowerType == EndpointType::CLOSED && a._upperType == EndpointType::OPEN)
					ret.AddInterval(ClosedOpenInterval(a._lower, a._upper));
				else
					ret.AddInterval(OpenInterval(a._lower, a._upper));
				return ret;
			}

			// Case: b completely contains a
			if (b._lower <= a._lower && b._upper >= a._upper)
				return ret; // Empty set

			// Case: b cuts left part of a
			if (b._lower <= a._lower && b._upper < a._upper) {
				EndpointType newLowerType = (b._upperType == EndpointType::CLOSED) ? EndpointType::OPEN : EndpointType::CLOSED;
				if (newLowerType == EndpointType::CLOSED && a._upperType == EndpointType::CLOSED)
					ret.AddInterval(ClosedInterval(b._upper, a._upper));
				else if (newLowerType == EndpointType::OPEN && a._upperType == EndpointType::CLOSED)
					ret.AddInterval(OpenClosedInterval(b._upper, a._upper));
				else if (newLowerType == EndpointType::CLOSED && a._upperType == EndpointType::OPEN)
					ret.AddInterval(ClosedOpenInterval(b._upper, a._upper));
				else
					ret.AddInterval(OpenInterval(b._upper, a._upper));
			}
			// Case: b cuts right part of a
			else if (b._lower > a._lower && b._upper >= a._upper) {
				EndpointType newUpperType = (b._lowerType == EndpointType::CLOSED) ? EndpointType::OPEN : EndpointType::CLOSED;
				if (a._lowerType == EndpointType::CLOSED && newUpperType == EndpointType::CLOSED)
					ret.AddInterval(ClosedInterval(a._lower, b._lower));
				else if (a._lowerType == EndpointType::OPEN && newUpperType == EndpointType::CLOSED)
					ret.AddInterval(OpenClosedInterval(a._lower, b._lower));
				else if (a._lowerType == EndpointType::CLOSED && newUpperType == EndpointType::OPEN)
					ret.AddInterval(ClosedOpenInterval(a._lower, b._lower));
				else
					ret.AddInterval(OpenInterval(a._lower, b._lower));
			}
			// Case: b is completely inside a (splits a into two intervals)
			else {
				// Left part: [a.lower, b.lower)
				EndpointType leftUpperType = (b._lowerType == EndpointType::CLOSED) ? EndpointType::OPEN : EndpointType::CLOSED;
				if (a._lowerType == EndpointType::CLOSED && leftUpperType == EndpointType::CLOSED)
					ret.AddInterval(ClosedInterval(a._lower, b._lower));
				else if (a._lowerType == EndpointType::OPEN && leftUpperType == EndpointType::CLOSED)
					ret.AddInterval(OpenClosedInterval(a._lower, b._lower));
				else if (a._lowerType == EndpointType::CLOSED && leftUpperType == EndpointType::OPEN)
					ret.AddInterval(ClosedOpenInterval(a._lower, b._lower));
				else
					ret.AddInterval(OpenInterval(a._lower, b._lower));

				// Right part: (b.upper, a.upper]
				EndpointType rightLowerType = (b._upperType == EndpointType::CLOSED) ? EndpointType::OPEN : EndpointType::CLOSED;
				if (rightLowerType == EndpointType::CLOSED && a._upperType == EndpointType::CLOSED)
					ret.AddInterval(ClosedInterval(b._upper, a._upper));
				else if (rightLowerType == EndpointType::OPEN && a._upperType == EndpointType::CLOSED)
					ret.AddInterval(OpenClosedInterval(b._upper, a._upper));
				else if (rightLowerType == EndpointType::CLOSED && a._upperType == EndpointType::OPEN)
					ret.AddInterval(ClosedOpenInterval(b._upper, a._upper));
				else
					ret.AddInterval(OpenInterval(b._upper, a._upper));
			}

			return ret;
		}

		/// @brief Computes the complement of an interval.
		/// Returns ℝ \ A, the set of all real numbers not in A.
		/// Results in zero, one, or two semi-infinite intervals.
		/// @param a The interval to complement.
		/// @return Composite interval representing ℝ \ A.
		/// @note Complement of (-∞, +∞) is the empty set.

		static Interval Complement(const BaseInterval& a) {
			Interval ret;

			// Handle special cases with infinity
			if (a._lowerType == EndpointType::NEG_INF && a._upperType == EndpointType::POS_INF)
				return ret; // Complement of R is empty

			// Add interval from -inf to lower bound
			if (a._lowerType != EndpointType::NEG_INF) {
				if (a._lowerType == EndpointType::CLOSED)
					ret.AddInterval(NegInfToOpenInterval(a._lower));
				else
					ret.AddInterval(NegInfToClosedInterval(a._lower));
			}

			// Add interval from upper bound to +inf
			if (a._upperType != EndpointType::POS_INF) {
				if (a._upperType == EndpointType::CLOSED)
					ret.AddInterval(OpenToInfInterval(a._upper));
				else
					ret.AddInterval(ClosedToInfInterval(a._upper));
			}

			return ret;
		}

		/// @brief Gets the minimum lower bound across all subintervals.
		Real getLowerBound() const {
			if (_intervals.empty())
				return 0;
			Real minLower = _intervals[0]->getLowerBound();
			for (const auto& interval : _intervals)
				minLower = std::min(minLower, interval->getLowerBound());
			return minLower;
		}

		/// @brief Gets the maximum upper bound across all subintervals.
		Real getUpperBound() const {
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
		Real getLength() const {
			return getUpperBound() - getLowerBound();
		}

		/// @brief Gets the total measure (sum of all subinterval lengths).
		/// @return Sum of individual subinterval lengths
		Real getMeasure() const {
			Real totalMeasure = 0;
			for (const auto& interval : _intervals)
				totalMeasure += interval->getLength();
			return totalMeasure;
		}

		/// @brief Returns false (composite intervals are not continuous by definition).
		bool isContinuous() const { return false; }

		/// @brief Checks if x is contained in any subinterval.
		/// @param x The point to test.
		/// @return True if x is in at least one subinterval.

		bool contains(Real x) const {
			// check for each interval if it contains x
			for (auto& interval : _intervals) {
				if (interval->contains(x))
					return true;
			}
			return false;
		}

		/// @brief Checks if another interval is fully contained.
		/// @param other The interval to test.
		/// @return True if other is entirely within one of the subintervals.
		/// @note This is a simplified check using only bounds, not full containment.

		bool contains(const BaseInterval& other) const {
			// Check if all points in 'other' are contained in at least one of our intervals
			// For simplicity, we check if other's bounds are contained
			for (const auto& interval : _intervals) {
				if (interval->contains(other.getLowerBound()) && interval->contains(other.getUpperBound()))
					return true;
			}
			return false;
		}

		/// @brief Checks if this composite intersects with another interval.
		/// @param other The interval to test for intersection.
		/// @return True if any subinterval overlaps with other.

		bool intersects(const BaseInterval& other) const {
			// Check if any of our intervals intersect with 'other'
			for (const auto& interval : _intervals) {
				// Two intervals intersect if one contains a point from the other
				if (interval->getLowerBound() <= other.getUpperBound() && interval->getUpperBound() >= other.getLowerBound())
					return true;
			}
			return false;
		}

		/// @brief Generates equidistant points across all subintervals.
		/// Points are distributed proportionally to subinterval lengths.
		/// @param numPoints Total number of points to generate.
		/// @param[out] points Vector to receive the generated points.

		void GetEquidistantCovering(int numPoints, std::vector<Real>& points) const {
			points.clear();
			if (_intervals.empty() || numPoints <= 0)
				return;

			// Distribute points across all intervals based on their relative lengths
			Real totalLength = getLength();
			if (totalLength <= 0)
				return;

			for (const auto& interval : _intervals) {
				int intervalPoints = std::max(1, (int)(numPoints * interval->getLength() / totalLength));
				std::vector<Real> intervalPointsVec;
				interval->GetEquidistantCovering(intervalPoints, intervalPointsVec);
				points.insert(points.end(), intervalPointsVec.begin(), intervalPointsVec.end());
			}
		}
	};
	/// /** @} */
	// End Composite Interval
} // namespace MML

#endif
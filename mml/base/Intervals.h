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

#include <mml/MMLBase.h>

#include <mml/interfaces/IInterval.h>

// Standard headers - include what we use
#include <algorithm>
#include <initializer_list>
#include <limits>
#include <memory>
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
			: BaseInterval(-std::numeric_limits<Real>::max(), EndpointType::NEG_INF, std::numeric_limits<Real>::max(),
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
			: BaseInterval(-std::numeric_limits<Real>::max(), EndpointType::NEG_INF, std::numeric_limits<Real>::max(),
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
		/// @note Periodic-hole intervals are treated as their continuous hull for set algebra.
		CompositeInterval(const BaseInterval& piece) {
			addPiece(piece._lower, piece._lowerType, piece._upper, piece._upperType);
		}

		/// @brief Adds a subinterval to the composite.
		/// @tparam _IntervalType Type derived from BaseInterval.
		/// @param interval The interval to add.
		/// @return Reference to this for chaining.

		template<class _IntervalType>
		CompositeInterval& AddInterval(const _IntervalType& interval) {
			_intervals.emplace_back(std::make_shared<_IntervalType>(interval));
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
			for (const auto& p : _intervals)
				segs.push_back({ p->_lower, p->_upper, p->_lowerType, p->_upperType });

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

			if (lo1 == lo2) {
				if (lo1T == EndpointType::NEG_INF && lo2T == EndpointType::NEG_INF) loT = EndpointType::NEG_INF;
				else loT = (lo1T == EndpointType::CLOSED && lo2T == EndpointType::CLOSED) ? EndpointType::CLOSED : EndpointType::OPEN;
			}
			else loT = (lo == lo1) ? lo1T : lo2T;

			if (hi1 == hi2) {
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
			CompositeInterval ret;
			ret._intervals = _intervals;
			ret._intervals.insert(ret._intervals.end(), other._intervals.begin(), other._intervals.end());
			ret.Normalize();
			return ret;
		}

		/// @brief Intersection with another interval set: this ∩ other (normalized).
		CompositeInterval intersectWith(const CompositeInterval& other) const {
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

	/// @brief Back-compatibility alias. `CompositeInterval` is the composite (union-of-pieces) set.
	using Interval = CompositeInterval;
} // namespace MML

#endif
///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        IInterval.h                                                         ///
///  Description: Interval interface for bounded ranges                               ///
///               Base class for open/closed interval representations                 ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @file IInterval.h
 * @brief Interface for mathematical intervals on the real line.
 * 
 * Defines the abstract interface for representing bounded ranges of real numbers.
 * Implementations can represent various interval types (open, closed, half-open)
 * and support operations like containment testing and point generation.
 */

#if !defined MML_IINTERVAL_H
#define MML_IINTERVAL_H

#include "MMLBase.h"
#include <vector>

namespace MML
{
	/**
	 * @brief Abstract interface for real-valued intervals.
	 * 
	 * This interface defines the contract for interval representations in MML.
	 * Intervals can be continuous or composed of multiple disjoint segments.
	 * 
	 * Implementations include:
	 * - OpenInterval: (a, b) - excludes endpoints
	 * - ClosedInterval: [a, b] - includes endpoints
	 * - HalfOpenInterval: [a, b) or (a, b]
	 * - UnionInterval: union of multiple intervals
	 * 
	 * @note Used extensively in numerical integration, function analysis,
	 *       and domain specification for solvers.
	 */
	class IInterval
	{
	public:
		/**
		 * @brief Get the lower bound of the interval.
		 * @return The infimum of the interval
		 */
		virtual Real getLowerBound() const = 0;
		
		/**
		 * @brief Get the upper bound of the interval.
		 * @return The supremum of the interval
		 */
		virtual Real getUpperBound() const = 0;
		
		/**
		 * @brief Calculate the length of the interval.
		 * 
		 * For continuous intervals, returns the measure: upperBound - lowerBound.
		 * 
		 * For non-continuous intervals (unions of disjoint segments), this returns
		 * the **hull length** (span from overall min to max), NOT the measure
		 * (sum of individual segment lengths). Use getMeasure() if you need the
		 * total length of all segments.
		 * 
		 * @return The hull length: getUpperBound() - getLowerBound()
		 * @see getMeasure() for total measure of non-continuous intervals
		 */
		virtual Real getLength() const = 0;

		/**
		 * @brief Calculate the measure (total length) of the interval.
		 * 
		 * For continuous intervals, this is identical to getLength().
		 * 
		 * For non-continuous intervals (unions of disjoint segments), this returns
		 * the sum of the lengths of all segments, which may be less than getLength().
		 * 
		 * Default implementation returns getLength() (suitable for continuous intervals).
		 * Override in union/compound interval implementations.
		 * 
		 * @return The total measure (sum of segment lengths)
		 */
		virtual Real getMeasure() const { return getLength(); }

		/**
		 * @brief Check if the interval is continuous (single connected segment).
		 * @return true if the interval has no gaps, false if it's a union of disjoint parts
		 */
		virtual bool isContinuous() const = 0;
		
		/**
		 * @brief Test if a point lies within the interval.
		 * @param x The point to test
		 * @return true if x is contained in the interval (respecting open/closed boundaries)
		 */
		virtual bool contains(Real x) const = 0;

		/**
		 * @brief Generate equidistant points covering the interval.
		 * 
		 * Creates a vector of numPoints equally-spaced points spanning
		 * from the lower to upper bound, suitable for sampling or plotting.
		 * 
		 * @param numPoints Number of points to generate
		 * @param[out] points Vector to store the generated points
		 */
		virtual void GetEquidistantCovering(int numPoints, std::vector<Real>& points) const = 0;

		virtual ~IInterval() {}
	};
}
#endif
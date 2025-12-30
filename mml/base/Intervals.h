///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Intervals.h                                                         ///
///  Description: Real and complex interval classes for interval arithmetic           ///
///               Bounds checking and range operations                                ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_INTERVALS_H
#define MML_INTERVALS_H

#include "MMLBase.h"

#include "interfaces/IInterval.h"

namespace MML
{
	// TODO - finalize intervals properly (and use them in test beds)
	///////////////////////////////////////////////   Interfaces    ///////////////////////////////////////////
	enum class EndpointType
	{
		OPEN,
		CLOSED,
		NEG_INF,
		POS_INF
	};

	class Interval;  // Forward declaration

	class BaseInterval : public IInterval
	{
		friend class Interval;  // Allow Interval to access protected members
		
	protected:
		Real _lower, _upper;
		EndpointType _lowerType, _upperType;

		BaseInterval(Real lower, EndpointType lowerType, Real upper, EndpointType upperType) : _lower(lower), _lowerType(lowerType), _upper(upper), _upperType(upperType) { }
	public:
		virtual ~BaseInterval() {}

		Real getLowerBound() const { return _lower; }
		Real getUpperBound() const { return _upper; }
		Real getLength()     const { return _upper - _lower; }

		virtual bool isContinuous()  const { return true; }     // we suppose continuous intervals by default

		void GetEquidistantCovering(int numPoints, std::vector<Real>& points) const
		{ 
			points.clear();
			if (numPoints <= 0) return;

			// Handle infinite endpoints
			Real lower = _lower;
			Real upper = _upper;

			// For infinite bounds, use practical limits
			if (_lowerType == EndpointType::NEG_INF)
				lower = -1e10;
			if (_upperType == EndpointType::POS_INF)
				upper = 1e10;

			if (numPoints == 1) {
				points.push_back((lower + upper) / 2.0);
				return;
			}

			Real delta = (upper - lower) / (numPoints - 1);
			
			for (int i = 0; i < numPoints; i++)
			{
				Real x = lower + i * delta;
				
				// Adjust for open endpoints
				if (i == 0 && _lowerType == EndpointType::OPEN)
					x += delta * 0.001;  // Small epsilon inside
				if (i == numPoints - 1 && _upperType == EndpointType::OPEN)
					x -= delta * 0.001;  // Small epsilon inside
					
				points.push_back(x);
			}
		}
	};

	class CompleteRInterval : public BaseInterval
	{
	public:
		CompleteRInterval() : BaseInterval(-std::numeric_limits<double>::max(), EndpointType::NEG_INF, std::numeric_limits<double>::max(), EndpointType::POS_INF) { }
		bool contains(Real x) const { return true; }
	};
	class CompleteRWithReccuringPointHoles : public BaseInterval
	{
		Real _hole0, _holeDelta;
	public:
		CompleteRWithReccuringPointHoles(Real hole0, Real holeDelta) : BaseInterval(-std::numeric_limits<double>::max(), EndpointType::NEG_INF, std::numeric_limits<double>::max(), EndpointType::POS_INF)
		{
			_hole0 = hole0;
			_holeDelta = holeDelta;
		}

		bool contains(Real x) const {


			if (x == _hole0)
				return false;

			Real diff = (x - _hole0) / _holeDelta;
			if (diff == (int)diff)
				return false;

			return true;
		}
		bool isContinuous() const { return false; }
	};
	class OpenInterval : public BaseInterval
	{
		// for equidistant covering
		Real _lowerRealDif = REAL(0.0000001);
	public:
		OpenInterval(Real lower, Real upper) : BaseInterval(lower, EndpointType::OPEN, upper, EndpointType::OPEN) { }
		bool contains(Real x) const { return (x > _lower) && (x < _upper); }
	};
	class OpenClosedInterval : public BaseInterval
	{
	public:
		OpenClosedInterval(Real lower, Real upper) : BaseInterval(lower, EndpointType::OPEN, upper, EndpointType::CLOSED) { }
		bool contains(Real x) const { return (x > _lower) && (x <= _upper); }
	};
	class ClosedInterval : public BaseInterval
	{
	public:
		ClosedInterval(Real lower, Real upper) : BaseInterval(lower, EndpointType::CLOSED, upper, EndpointType::CLOSED) { }

		bool contains(Real x) const {
			return (x >= _lower) && (x <= _upper);
		}
	};
	class ClosedIntervalWithReccuringPointHoles : public BaseInterval
	{
		Real _hole0, _holeDelta;
	public:
		ClosedIntervalWithReccuringPointHoles(Real lower, Real upper, Real hole0, Real holeDelta)
			: BaseInterval(lower, EndpointType::CLOSED, upper, EndpointType::CLOSED)
		{
			_hole0 = hole0;
			_holeDelta = holeDelta;
		}

		bool contains(Real x) const {
			// check for hole!
			if (x == _hole0)
				return false;

			Real diff = (x - _hole0) / _holeDelta;
			if (diff == (int)diff)
				return false;

			return (x >= _lower) && (x <= _upper);
		}
		bool isContinuous() const { return false; }
	};
	class ClosedOpenInterval : public BaseInterval
	{
	public:
		ClosedOpenInterval(Real lower, Real upper) : BaseInterval(lower, EndpointType::CLOSED, upper, EndpointType::OPEN) { }

		bool contains(Real x) const {
			return (x >= _lower) && (x < _upper);
		}
	};
	class NegInfToOpenInterval : public BaseInterval
	{
	public:
		NegInfToOpenInterval(Real upper) : BaseInterval(-std::numeric_limits<double>::max(), EndpointType::NEG_INF, upper, EndpointType::OPEN) { }

		bool contains(Real x) const {
			return x < _upper;
		}
	};
	class NegInfToClosedInterval : public BaseInterval
	{
	public:
		NegInfToClosedInterval(Real upper) : BaseInterval(-std::numeric_limits<double>::max(), EndpointType::NEG_INF, upper, EndpointType::CLOSED) { }

		bool contains(Real x) const {
			return x <= _upper;
		}
	};
	class OpenToInfInterval : public BaseInterval
	{
	public:
		OpenToInfInterval(Real lower) : BaseInterval(lower, EndpointType::OPEN, std::numeric_limits<double>::max(), EndpointType::POS_INF) { }

		bool contains(Real x) const {
			return x > _lower;
		}
	};
	class ClosedToInfInterval : public BaseInterval
	{
	public:
		ClosedToInfInterval(Real lower) : BaseInterval(lower, EndpointType::CLOSED, std::numeric_limits<double>::max(), EndpointType::POS_INF) { }

		bool contains(Real x) const {
			return x >= _lower;
		}
	};

	class Interval : public IInterval
	{
		Real _lower, _upper;
		std::vector<std::shared_ptr<BaseInterval>> _intervals;
	public:
		Interval() {}
		Interval(std::initializer_list<BaseInterval*> intervals)
		{
			for (BaseInterval* interval : intervals)
			{
				_intervals.emplace_back(std::shared_ptr<BaseInterval>(interval));
			}

			// sortirati po lower bound
			// i redom provjeriti presijecanja
		}

		template<class _IntervalType>
		Interval& AddInterval(const _IntervalType& interval)
		{
			_intervals.emplace_back(std::make_shared<_IntervalType>(interval));
			return *this;
		}

		static Interval Intersection(const BaseInterval& a, const BaseInterval& b)
		{
			Interval ret;
			
			// Find intersection bounds
			Real lower = std::max(a._lower, b._lower);
			Real upper = std::min(a._upper, b._upper);
			
			// Check if intervals actually intersect
			if (lower > upper)
				return ret;  // Empty intersection
				
			// Determine endpoint types
			EndpointType lowerType, upperType;
			
			if (lower == a._lower && lower == b._lower)
				lowerType = (a._lowerType == EndpointType::CLOSED && b._lowerType == EndpointType::CLOSED) ? EndpointType::CLOSED : EndpointType::OPEN;
			else if (lower == a._lower)
				lowerType = a._lowerType;
			else
				lowerType = b._lowerType;
				
			if (upper == a._upper && upper == b._upper)
				upperType = (a._upperType == EndpointType::CLOSED && b._upperType == EndpointType::CLOSED) ? EndpointType::CLOSED : EndpointType::OPEN;
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
		static Interval Difference(const BaseInterval& a, const BaseInterval& b)
		{
			Interval ret;
			
			// If b doesn't intersect a, return a
			if (b._upper <= a._lower || b._lower >= a._upper)
			{
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
				return ret;  // Empty set
				
			// Case: b cuts left part of a  
			if (b._lower <= a._lower && b._upper < a._upper)
			{
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
			else if (b._lower > a._lower && b._upper >= a._upper)
			{
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
			else
			{
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
		static Interval Complement(const BaseInterval& a)
		{
			Interval ret;
			
			// Handle special cases with infinity
			if (a._lowerType == EndpointType::NEG_INF && a._upperType == EndpointType::POS_INF)
				return ret;  // Complement of R is empty
				
			// Add interval from -inf to lower bound
			if (a._lowerType != EndpointType::NEG_INF)
			{
				if (a._lowerType == EndpointType::CLOSED)
					ret.AddInterval(NegInfToOpenInterval(a._lower));
				else
					ret.AddInterval(NegInfToClosedInterval(a._lower));
			}
			
			// Add interval from upper bound to +inf
			if (a._upperType != EndpointType::POS_INF)
			{
				if (a._upperType == EndpointType::CLOSED)
					ret.AddInterval(OpenToInfInterval(a._upper));
				else
					ret.AddInterval(ClosedToInfInterval(a._upper));
			}
			
			return ret;
		}
		Real getLowerBound() const 
		{ 
			if (_intervals.empty())
				return 0;
			Real minLower = _intervals[0]->getLowerBound();
			for (const auto& interval : _intervals)
				minLower = std::min(minLower, interval->getLowerBound());
			return minLower;
		}
		Real getUpperBound() const 
		{ 
			if (_intervals.empty())
				return 0;
			Real maxUpper = _intervals[0]->getUpperBound();
			for (const auto& interval : _intervals)
				maxUpper = std::max(maxUpper, interval->getUpperBound());
			return maxUpper;
		}
		Real getLength() const 
		{ 
			Real totalLength = 0;
			for (const auto& interval : _intervals)
				totalLength += interval->getLength();
			return totalLength;
		}

		bool isContinuous()  const { return false; }
		bool contains(Real x) const
		{
			// check for each interval if it contains x
			for (auto& interval : _intervals)
			{
				if (interval->contains(x))
					return true;
			}
			return false;
		}
		bool contains(const BaseInterval& other) const 
		{
			// Check if all points in 'other' are contained in at least one of our intervals
			// For simplicity, we check if other's bounds are contained
			for (const auto& interval : _intervals)
			{
				if (interval->contains(other.getLowerBound()) && 
				    interval->contains(other.getUpperBound()))
					return true;
			}
			return false;
		}
		
		bool intersects(const BaseInterval& other) const 
		{
			// Check if any of our intervals intersect with 'other'
			for (const auto& interval : _intervals)
			{
				// Two intervals intersect if one contains a point from the other
				if (interval->getLowerBound() <= other.getUpperBound() && 
				    interval->getUpperBound() >= other.getLowerBound())
					return true;
			}
			return false;
		}

		void GetEquidistantCovering(int numPoints, std::vector<Real>& points) const 
		{ 
			points.clear();
			if (_intervals.empty() || numPoints <= 0) return;
			
			// Distribute points across all intervals based on their relative lengths
			Real totalLength = getLength();
			if (totalLength <= 0) return;
			
			for (const auto& interval : _intervals)
			{
				int intervalPoints = std::max(1, (int)(numPoints * interval->getLength() / totalLength));
				std::vector<Real> intervalPointsVec;
				interval->GetEquidistantCovering(intervalPoints, intervalPointsVec);
				points.insert(points.end(), intervalPointsVec.begin(), intervalPointsVec.end());
			}
		}
	};
}

#endif
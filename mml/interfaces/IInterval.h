///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        IInterval.h                                                         ///
///  Description: Interval interface for bounded ranges                               ///
///               Base class for open/closed interval representations                 ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_IINTERVAL_H
#define MML_IINTERVAL_H

#include "MMLBase.h"
#include <vector>

namespace MML
{
	// group
	class IInterval
	{
	public:
		virtual Real getLowerBound() const = 0;
		virtual Real getUpperBound() const = 0;
		virtual Real getLength() const = 0;

		virtual bool isContinuous() const = 0;
		virtual bool contains(Real x) const = 0;

		virtual void GetEquidistantCovering(int numPoints, std::vector<Real>& points) const = 0;

		virtual ~IInterval() {}
	};
}
#endif
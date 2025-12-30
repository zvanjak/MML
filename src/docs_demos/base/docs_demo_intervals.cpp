#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Intervals.h"
#endif

using namespace MML;

void Docs_Demo_Intervals()
{
	std::cout << "***********************************************************" << std::endl;
	std::cout << "*****                     Intervals                  *******" << std::endl;

	// Simple interval types
	std::cout << "\n*** Simple Interval Types ***" << std::endl;
	
	ClosedInterval closed(0.0, 10.0);       // [0, 10]
	OpenInterval open(0.0, 10.0);           // (0, 10)
	OpenClosedInterval oc(0.0, 10.0);       // (0, 10]
	ClosedOpenInterval co(0.0, 10.0);       // [0, 10)
	
	std::cout << "ClosedInterval [0, 10]:" << std::endl;
	std::cout << "  contains(0) = " << (closed.contains(0.0) ? "true" : "false") << std::endl;
	std::cout << "  contains(5) = " << (closed.contains(5.0) ? "true" : "false") << std::endl;
	std::cout << "  contains(10) = " << (closed.contains(10.0) ? "true" : "false") << std::endl;
	
	std::cout << "\nOpenInterval (0, 10):" << std::endl;
	std::cout << "  contains(0) = " << (open.contains(0.0) ? "true" : "false") << std::endl;
	std::cout << "  contains(5) = " << (open.contains(5.0) ? "true" : "false") << std::endl;
	std::cout << "  contains(10) = " << (open.contains(10.0) ? "true" : "false") << std::endl;

	// Infinite bound intervals
	std::cout << "\n*** Infinite Bound Intervals ***" << std::endl;
	
	ClosedToInfInterval rightRay(5.0);      // [5, +inf)
	NegInfToClosedInterval leftRay(5.0);    // (-inf, 5]
	
	std::cout << "ClosedToInfInterval [5, +inf):" << std::endl;
	std::cout << "  contains(4) = " << (rightRay.contains(4.0) ? "true" : "false") << std::endl;
	std::cout << "  contains(5) = " << (rightRay.contains(5.0) ? "true" : "false") << std::endl;
	std::cout << "  contains(1000) = " << (rightRay.contains(1000.0) ? "true" : "false") << std::endl;

	// Interval properties
	std::cout << "\n*** Interval Properties ***" << std::endl;
	std::cout << "closed.getLowerBound() = " << closed.getLowerBound() << std::endl;
	std::cout << "closed.getUpperBound() = " << closed.getUpperBound() << std::endl;
	std::cout << "closed.getLength() = " << closed.getLength() << std::endl;

	// Equidistant covering
	std::cout << "\n*** Equidistant Covering ***" << std::endl;
	std::vector<Real> points;
	closed.GetEquidistantCovering(5, points);
	std::cout << "5 equidistant points in [0, 10]: ";
	for (const auto& p : points) {
		std::cout << p << " ";
	}
	std::cout << std::endl;

	// Interval set operations
	std::cout << "\n*** Interval Set Operations ***" << std::endl;
	
	ClosedInterval a(0.0, 10.0);
	ClosedInterval b(5.0, 15.0);
	
	Interval intersection = Interval::Intersection(a, b);
	std::cout << "Intersection of [0,10] and [5,15]: computed" << std::endl;
	
	ClosedInterval c(0.0, 10.0);
	ClosedInterval d(3.0, 7.0);
	Interval diff = Interval::Difference(c, d);
	std::cout << "Difference of [0,10] \\ [3,7]: computed (results in [0,3) U (7,10])" << std::endl;

	// Compound intervals
	std::cout << "\n*** Compound Intervals ***" << std::endl;
	Interval compound;
	compound.AddInterval(ClosedInterval(0.0, 5.0));
	compound.AddInterval(ClosedInterval(10.0, 15.0));
	std::cout << "Compound interval [0,5] U [10,15] created" << std::endl;
}
